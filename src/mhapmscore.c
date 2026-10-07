/*
 * mhapmscore: M-score statistics from an mHap file, per CpG or per region.
 *
 * Each mHap record (chr, start, end, haplotype, count[, strand]) is a read
 * covering consecutive CpGs of the CpG file. A read t seen c_t times
 * contributes N_t CpGs and Z_t = 1 if any of them is methylated, and
 *   reads = sum c_t, Nsum = sum c_t N_t, N2sum = sum c_t N_t^2,
 *   Sjd = sum c_t N_t Z_t, mscore = Sjd / Nsum, kappa = Nsum^2 / N2sum,
 *   Y_prime = Sjd * Nsum / N2sum
 * (the statistics of mhapasm and mHapDMR).
 *   - Per CpG (default): the reads whose haplotype spans the CpG; N_t and Z_t
 *     over all CpGs of the read, or with -w over those within INT bp of it.
 *   - Per region (-R): the reads with CpGs in the region; N_t and Z_t over
 *     those CpGs, as mHapDMR's .mscore_streaming.
 * As in mHapDMR and MscoreDMR, a record is used only if the CpG file has as
 * many CpGs in [start, end] as its haplotype has states.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <getopt.h>
#include <time.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>
#include "bsread.h"

#define MHAPMSCORE_VERSION "0.2.0"

KHASH_SET_INIT_STR(strset)

typedef struct {
    const char *in_fn, *cpg_fn, *out_fn, *region, *bed_fn, *regions_fn;
    int window, min_reads, threads, index;
} opt_t;

typedef struct {
    int64_t reads, nsum, n2sum, sjd;
} acc_t;

typedef struct {
    char *chr;
    hts_pos_t start, end;     /* 1-based positions of the first and last CpG */
    char *hap;
    int len;
    int64_t count;
} mrec_t;

typedef struct {              /* -r/-b: CpGs to report */
    int tid;                  /* contig in the mHap index */
    hts_pos_t beg, end;       /* 0-based, half-open */
} region_t;

typedef struct {              /* -R: one region of the BED file */
    char *chr;                /* contig as written in the BED file */
    int tid;                  /* contig in the mHap index, -1 if absent */
    hts_pos_t beg, end;       /* 0-based, half-open */
    acc_t a;
} rstat_t;

typedef struct {
    opt_t *o;
    htsFile *in;
    tbx_t *in_tbx;            /* mHap index, for -r/-b/-R */
    htsFile *cpg_fp;
    tbx_t *cpg_tbx;
    char *cpg_chr;            /* contig whose CpGs are loaded */
    int cpg_ok;               /* 0 if that contig is absent from the CpG file */
    hts_pos_t wbeg, wend;     /* the loaded CpGs are all those in [wbeg, wend) */
    cpgs_t cpg;
    acc_t *acc;               /* statistics of the CpGs [ia, ib) of the current block */
    int ia, ib;
    size_t macc;
    int *ones;                /* ones[k]: methylated states among the first k of a haplotype */
    int mones;
    BGZF *out;
    int out_sorted;           /* 0 if the -R output is not sorted by position */
    kstring_t ks, line;
    kstring_t no_cpg_chr;     /* contigs absent from the CpG file (first few) */
    int n_no_cpg_chr;
    /* run summary */
    int64_t n_records, n_bad, n_no_cpg, n_span, n_used, n_out, n_with_reads;
} ms_t;

/***************************************************************
 *                        input                                *
 ***************************************************************/

/* Split an mHap line in place; 0 for a well-formed record. */
static int parse_record(char *s, mrec_t *r)
{
    char *f[5], *e;
    int k = 1;
    f[0] = s;
    for (char *p = s; *p && k < 5; p++)
        if (*p == '\t') { *p = '\0'; f[k++] = p + 1; }
    if (k < 5) return -1;
    r->chr = f[0];
    r->start = strtoll(f[1], &e, 10);
    if (e == f[1] || *e) return -1;
    r->end = strtoll(f[2], &e, 10);
    if (e == f[2] || *e) return -1;
    r->count = strtoll(f[4], &e, 10);
    if (e == f[4] || (*e && *e != '\t' && *e != '\r')) return -1;
    r->hap = f[3];
    r->len = strlen(f[3]);
    if (!*r->chr || r->start < 1 || r->end < r->start || r->count < 1 || r->len == 0) return -1;
    for (int i = 0; i < r->len; i++)
        if (r->hap[i] != '0' && r->hap[i] != '1') return -1;
    return 0;
}

/* Split a BED line in place: 1 for a region, 0 for a comment or track line,
 * -1 if malformed (e.g. a column header). */
static int parse_bed(char *s, char **chr, hts_pos_t *beg, hts_pos_t *end)
{
    if (!*s || *s == '#' || !strncmp(s, "track", 5) || !strncmp(s, "browser", 7)) return 0;
    char *p = strchr(s, '\t'), *e;
    if (!p) return -1;
    *p = '\0';
    *chr = s;
    *beg = strtoll(p + 1, &e, 10);
    if (e == p + 1 || *e != '\t') return -1;
    *end = strtoll(e + 1, &p, 10);
    if (p == e + 1 || (*p && *p != '\t' && *p != '\r') || *beg < 0 || *end <= *beg) return -1;
    return 1;
}

#ifndef CPG_WINDOW
#define CPG_WINDOW 1000000    /* region queries load the CpGs of at least 1 Mb */
#endif
#ifndef CPG_PAD
#define CPG_PAD 10000         /* ... starting this far before the region */
#endif

static hts_pos_t add_pos(hts_pos_t a, hts_pos_t b) { return a > HTS_POS_MAX - b ? HTS_POS_MAX : a + b; }

/* Make sure the CpGs of chr in [beg, end) (0-based) are loaded; all CpGs of
 * the loaded window are kept, so a window loaded for one region serves its
 * neighbours. 0 if the contig is absent from the CpG file. */
static int load_cpgs(ms_t *c, const char *chr, hts_pos_t beg, hts_pos_t end)
{
    int same = c->cpg_chr && !strcmp(c->cpg_chr, chr);
    if (same && (!c->cpg_ok || (beg >= c->wbeg && end <= c->wend))) return c->cpg_ok;
    if (!same) {
        free(c->cpg_chr);
        c->cpg_chr = bs_strdup(chr);
    }
    if (beg < 0) beg = 0;
    if (end < add_pos(beg, CPG_WINDOW)) end = add_pos(beg, CPG_WINDOW);
    c->wbeg = beg;
    c->wend = end;
    c->cpg_ok = bs_cpgs_load(c->cpg_fp, c->cpg_tbx, chr, beg, end, &c->cpg, &c->ks) >= 0;
    if (!c->cpg_ok && c->n_no_cpg_chr++ < 5)
        ksprintf(&c->no_cpg_chr, "%s%s", c->no_cpg_chr.l ? ", " : "", chr);
    return c->cpg_ok;
}

/* Widen the loaded CpGs to the span of a record (region queries); 1 if they
 * were reloaded, which shifts all indexes into c->cpg. */
static int cover_record(ms_t *c, const mrec_t *r)
{
    if (r->start - 1 >= c->wbeg && r->end <= c->wend) return 0;
    load_cpgs(c, c->cpg_chr, r->start - 1 < c->wbeg ? r->start - 1 : c->wbeg, r->end > c->wend ? r->end : c->wend);
    return 1;
}

/* CpGs i0..i1 of a record; 0 if it does not match the CpG file. */
static int record_cpgs(ms_t *c, const mrec_t *r, int *i0, int *i1)
{
    *i0 = bs_lower_bound(c->cpg.pos, c->cpg.n, r->start - 1);  /* first CpG >= start */
    *i1 = bs_lower_bound(c->cpg.pos, c->cpg.n, r->end) - 1;    /* last CpG <= end */
    if (*i1 - *i0 + 1 == r->len) return 1;
    c->n_span++;
    return 0;
}

/***************************************************************
 *                        statistics                           *
 ***************************************************************/

static inline void acc_add(acc_t *a, int64_t cnt, int64_t n, int64_t z)
{
    a->reads += cnt;
    a->nsum += cnt * n;
    a->n2sum += cnt * n * n;
    a->sjd += cnt * n * z;
}

/* "\treads\tNsum\tN2sum\tSjd\tmscore\tkappa\tY_prime\n"; NA without reads */
static void put_stats(kstring_t *ks, const acc_t *a)
{
    ksprintf(ks, "\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t%" PRId64, a->reads, a->nsum, a->n2sum, a->sjd);
    if (a->nsum > 0) {        /* every read has N_t >= 1, so N2sum > 0 too */
        double nsum = a->nsum, n2sum = a->n2sum, sjd = a->sjd;
        ksprintf(ks, "\t%.6g\t%.6g\t%.6g\n", sjd / nsum, nsum * nsum / n2sum, sjd * nsum / n2sum);
    } else {
        kputs("\tNA\tNA\tNA\n", ks);
    }
}

static void write_line(ms_t *c)
{
    if (bgzf_write(c->out, c->line.s, c->line.l) < 0) {
        fprintf(stderr, "[mhapmscore] error writing %s\n", c->o->out_fn);
        exit(1);
    }
    c->n_out++;
}

/***************************************************************
 *                        per CpG                              *
 ***************************************************************/

/* Fresh statistics for the CpGs in [beg, end) of the loaded contig. */
static void block_start(ms_t *c, hts_pos_t beg, hts_pos_t end)
{
    c->ia = bs_lower_bound(c->cpg.pos, c->cpg.n, beg);
    c->ib = bs_lower_bound(c->cpg.pos, c->cpg.n, end);
    size_t n = c->ib - c->ia;
    if (n > c->macc) {
        free(c->acc);
        c->acc = bs_malloc(n * sizeof(acc_t));
        c->macc = n;
    }
    if (n) memset(c->acc, 0, n * sizeof(acc_t));
}

/* Add a read to the CpGs of the current block that it spans. */
static void add_record(ms_t *c, const mrec_t *r)
{
    const hts_pos_t *pos = c->cpg.pos;
    int i0, i1;
    if (!record_cpgs(c, r, &i0, &i1)) return;
    c->n_used++;
    int lo_j = i0 > c->ia ? i0 : c->ia;
    int hi_j = i1 < c->ib - 1 ? i1 : c->ib - 1;
    if (c->o->window == 0) {                 /* all CpGs of the read */
        int64_t z = memchr(r->hap, '1', r->len) != NULL;
        for (int j = lo_j; j <= hi_j; j++) acc_add(&c->acc[j - c->ia], r->count, r->len, z);
        return;
    }
    if (r->len + 1 > c->mones) {
        c->mones = r->len + 1;
        c->ones = bs_realloc(c->ones, c->mones * sizeof(int));
    }
    c->ones[0] = 0;
    for (int k = 0; k < r->len; k++) c->ones[k + 1] = c->ones[k] + (r->hap[k] == '1');
    /* CpGs lo..hi of the read lie within the window around CpG j; both ends
     * only move forward as j does */
    hts_pos_t w = c->o->window;
    int lo = i0, hi = i0;
    for (int j = lo_j; j <= hi_j; j++) {
        while (pos[lo] < pos[j] - w) lo++;
        if (hi < j) hi = j;
        while (hi < i1 && pos[hi + 1] <= pos[j] + w) hi++;
        acc_add(&c->acc[j - c->ia], r->count, hi - lo + 1, c->ones[hi - i0 + 1] > c->ones[lo - i0]);
    }
}

static void block_write(ms_t *c, const char *chr)
{
    for (int j = c->ia; j < c->ib; j++) {
        const acc_t *a = &c->acc[j - c->ia];
        if (a->reads < c->o->min_reads) continue;
        c->line.l = 0;
        ksprintf(&c->line, "%s\t%" PRIhts_pos, chr, c->cpg.pos[j] + 1);
        put_stats(&c->line, a);
        write_line(c);
    }
}

/* Whole file: the records of a contig must be together, in any order. */
static int run_file(ms_t *c)
{
    kstring_t s = {0, 0, NULL};
    khash_t(strset) *seen = kh_init(strset);
    char *chr = NULL;
    int ok = 0, ret, rc = 0;
    while ((ret = hts_getline(c->in, '\n', &s)) >= 0) {
        if (!s.l || s.s[0] == '#') continue;
        c->n_records++;
        mrec_t r;
        if (parse_record(s.s, &r) < 0) { c->n_bad++; continue; }
        if (!chr || strcmp(chr, r.chr)) {
            if (chr && ok) block_write(c, chr);
            int absent;
            khint_t k = kh_put(strset, seen, r.chr, &absent);
            if (!absent) {
                fprintf(stderr, "[mhapmscore] the records of %s are not together in %s (sort it by position)\n",
                        r.chr, c->o->in_fn);
                rc = -1;
                break;
            }
            kh_key(seen, k) = bs_strdup(r.chr);
            free(chr);
            chr = bs_strdup(r.chr);
            ok = load_cpgs(c, chr, 0, HTS_POS_MAX);
            if (ok) block_start(c, 0, HTS_POS_MAX);
        }
        if (!ok) { c->n_no_cpg++; continue; }
        add_record(c, &r);
    }
    if (rc == 0 && ret < -1) {
        fprintf(stderr, "[mhapmscore] error reading %s\n", c->o->in_fn);
        rc = -1;
    }
    if (rc == 0 && chr && ok) block_write(c, chr);
    for (khint_t k = kh_begin(seen); k != kh_end(seen); k++)
        if (kh_exist(seen, k)) free((char *) kh_key(seen, k));
    kh_destroy(strset, seen);
    free(chr);
    free(s.s);
    return rc;
}

static int region_cmp(const void *a, const void *b)
{
    const region_t *x = a, *y = b;
    if (x->tid != y->tid) return x->tid < y->tid ? -1 : 1;
    if (x->beg != y->beg) return x->beg < y->beg ? -1 : 1;
    return 0;
}

/* Regions from -r or -b on contigs of the mHap index, sorted and merged.
 * Returns their number, or -1 on error. */
static int read_regions(ms_t *c, region_t **out)
{
    const opt_t *o = c->o;
    region_t *reg = NULL;
    int n = 0, m = 0, n_skip = 0, err = 0;
    kstring_t s = {0, 0, NULL};
    htsFile *fp = NULL;
    if (o->bed_fn && !(fp = hts_open(o->bed_fn, "r"))) {
        fprintf(stderr, "[mhapmscore] cannot open %s\n", o->bed_fn);
        return -1;
    }
    for (int done = 0; !done; ) {
        char *name;
        hts_pos_t beg, end;
        if (fp) {
            if (hts_getline(fp, '\n', &s) < 0) break;
            int k = parse_bed(s.s, &name, &beg, &end);
            if (k == 0) continue;
            if (k < 0) { n_skip++; continue; }
        } else {
            const char *q = hts_parse_reg64(o->region, &beg, &end);
            if (!q || end <= beg) { fprintf(stderr, "[mhapmscore] invalid region %s\n", o->region); err = 1; break; }
            s.l = 0;
            kputsn(o->region, q - o->region, &s);
            name = s.s;
            done = 1;
        }
        int tid = bs_tbx_tid(c->in_tbx, name);
        if (tid < 0) { n_skip++; continue; }
        if (n == m) { m = m ? m * 2 : 64; reg = bs_realloc(reg, m * sizeof(region_t)); }
        reg[n].tid = tid;
        reg[n].beg = beg;
        reg[n++].end = end;
    }
    if (fp) hts_close(fp);
    free(s.s);
    if (err) { free(reg); return -1; }
    if (n_skip) fprintf(stderr, "[mhapmscore] skipped %d region(s): malformed or contig absent from the mHap file\n", n_skip);
    qsort(reg, n, sizeof(region_t), region_cmp);
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (k && reg[i].tid == reg[k - 1].tid && reg[i].beg <= reg[k - 1].end) {
            if (reg[i].end > reg[k - 1].end) reg[k - 1].end = reg[i].end;
        } else reg[k++] = reg[i];
    }
    *out = reg;
    return k;
}

/* -r/-b: only the CpGs inside the regions are reported, from all reads
 * spanning them (the tabix query also returns reads that start before). */
static int run_regions(ms_t *c)
{
    int nseq = 0, rc = 0;
    region_t *reg = NULL;
    int nreg = read_regions(c, &reg);
    if (nreg < 0) return -1;
    const char **names = tbx_seqnames(c->in_tbx, &nseq);
    kstring_t s = {0, 0, NULL};
    for (int i = 0; i < nreg && rc == 0; i++) {
        const char *chr = names[reg[i].tid];
        if (!load_cpgs(c, chr, reg[i].beg - CPG_PAD, add_pos(reg[i].end, CPG_PAD))) continue;
        block_start(c, reg[i].beg, reg[i].end);
        if (c->ib <= c->ia) continue;
        /* one base of slack for indexes made with another coordinate convention */
        hts_itr_t *itr = tbx_itr_queryi(c->in_tbx, reg[i].tid, reg[i].beg > 0 ? reg[i].beg - 1 : 0,
                                        reg[i].end < HTS_POS_MAX ? reg[i].end + 1 : reg[i].end);
        if (!itr) continue;
        int ret;
        while ((ret = tbx_itr_next(c->in, c->in_tbx, itr, &s)) >= 0) {
            c->n_records++;
            mrec_t r;
            if (parse_record(s.s, &r) < 0) { c->n_bad++; continue; }
            if (cover_record(c, &r)) {       /* the same CpGs at new indexes */
                c->ia = bs_lower_bound(c->cpg.pos, c->cpg.n, reg[i].beg);
                c->ib = bs_lower_bound(c->cpg.pos, c->cpg.n, reg[i].end);
            }
            add_record(c, &r);
        }
        if (ret < -1) { fprintf(stderr, "[mhapmscore] error reading %s\n", c->o->in_fn); rc = -1; }
        tbx_itr_destroy(itr);
        if (rc == 0) block_write(c, chr);
    }
    free(names);
    free(reg);
    free(s.s);
    return rc;
}

/***************************************************************
 *                        per region (-R)                      *
 ***************************************************************/

static int rstat_cmp(const void *a, const void *b)
{
    const rstat_t *x = *(rstat_t *const *) a, *y = *(rstat_t *const *) b;
    if (x->tid != y->tid) return x->tid < y->tid ? -1 : 1;
    if (x->beg != y->beg) return x->beg < y->beg ? -1 : 1;
    return 0;
}

/* The reads with CpGs in [beg, end) of a region; N_t and Z_t over those CpGs. */
static int region_stats(ms_t *c, rstat_t *g, kstring_t *s)
{
    int rlo = bs_lower_bound(c->cpg.pos, c->cpg.n, g->beg);
    int rhi = bs_lower_bound(c->cpg.pos, c->cpg.n, g->end) - 1;
    if (rlo > rhi) return 0;                 /* no CpG in the region */
    hts_itr_t *itr = tbx_itr_queryi(c->in_tbx, g->tid, g->beg > 0 ? g->beg - 1 : 0, g->end + 1);
    if (!itr) return 0;
    int ret;
    while ((ret = tbx_itr_next(c->in, c->in_tbx, itr, s)) >= 0) {
        c->n_records++;
        mrec_t r;
        int i0, i1;
        if (parse_record(s->s, &r) < 0) { c->n_bad++; continue; }
        if (cover_record(c, &r)) {
            rlo = bs_lower_bound(c->cpg.pos, c->cpg.n, g->beg);
            rhi = bs_lower_bound(c->cpg.pos, c->cpg.n, g->end) - 1;
        }
        if (!record_cpgs(c, &r, &i0, &i1)) continue;
        int lo = i0 > rlo ? i0 : rlo, hi = i1 < rhi ? i1 : rhi;
        if (lo > hi) continue;               /* no CpG of the read in the region */
        c->n_used++;
        acc_add(&g->a, r.count, hi - lo + 1, memchr(r.hap + (lo - i0), '1', hi - lo + 1) != NULL);
    }
    tbx_itr_destroy(itr);
    if (ret < -1) { fprintf(stderr, "[mhapmscore] error reading %s\n", c->o->in_fn); return -1; }
    return 0;
}

/* -R: one line per BED line, in input order. */
static int run_region_stats(ms_t *c)
{
    htsFile *fp = hts_open(c->o->regions_fn, "r");
    if (!fp) { fprintf(stderr, "[mhapmscore] cannot open %s\n", c->o->regions_fn); return -1; }
    kstring_t s = {0, 0, NULL};
    rstat_t *reg = NULL;
    int n = 0, m = 0, n_skip = 0, rc = 0;
    while (hts_getline(fp, '\n', &s) >= 0) {
        char *name;
        hts_pos_t beg, end;
        int k = parse_bed(s.s, &name, &beg, &end);
        if (k == 0) continue;
        if (k < 0) { n_skip++; continue; }
        if (n == m) { m = m ? m * 2 : 1024; reg = bs_realloc(reg, m * sizeof(rstat_t)); }
        memset(&reg[n], 0, sizeof(rstat_t));
        reg[n].chr = bs_strdup(name);
        reg[n].tid = bs_tbx_tid(c->in_tbx, name);
        reg[n].beg = beg;
        reg[n++].end = end;
    }
    hts_close(fp);
    if (n_skip) fprintf(stderr, "[mhapmscore] skipped %d malformed line(s) of %s\n", n_skip, c->o->regions_fn);

    /* processed by position, so that the CpGs of each contig are loaded once */
    rstat_t **ord = bs_malloc((n ? n : 1) * sizeof(rstat_t *));
    for (int i = 0; i < n; i++) ord[i] = &reg[i];
    qsort(ord, n, sizeof(rstat_t *), rstat_cmp);
    int nseq = 0;
    const char **names = tbx_seqnames(c->in_tbx, &nseq);
    for (int i = 0; i < n && rc == 0; i++) {
        rstat_t *g = ord[i];
        if (g->tid >= 0 && load_cpgs(c, names[g->tid], g->beg - CPG_PAD, add_pos(g->end, CPG_PAD)))
            rc = region_stats(c, g, &s);
    }
    free(names);
    free(ord);

    /* written in input order; the tabix index needs it sorted by position */
    khash_t(strset) *done = kh_init(strset);
    const char *prev = NULL;
    hts_pos_t prev_beg = 0;
    c->out_sorted = 1;
    for (int i = 0; i < n && rc == 0; i++) {
        rstat_t *g = &reg[i];
        if (!prev || strcmp(prev, g->chr)) {
            int absent;
            if (prev) kh_put(strset, done, prev, &absent);
            if (kh_get(strset, done, g->chr) != kh_end(done)) c->out_sorted = 0;
            prev = g->chr;
        } else if (g->beg < prev_beg) {
            c->out_sorted = 0;
        }
        prev_beg = g->beg;
        if (g->a.reads) c->n_with_reads++;
        c->line.l = 0;
        ksprintf(&c->line, "%s\t%" PRIhts_pos "\t%" PRIhts_pos, g->chr, g->beg, g->end);
        put_stats(&c->line, &g->a);
        write_line(c);
    }
    kh_destroy(strset, done);
    for (int i = 0; i < n; i++) free(reg[i].chr);
    free(reg);
    free(s.s);
    return rc;
}

/***************************************************************
 *                        main                                 *
 ***************************************************************/

static void usage(FILE *fp)
{
    fprintf(fp,
"mhapmscore " MHAPMSCORE_VERSION ": M-score statistics from an mHap file, per CpG or per region\n"
"\n"
"Usage: mhapmscore -i <in.mhap.gz> -c <CpG.gz> [-r chr:beg-end | -b regions.bed] [-o out.tsv.gz] [options]\n"
"       mhapmscore -i <in.mhap.gz> -c <CpG.gz> -R regions.bed [-o out.tsv] [options]\n"
"\n"
"Input:\n"
"  -i, --input FILE        mHap file (chr, start, end, haplotype, count[, strand]), plain or\n"
"                          bgzipped, the records of a contig together; -r/-b/-R need a tabix index\n"
"  -c, --cpg FILE          CpG positions (chr, 1-based C position, ...), bgzipped and tabix-indexed\n"
"Per CpG (default):\n"
"  -r, --region STR        only CpGs in this region (chr:beg-end)\n"
"  -b, --bed FILE          only CpGs in these BED regions\n"
"  -w, --window INT        N_t and Z_t from the CpGs of the read within INT bp of the CpG;\n"
"                          0 = all CpGs of the read [0]\n"
"  -m, --min-reads INT     report CpGs spanned by at least INT reads [1]\n"
"Per region:\n"
"  -R, --regions FILE      statistics of each region of this BED file instead of each CpG\n"
"Output:\n"
"  -o, --output FILE       [stdout]; bgzipped and tabix-indexed if it ends in .gz\n"
"      --no-index          do not write the tabix index\n"
"  -@, --threads INT       extra threads for decompression and compression [0]\n"
"  -h, --help              show this help\n"
"  -v, --version           show the version\n"
"\n"
"Output columns:\n"
"  per CpG:    chr pos reads Nsum N2sum Sjd mscore kappa Y_prime   (pos 1-based)\n"
"  per region: chr start end reads Nsum N2sum Sjd mscore kappa Y_prime   (BED coordinates,\n"
"              one line per BED line in input order; NA without reads)\n"
"  A read t, seen c_t times, contributes N_t CpGs and Z_t = 1 if any of them is methylated:\n"
"  reads = sum c_t, Nsum = sum c_t*N_t, N2sum = sum c_t*N_t^2, Sjd = sum c_t*N_t*Z_t,\n"
"  mscore = Sjd/Nsum, kappa = Nsum^2/N2sum, Y_prime = Sjd*Nsum/N2sum (as mhapasm and mHapDMR).\n"
"  Per CpG, the reads spanning the CpG count, N_t over all their CpGs (or within -w bp);\n"
"  per region, the reads with CpGs in the region count, N_t over those CpGs.\n");
}

int main(int argc, char **argv)
{
    bs_prog = "mhapmscore";
    opt_t o;
    memset(&o, 0, sizeof(o));
    o.out_fn = "-";
    o.min_reads = 1;
    o.index = 1;
    static const struct option lopts[] = {
        {"input", required_argument, NULL, 'i'},
        {"cpg", required_argument, NULL, 'c'},
        {"region", required_argument, NULL, 'r'},
        {"bed", required_argument, NULL, 'b'},
        {"regions", required_argument, NULL, 'R'},
        {"window", required_argument, NULL, 'w'},
        {"min-reads", required_argument, NULL, 'm'},
        {"output", required_argument, NULL, 'o'},
        {"no-index", no_argument, NULL, 1},
        {"threads", required_argument, NULL, '@'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}
    };
    int ch;
    while ((ch = getopt_long(argc, argv, "i:c:r:b:R:w:m:o:@:hv", lopts, NULL)) >= 0) {
        switch (ch) {
        case 'i': o.in_fn = optarg; break;
        case 'c': o.cpg_fn = optarg; break;
        case 'r': o.region = optarg; break;
        case 'b': o.bed_fn = optarg; break;
        case 'R': o.regions_fn = optarg; break;
        case 'w': o.window = atoi(optarg); break;
        case 'm': o.min_reads = atoi(optarg); break;
        case 'o': o.out_fn = optarg; break;
        case 1: o.index = 0; break;
        case '@': o.threads = atoi(optarg); break;
        case 'h': usage(stdout); return 0;
        case 'v': puts(MHAPMSCORE_VERSION); return 0;
        default: usage(stderr); return 1;
        }
    }
    if (!o.in_fn || !o.cpg_fn) { usage(stderr); return 1; }
    if (!!o.region + !!o.bed_fn + !!o.regions_fn > 1) { fprintf(stderr, "[mhapmscore] use only one of -r, -b and -R\n"); return 1; }
    if (o.regions_fn && o.window) { fprintf(stderr, "[mhapmscore] -w is for per-CpG statistics; -R uses the CpGs in each region\n"); return 1; }
    if (o.window < 0 || o.min_reads < 1) { fprintf(stderr, "[mhapmscore] -w must be >= 0 and -m >= 1\n"); return 1; }

    ms_t c;
    memset(&c, 0, sizeof(c));
    c.o = &o;
    clock_t t0 = clock();
    if (!(c.in = hts_open(o.in_fn, "r"))) { fprintf(stderr, "[mhapmscore] cannot open %s\n", o.in_fn); return 1; }
    if (o.threads > 0 && hts_get_format(c.in)->compression == bgzf) hts_set_threads(c.in, o.threads);
    if (o.region || o.bed_fn || o.regions_fn) {
        if (!(c.in_tbx = tbx_index_load(o.in_fn))) {
            fprintf(stderr, "[mhapmscore] -r/-b/-R need a tabix index of %s (tabix -b 2 -e 3)\n", o.in_fn);
            return 1;
        }
        hts_set_cache_size(c.in, 16 << 20);  /* neighbouring regions often share blocks */
    }
    if (!(c.cpg_fp = hts_open(o.cpg_fn, "r"))) { fprintf(stderr, "[mhapmscore] cannot open %s\n", o.cpg_fn); return 1; }
    if (!(c.cpg_tbx = tbx_index_load(o.cpg_fn))) {
        fprintf(stderr, "[mhapmscore] %s needs a tabix index (tabix -s 1 -b 2 -e 2)\n", o.cpg_fn);
        return 1;
    }
    int gz = strlen(o.out_fn) > 3 && !strcmp(o.out_fn + strlen(o.out_fn) - 3, ".gz");
    if (!(c.out = bgzf_open(o.out_fn, gz ? "w" : "wu"))) { fprintf(stderr, "[mhapmscore] cannot write %s\n", o.out_fn); return 1; }
    if (gz && o.threads > 0) bgzf_mt(c.out, o.threads, 256);
    const char *header = o.regions_fn ? "#chr\tstart\tend\treads\tNsum\tN2sum\tSjd\tmscore\tkappa\tY_prime\n"
                                      : "#chr\tpos\treads\tNsum\tN2sum\tSjd\tmscore\tkappa\tY_prime\n";
    if (bgzf_write(c.out, header, strlen(header)) < 0) { fprintf(stderr, "[mhapmscore] error writing %s\n", o.out_fn); return 1; }

    int rc = o.regions_fn ? run_region_stats(&c) : c.in_tbx ? run_regions(&c) : run_file(&c);

    if (bgzf_close(c.out) < 0) { fprintf(stderr, "[mhapmscore] error closing %s\n", o.out_fn); rc = -1; }
    if (rc == 0 && gz && o.index && strcmp(o.out_fn, "-")) {
        if (c.n_out == 0) {
            fprintf(stderr, "[mhapmscore] no output lines; index not written\n");
        } else if (o.regions_fn && !c.out_sorted) {
            fprintf(stderr, "[mhapmscore] the regions are not sorted by position; index not written\n");
        } else {
            tbx_conf_t conf_cpg = {TBX_GENERIC, 1, 2, 2, '#', 0}, conf_bed = {TBX_UCSC, 1, 2, 3, '#', 0};
            if (tbx_index_build(o.out_fn, 0, o.regions_fn ? &conf_bed : &conf_cpg) < 0) {
                fprintf(stderr, "[mhapmscore] failed to index %s\n", o.out_fn);
                rc = -1;
            }
        }
    }
    if (c.n_no_cpg_chr)
        fprintf(stderr, "[mhapmscore] %d contig(s) are absent from the CpG file, their records were skipped: %s%s\n",
                c.n_no_cpg_chr, c.no_cpg_chr.s, c.n_no_cpg_chr > 5 ? ", ..." : "");
    const char *out_name = strcmp(o.out_fn, "-") ? o.out_fn : "stdout";
    if (o.regions_fn)
        fprintf(stderr, "[mhapmscore] %" PRId64 " regions (%" PRId64 " with reads) written to %s; %" PRId64
                " records read, %" PRId64 " used (%" PRId64 " not matching the CpG file, %" PRId64 " malformed), %.2f s CPU\n",
                c.n_out, c.n_with_reads, out_name, c.n_records, c.n_used, c.n_span, c.n_bad,
                (double) (clock() - t0) / CLOCKS_PER_SEC);
    else
        fprintf(stderr, "[mhapmscore] %" PRId64 " records, %" PRId64 " used (%" PRId64 " not matching the CpG file, %" PRId64
                " malformed, %" PRId64 " on contigs absent from the CpG file), %" PRId64 " CpGs written to %s, %.2f s CPU\n",
                c.n_records, c.n_used, c.n_span, c.n_bad, c.n_no_cpg, c.n_out, out_name,
                (double) (clock() - t0) / CLOCKS_PER_SEC);
    if (c.n_span > 0.01 * (c.n_records - c.n_bad - c.n_no_cpg))
        fprintf(stderr, "[mhapmscore] warning: over 1%% of the records do not match the CpG file;"
                " is it the CpG file (genome build) used to make the mHap file?\n");

    hts_close(c.in);
    if (c.in_tbx) tbx_destroy(c.in_tbx);
    hts_close(c.cpg_fp);
    tbx_destroy(c.cpg_tbx);
    bs_cpgs_free(&c.cpg);
    free(c.cpg_chr);
    free(c.acc);
    free(c.ones);
    free(c.ks.s);
    free(c.line.s);
    free(c.no_cpg_chr.s);
    return rc == 0 ? 0 : 1;
}
