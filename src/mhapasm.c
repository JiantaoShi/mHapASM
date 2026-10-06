/*
 * mhapasm: allele-specific M-score statistics from bisulfite sequencing BAMs.
 *
 * For every SNP in a list, reads are paired into fragments, each fragment is
 * assigned to the REF or ALT allele, CpG methylation is called on the fragment,
 * and the per-allele sufficient statistics of the M-score are reported
 * (reads, Nsum, N2sum, Sjd, mscore, kappa, Y_prime; same definitions as
 * mHapDMR, and Sjd/Nsum/N2sum are the S/T/N_sq of MscoreDMR).
 *
 * The read-level rules follow wgbs_tools (nloyfer/wgbs_tools):
 *   - allele call:      snp_patter::compareSeqToRef + snp_patter::proc2lines
 *   - methylation call: patter::compareSeqToRef (CIGAR-cleaned read, CpG
 *                       context checked on the read itself)
 *   - mate merging:     patter_utils::merge_PE (disagreeing CpGs dropped)
 *   - strand:           patter_utils::is_bottom, with aligner strand tags
 *                       (XG, YD, ZS) taking precedence when present
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <getopt.h>
#include <time.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>

#define MHAPASM_VERSION "0.1.0"

KHASH_MAP_INIT_STR(pend, void *)
KHASH_MAP_INIT_STR(name, int)

/* allele call of one read or one fragment at a SNP */
enum {
    CALL_NOCOV = 0,   /* SNP not covered by the read */
    CALL_LOWQ,        /* base quality below the threshold */
    CALL_AMBIG,       /* bisulfite strand cannot tell the alleles apart */
    CALL_OTHER,       /* base matches neither allele (incl. deletions) */
    CALL_CONFLICT,    /* the two mates support different alleles */
    CALL_REF,
    CALL_ALT
};

typedef struct {
    const char *bam_fn, *cpg_fn, *snp_fn, *out_fn, *qc_fn, *dump_fn, *ref_fn;
    int min_mapq, min_snp_bq, min_cpg_bq, excl_flags, window, max_frag;
    int mask, threads, min_allele;
} opt_t;

typedef struct {
    int32_t chr;        /* index into the contig names of the SNP file */
    int32_t tid;        /* BAM target id, -1 if the contig is absent */
    hts_pos_t pos0;     /* 0-based SNP position */
    int64_t order;      /* input order, used to keep sorting stable */
    char ref, alt;
} snp_t;

typedef struct {
    hts_pos_t beg, end; /* 0-based reference span [beg, end) */
    int bottom;         /* 1 if the read comes from the bottom (OB/CTOB) strand */
    char *seq;          /* bases aligned to the reference, 'N' in deletions */
    uint8_t *qual;      /* base qualities aligned to the reference, 0 in deletions */
} mate_t;

typedef struct {
    char *qname;
    hts_pos_t beg, end; /* union of the mate spans */
    hts_pos_t mpos;     /* expected mate position while waiting for the mate */
    int heap_idx;       /* position in the pending-mate heap, -1 if not pending */
    int nmate;
    mate_t m[2];
    int ncall;
    int32_t *cpg;       /* called CpGs, as ascending indexes into the region CpG array */
    uint8_t *meth;      /* 1 = methylated, 0 = unmethylated */
} frag_t;

typedef struct {
    hts_pos_t *pos;     /* 0-based positions of the CpG C, ascending */
    uint8_t *mask;      /* 1 if the CpG overlaps a listed SNP */
    int n, m;
} cpgs_t;

typedef struct {
    int64_t reads, nsum, n2sum, sjd;
} acc_t;

typedef struct {
    opt_t *o;
    samFile *fp;
    sam_hdr_t *hdr;
    hts_idx_t *idx;
    htsFile *cpg_fp;
    tbx_t *tbx;
    FILE *out, *qc, *dump;
    bam1_t *b;
    khash_t(pend) *pend;  /* reads waiting for their mate, by QNAME */
    frag_t **heap;      /* the same reads as a min-heap on the expected mate position */
    int nheap, mheap;
    frag_t **act;       /* finished fragments that may still overlap a pending SNP */
    int nact, mact;
    cpgs_t cpg;
    int32_t *buf_i[2];  /* scratch buffers for the per-mate CpG calls */
    uint8_t *buf_m[2];
    int mbuf;
    kstring_t ks;
    kstring_t dks;      /* dump lines of the current SNP */
    char **names;       /* contig names of the SNP file */
    int nnames;
    char *warned_chr;   /* last contig reported as missing from the CpG file */
    /* run summary */
    int64_t n_reads, n_frags, n_ref, n_alt, n_snps, n_reported;
} ctx_t;

/***************************************************************
 *                        small helpers                        *
 ***************************************************************/

static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p) { fprintf(stderr, "[mhapasm] out of memory\n"); exit(1); }
    return p;
}

static void *xrealloc(void *p, size_t n)
{
    p = realloc(p, n ? n : 1);
    if (!p) { fprintf(stderr, "[mhapasm] out of memory\n"); exit(1); }
    return p;
}

static char *xstrdup(const char *s)
{
    size_t n = strlen(s) + 1;
    return memcpy(xmalloc(n), s, n);
}

static int is_pair(char a, char b, char x, char y)
{
    return (a == x && b == y) || (a == y && b == x);
}

static int is_acgt(char c)
{
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

/* first index i with v[i] >= x */
static int lower_bound(const hts_pos_t *v, int n, hts_pos_t x)
{
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (v[mid] < x) lo = mid + 1; else hi = mid;
    }
    return lo;
}

/* Look up a contig, accepting a missing or extra "chr" prefix and the usual
 * names of the mitochondrial genome. */
static int name2id_fallback(int (*lookup)(void *, const char *), void *h, const char *name)
{
    static const char *mito[] = {"chrM", "MT", "M", "chrMT"};
    int id = lookup(h, name);
    if (id >= 0) return id;
    for (int i = 0; i < 4; i++) {
        if (strcmp(name, mito[i])) continue;
        for (int j = 0; j < 4; j++)
            if (j != i && (id = lookup(h, mito[j])) >= 0) return id;
        return -1;
    }
    if (!strncmp(name, "chr", 3)) return lookup(h, name + 3);
    kstring_t ks = {0, 0, NULL};
    ksprintf(&ks, "chr%s", name);
    id = lookup(h, ks.s);
    free(ks.s);
    return id;
}

static int bam_lookup(void *h, const char *name) { return sam_hdr_name2tid((sam_hdr_t *) h, name); }
static int tbx_lookup(void *h, const char *name) { return tbx_name2id((tbx_t *) h, name); }

/***************************************************************
 *                        SNP list                             *
 ***************************************************************/

static int snp_cmp(const void *a, const void *b)
{
    const snp_t *x = a, *y = b;
    unsigned tx = (unsigned) x->tid, ty = (unsigned) y->tid;  /* -1 sorts last */
    if (tx != ty) return tx < ty ? -1 : 1;
    if (x->pos0 != y->pos0) return x->pos0 < y->pos0 ? -1 : 1;
    return x->order < y->order ? -1 : x->order > y->order;
}

/* Split on blanks into at most `max` fields; returns the number of fields. */
static int split_fields(char *p, char **f, int max)
{
    int nf = 0;
    while (nf < max && *p) {
        while (*p == ' ' || *p == '\t') p++;
        if (!*p) break;
        f[nf++] = p;
        while (*p && *p != ' ' && *p != '\t') p++;
        if (*p) *p++ = '\0';
    }
    return nf;
}

/* Load SNPs from "chr pos ref alt" (1-based pos) or from a VCF (plain or gz,
 * with or without header lines). Multi-allelic records are split into one
 * bi-allelic SNP per ALT; only single-base substitutions are kept. */
static int load_snps(ctx_t *c, const char *fn, snp_t **out, int64_t *n_out)
{
    htsFile *fp = hts_open(fn, "r");
    if (!fp) { fprintf(stderr, "[mhapasm] cannot open SNP file %s\n", fn); return -1; }
    khash_t(name) *nh = kh_init(name);
    int32_t *name_tid = NULL;
    kstring_t ks = {0, 0, NULL};
    snp_t *v = NULL;
    int64_t n = 0, m = 0, n_bad = 0, n_nocontig = 0;
    int layout = 0;   /* 0 = undecided, 1 = chr pos ref alt, 2 = VCF */
    while (hts_getline(fp, '\n', &ks) >= 0) {
        if (ks.l && ks.s[ks.l - 1] == '\r') ks.s[--ks.l] = '\0';
        if (ks.l == 0) continue;
        if (ks.s[0] == '#') {
            if (!strncmp(ks.s, "##fileformat=VCF", 16) || !strncmp(ks.s, "#CHROM", 6)) layout = 2;
            continue;
        }
        char *f[5];
        int nf = split_fields(ks.s, f, 5);
        if (!layout)   /* no header: a VCF holds the ID, not a base, in column 3 */
            layout = nf >= 5 && !(strlen(f[2]) == 1 && is_acgt(toupper((unsigned char) f[2][0]))) ? 2 : 1;
        if (nf < (layout == 2 ? 5 : 4)) { n_bad++; continue; }
        char *chr = f[0], *ref = f[layout == 2 ? 3 : 2], *alt = f[layout == 2 ? 4 : 3], *end;
        long long pos = strtoll(f[1], &end, 10);
        char r = toupper((unsigned char) ref[0]);
        if (*end || pos < 1 || ref[1] || !is_acgt(r)) { n_bad++; continue; }
        khint_t k = kh_get(name, nh, chr);
        int32_t ci;
        if (k == kh_end(nh)) {   /* new contig: keep one copy of its name */
            int ret;
            ci = c->nnames++;
            c->names = xrealloc(c->names, c->nnames * sizeof(char *));
            name_tid = xrealloc(name_tid, c->nnames * sizeof(int32_t));
            c->names[ci] = xstrdup(chr);
            name_tid[ci] = name2id_fallback(bam_lookup, c->hdr, chr);
            k = kh_put(name, nh, c->names[ci], &ret);
            kh_val(nh, k) = ci;
        } else {
            ci = kh_val(nh, k);
        }
        for (char *a = alt; ; ) {
            char *comma = strchr(a, ',');
            size_t len = comma ? (size_t) (comma - a) : strlen(a);
            char t = toupper((unsigned char) a[0]);
            if (len != 1 || !is_acgt(t) || t == r) {
                n_bad++;
            } else {
                if (n == m) { m = m ? m * 2 : 1024; v = xrealloc(v, m * sizeof(snp_t)); }
                snp_t *s = &v[n];
                s->chr = ci;
                s->tid = name_tid[ci];
                s->pos0 = pos - 1;
                s->ref = r;
                s->alt = t;
                s->order = n++;
                if (s->tid < 0) n_nocontig++;
            }
            if (!comma) break;
            a = comma + 1;
        }
    }
    free(ks.s);
    free(name_tid);
    kh_destroy(name, nh);
    hts_close(fp);
    if (n_bad)
        fprintf(stderr, "[mhapasm] skipped %" PRId64 " SNP record(s)/ALT allele(s) that are not single-base substitutions\n", n_bad);
    if (n_nocontig)
        fprintf(stderr, "[mhapasm] %" PRId64 " SNP(s) on contigs absent from the BAM header (reported with zero counts)\n", n_nocontig);
    if (n) qsort(v, n, sizeof(snp_t), snp_cmp);
    *out = v;
    *n_out = n;
    return 0;
}

/***************************************************************
 *                        CpG positions                        *
 ***************************************************************/

/* Load the CpGs of [beg, end) (0-based) from the tabix-indexed CpG file. */
static int load_cpgs(ctx_t *c, const char *chr, hts_pos_t beg, hts_pos_t end)
{
    cpgs_t *g = &c->cpg;
    g->n = 0;
    int tid = name2id_fallback(tbx_lookup, c->tbx, chr);
    if (tid < 0) {
        if (!c->warned_chr || strcmp(c->warned_chr, chr)) {
            fprintf(stderr, "[mhapasm] warning: contig %s is absent from the CpG file; its SNPs get zero counts\n", chr);
            free(c->warned_chr);
            c->warned_chr = xstrdup(chr);
        }
        return 0;
    }
    hts_itr_t *itr = tbx_itr_queryi(c->tbx, tid, beg < 0 ? 0 : beg, end);
    if (!itr) return 0;
    while (tbx_itr_next(c->cpg_fp, c->tbx, itr, &c->ks) >= 0) {
        char *p = strchr(c->ks.s, '\t'), *q;
        if (!p) continue;
        long long pos = strtoll(p + 1, &q, 10);
        if (q == p + 1 || pos < 1) continue;
        if (g->n && pos - 1 <= g->pos[g->n - 1]) continue;  /* duplicates */
        if (g->n == g->m) {
            g->m = g->m ? g->m * 2 : 4096;
            g->pos = xrealloc(g->pos, g->m * sizeof(hts_pos_t));
            g->mask = xrealloc(g->mask, g->m);
        }
        g->pos[g->n] = pos - 1;
        g->mask[g->n++] = 0;
    }
    tbx_itr_destroy(itr);
    return g->n;
}

/* Mask CpGs whose C or G sits on a listed SNP, so both alleles are compared
 * on the same CpG set. */
static void mask_snp_cpgs(cpgs_t *g, const snp_t *s, int64_t ns)
{
    for (int64_t k = 0; k < ns; k++) {
        int i = lower_bound(g->pos, g->n, s[k].pos0 - 1);
        for (; i < g->n && g->pos[i] <= s[k].pos0; i++) g->mask[i] = 1;
    }
}

/***************************************************************
 *                        reads and fragments                  *
 ***************************************************************/

/* Bisulfite strand of a read. Strand tags written by the aligner win
 * (Bismark XG, bwa-meth/BISCUIT YD, BSMAP ZS). Otherwise use the wgbs_tools
 * is_bottom rule (OB = read1 reverse / read2 forward), without requiring the
 * proper-pair bit. */
static int read_is_bottom(const bam1_t *b)
{
    uint8_t *t;
    if ((t = bam_aux_get(b, "XG")) && *t == 'Z') {
        const char *v = bam_aux2Z(t);
        if (v && !strcmp(v, "GA")) return 1;
        if (v && !strcmp(v, "CT")) return 0;
    }
    if ((t = bam_aux_get(b, "YD"))) {
        char v = 0;
        if (*t == 'A') v = bam_aux2A(t);
        else if (*t == 'Z') { const char *z = bam_aux2Z(t); v = z ? z[0] : 0; }
        if (v == 'r' || v == 'R') return 1;
        if (v == 'f' || v == 'F') return 0;
    }
    if ((t = bam_aux_get(b, "ZS")) && *t == 'Z') {
        const char *v = bam_aux2Z(t);
        if (v && v[0] == '-') return 1;
        if (v && v[0] == '+') return 0;
    }
    uint16_t fl = b->core.flag;
    if (fl & BAM_FPAIRED) {
        if (fl & BAM_FREAD1) return (fl & BAM_FREVERSE) != 0;
        if (fl & BAM_FREAD2) return (fl & BAM_FREVERSE) == 0;
    }
    return (fl & BAM_FREVERSE) != 0;
}

/* Lay the read out on the reference (patter_utils::clean_CIGAR): keep M/=/X,
 * drop I/S/H/P, fill D/N with 'N'. */
static int mate_init(mate_t *m, const bam1_t *b)
{
    const uint32_t *cig = bam_get_cigar(b);
    int ncig = b->core.n_cigar;
    hts_pos_t rlen = bam_cigar2rlen(ncig, cig);
    if (rlen <= 0) return -1;
    m->beg = b->core.pos;
    m->end = b->core.pos + rlen;
    m->bottom = read_is_bottom(b);
    m->seq = xmalloc(rlen);
    m->qual = xmalloc(rlen);
    const uint8_t *s = bam_get_seq(b), *q = bam_get_qual(b);
    hts_pos_t r = 0;
    int32_t qi = 0;
    for (int i = 0; i < ncig; i++) {
        int op = bam_cigar_op(cig[i]);
        int32_t len = bam_cigar_oplen(cig[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (int32_t j = 0; j < len; j++, r++, qi++) {
                m->seq[r] = seq_nt16_str[bam_seqi(s, qi)];
                m->qual[r] = q[qi];
            }
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            memset(m->seq + r, 'N', len);
            memset(m->qual + r, 0, len);
            r += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            qi += len;
        }
    }
    return 0;
}

static frag_t *frag_new(const bam1_t *b, const mate_t *m)
{
    frag_t *f = xmalloc(sizeof(frag_t));
    memset(f, 0, sizeof(frag_t));
    f->qname = xstrdup(bam_get_qname(b));
    f->m[0] = *m;
    f->nmate = 1;
    return f;
}

static void frag_free(frag_t *f)
{
    for (int i = 0; i < f->nmate; i++) { free(f->m[i].seq); free(f->m[i].qual); }
    free(f->cpg);
    free(f->meth);
    free(f->qname);
    free(f);
}

/* CpG calls of one read (patter::compareSeqToRef). The whole CpG must lie in
 * the aligned span and the read must show the CpG context itself: C/T + G on
 * the top strand, C + G/A on the bottom strand. Only informative calls are
 * returned. */
static int mate_calls(const mate_t *m, const cpgs_t *g, int min_bq, int32_t *idx, uint8_t *st)
{
    int n = 0;
    for (int i = lower_bound(g->pos, g->n, m->beg); i < g->n && g->pos[i] + 1 < m->end; i++) {
        if (g->mask[i]) continue;
        hts_pos_t c = g->pos[i] - m->beg;
        char base, ctx;
        uint8_t q;
        if (!m->bottom) { base = m->seq[c]; ctx = m->seq[c + 1]; q = m->qual[c]; }
        else            { base = m->seq[c + 1]; ctx = m->seq[c]; q = m->qual[c + 1]; }
        if (q < min_bq) continue;
        int s = -1;
        if (!m->bottom) { if (ctx == 'G') s = base == 'C' ? 1 : base == 'T' ? 0 : -1; }
        else            { if (ctx == 'C') s = base == 'G' ? 1 : base == 'A' ? 0 : -1; }
        if (s < 0) continue;
        idx[n] = i;
        st[n++] = (uint8_t) s;
    }
    return n;
}

static void grow_buffers(ctx_t *c, int need)
{
    if (need <= c->mbuf) return;
    c->mbuf = need * 2;
    for (int k = 0; k < 2; k++) {
        c->buf_i[k] = xrealloc(c->buf_i[k], c->mbuf * sizeof(int32_t));
        c->buf_m[k] = xrealloc(c->buf_m[k], c->mbuf);
    }
}

/* Upper bound on the number of CpGs inside [beg, end). */
static int cpgs_in_span(const cpgs_t *g, hts_pos_t beg, hts_pos_t end)
{
    return lower_bound(g->pos, g->n, end) - lower_bound(g->pos, g->n, beg);
}

/* Fragment is complete: merge the CpG calls of its mates (merge_PE: a CpG
 * seen by one mate is kept, a CpG on which the mates disagree is dropped) and
 * make it available to the SNP evaluation. */
static void frag_finish(ctx_t *c, frag_t *f)
{
    f->beg = f->m[0].beg;
    f->end = f->m[0].end;
    int n[2] = {0, 0};
    for (int k = 0; k < f->nmate; k++) {
        if (f->m[k].beg < f->beg) f->beg = f->m[k].beg;
        if (f->m[k].end > f->end) f->end = f->m[k].end;
        grow_buffers(c, cpgs_in_span(&c->cpg, f->m[k].beg, f->m[k].end) + 1);
        n[k] = mate_calls(&f->m[k], &c->cpg, c->o->min_cpg_bq, c->buf_i[k], c->buf_m[k]);
    }
    int tot = n[0] + n[1];
    f->cpg = xmalloc(tot * sizeof(int32_t));
    f->meth = xmalloc(tot);
    int i = 0, j = 0, k = 0;
    while (i < n[0] || j < n[1]) {
        if (j >= n[1] || (i < n[0] && c->buf_i[0][i] < c->buf_i[1][j])) {
            f->cpg[k] = c->buf_i[0][i]; f->meth[k++] = c->buf_m[0][i++];
        } else if (i >= n[0] || c->buf_i[1][j] < c->buf_i[0][i]) {
            f->cpg[k] = c->buf_i[1][j]; f->meth[k++] = c->buf_m[1][j++];
        } else {
            if (c->buf_m[0][i] == c->buf_m[1][j]) { f->cpg[k] = c->buf_i[0][i]; f->meth[k++] = c->buf_m[0][i]; }
            i++; j++;
        }
    }
    f->ncall = k;
    if (c->nact == c->mact) {
        c->mact = c->mact ? c->mact * 2 : 256;
        c->act = xrealloc(c->act, c->mact * sizeof(frag_t *));
    }
    c->act[c->nact++] = f;
    c->n_frags++;
}

/* Reads waiting for their mate are kept in a hash (by QNAME) and in a
 * min-heap on the expected mate position, so that reads whose mate never
 * shows up are found without scanning the hash. */
static void heap_swap(ctx_t *c, int i, int j)
{
    frag_t *t = c->heap[i];
    c->heap[i] = c->heap[j];
    c->heap[j] = t;
    c->heap[i]->heap_idx = i;
    c->heap[j]->heap_idx = j;
}

static void heap_up(ctx_t *c, int i)
{
    while (i > 0 && c->heap[(i - 1) / 2]->mpos > c->heap[i]->mpos) {
        heap_swap(c, i, (i - 1) / 2);
        i = (i - 1) / 2;
    }
}

static void heap_down(ctx_t *c, int i)
{
    for (;;) {
        int l = 2 * i + 1, r = l + 1, m = i;
        if (l < c->nheap && c->heap[l]->mpos < c->heap[m]->mpos) m = l;
        if (r < c->nheap && c->heap[r]->mpos < c->heap[m]->mpos) m = r;
        if (m == i) return;
        heap_swap(c, i, m);
        i = m;
    }
}

static void heap_push(ctx_t *c, frag_t *f)
{
    if (c->nheap == c->mheap) {
        c->mheap = c->mheap ? c->mheap * 2 : 256;
        c->heap = xrealloc(c->heap, c->mheap * sizeof(frag_t *));
    }
    c->heap[c->nheap] = f;
    f->heap_idx = c->nheap++;
    heap_up(c, f->heap_idx);
}

static void heap_remove(ctx_t *c, frag_t *f)
{
    int i = f->heap_idx;
    if (i != --c->nheap) {
        c->heap[i] = c->heap[c->nheap];
        c->heap[i]->heap_idx = i;
        heap_down(c, i);
        heap_up(c, i);
    }
    f->heap_idx = -1;
}

/* Add a filtered read: pair it with its mate (by QNAME) or keep it as a
 * single-read fragment, like wgbs_tools match_maker. */
static void add_read(ctx_t *c, const bam1_t *b)
{
    const bam1_core_t *co = &b->core;
    mate_t m;
    if (mate_init(&m, b) < 0) return;
    c->n_reads++;
    int expect = (co->flag & BAM_FPAIRED) && !(co->flag & BAM_FMUNMAP) && co->mtid == co->tid
                 && llabs((long long) (co->mpos - co->pos)) <= c->o->max_frag;
    if (expect) {
        khint_t k = kh_get(pend, c->pend, bam_get_qname(b));
        if (k != kh_end(c->pend)) {
            frag_t *f = kh_val(c->pend, k);
            kh_del(pend, c->pend, k);
            heap_remove(c, f);
            f->m[f->nmate++] = m;
            frag_finish(c, f);
            return;
        }
        if (co->mpos >= co->pos) {   /* mate still to come */
            frag_t *f = frag_new(b, &m);
            f->mpos = co->mpos;
            int ret;
            k = kh_put(pend, c->pend, f->qname, &ret);
            if (ret == 0) {          /* same QNAME already waiting: keep both as singles */
                frag_t *old = kh_val(c->pend, k);
                kh_key(c->pend, k) = f->qname;
                heap_remove(c, old);
                frag_finish(c, old);
            }
            kh_val(c->pend, k) = f;
            heap_push(c, f);
            return;
        }
    }
    frag_finish(c, frag_new(b, &m));
}

/* Reads whose mate was expected before position `pos` but never came (e.g.
 * the mate failed the filters) become single-read fragments. */
static void flush_pending(ctx_t *c, hts_pos_t pos)
{
    while (c->nheap && c->heap[0]->mpos < pos) {
        frag_t *f = c->heap[0];
        heap_remove(c, f);
        khint_t k = kh_get(pend, c->pend, f->qname);
        if (k != kh_end(c->pend)) kh_del(pend, c->pend, k);
        frag_finish(c, f);
    }
}

/* Drop fragments that end before `pos`; no later SNP can use them. */
static void evict(ctx_t *c, hts_pos_t pos)
{
    int j = 0;
    for (int i = 0; i < c->nact; i++) {
        if (c->act[i]->end <= pos) frag_free(c->act[i]);
        else c->act[j++] = c->act[i];
    }
    c->nact = j;
}

/***************************************************************
 *                        allele calls                         *
 ***************************************************************/

/* Bases allowed for allele `let` (snp_patter::compareSeqToRef): bisulfite can
 * turn C into T on the top strand and G into A on the bottom strand. */
static int base_allowed(char base, char let, char other, int bottom)
{
    if (let == 'C' && other != 'T' && !bottom) return base == 'C' || base == 'T';
    if (let == 'G' && other != 'A' && bottom) return base == 'G' || base == 'A';
    return base == let;
}

static int mate_allele(const mate_t *m, const snp_t *s, int min_bq)
{
    if (s->pos0 < m->beg || s->pos0 >= m->end) return CALL_NOCOV;
    hts_pos_t i = s->pos0 - m->beg;
    if (m->qual[i] < min_bq) return CALL_LOWQ;
    if (is_pair(s->ref, s->alt, 'C', 'T') && !m->bottom) return CALL_AMBIG;
    if (is_pair(s->ref, s->alt, 'G', 'A') && m->bottom) return CALL_AMBIG;
    if (base_allowed(m->seq[i], s->ref, s->alt, m->bottom)) return CALL_REF;
    if (base_allowed(m->seq[i], s->alt, s->ref, m->bottom)) return CALL_ALT;
    return CALL_OTHER;
}

/* Combine the mates (snp_patter::proc2lines): one informative mate decides,
 * two informative mates must agree. */
static int frag_allele(const frag_t *f, const snp_t *s, int min_bq)
{
    int c1 = mate_allele(&f->m[0], s, min_bq);
    int c2 = f->nmate > 1 ? mate_allele(&f->m[1], s, min_bq) : CALL_NOCOV;
    int i1 = c1 == CALL_REF || c1 == CALL_ALT;
    int i2 = c2 == CALL_REF || c2 == CALL_ALT;
    if (i1 && i2) return c1 == c2 ? c1 : CALL_CONFLICT;
    if (i1) return c1;
    if (i2) return c2;
    return c1 > c2 ? c1 : c2;
}

/***************************************************************
 *                        statistics                           *
 ***************************************************************/

static void print_acc(FILE *fp, const char *chr, hts_pos_t pos1, char allele, const acc_t *a)
{
    fprintf(fp, "%s\t%" PRIhts_pos "\t%c\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t%" PRId64,
            chr, pos1, allele, a->reads, a->nsum, a->n2sum, a->sjd);
    if (a->nsum > 0) fprintf(fp, "\t%.6g", (double) a->sjd / (double) a->nsum);
    else fputs("\tNA", fp);
    if (a->n2sum > 0)
        fprintf(fp, "\t%.6g\t%.6g", (double) a->nsum * (double) a->nsum / (double) a->n2sum,
                (double) a->sjd * (double) a->nsum / (double) a->n2sum);
    else fputs("\tNA\tNA", fp);
    fputc('\n', fp);
}

/* One line per allele-assigned fragment; fragments without a CpG call get
 * '.' as start, end and haplotype. */
static void dump_frag(ctx_t *c, const frag_t *f, int lo, int hi, const snp_t *s, char allele)
{
    const cpgs_t *g = &c->cpg;
    kstring_t *ks = &c->dks;
    const char *chr = c->names[s->chr];
    if (hi > lo) {
        ksprintf(ks, "%s\t%" PRIhts_pos "\t%" PRIhts_pos "\t", chr, g->pos[f->cpg[lo]] + 1, g->pos[f->cpg[hi - 1]] + 1);
        for (int j = lo, i = f->cpg[lo]; j < hi; i++) {
            if (f->cpg[j] == i) kputc(f->meth[j++] ? '1' : '0', ks);
            else kputc('.', ks);
        }
    } else {
        ksprintf(ks, "%s\t.\t.\t.", chr);
    }
    ksprintf(ks, "\t%c\t%s:%" PRIhts_pos ":%c:%c\t%c\t%s\n", f->m[0].bottom ? '-' : '+',
             chr, s->pos0 + 1, s->ref, s->alt, allele, f->qname);
}

/* Assign every fragment covering the SNP to an allele and accumulate the
 * M-score statistics over its CpG calls (restricted to the window if set). */
static void eval_snp(ctx_t *c, const snp_t *s)
{
    const opt_t *o = c->o;
    acc_t a[2];
    int64_t qc[CALL_ALT + 1];
    memset(a, 0, sizeof(a));
    memset(qc, 0, sizeof(qc));
    c->dks.l = 0;
    for (int i = 0; i < c->nact; i++) {
        const frag_t *f = c->act[i];
        if (s->pos0 < f->beg || s->pos0 >= f->end) continue;
        int call = frag_allele(f, s, o->min_snp_bq);
        qc[call]++;
        if (call != CALL_REF && call != CALL_ALT) continue;
        int lo = 0, hi = f->ncall;
        if (o->window > 0) {
            while (lo < hi && c->cpg.pos[f->cpg[lo]] < s->pos0 - o->window) lo++;
            while (hi > lo && c->cpg.pos[f->cpg[hi - 1]] > s->pos0 + o->window) hi--;
        }
        if (c->dump) dump_frag(c, f, lo, hi, s, call == CALL_ALT ? s->alt : s->ref);
        int64_t n = hi - lo, z = 0;
        if (n == 0) continue;
        for (int j = lo; j < hi; j++) if (f->meth[j]) { z = 1; break; }
        acc_t *x = &a[call == CALL_ALT];
        x->reads++;
        x->nsum += n;
        x->n2sum += n * n;
        x->sjd += n * z;
    }
    const char *chr = c->names[s->chr];
    if (qc[CALL_REF] >= o->min_allele && qc[CALL_ALT] >= o->min_allele) {
        print_acc(c->out, chr, s->pos0 + 1, s->ref, &a[0]);
        print_acc(c->out, chr, s->pos0 + 1, s->alt, &a[1]);
        if (c->dump && c->dks.l) fputs(c->dks.s, c->dump);
        c->n_reported++;
    }
    if (c->qc) {
        int64_t cov = 0;
        for (int k = CALL_LOWQ; k <= CALL_ALT; k++) cov += qc[k];
        fprintf(c->qc, "%s\t%" PRIhts_pos "\t%c\t%c\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\n",
                chr, s->pos0 + 1, s->ref, s->alt, cov, qc[CALL_REF], qc[CALL_ALT], qc[CALL_CONFLICT],
                qc[CALL_AMBIG], qc[CALL_OTHER] + qc[CALL_LOWQ]);
    }
    c->n_ref += qc[CALL_REF];
    c->n_alt += qc[CALL_ALT];
    c->n_snps++;
}

/***************************************************************
 *                        driver                               *
 ***************************************************************/

static int pass_filters(const bam1_t *b, const opt_t *o)
{
    if (b->core.flag & o->excl_flags) return 0;
    if (b->core.qual < o->min_mapq) return 0;
    return 1;
}

/* Stream the reads of one region and evaluate its SNPs in order. A SNP is
 * evaluated once the reads have moved more than max_frag past it, so every
 * fragment overlapping it is complete. */
static int process_region(ctx_t *c, int tid, hts_pos_t rbeg, hts_pos_t rend, const snp_t *s, int64_t ns)
{
    const opt_t *o = c->o;
    hts_pos_t L = o->max_frag;
    load_cpgs(c, sam_hdr_tid2name(c->hdr, tid), rbeg - L, rend + L);
    if (o->mask) mask_snp_cpgs(&c->cpg, s, ns);

    hts_itr_t *itr = sam_itr_queryi(c->idx, tid, rbeg, rend);
    if (!itr) { fprintf(stderr, "[mhapasm] failed to query %s\n", sam_hdr_tid2name(c->hdr, tid)); return -1; }
    int64_t k = 0;
    int ret;
    while ((ret = sam_itr_next(c->fp, itr, c->b)) >= 0) {
        hts_pos_t pos = c->b->core.pos;
        if (s[k].pos0 + L < pos) {
            flush_pending(c, pos);
            while (k < ns && s[k].pos0 + L < pos) {
                eval_snp(c, &s[k]);
                k++;
                evict(c, k < ns ? s[k].pos0 : HTS_POS_MAX);
            }
            if (k == ns) break;
        }
        if (pass_filters(c->b, o)) add_read(c, c->b);
    }
    hts_itr_destroy(itr);
    if (ret < -1) { fprintf(stderr, "[mhapasm] error reading %s\n", o->bam_fn); return -1; }
    flush_pending(c, HTS_POS_MAX);
    for (; k < ns; k++) eval_snp(c, &s[k]);
    evict(c, HTS_POS_MAX);
    return 0;
}

static void usage(FILE *fp)
{
    fprintf(fp,
"mhapasm " MHAPASM_VERSION ": allele-specific M-score statistics from bisulfite BAMs\n"
"\n"
"Usage: mhapasm [options] -b <in.bam> -c <CpG.gz> -s <snps>\n"
"\n"
"Input:\n"
"  -b, --bam FILE          coordinate-sorted and indexed BAM/CRAM\n"
"  -c, --cpg FILE          bgzipped, tabix-indexed CpG file (chr, 1-based C position, ...)\n"
"  -s, --snp FILE          SNPs as 'chr pos ref alt' (1-based) or VCF (header optional);\n"
"                          plain or gzipped; multi-allelic records are split\n"
"  -T, --reference FILE    reference FASTA (CRAM input only)\n"
"Output:\n"
"  -o, --out FILE          per-allele statistics [stdout]\n"
"      --qc FILE           per-SNP allele-calling counts\n"
"      --dump FILE         per-fragment CpG haplotypes assigned to each allele\n"
"Filters:\n"
"  -q, --min-mapq INT      minimum mapping quality [10]\n"
"  -F, --exclude-flags INT skip reads with any of these flags [0xF04]\n"
"  -Q, --min-snp-bq INT    minimum base quality at the SNP [0]\n"
"  -B, --min-cpg-bq INT    minimum base quality of a CpG call [0]\n"
"  -w, --window INT        only use CpGs within INT bp of the SNP; 0 = whole fragment [0]\n"
"  -L, --max-frag INT      maximum mate distance of a fragment, in bp [1000]\n"
"      --no-mask           keep CpGs that overlap a listed SNP\n"
"      --min-allele INT    report a SNP only if each allele is carried by >= INT\n"
"                          fragments (the --qc file still lists every SNP) [0]\n"
"Other:\n"
"  -@, --threads INT       extra BAM decompression threads [0]\n"
"  -h, --help              show this help\n"
"  -v, --version           show the version\n"
"\n"
"Output columns (two lines per SNP, REF first):\n"
"  chr pos allele reads Nsum N2sum Sjd mscore kappa Y_prime\n"
"  N_t = CpGs called on fragment t, Z_t = 1 if any of them is methylated;\n"
"  reads = #fragments, Nsum = sum N_t, N2sum = sum N_t^2, Sjd = sum N_t*Z_t,\n"
"  mscore = Sjd/Nsum, kappa = Nsum^2/N2sum, Y_prime = Sjd*Nsum/N2sum.\n");
}

int main(int argc, char **argv)
{
    opt_t o;
    memset(&o, 0, sizeof(o));
    o.min_mapq = 10;
    o.excl_flags = 0xF04;
    o.max_frag = 1000;
    o.mask = 1;

    static const struct option lopts[] = {
        {"bam", required_argument, NULL, 'b'},
        {"cpg", required_argument, NULL, 'c'},
        {"snp", required_argument, NULL, 's'},
        {"reference", required_argument, NULL, 'T'},
        {"out", required_argument, NULL, 'o'},
        {"qc", required_argument, NULL, 1},
        {"dump", required_argument, NULL, 2},
        {"no-mask", no_argument, NULL, 3},
        {"min-allele", required_argument, NULL, 4},
        {"min-mapq", required_argument, NULL, 'q'},
        {"exclude-flags", required_argument, NULL, 'F'},
        {"min-snp-bq", required_argument, NULL, 'Q'},
        {"min-cpg-bq", required_argument, NULL, 'B'},
        {"window", required_argument, NULL, 'w'},
        {"max-frag", required_argument, NULL, 'L'},
        {"threads", required_argument, NULL, '@'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}
    };
    int ch;
    while ((ch = getopt_long(argc, argv, "b:c:s:T:o:q:F:Q:B:w:L:@:hv", lopts, NULL)) >= 0) {
        switch (ch) {
        case 'b': o.bam_fn = optarg; break;
        case 'c': o.cpg_fn = optarg; break;
        case 's': o.snp_fn = optarg; break;
        case 'T': o.ref_fn = optarg; break;
        case 'o': o.out_fn = optarg; break;
        case 1: o.qc_fn = optarg; break;
        case 2: o.dump_fn = optarg; break;
        case 3: o.mask = 0; break;
        case 4: o.min_allele = atoi(optarg); break;
        case 'q': o.min_mapq = atoi(optarg); break;
        case 'F': o.excl_flags = (int) strtol(optarg, NULL, 0); break;
        case 'Q': o.min_snp_bq = atoi(optarg); break;
        case 'B': o.min_cpg_bq = atoi(optarg); break;
        case 'w': o.window = atoi(optarg); break;
        case 'L': o.max_frag = atoi(optarg); break;
        case '@': o.threads = atoi(optarg); break;
        case 'h': usage(stdout); return 0;
        case 'v': puts(MHAPASM_VERSION); return 0;
        default: usage(stderr); return 1;
        }
    }
    if (!o.bam_fn || !o.cpg_fn || !o.snp_fn || optind != argc) { usage(stderr); return 1; }
    if (o.max_frag < 1 || o.window < 0) { fprintf(stderr, "[mhapasm] -L must be >= 1 and -w >= 0\n"); return 1; }

    clock_t t0 = clock();
    ctx_t c;
    memset(&c, 0, sizeof(c));
    c.o = &o;
    if (!(c.fp = sam_open(o.bam_fn, "r"))) { fprintf(stderr, "[mhapasm] cannot open %s\n", o.bam_fn); return 1; }
    if (o.ref_fn && hts_set_fai_filename(c.fp, o.ref_fn) < 0) { fprintf(stderr, "[mhapasm] cannot use reference %s\n", o.ref_fn); return 1; }
    if (o.threads > 0) hts_set_threads(c.fp, o.threads);
    if (!(c.hdr = sam_hdr_read(c.fp))) { fprintf(stderr, "[mhapasm] cannot read the header of %s\n", o.bam_fn); return 1; }
    if (!(c.idx = sam_index_load(c.fp, o.bam_fn))) { fprintf(stderr, "[mhapasm] %s is not indexed (run samtools index)\n", o.bam_fn); return 1; }
    if (!(c.cpg_fp = hts_open(o.cpg_fn, "r"))) { fprintf(stderr, "[mhapasm] cannot open %s\n", o.cpg_fn); return 1; }
    if (!(c.tbx = tbx_index_load(o.cpg_fn))) { fprintf(stderr, "[mhapasm] %s has no tabix index\n", o.cpg_fn); return 1; }
    c.out = o.out_fn ? fopen(o.out_fn, "w") : stdout;
    if (!c.out) { fprintf(stderr, "[mhapasm] cannot write %s\n", o.out_fn); return 1; }
    if (o.qc_fn && !(c.qc = fopen(o.qc_fn, "w"))) { fprintf(stderr, "[mhapasm] cannot write %s\n", o.qc_fn); return 1; }
    if (o.dump_fn && !(c.dump = fopen(o.dump_fn, "w"))) { fprintf(stderr, "[mhapasm] cannot write %s\n", o.dump_fn); return 1; }
    c.b = bam_init1();
    c.pend = kh_init(pend);

    int64_t ns = 0;
    snp_t *s = NULL;
    if (load_snps(&c, o.snp_fn, &s, &ns) < 0) return 1;

    fputs("#chr\tpos\tallele\treads\tNsum\tN2sum\tSjd\tmscore\tkappa\tY_prime\n", c.out);
    if (c.qc) fputs("#chr\tpos\tref\talt\tn_cover\tn_ref\tn_alt\tn_conflict\tn_bs_ambiguous\tn_other\n", c.qc);
    if (c.dump) fputs("#chr\tstart\tend\thap\tstrand\tsnp\tallele\tqname\n", c.dump);

    /* SNPs closer than 4*max_frag share one read query */
    hts_pos_t pad = 2 * (hts_pos_t) o.max_frag;
    int64_t i = 0, rc = 0;
    while (i < ns && s[i].tid >= 0 && rc == 0) {
        int64_t j = i + 1;
        hts_pos_t rend = s[i].pos0 + pad + 1;
        while (j < ns && s[j].tid == s[i].tid && s[j].pos0 - pad <= rend) rend = s[j++].pos0 + pad + 1;
        hts_pos_t rbeg = s[i].pos0 - pad < 0 ? 0 : s[i].pos0 - pad;
        rc = process_region(&c, s[i].tid, rbeg, rend, s + i, j - i);
        i = j;
    }
    acc_t zero = {0, 0, 0, 0};
    for (; i < ns && rc == 0; i++) {   /* contigs absent from the BAM */
        if (o.min_allele == 0) {
            print_acc(c.out, c.names[s[i].chr], s[i].pos0 + 1, s[i].ref, &zero);
            print_acc(c.out, c.names[s[i].chr], s[i].pos0 + 1, s[i].alt, &zero);
            c.n_reported++;
        }
        if (c.qc)
            fprintf(c.qc, "%s\t%" PRIhts_pos "\t%c\t%c\t0\t0\t0\t0\t0\t0\n", c.names[s[i].chr], s[i].pos0 + 1, s[i].ref, s[i].alt);
        c.n_snps++;
    }

    fprintf(stderr, "[mhapasm] %" PRId64 " SNPs (%" PRId64 " reported), %" PRId64 " reads, %" PRId64 " fragments, %" PRId64 " REF / %" PRId64 " ALT assignments, %.2f s CPU\n",
            c.n_snps, c.n_reported, c.n_reads, c.n_frags, c.n_ref, c.n_alt, (double) (clock() - t0) / CLOCKS_PER_SEC);

    for (int k = 0; k < c.nnames; k++) free(c.names[k]);
    free(c.names);
    free(s);
    free(c.heap);
    free(c.dks.s);
    kh_destroy(pend, c.pend);
    free(c.act);
    free(c.cpg.pos);
    free(c.cpg.mask);
    for (int k = 0; k < 2; k++) { free(c.buf_i[k]); free(c.buf_m[k]); }
    free(c.ks.s);
    free(c.warned_chr);
    bam_destroy1(c.b);
    tbx_destroy(c.tbx);
    hts_close(c.cpg_fp);
    hts_idx_destroy(c.idx);
    sam_hdr_destroy(c.hdr);
    sam_close(c.fp);
    if (c.qc) fclose(c.qc);
    if (c.dump) fclose(c.dump);
    if (c.out != stdout && fclose(c.out) != 0) { fprintf(stderr, "[mhapasm] error closing %s\n", o.out_fn); return 1; }
    return rc ? 1 : 0;
}
