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
#include "bsread.h"

#define MHAPASM_VERSION "0.1.0"

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
    bs_pairer_t *pairer;
    frag_t **act;       /* finished fragments that may still overlap a pending SNP */
    int nact, mact;
    cpgs_t cpg;
    bs_buf_t buf;
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

static int is_pair(char a, char b, char x, char y)
{
    return (a == x && b == y) || (a == y && b == x);
}

static int is_acgt(char c)
{
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

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
            c->names = bs_realloc(c->names, c->nnames * sizeof(char *));
            name_tid = bs_realloc(name_tid, c->nnames * sizeof(int32_t));
            c->names[ci] = bs_strdup(chr);
            name_tid[ci] = bs_bam_tid(c->hdr, chr);
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
                if (n == m) { m = m ? m * 2 : 1024; v = bs_realloc(v, m * sizeof(snp_t)); }
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
    int n = bs_cpgs_load(c->cpg_fp, c->tbx, chr, beg, end, &c->cpg, &c->ks);
    if (n < 0) {
        if (!c->warned_chr || strcmp(c->warned_chr, chr)) {
            fprintf(stderr, "[mhapasm] warning: contig %s is absent from the CpG file; its SNPs get zero counts\n", chr);
            free(c->warned_chr);
            c->warned_chr = bs_strdup(chr);
        }
        return 0;
    }
    return n;
}

/* Mask CpGs whose C or G sits on a listed SNP, so both alleles are compared
 * on the same CpG set. */
static void mask_snp_cpgs(cpgs_t *g, const snp_t *s, int64_t ns)
{
    for (int64_t k = 0; k < ns; k++) {
        int i = bs_lower_bound(g->pos, g->n, s[k].pos0 - 1);
        for (; i < g->n && g->pos[i] <= s[k].pos0; i++) g->mask[i] = 1;
    }
}

/***************************************************************
 *                        reads and fragments                  *
 ***************************************************************/

/* A fragment is complete: call its CpGs and make it available to the SNP
 * evaluation. */
static void frag_done(frag_t *f, void *data)
{
    ctx_t *c = data;
    bs_frag_call(f, &c->cpg, c->o->min_cpg_bq, &c->buf);
    if (c->nact == c->mact) {
        c->mact = c->mact ? c->mact * 2 : 256;
        c->act = bs_realloc(c->act, c->mact * sizeof(frag_t *));
    }
    c->act[c->nact++] = f;
    c->n_frags++;
}

static void add_read(ctx_t *c, const bam1_t *b)
{
    mate_t m;
    if (bs_mate_init(&m, b) < 0) return;
    c->n_reads++;
    bs_pairer_add(c->pairer, b, &m);
}

/* Drop fragments that end before `pos`; no later SNP can use them. */
static void evict(ctx_t *c, hts_pos_t pos)
{
    int j = 0;
    for (int i = 0; i < c->nact; i++) {
        if (c->act[i]->end <= pos) bs_frag_free(c->act[i]);
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
            bs_pairer_flush(c->pairer, pos);
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
    bs_pairer_flush(c->pairer, HTS_POS_MAX);
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

    bs_prog = "mhapasm";
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
    c.pairer = bs_pairer_init(o.max_frag, frag_done, &c);

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
    free(c.dks.s);
    bs_pairer_destroy(c.pairer);
    free(c.act);
    free(c.cpg.pos);
    free(c.cpg.mask);
    bs_buf_free(&c.buf);
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
