# mhapasm

Allele-specific M-score statistics from bisulfite sequencing BAMs.

For every SNP in a list, `mhapasm` splits the reads (fragments) by allele and
computes, per allele, the statistics needed to compare M-scores: the sufficient
statistics used by [mHapDMR](https://github.com/JiantaoShi/mHapDMR) (DSS) and
[MscoreDMR](https://github.com/JiantaoShi/MscoreDMR) (GLS). It reads the BAM
directly with htslib; there is no intermediate allele BAM, mHap file or Java
step. The read-level rules (allele call, methylation call, mate merging) are
those of [wgbs_tools](https://github.com/nloyfer/wgbs_tools) `snp_patter` /
`patter`, re-implemented in C.

## Build

Requires a C compiler and htslib (>= 1.10).

```bash
make                          # htslib found in $CONDA_PREFIX, Homebrew, /usr/local or /usr
make HTSLIB=/opt/htslib-1.20  # any prefix holding include/htslib and lib/libhts.*
```

## Usage

```bash
mhapasm -b sample.bam -c hg38_CpG.gz -s snps.vcf -o sample.asm.tsv --qc sample.qc.tsv
```

| Option | Default | Meaning |
|---|---|---|
| `-b FILE` | | Coordinate-sorted, indexed BAM/CRAM (single- or paired-end, directional library). |
| `-c FILE` | | CpG file, bgzipped + tabix-indexed; column 1 = chr, column 2 = 1-based C position (the mHap annotation file, e.g. `hg38_CpG.gz`). |
| `-s FILE` | | SNPs: `chr pos ref alt` (1-based) or VCF (header optional), plain or gzipped. Multi-allelic records are split into one SNP per ALT; indels are skipped. |
| `-o FILE` | stdout | Per-allele statistics. |
| `--qc FILE` | | Per-SNP allele-calling counts (all SNPs). |
| `--dump FILE` | | Per-fragment CpG haplotypes assigned to each allele (for inspection). |
| `-q INT` | 10 | Minimum mapping quality. |
| `-F INT` | 0xF04 | Skip reads with any of these flags (unmapped, secondary, QC-fail, duplicate, supplementary). |
| `-Q INT` | 0 | Minimum base quality at the SNP. |
| `-B INT` | 0 | Minimum base quality of a CpG call. |
| `-w INT` | 0 | Use only CpGs within INT bp of the SNP; 0 = all CpGs of the fragment. |
| `-L INT` | 1000 | Maximum distance between mates of one fragment (bp). |
| `--no-mask` | | Keep CpGs that overlap a listed SNP (masked by default). |
| `--min-allele INT` | 0 | Report a SNP only if each allele is carried by >= INT fragments. |
| `-T FILE` | | Reference FASTA (CRAM only). |
| `-@ INT` | 0 | Extra BAM decompression threads. |

Contig names may differ between the BAM, the CpG file and the SNP file by a
`chr` prefix (`1` / `chr1`) and in the name of the mitochondrial genome
(`chrM` / `MT`). SNPs are processed and reported in BAM contig order and by
position; SNPs on contigs absent from the BAM are reported last with zero counts.

## Output

Two lines per SNP, REF first:

```
#chr    pos      allele  reads  Nsum  N2sum  Sjd  mscore  kappa    Y_prime
chr11   2019496  C       20     74    284    0    0       19.2817  0
chr11   2019496  A       21     89    389    89   1       20.3625  20.3625
```

For a fragment *t* assigned to the allele, *N_t* is the number of CpGs called
on it (within the window if `-w` is set) and *Z_t* = 1 if any of them is
methylated. Fragments with *N_t* = 0 are not counted.

| Column | Definition | mHapDMR | MscoreDMR GLS |
|---|---|---|---|
| `reads` | number of fragments | `reads` | |
| `Nsum` | sum of *N_t* | `Nsum` (*T_jd*) | `T_matrix` |
| `N2sum` | sum of *N_t*^2 | `N2sum` | `N_sq_matrix` |
| `Sjd` | sum of *N_t Z_t* | `Sjd` (*S_jd*) | `S_matrix` |
| `mscore` | `Sjd / Nsum` | `mscore` | |
| `kappa` | `Nsum^2 / N2sum` | `kappa` (DSS coverage) | effective coverage |
| `Y_prime` | `Sjd * Nsum / N2sum` | `Y_prime` (DSS methylated count) | |

`--qc` columns: `n_cover` (fragments with a mate on the SNP), `n_ref`, `n_alt`
(fragments assigned to each allele, with or without CpGs), `n_conflict` (mates
support different alleles), `n_bs_ambiguous` (C/T SNP on the top strand or G/A
SNP on the bottom strand, where bisulfite conversion hides the allele) and
`n_other` (base matching neither allele, deletion or low quality).

## How reads are processed

1. **Filters**: flags `-F`, MAPQ `-q`.
2. **Fragments**: mates are paired by QNAME (as `match_maker`). A read whose
   mate is unmapped, on another contig, farther than `-L`, or filtered out is a
   fragment on its own.
3. **Strand**: the aligner's strand tag when present (Bismark `XG`, bwa-meth /
   BISCUIT `YD`, BSMAP `ZS`); otherwise the `is_bottom` flag rule of wgbs_tools
   (bottom = read 1 reverse or read 2 forward), without requiring the
   proper-pair bit.
4. **Allele** (`snp_patter`): base at the SNP of the CIGAR-aligned read. C/T
   SNPs cannot be called on the top strand, G/A SNPs not on the bottom strand.
   On the top strand a C allele also accepts T (unless the other allele is T);
   on the bottom strand a G allele also accepts A (unless the other allele is A).
   One informative mate decides; mates supporting different alleles discard
   the fragment.
5. **Methylation** (`patter`): CpGs from the CpG file. A CpG is called only if
   its C and G both lie in the aligned part of the read and the read shows the
   CpG context itself (top strand: C/T followed by G; bottom strand: C followed
   by G/A). The calls of the two mates are merged; a CpG on which they disagree
   is dropped (`merge_PE`).
6. **Masking**: CpGs whose C or G is a listed SNP are ignored for both alleles,
   so the alleles are compared on the same CpGs.

Differences from the wgbs_tools pipeline (`samtools view | match_maker |
snp_patter`): the strand tags, the proper-pair bit, and supplementary
alignments being excluded by default. For directional libraries with proper
pairs the results are identical to wgbs_tools. Compared with mHapSuite /
mHapTools conversion, reads with soft clips or indels are laid out with their
CIGAR (mHapTools ignores it), and the CpG context check above is applied.

## Validation

* Read by read, the allele assignment and the CpG haplotypes are identical to
  the upstream wgbs_tools C++ chain (commit 6f24bed, April 2026) on the IGF2
  tutorial BAM (481 SNPs: every ordered base pair at 40 random positions), a
  simulated paired-end BAM with indels, soft clips, low MAPQ and missing mates
  (`test/simulate_pe.py`), and chr1 of a real BISCUIT paired-end sample (463
  SNPs; proper pairs, YD tags removed so that both tools use the flags). The
  allele assignment is also identical on the BISCUIT hg38 IGF2 BAM.
* The statistics equal those of mHapDMR's `.mscore_accumulate` /
  `.mscore_streaming` recomputed from the dumped haplotypes (windows 0 to 150 bp).
* Compare with a current wgbs_tools build. In older versions `snp_patter`
  writes out the next read in place of a paired-end read without its mate
  (fixed in 4de9ef5, 2024-11-17) and counts a base matching neither allele as
  the second allele when combining mates (fixed in a0524f4, 2025-03-12), and
  `patter` has no CpG context check (added in 6a018db, 2026-01-22).

```bash
make test                  # golden outputs (IGF2 tutorial BAM, simulated PE BAM), threads
WGBS_BIN=<dir> make test   # + read-by-read comparison with wgbs_tools (match_maker, snp_patter, patter)
```

Test data: `test/data` holds the wgbs_tools tutorial BAM (single-end, hg19,
IGF2 locus with the heterozygous SNP chr11:2019496 C/A) and the hg19 CpGs of
chr11:2,000,000-2,050,000; `test/sim` holds a simulated paired-end BAM made by
`test/simulate_pe.py` (phased SNPs, allele-specific methylation, indels, soft
clips, missing mates) and its truth table.

## Performance

| Data | SNPs | Time | Memory |
|---|---|---|---|
| simulated PE, 394k reads | 2,006 | 0.24 s | 8 MB |
| BISCUIT PE sample, 2.0M reads | 31.1M (dbSNP 146 common VCF, 6.4 GB) | 38 s | 1.9 GB |

Run time is dominated by reading the BAM and the SNP list; memory by the SNP
list (about 32 bytes per SNP).

## Comparing the alleles in R

```r
# S / T / N_sq matrices (rows = SNPs, columns = sample.allele) from mhapasm outputs
read_mhapasm <- function(files) {
  cols <- c("chr", "pos", "allele", "reads", "Nsum", "N2sum", "Sjd", "mscore", "kappa", "Y_prime")
  tabs <- lapply(files, function(f) {
    x <- read.table(f, sep = "\t", comment.char = "#", col.names = cols,
                    colClasses = c("character", "integer", "character", rep("numeric", 7)))
    r <- x[c(TRUE, FALSE), ]; a <- x[c(FALSE, TRUE), ]   # REF line, then ALT line
    list(id = paste(r$chr, r$pos, r$allele, a$allele, sep = ":"), ref = r, alt = a)
  })
  ids <- Reduce(intersect, lapply(tabs, `[[`, "id"))
  pick <- function(stat) {
    m <- do.call(cbind, lapply(tabs, function(t) {
      i <- match(ids, t$id); cbind(t$ref[[stat]][i], t$alt[[stat]][i])
    }))
    dimnames(m) <- list(ids, paste(rep(names(files), each = 2), c("ref", "alt"), sep = "."))
    m
  }
  list(S = pick("Sjd"), T = pick("Nsum"), N_sq = pick("N2sum"), reads = pick("reads"))
}

# One sample: Wald z of the allele difference on the arcsine scale of the
# MscoreDMR GLS, Var(asin(2p - 1)) ~ 1 / kappa (no overdispersion term).
asm_z <- function(file) {
  m <- read_mhapasm(c(x = file))
  ok <- m$T[, 1] > 0 & m$T[, 2] > 0
  z <- function(S, T) asin(2 * (S + 0.5) / (T + 1) - 1)
  d <- z(m$S[ok, 2], m$T[ok, 2]) - z(m$S[ok, 1], m$T[ok, 1])
  se <- sqrt(m$N_sq[ok, 1] / m$T[ok, 1]^2 + m$N_sq[ok, 2] / m$T[ok, 2]^2)
  data.frame(snp = rownames(m$S)[ok], mscore_ref = m$S[ok, 1] / m$T[ok, 1],
             mscore_alt = m$S[ok, 2] / m$T[ok, 2], z = d / se,
             p = 2 * pnorm(-abs(d / se)), row.names = NULL)
}

# Several samples: alleles paired within sample, MscoreDMR two-stage GLS
library(MscoreDMR)
m <- read_mhapasm(c(s1 = "s1.asm.tsv", s2 = "s2.asm.tsv", s3 = "s3.asm.tsv"))
keep <- rowSums(m$T > 0 & m$N_sq > 0) == ncol(m$T)     # coverage in every column
info <- data.frame(sample = factor(sub("\\..*$", "", colnames(m$S))),
                   allele = factor(sub("^.*\\.", "", colnames(m$S)), levels = c("ref", "alt")))
X <- model.matrix(~ sample + allele, info)
fit <- fit_mscore_gls(m$S[keep, ], m$T[keep, ], m$N_sq[keep, ], X,
                      contrast = as.numeric(colnames(X) == "allelealt"))
res <- data.frame(snp = rownames(m$S)[keep], stat = fit$statistic, p = fit$p_value)
```

`kappa` and `Y_prime` can likewise be used as `Cov` and `M` of a DSS `BSseq`
object (as in `mHapDMR::build_bsseq_mscore`) with a `~ sample + allele` design.

## Acknowledgements

The read-level rules (allele call, methylation call, mate merging) follow
`snp_patter`, `patter` and `match_maker` of
[wgbs_tools](https://github.com/nloyfer/wgbs_tools) (Loyfer *et al.*,
*Life Science Alliance* 9(4): e202503514, 2026), which also provides the
tutorial BAM in `test/data`. Aligner strand tags are handled as in
[mhapx](https://github.com/BestRivenWR/mhapx). The M-score statistics are
those of [mHapDMR](https://github.com/JiantaoShi/mHapDMR) and
[MscoreDMR](https://github.com/JiantaoShi/MscoreDMR).
