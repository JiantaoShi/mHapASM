# mHapASM

Read-level methylation tools for bisulfite sequencing BAMs, written in C on a
shared core (`src/bsread.c`) that lays reads out on the reference with their
CIGAR, takes the bisulfite strand from the aligner's tags, pairs mates into
fragments, and calls CpGs with the rules of
[wgbs_tools](https://github.com/nloyfer/wgbs_tools).

* **`mhapasm`**: allele-specific M-score statistics. For every SNP in a list,
  splits the fragments by allele and reports, per allele, the sufficient
  statistics used by [mHapDMR](https://github.com/JiantaoShi/mHapDMR) (DSS)
  and [MscoreDMR](https://github.com/JiantaoShi/MscoreDMR) (GLS).
* **`mhapconvert`**: BAM to [mHap](https://jiantaoshi.github.io/mHap/), a
  rewrite of `mhaptools convert` ([mHapTools](https://github.com/butyuhao/mHapTools))
  that handles soft clips and indels, the strand tags of Bismark, bwa-meth,
  BISCUIT and BSMAP, and non-directional libraries.
* **`mhapmscore`**: M-score statistics from an mHap file, for every CpG (the
  statistics of `mhapasm` over the reads that span it) or for given regions,
  as mHapDMR computes them.

## Build

Requires a C compiler and htslib (>= 1.10).

```bash
make                          # htslib found in $CONDA_PREFIX, Homebrew, /usr/local or /usr
make HTSLIB=/opt/htslib-1.20  # any prefix holding include/htslib and lib/libhts.*
```

Contig names may differ between the BAM, the CpG file and the SNP file by a
`chr` prefix (`1` / `chr1`) and in the name of the mitochondrial genome
(`chrM` / `MT`). The CpG file is the mHap annotation file (chr, 1-based C
position, ...), e.g. `hg38_CpG.gz`.

## mhapasm

```bash
mhapasm -b sample.bam -c hg38_CpG.gz -s snps.vcf -o sample.asm.tsv --qc sample.qc.tsv
```

| Option | Default | Meaning |
|---|---|---|
| `-b FILE` | | Coordinate-sorted, indexed BAM/CRAM (single- or paired-end). |
| `-c FILE` | | CpG file, bgzipped + tabix-indexed. |
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

SNPs are processed and reported in BAM contig order and by position; SNPs on
contigs absent from the BAM are reported last with zero counts.

### Output

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

### Allele calls

Base at the SNP of the CIGAR-aligned read (`snp_patter` of wgbs_tools). C/T
SNPs cannot be called on the top strand, G/A SNPs not on the bottom strand. On
the top strand a C allele also accepts T (unless the other allele is T); on the
bottom strand a G allele also accepts A (unless the other allele is A). One
informative mate decides; mates supporting different alleles discard the
fragment. CpGs whose C or G is a listed SNP are ignored for both alleles
(`--no-mask` keeps them), so the alleles are compared on the same CpGs.

### Comparing the alleles in R

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

## mhapconvert

```bash
mhapconvert -i sample.bam -c hg38_CpG.gz -o sample.mhap.gz
```

The short options of `mhaptools convert` (`-i -c -r -b -n -o -m`) are kept.

| Option | Default | Meaning |
|---|---|---|
| `-i FILE` | | Coordinate-sorted SAM/BAM/CRAM (`-` for stdin); `-r` and `-b` need an index. |
| `-c FILE` | | CpG file, plain or (b)gzipped; read per contig if tabix-indexed, otherwise loaded at once. |
| `-r STR` | | Only reads overlapping this region (`chr:beg-end`). |
| `-b FILE` | | Only reads overlapping these BED regions (overlapping regions are merged; every read is used once). |
| `-o FILE` | out.mhap.gz | Output; bgzipped and tabix-indexed (columns 1-3, 1-based) if the name ends in `.gz`, plain text otherwise. |
| `-n` | | Non-directional: strand `*`, haplotypes not grouped by strand. |
| `-m BS\|TAPS` | BS | BS (bisulfite, EM-seq): C = methylated; TAPS: T = methylated. |
| `-q INT` | 10 | Minimum mapping quality. |
| `-F INT` | 0xF04 | Skip reads with any of these flags (unmapped, secondary, QC-fail, duplicate, supplementary). |
| `-B INT` | 0 | Minimum base quality of a CpG call. |
| `-L INT` | 1000 | Maximum distance between mates of one fragment (bp). |
| `--max-unconv FLOAT` | 0.2 | Drop incompletely converted reads: at least 3 C outside CpGs, making up more than FLOAT of the C and T bases outside CpGs (G and A on the bottom strand); -1 = no filter. Not applied with `-m TAPS`. |
| `--max-ch INT` | -1 | Drop reads with more than INT unconverted cytosines outside CpGs, as Bismark's `filter_non_conversion`; -1 = no filter. |
| `--split` | | Cut a fragment at a sequenced CpG without a call instead of dropping the fragment. |
| `--qname` | | One line per record with the read name in column 7 (no collapsing). |
| `--no-index` | | Do not write the tabix index. |
| `-T FILE` | | Reference FASTA (CRAM only). |
| `-@ INT` | 0 | Extra threads for decompression and compression. |

Output: mHap records `chr start end haplotype count strand`, with start/end the
1-based positions of the first/last CpG and one `0`/`1` per consecutive CpG,
sorted by position and with identical records collapsed into `count`.

An mHap haplotype cannot hold a missing value, so a fragment's calls are cut
where a CpG has no call. CpGs between two mates that do not overlap were not
sequenced, and the fragment is written as two records (as in mHapTools). A
sequenced CpG without a call (sequencing error, CpG context not seen on the
read, deletion, mates disagreeing) drops the fragment by default, as in
mHapTools; cutting it instead (`--split`) keeps the data but turns one molecule
into several records, which lowers read-level statistics such as the M-score
that score a read as methylated when any of its CpGs is.

The conversion filter needs no reference. Outside CpGs, a C on a top-strand
read is an unconverted cytosine and a T is a converted cytosine or a genomic T,
so C makes up almost none of the C and T bases of a converted read and about
40% of those of an unconverted read (human genome). A count of unconverted
cytosines (`--max-ch`) also counts the bases filled in at fragment ends during
library preparation and genuine non-CpG methylation: on the BISCUIT sample
below, `--max-ch 2` dropped 37% of the reads and `--max-unconv 0.2` 0.5%.

### Changes from `mhaptools convert`

| | `mhaptools convert` | `mhapconvert` |
|---|---|---|
| Soft clips, indels | CIGAR ignored: bases taken by offset from the start of SEQ | Reads laid out on the reference with their CIGAR |
| Strand | Flags only; paired reads must be proper pairs | Aligner tag (Bismark `XG`, bwa-meth/BISCUIT `YD`, BSMAP `ZS`), else flags |
| Non-directional (CTOT/CTOB) | Read from the wrong strand | Correct with strand tags |
| CpG call | Base at the C (top) or G (bottom) only; any other base drops the read | Both bases must be on the read and show the CpG context (wgbs_tools `patter`) |
| Mates | Merged if their CpG spans overlap, higher base quality wins | Merged into one fragment; a CpG on which they disagree has no call |
| Incomplete conversion | Bismark `XM` only: any `X`/`H`/`U` drops the read | Share of unconverted C outside CpGs (`--max-unconv 0.2`), any aligner |
| MAPQ filter | None | `-q 10` |
| Contig names | Must match the CpG file exactly; otherwise an empty result without error | `chr` prefix and `chrM`/`MT` matched |
| Output | Whole contigs in memory, temporary files, no index | Streamed, bgzipped, tabix-indexed |

### Accuracy on simulated data

30,000 fragments per data set (`test/simulate_pe.py`, 0.5% sequencing errors),
CpG calls scored against the simulated methylation states (default options):

| Data | Fragments | mhaptools | mhapconvert |
|---|---|---|---|
| paired-end, directional | plain | 99.69% | 99.79% |
| | with indel / soft clip | 97.74% | 99.80% |
| single-end, directional | plain | 99.69% | 99.80% |
| | with indel / soft clip (fragments written) | 92.73% (1,604 of 4,521) | 99.80% (4,039 of 4,521) |
| paired-end, 30% non-directional, strand tags | CTOT/CTOB | 40.77% | 99.79% |
| paired-end, 5% unconverted fragments | all | 96.43% | 99.76% (96.90% with `--max-unconv -1`) |

Results are the same with the tags of Bismark (`XG:Z`), bwa-meth (`YD:Z`),
BISCUIT (`YD:A`) and BSMAP (`ZS:Z`). Without strand tags, CTOT/CTOB reads are
read on the wrong strand, where they look unconverted or fail the CpG context
check: 99.1% of these fragments are dropped rather than miscalled. To
reproduce a row (`test/mhaptools_emulate.py` reproduces `mhaptools convert`
read by read):

```bash
python3 test/simulate_pe.py --out-dir nd --length 100000 --n-frag 30000 --n-snp 100 --seed 7 \
    --non-directional 0.3 --aligner bismark --meth-truth nd/meth_truth.tsv
./mhapconvert -i nd/sim.bam -c nd/cpg.gz --qname -o nd/cvt.mhap.gz
python3 test/mhaptools_emulate.py --bam nd/sim.bam --cpg nd/cpg.gz > nd/mht.rec
python3 test/score_truth.py --records nd/cvt.mhap.gz --truth-dir nd --label mhapconvert
python3 test/score_truth.py --records nd/mht.rec --truth-dir nd --label mhaptools
```

### Real data

A BISCUIT paired-end sample (2.0M reads, hg38), fragment by fragment
(553,564 fragments with records in either tool; the reason for every
difference is taken from the reads):

| | Fragments |
|---|---|
| Identical records | 58.1% |
| A read has an indel or soft clip (misplaced by mhaptools) | 21.7% |
| MAPQ < 10 (kept by mhaptools) | 9.0% |
| CpG at a read end with only one of its bases on the read (called by mhaptools) | 4.3% |
| Mates without a CpG between them: one record instead of two | 3.1% |
| A sequenced CpG without a call: fragment dropped by mhapconvert | 2.0% |
| Combinations of the above | 1.2% |
| Improper pairs (dropped by mhaptools) | 0.4% |
| A non-C/T base at a CpG: read dropped by mhaptools | 0.2% |
| An incompletely converted read (dropped by mhapconvert) | 0.1% |

With the CpG file named as in the BAM (`1`, `2`, ...); with `chr1`-style names
`mhaptools convert` finds no CpG and writes nothing, without an error. The
table comes from `test/compare_mhaptools.py BAM mhaptools.rec
mhapconvert.qname.mhap.gz mhapconvert.nofilter.qname.mhap.gz` (the last made
with `--max-unconv -1`).

## mhapmscore

```bash
mhapmscore -i sample.mhap.gz -c hg38_CpG.gz -o sample.mscore.tsv.gz        # per CpG
mhapmscore -i sample.mhap.gz -c hg38_CpG.gz -R dmrs.bed -o sample.dmrs.tsv  # per region
```

M-score statistics from an mHap file (from `mhapconvert`, `mhaptools convert`
or any other source), per CpG or per region. Every mHap record is a read
covering consecutive CpGs. A read *t*, seen *c_t* times (the count column),
contributes *N_t* CpGs and *Z_t* = 1 if any of them is methylated, and the
reads are summarised with the statistics of `mhapasm`:

* **per CpG** (default): the reads whose haplotype spans the CpG, with *N_t*
  over all CpGs of the read, or with `-w` over those within INT bp of the CpG;
* **per region** (`-R`): the reads with CpGs in the region, with *N_t* over
  those CpGs, as in mHapDMR's `mscore_region_stats`.

| Option | Default | Meaning |
|---|---|---|
| `-i FILE` | | mHap file, plain or bgzipped. The records of a contig must be together, in any order; `-r` and `-b` need a tabix index. |
| `-c FILE` | | CpG file used to make the mHap file, bgzipped and tabix-indexed. |
| `-r STR` | | Only CpGs in this region (`chr:beg-end`). |
| `-b FILE` | | Only CpGs in these BED regions. |
| `-w INT` | 0 | *N_t* and *Z_t* from the CpGs of the read within INT bp of the CpG; 0 = all CpGs of the read. |
| `-m INT` | 1 | Report CpGs spanned by at least INT reads. |
| `-R FILE` | | Statistics of each region of this BED file instead of each CpG; needs a tabix index of the mHap file. |
| `-o FILE` | stdout | Output; bgzipped and tabix-indexed (`-s1 -b2 -e2`) if the name ends in `.gz`. |
| `--no-index` | | Do not write the tabix index. |
| `-@ INT` | 0 | Extra threads for decompression and compression. |

Output, one line per CpG (1-based position) after a `#` header line:

```
#chr    pos      reads  Nsum  N2sum  Sjd  mscore    kappa    Y_prime
chr11   2018724  4      15    57     11   0.733333  3.94737  2.89474
chr11   2018741  11     38    136    28   0.736842  10.6176  7.82353
```

`reads` = Σ *c_t*, `Nsum` = Σ *c_t N_t*, `N2sum` = Σ *c_t N_t*², `Sjd` =
Σ *c_t N_t Z_t*, `mscore` = `Sjd / Nsum`, `kappa` = `Nsum² / N2sum` and
`Y_prime` = `Sjd * Nsum / N2sum`, as in the `mhapasm` table above. With a
window, the statistics of a CpG are those mHapDMR computes for the region
[pos − w, pos + w] from the reads that span the CpG.

With `-R`, every BED line gets one output line, in input order: `chr start
end` as in the BED file, then the same seven columns. A region without reads
has zeros and `NA`; lines that are not regions (a column header, for
instance) are skipped with a message. A `.gz` output is tabix-indexed only if
the BED file is sorted by position.

```
#chr    start   end     reads  Nsum   N2sum   Sjd   mscore    kappa    Y_prime
chrS    3000    5000    1249   10278  112236  7523  0.731952  941.207  688.918
chrS    3014    3015    0      0      0       0     NA        NA       NA
```

A record is used only if the CpG file has as many CpGs in [start, end] as its
haplotype has states (the rule of mHapDMR and MscoreDMR). Other records are
counted and skipped, with a warning when they exceed 1% (a different CpG file
or genome build).

mHapDMR before version 0.1.1 looked CpGs up only within `margin` bp of the
region (150 by default), so `mscore_region_stats()` dropped every read with a
CpG farther outside, which then failed this check: 25 of 150 ESCC DMRs (reads
up to 233 bp) lost 2 to 4 reads each (M-score change up to 0.02), and 114 of
126 regions of the simulated paired-end data (records up to 340 bp) lost 1.3%
of their reads. mHapDMR 0.1.1 looks CpGs up over the full span of the reads,
like `mhapmscore`, and gives the same results (see Validation).

## How reads are processed

Shared by both tools (`src/bsread.c`):

1. **Filters**: flags `-F`, MAPQ `-q`; `mhapconvert` also drops incompletely
   converted reads (`--max-unconv`).
2. **Fragments**: mates are paired by QNAME (as `match_maker`). A read whose
   mate is unmapped, on another contig, farther than `-L`, or filtered out is a
   fragment on its own.
3. **Strand**: the aligner's strand tag when present (Bismark `XG`, bwa-meth /
   BISCUIT `YD`, BSMAP `ZS`); otherwise the `is_bottom` flag rule of wgbs_tools
   (bottom = read 1 reverse or read 2 forward), without requiring the
   proper-pair bit.
4. **Layout**: the read is placed on the reference with its CIGAR (`clean_CIGAR`):
   soft clips and insertions are skipped, deletions leave a gap.
5. **Methylation** (`patter`): CpGs from the CpG file. A CpG is called only if
   its C and G both lie in the aligned part of the read and the read shows the
   CpG context itself (top strand: C/T followed by G; bottom strand: C followed
   by G/A). The calls of the two mates are merged; a CpG on which they disagree
   is dropped (`merge_PE`).

Differences from the wgbs_tools pipeline (`samtools view | match_maker |
snp_patter` / `patter`): the strand tags, the proper-pair bit, and supplementary
alignments being excluded by default. For directional libraries with proper
pairs the results are identical to wgbs_tools.

## Validation

* `mhapasm`: read by read, the allele assignment and the CpG haplotypes are
  identical to the upstream wgbs_tools C++ chain (commit 6f24bed, April 2026)
  on the IGF2 tutorial BAM (481 SNPs: every ordered base pair at 40 random
  positions), a simulated paired-end BAM with indels, soft clips, low MAPQ and
  missing mates, and chr1 of a real BISCUIT paired-end sample (463 SNPs; proper
  pairs, YD tags removed so that both tools use the flags). The allele
  assignment is also identical on the BISCUIT hg38 IGF2 BAM. The statistics
  equal those of mHapDMR's `.mscore_accumulate` / `.mscore_streaming`
  recomputed from the dumped haplotypes (windows 0 to 150 bp).
* `mhapconvert`: with `--split`, the records of every fragment are exactly the
  call runs of wgbs_tools `patter --long` on the IGF2 BAM (2,576 fragments),
  simulated paired-end BAMs (7,799 and 195,408 fragments) and chr1 of the real
  BISCUIT sample (39,964 fragments); by default the same fragments are written
  identically and only fragments with such a gap are dropped. A per-read
  emulation of `mhaptools convert` used for the comparisons above reproduces
  its output exactly on all these data sets.
* `mhapmscore`: `test/check_mscore.R` recomputes every CpG in R from (read,
  CpG) pairs and passes the records spanning each CpG to mHapDMR's
  `.mscore_streaming` with the region [pos − w, pos + w] (the whole read for
  w = 0). The results are identical on the IGF2 and simulated mHap files
  (windows 0 to 500, also with records that do not match the CpG file) and on
  chr21 of an ESCC WGBS sample (1.23M records; all 373,647 CpGs, 2,000 of them
  also through mHapDMR; windows 0 and 150). With `-r` and `-b` the output
  equals the corresponding rows of the whole-file output, and `-@` does not
  change it. With `-R`, every region equals mHapDMR's `.mscore_streaming`
  (with all CpGs of the contig) and `mscore_region_stats()` (mHapDMR >= 0.1.1,
  or older versions with a margin as long as the longest read;
  `test/check_mscore_regions.R`): the edge cases of
  `test/regions_mscore.bed` (overlapping, unsorted, no CpG, absent contig)
  and the 150 DMRs of the ESCC samples. The region modes load CpGs in windows
  of at least 1 Mb, widened for reads that stick out; a build with 64-bp
  windows gives the same results (`make test`, and on the ESCC DMRs).
* Compare with a current wgbs_tools build. In older versions `snp_patter`
  writes out the next read in place of a paired-end read without its mate
  (fixed in 4de9ef5, 2024-11-17) and counts a base matching neither allele as
  the second allele when combining mates (fixed in a0524f4, 2025-03-12), and
  `patter` has no CpG context check (added in 6a018db, 2026-01-22).

```bash
make test                  # golden outputs, region modes, threads, index, mHapDMR checks (if installed)
WGBS_BIN=<dir> make test   # + read-by-read comparison with wgbs_tools (match_maker, snp_patter, patter)
```

Test data: `test/data` holds the wgbs_tools tutorial BAM (single-end, hg19,
IGF2 locus with the heterozygous SNP chr11:2019496 C/A) and the hg19 CpGs of
chr11:2,000,000-2,050,000; `test/sim` holds a simulated paired-end BAM made by
`test/simulate_pe.py` (phased SNPs, allele-specific methylation, indels, soft
clips, missing mates) and its truth table.

## Performance

| Tool | Data | Time | Memory |
|---|---|---|---|
| `mhapasm` | simulated PE, 394k reads, 2,006 SNPs | 0.24 s | 8 MB |
| `mhapasm` | BISCUIT PE sample, 2.0M reads, 31.1M SNPs (dbSNP 146 common VCF, 6.4 GB) | 38 s | 1.9 GB |
| `mhapconvert` | BISCUIT PE sample, 2.0M reads, hg38 | 13.9 s (10.0 s with `-@ 4`) | 61 MB |
| `mhaptools convert` 0.10 | same | 30.9 s | 293 MB |
| `mhapmscore` | ESCC WGBS sample, 88.5M mHap records, hg19 (27.4M CpGs written) | 2 min 41 s (1 min 6 s with `-@ 4`; `-w 150`: 2 min 43 s) | 126 MB |
| `mhapmscore -R` | same sample, 150 DMRs / 100k / 1M regions of 1 kb | 0.4 s / 20 s / 1 min 14 s | 52 MB / 144 MB (100k / 1M) |
| mHapDMR `mscore_region_stats()` | 150 DMRs, mHap records of the DMRs only | 0.5 s (`mhapmscore -R`: 0.02 s) | |

## Acknowledgements

The read-level rules (allele call, methylation call, mate merging) follow
`snp_patter`, `patter` and `match_maker` of
[wgbs_tools](https://github.com/nloyfer/wgbs_tools) (Loyfer *et al.*,
*Life Science Alliance* 9(4): e202503514, 2026). The tutorial BAM in
`test/data` comes from wgbs_tools; for these data cite Loyfer *et al.*,
[A DNA methylation atlas of normal human cell types](https://doi.org/10.1038/s41586-022-05580-6),
*Nature* 613: 355–364 (2023). Aligner strand tags are handled as in
[mhapx](https://github.com/BestRivenWR/mhapx). The M-score statistics are
those of [mHapDMR](https://github.com/JiantaoShi/mHapDMR) and
[MscoreDMR](https://github.com/JiantaoShi/MscoreDMR); `mhapconvert` replaces
`convert` of [mHapTools](https://github.com/butyuhao/mHapTools).

## License

MIT, see [LICENSE](LICENSE).
