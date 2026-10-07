#!/usr/bin/env Rscript
# Check mhapmscore -R against mHapDMR, region by region:
#   - mHapDMR:::.mscore_streaming on the records overlapping the region, with
#     all CpGs of the contig: the definition of the region statistics;
#   - mHapDMR::mscore_region_stats(), the function mHapDMR applies to given
#     intervals, with a CpG margin covering the longest read: must be equal.
#     With its default margin (150 bp), mHapDMR before 0.1.1 dropped reads
#     with CpGs farther than that from the region (they failed its CpG-count
#     check); the regions where the default differs are reported, not counted
#     as errors (none from mHapDMR 0.1.1 on).
#
# Usage: Rscript test/check_mscore_regions.R <mhapmscore> <mhap.gz> <cpg.gz> <regions.bed>
#   (mhap.gz and cpg.gz bgzipped and tabix-indexed)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  stopifnot(requireNamespace("mHapDMR", quietly = TRUE))
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: check_mscore_regions.R <mhapmscore> <mhap.gz> <cpg.gz> <regions.bed>")
exe <- args[1]; mhap <- args[2]; cpg <- args[3]; bed <- args[4]
stat_cols <- c("reads", "Nsum", "N2sum", "Sjd", "mscore", "kappa", "Y_prime")

out_fn <- tempfile("regions", fileext = ".tsv")
t0 <- proc.time()[["elapsed"]]
if (system2(exe, c("-i", mhap, "-c", cpg, "-R", bed, "-o", out_fn), stdout = FALSE, stderr = FALSE) != 0)
  stop("mhapmscore failed")
t_c <- proc.time()[["elapsed"]] - t0
out <- fread(out_fn, skip = 1, header = FALSE, col.names = c("chr", "start", "end", stat_cols),
             colClasses = list(character = 1))

# the BED lines mhapmscore reads as regions, in order
lines <- readLines(bed)
f <- strsplit(lines, "\t", fixed = TRUE)
is_reg <- vapply(f, function(x) length(x) >= 3 && grepl("^[0-9]+$", x[2]) && grepl("^[0-9]+$", x[3]) &&
                   as.numeric(x[3]) > as.numeric(x[2]), TRUE)
reg <- data.table(chr = vapply(f[is_reg], `[`, "", 1),
                  beg = as.integer(vapply(f[is_reg], `[`, "", 2)),
                  end = as.integer(vapply(f[is_reg], `[`, "", 3)))
stopifnot(nrow(reg) == nrow(out), reg$chr == out$chr, reg$beg == out$start, reg$end == out$end)

rd <- fread(cmd = paste("gzip -dc", shQuote(mhap)), header = FALSE, sep = "\t", fill = TRUE,
            colClasses = list(character = c(1, 4)))
setnames(rd, 1:5, c("chr", "start", "end", "hap", "count"))
cpg_cache <- list()
contig_cpgs <- function(ch) {
  if (is.null(cpg_cache[[ch]])) {
    l <- suppressWarnings(system2("tabix", c(cpg, ch), stdout = TRUE, stderr = FALSE))
    cpg_cache[[ch]] <<- sort(as.integer(vapply(strsplit(l, "\t", fixed = TRUE), `[`, "", 2)))
  }
  cpg_cache[[ch]]
}
same <- function(a, b) isTRUE(all.equal(unname(as.numeric(unlist(a))), unname(as.numeric(unlist(b))), tolerance = 1e-5))

# 1. the definition, region by region
n_bad <- 0L
for (i in seq_len(nrow(reg))) {
  g <- reg[i]
  r <- rd[chr == g$chr & start <= g$end & end >= g$beg + 1L]
  raw <- sprintf("%s\t%d\t%d\t%s\t%d", r$chr, r$start, r$end, r$hap, r$count)
  expect <- mHapDMR:::.mscore_streaming(raw, contig_cpgs(g$chr), g$beg + 1L, g$end)
  if (!same(out[i, ..stat_cols], expect)) {
    n_bad <- n_bad + 1L
    cat("[streaming]", g$chr, g$beg, g$end, "mhapmscore:", unlist(out[i, ..stat_cols]), " mHapDMR:", expect, "\n")
  }
}

# 2. mscore_region_stats, with a margin covering the longest read and with the default
gr <- GRanges(reg$chr, IRanges(reg$beg + 1L, reg$end))
span <- max(rd$end - rd$start) + 1L
t0 <- proc.time()[["elapsed"]]
full <- as.data.table(mcols(mHapDMR::mscore_region_stats(mhap, cpg, gr, margin = span, verbose = FALSE)))
t_r <- proc.time()[["elapsed"]] - t0
for (i in seq_len(nrow(reg))) {
  if (!same(out[i, ..stat_cols], full[i, ..stat_cols])) {
    n_bad <- n_bad + 1L
    cat("[region_stats]", reg$chr[i], reg$beg[i], reg$end[i], "mhapmscore:", unlist(out[i, ..stat_cols]),
        " mHapDMR:", unlist(full[i, ..stat_cols]), "\n")
  }
}
def <- as.data.table(mcols(mHapDMR::mscore_region_stats(mhap, cpg, gr, verbose = FALSE)))
diff_def <- which(!vapply(seq_len(nrow(reg)), function(i) same(out[i, ..stat_cols], def[i, ..stat_cols]), TRUE))
lost <- sum(out$reads[diff_def] - def$reads[diff_def], na.rm = TRUE)

cat(sprintf("%d regions (%d with reads): %d mismatches with .mscore_streaming and with mscore_region_stats(margin = %d)\n",
            nrow(reg), sum(out$reads > 0), n_bad, span))
cat(sprintf("mscore_region_stats with the default margin (150): %d regions differ, %g of %g reads there dropped\n",
            length(diff_def), lost, sum(out$reads[diff_def])))
cat(sprintf("time: mhapmscore %.2f s, mscore_region_stats %.2f s\n", t_c, t_r))
quit(status = if (n_bad) 1L else 0L)
