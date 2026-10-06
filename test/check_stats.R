#!/usr/bin/env Rscript
# Recompute the per-allele M-score statistics of mhapasm from its per-fragment
# dump and compare them with its main output.
#   - every fragment: N/Z counted in R, summed with mHapDMR:::.mscore_accumulate
#   - fragments without uncalled CpGs: additionally mHapDMR:::.mscore_streaming,
#     i.e. the exact code mHapDMR applies to mHap records of a region
# The window (-w) is applied in R to the full-fragment dump of a -w 0 run, so
# the trimming done in C is checked independently.
#
# Usage: Rscript test/check_stats.R <mhapasm> <bam> <cpg.gz> <snps.txt> <window> [extra mhapasm args]

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: check_stats.R <mhapasm> <bam> <cpg.gz> <snps.txt> <window> [args]")
exe <- args[1]; bam <- args[2]; cpg <- args[3]; snp_fn <- args[4]; w <- as.integer(args[5])
extra <- if (length(args) > 5) args[-(1:5)] else character(0)

tmp <- tempfile("chk"); dir.create(tmp)
dump0 <- file.path(tmp, "dump0.tsv"); out_w <- file.path(tmp, "out_w.tsv")
# mhapasm reports SNPs sorted by position (ties in input order), REF then ALT;
# feed it a stably sorted copy so that output rows 2i-1 / 2i belong to SNP i
snps <- fread(snp_fn, header = FALSE, col.names = c("chr", "pos", "ref", "alt"),
              colClasses = c("character", "integer", "character", "character"))
snps <- snps[order(pos)]
snp_fn <- file.path(tmp, "snps.sorted.txt")
fwrite(snps, snp_fn, sep = "\t", col.names = FALSE)
run <- function(a) {
  st <- system2(exe, a, stdout = FALSE, stderr = FALSE)
  if (st != 0) stop("mhapasm failed: ", paste(a, collapse = " "))
}
run(c("-b", bam, "-c", cpg, "-s", snp_fn, "--no-mask", "-w", "0", "--dump", dump0, "-o", file.path(tmp, "out0.tsv"), extra))
run(c("-b", bam, "-c", cpg, "-s", snp_fn, "--no-mask", "-w", w, "-o", out_w, extra))

cols <- c("chr", "start", "end", "hap", "strand", "snp", "allele", "qname")
dump <- fread(dump0, skip = 1, header = FALSE, col.names = cols, colClasses = "character")
dump <- dump[hap != "."]
out <- fread(out_w, skip = 1, header = FALSE,
             col.names = c("chr", "pos", "allele", "reads", "Nsum", "N2sum", "Sjd", "mscore", "kappa", "Y_prime"))

# CpG positions covering all fragments
chrom <- unique(dump$chr)
stopifnot(length(chrom) == 1)
rng <- range(as.integer(c(dump$start, dump$end)))
cpg_pos <- sort(as.integer(sapply(strsplit(system2("tabix", c(cpg, sprintf("%s:%d-%d", chrom, rng[1], rng[2])), stdout = TRUE), "\t"), `[`, 2)))

stat_cols <- c("reads", "Nsum", "N2sum", "Sjd", "mscore", "kappa", "Y_prime")
n_bad <- 0L; n_stream <- 0L
for (i in seq_len(nrow(snps))) {
  s <- snps[i]
  key <- sprintf("%s:%d:%s:%s", s$chr, s$pos, s$ref, s$alt)
  lo <- if (w > 0) s$pos - w else -Inf
  hi <- if (w > 0) s$pos + w else Inf
  for (k in 1:2) {
    a <- if (k == 1) s$ref else s$alt
    d <- dump[snp == key & allele == a]
    N <- integer(0); Z <- integer(0)
    for (j in seq_len(nrow(d))) {
      p <- cpg_pos[cpg_pos >= as.integer(d$start[j]) & cpg_pos <= as.integer(d$end[j])]
      ch <- strsplit(d$hap[j], "")[[1]]
      stopifnot(length(p) == length(ch))
      keep <- p >= lo & p <= hi & ch != "."
      N <- c(N, sum(keep)); Z <- c(Z, as.integer(any(ch[keep] == "1")))
    }
    ok <- N > 0
    expect <- mHapDMR:::.mscore_accumulate(N[ok], Z[ok], rep(1L, sum(ok)))
    if (nrow(d) > 0 && !any(grepl(".", d$hap, fixed = TRUE))) {
      raw <- sprintf("%s\t%s\t%s\t%s\t1\t%s", d$chr, d$start, d$end, d$hap, d$strand)
      rs <- if (w > 0) s$pos - w else 1L
      re <- if (w > 0) s$pos + w else .Machine$integer.max
      streamed <- mHapDMR:::.mscore_streaming(raw, cpg_pos, rs, re)
      if (!isTRUE(all.equal(streamed, expect))) { n_bad <- n_bad + 1L; cat("[stream]", key, a, "\n") }
      n_stream <- n_stream + 1L
    }
    row <- out[2L * i - 2L + k]
    stopifnot(row$pos == s$pos, row$allele == a)
    got <- unlist(as.list(row)[stat_cols])
    if (!isTRUE(all.equal(unname(as.numeric(got)), unname(expect), tolerance = 1e-5))) {
      n_bad <- n_bad + 1L
      cat("[mismatch]", key, a, "mhapasm:", got, " R:", expect, "\n")
    }
  }
}
cat(sprintf("window %d: %d SNP alleles checked, %d also via mHapDMR .mscore_streaming, %d mismatches\n",
            w, 2L * nrow(snps), n_stream, n_bad))
quit(status = if (n_bad) 1L else 0L)
