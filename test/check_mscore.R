#!/usr/bin/env Rscript
# Recompute the per-CpG statistics of mhapmscore in R and compare:
#   - every CpG: the valid records (as many CpGs in [start, end] as haplotype
#     states) are expanded into (read, CpG) pairs and summed in R, which also
#     checks which CpGs are reported;
#   - n_check CpGs (default: all): the records spanning the CpG are passed to
#     mHapDMR:::.mscore_streaming, the code mHapDMR applies to the mHap records
#     of a region, with the region [pos - w, pos + w] (the whole read for w = 0).
#
# Usage: Rscript test/check_mscore.R <mhapmscore> <mhap> <cpg.gz> <window> [n_check]

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: check_mscore.R <mhapmscore> <mhap> <cpg.gz> <window> [n_check]")
exe <- args[1]; mhap <- args[2]; cpg <- args[3]; w <- as.integer(args[4])
n_check <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_
stat_cols <- c("reads", "Nsum", "N2sum", "Sjd", "mscore", "kappa", "Y_prime")

out_fn <- tempfile("mscore", fileext = ".tsv")
if (system2(exe, c("-i", mhap, "-c", cpg, "-w", w, "-o", out_fn), stdout = FALSE, stderr = FALSE) != 0)
  stop("mhapmscore failed")
out <- fread(out_fn, skip = 1, header = FALSE, col.names = c("chr", "pos", stat_cols),
             colClasses = list(character = 1))

rd <- fread(if (grepl("\\.gz$", mhap)) cmd = paste("gzip -dc", shQuote(mhap)) else mhap,
            header = FALSE, sep = "\t", fill = TRUE, colClasses = list(character = c(1, 4)))
setnames(rd, 1:5, c("chr", "start", "end", "hap", "count"))
setorder(rd, chr, start)

cpg_cache <- list()
contig_cpgs <- function(ch) {
  if (is.null(cpg_cache[[ch]])) {
    l <- system2("tabix", c(cpg, ch), stdout = TRUE)
    cpg_cache[[ch]] <<- sort(as.integer(vapply(strsplit(l, "\t", fixed = TRUE), `[`, "", 2)))
  }
  cpg_cache[[ch]]
}

# all CpGs, from (read, CpG) pairs
expect <- rbindlist(lapply(unique(rd$chr), function(ch) {
  r <- rd[chr == ch]
  pos <- contig_cpgs(ch)
  if (!length(pos)) return(NULL)
  lo <- findInterval(r$start - 1L, pos) + 1L
  L <- nchar(r$hap)
  ok <- findInterval(r$end, pos) - lo + 1L == L
  r <- r[ok]; lo <- lo[ok]; L <- L[ok]
  if (!nrow(r)) return(NULL)
  t <- rep(seq_len(nrow(r)), L)
  j <- lo[t] + sequence(L) - 1L
  if (w == 0) {
    N <- L[t]
    Z <- grepl("1", r$hap, fixed = TRUE)[t]
  } else {
    a <- pmax(lo[t], findInterval(pos[j] - w - 1L, pos) + 1L)
    b <- pmin(lo[t] + L[t] - 1L, findInterval(pos[j] + w, pos))
    N <- b - a + 1L
    Z <- grepl("1", substr(r$hap[t], a - lo[t] + 1L, b - lo[t] + 1L), fixed = TRUE)
  }
  d <- data.table(j = j, cnt = as.double(r$count[t]), N = as.double(N), Z = as.double(Z))
  d <- d[, .(reads = sum(cnt), Nsum = sum(cnt * N), N2sum = sum(cnt * N * N), Sjd = sum(cnt * N * Z)), by = j]
  d[, `:=`(chr = ch, pos = pos[j], j = NULL)]
}))
expect[, `:=`(mscore = Sjd / Nsum, kappa = Nsum^2 / N2sum, Y_prime = Sjd * Nsum / N2sum)]

n_bad <- 0L
m <- merge(out, expect, by = c("chr", "pos"), all = TRUE, suffixes = c("", ".R"))
missing <- m[is.na(reads) | is.na(reads.R)]
if (nrow(missing)) {
  n_bad <- n_bad + nrow(missing)
  cat("[set] CpGs reported by only one side:", nrow(missing), "\n"); print(head(missing))
}
m <- m[!is.na(reads) & !is.na(reads.R)]
for (s in stat_cols) {
  bad <- if (s %in% c("reads", "Nsum", "N2sum", "Sjd")) m[[s]] != m[[paste0(s, ".R")]]
         else abs(m[[s]] - m[[paste0(s, ".R")]]) > 1e-5 * pmax(1, abs(m[[paste0(s, ".R")]]))
  if (any(bad)) {
    n_bad <- n_bad + sum(bad)
    cat("[all]", s, "differs at", sum(bad), "CpGs\n"); print(head(m[bad]))
  }
}

# mHapDMR's own code, CpG by CpG
n_stream <- 0L
if (requireNamespace("mHapDMR", quietly = TRUE)) {
  sel <- seq_len(nrow(out))
  if (!is.na(n_check) && n_check < length(sel)) { set.seed(1); sel <- sort(sample(sel, n_check)) }
  span <- max(rd$end - rd$start)
  for (i in sel) {
    o <- out[i]
    r <- rd[chr == o$chr]
    lo <- findInterval(o$pos - span - 1L, r$start) + 1L
    hi <- findInterval(o$pos, r$start)
    r <- r[seq_len(max(0L, hi - lo + 1L)) + lo - 1L][end >= o$pos]
    raw <- sprintf("%s\t%d\t%d\t%s\t%d", r$chr, r$start, r$end, r$hap, r$count)
    rs <- if (w > 0) o$pos - w else 1L
    re <- if (w > 0) o$pos + w else .Machine$integer.max
    got <- unlist(o[, ..stat_cols])
    exp <- mHapDMR:::.mscore_streaming(raw, contig_cpgs(o$chr), rs, re)
    if (!isTRUE(all.equal(unname(as.numeric(got)), unname(exp), tolerance = 1e-5))) {
      n_bad <- n_bad + 1L
      cat("[mHapDMR]", o$chr, o$pos, "mhapmscore:", got, " mHapDMR:", exp, "\n")
    }
    n_stream <- n_stream + 1L
  }
} else cat("mHapDMR not installed: CpG-by-CpG check skipped\n")

cat(sprintf("window %d: %d CpGs reported, %d expected in R, %d also checked with mHapDMR .mscore_streaming, %d mismatches\n",
            w, nrow(out), nrow(expect), n_stream, n_bad))
quit(status = if (n_bad) 1L else 0L)
