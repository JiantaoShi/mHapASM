#!/usr/bin/env python3
"""
Cross-check mhapasm against the wgbs_tools C++ chain, read by read.

  allele assignment : samtools view | match_maker | snp_patter   (run per SNP and allele)
  methylation calls : samtools view | match_maker | patter --long (run once)

SNPs are drawn at random positions inside the covered region of the BAM, with
every ordered pair of bases as REF/ALT, so all bisulfite-specific branches of
the allele-calling rules are exercised. mhapasm is run with --no-mask so that
the CpG calls are directly comparable with patter.

Example:
  python3 test/validate_vs_wgbs.py --bam BAM/Left_Ventricle_STL001.IGF2.bam \
      --cpg CpG/hg19_CpG.gz --wgbs-bin <dir with match_maker, snp_patter, patter> \
      --region chr11:2018713-2024739 --n-pos 40
"""
import argparse
import gzip
import itertools
import os
import random
import subprocess
import sys
import tempfile
from collections import defaultdict


def sh(cmd):
    res = subprocess.run(cmd, shell=True, executable="/bin/bash", capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"command failed ({res.returncode}): {cmd}\n{res.stderr}")
    return res.stdout


def wgbs_allele_reads(a, chrom, pos, let1, let2):
    """QNAMEs that snp_patter assigns to let1 (wgbs_tools split_by_allele)."""
    lo, hi = max(1, pos - 1000), pos + 1001
    cmd = (f"samtools view {a.bam} {chrom}:{lo}-{hi} -q {a.mapq} -F {a.flags}"
           f" | {a.wgbs_bin}/match_maker - 2>/dev/null"
           f" | {a.wgbs_bin}/snp_patter --snp_pos {pos} --snp_let1 {let1} --snp_let2 {let2}"
           f" --qual_filter {a.snp_bq} 2>/dev/null | cut -f1")
    return set(sh(cmd).split())


def make_dict(a, chrom, tmp):
    """CpG dictionary for patter: chr, 1-based C position, running CpG index."""
    path = os.path.join(tmp, "cpg.dict.gz")
    sh(f"tabix {a.cpg} {chrom} | awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,NR}}' | bgzip > {path}"
       f" && tabix -f -s1 -b2 -e2 {path}")
    pos2idx = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, p, i = line.split()
            pos2idx[int(p)] = int(i)
    return path, pos2idx


def patter_patterns(a, chrom, beg, end, dict_path):
    """QNAME -> (first CpG index, pattern) from patter --long."""
    lo, hi = max(1, beg - 2000), end + 2000
    cmd = (f"samtools view {a.bam} {chrom}:{lo}-{hi} -q {a.mapq} -F {a.flags}"
           f" | {a.wgbs_bin}/match_maker - 2>/dev/null"
           f" | {a.wgbs_bin}/patter {dict_path} {chrom}:{lo}-{hi} --long 2>/dev/null")
    out = {}
    for line in sh(cmd).splitlines():
        c, idx, pat, qname = line.split("\t")[:4]
        out[qname] = (int(idx), pat)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--cpg", required=True)
    ap.add_argument("--wgbs-bin", required=True, help="directory with match_maker, snp_patter, patter")
    ap.add_argument("--mhapasm", default="./mhapasm")
    ap.add_argument("--region", required=True, help="chr:beg-end to draw SNP positions from")
    ap.add_argument("--n-pos", type=int, default=40)
    ap.add_argument("--extra-snp", action="append", default=[], help="chr:pos:ref:alt, always tested")
    ap.add_argument("--snp-file", help="'chr pos ref alt' file, always tested")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mapq", type=int, default=10)
    ap.add_argument("--flags", type=int, default=1796)
    ap.add_argument("--snp-bq", type=int, default=0)
    a = ap.parse_args()

    chrom, rng = a.region.split(":")
    beg, end = map(int, rng.replace(",", "").split("-"))
    random.seed(a.seed)
    positions = sorted(random.sample(range(beg, end + 1), a.n_pos))
    snps = [(chrom, p, r, t) for p in positions for r, t in itertools.permutations("ACGT", 2)]
    for s in a.extra_snp:
        c, p, r, t = s.split(":")
        snps.append((c, int(p), r, t))
    if a.snp_file:
        with open(a.snp_file) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    c, p, r, t = line.split()[:4]
                    snps.append((c, int(p), r, t))

    with tempfile.TemporaryDirectory() as tmp:
        snp_fn = os.path.join(tmp, "snps.txt")
        with open(snp_fn, "w") as fh:
            for s in snps:
                fh.write("%s\t%d\t%s\t%s\n" % s)
        dump_fn = os.path.join(tmp, "dump.tsv")
        sh(f"{a.mhapasm} -b {a.bam} -c {a.cpg} -s {snp_fn} -q {a.mapq} -F {a.flags} -Q {a.snp_bq}"
           f" --no-mask --dump {dump_fn} -o {os.path.join(tmp, 'out.tsv')}")

        mine = defaultdict(set)        # (snp key, allele) -> QNAMEs
        mine_pat = {}                  # QNAME -> (first CpG position, haplotype)
        with open(dump_fn) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c, s, e, hap, strand, key, allele, qname = line.rstrip("\n").split("\t")
                mine[(key, allele)].add(qname)
                if hap != ".":
                    mine_pat[qname] = (int(s), hap)

        # allele assignment, SNP by SNP
        n_bad = n_reads = 0
        for c, p, r, t in snps:
            key = f"{c}:{p}:{r}:{t}"
            for let1, let2 in ((r, t), (t, r)):
                ref_set = wgbs_allele_reads(a, c, p, let1, let2)
                got = mine.get((key, let1), set())
                n_reads += len(ref_set)
                if got != ref_set:
                    n_bad += 1
                    print(f"[allele] {key} {let1}: wgbs-only {sorted(ref_set - got)[:5]} mhapasm-only {sorted(got - ref_set)[:5]}")
        print(f"allele assignment: {len(snps)} SNPs x 2 alleles, {n_reads} wgbs read assignments, {n_bad} mismatching sets")

        # methylation patterns of every fragment that mhapasm assigned anywhere
        dict_path, pos2idx = make_dict(a, chrom, tmp)
        ref_pat = patter_patterns(a, chrom, beg, end, dict_path)
        n_cmp = n_diff = 0
        for qname, (start, hap) in mine_pat.items():
            n_cmp += 1
            conv = (pos2idx[start], hap.replace("1", "C").replace("0", "T"))
            if ref_pat.get(qname) != conv:
                n_diff += 1
                print(f"[meth] {qname}: patter {ref_pat.get(qname)} mhapasm {conv}")
        print(f"methylation calls: {n_cmp} fragments compared with patter, {n_diff} differ")
        sys.exit(1 if (n_bad or n_diff) else 0)


if __name__ == "__main__":
    main()
