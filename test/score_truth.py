#!/usr/bin/env python3
"""
Score per-read mHap records against the simulated methylation truth.

records: mhapconvert --qname output (gz or plain), or chr start end hap strand
         qname (mhaptools_emulate.py); truth: simulate_pe.py --meth-truth (run
         with --out-dir equal to --truth-dir). Each call is correct,
wrong (state), at a CpG that the haplotype does not have (x), or outside the
fragment. Fragments are grouped by whether a read has an indel / soft clip and
whether they come from the complementary strands (CTOT/CTOB).
"""
import argparse
import bisect
import gzip
import subprocess
from collections import defaultdict, Counter

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", required=True)
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    d = a.truth_dir
    cpgs = [int(l.split("\t")[1]) for l in gzip.open(f"{d}/cpg.gz", "rt")]
    truth = {}
    for line in open(f"{d}/meth_truth.tsv"):
        q, s = line.rstrip("\n").split("\t")
        truth[q] = dict((int(p), v) for p, v in (x.split(":") for x in s.split(",") if x))
    strand = {l.split("\t")[0]: l.split("\t")[2] for l in open(f"{d}/truth.tsv") if not l.startswith("qname")}
    # fragment categories from the BAM
    cat = {}
    out = subprocess.run(["samtools", "view", f"{d}/sim.bam"], capture_output=True, text=True).stdout
    indel, comp = set(), set()
    for line in out.splitlines():
        f = line.split("\t")
        q, flag, cig = f[0], int(f[1]), f[5]
        if any(op in cig for op in "IDS"):
            indel.add(q)
        if not (flag & 1) or (flag & 0x40):          # read 1 or single-end read
            rev = bool(flag & 16)
            if (strand[q] == "+") == rev:
                comp.add(q)
    for q in strand:
        cat[q] = "complementary strand" if q in comp else "indel/soft clip" if q in indel else "plain"
    stats = defaultdict(Counter)
    frags = defaultdict(set)
    fh = gzip.open(a.records, "rt") if a.records.endswith(".gz") else open(a.records)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        c, s, e, hap, q = f[0], f[1], f[2], f[3], f[-1]
        lo = bisect.bisect_left(cpgs, int(s))
        pos = cpgs[lo:lo + len(hap)]
        assert pos[-1] == int(e), line
        k = cat[q]
        frags[k].add(q)
        for p, h in zip(pos, hap):
            t = truth[q].get(p)
            r = "outside" if t is None else "noCpG" if t == "x" else "ok" if t == h else "wrong"
            stats[k][r] += 1
            stats["all"][r] += 1
        frags["all"].add(q)
    n_cat = Counter(cat.values())
    n_cat["all"] = len(cat)
    for k in ["all", "plain", "indel/soft clip", "complementary strand"]:
        s = stats[k]
        tot = sum(s.values())
        if not n_cat[k]:
            continue
        print(f"{a.label:22s} {k:21s} fragments {len(frags[k]):6d}/{n_cat[k]:<6d} calls {tot:7d}  "
              f"accuracy {100 * s['ok'] / tot if tot else 0:6.2f}%  wrong {s['wrong']:5d}  "
              f"noCpG {s['noCpG']:4d}  outside {s['outside']:4d}")

if __name__ == "__main__":
    main()
