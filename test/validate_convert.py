#!/usr/bin/env python3
"""
Cross-check mhapconvert against wgbs_tools, read by read.

wgbs_tools (samtools view | match_maker | patter --long) gives one CpG pattern
per fragment, with '.' for CpGs without a call. mHap cannot hold missing
values, so with --split mhapconvert must write exactly the '.'-free runs of
that pattern, one record per run. Without --split (the default) a fragment is
either written exactly like that or, if it had a gap, dropped; the script
checks both.

Example:
  python3 test/validate_convert.py --bam test/sim/sim.bam --cpg test/sim/cpg.gz \
      --wgbs-bin <dir with match_maker and patter> --region chrS:1-20000
"""
import argparse
import gzip
import os
import subprocess
import sys
import tempfile
from collections import defaultdict


def sh(cmd):
    res = subprocess.run(cmd, shell=True, executable="/bin/bash", capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"command failed ({res.returncode}): {cmd}\n{res.stderr}")
    return res.stdout


def records(path):
    """QNAME -> sorted list of (start, hap) from mhapconvert --qname output."""
    out = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            chrom, start, end, hap, count, strand, qname = line.rstrip("\n").split("\t")
            out[qname].append((int(start), hap))
    return {q: sorted(v) for q, v in out.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--cpg", required=True)
    ap.add_argument("--wgbs-bin", required=True)
    ap.add_argument("--mhapconvert", default="./mhapconvert")
    ap.add_argument("--region", required=True, help="chr:beg-end (one contig)")
    ap.add_argument("--mapq", type=int, default=10)
    a = ap.parse_args()
    chrom = a.region.split(":")[0]

    with tempfile.TemporaryDirectory() as tmp:
        # CpG dictionary for patter: chr, position, CpG index
        dict_path = os.path.join(tmp, "cpg.dict.gz")
        sh(f"tabix {a.cpg} {chrom} | awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,NR}}' | bgzip > {dict_path}"
           f" && tabix -f -s1 -b2 -e2 {dict_path}")
        idx2pos = {}
        with gzip.open(dict_path, "rt") as fh:
            for line in fh:
                c, p, i = line.split()
                idx2pos[int(i)] = int(p)

        # wgbs_tools patterns
        cmd = (f"samtools view {a.bam} {a.region} -q {a.mapq} -F 1796"
               f" | {a.wgbs_bin}/match_maker - 2>/dev/null"
               f" | {a.wgbs_bin}/patter {dict_path} {a.region} --long 2>/dev/null")
        expect = {}
        for line in sh(cmd).splitlines():
            c, idx, pat, qname = line.split("\t")[:4]
            runs, i0 = [], int(idx)
            for k, ch in enumerate(pat + "."):
                if ch == ".":
                    continue
                if k == 0 or pat[k - 1] == ".":
                    runs.append([idx2pos[i0 + k], ""])
                runs[-1][1] += "1" if ch == "C" else "0"
            expect[qname] = sorted((s, h) for s, h in runs)

        # wgbs_tools has no conversion filter: compare without it
        common = f"{a.mhapconvert} -i {a.bam} -c {a.cpg} -r {a.region} -q {a.mapq} -F 1796 --max-unconv -1 --qname"
        split_fn, drop_fn = os.path.join(tmp, "split.mhap.gz"), os.path.join(tmp, "drop.mhap.gz")
        sh(f"{common} --split -o {split_fn} 2>/dev/null")
        sh(f"{common} -o {drop_fn} 2>/dev/null")
        got_split, got_drop = records(split_fn), records(drop_fn)

        # Reads that the region query cuts off from their mate are excluded
        # from neither tool: both see the same reads, so compare everything.
        bad = 0
        for q in set(expect) | set(got_split):
            if expect.get(q) != got_split.get(q):
                bad += 1
                if bad <= 10:
                    print(f"[split] {q}: wgbs {expect.get(q)} mhapconvert {got_split.get(q)}")
        print(f"--split: {len(expect)} fragments with calls in wgbs_tools, {bad} differ")

        bad_kept = bad_drop = 0
        for q, recs in got_drop.items():
            if recs != got_split.get(q):
                bad_kept += 1
        for q, recs in got_split.items():
            if q not in got_drop and len(recs) < 2:
                bad_drop += 1        # dropped although it had no gap
        dropped = len(set(got_split) - set(got_drop))
        print(f"default: {len(got_drop)} fragments written, identical to --split: {len(got_drop) - bad_kept}; "
              f"{dropped} dropped, all with a gap: {bad_drop == 0}")
        sys.exit(1 if (bad or bad_kept or bad_drop) else 0)


if __name__ == "__main__":
    main()
