#!/usr/bin/env python3
"""
Per-read emulation of `mhaptools convert` (butyuhao/mHapTools convert.cpp,
whole-file mode) so that its records can be traced back to reads. With
--collapse the output equals that of mhaptools convert (checked on the IGF2,
simulated and BISCUIT test data).

Rules reproduced: flag filter (QCFAIL/UNMAP/DUP/SECONDARY/SUPPLEMENTARY), no
MAPQ filter, Bismark XM filter (X/H/U), direction from flags only (paired reads
must be proper pairs), no CIGAR (read_end = pos + read length, bases indexed
from the start of SEQ), any non-C/T (+) or non-G/A (-) base at a CpG drops the
read, mates merged if their CpG spans overlap (higher base quality wins).

Output (stdout): chr start end hap strand qname   (one line per record)
--collapse writes the mhaptools format instead (chr start end hap count strand).
"""
import argparse
import bisect
import gzip
import subprocess
import sys
from collections import Counter

def load_cpgs(path):
    cpg = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            cpg.setdefault(f[0], []).append(int(f[1]))
    for c in cpg:
        cpg[c].sort()
    return cpg

def read_records(bam, region):
    cmd = ["samtools", "view", bam] + ([region] if region else [])
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        yield line.rstrip("\n").split("\t")
    p.wait()

def haplo(fields, cpgs):
    """(chr, cpg positions, haplotype, quals, direction) or None if the read is dropped."""
    qname, flag, chrom, pos = fields[0], int(fields[1]), fields[2], int(fields[3])
    seq, qual = fields[9], fields[10]
    if flag & (0x200 | 0x4 | 0x400 | 0x100 | 0x800):
        return None
    if seq == "*" or not seq:
        return None
    for t in fields[11:]:
        if t.startswith("XM:Z:") and any(ch in "XHU" for ch in t[5:]):
            return None
    paired, proper = flag & 1, flag & 2
    if paired and proper:
        wc = 0
        if flag & 0x40:
            if flag & 0x10: wc = -1
            elif flag & 0x20: wc = 1
        elif flag & 0x80:
            if flag & 0x10: wc = 1
            elif flag & 0x20: wc = -1
    elif not paired and not proper:
        wc = -1 if flag & 0x10 else 1
    else:
        return None
    rlen = len(seq)
    rstart, rend = pos, pos - 1 + rlen
    cl = cpgs.get(chrom, [])
    lo = bisect.bisect_left(cl, rstart - 500 - 1)
    hi = bisect.bisect_right(cl, rend + 500)
    pos_list, hap, quals = [], [], []
    for c in cl[lo:hi]:
        if wc >= 0:
            if c < rstart: continue
            if c > rend: break
            r = c - rstart
        else:
            if c < rstart - 1: continue
            if c > rend - 1: break
            r = c - rstart + 1
        if r >= rlen or r < 0:
            continue
        b = seq[r]
        if wc >= 0:
            if b == "C": hap.append("1")
            elif b == "T": hap.append("0")
            else: return None
        else:
            if b == "G": hap.append("1")
            elif b == "A": hap.append("0")
            else: return None
        pos_list.append(c)
        quals.append(ord(qual[r]) - 33 if qual != "*" else 255)
    if not pos_list:
        return None
    return chrom, pos_list, hap, quals, wc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--cpg", required=True)
    ap.add_argument("--region")
    ap.add_argument("--collapse", action="store_true")
    a = ap.parse_args()
    cpgs = load_cpgs(a.cpg)
    pending, out = {}, []
    for f in read_records(a.bam, a.region):
        h = haplo(f, cpgs)
        if h is None:
            continue
        q = f[0]
        if q not in pending:
            pending[q] = h
            continue
        F, R = pending.pop(q), h
        checkF = F[1][-1] >= R[1][0]
        checkR = R[1][-1] >= F[1][0]
        if F[0] == R[0] and ((checkF and checkR) or (not checkF and not checkR)):
            calls = {p: (s, ql) for p, s, ql in zip(F[1], F[2], F[3])}
            for p, s, ql in zip(R[1], R[2], R[3]):
                if p not in calls or ql > calls[p][1]:
                    calls[p] = (s, ql)
            ps = sorted(calls)
            out.append((F[0], ps[0], ps[-1], "".join(calls[p][0] for p in ps), F[4], q))
        else:
            for X in (F, R):
                out.append((X[0], X[1][0], X[1][-1], "".join(X[2]), X[4], q))
    for q, X in pending.items():
        out.append((X[0], X[1][0], X[1][-1], "".join(X[2]), X[4], q))
    sym = {1: "+", -1: "-", 0: "*"}
    if a.collapse:
        cnt = Counter((c, s, e, h, sym[w]) for c, s, e, h, w, q in out)
        for (c, s, e, h, w), n in sorted(cnt.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], kv[0][3], kv[0][4])):
            print(c, s, e, h, n, w, sep="\t")
    else:
        for c, s, e, h, w, q in out:
            print(c, s, e, h, sym[w], q, sep="\t")

if __name__ == "__main__":
    main()
