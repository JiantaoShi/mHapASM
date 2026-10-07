#!/usr/bin/env python3
"""
Fragment-by-fragment comparison of mhaptools (per-read emulation) and
mhapconvert (--qname) records, with the reason for every difference taken
from the reads: MAPQ, conversion filter, proper pair, aligner strand tag vs
flags, CIGAR.

usage: compare_mhaptools.py BAM MHT_RECORDS CVT_RECORDS [CVT_NO_FILTER]
  MHT_RECORDS:   chr start end hap strand qname   (mhaptools_emulate.py)
  CVT_RECORDS:   mhapconvert --qname output (gz or plain)
  CVT_NO_FILTER: the same with --max-unconv -1; fragments whose records it
                 changes are attributed to the conversion filter
"""
import gzip
import subprocess
import sys
from collections import Counter, defaultdict

bam, mht_fn, cvt_fn = sys.argv[1:4]
nofilter_fn = sys.argv[4] if len(sys.argv) > 4 else None

def opener(fn):
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)

def load_cvt(fn):
    out = defaultdict(set)
    for line in opener(fn):
        c, s, e, h, n, st, q = line.rstrip("\n").split("\t")
        out[q].add((c, s, e, h, st))
    return out

A = defaultdict(set)
for line in opener(mht_fn):
    c, s, e, h, st, q = line.rstrip("\n").split("\t")
    A[q].add((c, s, e, h, st))
B = load_cvt(cvt_fn)
unconv = set()
if nofilter_fn:
    C = load_cvt(nofilter_fn)
    unconv = {q for q in set(B) | set(C) if B.get(q) != C.get(q)}
    del C

attr = defaultdict(set)
p = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE, text=True)
for line in p.stdout:
    f = line.split("\t", 12)
    q, flag, mapq, cig = f[0], int(f[1]), int(f[4]), f[5]
    if flag & 0xF04:
        continue
    if mapq < 10:
        attr[q].add("mapq")
    if (flag & 1) and not (flag & 2):
        attr[q].add("improper")
    if any(op in cig for op in "IDS"):
        attr[q].add("cigar")
    yd = next((t[5] for t in f[11:] if t.startswith("YD:")), None) if len(f) > 11 else None
    if yd in ("f", "r"):
        rev = bool(flag & 16)
        bottom_flag = (rev if flag & 0x40 else not rev) if flag & 1 else rev
        if (yd == "r") != bottom_flag:
            attr[q].add("tag")
p.wait()

def subclass(a, b):
    a = sorted(a, key=lambda r: int(r[1]))
    b = sorted(b, key=lambda r: int(r[1]))
    if len(a) == 2 and len(b) == 1 and a[0][1] == b[0][1] and a[1][2] == b[0][2] and a[0][3] + a[1][3] == b[0][3]:
        return "mates joined into one record (no CpG between them)"
    if len(a) == len(b) and all(len(x[3]) > len(y[3]) and y[3] in x[3] and int(x[1]) <= int(y[1])
                                and int(x[2]) >= int(y[2]) for x, y in zip(a, b)):
        return "CpG at a read end called by mhaptools only"
    if len(a) == len(b) and all(x[1:3] == y[1:3] for x, y in zip(a, b)):
        return "same CpGs, different calls (mate conflict)"
    return "other"

examples = []
cat = Counter()
for q in set(A) | set(B):
    a, b = A.get(q, set()), B.get(q, set())
    if a == b:
        cat["identical"] += 1
        continue
    t = attr.get(q, set())
    if "mapq" in t:
        r = "MAPQ < 10 (filtered by mhapconvert)"
    elif q in unconv:
        r = "incompletely converted read (dropped by mhapconvert)"
    elif "improper" in t:
        r = "improper pair (dropped by mhaptools)"
    elif "tag" in t:
        r = "strand: aligner tag disagrees with flags"
    elif "cigar" in t:
        r = "CIGAR: indel or soft clip"
    elif a and not b:
        r = "only mhaptools: mhapconvert dropped (CpG without call / context)"
    elif b and not a:
        r = "only mhapconvert: mhaptools dropped (non-C/T base at a CpG)"
    else:
        r = "both, different: " + subclass(a, b)
        if r.endswith("other") and len(examples) < 8:
            examples.append((q, sorted(a), sorted(b)))
    cat[r] += 1
tot = sum(cat.values())
print(f"{'fragments':72s} {tot:9d}")
for k, v in cat.most_common():
    print(f"{k:72s} {v:9d}  {100 * v / tot:6.2f}%")
for e in examples:
    print("example:", e)
