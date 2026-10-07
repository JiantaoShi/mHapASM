#!/usr/bin/env python3
"""
Simulate a small bisulfite sequencing BAM for testing mhapasm and mhapconvert.

A random reference contig carries phased heterozygous SNPs (REF alleles on
haplotype 0, ALT alleles on haplotype 1, some of them CpG-destroying C>T / G>A).
Fragments are drawn from either haplotype, methylated with allele-specific
probabilities, bisulfite converted on the original top (OT) or bottom (OB)
strand, and sequenced as read pairs with the directional flags 99/147 (OT) and
83/163 (OB). Sequencing errors, soft clips, insertions, deletions, low MAPQ,
duplicates and missing mates are mixed in. Options add single-end reads,
the strand tags of an aligner, and non-directional (CTOT/CTOB) fragments,
whose flags look like OB/OT although they carry the OT/OB conversion. Without
these options the output is the same as before they existed.

Writes into --out-dir:
  sim.bam(.bai)       coordinate-sorted, indexed
  cpg.gz(.tbi)        chr, 1-based C position, 1-based G position (tabix -s1 -b2 -e3)
  snps.txt            chr pos ref alt
  truth.tsv           qname, haplotype, strand, fragment start/end (1-based)
  --meth-truth FILE   qname, then pos:state for every reference CpG of the
                      fragment (1 / 0, x if the haplotype has no CpG there)
"""
import argparse
import bisect
import os
import random
import subprocess

COMP = str.maketrans("ACGT", "TGCA")


def rand_base(rnd, weights):
    return rnd.choices("ACGT", weights=weights)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="test/sim")
    ap.add_argument("--seed", type=int, default=42, help="seed of the reference, SNPs and CpGs")
    ap.add_argument("--frag-seed", type=int, help="seed of the fragments [--seed]; vary it to get replicate samples")
    ap.add_argument("--length", type=int, default=20000)
    ap.add_argument("--n-frag", type=int, default=8000)
    ap.add_argument("--n-snp", type=int, default=40)
    ap.add_argument("--se", action="store_true", help="single-end reads")
    ap.add_argument("--aligner", choices=["none", "bismark", "bwameth", "biscuit", "bsmap"], default="none",
                    help="strand tags to write: XG (Bismark), YD:Z (bwa-meth), YD:A (BISCUIT), ZS (BSMAP)")
    ap.add_argument("--non-directional", type=float, default=0.0,
                    help="fraction of fragments read from the complementary strands (CTOT/CTOB)")
    ap.add_argument("--unconverted", type=float, default=0.0,
                    help="fraction of fragments left unconverted by bisulfite (they look fully methylated)")
    ap.add_argument("--meth-truth", help="write the true CpG states of every fragment")
    a = ap.parse_args()
    rnd = random.Random(a.seed)
    os.makedirs(a.out_dir, exist_ok=True)
    chrom, L = "chrS", a.length

    # reference: AT-rich background with three CpG-rich islands
    ref = [rand_base(rnd, (3, 2, 2, 3)) for _ in range(L)]
    for isl in (3000, 9000, 15000):
        for i in range(isl, isl + 1200):
            ref[i] = rand_base(rnd, (2, 3, 3, 2))
    ref = "".join(ref)
    cpgs = [i for i in range(L - 1) if ref[i] == "C" and ref[i + 1] == "G"]

    # SNPs: random substitutions plus CpG-destroying ones at C (C>T) and G (G>A)
    snps = {}
    cand = list(range(1500, L - 1500))
    for p in rnd.sample(cand, a.n_snp):
        snps[p] = rnd.choice([b for b in "ACGT" if b != ref[p]])
    for c in rnd.sample([c for c in cpgs if 1500 < c < L - 1500], 6):
        if rnd.random() < 0.5:
            snps[c] = "T"
        else:
            snps[c + 1] = "A"
    hap = [list(ref), list(ref)]
    for p, alt in snps.items():
        hap[1][p] = alt
    hap = ["".join(h) for h in hap]
    if a.frag_seed is not None:
        rnd = random.Random(a.frag_seed)

    # allele-specific methylation: first half strongly imprinted, second half
    # with all-or-none fragments at equal rates on both haplotypes
    def meth_prob(h, c, frag_state):
        if c < L // 2:
            return 0.92 if h == 1 else 0.06
        return 0.95 if frag_state else 0.05

    with open(os.path.join(a.out_dir, "snps.txt"), "w") as fh:
        for p in sorted(snps):
            fh.write(f"{chrom}\t{p + 1}\t{ref[p]}\t{snps[p]}\n")
    cpg_txt = os.path.join(a.out_dir, "cpg")
    with open(cpg_txt, "w") as fh:
        for c in cpgs:
            fh.write(f"{chrom}\t{c + 1}\t{c + 2}\n")
    subprocess.run(f"bgzip -f {cpg_txt} && tabix -f -s1 -b2 -e3 {cpg_txt}.gz", shell=True, check=True)

    sam_path = os.path.join(a.out_dir, "sim.sam")
    truth = open(os.path.join(a.out_dir, "truth.tsv"), "w")
    truth.write("qname\thap\tstrand\tstart\tend\n")
    sam = open(sam_path, "w")
    mtruth = open(a.meth_truth, "w") if a.meth_truth else None
    sam.write(f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{chrom}\tLN:{L}\n")

    def qual_str(n):
        return "".join(chr(33 + (rnd.randint(5, 15) if rnd.random() < 0.03 else rnd.randint(25, 40))) for _ in range(n))

    def layout(rl, kind):
        """CIGAR operations of one read as (op, length) pairs."""
        if kind == "del":
            d = rnd.randint(20, rl - 20)
            return [("M", d), ("D", rnd.randint(1, 3)), ("M", rl - d)]
        if kind == "ins":
            d, il = rnd.randint(20, rl - 20), rnd.randint(1, 2)
            return [("M", d), ("I", il), ("M", rl - il - d)]
        if kind == "sclip5":
            k = rnd.randint(5, 15)
            return [("S", k), ("M", rl - k)]
        if kind == "sclip3":
            k = rnd.randint(5, 15)
            return [("M", rl - k), ("S", k)]
        return [("M", rl)]

    def render(conv, pos, ops):
        """Read bases in reference orientation; I/S bases are random."""
        out, r = [], pos
        for op, ln in ops:
            if op == "M":
                out.append(conv[r:r + ln])
                r += ln
            elif op == "D":
                r += ln
            else:
                out.append("".join(rnd.choice("ACGT") for _ in range(ln)))
        return "".join(out)

    n = 0
    for i in range(a.n_frag):
        h = rnd.randint(0, 1)
        flen = max(110, min(700, int(rnd.gauss(320, 90))))
        fs = rnd.randint(0, L - flen - 5)
        fe = fs + flen
        bottom = rnd.random() < 0.5
        frag_state = rnd.random() < 0.3
        hs = hap[h]
        # bisulfite-converted fragment in reference orientation, offset fs
        # (a few extra bases so that a deletion near the end stays in range)
        seg_end = min(L, fe + 5)
        unconv = a.unconverted > 0 and rnd.random() < a.unconverted
        conv = list(hs[fs:seg_end])
        state = {}                          # CpG C position -> drawn methylation
        for j in range(fs, seg_end):
            if not bottom and hs[j] == "C":
                meth = j + 1 < L and hs[j + 1] == "G" and rnd.random() < meth_prob(h, j, frag_state)
                conv[j - fs] = "C" if meth or unconv else "T"
                if j + 1 < L and hs[j + 1] == "G":
                    state[j] = meth
            elif bottom and hs[j] == "G":
                meth = j > 0 and hs[j - 1] == "C" and rnd.random() < meth_prob(h, j - 1, frag_state)
                conv[j - fs] = "G" if meth or unconv else "A"
                if j > 0 and hs[j - 1] == "C":
                    state[j - 1] = meth
        conv = "".join(conv)
        rl = 150 if rnd.random() < 0.3 else 100
        rl = min(rl, flen)
        qname = f"frag{i}"
        dup = rnd.random() < 0.02
        # complementary strands (non-directional libraries): same conversion,
        # read 1 from the other end of the fragment
        comp = a.non_directional > 0 and rnd.random() < a.non_directional
        mates = []
        sides = ("L", "R") if not a.se else (("R",) if bottom != comp else ("L",))
        for side in sides:
            kind = rnd.choices(["plain", "del", "ins", "sclip5", "sclip3"], weights=(85, 4, 4, 3.5, 3.5))[0]
            if kind in ("del", "ins") and (rl < 60 or rl + 3 > flen):
                kind = "plain"
            ops = layout(rl, kind)
            span = sum(ln for op, ln in ops if op in "MD")
            # the left mate starts at the fragment start, the right mate ends at the fragment end
            pos = fs if side == "L" else fe - span
            seq = render(conv, pos - fs, ops)
            seq = "".join(rnd.choice([b for b in "ACGT" if b != x]) if rnd.random() < 0.005 else x for x in seq)
            cig = "".join(f"{ln}{op}" for op, ln in ops)
            mates.append([pos, cig, seq, qual_str(len(seq))])
        # directional library: OT = read1 forward on the left; OB = read1 reverse
        # on the right. CTOT/CTOB swap the ends, so their flags look like OB/OT.
        if a.se:
            (p1, c1, s1, q1), = mates
            recs = [(16 if bottom != comp else 0, p1, c1, s1, q1, -1, 0)]
        else:
            (lp, lc, ls, lq), (rp, rc, rs, rq) = mates
            if bottom == comp:
                recs = [(99, lp, lc, ls, lq, rp, 1), (147, rp, rc, rs, rq, lp, -1)]
            else:
                recs = [(163, lp, lc, ls, lq, rp, 1), (83, rp, rc, rs, rq, lp, -1)]
            drop = rnd.random() < 0.03      # one mate missing from the BAM
            if drop:
                recs = [recs[rnd.randint(0, 1)]]
        for flag, pos, cig, seq, qual, mpos, sign in recs:
            mapq = 3 if rnd.random() < 0.03 else 60
            if dup:
                flag |= 0x400
            tags = ""
            if a.aligner == "bismark":
                tags = "\tXG:Z:" + ("GA" if bottom else "CT")
            elif a.aligner == "bwameth":
                tags = "\tYD:Z:" + ("r" if bottom else "f")
            elif a.aligner == "biscuit":
                tags = "\tYD:A:" + ("r" if bottom else "f")
            elif a.aligner == "bsmap":
                tags = "\tZS:Z:" + ("-" if bottom else "+") + ("-" if flag & 16 else "+")
            mate = f"=\t{mpos + 1}\t{sign * flen}" if not a.se else "*\t0\t0"
            sam.write(f"{qname}\t{flag}\t{chrom}\t{pos + 1}\t{mapq}\t{cig}\t{mate}\t{seq}\t{qual}{tags}\n")
            n += 1
        truth.write(f"{qname}\t{h}\t{'-' if bottom else '+'}\t{fs + 1}\t{fe}\n")
        if mtruth:
            cs = cpgs[bisect.bisect_left(cpgs, fs - 1):bisect.bisect_left(cpgs, fe)]
            mtruth.write(qname + "\t" + ",".join(
                f"{c + 1}:{'x' if c not in state else int(state[c])}" for c in cs) + "\n")
    sam.close()
    if mtruth:
        mtruth.close()
    truth.close()
    bam = os.path.join(a.out_dir, "sim.bam")
    subprocess.run(f"samtools sort -o {bam} {sam_path} 2>/dev/null && samtools index {bam}", shell=True, check=True)
    os.remove(sam_path)
    print(f"{a.n_frag} fragments, {n} reads, {len(snps)} SNPs, {len(cpgs)} CpGs -> {a.out_dir}")


if __name__ == "__main__":
    main()
