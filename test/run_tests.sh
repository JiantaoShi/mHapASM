#!/bin/bash
# Regression and cross-validation tests for mhapasm and mhapconvert.
#
#   bash test/run_tests.sh                        # golden outputs + thread determinism
#   WGBS_BIN=<dir> bash test/run_tests.sh         # + read-by-read comparison with wgbs_tools
#                                                 #   (dir holding match_maker, snp_patter, patter)
# The IGF2 tests use the wgbs_tools tutorial BAM and the hg19 CpGs around IGF2
# in test/data (override with IGF2_BAM / HG19_CPG).
# The M-score arithmetic is also checked against mHapDMR when Rscript and the
# mHapDMR package are available. Outputs go to test/out/ (overwritten).
set -uo pipefail
cd "$(dirname "$0")/.."
EXE=./mhapasm
OUT=test/out
mkdir -p "$OUT"
fail=0

check() {   # name, command...
    local name=$1; shift
    if "$@" > "$OUT/$name.log" 2>&1; then
        echo "PASS  $name"
    else
        echo "FAIL  $name (see $OUT/$name.log)"
        fail=1
    fi
}

golden() {  # name, expected file, mhapasm args...
    local name=$1 expected=$2; shift 2
    if "$EXE" "$@" -o "$OUT/$name.tsv" > "$OUT/$name.log" 2>&1 &&
       diff "$expected" "$OUT/$name.tsv" >> "$OUT/$name.log" 2>&1; then
        echo "PASS  $name"
    else
        echo "FAIL  $name (see $OUT/$name.log)"
        fail=1
    fi
}

# golden test of mhapconvert: name, expected file, mhapconvert args...
golden_cvt() {
    local name=$1 expected=$2; shift 2
    if ./mhapconvert "$@" -o "$OUT/$name.mhap" > "$OUT/$name.log" 2>&1 &&
       diff "$expected" "$OUT/$name.mhap" >> "$OUT/$name.log" 2>&1; then
        echo "PASS  $name"
    else
        echo "FAIL  $name (see $OUT/$name.log)"
        fail=1
    fi
}

[ -x "$EXE" ] && [ -x ./mhapconvert ] || { echo "build first (make)"; exit 1; }

IGF2_BAM=${IGF2_BAM:-test/data/Left_Ventricle_STL001.IGF2.bam}
HG19_CPG=${HG19_CPG:-test/data/hg19_CpG.IGF2.gz}
IGF2=(-b "$IGF2_BAM" -c "$HG19_CPG" -s test/snp_igf2.txt)
SIM=(-b test/sim/sim.bam -c test/sim/cpg.gz -s test/sim/snps.txt)

if [ -f "$IGF2_BAM" ] && [ -f "$HG19_CPG" ]; then
    golden igf2 test/expected/igf2.tsv "${IGF2[@]}" --qc "$OUT/igf2.qc.tsv"
    check igf2_qc diff test/expected/igf2.qc.tsv "$OUT/igf2.qc.tsv"
else
    echo "SKIP  igf2 (set IGF2_BAM and HG19_CPG)"
fi

golden sim test/expected/sim.tsv "${SIM[@]}" --qc "$OUT/sim.qc.tsv"
check sim_qc diff test/expected/sim.qc.tsv "$OUT/sim.qc.tsv"
golden sim_w100 test/expected/sim.w100.tsv "${SIM[@]}" -w 100
golden sim_threads test/expected/sim.tsv "${SIM[@]}" -@ 2

# mhapconvert
SIMC=(-i test/sim/sim.bam -c test/sim/cpg.gz)
if [ -f "$IGF2_BAM" ] && [ -f "$HG19_CPG" ]; then
    golden_cvt convert_igf2 test/expected/convert_igf2.mhap -i "$IGF2_BAM" -c "$HG19_CPG"
fi
golden_cvt convert_sim test/expected/convert_sim.mhap "${SIMC[@]}"
golden_cvt convert_sim_split test/expected/convert_sim.split.mhap "${SIMC[@]}" --split
golden_cvt convert_sim_nondir test/expected/convert_sim.nondir.mhap "${SIMC[@]}" -n
golden_cvt convert_sim_taps test/expected/convert_sim.taps.mhap "${SIMC[@]}" -m TAPS
golden_cvt convert_sim_bed test/expected/convert_sim.bed.mhap "${SIMC[@]}" -b test/regions.bed
golden_cvt convert_sim_threads test/expected/convert_sim.mhap "${SIMC[@]}" -@ 2
check convert_sim_index bash -c './mhapconvert -i test/sim/sim.bam -c test/sim/cpg.gz -o test/out/convert_sim.mhap.gz &&
    cmp <(tabix test/out/convert_sim.mhap.gz chrS) test/expected/convert_sim.mhap'

if [ -n "${WGBS_BIN:-}" ]; then
    check wgbs_convert_sim python3 test/validate_convert.py --bam test/sim/sim.bam --cpg test/sim/cpg.gz \
        --wgbs-bin "$WGBS_BIN" --region chrS:1-20000
    if [ -f "$IGF2_BAM" ] && [ -f "$HG19_CPG" ]; then
        check wgbs_convert_igf2 python3 test/validate_convert.py --bam "$IGF2_BAM" --cpg "$HG19_CPG" \
            --wgbs-bin "$WGBS_BIN" --region chr11:2000000-2050000
    fi
    check wgbs_sim python3 test/validate_vs_wgbs.py --bam test/sim/sim.bam --cpg test/sim/cpg.gz \
        --wgbs-bin "$WGBS_BIN" --mhapasm "$EXE" --region chrS:1500-18500 --n-pos 10 --snp-file test/sim/snps.txt
    if [ -f "$IGF2_BAM" ] && [ -f "$HG19_CPG" ]; then
        check wgbs_igf2 python3 test/validate_vs_wgbs.py --bam "$IGF2_BAM" --cpg "$HG19_CPG" \
            --wgbs-bin "$WGBS_BIN" --mhapasm "$EXE" --region chr11:2018713-2024739 --n-pos 10 --extra-snp chr11:2019496:C:A
    fi
else
    echo "SKIP  wgbs_tools comparison (set WGBS_BIN)"
fi

if command -v Rscript > /dev/null && Rscript -e 'quit(status = !requireNamespace("mHapDMR", quietly = TRUE))' > /dev/null 2>&1; then
    check mhapdmr_sim_w0 Rscript test/check_stats.R "$EXE" test/sim/sim.bam test/sim/cpg.gz test/sim/snps.txt 0
    check mhapdmr_sim_w100 Rscript test/check_stats.R "$EXE" test/sim/sim.bam test/sim/cpg.gz test/sim/snps.txt 100
else
    echo "SKIP  mHapDMR arithmetic check (needs Rscript + mHapDMR)"
fi

exit $fail
