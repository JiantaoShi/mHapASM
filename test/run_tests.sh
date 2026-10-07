#!/bin/bash
# Regression and cross-validation tests for mhapasm, mhapconvert and mhapmscore.
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

# golden test of mhapmscore: name, expected file, mhapmscore args...
golden_ms() {
    local name=$1 expected=$2; shift 2
    if ./mhapmscore "$@" -o "$OUT/$name.tsv" > "$OUT/$name.log" 2>&1 &&
       diff "$expected" "$OUT/$name.tsv" >> "$OUT/$name.log" 2>&1; then
        echo "PASS  $name"
    else
        echo "FAIL  $name (see $OUT/$name.log)"
        fail=1
    fi
}

# mhapmscore with -r/-b must give the whole-file output restricted to the
# CpGs in the regions: awk condition on pos ($2), mhapmscore args...
ms_regions() {
    local cond=$1; shift
    cmp <(./mhapmscore -i "$OUT/convert_sim.mhap.gz" -c test/sim/cpg.gz "$@") \
        <(awk "NR == 1 || ($cond)" test/expected/mscore_sim.tsv)
}

# the records of a contig may come in any order
ms_shuffled() {
    awk 'BEGIN { srand(7) } { print rand() "\t" $0 }' test/expected/convert_sim.mhap |
        sort -k1,1g | cut -f2- > "$OUT/sim.shuffled.mhap" &&
    ./mhapmscore -i "$OUT/sim.shuffled.mhap" -c test/sim/cpg.gz | cmp - test/expected/mscore_sim.tsv
}

# ... but must be together: a contig seen again later is an error
ms_split_contig() {
    { cat test/expected/convert_sim.mhap; sed 's/^chrS/chrT/' test/expected/convert_sim.mhap | head -1
      head -1 test/expected/convert_sim.mhap; } > "$OUT/sim.split_contig.mhap"
    ! ./mhapmscore -i "$OUT/sim.split_contig.mhap" -c test/sim/cpg.gz > /dev/null
}

# -R: indexed output for a sorted BED file (0-based starts), none for an
# unsorted one; -w is refused
ms_regions_index() {
    rm -f "$OUT"/mscore_regions*.tsv.gz.tbi
    ./mhapmscore -i "$OUT/convert_sim.mhap.gz" -c test/sim/cpg.gz -R test/regions.bed -o "$OUT/mscore_regions.tsv.gz" &&
    [ "$(tabix "$OUT/mscore_regions.tsv.gz" chrS:4500-4600 | cut -f2,3 | tr '\n' ' ')" = "3000	5000 4000	7000 " ] &&
    ./mhapmscore -i "$OUT/convert_sim.mhap.gz" -c test/sim/cpg.gz -R test/regions_mscore.bed -o "$OUT/mscore_regions2.tsv.gz" &&
    [ ! -e "$OUT/mscore_regions2.tsv.gz.tbi" ] &&
    ! ./mhapmscore -i "$OUT/convert_sim.mhap.gz" -c test/sim/cpg.gz -R test/regions.bed -w 100 > /dev/null
}

# the region modes load CpGs in windows and widen them for reads that stick
# out; a build with 64-bp windows (make test) must give the same output
ms_tinywin() {
    local t=test/out/mhapmscore_tinywin m=(-i "$OUT/convert_sim.mhap.gz" -c test/sim/cpg.gz)
    cmp <("$t" "${m[@]}" -R test/regions_mscore.bed) test/expected/mscore_sim.regions.tsv &&
    cmp <("$t" "${m[@]}" -b test/regions.bed) <(./mhapmscore "${m[@]}" -b test/regions.bed) &&
    cmp <("$t" "${m[@]}" -r chrS:3001-8000 -w 150) <(./mhapmscore "${m[@]}" -r chrS:3001-8000 -w 150) &&
    cmp <("$t" "${m[@]}" -r chrS) test/expected/mscore_sim.tsv
}

ms_index() {
    ./mhapmscore -i "$OUT/convert_sim.mhap.gz" -c test/sim/cpg.gz -@ 2 -o "$OUT/mscore_sim.tsv.gz" &&
    cmp <(gzip -dc "$OUT/mscore_sim.tsv.gz") test/expected/mscore_sim.tsv &&
    cmp <(tabix "$OUT/mscore_sim.tsv.gz" chrS) <(tail -n +2 test/expected/mscore_sim.tsv)
}

# records that do not match the CpG file (haplotype too long, start off a
# CpG) must be skipped as mHapDMR skips them
ms_mismatch() {
    awk 'BEGIN { OFS = "\t" } NR % 50 == 0 { $4 = $4 "0" } NR % 77 == 0 { $2 = $2 + 1 } { print }' \
        test/expected/convert_sim.mhap > "$OUT/sim.mismatch.mhap" &&
    Rscript test/check_mscore.R ./mhapmscore "$OUT/sim.mismatch.mhap" test/sim/cpg.gz 100
}

[ -x "$EXE" ] && [ -x ./mhapconvert ] && [ -x ./mhapmscore ] || { echo "build first (make)"; exit 1; }

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

# mhapmscore, on the mHap files of mhapconvert
if [ -f "$HG19_CPG" ]; then
    golden_ms mscore_igf2 test/expected/mscore_igf2.tsv -i test/expected/convert_igf2.mhap -c "$HG19_CPG"
fi
golden_ms mscore_sim test/expected/mscore_sim.tsv -i test/expected/convert_sim.mhap -c test/sim/cpg.gz
golden_ms mscore_sim_w150 test/expected/mscore_sim.w150.tsv -i test/expected/convert_sim.mhap -c test/sim/cpg.gz -w 150
check mscore_sim_region ms_regions '$2 >= 3001 && $2 <= 8000' -r chrS:3001-8000
check mscore_sim_bed ms_regions '($2 > 3000 && $2 <= 7000) || ($2 > 15000 && $2 <= 16000)' -b test/regions.bed
check mscore_sim_shuffled ms_shuffled
check mscore_sim_split_contig ms_split_contig
check mscore_sim_index ms_index
golden_ms mscore_sim_R test/expected/mscore_sim.regions.tsv -i "$OUT/convert_sim.mhap.gz" -c test/sim/cpg.gz -R test/regions_mscore.bed
check mscore_sim_R_index ms_regions_index
if [ -x test/out/mhapmscore_tinywin ]; then
    check mscore_sim_tiny_windows ms_tinywin
else
    echo "SKIP  mscore_sim_tiny_windows (run make test)"
fi

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
    check mhapdmr_mscore_sim_w0 Rscript test/check_mscore.R ./mhapmscore test/expected/convert_sim.mhap test/sim/cpg.gz 0
    check mhapdmr_mscore_sim_w150 Rscript test/check_mscore.R ./mhapmscore test/expected/convert_sim.mhap test/sim/cpg.gz 150
    check mhapdmr_mscore_mismatch ms_mismatch
    check mhapdmr_mscore_regions Rscript test/check_mscore_regions.R ./mhapmscore "$OUT/convert_sim.mhap.gz" \
        test/sim/cpg.gz test/regions_mscore.bed
else
    echo "SKIP  mHapDMR arithmetic check (needs Rscript + mHapDMR)"
fi

exit $fail
