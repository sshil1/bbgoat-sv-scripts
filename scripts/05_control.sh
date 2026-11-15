#!/usr/bin/env bash
set -euo pipefail
E=/opt/conda/envs/bbg_cmp/bin; W=/data/bbg_pangenome/05_novelty
for S in 50000 120000; do
  export S
  $E/bcftools view $W/ours.vcf.gz \
    | perl -F'\t' -ane 'if(/^#/){print;next} $F[1]+=$ENV{S}; $F[7]=~s/END=(\d+)/"END=".($1+$ENV{S})/e; print join("\t",@F)' \
    | $E/bgzip -c > $W/shift$S.vcf.gz
  $E/tabix -f -p vcf $W/shift$S.vcf.gz
  rm -rf $W/bench_shift$S
  $E/truvari bench -b $W/db8.vcf.gz -c $W/shift$S.vcf.gz -o $W/bench_shift$S --sizemin 50 --sizefilt 50 \
    --sizemax 100000 --pctseq 0 --pctsize 0.7 --refdist 500 --dup-to-ins > $W/bench_shift$S.log 2>&1 \
    || { echo "truvari failed (shift $S):"; tail -5 $W/bench_shift$S.log; exit 1; }
  python3 - $S $W/bench_shift$S/summary.json <<'PY'
import json, sys
d = json.load(open(sys.argv[2])); n = d["TP-comp"] + d["FP"]
print(f"shift {sys.argv[1]}: ours={n}  matched_by_chance={d['TP-comp']} ({100*d['TP-comp']/n:.1f}%)")
PY
done
