#!/usr/bin/env bash
# Sensitivity of the "no match in the eight-breed catalogue" count to the Truvari matching settings.
# Reuses db8.vcf.gz and ours.vcf.gz written by 05_novelty.sh. Each setting runs in its own new directory;
# an existing directory is skipped, never deleted. Run inside a dedicated tmux window.
set -euo pipefail
E="${E:-/opt/conda/envs/bbg_cmp/bin}"
W="${W:-/data/bbg_pangenome/05_novelty}"
S="$W/sens"
LOG="${LOG:-/data/bbg_pangenome/34_novelty_sensitivity.log}"
exec > >(tee -a "$LOG") 2>&1

for f in "$W/db8.vcf.gz" "$W/ours.vcf.gz"; do
  [ -s "$f" ] || { echo "missing input: $f (run 05_novelty.sh first)"; exit 1; }
done
mkdir -p "$S"

# tag refdist pctsize pctseq   (baseline is the setting used for the reported 6,465 / 6,212)
CONFIGS="
baseline 500 0.7 0
refdist200 200 0.7 0
refdist1000 1000 0.7 0
refdist2000 2000 0.7 0
pctsize0.5 500 0.5 0
pctsize0.9 500 0.9 0
pctseq0.7 500 0.7 0.7
loose 2000 0.5 0
strict 200 0.9 0.7
"

count() { # $1 = vcf.gz ; prints "<DEL+INS count> <all count>"
  "$E/bcftools" query -f '%INFO/SVTYPE\n' "$1" | awk '{a++; if ($1=="DEL"||$1=="INS") d++} END {printf "%d %d\n", d+0, a+0}'
}

printf 'tag\trefdist\tpctsize\tpctseq\tmatched_delins\tnovel_delins\tpct_novel_delins\tmatched_all\tnovel_all\tpct_novel_all\n'
echo "$CONFIGS" | while read -r tag rd ps pq; do
  [ -n "${tag:-}" ] || continue
  out="$S/run_${tag}"
  if [ -e "$out" ]; then
    echo "reusing existing $out (not deleted, not rerun)" >&2
  else
    "$E/truvari" bench -b "$W/db8.vcf.gz" -c "$W/ours.vcf.gz" -o "$out" --sizemin 50 --sizefilt 50 \
      --sizemax 100000 --pctseq "$pq" --pctsize "$ps" --refdist "$rd" --chunksize "$(( rd > 1000 ? rd : 1000 ))" --dup-to-ins > "$S/${tag}.log" 2>&1 \
      || { echo "truvari failed for $tag:" >&2; tail -6 "$S/${tag}.log" >&2; continue; }
  fi
  read -r md ma < <(count "$out/tp-comp.vcf.gz")
  read -r nd na < <(count "$out/fp.vcf.gz")
  awk -v t="$tag" -v rd="$rd" -v ps="$ps" -v pq="$pq" -v md="$md" -v nd="$nd" -v ma="$ma" -v na="$na" 'BEGIN {
    printf "%s\t%s\t%s\t%s\t%d\t%d\t%.1f\t%d\t%d\t%.1f\n", t, rd, ps, pq, md, nd, 100*nd/(md+nd), ma, na, 100*na/(ma+na) }'
  if [ "$tag" = baseline ] && [ "$na" -ne 6465 ]; then
    echo "WARNING: baseline novel_all=$na but the reported value is 6465; check inputs before using this table" >&2
  fi
done
