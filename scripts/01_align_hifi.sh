#!/usr/bin/env bash
# Align PacBio HiFi reads to ARS1.2. Usage: 01_align_hifi.sh <label> <hifi.bam> <smoke|full>
set -euo pipefail
LABEL="$1"; READS="$2"; MODE="${3:-smoke}"
REF=/data/db/reference/goat_ARS1.fna
OUT=/data/bbg_pangenome/01_hifi_aln; THREADS=32
[[ "$LABEL" =~ ^[A-Za-z0-9_]+$ ]] || { echo "ERROR: bad label" >&2; exit 1; }
[[ "$MODE" == smoke || "$MODE" == full ]] || { echo "ERROR: mode must be smoke or full" >&2; exit 1; }
[[ -r "$READS" && -r "$REF" ]]   || { echo "ERROR: input not readable" >&2; exit 1; }
for t in minimap2 samtools; do command -v "$t" >/dev/null || { echo "ERROR: $t missing" >&2; exit 1; }; done
mkdir -p "$OUT"; BAM="$OUT/${LABEL}.${MODE}.ars1.bam"
[[ -e "$BAM" ]] && { echo "ERROR: $BAM exists, refusing to overwrite" >&2; exit 1; }

if [[ "$MODE" == smoke ]]; then
  SRC() { samtools fastq -@4 "$READS" 2>/dev/null | head -n 80000 || true; }
else
  SRC() { samtools fastq -@4 "$READS"; }
fi
SRC | minimap2 -ax map-hifi -t "$THREADS" --MD -Y -R "@RG\tID:${LABEL}\tSM:${LABEL}" "$REF" - \
    | samtools sort -@8 -m 2G -T "$OUT/${LABEL}.tmp" -o "$BAM" -
samtools index -@8 "$BAM"

# Checks: sorted, indexed, and at least 95% of reads mapped
SORTED=$(samtools view -H "$BAM" | grep -c 'SO:coordinate' || true); [[ "$SORTED" -ge 1 ]] || { echo "FAIL: not sorted"; exit 1; }
PCT=$(samtools flagstat "$BAM" | awk -F'[(%]' '/primary mapped/{print $2}')
awk -v p="$PCT" 'BEGIN{exit !(p>=95)}' || { echo "FAIL: only ${PCT}% mapped"; exit 1; }
echo "PASS: $BAM sorted, indexed, ${PCT}% primary mapped"
