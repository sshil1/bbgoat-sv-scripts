#!/usr/bin/env bash
# 32_asm_concordance.sh - assembly-based check of the two-caller SV set (DEL/INS, autosomes)
set -Eeuo pipefail
umask 027
trap 'echo "[$(date +%T)] FAILED at line $LINENO" >&2' ERR

OUT=/data/bbg_pangenome/32_asm
REF=/data/db/reference/goat_ARS1.fna
SET=/data/bbg_pangenome/03_sv/concordance/BBG_hiconf.collapsed.vcf
THREADS=${THREADS:-24}
CONDA=/opt/conda/bin/conda
SVENV=${SVENV:-/opt/conda/envs/svimasm}
DIRD=/data/gmbluser/pacbio_data/r84286_20260723_204918_1_D01
DIRC=/data/gmbluser/pacbio_data/r84286_20260723_204918_1_C01
mkdir -p "$OUT"; cd "$OUT"
log()   { echo "[$(date +%T)] $*"; }
done_() { [ -e "$OUT/.done_$1" ]; }
mark()  { touch "$OUT/.done_$1"; }

for t in minimap2 samtools bcftools bgzip tabix awk python3; do
  command -v "$t" >/dev/null || { echo "missing tool: $t"; exit 1; }
done
TRU=$(command -v truvari || true)
if [ -z "$TRU" ]; then TRU=$(ls /opt/conda/envs/*/bin/truvari 2>/dev/null | head -n1 || true); fi
[ -n "$TRU" ] || { echo "truvari not found"; exit 1; }
log "truvari: $TRU"

if [ ! -x "$SVENV/bin/svim-asm" ]; then
  log "creating env svimasm"
  "$CONDA" create -y -n svimasm python=3.11
  "$SVENV/bin/pip" install svim-asm
fi
"$SVENV/bin/svim-asm" --version || true

AUTOS=$(printf 'NC_0308%02d.1,' $(seq 8 36)); AUTOS=${AUTOS%,}

declare -A GFA=(
  [D01_hap1]="$DIRD/black_bengal_goat.bp.hap1.p_ctg.gfa"
  [D01_hap2]="$DIRD/black_bengal_goat.bp.hap2.p_ctg.gfa"
  [C01_hap1]="$DIRC/black_bengal_goat_C01.bp.hap1.p_ctg.gfa"
  [C01_hap2]="$DIRC/black_bengal_goat_C01.bp.hap2.p_ctg.gfa"
)
mkfa() {
  awk '$1=="S"{print ">"$2; print $3}' "$1" > "$2"
  [ "$(stat -c %s "$2")" -gt 1000000000 ] || { echo "fasta too small: $2"; exit 1; }
}

# 1. haplotype contigs -> ARS1.2
for k in D01_hap1 D01_hap2 C01_hap1 C01_hap2; do
  if done_ "align_$k"; then log "skip align $k"; continue; fi
  [ -s "${GFA[$k]}" ] || { echo "missing ${GFA[$k]}"; exit 1; }
  log "gfa->fasta $k"; mkfa "${GFA[$k]}" "$k.fa"
  log "align $k"
  minimap2 -t "$THREADS" -ax asm5 --cs -r2k "$REF" "$k.fa" 2> "$k.minimap2.log" \
    | samtools sort -@ 8 -T "$OUT/tmp_$k" -o "$k.bam" -
  samtools index "$k.bam"
  [ "$(samtools view -c -F 0x904 "$k.bam")" -gt 100 ] || { echo "too few alignments: $k"; exit 1; }
  rm -f "$k.fa"
  mark "align_$k"
done

# 2. our per-animal calls
if ! done_ setprep; then
  bcftools view -Oz -o set.vcf.gz "$SET"; tabix -f -p vcf set.vcf.gz
  bcftools query -l set.vcf.gz > set.samples.txt
  mark setprep
fi
log "samples in set: $(tr '\n' ' ' < set.samples.txt)"

HELP=$("$SVENV/bin/svim-asm" diploid --help 2>&1 || true)
EXTRA=()
if grep -q -- '--min_sv_size' <<<"$HELP"; then EXTRA+=(--min_sv_size 50); fi
if grep -q -- '--max_sv_size' <<<"$HELP"; then EXTRA+=(--max_sv_size 100000); fi

# 3. SVIM-asm diploid + truvari bench per animal
for A in D01 C01; do
  if done_ "bench_$A"; then log "skip bench $A"; continue; fi
  S=$(grep -m1 -i "$A" set.samples.txt || true)
  [ -n "$S" ] || { echo "no sample matching $A in set.samples.txt"; exit 1; }
  log "svim-asm $A (sample $S)"
  rm -rf "svim_$A"
  "$SVENV/bin/svim-asm" diploid ${EXTRA[@]+"${EXTRA[@]}"} "svim_$A" "${A}_hap1.bam" "${A}_hap2.bam" "$REF" > "svim_$A.log" 2>&1
  bcftools sort -Oz -o "svim_$A.all.vcf.gz" "svim_$A/variants.vcf"; tabix -f -p vcf "svim_$A.all.vcf.gz"
  bcftools view -r "$AUTOS" -i 'SVTYPE="DEL"||SVTYPE="INS"' -Oz -o "asm_$A.vcf.gz" "svim_$A.all.vcf.gz"; tabix -f -p vcf "asm_$A.vcf.gz"
  bcftools view -s "$S" -c1 -r "$AUTOS" -i 'SVTYPE="DEL"||SVTYPE="INS"' -Oz -o "ours_$A.vcf.gz" set.vcf.gz; tabix -f -p vcf "ours_$A.vcf.gz"
  [ "$(bcftools view -H "asm_$A.vcf.gz" | wc -l)" -gt 1000 ] || { echo "too few assembly calls: $A"; exit 1; }
  log "truvari bench $A"
  rm -rf "bench_$A"
  "$TRU" bench -b "asm_$A.vcf.gz" -c "ours_$A.vcf.gz" -o "bench_$A" \
    --sizemin 50 --sizemax 100000 --pctsize 0.7 --pctseq 0 --refdist 1000 --dup-to-ins > "bench_$A.log" 2>&1
  mark "bench_$A"
done

# 4. summaries
printf 'animal\tasm_calls\tour_calls\tTP_base\tTP_comp\tFP\tFN\tprecision\trecall\tF1\n' > summary.tsv
for A in D01 C01; do
python3 - "bench_$A/summary.json" "$A" <<'PYEOF' >> summary.tsv
import json, sys
s = json.load(open(sys.argv[1]))
assert s["TP-base"] + s["FN"] == s["base cnt"], "TP-base + FN != base cnt"
f = [sys.argv[2], s["base cnt"], s["comp cnt"], s["TP-base"], s["TP-comp"], s["FP"], s["FN"],
     round(s["precision"], 3), round(s["recall"], 3), round(s["f1"], 3)]
print("\t".join(str(x) for x in f))
PYEOF
done
strat() {
  bcftools query -f '%SVLEN\n' "$1" | awk -v L="$2" '{v=$1; if(v==".")next; if(v<0)v=-v;
    b=(v<100)?"50-99":(v<500)?"100-499":(v<5000)?"500-4999":"5000-100000"; n[b]++}
    END{for(b in n) printf "%s\t%s\t%d\n", L, b, n[b]}'
}
: > strata.tsv
for A in D01 C01; do
  for f in tp-base fn tp-comp fp; do strat "bench_$A/$f.vcf.gz" "${A}_$f" >> strata.tsv; done
done
sort -o strata.tsv strata.tsv
log "SUMMARY"; cat summary.tsv; log "STRATA (animal_class, size bin, count)"; cat strata.tsv
touch "$OUT/ALLDONE"; log "ALLDONE"
