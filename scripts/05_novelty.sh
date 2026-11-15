#!/usr/bin/env bash
set -euo pipefail
E=/opt/conda/envs/bbg_cmp/bin
REF=/data/db/reference/goat_ARS1.fna
INFO=/data/bbg_pangenome/04_public/SV_extracted/SV/non-redundantSV.info
HI=/data/bbg_pangenome/03_sv/concordance/BBG_hiconf.collapsed.vcf
W=/data/bbg_pangenome/05_novelty; mkdir -p $W
{
  echo '##fileformat=VCFv4.2'
  awk '{print "##contig=<ID="$1",length="$2">"}' ${REF}.fai
  echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">'
  echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">'
  echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End">'
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDB\n'
  awk -F'\t' -v OFS='\t' 'NR>1 {
    c=$1; if (c=="X1") ch="NW_017189516.1"; else if (c=="X2") ch="NW_017189517.1"; else ch=sprintf("NC_0%d.1", c+30807)
    if ($5 ~ /^DEL\./ || $5 ~ /-DEL-/) t="DEL"; else if ($5 ~ /^INS\./ || $5 ~ /-INS-/) t="INS"; else next
    l=length($4)-length($3); if (l<50 && l>-50) next
    print ch, $2, "db" NR, toupper($3), toupper($4), ".", "PASS", "SVTYPE=" t ";SVLEN=" l ";END=" ($2+length($3)-1), "GT", "0/1" }' $INFO
} > $W/db8.unsorted.vcf
$E/bcftools sort -O z -o $W/db8.vcf.gz $W/db8.unsorted.vcf && $E/tabix -f -p vcf $W/db8.vcf.gz
printf 'catalogue DEL/INS entries >=50bp\t'; $E/bcftools view -H $W/db8.vcf.gz | wc -l
$E/bcftools view -s D01 -t '^NW_017189516.1,NW_017189517.1,NW_017189518.1' -O z -o $W/ours.vcf.gz $HI
$E/tabix -f -p vcf $W/ours.vcf.gz
printf 'our autosomal high-confidence SVs\t'; $E/bcftools view -H $W/ours.vcf.gz | wc -l
rm -rf -- "${W:?}/bench"
$E/truvari bench -b $W/db8.vcf.gz -c $W/ours.vcf.gz -o $W/bench --sizemin 50 --sizefilt 50 \
  --sizemax 100000 --pctseq 0 --pctsize 0.7 --refdist 500 --dup-to-ins > $W/bench.log 2>&1 \
  || { echo "truvari failed:"; tail -6 $W/bench.log; exit 1; }
python3 - $W/bench/summary.json <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
n = d["TP-comp"] + d["FP"]
print(f"ours={n}\tin_8breed_catalogue={d['TP-comp']} ({100*d['TP-comp']/n:.1f}%)\tnot_in_catalogue={d['FP']} ({100*d['FP']/n:.1f}%)")
PY
for f in tp-comp fp; do echo "--- $f by type ---"; $E/bcftools query -f '%INFO/SVTYPE\n' $W/bench/$f.vcf.gz | sort | uniq -c; done
