#!/usr/bin/env bash
# HiFi large events (DEL/DUP >=10 kb, D01) vs 803M short-read CNV calls.
# RULES (stated before running): hit = same type and >=50% reciprocal overlap; sets: Sniffles PASS, and the subset also found by pbsv (same rule).
# Reported by type and size class (10-30, 30-100, >=100 kb) with Wilson 95% CI. Reverse direction is a descriptive extra.
set -euo pipefail
exec > >(tee -a /data/bbg_pangenome/28_hifi_vs_sr.log) 2>&1
CMP=/opt/conda/envs/bbg_cmp/bin
W=/data/bbg_pangenome/28_hifi_vs_sr; mkdir -p $W
SV=/data/bbg_pangenome/03_sv
ext() { # vcf -> bed (chrom start end type id): PASS or '.', DEL/DUP, END-POS >= 10 kb
  $CMP/bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%FILTER\t%ID\n' "$1" |
  awk -F'\t' -v OFS='\t' '($5=="PASS"||$5==".") && ($3 ~ /^[0-9]+$/) && ($4=="DEL"||$4=="DUP") { l=$3-$2; if (l>=10000) print $1,$2-1,$3,$4,$6 }' | sort -k1,1 -k2,2n
}
ext $SV/D01.vcf.gz > $W/sniffles.bed
ext $SV/concordance/D01.pbsv.vcf.gz > $W/pbsv.bed
awk -F'\t' -v OFS='\t' '$1=="BBGOWN-803-M" { t=($7=="deletion"?"DEL":($7=="duplication"?"DUP":"?")); if (t=="?") { print "bad type " $7 > "/dev/stderr"; exit 1 } print $4,$5-1,$6,t,$14 }' /data/bbg_comparative/cnvr_final/calls_annotated.tsv | sort -k1,1 -k2,2n > $W/sr803.bed
echo "sniffles PASS DEL/DUP>=10kb: $(wc -l < $W/sniffles.bed) ; pbsv: $(wc -l < $W/pbsv.bed) ; short-read 803M calls: $(wc -l < $W/sr803.bed)"
hit() { $CMP/bedtools intersect -a "$1" -b "$2" -f 0.5 -r -wa -wb | awk -F'\t' -v OFS='\t' '$4==$9 {print $1,$2,$3,$4,$5}' | sort -u; }
hit $W/sniffles.bed $W/pbsv.bed > $W/sniffles_with_pbsv.bed
hit $W/sniffles.bed $W/sr803.bed > $W/sniffles_with_sr.bed
hit $W/sr803.bed $W/sniffles.bed > $W/sr_with_sniffles.bed
echo "sniffles events also in pbsv: $(wc -l < $W/sniffles_with_pbsv.bed)"
$CMP/python - <<'PY'
import math
W='/data/bbg_pangenome/28_hifi_vs_sr'
def rd(f): return [tuple(l.rstrip('\n').split('\t')) for l in open(f)]
sn=rd(W+'/sniffles.bed'); two=set(rd(W+'/sniffles_with_pbsv.bed')); found=set(rd(W+'/sniffles_with_sr.bed'))
sr=rd(W+'/sr803.bed'); srf=set(rd(W+'/sr_with_sniffles.bed'))
def wil(k,n,z=1.96):
    if n==0: return (0,0)
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return 100*(c-h),100*(c+h)
def cls(e):
    l=int(e[2])-int(e[1]); return '10-30kb' if l<30000 else ('30-100kb' if l<100000 else '>=100kb')
def row(lab,ev,fd):
    n=len(ev); k=sum(1 for e in ev if e in fd); lo,hi=wil(k,n)
    print('%-34s n=%3d  found: %3d (%.1f%%, %.1f-%.1f)'%(lab,n,k,100*k/n if n else 0,lo,hi))
print('\nA) HiFi events found by 803M short-read calls')
for name,S in (('Sniffles PASS',sn),('Sniffles+pbsv',[e for e in sn if e in two])):
    row(name+' all',S,found)
    for t in ('DEL','DUP'): row(name+' '+t,[e for e in S if e[3]==t],found)
    for c in ('10-30kb','30-100kb','>=100kb'): row(name+' '+c,[e for e in S if cls(e)==c],found)
print('\nB) descriptive: short-read 803M calls with a same-type Sniffles event')
for t in ('DEL','DUP'): row('short-read '+t,[e for e in sr if e[3]==t],srf)
PY
