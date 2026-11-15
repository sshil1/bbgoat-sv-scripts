#!/usr/bin/env bash
# Follow-up to script 28 (added after seeing its output; post hoc).
# (i) relaxed overlap (any overlap, same type) to test whether fragmentation explains low concordance;
# (ii) HiFi depth support (script-23 rules) for each Sniffles event, by whether the short-read catalogue found it.
set -euo pipefail
exec > >(tee -a /data/bbg_pangenome/29_hifi_events_depth.log) 2>&1
S=/opt/conda/envs/bbg_pg/bin; CMP=/opt/conda/envs/bbg_cmp/bin
W=/data/bbg_pangenome/28_hifi_vs_sr; REF=/data/db/reference/goat_ARS1.fna
anyhit() { $CMP/bedtools intersect -a "$1" -b "$2" -wa -wb | awk -F'\t' -v OFS='\t' '$4==$9 {print $1,$2,$3,$4,$5}' | sort -u; }
anyhit $W/sniffles.bed $W/sr803.bed > $W/sniffles_with_sr_any.bed
anyhit $W/sr803.bed $W/sniffles.bed > $W/sr_with_sniffles_any.bed
awk -F'\t' -v OFS='\t' 'FNR==NR {len[$1]=$2; next} { i=FNR-1; if ($2-12000<0 || $3+12000>len[$1]) next; print $1,$2,$3,i"|int"; print $1,$2-12000,$2-2000,i"|L"; print $1,$3+2000,$3+12000,i"|R" }' $REF.fai $W/sniffles.bed > $W/ev_regions.bed
$S/samtools bedcov $W/ev_regions.bed /data/bbg_pangenome/01_hifi_aln/D01.full.ars1.bam > $W/ev.bedcov
$CMP/python - <<'PY'
import math, statistics as st, csv
W='/data/bbg_pangenome/28_hifi_vs_sr'
def rd(f): return [tuple(l.rstrip('\n').split('\t')) for l in open(f)]
sn=rd(W+'/sniffles.bed'); two=set(rd(W+'/sniffles_with_pbsv.bed')); strict=set(rd(W+'/sniffles_with_sr.bed'))
anyf=set(rd(W+'/sniffles_with_sr_any.bed')); sr=rd(W+'/sr803.bed'); sra=set(rd(W+'/sr_with_sniffles_any.bed'))
cov={}
for l in open(W+'/ev.bedcov'):
    c,s,e,n,t=l.rstrip('\n').split('\t'); i,p=n.split('|'); cov.setdefault(int(i),{})[p]=float(t)/(int(e)-int(s))
fl=[]
with open('/data/bbg_pangenome/23_cnvr_hifi/per_cnvr.tsv') as f:
    r=csv.reader(f,delimiter='\t'); next(r)
    for row in r: fl.append(float(row[5]))
med=st.median(fl)
def wil(k,n,z=1.96):
    if n==0: return (0,0)
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return 100*(c-h),100*(c+h)
def cls(e):
    l=int(e[2])-int(e[1]); return '10-30kb' if l<30000 else ('30-100kb' if l<100000 else '>=100kb')
def line(lab,k,n):
    lo,hi=wil(k,n); print('%-40s %3d/%3d (%.1f%%, %.1f-%.1f)'%(lab,k,n,100*k/n if n else 0,lo,hi))
print('A) relaxed overlap (any overlap, same type): HiFi events found by short-read calls')
line('Sniffles PASS all',sum(e in anyf for e in sn),len(sn))
for c in ('10-30kb','30-100kb','>=100kb'):
    S=[e for e in sn if cls(e)==c]; line('  '+c,sum(e in anyf for e in S),len(S))
line('short-read 803M calls with Sniffles event',sum(e in sra for e in sr),len(sr))
print('\nB) HiFi depth support of Sniffles events (informative only)')
rows=[]
for i,e in enumerate(sn):
    v=cov.get(i)
    if not v or len(v)<3: continue
    f=(v['L']+v['R'])/2
    if f<0.4*med: continue
    ratio=v['int']/f; ok=ratio<=0.75 if e[3]=='DEL' else ratio>=1.25
    rows.append((e,ok,ratio))
print('informative events: %d of %d'%(len(rows),len(sn)))
def grp(lab,sel): line(lab,sum(1 for r in sel if r[1]),len(sel))
grp('all events',rows)
grp('  found by short-read (strict)',[r for r in rows if r[0] in strict])
grp('  NOT found by short-read (strict)',[r for r in rows if r[0] not in strict])
grp('  two-caller, not found (strict)',[r for r in rows if r[0] in two and r[0] not in strict])
grp('  Sniffles-only, not found (strict)',[r for r in rows if r[0] not in two and r[0] not in strict])
for c in ('10-30kb','30-100kb','>=100kb'):
    grp('  '+c+' (all)',[r for r in rows if cls(r[0])==c])
PY
