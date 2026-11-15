#!/usr/bin/env bash
# Last slice for item 1 (post hoc): HiFi depth support of Sniffles events by size class x found-by-short-read (strict). Reads outputs of scripts 28/29.
set -euo pipefail
exec > >(tee -a /data/bbg_pangenome/30_crosstab.log) 2>&1
/opt/conda/envs/bbg_cmp/bin/python - <<'PY'
import math, statistics as st, csv
W='/data/bbg_pangenome/28_hifi_vs_sr'
def rd(f): return [tuple(l.rstrip('\n').split('\t')) for l in open(f)]
sn=rd(W+'/sniffles.bed'); two=set(rd(W+'/sniffles_with_pbsv.bed')); strict=set(rd(W+'/sniffles_with_sr.bed'))
cov={}
for l in open(W+'/ev.bedcov'):
    c,s,e,n,t=l.rstrip('\n').split('\t'); i,p=n.split('|'); cov.setdefault(int(i),{})[p]=float(t)/(int(e)-int(s))
fl=[]
with open('/data/bbg_pangenome/23_cnvr_hifi/per_cnvr.tsv') as f:
    r=csv.reader(f,delimiter='\t'); next(r)
    for row in r: fl.append(float(row[5]))
med=st.median(fl)
rows=[]
for i,e in enumerate(sn):
    v=cov.get(i)
    if not v or len(v)<3: continue
    f=(v['L']+v['R'])/2
    if f<0.4*med: continue
    ratio=v['int']/f
    ok=ratio<=0.75 if e[3]=='DEL' else ratio>=1.25
    rows.append((e,ok))
def wil(k,n,z=1.96):
    if n==0: return (0,0)
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return 100*(c-h),100*(c+h)
def line(lab,sel):
    k=sum(1 for r in sel if r[1]); n=len(sel); lo,hi=wil(k,n)
    print('%-44s depth-supported %3d/%3d (%.1f%%, %.1f-%.1f)'%(lab,k,n,100*k/n if n else 0,lo,hi))
size=lambda e: '10-100kb' if int(e[2])-int(e[1])<100000 else '>=100kb'
for g in ('10-100kb','>=100kb'):
    G=[r for r in rows if size(r[0])==g]
    print('\n'+g+' (informative events: %d)'%len(G))
    line('  found by short-read (strict)',[r for r in G if r[0] in strict])
    line('  NOT found by short-read (strict)',[r for r in G if r[0] not in strict])
    line('    ... unfound, two-caller',[r for r in G if r[0] not in strict and r[0] in two])
    line('    ... unfound, Sniffles-only',[r for r in G if r[0] not in strict and r[0] not in two])
PY
