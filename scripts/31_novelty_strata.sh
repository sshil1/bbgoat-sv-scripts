#!/usr/bin/env bash
# Novelty vs the 8-breed catalogue by SV type and size class (descriptive). Strict bench + 50-kb-shift control.
set -euo pipefail
exec > >(tee -a /data/bbg_pangenome/31_novelty_strata.log) 2>&1
/opt/conda/envs/bbg_cmp/bin/python - <<'PY'
import subprocess, math, collections
B='/opt/conda/envs/bbg_cmp/bin/bcftools'; N='/data/bbg_pangenome/05_novelty'
def load(vcf):
    out=subprocess.check_output([B,'query','-f','%INFO/SVTYPE\t%INFO/SVLEN\n',vcf]).decode().splitlines()
    res=[]; miss=0
    for l in out:
        t,s=l.split('\t')
        try: n=abs(int(float(s)))
        except ValueError: miss+=1; continue
        res.append((t,n))
    return res,miss
bins=[(0,50,'<50'),(50,100,'50-99'),(100,500,'100-499'),(500,1000,'500-999'),(1000,5000,'1-5kb'),(5000,50000,'5-50kb'),(50000,10**12,'>=50kb')]
def bn(n):
    for lo,hi,name in bins:
        if lo<=n<hi: return name
stp,m1=load(N+'/bench/tp-comp.vcf.gz'); sfp,m2=load(N+'/bench/fp.vcf.gz')
htp,m3=load(N+'/bench_shift50000/tp-comp.vcf.gz'); hfp,m4=load(N+'/bench_shift50000/fp.vcf.gz')
print('strict: matched %d, novel %d, total %d (records without SVLEN: %d, %d)'%(len(stp),len(sfp),len(stp)+len(sfp),m1,m2))
print('shift50kb: matched %d, novel %d, total %d (without SVLEN: %d, %d)'%(len(htp),len(hfp),len(htp)+len(hfp),m3,m4))
cnt=collections.defaultdict(lambda:[0,0,0,0])
for t,n in stp: cnt[(t,bn(n))][0]+=1
for t,n in sfp: cnt[(t,bn(n))][1]+=1
for t,n in htp: cnt[(t,bn(n))][2]+=1
for t,n in hfp: cnt[(t,bn(n))][3]+=1
def wil(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return 100*(c-h),100*(c+h)
order={name:i for i,(lo,hi,name) in enumerate(bins)}
types=sorted({k[0] for k in cnt})
print('\n%-6s %-8s %6s %6s  %-26s %s'%('type','size','n','novel','% novel (95% CI)','shift50kb matched %'))
for t in types:
    tot=[sum(cnt[k][i] for k in cnt if k[0]==t) for i in range(4)]
    rows=[(k[1],cnt[k]) for k in cnt if k[0]==t]+[('ALL',tot)]
    for name,c in sorted(rows,key=lambda r:(order.get(r[0],99))):
        n=c[0]+c[1]
        if n<30 and name!='ALL': continue
        lo,hi=wil(c[1],n); sn=c[2]+c[3]
        print('%-6s %-8s %6d %6d  %5.1f (%4.1f-%4.1f)           %5.1f'%(t,name,n,c[1],100*c[1]/n,lo,hi,100*c[2]/sn if sn else 0))
PY
