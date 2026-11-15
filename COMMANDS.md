# Steps run by hand

Variables: `REF` = ARS1.2 FASTA (GCF_001704415.2); `A` = HiFi alignments from scripts/01_align_hifi.sh; `V` = Sniffles2 output directory; `C` = concordance directory; `s` = D01 (803M) or C01 (D863F).

## pbmm2 and pbsv (pbmm2 26.2.99, pbsv 2.11.0)
    pbmm2 align "$REF" <hifi_reads.bam> ${s}_vs_ARS1.2.aligned.bam --preset CCS --sort -j 16
    pbsv discover ${s}_vs_ARS1.2.aligned.bam ${s}_vs_ARS1.2.svsig.gz
    pbsv call "$REF" ${s}_vs_ARS1.2.svsig.gz ${s}_vs_ARS1.2.pbsv.vcf -j 16
The 803M files carry the prefix `803m_`. The VCFs were bgzipped, indexed and stored as `$C/D01.pbsv.vcf.gz` and `$C/C01.pbsv.vcf.gz`.

## Sniffles2 (2.8.1)
    sniffles --input $A/$s.full.ars1.bam --snf $V/$s.snf --vcf $V/$s.vcf.gz \
             --reference $REF --threads 20 --minsvlen 50 --sample-id $s
    sniffles --input $V/D01.snf $V/C01.snf --vcf $V/BBG_joint.vcf.gz --reference $REF --threads 20

## Two-caller consensus (Truvari 5.5.0, bcftools 1.25)
    truvari bench -b $C/$s.pbsv.vcf.gz -c $V/$s.vcf.gz -o $C/$s.bench \
        --passonly --sizemin 50 --sizefilt 50 --sizemax 100000 \
        --pctseq 0 --pctsize 0.7 --refdist 1000 --dup-to-ins
    bcftools merge -m none --force-samples -O z -o $C/BBG_hiconf.vcf.gz \
        $C/D01.bench/tp-comp.vcf.gz $C/C01.bench/tp-comp.vcf.gz
    truvari collapse -i $C/BBG_hiconf.vcf.gz -o $C/BBG_hiconf.collapsed.vcf --sizemin 50 \
        --refdist 1000 --pctsize 0.7 --pctseq 0 -c $C/BBG_hiconf.removed.vcf
These collapse settings reproduce the released set (38,632 retained calls).
