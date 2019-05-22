#!/bin/bash

## Build BWA indexes and reformat VCFs.

## bwa index
bwa index ref.fa

## VCF for svtyper and delly
for vcf in calls truth
do
    grep -v "SVTYPE=INS" $vcf.symb.vcf | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8}' > $vcf.symb.svtyper.vcf
    grep -v "SVTYPE=DEL" $vcf.symb.vcf | grep -v "SVTYPE=INS" > $vcf.symb.inv.vcf
done
