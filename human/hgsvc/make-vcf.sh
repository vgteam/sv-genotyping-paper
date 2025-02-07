#!/bin/bash

if [ ! -f hg38.fa.gz ]
then
	 wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	 gzip -d hg38.fa
	 bgzip -@ 4 hg38.fa
fi

# Make individual VCF for each sample (required for do-by-add.sh)
for SAMPLE in HG00514 HG00733 NA19240
do
	 # Make diploid VCF out of haploid vcfs:
	 echo "Making diploid VCF for ${SAMPLE}"
	 #	 ./hap-merge.py ${SAMPLE}.hap0.vcf.gz ${SAMPLE}.hap1.vcf.gz | vcfuniq | vcfkeepinfo - NA | vcffixup - | bcftools norm - -f hg38.fa.gz --multiallelic +both | vcfuniq | bgzip > HGSVC.${SAMPLE}.vcf.gz
	 # This is the line used for the analysis for the paper.  If repeating, the above would be preferred to have a slightly cleaner VCF (it adds bcftools norm)
	 ./hap-merge.py ${SAMPLE}.hap0.vcf.gz ${SAMPLE}.hap1.vcf.gz | sed 's/HG005733/HG00733/g' | vcfuniq | vcfkeepinfo - NA | vcffixup - | vcfuniq | bgzip > HGSVC.${SAMPLE}.vcf.gz
	 tabix -f -p vcf HGSVC.${SAMPLE}.vcf.gz
done

# Merge the samples
echo "Making merged VCF for all samples"
#bcftools merge -0 HGSVC.HG00514.vcf.gz HGSVC.HG00733.vcf.gz HGSVC.NA19240.vcf.gz | sed -e 's/0\/0/0\|0/g' | bcftools norm - -f hg38.fa.gz --multiallelic +both | bgzip > HGSVC.haps.vcf.gz
# This is the line used for the analysis for the paper.  If repeating, the above would be preferred to have a slightly cleaner VCF (it adds bcftools norm)
bcftools merge -0 HGSVC.HG00514.vcf.gz HGSVC.HG00733.vcf.gz HGSVC.NA19240.vcf.gz | sed -e 's/0\/0/0\|0/g' | bgzip > HGSVC.haps.vcf.gz


tabix -f -p vcf HGSVC.haps.vcf.gz



