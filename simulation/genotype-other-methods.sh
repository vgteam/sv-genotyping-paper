#!/bin/bash

## Genotypes each sample with Delly, SVTyper, and BayesTyper.

COV=$1
cd cov$COV
WD=`pwd`

THREADS=2

for samp in s0 s1 s2
do
    echo `date` $samp
    ## kmers from reads for BayesTyper
    kmc -k55 -ci1 -m8 $samp.fastq.gz $samp .
    ~/build/bayesTyper_v1.5_beta/bin/bayesTyperTools makeBloom -k $samp -p $THREADS
    echo -e "$samp\tF\t$WD/$samp" > samples.$samp.tsv
    ## Input vcf with/without errors
    for vcf in truth calls
    do
	echo `date` '  ' $vcf
	## Map reads with bwa
	bgzip -dc $samp.fastq.gz - | sed 's/fragment_\([0-9]*\)_[0-9]/fragment_\1/g' | bgzip > ${samp}_bwa.fastq.gz
	bwa mem ../ref.fa ${samp}_bwa.fastq.gz -p -R "@RG\tID:$samp\tSM:$samp" -t $THREADS | samtools sort - > $samp.bam
	samtools index $samp.bam
	## svtyper
	svtyper -i ../$vcf.symb.svtyper.vcf -B $samp.bam -l $samp.bam.json > $vcf.$samp.svtyper.vcf
	## Delly
	delly call -g ../ref.fa -v ../$vcf.vcf -o temp.bcf $samp.bam
	bcftools convert temp.bcf > $vcf.$samp.delly.vcf
	## Delly inversions
	delly call -g ../ref.fa -v ../$vcf.symb.inv.vcf -o temp.bcf $samp.bam
	bcftools convert temp.bcf > $vcf.$samp.delly.inv.vcf
	## BayesTyper
	# rm -rf bayestyper_* *kmc* *bloom* $vcf.$samp.bayestyper.vcf $samp.snv.vcf
	~/build/bayesTyper_v1.5_beta/bin/bayesTyper cluster -v ../$vcf.vcf -s samples.$samp.tsv -g ../ref.fa -p $THREADS
	~/build/bayesTyper_v1.5_beta/bin/bayesTyper genotype -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_cluster_data -s samples.$samp.tsv -g ../ref.fa -o bayestyper_unit_1/bayestyper -z -p $THREADS --min-genotype-posterior 0
	zcat bayestyper_unit_1/bayestyper.vcf.gz > $vcf.$samp.bayestyper.vcf
	mv bayestyper_unit_1 bayestyper_unit_1.$vcf.$samp
	mv bayestyper_cluster_data bayestyper_cluster_data.$vcf.$samp
    done
done
