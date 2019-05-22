#!/bin/bash

## Genotypes each sample using toil-vg.

COV=$1
cd cov$COV

THREADS=2

for samp in s0 s1 s2
do
    echo `date` '  ' $samp
    for graph in calls truth calls.symb truth.symb
    do
	echo `date` '    ' $graph
	## Map reads
	vg map -x ../$graph.xg -g ../$graph.gcsa -f $samp.fastq.gz -i -N $samp -t $THREADS > $graph.$samp.gam
	rm -fr $graph.$samp.gam.index
	vg index -t $THREADS -d $graph.$samp.gam.index -N $graph.$samp.gam
	## toil-vg recall
	toil-vg call js ../$graph.xg $samp . --gams $graph.$samp.gam --chroms chr1 --recall --logFile $graph.$samp.log --defaultCores $THREADS
	mv $samp.vcf.gz $graph.$samp.toilvg.vcf.gz
	mv $samp.vcf.gz.tbi $graph.$samp.toilvg.vcf.gz.tbi
    done
done
