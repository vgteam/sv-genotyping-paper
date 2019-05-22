#!/bin/bash

## Simulates reads for each sample using the graph containing the true SVs.

COV=$1
THREADS=2

## Cleanup separate folder for the results
# rm -rf cov$COV
mkdir -p cov$COV
cd cov$COV

## Number of 2x150bp reads to simulate PER HAPLOTYPE
NREADS=`sed 1d ../ref.fa | awk -v cov=$COV '{tot=tot+length($0)}END{print int(cov*tot/(2*2*150))}'`
FASTQ=../U5a_AGTCAA_L002_R1_007.fastq.gz

## Create graphs with sample threads and simulate reads
for samp in s0 s1 s2
do
    echo `date` '  ' $samp
    for hap in 0 1
    do
	vg mod -D ../truth.vg > thread.merge.vg
	vg paths --gbwt ../truth.gbwt --extract-vg -x ../truth.xg -q _thread_${samp}_chr1_${hap} >> thread.merge.vg
	vg mod -N thread.merge.vg > thread_${samp}_${hap}.vg
	vg index -t $THREADS -x thread_${samp}_${hap}.xg thread_${samp}_${hap}.vg
    done
    ## Simulate fastq files
    vg sim -x thread_${samp}_0.xg -l 100 -n $NREADS -F $FASTQ -a -p 400 -v 50 | vg view -X - > $samp.fastq
    vg sim -x thread_${samp}_1.xg -l 100 -n $NREADS -F $FASTQ -a -p 400 -v 50 | vg view -X - >> $samp.fastq
    gzip -f $samp.fastq
done

rm -f thread*
