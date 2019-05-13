# HGSVC

##  Commands used to download the data

```
for name in HG00514 HG00733 NA19240 HG00514-sim
do
aws s3 sync s3://${OUTSTORE}/HGSVC/eval-${name} ./HGSVC-jan5-eval-${name}
done
```

The merged VCF with SVs from the three samples was created using [make-vcf.sh](make-vcf.sh).
It creates a `HGSVC.haps.vcf.gz` file.

## toil-vg

### Setup

The following environment variables need to be set

```
# jobstore prefix: EX my-jobstore
export JOBSTORE=
# outstore: EX my-outstore/HGSVC
export OUTSTORE=
# cluster name prefix EX: my-cluster where my-cluster1 and my-cluster2 were already created with toil-vg/scripts/create-ec2-leader.sh
export CLUSTER=
# BED file for evaluation
export COMPARE_REGIONS_BED=
# interleaved fq template for simulation error model
export TEMPLATE_FQ=
# input reads root location (will look for reads in subdirs)
export FQBASE=
```

### Running

Using the helper scripts from `../toil-scripts`.

```
# Construct all graphs and indexes.
./construct.sh -c ${CLUSTER}1 ${JOBSTORE}1 ${OUTSTORE}

# Simulate a GAM and fastq from the two HG00514 haplotypes
./simulate.sh -c ${CLUSTER}1 ${JOBSTORE}1 ${OUTSTORE}/sim s3://${OUTSTORE}/HGSVC_HG00514_haplo_thread_0.xg s3://${OUTSTORE}/HGSVC_HG00514_haplo_thread_1.xg ${TEMPLATE_FQ}

## Mapping, Calling and Evaluation.
## Use -M SKIP, -C SKIP, -E SKIP to bypass mapping, calling, evaluation respectively (ex, if rerunning a step)

# Simulation
./mce.sh -c ${CLUSTER}1  ${JOBSTORE}1 ${OUTSTORE}/HGSVC s3://${OUTSTORE}/HGSVC/HGSVC HG00514 HG00514-sim s3://${OUTSTORE}/HGSVC/HGSVC.haps.vcf.gz ${COMPARE_REGIONS_BED} s3://${OUTSTORE}/HGSVC-chroms-dec5/sim/sim-HG00514-30x.fq.gz

# Three HGSVC Samples (the reads can be found publicly by looking up the run names on EBI's ENA.  I had them mirrored in FQBASE for faster access)
./mce.sh -c ${CLUSTER}1  ${JOBSTORE}1 ${OUTSTORE}/HGSVC s3://${OUTSTORE}/HGSVC/HGSVC HG00514 HG00514 s3://${OUTSTORE}/HGSVC/HGSVC.haps.vcf.gz ${COMPARE_REGIONS_BED} ${FQBASE}/HG00514/ERR903030_1.fastq.gz ${FQBASE}/HG00514/ERR903030_2.fastq.gz 

./mce.sh -c ${CLUSTER}2  ${JOBSTORE}2 ${OUTSTORE}/HGSVC s3://${OUTSTORE}/HGSVC/HGSVC HG00733 HG00733 s3://${OUTSTORE}/HGSVC/HGSVC.haps.vcf.gz ${COMPARE_REGIONS_BED} ${FQBASE}/HG00733/ERR895347_1.fastq.gz ${FQBASE}/HG00733/ERR895347_2.fastq.gz 

./mce.sh -c ${CLUSTER}3  ${JOBSTORE}3 ${OUTSTORE}/HGSVC s3://${OUTSTORE}/HGSVC/HGSVC NA19240 NA19240 s3://${OUTSTORE}/HGSVC/HGSVC.haps.vcf.gz ${COMPARE_REGIONS_BED} ${FQBASE}/NA19240/ERR894724_1.fastq.gz ${FQBASE}/NA19240/ERR894724_2.fastq.gz 
```

## Other methods

### Mapping reads to GRCh38

```
bwa mem -t 12 hg38_noalt.fa ERR894724_1.fq.gz ERR894724_2.fq.gz | samtools view -bh - | samtools sort -o ERR894724_bwa_mem_sort.bam -
samtools index ERR894724_bwa_mem_sort.bam
picard AddOrReplaceReadGroups I=ERR894724_bwa_mem_sort.bam O=ERR894724_bwa_mem_sort_rg.bam LB=ERR894724 PL=ILLUMINA PU=ILLUMINA SM=ERR894724
picard MarkDuplicates I=ERR894724_bwa_mem_sort_rg.bam O=ERR894724_bwa_mem_sort_rg_md.bam M=md_metrics.txt CREATE_INDEX=true TMP_DIR=.
rm ERR894724_bwa_mem_sort_rg.*
```

### Delly

For each sample:

```
VCF=HGSVC.haps.vcf.gz
BAM=XXX.bam
REF=hg38.fa
SAMP=XXX

delly call -g $REF -v $VCF -o temp.bcf $BAM
bcftools convert temp.bcf > $SAMP.delly.vcf
bgzip $SAMP.delly.vcf
tabix $SAMP.delly.vcf.gz
```

### SVTyper

The explicit VCF was first converted into a symbolic VCF, for deletions only because SVTyper cannot genotype insertions.

```
zcat HGSVC.haps.vcf.gz | python VCFtoSymbolicDEL.py > HGSVC.haps.symbolic.dels.vcf
```

Then for each sample:

```
BAM=XXX.bam
SAMP=XXX

svtyper -i HGSVC.haps.symbolic.dels.vcf -B $BAM -l libinfo.json > $SAMP.svtyper.vcf
vcf-sort $SAMP.svtyper.vcf | vcfkeepinfo - END SVLEN SVTYPE NS AN | grep -v 'contig' > $SAMP.svtyper.clean.vcf
bayesTyperTools convertAllele -v $SAMP.svtyper.clean.vcf -g ../../haps/hg38.fa --keep-imprecise -o $SAMP.svtyper.explicit
vcfkeepinfo $SAMP.svtyper.explicit.vcf NS AN | bgzip > $SAMP.svtyper.explicit.vcf.gz
tabix -f $SAMP.svtyper.explicit.vcf.gz
```

### BayesTyper

```
FASTQ=XXX.fastq
kmc -v -k55 -m128 -sm -ci1 -t12 @reads.txt ERR894724_k55 .
bayesTyperTools makeBloom -p 12 -k ERR894724_k55
```

With `reads.txt` containing two lines with the paths to the FASQ files for the sample.

#### Calling SNVs

Haplotype caller

```
gatk3 -Xmx256G -T HaplotypeCaller -R hg38_noalt.fa -nct 12 -I ERR894724_bwa_mem_sort_rg_md.bam -o ERR894724_bwa_mem_sort_md_haplotypecaller.vcf
gzip ERR894724_bwa_mem_sort_md_haplotypecaller.vcf

# Norm
VCF_PREFIX=ERR894724_bwa_mem_sort_md_haplotypecaller
GENOME_DIR=XXX
bcftools norm -c x -m -any -f ${GENOME_DIR}/hg38_noalt.fa ${VCF_PREFIX}.vcf.gz | cut -f1-8 | bcftools annotate -x FORMAT,INFO > ${VCF_PREFIX}_norm.vcf
gzip ${VCF_PREFIX}_norm.vcf

# Filter
filterStructuralVariants ERR894724_bwa_mem_sort_md_haplotypecaller_norm.vcf.gz ERR894724_bwa_mem_sort_md_haplotypecaller_norm_sv49 49 -49 49
```

Platypus

```
platypus callVariants --bamFiles=ERR894724_bwa_mem_sort.bam --refFile=hg38_noalt.fa --output=ERR894724_bwa_mem_sort_platypus.vcf --logFileName=platypus_log.txt --nCPU 12 --assemble=1 --assembleBrokenPairs=1
bgzip ERR894724_bwa_mem_sort_platypus.vcf
tabix ERR894724_bwa_mem_sort_platypus.vcf.gz

# Norm
VCF_PREFIX=ERR894724_bwa_mem_sort_platypus
GENOME_DIR=XXX
bcftools norm -c x -m -any -f ${GENOME_DIR}/hg38_noalt.fa ${VCF_PREFIX}.vcf.gz | cut -f1-8 | bcftools annotate -x FORMAT,INFO > ${VCF_PREFIX}_norm.vcf
picard UpdateVcfSequenceDictionary I=${VCF_PREFIX}_norm.vcf O=${VCF_PREFIX}_norm_header.vcf SD=${GENOME_DIR}/hg38_noalt.dict
picard SortVcf I=${VCF_PREFIX}_norm_header.vcf O=${VCF_PREFIX}_norm_sort.vcf SD=${GENOME_DIR}/hg38_noalt.dict
rm ${VCF_PREFIX}_norm.vcf
rm ${VCF_PREFIX}_norm_header.vcf
rm ${VCF_PREFIX}_norm_sort.vcf.idx
gzip ${VCF_PREFIX}_norm_sort.vcf

# Filter
filterStructuralVariants ERR894724_bwa_mem_sort_platypus_norm_sort.vcf.gz ERR894724_bwa_mem_sort_platypus_norm_sort_sv49 49 -49 49
```

Variants were then combined across all samples:

```
bayesTyperTools combine -v hgsvc:HGSVC.haps_norm_sort.vcf.gz,pp:ERR903030_bwa_mem_sort_platypus_norm_sort_sv49.vcf.gz,pp:ERR895347_bwa_mem_sort_platypus_norm_sort_sv49.vcf.gz,pp:ERR894724_bwa_mem_sort_platypus_norm_sort_sv49.vcf.gz,hc:ERR903030_bwa_mem_sort_md_haplotypecaller_norm_sv49.vcf.gz,hc:ERR895347_bwa_mem_sort_md_haplotypecaller_norm_sv49.vcf.gz,hc:ERR894724_bwa_mem_sort_md_haplotypecaller_norm_sv49.vcf.gz -o all_hgsvc_pp49_hc49 -z
```

BayesTyper genotyping

```
MAIN_DIR=XXX
OUT_PREFIX=bayestyper_unit_1/all_hgsvc_pp49_hc49_bayestyper

bayesTyper cluster -v all_hgsvc_pp49_hc49.vcf.gz -s samples.txt -g ${MAIN_DIR}/genome/hg38_canon.fa -d ${MAIN_DIR}/genome/hg38_decoy.fa -p 24 -r 12345678 --min-number-of-unit-variants 10000000
bayesTyper genotype -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_cluster_data -s samples.txt -g ${MAIN_DIR}/genome/hg38_canon.fa -d ${MAIN_DIR}/genome/hg38_decoy.fa -o ${OUT_PREFIX} -z -p 24 -r 12345678 --min-genotype-posterior 0

bcftools filter -i 'FILTER=="PASS"' ${OUT_PREFIX}.vcf.gz | gzip -c 1> ${OUT_PREFIX}_pass.vcf.gz
~/tools/filterAlleleCallsetOrigin ${OUT_PREFIX}_pass.vcf.gz ${OUT_PREFIX}_pass_hgsvc hc,pp,. 1
```

`samples.txt` contains:

```
HG00514	F	/public/groups/cgl/graph-genomes/jsibbese/hgsvc/real/HG00514/kmers/ERR903030_k55
HG00733	F	/public/groups/cgl/graph-genomes/jsibbese/hgsvc/real/HG00733/kmers/ERR895347_k55
NA19240	F	/public/groups/cgl/graph-genomes/jsibbese/hgsvc/real/NA19240/kmers/ERR894724_k55
```
