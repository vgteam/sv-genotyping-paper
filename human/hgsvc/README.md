# HGSVC

## Data

- [GRCh38 reference genome](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).
- Phased VCFs for the three Human Genome Structural Variation Consortium (HGSVC) samples from [Chaisson et al. Nat Comms 2019](https://www.nature.com/articles/s41467-018-08148-z) are in this folder (names in *SAMP.hap[01].vcf.gz*).

The [make-vcf.sh](make-vcf.sh) script creates a combined VCF file `HGSVC.haps.vcf.gz` with SVs from the three samples.

### Reads

- HG00514: [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR903/ERR903030/ERR903030_1.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR903/ERR903030/ERR903030_1.fastq.gz) [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR903/ERR903030/ERR903030_2.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR903/ERR903030/ERR903030_2.fastq.gz)
- HG00733: [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_1.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_1.fastq.gz) [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_2.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_2.fastq.gz)
- NA19240: [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_1.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_1.fastq.gz) [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_2.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_2.fastq.gz)

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
./construct.sh -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}
# Simulate a GAM and fastq from the two HG00514 haplotypes
./simulate.sh -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}/sim s3://${OUTSTORE}/HGSVC_HG00514_haplo_thread_0.xg s3://${OUTSTORE}/HGSVC_HG00514_haplo_thread_1.xg ${TEMPLATE_FQ}
## Mapping, Calling and Evaluation.
## Use -M SKIP, -C SKIP, -E SKIP to bypass mapping, calling, evaluation respectively (ex, if rerunning a step)
# Simulation
./map.sh -c ${CLUSTER}  ${JOBSTORE} ${OUTSTORE}/HGSVC/map-HG00514-sim s3://${OUTSTORE}/HGSVC HG00514 s3://${OUTSTORE}/HGSVC-chroms-dec5/sim/sim-HG00514-30x.fq.gz
./call.sh -c ${CLUSTER} -v s3://${OUTSTORE}/HGSVC/HGSVC.haps.vcf.gz -l s3://${OUTSTORE}/HGSVC/HGSVC_alts.gam ${JOBSTORE} ${OUTSTORE}/HGSVC/call-HG00514-sim s3://${OUTSTORE}/HGSVC/HGSVC.xg HG00514 s3://${OUTSTORE}/HGSVC/map-HG00514-sim/HG00514_chr

# Three HGSVC Samples (the reads can be found publicly by looking up the run names on EBI's ENA.  I had them mirrored in FQBASE for faster access)
./map.sh -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}/HGSVC/map-HG00514 s3://${OUTSTORE}/HGSVC/HGSVC HG00514 ${FQBASE}/HG00514/ERR903030_1.fastq.gz ${FQBASE}/HG00514/ERR903030_2.fastq.gz 
./map.sh -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}/HGSVC/map-HG00733 s3://${OUTSTORE}/HGSVC/HGSVC HG00733 ${FQBASE}/HG00733/ERR895347_1.fastq.gz ${FQBASE}/HG00733/ERR895347_2.fastq.gz
./map.sh -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}/HGSVC/map-NA19240 s3://${OUTSTORE}/HGSVC/HGSVC NA19240 ${FQBASE}/NA19240/ERR894724_1.fastq.gz ${FQBASE}/NA19240/ERR894724_2.fastq.gz

./call.sh -c ${CLUSTER} -v s3://${OUTSTORE}/HGSVC/HGSVC.haps.vcf.gz -l s3://${OUTSTORE}/HGSVC/HGSVC_alts.gam ${JOBSTORE} ${OUTSTORE}/HGSVC/call-HG00514 s3://${OUTSTORE}/HGSVC/HGSVC.xg HG00514 s3://${OUTSTORE}/HGSVC/map-HG00514/HG00514_chr
./call.sh -c ${CLUSTER} -v s3://${OUTSTORE}/HGSVC/HGSVC.haps.vcf.gz -l s3://${OUTSTORE}/HGSVC/HGSVC_alts.gam ${JOBSTORE} ${OUTSTORE}/HGSVC/call-HG00733 s3://${OUTSTORE}/HGSVC/HGSVC.xg HG00733 s3://${OUTSTORE}/HGSVC/map-HG00733/HG00733_chr
./call.sh -c ${CLUSTER} -v s3://${OUTSTORE}/HGSVC/HGSVC.haps.vcf.gz -l s3://${OUTSTORE}/HGSVC/HGSVC_alts.gam ${JOBSTORE} ${OUTSTORE}/HGSVC/call-NA19240 s3://${OUTSTORE}/HGSVC/HGSVC.xg NA19240 s3://${OUTSTORE}/HGSVC/map-NA19240/NA19240_chr

# vg call no longer needs toil-vg
# for the timing benchmarks in the paper, it was run directly on the whole genome output (GAM created by catting chromosome gams together)
# vg pack -x HGSVC.xg -Q 5 -t 10 -g HG00514.gam -o HGSVC.HG00514.pack
# vg call HGSVC.xg -k HGSVC.HG00514.pack -v HGSVC.haps.vcf.gz -t 10 -r HGSVC.snarls | bgzip > HG00514.vg.vcf.gz
```

## Mapping reads to GRCh38

Other methods use reads mapped to the linear reference genome, either for genotyping (Delly, SVTyper) or to identify SNVs/indels (BayesTyper).

Reads were mapped with bwa. 
For example for NA19240:

```
bwa mem -t 12 hg38_noalt.fa ERR894724_1.fq.gz ERR894724_2.fq.gz | samtools view -bh - | samtools sort -o ERR894724_bwa_mem_sort.bam -
samtools index ERR894724_bwa_mem_sort.bam
picard AddOrReplaceReadGroups I=ERR894724_bwa_mem_sort.bam O=ERR894724_bwa_mem_sort_rg.bam LB=ERR894724 PL=ILLUMINA PU=ILLUMINA SM=ERR894724
picard MarkDuplicates I=ERR894724_bwa_mem_sort_rg.bam O=ERR894724_bwa_mem_sort_rg_md.bam M=md_metrics.txt CREATE_INDEX=true TMP_DIR=.
rm ERR894724_bwa_mem_sort_rg.*
```

## Paragraph

Each sample was genotyped individually.
For example for HG00514:

```
echo -e "id\tpath\tdepth\tread length\nHG00514\tHG00514.bam\t30\t150" > samples_for_paragraph_HG00514.txt
multigrmpy.py -m samples_for_paragraph_HG00514.txt -i HGSVC.haps.vcf.gz -r hg38.fa -t 10 -o paragraph_out_HG00514
```

The output genotypes are in `paragraph_out_HG00514/genotypes.vcf.gz`.

## Delly

Each sample was genotyped individually.
For example for NA19240:

```
delly call -g hg38.fa -v HGSVC.haps.vcf.gz -o temp.bcf ERR894724_bwa_mem_sort_rg_md.bam
bcftools convert temp.bcf > NA19240.delly.vcf
bgzip NA19240.delly.vcf
tabix NA19240.delly.vcf.gz
```

## SVTyper

The explicit VCF was first converted into a symbolic VCF, for deletions only because SVTyper cannot genotype insertions.
See *VCFtoSymbolicDEL.py* in [../misc-scripts](../misc-scripts).

```
zcat HGSVC.haps.vcf.gz | python VCFtoSymbolicDEL.py > HGSVC.haps.symbolic.dels.vcf
```

Each sample was then genotyped individually.
The output VCF was converted back to explicit format using a helper scripts from the BayesTyper tools.
For example for NA19240:

```
svtyper -i HGSVC.haps.symbolic.dels.vcf -B ERR894724_bwa_mem_sort_rg_md.bam -l libinfo.json > NA19240.svtyper.vcf
vcf-sort NA19240.svtyper.vcf | vcfkeepinfo - END SVLEN SVTYPE NS AN | grep -v 'contig' > NA19240.svtyper.clean.vcf
bayesTyperTools convertAllele -v NA19240.svtyper.clean.vcf -g ../../haps/hg38.fa --keep-imprecise -o NA19240.svtyper.explicit
vcfkeepinfo NA19240.svtyper.explicit.vcf NS AN | bgzip > NA19240.svtyper.explicit.vcf.gz
tabix -f NA19240.svtyper.explicit.vcf.gz
```

## BayesTyper

### kmer counting

For each sample (here for NA19240):

```
kmc -v -k55 -m128 -sm -ci1 -t12 @reads.txt ERR894724_k55 .
bayesTyperTools makeBloom -p 12 -k ERR894724_k55
```

With `reads.txt` containing two lines with the paths to the FASTQ files for the sample, for example:

```
ERR894724_1.fq.gz
ERR894724_2.fq.gz
```

### SNVs

SNVs were called using GATK Haplotype Caller and Platypus.
Variants were normalized and then filtered to keep only variants smaller than 50 bp.

#### GATK Haplotype Caller

For each sample (here for NA19240):

```
gatk3 -Xmx256G -T HaplotypeCaller -R hg38_noalt.fa -nct 12 -I ERR894724_bwa_mem_sort_rg_md.bam -o ERR894724_bwa_mem_sort_md_haplotypecaller.vcf
gzip ERR894724_bwa_mem_sort_md_haplotypecaller.vcf
VCF_PREFIX=ERR894724_bwa_mem_sort_md_haplotypecaller
bcftools norm -c x -m -any -f hg38_noalt.fa ${VCF_PREFIX}.vcf.gz | cut -f1-8 | bcftools annotate -x FORMAT,INFO > ${VCF_PREFIX}_norm.vcf
gzip ${VCF_PREFIX}_norm.vcf
filterStructuralVariants ERR894724_bwa_mem_sort_md_haplotypecaller_norm.vcf.gz ERR894724_bwa_mem_sort_md_haplotypecaller_norm_sv49 49 -49 49
```

#### Platypus

For each sample (here for NA19240):

```
platypus callVariants --bamFiles=ERR894724_bwa_mem_sort.bam --refFile=hg38_noalt.fa --output=ERR894724_bwa_mem_sort_platypus.vcf --logFileName=platypus_log.txt --nCPU 12 --assemble=1 --assembleBrokenPairs=1
bgzip ERR894724_bwa_mem_sort_platypus.vcf
tabix ERR894724_bwa_mem_sort_platypus.vcf.gz
VCF_PREFIX=ERR894724_bwa_mem_sort_platypus
bcftools norm -c x -m -any -f hg38_noalt.fa ${VCF_PREFIX}.vcf.gz | cut -f1-8 | bcftools annotate -x FORMAT,INFO > ${VCF_PREFIX}_norm.vcf
picard UpdateVcfSequenceDictionary I=${VCF_PREFIX}_norm.vcf O=${VCF_PREFIX}_norm_header.vcf SD=hg38_noalt.dict
picard SortVcf I=${VCF_PREFIX}_norm_header.vcf O=${VCF_PREFIX}_norm_sort.vcf SD=hg38_noalt.dict
rm ${VCF_PREFIX}_norm.vcf
rm ${VCF_PREFIX}_norm_header.vcf
rm ${VCF_PREFIX}_norm_sort.vcf.idx
gzip ${VCF_PREFIX}_norm_sort.vcf
filterStructuralVariants ERR894724_bwa_mem_sort_platypus_norm_sort.vcf.gz ERR894724_bwa_mem_sort_platypus_norm_sort_sv49 49 -49 49
```

### Combining variants

SNVs/indels across all samples were combined with the SV catalog.

```
bayesTyperTools combine -v hgsvc:HGSVC.haps_norm_sort.vcf.gz,pp:ERR903030_bwa_mem_sort_platypus_norm_sort_sv49.vcf.gz,pp:ERR895347_bwa_mem_sort_platypus_norm_sort_sv49.vcf.gz,pp:ERR894724_bwa_mem_sort_platypus_norm_sort_sv49.vcf.gz,hc:ERR903030_bwa_mem_sort_md_haplotypecaller_norm_sv49.vcf.gz,hc:ERR895347_bwa_mem_sort_md_haplotypecaller_norm_sv49.vcf.gz,hc:ERR894724_bwa_mem_sort_md_haplotypecaller_norm_sv49.vcf.gz -o all_hgsvc_pp49_hc49 -z
```

### Genotyping

The three samples were genotyped jointly.

```
OUT_PREFIX=bayestyper_unit_1/all_hgsvc_pp49_hc49_bayestyper

bayesTyper cluster -v all_hgsvc_pp49_hc49.vcf.gz -s samples.txt -g hg38_canon.fa -d hg38_decoy.fa -p 24 -r 12345678 --min-number-of-unit-variants 10000000
bayesTyper genotype -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_cluster_data -s samples.txt -g hg38_canon.fa -d hg38_decoy.fa -o ${OUT_PREFIX} -z -p 24 -r 12345678 --min-genotype-posterior 0

bcftools filter -i 'FILTER=="PASS"' ${OUT_PREFIX}.vcf.gz | gzip -c 1> ${OUT_PREFIX}_pass.vcf.gz
filterAlleleCallsetOrigin ${OUT_PREFIX}_pass.vcf.gz ${OUT_PREFIX}_pass_hgsvc hc,pp,. 1
```

`samples.txt` contains paths to the kmers counts:

```
HG00514	F	ERR903030_k55
HG00733	F	ERR895347_k55
NA19240	F	ERR894724_k55
```
