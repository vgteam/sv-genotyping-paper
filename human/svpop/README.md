# SVPOP

## Data

- [GRCh38 reference genome](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).

The VCF file with the SV catalog (`sv-pop-explicit.vcf.gz`) was prepared with:

```
wget -nc http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181025_EEE_SV-Pop_1/VariantCalls_EEE_SV-Pop_1/EEE_SV-Pop_1.ALL.sites.20181204.bed.gz
wget -nc http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181025_EEE_SV-Pop_1/VariantCalls_EEE_SV-Pop_1/EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz
wget -nc http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181025_EEE_SV-Pop_1/VariantCalls_EEE_SV-Pop_1/EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz.tbi

# fixes up some basic info tags like AC and AF
vcfkeepinfo EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz SVTYPE | vcffixup - | bgzip > EEE_SV-Pop_1.ALL.sites.20181204.fix.vcf.gz
./vcf-add-bed-seqs.py --inv leave EEE_SV-Pop_1.ALL.sites.20181204.fix.vcf.gz EEE_SV-Pop_1.ALL.sites.20181204.bed.gz | bgzip > sv-pop-explicit.vcf.gz
tabix -f -p vcf sv-pop-explicit.vcf.gz
```

### Reads

- HG00514: [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR903/ERR903030/ERR903030_1.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR903/ERR903030/ERR903030_1.fastq.gz) [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR903/ERR903030/ERR903030_2.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR903/ERR903030/ERR903030_2.fastq.gz)
- HG00733: [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_1.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_1.fastq.gz) [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_2.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_2.fastq.gz)
- NA19240: [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_1.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_1.fastq.gz) [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_2.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_2.fastq.gz)


## toil-vg

Using the helper scripts from `../toil-scripts`.

```
# Construct the graph and index for svpop (15 sv samples) including inversions
./construct.sh -p -i -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}/SVPOP

# Mapping and calling

./map.sh -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}/SVPOP/map-HG00514 s3://${OUTSTORE}/SVPOP/SVPOP HG00514 ${FQBASE}/HG00514/ERR903030_1.fastq.gz ${FQBASE}/HG00514/ERR903030_2.fastq.gz 
./map.sh -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}/SVPOP/map-HG00733 s3://${OUTSTORE}/SVPOP/SVPOP HG00733 ${FQBASE}/HG00733/ERR895347_1.fastq.gz ${FQBASE}/HG00733/ERR895347_2.fastq.gz
./map.sh -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}/SVPOP/map-NA19240 s3://${OUTSTORE}/SVPOP/SVPOP NA19240 ${FQBASE}/NA19240/ERR894724_1.fastq.gz ${FQBASE}/NA19240/ERR894724_2.fastq.gz

./call.sh -c ${CLUSTER} -v s3://${OUTSTORE}/SVPOP/sv-pop-explicit.vcf.gz -l s3://${OUTSTORE}/SVPOP/SVPOP_alts.gam ${JOBSTORE} ${OUTSTORE}/SVPOP/call-HG00514 s3://${OUTSTORE}/SVPOP/SVPOP.xg HG00514 s3://${OUTSTORE}/SVPOP/map-HG00514/HG00514_chr
./call.sh -c ${CLUSTER} -v s3://${OUTSTORE}/SVPOP/sv-pop-explicit.vcf.gz -l s3://${OUTSTORE}/SVPOP/SVPOP_alts.gam ${JOBSTORE} ${OUTSTORE}/SVPOP/call-HG00733 s3://${OUTSTORE}/SVPOP/SVPOP.xg HG00733 s3://${OUTSTORE}/SVPOP/map-HG00733/HG00733_chr
./call.sh -c ${CLUSTER} -v s3://${OUTSTORE}/SVPOP/sv-pop-explicit.vcf.gz -l s3://${OUTSTORE}/SVPOP/SVPOP_alts.gam ${JOBSTORE} ${OUTSTORE}/SVPOP/call-NA19240 s3://${OUTSTORE}/SVPOP/SVPOP.xg NA19240 s3://${OUTSTORE}/SVPOP/map-NA19240/NA19240_chr

# vg call no longer needs toil-vg
# for the timing benchmarks in the paper, it was run directly on the whole genome output (GAM created by catting chromosome gams together)
# vg pack -x SVPOP.xg -Q 5 -t 10 -g HG00514.gam -o SVPOP.HG00514.pack
# vg call SVPOP.xg -k SVPOP.HG00514.pack -v SVPOP.haps.vcf.gz -r SVPOP.snarls -t 10 | bgzip > HG00514.vg.vcf.gz
```

## SMRT-SV v2

SMRT-SV was obtained from [github](https://github.com/EichlerLab/smrtsv2) using commit 43d28cbd4e76ed2a4e820ded11a80592675db6c9 (the then-current master branch).
It was then run as described in its [README](https://github.com/EichlerLab/smrtsv2/blob/master/GENOTYPE.md), on a single server
```
${SMRTSV_DIR}/smrtsv --jobs 30 genotype genotyper.json variants.vcf.gz

```
Where genotyper.json contained (the VCF and BAM were downloaded from the FTP site listed above (ftp.1000genomes.ebi.ac.uk))
```
{
  "sample_manifest": "svpop-samples.tab",
  "sv_reference": "hg38.fa",
  "sv_calls": "EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz",
  "sv_contigs": "EEE_SV-Pop_1.ALL.sites.20181204.bam",
  "model": "30x-4",
  "min_call_depth": 8
}

```
and hsvc-samples.tab contained:
```
SAMPLE	SEX	DATA
HG00514	F	hgsvc-bams/ERR903030_bwa_mem_sort_rg_md.bam
HG00733	F	hgsvc-bams/ERR895347_bwa_mem_sort_rg_md.bam
NA19240	F	hgsvc-bams/ERR894724_bwa_mem_sort_rg_md.bam
```

The BAMs were made as described [here](https://github.com/vgteam/sv-genotyping-paper/tree/master/human/hgsvc)

The VCFs output by this version of SMRT-SV v2 are invalid due to a bug which outputs the wrong values in the REF column.  This was corrected by re-extracting the REF values from the FASTA using [this script](https://github.com/vgteam/sv-genotyping-paper/blob/master/human/chmpd/make-explicit.py) on variants.vcf.gz.  At the same time, the symbolic alleles were converted to explict alleles.  

```
./make-explicit.py variants.vcf.gz --fasta hg38.fa.gz  | vcfkeepinfo - NA | vcffixup - | bgzip > svpop-smrtsv-explicit.vcf.gz
```
