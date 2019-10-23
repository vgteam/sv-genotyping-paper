# Pseudo-diploid CHM genome

## Data

- [GRCh38 reference genome](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).

This dataset was generated and used in [Audano et al. Cell 2019](https://doi.org/10.1016/j.cell.2018.12.019).
It is not, to our knowledge, publicly available. 
We obtained the following files (referred to below) directly from SMRT-SV v2's author, Peter Audano.

```
# VCF to genotype
variants_reduced.vcf.gz
# Information needed for evaluaiont
reduced.tab
# Corresponding contig alignment
contigs_reduced.bam
# 30X short reads
pseudo_30x.bam
```

VCF `pseudo_diploid-explicit.vcf.gz` was created using:

```
# remove the END tag, which was interfering with a later script
vt rminfo variants_reduced.vcf.gz -t END -o pseudo_diploid.vcf
bgzip -f pseudo_diploid.vcf
tabix -f -p vcf pseudo_diploid.vcf.gz
# extract the genotype information from the table and put it into the VCF
./add-genotypes.py pseudo_diploid.vcf.gz reduced.tab | bgzip > pseudo_diploid_gt.vcf.gz
# fix bugs in REF columns and convert to explicit variants which are more robustly-supported by vg
# (and to be more consistent through evaluation down the road)
# Note: HG38 should point to a reference fasta that has a .fai index
./make-explicit.py pseudo_diploid_gt.vcf.gz --fasta ${HG38} | vcfkeepinfo - NA | vcffixup - | bgzip > pseudo_diploid-explicit.vcf.gz
tabix -f -p vcf pseudo_diploid-explicit.vcf.gz
```
## toil-vg

Using the helper scripts from `../toil-scripts`.

```
# make the CHM-PSEUDODIPLOID graph
./construct.sh -s -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}/CHMPD

# map the 30x reads and call variants
./map.sh -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}/CHMPD/map-PSEUDOSET s3://${OUTSTORE}/CHMPD PSEUDOSET ${FQBASE}/aln_30x.gam

./call.sh -c ${CLUSTER} -v s3://${OUTSTORE}/CHMPD/pseudo_diploid.vcf.gz  -l s3://${OUTSTORE}/CHMPD/CHMPD_alts.gam ${JOBSTORE} ${OUTSTORE}/CHMPD/call-PSEUDOSET s3://${OUTSTORE}/CHMPD/CHMPD.xg PSEUDOSET s3://${OUTSTORE}/CHMPD/map-PSEUDOSET/PSEUDOSET_chr
```

## SMRT-SV v2

SMRT-SV was obtained from [github](https://github.com/EichlerLab/smrtsv2) using commit 43d28cbd4e76ed2a4e820ded11a80592675db6c9 (the then-current master branch).
It was then run as described in its [README](https://github.com/EichlerLab/smrtsv2/blob/master/GENOTYPE.md), on a single server
```
${SMRTSV_DIR}/smrtsv --jobs 20 genotype genotyper.json variants.vcf.gz

```
Where genotyper.json contained (the VCF and BAM were downloaded from the FTP site listed above (ftp.1000genomes.ebi.ac.uk))
```
{
  "sample_manifest": "chmpd-samples.tab",
  "sv_reference": "hg38.fa",
  "sv_calls": "variants_reduced.vcf.gz",
  "sv_contigs": "contigs_reduced.bam",
  "model": "30x-4",
  "min_call_depth": 8
}
```

and chmpd-samples.tab contained:
```
SAMPLE	SEX	DATA
PSEUDOSET	U	pseudo_30x.bam
```

The VCFs output by this version of SMRT-SV v2 are invalid due to a bug which outputs the wrong values in the REF column.  This was corrected by re-extracting the REF values from the FASTA using [this script](https://github.com/vgteam/sv-genotyping-paper/blob/master/human/chmpd/make-explicit.py) on variants.vcf.gz.  At the same time, the symbolic alleles were converted to explict alleles.  

```
./make-explicit.py variants.vcf.gz --fasta hg38.fa.gz  | vcfkeepinfo - NA | vcffixup - | bgzip > PSEUDOSET-smrtsv-explicit.vcf.gz
```
