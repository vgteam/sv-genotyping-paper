# Pseudo-diploid CHM genome

## Data

- [GRCh38 reference genome](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).

VCF `pseudo_diploid-explicit.vcf.gz` was created using:

```
./download.sh
./add-genotypes.py pseudo_diploid.vcf.gz reduced.tab | bgzip > pseudo_diploid_gt.vcf.gz
# Note: HG38 should point to a reference fasta that has a .fai index
./make-explicit.py pseudo_diploid_gt.vcf.gz --fasta ${HG38} | vcfkeepinfo - NA | vcffixup - | bgzip > pseudo_diploid-explicit.vcf.gz
tabix -f -p vcf pseudo_diploid-explicit.vcf.gz
```

Reads?


## toil-vg

Using the helper scripts from `../toil-scripts`.

```
# make the CHM-PSEUDODIPLOID graph
./construct.sh -s -c ${CLUSTER}2 ${JOBSTORE}2 ${OUTSTORE}/CHMPD-feb12

# map the 30x reads and call variants
./mapcall.sh -c ${CLUSTER}3 ${JOBSTORE}3 ${OUTSTORE}/CHMPD-feb12 s3://${OUTSTORE}/CHMPD-feb12/CHMPD PSEUDOSET PSEUDOSET-30 s3://majorsv-ucsc/gt/gam/aln_30x.gam

# download VCF file
rm -rf ./chmpd-vg-chmpd.vcf.gz ; aws s3 sync s3://${OUTSTORE}/CHMPD-feb12/CHMPD/call-PSEUDOSET-30/PSEUDOSET.vcf.gz ./chmpd-vg-chmpd.vcf.gz
```

## SMRTSV2

```
aws s3 sync s3://${OUTSTORE}/CHMPD-feb12/eval-PSEUDOSET-30 ./CHMPD-feb12-eval-PSEUDOSET
# then do the SMRTSV that we copied from courtyard
# note that we made it explicit with
# ./make-explicit.py PSEUDOSET-smrtsv.vcf.gz --fasta ~/dev/work/references/hg38.fa.gz  | vcfkeepinfo - NA | vcffixup - | bgzip > PSEUDOSET-smrtsv-explicit.vcf.gz
```
