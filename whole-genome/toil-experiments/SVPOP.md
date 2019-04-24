# SVPOP

```
# Construct the graph and index for svpop (15 sv samples) including inversions
./construct-hgsvc.sh -p -i -c ${CLUSTER}2 ${JOBSTORE}2 ${OUTSTORE}/SVPOP-jan10

# Mapping, Calling, Evaluation

./mce-hgsvc.sh -c ${CLUSTER}1 ${JOBSTORE}1 ${OUTSTORE}/SVPOP-jan10 s3://${OUTSTORE}/SVPOP-jan10/SVPOP HG00514 HG00514 s3://${OUTSTORE}/SVPOP-jan10/sv-pop-explicit.vcf.gz ${COMPARE_REGIONS_BED} ${FQBASE}/HG00514/ERR903030_1.fastq.gz ${FQBASE}/HG00514/ERR903030_2.fastq.gz

./mce-hgsvc.sh -c ${CLUSTER}2 ${JOBSTORE}2 ${OUTSTORE}/SVPOP-jan10 s3://${OUTSTORE}/SVPOP-jan10/SVPOP HG00733 HG00733 s3://${OUTSTORE}/SVPOP-jan10/sv-pop-explicit.vcf.gz ${COMPARE_REGIONS_BED} ${FQBASE}/HG00733/ERR895347_1.fastq.gz ${FQBASE}/HG00733/ERR895347_2.fastq.gz 

./mce-hgsvc.sh -c ${CLUSTER}3 ${JOBSTORE}3 ${OUTSTORE}/SVPOP-jan10 s3://${OUTSTORE}/SVPOP-jan10/SVPOP NA19240 NA19240 s3://${OUTSTORE}/SVPOP-jan10/sv-pop-explicit.vcf.gz ${COMPARE_REGIONS_BED} ${FQBASE}/NA19240/ERR894724_1.fastq.gz ${FQBASE}/NA19240/ERR894724_2.fastq.gz

# Download

rm -rf ./SVPOP-jan10-eval-HG00514 ; aws s3 sync s3://${OUTSTORE}/SVPOP-jan10/eval-HG00514 ./SVPOP-jan10-eval-HG00514
rm -rf ./SVPOP-jan10-eval-HG00733 ; aws s3 sync s3://${OUTSTORE}/SVPOP-jan10/eval-HG00733 ./SVPOP-jan10-eval-HG00733
rm -rf ./SVPOP-jan10-eval-NA19240 ; aws s3 sync s3://${OUTSTORE}/SVPOP-jan10/eval-NA19240 ./SVPOP-jan10-eval-NA19240

#### SMRTSV

# Smrtsv2 was run on courtyard, then made explicity using the same process as the original svpop graph (see construg-hgsvc.sh)
# (all three samples are in one VCF: svpop-smrtsv-hgsvc-explicit.vcf.gz)

./eval-hgsvc.sh -c ${CLUSTER}1 ${JOBSTORE}1 ${OUTSTORE}/SVPOP-jan10/eval-HG00514-smrtsv s3://${OUTSTORE}/SVPOP-jan10/sv-pop-explicit.vcf.gz  s3://glennhickey/outstore/SVPOP-jan10/call-all-smrtsv/svpop-smrtsv-hgsvc-explicit.vcf.gz ${COMPARE_REGIONS_BED} HG00514

./eval-hgsvc.sh -c ${CLUSTER}1 ${JOBSTORE}1 ${OUTSTORE}/SVPOP-jan10/eval-HG00733-smrtsv s3://${OUTSTORE}/SVPOP-jan10/sv-pop-explicit.vcf.gz  s3://glennhickey/outstore/SVPOP-jan10/call-all-smrtsv/svpop-smrtsv-hgsvc-explicit.vcf.gz ${COMPARE_REGIONS_BED} HG00733

./eval-hgsvc.sh -c ${CLUSTER}1 ${JOBSTORE}1 ${OUTSTORE}/SVPOP-jan10/eval-NA19240-smrtsv s3://${OUTSTORE}/SVPOP-jan10/sv-pop-explicit.vcf.gz  s3://glennhickey/outstore/SVPOP-jan10/call-all-smrtsv/svpop-smrtsv-hgsvc-explicit.vcf.gz ${COMPARE_REGIONS_BED} NA19240

rm -rf ./SVPOP-jan10-eval-HG00514-smrtsv ; aws s3 sync s3://${OUTSTORE}/SVPOP-jan10/eval-HG00514-smrtsv ./SVPOP-jan10-eval-HG00514-smrtsv
rm -rf ./SVPOP-jan10-eval-HG00733-smrtsv ; aws s3 sync s3://${OUTSTORE}/SVPOP-jan10/eval-HG00733-smrtsv ./SVPOP-jan10-eval-HG00733-smrtsv
rm -rf ./SVPOP-jan10-eval-NA19240-smrtsv ; aws s3 sync s3://${OUTSTORE}/SVPOP-jan10/eval-NA19240-smrtsv ./SVPOP-jan10-eval-NA19240-smrtsv
