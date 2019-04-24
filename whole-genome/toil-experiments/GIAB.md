# this is inspired by the analysis from https://www.biorxiv.org/content/biorxiv/early/2019/02/24/558247.full.pdf
# they used
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_UnionSVs_12122017/

# a newer one with bed files is here
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/

# construct a 0.60 graph with HG002 control
./construct-hgsvc.sh  -H -c ${CLUSTER}4 ${JOBSTORE}4 ${OUTSTORE}/GIAB-FEB26

# construct a 0.50 graph with HG002 control
./construct-hgsvc.sh  -G -c ${CLUSTER}1 ${JOBSTORE}1x ${OUTSTORE}/GIAB-0.5-FEB26

# map our downsampled 50X reads from the paper

./mce-hgsvc.sh -p -E " -p -n -m 20" -C " -p" -c ${CLUSTER}4  ${JOBSTORE}4 ${OUTSTORE}/GIAB-FEB26 s3://${OUTSTORE}/GIAB-FEB26/GIAB HG002 HG002 s3://${OUTSTORE}/GIAB-FEB26/giab-0.6.vcf.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed s3://glennhickey/mike-hs37d5-hg002-nov18/bams/HG002-NA24385-50x.bam

./mce-hgsvc.sh -p -E " -p -n -m 20" -C " -p" -c ${CLUSTER}1  ${JOBSTORE}1 ${OUTSTORE}/GIAB-0.5-FEB26 s3://${OUTSTORE}/GIAB-0.5-FEB26/GIAB HG002 HG002 s3://${OUTSTORE}/GIAB-0.5-FEB26/giab-0.5.vcf.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed s3://glennhickey/mike-hs37d5-hg002-nov18/bams/HG002-NA24385-50x.bam

aws s3 sync s3://${OUTSTORE}/GIAB-FEB26/eval-HG002 ./GIAB-FEB26-eval-HG002
aws s3 sync s3://${OUTSTORE}/GIAB-0.5-FEB26/eval-HG002 ./GIAB-0.5-FEB26-eval-HG002

