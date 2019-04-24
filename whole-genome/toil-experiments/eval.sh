# Evaluate SV calls on a hgsvc graph

#!/bin/bash

BID=0.53
RESUME=0
REGION="us-west-2"
HEAD_NODE_OPTS=""
DOWNLOAD=0
MIN_LEN=50
HG38=1
NORM_TRUTH=1

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [OPTIONS] <JOBSTORE-NAME> <OUTSTORE-NAME> <TRUTH-VCF> <CALLS-VCF> <REGION-BED> <SAMPLE>\n"
	 printf "Arguments:\n"
	 printf "   JOBSTORE-NAME: Name of Toil S3 Jobstore (without any prefix). EX: my-job-store \n"
	 printf "   OUTSTORE-NAME: Name of output bucket (without prefix or trailing slash). EX my-bucket/hgsvc\n"
	 printf "   TRUTH-VCF:     Truth VCF\n"
	 printf "   CALLS-VCF:     Calls VCF\n"
	 printf "   REGIONS-BED:   Regions to clip to\n"
	 printf "   SAMPLE:        Sample name\n"
	 printf "Options:\n"
	 printf "   -b BID  Spot bid in dollars for r3.8xlarge nodes [${BID}]\n"
	 printf "   -r      Resume existing job\n"
	 printf "   -g      Aws region [${REGION}]\n"
	 printf "   -c      Toil Cluster Name (created with https://github.com/vgteam/toil-vg/blob/master/scripts/create-ec2-leader.sh).  Only use if not running from head node.\n"
	 printf "   -d      Download results locally\n"
	 printf "   -m      Minimum length [${MIN_LEN}]\n"
	 printf "   -p      hs37d5 reference\n"
	 printf "   -n      only normalize calls\n"
    exit 1
}

while getopts "b:re:c:dm:pn" o; do
    case "${o}" in
        b)
            BID=${OPTARG}
            ;;
        r)
            RESUME=1
            ;;
		  e)
            REGION=${OPTARG}
            ;;
		  c)
				HEAD_NODE_OPTS="-l ${OPTARG}"
				;;
		  d)
				DOWNLOAD=1
				;;
		  m)
				MIN_LEN=${OPTARG}
				;;
		  p)
				HG38=0
				;;
		  n)
				NORM_TRUTH=0
				;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [[ "$#" -lt "6" ]]; then
    # Too few arguments
    usage
fi

# of the form aws:us-west:name
JOBSTORE_NAME="${1}"
shift
OUTSTORE_NAME="${1}"
shift
TRUTH_VCF="${1}"
shift
CALLS_VCF="${1}"
shift
REGIONS_BED="${1}"
shift
SAMPLE_NAME="${1}"

# assume interleaved if READS2 not given
if [ -z ${READS2} ]
then
	 READS_OPTS="--fastq ${READS1} --interleaved"
else
	 READS_OPTS="--fastq ${READS1} ${READS2}"
fi

# pull in ec2-run from git if not found in current dir
wget -nc https://raw.githubusercontent.com/vgteam/toil-vg/master/scripts/ec2-run.sh
chmod 777 ec2-run.sh

# without -r we start from scratch! 
if [ $RESUME == 0 ]
then
	 toil clean aws:${REGION}:${JOBSTORE_NAME}
fi

if [ $HG38 == 1 ]
then
	 FASTA="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
else
	 FASTA="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz"
fi

SVEVAL_OPTS="--sveval --min_sv_len $MIN_LEN --vcfeval_sample $SAMPLE_NAME --check_inv"

if [ $NORM_TRUTH == 1 ]
then
	 NORM_OPTS="--normalize --vcfeval_fasta ${FASTA}"
else
	 NORM_OPTS="--normalize_calls --vcfeval_fasta ${FASTA}"
fi

set -e

# run the job
./ec2-run.sh ${HEAD_NODE_OPTS} -n r3.8xlarge:${BID},r3.8xlarge "vcfeval aws:${REGION}:${JOBSTORE_NAME} aws:${REGION}:${OUTSTORE_NAME}/sveval --whole_genome_config --vcfeval_baseline ${TRUTH_VCF} --call_vcf ${CALLS_VCF} ${SVEVAL_OPTS}"

./ec2-run.sh ${HEAD_NODE_OPTS} -n r3.8xlarge:${BID},r3.8xlarge "vcfeval aws:${REGION}:${JOBSTORE_NAME} aws:${REGION}:${OUTSTORE_NAME}/sveval-clip --whole_genome_config --vcfeval_baseline ${TRUTH_VCF} --call_vcf ${CALLS_VCF} ${SVEVAL_OPTS} --vcfeval_bed_regions ${REGIONS_BED}"

./ec2-run.sh ${HEAD_NODE_OPTS} -n r3.8xlarge:${BID},r3.8xlarge "vcfeval aws:${REGION}:${JOBSTORE_NAME} aws:${REGION}:${OUTSTORE_NAME}/sveval-norm --whole_genome_config --vcfeval_baseline ${TRUTH_VCF} --call_vcf ${CALLS_VCF} ${SVEVAL_OPTS} ${NORM_OPTS}"

./ec2-run.sh ${HEAD_NODE_OPTS} -n r3.8xlarge:${BID},r3.8xlarge "vcfeval aws:${REGION}:${JOBSTORE_NAME} aws:${REGION}:${OUTSTORE_NAME}/sveval-clip-norm --whole_genome_config --vcfeval_baseline ${TRUTH_VCF} --call_vcf ${CALLS_VCF} ${SVEVAL_OPTS} --vcfeval_bed_regions ${REGIONS_BED} ${NORM_OPTS}"


###  Same thing againt with genotype comparisons.

./ec2-run.sh ${HEAD_NODE_OPTS} -n r3.8xlarge:${BID},r3.8xlarge "vcfeval aws:${REGION}:${JOBSTORE_NAME} aws:${REGION}:${OUTSTORE_NAME}/sveval-gt --whole_genome_config --vcfeval_baseline ${TRUTH_VCF} --call_vcf ${CALLS_VCF} ${SVEVAL_OPTS} --genotype_eval"

./ec2-run.sh ${HEAD_NODE_OPTS} -n r3.8xlarge:${BID},r3.8xlarge "vcfeval aws:${REGION}:${JOBSTORE_NAME} aws:${REGION}:${OUTSTORE_NAME}/sveval-clip-gt --whole_genome_config --vcfeval_baseline ${TRUTH_VCF} --call_vcf ${CALLS_VCF} ${SVEVAL_OPTS} --vcfeval_bed_regions ${REGIONS_BED} --genotype_eval"

./ec2-run.sh ${HEAD_NODE_OPTS} -n r3.8xlarge:${BID},r3.8xlarge "vcfeval aws:${REGION}:${JOBSTORE_NAME} aws:${REGION}:${OUTSTORE_NAME}/sveval-norm-gt --whole_genome_config --vcfeval_baseline ${TRUTH_VCF} --call_vcf ${CALLS_VCF} ${SVEVAL_OPTS} ${NORM_OPTS} --genotype_eval"

./ec2-run.sh ${HEAD_NODE_OPTS} -n r3.8xlarge:${BID},r3.8xlarge "vcfeval aws:${REGION}:${JOBSTORE_NAME} aws:${REGION}:${OUTSTORE_NAME}/sveval-clip-norm-gt --whole_genome_config --vcfeval_baseline ${TRUTH_VCF} --call_vcf ${CALLS_VCF} ${SVEVAL_OPTS} --vcfeval_bed_regions ${REGIONS_BED} ${NORM_OPTS} --genotype_eval"

TOIL_ERROR=!$

#./ec2-run.sh ${HEAD_NODE_OPTS} -n r3.8xlarge:${BID},r3.8xlarge "vcfeval aws:${REGION}:${JOBSTORE_NAME} aws:${REGION}:${OUTSTORE_NAME}/vcfeval --whole_genome_config --vcfeval_baseline ${TRUTH_VCF} --call_vcf ${CALLS_VCF}  --vcfeval_fasta ${FASTA} --vcfeval_opts \" --squash-ploidy --Xmax-length 15000\""

#./ec2-run.sh ${HEAD_NODE_OPTS} -n r3.8xlarge:${BID},r3.8xlarge "vcfeval aws:${REGION}:${JOBSTORE_NAME} aws:${REGION}:${OUTSTORE_NAME}/vcfeval-clip --whole_genome_config --vcfeval_baseline ${TRUTH_VCF} --call_vcf ${CALLS_VCF}  --vcfeval_bed_regions ${REGIONS_BED} ${NORM_OPTS}   --vcfeval_fasta ${FASTA} --vcfeval_opts \" --squash-ploidy --Xmax-length 15000\""


if [ $DOWNLOAD == 1 ]
then
	 aws s3 sync s3://${OUTSTORE_NAME} ./`basename $OUTSTORE_NAME`
fi

exit $TOIL_ERROR
