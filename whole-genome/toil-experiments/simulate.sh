# Simulate 30X coverage for the Construct a graph for hg38 (chromosomes only) from the HGSVC vcf
# EX ./simulate.sh -c my-cluster my-jobstore my-bucket/hgsvc/sim/ s3://my-bucket/hgsvc/HGSVC.chroms_HG00514_haplo_thread_0.xg  s3://my-bucket/hgsvc/HGSVC.chroms_HG00514_haplo_thread_1.xg s3://my-bucket/interleaved.fq

#!/bin/bash

BID=0.78
RESUME=0
REGION="us-west-2"
HEAD_NODE_OPTS=""
PATH_OPTS=""
OUT_NAME="sim-HG00514-30x"
NUM_READS="256000000"

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [OPTIONS] <JOBSTORE-NAME> <OUTSTORE-NAME> <HAPLO-XG0> <HAPLO-XG1> <TEMPLATE-FQ>\n"
	 printf "Arguments:\n"
	 printf "   JOBSTORE-NAME: Name of Toil S3 Jobstore (without any prefix). EX: my-job-store \n"
	 printf "   OUTSTORE-NAME: Name of output bucket (without prefix or trailing slash). EX my-bucket/hgsvc\n"
	 printf "   HAPLO-XG0:     Full path to xg for haplotype 0\n"
	 printf "   HAPLO-XG2:     Full path to xg for haplotype 1\n"
	 printf "   TEMPLATE-FQ:   Full path to interleaved fastq file to use as template\n"
	 printf "Options:\n"
	 printf "   -b BID  Spot bid in dollars for i3.8xlarge nodes [${BID}]\n"
	 printf "   -r      Resume existing job\n"
	 printf "   -g      Aws region [${REGION}]\n"
	 printf "   -c      Toil Cluster Name (created with https://github.com/vgteam/toil-vg/blob/master/scripts/create-ec2-leader.sh).  Only use if not running from head node.\n"
	 printf "   -p PATH Simulate from this path\n"
	 printf "   -n N    Simulate N reads\n"
    exit 1
}

while getopts "b:re:c:p:n:" o; do
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
		  p)
				PATH_OPTS="--path ${OPTARG}"
				OUT_NAME="${OUT_NAME}-${OPTARG}"
				;;
		  n)
				NUM_READS=${OPTARG}
				;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [[ "$#" -lt "5" ]]; then
    # Too few arguments
    usage
fi

# of the form aws:us-west:name
JOBSTORE_NAME="${1}"
shift
OUTSTORE_NAME="${1}"
shift
XG0_PATH="${1}"
shift
XG1_PATH="${1}"
shift
TEMPLATE_PATH="${1}"
shift

# pull in ec2-run from git if not found in current dir
wget -nc https://raw.githubusercontent.com/vgteam/toil-vg/master/scripts/ec2-run.sh
chmod 777 ec2-run.sh

# without -r we start from scratch!
RESTART_FLAG=""
if [ $RESUME == 0 ]
then
	 toil clean aws:${REGION}:${JOBSTORE_NAME}
else
	 RESTART_FLAG="--restart"
fi

CMD="sim aws:${REGION}:${JOBSTORE_NAME} ${XG0_PATH} ${XG1_PATH} ${NUM_READS} aws:${REGION}:${OUTSTORE_NAME} --out_name ${OUT_NAME} --gam --fastq_out --fastq  ${TEMPLATE_PATH} --sim_opts \"-p 570 -v 165 -i 0.002 -I\" --sim_chunks 20 --seed 23 --validate --whole_genome_config --logFile simulate.hgsvc.log ${RESTART_FLAG} ${PATH_OPTS}"

# run the job
./ec2-run.sh ${HEAD_NODE_OPTS} -m 50 -n i3.8xlarge:${BID},i3.8xlarge "${CMD}" | tee sim.hgsvc.$(basename ${OUTSTORE_NAME}).stdout

aws s3 cp sim.hgsvc.$(basename ${OUTSTORE_NAME}).stdout s3://${OUTSTORE_NAME}/
