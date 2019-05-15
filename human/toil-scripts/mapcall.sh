# Map reads to a graph and calls variants

#!/bin/bash

BID=0.53
RESUME=0
REGION="us-west-2"
HEAD_NODE_OPTS=""
MAP_OPTS=""
CALL_OPTS=""
CHR_PREFIX="chr"

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [OPTIONS] <JOBSTORE-NAME> <OUTSTORE-NAME> <INDEX-BASE> <SAMPLE> <NAME> <FASTQ1> [FASTQ2]\n"
	 printf "Arguments:\n"
	 printf "   JOBSTORE-NAME: Name of Toil S3 Jobstore (without any prefix). EX: my-job-store \n"
	 printf "   OUTSTORE-NAME: Name of output bucket (without prefix or trailing slash). EX my-bucket/hgsvc\n"
	 printf "   INDEX-BASE:    Full path of indexes minus the file extension\n"
	 printf "   SAMPLE:        Name of sample\n"
	 printf "   NAME:          Name of output \n"
	 printf "   FASTQ1:        Path of fastq reads (assume interleaved if FASTQ2 not given)\n"
	 printf "   FASTQ2:        Path of fastq reads (optional)\n"
	 printf "Options:\n"
	 printf "   -b BID  Spot bid in dollars for r3.8xlarge nodes [${BID}]\n"
	 printf "   -c      Toil Cluster Name (created with https://github.com/vgteam/toil-vg/blob/master/scripts/create-ec2-leader.sh).  Only use if not running from head node.\n"
	 printf "   -M      options for map.sh in quotes (SKIP to skip)\n"
	 printf "   -C      options for call.sh in quotes (SKIP to skip)\n"
	 printf "   -p      No chr prefix in chromosome names\n"
    exit 1
}

while getopts "b:c:M:C:E:p" o; do
    case "${o}" in
        b)
            HEAD_NODE_OPTS="${HEAD_NODE_OPTS} -b ${OPTARG}"
            ;;
		  c)
				HEAD_NODE_OPTS="${HEAD_NODE_OPTS} -c ${OPTARG}"
				;;
		  M)
				MAP_OPTS="${OPTARG}"
				;;
		  C)
				CALL_OPTS="${OPTARG}"
				;;
		  p)
				CHR_PREFIX=""
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
INDEX_BASE="${1}"
shift
SAMPLE="${1}"
shift
NAME="${1}"
shift
READS1="${1}"
shift
READS2="${1}"
shift

set -ex

if [ "${MAP_OPTS}" != "SKIP" ]
then
	 ./map.sh $HEAD_NODE_OPTS $MAP_OPTS $JOBSTORE_NAME "${OUTSTORE_NAME}/map-${NAME}" $INDEX_BASE $SAMPLE $READS1 $READS2
fi

if [ "${CALL_OPTS}" != "SKIP" ]
then
	 ./call.sh $HEAD_NODE_OPTS $CALL_OPTS $JOBSTORE_NAME "${OUTSTORE_NAME}/call-${NAME}" ${INDEX_BASE}.xg $SAMPLE "s3://${OUTSTORE_NAME}/map-${NAME}/${SAMPLE}_${CHR_PREFIX}"
fi
