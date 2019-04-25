# Construct a graph for hg38 (chromosomes only) from the HGSVC vcf
# EX ./construct.sh -c my-cluster my-jobstore my-bucket/hgsvc

#!/bin/bash

BID=0.83
RESUME=0
REGION="us-west-2"
HEAD_NODE_OPTS=""
HEAD_NODE=""
INCLUDE_1KG=0
GRAPH="HGSVC"
INVERSIONS=0
MIN_AF=0
FORCE_OUTSTORE=0
NO_UNFOLD=0
CONFIG_PATH=""
#ALT_REGIONS="--alt_regions https://raw.githubusercontent.com/glennhickey/toil-vg/construct/data/grch38-alt-positions.bed"
#ALT_REGIONS="--alt_regions s3://glennhickey/grch38/grch38-alt-positions-no-hla.bed"
ALT_REGIONS=""
NODECOY=1
HG38=1

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [OPTIONS] <JOBSTORE-NAME> <OUTSTORE-NAME>\n"
	 printf "Arguments:\n"
	 printf "   JOBSTORE-NAME: Name of Toil S3 Jobstore (without any prefix). EX: my-job-store \n"
	 printf "   OUTSTORE-NAME: Name of output bucket (without prefix or trailing slash). EX my-bucket/hgsvc\n"
	 printf "Options:\n"
	 printf "   -b BID     Spot bid in dollars for i3.8xlarge nodes [${BID}]\n"
	 printf "   -r         Resume existing job\n"
	 printf "   -g REGION  Aws region [${REGION}]\n"
	 printf "   -c NAME    Toil Cluster Name (created with https://github.com/vgteam/toil-vg/blob/master/scripts/create-ec2-leader.sh).  Only use if not running from head node.\n"
	 printf "   -k         include thousand genomes VCFs\n"
	 printf "   -p         use sv-pop instead of HGSVC vcf\n"
	 printf "   -s         use pseudo-diploid instead of HGSVC vcf\n"
	 printf "   -G         use giab version 0.5.0 12122017 instead of HGSVC vcf\n"
	 printf "   -H         use giab version 0.6.0 instead of HGSVC vcf\n"
	 printf "   -i         include inversions in vcf\n"
	 printf "   -a AF      min af threshold\n"
	 printf "   -f         run with --force_outstore to preserve all intermediate files\n"
	 printf "   -n         no unfolding for gcsa prunding\n"
	 printf "   -l BED     alt regions\n"
	 printf "   -o         (local) Path of config file\n"
	 printf "   -D         Include HS38d1 decoy sequences\n"
    exit 1
}

while getopts "b:re:c:kpsia:fnl:o:DGH" o; do
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
				HEAD_NODE=${OPTARG}
				HEAD_NODE_OPTS="-l ${OPTARG}"
				;;
		  k)
				INCLUDE_1KG=1
				;;
		  p)
				GRAPH="SVPOP"
				;;
		  s)
				GRAPH="CHMPD"
				;;
		  G)
				GRAPH="GIAB-0.5"
				HG38=0
				;;
		  H)
				GRAPH="GIAB-0.6"
				HG38=0
				;;
		  i)
				INVERSIONS=1
				;;
		  a)
				MIN_AF=${OPTARG}
				;;
		  f)
				FORCE_OUTSTORE=1
				;;
		  n)
				NO_UNFOLD=1
				;;
		  l)
				ALT_REGIONS="--alt_regions_bed ${OPTARG}"
				;;
		  o)
				CONFIG_PATH=${OPTARG}
				;;
		  D)
				NODECOY=0
				;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [[ "$#" -lt "2" ]]; then
    # Too few arguments
    usage
fi

# of the form aws:us-west:name
JOBSTORE_NAME="${1}"
shift
OUTSTORE_NAME="${1}"
shift

# pull in ec2-run from git if not found in current dir
wget -nc https://raw.githubusercontent.com/vgteam/toil-vg/master/scripts/ec2-run.sh
chmod 777 ec2-run.sh

# make our vcf
if [ $GRAPH == "HGSVC" ]
then
	 pushd ../hgsvc		  
	 if [ $INVERSIONS == 0 ]
	 then
		  ./make-vcf.sh
		  popd
		  VCF=../hgsvc/HGSVC.haps.vcf.gz
	 else
		  ./make-inv-vcf.sh
		  popd
		  VCF=../hgsvc/HGSVC.haps.inv.vcf.gz
	 fi
	 NAME=HGSVC	  
elif [ $GRAPH == "SVPOP" ]
then
	 pushd ../svpop
	 ./download.sh
	 if [ $INVERSIONS == 0 ]
	 then
		  INV_OPTS="--inv drop"
	 else
		  INV_OPTS="--inv leave"
	 fi
	 # fixes up some basic info tags like AC and AF
	 vcfkeepinfo EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz SVTYPE | vcffixup - | bgzip > EEE_SV-Pop_1.ALL.sites.20181204.fix.vcf.gz
	 ./vcf-add-bed-seqs.py ${INV_OPTS} EEE_SV-Pop_1.ALL.sites.20181204.fix.vcf.gz EEE_SV-Pop_1.ALL.sites.20181204.bed.gz | bgzip > sv-pop.vcf.gz
	 tabix -f -p vcf sv-pop.vcf.gz
	 popd
	 VCF=../svpop/sv-pop.vcf.gz
	 NAME=SVPOP
elif [ $GRAPH == "CHMPD" ]
then
	 pushd ../chmpd
	 ./download.sh
	 ./add-genotypes.py pseudo_diploid.vcf.gz reduced.tab | bgzip > pseudo_diploid_gt.vcf.gz
	 # Note: HG38 should point to a reference fasta that has a .fai index
	 ./make-explicit.py pseudo_diploid_gt.vcf.gz --fasta ${HG38} | vcfkeepinfo - NA | vcffixup - | bgzip > pseudo_diploid-explicit.vcf.gz
	 tabix -f -p vcf pseudo_diploid_gt.vcf.gz
	 tabix -f -p vcf pseudo_diploid-explicit.vcf.gz
	 popd
	 VCF=../chmpd/pseudo_diploid-explicit.vcf.gz
	 NAME=CHMPD
elif [ $GRAPH == "GIAB-0.5" ]
then
	 pushd ../giab
	 wget -nc ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_UnionSVs_12122017/svanalyzer_union_171212_v0.5.0_annotated.vcf.gz
	 wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_UnionSVs_12122017/svanalyzer_union_171212_v0.5.0_annotated.vcf.gz.tbi
	 vcfkeepinfo svanalyzer_union_171212_v0.5.0_annotated.vcf.gz NA | vcffixup - | bgzip > giab-0.5.vcf.gz
	 tabix -f -p vcf giab-0.5.vcf.gz
	 popd
	 VCF=../giab/giab-0.5.vcf.gz
	 NAME=GIAB
elif [ $GRAPH == "GIAB-0.6" ]
then
	 pushd ../giab
	 wget -nc ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
	 wget -nc ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
	 vcfkeepinfo HG002_SVs_Tier1_v0.6.vcf.gz NA | vcffixup - | bgzip > giab-0.6.vcf.gz
	 tabix -f -p vcf giab-0.6.vcf.gz
	 popd
	 VCF=../giab/giab-0.6.vcf.gz
	 NAME=GIAB	 
fi

# Get our vcf on S3 in our outstore
aws s3 mb s3://${OUTSTORE_NAME} --region ${REGION}
sleep 5
aws s3 cp ${VCF} s3://${OUTSTORE_NAME}/
aws s3 cp ${VCF}.tbi s3://${OUTSTORE_NAME}/
S3VCF="s3://${OUTSTORE_NAME}/$(basename $VCF)"

# without -r we start from scratch!
RESTART_FLAG=""
if [ $RESUME == 0 ]
then
	 toil clean aws:${REGION}:${JOBSTORE_NAME}
else
	 RESTART_FLAG="--restart"
fi

if [ -z ${CONFIG_PATH} ]
then
	 CONFIG_OPTS="--whole_genome_config"
else
	 # pass a local config to our job by way of the S3 outstore
	 CONF_NAME=`basename $CONFIG_PATH`
	 aws s3 cp $CONFIG_PATH s3://${OUTSTORE_NAME}/${CONF_NAME}
	 toil ssh-cluster --insecure --logOff --zone=us-west-2a ${CLUSTER_NAME} ${HEAD_NODE} /venv/bin/aws s3 cp s3://${OUTSTORE_NAME}/${CONF_NAME} .
	 CONFIG_OPTS="--config $CONF_NAME"
fi

if [ $NODECOY == 0 ]
then
	 REGION_REGEX="--regions_regex \"chr.*_random\" \"chrUn_[a-zA-Z0-9]*\" \"chr.*decoy\" \"chrEBV\""
else
	 REGION_REGEX="--regions_regex \"chr.*_random\" \"chrUn_[a-zA-Z0-9]*\""
fi

if [ $HG38 == 1 ]
then
	 REGIONS="--regions $(for i in $(seq 1 22; echo X; echo Y; echo M); do echo chr${i}; done) --fasta_regions ${ALT_REGIONS} ${REGION_REGEX} --add_chr_prefix --mask_ambiguous"
	 FASTA="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
else
	 REGIONS="--regions $(for i in $(seq 1 22; echo X; echo Y; echo MT); do echo ${i}; done) --fasta_regions --remove_chr_prefix --mask_ambiguous"
	 FASTA="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz"
fi	 

CONTROLS="--pangenome"

if [ $INCLUDE_1KG == 1 ]
then
	 # Pass in a mix of our HGSVC and 1KG vcfs
	 VCFS="$(for i in $(seq 1 22; echo X; echo Y); do echo ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr${i}_GRCh38.genotypes.20170504.vcf.gz,${S3VCF}; done)"
	 OUT_NAME="${NAME}_1KG"
else
	 # just the HGSVC SVs (we pass in the same vcf for each region so toil-vg doesn't use it for the decoy sequences.
	 # this shouldn't matter, but the the job-chain can get too big for toil when doing unfold-pruning)
	 VCFS="$(for i in $(seq 1 22; echo X; echo Y); do echo ${S3VCF}; done)"
	 OUT_NAME="${NAME}"
	 if [ $GRAPH == "HGSVC" ]
	 then
		  CONTROLS="--pos_control HG00514 --haplo_sample HG00514 --neg_control HG00514 --pangenome"
	 elif  [ $GRAPH == "CHMPD" ]
	 then
		  CONTROLS="--pos_control PSEUDOSET --haplo_sample PSEUDOSET --neg_control PSEUDOSET --pangenome"
	 elif [[ $GRAPH =~ ^GIAB ]]
	 then
		  CONTROLS="--pos_control HG002 --pangenome"
	 fi
fi

if [ $GRAPH != "SVPOP" ] && [ $GRAPH != "GIAB-0.5" ] && [ $GRAPH != "GIAB-0.6" ]
then
	 INDEX_OPTS="--all_index"
	 if [ $NO_UNFOLD == 0 ]
	 then
		  INDEX_OPTS="${INDEX_OPTS} --gbwt_prune"
	 fi
else
	 INDEX_OPTS="--xg_index --gcsa_index --id_ranges_index --snarls_index"
fi

if [ $MIN_AF != 0 ]
then
	 INDEX_OPTS="${INDEX_OPTS} --pre_min_af ${MIN_AF}"
fi

if [ $FORCE_OUTSTORE == 1 ]
then
	 INDEX_OPTS="${INDEX_OPTS} --force_outstore"
fi

INDEX_OPTS="${INDEX_OPTS} --alt_paths --alt_path_gam_index"

# run the job
./ec2-run.sh ${HEAD_NODE_OPTS} -m 50 -n i3.8xlarge:${BID},i3.8xlarge "construct aws:${REGION}:${JOBSTORE_NAME} aws:${REGION}:${OUTSTORE_NAME} --fasta ${FASTA} --vcf ${VCFS}  --out_name ${OUT_NAME} --flat_alts ${ALL_INDEX} ${CONTROLS} --normalize ${REGIONS} ${INDEX_OPTS} --merge_graphs --keep_vcfs --validate --handle_svs ${CONFIG_OPTS} --logFile construct.${OUT_NAME}.log ${RESTART_FLAG}" | tee "construct.${OUT_NAME}-$(date).stdout"

aws s3 cp construct.${OUT_NAME}-$(date).stdout s3://${OUTSTORE_NAME}/
