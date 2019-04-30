#!/usr/bin/env bash
# create-leader.sh: launch a Toil leader node with everything needed to run toil-vg
# (basically a subset of  vg/scripts/compare-graphs.sh)
# Tested with Toil 3.13.0.  Will not work on earlier versions

# Need Toil installed
if ! [ -x "$(command -v toil)" ]; then
	 printf "Toil must be installed in order to create leader nodes\n"
	 printf "ex:\n"
	 printf "virtualenv venv\n"
	 printf ". venv/bin/activate\n"
	 printf "pip install -U toil[aws,mesos]==3.16.0\n"
	 exit 1
fi

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [Options] CLUSTER_NAME KEYPAIR_NAME \n"
    printf "Options:\n\n"
    printf "\t-p PACKAGE\tUse the given Python package specifier to install toil-vg (default to current master).\n"
    printf "\t-t CONTAINER\tUse the given Toil container in the cluster (default: ${TOIL_APPLIANCE_SELF}).\n"
    exit 1
}

while getopts "hp:t:" o; do
    case "${o}" in
        p)
            TOIL_VG_PACKAGE="${OPTARG}"
            ;;
        t)
            TOIL_APPLIANCE_SELF="${OPTARG}"
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

CLUSTER_NAME="${1}"
shift
KEYPAIR_NAME="${1}"
shift

# default toil-vg to current master
if [ -z ${TOIL_VG_PACKAGE} ]; then
	 TOIL_VG_PACKAGE="git+https://github.com/vgteam/toil-vg.git"
fi

set -x

if [ ! -z "${TOIL_APPLIANCE_SELF}" ]; then
	 TOIL_APPLIANCE_SELF="${TOIL_APPLIANCE_SELF}" $PREFIX toil launch-cluster "${CLUSTER_NAME}" --leaderNodeType=t2.medium -z us-west-2a "--keyPairName=${KEYPAIR_NAME}"
else
	 $PREFIX toil launch-cluster "${CLUSTER_NAME}" --leaderNodeType=t2.medium -z us-west-2a "--keyPairName=${KEYPAIR_NAME}"
fi

# We need to manually install git to make pip + git work...
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" apt update
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" apt install git -y


# For hot deployment to work, toil-vg needs to be in a virtualenv that can see the system Toil
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" virtualenv --system-site-packages venv

$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install pyyaml
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install aws
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install awscli
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install numpy
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install scipy==1.0.0rc2
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install scikit-learn
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install "${TOIL_VG_PACKAGE}"

# Destroy cluster if toil-vg not running
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg --help 2> /dev/null
if [ $? -ne 0 ]; then
	 printf "toil-vg was not successfully installed.  destroying cluster\n"
    $PREFIX toil destroy-cluster "${CLUSTER_NAME}" -z us-west-2a
	 exit 1
fi

set +x

printf "Cluster ready:\n"
printf "Log in with:\n\n"
printf "toil ssh-cluster --insecure --zone=us-west-2a ${CLUSTER_NAME}\n\n"
printf "Launch toil-vg commands as follows:\n\n"
printf ". venv/bin/activate"
printf "script"
printf "screen #(or screen -r if screen already running)\n"
printf "toil-vg --provisioner=aws --mesosMaster=\$(hostname -i):5050 ...\n\n"
printf "Destroy the cluster with:\n\n"
printf "toil destroy-cluster ${CLUSTER_NAME}\n\n"
