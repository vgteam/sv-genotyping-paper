#Activate virtual environment if not already done
source cactus_env/bin/activate
ssh-add

#Launch cluster for cactus
toil launch-cluster -z us-west-2a cactus --keyPairName <KEYPAIRNAME> --leaderNodeType t2.medium

#Copy repeat-masked assemblies and seqfile onto leader node
toil rsync-cluster -z us-west-2a cactus -avP yeast/cactus/all/seqfile_yeast.txt yeast/assemblies/assemblies_repeatmasked/CBS432.genome.fa.masked yeast/assemblies/assemblies_repeatmasked/DBVPG6044.genome.fa.masked yeast/assemblies/assemblies_repeatmasked/DBVPG6765.genome.fa.masked yeast/assemblies/assemblies_repeatmasked/N44.genome.fa.masked yeast/assemblies/assemblies_repeatmasked/S288C.genome.fa.masked yeast/assemblies/assemblies_repeatmasked/SK1.genome.fa.masked yeast/assemblies/assemblies_repeatmasked/UFRJ50816.genome.fa.masked yeast/assemblies/assemblies_repeatmasked/UWOPS034614.genome.fa.masked yeast/assemblies/assemblies_repeatmasked/UWOPS919171.genome.fa.masked yeast/assemblies/assemblies_repeatmasked/Y12.genome.fa.masked yeast/assemblies/assemblies_repeatmasked/YPS128.genome.fa.masked yeast/assemblies/assemblies_repeatmasked/YPS138.genome.fa.masked :/

#Connect to leader node
toil ssh-cluster -z us-west-2a cactus

#On leader node:
apt update
apt install -y git tmux
virtualenv --system-site-packages venv
source venv/bin/activate
git clone https://github.com/comparativegenomicstoolkit/cactus.git
cd cactus
pip install --upgrade .
cd /
screen
cactus --nodeTypes c5.4xlarge:0.4,r4.2xlarge --minNodes 0,0 --maxNodes 1,1 --provisioner aws --batchSystem mesos --metrics aws:us-west-2:cactus seqfile_yeast.txt cactusoutput.hal
#Disconnect from leader node

#Copy cactus alignment result from leader node to the vg mapping pipeline directory
toil rsync-cluster -z us-west-2a cactus -avP :/cactusoutput.hal yeast/graphs/cactus_all/
#Destroy cluster
toil destroy-cluster -z us-west-2a cactus