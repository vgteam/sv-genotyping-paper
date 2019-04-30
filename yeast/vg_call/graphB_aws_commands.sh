#Activate virtual environment if not already done
source toilvenv/bin/activate
ssh-add

#Create leader node
yeast_sv/vg_call/create-ec2-leader.sh vgtoil <KEYPAIRNAME>

#Move graph, alignment GAMs and bash script to leader node
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/AWS/callmultiple.sh :/
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/exploded/component0.xg :/
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/mappings/SRR4074257.mapped.sorted.gam :/
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/mappings/SRR4074358.mapped.sorted.gam :/
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/mappings/SRR4074383.mapped.sorted.gam :/
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/mappings/SRR4074384.mapped.sorted.gam :/
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/mappings/SRR4074412.mapped.sorted.gam :/
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/mappings/SRR4074413.mapped.sorted.gam :/
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/mappings/SRR4074256.mapped.sorted.gam :/
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/mappings/SRR4074258.mapped.sorted.gam :/
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/mappings/SRR4074385.mapped.sorted.gam :/
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/mappings/SRR4074394.mapped.sorted.gam :/
toil rsync-cluster -z us-west-2a vgtoil -avP graphs/cactus_all/mappings/SRR4074411.mapped.sorted.gam :/

#Connect to leader node
toil ssh-cluster --insecure -z us-west-2a vgtoil

#On leader node:
screen
. venv/bin/activate
#Get latest toil-vg
#pip install -U git+https://github.com/vgteam/toil-vg.git
bash callmultiple.sh

#Destroy cluster
toil destroy-cluster -z us-west-2a vgtoil
