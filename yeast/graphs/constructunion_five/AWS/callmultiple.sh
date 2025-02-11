declare -a arr=("SRR4074256" "SRR4074257" "SRR4074258" "SRR4074358" "SRR4074383" "SRR4074384" "SRR4074385" "SRR4074394" "SRR4074411" "SRR4074412" "SRR4074413")

## now loop through the above array
for i in "${arr[@]}"
do
    echo "$i"
    MASTER_IP=`ifconfig eth0 |grep "inet addr" |awk '{print $2}' |awk -F: '{print $2}'`
    toil clean aws:us-west-2:vgcall-yeast-constructunion-five-jobstore
    toil-vg call --realTimeLogging --realTimeStderr --whole_genome_config --nodeTypes r3.8xlarge:0.53 --maxNodes 20 --defaultPreemptable --provisioner aws --batchSystem mesos --mesosMaster=${MASTER_IP}:5050 --metrics aws:us-west-2:vgcall-yeast-constructunion-five-jobstore construct.xg $i.recall.constructunion.five aws:us-west-2:vgcall-yeast-constructunion-five-outstore --gams $i.mapped.sorted.gam --recall --call_chunk_size 0 --chroms chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI 2> $i.recall.constructunion.five.log
done
