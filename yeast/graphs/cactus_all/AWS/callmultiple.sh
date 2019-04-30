declare -a arr=("SRR4074413" "SRR4074412" "SRR4074411" "SRR4074394" "SRR4074385" "SRR4074384" "SRR4074383" "SRR4074358" "SRR4074258" "SRR4074257" "SRR4074256")

## now loop through the above array
for i in "${arr[@]}"
do
    echo "$i"
    MASTER_IP=`ifconfig eth0 |grep "inet addr" |awk '{print $2}' |awk -F: '{print $2}'`
    toil clean aws:us-west-2:vgcall-yeast-cactus-all-jobstore
    toil-vg call --realTimeLogging --realTimeStderr --whole_genome_config --nodeTypes r3.8xlarge:0.53 --maxNodes 20 --defaultPreemptable --provisioner aws --batchSystem mesos --mesosMaster=${MASTER_IP}:5050 --metrics aws:us-west-2:vgcall-yeast-cactus-all-jobstore component0.xg $i.recall.cactus.all aws:us-west-2:vgcall-yeast-cactus-all-outstore --gams $i.mapped.sorted.gam --recall --call_chunk_size 0 --chroms S288C.chrI S288C.chrII S288C.chrIII S288C.chrIV S288C.chrV S288C.chrVI S288C.chrVII S288C.chrVIII S288C.chrIX S288C.chrX S288C.chrXI S288C.chrXII S288C.chrXIII S288C.chrXIV S288C.chrXV S288C.chrXVI 2> $i.recall.cactus.all.log
done

