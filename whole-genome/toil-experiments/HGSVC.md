### HGSVG Whole Genome Experiment

## Setup
The following environment variables need to be set
```
# jobstore prefix: EX my-jobstore
export JOBSTORE=
# outstore: EX my-outstore/HGSVC
export OUTSTORE=
# cluster name prefix EX: my-cluster where my-cluster1 and my-cluster2 were already created with toil-vg/scripts/create-ec2-leader.sh
export CLUSTER=
# BED file for evaluation
export COMPARE_REGIONS_BED=
# interleaved fq template for simulation error model
export TEMPLATE_FQ=
# input reads root location (will look for reads in subdirs)
export FQBASE=
```

## Running
```
# Construct all graphs and indexes.
./construct.sh -c ${CLUSTER}1 ${JOBSTORE}1 ${OUTSTORE}

# Simulate a GAM and fastq from the two HG00514 haplotypes
./simulate.sh -c ${CLUSTER}1 ${JOBSTORE}1 ${OUTSTORE}/sim s3://${OUTSTORE}/HGSVC_HG00514_haplo_thread_0.xg s3://${OUTSTORE}/HGSVC_HG00514_haplo_thread_1.xg ${TEMPLATE_FQ}

## Mapping, Calling and Evaluation.
## Use -M SKIP, -C SKIP, -E SKIP to bypass mapping, calling, evaluation respectively (ex, if rerunning a step)

# Simulation
./mce.sh -c ${CLUSTER}1  ${JOBSTORE}1 ${OUTSTORE}/HGSVC s3://${OUTSTORE}/HGSVC/HGSVC HG00514 HG00514-sim s3://${OUTSTORE}/HGSVC/HGSVC.haps.vcf.gz ${COMPARE_REGIONS_BED} s3://${OUTSTORE}/HGSVC-chroms-dec5/sim/sim-HG00514-30x.fq.gz

# Three HGSVC Samples (the reads can be found publicly by looking up the run names on EBI's ENA.  I had them mirrored in FQBASE for faster access)
./mce.sh -c ${CLUSTER}1  ${JOBSTORE}1 ${OUTSTORE}/HGSVC s3://${OUTSTORE}/HGSVC/HGSVC HG00514 HG00514 s3://${OUTSTORE}/HGSVC/HGSVC.haps.vcf.gz ${COMPARE_REGIONS_BED} ${FQBASE}/HG00514/ERR903030_1.fastq.gz ${FQBASE}/HG00514/ERR903030_2.fastq.gz 

./mce.sh -c ${CLUSTER}2  ${JOBSTORE}2 ${OUTSTORE}/HGSVC s3://${OUTSTORE}/HGSVC/HGSVC HG00733 HG00733 s3://${OUTSTORE}/HGSVC/HGSVC.haps.vcf.gz ${COMPARE_REGIONS_BED} ${FQBASE}/HG00733/ERR895347_1.fastq.gz ${FQBASE}/HG00733/ERR895347_2.fastq.gz 

./mce.sh -c ${CLUSTER}3  ${JOBSTORE}3 ${OUTSTORE}/HGSVC s3://${OUTSTORE}/HGSVC/HGSVC NA19240 NA19240 s3://${OUTSTORE}/HGSVC/HGSVC.haps.vcf.gz ${COMPARE_REGIONS_BED} ${FQBASE}/NA19240/ERR894724_1.fastq.gz ${FQBASE}/NA19240/ERR894724_2.fastq.gz 
```

#commands used to download the data

```
for name in HG00514 HG00733 NA19240 HG00514-sim
do
aws s3 sync s3://${OUTSTORE}/HGSVC/eval-${name} ./HGSVC-jan5-eval-${name}
done
```


