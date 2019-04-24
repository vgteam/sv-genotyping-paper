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
# Construct all graphs and indexes.  Haploid vcfs from ../haps/haps.urls must have been downloaded into ../haps
./construct-hgsvc.sh -c ${CLUSTER}1 ${JOBSTORE}1 ${OUTSTORE}

# Simulate a GAM and fastq from the two HG00514 haplotypes
./simulate-hgsvc.sh -c ${CLUSTER}1 ${JOBSTORE}1 ${OUTSTORE}/sim s3://${OUTSTORE}/HGSVC.chroms_HG00514_haplo_thread_0.xg s3://${OUTSTORE}/HGSVC.chroms_HG00514_haplo_thread_1.xg ${TEMPLATE_FQ}

## Mapping, Calling and Evaluation.
## Use -M SKIP, -C SKIP, -E SKIP to bypass mapping, calling, evaluation respectively (ex, if rerunning a step)

# Simulation
./mce-hgsvc.sh -c ${CLUSTER}1  ${JOBSTORE}1 ${OUTSTORE}/HGSVC-jan5 s3://${OUTSTORE}/HGSVC-jan5/HGSVC HG00514 HG00514-sim s3://${OUTSTORE}/HGSVC-jan5/HGSVC.haps.vcf.gz ${COMPARE_REGIONS_BED} s3://${OUTSTORE}/HGSVC-chroms-dec5/sim/sim-HG00514-30x.fq.gz

# Three HGSVC Samples
./mce-hgsvc.sh -c ${CLUSTER}1  ${JOBSTORE}1 ${OUTSTORE}/HGSVC-jan5 s3://${OUTSTORE}/HGSVC-jan5/HGSVC HG00514 HG00514 s3://${OUTSTORE}/HGSVC-jan5/HGSVC.haps.vcf.gz ${COMPARE_REGIONS_BED} ${FQBASE}/HG00514/ERR903030_1.fastq.gz ${FQBASE}/HG00514/ERR903030_2.fastq.gz 

./mce-hgsvc.sh -c ${CLUSTER}2  ${JOBSTORE}2 ${OUTSTORE}/HGSVC-jan5 s3://${OUTSTORE}/HGSVC-jan5/HGSVC HG00733 HG00733 s3://${OUTSTORE}/HGSVC-jan5/HGSVC.haps.vcf.gz ${COMPARE_REGIONS_BED} ${FQBASE}/HG00733/ERR895347_1.fastq.gz ${FQBASE}/HG00733/ERR895347_2.fastq.gz 

./mce-hgsvc.sh -c ${CLUSTER}3  ${JOBSTORE}3 ${OUTSTORE}/HGSVC-jan5 s3://${OUTSTORE}/HGSVC-jan5/HGSVC NA19240 NA19240 s3://${OUTSTORE}/HGSVC-jan5/HGSVC.haps.vcf.gz ${COMPARE_REGIONS_BED} ${FQBASE}/NA19240/ERR894724_1.fastq.gz ${FQBASE}/NA19240/ERR894724_2.fastq.gz 
```

#commands used to download the data

```
for name in HG00514 HG00733 NA19240 HG00514-sim
do
aws s3 sync s3://${OUTSTORE}/HGSVC-jan5/eval-${name} ./HGSVC-jan5-eval-${name}
done
```

#commands for bayestyper evaluation

```
./eval-hgsvc.sh  -c ${CLUSTER}1 ${JOBSTORE}1  ${OUTSTORE}/HGSVC-jan5/eval-HG00514-bayestyper-full s3://${OUTSTORE}/HGSVC-jan5/HGSVC.haps_HG00514.vcf.gz s3://${OUTSTORE}/HGSVC-Bayestyper/ERR903030_hgsvc_platypus_bayestyper_pass_nomis_maxgpp.vcf.gz ${COMPARE_REGIONS_BED}

./eval-hgsvc.sh  -c ${CLUSTER}2 ${JOBSTORE}2  ${OUTSTORE}/HGSVC-jan5/eval-HG00514-sim-bayestyper-full s3://${OUTSTORE}/HGSVC-jan5/HGSVC.haps_HG00514.vcf.gz s3://${OUTSTORE}/HGSVC-Bayestyper/HG00514_sim30x_hgsvc_platypus_bayestyper_pass_nomis_maxgpp.vcf.gz ${COMPARE_REGIONS_BED}

./eval-hgsvc.sh  -c ${CLUSTER}2 ${JOBSTORE}2  ${OUTSTORE}/HGSVC-jan5/eval-HG00514-sim-bayestyper-feb20 s3://${OUTSTORE}/HGSVC-jan5/HGSVC.haps_HG00514.vcf.gz s3://${OUTSTORE}/HGSVC-Bayestyper/HG00514_sim30x_hgsvc_bayestyper_pass_nomis-feb20.vcf.gz ${COMPARE_REGIONS_BED}

./eval-hgsvc.sh  -c ${CLUSTER}2 ${JOBSTORE}2  ${OUTSTORE}/HGSVC-jan5/eval-HG00514-bayestyper-manta-feb20 s3://${OUTSTORE}/HGSVC-jan5/HGSVC.haps_HG00514.vcf.gz s3://${OUTSTORE}/HGSVC-Bayestyper/ERR903030_hgsvc_pp_hc_mt_bayestyper_pass_nomis_feb20.vcf.gz ${COMPARE_REGIONS_BED}

./eval-hgsvc.sh  -c ${CLUSTER}2 ${JOBSTORE}2  ${OUTSTORE}/HGSVC-jan5/eval-HG00514-bayestyper-feb20 s3://${OUTSTORE}/HGSVC-jan5/HGSVC.haps_HG00514.vcf.gz s3://${OUTSTORE}/HGSVC-Bayestyper/ERR903030_hgsvc_pp_bayestyper_pass_nomis_feb20.vcf.gz ${COMPARE_REGIONS_BED}

aws s3 sync s3://${OUTSTORE}/HGSVC-jan5/eval-HG00514-bayestyper-full ./HGSVC-jan5-eval-HG00514-bayestyper
aws s3 sync s3://${OUTSTORE}/HGSVC-jan5/eval-HG00514-sim-bayestyper-full ./HGSVC-jan5-eval-HG00514-sim-bayestyper
aws s3 sync s3://${OUTSTORE}/HGSVC-jan5/eval-HG00514-sim-bayestyper-feb20 ./HGSVC-jan5-eval-HG00514-sim-bayestyper-feb20
aws s3 sync s3://${OUTSTORE}/HGSVC-jan5/eval-HG00514-bayestyper-feb20 ./HGSVC-jan5-eval-HG00514-bayestyper-feb20
aws s3 sync s3://${OUTSTORE}/HGSVC-jan5/eval-HG00514-bayestyper-manta-feb20 ./HGSVC-jan5-eval-HG00514-bayestyper-manta-feb20
```

