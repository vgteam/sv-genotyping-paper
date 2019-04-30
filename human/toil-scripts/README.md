These scripts are used to automate the analysis of each dataset on AWS.
The commands for each dataset are located in their dedicated folder (e.g. `../hgsvc`).

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
# input reads 1
export FQ1=
# input reads 2
export FQ2=
```

## Scripts

- `construct.sh` constructs a graph from a VCF.
- `simulate.sh` simulates 30X coverage from the graph.
- `mce.sh` maps reads to a graph.
- `call.sh` call SVs on a graph.
