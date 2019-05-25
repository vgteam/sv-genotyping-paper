# Whole genome experiments in human

There is one folder for each dataset with the commands to download/prepare the data and genotype SV with vg and the other methods.

* [Human Genome Structural Variation Consortium (HGSVC)](hgsvc)
* [Genome in a Bottle (GiaB)](giab)
* [Pseudo-diploid CHM genome (CHMPD)](chmpd)
* [SV catalog from Audano et al. Cell 2019 (SVPOP)](svpop)

There is also a [`toil-scripts`](toil-scripts) folder with helper scripts to were used to run the analysis on AWS.

The commands for the **evaluation**, using Snakemake, are available in the [`sveval`](sveval) folder.

The VCFs produced produced by vg and the other methods across these datasets are available at [https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=vgsv2019/vcfs/](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=vgsv2019/vcfs/). 
