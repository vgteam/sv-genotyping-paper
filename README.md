# Genotyping structural variants in pangenome graphs using the vg toolkit


This repository contains the commands and scripts used for [*Genotyping structural variants in pangenome graphs using the vg toolkit, 2019, in press*](https://www.biorxiv.org/content/10.1101/654566v1).
They are primarily dependent on [toil-vg](https://github.com/vgteam/toil-vg), which can run most other dependencies via Docker. 
There is a [WIKI here](https://github.com/vgteam/toil-vg/wiki/Genotyping-Structural-Variants) for genotyping SVs with toil-vg. 
[Github issues](https://github.com/vgteam/sv-genotyping-paper/issues/new) is the best place to raise questions or concerns.

Links to necessary data are also listed for each analysis.

Of note the code to reproduce the figures and tables in the manuscript are available in [the manuscript's repo](https://github.com/jmonlong/manu-vgsv/tree/master/figures).

## Simulation experiment

The different methods were compared on simulated sequence and SVs.
Different depth were tested. 
We also tested the effect of errors in the breakpoint location of the SVs.
Scripts are available and described in the [`simulation`](simulation) folder.

## Whole genome experiments in human

These were run on AWS via [Toil](http://toil.ucsc-cgl.org/). 
In theory, they could use any other framework that Toil supports, though the scripts will have to be modified accordingly. 

In the [`human`](human) directory, there is one folder for each dataset with the commands to download/prepare the data and genotype SV with vg and the other methods.

* [Human Genome Structural Variation Consortium (HGSVC)](human/hgsvc)
* [Genome in a Bottle (GiaB)](human/giab)
* [Pseudo-diploid CHM genome (CHMPD)](human/chmpd)
* [SV catalog from Audano et al. Cell 2019 (SVPOP)](human/svpop)

There is also a [`toil-scripts`](human/toil-scripts) folder with helper scripts that were used to run the analysis on AWS.

The commands for the evaluation, using [Snakemake](https://snakemake.readthedocs.io/en/stable/), are available in the [`sveval`](human/sveval) folder.

The VCFs produced produced by vg and the other methods across these datasets are available at [https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=vgsv2019/vcfs/](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=vgsv2019/vcfs/). 

## De-novo assembly graph experiments in yeast

The yeast experiments are written as [snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines.
Each pipeline consists of a set of rules that process a set of input files into a set of output files.

In the [`yeast`](yeast) directory, there are several folders for the different phases of the experiment as well as detailed descriptions on how to re-run it.

