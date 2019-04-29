# Genotyping structural variation in variation graphs with the vg toolkit

This repositiory contains the commands and scripts used for *Genotyping structural variation in variation graphs with the vg toolkit, 2019, in press*. 
They are primarily dependent on [toil-vg](https://github.com/vgteam/toil-vg), which can run most other dependencies via Docker. 
[Github issues](https://github.com/vgteam/sv-genotyping-paper/issues/new) is the best place to raise questions or concerns.

## Whole genome experiments in human

These were run on AWS via [Toil](http://toil.ucsc-cgl.org/). 
In theory, they could use any other framework that Toil supports, though the scripts will have to be modified accordingly. 

In the [`human`](human) directory, there is one folder for each dataset with the commands to download/prepare the data and genotype SV with vg and the other methods.

* [Human Genome Structural Variation Consortium (HGSVC)](human/hgsvc)
* [Genome in a Bottle (GiaB)](human/giab)
* [Pseudo-diploid CHM genome (CHMPD)](human/chmpd)
* [SV catalog from Audano et al. Cell 2019 (SVPOP)](human/svpop)

There is also a [`toil-scripts`](human/toil-scripts) folder with helper scripts to were used to run the analysis on AWS.

## From de novo assembly alignment in yeast 

*Add scripts in the [`yeast`](yeast) folder.*
