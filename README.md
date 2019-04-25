## Genotyping structural variation in variation graphs with the vg toolkit

This repositiory contains the commands and scripts used for *Genotyping structural variation in variation graphs with the vg toolkit, 2019, in press*.  They are primarily dependent on [toil-vg](https://github.com/vgteam/toil-vg), which can run most other dependencies via Docker.  [Github issues](https://github.com/vgteam/sv-genotyping-paper/issues/new) is the best place to raise questions or concerns.

### Whole genome experiments

These were run on AWS via [Toil](http://toil.ucsc-cgl.org/).  In theory, they could use any other framework that Toil supports, though the scripts will have to be modified accordingly. 

* [HGSVC](whole-genome/toil-experiments/HGSVC.md)
* [GIAB](whole-genome/toil-experiments/GIAB.md)
* [CHMPD](whole-genome/toil-experiments/CHMPD.md)
* [SVPOP](whole-genome/toil-experiments/SVPOP.md)

