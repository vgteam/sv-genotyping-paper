# De-novo assembly graph experiments in yeast

This directory contains the snakemake pipelines and other scripts for SV genotyping using [vg](vgteam/vg) on 12 yeast strains.
Primarily, the SV genotyping performance of two different vg graphs is compared:

- **VCF graph**: Created with `vg construct` using S288C as reference strain and VCFs with SVs from Assemblytics (Nattestad et al.), AsmVar (Huang et al.) and paftools (Li) as variants
- **Cactus graph**: Created from multiple genome alignment of yeast strains using [cactus](https://github.com/ComparativeGenomicsToolkit/cactus) and [hal2vg](https://github.com/ComparativeGenomicsToolkit/hal2vg)

The experiments comprise several phases which are explained in the following:

## 1. Prepare yeast assemblies and call SVs relative to the reference strain S288C

```
git clone https://github.com/vgteam/sv-genotyping-paper.git
cd yeast/assemblies
snakemake
```

This will clone this repo including the contained PacBio assemblies from the [Yeast Population Reference Panel (YPRP)](https://yjx1217.github.io/Yeast_PacBio_2016/welcome/).
Subsequently, it will run 3 SV calling pipelines on the 11 non-reference strains:
- [Minimap2](https://github.com/lh3/minimap2) + paftools call
- [LAST](http://last.cbrc.jp/) + [AsmVar](https://github.com/bioinformatics-centre/AsmVar)
- [nucmer](https://github.com/mummer4/mummer) + [Assemblytics](https://github.com/marianattestad/assemblytics)

Because the resulting callsets differ considerably, the union of all 3 callsets is produced (high-sensitivy callset).
Finally, the high-sensitivity callsets of all 11 non-reference strains are merged.

## 2. Download Illumina reads for the 12 yeast strains

```
cd yeast/illumina_reads
snakemake
```

## 3. Create, index and map Illumina reads to the VCF graph

There are 2 different subdirectories for the VCF graph in `yeast/graphs`:
- `yeast/graphs/constructunion_all` - Integrate variants from all strains into graph
- `yeast/graphs/constructunion_five` - Integrate variants from 3 S. cerevisiae and 2 S. paradoxus strains into the graph

Let's call them the two strain sets. For the first strain set (`all`), run:
```
cd yeast/graphs/constructunion_all
snakemake
```

For the second strain set (`five`), run:
```
cd yeast/graphs/constructunion_five
snakemake
```

This will create the VCF graphs for the two strain sets using S288C as reference strain and the high-sensitivity variant callsets produced before. 

As a result `yeast/graphs/constructunion_*/mappings` will contain the sorted GAM alignments of the Illumina reads against the graph.


## 4. Run cactus to produce an assembly-derived graph

The cactus genome-genome aligner is run using toil.
To install toil in a virtual environment:
```
pip install virtualenv
virtualenv cactus_env
source cactus_env/bin/activate
pip install --upgrade toil
```

Like for the VCF graph, there are 2 different strain sets for the cactus graph: `all` and `five`. To create cactus graphs for both of them, now follow the steps in `yeast/cactus/all/aws_commands.sh` and `yeast/cactus/five/aws_commands.sh`.


## 5. Create, index and map Illumina reads to the cactus graph

```
cd yeast/graphs/cactus_all
snakemake
```

```
cd yeast/graphs/cactus_five
snakemake
```

This will a) create each of the two cactus graphs from the cactus alignments produced in the previous step and b) map the Illumina reads to them.


As a result `yeast/graphs/cactus_*/mappings` will contain the sorted GAM alignments of the Illumina reads against the graphs.
Additionally, it will produce png plots of mapping quality and identity in `yeast/graphs/cactus_*/mappings/stats`.


## 6. Genotype SVs on the VCF graph and the cactus graph

The SV genotyping is run using toil-vg.
First, install toil in a separate virtual environment:
```
virtualenv toilvenv
source toilvenv/bin/activate
pip install toil[aws,mesos]==3.18.0
pip install toil-vg
```

For the VCF graph, follow the steps in `yeast/vg_call/graphA_aws_commands.sh`. It contains commands for the `all` strain set but can be easily modified for the `five` strain set.
For the cactus graph, follow the steps in `yeast/vg_call/graphB_aws_commands.sh`. It contains commands for the `all` strain set but can be easily modified for the `five` strain set.

The compressed vcf results are written into an outstore on S3 of the form `aws:us-west-2:vgcall-yeast-<graph>-<strain-set>-outstore`


## 7. Evaluate genotyping performance

To evaluate the SV genotype calls made by `vg call` you need to download them from the S3 outstore to the respective local directory `yeast/evaluation/calls/<graph>_<strain-set>/vcf`. This can be done via the S3 web interface or commands of the form `aws s3 cp s3://vgcall-yeast-<graph>-<strain-set>-outstore/<sample>*.vcf.gz yeast/evaluation/calls/<graph>_<strain-set>/vcf`.

Now, the evaluation pipeline can be started with:
```
source toilvenv/bin/activate
cd yeast/evaluation/
snakemake
```

It will produce pdf plots showing average deltas in mapping identity and mapping quality of reads aligned to the sample graphs in `yeast/evaluation/results/assemblyeval_svonly`.
