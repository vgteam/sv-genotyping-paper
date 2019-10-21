The evaluation pipeline uses the [sveval R package](https://github.com/jmonlong/sveval) (and bcftools for variant normalization).
[Snakemake](https://snakemake.readthedocs.io/en/stable/) was used to automate evaluation across the different samples, methods, and parameters (SV presence vs genotype, whole-genome vs non-repeat regions).
For more details see the [Snakefile](Snakefile).

To run the evaluation

```
## folders containing input and output files
mkdir -p vcf rdata tsv bed

## evaluate each dataset
snakemake -p hgsvc
snakemake -p giab5
snakemake -p svpop
snakemake -p chmpd
snakemake -p svpopregions
## or all at once
snakemake -p all

## same but with more stringent matching criterion
snakemake -p stringent

## merge all the TSV files into two TSV files: human-merged-prcurve.tsv and human-merged-persize.tsv 
Rscript mergeTSVs.R

## compute size distribution of SVs in the input catalogs
Rscript refsize.R tsv/human-ref-size.tsv HGSVC.haps.vcf.gz sv-pop-explicit.vcf.gz pseudo_diploid-explicit.vcf.gz giab-0.5.vcf.gz
```

After running the evaluation, we also annotated the deleted/inserted sequence of the SVs to explore the performance across repeat classes.
The results of the evaluation (in the `rdata` folder after running the analysis) were combined with annotated VCFs.
The RepeatMasker annotation code is available in the [`repeatmasked`](repeatmasked) folder.

## Data

For HGSVS, SVPOP and CHMPD datasets:

- [GRCh38 reference genome](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).
- Non-repeat regions in GRCh38: *hg38_non_repeats.bed.gz*.

For the Genome in a bottle dataset:

- Hg19 reference genome: [ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz](ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz)
- Genome in a bottle high-confidence regions BED file (Hg19): [ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed)


## Preparing the VCF files

The VCFs containing genotype prediction for each method on the different datasets were produced with the commands in [each dataset's folder](..).
Those VCFs are also available at [https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=vgsv2019/vcfs/](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=vgsv2019/vcfs/). 

Before evaluation, some VCF files were reformatted or renamed.

Notably, some methods don't output only the genotyped samples but also other samples from the truth sets with "./.". 
These samples disturbed the pipeline, especially variant normalization.
We used the XXX script to select the appropriate sample in the output of several methods.

### BayesTyper

BayesTyper outputs gzipped VCFs. 
Before evaluation, they were bgzipped.

For HGSVC, the VCF contains genotypes for the three samples so we create links with appropriate file names that match the Snakemake file.

```
cd vcf

# HGSVC
zcat HG00514_sim30x_hgsvc_pp_bayestyper_pass_nomis.vcf.gz | bgzip > hgsvcsim-bayestyper-HG00514.vcf.gz
zcat all_hgsvc_pp49_hc49_bayestyper_pass_hgsvc.vcf.gz | bgzip > hgsvc-bayestyper-HG00514.vcf.gz
ln -s hgsvc-bayestyper-HG00514.vcf.gz hgsvc-bayestyper-HG00733.vcf.gz
ln -s hgsvc-bayestyper-HG00514.vcf.gz hgsvc-bayestyper-NA19240.vcf.gz

# GIAB
zcat giab_AJ_sv_HG002_all19_pass_giab.vcf.gz | bgzip > vcf/giab5-bayestyper-HG002.vcf.gz
```

### SMRT-SV2

The *SVTYPE* and *END* fields were filtered out from the header.
Otherwise the sveval package expects variants in symbolic form.
We used explicit VCFs to allow for normalization with `bcftools norm`.

The VCF for SVPOP contains genotypes for multiple samples so we create links with appropriate file names that match the Snakemake file.

```
cd vcf
zgrep -v "##INFO=<ID=SVTYPE" svpop-smrtsv-hgsvc-explicit.vcf.gz | grep -v "##INFO=<ID=END" | bgzip > temp.vcf.gz
mv temp.vcf.gz svpop-smrtsv-HG00514.vcf.gz
ln -s svpop-smrtsv-HG00514.vcf.gz vcf/svpop-smrtsv-HG00733.vcf.gz
ln -s svpop-smrtsv-HG00514.vcf.gz vcf/svpop-smrtsv-NA19240.vcf.gz
```
