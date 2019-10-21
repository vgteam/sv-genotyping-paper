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

Notably, some methods outputs the genotypes of all samples in one VCF.
Others sometimes output only the genotyped sample but keeps the other samples from the truth sets (and fill genotypes with "./."). 
These samples disturbed the pipeline, especially variant normalization.
The way we ran some of the tools resulted in some VCF having the filename of the BAM file as sample name and we want to replace them with the actual sample name.
We used the [selectSample.py](../misc-scripts/selectSample.py) script to select the appropriate sample or rename the sample in the output VCF of several methods.


### BayesTyper

```
# HGSVC
zcat HG00514_sim30x_hgsvc_pp_bayestyper_pass_nomis.vcf.gz | bgzip > hgsvcsim-bayestyper-HG00514.vcf.gz
zcat all_hgsvc_pp49_hc49_bayestyper_pass_hgsvc.vcf.gz | bgzip > vcf/hgsvc-bayestyper-HG00514-HG00733-NA19240.vcf.gz
python ../misc-scripts/selectSample.py -s HG00514 -v vcf/hgsvc-bayestyper-HG00514-HG00733-NA19240.vcf.gz | bgzip > vcf/hgsvc-bayestyper-HG00514.vcf.gz
python ../misc-scripts/selectSample.py -s HG00733 -v vcf/hgsvc-bayestyper-HG00514-HG00733-NA19240.vcf.gz | bgzip > vcf/hgsvc-bayestyper-HG00733.vcf.gz
python ../misc-scripts/selectSample.py -s NA19240 -v vcf/hgsvc-bayestyper-HG00514-HG00733-NA19240.vcf.gz | bgzip > vcf/hgsvc-bayestyper-NA19240.vcf.gz

# GIAB
zcat giab_AJ_sv_HG002_all19_pass_giab.vcf.gz | bgzip > vcf/giab5-bayestyper-HG002.vcf.gz
```

### Delly

```
## HGSVC
python ../misc-scripts/selectSample.py -s HG00514 -v HG00514_sim30x.delly.vcf.gz | bgzip > vcf/hgsvcsim-delly-HG00514.vcf.gz
python ../misc-scripts/selectSample.py -s HG00514 -v ERR903030.delly.vcf.gz | bgzip > vcf/hgsvc-delly-HG00514.vcf.gz
python ../misc-scripts/selectSample.py -s HG00733 -v HG00733.delly.vcf.gz | bgzip > vcf/hgsvc-delly-HG00733.vcf.gz
python ../misc-scripts/selectSample.py -s NA19240 -v NA19240.delly.vcf.gz | bgzip > vcf/hgsvc-delly-NA19240.vcf.gz

## GIAB
python ../misc-scripts/selectSample.py -s HG002 -v HG002-giab-0.5.delly.vcf.gz | bgzip > vcf/giab5-delly-HG002.vcf.gz
```

## SVTyper

```
## HGSVC
python ../misc-scripts/selectSample.py -s HG00514 -v HG00514_sim30x.svtyper.explicit.vcf.gz | bgzip > vcf/hgsvcsim-svtyper-HG00514.vcf.gz
python ../misc-scripts/selectSample.py -s HG00514 -v HG00514.svtyper.explicit.vcf.gz | bgzip > vcf/hgsvc-svtyper-HG00514.vcf.gz
python ../misc-scripts/selectSample.py -s HG00733 -v HG00733.svtyper.explicit.vcf.gz | bgzip > vcf/hgsvc-svtyper-HG00733.vcf.gz
python ../misc-scripts/selectSample.py -s NA19240 -v NA19240.svtyper.explicit.vcf.gz | bgzip > vcf/hgsvc-svtyper-NA19240.vcf.gz

## GIAB
python ../misc-scripts/selectSample.py -s HG002 -v ../svtyper/HG002/HG002.svtyper.explicit.vcf.gz | bgzip > vcf/giab5-svtyper-HG002.vcf.gz
```

## Paragraph

```
## HGSVC
python ../misc-scripts/selectSample.py -s HG00514 -v HG00514_sim30x/paragraph_out/genotypes.vcf.gz | bgzip > vcf/hgsvcsim-paragraph-HG00514.vcf.gz
python ../misc-scripts/selectSample.py -s HG00514 -v HG00514/paragraph_out/genotypes.vcf.gz | bgzip > vcf/hgsvc-paragraph-HG00514.vcf.gz
python ../misc-scripts/selectSample.py -s HG00733 -v HG00733/paragraph_out/genotypes.vcf.gz | bgzip > vcf/hgsvc-paragraph-HG00733.vcf.gz
python ../misc-scripts/selectSample.py -s NA19240 -v NA19240/paragraph_out/genotypes.vcf.gz | bgzip > vcf/hgsvc-paragraph-NA19240.vcf.gz

## GIAB
vcf-sort ../paragraph/HG002/paragraph_out/genotypes.vcf.gz | vcfkeepsamples - HG002 | bgzip -c > vcf/giab5-paragraph-HG002.vcf.gz

## SVPOP
python ../misc-scripts/selectSample.py -s HG00514 -v HG00514/paragraph_out/genotypes.vcf.gz | bgzip > vcf/svpop-paragraph-HG00514.vcf.gz
python ../misc-scripts/selectSample.py -s HG00733 -v HG00733/paragraph_out/genotypes.vcf.gz | bgzip > vcf/svpop-paragraph-HG00733.vcf.gz
python ../misc-scripts/selectSample.py -s NA19240 -v NA19240/paragraph_out/genotypes.vcf.gz | bgzip > vcf/svpop-paragraph-NA19240.vcf.gz
```

### SMRT-SV2

The *SVTYPE* and *END* fields were filtered out from the header.
Otherwise the sveval package expects variants in symbolic form.
We used explicit VCFs to allow for normalization with `bcftools norm`.

The VCF for SVPOP contains genotypes for multiple samples so we create links with appropriate file names that match the Snakemake file.

```
## SVPOP
zgrep -v "##INFO=<ID=SVTYPE" svpop-smrtsv-hgsvc-explicit.vcf.gz | grep -v "##INFO=<ID=END" | bgzip > svpop-smrtsv-HG00514-HG00733-NA19240.vcf.gztemp.vcf.gz
python ../misc-scripts/selectSample.py -s HG00514 -v svpop-smrtsv-HG00514-HG00733-NA19240.vcf.gz | bgzip > vcf/svpop-smrtsv-HG00514.vcf.gz
python ../misc-scripts/selectSample.py -s HG00733 -v svpop-smrtsv-HG00514-HG00733-NA19240.vcf.gz | bgzip > vcf/svpop-smrtsv-HG00733.vcf.gz
python ../misc-scripts/selectSample.py -s NA19240 -v svpop-smrtsv-HG00514-HG00733-NA19240.vcf.gz | bgzip > vcf/svpop-smrtsv-NA19240.vcf.gz
```
