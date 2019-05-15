The evaluation pipeline uses `toil-vg vcfeval` and, internally, the [sveval R package](https://github.com/jmonlong/sveval).
[Snakemake](https://snakemake.readthedocs.io/en/stable/) was used to automate evaluation across the different samples, methods, and parameters (SV presence vs genotype, whole-genome vs non-repeat regions).
For more details see the [Snakemake file](Snakefile).

To run the evaluation

```
# activate toil-vg environment
source toilvenv/bin/activate

snakemake -p hgsvc
snakemake -p giab5
snakemake -p svpop
snakemake -p chmpd
```

## Data

For HGSVS, SVPOP and CHMPD datasets:

- [GRCh38 reference genome](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).
- Non-repeat regions in GRCh38: *hg38_non_repeats.bed.gz*.

For the Genome in a bottle dataset:

- [Hg19 reference genome](ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz).
- Genome in a bottle [high-confidence regions BED file](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed) (Hg19).


## Preparing the VCF files

Before evaluation, some VCF files were reformatted or renamed.

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
