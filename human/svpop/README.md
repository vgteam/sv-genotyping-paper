# SVPOP

## Data

- [GRCh38 reference genome](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).
- SVPOP catalog from [Audano et al. Cell 2019](https://doi.org/10.1016/j.cell.2018.12.019).

The VCF file with the SV catalog (`sv-pop-explicit.vcf.gz`) was prepared with:

```
wget -nc http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181025_EEE_SV-Pop_1/VariantCalls_EEE_SV-Pop_1/EEE_SV-Pop_1.ALL.sites.20181204.bed.gz
wget -nc http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181025_EEE_SV-Pop_1/VariantCalls_EEE_SV-Pop_1/EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz
wget -nc http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181025_EEE_SV-Pop_1/VariantCalls_EEE_SV-Pop_1/EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz.tbi

# fixes up some basic info tags like AC and AF
vcfkeepinfo EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz SVTYPE | vcffixup - | bgzip > EEE_SV-Pop_1.ALL.sites.20181204.fix.vcf.gz
./vcf-add-bed-seqs.py --inv leave EEE_SV-Pop_1.ALL.sites.20181204.fix.vcf.gz EEE_SV-Pop_1.ALL.sites.20181204.bed.gz | bgzip > sv-pop-explicit.vcf.gz
tabix -f -p vcf sv-pop-explicit.vcf.gz
```

### Reads

- HG00514: [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR903/ERR903030/ERR903030_1.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR903/ERR903030/ERR903030_1.fastq.gz) [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR903/ERR903030/ERR903030_2.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR903/ERR903030/ERR903030_2.fastq.gz)
- HG00733: [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_1.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_1.fastq.gz) [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_2.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_2.fastq.gz)
- NA19240: [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_1.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_1.fastq.gz) [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_2.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR894/ERR894724/ERR894724_2.fastq.gz)


## toil-vg

Using the helper scripts from `../toil-scripts`.

```
# Construct the graph and index for svpop (15 sv samples) including inversions
./construct.sh -p -i -c ${CLUSTER}2 ${JOBSTORE}2 ${OUTSTORE}/SVPOP-jan10

# Mapping and calling

./mapcall.sh -c ${CLUSTER}1 ${JOBSTORE}1 ${OUTSTORE}/SVPOP-jan10 s3://${OUTSTORE}/SVPOP-jan10/SVPOP HG00514 HG00514 ${FQBASE}/HG00514/ERR903030_1.fastq.gz ${FQBASE}/HG00514/ERR903030_2.fastq.gz

./mapcall.sh -c ${CLUSTER}2 ${JOBSTORE}2 ${OUTSTORE}/SVPOP-jan10 s3://${OUTSTORE}/SVPOP-jan10/SVPOP HG00733 HG00733 ${FQBASE}/HG00733/ERR895347_1.fastq.gz ${FQBASE}/HG00733/ERR895347_2.fastq.gz 

./mapcall.sh -c ${CLUSTER}3 ${JOBSTORE}3 ${OUTSTORE}/SVPOP-jan10 s3://${OUTSTORE}/SVPOP-jan10/SVPOP NA19240 NA19240 ${FQBASE}/NA19240/ERR894724_1.fastq.gz ${FQBASE}/NA19240/ERR894724_2.fastq.gz

# Download VCF files
rm -rf ./svpop-vg-HG00514.vcf.gz ; aws s3 sync s3://${OUTSTORE}/SVPOP-jan10/call-HG00514/HG00514.vcf.gz ./svpop-vg-HG00514.vcf.gz
rm -rf ./svpop-vg-HG00733.vcf.gz ; aws s3 sync s3://${OUTSTORE}/SVPOP-jan10/call-HG00733/HG00733.vcf.gz ./svpop-vg-HG00733.vcf.gz
rm -rf ./svpop-vg-NA19240.vcf.gz ; aws s3 sync s3://${OUTSTORE}/SVPOP-jan10/call-NA19240/NA19240.vcf.gz ./svpop-vg-NA19240.vcf.gz
```

## Paragraph

Padding sequence was added when needed to avoid errors:

```
zcat sv-pop-explicit.vcf.gz | python ../misc-scripts/addMissingPaddingHg38.py | bgzip -c > sv-pop-explicit-padded.vcf.gz
```

Each sample was genotyped individually.
For example for HG00514:

```
echo -e "id\tpath\tdepth\tread length\nHG00514\tHG00514.bam\t30\t150" > samples_for_paragraph_HG00514.txt
multigrmpy.py -m samples_for_paragraph_HG00514.txt -i sv-pop-explicit-padded.vcf.gz -r hg38.fa -t 10 -o paragraph_out_HG00514
```

The output genotypes are in `paragraph_out_HG00514/genotypes.vcf.gz`.

## SMRT-SV v2

SMRT-SV was obtained from [github](https://github.com/EichlerLab/smrtsv2) using commit 43d28cbd4e76ed2a4e820ded11a80592675db6c9 (the then-current master branch).
It was then run as described in its [README](https://github.com/EichlerLab/smrtsv2/blob/master/GENOTYPE.md), on a single server
```
${SMRTSV_DIR}/smrtsv --jobs 20 genotype genotyper.json variants.vcf.gz

```
Where genotyper.json contained (the VCF and BAM were downloaded from the FTP site listed above (ftp.1000genomes.ebi.ac.uk))
```
{
  "sample_manifest": "hgsvc-samples.tab",
  "sv_reference": "hg38.fa",
  "sv_calls": "EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz",
  "sv_contigs": "EEE_SV-Pop_1.ALL.sites.20181204.bam",
  "model": "30x-4",
  "min_call_depth": 8
}

```
and hsvc-samples.tab contained:
```
SAMPLE	SEX	DATA
HG00514	F	hgsvc-bams/ERR903030_bwa_mem_sort_rg_md.bam
HG00733	F	hgsvc-bams/ERR895347_bwa_mem_sort_rg_md.bam
NA19240	F	hgsvc-bams/ERR894724_bwa_mem_sort_rg_md.bam
```

The BAMs were made as described [here](https://github.com/vgteam/sv-genotyping-paper/tree/master/human/hgsvc)

The VCFs output by this version of SMRT-SV v2 are invalid due to a bug which outputs the wrong values in the REF column.  This was corrected by re-extracting the REF values from the FASTA using [this script](https://github.com/vgteam/sv-genotyping-paper/blob/master/human/chmpd/make-explicit.py) on variants.vcf.gz.  At the same time, the symbolic alleles were converted to explict alleles.  

```
./make-explicit.py variants.vcf.gz --fasta hg38.fa.gz  | vcfkeepinfo - NA | vcffixup - | bgzip > svpop-smrtsv-explicit.vcf.gz
```
