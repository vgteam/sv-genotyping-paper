# SVPOP

## Data

- [GRCh38 reference genome](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).

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

## SMRTSV v2

```
# Smrtsv2 was run on courtyard, then made explicity using the same process as the original svpop graph (see construg-hgsvc.sh)
# (all three samples are in one VCF: svpop-smrtsv-hgsvc-explicit.vcf.gz)
```
