# GIAB

## Data

- Hg19 reference genome: [ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz](ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz)
- GIAB v0.5.0 VCF: [ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_UnionSVs_12122017/svanalyzer_union_171212_v0.5.0_annotated.vcf.gz](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_UnionSVs_12122017/svanalyzer_union_171212_v0.5.0_annotated.vcf.gz)

The SV catalog VCF file was prepared with:

```
wget -nc ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_UnionSVs_12122017/svanalyzer_union_171212_v0.5.0_annotated.vcf.gz
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_UnionSVs_12122017/svanalyzer_union_171212_v0.5.0_annotated.vcf.gz.tbi
vcfkeepinfo svanalyzer_union_171212_v0.5.0_annotated.vcf.gz NA | vcffixup - | bgzip > giab-0.5.vcf.gz
tabix -f -p vcf giab-0.5.vcf.gz
```

## Reads

- [Original reads mapped to hg19](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam)

The reads were downsampled with samtools to 50x to create `HG002-NA24385-50x.bam`:

## toil-vg

Using the helper scripts from `../toil-scripts`.

```
# construct a 0.50 graph with HG002 control
./construct.sh  -G -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}/GIAB-0.5

# map our downsampled 50X reads and call variants
./map.sh -c ${CLUSTER} ${JOBSTORE} ${OUTSTORE}/GIAB-0.5/map-HG002 s3://${OUTSTORE}/GIAB-0.5/GIAB HG002 ${FQBASE}/HG002-NA24385-50x.bam

./call.sh -c ${CLUSTER} -v s3://${OUTSTORE}/GIAB-0.5/giab-0.5.vcf.gz -s s3://${OUTSTORE}/GIAB-0.5/GIAB.snarls -l s3://${OUTSTORE}/GIAB-0.5/GIAB_alts.gam ${JOBSTORE} ${OUTSTORE}/GIAB-0.5/call-HG002 s3://${OUTSTORE}/GIAB-0.5/GIAB.xg HG002 s3://${OUTSTORE}/GIAB-0.5/map-HG002/HG002_chr

```

## Delly

```
delly call -g hs37d5.fa -v giab-0.5.vcf.gz -o temp.bcf HG002-NA24385-50x.bam
bcftools convert temp.bcf > giab5-delly-HG002.vcf
bgzip giab5-delly-HG002.vcf
tabix giab5-delly-HG002.vcf.gz
```

## SVTyper

The explicit VCF was first converted into a symbolic VCF, for deletions only because SVTyper cannot genotype insertions.
See *VCFtoSymbolicDEL.py* in [../misc-scripts](../misc-scripts).

```
zcat giab-0.5.vcf.gz | python VCFtoSymbolicDEL.py > giab-0.5.symbolic.dels.vcf
```

Then:

```
svtyper -i giab-0.5.symbolic.dels.vcf -B HG002-NA24385-50x.bam -l libinfo.json > HG002.svtyper.vcf
vcf-sort HG002.svtyper.vcf | vcfkeepinfo - END SVLEN SVTYPE NS AN | grep -v 'contig' > HG002.svtyper.clean.vcf
bayesTyperTools convertAllele -v HG002.svtyper.clean.vcf -g ../../haps/hg38.fa --keep-imprecise -o giab5-svtyper-HG002
vcfkeepinfo giab5-svtyper-HG002.vcf NS AN | bgzip > giab5-svtyper-HG002.vcf.gz
tabix -f giab5-svtyper-HG002.vcf.gz
```

## BayesTyper

### kmer counting

```
kmc -v -k55 -ci1 -fbam -t12 HG002_NA24385_50x.bam HG002_NA24385_50x_k55 .
bayesTyperTools makeBloom -p 12 -k HG002_NA24385_50x_k55
```

### SNVs/indels 

SNVs downloaded from GIAB releases, normalized and filtered to keep only variants smaller than 50 bp.

```
GIAB_FTP="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh37/supplementaryFiles"
wget ${GIAB_FTP}/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_all.vcf.gz
bcftools view -r X,Y -O z ${GIAB_FTP}/inputvcfsandbeds/HG002_GRCh37_CHROM1-MT_novoalign_Ilmn250x250_FB.vcf.gz > HG002_GRCh37_CHROM1-MT_novoalign_Ilmn250x250_FB_xy.vcf.gz
bcftools view -r X,Y -O z ${GIAB_FTP}/inputvcfsandbeds/HG002_GRCh37_CHROM1-MT_novoalign_Ilmn250x250_GATKHC.vcf.gz > HG002_GRCh37_CHROM1-MT_novoalign_Ilmn250x250_GATKHC_xy.vcf.gz
cp giab-0.5.vcf.gz AJ_GRCh37_GIAB_sv_0.5.vcf.gz
rm *.tbi

# Variant normalization
for VCF_FILE in *.vcf.gz; do
    VCF_FILENAME=$(basename ${VCF_FILE})
    VCF_PREFIX=${VCF_FILENAME%.vcf.gz}
    echo ${VCF_PREFIX}
    zcat ${VCF_PREFIX}.vcf.gz | cut -f1-8 | bcftools norm -c x -m -any -f hs37d5.fa | bcftools annotate -x FORMAT,INFO > ${VCF_PREFIX}_norm.vcf
    picard UpdateVcfSequenceDictionary I=${VCF_PREFIX}_norm.vcf O=${VCF_PREFIX}_norm_header.vcf SD=hs37d5.dict
    picard SortVcf I=${VCF_PREFIX}_norm_header.vcf O=${VCF_PREFIX}_norm_sort.vcf SD=hs37d5.dict
    gzip ${VCF_PREFIX}_norm_sort.vcf
done

rm *_norm.vcf
rm *_norm_header.vcf
rm *_norm_sort.vcf.idx

# Filter by size
for VCF_FILE in HG002_GRCh37*_sort.vcf.gz; do
    VCF_FILENAME=$(basename ${VCF_FILE})
    VCF_PREFIX=${VCF_FILENAME%.vcf.gz}
    echo ${VCF_PREFIX}
    filterStructuralVariants ${VCF_PREFIX}.vcf.gz ${VCF_PREFIX}_sv19 19 -19 19
done
```

### Combining variants

SNVs/indels were combined with the SV catalog.

```
bayesTyperTools combine -v sv:AJ_GRCh37_GIAB_sv_0.5_norm_sort.vcf.gz,all19:HG002_GRCh37_CHROM1-MT_novoalign_Ilmn250x250_FB_xy_norm_sort_sv19.vcf.gz,all19:HG002_GRCh37_CHROM1-MT_novoalign_Ilmn250x250_GATKHC_xy_norm_sort_sv19.vcf.gz,all19:HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_all_norm_sort_sv19.vcf.gz -o giab_AJ_sv_HG002_all19 -z
```

### Genotyping

```
OUT_PREFIX=bayestyper_unit_1/giab_AJ_sv_HG002_all19

bayesTyper cluster -v giab_AJ_sv_HG002_all19.vcf.gz -s sample.txt -g hs37d5_canon.fa -d hs37d5_decoy.fa -p 24 -r 12345678
bayesTyper genotype -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_cluster_data -s sample.txt -g hs37d5_canon.fa -d hs37d5_decoy.fa -o ${OUT_PREFIX} -z -p 24 -r 12345678 --min-genotype-posterior 0

bcftools filter -i 'FILTER=="PASS"' ${OUT_PREFIX}.vcf.gz | gzip -c 1> ${OUT_PREFIX}_pass.vcf.gz
filterAlleleCallsetOrigin ${OUT_PREFIX}_pass.vcf.gz ${OUT_PREFIX}_pass_hgsvc all19,. 1
```

`samples.txt` contains paths to the kmers counts:

```
HG002	F	HG002_NA24385_50x_k55
```
