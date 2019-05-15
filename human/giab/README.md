# GIAB

## Data

- [Hg19 reference genome](ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz).
- [GIAB v0.5.0 VCF](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_UnionSVs_12122017/svanalyzer_union_171212_v0.5.0_annotated.vcf.gz).
- [High-confidence regions BED file](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed).

## Reads

- Original reads mapped to hg19: []()

The reads were downsampled to 50x using to create `HG002-NA24385-50x.bam`:

```
???
```

## toil-vg

Using the helper scripts from `../toil-scripts`.

```
# construct a 0.50 graph with HG002 control
./construct-hgsvc.sh  -G -c ${CLUSTER}1 ${JOBSTORE}1x ${OUTSTORE}/GIAB-0.5-FEB26

# map our downsampled 50X reads from the paper
./mce-hgsvc.sh -p -E " -p -n -m 20" -C " -p" -c ${CLUSTER}1  ${JOBSTORE}1 ${OUTSTORE}/GIAB-0.5-FEB26 s3://${OUTSTORE}/GIAB-0.5-FEB26/GIAB HG002 HG002 s3://${OUTSTORE}/GIAB-0.5-FEB26/giab-0.5.vcf.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed s3://glennhickey/mike-hs37d5-hg002-nov18/bams/HG002-NA24385-50x.bam

aws s3 sync s3://${OUTSTORE}/GIAB-FEB26/eval-HG002 ./GIAB-FEB26-eval-HG002
aws s3 sync s3://${OUTSTORE}/GIAB-0.5-FEB26/eval-HG002 ./GIAB-0.5-FEB26-eval-HG002
```


## Delly

```
VCF=svanalyzer_union_171212_v0.5.0_annotated.vcf.gz
BAM=HG002-NA24385-50x.bam
REF=hs37d5.fa

delly call -g $REF -v $VCF -o temp.bcf $BAM
bcftools convert temp.bcf > HG002.delly.vcf
bgzip HG002.delly.vcf
tabix HG002.delly.vcf.gz
```

## SVTyper

The explicit VCF was first converted into a symbolic VCF, for deletions only because SVTyper cannot genotype insertions.
See *VCFtoSymbolicDEL.py* in [../misc-scripts](../misc-scripts).

```
zcat svanalyzer_union_171212_v0.5.0_annotated.vcf.gz | python VCFtoSymbolicDEL.py > giab-0.5.symbolic.dels.vcf
```

Then:

```
BAM=HG002-NA24385-50x.bam

svtyper -i HGSVC.haps.symbolic.dels.vcf -B $BAM -l libinfo.json > HG002.svtyper.vcf
vcf-sort HG002.svtyper.vcf | vcfkeepinfo - END SVLEN SVTYPE NS AN | grep -v 'contig' > HG002.svtyper.clean.vcf
bayesTyperTools convertAllele -v HG002.svtyper.clean.vcf -g ../../haps/hg38.fa --keep-imprecise -o HG002.svtyper.explicit
vcfkeepinfo HG002.svtyper.explicit.vcf NS AN | bgzip > HG002.svtyper.explicit.vcf.gz
tabix -f HG002.svtyper.explicit.vcf.gz
```

## BayesTyper

### kmer counting

```
BAM=HG002_NA24385_50x.bam
kmc -v -k55 -ci1 -fbam -t12 $BAM HG002_NA24385_50x_k55 .
bayesTyperTools makeBloom -p 12 -k HG002_NA24385_50x_k55
```

### SNVs/indels 

SNVs downloaded from GIAB releases, normalized and filtered to keep only variants smaller than 50 bp.

```
GIAB_FTP="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh37/supplementaryFiles"
wget ${GIAB_FTP}/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_all.vcf.gz
bcftools view -r X,Y -O z ${GIAB_FTP}/inputvcfsandbeds/HG002_GRCh37_CHROM1-MT_novoalign_Ilmn250x250_FB.vcf.gz > HG002_GRCh37_CHROM1-MT_novoalign_Ilmn250x250_FB_xy.vcf.gz
bcftools view -r X,Y -O z ${GIAB_FTP}/inputvcfsandbeds/HG002_GRCh37_CHROM1-MT_novoalign_Ilmn250x250_GATKHC.vcf.gz > HG002_GRCh37_CHROM1-MT_novoalign_Ilmn250x250_GATKHC_xy.vcf.gz
cp svanalyzer_union_171212_v0.5.0_annotated.vcf.gz giab-0.5.vcf.gz AJ_GRCh37_GIAB_sv_0.5.vcf.gz
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
~/tools/filterAlleleCallsetOrigin ${OUT_PREFIX}_pass.vcf.gz ${OUT_PREFIX}_pass_hgsvc all19,. 1
```

`samples.txt` contains paths to the kmers counts:

```
HG002	F	HG002_NA24385_50x_k55
```
