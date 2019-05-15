#### Copy vcfs
mkdir -p vcf

source ~/build/toil-vg/toilvenv/bin/activate

## One method/sample pair
snakemake -np single-hgsvc-bayestypersv-HG00514


#### HGSVC
## vg
aws s3 cp s3://glennhickey/outstore/HGSVC-jan5/call-HG00514-sim/HG00514.vcf.gz vcf/hgsvcsim-vg-HG00514.vcf.gz
aws s3 cp s3://glennhickey/outstore/HGSVC-jan5/call-HG00514/HG00514.vcf.gz vcf/hgsvc-vg-HG00514.vcf.gz
aws s3 cp s3://glennhickey/outstore/HGSVC-jan5/call-HG00733/HG00733.vcf.gz vcf/hgsvc-vg-HG00733.vcf.gz
aws s3 cp s3://glennhickey/outstore/HGSVC-jan5/call-NA19240/NA19240.vcf.gz vcf/hgsvc-vg-NA19240.vcf.gz

## Delly
cp ../delly/HG00514_sim30x/HG00514_sim30x.delly.vcf.gz vcf/hgsvcsim-delly-HG00514.vcf.gz
cp ../delly/ERR903030/ERR903030.delly.vcf.gz vcf/hgsvc-delly-HG00514.vcf.gz
cp ../delly/HG00733/HG00733.delly.vcf.gz vcf/hgsvc-delly-HG00733.vcf.gz
cp ../delly/NA19240/NA19240.delly.vcf.gz vcf/hgsvc-delly-NA19240.vcf.gz

## svtyper
cp ../svtyper/HG00514_sim30x/HG00514_sim30x.svtyper.explicit.vcf.gz vcf/hgsvcsim-svtyper-HG00514.vcf.gz
cp ../svtyper/HG00514/HG00514.svtyper.explicit.vcf.gz vcf/hgsvc-svtyper-HG00514.vcf.gz
cp ../svtyper/HG00733/HG00733.svtyper.explicit.vcf.gz vcf/hgsvc-svtyper-HG00733.vcf.gz
cp ../svtyper/NA19240/NA19240.svtyper.explicit.vcf.gz vcf/hgsvc-svtyper-NA19240.vcf.gz

## bayestyper
cp ../bayestyper/HG00514_sim30x_hgsvc_pp_bayestyper_pass_nomis.vcf.gz vcf/hgsvcsim-bayestyper-HG00514.vcf.gz
zcat ../../../jsibbese/hgsvc/real/all/bayestyper_v1.5beta/hgsvc_pp49_hc49/bayestyper_unit_1/all_hgsvc_pp49_hc49_bayestyper_pass_hgsvc.vcf.gz | bgzip > vcf/hgsvc-bayestyper-HG00514.vcf.gz
ln -s hgsvc-bayestyper-HG00514.vcf.gz vcf/hgsvc-bayestyper-HG00733.vcf.gz
ln -s hgsvc-bayestyper-HG00514.vcf.gz vcf/hgsvc-bayestyper-NA19240.vcf.gz
zcat ../../../jsibbese/hgsvc/real/all/bayestyper_v1.5beta/hgsvc/bayestyper_unit_1/all_hgsvc_bayestyper_pass_nomis.vcf.gz | bgzip > vcf/hgsvc-bayestypersv-HG00514.vcf.gz

snakemake -np hgsvc


#### GiaB
## vg
aws s3 cp s3://glennhickey/outstore/GIAB-0.5-FEB26/call-HG002/HG002.vcf.gz vcf/giab5-vg-HG002.vcf.gz

## Delly
cp ../delly/HG002/HG002-giab-0.5.delly.vcf.gz vcf/giab5-delly-HG002.vcf.gz

## svtyper
cp ../svtyper/HG002/HG002.svtyper.explicit.vcf.gz vcf/giab5-svtyper-HG002.vcf.gz

## bayestyper
zcat /public/groups/cgl/graph-genomes/jsibbese/giab/bayestyper_v1.5beta/giab_AJ_sv_HG002_all19/bayestyper_unit_1/giab_AJ_sv_HG002_all19_pass_giab.vcf.gz | bgzip > vcf/giab5-bayestyper-HG002.vcf.gz

snakemake -np giab5


#### SVPOP
## vg
aws s3 cp s3://glennhickey/outstore/SVPOP-jan10/call-HG00514/HG00514.vcf.gz vcf/svpop-vg-HG00514.vcf.gz
aws s3 cp s3://glennhickey/outstore/SVPOP-jan10/call-HG00733/HG00733.vcf.gz vcf/svpop-vg-HG00733.vcf.gz
aws s3 cp s3://glennhickey/outstore/SVPOP-jan10/call-NA19240/NA19240.vcf.gz vcf/svpop-vg-NA19240.vcf.gz

## SMRT-SV (fix header to ensure explicit mode)
aws s3 cp s3://glennhickey/outstore/SVPOP-jan10/call-all-smrtsv/svpop-smrtsv-hgsvc-explicit.vcf.gz temp.vcf.gz
zgrep -v "##INFO=<ID=SVTYPE" temp.vcf.gz | grep -v "##INFO=<ID=END" | bgzip > vcf/svpop-smrtsv-HG00514.vcf.gz
rm temp.vcf.gz
ln -s svpop-smrtsv-HG00514.vcf.gz vcf/svpop-smrtsv-HG00733.vcf.gz
ln -s svpop-smrtsv-HG00514.vcf.gz vcf/svpop-smrtsv-NA19240.vcf.gz

snakemake -np svpop

## No-call experiment
aws s3 cp s3://glennhickey/outstore/SVPOP-jan10/call-HG00514/HG00514.vcf.gz vcf/svpop-vgref-HG00514.vcf.gz

snakemake vcf/svpop-vg-HG00514-norm.vcf.gz
snakemake vcf/svpop-smrtsv-HG00514-norm.vcf.gz
snakemake ../glenndata/sv-pop-explicit-norm.vcf.gz

Rscript findNocallsRegions.R ../glenndata/sv-pop-explicit-norm.vcf.gz HG00514 nocalls-HG00514-vg-smrtsv.bed vcf/svpop-vg-HG00514-norm.vcf.gz vcf/svpop-smrtsv-HG00514-norm.vcf.gz
Rscript invertRegions.R nocalls-HG00514-vg-smrtsv.bed called-HG00514-vg-smrtsv.bed

snakemake singlenocalls-svpop-vg-HG00514
snakemake singlenocalls-svpop-smrtsv-HG00514

#### CHM-PD
## vg
aws s3 cp s3://glennhickey/outstore/CHMPD-feb12/call-PSEUDOSET-30/PSEUDOSET.vcf.gz vcf/chmpd-vg-chmpd.vcf.gz

## SMRT-SV
aws s3 cp s3://glennhickey/outstore/CHMPD-feb12/call-PSEUDOSET-30-smrtsv/PSEUDOSET-smrtsv.explicit.vcf.gz vcf/chmpd-smrtsv-chmpd.vcf.gz

snakemake -np chmpd


