```
# make the CHM-PSEUDODIPLOID graph
./construct-hgsvc.sh -s -c ${CLUSTER}2 ${JOBSTORE}2 ${OUTSTORE}/CHMPD-feb12

# map the 30x reads

./mce-hgsvc.sh -c ${CLUSTER}3 ${JOBSTORE}3 ${OUTSTORE}/CHMPD-feb12 s3://${OUTSTORE}/CHMPD-feb12/CHMPD PSEUDOSET PSEUDOSET-30 s3://${OUTSTORE}/CHMPD-feb12/pseudo_diploid-explicit.vcf.gz ${COMPARE_REGIONS_BED} s3://majorsv-ucsc/gt/gam/aln_30x.gam
#done

aws s3 sync s3://${OUTSTORE}/CHMPD-feb12/eval-PSEUDOSET-30 ./CHMPD-feb12-eval-PSEUDOSET


# then do the SMRTSV that we copied from courtyard
# note that we made it explicit with
# ./make-explicit.py PSEUDOSET-smrtsv.vcf.gz --fasta ~/dev/work/references/hg38.fa.gz  | vcfkeepinfo - NA | vcffixup - | bgzip > PSEUDOSET-smrtsv-explicit.vcf.gz

./eval-hgsvc.sh -c ${CLUSTER}1 ${JOBSTORE}1x ${OUTSTORE}/CHMPD-feb12/eval-PSEUDOSET-30-smrtsv s3://${OUTSTORE}/CHMPD-feb12/pseudo_diploid-explicit.vcf.gz  s3://glennhickey/outstore/CHMPD-feb12/call-PSEUDOSET-30-smrtsv/PSEUDOSET-smrtsv.explicit.vcf.gz ${COMPARE_REGIONS_BED} PSEUDOSET

aws s3 sync s3://${OUTSTORE}/CHMPD-feb12/eval-PSEUDOSET-30-smrtsv ./CHMPD-feb12-eval-PSEUDOSET-smrtsv



### move to manuscript sv
pushd CHMPD-feb12-eval-PSEUDOSET-smrtsv
for i in sveval* ; do cd $i; tar zxf sv_evaluation.tar.gz; cd ..; done
cp sveval-norm/sv_evaluation/prcurve.tsv ~/Documents/Research/manu-vgsv/figures/data/chmpd-smrtsv2-prcurve.tsv 
cp sveval-clip-norm/sv_evaluation/prcurve.tsv ~/Documents/Research/manu-vgsv/figures/data/chmpd-smrtsv2-clip-prcurve.tsv 
cp sveval-clip-norm-gt/sv_evaluation/prcurve.tsv ~/Documents/Research/manu-vgsv/figures/data/chmpd-smrtsv2-clip-prcurve-geno.tsv
cp sveval-norm-gt/sv_evaluation/prcurve.tsv ~/Documents/Research/manu-vgsv/figures/data/chmpd-smrtsv2-prcurve-geno.tsv
popd

pushd CHMPD-feb12-eval-PSEUDOSET
for i in sveval* ; do cd $i; tar zxf sv_evaluation.tar.gz; cd ..; done
cp sveval-norm/sv_evaluation/prcurve.tsv ~/Documents/Research/manu-vgsv/figures/data/chmpd-construct-prcurve.tsv 
cp sveval-clip-norm/sv_evaluation/prcurve.tsv ~/Documents/Research/manu-vgsv/figures/data/chmpd-construct-clip-prcurve.tsv 
cp sveval-clip-norm-gt/sv_evaluation/prcurve.tsv ~/Documents/Research/manu-vgsv/figures/data/chmpd-construct-clip-prcurve-geno.tsv
cp sveval-norm-gt/sv_evaluation/prcurve.tsv ~/Documents/Research/manu-vgsv/figures/data/chmpd-construct-prcurve-geno.tsv
popd
```
