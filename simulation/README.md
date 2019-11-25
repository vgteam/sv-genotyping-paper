## Simulation of a region and SVs

A buffer of sequence is simulated between each SV (500 bp) to avoid overlap.
SVs are insertions, deletions and inversions following the size distribution of the [HGSVC dataset](../human/hgsvc).
Genotypes across three samples are simulated, (hom ref, het, hom alt randomly assigned for each variant).

To mimick errors during SV calling the python script output a *truth.vcf* with the original SVs (used for example to simulate reads), and a *calls.vcf* with slightly different breakpoint location (or sequence for insertions).
It's called *calls.vcf* because it's what a SV caller might generate.
Two types of errors are simulated: 

- **shift** the coordinates or deletion or insertions by 1-10 bp
- **del**ete 1-10 bp from the start or end of the inserted sequence.

The simulated reference sequence and reads that we used in the manuscript were deposited at 
[`https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=vgsv2019/simulation`](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=vgsv2019/simulation).


### 1. Prepare the environment

- Create folders.
- Download real reads to provide an error model for the read simulation.
- Activate the toil-vg environment.

```
mkdir -p logs reports
curl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_005_BH814YADXX/Project_RM8398/Sample_U5a/U5a_AGTCAA_L002_R1_007.fastq.gz > U5a_AGTCAA_L002_R1_007.fastq.gz

source ~/build/toil-vg/toilvenv/bin/activate

## Eventually clean up previous run
rm -rf *.vg *.xg *.threads *.gcsa* *.gbwt *.fai cov* s*.fastq* *.bam *.gam *.gam *.index *.vcf* ref.fa* logs/*
```

In several of the following scripts there is a *THREADS* variable that can be changed to adjust the number of cores used by the tools.
It is now set to `THREADS=2`.

### 2. Simulate sequence and SVs

- `HGSVC.haps.vcf.gz` is the VCF with the HGSVC SV catalog created [here](../human/hgsvc).
- `vcf-sim-errors-buffer.py` simulates the sequence and SVs, creating fasta and vcf files.
- `make-graphs.sh` constructs graph with vg using the true SVs and the "calls" SVs.

```
python vcf-sim-errors-buffer.py -n 1000 -r HGSVC.haps.vcf.gz
sh make-graphs.sh &> logs/make-graphs.log
```

### 3. Simulate reads

- `sim-reads.sh` simulates reads for each sample using the graph containing the true SVs.
- Its parameter is the desired depth: 1x, 3x, 7x, 13x, 20x.
- This script creates folders staring with `cov` for each depth/coverage, that will contains reads and genotype calls.

```
for COV in 1 3 7 13 20
do
	sh sim-reads.sh $COV &> logs/sim-reads-$COV.log
done
```

### 4. Genotype SVs with vg

- `genotype-vg.sh` genotypes each sample using toil-vg.
- Like in the previous section, the parameter is the depth/coverage of the experiment. 

```
for COV in 1 3 7 13 20
do
	sh genotype-vg.sh $COV &> logs/genotype-vg-$COV.log
done
```

### 5. Genotype SVs with other methods

- `prepare-other-methods.sh` build BWA indexes and reformat VCFs.
- `genotype-other-methods.sh` genotypes each sample with Delly, SVTyper, and BayesTyper.
- Like in the previous sections, the parameter is the depth/coverage of the experiment. 

```
sh prepare-other-methods.sh &> logs/prepare-other-methods.log

for COV in 1 3 7 13 20
do
	sh genotype-other-methods.sh $COV &> logs/genotype-other-methods-$COV.log
done
```

### 6. Evaluation

- `eval.R` reads all the VCFs and evaluate the genotyping accuracy using the [sveval package](https://github.com/jmonlong/sveval).
- It produces *simerror-prcurve.tsv* and *simerror-geno-prcurve.tsv*.
- The figures in the manuscript are generated using these files. Find the code on the [manuscript's repo](https://github.com/jmonlong/manu-vgsv/tree/master/figures)

```
Rscript eval.R &> logs/eval.log
```
