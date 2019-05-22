#!/bin/bash

## Constructs graph with vg using the true SVs and the "calls" SVs.

#### Primary
vg construct -r ref.fa -m 32 > primary.vg
vg index -x primary.xg -g primary.gcsa primary.vg

#### Thruth
## Prepare truth.vcf
vcf-sort truth.vcf | bgzip -c > truth.vcf.gz
tabix -f truth.vcf.gz
## trugh graph
vg construct -r ref.fa -v truth.vcf.gz -a -f -m 32 > truth.vg
vg index -x truth.xg -G truth.gbwt -v truth.vcf.gz -g truth.gcsa -F truth.threads truth.vg

#### Calls
## Prepare inexact vcf
vcf-sort calls.vcf | bgzip -c > calls.vcf.gz
tabix -f calls.vcf.gz
## calls graph
vg construct -r ref.fa -v calls.vcf.gz -m 32 -a -f > calls.vg
vg index -x calls.xg -G calls.gbwt -v calls.vcf.gz -g calls.gcsa -F calls.threads calls.vg


#### Thruth with symbolic VCF
## Prepare truth.vcf
vcf-sort truth.symb.vcf | bgzip -c > truth.symb.vcf.gz
tabix -f truth.symb.vcf.gz
## trugh graph
vg construct -r ref.fa -v truth.symb.vcf.gz -a -f -m 32 -S > truth.symb.vg
vg index -x truth.symb.xg -G truth.symb.gbwt -v truth.symb.vcf.gz -g truth.symb.gcsa -F truth.symb.threads truth.symb.vg

#### Calls with symbolic VCF
## Prepare inexact vcf
vcf-sort calls.symb.vcf | bgzip -c > calls.symb.vcf.gz
tabix -f calls.symb.vcf.gz
## calls graph
vg construct -r ref.fa -v calls.symb.vcf.gz -m 32 -a -f -S > calls.symb.vg
vg index -x calls.symb.xg -G calls.symb.gbwt -v calls.symb.vcf.gz -g calls.symb.gcsa -F calls.symb.threads calls.symb.vg
