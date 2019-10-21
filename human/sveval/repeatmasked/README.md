The Snakemake pipeline runs the python script that annotate a VCF file with RepeatMasker ([`rmskAnnotateVcf.py`](rmskAnnotateVcf.py)).
It was run inside a Docker container built from the [Dockerfile](Dockerfile).

```
snakemake -p all
```

Three new fields are added to the VCFs: 

- *RMSKNAME* with the repeat element's name.
- *RMSKCLASS* with the repeat's class and family.
- *RMSKCOV* with the proportion of the variant covered by the repeat.


The [`eval.R`](eval.R) script retrieves variants covered at at 80% by a repeat element and  produces a *eval-rmsk-hgsvc-vg-HG00514-call-geno.tsv* TSV file with the precision/recall/F1 for the most frequent repeat classes.
This file is used in to make the manuscript's figure (see [manuscript repo](https://github.com/jmonlong/manu-vgsv/tree/master/figures)).
