#### `VCFtoSymbolicDEL.py`

Reads through a VCF file in explicit format, identifies deletions and produces a symbolic output. 
Multiple ALT sequences are split into multiple symbolic records. 
Common prefixes/suffixes in the REF and ALT sequences are removed.

#### Scripts from the BayesTyper tool

Some helper scripts provided by BayesTyper were used to manipulate the VCFs.
In particular:

- `filterStructuralVariants` to filter SNVs/indels/SVs by size.
- `bayesTyperTools combine` to combine multiple VCF files.
- `bayesTyperTools convertAllele` to convert a symnolic VCF into a VCF with explicit alleles.
