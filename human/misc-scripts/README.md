#### `VCFtoSymbolicDEL.py`

Reads through a VCF file in explicit format, identifies deletions and produces a symbolic output. 
Multiple ALT sequences are split into multiple symbolic records. 
Common prefixes/suffixes in the REF and ALT sequences are removed.

#### `addMissingPaddingHg37.py`

Simple script to add padding sequence when needed. 
Used to avoid Paragraph throwing errors on the GIAB VCF file.
The path to the reference path is hard-coded.

#### `addMissingPaddingHg38.py`

Simple script to add padding sequence when needed. 
Used to avoid Paragraph throwing errors on the SVPOP VCF file.
The path to the reference path is hard-coded.

#### Scripts from the BayesTyper tool

Some helper scripts provided by [BayesTyper](https://github.com/bioinformatics-centre/BayesTyper) were used to manipulate the VCFs.
In particular:

- `filterStructuralVariants` to filter SNVs/indels/SVs by size.
- `filterAlleleCallsetOrigin` to filter alleles based on their origin.
- `bayesTyperTools combine` to combine multiple VCF files.
- `bayesTyperTools convertAllele` to convert a symnolic VCF into a VCF with explicit alleles.
