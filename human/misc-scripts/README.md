#### `VCFtoSymbolicDEL.py`

Reads through a VCF file in explicit format, identifies deletions and produces a symbolic output. 
Multiple ALT sequences are split into multiple symbolic records. 
Common prefixes/suffixes in the REF and ALT sequences are removed.
