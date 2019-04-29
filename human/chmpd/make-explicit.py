#!/usr/bin/env python2.7

"""
Assume:
POS coordinates are correct
REF entries are wrong
SEQ elements are correct (but don't contain match reference base at beginning) but shifted for deletions

-> pull REF from Fasta
-> deleteion: REF = REF + SEQ, ALT = REF
-> insertion: REF = REF, ALT = REF + SEQ

don't touch inversions. 

throw error if more than one alt
(hacked together from vcf-add-bed-seqs.py)
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools, math
import vcf, collections, gzip, re
import pysam
from Bio.Seq import Seq

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("vcf", type=str,
                        help="VCF whose SV sequences we want to fill out")
    parser.add_argument("--fasta",
                        default='ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/'
                        'GRCh38_full_analysis_set_plus_decoy_hla.fa',
                        help="Fasta file for reference.  Needs .fai as well")
                        
    args = args[1:]
    options = parser.parse_args(args)
    return options


def open_input(file_path):
    open_fn = gzip.open if file_path.endswith('.gz') else open
    return open_fn(file_path, 'r')

def main(args):
    options = parse_args(args)


    # fasta index needed for reference bases
    faidx = pysam.FastaFile(options.fasta)

    # print the edited vcf
    with open_input(options.vcf) as vcf_file:
        add_span_info = True
        for line in vcf_file:
            if line.startswith('#'):
                sys.stdout.write(line)
            elif line:
                vcf_toks = line.split('\t')
                vcf_chrom = vcf_toks[0]
                vcf_pos = int(vcf_toks[1])
                vcf_sv_type = vcf_toks[4][1:-1]

                sv_seq = [i for i in vcf_toks[7].split(';') if i.startswith('SEQ=')][0][4:]
                ref_seq = faidx.fetch(vcf_chrom, vcf_pos - 1, vcf_pos)

                if vcf_sv_type == 'INS':
                    vcf_toks[3] = ref_seq;
                    vcf_toks[4] = ref_seq + sv_seq
                elif vcf_sv_type == 'DEL':
                    vcf_toks[3] = ref_seq + sv_seq
                    # it looks like the vcf_seqs are shifted.  we assume POS is gospel an reload from Fasta
                    ref_del_seq = faidx.fetch(vcf_chrom, vcf_pos - 1, vcf_pos - 1 + len(vcf_toks[3]))
                    vcf_toks[3] = ref_del_seq
                    vcf_toks[4] = ref_seq

                sys.stdout.write('\t'.join(vcf_toks))

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
