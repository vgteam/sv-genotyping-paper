#!/usr/bin/env python2.7

"""
We need a VCF for evaluation.  This script takes the pseudodiploid VCF and genotype table and uses it to make a 
phased pseudo diploid VCF where CHM1 is haplotype 0 and CHM13 is haplotype 1
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools, math
import vcf, collections, gzip, re


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("vcf", type=str,
                        help="VCF from which to extract samples")
    parser.add_argument("tab", type=str,
                        help="table to look up genotypes")
    parser.add_argument("--no-hom-ref", action="store_true",
                        help="don't include 0|0 calls")
    parser.add_argument("--ref-only", action="store_true",
                        help="only fix REF column")
                        
    args = args[1:]
    options = parser.parse_args(args)
    return options

def tab_header(line):
    """ dict mapping column title to column number
    """
    return dict([(col_name, col_number) for col_number, col_name in enumerate(line.strip().split('\t'))])

def open_input(file_path):
    open_fn = gzip.open if file_path.endswith('.gz') else open
    return open_fn(file_path, 'r')

def main(args):
    options = parse_args(args)

    # read the table into memory
    tab_map = {}
    with open_input(options.tab) as tab_file:
        header_map = tab_header(tab_file.readline())
        for line in tab_file:
            toks = line.strip().split('\t')
            if toks and not toks[0].startswith('#'):
                tab_map[toks[header_map['ID']]] = toks

    # print the exctracted vcf
    with open_input(options.vcf) as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                line = line.replace('##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">',
                                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
                sys.stdout.write(line)
            elif line:
                vcf_toks = line.split('\t')
                tab_toks = tab_map[vcf_toks[2]]

                if not options.ref_only:
                    if tab_toks[header_map['SAMPLE']] in ['CHM1', 'CHM13']:
                        if tab_toks[header_map['CALL']] == 'HOM_ALT':
                            gt = '1|1'
                        elif tab_toks[header_map['CALL']] == 'HET':
                            gt = '1|0' if tab_toks[header_map['SAMPLE']] == 'CHM1' else '0|1'
                        else:
                            assert False
                    elif not options.no_hom_ref:
                        gt = '0|0'
                    else:
                        gt = None

                    if gt:
                        vcf_toks[-1] = gt + vcf_toks[-1][1:]

                # we add one in an effort to make the reference base match the fasta
                # vg doesn't care for symbolic alleles, but the vcf comparators do
                vcf_toks[1] = str(int(vcf_toks[1]) + 1)

                sys.stdout.write('\t'.join(vcf_toks))

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
