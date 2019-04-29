#!/usr/bin/env python2.7

"""
We construct HGSVC VCFs from pairs of haplotype VCFs.  There's no off the shelf
tool to merge them together that I'm aware of, so we do it here.
Example: 
File 1: chr1 100 A TT  .
and
File 2: chr1 100 A TTT .
goes to 
chr1 100 A TT/TTT 1/2
in the merged file.  
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools, math
import vcf, collections

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("hap1_vcf", type=str,
                        help="Input vcf file for haplotype 1")
    parser.add_argument("hap2_vcf", type=str,
                        help="Input vcf file for haplotype 2")
                        
    args = args[1:]
    options = parser.parse_args(args)
    return options


def merge_records(record1, record2):
    """ merge two haplotype records into a diploid record
    """
    # normalize to have same reference
    if len(record1.REF) > len(record2.REF):
        record2.REF += record1.REF[len(record2.REF):]
        record2.ALT += record1.REF[len(record2.REF):]
    elif len(record1.REF) < len(record2.REF):
        record1.REF += record2.REF[len(record1.REF):]
        record1.ALT += record2.REF[len(record1.REF):]
    
    # merge into record
    assert record1.REF == record2.REF
    if record2.ALT != record1.ALT:
        record1.ALT += record2.ALT
        biallele = True
    else:
        biallele = False

    # convert from namedtuple to dict so we can modify
    #data = record1.samples[0].data
    #data_dict = data._asdict()
    data_dict = {}
    data_dict['GT'] = "1|2" if biallele else "1|1"
    nt = collections.namedtuple('CallData', ' '.join(data_dict.keys()))(**data_dict)
    record1.samples[0].data = nt

    return record1

def add_genotype(record, g):
    """ make a het genotype for unmerged records """
    data_dict = {'GT': g}
    nt = collections.namedtuple('CallData', ' '.join(data_dict.keys()))(**data_dict)
    record.samples[0].data = nt
    return record

def main(args):
    options = parse_args(args)

    # load one file into memory to simplify iteration
    with open(options.hap2_vcf) as hap2_file:
        records2 = [record for record in vcf.Reader(hap2_file)]
    
    with open(options.hap1_vcf) as hap1_file:
        reader1 = vcf.Reader(hap1_file)
        vcf_writer = vcf.Writer(sys.stdout, reader1)
        i = 0

        for record1 in reader1:
            # catch up on second stream
            while i < len(records2) and (records2[i].CHROM, records2[i].POS) < (record1.CHROM, record1.POS):
                add_genotype(records2[i], '0|1')
                vcf_writer.write_record(records2[i])
                i += 1

            # output merged record if variant at same position in both streams
            if i < len(records2) and (records2[i].CHROM, records2[i].POS) == (record1.CHROM, record1.POS):
                record = merge_records(record1, records2[i])
                i += 1
            else:
                record = add_genotype(record1, '1|0')
            vcf_writer.write_record(record)            

    # print everything left in second stream
    while i < len(records2):
        add_genotype(records2[i], '0|1')
        vcf_writer.write_record(records2[i])
        i += 1

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
