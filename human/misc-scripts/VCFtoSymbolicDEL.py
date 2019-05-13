"""
Reads through a VCF file in explicit format, identifies deletions and produces a 
symbolic output. Multiple ALT sequences are split into multiple symbolic records.
Common prefixes/suffixes in the REF and ALT sequences are removed.
"""

import fileinput


def findSuffix(ref, alt):
    suff = ''
    cpt = 1
    while(cpt < len(ref) and cpt < len(alt) and ref[-cpt] == alt[-cpt]):
        suff = ref[-cpt] + suff
        cpt += 1
    return suff


def findPrefix(ref, alt):
    pref = ''
    cpt = 0
    while(cpt < len(ref) and cpt < len(alt) and ref[cpt] == alt[cpt]):
        pref = pref + ref[-cpt]
        cpt += 1
    return pref


in_headers = True
for line in fileinput.input():
    if line[:2] == '##':
        line = line.rstrip()
        print line
    elif in_headers:
        new_heads = '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n'
        new_heads += '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'
        new_heads += '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
        new_heads += '##INFO=<ID=REFSEQ,Number=1,Type=String,Description="The original REF sequence">\n'
        new_heads += '##INFO=<ID=ALTSEQ,Number=1,Type=String,Description="The original ALT sequence">\n'
        new_heads += '##ALT=<ID=DEL,Description="Deletion">'
        print new_heads
        in_headers = False
        # skip format and genotypes
        line = line.rstrip().split('\t')
        line = line[:8]
        print '\t'.join(line)
    else:
        line = line.rstrip().split('\t')
        ref = line[3]
        alts = line[4].split(',')
        for alt in alts:
            suff = findSuffix(ref, alt)
            pref = findPrefix(ref, alt)
            ref_short = ref[len(pref):(len(ref)-len(suff))]
            ref_padd = ref[len(pref)-1]
            alt = alt[len(pref):(len(alt)-len(suff))]
            if len(ref_short) <= len(alt):
                # Skip insertion and SNVs
                continue
            else:
                linealt = list(line)
                # Deletion
                # shift position by prefix length
                linealt[1] = str(int(linealt[1]) + len(pref) - 1)
                # update info
                del_size = len(ref_short)
                end_pos = int(linealt[1]) + del_size
                linealt[7] += ';SVLEN=-{};SVTYPE=DEL;CIPOS=-100,100;CIEND=-100,100;END={};REFSEQ={};ALTSEQ={}'.format(del_size, end_pos, ref_padd + ref_short, ref_padd + alt)
                # update ALT
                linealt[3] = ref_padd
                linealt[4] = '<DEL>'
                # skip format and genotypes
                linealt = linealt[:8]
                print '\t'.join(linealt)
