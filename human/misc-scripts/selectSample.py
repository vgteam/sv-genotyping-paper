import argparse
import gzip

parser = argparse.ArgumentParser(description='Filter/edit sample name in VCF.')
parser.add_argument('-v', help='the VCF file', required=True)
parser.add_argument('-s', help='the sample', required=True)
args = parser.parse_args()

sampcol = -1
vcff = gzip.open(args.v, 'r')
for line in vcff:
    line = line.rstrip()
    # write headers
    if line[:2] == '##':
        if '##SAMPLE=' in line:
            continue
        print line
        continue
    line = line.split('\t')
    if '#CHROM'in line[0]:
        print '##SAMPLE=<ID={}>'.format(args.s)
        for ii in range(len(line)):
            if line[ii] == args.s:
                sampcol = ii
        if sampcol == -1:
            sampcol = 9
        line[sampcol] = args.s
    print '\t'.join(line[:8]) + '\t' + line[8] + '\t' + line[sampcol]
vcff.close()
