from pyfaidx import Fasta
import fileinput


# Load reference genome
fa = Fasta('hs37d5.fa')

# stream VCF
for line in fileinput.input():
    # if header, print and got to next line
    if line[0] == '#':
        line = line.rstrip()
        print line
        continue
    # else parse variant record
    line = line.rstrip().split('\t')
    chrom = line[0]
    pos = int(line[1])
    ref = line[3]
    alt = line[4]
    if ref[0] == alt[0]:
        # if variant and padding present, print
        print '\t'.join(line)
    else:
        # else missing padding
        # find padding base
        nuc = fa[str(chrom)][pos - 2]
        line[3] = nuc.seq + line[3]
        line[4] = nuc.seq + line[4]
        # shift position
        line[1] = str(pos - 1)
        # print padded record
        print '\t'.join(line)
