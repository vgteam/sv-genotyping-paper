import random as rd
import numpy as np
import argparse
import vcf
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

nuc = np.array(["A", "T", "C", "G"])


class Insertion:
    def __init__(self, region, ins_size, gt=[], label=''):
        self.region = region
        self.label = label
        self.pos = len(self.region.seq)
        seqArray = nuc[[int(rd.random()*4) for i in xrange(ins_size)]]
        self.seq = MutableSeq("".join(seqArray), generic_dna)
        self.gt = list(gt)
        self.region.addBuffer()

    def shiftError(self, bp):
        self.label += 'sE{}'.format(bp)
        self.pos += bp

    def delError(self, bp, pos='end'):
        self.label += 'dE{}{}'.format(bp, pos)
        if(pos == 'start'):
            self.seq = self.seq[bp:]
        if(pos == 'end'):
            self.seq = self.seq[:(len(self.seq)-bp)]

    def size(self):
        return len(self.seq)

    def nbSamp(self):
        return len(self.gt)

    def toVcf(self):
        if(self.label == ''):
            var_id = 'ins' + str(self.pos) + '' + str(len(self.seq))
        else:
            var_id = self.label
        vcf = self.region.chrom + '\t' + str(self.pos) + '\t' + var_id + '\t'
        vcf += str(self.region.seq[self.pos-1])
        vcf += '\t' + str(self.region.seq[self.pos-1]) + str(self.seq)
        vcf += '\t99\t.\tAC=1\tGT'
        for gt in self.gt:
            vcf += '\t' + gt
        return vcf

    def toSymbVcf(self):
        if(self.label == ''):
            var_id = 'ins' + str(self.pos) + '' + str(len(self.seq))
        else:
            var_id = self.label
        vcf = self.region.chrom + '\t' + str(self.pos+1) + '\t' + var_id + '\t'
        vcf += str(self.region.seq[self.pos])
        vcf += '\t<INS>\t99\t.\tAC=1;SVLEN={};CIPOS=-10,10;CIEND=-10,10;SVTYPE=INS;END={};SEQ={}\tGT'.format(self.size(), self.pos+1, str(self.seq))
        for gt in self.gt:
            vcf += '\t' + gt
        return vcf


class Deletion:
    def __init__(self, region, del_size, gt=[], label=''):
        self.region = region
        self.label = label
        self.pos = len(self.region.seq)
        self.size = del_size
        self.gt = list(gt)
        self.region.addBuffer(del_size)

    def shiftError(self, bp, bkpt='both'):
        self.label += 'sE{}{}'.format(bp, bkpt)
        if(bkpt == 'both'):
            self.pos += bp
        elif(bkpt == 'start'):
            self.pos += bp
            self.size += -1 * bp
        elif(bkpt == 'end'):
            self.size += bp

    def nbSamp(self):
        return len(self.gt)

    def toVcf(self):
        if(self.label == ''):
            var_id = 'del' + str(self.pos) + '' + str(self.size)
        else:
            var_id = self.label
        vcf = self.region.chrom + '\t' + str(self.pos) + '\t' + var_id + '\t'
        vcf += str(self.region.seq[(self.pos-1):(self.pos+self.size)])
        vcf += '\t' + str(self.region.seq[self.pos-1]) + '\t99\t.\tAC=1\tGT'
        for gt in self.gt:
            vcf += '\t' + gt
        return vcf

    def toSymbVcf(self):
        if(self.label == ''):
            var_id = 'del' + str(self.pos) + '' + str(self.size)
        else:
            var_id = self.label
        vcf = self.region.chrom + '\t' + str(self.pos) + '\t' + var_id + '\t'
        vcf += str(self.region.seq[self.pos-1])
        vcf += '\t<DEL>\t99\t.\tAC=1;SVLEN=-{};CIPOS=-10,10;CIEND=-10,10;SVTYPE=DEL;END={}\tGT'.format(self.size, self.pos+self.size)
        for gt in self.gt:
            vcf += '\t' + gt
        return vcf


class Inversion:
    def __init__(self, region, inv_size, gt=[], label=''):
        self.region = region
        self.label = label
        self.pos = len(self.region.seq)
        self.size = inv_size
        self.gt = list(gt)
        self.region.addBuffer(inv_size)

    def shiftError(self, bp, bkpt='both'):
        self.label += 'sE{}{}'.format(bp, bkpt)
        if(bkpt == 'both'):
            self.pos += bp
        elif(bkpt == 'start'):
            self.pos += bp
            self.size += -1 * bp
        elif(bkpt == 'end'):
            self.size += bp

    def nbSamp(self):
        return len(self.gt)

    def toVcf(self):
        if(self.label == ''):
            var_id = 'inv' + str(self.pos) + '' + str(self.size)
        else:
            var_id = self.label
        vcf = self.region.chrom + '\t' + str(self.pos) + '\t' + var_id + '\t'
        invseq = self.region.seq[self.pos:(self.pos+self.size)]
        invseq.reverse_complement()
        vcf += str(self.region.seq[(self.pos-1):(self.pos+self.size)])
        vcf += '\t' + str(self.region.seq[self.pos-1]) + str(invseq)
        vcf += '\t99\t.\tAC=1\tGT'
        for gt in self.gt:
            vcf += '\t' + gt
        return vcf

    def toSymbVcf(self):
        if(self.label == ''):
            var_id = 'inv' + str(self.pos) + '' + str(self.size)
        else:
            var_id = self.label
        # In vg, POS is not included in the inversion
        vcf = self.region.chrom + '\t' + str(self.pos) + '\t' + var_id + '\t'
        vcf += str(self.region.seq[self.pos-1])
        vcf += '\t<INV>\t99\t.\tAC=1;SVLEN=0;SVMETHOD=EMBL.DELLYv0.7.9;CT=5to5;CIPOS=-10,10;CIEND=-10,10;SVTYPE=INV;END={}\tGT'.format(self.pos+self.size)
        for gt in self.gt:
            vcf += '\t' + gt
        return vcf


class SNV:
    def __init__(self, region, pos, gt=[], label=''):
        self.region = region
        self.label = label
        self.pos = pos
        self.gt = list(gt)
        alts = ['A', 'T', 'C', 'G']
        rd.shuffle(alts)
        if alts[0] == str(self.region.seq[(self.pos-1)]):
            self.alt = alts[1]
        else:
            self.alt = alts[0]

    def nbSamp(self):
        return len(self.gt)

    def toVcf(self):
        if(self.label == ''):
            var_id = 'snv_' + str(self.pos)
        else:
            var_id = self.label
        vcf = self.region.chrom + '\t' + str(self.pos) + '\t' + var_id + '\t'
        vcf += str(self.region.seq[(self.pos-1)])
        vcf += '\t' + self.alt + '\t99\t.\tAC=1\tGT'
        for gt in self.gt:
            vcf += '\t' + gt
        return vcf


class Region:
    def __init__(self, size, chrom='chr1'):
        self.buff_size = size
        self.chrom = chrom
        seqArray = nuc[[int(rd.random()*4) for i in xrange(size)]]
        self.seq = MutableSeq("".join(seqArray), generic_dna)

    def addBuffer(self, extra_seq=0):
        buff_size = self.buff_size + extra_seq
        seqArray = nuc[[int(rd.random()*4) for i in xrange(buff_size)]]
        self.seq += MutableSeq("".join(seqArray), generic_dna)

    def writeRefFasta(self, fasta_file, desc=""):
        recs = []
        if(desc == ''):
            desc = 'length={}'.format(len(self.seq))
        recs.append(SeqRecord(self.seq, id=self.chrom, description=desc))
        SeqIO.write(recs, fasta_file, "fasta")


def writeVcf(svs, vcf_file, symbolic=False):
        vcff = open(vcf_file, 'w')
        vcff.write('##fileformat=VCFv4.1\n')
        if symbolic:
            vcff.write('##ALT=<ID=DEL,Description="Deletion">\n')
            vcff.write('##ALT=<ID=INS,Description="Insertion">\n')
            vcff.write('##ALT=<ID=INV,Description="Inversion">\n')
            vcff.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
            vcff.write('##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">\n')
            vcff.write('##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">\n')
            vcff.write('##INFO=<ID=SEQ,Number=1,Type=String,Description="Inserted sequence">\n')
            vcff.write('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
            vcff.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
            vcff.write('##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">\n')
            vcff.write('##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">\n')

        vcff.write('##INFO=<ID=AC,Number=A,Type=Integer,Description='
                   '"Alternate allele count">\n##FORMAT=<ID=GT,Number=1,'
                   'Type=String,Description="Genotype">\n')
        for samp in xrange(svs[0].nbSamp()):
            vcff.write('##SAMPLE=<ID=s' + str(samp) + '>\n')
        vcff.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT')
        for samp in xrange(svs[0].nbSamp()):
            vcff.write('\ts' + str(samp))
        vcff.write('\n')
        for var in svs:
            if symbolic:
                vcff.write(var.toSymbVcf() + '\n')
            else:
                vcff.write(var.toVcf() + '\n')
        vcff.close()
        return vcf_file


parser = argparse.ArgumentParser(description='Simulate ref and VCFs with ins/dels.')
parser.add_argument('-n', type=int, help='the number of del/ins to simulate',
                    required=True)
parser.add_argument('-r', help='a VCF with SVs (for size sampling)',
                    required=True)
parser.add_argument('-b', type=int, default=500,
                    help='the buffer size. Default 500 (bp).')

args = parser.parse_args()

# For each variant, 3 genotypes randomly attributed to a sample each
gts = ['1|1', '0|1', '0/0']

# Read VCF and save sizes for deletions and insertions
vcf_reader = vcf.Reader(filename=args.r)
del_s = []
ins_s = []
for record in vcf_reader:
    max_alt = 0
    for alt in record.ALT:
        max_alt = max(max_alt, len(alt))
    var_s = max_alt - len(record.REF)
    if abs(var_s) > 20: # Minimum 20 bp
        if var_s > 0:
            ins_s.append(var_s)
        else:
            del_s.append(abs(var_s))

# Regions and buffer size
reg = Region(size=args.b)
# Insertions
ins_l = []
for ii in range(args.n):
    ins_size = ins_s[int(rd.random()*len(ins_s))]
    ins_id = 'ins' + str(ii)
    rd.shuffle(gts)
    ins_l.append(Insertion(reg, ins_size, gt=gts, label=ins_id))
# Deletions
del_l = []
for ii in range(args.n):
    del_size = del_s[int(rd.random()*len(del_s))]
    del_id = 'del' + str(ii)
    rd.shuffle(gts)
    del_l.append(Deletion(reg, del_size, gt=gts, label=del_id))
# Inversions
inv_l = []
for ii in range(args.n):
    inv_size = del_s[int(rd.random()*len(del_s))]
    inv_id = 'inv' + str(ii)
    rd.shuffle(gts)
    inv_l.append(Inversion(reg, inv_size, gt=gts, label=inv_id))
# Write reference sequence fasta
reg.writeRefFasta(fasta_file='ref.fa')
# Write true VCF
writeVcf(ins_l + del_l + inv_l, vcf_file='truth.vcf')
writeVcf(ins_l + del_l + inv_l, vcf_file='truth.symb.vcf', symbolic=True)
# Redo with errors in calls
for ins in ins_l:
    shift_error = int(rd.random()*10) + 1
    ins.shiftError(shift_error)
    del_error = int(rd.random()*3)
    ins.delError(del_error)
error_type = ['both', 'start']
for de in del_l:
    shift_error = int(rd.random()*10) + 1
    type_error = int(rd.random()*2)
    de.shiftError(shift_error, error_type[type_error])
writeVcf(ins_l + del_l, vcf_file='calls.vcf')
for inv in inv_l:
    shift_error = int(rd.random()*10) + 1
    type_error = int(rd.random()*2)
    inv.shiftError(shift_error, error_type[type_error])
writeVcf(ins_l + del_l + inv_l, vcf_file='calls.vcf')
writeVcf(ins_l + del_l + inv_l, vcf_file='calls.symb.vcf', symbolic=True)
