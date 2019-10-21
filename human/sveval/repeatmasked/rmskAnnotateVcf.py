import argparse
import subprocess


parser = argparse.ArgumentParser(description='Annotate VCF with RepeatMasker')
parser.add_argument('-i', help='the input VCF', required=True)
parser.add_argument('-o', help='the output VCF', required=True)
parser.add_argument('-t', help='number of threads to use', default=1, type=int)
parser.add_argument('-s', help='the minimum size', default=20, type=int)
args = parser.parse_args()

# Read input VCF and make FASTA file
temp_fa = args.i + '.fa'
tfile = open(temp_fa, 'w')
svids = []
seq_lens = {}
for line in open(args.i, 'r'):
    if line[0] == '#':
        continue
    line = line.rstrip().split('\t')
    ref_seq = line[3]
    alt_seq = line[4]
    if abs(len(ref_seq)-len(alt_seq)) < args.s:
        continue
    svid = '{}_{}_{}'.format(line[0], line[1],
                             len(ref_seq)-len(alt_seq))
    svids.append(svid)
    tfile.write('>' + svid + '\n')
    if len(ref_seq) > len(alt_seq):
        seq = ref_seq
    else:
        seq = alt_seq
    tfile.write(seq + '\n')
    seq_lens[svid] = len(seq)

# Run RepeatMasker
rm_cmd = ['RepeatMasker', temp_fa, '--species', 'human', '-pa', str(args.t)]
dump = open('/dev/null')
rm_out = subprocess.check_output(rm_cmd, stderr=dump)
dump.close()

# Read RepeatMasker output
repsv = {}
for line in open(temp_fa + '.out', 'r'):
    line = line.rstrip().split()
    if len(line) < 4:
        continue
    svid = line[4]
    if svid in svids:
        repw = int(line[6]) - int(line[5])
        if svid not in repsv or repsv[svid]['w'] < repw:
            repsv[svid] = {}
            repsv[svid]['w'] = repw
            repsv[svid]['name'] = line[9]
            repsv[svid]['classfam'] = line[10]

# Write output VCF
outf = open(args.o, 'w')
info_added = False
new_info = '''##INFO=<ID=RMSKNAME,Number=1,Type=String,Description="Repeat name">
##INFO=<ID=RMSKCLASS,Number=1,Type=String,Description="Repeat class and family">
##INFO=<ID=RMSKCOV,Number=1,Type=Float,Description="Coverage of the repeat">
'''
for line in open(args.i, 'r'):
    line = line.rstrip()
    if line[0] == '#':
        if '##INFO=' in line and not info_added:
            outf.write(new_info)
            info_added = True
        outf.write(line + '\n')
        continue
    line = line.split('\t')
    ref_seq = line[3]
    alt_seq = line[4]
    svid = '{}_{}_{}'.format(line[0], line[1],
                             len(ref_seq)-len(alt_seq))
    if abs(len(ref_seq)-len(alt_seq)) < args.s or svid not in repsv:
        outf.write('\t'.join(line) + '\n')
        continue
    newf = 'RMSKNAME={};RMSKCLASS={};RMSKCOV={}'
    newf = newf.format(repsv[svid]['name'],
                       repsv[svid]['classfam'],
                       round(float(repsv[svid]['w']) / seq_lens[svid], 3))
    sep = ''
    if len(line[7]) > 0:
        sep = ';'
    line[7] = newf + sep + line[7]
    outf.write('\t'.join(line) + '\n')
outf.close()
