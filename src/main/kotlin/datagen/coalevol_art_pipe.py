from sys import argv
from random import shuffle
import os
import argparse

os.chdir('../../resources')
coalevol_path = '../../../side-packages/CoalEvol-7.3.5/exe/'
art_path = '../../../side-packages/art_bin_MountRainier/'
align_exe = '../../../side-packages/align'

parser = argparse.ArgumentParser()
parser.add_argument('dirs', type=str)
parser.add_argument('-n', action='store', type=int)
parser.add_argument('-s', action='store', type=int, default=0)
parser.add_argument('-l', action='store', type=int, default=150)
parser.add_argument('-u', action='store', type=float)
parser.add_argument('-r', action='store', type=int, default=0)

parser.add_argument('--errfree', action='store_true')
parser.add_argument('--noclean', action='store_true')

args = parser.parse_args()
N, S, u, R = args.n, args.s, args.u, args.r
readlen = args.l
if S == 0:
    S = N

for d in argv[1].split('/'):
    os.system(f'rm -rf {d}/reads_*')
    os.system(f'rm -rf {d}/_aligned')

    with open(f'{d}/_reference') as f:
        lines = f.readlines()
        reference = ''.join(list(map(str.strip, lines))).upper()
    size = len(reference)
    if R == 0:
        R = size // 5

    with open('seqGMRCA', 'w') as seq:
        print(reference, file=seq)
    os.system(f'{coalevol_path}/CoalEvol7.3.5 -n1 -s{S} {size} -e{N} 1 -u{u} -xseqGMRCA >/dev/null')

    with open('Results/sequences') as res:
        lines = list(filter(lambda s: len(s) > 0, map(str.strip, res.readlines())))[1:]
        sequences = list(map(lambda s: s.split()[-1], lines))
        C = len(sequences)
    data = {}
    for s in sequences:
        if s not in data.keys():
            data[s] = 0
        data[s] += 1
    for s in data.keys():
        data[s] /= C

    with open(f'{d}/_origin', 'w') as origin:
        for k, v in data.items():
            print(k, v, file=origin)
    for i, (k, v) in enumerate(data.items()):
        with open('_seq.fa', 'w') as haplo:
            print(f'>{d}', file=haplo)
            idx = 0
            while idx < len(k):
                print(k[idx:idx + 70], file=haplo)
                idx += 70
        reads = int(R * v)
        cmd = f'{art_path}/art_illumina -ss MSv3 -i _seq.fa -l {readlen} -c {reads} -o {d}/reads_{i}'
        if args.errfree:
            cmd += ' --errfree'
        os.system(f'{cmd} >/dev/null')

    total = []
    if not args.errfree:
        for i in range(len(data)):
            with open(f'{d}/reads_{i}.fq') as reads:
                lines = list(map(str.strip, reads.readlines()))
            idx = 1
            while idx < len(lines):
                total.append(lines[idx])
                idx += 4
    else:
        for i in range(len(data)):
            with open(f'{d}/reads_{i}_errFree.sam') as reads:
                lines = list(map(str.strip, reads.readlines()))
            for line in lines:
                if line[0] == '@':
                    continue
                total.append(line.split('\t')[-2])

    shuffle(total)
    with open(f'{d}/_reads', 'w') as all_reads:
        for read in total:
            print(read, file=all_reads)

    os.system(f'g++ ../kotlin/aligner/edlib_aligner.cpp -std=c++17 -ledlib -lpthread -o {align_exe}')
    os.system(f'{align_exe} {d}')

    if not args.noclean:
        os.remove('seqGMRCA')
        os.remove('_seq.fa')
        os.system(f'rm -rf {d}/reads_*')

    print(f'Sample {d} done')
