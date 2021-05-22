from sys import argv
from random import shuffle
import os

os.chdir('../../resources')
coalevol_path = '../../../side-packages/CoalEvol-7.3.5/exe/'
art_path = '../../../side-packages/art_bin_MountRainier/'

N = 5000
# S = 200
S = N
R = 10000
u = 1e-8

for d in argv[1:]:
    os.system(f'rm -rf {d}/reads_*')

    with open(f'{d}/_reference') as f:
        lines = f.readlines()
        reference = ''.join(list(map(str.strip, lines))).upper()
    size = len(reference)

    with open('seqGMRCA', 'w') as seq:
        print(reference, file=seq)
    os.system(f'{coalevol_path}/CoalEvol7.3.5 -n1 -s{S} {size} -e{N} 1 -u{u} -xseqGMRCA')

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
        os.system(f'{art_path}/art_illumina -ss MSv3 -i _seq.fa -l 250 -c {reads} -o {d}/reads_{i}')

    total = []
    for i in range(len(data)):
        with open(f'{d}/reads_{i}.fq') as reads:
            lines = list(map(str.strip, reads.readlines()))
        idx = 1
        while idx < len(lines):
            total.append(lines[idx])
            idx += 4
    shuffle(total)
    with open(f'{d}/_reads', 'w') as all_reads:
        for read in total:
            print(read, file=all_reads)

    os.remove('seqGMRCA')
    os.remove('_seq.fa')
    os.system(f'rm -rf {d}/reads_*')

    print(f'Sample {d} done')
