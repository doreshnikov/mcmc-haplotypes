from sys import argv
import os

os.chdir('../../resources')
coalevol_path = '../../../side-packages/CoalEvol-7.3.5/exe/'
art_path = '../../../side-packages/art_bin_MountRainier/'

S = 200
N = 5000
R = 1000
u = 1e-8

for d in argv[1:]:
    with open(f'{d}/_reference') as f:
        lines = f.readlines()
        reference = ''.join(list(map(str.strip, lines))).upper()
    size = len(reference)

    with open('_seq', 'w') as seq:
        print(reference, file=seq)
    os.system(f'{coalevol_path}/CoalEvol7.3.5 -n1 -s{S} {size} -e{N} 1 -u{u} -x_seq')

    with open('Results/sequences') as res:
        lines = list(filter(lambda s: len(s) > 0, map(str.strip, res.readlines())))[1:]
        sequences = list(map(lambda s: s.split()[-1], lines))
        C = len(sequences)
    data = {}
    for s in sequences:
        if s not in data.keys():
            data[s] = 1
        data[s] += 1
    for s in data.keys():
        data[s] /= C

    with open(f'{d}/_origin', 'w') as origin:
        for k, v in data.items():
            print(k, v, file=origin)
    for i, (k, v) in enumerate(data.items()):
        with open('_seq.fa') as haplo:
            print(k, file=haplo)
        reads = int(C * v)
        os.system(f'{art_path}/art_illumina -ss MSv3 -sam -i _seq.fa -l 300 -c {reads} -o Results/reads_{i}')

    os.remove('_seq')
    os.remove('_seq.fa')
    print(f'Sample {d} done')