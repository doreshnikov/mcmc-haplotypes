from sys import argv
import os

os.chdir('../../resources')
exe_path = 'side-packages/CoalEvol-7.3.5/exe/'

S = 200
N = 5000
u = 1e-8

for d in argv[1:]:
    with open(f'{d}/_reference') as f:
        lines = f.readlines()
        reference = ''.join(list(map(str.strip, lines))).upper()
    size = len(reference)

    os.chdir(f'../../../{exe_path}')
    with open('seqGMRCA', 'w') as seq:
        print(reference, file=seq)
    os.system(f'./CoalEvol7.3.5 -n1 -s{S} {size} -e{N} 1 -u{u} -xseqGMRCA')

    with open('Results/sequences') as res:
        lines = list(filter(lambda s: len(s) > 0, map(str.strip, res.readlines())))[1:]
        sequences = list(map(lambda s: s.split()[-1], lines))
    data = {}
    for s in sequences:
        if s not in data.keys():
            data[s] = 1
        data[s] += 1
    for s in data.keys():
        data[s] /= len(sequences)

    os.chdir('../../../src/main/resources')
    with open(f'{d}/_origin', 'w') as origin:
        for k, v in data.items():
            print(k, v, file=origin)

    print(f'Sample {d} done')