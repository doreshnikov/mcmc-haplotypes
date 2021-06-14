import argparse
import os

coalevol_path = '../side-packages/CoalEvol-7.3.5/exe/'
art_path = '../side-packages/art_bin_MountRainier/'
align_exe = '../side-packages/align'

parser = argparse.ArgumentParser()
parser.add_argument('-e', action='store', type=int)
parser.add_argument('-s', action='store', type=int, default=0)
parser.add_argument('-u', action='store', type=str)
parser.add_argument('-r', action='store', type=int)
parser.add_argument('-f', action='store', type=int, default=10)

parser.add_argument('-o', action='store', type=str, default='test')

args = parser.parse_args()
s = args.s if args.s > 0 else 1000

cmd = f'{coalevol_path}/CoalEvol7.3.5 -n1 -s{s} {args.r} -e{args.e} 1 -u{args.u} -$ -bsequences'
print(cmd)
os.system(cmd)
sequences = list(map(str.strip, open('Results/sequences').readlines()))[1:-1]
data = list(map(str.split, sequences))

d = {}
if args.o not in os.listdir():
    os.mkdir(args.o)

with open(f'{args.o}/sequences.fa', 'w') as fa:
    for item in data:
        seq = item[-1]
        if seq in d.keys():
            d[seq] += 1
        else:
            d[seq] = 1
        print(f'>{item[0]}\n{seq}', file=fa)

ref = open('Results/GMRCA00001').readlines()[1]
sequence = ref.strip().split()[-1]
with open(f'{args.o}/_reference', 'w') as reference:
    print(sequence, file=reference)
with open(f'{args.o}/_reference.fa', 'w') as referenceF:
    lines = []
    i = 0
    while i < len(sequence):
        lines.append(sequence[i: i + 70])
        i += 70
    print(f'>{args.o}', file=referenceF)
    for line in lines:
        print(line, file=referenceF)

total = len(data)
for k in d.keys():
    d[k] /= total
with open(f'{args.o}/_origin', 'w') as origin:
    for k, v in d.items():
        print(k, v, file=origin)

cmd = f'{art_path}/art_illumina -ss MSv3 -i {args.o}/sequences.fa -l 250 -f {args.f} -o {args.o}/reads -ef'
os.system(cmd)
sam = list(map(str.strip, open(f'{args.o}/reads_errFree.sam').readlines()))
i = 0
while sam[i][0] == '@':
    i += 1
sam = list(map(str.split, sam[i:]))
with open(f'{args.o}/_aligned', 'w') as aligned:
    for item in sam:
        print(item[9], int(item[3]) - 1, file=aligned)