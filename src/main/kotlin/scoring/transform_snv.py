import os
from sys import argv
import os

d = argv[1]
lines = list(map(str.strip, open(f'{d}/snv.fasta').readlines()))

with open(f'{d}/snv', 'w') as snv:
    for i in range(len(lines)):
        if not lines[i].startswith('>'):
            continue
        p = float(lines[i].split('_')[-1])
        s = lines[i + 1]
        print(s, p, file=snv)

os.system(f'rm -f {d}/snv.fasta')