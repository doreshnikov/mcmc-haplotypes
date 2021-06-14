import os
from sys import argv
import os

d = argv[1]
ref = open(f'{d}/_reference').readlines()[0].strip()
lines = list(map(str.strip, open(f'{d}/abqr.txt').readlines()))

with open(f'{d}/abqr', 'w') as snv:
    for i in range(len(lines)):
        if not lines[i].startswith('V'):
            continue
        p = float(lines[i].split()[-1])
        s = lines[i + 1]
        for i in range(len(s)):
                if s[i] == '-':
                        s[i] = ref[i]
        print(s, p, file=snv)

os.system(f'rm -f {d}/abqr.txt')
