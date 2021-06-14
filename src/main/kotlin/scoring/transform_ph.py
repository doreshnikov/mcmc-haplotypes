import os
from sys import argv

d = argv[1]
file = f'{argv[2]}.fas' if len(argv) > 2 else 'ph.fas'

lines = list(map(str.strip, open(f'{d}/{file}').readlines()))
res = {}

i = 0
while i < len(lines):
    if not lines[i].startswith('>'):
        i += 1
        continue
    p = float(lines[i].split()[1].split(':')[-1][:-1])
    i += 1
    s = ''
    while i < len(lines) and lines[i][0] != '>':
        s += lines[i]
        i += 1

    if s not in res.keys():
        res[s] = 0.0
    res[s] += p

excess = sum(res.values())
for s in res.keys():
    res[s] /= excess

with open(f'{d}/{file.split(".")[0]}', 'w') as ph:
    for s, p in res.items():
        print(s, p, file=ph)

print(f'total weight: {sum(res.values())}')
os.system(f'rm -f {d}/{file}')
