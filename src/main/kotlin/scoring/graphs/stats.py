import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-o', action='store', type=str, default='stats')
args = parser.parse_args()

f = open(f'../stats')
history = list(map(str.strip, f.readlines()))
res = {}
for l in history:
    items = l.split()
    key = items[0][:-1]
    if key not in res.keys():
        res[key] = []
    res[key].append(float(items[-1]))

for k, v in res.items():
    plt.hist(v)
    plt.title(f'log g/g on \'{k}\'')
    plt.xlabel(f'avg = {sum(v) / len(v)}')
    plt.savefig(f'{args.o}_{k}')
    plt.clf()