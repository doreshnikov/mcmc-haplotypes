from sys import argv


def read(fname):
    return list(map(str.strip, open(fname).readlines()))


d = argv[1]
rbwa = read(f'{d}/reads_bwa.sam')
rerf = read(f'{d}/reads_errFree.sam')

idx = 0
while rerf[idx].startswith('@'):
    idx += 1
rerf = list(map(str.split, rerf[idx:]))
idx = 0

with open(f'{d}/reads_comb.sam', 'w') as out:
    for line in rbwa:
        if line.startswith('@'):
            print(line, file=out)
        else:
            items = line.split()
            while rerf[idx][0] != items[0]:
                idx += 1
            items[9] = rerf[idx][9]
            print('\t'.join(items), file=out)