from sys import argv
import os


def read_data(filename):
    def parse(line: str):
        items = line.strip().split()
        return items[0], float(items[1])

    with open(filename) as f:
        lines = f.readlines()
        return list(map(parse, lines))


def dist(s1, s2):
    c = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            c += 1
    return c


def stats(exp, act, t=0):
    tp, fp = 0.0, 0.0
    for h, p in exp:
        total = 0.0
        for r, q in act:
            if dist(h, r) < t:
                total += q
        if total > 0.0:
            tp += p
    for r, q in act:
        total = 0.0
        for h, p in exp:
            if dist(h, r) < t:
                total += p
        if total == 0.0:
            fp += q

    return tp / (tp + fp), tp


if __name__ == '__main__':
    dd = argv[1]
    for d in dd.split('/'):
        print(f'Scoring {d}:')
        expect = read_data(f'{d}/_origin')
        for fname in os.listdir(d):
            # # only new results
            # if not fname.startswith('ph'):
            #     continue
            if fname[0] == '_' and fname != '_reference':
                continue
            fail = map(fname.startswith, ['@', 'reads', 'sequences'])
            if any(fail) or fname.endswith('.png'):
                continue
            if fname == '_reference':
                data = ''.join(list(map(str.strip, open(f'{d}/{fname}').readlines())))
                data = [(data, 1.0)]
            else:
                data = read_data(f'{d}/{fname}')
            # print(f' - {fname}: {stats(expect, data, len(expect[0][0]) // 100)}')
            print(f' - {fname}: {stats(expect, data, 9)}')
