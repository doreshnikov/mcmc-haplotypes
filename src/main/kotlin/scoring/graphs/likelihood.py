import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('dirs', type=str)
parser.add_argument('-o', action='store', type=str, default='history')
parser.add_argument('-i', action='store', type=str, default='')

args = parser.parse_args()

for d in args.dirs.split('/'):
    fname = '_loglhistory'
    if len(args.i) > 0:
        fname = fname + '-' + args.i
    f = open(f'{d}/{fname}')
    history = list(map(str.strip, f.readlines()))
    if history[0][0] == '(':
        history = list(map(lambda x: x[:-1].split()[-1], history))
    history = list(map(float, history))

    plt.plot(history)
    plt.xlabel('Итерация')
    plt.ylabel('Правдоподобие')
    plt.savefig(f'{d}/{args.o}.png')