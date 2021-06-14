import matplotlib.pyplot as plt
from math import log2


def f(x):
    return list(map(lambda t: log2(t) if t is not None else None, x[1:]))


names = ['Reference', 'CliqueSNV', 'PredictHaplo', 'MCMC']
styles = ['-.b', '-or', '-xy' ,'-^g']
unifrac = {
    500: [
        [0.0, 0.008, 0.017, 1.464, 2.208, 8.941, 21.947, 86.229],
        [0.0, 0.008, 0.017, 1.272, 2.118, 4.889, None, None],
        [0.0, 0.0082, 0.0189, 1.464, 2.118, 5.935, 18.598, 18.878],
        [0.0, 0.0006, 0.0076, 0.202, 1.659, 5.458, 18.159, 59.621],
    ],
    1000: [
        [0.0, 0.006, 1.350, 0.929, 3.255, 14.918, 57.114, 130.574],
        [0.0, 0.006, 1.350, 0.812, 2.233, 4.007, None, None],
        [0.0, 0.0076, 1.349, 0.939, 2.792, 5.528, 8.241, 18.437],
        [0.0, 0.00064, 0.0127, 0.632, 3.638, 5.683, 48.551, 145.358]
    ],
    2500: [
        [0.0, 0.004, 0.604, 8.152, 13.82, 43.07, 100.627, None],
        [0.0, 0.004, 0.604, 3.907, 5.791, None, None, None],
        [0.0, 0.0071, 0.599, 6.228, 10.122, 11.555, 14.639, None],
        [0.0, 0.00027, 0.537, 3.954, 9.143, 65.344, None, None]
    ]
}

for n, v in unifrac.items():
    plt.clf()

    ticks = [8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5]
    if v[0][0] == 0.0:
        ticks = ticks[1:]
    plt.xticks(range(len(ticks)), ticks)

    for i, score in enumerate(v):
        plt.plot(f(score), styles[i], label=names[i])

    plt.xlabel('log10 μ')
    plt.ylabel('log2 UniFrac')
    plt.title(f'Ne = {n}')
    plt.legend(loc='right')

    plt.savefig(f'compare-unifrac-{n}.png')
