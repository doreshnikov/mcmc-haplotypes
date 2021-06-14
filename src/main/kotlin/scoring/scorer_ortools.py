import editdistance
from ortools.linear_solver import pywraplp
from collections import defaultdict
from sys import argv
import os


def read_data(filename):
    def parse(line: str):
        items = line.strip().split()
        return items[0], float(items[1])

    with open(filename) as f:
        lines = f.readlines()
        return list(map(parse, lines))


def earth_mover_distance(dist1, dist2):
    solver = pywraplp.Solver(
        'earth_mover_distance',
        pywraplp.Solver.GLOP_LINEAR_PROGRAMMING
    )

    variables = dict()

    dirt_leaving_constraints = defaultdict(lambda: 0)
    dirt_filling_constraints = defaultdict(lambda: 0)

    objective = solver.Objective()
    objective.SetMinimization()

    for x, dirt_at_x in dist1:
        for y, capacity_of_y in dist2:
            amount_to_move_x_y = solver.NumVar(0, solver.infinity(),
                                               'z_{%s, %s}' % (x, y))
            variables[(x, y)] = amount_to_move_x_y
            dirt_leaving_constraints[x] += amount_to_move_x_y
            dirt_filling_constraints[y] += amount_to_move_x_y
            objective.SetCoefficient(amount_to_move_x_y,
                                     editdistance.eval(x, y))

    dist1 = {k: v for k, v in dist1}
    dist2 = {k: v for k, v in dist2}
    for x, linear_combination in dirt_leaving_constraints.items():
        solver.Add(linear_combination == dist1[x])

    for y, linear_combination in dirt_filling_constraints.items():
        solver.Add(linear_combination == dist2[y])

    status = solver.Solve()
    if status not in [solver.OPTIMAL, solver.FEASIBLE]:
        raise Exception('Unable to find feasible solution')

    return objective.Value()


if __name__ == '__main__':
    dd = argv[1]
    for d in dd.split('/'):
        print(f'Scoring {d}:')
        expect = read_data(f'{d}/_origin')
        for fname in os.listdir(d):
            # only new results
            if not fname.startswith('ph') and not fname.startswith('_'):
                continue
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
            try:
                print(f' - {fname}: {earth_mover_distance(expect, data)}')
            except Exception:
                print(f' - {fname}: No solution found')
