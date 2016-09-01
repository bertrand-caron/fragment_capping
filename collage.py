COLUMNS, LINES = range(2)

def best_grid(l):
    '''Returns (n, m) such as n * m >= l and the grid is as dense as posible'''
    from math import ceil, sqrt
    from itertools import product

    highest_dim = int(ceil(sqrt(l)))
    possible_grids = product(range(1, highest_dim + 1), repeat=2)
    sorted_grids = sorted(
        filter(
            lambda grid: grid[LINES] * grid[COLUMNS] >= l,
            possible_grids,
        ),
        key=lambda grid: (grid[LINES] * grid[COLUMNS]) - l,
    )
    return sorted_grids[0]

def indices_for_subplot(n, subplot_dims):
    return ((n // subplot_dims[LINES]), n - (n // subplot_dims[LINES]) * subplot_dims[LINES])
