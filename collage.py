COLUMNS, LINES = range(2)

MEAN_IMAGE_ASPECT_RATIO = 4. / 3.

def best_grid(l, aspect_ratio=(1, 1)):
    '''Returns (n, m) such as n * m >= l and the grid is as dense as posible'''
    from math import ceil, sqrt
    from itertools import product

    target_aspect_ratio = (aspect_ratio[1] / aspect_ratio[0]) / MEAN_IMAGE_ASPECT_RATIO

    possible_grids = product(range(1, l), repeat=2)

    def grid_aspect_ratio(grid):
        return grid[COLUMNS] / grid[LINES]

    def grid_fitness(grid):
        fitness = (
            abs(
                grid_aspect_ratio(grid) - target_aspect_ratio,
            )
            +
            0.1 * ((grid[LINES] * grid[COLUMNS]) - l),
        )
        return fitness

    sorted_grids = sorted(
        filter(
            lambda grid: grid[LINES] * grid[COLUMNS] >= l,
            possible_grids,
        ),
        key=grid_fitness,
    )

    return sorted_grids[0]

def indices_for_subplot(n, subplot_dims):
    return ((n // subplot_dims[LINES]), n - (n // subplot_dims[LINES]) * subplot_dims[LINES])
