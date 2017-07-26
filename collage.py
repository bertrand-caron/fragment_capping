from typing import List, Tuple, Union

Numeric = Union[int, float]

COLUMNS, LINES = range(2)

MEAN_IMAGE_ASPECT_RATIO = 4. / 3.

def best_grid(l: int, aspect_ratio: Tuple[Numeric, Numeric] = (1, 1)) -> Tuple[int, int]:
    '''Returns (n, m) such as n * m >= l and the grid is as dense as posible'''
    from math import ceil, sqrt
    from itertools import product

    target_aspect_ratio = (aspect_ratio[1] / aspect_ratio[0]) / MEAN_IMAGE_ASPECT_RATIO

    possible_grids = product(range(1, l + 1), repeat=2)

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

def indices_for_subplot(n: int, subplot_dims: Tuple[int, int], always_return_tuple: bool = False) -> Union[Tuple[int, int], int]:
    indices = (
        n // subplot_dims[LINES],
        n - (n // subplot_dims[LINES]) * subplot_dims[LINES],
    )

    assert 0 <= indices[COLUMNS] < subplot_dims[COLUMNS]
    assert 0 <= indices[LINES] < subplot_dims[LINES]

    if subplot_dims[COLUMNS] == 1 and not always_return_tuple:
        return indices[LINES]
    elif subplot_dims[LINES] == 1 and not always_return_tuple:
        return indices[COLUMNS]
    else:
        return indices
