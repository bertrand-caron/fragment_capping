from typing import Sequence, Any
from functools import reduce

from fragment_capping.helpers.exceptions import Too_Many_Permutations

def concat(list_of_lists: Sequence[Sequence[Any]]):
    if len(list_of_lists) == 0:
        return list_of_lists
    else:
        return reduce(
            lambda acc, e: acc + e,
            list_of_lists,
        )

assert concat([[1], [2], [3]]) == [1, 2, 3], concat([[1], [2], [3]])

MAXIMUM_PERMUTATION_NUMBER = 600_000

def product_len(list_of_lists: Sequence[Sequence[Any]]) -> int:
    _product_len = reduce(lambda acc, e: acc * len(e), list_of_lists, 1)
    if _product_len > MAXIMUM_PERMUTATION_NUMBER and False:
        raise Too_Many_Permutations(_product_len)
    return _product_len
