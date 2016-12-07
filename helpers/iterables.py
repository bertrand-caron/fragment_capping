from typing import Sequence, Any
from functools import reduce

def concat(list_of_lists: Sequence[Sequence[Any]]):
    if len(list_of_lists) == 0:
        return list_of_lists
    else:
        return reduce(
            lambda acc, e: acc + e,
            list_of_lists,
        )

assert concat([[1], [2], [3]]) == [1, 2, 3], concat([[1], [2], [3]])
