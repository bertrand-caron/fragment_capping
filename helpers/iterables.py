from typing import Sequence, Any

def concat(list_of_lists: Sequence[Sequence[Any]]):
    return reduce(
        lambda acc, e: acc + e,
        list_of_lists,
        type(list_of_lists)(),
    )
