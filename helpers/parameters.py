from typing import NamedTuple, Any, List, Sequence, Tuple, Union
from itertools import groupby
from re import match, search

from fragment_capping.helpers.iterables import concat

FULL_VALENCES = {
    'C': 4,
    'N': 3,
    'O': 2,
    'H': 1,
    'S': 2,
    'P': 5,
    'CL': 1,
}

POSSIBLE_BOND_ORDERS = {
    'S': (1, 2,),
    'C': (1, 2,),
    'H': (1,),
    'O': (1, 2,),
    'N': (1, 2,),
    'P': (1, 2,),
    'CL': (1,),
}

POSSIBLE_CHARGES = {
    'S': (0,),
    'C': (0,),
    'H': (0,),
    'O': (0, -1,),
    'N': (0, +1,),
    'P': (0,),
    'CL': (0,),
}

Capping_Strategy = NamedTuple('Capping_Strategy', [('new_atoms', Sequence[str]), ('new_bonds', Sequence[Tuple[int, int]]), ('new_valences', Sequence[int])])

NO_CAP = Capping_Strategy((), (), ())
H_CAP = Capping_Strategy(('H',), ((0, 1),), (1,))
H2_CAP = Capping_Strategy(('H', 'H'), ((0, 1), (0, 2)), (1, 1))
H3_CAP = Capping_Strategy(('H', 'H', 'H'), ((0, 1), (0, 2), (0, 3)), (1, 1, 1))
H4_CAP = Capping_Strategy(('H', 'H', 'H', 'H'), ((0, 1), (0, 2), (0, 3), (0, 4)), (1, 1, 1, 1))
H_CH2_CAP = Capping_Strategy(('H', 'C', 'H', 'H'), ((0, 1), (0, 2), (2, 3), (2, 4)), (1, 3, 1, 1))
CH3_CAP = Capping_Strategy(('C', 'H', 'H', 'H'), ((0, 1), (1, 2), (1, 3), (1, 4)), (4, 1, 1, 1))

INDIVIDUAL_CAPPING_OPTIONS = {
    'H1': (NO_CAP,),
    'O1': (NO_CAP,),
    'O2': (H_CAP,),
    'S1': (NO_CAP,),
    'S2': (H_CAP,),
    'C4': (H3_CAP,),
    'C3': (H2_CAP, H_CH2_CAP),
    'N2': (H_CAP, CH3_CAP,),
    'N3': (H2_CAP,),
    'N4': (H3_CAP,),
    'P5': (H4_CAP,),
    'CL1': (NO_CAP,),
}

def group_if_not_none(match: Any) -> Union[str, None]:
    return match.group() if match is not None else None

def on_letters_of_dict_key(item: List[Any]) -> str:
    return group_if_not_none(match('[A-Z]+', item[0]))

def get_capping_options(use_neighbour_valences: bool, debug: bool = False) -> List[Capping_Strategy]:
    if not use_neighbour_valences:
        'Aggregate capping strategies for a given element.'
        capping_options = dict(
            [
                (element, concat([x[1] for x in list(group)]))
                for (element, group) in
                groupby(
                    sorted(
                        INDIVIDUAL_CAPPING_OPTIONS.items(),
                        key=on_letters_of_dict_key,
                    ),
                    key=on_letters_of_dict_key,
                )
            ]
        )
    else:
        capping_options = INDIVIDUAL_CAPPING_OPTIONS

    if debug:
        print([(key, value) for (key, value) in list(INDIVIDUAL_CAPPING_OPTIONS.items())])
        print([(key, value) for (key, value) in list(capping_options.items())])

    return capping_options
