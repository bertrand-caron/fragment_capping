from typing import NamedTuple, Any, List, Sequence, Tuple, Union, Dict, Set, Optional
from itertools import groupby, product, combinations_with_replacement
from functools import reduce
from re import match, search
from numpy import array as vector
from numpy.random import uniform
from numpy.linalg import norm

from fragment_capping.helpers.iterables import concat
from fragment_capping.helpers.types_helpers import Atom

FULL_VALENCES = {
    'C': {4},
    'N': {3, 5},
    'O': {2},
    'H': {1},
    'S': {2, 4, 6},
    'P': {3, 5},
    'CL': {1},
    'BR': {1},
    'F': {1},
    'I': {1},
    'B': {3},
    'SI': {4},
}

POSSIBLE_BOND_ORDERS = {
    'S': {1, 2},
    'C': {1, 2, 3},
    'H': {1},
    'O': {1, 2},
    'N': {1, 2, 3},
    'P': {1, 2},
    'CL': {1},
    'BR': {1},
    'F': {1},
    'I': {1},
    'B': {1},
    'SI': {1, 2, 3},
}

POSSIBLE_CHARGES = {
    'S': {0},
    'C': {0},
    'H': {0},
    'O': {0, -1},
    'N': {0, +1},
    'P': {0},
    'CL': {0},
    'BR': {0},
    'F': {0},
    'I': {0},
    'B': {0},
    'SI': {0},
}

assert set(FULL_VALENCES.keys()) == set(POSSIBLE_BOND_ORDERS.keys()) == set(POSSIBLE_CHARGES.keys())

ALL_ELEMENTS = set(FULL_VALENCES.keys())

Capping_Strategy = NamedTuple(
    'Capping_Strategy',
    [
        ('new_atoms', Sequence[str]),
        ('new_bonds', Sequence[Tuple[int, int]]),
        ('new_valences', Sequence[int]),
    ],
)

NO_CAP = Capping_Strategy((), (), ())
H_CAP = Capping_Strategy(('H',), ((0, 1),), (1,))
H2_CAP = Capping_Strategy(('H', 'H'), ((0, 1), (0, 2)), (1, 1))
H3_CAP = Capping_Strategy(('H', 'H', 'H'), ((0, 1), (0, 2), (0, 3)), (1, 1, 1))
H4_CAP = Capping_Strategy(('H', 'H', 'H', 'H'), ((0, 1), (0, 2), (0, 3), (0, 4)), (1, 1, 1, 1))
H_CH2_CAP = Capping_Strategy(('H', 'C', 'H', 'H'), ((0, 1), (0, 2), (2, 3), (2, 4)), (1, 3, 1, 1))
CH3_CAP = Capping_Strategy(('C', 'H', 'H', 'H'), ((0, 1), (1, 2), (1, 3), (1, 4)), (4, 1, 1, 1))
O3_CAP = Capping_Strategy(('O', 'O', 'O'), ((0, 1), (0, 2), (0, 3)), (2, 1, 1))

INDIVIDUAL_CAPPING_OPTIONS = {
    'H1': [NO_CAP],
    'O1': [NO_CAP],
    'O2': [H_CAP],
    'S1': [NO_CAP],
    'S2': [H_CAP],
    'S4': [NO_CAP],
    'S6': [NO_CAP],
    'C4': [H_CAP, H2_CAP, H3_CAP],
    'C3': [H_CAP, H2_CAP, H_CH2_CAP],
    'C2': [H_CAP],
    'N1': [NO_CAP],
    'N2': [H_CAP, CH3_CAP],
    'N3': [H2_CAP],
    'N4': [H3_CAP],
    'P5': [O3_CAP],
    'CL1': [NO_CAP],
    'BR1': [NO_CAP],
    'F1': [NO_CAP],
    'I1': [NO_CAP],
    'B3': [H2_CAP, H_CAP],
    'SI4': [NO_CAP],
    'SI3': [NO_CAP],
    'SI2': [NO_CAP],
}

ALL_ELEMENT_VALENCE_COMBINATIONS = [
    '{element}{neighbour_count}'.format(element=element, neighbour_count=valence)
    for (element, valence) in reduce(
        lambda acc, e: acc + e,
        [list(product([element], valences))for (element, valences) in FULL_VALENCES.items()],
        [],
    )
]

POSSIBLE_BOND_ORDER_FOR_PAIR = {
    (element_1, element_2): POSSIBLE_BOND_ORDERS[element_1] & POSSIBLE_BOND_ORDERS[element_2]
    for (element_1, element_2) in product(POSSIBLE_BOND_ORDERS.keys(), repeat=2)
}

def max_valence(atom: Atom) -> int:
    return max(FULL_VALENCES[atom.element]) - min(POSSIBLE_CHARGES[atom.element])

def min_valence(atom: Atom) -> int:
    return (min(FULL_VALENCES[atom.element]) + max(POSSIBLE_CHARGES[atom.element])) // max(POSSIBLE_BOND_ORDERS[atom.element])

def possible_sets_of_bond_orders_for_atom(atom: Atom) -> List[int]:
    if atom.valence is None:
        return POSSIBLE_BOND_ORDERS[atom.element]
    else:
        if POSSIBLE_CHARGES[atom.element] == {0}:
            return reduce(
                lambda acc, e: acc | e,
                [
                    set(bond_orders)
                    for bond_orders in combinations_with_replacement(POSSIBLE_BOND_ORDERS[atom.element], r=atom.valence)
                    if sum(bond_orders) in FULL_VALENCES[atom.element]
                ],
                set(),
            )
        else:
            return POSSIBLE_BOND_ORDERS[atom.element]

def possible_bond_order_for_atom_pair(atoms: Tuple[Atom, Atom]) -> List[int]:
    return possible_sets_of_bond_orders_for_atom(atoms[0]) & possible_sets_of_bond_orders_for_atom(atoms[1])

def possible_charge_for_atom(atom: Atom) -> Set[int]:
    if atom.valence is None:
        return POSSIBLE_CHARGES[atom.element]
    else:
        if len(FULL_VALENCES[atom.element]) == 1 and atom.valence in FULL_VALENCES[atom.element]:
            return {0}
        else:
            return POSSIBLE_CHARGES[atom.element]

def new_atom_for_capping_strategy(capping_strategy: Capping_Strategy) -> int:
    return len([1 for (atom_id_1, atom_id_2) in capping_strategy.new_bonds if atom_id_1 == 0 or atom_id_2 == 0])

def group_if_not_none(match: Any) -> Union[str, None]:
    return match.group() if match is not None else None

def on_letters_of_dict_key(item: List[Any]) -> str:
    return group_if_not_none(match('[A-Z]+', item[0]))

def get_capping_options(use_neighbour_valences: bool, debug: bool = False) -> Dict[str, List[Capping_Strategy]]:
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

DOUBLE_BOND_ELEMENTS = [
    element
    for element in ALL_ELEMENTS
    if 2 in POSSIBLE_BOND_ORDERS[element]
]

ALL_POSSIBLE_DOUBLE_BONDS = {
    frozenset([element_1, element_2])
    for (element_1, element_2) in product(DOUBLE_BOND_ELEMENTS, DOUBLE_BOND_ELEMENTS)
}

ELECTRONEGATIVITIES = {
    # Source: https://en.wikipedia.org/wiki/Electronegativity
    'C': 2.55,
    'N': 3.04,
    'O': 3.44,
    'H': 2.20,
    'S': 2.58,
    'P': 2.19,
    'CL': 3.16,
    'BR': 2.96,
    'F': 3.98,
    'I': 2.66,
    'B': 2.04,
    'SI': 1.90,
}

assert set(ELECTRONEGATIVITIES.keys()) >= ALL_ELEMENTS, ALL_ELEMENTS - set(ELECTRONEGATIVITIES.keys())

def electronegativity_spread(elements: Tuple[str, str]) -> float:
    assert len(elements) > 0, elements

    electronegativites = list(map(
        lambda element: ELECTRONEGATIVITIES[element],
        elements,
    ))

    return max(electronegativites) - min(electronegativites)

def coordinates_n_angstroms_away_from(atom: Atom, n: float) -> Optional[Tuple[float, float, float]]:
    if atom.coordinates is None:
        return None
    else:
        random_vector = uniform(-1, 1, 3)
        return tuple((vector(atom.coordinates) + random_vector / norm(random_vector) * n).tolist())

if __name__ == '__main__':
    print(
        sorted(
            [
                (
                    element + str(valence),
                    possible_sets_of_bond_orders_for_atom(
                        Atom(index=None, element=element, valence=valence, capped=True, coordinates=None),
                    ),
                )
                for (element, valence) in map(lambda x: (x[:-1], int(x[-1])), INDIVIDUAL_CAPPING_OPTIONS.keys())
            ],
        ),
    )

    print(coordinates_n_angstroms_away_from(Atom(index=None, element=None, valence=None, capped=False, coordinates=(1. ,2. ,3.)), 1.2))
