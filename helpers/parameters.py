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
    'SE': {2, 4, 6},
    'P': {3, 5},
    'CL': {1},
    'BR': {1},
    'F': {1},
    'I': {1, 3, 5},
    'B': {3, 5},
    'SI': {4},
}

VALENCE_ELECTRONS = {
    'H':  1, 'HE': 2,
    'LI': 1, 'BE': 2, 'B':  3, 'C':  4, 'N':  5, 'O':  6, 'F':  7, 'NE': 8,
    'NA': 1, 'MG': 2, 'AL': 3, 'SI': 4, 'P':  5, 'S':  6, 'CL': 7, 'AR': 8,
    'K':  1, 'CA': 2, 'GA': 3, 'GE': 4, 'AS': 5, 'SE': 6, 'BR': 7, 'KR': 8,
    'RB': 1, 'SR': 2, 'IN': 3, 'SN': 4, 'SB': 5, 'TE': 6, 'I':  7, 'XE': 8,
}

ELECTRONEGATIVITIES = {
    # Source: https://en.wikipedia.org/wiki/Electronegativity
    'H':  2.20, 'HE': None,
    'LI': 0.98, 'BE': 1.57, 'B':  2.04, 'C':  2.55, 'N':  3.04, 'O':  3.44, 'F':  3.98, 'NE': None,
    'NA': 0.93, 'MG': 1.31, 'AL': 1.61, 'SI': 1.90, 'P':  2.19, 'S':  2.58, 'CL': 3.16, 'AR': None,
    'K':  0.82, 'CA': 1.00, 'GA': 1.81, 'GE': 2.01, 'AS': 2.18, 'SE': 2.55, 'BR': 2.96, 'KR': 3.00,
    'RB': 0.82, 'SR': 0.95, 'IN': 1.78, 'SN': 1.96, 'SB': 2.05, 'TE': 2.10, 'I':  2.66, 'XE': 2.60,
}

MIN_ABSOLUTE_CHARGE, MAX_ABSOLUTE_CHARGE = 0, 9
MIN_BOND_ORDER, MAX_BOND_ORDER = 1, 3
MAX_NONBONDED_ELECTRONS = 18
ELECTRONS_PER_BOND = 2

MUST_BE_INT = lambda x: round(x)

POSSIBLE_BOND_ORDERS = {
    'S': {1, 2},
    'SE': {1, 2},
    'C': {1, 2, 3},
    'H': {1},
    'O': {1, 2},
    'N': {1, 2, 3},
    'P': {1, 2},
    'CL': {1},
    'BR': {1},
    'F': {1},
    'I': {1, 2},
    'B': {1, 2},
    'SI': {1, 2, 3},
}

POSSIBLE_CHARGES = {
    'S': {0},
    'SE': {0},
    'C': {0},
    'H': {0},
    'O': {-1, 0},
    'N': {-1, 0, +1},
    'P': {0},
    'CL': {0},
    'BR': {0},
    'F': {0},
    'I': {0},
    'B': {-1, 0},
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

NO_CAP = Capping_Strategy([], [], [])

def merge_caps(*list_of_capping_strategies: List[Capping_Strategy]) -> Capping_Strategy:
    def renumber_atom_ids(atom_ids: Tuple[int, int], acc: Capping_Strategy) -> Tuple[int, int]:
        return tuple(
            [
                0 if atom_id == 0 else atom_id + len(acc.new_atoms)
                for atom_id in atom_ids
            ],
        )

    return reduce(
        lambda acc, capping_strategy: Capping_Strategy(
            new_atoms=acc.new_atoms + capping_strategy.new_atoms,
            new_bonds=acc.new_bonds + [renumber_atom_ids(atom_ids, acc) for atom_ids in capping_strategy.new_bonds],
            new_valences=acc.new_valences + capping_strategy.new_valences,
        ),
        list_of_capping_strategies,
        NO_CAP,
    )

H_CAP = Capping_Strategy(['H'], [(0, 1)], [1])
H2_CAP = merge_caps(*[H_CAP] * 2)
H3_CAP = merge_caps(*[H_CAP] * 3)
H4_CAP = merge_caps(*[H_CAP] * 4)
CH2_CAP = Capping_Strategy(['C', 'H', 'H'], [(0, 1), (1, 2), (1, 3)], [3, 1, 1])
H_CH2_CAP = merge_caps(H_CAP, CH2_CAP)
C_H_CAP = Capping_Strategy(['C', 'H'], [(0, 1), (1, 2)], [2, 1])
CH3_CAP = Capping_Strategy(['C', 'H', 'H', 'H'], [(0, 1), (1, 2), (1, 3), (1, 4)], [4, 1, 1, 1])
O_CAP = Capping_Strategy(['O'], [(0, 1)], [1])
O3_CAP = merge_caps(*[O_CAP] * 3)
OH_CAP = Capping_Strategy(['O', 'H'], [(0, 1), (1, 2)], [2, 1])
O_OH_OH_CAP = merge_caps(O_CAP, OH_CAP, OH_CAP)

INDIVIDUAL_CAPPING_OPTIONS = {
    'H1': [NO_CAP, H_CAP],
    'O1': [NO_CAP],
    'O2': [H_CAP, H2_CAP],
    'S1': [NO_CAP, H_CAP],
    'S2': [H_CAP, H2_CAP],
    'S4': [NO_CAP],
    'S6': [NO_CAP],
    'SE1': [NO_CAP, H_CAP],
    'SE2': [H_CAP, H2_CAP],
    'SE4': [NO_CAP],
    'SE6': [NO_CAP],
    'C4': [NO_CAP, H_CAP, H2_CAP, H3_CAP, H4_CAP],
    'C3': [H_CAP, H2_CAP, H_CH2_CAP],
    'C2': [H_CAP],#, C_H_CAP],
    'N1': [NO_CAP],
    'N2': [H_CAP], #CH3_CAP],
    'N3': [H_CAP, H2_CAP],
    'N4': [H3_CAP, H4_CAP],
    'P5': [H_CAP, H2_CAP, O_OH_OH_CAP],
    'CL1': [NO_CAP, H_CAP],
    'BR1': [NO_CAP, H_CAP],
    'F1': [NO_CAP, H_CAP],
    'I1': [NO_CAP, H_CAP],
    'B3': [H_CAP, H2_CAP, H3_CAP],
    'SI4': [H_CAP, H2_CAP, H3_CAP, H4_CAP],
    'SI3': [H_CAP, H2_CAP, H_CH2_CAP],
    'SI2': [H_CAP],
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

def max_valence_for(atom: Atom) -> int:
    return max(FULL_VALENCES[atom.element]) - min(POSSIBLE_CHARGES[atom.element])

def min_valence_for(atom: Atom) -> int:
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

def make_capping_strategy_hashable(capping_strategy: Capping_Strategy) -> Capping_Strategy:
    return Capping_Strategy(
        *map(
            tuple,
            capping_strategy,
        ),
    )

GROUPED_CAPPING_OPTIONS = {
    element: reduce(
        lambda acc, e: acc | set(e),
        [x[1] for x in list(group)],
        set(),
    )
    for (element, group) in
    groupby(
        sorted(
            [(atom_desc, [make_capping_strategy_hashable(c) for c in C]) for (atom_desc, C) in INDIVIDUAL_CAPPING_OPTIONS.items()],
            key=on_letters_of_dict_key,
        ),
        key=on_letters_of_dict_key,
    )
}

def validate_capping_strategies() -> None:
    '''Ensures there exist a capping strategy for any valence of any atom.'''
    for (element, capping_strategies) in GROUPED_CAPPING_OPTIONS.items():
        capping_atom_valences = reduce(lambda acc, e: acc | {new_atom_for_capping_strategy(e)}, capping_strategies, set())
        try:
            assert capping_atom_valences >= set(range(1, min(FULL_VALENCES[element]) + 1)), (
                element,
                capping_atom_valences,
                set(range(1, min(FULL_VALENCES[element]) + 1)),
            )
        except AssertionError:
            from traceback import print_exc
            print_exc()
            print()

validate_capping_strategies()

ALL_CAPPING_OPTIONS = {
        **INDIVIDUAL_CAPPING_OPTIONS,
        **GROUPED_CAPPING_OPTIONS,
}

DOUBLE_BOND_ELEMENTS = [
    element
    for element in ALL_ELEMENTS
    if 2 in POSSIBLE_BOND_ORDERS[element]
]

ALL_POSSIBLE_DOUBLE_BONDS = {
    frozenset([element_1, element_2])
    for (element_1, element_2) in product(DOUBLE_BOND_ELEMENTS, DOUBLE_BOND_ELEMENTS)
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
