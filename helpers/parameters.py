from typing import NamedTuple, Any

FULL_VALENCES = {
    'C': 4,
    'N': 3,
    'O': 2,
    'H': 1,
    'S': 2,
    'P': 5,
}

POSSIBLE_BOND_ORDERS = {
    'S': (1, 2,),
    'C': (1, 2,),
    'H': (1,),
    'O': (1, 2,),
    'N': (1, 2,),
    'P': (1, 2,),
}

POSSIBLE_CHARGES = {
    'S': (0,),
    'C': (0,),
    'H': (0,),
    'O': (0, -1,),
    'N': (0, +1,),
    'P': (0,),
}

Capping_Strategy = NamedTuple('Capping_Strategy', [('new_atoms', Any), ('new_bonds', Any), ('new_valences', Any)])

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
}
