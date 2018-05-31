from sys import stderr
from traceback import format_exc
import unittest

from fragment_capping.helpers.types_helpers import Atom
from fragment_capping.helpers.molecule import Molecule

CAN_FAIL, CAN_NOT_FAIL = True, False

USE_OCTET_RULE, DO_NOT_USE_OCTET_RULE = True, False

ALL_EXAMPLES = {
    Molecule(
        [
            Atom(index=1, element='C', valence=4, capped=True, coordinates=None),
            Atom(index=2, element='C', valence=3, capped=True, coordinates=None),
            Atom(index=3, element='H', valence=None, capped=False, coordinates=None),
            Atom(index=4, element='H', valence=None, capped=False, coordinates=None),
            Atom(index=5, element='H', valence=None, capped=False, coordinates=None),
            Atom(index=6, element='N', valence=None, capped=False, coordinates=None),
            Atom(index=7, element='C', valence=None, capped=False, coordinates=None),
        ],
        [
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5),
            (2, 6),
            (2, 7),
        ],
        name='H,H,H_C_C_N,C',
    ): (0, 0, 'C3H7N', USE_OCTET_RULE, CAN_NOT_FAIL),
    Molecule(
        [
            Atom(index=1, element='C', valence=4, capped=True, coordinates=None),
            Atom(index=2, element='C', valence=4, capped=True, coordinates=None),
            Atom(index=3, element='H', valence=None, capped=False, coordinates=None),
            Atom(index=4, element='H', valence=None, capped=False, coordinates=None),
            Atom(index=5, element='H', valence=None, capped=False, coordinates=None),
            Atom(index=6, element='H', valence=None, capped=False, coordinates=None),
            Atom(index=7, element='H', valence=None, capped=False, coordinates=None),
            Atom(index=8, element='H', valence=None, capped=False, coordinates=None),
        ],
        [
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5),
            (2, 6),
            (2, 7),
            (2, 8),
        ],
        name='H,H,H_C_C_H,H,H',
    ): (0, 0, 'C2H6', USE_OCTET_RULE, CAN_NOT_FAIL),
    Molecule(
        [
            Atom(index=1, element='C', valence=4, capped=True, coordinates=None),
            Atom(index=2, element='C', valence=4, capped=True, coordinates=None),
            Atom(index=3, element='C', valence=None, capped=False, coordinates=None),
            Atom(index=4, element='C', valence=None, capped=False, coordinates=None),
            Atom(index=5, element='C', valence=None, capped=False, coordinates=None),
            Atom(index=6, element='C', valence=None, capped=False, coordinates=None),
            Atom(index=7, element='C', valence=None, capped=False, coordinates=None),
            Atom(index=8, element='C', valence=None, capped=False, coordinates=None),
        ],
        [
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5),
            (2, 6),
            (2, 7),
            (2, 8),
        ],
        name='C,C,C_C_C_C,C,C',
    ): (0, 0, 'C8H18', USE_OCTET_RULE, CAN_NOT_FAIL),
    Molecule(
        [
            Atom(index=1, element='C', valence=3, capped=True, coordinates=None),
            Atom(index=2, element='C', valence=3, capped=True, coordinates=None),
            Atom(index=3, element='C', valence=None, capped=False, coordinates=None),
            Atom(index=4, element='C', valence=None, capped=False, coordinates=None),
            Atom(index=5, element='C', valence=None, capped=False, coordinates=None),
            Atom(index=6, element='C', valence=None, capped=False, coordinates=None),
        ],
        [
            (1, 2),
            (1, 3),
            (1, 4),
            (2, 5),
            (2, 6),
        ],
        name='C,C_C_C_C,C',
    ): (0, 0, 'C6H12', USE_OCTET_RULE, CAN_FAIL),
    Molecule(
        [
            Atom(index=1, element='O', valence=2, capped=False, coordinates=None),
            Atom(index=2, element='N', valence=2, capped=False, coordinates=None),
            Atom(index=3, element='C', valence=None, capped=True, coordinates=None),
            Atom(index=4, element='C', valence=None, capped=True, coordinates=None),
            Atom(index=5, element='N', valence=None, capped=False, coordinates=None),
            Atom(index=6, element='C', valence=3, capped=False, coordinates=None),
        ],
        [
            (1, 3),
            (2, 3),
            (3, 4),
            (4, 5),
            (4, 6),
        ],
        name='O,N_C_C_N,C',
    ): (0, 0, 'C3H6N2O', USE_OCTET_RULE, CAN_NOT_FAIL),
    Molecule(
        [
            Atom(index=1, element='O', valence=None, capped=False, coordinates=None),
            Atom(index=2, element='N', valence=None, capped=False, coordinates=None),
            Atom(index=3, element='C', valence=3, capped=True, coordinates=None),
            Atom(index=4, element='C', valence=3, capped=True, coordinates=None),
            Atom(index=5, element='H', valence=None, capped=False, coordinates=None),
            Atom(index=6, element='H', valence=None, capped=False, coordinates=None),
        ],
        [
            (1, 3),
            (2, 3),
            (3, 4),
            (4, 5),
            (4, 6),
        ],
        name='O,N_C_C_H,H',
    ): (0, 0, 'C2H5NO', USE_OCTET_RULE, CAN_NOT_FAIL),
    Molecule(
        [
            Atom(index=1, element='O', valence=None, capped=False, coordinates=None),
            Atom(index=2, element='N', valence=None, capped=False, coordinates=None),
            Atom(index=3, element='C', valence=3, capped=True, coordinates=None),
            Atom(index=4, element='C', valence=3, capped=True, coordinates=None),
            Atom(index=5, element='C', valence=None, capped=False, coordinates=None),
            Atom(index=6, element='H', valence=None, capped=False, coordinates=None),
        ],
        [
            (1, 3),
            (2, 3),
            (3, 4),
            (4, 5),
            (4, 6),
        ],
        name='O,N_C_C_C,H',
    ): (0, 0, 'C3H7NO', USE_OCTET_RULE, CAN_NOT_FAIL),
}

class Test_Dihedral_Fragment_Capping(unittest.TestCase):
    pass

for (molecule, (expected_netcharge, expected_abs_netcharge, expected_formula, use_octet_rule, can_fail)) in sorted(ALL_EXAMPLES.items()):
    def dynamic_test(self, molecule=molecule, expected_netcharge=expected_netcharge, expected_abs_netcharge=expected_abs_netcharge, expected_formula=expected_formula, use_octet_rule=use_octet_rule, can_fail=can_fail):
        molecule.write_graph('uncapped')
        print(molecule.name + '...', end='')
        molecule.get_best_capped_molecule_with_ILP(enforce_octet_rule=use_octet_rule)
        molecule.write_graph('capped')
        print('OK')

        try:
            assert molecule.netcharge() == expected_netcharge, (molecule.netcharge(), molecule.formal_charges)
            assert molecule.net_abs_charges() == expected_abs_netcharge, (molecule.net_abs_charges(), molecule.formal_charges)
            assert molecule.formula() == expected_formula, molecule.formula()
        except AssertionError as e:
            print(molecule.name)
            if can_fail:
                print(str(e))
                print(format_exc())
            else:
                raise

    setattr(
        Test_Dihedral_Fragment_Capping,
        'test_' + molecule.name,
        dynamic_test,
    )

if __name__ == '__main__':
    unittest.main()
