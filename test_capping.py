from sys import stderr

from fragment_capping.helpers.types_helpers import Atom
from fragment_capping.helpers.molecule import Molecule

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
    ): (0, 'C3H7N'),
    Molecule(
        [
            Atom(index=1, element='C', valence=4, capped=True, coordinates=None),
            Atom(index=2, element='C', valence=3, capped=True, coordinates=None),
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
    ): (0, 'C2H6'),
}

if __name__ == '__main__':
    for (molecule, (expected_netcharge, expected_formula)) in ALL_EXAMPLES.items():
        molecule.write_graph('uncapped')
        molecule.get_best_capped_molecule_with_ILP(enforce_octet_rule=True)
        molecule.write_graph('capped')

        try:
            assert molecule.netcharge() == expected_netcharge, (molecule.netcharge(), molecule.formal_charges)
            assert molecule.formula() == expected_formula, molecule.formula()
        except AssertionError:
            print(molecule.name)
            raise
