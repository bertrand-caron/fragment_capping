from sys import stderr

from fragment_capping.helpers.types_helpers import Atom
from fragment_capping.helpers.molecule import Molecule

CAPPING_FUNCTION_NAME = 'get_best_capped_molecule'

def example_1() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=3, capped=False, coordinates=None),
            2: Atom(index=2, element='C', valence=3, capped=False, coordinates=None),
        },
        [
            (1, 2),
        ],
        name='example_1',
    )

    uncapped_molecule.write_graph('uncapped_molecule')

    for use_ILP in (True,):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C2H4', capped_molecule.formula(charge=True)

def example_2() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=3, capped=False, coordinates=None),
            2: Atom(index=2, element='O', valence=1, capped=False, coordinates=None),
        },
        [
            (1, 2),
        ],
        name='example_2',
    )

    uncapped_molecule.write_graph('uncapped_molecule')

    for use_ILP in (True,):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'CH2O', capped_molecule.formula(charge=True)

def example_3() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=3, capped=False, coordinates=None),
            2: Atom(index=2, element='C', valence=3, capped=False, coordinates=None),
            3: Atom(index=3, element='C', valence=3, capped=False, coordinates=None),
            4: Atom(index=4, element='C', valence=3, capped=False, coordinates=None),
            5: Atom(index=5, element='C', valence=3, capped=False, coordinates=None),
            6: Atom(index=6, element='C', valence=3, capped=False, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 1),
        ]
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C6H6', capped_molecule.formula(charge=True)

def example_4() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=4, capped=False, coordinates=None),
        },
        [
        ],
        name='example_4',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'CH4', capped_molecule.formula(charge=True)

def example_5() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=2, capped=True, coordinates=None),
            2: Atom(index=2, element='O', valence=1, capped=True, coordinates=None),
            3: Atom(index=3, element='O', valence=1, capped=True, coordinates=None),
        },
        [
            (1, 2),
            (1, 3),
        ],
        name='example_5',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'CO2', capped_molecule.formula(charge=True)

def example_6() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=4, capped=False, coordinates=None),
        },
        [
        ],
        name='example_6',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'CH4', capped_molecule.formula(charge=True)

def example_7() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=3, capped=False, coordinates=None),
            2: Atom(index=2, element='O', valence=2, capped=True, coordinates=None),
            3: Atom(index=3, element='O', valence=1, capped=True, coordinates=None),
            4: Atom(index=4, element='H', valence=1, capped=True, coordinates=None),
        },
        [
            (1, 2),
            (1, 3),
            (2, 4),
        ],
        name='example_7',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'CH2O2', capped_molecule.formula(charge=True)

def example_8() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='N', valence=None, capped=False, coordinates=None),
            2: Atom(index=2, element='H', valence=1, capped=True, coordinates=None),
        },
        [
            (1, 2),
        ],
        name='example_0',
    )

    uncapped_molecule.write_graph('uncapped_molecule')

    for use_ILP in (True,):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'H3N', capped_molecule.formula(charge=True)

def example_9() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='N', valence=None, capped=False, coordinates=None),
            2: Atom(index=2, element='H', valence=1, capped=True, coordinates=None),
        },
        [
            (1, 2),
        ],
        name='example_0',
    )

    uncapped_molecule.write_graph('uncapped_molecule')

    for use_ILP in (True,):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'H3N', capped_molecule.formula(charge=True)


# Source: https://doi.org/10.1016/j.jmgm.2005.12.005

def example_wang_1() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=3, capped=True, coordinates=None),
            2: Atom(index=2, element='C', valence=3, capped=True, coordinates=None),
            3: Atom(index=3, element='C', valence=3, capped=True, coordinates=None),
            4: Atom(index=4, element='C', valence=3, capped=True, coordinates=None),
            5: Atom(index=5, element='C', valence=3, capped=True, coordinates=None),
            6: Atom(index=6, element='C', valence=3, capped=True, coordinates=None),
            7: Atom(index=7, element='H', valence=1, capped=True, coordinates=None),
            8: Atom(index=8, element='H', valence=1, capped=True, coordinates=None),
            9: Atom(index=9, element='H', valence=1, capped=True, coordinates=None),
            10: Atom(index=10, element='H', valence=1, capped=True, coordinates=None),
            11: Atom(index=11, element='H', valence=1, capped=True, coordinates=None),
            12: Atom(index=12, element='O', valence=1, capped=True, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 1),
            (7, 1),
            (8, 2),
            (9, 3),
            (10, 4),
            (11, 5),
            (12, 6),
        ],
        name='example_wang_1',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C6H5O 1-', capped_molecule.formula(charge=True)

def example_wang_2() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='O', valence=1, capped=True, coordinates=None),
            2: Atom(index=2, element='C', valence=2, capped=True, coordinates=None),
            3: Atom(index=3, element='S', valence=3, capped=True, coordinates=None),
            4: Atom(index=4, element='C', valence=4, capped=False, coordinates=None),
            5: Atom(index=5, element='C', valence=4, capped=False, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (3, 4),
            (3, 5),
        ],
        name='example_wang_2',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C3H6OS', capped_molecule.formula(charge=True)

def example_wang_3() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=4, capped=False, coordinates=None),
            2: Atom(index=2, element='C', valence=4, capped=False, coordinates=None),
            3: Atom(index=3, element='C', valence=4, capped=False, coordinates=None),
            4: Atom(index=4, element='P', valence=4, capped=True, coordinates=None),
            5: Atom(index=5, element='O', valence=1, capped=True, coordinates=None),
            6: Atom(index=6, element='O', valence=2, capped=True, coordinates=None),
            7: Atom(index=7, element='O', valence=2, capped=True, coordinates=None),
            8: Atom(index=8, element='C', valence=4, capped=False, coordinates=None),
            9: Atom(index=9, element='C', valence=4, capped=False, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (2, 4),
            (4, 5),
            (4, 6),
            (4, 7),
            (6, 8),
            (7, 9),
        ],
        name='example_wang_3',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C5H13O3P', capped_molecule.formula(charge=True)

def example_wang_4() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='O', valence=1, capped=True, coordinates=None),
            2: Atom(index=2, element='C', valence=3, capped=True, coordinates=None),
            3: Atom(index=3, element='C', valence=3, capped=False, coordinates=None),
            4: Atom(index=4, element='C', valence=3, capped=False, coordinates=None),
            5: Atom(index=5, element='C', valence=3, capped=False, coordinates=None),
            6: Atom(index=6, element='N', valence=2, capped=True, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 2),
        ],
        name='example_wang_4',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C4H3NO', capped_molecule.formula(charge=True)

def example_wang_5() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=4, capped=False, coordinates=None),
            2: Atom(index=2, element='C', valence=4, capped=False, coordinates=None),
            3: Atom(index=3, element='C', valence=4, capped=False, coordinates=None),
            4: Atom(index=4, element='C', valence=2, capped=True, coordinates=None),
            5: Atom(index=5, element='N', valence=2, capped=True, coordinates=None),
            6: Atom(index=6, element='C', valence=4, capped=False, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (2, 4),
            (4, 5),
            (5, 6),
        ],
        name='example_wang_5',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C5H10N 1+', capped_molecule.formula(charge=True)

def example_wang_6() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=3, capped=False, coordinates=None),
            2: Atom(index=2, element='C', valence=3, capped=False, coordinates=None),
            3: Atom(index=3, element='C', valence=3, capped=False, coordinates=None),
            4: Atom(index=4, element='C', valence=3, capped=False, coordinates=None),
            5: Atom(index=5, element='C', valence=3, capped=True, coordinates=None),
            6: Atom(index=6, element='C', valence=3, capped=False, coordinates=None),
            7: Atom(index=7, element='C', valence=3, capped=False, coordinates=None),
            8: Atom(index=8, element='C', valence=3, capped=False, coordinates=None),
            9: Atom(index=9, element='C', valence=3, capped=False, coordinates=None),
            10: Atom(index=10, element='C', valence=3, capped=True, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 1),
            (10, 5),
        ],
        name='example_wang_6',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C10H8', capped_molecule.formula(charge=True)

def example_wang_7() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=4, capped=False, coordinates=None),
            2: Atom(index=2, element='C', valence=2, capped=True, coordinates=None),
            3: Atom(index=3, element='N', valence=2, capped=True, coordinates=None),
            4: Atom(index=4, element='S', valence=1, capped=True, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (3, 4),
        ],
        name='example_wang_7',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C2H3NS', capped_molecule.formula(charge=True)

def example_wang_8() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='N', valence=1, capped=True, coordinates=None),
            2: Atom(index=2, element='N', valence=2, capped=True, coordinates=None),
            3: Atom(index=3, element='C', valence=3, capped=True, coordinates=None),
            4: Atom(index=4, element='C', valence=3, capped=False, coordinates=None),
            5: Atom(index=5, element='P', valence=2, capped=True, coordinates=None),
            6: Atom(index=6, element='N', valence=2, capped=True, coordinates=None),
            7: Atom(index=7, element='C', valence=3, capped=False, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 3),
        ],
        name='example_wang_8',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C3H2N3P', capped_molecule.formula(charge=True)

def example_wang_9() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='H', valence=1, capped=True, coordinates=None),
            2: Atom(index=2, element='N', valence=2, capped=True, coordinates=None),
            3: Atom(index=3, element='C', valence=3, capped=False, coordinates=None),
            4: Atom(index=4, element='C', valence=3, capped=True, coordinates=None),
            5: Atom(index=5, element='C', valence=3, capped=False, coordinates=None),
            6: Atom(index=6, element='P', valence=2, capped=True, coordinates=None),
            7: Atom(index=7, element='H', valence=1, capped=True, coordinates=None),
            8: Atom(index=8, element='N', valence=2, capped=True, coordinates=None),
            9: Atom(index=9, element='N', valence=1, capped=True, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (4, 8),
            (8, 9),
        ],
        name='example_wang_9',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C3H4N3P', capped_molecule.formula(charge=True)

def example_wang_10() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=3, capped=True, coordinates=None),
            2: Atom(index=2, element='C', valence=3, capped=True, coordinates=None),
            3: Atom(index=3, element='C', valence=3, capped=True, coordinates=None),
            4: Atom(index=4, element='C', valence=3, capped=True, coordinates=None),
            5: Atom(index=5, element='C', valence=3, capped=True, coordinates=None),
            6: Atom(index=6, element='C', valence=3, capped=True, coordinates=None),
            7: Atom(index=7, element='O', valence=1, capped=True, coordinates=None),
            8: Atom(index=8, element='H', valence=1, capped=True, coordinates=None),
            9: Atom(index=9, element='H', valence=1, capped=True, coordinates=None),
            10: Atom(index=10, element='O', valence=1, capped=True, coordinates=None),
            11: Atom(index=11, element='H', valence=1, capped=True, coordinates=None),
            12: Atom(index=12, element='H', valence=1, capped=True, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 1),
            (7, 1),
            (8, 2),
            (9, 3),
            (10, 4),
            (11, 5),
            (12, 6),
        ],
        name='example_wang_10',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C6H4O2', capped_molecule.formula(charge=True)

def example_wang_11() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=3, capped=False, coordinates=None),
            2: Atom(index=2, element='N', valence=3, capped=True, coordinates=None),
            3: Atom(index=3, element='H', valence=1, capped=True, coordinates=None),
            4: Atom(index=4, element='C', valence=3, capped=True, coordinates=None),
            5: Atom(index=5, element='C', valence=3, capped=False, coordinates=None),
            6: Atom(index=6, element='N', valence=2, capped=True, coordinates=None),
            7: Atom(index=7, element='S', valence=2, capped=True, coordinates=None),
            8: Atom(index=8, element='C', valence=3, capped=True, coordinates=None),
            9: Atom(index=9, element='N', valence=3, capped=True, coordinates=None),
            10: Atom(index=10, element='H', valence=1, capped=True, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (2, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (9, 1),
            (4, 8),
        ],
        name='example_wang_11',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C4H4N3S 1+', capped_molecule.formula(charge=True)

def example_wang_12() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=3, capped=True, coordinates=None),
            2: Atom(index=2, element='C', valence=3, capped=False, coordinates=None),
            3: Atom(index=3, element='C', valence=3, capped=False, coordinates=None),
            4: Atom(index=4, element='C', valence=3, capped=False, coordinates=None),
            5: Atom(index=5, element='C', valence=3, capped=False, coordinates=None),
            6: Atom(index=6, element='P', valence=4, capped=True, coordinates=None),
            7: Atom(index=7, element='O', valence=2, capped=True, coordinates=None),
            8: Atom(index=8, element='H', valence=1, capped=True, coordinates=None),
            9: Atom(index=9, element='O', valence=2, capped=True, coordinates=None),
            10: Atom(index=10, element='H', valence=1, capped=True, coordinates=None),
            11: Atom(index=11, element='O', valence=1, capped=True, coordinates=None),
        },
        [
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 1),
            (1, 6),
            (6, 7),
            (7, 8),
            (6, 9),
            (9, 10),
            (6, 11),
        ],
        name='example_wang_12',
    )

    uncapped_molecule.write_graph('uncapped_molecule')
    for use_ILP in (True, False):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C5H6O3P 1-', capped_molecule.formula(charge=True)

def example_taxol_core() -> None:
    uncapped_molecule = Molecule(
        {5: Atom(index=5, element="C", valence=3, capped=False, coordinates=(0.746, 3.138, 0.794)), 7: Atom(index=7, element="O", valence=2, capped=True, coordinates=(1.175, 1.853, 0.61)), 8: Atom(index=8, element="C", valence=4, capped=True, coordinates=(1.672, 1.0, 1.641)), 9: Atom(index=9, element="C", valence=4, capped=True, coordinates=(2.696, 1.617, 2.644)), 10: Atom(index=10, element="H", valence=1, capped=True, coordinates=(2.871, 2.705, 2.545)), 11: Atom(index=11, element="H", valence=1, capped=True, coordinates=(3.667, 1.113, 2.766)), 12: Atom(index=12, element="O", valence=2, capped=True, coordinates=(1.805, 1.347, 3.756)), 13: Atom(index=13, element="C", valence=4, capped=True, coordinates=(0.741, 0.845, 2.897)), 14: Atom(index=14, element="H", valence=1, capped=True, coordinates=(-0.111, 1.567, 2.942)), 15: Atom(index=15, element="C", valence=4, capped=False, coordinates=(0.262, -0.558, 3.171)), 18: Atom(index=18, element="C", valence=4, capped=False, coordinates=(1.353, -1.632, 3.131)), 20: Atom(index=20, element="O", valence=2, capped=True, coordinates=(1.998, -1.605, 4.383)), 21: Atom(index=21, element="H", valence=1, capped=True, coordinates=(2.798, -2.116, 4.304)), 22: Atom(index=22, element="C", valence=4, capped=False, coordinates=(2.315, -1.541, 1.879)), 29: Atom(index=29, element="C", valence=4, capped=True, coordinates=(1.951, -0.352, 0.936)), 30: Atom(index=30, element="H", valence=1, capped=True, coordinates=(0.931, -0.608, 0.548)), 31: Atom(index=31, element="C", valence=4, capped=True, coordinates=(2.911, -0.285, -0.293)), 32: Atom(index=32, element="H", valence=1, capped=True, coordinates=(3.702, -1.082, -0.163)), 33: Atom(index=33, element="O", valence=2, capped=False, coordinates=(3.52, 1.006, -0.367)), 101: Atom(index=101, element="C", valence=4, capped=False, coordinates=(0.908, -3.278, 0.308)), 47: Atom(index=47, element="C", valence=4, capped=True, coordinates=(2.191, -0.548, -1.701)), 48: Atom(index=48, element="O", valence=2, capped=True, coordinates=(3.132, -0.275, -2.729)), 49: Atom(index=49, element="H", valence=1, capped=True, coordinates=(3.241, 0.668, -2.768)), 50: Atom(index=50, element="C", valence=4, capped=False, coordinates=(0.937, 0.345, -1.879)), 53: Atom(index=53, element="C", valence=4, capped=False, coordinates=(1.845, -2.08, -1.847)), 54: Atom(index=54, element="C", valence=3, capped=True, coordinates=(0.733, -2.32, -0.843)), 59: Atom(index=59, element="C", valence=3, capped=False, coordinates=(-0.388, -1.56, -0.899))},
        {frozenset({5, 7}), frozenset({13, 15}), frozenset({18, 20}), frozenset({18, 15}), frozenset({50, 47}), frozenset({18, 22}), frozenset({8, 29}), frozenset({8, 9}), frozenset({32, 31}), frozenset({48, 47}), frozenset({9, 10}), frozenset({29, 30}), frozenset({29, 31}), frozenset({9, 11}), frozenset({53, 54}), frozenset({9, 12}), frozenset({20, 21}), frozenset({48, 49}), frozenset({29, 22}), frozenset({33, 31}), frozenset({31, 47}), frozenset({12, 13}), frozenset({13, 14}), frozenset({59, 54}), frozenset({8, 7}), frozenset({53, 47}), frozenset({101, 54}), frozenset({8, 13})},
        name='example_taxol_core',
    )

    uncapped_molecule.write_graph('uncapped_molecule', output_size=(1200, 1200))
    for use_ILP in (True,):
        capped_molecule = getattr(uncapped_molecule, CAPPING_FUNCTION_NAME)(debug=None)
        capped_molecule.write_graph('capped_molecule_with_{0}'.format('ILP' if use_ILP else 'bruteforce'), output_size=(1200, 1200))
        print(capped_molecule.dummy_pdb())

    assert capped_molecule.formula(charge=True) == 'C16H26O5', capped_molecule.formula(charge=True)

ALL_EXAMPLES = [
    example_1,
    example_2,
    example_3,
    example_4,
    example_5,
    example_6,
    example_7,
    example_8,
    example_9,
    example_wang_1,
    example_wang_2,
    example_wang_3,
    example_wang_4,
    example_wang_5,
    example_wang_6,
    example_wang_7,
    example_wang_8,
    example_wang_9,
    example_wang_10,
    example_wang_11,
    example_wang_12,
    example_taxol_core,
]

if __name__ == '__main__':
    for example in ALL_EXAMPLES:
        try:
            example()
        except AssertionError as e:
            print(str(e))
            raise
