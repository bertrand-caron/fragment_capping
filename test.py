from sys import stderr

from fragment_capping.helpers.types_helpers import Atom
from fragment_capping.helpers.molecule import Molecule

def example_1() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=3, capped=False, coordinates=None),
            2: Atom(index=2, element='C', valence=3, capped=False, coordinates=None),
        },
        [
            (1, 2),
        ],
    )

    uncapped_molecule.write_graph('uncapped_molecule')

    capped_molecule = uncapped_molecule.get_best_capped_molecule(debug=None)

    capped_molecule.write_graph('capped_molecule')

    print(capped_molecule.dummy_pdb())

def example_2() -> None:
    uncapped_molecule = Molecule(
        {
            1: Atom(index=1, element='C', valence=3, capped=False, coordinates=None),
            2: Atom(index=2, element='O', valence=1, capped=False, coordinates=None),
        },
        [
            (1, 2),
        ]
    )

    uncapped_molecule.write_graph('uncapped_molecule')

    capped_molecule = uncapped_molecule.get_best_capped_molecule()

    capped_molecule.write_graph('capped_molecule')

    print(capped_molecule.dummy_pdb())

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

    capped_molecule = uncapped_molecule.get_best_capped_molecule()

    capped_molecule.write_graph('capped_molecule')

    print(capped_molecule.dummy_pdb())

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
    capped_molecule = uncapped_molecule.get_best_capped_molecule()
    capped_molecule.write_graph('capped_molecule')

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
    capped_molecule = uncapped_molecule.get_best_capped_molecule()
    capped_molecule.write_graph('capped_molecule')

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
    capped_molecule = uncapped_molecule.get_best_capped_molecule()
    capped_molecule.write_graph('capped_molecule')

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
    capped_molecule = uncapped_molecule.get_best_capped_molecule()
    capped_molecule.write_graph('capped_molecule')

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
    capped_molecule = uncapped_molecule.get_best_capped_molecule()
    capped_molecule.write_graph('capped_molecule')

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
    capped_molecule = uncapped_molecule.get_best_capped_molecule()
    capped_molecule.write_graph('capped_molecule')

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
    capped_molecule = uncapped_molecule.get_best_capped_molecule()
    capped_molecule.write_graph('capped_molecule')

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
    capped_molecule = uncapped_molecule.get_best_capped_molecule()
    capped_molecule.write_graph('capped_molecule')

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
    capped_molecule = uncapped_molecule.get_best_capped_molecule()
    capped_molecule.write_graph('capped_molecule')

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
    capped_molecule = uncapped_molecule.get_best_capped_molecule()
    capped_molecule.write_graph('capped_molecule')

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
    capped_molecule = uncapped_molecule.get_best_capped_molecule()
    capped_molecule.write_graph('capped_molecule')

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
    capped_molecule = uncapped_molecule.get_best_capped_molecule()
    capped_molecule.write_graph('capped_molecule')

ALL_EXAMPLES = [
    example_1,
    example_2,
    example_3,
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
]

if __name__ == '__main__':
    for example in ALL_EXAMPLES:
        example()
