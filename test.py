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
        ]
    )

    uncapped_molecule.write_graph('uncapped_molecule')

    capped_molecule = uncapped_molecule.get_best_capped_molecule()

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

if __name__ == '__main__':
    for example in [example_1, example_2, example_3]:
        example()
