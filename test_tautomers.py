from sys import stderr

from fragment_capping.helpers.types_helpers import Atom
from fragment_capping.helpers.molecule import Molecule, molecule_from_pdb_str

def example_porphyrin(use_ILP: bool = True) -> None:
    with open('pdbs/tetraphenylporphyrin.pdb') as fh:
        molecule = molecule_from_pdb_str(
            fh.read(),
            name='porphyrin',
        )
    molecule.remove_all_hydrogens()
    molecule.write_graph(
        'input',
        output_size=(1200, 1200),
        graph_kwargs={'include_atom_index': False},
    )
    molecule.get_all_tautomers(
        net_charge=0,
        total_number_hydrogens=30,
        enforce_octet_rule=True,
        allow_radicals=False,
    )
    molecule.write_graph(
        'capped',
        output_size=(1200, 1200),
        graph_kwargs={'include_atom_index': False},
    )


ALL_EXAMPLES = [
    example_porphyrin,
]

def main() -> None:
    for example in ALL_EXAMPLES:
        try:
            example()
        except AssertionError as e:
            print(str(e))
            raise

if __name__ == '__main__':
    main()
