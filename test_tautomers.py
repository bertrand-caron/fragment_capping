from sys import stderr
import unittest

from fragment_capping.helpers.types_helpers import Atom
from fragment_capping.helpers.molecule import Molecule, molecule_from_pdb_str

class Test_Capping(unittest.TestCase):
    def test_porphyrin(self, use_ILP: bool = True) -> None:
        with open('pdbs/tetraphenylporphyrin.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='porphyrin',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(1200, 1200),
            graph_kwargs={'include_atom_index': False},
        )
        tautomers = input_molecule.get_all_tautomers(
            net_charge=0,
            total_number_hydrogens=30,
            enforce_octet_rule=True,
            allow_radicals=False,
        )

        assert len(tautomers) == 46, tautomers

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(1200, 1200),
                graph_kwargs={'include_atom_index': False},
            )

    def test_methylimidazole(self, use_ILP: bool = True) -> None:
        with open('pdbs/methylimidazole.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='methylimidazole',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(600, 600),
            graph_kwargs={'include_atom_index': False},
        )
        tautomers = input_molecule.get_all_tautomers(
            net_charge=0,
            total_number_hydrogens=6,
            enforce_octet_rule=True,
            allow_radicals=False,
        )
        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(600, 600),
                graph_kwargs={'include_atom_index': False},
            )

    def test_benzene(self):
        with open('pdbs/benzene.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='benzene',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(600, 600),
            graph_kwargs={'include_atom_index': True},
        )
        tautomers = input_molecule.get_all_tautomers(
            net_charge=0,
            total_number_hydrogens=6,
            enforce_octet_rule=True,
            allow_radicals=False,
        )

        assert len(tautomers) == 1, tautomers

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(600, 600),
                graph_kwargs={'include_atom_index': False},
            )

    def test_ethanal(self):
        with open('pdbs/ethanal.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='ethanal',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(600, 600),
            graph_kwargs={'include_atom_index': False},
        )
        tautomers = input_molecule.get_all_tautomers(
            net_charge=0,
            total_number_hydrogens=4,
            enforce_octet_rule=True,
            allow_radicals=False,
        )

        assert len(tautomers) == 2, tautomers

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(600, 600),
                graph_kwargs={'include_atom_index': True},
            )

    def test_propadiene(self):
        with open('pdbs/propadiene.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='propadiene',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(600, 600),
            graph_kwargs={'include_atom_index': False},
        )
        tautomers = input_molecule.get_all_tautomers(
            net_charge=0,
            total_number_hydrogens=4,
            enforce_octet_rule=True,
            allow_radicals=False,
        )

        assert len(tautomers) == 2, tautomers

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(600, 600),
                graph_kwargs={'include_atom_index': True},
            )

    def test_cyclooctyne(self):
        with open('pdbs/cyclooctyne.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='cyclooctyne',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(600, 600),
            graph_kwargs={'include_atom_index': False},
        )
        tautomers = input_molecule.get_all_tautomers(
            net_charge=0,
            total_number_hydrogens=12,
            enforce_octet_rule=True,
            allow_radicals=False,
        )

        assert len(tautomers) == 5, tautomers

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(600, 600),
                graph_kwargs={'include_atom_index': True},
            )

if __name__ == '__main__':
    print('Note: This test is very long (~15 minutes)')
    unittest.main(warnings='ignore')
