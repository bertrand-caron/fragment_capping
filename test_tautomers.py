from sys import stderr
import unittest

from fragment_capping.helpers.types_helpers import Atom
from fragment_capping.helpers.molecule import Molecule, molecule_from_pdb_str

class Test_Capping(unittest.TestCase):
    def test_violuric_acid(self, use_ILP: bool = True) -> None:
        with open('pdbs/violuric_acid.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='violuric_acid',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(1200, 1200),
            graph_kwargs={'include_atom_index': False},
        )
        tautomers = input_molecule.get_all_tautomers(
            #net_charge=0,
            total_number_hydrogens=3,
            enforce_octet_rule=True,
            allow_radicals=False,
        )

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(1200, 1200),
                graph_kwargs={'include_atom_index': False, 'vertex_color_scheme': 'elements', 'vertex_label_template': ''},
            )

        assert len(tautomers) == 15, len(tautomers)

    def test_warfarin(self, use_ILP: bool = True) -> None:
        with open('pdbs/warfarin.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='warfarin',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(1200, 1200),
            graph_kwargs={'include_atom_index': False},
        )
        tautomers = input_molecule.get_all_tautomers(
            net_charge=0,
            total_number_hydrogens=16,
            enforce_octet_rule=True,
            allow_radicals=False,
        )

        assert len(tautomers) == 43, len(tautomers)

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(1200, 1200),
                graph_kwargs={'include_atom_index': False, 'vertex_color_scheme': 'elements', 'vertex_label_template': ''},
            )

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

        assert len(tautomers) == 21, len(tautomers)

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

        assert len(tautomers) == 1, len(tautomers)

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(600, 600),
                graph_kwargs={'include_atom_index': False},
            )

    def test_benzene_fail(self):
        with open('pdbs/benzene.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='benzene_fail',
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
            disallow_triple_bond_in_small_rings=False,
            disallow_allenes_in_small_rings=False,
            disallow_allenes_completely=False,
            lock_phenyl_rings=False,
        )

        assert len(tautomers) == 5, len(tautomers)

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

        assert len(tautomers) == 2, len(tautomers)

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
            disallow_allenes_completely=False,
        )

        assert len(tautomers) == 2, len(tautomers)

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

        assert len(tautomers) == 4, len(tautomers)

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(600, 600),
                graph_kwargs={'include_atom_index': True},
            )

    def test_methane(self):
        with open('pdbs/methane.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='methane',
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

        assert len(tautomers) == 1, len(tautomers)

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(600, 600),
                graph_kwargs={'include_atom_index': True},
            )

    def test_ammonium(self):
        with open('pdbs/ammonium.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='ammonium',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(600, 600),
            graph_kwargs={'include_atom_index': False},
        )
        tautomers = input_molecule.get_all_tautomers(
            net_charge=1,
            total_number_hydrogens=4,
            enforce_octet_rule=True,
            allow_radicals=False,
        )

        assert len(tautomers) == 1, len(tautomers)

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(600, 600),
                graph_kwargs={'include_atom_index': True},
            )

    def test_ammonium_fail(self):
        with open('pdbs/ammonium.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='ammonium',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(600, 600),
            graph_kwargs={'include_atom_index': False},
        )
        with self.assertRaises(Exception):
            tautomers = input_molecule.get_all_tautomers(
                net_charge=1,
                total_number_hydrogens=4,
                maximum_number_hydrogens_per_atom=2,
                enforce_octet_rule=True,
                allow_radicals=False,
            )
        with self.assertRaises(Exception):
            tautomers = input_molecule.get_all_tautomers(
                net_charge=0,
                total_number_hydrogens=4,
                maximum_number_hydrogens_per_atom=3,
                enforce_octet_rule=True,
                allow_radicals=False,
            )

    def test_radicals(self):
        with open('pdbs/nitric_oxide.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='nitric_oxide',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(600, 600),
            graph_kwargs={'include_atom_index': False},
        )

        with self.assertRaises(Exception):
            tautomers = input_molecule.get_all_tautomers(
                net_charge=0,
                total_number_hydrogens=0,
                enforce_octet_rule=True,
                allow_radicals=False,
            )

        tautomers = input_molecule.get_all_tautomers(
            net_charge=0,
            total_number_hydrogens=0,
            enforce_octet_rule=False,
            allow_radicals=True,
        )

        assert len(tautomers) == 1, len(tautomers)

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(600, 600),
                graph_kwargs={'include_atom_index': True},
            )

        tautomer = tautomers[0]

        assert list(tautomer.bond_orders.items()) == [2]

    def test_lewis_octet_exceptions(self):
        with open('pdbs/methylene.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='methylene',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(600, 600),
            graph_kwargs={'include_atom_index': False},
        )

        with self.assertRaises(Exception):
            tautomers = input_molecule.get_all_tautomers(
                net_charge=0,
                total_number_hydrogens=2,
                enforce_octet_rule=True,
                allow_radicals=False,
            )

        tautomers = input_molecule.get_all_tautomers(
            net_charge=0,
            total_number_hydrogens=2,
            enforce_octet_rule=False,
            allow_radicals=False,
        )

        assert len(tautomers) == 1, len(tautomers)

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(600, 600),
                graph_kwargs={'include_atom_index': True},
            )

    def test_skeletal_rearrangement(self):
        with open('pdbs/warfarin.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='warfarin_extended',
            )
        input_molecule.remove_all_hydrogens(mark_all_uncapped=True)
        input_molecule.write_graph(
            'input',
            output_size=(1200, 1200),
            graph_kwargs={'include_atom_index': False},
        )

        with self.assertRaises(TypeError):
            tautomers = input_molecule.get_all_tautomers(
                net_charge=0,
                total_number_hydrogens=16,
                enforce_octet_rule=True,
                allow_radicals=False,
                skeletal_rearrangements=[(1, 2, 3 )],
            )

        tautomers = input_molecule.get_all_tautomers(
            net_charge=0,
            total_number_hydrogens=16,
            enforce_octet_rule=True,
            allow_radicals=False,
            skeletal_rearrangements=[(2, 10), (2, 8)], # Bonds ('C13', 'O1') and ('C13', 'O3')
        )

        assert len(tautomers) == 53, len(tautomers)

        for (n, molecule) in enumerate(tautomers, start=1):
            molecule.name = input_molecule.name
            molecule.write_graph(
                '_tautomer_{0}'.format(n),
                output_size=(1200, 1200),
                graph_kwargs={'include_atom_index': False, 'vertex_color_scheme': 'elements', 'vertex_label_template': ''},
            )

if __name__ == '__main__':
    print('Note: This test is very long (~15 minutes)')
    unittest.main(warnings='ignore', verbosity=2)
