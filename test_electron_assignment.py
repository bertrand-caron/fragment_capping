from sys import stderr
import unittest

from fragment_capping.helpers.types_helpers import Atom
from fragment_capping.helpers.molecule import Molecule, molecule_from_pdb_str

class Test_Capping(unittest.TestCase):
    def test_carbon_dioxide(self) -> None:
        with open('pdbs/carbon_dioxide.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='carbon_dioxide',
            )

        input_molecule.assign_bond_orders_and_charges_with_ILP(
            enforce_octet_rule=True,
        )

        input_molecule.write_graph(
            '',
            output_size=(400, 400),
            graph_kwargs={'include_atom_index': True},
        )

    def test_cyano(self) -> None:
        with open('pdbs/cyano.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='cyano',
            )

        input_molecule.assign_bond_orders_and_charges_with_ILP(
            enforce_octet_rule=True,
        )

        input_molecule.write_graph(
            '',
            output_size=(400, 400),
            graph_kwargs={'include_atom_index': True},
        )

    def test_guanidinium(self) -> None:
        with open('pdbs/guanidinium.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='guanidinium_fail',
            )

        input_molecule.assign_bond_orders_and_charges_with_ILP(
            enforce_octet_rule=True,
            debug=stderr,
        )

        input_molecule.write_graph(
            '',
            output_size=(400, 400),
            graph_kwargs={
                'include_atom_index': False,
                'vertex_color_scheme': 'elements',
                'vertex_label_template': '{charge_str}',
            },
        )

    def test_guanidinium_with_charge(self) -> None:
        with open('pdbs/guanidinium.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='guanidinium_charge',
            )

        input_molecule.assign_bond_orders_and_charges_with_ILP(
            net_charge=+1,
            enforce_octet_rule=True,
            debug=stderr,
        )

        input_molecule.write_graph(
            '',
            output_size=(400, 400),
            graph_kwargs={
                'include_atom_index': False,
                'vertex_color_scheme': 'elements',
                'vertex_label_template': '{charge_str}',
            },
        )

    def test_guanidinium_with_total_electrons(self) -> None:
        with open('pdbs/guanidinium.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='guanidinium_electrons',
            )

        input_molecule.assign_bond_orders_and_charges_with_ILP(
            total_electrons=24,
            enforce_octet_rule=True,
            debug=stderr,
        )

        input_molecule.write_graph(
            '',
            output_size=(400, 400),
            graph_kwargs={
                'include_atom_index': False,
                'vertex_color_scheme': 'elements',
                'vertex_label_template': '{charge_str}',
            },
        )

if __name__ == '__main__':
    unittest.main(warnings='ignore', verbosity=2)
