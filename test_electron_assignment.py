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

if __name__ == '__main__':
    unittest.main(warnings='ignore', verbosity=2)
