from itertools import groupby
from functools import reduce
from datetime import datetime
from copy import deepcopy
import unittest
from sys import stdout

from fragment_capping.helpers.molecule import Atom, Molecule, molecule_from_pdb_str

UNKNOWN_VALENCE = None

def atom_index_for_line(line: str) -> str:
    return int(line[6:11].replace(' ', ''))

def element_for_line(line: str) -> str:
    return line[76:78].replace(' ', '')

def atom_name_for_line(line: str) -> str:
    return line[13:16].replace(' ', '')

class Test_Peptide_Capping(unittest.TestCase):
    def test_capping_all(self) -> None:
        with open('pdbs/2OVN_with_connects.pdb') as fh:
            input_molecule = molecule_from_pdb_str(
                fh.read(),
                name='2OVN_complete',
            )

        old_formula, old_charge = 'C89H144N26O27 4-', -4

        input_molecule.remove_atoms_with_predicate(
            lambda atom: atom.element == 'H',
        )

        now = datetime.now()
        new_molecule = input_molecule.get_best_capped_molecule_with_ILP(
            enforce_octet_rule=True,
            net_charge=old_charge,
            #debug=stdout,
        )
        print((datetime.now() - now).total_seconds())
        new_formula = new_molecule.formula(charge=True)

        assert old_formula == new_formula, (old_formula, new_formula)

    def test_capping_backbone(self) -> None:
        import faulthandler, signal
        faulthandler.register(signal.SIGTERM)

        model_1_atom_lines = []
        model_1_connect_lines = []

        BACKBONE_CARBONS = ['CA', 'N', 'CB']
        BACKBONE_HYDROGENS = ['HA', 'H', 'HB2', 'HB3']

        with open('pdbs/2OVN_with_connects.pdb') as fh:
            for line in fh.read().splitlines():
                if line.startswith('HETATM') or line.startswith('ATOM '):
                    model_1_atom_lines.append(line)
                if line.startswith('ENDMDL'):
                    break
                if line.startswith('CONECT'):
                    model_1_connect_lines.append(line)

        all_bonds = [
            (int(bond_index_1), int(bond_index_2))
            for (_, bond_index_1, *bond_indices) in map(lambda line: line.split(), model_1_connect_lines)
            for bond_index_2 in bond_indices
        ]

        uncapped_molecule = Molecule(
            [
                Atom(
                    index=atom_index_for_line(atom_line),
                    element=element_for_line(atom_line),
                    valence=None,
                    capped=False if atom_name_for_line(atom_line) in BACKBONE_CARBONS else True,
                    coordinates=None,
                )
                for atom_line in model_1_atom_lines
            ],
            all_bonds,
            name='2OVN',
        )

        molecule_copy = deepcopy(uncapped_molecule)
        molecule_copy.assign_bond_orders_and_charges_with_ILP()
        old_formula, old_charge = molecule_copy.formula(charge=True), molecule_copy.netcharge()

        uncapped_molecule.remove_atoms_with_predicate(
            lambda atom: atom.index in {
                atom_index_for_line(line)
                for line in model_1_atom_lines
                if atom_name_for_line(line) in BACKBONE_HYDROGENS
            }
        )

        print(uncapped_molecule.write_graph('uncapped', graph_kwargs={'include_atom_index': False}))
        now = datetime.now()
        new_molecule = uncapped_molecule.get_best_capped_molecule_with_ILP(
            enforce_octet_rule=True,
            net_charge=old_charge,
            #debug=stdout,
        )
        print((datetime.now() - now).total_seconds())
        print(uncapped_molecule.write_graph('capped', graph_kwargs={'include_atom_index': False}))

        new_formula = new_molecule.formula(charge=True)

        assert old_formula == new_formula, (old_formula, new_formula)

if __name__ == '__main__':
    unittest.main(warnings='ignore', verbosity=2)
