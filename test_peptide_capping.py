from itertools import groupby
from functools import reduce
from datetime import datetime
from copy import deepcopy
import unittest
from sys import stdout
from os.path import basename
from typing import Optional

from fragment_capping.helpers.molecule import Atom, Molecule, molecule_from_pdb_str
from fragment_capping.helpers.compare import compare_capped_molecules

UNKNOWN_VALENCE = None

def atom_index_for_line(line: str) -> str:
    return int(line[6:11].replace(' ', ''))

def element_for_line(line: str) -> str:
    return line[76:78].replace(' ', '')

def atom_name_for_line(line: str) -> str:
    return line[13:16].replace(' ', '')

def test_capping_peptide(peptide_pdb_filepath: str, peptide_net_charge: int, peptide_number_hydrogens: Optional[int] = None) -> None:
    peptide_name, _ = basename(peptide_pdb_filepath).split('.')

    with open(peptide_pdb_filepath) as fh:
        input_molecule = molecule_from_pdb_str(
            fh.read(),
            name=peptide_name,
            net_charge=peptide_net_charge,
        )

    molecule_copy = deepcopy(input_molecule)
    molecule_copy.assign_bond_orders_and_charges_with_ILP()
    old_formula, old_charge = molecule_copy.formula(charge=True), molecule_copy.netcharge()

    molecule_copy.write_graph(
        'input',
        graph_kwargs={
            'include_atom_index': False,
            'vertex_color_scheme': 'elements',
            'vertex_label_template': '{charge_str}',
        }
    )

    input_molecule.remove_atoms_with_predicate(
        lambda atom: atom.element == 'H',
        reset_valences=False,
        renumber_atoms=False,
    )
    kept_atom_ids = [atom.index for atom in input_molecule.atoms.values()]

    def enforce_valence_aromatic_rings() -> None:
        molecule_copy.assign_aromatic_bonds()
        aromatic_atom_ids = {
            atom_id
            for atom_id in reduce(
                lambda acc, e: acc | set(e),
                molecule_copy.aromatic_rings,
                set(),
            )
        }

        input_molecule.atoms = {
            atom_id: atom if (not atom.capped and atom_id in aromatic_atom_ids) else atom._replace(valence=None)
            for (atom_id, atom) in input_molecule.atoms.items()
        }

    enforce_valence_aromatic_rings()

    now = datetime.now()
    new_molecule = input_molecule.get_best_capped_molecule_with_ILP(
        enforce_octet_rule=True,
        net_charge=old_charge,
        number_hydrogens=peptide_number_hydrogens,
    )
    print((datetime.now() - now).total_seconds())
    new_formula = new_molecule.formula(charge=True)

    print(
        new_molecule.write_graph(
            'capped',
            graph_kwargs={
                'include_atom_index': False,
                'vertex_color_scheme': 'elements',
                'vertex_label_template': '{charge_str}',
            }
        ,)
    )

    with open(peptide_pdb_filepath.replace('.pdb', '_capped.pdb'), 'wt') as fh:
        fh.write(new_molecule.dummy_pdb())

    print(compare_capped_molecules(molecule_copy, new_molecule, kept_atom_ids))

    assert old_formula == new_formula, (
        old_formula,
        new_formula,
        compare_capped_molecules(molecule_copy, new_molecule, kept_atom_ids)
    )

class Test_Peptide_Capping(unittest.TestCase):
    def test_capping_AAA(self):
        test_capping_peptide('pdbs/AAA.pdb', 0, None)

    def test_capping_ARNDCEQGHILKMFPSTWYV(self):
        test_capping_peptide('pdbs/ARNDCEQGHILKMFPSTWYV.pdb', 0, None)

    def test_capping_2OVN_all(self) -> None:
        test_capping_peptide('pdbs/2OVN_with_connects.pdb', -4, 144)

if __name__ == '__main__':
    unittest.main(warnings='ignore', verbosity=2)
