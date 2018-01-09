from itertools import groupby
from functools import reduce
from fragment_capping.helpers.molecule import Atom, Molecule

UNKNOWN_VALENCE = None

def atom_index_for_line(line: str) -> str:
    return int(line[6:11].replace(' ', ''))

def element_for_line(line: str) -> str:
    return line[76:78].replace(' ', '')

def atom_name_for_line(line: str) -> str:
    return line[13:16].replace(' ', '')

if __name__ == '__main__':
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

    uncapped_molecule.remove_atoms_with_predicate(
        lambda atom: atom.index in {
            atom_index_for_line(line)
            for line in model_1_atom_lines
            if atom_name_for_line(line) in BACKBONE_HYDROGENS
        }
    )

    from sys import stdout
    print(uncapped_molecule.write_graph('uncapped', graph_kwargs={'include_atom_index': False}))
    uncapped_molecule.get_best_capped_molecule_with_ILP(
        enforce_octet_rule=True,
        #debug=stdout,
    )
    print(uncapped_molecule.write_graph('capped', graph_kwargs={'include_atom_index': False}))
