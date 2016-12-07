from typing import Any, List, Optional

from fragment_capping.helpers.molecule import Uncapped_Molecule, Molecule
from fragment_capping.helpers.types_helpers import Fragment

from fragment_dihedrals.fragment_dihedral import element_valence_for_atom

DEBUG = False

def best_capped_molecule_for_dihedral_fragment(fragment_str: Fragment) -> Molecule:
    molecule = uncapped_molecule_for_dihedral_fragment(fragment_str).get_best_capped_molecule()

    if DEBUG:
        print(molecule)

    return molecule

def uncapped_molecule_for_dihedral_fragment(dihedral_fragment: Fragment) -> Uncapped_Molecule:
    if dihedral_fragment.count('|') == 3:
        neighbours_1, atom_2, atom_3, neighbours_4 = dihedral_fragment.split('|')
        cycles = []
        neighbours_1, neighbours_4 = neighbours_1.split(','), neighbours_4.split(',')
    elif dihedral_fragment.count('|') == 4:
        neighbours_1, atom_2, atom_3, neighbours_4, cycles = dihedral_fragment.split('|')
        neighbours_1, neighbours_4, cycles = neighbours_1.split(','), neighbours_4.split(','), cycles.split(',')
    else:
        raise Exception('Invalid dihedral_fragment: "{0}"'.format(dihedral_fragment))

    ids = [n for (n, _) in enumerate(neighbours_1 + [atom_2, atom_3] + neighbours_4)]

    neighbours_id_1, atom_id_2, atom_id_3, neighbours_id_4 = ids[0:len(neighbours_1)], ids[len(neighbours_1)], ids[len(neighbours_1) + 1], ids[len(neighbours_1) + 2:]
    CENTRAL_BOND = (atom_id_2, atom_id_3)

    elements = dict(
        list(zip(
            ids,
            [element_valence_for_atom(neighbour)[0] for neighbour in neighbours_1] + [atom_2, atom_3] + [element_valence_for_atom(neighbour)[0] for neighbour in neighbours_4],
        )),
    )

    valences = dict(
        list(zip(
            ids,
            [element_valence_for_atom(neighbour)[1] for neighbour in neighbours_1] + [len(neighbours_1) + 1, len(neighbours_4) + 1] + [element_valence_for_atom(neighbour)[1] for neighbour in neighbours_4],
        )),
    )

    bonds = (
        [
            (neighbour_id, atom_id_2)
            for neighbour_id in neighbours_id_1
        ]
        +
        [CENTRAL_BOND]
        +
        [
            (atom_id_3, neighbour_id)
            for neighbour_id in neighbours_id_4
        ]
    )

    molecule = Molecule(
        dict(
            list(zip(
                ids,
                [
                    {
                        'valence': valences[atom_id],
                        'element': elements[atom_id],
                        'index':atom_id,
                        'capped': (atom_id not in (neighbours_id_1 + neighbours_id_4)),
                    }
                    for atom_id in ids],

            ))
        ),
        bonds,
        name=dihedral_fragment.replace('|', '_'),
    )

    if DEBUG:
        print(molecule)

    for (i, n, j) in map(lambda cycle: map(int, cycle), cycles):
        i_id, j_id = neighbours_id_1[i], neighbours_id_4[j]
        if n == 0:
            # i and j are actually the same atoms
            del molecule.atoms[j_id]
            replace_j_by_i = lambda x: i_id if x == j_id else x
            molecule.bonds = [
                list(map(replace_j_by_i, bond))
                for bond in molecule.bonds
            ]
        else:
            NEW_ATOM_ID = -1
            NEW_ATOM = {
                'valence': NO_VALENCE,
                'element': 'C',
                'index': NEW_ATOM_ID, # This will get overwritten by Molecule.add_atom
                'capped': False,
            }
            atom_chain_id = [i_id] + [molecule.add_atom(NEW_ATOM) for i in range(n - 1)] + [j_id]
            new_bonds = zip(atom_chain_id[:-1], atom_chain_id[1:])
            molecule.bonds += new_bonds

    return molecule


