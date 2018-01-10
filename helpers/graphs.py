from typing import List, Tuple

from networkx import Graph
from networkx.algorithms import is_isomorphic

def lewis_graph(molecule: 'Molecule') -> Graph:
    '''
    Return the Lewis graph of a molecule (Using networkx' Graph() class).
    '''
    G = Graph()

    for atom in molecule.atoms.values():
        G.add_node(atom.index, element=atom.element, non_bonded_electrons=molecule.non_bonded_electrons[atom.index])

    for bond in molecule.bonds:
        G.add_edge(*bond, order=molecule.bond_orders[bond])

    return G

def unique_molecules(molecules: List['Molecule']) -> List['Molecule']:
    '''
    Return list of one-by-one non-isomorphic graphs (including bond order, chemical elements and number of lone pairs)
    '''
    lewis_graphs = [
        lewis_graph(molecule)
        for molecule in molecules
    ]

    unique_molecules: List[Tuple['Molecules', Graph]] = []

    for (molecule_1, graph_1) in zip(molecules, lewis_graphs):
        print('Assessing {0} (unique_molecules={1})'.format(molecule_1.name, [m.name for (m, _) in unique_molecules]))
        if all(
            not is_isomorphic(
                graph_1,
                graph_2,
                node_match=lambda node_1, node_2: node_1['element'] == node_2['element'] and node_1['non_bonded_electrons'] == node_2['non_bonded_electrons'],
                edge_match=lambda edge_1, edge_2: edge_1['order'] == edge_2['order'],
            )
            for (_, graph_2) in unique_molecules
        ):
            print('UNIQUE: {0}'.format(molecule_1.name))
            unique_molecules.append(
                (molecule_1, graph_1),
            )
        else:
            print('NOT UNIQUE: {0}'.format(molecule_1.name))

    return [
        molecule
        for (molecule, _) in unique_molecules
    ]
