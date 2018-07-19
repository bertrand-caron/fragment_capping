from typing import List, Tuple, Callable, Sequence

from networkx import Graph
from networkx.algorithms import is_isomorphic

def lewis_graph(molecule: 'Molecule', use_non_bonded_electrons: bool = True) -> Graph:
    '''
    Return the Lewis graph of a molecule (Using networkx' Graph() class).
    '''
    G = Graph()

    for atom in molecule.atoms.values():
        G.add_node(
            atom.index,
            element=atom.element,
            non_bonded_electrons=molecule.non_bonded_electrons[atom.index] if use_non_bonded_electrons else None,
        )

    for bond in molecule.bonds:
        G.add_edge(*bond, order=molecule.bond_orders[bond])

    return G

def are_nodes_equivalent(node_1: 'Node', node_2: 'Node') -> bool:
    return node_1['element'] == node_2['element'] and node_1['non_bonded_electrons'] == node_2['non_bonded_electrons']

def are_edges_equivalent(edge_1: 'Edge', edge_2: 'Edge') -> bool:
    return edge_1['order'] == edge_2['order']

def are_graphs_isomorphic(
    graphs: Sequence[Graph],
    node_match: Callable[['Node', 'Node'], bool] = are_nodes_equivalent,
    edge_match: Callable[['Edge', 'Edge'], bool] = are_edges_equivalent,
) -> bool:
    return is_isomorphic(
        *graphs,
        node_match=node_match,
        edge_match=edge_match,
    )

def unique_molecules(molecules: List['Molecule'], debug: bool = False) -> List['Molecule']:
    '''
    Return list of one-by-one non-isomorphic graphs (including bond order, chemical elements and number of lone pairs)
    '''
    lewis_graphs = [
        lewis_graph(molecule)
        for molecule in molecules
    ]

    unique_molecules: List[Tuple['Molecules', Graph]] = []

    print_if_debug = lambda *args, **kwargs: print(*args, **kwargs) if debug else None

    for (molecule_1, graph_1) in zip(molecules, lewis_graphs):
        print_if_debug('Assessing {0} (unique_molecules={1})'.format(molecule_1.name, [m.name for (m, _) in unique_molecules]))
        if all(
            not are_graphs_isomorphic(
                (graph_1, graph_2),
                node_match=are_nodes_equivalent,
                edge_match=are_edges_equivalent,
            )
            for (_, graph_2) in unique_molecules
        ):
            print_if_debug('UNIQUE: {0}'.format(molecule_1.name))
            unique_molecules.append(
                (molecule_1, graph_1),
            )
        else:
            print_if_debug('NOT UNIQUE: {0}'.format(molecule_1.name))

    return [
        molecule
        for (molecule, _) in unique_molecules
    ]
