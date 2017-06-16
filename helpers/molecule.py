from typing import Any, Optional, List, Tuple, Dict, Union
from copy import deepcopy
from operator import itemgetter
from itertools import groupby, product
from io import StringIO
from functools import reduce
from os.path import join

from fragment_capping.helpers.types_helpers import Fragment, ATB_Molid, Atom
from fragment_capping.helpers.parameters import FULL_VALENCES, POSSIBLE_BOND_ORDERS, POSSIBLE_CHARGES, get_capping_options, new_atom_for_capping_strategy, POSSIBLE_BOND_ORDER_FOR_PAIR, BEST_DOUBLE_BONDS

from dihedral_fragments.dihedral_fragment import element_valence_for_atom, on_asc_number_electron_then_asc_valence, NO_VALENCE

DRAW_ALL_GRAPHS = True

DEBUG = False

DEBUG_MOLECULE = lambda molecule: (molecule.formula() == 'C6H12')

MAXIMUM_PERMUTATION_NUMBER = 10000

class Too_Many_Permutations(Exception):
    pass

class Molecule:
    def __init__(self, atoms: Any, bonds: Any, name: Optional[str] = None) -> None:
        self.atoms = atoms
        self.bonds = bonds
        self.name = name

        self.use_neighbour_valences = (
            True
            if all([atom.valence is not NO_VALENCE for atom in list(self.atoms.values())])
            else False
        )

    def atom_desc(self, atom: Atom):
        if self.use_neighbour_valences:
            return atom.element + str(atom.valence)
        else:
            return atom.element

    def assert_use_neighbour_valences(self) -> None:
        assert self.use_neighbour_valences, 'ERROR: self.use_neighbour_valences is set to False'

    def valence(self, atom_id: int):
        self.assert_use_neighbour_valences()
        return self.atoms[atom_id].valence

    def element(self, atom_id: int):
        return self.atoms[atom_id].element

    def ids(self):
        return list(self.atoms.keys())

    def __str__(self) -> str:
        return 'Molecule; atoms={0}; bonds={1}; formula={2}'.format(self.atoms, self.bonds, self.formula(charge=False))

    def add_atom(self, atom: Atom) -> Any:
        highest_id = max(self.atoms.keys())
        atom_id = highest_id + 1
        self.atoms[atom_id] = Atom(atom_id, *atom[1:])
        return atom_id

    def capped_molecule_with(self, capping_strategies: List[Any], atoms_need_capping: Any) -> Any:
        capped_molecule = deepcopy(self)

        for (atom, capping_strategy) in zip(atoms_need_capping, capping_strategies):
            atom_id = atom.index
            new_atoms, fragment_bonds, new_valences = capping_strategy

            last_used_id = sorted(capped_molecule.ids())[-1]
            new_ids = list(map(
                lambda id__: id__[0] + last_used_id + 1,
                enumerate(new_atoms),
            ))
            new_bonds = [
                tuple(
                    list(map(
                        lambda id: atom_id if id == 0 else id + last_used_id,
                        bond,
                    )),
                )
                for bond in fragment_bonds
            ]

            capped_molecule.bonds += new_bonds

            assert len(new_ids) == len(new_atoms) == len(new_valences), 'Wrong dimensions: {0}, {1}, {2}'.format(
                new_ids,
                new_atoms,
                new_valences,
            )

            for (new_id, new_atom, new_valence) in zip(new_ids, new_atoms, new_valences):
                capped_molecule.atoms[new_id] = Atom(
                    new_id,
                    new_atom,
                    new_valence if self.use_neighbour_valences else NO_VALENCE,
                    True,
                )
            capped_molecule.atoms[atom_id] = Atom(*atom[:-1], True)

        assert all([atom.capped for atom in list(capped_molecule.atoms.values())]), 'Some atoms were not capped: {0}'.format(
            [atom for atom in list(capped_molecule.atoms.values()) if not atom.capped],
        )

        if self.use_neighbour_valences:
            capped_molecule.check_valence()

        try:
            capped_molecule.assign_bond_orders_and_charges()
            return capped_molecule
        except AssertionError as e:
            if DEBUG and DEBUG_MOLECULE(capped_molecule):
                print('AssertionError for capped molecule {0}:\n{1}'.format(
                    capped_molecule,
                    str(e),
                ))
            return None

    def check_valence(self) -> None:
        self.assert_use_neighbour_valences()

        try:
            for atom in list(self.atoms.values()):
                atom_id = atom.index
                assert atom.valence == sum([1 for bond in self.bonds if atom_id in bond]), 'Atom {2}: {0} != {1} (bonds={3})'.format(
                    atom.valence,
                    sum([1 for bond in self.bonds if atom_id in bond]),
                    atom,
                    [bond for bond in self.bonds if atom_id in bond],
                )
        except:
            print('ERROR')
            print('Atoms are: {0}'.format(self.atoms))
            print('Bonds are: {0}'.format(self.bonds))
            raise

    def get_best_capped_molecule(self, debug: bool = False):
        capping_options = get_capping_options(self.use_neighbour_valences)

        neighbour_counts = self.neighbours_for_atoms()

        def possible_capping_strategies_for_atom(atom: Atom):
            if DEBUG:
                for capping_strategy in capping_options[self.atom_desc(atom)]:
                    print(
                        atom,
                        capping_strategy,
                        new_atom_for_capping_strategy(capping_strategy),
                        FULL_VALENCES[self.atom_desc(atom)] - min(POSSIBLE_CHARGES[atom.element]) - neighbour_counts[atom.index],
                    )
            return [
                capping_strategy
                for capping_strategy in capping_options[self.atom_desc(atom)]
                if new_atom_for_capping_strategy(capping_strategy) <= max(FULL_VALENCES[self.atom_desc(atom)]) - min(POSSIBLE_CHARGES[atom.element]) - neighbour_counts[atom.index]
            ]

        atoms_need_capping = [atom for atom in self.sorted_atoms() if not atom.capped]
        capping_schemes = list(
            product(
                *[
                    possible_capping_strategies_for_atom(atom)
                    for atom in atoms_need_capping
                ]
            ),
        )

        if len(capping_schemes) >= MAXIMUM_PERMUTATION_NUMBER:
            raise Too_Many_Permutations(len(capping_schemes))

        if debug:
            print('atoms_need_capping: {0}'.format(atoms_need_capping))
            print('capping_schemes: {0}'.format(capping_schemes))
            print('capping_options: {0}'.format([
                len(capping_options[self.atom_desc(atom)])
                for atom in atoms_need_capping
            ]))

        possible_capped_molecules = sorted(
            filter(
                lambda mol: mol is not None,
                [
                    self.capped_molecule_with(capping_strategies, atoms_need_capping)
                    for capping_strategies in capping_schemes
                ],
            ),
            key=lambda mol: (mol.net_abs_charges(), mol.n_atoms(), mol.double_bonds_fitness()),
        )

        print('Possible capped molecules: {0} ({1}/{2})'.format(
            [(mol.formula(charge=True), mol.net_abs_charges(), mol.double_bonds_fitness()) for mol in possible_capped_molecules],
            len(possible_capped_molecules),
            len(capping_schemes),
        ))

        if DRAW_ALL_GRAPHS:
            try:
                for (i, molecule) in enumerate(possible_capped_molecules):
                    molecule.write_graph(i)
            except Exception as e:
                print(
                    'ERROR: Could not plot graphs (error was: {0})'.format(
                        str(e),
                    ),
                )
                raise

        best_molecule = possible_capped_molecules[0]
        return best_molecule

    def write_graph(self, unique_id: Union[str, int]) -> str:
        graph_filepath = join('graphs' ,'_'.join((self.name, str(unique_id))) + '.png')
        try:
            from py_graphs.pdb import graph_from_pdb
            from py_graphs.moieties import draw_graph
            graph = self.graph()
            draw_graph(
                graph,
                fnme=graph_filepath,
                force_regen=True,
            )
        except Exception as e:
            print(
                'ERROR: Could not plot graphs (error was: {0})'.format(
                    str(e),
                ),
            )
            raise
        return graph_filepath

    def formula(self, charge: bool = False) -> str:
        elements =  [atom.element for atom in list(self.atoms.values())]

        return ''.join(
            list(map(
                lambda element_number: element_number[0] + (str(element_number[1]) if element_number[1] > 1 else ''),
                [
                    (key, len(list(group)))
                    for (key, group) in
                    groupby(
                        sorted(
                            elements,
                            key=lambda element: (element != 'C', element),
                        ),
                    )
                ],
            ))
            +
            (
                    [
                        (' ' + str(abs(self.netcharge()))) if self.netcharge() != 0 else '',
                        '+' if self.netcharge() > 0 else ('-' if self.netcharge() < 0 else ''),
                    ]
                    if charge == True
                    else []
            ),
        )

    def n_atoms(self) -> int:
        return len(self.atoms)

    def dummy_pdb(self) -> str:
        from atb_helpers.pdb import PDB_TEMPLATE, pdb_conect_line
        io = StringIO()

        ordered_atoms = sorted(
            self.atoms.values(),
            key=lambda atom: atom.index,
        )

        pdb_ids = dict(
            zip(
                [atom.index for atom in ordered_atoms],
                range(1, len(ordered_atoms) + 1),
            ),
        )

        for (atom_index, pdb_id) in sorted(pdb_ids.items(), key=itemgetter(1)):
            print(PDB_TEMPLATE.format(
                'HETATM',
                pdb_id,
                'D',
                'R',
                '',
                pdb_id,
                1.5 * pdb_id,
                0.1 * (-1 if pdb_id % 2 == 0 else +1),
                0.,
                '',
                '',
                self.atoms[atom_index].element.title(),
                '',
            ), file=io)

        for (atom_index, pdb_id) in sorted(pdb_ids.items(), key=itemgetter(1)):
            print(
                pdb_conect_line(
                    [pdb_id]
                    +
                    [pdb_ids[bond[0] if bond[1] == atom_index else bond[1]] for bond in self.bonds if atom_index in bond]
                ),
                file=io,
            )

        return io.getvalue()

    def representation(self, out_format: str) -> str:
        assert out_format in ('smiles', 'inchi'), 'Wrong representation format: {0}'.format(out_format)

        from atb_helpers.babel import babel_output
        return babel_output(
            self.dummy_pdb(),
            in_format='pdb',
            out_format=out_format,
            dont_add_H=True,
        )

    def smiles(self) -> str:
        return self.representation('smiles')

    def inchi(self) -> str:
        return self.representation('inchi')

    def graph(self) -> Any:
        try:
            from graph_tool.all import Graph
        except:
            return None

        g = Graph(directed=False)

        vertex_types = g.new_vertex_property("string")
        g.vertex_properties['type'] = vertex_types

        vertices = {}
        for atom_index in sorted(self.atoms.keys()):
            v = g.add_vertex()
            vertex_types[v] = '{element}{valence}'.format(
                element=self.atoms[atom_index].element,
                valence=self.atoms[atom_index].valence if self.use_neighbour_valences else '',
            )
            vertices[atom_index] = v

        for (i, j) in self.bonds:
            g.add_edge(vertices[i], vertices[j])

        return g

    def sorted_atoms(self) -> List[Atom]:
        return [atom for (atom_id, atom) in sorted(self.atoms.items())]

    def sorted_atom_ids(self) -> List[int]:
        return [atom_id for (atom_id, atom) in sorted(self.atoms.items())]

    def assign_bond_orders_and_charges(self) -> None:
        possible_bond_orders_lists = list(
            product(
                *[
                    POSSIBLE_BOND_ORDER_FOR_PAIR[element_1, element_2]
                    for (element_1, element_2) in
                    list(
                        map(
                            lambda bond: (self.atoms[bond[0]].element, self.atoms[bond[1]].element),
                            self.bonds,
                        )
                    )
                ]
            ),
        )

        assert len(possible_bond_orders_lists) >= 1, 'No possible bond orders found'

        possible_charges_dicts = list(map(
            lambda charges: dict(list(zip(self.sorted_atom_ids(), charges))),
            product(
                *[
                    POSSIBLE_CHARGES[atom.element]
                    for atom in self.sorted_atoms()
                ]
            ),
        ))

        assert len(possible_charges_dicts) >= 1, 'No possible charges assignment found'

        possible_bond_orders_and_charges = list(product(possible_bond_orders_lists, possible_charges_dicts))

        print('INFO: Found {0} posible charge and bond order assignments'.format(len(possible_bond_orders_and_charges)))

        acceptable_bond_orders_and_charges = sorted(
            [
                (
                    list(zip(self.bonds, bond_orders)),
                    charges,
                )
                for (bond_orders, charges) in possible_bond_orders_and_charges
                if self.is_valid(bond_orders, charges)
            ],
            key=lambda __charges: sum(map(abs, list(__charges[1].values()))),
        )

        if DEBUG and DEBUG_MOLECULE(self):
            if len(acceptable_bond_orders_and_charges) != 1:
                print('acceptable_bond_orders_and_charges: {0}'.format(acceptable_bond_orders_and_charges))

        assert len(acceptable_bond_orders_and_charges) >= 1, 'No valid bond_orders and charges found amongst {0} tried.'.format(len(possible_bond_orders_and_charges))

        self.bond_orders, self.charges = acceptable_bond_orders_and_charges[0]

    def netcharge(self) -> int:
        try:
            return sum(self.charges.values())
        except:
            raise Exception('Assign charges and bond_orders first.')

    def net_abs_charges(self) -> int:
        try:
            return sum(map(abs, list(self.charges.values())))
        except:
            raise Exception('Assign charges and bond_orders first.')

    def is_valid(self, bond_orders: List[int], charges: Dict[int, int]) -> bool:
        assert len(self.bonds) == len(bond_orders), 'Unmatched bonds and bond_orders: {0} != {1}'.format(
            len(self.bonds),
            len(bond_orders),
        )

        on_atom_id = lambda atom_id_bond_order: atom_id_bond_order[0]
        on_bond_order = lambda atom_id_bond_order1: atom_id_bond_order1[1]

        if DEBUG and DEBUG_MOLECULE(self):
            print(
                'is_valid',
                [
                    (atom_id, sum(map(on_bond_order, group)), FULL_VALENCES[self.atoms[atom_id].element] + charges[atom_id])
                    for (atom_id, group) in
                    groupby(
                        sorted(
                            reduce(
                                lambda acc, e: acc + e,
                                [
                                    ((atom_id_1, bond_order), (atom_id_2, bond_order))
                                    for ((atom_id_1, atom_id_2), bond_order) in
                                    zip(self.bonds, bond_orders)
                                ],
                                (),
                            ),
                            key=on_atom_id,
                        ),
                        key=on_atom_id,
                    )
                ],
            )

        valid = all(
            [
                sum(map(on_bond_order, group)) - charges[atom_id] in FULL_VALENCES[self.atoms[atom_id].element]
                for (atom_id, group) in
                groupby(
                    sorted(
                        reduce(
                            lambda acc, e: acc + e,
                            [
                                ((atom_id_1, bond_order), (atom_id_2, bond_order))
                                for ((atom_id_1, atom_id_2), bond_order) in
                                zip(self.bonds, bond_orders)
                            ],
                            (),
                        ),
                        key=on_atom_id,
                    ),
                    key=on_atom_id,
                )
            ],
        )

        return valid

    def double_bonds_fitness(self) -> Tuple[int, int, int]:
        '''Sorted ASC (low fitness is better)'''

        grouped_double_bonds = dict([
            (key, len(list(group)))
            for (key, group) in
            groupby(
                sorted(
                    [
                        frozenset([self.atoms[atom_id].element for atom_id in bond])
                        for (bond, bond_order) in self.bond_orders
                        if bond_order == 2
                    ]
                ),
            )
        ])

        assert set(grouped_double_bonds.keys()) <= set(BEST_DOUBLE_BONDS), 'Unexpected double bonds: {0}'.format(set(grouped_double_bonds.keys()) - set(BEST_DOUBLE_BONDS))

        return tuple(
            [
                (- grouped_double_bonds[double_bond_type] if double_bond_type in grouped_double_bonds else 0)
                for double_bond_type in BEST_DOUBLE_BONDS
            ]
        )

    def neighbours_for_atoms(self) -> Dict[int, int]:
        bond_atoms = reduce(
            lambda acc, e: acc + e,
            self.bonds,
            (),
        )

        return {
            key: len(list(group))
            for (key, group) in groupby(
                sorted(
                    bond_atoms,
                )
            )
        }

Uncapped_Molecule = Molecule
