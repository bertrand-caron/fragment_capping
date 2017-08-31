from typing import Any, Optional, List, Tuple, Dict, Union, Set, FrozenSet, Sequence, TextIO
from copy import deepcopy
from operator import itemgetter
from itertools import groupby, product
from io import StringIO
from functools import reduce, lru_cache
from os.path import join
from hashlib import md5
from sys import stderr

from fragment_capping.helpers.types_helpers import Fragment, ATB_Molid, Atom, FRAGMENT_CAPPING_DIR
from fragment_capping.helpers.parameters import FULL_VALENCES, POSSIBLE_CHARGES, get_capping_options, new_atom_for_capping_strategy, Capping_Strategy, possible_bond_order_for_atom_pair, min_valence, max_valence, coordinates_n_angstroms_away_from, possible_charge_for_atom, ALL_ELEMENTS, electronegativity_spread, ELECTRONEGATIVITIES

DRAW_ALL_POSSIBLE_GRAPHS = False

MAXIMUM_PERMUTATION_NUMBER = 600000

DESC = lambda x: -x

class No_Charges_And_Bond_Orders(Exception):
    pass

def product_len(list_of_lists: List[List[Any]]) -> int:
    _product_len = reduce(lambda acc, e: acc * len(e), list_of_lists, 1)
    if _product_len > MAXIMUM_PERMUTATION_NUMBER and False:
        raise Too_Many_Permutations(_product_len)
    return _product_len

class Too_Many_Permutations(Exception):
    pass

class Not_Capped_Error(Exception):
    pass

def write_to_debug(debug: Optional[TextIO], *objects: Sequence[Any]) -> None:
     if debug is not None:
        debug.write(' '.join(map(str, objects)) + '\n')

class Molecule:
    def __init__(self, atoms: Dict[int, Atom], bonds: List[Tuple[int, int]], name: Optional[str] = None, **kwargs: Dict[str, Any]) -> None:
        self.atoms = validated_atoms_dict(atoms)
        self.bonds = set(map(frozenset, bonds))
        validate_bond_dict(self.atoms, self.bonds)
        self.name = name if name is not None else md5((str(sorted(atoms.values())) + str(sorted(bonds))).encode()).hexdigest()

        self.use_neighbour_valences = (
            True
            if all([atom.valence is not None for atom in self.atoms.values()])
            else False
        )

        self.bond_orders, self.charges = None, None
        self.previously_uncapped = set()

    def atom_desc(self, atom: Atom):
        if self.use_neighbour_valences:
            return atom.element + str(atom.valence)
        else:
            return atom.element

    def assert_use_neighbour_valences(self) -> None:
        assert self.use_neighbour_valences, 'ERROR: self.use_neighbour_valences is set to False'

    def valence(self, atom_id: int) -> int:
        self.assert_use_neighbour_valences()
        return self.atoms[atom_id].valence

    def element(self, atom_id: int) -> str:
        return self.atoms[atom_id].element

    def ids(self) -> List[int]:
        return list(self.atoms.keys())

    def __str__(self) -> str:
        try:
            netcharge = self.netcharge()
        except No_Charges_And_Bond_Orders:
            netcharge = None
        return 'Molecule(atoms={0}, bonds={1}, formula="{2}", netcharge={3})'.format(self.atoms, self.bonds, self.formula(charge=False), netcharge)

    def add_atom(self, atom: Atom, bonded_to: Optional[List[int]] = None) -> int:
        highest_id = max(self.atoms.keys())
        atom_id = highest_id + 1
        self.atoms[atom_id] = Atom(atom_id, *atom[1:])

        if bonded_to is not None:
            self.add_bonds({frozenset([bonded_atom_id, atom_id]) for bonded_atom_id in bonded_to})

        return atom_id

    def add_bond(self, bond: Tuple[int, int]) -> None:
        self.bonds.add(frozenset(bond))

    def add_bonds(self, bonds: List[Tuple[int, int]]) -> None:
        self.bonds |= set(map(frozenset, bonds))

    def capped_molecule_with(self, capping_strategies: List[Any], atoms_need_capping: Any, debug: Optional[TextIO] = None, debug_line: Optional[Any] = None, use_ILP: bool = True) -> Any:
        capped_molecule = deepcopy(self)

        for (atom, capping_strategy) in zip(atoms_need_capping, capping_strategies):
            atom_id = atom.index
            new_atoms, fragment_bonds, new_valences = capping_strategy

            last_used_id = sorted(capped_molecule.ids())[-1]
            new_ids = list(map(
                lambda id__: id__[0] + last_used_id + 1,
                enumerate(new_atoms),
            ))
            new_bonds = {
                frozenset(
                    map(
                        lambda id: atom_id if id == 0 else id + last_used_id,
                        bond,
                    ),
                )
                for bond in fragment_bonds
            }

            capped_molecule.add_bonds(new_bonds)

            assert len(new_ids) == len(new_atoms) == len(new_valences), 'Wrong dimensions: {0}, {1}, {2}'.format(
                new_ids,
                new_atoms,
                new_valences,
            )

            for (new_id, new_atom, new_valence) in zip(new_ids, new_atoms, new_valences):
                capped_molecule.atoms[new_id] = Atom(
                    index=new_id,
                    element=new_atom,
                    valence=new_valence if self.use_neighbour_valences else None,
                    capped=True,
                    coordinates=coordinates_n_angstroms_away_from(atom, 1.2),
                )
            capped_molecule.atoms[atom_id] = atom._replace(capped=True)
            self.previously_uncapped.add(atom_id)

        assert all([atom.capped for atom in capped_molecule.atoms.values()]), 'Some atoms were not capped: {0}'.format(
            [atom for atom in capped_molecule.atoms.values() if not atom.capped],
        )

        if self.use_neighbour_valences:
            capped_molecule.check_valence()

        try:
            write_to_debug(debug, debug_line)
            if use_ILP:
                capped_molecule.assign_bond_orders_and_charges_with_ILP(debug=debug)
            else:
                capped_molecule.assign_bond_orders_and_charges(debug=debug)
            return capped_molecule
        except AssertionError as e:
            write_to_debug(
                debug,
                'AssertionError for capped molecule {0}:\n{1}'.format(
                    capped_molecule,
                    str(e),
                ),
            )
            return None

    def check_valence(self) -> None:
        '''Check that the valence of each individual atom matches its current number of bonded atoms'''
        self.assert_use_neighbour_valences()

        try:
            for atom in self.atoms.values():
                atom_id = atom.index
                assert atom.valence == sum([1 for bond in self.bonds if atom_id in bond]), '{2}: expected_valence={0} != number_neighbours={1} (bonds={3})'.format(
                    atom.valence,
                    sum([1 for bond in self.bonds if atom_id in bond]),
                    atom,
                    [bond for bond in self.bonds if atom_id in bond],
                )
        except AssertionError:
            print('ERROR')
            print('Atoms were: {0}'.format(self.atoms))
            print('Bonds were: {0}'.format(self.bonds))
            raise

    def get_best_capped_molecule(self, draw_all_possible_graphs: bool = DRAW_ALL_POSSIBLE_GRAPHS, debug: Optional[TextIO] = None, use_ILP: bool = True):
        capping_options = get_capping_options(self.use_neighbour_valences)

        neighbour_counts = self.neighbours_for_atoms()

        def keep_capping_strategy_for_atom(capping_strategy: Capping_Strategy, atom: Atom):
            if self.use_neighbour_valences:
                return neighbour_counts[atom.index] + new_atom_for_capping_strategy(capping_strategy) == atom.valence
            else:
                return min_valence(atom) <= neighbour_counts[atom.index] + new_atom_for_capping_strategy(capping_strategy) <= max_valence(atom)

        def possible_capping_strategies_for_atom(atom: Atom) -> List[Capping_Strategy]:
            if debug is not None:
                write_to_debug(debug, atom)
                for capping_strategy in capping_options[self.atom_desc(atom)]:
                    write_to_debug(
                        debug, 
                        capping_strategy,
                        new_atom_for_capping_strategy(capping_strategy),
                        keep_capping_strategy_for_atom(capping_strategy, atom)
                    )
            return [
                capping_strategy
                for capping_strategy in capping_options[self.atom_desc(atom)]
                if keep_capping_strategy_for_atom(capping_strategy, atom)
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

        assert len(capping_schemes) > 0, [
            (
                atom,
                possible_capping_strategies_for_atom(atom),
            )
            for atom in atoms_need_capping
            if len(possible_capping_strategies_for_atom(atom)) == 0
        ]

        write_to_debug(
            debug,
            [
                possible_capping_strategies_for_atom(atom)
                for atom in atoms_need_capping
            ],
        )

        if len(capping_schemes) >= MAXIMUM_PERMUTATION_NUMBER:
            raise Too_Many_Permutations(len(capping_schemes))

        write_to_debug(debug, 'atoms_need_capping: {0}'.format(atoms_need_capping))
        write_to_debug(debug, 'capping_schemes: {0}'.format(capping_schemes))
        write_to_debug(debug, 'capping_options: {0}'.format([
            len(capping_options[self.atom_desc(atom)])
            for atom in atoms_need_capping
        ]))

        write_to_debug(debug, 'atoms_need_capping: {0}'.format(atoms_need_capping))
        write_to_debug(debug, 'INFO: Will try all {0} possible capped molecules'.format(len(capping_schemes)))

        possible_capped_molecules = sorted(
            filter(
                lambda mol: mol is not None,
                [
                    self.capped_molecule_with(capping_strategies, atoms_need_capping, debug=debug, debug_line='molecule {0}/{1}'.format(i, len(capping_schemes)), use_ILP=use_ILP)
                    for (i, capping_strategies) in enumerate(capping_schemes, start=1)
                ],
            ),
            key=lambda mol: (mol.net_abs_charges(), mol.n_atoms(), mol.double_bonds_fitness()),
        )

        write_to_debug(debug, 'Possible capped molecules: {0} ({1}/{2})'.format(
            [(mol.formula(charge=True), mol.net_abs_charges(), mol.double_bonds_fitness()) for mol in possible_capped_molecules],
            len(possible_capped_molecules),
            len(capping_schemes),
        ))

        if draw_all_possible_graphs:
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

        if len(possible_capped_molecules) == 0:
            raise Not_Capped_Error(self)

        best_molecule = possible_capped_molecules[0]
        return best_molecule

    def write_graph(self, unique_id: Union[str, int], **kwargs: Dict[str, Any]) -> str:
        graph_filepath = join(FRAGMENT_CAPPING_DIR, 'graphs' ,'_'.join((self.name, str(unique_id))) + '.png')
        try:
            from chem_graph_tool.pdb import graph_from_pdb
            from chem_graph_tool.moieties import draw_graph
            graph = self.graph()
            draw_graph(
                graph,
                fnme=graph_filepath,
                force_regen=True,
                **kwargs,
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
        elements = [atom.element for atom in self.atoms.values()]

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
        from chemistry_helpers.pdb import PDB_TEMPLATE, pdb_conect_line
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
            atom = self.atoms[atom_index]
            coordinates = atom.coordinates or (
                1.3 * pdb_id,
                0.1 * (-1 if pdb_id % 2 == 0 else +1),
                0.1 * (pdb_id % 5),
            )

            print(PDB_TEMPLATE.format(
                'HETATM',
                pdb_id,
                (atom.element.title() + str(atom_index))[:4],
                'R',
                '',
                pdb_id,
                *coordinates,
                '',
                '',
                atom.element.title(),
                '',
            ), file=io)

        for (atom_index, pdb_id) in sorted(pdb_ids.items(), key=itemgetter(1)):
            print(
                pdb_conect_line(
                    [pdb_id]
                    +
                    [pdb_ids[list(bond - frozenset([atom_index]))[0]] for bond in self.bonds if atom_index in bond]
                ),
                file=io,
            )

        return io.getvalue()

    def representation(self, out_format: str) -> str:
        assert out_format in ('smiles', 'inchi'), 'Wrong representation format: {0}'.format(out_format)

        from chemistry_helpers.babel import babel_output
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

        edge_types = g.new_edge_property("string")
        g.edge_properties['type'] = edge_types

        vertex_colors = g.new_vertex_property("string")
        g.vertex_properties['color'] = vertex_colors

        vertices = {}
        for (atom_index, atom) in sorted(self.atoms.items()):
            v = g.add_vertex()
            if self.charges is None:
                possible_charges = possible_charge_for_atom(atom)
                charge_str = str(possible_charges).replace(' ', '') if len(possible_charges) > 1 else ''
            else:
                charge = self.charges[atom_index]
                charge_str = (str(abs(charge)) + ('-' if charge < 0 else '+')) if charge != 0 else ''
            vertex_types[v] = '{element}{valence}{charge_str}'.format(
                element=atom.element,
                valence=atom.valence if self.use_neighbour_valences else '',
                charge_str=(' ' if charge_str else '') + charge_str,
            )
            if atom.index in self.previously_uncapped:
                vertex_colors[v] = '#90EE90'
            elif atom.capped:
                vertex_colors[v] = '#EEEEEE'
            else:
                vertex_colors[v] = '#FF91A4'
            vertices[atom_index] = v

        for bond in self.bonds:
            (i, j) = bond
            e = g.add_edge(vertices[i], vertices[j])
            if self.bond_orders is None:
                possible_bond_orders = possible_bond_order_for_atom_pair((self.atoms[i], self.atoms[j]))
                edge_text = str(possible_bond_orders)[1:-1] if len(possible_bond_orders) > 1 else ''
            else:
                edge_text = str(self.bond_orders[bond])
            edge_types[e] = edge_text

        return g

    def sorted_atoms(self) -> List[Atom]:
        return [atom for (atom_id, atom) in sorted(self.atoms.items())]

    def sorted_atom_ids(self) -> List[int]:
        return [atom_id for (atom_id, atom) in sorted(self.atoms.items())]

    def assign_bond_orders_and_charges(self, debug: Optional[TextIO] = None) -> None:
        list_of_possible_bond_orders_per_bond = [
            possible_bond_order_for_atom_pair((self.atoms[atom_id_1], self.atoms[atom_id_2]))
            for (atom_id_1, atom_id_2) in self.bonds
        ]

        number_bond_order_permutations = product_len(list_of_possible_bond_orders_per_bond)
        assert number_bond_order_permutations > 0, 'No possible bond orders found'

        possible_bond_orders_lists = product(
            *list_of_possible_bond_orders_per_bond,
        )

        list_of_possible_charges_per_atom = [
            possible_charge_for_atom(atom)
            for atom in self.sorted_atoms()
        ]

        number_charges_permutations = product_len(list_of_possible_charges_per_atom)
        assert number_charges_permutations > 0, 'No possible charges assignment found'

        possible_charges_dicts = map(
            lambda charges: dict(list(zip(self.sorted_atom_ids(), charges))),
            product(
                *list_of_possible_charges_per_atom,
            ),
        )

        len_possible_bond_orders_and_charges = number_bond_order_permutations * number_charges_permutations

        if len_possible_bond_orders_and_charges > MAXIMUM_PERMUTATION_NUMBER:
            raise Too_Many_Permutations([number_bond_order_permutations, number_charges_permutations, len_possible_bond_orders_and_charges, list_of_possible_bond_orders_per_bond, list_of_possible_charges_per_atom])

        write_to_debug(debug, 'INFO: Found {0} possible charge and bond order assignments'.format(len_possible_bond_orders_and_charges))

        possible_bond_orders_and_charges = product(possible_bond_orders_lists, possible_charges_dicts)

        acceptable_bond_orders_and_charges = sorted(
            [
                (
                    list(zip(self.bonds, bond_orders)),
                    charges,
                )
                for (i, (bond_orders, charges)) in enumerate(possible_bond_orders_and_charges)
                if self.is_valid(bond_orders, charges, debug=debug, debug_line='{0}/{1}'.format(i, len_possible_bond_orders_and_charges))
            ],
            key=lambda __charges: sum(map(abs, list(__charges[1].values()))),
        )

        if len(acceptable_bond_orders_and_charges) != 1:
            write_to_debug(debug, 'acceptable_bond_orders_and_charges: {0}'.format(acceptable_bond_orders_and_charges))

        assert len(acceptable_bond_orders_and_charges) >= 1, 'No valid bond_orders and charges found amongst {0} tried.'.format(len_possible_bond_orders_and_charges)

        best_bond_orders, self.charges = acceptable_bond_orders_and_charges[0]

        self.bond_orders = {
            bond: bond_order
            for (bond, bond_order) in best_bond_orders
        }

    def netcharge(self) -> int:
        try:
            return sum(self.charges.values())
        except:
            raise No_Charges_And_Bond_Orders('Assign charges and bond_orders first.')

    def net_abs_charges(self) -> int:
        try:
            return sum(map(abs, list(self.charges.values())))
        except:
            raise No_Charges_And_Bond_Orders('Assign charges and bond_orders first.')

    def is_valid(self, bond_orders: List[int], charges: Dict[int, int], debug: Optional[TextIO] = None, debug_line: Optional[Any] = None) -> bool:
        assert len(self.bonds) == len(bond_orders), 'Unmatched bonds and bond_orders: {0} != {1}'.format(
            len(self.bonds),
            len(bond_orders),
        )

        on_atom_id = lambda atom_id_bond_order: atom_id_bond_order[0]
        on_bond_order = lambda atom_id_bond_order1: atom_id_bond_order1[1]

        if debug is not None:
            write_to_debug(
                debug,
                'is_valid [(atom_id, element, bonded_electrons, charge, allowed_valences, is_valid) for atom_id]',
            )
            write_to_debug(
                debug,
                'is_valid',
                [
                    (
                        atom_id,
                        self.atoms[atom_id].element,
                        sum(map(on_bond_order, group)),
                        charges[atom_id],
                        FULL_VALENCES[self.atoms[atom_id].element],
                        sum(map(on_bond_order, group)) - charges[atom_id] in FULL_VALENCES[self.atoms[atom_id].element],
                    )
                    for (atom_id, group) in
                    map( # Cast groupby's iterators to lists to be able to read them several times
                        lambda T: (T[0], list(T[1])),
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
                        ),
                    )
                    if not sum(map(on_bond_order, group)) - charges[atom_id] in FULL_VALENCES[self.atoms[atom_id].element]
                ],
            )

        valid = all(
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
        )

        write_to_debug(debug, debug_line)

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
                        for (bond, bond_order) in self.bond_orders.items()
                        if bond_order == 2
                    ]
                ),
            )
        ])

        return DESC(
            sum(
                [
                    electronegativity_spread(double_bond_elements)
                    for double_bond_elements in grouped_double_bonds
                ]
            )
        )

    def neighbours_for_atoms(self) -> Dict[int, int]:
        bond_atoms = reduce(
            lambda acc, e: acc + list(e),
            self.bonds,
            [],
        )

        if len(self.atoms) == 1:
            return {
                atom_id: 0
                for atom_id in self.atoms.keys()
            }
        else:
            return {
                key: len(list(group))
                for (key, group) in groupby(
                    sorted(
                        bond_atoms,
                    )
                )
            }

    def assign_bond_orders_and_charges_with_ILP(self, debug: Optional[TextIO] = None) -> None:
        from pulp import LpProblem, LpMinimize, LpInteger, LpVariable, LpBinary

        problem = LpProblem("Lewis problem (bond order and charge assignment) for molecule {0}".format(self.name), LpMinimize)

        MIN_ABSOLUTE_CHARGE, MAX_ABSOLUTE_CHARGE = 0, 9
        MIN_BOND_ORDER, MAX_BOND_ORDER = 1, 3

        charges = {
            atom.index: LpVariable("C_{i}".format(i=atom.index), -MAX_ABSOLUTE_CHARGE, MAX_ABSOLUTE_CHARGE, LpInteger)
            for atom in self.atoms.values()
        }

        # Extra variable use to bind charges
        absolute_charges = {
            atom_id: LpVariable("Z_{i}".format(i=atom_id), MIN_ABSOLUTE_CHARGE, MAX_ABSOLUTE_CHARGE, LpInteger)
            for atom_id in charges.keys()
        }

        def can_atom_have_lone_pairs(atom: Atom) -> Atom:
            if atom.element in {'C'}:
                return False
            else:
                return True

        non_bonded_pairs = {
            atom_id: LpVariable("N_{i}".format(i=atom_id), 0, (18 / 2) if can_atom_have_lone_pairs(atom) else 0, LpInteger)
            for (atom_id, atom) in self.atoms.items()
        }

        # Maps a bond to an integer
        bond_mapping = {
            bond: i
            for (i, bond) in enumerate(self.bonds)
        }

        # Maps an integer to a bond
        bond_reverse_mapping = {v: k for (k, v) in bond_mapping.items()}

        bond_orders = {
            bond: LpVariable("B_{i}".format(i=bond_mapping[bond]), MIN_BOND_ORDER, MAX_BOND_ORDER, LpInteger)
            for bond in self.bonds
        }

        max_electronegativity_score = sum(MAX_ABSOLUTE_CHARGE * ELECTRONEGATIVITIES[self.atoms[atom_id].element] for (atom_id, atom) in charges.items())

        problem += sum(absolute_charges.values()) + (1 / max_electronegativity_score) * sum([charge * ELECTRONEGATIVITIES[self.atoms[atom_id].element] for (atom_id, charge) in charges.items()]), 'Minimise total sum of absolute charges'

        VALENCE_ELECTRONS = {
            'H': 1,
            'HE': 2,
            'LI': 1,
            'BE': 2,
            'B': 3,
            'C': 4,
            'N': 5,
            'O': 6,
            'F': 7,
            'NE': 8,
            'NA': 1,
            'MG': 2,
            'AL': 3,
            'SI': 4,
            'P': 5,
            'S': 6,
            'CL': 7,
            'AR': 8,
        }

        for atom in self.atoms.values():
            problem += charges[atom.index] == VALENCE_ELECTRONS[atom.element] - sum([bond_orders[bond] for bond in self.bonds if atom.index in bond]) - 2 * non_bonded_pairs[atom.index], '{element}_{index}'.format(element=atom.element, index=atom.index)

        # Deal with absolute values
        for atom in self.atoms.values():
            problem += charges[atom.index] <= absolute_charges[atom.index], 'Absolute charge contraint 1 {i}'.format(i=atom.index)
            problem += -charges[atom.index] <= absolute_charges[atom.index], 'Absolute charge contraint 2 {i}'.format(i=atom.index)

        problem.writeLP("lewis_{0}.lp".format(self.name))
        problem.solve()

        self.charges, self.bond_orders, self.lone_pairs = {}, {}, {}
        must_be_int = lambda x: round(x)

        for v in problem.variables():
            variable_type, atom_or_bond_index_str = v.name.split('_')
            atom_or_bond_index = int(atom_or_bond_index_str)
            if variable_type == 'C':
                self.charges[atom_or_bond_index] = must_be_int(v.varValue)
            elif variable_type == 'B':
                self.bond_orders[bond_reverse_mapping[atom_or_bond_index]] = must_be_int(v.varValue)
                pass
            elif variable_type == 'Z':
                pass
            elif variable_type == 'N':
                self.lone_pairs[atom_or_bond_index] = must_be_int(v.varValue)
            else:
                raise Exception('Unknown variable type: {0}'.format(variable_type))
        write_to_debug(debug, 'molecule_name', self.name)
        write_to_debug(debug, 'bond_orders:', self.bond_orders)
        write_to_debug(debug, 'charges', self.charges)
        write_to_debug(debug, 'lone_pairs', self.lone_pairs)

Uncapped_Molecule = Molecule

def validated_atoms_dict(atoms: Dict[int, Atom]) -> Dict[int, Atom]:
    assert {atom.element for atom in atoms.values()} <= ALL_ELEMENTS, 'Unsupported elements: {0}'.format({atom.element for atom in atoms.values()} - ALL_ELEMENTS)

    for (atom_index, atom) in atoms.items():
        assert atom_index == atom.index, 'Invalid index={0} for atom={1}'.format(atom_index, atom)

    return atoms

def validate_bond_dict(atoms: Dict[int, Atom], bonds: Set[FrozenSet[int]]) -> None:
    all_bond_indices = reduce(
        lambda acc, e: acc | e,
        bonds,
        set(),
    )

    all_atom_indices = set(atoms.keys())

    if all_atom_indices - all_bond_indices != set() and len(atoms) > 1:
        raise AssertionError('The atoms with the following indices have no bonds: {0}'.format(all_atom_indices - all_bond_indices))

    if all_bond_indices - all_atom_indices != set():
        raise AssertionError('The following atoms indices in the bonds reference non-existing atoms: {0}'.format(all_bond_indices - all_atom_indices))
