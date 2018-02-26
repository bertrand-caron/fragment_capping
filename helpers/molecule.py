from typing import Any, Optional, List, Tuple, Dict, Union, Set, Sequence, TextIO, FrozenSet, Callable
from copy import deepcopy
from operator import itemgetter
from itertools import groupby, product, count
from io import StringIO
from functools import reduce, lru_cache
from os.path import join
from hashlib import md5
from sys import stderr
from math import sqrt
from warnings import warn

from fragment_capping.helpers.types_helpers import Fragment, ATB_Molid, Atom, FRAGMENT_CAPPING_DIR, Bond, ATOM_INDEX, MIN, MAX
from fragment_capping.helpers.parameters import FULL_VALENCES, POSSIBLE_CHARGES,  Capping_Strategy, possible_bond_order_for_atom_pair, coordinates_n_angstroms_away_from, possible_charge_for_atom, ALL_ELEMENTS, electronegativity_spread, ELECTRONEGATIVITIES, VALENCE_ELECTRONS, MIN_ABSOLUTE_CHARGE, MAX_ABSOLUTE_CHARGE, MIN_BOND_ORDER, MAX_BOND_ORDER, MUST_BE_INT, MAX_NONBONDED_ELECTRONS, NO_CAP, ELECTRONS_PER_BOND, ALL_CAPPING_OPTIONS
from fragment_capping.helpers.babel import energy_minimised_pdb
from fragment_capping.helpers.rings import bonds_for_ring, atom_indices_in_phenyl_rings_for
from fragment_capping.helpers.graphs import unique_molecules
from fragment_capping.helpers.tautomers import get_all_tautomers, get_all_tautomers_naive
from fragment_capping.helpers.capping import get_best_capped_molecule_with_ILP, get_best_capped_molecule
from fragment_capping.helpers.misc import write_to_debug
from fragment_capping.helpers.exceptions import No_Charges_And_Bond_Orders, Too_Many_Permutations, Not_Capped_Error

PropertyMap = Any

class Molecule:
    def __init__(
        self,
        atoms: Union[Dict[ATOM_INDEX, Atom], Sequence[Atom]],
        bonds: List[Tuple[int, int]],
        name: Optional[str] = None,
        net_charge: Optional[int] = None,
        formal_charges: Optional[Dict[ATOM_INDEX, int]] = None,
        bond_orders: Optional[Dict[Bond, int]] = None,
        **kwargs: Dict[str, Any]
    ) -> None:
        if type(atoms) in (list, set, frozenset):
            assert all(isinstance(atom, Atom) for atom in atoms), 'Invalid atom types: {0}'.format([type(atom) for atom in atoms])
            atoms_dict = {atom.index: atom for atom in atoms}
        elif isinstance(atoms, dict):
            assert all(isinstance(atom, Atom) for atom in atoms.values()), 'Invalid atom types: {0}'.format([type(atom) for atom in atoms])
            atoms_dict = atoms
        else:
            raise Exception('Invalide type(atoms): {0}'.format(type(atoms)))

        self.atoms = validated_atoms_dict(atoms_dict)
        self.bonds = set(map(frozenset, bonds))
        validate_bond_dict(self.atoms, self.bonds)
        self.name = name if name is not None else md5((str(sorted(self.atoms.values())) + str(sorted(bonds))).encode()).hexdigest()
        self.net_charge = net_charge

        if formal_charges is not None:
            assert all(isinstance(formal_charge, int) for (_, formal_charge) in formal_charges.items()), 'Non-integer formal charges: {0}'.format(
                {atom: formal_charge for (atom, formal_charge) in formal_charges.items() if not isinstance(formal_charge, int)}
            )

        if bond_orders is not None:
            assert all(isinstance(bond_order, int) for (_, bond_order) in bond_orders.items()), 'Non-integer bond orders: {0}'.format(
                {bond: bond_order for (bond, bond_order) in bond_orders.items() if not isinstance(bond_order, int)}
            )

        self.bond_orders, self.formal_charges = bond_orders, formal_charges
        self.previously_uncapped = set()
        self.aromatic_bonds, self.aromatic_rings = (None, None)
        self.non_bonded_electrons = None
        self.phenyl_atoms = set()

    def use_neighbour_valences(self) -> bool:
        return all([atom.valence is not None for atom in self.atoms.values()])

    def atom_desc(self, atom: Atom):
        if self.use_neighbour_valences():
            return atom.element + str(atom.valence)
        else:
            return atom.element

    def assert_use_neighbour_valences(self) -> None:
        assert self.use_neighbour_valences(), 'ERROR: self.use_neighbour_valences() returned False'

    def valence(self, atom_id: ATOM_INDEX) -> int:
        self.assert_use_neighbour_valences()
        return self.atoms[atom_id].valence

    def element(self, atom_id: ATOM_INDEX) -> str:
        return self.atoms[atom_id].element

    def ids(self) -> List[ATOM_INDEX]:
        return list(self.atoms.keys())

    def assert_no_floating_atoms(self) -> None:
        if len(self.atoms) > 1:
            assert all([len([bond for bond in self.bonds if atom_index in bond]) > 0 for atom_index in self.atoms.keys()])
        return None

    def __str__(self) -> str:
        try:
            netcharge = self.netcharge()
        except No_Charges_And_Bond_Orders:
            netcharge = None

        return 'Molecule(atoms={0}, bonds={1}, formula="{2}", netcharge={3}, charges={4}, bond_orders={5}, non_bonded_electrons={6}, name={7})'.format(
            self.atoms,
            self.bonds,
            self.formula(charge=False),
            netcharge,
            self.formal_charges,
            self.bond_orders,
            self.non_bonded_electrons,
            self.name,
        )

    def __lt__(self, other) -> bool:
        return self.name < other.name

    def add_atom(self, atom: Atom, bonded_to: Optional[List[int]] = None) -> int:
        highest_id = max(self.atoms.keys())
        atom_id = highest_id + 1
        self.atoms[atom_id] = Atom(atom_id, *atom[1:])

        if bonded_to is not None:
            self.add_bonds({frozenset([bonded_atom_id, atom_id]) for bonded_atom_id in bonded_to})

        return atom_id

    def remove_atom_with_index(self, atom_index: int) -> None:
        del self.atoms[atom_index]
        old_bonds = deepcopy(self.bonds)
        self.bonds = {bond for bond in self.bonds if atom_index not in bond}
        if self.formal_charges:
            del self.formal_charges[atom_index]
        if self.non_bonded_electrons:
            del self.non_bonded_electrons[atom_index]
        for bond in old_bonds - self.bonds:
            del self.bond_orders[bond]

    def remove_atom(self, atom: Atom) -> None:
        return remove_atom_with_index(atom.index)

    def remove_atoms(self, atoms: List[Atom]) -> None:
        [self.remove_atom_with_index(atom.index) for atom in atoms]

    def add_bond(self, bond: Tuple[int, int]) -> None:
        self.bonds.add(frozenset(bond))

    def add_bonds(self, bonds: List[Tuple[int, int]]) -> None:
        self.bonds |= set(map(frozenset, bonds))

    def extend_molecule_with(self, atom: Atom, capping_strategy: Capping_Strategy) -> Tuple[List[Atom], List[Bond]]:
        atom_id = atom.index
        new_atoms, fragment_bonds, new_valences = capping_strategy

        last_used_id = sorted(self.ids())[-1]
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

        self.add_bonds(new_bonds)

        assert len(new_ids) == len(new_atoms) == len(new_valences), 'Wrong dimensions: {0}, {1}, {2}'.format(
            new_ids,
            new_atoms,
            new_valences,
        )

        for (new_id, new_atom, new_valence) in zip(new_ids, new_atoms, new_valences):
            self.atoms[new_id] = Atom(
                index=new_id,
                element=new_atom,
                valence=new_valence if self.use_neighbour_valences() else None,
                capped=True,
                coordinates=coordinates_n_angstroms_away_from(atom, 1.2),
            )
        self.atoms[atom_id] = atom._replace(capped=True)
        self.previously_uncapped.add(atom_id)

        return ([atom for (atom_id, atom) in self.atoms.items() if atom_id in new_ids], new_bonds)

    def capped_molecule_with(self, capping_strategies: List[Any], atoms_need_capping: Any, debug: Optional[TextIO] = None, debug_line: Optional[Any] = None, use_ILP: bool = True) -> Any:
        capped_molecule = deepcopy(self)

        for (atom, capping_strategy) in zip(atoms_need_capping, capping_strategies):
            capped_molecule.extend_molecule_with(atom, capping_strategy)

        assert all([atom.capped for atom in capped_molecule.atoms.values()]), 'Some atoms were not capped: {0}'.format(
            [atom for atom in capped_molecule.atoms.values() if not atom.capped],
        )

        if capped_molecule.use_neighbour_valences():
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

    def write_graph(
        self,
        unique_id: Union[str, int],
        graph_kwargs: Dict[str, Any] = {},
        sfdp_layout_kwargs: Dict[str, Any] = {},
        output_size: Optional[Tuple[int, int]] = None,
        g: Optional[Any] = None,
        pos: Optional[PropertyMap] = None,
        **kwargs: Dict[str, Any]
    ) -> Optional[Tuple[str, PropertyMap, Any]]:
        if output_size is None:
            output_size = [int(sqrt(len(self.atoms) / 3) * 300)] * 2

        graph_filepath = join(FRAGMENT_CAPPING_DIR, 'graphs' ,'_'.join((self.name, str(unique_id))) + '.pdf')
        try:
            SOFT_IMPORT_FAIL = True
            try:
                from chem_graph_tool.moieties import draw_graph, sfdp_layout
            except ImportError:
                if SOFT_IMPORT_FAIL:
                    stderr.write('\n' + "ERROR: Please download chem_graph_tool (https://github.com/bertrand-caron/chem_graph_tool) and make sure it can be found in your PYTHONPATH. You'll also need graph-tool (https://graph-tool.skewed.de)")
                    return None
                else:
                    raise ImportError("Please download chem_graph_tool (https://github.com/bertrand-caron/chem_graph_tool) and make sure it can be found in your PYTHONPATH. You'll also need graph-tool (https://graph-tool.skewed.de)")
            graph = self.graph(g=g, **graph_kwargs)

            if pos is None:
                pos = sfdp_layout(
                    graph,
                    **{
                        **sfdp_layout_kwargs,
                        #**dict(p=3, epsilon=0.001, weighted_coarse=True),
                    }
                )
            assert pos is not None, pos

            draw_graph(
                graph,
                pos=pos,
                fnme=graph_filepath,
                force_regen=True,
                **{**{'output_size': output_size}, **kwargs},
            )
        except Exception as e:
            print(
                'ERROR: Could not plot graphs (error was: {0})'.format(
                    str(e),
                ),
            )
            raise
        return (graph_filepath, graph, pos)

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
                            key=lambda element: (element != 'C', element != 'H', element), # Hill ordering
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

    def n_electrons(self) -> int:
        return 2 * sum(self.bond_orders.values()) + sum(self.non_bonded_electrons.values())

    def energy_minimised_pdb(self, **kwargs: Dict[str, Any]) -> str:
        return energy_minimised_pdb(
            pdb_str=self.dummy_pdb(),
            **kwargs,
        )

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

    def mol2(self) -> str:
        from chemistry_helpers.pdb import PDB_TEMPLATE, pdb_conect_line
        io = StringIO()

        def print_to_io(*args):
            print(*args, file=io)

        def sibyl_atom_type(atom: Atom) -> str:
            if atom.element in {'C', 'N', 'O', 'S', 'P'}:
                if atom.element == 'C':
                    valence = atom.valence - 1
                elif atom.element == 'O':
                    valence = atom.valence + 1
                elif atom.element == 'S':
                    if atom.valence == 3:
                        valence = 'o'
                    elif atom.valence == 4:
                        valence = 'o2'
                    elif atom.valence == 2:
                        valence = 3
                    else:
                        valence = 2
                else:
                    valence = atom.valence
                return '{element}.{valence}'.format(
                    element=atom.element,
                    valence=valence,
                )
            else:
                return atom.element

        print_to_io('@<TRIPOS>MOLECULE')
        print_to_io(self.name)
        print_to_io(
            '{num_atoms} {num_bonds}'.format(
                num_atoms=len(self.atoms),
                num_bonds=len(self.bonds)
            ),
        )
        print_to_io('SMALL')
        print_to_io('USER_CHARGES')

        print_to_io('@<TRIPOS>ATOM')
        for atom in self.sorted_atoms():
            print_to_io(
                '{index} {name} {coordinates} {sibyl_atom_type} {subst_id} {subst_name} {charge}'.format(
                    index=atom.index,
                    name='A' + str(atom.index),
                    coordinates=' '.join(map(str, atom.coordinates)),
                    sibyl_atom_type=sibyl_atom_type(atom),
                    subst_id=1,
                    subst_name='<1>',
                    charge=float(self.formal_charges[atom.index]),
                ),
            )

        assert self.aromatic_bonds is not None, 'Aromatic bonds were not assigned. Use Molecule.assign_aromatic_bonds()'
        print_to_io('@<TRIPOS>BOND')
        for (bond_id, bond) in enumerate(self.bonds):
            print_to_io(
                '{0} {1} {2} {3}'.format(
                    bond_id,
                    *list(bond),
                    self.bond_orders[bond] if bond not in self.aromatic_bonds else 'ar',
                ),
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

    def graph(self, g: Optional[Any] = None, include_atom_index: bool = True, include_atom_valence: bool = False) -> Any:
        try:
            from graph_tool.all import Graph
        except:
            return None

        if g is None:
            g = Graph(directed=False)

            vertex_types = g.new_vertex_property("string")
            g.vertex_properties['type'] = vertex_types

            edge_types = g.new_edge_property("string")
            g.edge_properties['type'] = edge_types

            vertex_colors = g.new_vertex_property("string")
            g.vertex_properties['color'] = vertex_colors

            _vertices = {}
            for (atom_index, atom) in sorted(self.atoms.items()):
                v = g.add_vertex()
                _vertices[atom_index] = v

            _edges = {}
            for bond in self.bonds:
                (i, j) = bond
                e = g.add_edge(_vertices[i], _vertices[j])
                _edges[bond] = e
        else:
            (vertex_types, vertex_colors) = map(lambda type_str: g.vertex_properties[type_str], ['type', 'color'])
            (edge_types,) = map(lambda type_str: g.edge_properties[type_str], ['type'])

        def non_bonded_str(atom_index: int) -> str:
            if self.non_bonded_electrons is not None:
                return '*' * (self.non_bonded_electrons[atom_index] // 2) + '.' * (self.non_bonded_electrons[atom_index] % 2)
            else:
                return ''

        for (atom_index, atom) in sorted(self.atoms.items()):
            v = _vertices[atom_index]
            if self.formal_charges is None:
                possible_charges = possible_charge_for_atom(atom)
                charge_str = '?'
            else:
                charge = self.formal_charges[atom_index]
                charge_str = ((str(abs(charge)) if abs(charge) != 1 else '') + ('-' if charge < 0 else '+')) if charge != 0 else ''
            vertex_types[v] = '{element}{valence}{charge_str}{non_bonded_str}{index}'.format(
                element=atom.element,
                valence=atom.valence if self.use_neighbour_valences() and include_atom_valence else '',
                charge_str=('' if charge_str else '') + charge_str,
                non_bonded_str=non_bonded_str(atom_index),
                index=' ({0})'.format(atom.index) if include_atom_index else '',
            )
            if atom.index in self.previously_uncapped:
                vertex_colors[v] = '#90EE90' # Green
            elif atom.capped:
                vertex_colors[v] = '#EEEEEE' # Grey
            else:
                vertex_colors[v] = '#FF91A4' # Red

        for bond in self.bonds:
            e = _edges[bond]
            if self.bond_orders is None:
                possible_bond_orders = possible_bond_order_for_atom_pair((self.atoms[i], self.atoms[j]))
                edge_text = '?'
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

        best_bond_orders, self.formal_charges = acceptable_bond_orders_and_charges[0]

        self.bond_orders = {
            bond: bond_order
            for (bond, bond_order) in best_bond_orders
        }
        self.assign_aromatic_bonds()

    def assert_molecule_coherence(self) -> None:
        bonded_atoms = reduce(lambda acc, e: acc | e, self.bonds, frozenset())
        try:
            if len(self.atoms) > 1:
                assert set(bonded_atoms) == set(self.atoms.keys()), set(bonded_atoms) ^ set(self.atoms.keys())
            assert set(self.bonds) == set(self.bond_orders.keys()), set(self.bonds) ^ set(self.bond_orders.keys())
            assert set(self.atoms.keys()) == set(self.formal_charges.keys()), set(self.atoms.keys()) ^ set(self.formal_charges.keys())
        except AssertionError:
            stderr.write('\n' + str(self))
            raise

    def netcharge(self) -> int:
        try:
            return sum(self.formal_charges.values())
        except:
            raise No_Charges_And_Bond_Orders('Assign charges and bond_orders first.')

    def net_abs_charges(self) -> int:
        try:
            return sum(map(abs, list(self.formal_charges.values())))
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

    def assign_bond_orders_and_charges_with_ILP(
        self,
        enforce_octet_rule: bool = True,
        allow_radicals: bool = False,
        debug: Optional[TextIO] = None,
        bond_order_constraints: List[Tuple[Bond, int]] = [],
    ) -> None:
        from pulp import LpProblem, LpMinimize, LpInteger, LpVariable, LpBinary, LpStatus
        from pulp.solvers import PulpSolverError

        assert not (allow_radicals and enforce_octet_rule), "Can't simultaneously allow_radicals and enforce octet rule."

        problem = LpProblem("Lewis problem (bond order and charge assignment) for molecule {0}".format(self.name), LpMinimize)

        ELECTRON_MULTIPLIER = (2 if not allow_radicals else 1)

        charges = {
            atom.index: LpVariable("C_{i}".format(i=atom.index), -MAX_ABSOLUTE_CHARGE, MAX_ABSOLUTE_CHARGE, LpInteger)
            for atom in self.atoms.values()
        }

        # Extra variable use to bind charges
        absolute_charges = {
            atom_id: LpVariable("Z_{i}".format(i=atom_id), MIN_ABSOLUTE_CHARGE, MAX_ABSOLUTE_CHARGE, LpInteger)
            for atom_id in charges.keys()
        }

        non_bonded_electrons = {
            atom_id: LpVariable("N_{i}".format(i=atom_id), 0, MAX_NONBONDED_ELECTRONS // ELECTRON_MULTIPLIER, LpInteger)
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

        OBJECTIVES = [
            MIN(sum(absolute_charges.values())),
            #FIXME: sum(charges.values()) as close to zero as possible (cf example_wang_8 and example_wang_9)
            MIN(sum([charge * ELECTRONEGATIVITIES[self.atoms[atom_id].element] for (atom_id, charge) in charges.items()])),
        ] + (
            [MIN(sum([bond_order * ELECTRONEGATIVITIES[self.atoms[atom_id].element] for (bond, bond_order) in bond_orders.items() for atom_id in bond]))]
            if len(bond_orders) > 0
            else []
        )

        if self.net_charge is not None:
            problem += sum(charges.values()) == self.net_charge, 'Total net charge'

        for atom in self.atoms.values():
            problem += charges[atom.index] == VALENCE_ELECTRONS[atom.element] - sum([bond_orders[bond] for bond in self.bonds if atom.index in bond]) - ELECTRON_MULTIPLIER * non_bonded_electrons[atom.index], '{element}_{index}'.format(element=atom.element, index=atom.index)

        # Deal with absolute values
        for atom in self.atoms.values():
            problem += charges[atom.index] <= absolute_charges[atom.index], 'Absolute charge contraint 1 {i}'.format(i=atom.index)
            problem += -charges[atom.index] <= absolute_charges[atom.index], 'Absolute charge contraint 2 {i}'.format(i=atom.index)

            if enforce_octet_rule:
                if atom.element not in {'B', 'BE', 'P', 'S'}:
                    problem += (
                        ELECTRONS_PER_BOND * sum([bond_orders[bond] for bond in self.bonds if atom.index in bond]) + ELECTRON_MULTIPLIER * non_bonded_electrons[atom.index] == (2 if atom.element in {'H', 'HE'} else 8),
                        'Octet for atom {element}_{index}'.format(element=atom.element, index=atom.index),
                    )

        for (bond, bond_order) in bond_order_constraints:
            problem += bond_orders[frozenset(bond)] == bond_order, 'Constraint bond {0} == {1}'.format(bond, bond_order)

        try:
            problem.sequentialSolve(OBJECTIVES)
            assert problem.status == 1, (self.name, LpStatus[problem.status])
        except (AssertionError, PulpSolverError) as e:
            args_id = ','.join(map(str, [enforce_octet_rule, allow_radicals, bond_order_constraints]))
            debug_file = '{0}_{1}_debug.lp'.format(self.name, args_id)
            problem.writeLP(debug_file)
            self.write_graph('DEBUG', output_size=(1000, 1000))
            stderr.write('\n' + 'Failed LP written to "{0}"'.format(debug_file))
            raise

        self.formal_charges, self.bond_orders, self.non_bonded_electrons = {}, {}, {}

        for v in problem.variables():
            variable_type, variable_substr = v.name.split('_')
            if variable_type == 'C':
                atom_index = int(variable_substr)
                self.formal_charges[atom_index] = MUST_BE_INT(v.varValue)
            elif variable_type == 'B':
                bond_index = int(variable_substr)
                self.bond_orders[bond_reverse_mapping[bond_index]] = MUST_BE_INT(v.varValue)
                pass
            elif variable_type == 'Z':
                pass
            elif variable_type == 'N':
                atom_index = int(variable_substr)
                self.non_bonded_electrons[atom_index] = MUST_BE_INT(v.varValue) * ELECTRON_MULTIPLIER
                if allow_radicals and self.non_bonded_electrons[atom_index] % 2 == 1:
                    stderr.write('Warning: Radical molecule...')
            else:
                raise Exception('Unknown variable type: {0}'.format(variable_type))
        write_to_debug(debug, 'molecule_name', self.name)
        write_to_debug(debug, 'bond_orders:', self.bond_orders)
        write_to_debug(debug, 'formal_charges', self.formal_charges)
        write_to_debug(debug, 'non_bonded_electrons', self.non_bonded_electrons)
        self.assign_aromatic_bonds()

    def rings(self) -> List[Any]:
        from networkx import Graph, cycle_basis

        graph = Graph()
        graph.add_nodes_from([atom.index for atom in self.sorted_atoms()])
        graph.add_edges_from([tuple(bond) for bond in self.bonds])

        return list(map(tuple, cycle_basis(graph)))

    def assign_aromatic_bonds(self):
        rings = self.rings()

        bonds_for_rings = {
            ring: bonds_for_ring(ring)
            for ring in rings
        }

        bond_orders_for_rings = {
            ring: [self.bond_orders[bond] for bond in bonds]
            for (ring, bonds) in bonds_for_rings.items()
        }

        def is_hucklel_compatible(bond_orders: List[int]):
            return (sum([2 for bond_order in bond_orders if bond_order == 2]) - 2) % 4 == 0

        def is_bond_sequence_aromatic(bond_orders: List[int]) -> bool:
            cyclic_bond_orders = [bond_orders[-1]] + list(bond_orders) + [bond_orders[0]]
            for pair_of_bond_orders in zip(cyclic_bond_orders, [cyclic_bond_orders[-1]] + cyclic_bond_orders[:-1]):
                if set(pair_of_bond_orders) == {1, 2}:
                    continue
                else:
                    return False
            else:
                return is_hucklel_compatible(bond_orders)

        self.aromatic_rings, self.aromatic_bonds = set(), set()
        for ring in rings:
            if is_bond_sequence_aromatic(bond_orders_for_rings[ring]):
                self.aromatic_rings.add(ring)
                for bond in bonds_for_rings[ring]:
                    self.aromatic_bonds.add(bond)

    def renumber_atoms(self) -> None:
        atom_mapping = {
            atom_index: i
            for (i, atom_index) in enumerate(self.atoms.keys(), start=1)
        }

        remap_atom = lambda atom_index: atom_mapping[atom_index]
        remap_bond = lambda bond: frozenset(map(remap_atom, bond))

        self.atoms = {
            remap_atom(atom.index): atom._replace(index=remap_atom(atom.index))
            for atom in self.atoms.values()
        }

        self.previously_uncapped = {
            remap_atom(atom_index)
            for atom_index in self.previously_uncapped
        }

        self.bonds = {
            remap_bond(bond)
            for bond in self.bonds
        }

        if self.formal_charges is not None:
            self.formal_charges = {
                remap_atom(atom_index): charge
                for (atom_index, charge) in self.formal_charges.items()
            }

        if self.bond_orders is not None:
            self.bond_orders = {
                remap_bond(bond): bond_order
                for (bond, bond_order) in self.bond_orders.items()
            }

        if self.aromatic_bonds is not None:
            self.aromatic_bonds = {
                remap_bond(bond)
                for bond in self.aromatic_bonds
            }

        if self.phenyl_atoms is not None:
            self.phenyl_atoms = {
                remap_atom(atom_index)
                for atom_index in self.phenyl_atoms
            }

        return None

    def remove_atoms_with_predicate(self, predicate: Callable[[Atom], bool], reset_valences: bool = True, reset_capped: bool = True) -> None:
        deleted_atom_ids = {
            atom.index
            for atom in self.atoms.values()
            if predicate(atom)
        }

        atoms_connected_to_deleted_atoms = reduce(
            lambda acc, e: acc | e,
            [
                bond
                for bond in self.bonds
                if len(bond & deleted_atom_ids) != 0
            ],
            set(),
        ) - deleted_atom_ids

        self.atoms = {
            atom.index: atom if atom.index not in atoms_connected_to_deleted_atoms else atom._replace(**{'valence':None} if reset_valences else {}, **{'capped': False} if reset_capped else {})
            for atom in self.atoms.values()
            if atom.index not in deleted_atom_ids
        }

        self.bonds = {
            bond
            for bond in self.bonds
            if bond & deleted_atom_ids == set()
        }

        if self.formal_charges is not None:
            self.formal_charges = {
                atom_index: charge
                for (atom_index, charge) in self.formal_charges.items()
                if atom_index not in deleted_atom_ids
            }

        if self.bond_orders is not None:
            self.bond_orders = {
                bond: bond_order
                for (bond, bond_order) in self.bond_orders.items()
                if bond in self.bonds
            }

        if self.aromatic_bonds is not None:
            self.aromatic_bonds = {
                bond
                for bond in self.aromatic_bonds
                if bond in self.bonds
            }

        if self.non_bonded_electrons is not None:
            self.non_bonded_electrons = {
                atom_index: non_bonded_electrons
                for (atom_index, non_bonded_electrons) in self.non_bonded_electrons.items()
                if atom_index not in deleted_atom_ids
            }

        return self.renumber_atoms()

    def remove_all_hydrogens(self, mark_all_uncapped: bool = False) -> None:
        # Store phenyl rings atoms before removing their hydrogens
        self.phenyl_atoms = atom_indices_in_phenyl_rings_for(self)

        self.remove_atoms_with_predicate(
            lambda atom: atom.element == 'H',
        )

        if mark_all_uncapped:
            self.atoms = {
                atom_index: atom._replace(capped=False)
                for (atom_index, atom) in self.atoms.items()
            }

    def remove_united_hydrogens(self, **kwargs: Dict[str, Any]) -> None:
        first_neighbours_ids = {
            atom.index: reduce(
                lambda acc, e: acc | e,
                [bond for bond in self.bonds if atom.index in bond],
                frozenset(),
            ) - {atom.index}
            for atom in self.atoms.values()
        }

        return self.remove_atoms_with_predicate(
            lambda atom: all(
                [
                    atom.element == 'H',
                    atom.valence == 1,
                    {self.atoms[neighbour_index].element for neighbour_index in first_neighbours_ids[atom.index]} == {'C'},
                    {self.formal_charges[neighbour_index] for neighbour_index in first_neighbours_ids[atom.index]} == {0},
                ]
            ),
            **kwargs,
        )

    def get_all_tautomers(self, *args, **kwargs):
        return get_all_tautomers(self, *args, **kwargs)

    def get_all_tautomers_naive(self, *args, **kwargs):
        return get_all_tautomers_naive(self, *args, **kwargs)

    def get_best_capped_molecule_with_ILP(self, *args, **kwargs):
        return get_best_capped_molecule_with_ILP(self, *args, **kwargs)

    def get_best_capped_molecule(self, *args, **kwargs):
        return get_best_capped_molecule(self, *args, **kwargs)

Uncapped_Molecule = Molecule

def validated_atoms_dict(atoms: Dict[int, Atom]) -> Dict[int, Atom]:
    assert {atom.element for atom in atoms.values()} <= ALL_ELEMENTS, 'Unsupported elements: {0}'.format({atom.element for atom in atoms.values()} - ALL_ELEMENTS)

    for (atom_index, atom) in atoms.items():
        assert atom_index == atom.index, 'Invalid index={0} for atom={1}'.format(atom_index, atom)

    return atoms

def validate_bond_dict(atoms: Dict[int, Atom], bonds: Set[Bond]) -> None:
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

def bonds_for_pdb_line(pdb_line: str) -> List[FrozenSet[int]]:
    atom_index, *other_atom_indices = map(int, pdb_line.split()[1:])
    return {
        frozenset((atom_index, other_atom_index))
        for other_atom_index in other_atom_indices
    }

def molecule_from_pdb_str(pdb_str: str, **kwargs: Dict[str, Any]) -> Molecule:
    from chemistry_helpers.pdb import pdb_atoms_in, is_pdb_connect_line

    bonds = reduce(
        lambda acc, e: acc | e,
        [
            bonds_for_pdb_line(line)
            for line in pdb_str.splitlines()
            if is_pdb_connect_line(line)
        ],
        set(),
    )

    return Molecule(
        {
            atom.index: Atom(
                index=atom.index,
                element=atom.element.upper(),
                valence=len([1 for bond in bonds if atom.index in bond]),
                capped=True,
                coordinates=atom.coordinates,
            )
            for atom in pdb_atoms_in(pdb_str)
        },
        bonds,
        **kwargs,
    )

def molecule_from_mol2_str(mol2_str: str, **kwargs: Dict[str, Any]) -> Molecule:
    assert mol2_str.startswith('@<TRIPOS>MOLECULE'), 'Error: MOL2 file does not start with "@<TRIPOS>MOLECULE"'
    assert mol2_str.count('@<TRIPOS>MOLECULE') == 1, 'Only one molecule at a time'

    def atom_for_atom_line(line: str) -> Tuple[Atom, float]:
        index_str, name_str, x, y, z, sybil_atom_type, _, _, partial_charge = line.split()
        element, *_ = sybil_atom_type.split('.')

        return (
            Atom(
                index=int(index_str),
                element=element,
                valence=None,
                coordinates=(x, y, z),
                capped=True,
            ),
            float(partial_charge),
        )

    AROMATIC_BOND, AMIDE_BOND = 'ar', 'am'
    def bond_for_atom_line(line: str) -> Tuple[Bond, int]:
        bond_label, atom_id_1, atom_id_2, bond_order_str = line.split()
        if bond_order_str == AROMATIC_BOND:
            bond_order = 1.5
        elif bond_order_str == AMIDE_BOND:
            bond_order = 1
        else:
            bond_order = int(bond_order_str)
        return (frozenset([int(atom_id_1), int(atom_id_2)]), bond_order)

    read_lines, atoms, bonds = False, [], []
    for (i, line) in enumerate(mol2_str.splitlines()):
        if i == 1:
            molecule_name = line
        elif line.startswith('@<TRIPOS>ATOM'):
            container, line_reading_fct, read_lines = atoms, atom_for_atom_line, True
        elif line.startswith('@<TRIPOS>BOND'):
            container, line_reading_fct, read_lines = bonds, bond_for_atom_line, True
        elif line.startswith('@'):
            read_lines = False
        else:
            if read_lines:
                container.append(line_reading_fct(line))

    total_net_charge = sum(partial_charge for (atom, partial_charge) in atoms)
    assert abs(total_net_charge - round(total_net_charge)) <= 0.01, total_net_charge

    return Molecule(
        [atom for (atom, _) in atoms],
        [bond for (bond, _) in bonds],
        formal_charges={atom.index: round(partial_charge) for (atom, partial_charge) in atoms},
        bond_orders={bond: bond_order for (bond, bond_order) in bonds},
        name=molecule_name,
        net_charge=round(total_net_charge),
    )

if __name__ == '__main__':
    TEST_PDB = '''HEADER    UNCLASSIFIED                            10-Sep-17
    TITLE     ALL ATOM STRUCTURE FOR MOLECULE UNK                                   
    AUTHOR    GROMOS AUTOMATIC TOPOLOGY BUILDER REVISION 2017-07-03 14:53:07
    AUTHOR   2  http://compbio.biosci.uq.edu.au/atb
    HETATM    1   H9 G223    0       3.249   0.792  -0.673  1.00  0.00           H
    HETATM    2   C7 G223    0       2.882   0.023   0.017  1.00  0.00           C
    HETATM    3   H7 G223    0       3.252   0.283   1.017  1.00  0.00           H
    HETATM    4   H8 G223    0       3.319  -0.940  -0.264  1.00  0.00           H
    HETATM    5   C4 G223    0       1.373  -0.026   0.001  1.00  0.00           C
    HETATM    6   N1 G223    0       0.725   1.151   0.064  1.00  0.00           N
    HETATM    7   C3 G223    0      -0.615   1.138   0.067  1.00  0.00           C
    HETATM    8   H4 G223    0      -1.098   2.114   0.119  1.00  0.00           H
    HETATM    9   C2 G223    0      -1.404  -0.016   0.006  1.00  0.00           C
    HETATM   10   C1 G223    0      -2.911   0.059  -0.003  1.00  0.00           C
    HETATM   11   H1 G223    0      -3.319  -0.237  -0.978  1.00  0.00           H
    HETATM   12   H2 G223    0      -3.350  -0.611   0.746  1.00  0.00           H
    HETATM   13   H3 G223    0      -3.259   1.075   0.208  1.00  0.00           H
    HETATM   14   C6 G223    0      -0.716  -1.234  -0.056  1.00  0.00           C
    HETATM   15   H6 G223    0      -1.266  -2.171  -0.103  1.00  0.00           H
    HETATM   16   C5 G223    0       0.677  -1.241  -0.059  1.00  0.00           C
    HETATM   17   H5 G223    0       1.223  -2.179  -0.107  1.00  0.00           H
    CONECT    1    2
    CONECT    2    1    3    4    5
    CONECT    3    2
    CONECT    4    2
    CONECT    5    2    6   16
    CONECT    6    5    7
    CONECT    7    6    8    9
    CONECT    8    7
    CONECT    9    7   10   14
    CONECT   10    9   11   12   13
    CONECT   11   10
    CONECT   12   10
    CONECT   13   10
    CONECT   14    9   15   16
    CONECT   15   14
    CONECT   16    5   14   17
    CONECT   17   16
    END'''
    molecule = molecule_from_pdb_str(TEST_PDB)
    #molecule.assign_bond_orders_and_charges_with_ILP()
    #print(molecule.charges, molecule.non_bonded_electrons, molecule.bond_orders)
    #print(molecule.mol2())

    TEST_MOL2 = '''@<TRIPOS>MOLECULE
AGLYSL01
   10     9    1     0    0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
   1 C1      -1.6234     1.6965     8.8431 C.3      1  AGLY  0.3310
   2 C2      -1.5438     0.1710     8.8960 C.2      1  AGLY  0.6590
   3 H1      -1.5827     1.7640     6.8094 H        1  AGLY  0.3600
   4 H3      -0.1271     1.8630     7.4736 H        1  AGLY  0.3600
   5 H5      -2.6707     1.9987     8.9343 H        1  AGLY  0.0000
   6 H6      -1.0462     2.1092     9.6756 H        1  AGLY  0.0000
   7 H7      -2.3655     0.3289    10.6732 H        1  AGLY  0.5000
   8 N1      -1.0818     2.2182     7.5784 N.3      1  AGLY -0.9900
   9 O5      -2.0344    -0.3437    10.0478 O.3      1  AGLY -0.6500
  10 O6      -1.0893    -0.5326     8.0048 O.2      1  AGLY -0.5700
@<TRIPOS>BOND
   1    1    2 1 
   2    1    5 1 
   3    1    6 1 
   4    1    8 1 
   5    2    9 1 
   6    2   10 2 
   7    3    8 1 
   8    4    8 1 
   9    7    9 1 
@<TRIPOS>SUBSTRUCTURE
   1  AGLY    1
@<TRIPOS>COMMENT
COMMENT AMMONIUM GLYCINIUM SULFATE (NEUTRON STUDY) PEPSEQ A=1 GLY'''

    molecule = molecule_from_mol2_str(TEST_MOL2)
    print(molecule)
    molecule.write_graph('RAW')
    molecule.assign_bond_orders_and_charges_with_ILP()
    molecule.write_graph('ILP')
