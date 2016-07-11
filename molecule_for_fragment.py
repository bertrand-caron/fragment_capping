from pickle import load
from itertools import groupby, product

from API_client.api import API
from fragment_dihedrals.fragment_dihedral import element_valence_for_atom, on_asc_number_electron_then_asc_valence, NO_VALENCE
from collections import namedtuple
from StringIO import StringIO

DRAW_GRAPHS = True

DEBUG = False

def concat(list_of_lists):
    return reduce(
        lambda acc, e: acc + e,
        list_of_lists,
    )

class Molecule:
    FULL_VALENCES = {
        'C': 4,
        'N': 3,
        'O': 2,
        'H': 1,
        'S': 2,
    }

    def __init__(self, atoms, bonds, name=None):
        self.atoms = atoms
        self.bonds = bonds
        self.name = name

        self.use_neighbour_valences = (True if all([atom['valence'] is not NO_VALENCE for atom in self.atoms.values()]) else False)

    def atom_desc(self, atom):
        if self.use_neighbour_valences:
            return atom['element'] + str(atom['valence'])
        else:
            return atom['element']

    def assert_use_neighbour_valences(self):
        assert self.use_neighbour_valences, 'ERROR: self.use_neighbour_valences is set to False'

    def valence(self, atom_id):
        self.assert_use_neighbour_valences()
        return self.atoms[atom_id]['valence']

    def element(self, atom_id):
        return self.atoms[atom_id]['element']

    def ids(self):
        return self.atoms.keys()

    def __str__(self):
        return 'Molecule; atoms={0}; bonds={1}'.format(self.atoms, self.bonds)

    def add_atom(self, atom):
        highest_id = max(self.atoms.keys())
        atom_id = highest_id + 1
        self.atoms[atom_id] = dict(atom.items() + dict(index=atom_id).items())
        return atom_id

    def capped_molecule_with(self, capping_strategies, atoms_need_capping):
        from copy import deepcopy
        capped_molecule = deepcopy(self)

        for (atom, capping_strategy) in zip(atoms_need_capping, capping_strategies):
            atom_id = atom['index']
            new_atoms, fragment_bonds, new_valences = capping_strategy

            last_used_id = sorted(capped_molecule.ids())[-1]
            new_ids = map(
                lambda (id, _): id + last_used_id + 1,
                enumerate(new_atoms),
            )
            new_bonds = [
                tuple(
                    map(
                        lambda id: atom_id if id == 0 else id + last_used_id,
                        bond,
                    ),
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
                capped_molecule.atoms[new_id] = {
                    'element': new_atom,
                    'valence': new_valence if self.use_neighbour_valences else NO_VALENCE,
                    'index': new_id,
                    'capped': True,
                }
            capped_molecule.atoms[atom_id]['capped'] = True

        assert all([atom['capped'] for atom in capped_molecule.atoms.values()]), 'Some atoms were not capped: {0}'.format(
            [atom for atom in capped_molecule.atoms.values() if not atom['capped']],
        )

        if self.use_neighbour_valences:
            capped_molecule.check_valence()

        try:
            capped_molecule.assign_bond_orders_and_charges()
            return capped_molecule
        except AssertionError as e:
            if DEBUG:
                print 'AssertionError for capped molecule {0}:\n{1}'.format(
                    capped_molecule,
                    str(e),
                )
            return None

    def check_valence(self):
        self.assert_use_neighbour_valences()

        try:
            for atom in self.atoms.values():
                atom_id = atom['index']
                assert atom['valence'] == sum([1 for bond in self.bonds if atom_id in bond]), 'Atom {2}: {0} != {1} (bonds={3})'.format(
                    atom['valence'],
                    sum([1 for bond in self.bonds if atom_id in bond]),
                    atom,
                    [bond for bond in self.bonds if atom_id in bond],
                )
        except:
            print 'ERROR'
            print 'Atoms are: {0}'.format(self.atoms)
            print 'Bonds are: {0}'.format(self.bonds)
            raise

    def get_best_capped_molecule(self):
        Capping_Strategy = namedtuple('Capping_Strategy', 'new_atoms, new_bonds, new_valences')

        NO_CAP = Capping_Strategy((), (), ())
        H_CAP = Capping_Strategy(('H',), ((0, 1),), (1,))
        H2_CAP = Capping_Strategy(('H', 'H'), ((0, 1), (0, 2)), (1, 1))
        H3_CAP = Capping_Strategy(('H', 'H', 'H'), ((0, 1), (0, 2), (0, 3)), (1, 1, 1))
        H_CH2_CAP = Capping_Strategy(('H', 'C', 'H', 'H'), ((0, 1), (0, 2), (2, 3), (2, 4)), (1, 3, 1, 1))
        CH3_CAP = Capping_Strategy(('C', 'H', 'H', 'H'), ((0, 1), (1, 2), (1, 3), (1, 4)), (4, 1, 1, 1))

        INDIVIDUAL_CAPPING_OPTIONS = {
            'H1': (NO_CAP,),
            'O1': (NO_CAP,),
            'O2': (H_CAP,),
            'S1': (NO_CAP,),
            'S2': (H_CAP,),
            'C4': (H3_CAP,),
            'C3': (H2_CAP, H_CH2_CAP),
            'N2': (H_CAP, CH3_CAP,),
            'N3': (H2_CAP,),
            'N4': (H3_CAP,),
        }

        on_first_atom_desc_letter = lambda (atom_desc, capping_strategies): atom_desc[0]

        if not self.use_neighbour_valences:
            'Aggregate capping strategies for a given element.'
            CAPPING_OPTIONS = dict(
                [
                    (element, concat(map(lambda x: x[1], group)))
                    for (element, group) in
                    groupby(
                        sorted(INDIVIDUAL_CAPPING_OPTIONS.items(), key=on_first_atom_desc_letter),
                        key=on_first_atom_desc_letter,
                    )
                ]
            )
        else:
            CAPPING_OPTIONS = INDIVIDUAL_CAPPING_OPTIONS

        if DEBUG:
            print [(key, len(value)) for (key, value) in INDIVIDUAL_CAPPING_OPTIONS.items()]
            print [(key, len(value)) for (key, value) in CAPPING_OPTIONS.items()]

        atoms_need_capping = [atom for atom in self.sorted_atoms() if not atom['capped']]
        capping_schemes = list(
            product(
                *[
                    CAPPING_OPTIONS[self.atom_desc(atom)]
                    for atom in atoms_need_capping
                ]
            ),
        )

        if DEBUG:
            print 'atoms_need_capping: {0}'.format(atoms_need_capping)
            print 'capping_schemes: {0}'.format(capping_schemes)
            print 'capping_options: {0}'.format([
                len(CAPPING_OPTIONS[self.atom_desc(atom)])
                for atom in atoms_need_capping
            ])

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

        print 'Possible capped molecules: {0} ({1}/{2})'.format(
            [(mol.formula(charge=True), mol.net_abs_charges(), mol.double_bonds_fitness()) for mol in possible_capped_molecules],
            len(possible_capped_molecules),
            len(capping_schemes),
        )

        if DRAW_GRAPHS:
            try:
                from py_graphs.pdb import graph_from_pdb
                from py_graphs.moieties import draw_graph
                for (i, molecule) in enumerate(possible_capped_molecules):
                    graph = molecule.graph()
                    draw_graph(
                        graph,
                        fnme=join('graphs' ,'_'.join((self.name, str(i))) + '.png'),
                    )
            except Exception as e:
                print 'ERROR: Could not plot graphs (error was: {0})'.format(
                    str(e),
                )

        best_molecule = possible_capped_molecules[0]
        return best_molecule

    def formula(self, charge=False):
        elements =  [atom['element'] for atom in self.atoms.values()]

        return ''.join(
            map(
                lambda (element, number): element + (str(number) if number > 1 else ''),
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
            )
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

    def n_atoms(self):
        return len(self.atoms)

    def dummy_pdb(self):
        from atb_helpers.pdb import PDB_TEMPLATE, CONECT_TEMPLATE, pdb_conect_line
        io = StringIO()

        def shift_one(x):
            return x + 1

        for (i, atom) in enumerate(sorted(self.atoms.values(), key=lambda atom: atom['index'])):
            print >> io, PDB_TEMPLATE.format(
                'HETATM',
                shift_one(i),
                'D',
                'R',
                '',
                shift_one(i),
                1. * i,
                0.,
                0.,
                '',
                '',
                atom['element'].title(),
                '',
            )

        for (i, atom) in enumerate(sorted(self.atoms.values(), key=lambda atom: atom['index'])):
            print >> io, pdb_conect_line(
                [shift_one(atom['index'])]
                +
                [shift_one(bond[0] if bond[1] == atom['index'] else bond[1]) for bond in self.bonds if atom['index'] in bond]
            )

        return io.getvalue()

    def representation(self, out_format):
        assert out_format in ('smiles', 'inchi'), 'Wrong representation format: {0}'.format(out_format)

        from atb_helpers.babel import babel_output
        return babel_output(
            self.dummy_pdb(),
            in_format='pdb',
            out_format=out_format,
            dont_add_H=True,
        )

    def smiles(self):
        return self.representation('smiles')

    def inchi(self):
        return self.representation('inchi')

    def graph(self):
        try:
            from graph_tool.all import Graph
        except:
            return None

        g = Graph(directed=False)

        vertex_types = g.new_vertex_property("string")
        g.vertex_properties['type'] = vertex_types

        vertices = []
        for atom_index in sorted(self.atoms.keys()):
            v = g.add_vertex()
            vertex_types[v] = '{element}{valence}'.format(
                element=self.atoms[atom_index]['element'],
                valence=self.atoms[atom_index]['valence'] if self.use_neighbour_valences else '',
            )
            vertices.append(v)

        for (i, j) in self.bonds:
            g.add_edge(vertices[i], vertices[j])

        return g

    def sorted_atoms(self):
        return [atom  for (atom_id, atom) in sorted(self.atoms.items())]

    def sorted_atom_ids(self):
        return [atom_id  for (atom_id, atom) in sorted(self.atoms.items())]

    def assign_bond_orders_and_charges(self):
        POSSIBLE_BOND_ORDERS = {
            'S': (1, 2,),
            'C': (1, 2,),
            'H': (1,),
            'O': (1, 2,),
            'N': (1, 2,),
        }

        POSSIBLE_CHARGES = {
            'S': (0,),
            'C': (0,),
            'H': (0,),
            'O': (0, -1,),
            'N': (0, +1,),
        }

        possible_bond_orders_lists = list(
            product(
                *[
                    set.intersection(set(POSSIBLE_BOND_ORDERS[element_1]), set(POSSIBLE_BOND_ORDERS[element_2]))
                    for (element_1, element_2) in
                    map(
                        lambda bond: (self.atoms[bond[0]]['element'], self.atoms[bond[1]]['element']),
                        self.bonds,
                    )
                ]
            ),
        )

        assert len(possible_bond_orders_lists) >= 1, 'No possible bond orders found'

        possible_charges_dicts = map(
            lambda charges: dict(zip(self.sorted_atom_ids(), charges)),
            product(
                *[
                    POSSIBLE_CHARGES[atom['element']]
                    for atom in self.sorted_atoms()
                ]
            ),
        )

        assert len(possible_charges_dicts) >= 1, 'No possible charges assignment found'

        possible_bond_orders_and_charges = list(product(possible_bond_orders_lists, possible_charges_dicts))

        acceptable_bond_orders_and_charges = sorted(
            [
                (
                    zip(self.bonds, bond_orders),
                    charges,
                )
                for (bond_orders, charges) in possible_bond_orders_and_charges
                if self.is_valid(bond_orders, charges)
            ],
            key=lambda (_, charges): sum(map(abs, charges.values())),
        )

        if DEBUG:
            if len(acceptable_bond_orders_and_charges) != 1:
                print 'acceptable_bond_orders_and_charges: {0}'.format(acceptable_bond_orders_and_charges)

        assert len(acceptable_bond_orders_and_charges) >= 1, 'No valid bond_orders and charges found amongst {0} tried.'.format(len(possible_bond_orders_and_charges))

        self.bond_orders, self.charges = acceptable_bond_orders_and_charges[0]

    def netcharge(self):
        try:
            return sum(self.charges.values())
        except:
            raise Exception('Assign charges and bond_orders first.')

    def net_abs_charges(self):
        try:
            return sum(map(abs, self.charges.values()))
        except:
            raise Exception('Assign charges and bond_orders first.')

    def is_valid(self, bond_orders, charges):
        assert len(self.bonds) == len(bond_orders), 'Unmatched bonds and bond_orders: {0} != {1}'.format(
            len(self.bonds),
            len(bond_orders),
        )

        on_atom_id = lambda (atom_id, bond_order): atom_id
        on_bond_order = lambda (atom_id, bond_order): bond_order

        valid = all(
            [
                sum(map(on_bond_order, group)) == Molecule.FULL_VALENCES[self.atoms[atom_id]['element']] + charges[atom_id]
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

    def double_bonds_fitness(self):
        '''Sorted ASC (low fitness is better)'''

        BEST_DOUBLE_BONDS = (
            # From best, to worst
            'CO',
            'CN',
            'CC',
        )

        grouped_double_bonds = dict([
            (key, len(list(group)))
            for (key, group) in
            groupby(
                sorted(
                    [
                        ''.join(
                            sorted([self.atoms[atom_id]['element'] for atom_id in bond]),
                        )
                        for (bond, bond_order) in self.bond_orders
                        if bond_order == 2
                    ]
                ),
            )
        ])
        return tuple([(- grouped_double_bonds[double_bond_type] if double_bond_type in grouped_double_bonds else 0) for double_bond_type in BEST_DOUBLE_BONDS])


api = API(
    host='http://scmb-atb.biosci.uq.edu.au/atb-uqbcaron', #'https://atb.uq.edu.au',
    debug=False,
    api_format='pickle',
)

def truncated_molecule(molecule):
    return dict(
        n_atoms=molecule.n_atoms,
        num_dihedral_fragments=len(molecule.dihedral_fragments),
        molid=molecule.molid,
        formula=molecule.formula,
    )

def reduce_iterables(iterables):
    reduced_iterable = reduce(
        lambda acc, e: acc + e,
        iterables
    )
    return reduced_iterable

assert reduce_iterables([[1], [2], [3]]) == [1, 2, 3], reduce_iterables([[1]], [[2]], [[3]])

def best_capped_molecule_for_dihedral_fragment(fragment_str):
    if fragment_str.count('|') == 3:
        neighbours_1, atom_2, atom_3, neighbours_4 = fragment_str.split('|')
        cycles = []
        neighbours_1, neighbours_4 = neighbours_1.split(','), neighbours_4.split(',')
    elif fragment_str.count('|') == 4:
        neighbours_1, atom_2, atom_3, neighbours_4, cycles = fragment_str.split('|')
        neighbours_1, neighbours_4, cycles = neighbours_1.split(','), neighbours_4.split(','), cycles.split(',')
    else:
        raise Exception('Invalid fragment_str: "{0}"'.format(fragment_str))

    print cycles

    ids = [n for (n, _) in enumerate(neighbours_1 + [atom_2, atom_3] + neighbours_4)]

    neighbours_id_1, atom_id_2, atom_id_3, neighbours_id_4 = ids[0:len(neighbours_1)], ids[len(neighbours_1)], ids[len(neighbours_1) + 1], ids[len(neighbours_1) + 2:]
    #print ids
    #print neighbours_id_1, atom_2, atom_3, neighbours_id_4
    CENTRAL_BOND = (atom_id_2, atom_id_3)

    elements = dict(
        zip(
            ids,
            [element_valence_for_atom(neighbour)[0] for neighbour in neighbours_1] + [atom_2, atom_3] + [element_valence_for_atom(neighbour)[0] for neighbour in neighbours_4],
        ),
    )

    valences = dict(
        zip(
            ids,
            [element_valence_for_atom(neighbour)[1] for neighbour in neighbours_1] + [len(neighbours_1) + 1, len(neighbours_4) + 1] + [element_valence_for_atom(neighbour)[1] for neighbour in neighbours_4],
        ),
    )

    bonds = [(neighbour_id, atom_id_2) for neighbour_id in neighbours_id_1] + [CENTRAL_BOND] + [(atom_id_3, neighbour_id) for neighbour_id in neighbours_id_4]

    m = Molecule(
        dict(
            zip(
                ids,
                [
                    {
                        'valence': valences[atom_id],
                        'element': elements[atom_id],
                        'index':atom_id,
                        'capped': (atom_id not in (neighbours_id_1 + neighbours_id_4)),
                    }
                    for atom_id in ids],

            )
        ),
        bonds,
        name=fragment_str.replace('|', '_'),
    )

    print m
    for (i, n, j) in map(lambda cycle: map(int, cycle), cycles):
        i_id, j_id = neighbours_id_1[i], neighbours_id_4[j]
        if n == 0:
            # i and j are actually the same atoms
            del m.atoms[j_id]
            replace_j_by_i = lambda x: i_id if x == j_id else x
            m.bonds = [
                map(replace_j_by_i, bond)
                for bond in m.bonds
            ]
        else:
            NEW_ATOM_ID = -1
            NEW_ATOM = {
                'valence': NO_VALENCE,
                'element': 'C',
                'index': NEW_ATOM_ID, # This will get overwritten by Molecule.add_atom
                'capped': False,
            }
            atom_chain_id = [i_id] + [m.add_atom(NEW_ATOM) for i in range(n - 1)] + [j_id]
            new_bonds = zip(atom_chain_id[:-1], atom_chain_id[1:])
            m.bonds += new_bonds

    m = m.get_best_capped_molecule()
    print m
    return m

def cap_fragment(fragment, count=None, i=None, fragments=None):
    if all([x is not None for x in (count, i, fragments)]):
        print 'Running fragment {0}/{1} (count={2}): "{3}"'.format(
            i + 1,
            len(fragments),
            count,
            fragment,
        )

    molecule = best_capped_molecule_for_dihedral_fragment(fragment)
    #print molecule.inchi()
    #print molecule.dummy_pdb()

    api_response = api.Molecules.structure_search(
        netcharge='*',
        structure_format='pdb',
        structure=molecule.dummy_pdb(),
        return_type='molecules',
    )

    molecules = api_response['matches']

    if molecules:
        print [atb_molecule['molid'] for atb_molecule in molecules]
        best_molid = sorted(
            molecules,
            key=lambda atb_molecule: int(atb_molecule['molid']),
        )[0]['molid']

        best_molecule = api.Molecules.molid(
            molid=best_molid,
        )

        try:
            assert best_molecule.is_finished, 'Molecule is still running'
            #assert fragment in best_molecule.dihedral_fragments, 'Dihedral fragment not found in molecule. Maybe it is still running ?'
            if fragment != 'N,C,H|C|C|O,O':
                assert set([atb_molecule['InChI'] for atb_molecule in molecules]) == set([best_molecule.inchi]), 'Several molecules with different InChI have been found: {0}'.format(
                    set([atb_molecule['InChI'] for atb_molecule in molecules]),
                )

        except AssertionError as e:
            print e
            best_molid = None

    else:
        print 'Capped fragment not found in ATB.'
        print molecule.formula()
        print molecule.dummy_pdb()
        best_molid = None

    print

    safe_fragment_name = fragment.replace('|', '_')

    with open('pdbs/{fragment}.pdb'.format(fragment=safe_fragment_name), 'w') as fh:
        fh.write(molecule.dummy_pdb())

    return best_molid

def get_matches():
    matches = [
        (fragment, cap_fragment(fragment, count=count, i=i, fragments=protein_fragments))
        for (i, (fragment, count)) in
        enumerate(protein_fragments)
    ]
    return matches

def parse_args():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--only-id', type=int, help='Rerun a single fragment')

    return parser.parse_args()

def generate_collage():

    matches = cached(get_matches, (), {})
    counts = dict(protein_fragments)

    for (fragment, molid) in matches:
        if molid:
            print 'python3 test.py --download {molid} --submit --dihedral-fragment "{dihedral_fragment}"'.format(
                molid=molid,
                dihedral_fragment=fragment,
            )

    print matches
    print 'INFO: Assigned {0}/{1} molecules (missing_ids = {2})'.format(
        len(filter(lambda (_, molid): molid is not None, matches)),
        len(matches),
        [i for (i, (fragment, molid)) in enumerate(matches) if molid is None],
    )

    def png_file_for(molid):
        PNG_DIR = 'pngs'

        png_file = join(
            PNG_DIR,
            '{molid}.png'.format(molid=molid),
        )

        if not exists(png_file):
            svg2png(
                url='https://atb.uq.edu.au/img2D/{molid}_thumb.svg'.format(molid=molid),
                write_to=png_file,
            )
        return png_file

    png_files = dict([(molid, png_file_for(molid)) for (_, molid) in matches if molid is not None])

    def figure_collage():
        import matplotlib.pyplot as p
        from PIL import Image
        subplot_dim = int(ceil(sqrt(len(matches))))

        fig, axarr = p.subplots(*[subplot_dim]*2, figsize=(30, 15))

        def indices_for_fig(n):
            return ((n // subplot_dim), n - (n // subplot_dim) * subplot_dim)

        for (n, (fragment, molid)) in enumerate(matches):
            if molid in png_files:
                image = Image.open(png_files[molid])
                axarr[indices_for_fig(n)].imshow(image)
            axarr[indices_for_fig(n)].set_title(
                fragment + (' (id={1}, molid={0})'.format(molid, n)),
                fontsize=11,
                fontname='Andale Mono',
            )
            axarr[indices_for_fig(n)].set_axis_off()

        for n in range(len(matches), subplot_dim**2):
            axarr[indices_for_fig(n)].set_axis_off()

        p.tight_layout()
        p.show()
        fig.savefig('collage.png')

    figure_collage()

REMOVE_VALENCES = True

def get_protein_fragments():
    with open('cache/protein_fragments.pickle') as fh:
        protein_fragments = load(fh)

    if REMOVE_VALENCES:
        from re import sub
        protein_fragments = [(sub('[0-9]', '', fragment), count) for (fragment, count) in protein_fragments]

    return protein_fragments

def test_cyclic_fragments():
    #print cap_fragment('C,H,H|C|C|C,H,H|000')
    print cap_fragment('C,H,H|C|C|C,H,H|010')
    print cap_fragment('C,H,H|C|C|C,H,H|020')
    print cap_fragment('C,H,H|C|C|C,H,H|030')

if __name__ == '__main__':
    from cache import cached
    from cairosvg import svg2png
    from os.path import join, exists
    from math import sqrt, ceil

    args = parse_args()
    protein_fragments = get_protein_fragments()

    if args.only_id:
        print cap_fragment(
            protein_fragments[args.only_id][0],
            i=0,
            count='unknown',
            fragments=(True,),
        )
    else:
        generate_collage()

