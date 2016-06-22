from pickle import load
from itertools import groupby, product

from API_client.api import API
from fragment_dihedrals.fragment_dihedral import element_valence_for_atom, on_asc_number_electron_then_asc_valence
from collections import namedtuple
from StringIO import StringIO

DRAW_GRAPHS = False

def atom_desc(atom):
    return atom['element'] + str(atom['valence'])

class Molecule:
    FULL_VALENCES = {
        'C': 4,
        'N': 3,
        'O': 2,
        'H': 1,
        'S': 2,
    }

    def __init__(self, atoms, bonds):
        self.atoms = atoms
        self.bonds = bonds

    def cap_atom(atom_id):
        assert atom_id in self.atoms

    def valence(self, atom_id):
        return self.atoms[atom_id]['valence']

    def element(self, atom_id):
        return self.atoms[atom_id]['element']

    def ids(self):
        return self.atoms.keys()

    def __str__(self):
        return 'Molecule; atoms={0}; bonds={1}'.format(self.atoms, self.bonds)

    def capped_molecule_with(self, capping_scheme, atoms_need_capping):
        from copy import deepcopy
        capped_molecule = deepcopy(self)

        for (atom, capping_strategy) in zip(atoms_need_capping, capping_scheme):
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

            assert len(new_ids) == len(new_atoms) == len(new_valences)
            for (new_id, new_atom, new_valence) in zip(new_ids, new_atoms, new_valences):
                capped_molecule.atoms[new_id] = {
                    'element': new_atom,
                    'valence': new_valence,
                    'index': new_id,
                    'capped': True,
                }

        try:
            capped_molecule.check_valence()
            return capped_molecule
        except:
            return None

    def check_valence(self):
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

    def get_capped_molecule(self):
        Capping_Strategy = namedtuple('Capping_Strategy', 'new_atoms, new_bonds, new_valences')

        NO_CAP = Capping_Strategy((), (), ())
        H_CAP = Capping_Strategy(('H',), ((0, 1),), (1,))
        H2_CAP = Capping_Strategy(('H', 'H'), ((0, 1), (0, 2)), (1, 1))
        H3_CAP = Capping_Strategy(('H', 'H', 'H'), ((0, 1), (0, 2), (0, 3)), (1, 1, 1))
        H_CH2_CAP = Capping_Strategy(('H', 'C', 'H', 'H'), ((0, 1), (0, 2), (2, 3), (2, 4)), (1, 4, 1, 1))
        CH3_CAP = Capping_Strategy(('C', 'H', 'H', 'H'), ((0, 1), (1, 2), (1, 3), (1, 4)), (4, 1, 1, 1))

        CAPPING_OPTIONS = {
            'H1': (NO_CAP,),
            'O1': (NO_CAP,),
            'O2': (H_CAP,),
            'S1': (NO_CAP,),
            'S2': (H_CAP,),
            'C4': (H3_CAP,),
            'C3': (H2_CAP, H_CH2_CAP),
            'N2': (CH3_CAP,),
            'N3': (H2_CAP,),
            'N4': (H3_CAP,),
        }

        atoms_need_capping = [atom for atom in self.sorted_atoms() if not atom['capped']]
        capping_schemes = product(
            *[
                CAPPING_OPTIONS[atom_desc(atom)]
                for atom in atoms_need_capping
            ]
        )

        possible_capped_molecules = sorted(
            filter(
                lambda mol: mol is not None,
                [
                    self.capped_molecule_with(capping_scheme, atoms_need_capping)
                    for capping_scheme in capping_schemes
                ],
            ),
            key=lambda mol: mol.n_atoms(),
        )

        print 'Possible capped molecules: {0}'.format([mol.formula() for mol in possible_capped_molecules])

        best_molecule = possible_capped_molecules[0]
        return best_molecule

    def formula(self):
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
            ),
        )

    def n_atoms(self):
        return len(self.atoms)

    def dummy_pdb(self):
        from atb_helpers.pdb import PDB_TEMPLATE
        io = StringIO()

        for (i, atom) in enumerate(sorted(self.atoms.values(), key=lambda atom: atom['index'])):
            print >> io, PDB_TEMPLATE.format(
                'HETATM',
                i,
                'D',
                'R',
                '',
                i,
                0.,
                0.,
                0.,
                '',
                '',
                atom['element'].title(),
                '',
            )

        for bond in self.bonds:
            print >> io, ' '.join(['CONECT'] + [str(id) for id in bond])

        return io.getvalue()

    def representation(self, out_format):
        assert out_format in ('smiles', 'inchi')

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
            vertex_types[v] = '{element}{valence}'.format(**self.atoms[atom_index])
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
            'N': (0,),
        }

        possible_bond_orders_lists = product(
            *[
                set.intersection(set(POSSIBLE_BOND_ORDERS[element_1]), set(POSSIBLE_BOND_ORDERS[element_2]))
                for (element_1, element_2) in
                map(
                    lambda bond: (self.atoms[bond[0]]['element'], self.atoms[bond[1]]['element']),
                    self.bonds,
                )
            ]
        )

        possible_charges_dicts = map(
            lambda charges: dict(zip(self.sorted_atom_ids(), charges)),
            product(
                *[
                    POSSIBLE_CHARGES[atom['element']]
                    for atom in self.sorted_atoms()
                ]
            ),
        )

        acceptable_bond_orders_and_charges = [
            (
                zip(self.bonds, bond_orders),
                charges,
            )
            for (bond_orders, charges) in product(possible_bond_orders_lists, possible_charges_dicts)
            if self.is_valid(bond_orders, charges)
        ]

        if len(acceptable_bond_orders_and_charges) != 1:
            print acceptable_bond_orders_and_charges

    def is_valid(self, bond_orders, charges):
        assert len(self.bonds) == len(bond_orders)

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

api = API(
    host='http://scmb-atb.biosci.uq.edu.au/atb-uqbcaron', #'https://atb.uq.edu.au',
    debug=False,
    api_format='pickle',
)

with open('cache/protein_fragments.pickle') as fh:
    protein_fragments = load(fh)

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

def molecule_for_capped_dihedral_fragment(fragment):
    assert fragment.count('|') == 3
    neighbours_1, atom_2, atom_3, neighbours_4 = fragment.split('|')
    neighbours_1, neighbours_4 = neighbours_1.split(','), neighbours_4.split(',')

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
    )

    m = m.get_capped_molecule()
    m.assign_bond_orders_and_charges()
    return m

def get_matches():
    matches = {}

    for (i, (fragment, count)) in enumerate(protein_fragments):
        print 'Running fragment {0}/{1} (count={2}): "{3}"'.format(i + 1, len(protein_fragments), count, fragment)
        molecule = molecule_for_capped_dihedral_fragment(fragment)
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
                assert fragment in best_molecule.dihedral_fragments, 'Dihedral fragment not found in molecule. Maybe it is still running ?'
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

        matches[fragment] = best_molid

        safe_fragment_name = fragment=fragment.replace('|', '_')

        if DRAW_GRAPHS:
            from py_graphs.pdb import graph_from_pdb
            from py_graphs.moieties import draw_graph
            graph = molecule.graph()
            draw_graph(
                graph,
                fnme='graphs/{fragment}.png'.format(
                    fragment=safe_fragment_name,
                ),
            )

        with open('pdbs/{fragment}.pdb'.format(fragment=safe_fragment_name), 'w') as fh:
            fh.write(molecule.dummy_pdb())

        print


    print dict(
        [
            (fragment, matches[fragment])
            for (fragment, count) in
            sorted(
                protein_fragments,
                key=lambda (fragment, count): count,
                reverse=True,
            )
        ],
    )

    for (fragment, molid) in matches.items():
        if molid:
            print 'python3 test.py --download {molid} --submit --dihedral-fragment "{dihedral_fragment}"'.format(
                molid=molid,
                dihedral_fragment=fragment,
            )
    return matches

if __name__ == '__main__':
    from cache import cached
    from cairosvg import svg2png
    from os.path import join, exists
    from math import sqrt, ceil

    matches = cached(get_matches, (), {})
    print matches

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

    png_files = dict([(molid, png_file_for(molid)) for molid in matches.values() if molid is not None])

    def figure_collage():
        import matplotlib.pyplot as p
        from PIL import Image
        subplot_dim = int(ceil(sqrt(len(matches))))

        fig, axarr = p.subplots(*[subplot_dim]*2, figsize=(30, 15))

        def indices_for_fig(n):
            return ((n // subplot_dim), n - (n // subplot_dim) * subplot_dim)

        for (n, (fragment, molid)) in enumerate(sorted(matches.items())):
            if molid in png_files:
                image = Image.open(png_files[molid])
                axarr[indices_for_fig(n)].imshow(image)
            axarr[indices_for_fig(n)].set_title(
                fragment + (' (molid={0})'.format(molid) if molid in png_files else ''),
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
