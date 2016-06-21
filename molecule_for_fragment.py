from pickle import load
from itertools import groupby, product

from API_client.api import API
from fragment_dihedrals.fragment_dihedral import element_valence_for_atom, on_asc_number_electron_then_asc_valence
from collections import namedtuple
from StringIO import StringIO

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

    def cap_neighbour(self, atom_id):
        atom_desc = self.element(atom_id) + str(self.valence(atom_id))

        if atom_desc == 'H1':
            new_atoms = []
            fragment_bonds = []
            new_valences = []
        elif atom_desc == 'C4':
            new_atoms = ['H'] * 3
            fragment_bonds = [(0, 1), (0, 2), (0, 3)]
            new_valences = [1] * 3
        elif atom_desc == 'C3':
            new_atoms = ['H']* 1 + ['C'] * 1 + ['H'] * 2
            fragment_bonds = [(0, 1), (0, 2), (2, 3), (2, 4)]
            new_valences = [1, 3, 1, 1]
        elif atom_desc == 'N4':
            new_atoms = ['H'] * 3
            fragment_bonds = [(0, 1), (0, 2), (0, 3)]
            new_valences = [1] * 3
        elif atom_desc == 'N3':
            new_atoms = ['H'] * 2
            fragment_bonds = [(0, 1), (0, 2)]
            new_valences = [1] * 2
        elif atom_desc == 'N2':
            new_atoms = ['C'] * 1 + ['H'] * 3
            fragment_bonds = [(0, 1), (1, 2), (1, 3), (1, 4)]
            new_valences = [4, 1, 1, 1]
        elif atom_desc == 'O2':
            new_atoms = ['H'] * 1
            fragment_bonds = [(0, 1)]
            new_valences = [1] * 1
        elif atom_desc == 'O1':
            new_atoms = []
            fragment_bonds = []
            new_valences = []
        elif atom_desc == 'S2':
            new_atoms = ['H'] * 1
            fragment_bonds = [(0, 1)]
            new_valences = [1] * 1
        else:
            raise Exception('No rules to cap fragment {0}'.format(atom_desc))

        last_used_id = sorted(self.ids())[-1]
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

        self.bonds += new_bonds

        assert len(new_ids) == len(new_atoms) == len(new_valences)
        for (new_id, new_atom, new_valence) in zip(new_ids, new_atoms, new_valences):
            self.atoms[new_id] = {
                'element': new_atom,
                'valence': new_valence,
                'index': new_id,
            }

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
            print self.atoms
            print self.bonds
            raise

    def cap_molecule(self, neighbours_id_1, neighbours_id_4):
        for neighbour_id in (neighbours_id_1 + neighbours_id_4):
            self.cap_neighbour(neighbour_id)


        if False:
            print map(
                lambda (id_1, id_2): (self.element(id_1), self.element(id_2)),
                self.bonds
            )
        self.check_valence()
        return self.formula()

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

    def assign_bond_orders(self):
        POSSIBLE_BOND_ORDER = (1, 2)

        possible_bond_orders_lists = product(POSSIBLE_BOND_ORDER, repeat=len(self.bonds))
        acceptable_bond_orders = [zip(self.bonds, bond_orders) for bond_orders in possible_bond_orders_lists if self.is_valid(bond_orders)]

        assert len(acceptable_bond_orders) == 1

    def is_valid(self, bond_orders):
        assert len(self.bonds) == len(bond_orders)

        on_atom_id = lambda (atom_id, bond_order): atom_id
        on_bond_order = lambda (atom_id, bond_order): bond_order

        valid = all(
            [
                sum(map(on_bond_order, group)) == self.atoms[atom_id]['valence']
                for (atom_id, group) in
                groupby(
                    sorted(
                        reduce(
                            lambda acc, e: acc + e,
                            [((atom_id_1, bond_order), (atom_id_2, bond_order)) for ((atom_id_1, atom_id_2), bond_order) in zip(self.bonds, bond_orders)],
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
                [{'valence': valences[atom_id], 'element': elements[atom_id], 'index':atom_id,} for atom_id in ids],

            )
        ),
        bonds,
    )

    formula = m.cap_molecule(neighbours_id_1, neighbours_id_4)
    m.assign_bond_orders()
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
