from pickle import load
from os.path import join, exists
from typing import Any, List, Optional, Tuple

from fragment_capping.cache import cached
from fragment_capping.helpers.types_helpers import Fragment, ATB_Molid, Atom
from fragment_capping.helpers.molecule import Molecule

from API_client.api import API, HTTPError
from fragment_dihedrals.fragment_dihedral import element_valence_for_atom, on_asc_number_electron_then_asc_valence, NO_VALENCE

from cairosvg import svg2png # pylint: disable=no-name-in-module

FIGSIZE = (30, 15) # inches ?

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

    bonds = [(neighbour_id, atom_id_2) for neighbour_id in neighbours_id_1] + [CENTRAL_BOND] + [(atom_id_3, neighbour_id) for neighbour_id in neighbours_id_4]

    m = Molecule(
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
        name=fragment_str.replace('|', '_'),
    )

    print(m)
    for (i, n, j) in map(lambda cycle: map(int, cycle), cycles):
        i_id, j_id = neighbours_id_1[i], neighbours_id_4[j]
        if n == 0:
            # i and j are actually the same atoms
            del m.atoms[j_id]
            replace_j_by_i = lambda x: i_id if x == j_id else x
            m.bonds = [
                list(map(replace_j_by_i, bond))
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
    print(m)
    return m

def cap_fragment(fragment: Fragment, count: Optional[int] = None, i: Optional[int] = None, fragments: Optional[List[Any]] = None):
    if all([x is not None for x in (count, i, fragments)]):
        print('Running fragment {0}/{1} (count={2}): "{3}"'.format(
            i + 1,
            len(fragments),
            count,
            fragment,
        ))

    molecule = best_capped_molecule_for_dihedral_fragment(fragment)

    try:
        api_response = api.Molecules.structure_search(
            netcharge='*',
            structure_format='pdb',
            structure=molecule.dummy_pdb(),
            return_type='molecules',
        )
    except HTTPError:
        print(molecule.dummy_pdb())
        raise

    molecules = api_response['matches']

    if molecules:
        print([atb_molecule['molid'] for atb_molecule in molecules])
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
            print(e)
            best_molid = None

    else:
        print('Capped fragment not found in ATB.')
        print(molecule.formula())
        print(molecule.dummy_pdb())
        best_molid = None

    print()

    safe_fragment_name = fragment.replace('|', '_')

    with open('pdbs/{fragment}.pdb'.format(fragment=safe_fragment_name), 'w') as fh:
        fh.write(molecule.dummy_pdb())

    return best_molid

def get_matches(protein_fragments: List[Tuple[Fragment, int]]) -> List[Tuple[Fragment, ATB_Molid]]:
    matches = [
        (fragment, cap_fragment(fragment, count=count, i=i, fragments=protein_fragments))
        for (i, (fragment, count)) in
        enumerate(protein_fragments)
    ]

    for (fragment, molid) in matches:
        if molid:
            print('python3 test.py --download {molid} --submit --dihedral-fragment "{dihedral_fragment}"'.format(
                molid=molid,
                dihedral_fragment=fragment,
            ))
    return matches

def parse_args() -> Any:
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--only-id', type=int, help='Rerun a single fragment')
    parser.add_argument('--figsize', nargs=2, type=int, default=FIGSIZE, help='Figure dimensions (in inches)')

    return parser.parse_args()

def generate_collage(protein_fragments, figsize=FIGSIZE) -> bool:
    matches = cached(get_matches, (protein_fragments,), {}, hashed=True)
    counts = dict(protein_fragments)

    for (fragment, molid) in matches:
        if molid:
            print(
                'python3 test.py --download {molid} --submit --dihedral-fragment "{dihedral_fragment}"'.format(
                    molid=molid,
                    dihedral_fragment=fragment,
                ),
            )

    print(matches)
    print('INFO: Assigned {0}/{1} molecules (missing_ids = {2})'.format(
        len([__molid for __molid in matches if __molid[1] is not None]),
        len(matches),
        [i for (i, (fragment, molid)) in enumerate(matches) if molid is None],
    ))

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

    def figure_collage(figsize=FIGSIZE):
        import matplotlib.pyplot as p
        from PIL import Image

        from fragment_capping.collage import best_grid, COLUMNS, LINES, indices_for_subplot

        subplot_dims = best_grid(
            len(matches),
            aspect_ratio=figsize,
        )

        fig, axarr = p.subplots(*subplot_dims, figsize=figsize)

        for (n, (fragment, molid)) in enumerate(matches):
            indices = indices_for_subplot(n, subplot_dims)
            if molid in png_files:
                image = Image.open(png_files[molid])
                axarr[indices].imshow(image)
            axarr[indices].set_title(
                fragment + (' (id={1}, molid={0})'.format(molid, n)),
                fontsize=11,
                fontname='Andale Mono',
            )
            axarr[indices].set_axis_off()

        for n in range(len(matches), subplot_dims[LINES] * subplot_dims[COLUMNS]):
            indices = indices_for_subplot(n, subplot_dims)
            axarr[indices].set_axis_off()

        p.tight_layout()
        #p.show()
        fig.savefig('collage.png', dpi=400)

    figure_collage(figsize=figsize)

    return True

REMOVE_VALENCES = False

EXCLUDE_CYCLIC_FRAGMENTS = True

DUMP_NUMBERED_FRAGMENTS = True

def get_protein_fragments() -> Any:
    with open('cache/protein_fragments.pickle', 'rb') as fh:
        protein_fragments = load(fh)

    if REMOVE_VALENCES:
        from re import sub
        protein_fragments = [(sub('[0-9]', '', fragment), count) for (fragment, count) in protein_fragments]

    if EXCLUDE_CYCLIC_FRAGMENTS:
        print(protein_fragments)
        protein_fragments = [(fragment, count) for (fragment, count) in protein_fragments if fragment.count('|') == 3]

    if DUMP_NUMBERED_FRAGMENTS:
        print()
        def numbered_fragments():
            return {fragment: n for (n, (fragment, count)) in enumerate(protein_fragments)}

        cached(
            numbered_fragments,
            (),
            {},
        )

    return protein_fragments

def main(only_id: Optional[int] = None, figsize: Tuple[int, int] = FIGSIZE):
    protein_fragments = get_protein_fragments()

    if only_id:
        print(cap_fragment(
            protein_fragments[only_id][0],
            i=0,
            count='unknown',
            fragments=(True,),
        ))
    else:
        generate_collage(protein_fragments, figsize=figsize)

if __name__ == '__main__':
    args = parse_args()
    main(
        only_id=args.only_id,
        figsize=tuple(args.figsize),
    )
