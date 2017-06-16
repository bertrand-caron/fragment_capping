from pickle import load
from itertools import groupby, product
from math import sqrt, ceil
from os.path import join, exists, basename, dirname, abspath
from io import StringIO
from functools import reduce
from operator import itemgetter
from typing import Any, List, Optional, Tuple
from re import sub, search
from urllib.request import urlopen
from os.path import dirname, abspath, join

from fragment_capping.cache import cached
from fragment_capping.helpers.types_helpers import Fragment, ATB_Molid, Atom, FRAGMENT_CAPPING_DIR
from fragment_capping.helpers.molecule import Molecule, Too_Many_Permutations
from fragment_capping.helpers.capping import best_capped_molecule_for_dihedral_fragment
from fragment_capping.helpers.babel import energy_minimised_pdb

from API_client.api import API, HTTPError
from dihedral_fragments.dihedral_fragment import element_valence_for_atom, on_asc_number_electron_then_asc_valence, NO_VALENCE

from cairosvg import svg2png # pylint: disable=no-name-in-module

api = API(
    host='http://scmb-atb.biosci.uq.edu.au/atb-uqbcaron', #'https://atb.uq.edu.au',
    debug=False,
    api_format='pickle',
)

def molid_after_capping_fragment(fragment: Fragment, count: Optional[int] = None, i: Optional[int] = None, fragments: Optional[List[Any]] = None) -> ATB_Molid:
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
        print('Formula', molecule.formula())
        print('Smiles', molecule.smiles())
        print('Dummy PDB')
        print(molecule.dummy_pdb())
        print('Energy Minimised Dummy PDB')
        print(energy_minimised_pdb(pdb_str=molecule.dummy_pdb()))
        molecule.write_graph('BEST')
        best_molid = None

    print()

    safe_fragment_name = fragment.replace('|', '_')

    with open(join(FRAGMENT_CAPPING_DIR, 'pdbs/{fragment}.pdb'.format(fragment=safe_fragment_name)), 'w') as fh:
        fh.write(molecule.dummy_pdb())

    return best_molid

def get_matches(protein_fragments: List[Tuple[Fragment, int]]) -> List[Tuple[Fragment, ATB_Molid]]:
    matches = [
        (fragment, molid_after_capping_fragment(fragment, count=count, i=i, fragments=protein_fragments))
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

def truncated_molecule(molecule: Molecule):
    return dict(
        n_atoms=molecule.n_atoms,
        num_dihedral_fragments=len(molecule.dihedral_fragments),
        molid=molecule.molid,
        formula=molecule.formula,
    )

FIGSIZE = (15, 15 * 1.4) # inches ?

def parse_args() -> Any:
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--only-id', type=int, help='Rerun a single fragment')
    parser.add_argument('--figsize', nargs=2, type=int, default=FIGSIZE, help='Figure dimensions (in inches)')

    return parser.parse_args()

def png_file_for(molid: int, force_regen: bool = False, remove_background: bool = True, pixel_height: int = 600, dpi: int = 400):
    PNG_DIR = 'pngs'

    png_file = join(
        dirname(abspath(__file__)),
        PNG_DIR,
        '{molid}.png'.format(molid=molid),
    )

    if (not exists(png_file)) or force_regen:
        try:
            urlopen('https://atb.uq.edu.au/outputs_babel_img.py?molid={molid}'.format(molid=molid))

            vanilla_svg_bytes = urlopen('https://atb.uq.edu.au/cache/img2D/{molid}_thumb.svg'.format(molid=molid)).read()

            modified_svg_bytes = (
                sub(
                    b'<rect.*?/>',
                    b'',
                    vanilla_svg_bytes,
                )
                if remove_background
                else vanilla_svg_bytes
            )


            if remove_background:
                assert b'rect' not in modified_svg_bytes, modified_svg_bytes.decode()
        except:
            raise

        if False:
            match = search(
                b'<svg.* width="([^"])" height="[^"]")',
                modified_svg_bytes,
            )

            svg_width, svg_height = match.group(1), match.group()

        svg2png(
            bytestring=modified_svg_bytes,
            write_to=png_file,
            dpi=dpi,
            height=pixel_height,
        )
    return png_file

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
                fragment + ('\n(id={1}, molid={0})'.format(molid, n)),
                fontsize=11,
                fontname='Andale Mono',
            )
            axarr[indices].set_axis_off()

        for n in range(len(matches), subplot_dims[LINES] * subplot_dims[COLUMNS]):
            indices = indices_for_subplot(n, subplot_dims)
            axarr[indices].set_axis_off()

        p.tight_layout()
        #p.show()
        fig.savefig('protein_fragment_molecules.png', dpi=400)

    figure_collage(figsize=figsize)

    return True

REMOVE_VALENCES = True

EXCLUDE_CYCLIC_FRAGMENTS = True

DUMP_NUMBERED_FRAGMENTS = True

def get_protein_fragments() -> Any:
    with open('data/protein_fragments_with_count.pickle', 'rb') as fh:
        protein_fragments = load(fh)

    if REMOVE_VALENCES:
        from re import sub
        protein_fragments = [(sub('[0-9]', '', fragment).replace('*', ''), count) for (fragment, count) in protein_fragments]

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
        print(molid_after_capping_fragment(
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
