from typing import List
from pprint import pprint

def compare_capped_molecules(molecule: 'Molecule', _molecule: 'Molecule', base_atoms_ids: List[int]) -> None:
    neighbour_counts, _neighbour_counts = map(
        lambda __molecule: __molecule.neighbours_for_atoms(),
        (molecule, _molecule),
    )

    pprint(
        list(
            filter(
                lambda T: len(set(T[1])) == 2,
                [
                    (atom_id, tuple(neighbours[atom_id] for neighbours in (neighbour_counts, _neighbour_counts)))
                    for atom_id in base_atoms_ids
                ]
            ),
        ),
    )
