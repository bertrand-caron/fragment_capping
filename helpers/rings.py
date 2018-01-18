from typing import List, Set
from functools import reduce

from fragment_capping.helpers.types_helpers import ATOM_INDEX, Bond, Atom

SMALL_RING = 7

def bonds_for_ring(ring: List[ATOM_INDEX]) -> List[Bond]:
    return [
        frozenset(atom_indices)
        for atom_indices in zip(ring, ring[1:] + (ring[0],))
    ]

def atoms_in_small_rings_for(molecule: 'Molecule') -> List[Atom]:
    return [
        molecule.atoms[atom_index]
        for atom_index in reduce(
            lambda acc, e: acc | set(e),
            [ring for ring in molecule.rings() if 0 < len(ring) < SMALL_RING],
            set(),
        )
    ]

def bonds_in_small_rings_for(molecule: 'Molecule') -> Set[Bond]:
    return reduce(
        lambda acc, e: acc | set(e),
        [bonds_for_ring(ring) for ring in molecule.rings() if 0 < len(ring) < SMALL_RING],
        set(),
    )

def atom_indices_in_phenyl_rings_for(molecule: 'Molecule') -> Set[ATOM_INDEX]:
    return reduce(
        lambda acc, e: acc | set(e),
        [
            ring
            for ring in molecule.rings()
            if (
                len(ring) == 6
                and
                all(
                    len([bond for bond in molecule.bonds if atom_index in bond]) == 3
                    and
                    molecule.atoms[atom_index].element == 'C'
                    for atom_index in ring
                )
            )
        ],
        set(),
    )
