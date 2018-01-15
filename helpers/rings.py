from typing import List, Set
from functools import reduce

from fragment_capping.helpers.types_helpers import ATOM_INDEX, Bond

SMALL_RING = 7

def bonds_for_ring(ring: List[ATOM_INDEX]) -> List[Bond]:
    return [
        frozenset(atom_indices)
        for atom_indices in zip(ring, ring[1:] + (ring[0],))
    ]

def bonds_in_small_rings_for(molecule: 'Molecule') -> Set[Bond]:
    return reduce(
        lambda acc, e: acc | set(e),
        [bonds_for_ring(ring) for ring in molecule.rings() if 0 < len(ring) < SMALL_RING],
        set(),
    )
