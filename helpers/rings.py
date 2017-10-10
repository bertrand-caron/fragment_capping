from typing import List

from fragment_capping.helpers.types_helpers import ATOM_INDEX, Bond

def bonds_for_ring(ring: List[ATOM_INDEX]) -> List[Bond]:
    return [
        frozenset(atom_indices)
        for atom_indices in zip(ring, ring[1:] + (ring[0],))
    ]
