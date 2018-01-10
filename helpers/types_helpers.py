from typing import Dict, Any, NamedTuple, Optional, Tuple, FrozenSet
from os.path import dirname, abspath

Fragment = str

ATB_Molid = int

ATOM_INDEX = int

Atom = NamedTuple(
    'Atom',
    [
        ('index', ATOM_INDEX),
        ('element', str),
        ('valence', Optional[int]),
        ('capped', bool),
        ('coordinates', Optional[Tuple[float, float, float]],)
    ],
)

Bond = FrozenSet[ATOM_INDEX]

FRAGMENT_CAPPING_DIR = dirname(dirname(abspath(__file__)))

DESC = lambda x: -x

MAX, MIN = (lambda x: DESC(x), lambda x: x)
