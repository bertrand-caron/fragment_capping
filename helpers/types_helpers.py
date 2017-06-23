from typing import Dict, Any, NamedTuple, Optional, Tuple
from os.path import dirname, abspath

Fragment = str

ATB_Molid = int

Atom = NamedTuple(
    'Atom',
    [
        ('index', int),
        ('element', str),
        ('valence', Optional[int]),
        ('capped', bool),
        ('coordinates', Optional[Tuple[float, float, float]],)
    ],
)


FRAGMENT_CAPPING_DIR = dirname(dirname(abspath(__file__)))
