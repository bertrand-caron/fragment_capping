from typing import Dict, Any, NamedTuple, Optional
from os.path import dirname, abspath

Fragment = str

ATB_Molid = int

Atom = NamedTuple(
    'Atom',
    [('index', int), ('element', str), ('valence', Optional[int]), ('capped', bool)],
)


FRAGMENT_CAPPING_DIR = dirname(dirname(abspath(__file__)))
