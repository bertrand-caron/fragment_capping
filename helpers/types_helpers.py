from typing import Dict, Any, NamedTuple, Optional

Fragment = str

ATB_Molid = int

Atom = NamedTuple(
    'Atom',
    [('index', int), ('element', str), ('valence', Optional[int]), ('capped', bool)],
)
