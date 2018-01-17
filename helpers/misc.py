from typing import Optional, TextIO, Sequence, Any

from fragment_capping.helpers.types_helpers import Atom

def write_to_debug(debug: Optional[TextIO], *objects: Sequence[Any]) -> None:
     if debug is not None:
        debug.write(' '.join(map(str, objects)) + '\n')

def atom_short_desc(atom: Atom) -> str:
    return '{element}_{index}'.format(element=atom.element, index=atom.index)
