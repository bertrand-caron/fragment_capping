from typing import Optional, TextIO, Sequence, Any

def write_to_debug(debug: Optional[TextIO], *objects: Sequence[Any]) -> None:
     if debug is not None:
        debug.write(' '.join(map(str, objects)) + '\n')
