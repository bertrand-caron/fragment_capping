from typing import Optional
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from os import remove

BABEL_MINIMISE_EXECUTABLE = '/usr/local/bin/obminimize'

def energy_minimised_pdb(pdb_filepath: Optional[str] = None, pdb_str: Optional[str] = None, n_steps: int = 5000, convergence: float = 1E-1, debug: bool = False) -> str:
    assert pdb_filepath is not None or pdb_str is not None, [pdb_filepath, pdb_str]

    if pdb_filepath is None:
        with NamedTemporaryFile(mode='w+t', suffix='.pdb', delete=False) as pdb_file:
            pdb_file.write(pdb_str)
            pdb_filepath = pdb_file.name
        should_remove_filepath = True
    else:
        assert pdb_filepath.startswith('/'), pdb_filepath
        should_remove_filepath = False

    subprocess_call = '{BABEL_MINIMISE_EXECUTABLE} -n {n_steps} -c {convergence} {input_pdb_filepath}'.format(
            BABEL_MINIMISE_EXECUTABLE=BABEL_MINIMISE_EXECUTABLE,
            input_pdb_filepath=pdb_filepath,
            n_steps=n_steps,
            convergence=convergence,
    )

    if debug:
        print(subprocess_call)

    try:
        stdout, stderr = map(
            lambda stream: stream.decode(),
            Popen(
                subprocess_call.split(),
                stdout=PIPE,
                stderr=PIPE,
            ).communicate(),
        )
    except FileNotFoundError:
        raise Exception('Could not find "{0}" in $PATH'.format(BABEL_MINIMISE_EXECUTABLE))
    finally:
        if should_remove_filepath:
            remove(pdb_filepath)

    if not stdout:
        raise Exception('Shell command "{0}" failed with error message: "{1}"'.format(subprocess_call, stderr[:-1]))

    return stdout

if __name__ == '__main__':
    print(energy_minimised_pdb(pdb_filepath='pdbs/O,O,C_P_O_C.pdb'))

    with open('pdbs/O,O,C_P_O_C.pdb') as fh:
        print(energy_minimised_pdb(pdb_str=fh.read()))
