from typing import Optional

def sybyl_atom_type(atom_element: Optional[str], atom_valence: int) -> str:
    '''
    Source: http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

    Args:
        ``atom_element``: String (capitalised or not) describing the chemical element of the atom, or None.
        ``atom_valence``: Integer describing the number of bonded atoms.

    Returns:
        Sybyl atom type (String)
    '''
    assert atom_valence is not None, 'Missing valence: {0}'.format(atom_valence)

    if atom_element is None:
        return 'Du' #1.1
    elif atom_element == 'D':
        return 'H' # 1.3
    elif atom_element == 'P':
        return 'P.3' # 1.4
    elif atom_element.upper() in {'RU', 'CO'}:
        return atom_element.title() + '.oh' # 1.5
    elif atom_element in {'C', 'O', 'N', 'S'}: # 1.6, 1.7, 1.8, 1.9,
        if atom_element == 'C': # 1.6
            if atom_valence >= 4:
                valence_suffix = 3
            elif atom_valence == 1:
                valence_suffix = 1
            else:
                valence_suffix = atom_valence - 1
        elif atom_element == 'O': # 1.7
            valence_suffix = atom_valence + 1
        elif atom_element == 'S': # 1.9
            if atom_valence == 3:
                valence_suffix = 'o'
            elif atom_valence == 4:
                valence_suffix = 'o2'
            elif atom_valence == 2:
                valence_suffix = 3
            else:
                valence_suffix = 2
        else: # 1.8
            valence_suffix = atom_valence

        if isinstance(valence_suffix, int):
            assert valence_suffix > 0, (atom, valence_suffix)

        return '{element}.{valence_suffix}'.format(
            element=atom_element.title(),
            valence_suffix=valence_suffix,
        )
    elif atom_element.upper() in {'TI', 'CR'}: # 1.10
        if atom_valence <= 4:
            return atom_element.title() + '.th' # 1.10.1
        else:
            return atom_element.title() + '.oh' # 1.10.2
    else:
        return atom_element.title() # 1.11

if __name__ == '__main__':
    print(sybyl_atom_type('C', 4))
    print(sybyl_atom_type('RU', 4))
