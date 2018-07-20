from typing import Optional, TextIO, Dict, Any, List
from sys import stdout

from fragment_capping.config import ILP_SOLVER_TIMEOUT
from fragment_capping.helpers.types_helpers import Atom, MIN, MAX
from fragment_capping.helpers.parameters import MAX_ABSOLUTE_CHARGE, MIN_ABSOLUTE_CHARGE, MAX_NONBONDED_ELECTRONS, MAX_BOND_ORDER, MIN_BOND_ORDER, VALENCE_ELECTRONS, ELECTRONS_PER_BOND, MUST_BE_INT, ALL_CAPPING_OPTIONS, ELECTRONEGATIVITIES, Capping_Strategy, NO_CAP, new_atom_for_capping_strategy, max_valence_for, min_valence_for
from fragment_capping.helpers.misc import write_to_debug
from fragment_capping.helpers.exceptions import Too_Many_Permutations, Not_Capped_Error
from fragment_capping.helpers.iterables import MAXIMUM_PERMUTATION_NUMBER

def get_best_capped_molecule_with_ILP(
    molecule: 'Molecule',
    net_charge: Optional[int] = None,
    number_hydrogens: Optional[int] = None,
    enforce_octet_rule: bool = True,
    allow_radicals: bool = False,
    debug: Optional[TextIO] = None,
) -> 'Molecule':
    '''
    Use an ILP to cap (complete the valence) of an uncapped molecule using a library a capping fragments.

    Args:
        ``molecule``: Molecule to be capped. Some of its atoms should have the ``capped`` attribute set to False.
        ``net_charge``: (Optional) Constraint the total net charge for the capped molecule.
        ``number_hydrogens``: (Optional) Constraint the total number of hydrogens for the capped molecule (including hydrogens already present in the molecule).
        ``enforce_octet_rule``: (Optional) Constraint organic elements (H, C, N, O) to satisfy the octet rule.
        ``allow_radicals``: (Optional) Allow unpaired non-bonded electrons.
        ``debug``: (Optional) Print very detailed debugging information

    Returns:
        The modified, capped ``molecule``.
    '''
    neighbour_counts = molecule.neighbours_for_atoms()

    def keep_capping_strategy_for_atom(capping_strategy: Capping_Strategy, atom: Atom) -> bool:
        if atom.valence is not None:
            if False:
                if neighbour_counts[atom.index] + new_atom_for_capping_strategy(capping_strategy) == atom.valence:
                    print(atom, capping_strategy)
            return neighbour_counts[atom.index] + new_atom_for_capping_strategy(capping_strategy) == atom.valence
        else:
            return min_valence_for(atom) <= neighbour_counts[atom.index] + new_atom_for_capping_strategy(capping_strategy) <= max_valence_for(atom)

    def possible_capping_strategies_for_atom(atom: Atom) -> List[Capping_Strategy]:
        if debug is not None:
            write_to_debug(debug, '')
            write_to_debug(debug, atom)
            write_to_debug(debug, 'capping_strategy, new_atom_for_capping_strategy(), keep_capping_strategy_for_atom()')
            for capping_strategy in ALL_CAPPING_OPTIONS[molecule.atom_desc(atom)]:
                write_to_debug(
                    debug,
                    capping_strategy,
                    new_atom_for_capping_strategy(capping_strategy),
                    keep_capping_strategy_for_atom(capping_strategy, atom)
                )
        return [
            capping_strategy
            for capping_strategy in ALL_CAPPING_OPTIONS[molecule.atom_desc(atom)]
            if keep_capping_strategy_for_atom(capping_strategy, atom)
        ]

    atoms_need_capping = [atom for atom in molecule.sorted_atoms() if not atom.capped]

    assert len(atoms_need_capping) > 0, 'Error: There are no uncapped atoms in the molecule.'

    if False:
        capping_schemes = list(
            product(
                *[
                    possible_capping_strategies_for_atom(atom)
                    for atom in atoms_need_capping
                ]
            ),
        )

        assert len(capping_schemes) > 0, [
            (
                atom,
                possible_capping_strategies_for_atom(atom),
            )
            for atom in atoms_need_capping
            if len(possible_capping_strategies_for_atom(atom)) == 0
        ]

    write_to_debug(
        debug,
        [
            possible_capping_strategies_for_atom(atom)
            for atom in atoms_need_capping
        ],
    )

    from pulp import LpProblem, LpMinimize, LpInteger, LpVariable, LpBinary, LpStatus, lpSum

    problem = LpProblem("Capping problem for molecule {0}".format(molecule.name), LpMinimize)

    ELECTRON_MULTIPLIER = (2 if not allow_radicals else 1)

    hydrogens_before_capping = len([1 for atom in molecule.atoms.values() if atom.element == 'H'])

    counter_charges = {}
    fragment_switches, fragment_scores, fragment_H_scores = {}, {}, {}
    capping_atoms_for = {}
    new_bonds_sets = {}
    for uncapped_atom in atoms_need_capping:
        possible_capping_strategies = possible_capping_strategies_for_atom(uncapped_atom)
        if len(possible_capping_strategies) == 0 or len(possible_capping_strategies) == 1 and possible_capping_strategies[0] == NO_CAP:
            pass
            stdout.write('\nWarning: No capping strategy for atom: {0}. The ILP will be infeasible if a suitable cap is not available.'.format(uncapped_atom))
        else:
            for (i, capping_strategy) in enumerate(sorted(possible_capping_strategies), start=1):
                write_to_debug(debug, uncapped_atom, capping_strategy, i)
                # Add switch variable
                fragment_switches[uncapped_atom.index, i] = LpVariable(
                    'F_{i},{j}'.format(i=uncapped_atom.index, j=i),
                    0,
                    1,
                    LpBinary,
                )

                new_atoms, new_bonds = molecule.extend_molecule_with(uncapped_atom, capping_strategy)
                write_to_debug(debug, i, [atom for atom in new_atoms])
                capping_atoms_for[uncapped_atom.index, i] = new_atoms
                new_bonds_sets[uncapped_atom.index, i] = [bond for bond in new_bonds if uncapped_atom.index in bond]
                fragment_scores[uncapped_atom.index, i] = len(capping_atoms_for[uncapped_atom.index, i])
                fragment_H_scores[uncapped_atom.index, i] = len([atom for atom in capping_atoms_for[uncapped_atom.index, i] if atom.element == 'H'])

                for capping_atom in new_atoms:
                    # Add counter-charge variable S_i for every atom of the capping strategy
                    counter_charges[capping_atom.index] = LpVariable(
                        "S_{i}".format(i=capping_atom.index),
                        -MAX_ABSOLUTE_CHARGE,
                        MAX_ABSOLUTE_CHARGE,
                        LpInteger,
                    )
                    problem += (counter_charges[capping_atom.index] <=  (1 - fragment_switches[uncapped_atom.index, i]) * MAX_ABSOLUTE_CHARGE, 'Maximum counter charge for capping atom {element}_{index}'.format(element=capping_atom.element, index=capping_atom.index))
                    problem += (counter_charges[capping_atom.index] >= -(1 - fragment_switches[uncapped_atom.index, i]) * MAX_ABSOLUTE_CHARGE, 'Minimum counter charge for capping atom {element}_{index}'.format(element=capping_atom.element, index=capping_atom.index))

            # Only choose one capping strategy at a time
            problem += (lpSum(F_i for ((atom_id, _), F_i) in fragment_switches.items() if atom_id == uncapped_atom.index) == 1, 'Single capping strategy for atom {element}_{index}'.format(element=uncapped_atom.element, index=uncapped_atom.index))

    all_capping_atoms = {atom for atoms in capping_atoms_for.values() for atom in atoms}

    if True:
        molecule.write_graph('debug')

    charges = {
        atom.index: LpVariable("C_{i}".format(i=atom.index), -MAX_ABSOLUTE_CHARGE, MAX_ABSOLUTE_CHARGE, LpInteger)
        for atom in molecule.atoms.values()
    }
    original_charges = list(charges.values())

    # Extra variable use to bind charges
    absolute_charges = {
        atom_id: LpVariable("Z_{i}".format(i=atom_id), MIN_ABSOLUTE_CHARGE, MAX_ABSOLUTE_CHARGE, LpInteger)
        for atom_id in charges.keys()
    }

    non_bonded_electrons = {
        atom_id: LpVariable("N_{i}".format(i=atom_id), 0, MAX_NONBONDED_ELECTRONS // ELECTRON_MULTIPLIER, LpInteger)
        for (atom_id, atom) in molecule.atoms.items()
    }

    # Maps a bond to an integer
    bond_mapping = {
        bond: i
        for (i, bond) in enumerate(molecule.bonds)
    }

    # Maps an integer to a bond
    bond_reverse_mapping = {v: k for (k, v) in bond_mapping.items()}
    bond_key = lambda bond: ','.join(map(str, sorted(bond)))

    bond_orders = {
        bond: LpVariable(
            "B_{bond_key}".format(bond_key=bond_key(bond)),
            0 if any(bond in new_bonds for new_bonds in new_bonds_sets.values()) else MIN_BOND_ORDER,
            MAX_BOND_ORDER,
            LpInteger,
        )
        for bond in molecule.bonds
    }

    for ((uncapped_atom_id, i), new_bonds) in new_bonds_sets.items():
        for new_bond in new_bonds:
            problem += (bond_orders[new_bond] >= fragment_switches[uncapped_atom_id, i], 'Minimum bond order for fragment bond {bond_key}'.format(bond_key=bond_key(new_bond)))
            problem += (bond_orders[new_bond] <= MAX_BOND_ORDER * fragment_switches[uncapped_atom_id, i], 'Maximum bond order for fragment bond {bond_key}'.format(bond_key=bond_key(new_bond)))

    OBJECTIVES = [
        MIN(lpSum(absolute_charges.values())),
    ]

    H_size_objective = MAX(lpSum([F_i * fragment_H_scores[uncapped_atom_id, i] for ((uncapped_atom_id, i), F_i) in fragment_switches.items()]))
    has_non_null_H_size_objective = (sum([fragment_H_scores[uncapped_atom_id, i] for ((uncapped_atom_id, i), F_i) in fragment_switches.items()]) != 0)
    if has_non_null_H_size_objective:
        OBJECTIVES.append(H_size_objective)

    total_size_objective = MIN(lpSum([F_i * fragment_scores[uncapped_atom_id, i] for ((uncapped_atom_id, i), F_i) in fragment_switches.items()]))
    if sum([fragment_scores[uncapped_atom_id, i] for ((uncapped_atom_id, i), F_i) in fragment_switches.items()]) != 0:
        OBJECTIVES.append(total_size_objective)

    OBJECTIVES.extend([
        MIN(lpSum([charge * ELECTRONEGATIVITIES[molecule.atoms[atom_id].element] for (atom_id, charge) in charges.items()])),
        MIN(lpSum([bond_order * ELECTRONEGATIVITIES[molecule.atoms[atom_id].element] for (bond, bond_order) in bond_orders.items() for atom_id in bond])),
    ])

    if net_charge is not None:
        problem += (lpSum(charges.values()) == net_charge, 'Known net charge')

    if number_hydrogens is not None:
        problem += (
            lpSum(
                [
                    F_i * fragment_H_scores[uncapped_atom_id, i]
                    for ((uncapped_atom_id, i), F_i) in fragment_switches.items()
                ],
            ) + hydrogens_before_capping == number_hydrogens,
            'Total number of hydrogens: {0}'.format(number_hydrogens)
        )

    for atom in molecule.atoms.values():
        problem += (
            charges[atom.index]
            ==
            VALENCE_ELECTRONS[atom.element]
            -
            lpSum([bond_orders[bond] for bond in molecule.bonds if atom.index in bond])
            -
            ELECTRON_MULTIPLIER * non_bonded_electrons[atom.index]
            -
            (counter_charges[atom.index] if atom.index in counter_charges else 0)
            ,
            'Electron balance for atom {element}_{index}'.format(element=atom.element, index=atom.index),
        )

    # Deal with absolute values
    for atom in molecule.atoms.values():
        problem += charges[atom.index] <= absolute_charges[atom.index], 'Absolute charge contraint 1 {i}'.format(i=atom.index)
        problem += -charges[atom.index] <= absolute_charges[atom.index], 'Absolute charge contraint 2 {i}'.format(i=atom.index)

        if enforce_octet_rule and atom not in all_capping_atoms:
            if atom.element not in {'B', 'BE', 'P', 'S'}:
                problem += (
                    ELECTRONS_PER_BOND * lpSum([bond_orders[bond] for bond in molecule.bonds if atom.index in bond]) + ELECTRON_MULTIPLIER * non_bonded_electrons[atom.index] == (2 if atom.element in {'H', 'HE'} else 8),
                    'Octet for atom {element}_{index}'.format(element=atom.element, index=atom.index),
                )

    try:
        problem.sequentialSolve(OBJECTIVES, timeout=ILP_SOLVER_TIMEOUT)
        assert problem.status == 1, (molecule.name, LpStatus[problem.status])
        #assert False
    except Exception as e:
        problem.writeLP('debug.lp')
        molecule.write_graph('DEBUG', output_size=(1000, 1000))
        print('Failed LP written to "debug.lp"')
        raise

    DELETE_FAILED_CAPS = True

    molecule.formal_charges, molecule.bond_orders, molecule.non_bonded_electrons = {}, {}, {}
    atoms_to_remove = set()
    for v in problem.variables():
        variable_type, variable_substr = v.name.split('_')
        if variable_type == 'C':
            atom_index = int(variable_substr)
            molecule.formal_charges[atom_index] = MUST_BE_INT(v.varValue)
        elif variable_type == 'B':
            if False:
                bond_index = int(variable_substr)
                molecule.bond_orders[bond_reverse_mapping[bond_index]] = MUST_BE_INT(v.varValue)
            else:
                bond = frozenset(map(int, variable_substr.split(',')))
                molecule.bond_orders[bond] = MUST_BE_INT(v.varValue)
        elif variable_type == 'Z':
            pass
        elif variable_type == 'N':
            atom_index = int(variable_substr)
            molecule.non_bonded_electrons[atom_index] = MUST_BE_INT(v.varValue) * ELECTRON_MULTIPLIER
        elif variable_type == 'F':
            uncapped_atom_id, capping_strategy_id = map(int, variable_substr.split(','))
            if MUST_BE_INT(v.varValue) == 0 and DELETE_FAILED_CAPS:
                atoms_to_remove.add((uncapped_atom_id, capping_strategy_id))
        elif variable_type == 'S':
            capping_atom_id = int(variable_substr)
        else:
            raise Exception('Unknown variable type: {0}'.format(variable_type))

    if DELETE_FAILED_CAPS:
        for (uncapped_atom_id, capping_strategy_id) in atoms_to_remove:
            molecule.remove_atoms(atom for atom in capping_atoms_for[uncapped_atom_id, capping_strategy_id])

    if not allow_radicals and False:
        assert all([nonbonded_electrons % 2 == 0 for nonbonded_electrons in molecule.non_bonded_electrons.values()]), {
            molecule.atoms[atom_index]: electrons
            for (atom_index, electrons) in molecule.non_bonded_electrons.items()
            if electrons % 2 == 1
        }

    molecule.update_valences()
    molecule.assign_aromatic_bonds()
    molecule.assert_molecule_coherence()
    return molecule

from itertools import product
from functools import reduce

def get_best_capped_molecule(
    molecule: 'Molecule',
    draw_all_possible_graphs: bool = False,
    debug: Optional[TextIO] = None,
    use_ILP: bool = True,
    **kwargs: Dict[str, Any],
) -> 'Molecule':
    '''
    Use a brute-force combinatorial enumeration to cap (complete the valence) an uncapped molecule using a library a capping fragments.
    Only meant to be used on very small molecules to compare against the ILP version (``get_best_capped_molecule_with_ILP``).

    Args:
        ``molecule``: Molecule to be capped. Some of its atoms should have the ``capped`` attribute set to False.
        ``draw_all_possible_graphs``: (Optional) Draw all possible graphs for debug purposes.
        ``debug``: (Optional) Print very detailed debugging information
        ``use_ILP``: (Optional) Use an ILP for solving the electron assignment problem.
            If set to ``False``, will use a combinatorial enumeration (very slow and error-prone!).
        ``**kwargs``: (Optional) Additional keyword arguments, which are ignored.

    Returns:
        The modified, capped ``molecule``.
    '''
    neighbour_counts = molecule.neighbours_for_atoms()

    def keep_capping_strategy_for_atom(capping_strategy: Capping_Strategy, atom: Atom):
        if molecule.use_neighbour_valences():
            return neighbour_counts[atom.index] + new_atom_for_capping_strategy(capping_strategy) == atom.valence
        else:
            return min_valence_for(atom) <= neighbour_counts[atom.index] + new_atom_for_capping_strategy(capping_strategy) <= max_valence_for(atom)

    def possible_capping_strategies_for_atom(atom: Atom) -> List[Capping_Strategy]:
        if debug is not None:
            write_to_debug(debug, atom)
            for capping_strategy in ALL_CAPPING_OPTIONS[molecule.atom_desc(atom)]:
                write_to_debug(
                    debug,
                    capping_strategy,
                    new_atom_for_capping_strategy(capping_strategy),
                    keep_capping_strategy_for_atom(capping_strategy, atom)
                )
        return [
            capping_strategy
            for capping_strategy in ALL_CAPPING_OPTIONS[molecule.atom_desc(atom)]
            if keep_capping_strategy_for_atom(capping_strategy, atom)
        ]

    atoms_need_capping = [atom for atom in molecule.sorted_atoms() if not atom.capped]
    capping_schemes = list(
        product(
            *[
                possible_capping_strategies_for_atom(atom)
                for atom in atoms_need_capping
            ]
        ),
    )

    assert len(capping_schemes) > 0, [
        (
            atom,
            possible_capping_strategies_for_atom(atom),
        )
        for atom in atoms_need_capping
        if len(possible_capping_strategies_for_atom(atom)) == 0
    ]

    write_to_debug(
        debug,
        [
            possible_capping_strategies_for_atom(atom)
            for atom in atoms_need_capping
        ],
    )

    if len(capping_schemes) >= MAXIMUM_PERMUTATION_NUMBER:
        raise Too_Many_Permutations(len(capping_schemes))

    write_to_debug(debug, 'atoms_need_capping: {0}'.format(atoms_need_capping))
    write_to_debug(debug, 'capping_schemes: {0}'.format(capping_schemes))
    write_to_debug(debug, 'capping_options: {0}'.format([
        len(ALL_CAPPING_OPTIONS[molecule.atom_desc(atom)])
        for atom in atoms_need_capping
    ]))

    write_to_debug(debug, 'atoms_need_capping: {0}'.format(atoms_need_capping))
    write_to_debug(debug, 'INFO: Will try all {0} possible capped molecules'.format(len(capping_schemes)))

    possible_capped_molecules = sorted(
        filter(
            lambda mol: mol is not None,
            [
                molecule.capped_molecule_with(capping_strategies, atoms_need_capping, debug=debug, debug_line='molecule {0}/{1}'.format(i, len(capping_schemes)), use_ILP=use_ILP)
                for (i, capping_strategies) in enumerate(capping_schemes, start=1)
            ],
        ),
        key=lambda mol: (mol.net_abs_charges(), mol.n_atoms(), mol.double_bonds_fitness()),
    )

    write_to_debug(debug, 'Possible capped molecules: {0} ({1}/{2})'.format(
        [(mol.formula(charge=True), mol.net_abs_charges(), mol.double_bonds_fitness()) for mol in possible_capped_molecules],
        len(possible_capped_molecules),
        len(capping_schemes),
    ))

    if draw_all_possible_graphs:
        try:
            for (i, molecule) in enumerate(possible_capped_molecules):
                molecule.write_graph(i)
        except Exception as e:
            print(
                'ERROR: Could not plot graphs (error was: {0})'.format(
                    str(e),
                ),
            )
            raise

    if len(possible_capped_molecules) == 0:
        raise Not_Capped_Error(molecule)

    best_molecule = possible_capped_molecules[0]
    return best_molecule
