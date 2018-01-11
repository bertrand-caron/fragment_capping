from typing import List, Optional, TextIO
from copy import deepcopy
from operator import itemgetter
from itertools import groupby, chain, count
from pulp import LpProblem, LpMinimize, LpInteger, LpVariable, LpBinary, LpStatus, lpSum
from pulp.solvers import GUROBI_CMD

from fragment_capping.helpers.types_helpers import Atom, MIN, MAX
from fragment_capping.helpers.parameters import MAX_ABSOLUTE_CHARGE, MIN_ABSOLUTE_CHARGE, MAX_NONBONDED_ELECTRONS, MAX_BOND_ORDER, MIN_BOND_ORDER, VALENCE_ELECTRONS, ELECTRONS_PER_BOND, MUST_BE_INT, Capping_Strategy, NO_CAP, H_CAP, H2_CAP, H3_CAP, H4_CAP, ELECTRONEGATIVITIES
from fragment_capping.helpers.misc import write_to_debug
from fragment_capping.helpers.graphs import unique_molecules

def get_all_tautomers_naive(
    molecule: 'Molecule',
    total_number_hydrogens: Optional[int] = None,
    net_charge: Optional[int] = None,
    enforce_octet_rule: bool = True,
    allow_radicals: bool = False,
    max_tautomers: Optional[int] = None,
) -> List['Molecule']:
    ELECTRON_MULTIPLIER = (2 if not allow_radicals else 1)

    molecule.remove_all_hydrogens()

    problem = LpProblem("Lewis problem (tautomers) for molecule {0}".format(molecule.name), LpMinimize)

    def maximum_number_hydrogens_for(atom: Atom) -> int:
        if atom.element == 'C' and atom.valence == 3:
            return (3 - len([1 for bond in molecule.bonds if atom.index in bond]))
        else:
            return 3

    capping_atom_ids = set()
    core_atoms = list(molecule.atoms.values())
    for atom in filter(lambda atom: atom.element != 'H', core_atoms):
        # Add up to 3 Hydrogens FIXME: Might use different number of different elements
        # FIXME: Currently cause any unsaturation to get saturated ...
        for _ in range(maximum_number_hydrogens_for(atom)):
            capping_atom_ids.add(
                molecule.add_atom(
                    Atom(index=None, element='H', valence=1, capped=True, coordinates=None),
                    bonded_to=[atom.index]
                ),
            )

    is_capping_bond = lambda bond: len(bond & capping_atom_ids) > 0

    capping_bonds = {
        bond
        for bond in molecule.bonds
        if len(bond & capping_atom_ids) != 0
    }

    charges = {
        atom.index: LpVariable("C_{i}".format(i=atom.index), -MAX_ABSOLUTE_CHARGE, MAX_ABSOLUTE_CHARGE, LpInteger)
        for atom in molecule.atoms.values()
        if atom.index not in capping_atom_ids
    }

    # Extra variable use to bind charges
    absolute_charges = {
        atom_id: LpVariable("Z_{i}".format(i=atom_id), MIN_ABSOLUTE_CHARGE, MAX_ABSOLUTE_CHARGE, LpInteger)
        for atom_id in charges.keys()
    }

    non_bonded_electrons = {
        atom_id: LpVariable("N_{i}".format(i=atom_id), 0, MAX_NONBONDED_ELECTRONS // ELECTRON_MULTIPLIER, LpInteger)
        for atom_id in charges.keys()
    }

    # Maps a bond to an integer
    bond_mapping = {
        bond: i
        for (i, bond) in enumerate(molecule.bonds)
    }

    # Maps an integer to a bond
    bond_reverse_mapping = {v: k for (k, v) in bond_mapping.items()}

    bond_orders = {
        bond: (
            LpVariable("B_{i}".format(i=bond_mapping[bond]), 0, 1, LpBinary)
            if bond in capping_bonds
            else LpVariable("B_{i}".format(i=bond_mapping[bond]), MIN_BOND_ORDER, MAX_BOND_ORDER, LpInteger)
        )
        for bond in molecule.bonds
    }

    OBJECTIVES = [
        MIN(sum(absolute_charges.values())),
        #MIN(sum([charge * ELECTRONEGATIVITIES[molecule.atoms[atom_id].element] for (atom_id, charge) in charges.items()])),
        #MIN(sum([bond_order * ELECTRONEGATIVITIES[molecule.atoms[atom_id].element] for (bond, bond_order) in bond_orders.items() for atom_id in bond])),
    ]

    if net_charge is not None:
        problem += sum(charges.values()) == net_charge, 'Total net charge'

    if total_number_hydrogens is not None:
        problem += sum(bond_orders[bond] for bond in capping_bonds) == total_number_hydrogens, 'Total number of hydrogens'
    else:
        OBJECTIVES.append(MAX(sum(bond_orders[bond] for bond in capping_bonds))) # Add as many hydrogens as possible

    #OBJECTIVES.append(MAX(sum(bond_orders.values())))
    OBJECTIVES.append(MAX(sum(non_bonded_electrons.values())))

    for atom in map(lambda atom_index: molecule.atoms[atom_index], charges.keys()):
        problem += charges[atom.index] == VALENCE_ELECTRONS[atom.element] - sum([bond_orders[bond] for bond in molecule.bonds if atom.index in bond]) - ELECTRON_MULTIPLIER * non_bonded_electrons[atom.index], '{element}_{index}'.format(element=atom.element, index=atom.index)

        # Deal with absolute values
        problem += charges[atom.index] <= absolute_charges[atom.index], 'Absolute charge contraint left {i}'.format(i=atom.index)
        problem += -charges[atom.index] <= absolute_charges[atom.index], 'Absolute charge contraint right {i}'.format(i=atom.index)

        if enforce_octet_rule:
            if atom.element not in {'B', 'BE', 'P', 'S'}:
                problem += (
                    ELECTRONS_PER_BOND * sum([bond_orders[bond] for bond in molecule.bonds if atom.index in bond]) + ELECTRON_MULTIPLIER * non_bonded_electrons[atom.index] == (2 if atom.element in {'H', 'HE'} else 8),
                    'Octet for atom {element}_{index}'.format(element=atom.element, index=atom.index),
                )

    capping_atom_for_bond = lambda atom_index, bond: list({atom.index} & bond)[0]
    bond_for_capping_atom_id = lambda atom_index: [bond for bond in molecule.bonds if atom_index in bond][0]

    unique_solution_equation = sum(2**i * v for (i, v) in enumerate([bond_orders[bond] for bond in capping_bonds]))

    def encode_solution() -> int:
        return sum(2**i * int(v.varValue) for (i, v) in enumerate([bond_orders[bond] for bond in capping_bonds]))

    def new_molecule_for_current_solution(n: Optional[int] = None) -> 'Molecule':
        new_molecule = deepcopy(molecule)
        if n is not None:
            new_molecule.name += '-tautomer_{0}'.format(n)

        new_molecule.formal_charges, new_molecule.bond_orders, new_molecule.non_bonded_electrons = {}, {}, {}

        atom_indices_to_delete = set()
        for v in problem.variables():
            variable_type, variable_substr = v.name.split('_')
            if variable_type == 'C':
                atom_index = int(variable_substr)
                new_molecule.formal_charges[atom_index] = MUST_BE_INT(v.varValue)
            elif variable_type == 'B':
                bond_index = int(variable_substr)
                if MUST_BE_INT(v.varValue) == 0:
                    atom_indices_to_delete |= (capping_atom_ids & set([atom_index for atom_index in new_molecule.atoms.keys() if atom_index in bond_reverse_mapping[bond_index]]))
                new_molecule.bond_orders[bond_reverse_mapping[bond_index]] = MUST_BE_INT(v.varValue)
                pass
            elif variable_type == 'Z':
                pass
            elif variable_type == 'N':
                atom_index = int(variable_substr)
                new_molecule.non_bonded_electrons[atom_index] = ELECTRON_MULTIPLIER * MUST_BE_INT(v.varValue)
            elif variable_type == 'K':
                pass
            else:
                raise Exception('Unknown variable type: {0}'.format(variable_type))

        for atom_index in capping_atom_ids:
            # Manually add electronic properties of hydrogens: charge=0, non_bonded_electrons=0
            new_molecule.non_bonded_electrons[atom_index] = 0
            new_molecule.formal_charges[atom_index] = 0

        REMOVE_UNUSED_HYDROGENS = True
        if REMOVE_UNUSED_HYDROGENS:
            [new_molecule.remove_atom_with_index(atom_index) for atom_index in atom_indices_to_delete]

        if net_charge is not None:
            assert new_molecule.netcharge() == net_charge

        if total_number_hydrogens is not None and REMOVE_UNUSED_HYDROGENS:
            assert len([1 for atom in new_molecule.atoms.values() if atom.element == 'H']) == total_number_hydrogens

        return new_molecule

    def debug_failed_ILP(n: Optional[int] = None) -> None:
        debug_filename = 'debug{0}.lp'.format('' if n is None else '_{n}'.format(n=n))
        problem.writeLP(debug_filename)
        print('Failed LP written to "{0}"'.format(debug_filename))

    # Solve once to find optimal solution with lowest encode_solution()
    try:
        print('Solving')
        problem.sequentialSolve(OBJECTIVES + [MIN(unique_solution_equation)])
    except Exception as e:
        problem.writeLP('debug.lp')
        new_molecule.write_graph('DEBUG', output_size=(2000, 2000))
        print('Failed LP written to "debug.lp"')
        raise

    all_tautomers = [new_molecule_for_current_solution(n=0)]

    # Remove redundant constraints from previous multi-objective optimistion
    for constraint_name in ['1_Sequence_Objective', '2_Sequence_Objective']:
        del problem.constraints[constraint_name]

    # Iterate until no more tautomers are found
    for n in (count(1) if max_tautomers is None else range(1, max_tautomers)):
        problem += unique_solution_equation >= encode_solution() + 1, 'Solution {n}'.format(n=n)
        print('Excluding all solution below {0}'.format(encode_solution()))

        try:
            problem.setObjective(MIN(unique_solution_equation))
            problem.solve()
        except Exception as e:
            debug_failed_ILP()
            raise

        if problem.status == 1:
            pass
        else:
            print(problem.status, LpStatus[problem.status])
            break

        all_tautomers.append(new_molecule_for_current_solution(n=n))

    return unique_molecules(all_tautomers)

def get_all_tautomers(
    molecule: 'Molecule',
    total_number_hydrogens: Optional[int] = None,
    net_charge: Optional[int] = None,
    enforce_octet_rule: bool = True,
    allow_radicals: bool = False,
    max_tautomers: Optional[int] = 100,
    use_gurobi: bool = True,
    debug: Optional[TextIO] = None,
) -> 'Molecule':
    if len([1 for atom in molecule.atoms.values() if atom.element == 'H']) > 0:
        molecule.remove_all_hydrogens(mark_all_uncapped=True)

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
        return [NO_CAP, H_CAP, H2_CAP, H3_CAP] + ([H4_CAP] if False else [])
#        if debug is not None:
#            write_to_debug(debug, '')
#            write_to_debug(debug, atom)
#            write_to_debug(debug, 'capping_strategy, new_atom_for_capping_strategy(), keep_capping_strategy_for_atom()')
#            for capping_strategy in ALL_CAPPING_OPTIONS[molecule.atom_desc(atom)]:
#                write_to_debug(
#                    debug,
#                    capping_strategy,
#                    new_atom_for_capping_strategy(capping_strategy),
#                    keep_capping_strategy_for_atom(capping_strategy, atom)
#                )
#        return [
#            capping_strategy
#            for capping_strategy in ALL_CAPPING_OPTIONS[molecule.atom_desc(atom)]
#            if keep_capping_strategy_for_atom(capping_strategy, atom)
#        ]

    atoms_need_capping = [atom for atom in molecule.sorted_atoms() if not atom.capped]

    write_to_debug(
        debug,
        [
            possible_capping_strategies_for_atom(atom)
            for atom in atoms_need_capping
        ],
    )

    problem = LpProblem("Tautomer enumeration problem for molecule {0}".format(molecule.name), LpMinimize)

    ELECTRON_MULTIPLIER = (2 if not allow_radicals else 1)

    fragment_switches, fragment_scores, fragment_H_scores = {}, {}, {}
    capping_atoms_for = {}
    new_bonds_sets = {}
    for uncapped_atom in atoms_need_capping:
        possible_capping_strategies = possible_capping_strategies_for_atom(uncapped_atom)
        if len(possible_capping_strategies) == 0 or len(possible_capping_strategies) == 1 and possible_capping_strategies[0] == NO_CAP:
            pass
            print('PASS')
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

            # Only choose one capping strategy at a time
            problem += (lpSum(F_i for ((atom_id, _), F_i) in fragment_switches.items() if atom_id == uncapped_atom.index) == 1, 'Single capping strategy for atom {element}_{index}'.format(element=uncapped_atom.element, index=uncapped_atom.index))

    all_capping_atoms = {atom for atoms in capping_atoms_for.values() for atom in atoms}
    all_capping_atom_ids = {atom.index for atom in all_capping_atoms}
    non_capping_atoms = [atom for atom in molecule.atoms.values() if atom.index not in all_capping_atom_ids]
    fragment_switch_for_bond = {
        bond: fragment_switch
        for ((uncapped_atom_id, fragment_id), fragment_switch) in fragment_switches.items()
        for bond in new_bonds_sets[uncapped_atom_id, fragment_id]
    }
    by_switch_name = lambda T: T[1].name
    get_bond = lambda T: T[0]
    bonds_for_fragment_switch = {
        key: [get_bond(T) for T in group]
        for (key, group) in chain(
            groupby(
                sorted(
                    fragment_switch_for_bond.items(),
                    key=by_switch_name,
                ),
                key=by_switch_name,
            ),
            # Switches with no bonds, which are being "ignored" by groupby
            [
                (fragment_switches[uncapped_atom_id, i].name, [])
                for ((uncapped_atom_id, i), bonds) in new_bonds_sets.items()
                if len(bonds) == 0
            ]
        )
    }

    if True:
        molecule.write_graph('debug')

    charges = {
        atom.index: LpVariable("C_{i}".format(i=atom.index), -MAX_ABSOLUTE_CHARGE, MAX_ABSOLUTE_CHARGE, LpInteger)
        for atom in non_capping_atoms
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
            MIN_BOND_ORDER,
            MAX_BOND_ORDER,
            LpInteger,
        )
        if bond not in fragment_switch_for_bond
        else fragment_switch_for_bond[bond]
        for bond in molecule.bonds
    }

    #for ((uncapped_atom_id, i), new_bonds) in new_bonds_sets.items():
    #    for new_bond in new_bonds:
    #        problem += (bond_orders[new_bond] >= fragment_switches[uncapped_atom_id, i], 'Minimum bond order for fragment bond {bond_key}'.format(bond_key=bond_key(new_bond)))
    #        problem += (bond_orders[new_bond] <= MAX_BOND_ORDER * fragment_switches[uncapped_atom_id, i], 'Maximum bond order for fragment bond {bond_key}'.format(bond_key=bond_key(new_bond)))

    OBJECTIVES = [
        MIN(lpSum(absolute_charges.values())),
    ]

    H_size = lpSum([F_i * fragment_H_scores[uncapped_atom_id, i] for ((uncapped_atom_id, i), F_i) in fragment_switches.items()])

    if total_number_hydrogens is not None:
        problem += H_size == total_number_hydrogens, 'Total number of hydrogens={0}'.format(total_number_hydrogens)

    total_size_objective = MIN(lpSum([F_i * fragment_scores[uncapped_atom_id, i] for ((uncapped_atom_id, i), F_i) in fragment_switches.items()]))
    if sum([fragment_scores[uncapped_atom_id, i] for ((uncapped_atom_id, i), F_i) in fragment_switches.items()]) != 0:
        OBJECTIVES.append(total_size_objective)

    OBJECTIVES.extend([
        MIN(lpSum([charge * ELECTRONEGATIVITIES[molecule.atoms[atom_id].element] for (atom_id, charge) in charges.items()])),
        MIN(lpSum([bond_order * ELECTRONEGATIVITIES[molecule.atoms[atom_id].element] for (bond, bond_order) in bond_orders.items() for atom_id in bond])),
    ])

    if net_charge is not None:
        problem += (lpSum(charges.values()) == molecule.net_charge, 'Known net charge={0}'.format(net_charge))

    for atom in non_capping_atoms:
        problem += (
            charges[atom.index]
            ==
            VALENCE_ELECTRONS[atom.element]
            -
            lpSum([bond_orders[bond] for bond in molecule.bonds if atom.index in bond])
            -
            ELECTRON_MULTIPLIER * non_bonded_electrons[atom.index]
            ,
            'Electron balance for atom {element}_{index}'.format(element=atom.element, index=atom.index),
        )

    # Deal with absolute values
    for atom in non_capping_atoms:
        problem += charges[atom.index] <= absolute_charges[atom.index], 'Absolute charge contraint left {i}'.format(i=atom.index)
        problem += -charges[atom.index] <= absolute_charges[atom.index], 'Absolute charge contraint right {i}'.format(i=atom.index)

        if enforce_octet_rule and atom not in all_capping_atoms:
            if atom.element not in {'B', 'BE', 'P', 'S'}:
                problem += (
                    ELECTRONS_PER_BOND * lpSum([bond_orders[bond] for bond in molecule.bonds if atom.index in bond]) + ELECTRON_MULTIPLIER * non_bonded_electrons[atom.index] == (2 if atom.element in {'H', 'HE'} else 8),
                    'Octet for atom {element}_{index}'.format(element=atom.element, index=atom.index),
                )

    def debug_failed_ILP(n: Optional[int] = None) -> None:
        debug_filename = 'debug{0}.lp'.format('' if n is None else '_{n}'.format(n=n))
        problem.writeLP(debug_filename)
        print('Failed LP written to "{0}"'.format(debug_filename))

    unique_solution_equation = sum(2**i * v for (i, v) in enumerate(fragment_switches.values()))

    def encode_solution() -> int:
        return sum(2**i * int(v.varValue) for (i, v) in enumerate(fragment_switches.values()))

    def new_molecule_for_current_solution(n: Optional[int] = None) -> 'Molecule':
        new_molecule = deepcopy(molecule)
        if n is not None:
            new_molecule.name += '_tautomer_{0}'.format(n)

        DELETE_FAILED_CAPS = True

        new_molecule.formal_charges, new_molecule.bond_orders, new_molecule.non_bonded_electrons = {}, {}, {}
        atoms_to_remove = set()
        for v in problem.variables():
            variable_type, variable_substr = v.name.split('_')
            if variable_type == 'C':
                atom_index = int(variable_substr)
                new_molecule.formal_charges[atom_index] = MUST_BE_INT(v.varValue)
            elif variable_type == 'B':
                if False:
                    bond_index = int(variable_substr)
                    new_molecule.bond_orders[bond_reverse_mapping[bond_index]] = MUST_BE_INT(v.varValue)
                else:
                    bond = frozenset(map(int, variable_substr.split(',')))
                    new_molecule.bond_orders[bond] = MUST_BE_INT(v.varValue)
            elif variable_type == 'Z':
                pass
            elif variable_type == 'N':
                atom_index = int(variable_substr)
                new_molecule.non_bonded_electrons[atom_index] = MUST_BE_INT(v.varValue) * ELECTRON_MULTIPLIER
            elif variable_type == 'F':
                uncapped_atom_id, capping_strategy_id = map(int, variable_substr.split(','))
                if MUST_BE_INT(v.varValue) == 0 and DELETE_FAILED_CAPS:
                    atoms_to_remove.add((uncapped_atom_id, capping_strategy_id))
                for bond in bonds_for_fragment_switch[v.name]:
                    new_molecule.bond_orders[bond] = MUST_BE_INT(v.varValue)
            elif variable_type == 'S':
                capping_atom_id = int(variable_substr)
            else:
                raise Exception('Unknown variable type: {0}'.format(variable_type))

        for atom_index in all_capping_atom_ids:
            # Manually add electronic properties of hydrogens: charge=0, non_bonded_electrons=0
            new_molecule.formal_charges[atom_index] = 0
            new_molecule.non_bonded_electrons[atom_index] = 0

        if DELETE_FAILED_CAPS:
            for (uncapped_atom_id, capping_strategy_id) in atoms_to_remove:
                new_molecule.remove_atoms(atom for atom in capping_atoms_for[uncapped_atom_id, capping_strategy_id])

        if not allow_radicals and False:
            assert all([nonbonded_electrons % 2 == 0 for nonbonded_electrons in new_molecule.non_bonded_electrons.values()]), {
                new_molecule.atoms[atom_index]: electrons
                for (atom_index, electrons) in new_molecule.non_bonded_electrons.items()
                if electrons % 2 == 1
            }

        new_molecule.assert_molecule_coherence()
        return new_molecule

    # Solve once to find optimal solution with lowest encode_solution()
    try:
        problem.sequentialSolve(OBJECTIVES)
        assert problem.status == 1, (molecule.name, LpStatus[problem.status])
    except Exception as e:
        debug_failed_ILP(0)
        raise

    debug_failed_ILP(0)
    all_tautomers = [new_molecule_for_current_solution(n=0)]

    # Remove redundant constraints from previous multi-objective optimistion
    for constraint_name in ['1_Sequence_Objective', '2_Sequence_Objective', '3_Sequence_Objective']:
        del problem.constraints[constraint_name]

    # Iterate until no more tautomers are found
    for n in (count(1) if max_tautomers is None else range(1, max_tautomers)):
        problem += unique_solution_equation >= encode_solution() + 1, 'Solution {n}'.format(n=n)
        print('Excluding all solution below {0}'.format(encode_solution()))

        try:
            problem.setObjective(MIN(unique_solution_equation))
            problem.solve(
                solver=GUROBI_CMD() if use_gurobi else None,
            )
        except Exception as e:
            debug_failed_ILP(n)
            raise

        if problem.status == 1:
            pass
        else:
            print(problem.status, LpStatus[problem.status])
            debug_failed_ILP(n)
            break

        all_tautomers.append(new_molecule_for_current_solution(n=n))

    return unique_molecules(all_tautomers)
