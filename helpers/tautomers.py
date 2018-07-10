from sys import stderr
from typing import List, Optional, TextIO, Tuple
from copy import deepcopy
from operator import itemgetter
from itertools import groupby, chain, count
from pulp import LpProblem, LpMinimize, LpInteger, LpVariable, LpBinary, LpStatus, lpSum
from pulp.solvers import GUROBI_CMD

from fragment_capping.config import ILP_SOLVER_TIMEOUT
from fragment_capping.helpers.types_helpers import Atom, MIN, MAX
from fragment_capping.helpers.parameters import MAX_ABSOLUTE_CHARGE, MIN_ABSOLUTE_CHARGE, MAX_NONBONDED_ELECTRONS, \
    MAX_BOND_ORDER, MIN_BOND_ORDER, VALENCE_ELECTRONS, ELECTRONS_PER_BOND, MUST_BE_INT, Capping_Strategy, NO_CAP, H_CAP, \
    ELECTRONEGATIVITIES, new_atom_for_capping_strategy, min_valence_for, max_valence_for, merge_caps
from fragment_capping.helpers.misc import write_to_debug, atom_short_desc
from fragment_capping.helpers.graphs import unique_molecules
from fragment_capping.helpers.rings import bonds_in_small_rings_for, SMALL_RING, atoms_in_small_rings_for

def get_all_tautomers_naive(
    molecule: 'Molecule',
    total_number_hydrogens: Optional[int] = None,
    net_charge: Optional[int] = None,
    enforce_octet_rule: bool = True,
    allow_radicals: bool = False,
    max_tautomers: Optional[int] = None,
    debug: Optional[TextIO] = stderr,
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

    non_capping_bonds = molecule.bonds - capping_bonds

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
        problem += charges[atom.index] == VALENCE_ELECTRONS[atom.element] - sum([bond_orders[bond] for bond in molecule.bonds if atom.index in bond]) - ELECTRON_MULTIPLIER * non_bonded_electrons[atom.index], atom_short_desc(atom)

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

    def encode_solution() -> int:
        # bitshift is faster than multiplication by 2**i
        return sum(int(v.varValue) << i for i, v in enumerate([bond_orders[bond] for bond in capping_bonds]))

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

    # Solve once to find optimal solution
    try:
        write_to_debug(debug, 'Solving')
        problem.sequentialSolve(OBJECTIVES, timeout=ILP_SOLVER_TIMEOUT)
    except Exception as e:
        problem.writeLP('debug.lp')
        molecule.write_graph('DEBUG', output_size=(2000, 2000))
        print('Failed LP written to "debug.lp"')
        raise

    all_tautomers = [new_molecule_for_current_solution(n=0)]

    # Remove redundant constraints from previous multi-objective optimistion
    for constraint_name in ['1_Sequence_Objective', '2_Sequence_Objective']:
        if constraint_name in problem.constraints:
            del problem.constraints[constraint_name]

    # Iterate until no more tautomers are found
    solutions = []
    for n in (count(1) if max_tautomers is None else range(1, max_tautomers)):
        # exclude current optimal solution F*: \sum_{i}|F_i - F*_i| >= 1
        problem += \
            sum(v for v in [bond_orders[bond] for bond in capping_bonds] if not v.varValue) + \
            sum(1 - v for v in [bond_orders[bond] for bond in capping_bonds] if v.varValue) >= 1, \
            'Solution {n}'.format(n=n)

        s = encode_solution()
        write_to_debug(debug, 'Excluding solution {0}'.format(s))
        solutions.append(s)

        try:
            problem.solve(timeout=ILP_SOLVER_TIMEOUT)
        except Exception as e:
            debug_failed_ILP()
            raise

        if problem.status == 1:
            pass
        else:
            write_to_debug(debug, problem.status, LpStatus[problem.status])
            break

        all_tautomers.append(new_molecule_for_current_solution(n=n))

    # check that all solutions are unique
    assert len(solutions) == len(set(solutions))
    return unique_molecules(all_tautomers)

def get_all_tautomers(
    molecule: 'Molecule',
    total_number_hydrogens: Optional[int] = None,
    net_charge: Optional[int] = None,
    enforce_octet_rule: bool = True,
    allow_radicals: bool = False,
    max_tautomers: Optional[int] = 2000,
    disallow_triple_bond_in_small_rings: bool = True,
    disallow_allenes_in_small_rings: bool = True,
    disallow_allenes_completely: bool = True,
    lock_phenyl_rings: bool = True,
    use_gurobi: bool = False,
    maximum_number_hydrogens_per_atom: Optional[int] = 3,
    skeletal_rearrangements: Optional[List[Tuple[int, int]]] = None,
    debug: Optional[TextIO] = None,
) -> 'Molecule':
    '''
    Args:
        ``disallow_triple_bond_in_small_rings``: Disallow triple bonds in small rings (rings with size <= ``SMALL_RING``).
        ``disallow_allenes_in_small_rings``: Disallow allenes (=C=) in small rings (rings with size <= ``SMALL_RING``).
        ``lock_phenyl_rings``: Prevent phenyl rings (6 membered rings with all sp3 carbons) from being hydrogenised.
        ``disallow_allenes_completely``: Disallow allenes (=C=) completely.
        ``maximum_number_hydrogens_per_atom``: Maximum number of hydrogen atoms carried by a single heavy atom. Single atom molecules will be allowed ``maximum_number_hydrogens_per_atom + 1`` hydrogens.
        ``skeletal_rearrangements``: Optional list of bonds (tuple of two atom indices) that can potentially be formed or destroyed during the tautomerisation process.

    Returns:
        A list of tautomers (Molecule).
    '''
    print_if_debug = lambda *args: write_to_debug(debug, *args)

    try:
        optional_bonds = set() if skeletal_rearrangements is None else {frozenset([atom_1, atom_2]) for (atom_1, atom_2) in skeletal_rearrangements}
    except:
        raise TypeError('Invalid skeletal_rearrangements: {0}. Example: [(1, 2), (1, 3)]'.format(skeletal_rearrangements))

    molecule.bonds |= optional_bonds

    if len([1 for atom in molecule.atoms.values() if atom.element == 'H']) > 0:
        molecule.remove_all_hydrogens(mark_all_uncapped=True)

    neighbour_counts = molecule.neighbours_for_atoms()

    def keep_capping_strategy_for_atom(capping_strategy: Capping_Strategy, atom: Atom) -> bool:
        if atom.valence is not None:
            if False:
                if neighbour_counts[atom.index] + new_atom_for_capping_strategy(capping_strategy) == atom.valence:
                    print_if_debug(atom, capping_strategy)
            return neighbour_counts[atom.index] + new_atom_for_capping_strategy(capping_strategy) == atom.valence
        else:
            return min_valence_for(atom) <= neighbour_counts[atom.index] + new_atom_for_capping_strategy(capping_strategy) <= max_valence_for(atom)

    assert maximum_number_hydrogens_per_atom >= 0, 'Maximum number of hydrogens should be greater or equal than 0'
    def possible_capping_strategies_for_atom(atom: Atom, is_phenyl_atom: bool = False) -> List[Capping_Strategy]:
        if is_phenyl_atom and lock_phenyl_rings:
            return [NO_CAP, H_CAP]
        else:
            return (
                [NO_CAP]
                + [merge_caps(* [H_CAP] * i) for i in range(1, maximum_number_hydrogens_per_atom + 1 + (1 if len(molecule.atoms) == 1 else 0))]
            )

    atoms_need_capping = [atom for atom in molecule.sorted_atoms() if not atom.capped]

    print_if_debug(
        [
            possible_capping_strategies_for_atom(atom, is_phenyl_atom=(atom.index in molecule.phenyl_atoms))
            for atom in atoms_need_capping
        ],
    )

    problem = LpProblem("Tautomer enumeration problem for molecule {0}".format(molecule.name), LpMinimize)

    ELECTRON_MULTIPLIER = (2 if not allow_radicals else 1)

    non_allene_atoms = {}
    for atom in molecule.atoms.values():
        if disallow_allenes_completely:
            atom_bonds = [bond for bond in molecule.bonds if atom.index in bond]
            if atom.element == 'C' and len(atom_bonds) == 2:
                non_allene_atoms[atom] = atom_bonds

    fragment_switches, fragment_scores, fragment_H_scores = {}, {}, {}
    capping_atoms_for = {}
    new_bonds_sets = {}
    for uncapped_atom in atoms_need_capping:
        possible_capping_strategies = possible_capping_strategies_for_atom(uncapped_atom, is_phenyl_atom=(uncapped_atom.index in molecule.phenyl_atoms))
        if len(possible_capping_strategies) == 0 or len(possible_capping_strategies) == 1 and possible_capping_strategies[0] == NO_CAP:
            pass
        else:
            for (i, capping_strategy) in enumerate(sorted(possible_capping_strategies), start=1):
                print_if_debug(uncapped_atom, capping_strategy, i)
                # Add switch variable
                fragment_switches[uncapped_atom.index, i] = LpVariable(
                    'F_{i},{j}'.format(i=uncapped_atom.index, j=i),
                    0,
                    1,
                    LpBinary,
                )

                new_atoms, new_bonds = molecule.extend_molecule_with(uncapped_atom, capping_strategy)
                print_if_debug(i, [atom for atom in new_atoms])
                capping_atoms_for[uncapped_atom.index, i] = new_atoms
                new_bonds_sets[uncapped_atom.index, i] = [bond for bond in new_bonds if uncapped_atom.index in bond]
                fragment_scores[uncapped_atom.index, i] = len(capping_atoms_for[uncapped_atom.index, i])
                fragment_H_scores[uncapped_atom.index, i] = len([atom for atom in capping_atoms_for[uncapped_atom.index, i] if atom.element == 'H'])

            # Only choose one capping strategy at a time
            problem += (lpSum(F_i for ((atom_id, _), F_i) in fragment_switches.items() if atom_id == uncapped_atom.index) == 1, 'Single capping strategy for atom {atom_desc}'.format(atom_desc=atom_short_desc(uncapped_atom)))

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

    non_capping_bonds = {
        bond for bond in molecule.bonds
        if len(bond & all_capping_atom_ids) == 0
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

    if disallow_triple_bond_in_small_rings:
        print_if_debug('Note: Excluding triple bonds in small rings (<= {0})'.format(SMALL_RING))
        bonds_in_small_rings = bonds_in_small_rings_for(molecule)
    else:
        bonds_in_small_rings = set()

    bond_orders = {
        bond: LpVariable(
            "B_{bond_key}".format(bond_key=bond_key(bond)),
            MIN_BOND_ORDER if bond not in optional_bonds else 0,
            MAX_BOND_ORDER if bond not in bonds_in_small_rings else 2,
            LpInteger,
        )
        if bond not in fragment_switch_for_bond
        else fragment_switch_for_bond[bond]
        for bond in molecule.bonds
    }

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
        problem += (lpSum(charges.values()) == net_charge, 'Known net charge={0}'.format(net_charge))

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

    for (atom, (bond_1, bond_2)) in non_allene_atoms.items():
        new_allene_switch = LpVariable('A_{i}'.format(i=atom.index), 0, 1, LpBinary)
        problem += 2 * bond_orders[bond_1] - bond_orders[bond_2] + 4 * new_allene_switch >= 3
        problem += 2 * bond_orders[bond_1] - bond_orders[bond_2] + 4 * new_allene_switch <= 5

    for atom in atoms_in_small_rings_for(molecule):
        if disallow_allenes_in_small_rings:
            if atom.element in {'C', 'N'}:
                adjacent_non_hydrogen_bonds = [bond for bond in non_capping_bonds if atom.index in bond]
                if len(adjacent_non_hydrogen_bonds) == 2:
                    problem += sum(bond_orders[bond] for bond in adjacent_non_hydrogen_bonds) <= 3, 'No allenes for atom {atom_desc} in short ring'.format(atom_desc=atom_short_desc(atom))

    def debug_failed_ILP(n: Optional[int] = None) -> None:
        debug_filename = 'debug{0}.lp'.format('' if n is None else '_{n}'.format(n=n))
        problem.writeLP(debug_filename)
        print('Failed LP written to "{0}"'.format(debug_filename))

    def encode_solution() -> int:
        # bitshift is faster than multiplication by 2**i
        return sum(int(v.varValue) << i for i, v in enumerate(fragment_switches.values()))

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
                    bond, bond_order = frozenset(map(int, variable_substr.split(','))), MUST_BE_INT(v.varValue)
                    if bond in optional_bonds:
                        if bond_order > 0:
                            print_if_debug('new_bond', bond, bond_order, [new_molecule.atoms[atom_index] for atom_index in bond])
                            new_molecule.bond_orders[bond] = bond_order
                        else:
                            print_if_debug('no_new_bond', bond, bond_order, [new_molecule.atoms[atom_index] for atom_index in bond])
                            new_molecule.bonds.remove(bond)
                    else:
                        if bond_order == 0:
                            print_if_debug('removing_bond', bond, bond_order, [new_molecule.atoms[atom_index] for atom_index in bond])
                            new_molecule.bonds.remove(bond)
                        else:
                            new_molecule.bond_orders[bond] = bond_order
            elif variable_type == 'Z':
                pass
            elif variable_type == 'A':
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

        new_molecule.update_valences()
        new_molecule.assign_aromatic_bonds()
        new_molecule.assert_molecule_coherence()
        return new_molecule

    # Solve once to find optimal solution with lowest encode_solution()
    try:
        problem.sequentialSolve(OBJECTIVES, timeout=ILP_SOLVER_TIMEOUT)
        assert problem.status == 1, (molecule.name, LpStatus[problem.status])
    except Exception as e:
        debug_failed_ILP(0)
        raise

    if debug is not None:
        debug_failed_ILP(0)

    all_tautomers = [new_molecule_for_current_solution(n=0)]

    # Remove redundant constraints from previous multi-objective optimistion
    for constraint_name in ['1_Sequence_Objective', '2_Sequence_Objective', '3_Sequence_Objective']:
        del problem.constraints[constraint_name]

    solutions = []
    # Iterate until no more tautomers are found
    for n in (count(1) if max_tautomers is None else range(1, max_tautomers)):
        # exclude current optimal solution F*: \sum_{i}|F_i - F*_i| >= 1
        problem += \
            sum(v for v in fragment_switches.values() if not v.varValue) + \
            sum(1 - v for v in fragment_switches.values() if v.varValue) >= 1,\
            'Solution {n}'.format(n=n)

        s = encode_solution()
        print_if_debug('Excluding solution {0}'.format(s))
        solutions.append(s)

        try:
            problem.solve(
                solver=GUROBI_CMD() if use_gurobi else None,
                timeout=ILP_SOLVER_TIMEOUT,
            )
        except Exception as e:
            debug_failed_ILP(n)
            raise

        if problem.status == 1:
            pass
        else:
            print_if_debug(problem.status, LpStatus[problem.status])
            if debug is not None:
                debug_failed_ILP(n)
            break

        all_tautomers.append(new_molecule_for_current_solution(n=n))

    # check that all solutions are unique
    assert len(solutions) == len(set(solutions))
    return unique_molecules(all_tautomers)
