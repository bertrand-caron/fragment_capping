from typing import List, Optional
from copy import deepcopy
from pulp import LpProblem, LpMinimize, LpInteger, LpVariable, LpBinary, LpStatus

from fragment_capping.helpers.types_helpers import Atom, MIN, MAX
from fragment_capping.helpers.parameters import MAX_ABSOLUTE_CHARGE, MIN_ABSOLUTE_CHARGE, MAX_NONBONDED_ELECTRONS, MAX_BOND_ORDER, MIN_BOND_ORDER, VALENCE_ELECTRONS, ELECTRONS_PER_BOND, MUST_BE_INT

def get_all_tautomers(
    molecule: 'Molecule',
    total_number_hydrogens: Optional[int] = None,
    net_charge: Optional[int] = None,
    enforce_octet_rule: bool = True,
    allow_radicals: bool = False,
    max_tautomers: Optional[int] = 100,
) -> List['Molecule']:
    ELECTRON_MULTIPLIER = (2 if not allow_radicals else 1)

    molecule.remove_all_hydrogens()

    problem = LpProblem("Lewis problem (tautomers) for molecule {0}".format(molecule.name), LpMinimize)

    def maximum_number_hydrogens_for(atom: Atom) -> int:
        print(atom)
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
        problem += charges[atom.index] <= absolute_charges[atom.index], 'Absolute charge contraint 1 {i}'.format(i=atom.index)
        problem += -charges[atom.index] <= absolute_charges[atom.index], 'Absolute charge contraint 2 {i}'.format(i=atom.index)

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
            new_molecule.name += '_tautomer_{0}'.format(n)

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

        REMOVE_UNUSED_HYDROGENS = False
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

    return all_tautomers
    #return unique_molecules(all_tautomers)
