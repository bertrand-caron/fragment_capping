from copy import copy, deepcopy
from argparse import ArgumentParser, Namespace

from fragment_capping.helpers.molecule import molecule_from_mol2_str
from fragment_capping.helpers.rings import bonds_for_ring

def parse_args() -> Namespace:
    parser = ArgumentParser()

    parser.add_argument('--all', action='store_true')
    parser.add_argument('--names', nargs='*', default=[])

    return parser.parse_args()

if __name__ == '__main__':
    with open('MMFF94/MMFF94_hypervalent.mol2') as fh:
        all_mol2_str = [('@<TRIPOS>MOLECULE' + mol2_str) for mol2_str in fh.read().split('@<TRIPOS>MOLECULE')[1:]]

    args = parse_args()

    if args.all:
        test_mol2_str = all_mol2_str
    elif len(args.names) > 0:
        test_mol2_str = list(filter(
            lambda mol2_str: mol2_str.splitlines()[1] in args.names,
            all_mol2_str,
        ))
    else:
        raise Exception()

    for (i, mol2_str) in enumerate(test_mol2_str, start=1):
        print(i, '/', len(test_mol2_str))
        try:
            molecule = molecule_from_mol2_str(mol2_str)
        except Exception as e:
            print(mol2_str)
            print(str(e))
        old_bond_orders = deepcopy(molecule.bond_orders)
        _, graph, pos = molecule.write_graph('RAW', graph_kwargs=dict(include_atom_index=True))

        try:
            molecule.assign_bond_orders_and_charges_with_ILP()
        except AssertionError as e:
            print('--SOLVER ERROR--')
            print(molecule, molecule.net_charge)
            print(str(e))

        new_bond_orders = deepcopy(molecule.bond_orders)
        molecule.write_graph('ILP', g=graph, pos=pos, graph_kwargs=dict(include_atom_index=True))

        def assert_bond_orders_match() -> None:
            assert old_bond_orders == new_bond_orders, (
                molecule.name,
                {bond: (old_bond_orders[bond], new_bond_orders[bond]) for bond in old_bond_orders.keys() if old_bond_orders[bond] != new_bond_orders[bond]},
            )

        try:
            assert_bond_orders_match()
        except AssertionError as e:
            print('SOFT FAIL (AROMATIC EQUIVALENCE)')
            first_bond = lambda ring: bonds_for_ring(ring)[0]
            try:
                molecule.assign_aromatic_bonds()
                # Rerun solver with constraining bond order of one bond of each aromatic ring
                molecule.assign_bond_orders_and_charges_with_ILP(
                    bond_order_constraints=[
                        (first_bond(ring), old_bond_orders[first_bond(ring)])
                        for ring in molecule.aromatic_rings
                    ],
                )
                new_bond_orders = deepcopy(molecule.bond_orders)
                assert_bond_orders_match()
            except AssertionError as e:
                print(str(e))
                print(
                    [
                        [
                            (old_bond_orders[bond], new_bond_orders[bond])
                            for bond in bonds_for_ring(ring)
                        ]
                        for ring in molecule.rings()
                    ]
                )
                molecule.assign_aromatic_bonds()
                print(molecule.aromatic_rings)

