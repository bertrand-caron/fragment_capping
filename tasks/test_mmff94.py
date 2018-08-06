from copy import copy, deepcopy
from argparse import ArgumentParser, Namespace
from itertools import groupby
from operator import itemgetter
from datetime import datetime
from typing import Any, Dict, Tuple, Optional
from os.path import join, exists
from re import match

from fragment_capping.helpers.molecule import molecule_from_mol2_str
from fragment_capping.helpers.rings import bonds_for_ring
from fragment_capping.helpers.graphs import lewis_graph, are_graphs_isomorphic

CONSTRAINT_TOTAL_CHARGE = True

def parse_args() -> Namespace:
    parser = ArgumentParser()

    parser.add_argument('--all', action='store_true')
    parser.add_argument('--names', nargs='*', default=[])

    return parser.parse_args()

def molecule_name(mol2_str: str) -> str:
    return mol2_str.splitlines()[1]

def plot_timings(timings: Dict[str, Tuple[int, int, Any]]) -> None:
    DATA = [
        (atom_number, bond_number, timing.total_seconds()) for (atom_number, bond_number, timing) in timings.values()
    ]

    xs, ys, zs = zip(*DATA)

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from numpy import linspace, mean
    from scipy.stats import linregress

    WRITE_TO_FILE = True

    def plot_3D():
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(xs, ys, zs)

        ax.set_xlabel('Atom Number')
        ax.set_ylabel('Bond Number')
        ax.set_zlabel('Timing (s)')

        if not WRITE_TO_FILE:
            plt.show()
        else:
            plt.savefig('graphs/timings_3d.pdf')

    def plot_2D():
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(xs, zs, alpha=0.5, s=5.0)
        ax.set_xlabel('Atom Number')
        ax.set_ylabel('Timing (s)')

        slope, intercept, r_value, *_ = linregress(xs, zs)
        x_range = min(xs), max(xs)
        fine_xs = linspace(*x_range, 100)
        ax.plot(fine_xs, intercept + fine_xs * slope, color='red')

        if not WRITE_TO_FILE:
            plt.show()
        else:
            plt.savefig('graphs/timings_2d.pdf')

    plot_3D()
    plot_2D()

    print(mean(zs))

def try_plotting_molecule(mol2_str: str) -> None:
    try:
        from chemistry_helpers.babel import babel_output
        svg_filepath = join('graphs', '{0}.svg'.format(molecule_name(mol2_str)))
        if not exists(svg_filepath):
            with open(svg_filepath, 'wt') as fh:
                fh.write(babel_output(mol2_str, in_format='mol2', out_format='svg'))
    except:
        raise
        pass

if __name__ == '__main__':
    with open('MMFF94/MMFF94_hypervalent.mol2') as fh:
        all_mol2_str = [('@<TRIPOS>MOLECULE' + mol2_str) for mol2_str in fh.read().split('@<TRIPOS>MOLECULE')[1:]]

    assert len(all_mol2_str) == 761, len(all_mol2_str)

    args = parse_args()

    if args.all:
        test_mol2_str = all_mol2_str
    elif len(args.names) > 0:
        test_mol2_str = list(filter(
            lambda mol2_str: molecule_name(mol2_str) in args.names,
            all_mol2_str,
        ))
    else:
        raise Exception()

    SAME_BOND_ORDER, DIFFERENT_BOND_ORDER, WRONG_MOL2, SOLVER_FAILURE = 'S', 'F', 'F_MOL2', 'F_SOLVER'
    NO_SECOND_RUN = '-'
    SAME_NET_CHARGE, DIFFERENT_NET_CHARGE = 'C', 'NC'
    status, timings = {}, {}
    try:
        for (i, mol2_str) in enumerate(test_mol2_str, start=1):
            assert len(status) == i - 1, (len(status), i - 1)
            print(i, '/', len(test_mol2_str))
            try_plotting_molecule(mol2_str)
            try:
                molecule = molecule_from_mol2_str(mol2_str)
            except KeyboardInterrupt:
                raise
            except Exception as e:
                print(mol2_str)
                print(str(e))
                status[molecule_name(mol2_str)] = (WRONG_MOL2,) + ((None, None),) * 2
                continue

            net_charge = int(
                round(
                    sum(
                        float(line.split()[8])
                        for line in mol2_str.splitlines()
                        if len(line.split()) == 9 and match('^ *[0-9]+ ', line) is not None
                    ),
                )
            )

            molecule.net_charge = net_charge if CONSTRAINT_TOTAL_CHARGE else None
            old_molecule = deepcopy(molecule)
            old_bond_orders = deepcopy(molecule.bond_orders)
            _, graph, pos = molecule.write_graph(
                'RAW',
                graph_kwargs={'include_atom_index': False, 'vertex_color_scheme': 'elements', 'vertex_label_template': ''},
            )


            try:
                t1 = datetime.now()
                molecule.assign_bond_orders_and_charges_with_ILP()
                timing = datetime.now() - t1
                timings[molecule.name] = (len(molecule.atoms), len(molecule.bonds), timing)
            except KeyboardInterrupt:
                raise
            except Exception as e:
                print('--SOLVER ERROR--')
                print(molecule, molecule.net_charge)
                print(str(e))
                print('----')
                status[molecule.name] = (SOLVER_FAILURE,) + ((None, None),) * 2
                continue

            new_bond_orders = deepcopy(molecule.bond_orders)
            molecule.write_graph(
                'ILP',
                g=graph,
                pos=pos,
                graph_kwargs={'include_atom_index': False, 'vertex_color_scheme': 'elements', 'vertex_label_template': '{charge_str}'},
            )

            def assert_bond_orders_match(e: Optional['AssertionError']) -> None:
                assert are_graphs_isomorphic(
                    [
                        lewis_graph(
                            molecule,
                            # Remove non-bonded electrons since the reference do not have them
                            use_non_bonded_electrons=False,
                        )
                        for molecule in  (old_molecule, molecule)
                    ],
                ), (
                    molecule.name,
                    {bond: (old_bond_orders[bond], new_bond_orders[bond]) for bond in old_bond_orders.keys() if old_bond_orders[bond] != new_bond_orders[bond]},
                    SAME_NET_CHARGE if molecule.netcharge() == net_charge else DIFFERENT_NET_CHARGE,
                )
                if e is None:
                    status[molecule.name] = (
                        SAME_BOND_ORDER,
                        (0, (SAME_NET_CHARGE if molecule.netcharge() == net_charge else DIFFERENT_NET_CHARGE)),
                        (NO_SECOND_RUN, None),
                    )
                else:
                    status[molecule.name] = (
                        SAME_BOND_ORDER,
                        (len(e.args[0][1]), e.args[0][2]),
                        (0, (SAME_NET_CHARGE if molecule.netcharge() == net_charge else DIFFERENT_NET_CHARGE)),
                    )

            try:
                assert_bond_orders_match(None)
            except AssertionError as e:
                print('SOFT FAIL (AROMATIC EQUIVALENCE)')
                first_bond = lambda ring: bonds_for_ring(ring)[0]
                try:
                    molecule.assign_aromatic_bonds()
                    bond_order_constraints = [
                        (first_bond(ring), old_bond_orders[first_bond(ring)])
                        for ring in molecule.aromatic_rings
                    ]
                    if False:
                        print(f'PDB: {molecule.pdb()}')
                        print(f'aromatic_rings: {molecule.aromatic_rings}')
                        print(f'bond_order_constraints: {bond_order_constraints}')
                    # Rerun solver with constraining bond order of one bond of each aromatic ring
                    molecule.assign_bond_orders_and_charges_with_ILP(
                        bond_order_constraints=bond_order_constraints,
                    )
                    new_bond_orders = deepcopy(molecule.bond_orders)
                    assert_bond_orders_match(e)
                except AssertionError as f:
                    print(str(e))
                    status[molecule.name] = (
                        DIFFERENT_BOND_ORDER,
                        (len(e.args[0][1]), e.args[0][2]),
                        (len(f.args[0][1]), f.args[0][2]),
                    )
    finally:
        get_molecule_name, on_status = itemgetter(0), itemgetter(1)
        grouped = {
            status: len(list(group))
            for (status, group) in
            groupby(
                sorted(
                    status.items(),
                    key=on_status,
                ),
                key=on_status,
            )
        }
        print(grouped)
        assert sum(grouped.values()) == len(test_mol2_str), (sum(grouped.values()), len(test_mol2_str))
        plot_timings(timings)
        print(status)
