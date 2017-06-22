from argparse import ArgumentParser
from typing import Any, List, Optional

from fragment_capping.molecule_for_fragment import molid_after_capping_fragment, Fragment, Too_Many_Permutations

def parse_args():
    parser = ArgumentParser()

    parser.add_argument('--fragment', type=Fragment, help='')
    parser.add_argument('--profile', action='store_true', help='')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    if not args.profile:
        print(
            molid_after_capping_fragment(args.fragment, quick_run=False),
        )
    else:
        from cProfile import runctx as profile_run
        from pstats import Stats

        time_file='.time'
        profile_run(
            'molid_after_capping_fragment(fragment, quick_run=True)',
            {},
            dict(molid_after_capping_fragment=molid_after_capping_fragment, fragment=args.fragment),
            filename=time_file,
        )

        stats = Stats(time_file).sort_stats('tottime', 'cumtime')
        stats.print_stats(500)
