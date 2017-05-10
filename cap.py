from argparse import ArgumentParser
from typing import Any, List, Optional

from fragment_capping.molecule_for_fragment import molid_after_capping_fragment, Fragment, Too_Many_Permutations

def parse_args():
    parser = ArgumentParser()

    parser.add_argument('--fragment', type=Fragment, help='')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    print(
        molid_after_capping_fragment(args.fragment),
    )
