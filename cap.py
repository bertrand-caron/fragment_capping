from argparse import ArgumentParser
from typing import Any, List, Optional

from fragment_capping.molecule_for_fragment import cap_fragment, Fragment

def parse_args():
    parser = ArgumentParser()

    parser.add_argument('--fragment', type=Fragment, help='')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    print(
        cap_fragment(args.fragment),
    )
