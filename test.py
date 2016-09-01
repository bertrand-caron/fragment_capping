from molecule_for_fragment import cap_fragment

def test_cyclic_fragments():
    print(cap_fragment('C,H,H|C|C|C,H,H|000'))
    print(cap_fragment('C,H,H|C|C|C,H,H|010'))
    print(cap_fragment('C,H,H|C|C|C,H,H|020'))
    print(cap_fragment('C,H,H|C|C|C,H,H|030'))

if __name__ == '__main__':
    test_cyclic_fragments()
