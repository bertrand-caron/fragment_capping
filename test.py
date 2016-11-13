from fragment_capping.molecule_for_fragment import molid_after_capping_fragment

def test_cyclic_fragments():
    print(molid_after_capping_fragment('C,H,H|C|C|C,H,H|000'))
    print(molid_after_capping_fragment('C,H,H|C|C|C,H,H|010'))
    print(molid_after_capping_fragment('C,H,H|C|C|C,H,H|020'))
    print(molid_after_capping_fragment('C,H,H|C|C|C,H,H|030'))

if __name__ == '__main__':
    test_cyclic_fragments()
