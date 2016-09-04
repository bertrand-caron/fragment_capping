fragments: cache/protein_fragments.pickle
	python3 molecule_for_fragment.py
.PHONY: fragments

cache/protein_fragments.pickle:
	scp scmb-atb.biosci.uq.edu.au:/home/uqbcaron/ATB/ivan_dihedrals/protein_fragments.pickle $@
