PYTHON_BIN_DIR = /usr/local/python35/bin

fragments: cache/protein_fragments.pickle
	python3 molecule_for_fragment.py
.PHONY: fragments

cache/protein_fragments.pickle:
	scp scmb-atb.biosci.uq.edu.au:/home/uqbcaron/ATB/ivan_dihedrals/protein_fragments.pickle $@

mypy: $(PYTHON_BIN_DIR)/mypy
	MYPYPATH=$(PYTHONPATH) $(PYTHON_BIN_DIR)/mypy *.py helpers/*.py --ignore-missing-imports
.PHONY: mypy

errors:
	PYTHONPATH=$(PYTHONPATH) $(PYTHON_BIN_DIR)/pylint -j 4 -E $$(find . -name '*.py')
.PHONY: errors

protein_fragment_molecules.png: molecule_for_fragment.py
	python3 $<
