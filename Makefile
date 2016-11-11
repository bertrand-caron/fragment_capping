PYTHON_BIN_DIR = /home/uqbcaron/.local/bin

fragments: cache/protein_fragments.pickle
	python3 molecule_for_fragment.py
.PHONY: fragments

cache/protein_fragments.pickle:
	scp scmb-atb.biosci.uq.edu.au:/home/uqbcaron/ATB/ivan_dihedrals/protein_fragments.pickle $@

mypy: $(PYTHON_BIN_DIR)/mypy
	MYPYPATH=$(PYTHONPATH) $(PYTHON_BIN_DIR)/mypy *.py helpers/*.py
.PHONY: mypy

errors:
	PYTHONPATH=$(PYTHONPATH) $(PYTHON_BIN_DIR)/pylint -j 4 -E *.py
.PHONY: errors
