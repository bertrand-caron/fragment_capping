PYTHON_BIN_DIR = /opt/modules/Anaconda3/5.0.1/bin

test_tautomers:
	-@rm graphs/*tautomer*.pdf
	python3 test_tautomers.py
.PHONY: test_tautomers.py

test:
	python3 test_capping_wang.py
	python3 test_tautomers.py
	python3 test_capping_dihedral_fragments.py
.PHONY: test

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

MMFF94.tar.gz:
	wget http://server.ccl.net/cca/data/MMFF94/MMFF94.tar.gz -O $@

MMFF94: MMFF94.tar.gz
	mkdir $@
	gunzip -c MMFF94.tar.gz | tar xvo --directory MMFF94
.PHONY: MMFF94

setup: requirements.txt
	pip install -r requirements.txt
.PHONY: setup

test_MMFF94:
	python3 tasks/test_mmff94.py --all
.PHONY: test_MMFF94

test_SIDRUS:
	python3 tasks/test_mmff94.py --names SIDRUS JABGAU
.PHONY: test_SIDRUS
