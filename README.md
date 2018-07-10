[![DOI](https://zenodo.org/badge/96262215.svg)](https://zenodo.org/badge/latestdoi/96262215)

# `fragment_capping`: An Integer Linear Program (ILP) for predicting Lewis structures and applications to capping molecular fragments and enumerating tautomers

## Introduction

This `python` module implements several graph-based Integer Linear Programs (ILPs)
which are described in the paper:

```
Caron, Bertrand; Engler, Martin; Mark, Alan E. and Klau, Gunnar W.
An Integer Linear Program (ILP) for predicting Lewis structures and applications to capping molecular fragments and enumerating tautomers
```

## Standard-dependencies

See `requirements.txt`. Install dependencies by running `make setup`.

* `Python >=3.6`
* `pulp`
* `networkx`
* `numpy`

## Non-standard dependencies

* [graph_tool](https://graph-tool.skewed.de): Efficient network-analysis with Python
* [chem_graph_tool](https://github.com/bertrand-caron/chem_graph_tool): Helper library to use `graph_tool` with chemical structure files (PDB format).
* [chemistry_helpers](https://github.com/bertrand-caron/chemistry_helpers): Helper library for interacting with chemical structure files (PDB files, etc.)
* [OpenBabel](http://openbabel.org): The Open Source Chemistry Toolbox

All Python module should be located in your `PYTHONPATH`.

## Algorithms and examples

### Lewis structure prediction
Given a complete molecular graph, returns the best Lewis structure of a molecule by distributing valence electrons to minimise the sum of absolute charge of the molecule and respecting relative atomic electronegativities.

```
from fragment_capping.helpers.molecule import Molecule, Atom

molecule = Molecule(
	atoms=[
		Atom(element='C', capped=True, index=1, valence=2, coordinates=(0, 0, 0)),
		Atom(element='O', capped=True, index=2, valence=1, coordinates=(1, 0, 0)),
		Atom(element='O', capped=True, index=3, valence=1, coordinates=(-1, 0, 0)),
	],
	bonds=[(1, 2), (1, 3)],
	name='CO2',
)

molecule.assign_bond_orders_and_charges_with_ILP()

print(molecule.mol2())
```

### Fragment capping
Given a partial molecule graph (**fragment**), find the smallest, most aliphatic, least charged molecule containing this fragment.

```
from fragment_capping.helpers.molecule import Molecule, Atom

fragment = Molecule(
	atoms=[
		Atom(element='C', capped=False, index=1, valence=None, coordinates=[0, 0, 0]),
		Atom(element='O', capped=True, index=2, valence=1, coordinates=(1, 0, 0)),
	],
	bonds=[(1, 2)],
	name='CO',
)

fragment.get_best_capped_molecule_with_ILP()

print(fragment.mol2())
```

### Tautomer enumeration
Tautomers are compounds that readily interconvert by the movement of an atom (usually hydrogen) or group of atoms from one site to another within the molecular structure (Katritzky **et al.**, 2010).
Given a molecule, enumerate all possible tautomers.

```
from fragment_capping.helpers.molecule import Molecule, Atom

molecule = Molecule(
	atoms=[
		Atom(element='C', capped=True, index=1, valence=3, coordinates=(0, 0, 0)),
		Atom(element='O', capped=True, index=2, valence=1, coordinates=(1, 0, 0)),
		Atom(element='H', capped=True, index=3, valence=1, coordinates=(0, 1, 0)),
		Atom(element='H', capped=True, index=4, valence=1, coordinates=(0, -1, 0)),
	],
	bonds=[(1, 2), (1, 3), (1, 4)],
	name='CH2O',
)

for tautomer in molecule.get_all_tautomers():
	print(tautomer.mol2())
```
