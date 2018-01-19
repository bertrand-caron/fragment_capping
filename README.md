[![DOI](https://zenodo.org/badge/96262215.svg)](https://zenodo.org/badge/latestdoi/96262215)

# Exhaustive fragment capping

## Introduction

Aim: Finding appropriate molecules to parametrise dihedral fragments on.

## Standard-dependencies

See `requirements.txt`. Install dependencies by running `make setup`.

* `Python >=3.5`
* `pulp`
* `networkx`
* `numpy`

## Non-standard dependencies

* [graph_tool](https://graph-tool.skewed.de): Efficient network-analysis with Python
* [chem_graph_tool](https://github.com/bertrand-caron/chem_graph_tool): Helper library to use `graph_tool` with chemical structure files (PDB format).
* [chemistry_helpers](https://github.com/bertrand-caron/chemistry_helpers): Helper library for interacting with chemical structure files (PDB files, etc.)
* [OpenBabel](http://openbabel.org): The Open Source Chemistry Toolbox

All Python module should be located in your `PYTHONPATH`.


## Problem

Given a partial molecule graph, construct molecules matching constraints
(size constraints e.g. smallest, topological constraints e.g. acyclic)
and order them based on 'chemical relevance'.

## Methods

Heuristics:
* Allowed partial charges per element
* Allowed bond order per element

## Algorithmic

Possible orders of a bond is the set intersection of the possible bond orders of both partners.

Marking uncapped atoms increases comptational efficiency (drastically decreases the number of combinations to try).
Library of capping strategies for each element
(e.g. H3, H2, H-CH2 and potentially nothing for a carbon atom, depending on whether or not uncapped atoms were marked).
Cartesian product of all possible capping strategy.
For each capped molecule, a list of all possible bond orders and charges assignment is constructed.

## Chemical relevance sorting heuristic

Valid capped molecules are sorted with the following heuristic:
* Sum of absolute values of atomic partial charges (limits charge separation)
* Number of atoms (smallest molecules)
* Number of types of double bonds (to respect electronegativity): ('CO', 'CN', 'CC')
