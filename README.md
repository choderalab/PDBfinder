# PDBfinder
This repository contains a collection of scripts and files to help in PDB retrieval for various projects. These scripts have been used to retrieve kinase-inhibitor complexes.

## Kinase structures
The PDB structures of kinases bound to FDA approved inhibitors can be found in `pdbs`. This folder not only contains the raw PDB structures, but also 'refined' structures that are ready for molecular simulation. The structure refinement method is detailed in `kinase-refinement/README.md`. The refinement procedure included, but was not limitied to protonation and missing residues modelling
* Automatic structure preparation is an imperfect science---use these structures with caution!

## PDB retrieval and refinement scripts
```
PDBfinder.py
```
```
protprep.py
```