# PDBfinder
This repository contains a collection of scripts and files to help in PDB retrieval for various projects. These scripts have been used to retrieve kinase-inhibitor complexes.

## Kinase structures
The PDB structures of kinases bound to FDA approved inhibitors can be found in `kinase-pdbs`. This folder not only contains the raw PDB structures, but also 'refined' structures that are ready for molecular simulation. The structure refinement method is detailed in `kinase-refinement/README.md`. The refinement procedure included, but was not limitied to protonation, missing residues modelling, and energy minimization in explicit solvent. 

* The cleanest structures can be found in `kinase-pdbs/*/explicit_water_minimized/`
* Automatic structure preparation is an imperfect science---use these structures with caution!

## PDB retrieval and refinement scripts
```
PDBfinder.py
```
```
protprep.py
```
```
kinase-refinement/minimize_explicit_solvent.py
```