## Large-scale kinase structure preparation

#### Initial preparation with Schrodinger
* Kinase-inhibitor complexes were retrieved using `PDBFinder.py`.
* Retrieved structures automatically prepared in Schrodinger's Maestro using
`protprep.py`. 
    * `protprep.py` modelled missing protein loops and protonates the inhibitor
    using Schroginger's 'Epic'. Schrodinger defaults used throughout [CHECK!].
    * The structures created from this procedure are named with the 
     convention `XXXX-fixed.pdb`, where `XXXX` is the PDB code.
   
#### Refinement with OpenMM
* Following visual inspection and preliminary simulations, the following
 `XXXX-fixed.pdb` structures were manually edited to create `XXXX-manual.pdb`:
    * Removed duplicate inhibitors from 4G5J such that only one inhibitor is present in the binding site
    * Removed extra copy of the inhibitor from 3ZOS such that only one inhibitor is present in the binding site
    * Renamed duplicate hydrogen atom names (H11 and H12 --> H111 and H112) of P06 in 3CJG
    * Structure A was taken of multiply resolved residue CSO in 3CJG
    * Renamed duplicate hydrogen name (H11 --> H111) of LEV in 3WZD
    * Renamed duplicate hydrogen name (H11 --> H111) of FMM in 1XKK
    * Removed odd spacing in atom names 'H9 1' and 'H9 1' in P06 in 5CSW
    * Removed odd spacing in of atom names 'H9 1' and 'H9 1' in P06 in 5HIE
    * Removed odd spacing in of atom names 'H9 1' and 'H9 1' in P06 in 4XV2 
    * Removed chain B from 4XV2
    * Removed ions ACT and CA from 2EUF
* `XXXX-fixed.pdb` or `XXXX-manual.pdb` structures were minimized with `OpenMM`
in explicit solvent using `minimize_explicit_solvent.py`. Details can be found
in that script. Briefly:
    * 4G5P, 4WKQ, 4XE0, 4ZAU had missing loops longer than 20 residues, so were ommited
    * Small organic molecules, such as DMSO, were automatically deleted
    * The `amber99sbildn`, `gaff`, and `TIP3P` forcefields were used
* The minimized structures are located in the `pdbs/*/explicit_water_minimized/` directories.

##### Structures that did not surive the pipeline
The following structures did not work with `OpenMM` refinement script, or were omitted due to tricky structures:
* Cobimetinib-MEK, 4AN2 and 4LMN. Ligand bound with ATP and a metal ion, no parameters for residue ACP
* Afatinib-EGFR, 4G5P. Loop of over 20 residues missing. Appears to be missing capping residue NME.
* Osimertinib-EGFR: 4ZAU. Loop of over 20 residues missing. Appears to be missing capping reside NME.
* Idelalisib-PI3K: 4XEO. Loop of over 20 residues missing. Appears to have minimized anyway.
* Gefitinib-EGFR: 4WKQ. Loop of over 20 residues missing. Also include post-translational modification CSX.
* Pazopanib-VEGFR1: 3CJG. Missing template for post-translational modification CSO.


##### Notes
* Structure `Imatinib-BCR-ABL/fixed/3PYY-fixed.pdb` has an activating cofactor bound, so may also be 
inappropriate with `MCCE`.
* The first attempt to refine the structures with `OpenMM` was made using `old_vacuum_refinement.ipynb`. Minimizations
were performed in vacuum, and no effort was made to fix the simulations that failed. The results of these runs are located in the `minimized/` directories.
