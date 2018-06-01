
#!/bin/bash
#BSUB -J pdbfinder
#BSUB -n 1
#BSUB -R span[ptile=1]
#BSUB -R rusage[mem=10]
#BSUB -W 144:00
#BSUB -o %J.stout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

cd $LS_SUBCWD
python PDBfinder.py --mode LigAll -l Imatinib --fix --ph 7.4 --biological_unit
