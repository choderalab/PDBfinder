
#!/bin/bash
#BSUB -J pdbfinder
#BSUB -n 1
#BSUB -R span[ptile=1]
#BSUB -R rusage[mem=10]
#BSUB -W 10:00
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash
cd $LS_SUBCWD


# Launch job 
python PDBfinder.py -l Imatinib --mode LigAll --biological --fix --ph 7.4
