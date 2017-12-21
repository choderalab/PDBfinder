""""

A short script to run Schrodinger's Protein Prep, relies on Schrodinger 2016-3

Written by Steven Albanese, with gratuitious borrowing from Openmoltools' Schrodinger package

"""


#################
#     Import    #
#################

from openmoltools import utils
from openmoltools import schrodinger
import os
import sys
import logging
import shutil
import csv
import subprocess
import mdtraj
from openmoltools.schrodinger import need_schrodinger

def write_file(filename, contents):
    """
    Little helper function to write the pdb files

    Args:
        filename: String, 4-letter PDB ID
        contents: string that will be written to the file

    Returns: Nothing, just writes the file

    """

    with open(filename, 'w') as outfile:
        outfile.write(contents)

@need_schrodinger
def protein_prep(input_file_path, output_file_path, pdbid, pH=7.4, fillsidechains=True, fillloops=True,
                 noepik=False, rehtreat=True, max_states=32, tolerance=0):

    # Locate PrepWizard executable
    prepwiz_path = os.path.join(os.environ['SCHRODINGER'], 'utilities', 'prepwizard')

    # Normalize paths
    input_file_path = os.path.abspath(input_file_path)
    output_file_path = os.path.abspath(output_file_path)
    output_dir = os.path.join(output_file_path, '%s-fixed' % pdbid)

    output_file_name = '../%s-fixed.pdb' % pdbid

    # Check for output file pathway
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Format arguments for PrepWizard command

    wiz_args = dict(ms=max_states, ph=pH)
    wiz_args['fillsidechains'] = '-fillsidechains' if fillsidechains else ''
    wiz_args['fillloops'] = '-fillloops' if fillloops else ''
    wiz_args['pht'] = tolerance
    wiz_args['rehtreat'] = '-rehtreat' if rehtreat else ''
    wiz_args['water_hbond_cutoff'] = 3
    wiz_args['noepik'] = '-noepik' if noepik else ''

    cmd = [prepwiz_path]
    cmd += '-mse -propka_pH {ph} {fillsidechains} {fillloops} {rehtreat} {noepik} -delwater_hbond_cutoff {water_hbond_cutoff} ' \
           '-keepfarwat -disulfides -ms {ms} -minimize_adj_h -epik_pH {ph} -f 3 -epik_pHt {pht} -fix -NOJOBID'.format(**wiz_args).split()

    cmd.append(input_file_path)
    cmd.append(output_file_name)

    with utils.temporary_cd(output_dir):
        log = schrodinger.run_and_log_error(cmd)
        write_file('%s.log' % pdbid, log)
