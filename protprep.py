""""

A short script to run Schrodinger's Protein Prep

Written by Steven Albanese, with gratuitious borrowing from Openmoltools' Schrodinger package

"""

"""
/opt/schrodinger/suites2015-3/utilities/prepwizard -c -mse -fillsidechains -fillloops -propka_pH 7.4 3GCS.pdb 3GCS-sch.pdb

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


@need_schrodinger
def protein_prep(input_file_path, output_file_path, pdbid, pH=7.4, fillsidechains=True, fillloops=True, max_states=32):

    # Locate PrepWizard executable
    prepwiz_path = os.path.join(os.environ['SCHRODINGER'], 'utilities', 'prepwizard')

    # Normalize paths
    input_file_path = os.path.abspath(input_file_path)
    output_file_path = os.path.abspath(output_file_path)
    output_dir = os.path.join(os.path.dirname(output_file_path), 'fixed')

    output_file_name = '%s-fixed.pdb' % pdbid

    # Check for output file pathway
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Format arguments for PrepWizard command

    wiz_args = dict(ms=max_states, ph=pH)
    wiz_args['fillsidechains'] = '-fillsidechains' if fillsidechains else ''
    wiz_args['fillloops'] = '-fillloops' if fillloops else ''

    cmd = [prepwiz_path]
    cmd += '-c -mse -propka_pH {ph} {fillsidechains} {fillloops}'.format(**wiz_args).split()
    cmd.append(input_file_path)
    cmd.append(output_file_name)

    with utils.temporary_cd(output_dir):
        schrodinger.run_and_log_error(cmd)
