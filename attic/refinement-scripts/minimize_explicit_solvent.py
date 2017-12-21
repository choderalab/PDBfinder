from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmoltools.forcefield_generators import gaffTemplateGenerator
from glob import glob
import os

# Written by G. A. Ross

#------Functions used to clean and minimizing with openmm------#
def discard_organic(model, verbose=True):
    """
    Return an OpenMM modeller object that doesn't contain small organic molecules from the mother liquor

    Parameter
    ---------
    model: simtk.openmm.app.modeller.Modeller
        The modeller object with the PDB file of interest

    Returns
    -------
    model: simtk.openmm.app.modeller.Modeller
        The same object as the input with certain residues discarded
    """
    unwanted = ['DTT', 'EDO', 'GOL', 'SO4', 'PO4', 'DMS', 'ACT']
    for junk in unwanted:
        atms = [atm for atm in model.topology.atoms() if atm.residue.name == junk]
        if len(atms) > 0:
            if verbose == True: print('Deleting {0} from topology'.format(junk))
            model.delete(atms)
    return model

def get_upper_location(pdb_location, back=2):
    """
    Get the path of the directory that's up from the path specified from pdb_location.
    """
    path = pdb_location.split('/')[0:-back]
    upper_path = ''
    for location in path:
        upper_path += location + '/'
    return upper_path

def openmm_clean(pdb_filename, pdbname, out_folder='minimized', gpu=False, solvate=False):
    """
    Minimize a supplied system with openmm. Can handle small molecules as long as CONECT
    records are supplied.
    """
    # Initialize forcefield with small molecule capabilities
    forcefield = ForceField('gaff.xml','tip3p.xml','amber99sbildn.xml')
    forcefield.registerTemplateGenerator(gaffTemplateGenerator)

    # Use modeller to remove unwanted residues
    pdb = PDBFile(pdb_filename)
    model = Modeller(pdb.topology, pdb.positions)

    # Remove unwanted molecules
    model = discard_organic(model, verbose=False)

    # Add waters in a cubic box
    if solvate == True:
        model.addSolvent(forcefield, padding=1.0*nanometers)

    # Create the system with a cheap electrostatic cutoff
    system = forcefield.createSystem(model.topology, nonbondedMethod=CutoffNonPeriodic)

    # Minimize system with a placeholder integrator
    integrator = VerletIntegrator(0.001*picoseconds)
    if gpu == True:
        platform = Platform.getPlatformByName('OpenCL')
        properties = {'OpenCLPrecision': 'mixed'}
        simulation = Simulation(model.topology, system, integrator, platform, properties)
    else:
        simulation = Simulation(model.topology, system, integrator)
    simulation.context.setPositions(model.positions)
    simulation.minimizeEnergy()

    # Print PDB
    positions = simulation.context.getState(getPositions=True).getPositions()
    out_directory = get_upper_location(pdb_filename)
    try:
        os.mkdir(out_directory + out_folder)
    except OSError:
        pass
    PDBFile.writeFile(simulation.topology, positions, open(out_directory + out_folder + '/' + pdbname + '-minimized.pdb', 'w'))

if __name__ == "__main__":

    #------Generating the list of useable structures------#

    # Structures that will be omited for now due to tricky structure :
    # - 4AN2 and 4LMN have cobimetinib bound with ATP and an ion
    # - 3CJG has (presumably) post translation modifications to protein without CONECT records
    omit_pdbs = ['4AN2', '4LMN', '3CJG']

    # Structures that needed manual refinement, and are called XXXX-manual.pdb
    manual_edits = ['4G5J', '4XV2', '3ZOS', '3CJG', '3WZD', '1XKK', '5CSW', '5HIE', ' 2EUF']

    # Pre-assigment
    complete_pdbs = []     # The X-ray structures that are complete, with no modelled in loops
    missing_pdbs = []      # The X-ray structures that had missing loops longer than 20 residues, and were not modelled
    refine_pdbs = []       # All the structures that will be further refined with OpenMM
    refine_pdbs_location = [] # The relative location of the protein-ligand files that will be minimized

    # Looping through all the directories and locating the PDB structures that will be minimized with openmm
    for folder in glob('../pdbs/*'):
        for subfolder in glob(folder + '/fixed/*'):
            filename = subfolder.split('/')[-1]
            if filename.split('.')[-1] == 'pdb':
                pdb_code = filename.split('.')[0].split('-')[0]
                if pdb_code not in omit_pdbs:
                    logfile = folder + '/fixed/' +  pdb_code + '-fixed' + '/' + pdb_code.lower() + '-missing-loops.log'
                    logfile = open(logfile, 'r').read()
                    finder = logfile.find('Found loops without PDB coordinates')
                    if finder == -1:
                        # If no loops have been modeled, the structure was complete from the start!
                        complete_pdbs.append(pdb_code)
                    finder = logfile.find('is too long, skipping')
                    if finder != -1:
                        # Disregard structurs whose missing loops were too long to model in
                        missing_pdbs.append(pdb_code)
                    else:
                        refine_pdbs.append(pdb_code)
                        refine_pdbs_location.append(folder + '/fixed/' + filename)

    # Now removing the duplicate structures that have arisen due to the manaully edited PDB structures.
    index = 0
    for filename in refine_pdbs_location:
        for pdb in manual_edits:
            if  filename.split('/')[-1] == pdb + '-fixed.pdb':
                refine_pdbs.pop(index)
                refine_pdbs_location.pop(index)
        index += 1

    import logging
    logging.basicConfig(filename='explicit_solvent_minimization.log',level=logging.DEBUG)
    logging.captureWarnings(True)

    #------Cleaning and minimizing with openmm------#
    # Pre-assignment
    broken_pdbs = []
    minimized_pdbs = []
    minimized_files = []

    # Running through refined structures and seeing if they work
    for filename, pdb in zip(refine_pdbs_location, refine_pdbs):
        logging.info('Refining structure {0}'.format(pdb))
        try:
            openmm_clean(filename, pdb, out_folder='explicit_water_minimized', gpu=True, solvate=True)
            minimized_pdbs.append(pdb)
            minimized_files.append(filename)
            logging.info('Structure {0} COMPLETED.'.format(pdb))
            logging.info('Original file located in {0}'.format(pdb, filename))
            logging.info('-----------------------------')
        except Exception as detail:
            logging.error(detail)
            logging.info('Structure {0} BROKEN'.format(pdb))
            logging.info('Original file located in {0}'.format(pdb, filename))
            logging.info('-----------------------------')
            broken_pdbs.append(filename)

    with open('completed.txt', 'w') as the_file:
        the_file.writelines(["%s\n" % item  for item in minimized_pdbs])

    with open('broken.txt', 'w') as the_file:
        the_file.writelines(["%s\n" % item  for item in broken_pdbs])