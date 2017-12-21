# Implemented by SKA, Chodera Lab
# Last edited: 9/15/17
# TO DO: If we ever publish something with this script, check how to cite PyPDB, the query and search functions
#        are based on (but very much altered from) the code in that repo

##################################################
import csv
import argparse
import xmltodict
import urllib
import urllib.request
import warnings
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile
import os
from natsort import natsorted
import pypdb
import protprep

#################################################

parser = argparse.ArgumentParser(description="Automated script to search PDB by chemical ID")
parser.add_argument('-l', required=False, dest='lig', default='None',
                    help='The ligand that you want all PDBS for. Specify when using Lig or LigAndTarget modes')

parser.add_argument('--mode', required=False, default='LigAll', dest='mode',
                    help='What mode query would you like to make. Currently implemented modes are Lig, LigAndTarget, LigAll and Apo modes')

parser.add_argument('--fix', required=False, action='store_true', dest='fix',
                    help='Setting flag fixes problems with the PDB')

parser.add_argument('--ph', required=False, default=7.4, type=float, dest='ph',
                    help='Use to set pH to something other than 7.0')

parser.add_argument('--renumber', required=False, action='store_false', dest='keepNumbers',
                    help='Set flag to renumber PDB and not use original numbering')

parser.add_argument('--biological_unit', required=False, action='store_true', dest='biological_unit',
                    help='Set flag to retrieve biological unit for all structures')

args = parser.parse_args()

ligand = args.lig
query_mode = args.mode
fix = args.fix
ph = args.ph
keepNumbers = args.keepNumbers
bunit = args.biological_unit


#########################################
#        Helper Functions               #
#########################################

def get_pdb_biological_unit(pdb_id):
    """

    Args:
        pdb_id: 4-letter PDB code

    Returns: a string with the full PDB file in it

    """

    fullurl = 'https://files.rcsb.org/download/' + pdb_id + '.pdb1'
    req = urllib.request.Request(fullurl)
    f = urllib.request.urlopen(req)
    result = f.read()
    result = result.decode('unicode_escape')
    result = result.replace('XXXX', pdb_id)

    return result


def write_file(filename, contents):
    """
    Little helper function to write the pdb files

    Args:
        filename: String, name of the file to be written
        contents: String, what will be written to the file

    Returns: Nothing, just writes the file

    """

    with open(filename, 'w') as outfile:
        outfile.write(contents)


def gen_query(search_ligand, search_protein=None, querymode=query_mode):
    """
    Args:
        search_ligand: the Chem_ID of the FDA-approved inhibitor of interest
        search_protein: UniProt Accession ID for the approved target for the inhibitor
        querymode: string indicating whether searching by ligand or for approved target/inhibitor

    Returns: xml parsed dictionary to be used with pypdb.do_search

    """

    if querymode == 'Lig':
        xml = """
        <orgPdbQuery>
         <queryType>org.pdb.query.simple.ChemCompIdQuery</queryType>
         <description>Chemical ID(s): searching for FDA approved kinase inhibitors</description>
            <chemCompId>%s</chemCompId>
            <polymericType>Any</polymericType>
        </orgPdbQuery>
        """ % search_ligand

        final_params = xmltodict.parse(xml)

    elif querymode == 'Apo':
        xml = """
        <orgPdbCompositeQuery>
         <queryRefinement>
          <queryRefinementLevel>0</queryRefinementLevel>
          <orgPdbQuery>
             <queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>
             <accessionIdList>%s</accessionIdList>
          </orgPdbQuery>
         </queryRefinement>
         <queryRefinement>
          <queryRefinementLevel>1</queryRefinementLevel>
          <conjunctionType>and</conjunctionType>
          <orgPdbQuery>
            <queryType>org.pdb.query.simple.NoLigandQuery</queryType>
            <haveLigands>no</haveLigands>
          </orgPdbQuery>
         </queryRefinement>
         <queryRefinement>
          <queryRefinementLevel>2</queryRefinementLevel>
          <conjunctionType>and</conjunctionType>
          <orgPdbQuery>
            <version>head</version>
            <queryType>org.pdb.query.simple.EnzymeClassificationQuery</queryType>
            <Enzyme_Classification>2.7.11.*</Enzyme_Classification>
          </orgPdbQuery>
          <conjunctionType>or</conjunctionType>
          <orgPdbQuery>
            <version>head</version>
            <queryType>org.pdb.query.simple.EnzymeClassificationQuery</queryType>
            <Enzyme_Classification>2.7.10.*</Enzyme_Classification>
          </orgPdbQuery>
           </queryRefinement>
        </orgPdbCompositeQuery>

        """"" % (search_protein)

        final_params = xmltodict.parse(xml)

    else:
        xml = """
        <orgPdbCompositeQuery>
         <queryRefinement>
          <queryRefinementLevel>0</queryRefinementLevel>
          <orgPdbQuery>
            <queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>
            <accessionIdList>%s</accessionIdList>
          </orgPdbQuery>
         </queryRefinement>
         <queryRefinement>
          <queryRefinementLevel>1</queryRefinementLevel>
          <conjunctionType>and</conjunctionType>
          <orgPdbQuery>
            <queryType>org.pdb.query.simple.ChemCompIdQuery</queryType>
            <chemCompId>%s</chemCompId>
          </orgPdbQuery>
         </queryRefinement>
        </orgPdbCompositeQuery>
        """ % (search_protein, search_ligand)

        final_params = xmltodict.parse(xml)

    return final_params


def search(scan_params):
    """Converts a dict() into an XML query to send to the RCSB search url. Returns a list of PDBS that needs to be
    cleaned by clean_pdb function

    Args:
        scan_params: a valid dict(), preferably generated by the gen_query function

    Returns: list of pdbs that should be cleaned up using clean_pdb function

    """

    url = 'http://www.rcsb.org/pdb/rest/search'

    queryText = xmltodict.unparse(scan_params, pretty=False)
    queryText = queryText.encode()

    req = urllib.request.Request(url, data=queryText)
    f = urllib.request.urlopen(req)
    result = f.read()

    if not result:
        warnings.warn('No results were obtained for this search')

    idlist = str(result)

    return idlist


def clean_pdb(list_pdb_ids, querymode=query_mode):
    """Cleans up the string returned by search()

    Args:
        list_pdb_ids: string, output from search()
        querymode: Type of search being performed

    Returns: A list of 4-letter PDB codes

    """
    if querymode == 'Lig':
        list_pdb_ids = list_pdb_ids.split('\\n')
        list_pdb_ids[0] = list_pdb_ids[0][-4:]
        popout = list_pdb_ids.pop(-1)

    else:
        list_pdb_ids = list_pdb_ids.split('\\n')
        list_pdb_ids[0] = list_pdb_ids[0][2:]
        popout = list_pdb_ids.pop(-1)
        for i, entry in enumerate(list_pdb_ids):
            list_pdb_ids[i] = list_pdb_ids[i][:4]

    return list_pdb_ids


def pdb_fix_pdbfixer(pdbid, file_pathway, ph, chains_to_remove):
    """

    Args:
        pdbid: 4 letter string specifying the PDB ID of the file yoou want to fix
        file_pathway: a string containing the pathway specifying how you want to organize the PDB files once written
        ph: the pH at which hydrogens will be determined and added
        chains_to_remove: dictionary containing pdbs with chains to remove
    Returns: nothing, but it does right PDB files

    """
    print(pdbid)

    # Download the topology from rcsb based on pdbod
    fixer = PDBFixer(pdbid=pdbid)

    # Remove chains based on hand curated .csv file
    if pdbid in chains_to_remove['pdbid']:
        chains = chains_to_remove['chain_to_remove'][chain_to_remove['pdbid'].index(pdbid)]
        chains_list = chains.split()
        fixer.removeChains(chainIds=chains_list)

    # Determine the first and last residue resolved in chain 0
    chains = [chain for chain in fixer.topology.chains()]
    resindices = [residue.index for residue in chains[0].residues()]
    resindices = natsorted(resindices)
    first_resindex = resindices[0]
    last_resindex = resindices[-1]

    # Find Missing residues and determine if they are C or N terminal fragments (which will be removed)

    fixer.findMissingResidues()
    if len(fixer.missingResidues) > 0:
        if sorted(fixer.missingResidues.keys())[0][-1] <= first_resindex:
            fixer.missingResidues.pop((sorted(fixer.missingResidues.keys())[0]))

        if sorted(fixer.missingResidues.keys())[-1][-1] >= last_resindex:
            fixer.missingResidues.pop((sorted(fixer.missingResidues.keys())[-1]))

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    # Write fixed PDB file, with all of the waters and ligands
    PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(file_pathway,
                                                                         '%s_fixed_ph%s.pdb' % (pdbid, ph)), 'w'),
                      keepIds=keepNumbers)

    # Remove the ligand and write a pdb file
    fixer.removeHeterogens(True)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(file_pathway,
                                                                         '%s_fixed_ph%s_apo.pdb' % (pdbid, ph)), 'w'),
                      keepIds=keepNumbers)
    # Remove the waters and write a pdb file
    fixer.removeHeterogens(False)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(file_pathway,
                                                                         '%s_fixed_ph%s_apo_nowater.pdb' % (pdbid, ph)),
                                                            'w'), keepIds=keepNumbers)


def pdb_fix_schrodinger(pdbid, file_pathway, ph):

    print(pdbid)
    input_file = os.path.join(file_pathway, '%s.pdb' % pdbid)
    output_file_pathway = os.path.join(file_pathway, 'fixed')
    protprep.protein_prep(input_file, output_file_pathway, pdbid, pH=ph)


def download_pdb(pdbid, file_pathway):
    """

    Args:
        pdbid: 4 letter string specifying the PDB ID of the file yoou want to fix
        file_pathway: a string containing the pathway specifying how you want to organize the PDB files once written

    Returns: nothing, but it does write the PDB file

    ***Note: this function does NOT fix any mistakes with the PDB file

    """

    if not os.path.exists(file_pathway):
        os.makedirs(file_pathway)

    if bunit == True:
        pdb = get_pdb_biological_unit(pdbid)

    else:
        pdb = pypdb.get_pdb_file(pdbid)

    write_file(os.path.join(file_pathway, '%s.pdb' % pdbid), pdb)


def ligand_search_mode(inhibitor_list, ligname, pH, fixpdb):
    pdb_list = []
    pathway = 'pdbs/%s' % ligname
    for id in inhibitor_list:
        print('Searching for PDBS containing %s' % id)
        query = gen_query(search_ligand=id, querymode=query_mode)
        found_pdb = search(query)
        found_pdb = clean_pdb(found_pdb)
        if len(found_pdb) > 0:
            print('found %s PDBS for %s' % (len(found_pdb), id))
            for s in found_pdb:
                download_pdb(s, pathway)
                if fixpdb is True:
                    pdb_fix_schrodinger(s, pathway, pH)


def ligand_target_search_mode(inhibitor_list, dictionary, ligname, pH, fixpdb):
    pdb_list = []
    accessions = dictionary['Accession_ID'][dictionary['inhibitor'].index(ligand)]
    accessions_list = accessions.split()
    targets = dictionary['approved_target'][dictionary['inhibitor'].index(ligand)]
    targets_list = targets.split()

    print('The FDA approved targets for %s are:' % args.lig)
    for i, ac_id in enumerate(accessions_list):  # loop through all of the ids in the human target list
        print('(%s)  %s: %s' % (i + 1, targets_list[i], ac_id))
        for id in inhibitor_list:  # Loop through all of the chem_ids for a given ligand
            print('Searching for PDBs containing %s and %s' % (id, targets_list[i]))
            query = gen_query(search_ligand=id,
                              search_protein=ac_id)  # Generates and returns a dict() for queyring RCSB
            found_pdb = search(query)  # Converts dict() to XML and query's PDB. Returns string with PDBIDs
            found_pdb = clean_pdb(found_pdb)  # Cleans up extra characters and returns a list of PDBIDs
            if len(found_pdb) > 0:
                print('found %s PDB(s) for %s/%s' % (len(found_pdb), id, targets_list[i]))
                for s in found_pdb:
                    pathway = 'pdbs/%s-%s' % (ligname, targets_list[i])
                    download_pdb(s, pathway)
                    if fixpdb is True:
                        pdb_fix_schrodinger(s, pathway, pH)


def all_ligand_search_mode(dictionary, pH, fixpdb):
    for lig in dictionary['inhibitor']:
        # Make list of ChemIDs for ligand
        chem_id_list = make_chem_id_list(dictionary, lig)

        # Make a list of accession ids associated as targets of the ligand
        accessions = dictionary['Accession_ID'][dictionary['inhibitor'].index(lig)]
        accessions_list = accessions.split()

        # Make a list of target names associated with ligand, this is for readability when saving files
        targets = dictionary['approved_target'][dictionary['inhibitor'].index(lig)]
        targets_list = targets.split()
        print('The FDA approved targets for %s are:' % lig)
        for i, ac_id in enumerate(accessions_list):
            print('(%s)  %s: %s' % (i + 1, targets_list[i], ac_id))
            for chem_id in chem_id_list:
                print('Searching for PDBs containing %s and %s' % (chem_id, targets_list[i]))
                query = gen_query(search_ligand=chem_id, search_protein=ac_id)
                found_pdb = search(query)
                found_pdb = clean_pdb(found_pdb)
                if len(found_pdb) > 0:
                    print('found %s PDB(s) for %s/%s' % (len(found_pdb), chem_id, targets_list[i]))
                    for s in found_pdb:
                        pathway = 'pdbs/%s-%s' % (lig, targets_list[i])
                        download_pdb(s, pathway)
                        if fixpdb is True:
                            pdb_fix_schrodinger(s, pathway, pH)


def apo_search_mode(dictionary, pH, fixpdb):
    accessions = dictionary['Accession_ID']
    accessions_list = set()
    for accession in range(len(accessions)):
        newlist = accessions[accession].split()
        for new_accession in newlist:
            accessions_list.add(new_accession)
    for accession_id in accessions_list:
        print('Searching for PDBs containing %s and no ligands' % accession_id)
        query = gen_query(search_ligand=ligand, search_protein=accession_id)
        found_pdb = search(query)
        found_pdb = clean_pdb(found_pdb)
        if len(found_pdb) > 0:
            print('found %s PDB(s) for %s' % (len(found_pdb), accession_id))
            for s in found_pdb:
                pathway = 'pdbs/apo/%s' % accession_id
                download_pdb(s, pathway)
                if fixpdb is True:
                    pdb_fix_schrodinger(s, pathway, pH)


def make_chem_id_list(dictionary, ligname):
    # Create list of Chem_IDs for the FDA-Approved drug name from CSV file
    ids = dictionary['Chem_ID'][dictionary['inhibitor'].index(ligname)]
    id_list = ids.split()
    if id_list == 0:
        print('Missing ChemID for %s' % ligname)

    return id_list


def convert_csv_to_dict(filename):
    reader = csv.DictReader(open(filename))
    dictionary = {}
    for row in reader:
        for column, value in row.items():
            dictionary.setdefault(column, []).append(value)

    return dictionary


if __name__ == '__main__':

    # Assert that query_mode is an implemented search typ
    assert query_mode in {'Lig', 'LigAndTarget', 'LigAll', 'Apo'}

    # Open and read csv file containing list of approved inhibitors and targets
    main_dictionary = convert_csv_to_dict('approved/clinical-kinase-inhibitors.csv')

    # Create a dictionary containing the curated PDBs that must have chains removed
    chain_to_remove = convert_csv_to_dict('remove_chains.csv')

    # Query mode Lig searches for all PDBs with a given FDA-approved kinase inhibitors in them
    if query_mode == 'Lig':
        chem_id_list = make_chem_id_list(main_dictionary, ligand)
        ligand_search_mode(chem_id_list, ligand, ph, fix)

    # Query Mode LigAndTarget searches for all inhibitor:approved target PDB files
    elif query_mode == 'LigAndTarget':
        chem_id_list = make_chem_id_list(main_dictionary, ligand)
        list_of_PDBS = ligand_target_search_mode(chem_id_list, main_dictionary, ligand, ph, fix)

    # LigAll downloads all ligands and their HUMAN targets
    elif query_mode == 'LigAll':
        all_ligand_search_mode(main_dictionary, ph, fix)

    elif query_mode == 'Apo':
        apo_search_mode(main_dictionary, ph, fix)

    else:
        warnings.warn("I think you've specified a search mode that isn't supported yet! Check --mode")
