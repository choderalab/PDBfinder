# Script to renumber PDBs in the Hauser benchmark set
# Steven K. Albanese, with help from Stack overflow posts by Gordon Wells and Maximillian Peters

import argparse
import os
from glob import glob

from Bio import AlignIO, SeqIO, ExPASy, SwissProt
from Bio import PDB
from Bio import pairwise2
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat.MatrixInfo import blosum62


def renumber(file_name, output_dir, accession_id, cap=True, lig_name=None):
    # read in PDB
    pdb_io = PDB.PDBIO()

    # File name will point to the name of the file without .pdb at the end
    filename = os.path.basename(file_name)
    print('Processing %s' % filename)
    id_name = filename[0:-4]

    # PDB name will point to the actual file
    pdb_name = file_name
    pdbparser = PDB.PDBParser(PERMISSIVE=1)
    structure = pdbparser.get_structure(id_name, pdb_name)

    # Extract sequence from PDB file and skip cap residues (which won't align properly)
    oneletter = {
        'ASP': 'D', 'GLU': 'E', 'ASN': 'N', 'GLN': 'Q',
        'ARG': 'R', 'LYS': 'K', 'PRO': 'P', 'GLY': 'G',
        'CYS': 'C', 'THR': 'T', 'SER': 'S', 'MET': 'M',
        'TRP': 'W', 'PHE': 'F', 'TYR': 'Y', 'HIS': 'H',
        'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I',
        'PTR': 'T'
    }

    cap_res = ['NME', 'NMA', 'ACE']

    pdbseq_list = []
    pdb_len = 0
    for residue in structure.get_residues():
        three_letter = residue.get_resname()
        if three_letter not in oneletter:
            pass
        else:
            pdb_len += 1
            one_letter_name = oneletter[three_letter]
            pdbseq_list.append(one_letter_name)

    pdbseq_str = ''.join(pdbseq_list)

    # Write out FASTA file for PDB structure
    alnPDBseq = SeqRecord(Seq(pdbseq_str, IUPAC.protein), id=file_name)
    SeqIO.write(alnPDBseq, "%s.fasta" % file_name, "fasta")

    # Retrieve reference sequence
    accession = accession_id  # eventually this can be specified by user somehow (or read in from the csv file)
    handle = ExPASy.get_sprot_raw(accession)
    swissseq = SwissProt.read(handle)
    refseq = SeqRecord(Seq(swissseq.sequence, IUPAC.protein), id=accession)
    SeqIO.write(refseq, "%s.fasta" % accession, "fasta")

    # Align using Biopython's EMBOSS needle wrapper
    seq1 = SeqIO.read("%s.fasta" % accession, "fasta")  # Uniprot  sequence
    seq2 = SeqIO.read("%s.fasta" % file_name, "fasta")  # PDB sequence

    alignments = pairwise2.align.localds(seq1.seq, seq2.seq, blosum62, -10, -0.5)

    if len(alignments) == 0:
        print('%s does not align to %s. Move to the next structure' % (file_name, accession_id))
        return None
    else:
        alignment_start = alignments[0][3] + 1
        alignment_end = alignments[0][4] + 1

    # Clean up FASTA files
    os.remove("%s.fasta" % file_name)
    os.remove("%s.fasta" % accession)

    # Create the list of new residue numbers from the alignment
    new_resnums = list(range(alignment_start, alignment_start + pdb_len))
    # Check if the renumbering is necessary
    if cap:
        first_in_structure = list(structure.get_residues())[1].get_id()[1]
    else:
        first_in_structure = list(structure.get_residues())[0].get_id()[1]
    if first_in_structure == alignment_start:
        print("%s does not need to be renumbered! Just checking the cap residues" % file_name)
        for i, residue in enumerate(structure.get_residues()):
            three_letter = residue.get_resname()
            if three_letter == 'ACE':  # Set resid for ACE to 1 before the start of the alignment
                res_id = list(residue.id)
                res_id[1] = new_resnums[0] - 1
                if residue.id != tuple(res_id):
                    residue.id = tuple(res_id)
                residue.resname = 'ACE'
            elif three_letter == 'NME' or three_letter == 'NMA':  # Set resid for NME or NMA to last residue number
                res_id = list(residue.id)
                res_id[1] = new_resnums[-1]
                if residue.id != tuple(res_id):
                    residue.id = tuple(res_id)
                residue.resname = 'NME'
            elif three_letter == lig_name:
                res_id = list(residue.id)
                res_id[1] = new_resnums[-1] + 1
                if residue.id != tuple(res_id):
                    residue.id = tuple(res_id)
            else:
                pass
        pdb_io.set_structure(structure)
        output_filename = os.path.join(output_dir, id_name + '-uniprot.pdb')
        pdb_io.save(output_filename)

    elif first_in_structure > alignment_start:
        print("%s does need to be renumbered" % file_name)
        # Renumber the residue IDs in the structure
        for i, residue in enumerate(structure.get_residues()):
            three_letter = residue.get_resname()
            if three_letter == 'ACE':  # Set resid for ACE to 1 before the start of the alignment
                res_id = list(residue.id)
                res_id[1] = new_resnums[0] - 1
                if residue.id != tuple(res_id):
                    residue.id = tuple(res_id)
                residue.resname = 'ACE'
            elif three_letter == 'NME' or three_letter == 'NMA':  # Set resid for NME or NMA to last residue number
                res_id = list(residue.id)
                res_id[1] = new_resnums[-1]
                if residue.id != tuple(res_id):
                    residue.id = tuple(res_id)
                residue.resname = 'NME'
            elif three_letter == lig_name:
                res_id = list(residue.id)
                res_id[1] = new_resnums[-1] + 1
                if residue.id != tuple(res_id):
                    residue.id = tuple(res_id)
            elif three_letter not in cap_res and three_letter not in oneletter:
                pass
            else:
                if cap:
                    index = i - 1
                else:
                    index = i
                res_id = list(residue.id)
                res_id[1] = new_resnums[index]
                residue.id = tuple(res_id)
        #  Write the renumbered PDB file

        pdb_io.set_structure(structure)
        output_filename = os.path.join(output_dir, id_name + '-uniprot.pdb')
        pdb_io.save(output_filename)


    elif first_in_structure < alignment_start:
        print("%s does need to be renumbered" % file_name)
        temp_resnums = list(range(alignment_end, alignment_end + pdb_len))
        print(len(new_resnums))
        for i, residue in enumerate(structure.get_residues()):
            three_letter = residue.get_resname()
            if cap:
                index = i - 1
            else:
                index = i

            if three_letter not in cap_res and three_letter not in oneletter:
                pass

            else:
                res_id = list(residue.id)
                res_id[1] = temp_resnums[index]
                residue.id = tuple(res_id)

        for i, residue in enumerate(structure.get_residues()):
            three_letter = residue.get_resname()
            if three_letter == 'ACE':  # Set resid for ACE to 1 before the start of the alignment
                res_id = list(residue.id)
                res_id[1] = new_resnums[0] - 1
                if residue.id != tuple(res_id):
                    residue.id = tuple(res_id)
                residue.resname = 'ACE'
            elif three_letter == 'NME' or three_letter == 'NMA':  # Set resid for NME or NMA to last residue number
                res_id = list(residue.id)
                res_id[1] = new_resnums[-1]
                if residue.id != tuple(res_id):
                    residue.id = tuple(res_id)
            elif three_letter == lig_name:
                res_id = list(residue.id)
                res_id[1] = new_resnums[-1] + 1
                if residue.id != tuple(res_id):
                    residue.id = tuple(res_id)
            elif three_letter not in cap_res and three_letter not in oneletter:
                pass
            else:
                if cap:
                    index = i - 1
                else:
                    index = i
                res_id = list(residue.id)
                res_id[1] = new_resnums[index]
                residue.id = tuple(res_id)


        #  Write the renumbered PDB file

        pdb_io.set_structure(structure)
        output_filename = os.path.join(output_dir, id_name + '-uniprot.pdb')
        pdb_io.save(output_filename)

    return None


