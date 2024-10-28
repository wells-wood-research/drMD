from typing import List, Dict


def get_residue_heavy_atom_counts():
    # Dictionary with the number of heavy side chain atoms for each amino acid residue
    heavySideChainAtomCounts: Dict[str, int] = {
        "ALA": 1,  # Alanine (CH3)
        "ARG": 7,  # Arginine (C3H6N3)
        "ASN": 4,  # Asparagine (C2H4ON)
        "ASP": 4,  # Aspartic acid (C2H4O2)
        "ASH": 4,   # Aspartic Acid protonated (C2H4O2H)
        "CYS": 2,  # Cysteine (CH2S)
        "CYX": 2,  # Cysteine in a disulphided group (CH2S)
        "CYM": 2,  # Cysteine in a disulphide group (CH2S)
        "GLN": 5,  # Glutamine (C3H6ON)
        "GLU": 5,  # Glutamic acid (C3H6O2)
        "GLH": 5,  # Glutamic acid protonated (C3H6O2H)
        "GLY": 0,  # Glycine (no side chain)
        "HIS": 6,  # Histidine (C4H5N2)
        "HIP": 6,  # Histidine 2 x protonated (C4H5N2H2)
        "HIE": 6,  # Histidine epsilon protonated (C4H5N2H)
        "HID": 6,  # Histidine 4 delta protonated (C4H5N2H)
        "ILE": 4,  # Isoleucine (C4H9)
        "LEU": 4,  # Leucine (C4H9)
        "LYS": 5,  # Lysine (C4H8N)
        "LYN": 5,  # Lysine deprotonated (C4H7N)
        "MET": 4,  # Methionine (C3H7S)
        "PHE": 7,  # Phenylalanine (C7H7)
        "PRO": 3,  # Proline (C3H6)
        "SER": 2,  # Serine (CH2O)
        "THR": 3,  # Threonine (C2H5O)
        "TRP": 10, # Tryptophan (C9H8N)
        "TYR": 8,  # Tyrosine (C7H7O)
        "VAL": 3,   # Valine (C3H7)
        ## capping groups
        "ACE": 3,  # Acetylated (COCH3)
        "NME": 2,  # N-Methylated (NCH3)
        "NHE": 1,  # N-Hydroxylated (NH2)
    }
    return heavySideChainAtomCounts

##################################################################################
def get_amino_acid_residue_names() -> set:
    """
    Returns a list of the amino acid names.
    """
    return  {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN',
            'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS',
            'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            'ASH', 'GLH', 'HIP', 'HIE', 'HID', 'CYX', 'CYM', 'LYN',   ## oddball protonations
            'ACE','NME','NHE'}              ## caps
##################################################################################
def get_ion_residue_names() -> set:
    """ 
    Returns a list of the ion atom names.
    """
    return {"Cl-", "Na+", ## counter-ions
        "LI", "NA", "K", "RB", "CS", "TL", "CU", "AG", "NH4", "H3O", "F", "CL", "BR", "I",
        "BE2", "CU2", "NI2", "PT2", "ZN2", "CO2", "PD2", "AG2", "CR2", "FE2", 
        "MG2", "V2", "MN2", "HG2", "CD2", "YB2", "CA2", "SN2", "PB2", "EU2", 
        "SR2", "SM2", "BA2", "RA2", "AL3", "FE3", "CR3", "IN3", "TL3", "Y3", 
        "LA3", "CE3", "PR3", "ND3", "SM3", "EU3", "GD3", "TB3", "DY3", "ER3", 
        "TM3", "LU3", "HF4", "ZR4", "CE4", "U4", "PU4", "TH4"}

##################################################################################
def get_backbone_atom_names() -> set:
    """
    Returns a list of the backbone atom names.
    """
    return {"N","CA","C","O"}
##################################################################################
def get_solvent_residue_names() -> set:
    """
    Returns a list of the solvent residue names.
    """
    return {"HOH", "WAT"}
##################################################################################
def get_not_a_run_dir() -> set:
    """
    Returns a list of the not a run directory names.
    """
    return {"00_AutoMethods",
             "00_clustered_pdbs",
               "00_configs",
                 "00_drMD_logs",
                   "00_vitals_reports",
                     "01_ligand_parameters",
                       "00_collated_pdbs"}