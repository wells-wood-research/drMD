## BASIC PYTHON LIBRARIES
import os
from os import path as p
import pandas as pd
import logging
from collections import Counter

## drMD LIBRARIES
from ExaminationRoom import drLogger
from UtilitiesCloset import drSplash, drListInitiator

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

##  CLEAN CODE
from typing import Dict, Callable, Optional, List, Tuple, Set
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

"""
        PDB FILES MUST OBEY THE FOLLOWING RULES

## NO MISSING ATOMS

1. Protein chains must be unbroken
2. Amino-acid residues must have all their required atoms

## NO DUPLICATE ATOMS
3. No residues can have duplicate atoms (i.e. same ATOM_NAME field)

## OBEY RULES REGARDING CHAIN IDENTIFIERS
4. All atoms must have CHAIN_ID field
5. Each protein chain must have a unique CHAIN_ID
6. Each ligand chain must have a unique CHAIN_ID

## BE COMPATIBLE WITH AMBER PARAMETERS
7. No organometallic ligands (antechamber and tleap don't like them)
8. No non-canonical amino acids (working on this now!)
9. Ions must use the correct ATOM_NAME and RES_NAME field to be compatible with AMBER  

"""




#################################################################################################
def pdb_triage(pdbDir: DirectoryPath, config: dict) -> None:
    """
    This function iterates through all PDB files in a directory, checks for common problems,
    and logging.infos the results to the terminal.

    Args:
        pdbDir (DirectoryPath): The directory containing the PDB files.

    Returns:
        None
    """
    ## set up logging paths
    outDir: DirectoryPath = config["pathInfo"]["outputDir"]
    logDir: DirectoryPath = p.join(outDir, "00_drMD_logs")
    os.makedirs(logDir, exist_ok=True)
    pdbTriageLog: FilePath = p.join(logDir, "01_pdbTriage.log")
    ## set up logging
    drLogger.setup_logging(pdbTriageLog)


    ## get list of pdb files
    pdbNames = [p.splitext(file)[0] for file in os.listdir(pdbDir) if file.endswith(".pdb")]
    inputPdbs = [p.join(pdbDir, file) for file in os.listdir(pdbDir) if file.endswith(".pdb")]
    ## convert to dataframes
    pdbDfs = [pdbUtils.pdb2df(pdbFile) for pdbFile in inputPdbs]

    pdbDisorders = {}
    ## check for pdb files with problems
    pdbDisorders["01_broken_protein_chains"] = check_for_broken_chains(pdbDfs, pdbNames)
    pdbDisorders["02_residues_missing_atoms"] = check_for_missing_sidechains(pdbDfs, pdbNames)
    pdbDisorders["03_residues_with_duplicate_atoms"] = check_for_duplicate_atoms(pdbDfs, pdbNames)
    pdbDisorders["04_atoms_with_no_chain_id"] = check_for_missing_chain_ids(pdbDfs, pdbNames)  
    pdbDisorders["05_protein_chains_unique_chain_ids"] = check_for_termini_in_chain_middles(pdbDfs, pdbNames)
    pdbDisorders["06_ligands_and_protein_sharing_chain_ids"] = check_for_shared_chains(pdbDfs, pdbNames)
    pdbDisorders["07_organometallic_ligands"] = check_for_organometallic_ligand(pdbDfs, pdbNames)
    pdbDisorders["08_non-canonical_amino_acids"] = check_for_non_canonical_amino_acids(pdbDfs, pdbNames, pdbDir)
    pdbDisorders["09_ions_with_incorrect_names"] = check_for_ions_with_incorrect_names(pdbDfs, pdbNames)
    
    if any([len(problemPdbs) > 0 for problemPdbs in pdbDisorders.values()]):
        drSplash.print_pdb_error(pdbDisorders)



    ## deactivate logging
    drLogger.close_logging()
#################################################################################################
def check_for_ions_with_incorrect_names(pdbDfs: List[pd.DataFrame], pdbNames: List[str]) -> List[str]:
    ionResidueNames = drListInitiator.get_ion_residue_names()
    problemPdbs = []
    for pdbDf, pdbName in zip(pdbDfs, pdbNames):
        for chainId, chainDf in pdbDf.groupby(f"CHAIN_ID"):
            for resId, resDf in chainDf.groupby("RES_ID"):
                ## skip non-ion residues
                if not len(resDf) == 1:
                    continue
                resName = resDf["RES_NAME"].values[0]
                if not resName in ionResidueNames:
                    problemPdbs.append(pdbName)
                    break
                break
    return problemPdbs
#################################################################################################
def  check_for_missing_chain_ids(pdbDfs: List[pd.DataFrame], pdbNames: List[str]) -> List[str]:
    problemPdbs = []
    for pdbDf, pdbName in zip(pdbDfs, pdbNames):
        chainIds = pdbDf['CHAIN_ID'].to_list()
        if "" in chainIds:
            problemPdbs.append(pdbName)
        elif None in chainIds:
            problemPdbs.append(pdbName)

    return problemPdbs


#################################################################################################
def check_for_shared_chains(pdbDfs: List[pd.DataFrame], pdbNames: List[str]) -> List[str]:
    """
    Check to see if ligand and proteins are in the same chain

    Args:   
        pdbDf (pd.DataFrame): The pdb dataframe.

    Returns:
        isSharedChains  (bool): A boolean indicating if ligand and proteins are in the same chain
        sharedChains (Optional[Dict[str, int]]): A dictionary with the residue IDs of the shared chains and the number of non-organic atoms in each
    """
    aminoAcidResNames = drListInitiator.get_amino_acid_residue_names()
    problemPdbs = []
    for pdbDf, pdbName in zip(pdbDfs, pdbNames):
        for chainId, chainDf in pdbDf.groupby(f"CHAIN_ID"):
            chainProteinResidues = chainDf[chainDf["RES_NAME"].isin(aminoAcidResNames)].RES_ID.tolist()
            chainLigandResidues = chainDf[~chainDf["RES_NAME"].isin(aminoAcidResNames)].RES_ID.tolist()
            if len(chainProteinResidues) > 0 and len(chainLigandResidues) > 0:
                problemPdbs.append(pdbName)
   
    return problemPdbs



#################################################################################################
def check_for_organometallic_ligand(pdbDfs: List[pd.DataFrame], pdbNames: List[str]) -> List[str]:
    """
    Check for organometallic ligand in a pdb dataframe.

    Args:
        pdbDf (pd.DataFrame): The pdb dataframe.

    Returns:
        Tuple[bool, Optional[Dict[str, int]]]: A tuple containing a boolean indicating if organometallic ligand were found,
        and a dictionary with the residue IDs of the organometallic ligand and the number of non-organic atoms in each.
    """
    # Initialize lists
    aminoAcidResNames = drListInitiator.get_amino_acid_residue_names()
    
    organicElements = {"C", "N", "H", "O", "S", "P", "F", "CL",
                        "BR", "I", "SE", "B", "SI"}

    # Dictionary to store residue IDs and number of non-organic atoms
    problemPdbs = []
    for pdbDf, pdbName in zip(pdbDfs, pdbNames):
        # Loop through chains and residues in the pdb dataframe
        for chainId, chainDf in pdbDf.groupby(f"CHAIN_ID"):
            for resId, resDf in chainDf.groupby(f"RES_ID"):
                ## skip single-atom ions
                if len(resDf) == 1:
                    continue

                # Skip if amino acid residue
                resName: str = resDf["RES_NAME"].iloc[0]
                if resName in aminoAcidResNames:
                    continue
                try: 
                    resElements: Set[str] = set(resDf["ELEMENT"])
                except:
                    continue

                # If there are non-organic atoms, add the residue ID and atom count to the dictionary
                inorganicElements = [ele for ele in resElements if ele.upper() not in organicElements]
                
                if len(inorganicElements) > 0:
                    problemPdbs.append(pdbName)

    return problemPdbs


    
#################################################################################################
def check_for_non_canonical_amino_acids(pdbDfs: List[pd.DataFrame], pdbNames: List[str], inputDir: DirectoryPath) -> List[str]:
    """
    Check for non-canonical amino acids in the pdb dataframe.

    Args:
        pdbDf (pd.DataFrame): The pdb dataframe.
        inputDir (DirectoryPath): The directory path where the frcmod and mol2 files are located.

    Returns:
        Tuple[bool, Optional[Dict[str, int]]]: A tuple containing a boolean indicating if non-canonical amino acids were found,
        and a dictionary with the residue IDs of the non-canonical amino acids and a message indicating what files are missing.
    """
    # Initialize lists
    aminoAcidResNames = drListInitiator.get_amino_acid_residue_names()

    backboneAtoms: set  = {"N", "CA", "C", "O"}
    problemPdbs: List = []
    # Dictionary to store residue IDs and messages
    for pdbDf, pdbName in zip(pdbDfs, pdbNames):
        # Loop through chains and residues in the pdb dataframe
        for chainId, chainDf in pdbDf.groupby(f"CHAIN_ID"):
            for resId, resDf in chainDf.groupby(f"RES_ID"):
                resName: str = resDf["RES_NAME"].iloc[0]
                # Skip if cannonical amio acid residue, water,
                if resName in aminoAcidResNames or resName == "HOH":
                    continue
                # Skip residues with no backbone residues (i.e. ligand)
                if  not  backboneAtoms.issubset(resDf["ATOM_NAME"].unique()):
                    continue
                if look_for_ncaa_params(resName, inputDir):
                    continue
                problemPdbs.append(pdbName)

    # Return boolean indicating if non-canonical amino acids were found and the dictionary
    return problemPdbs

#################################################################################################
def look_for_ncaa_params(resName, inputDir):
    """
    Looks for AMBER parameters for a non-canonical amino acid [MOL2, FRCMOD, LIB]

    Args:
        resName (str): Residue name
        inputDir (DirectoryPath): Path to the input directory

    Returns:
        bool: True if parameters are found
    
    """
    frcmod = p.join(inputDir, f"{resName}.frcmod")
    mol2 = p.join(inputDir, f"{resName}.mol2")
    lib = p.join(inputDir, f"{resName}.lib")
    if p.exists(frcmod) and p.exists(mol2) and p.exists(lib):
        True
    else:
        False
#################################################################################################
def check_for_missing_sidechains(pdbDfs: List[pd.DataFrame], pdbNames: List[str]) -> List[str]:
    """
    Checks if a pdb dataframe contains missing sidechains.
    
    Args:
        pdbDfs (pd.DataFrame): The pdb dataframe.
        
    Returns:
        Tuple[bool, Optional[Dict[str, List[str]]]]: A tuple containing a boolean indicating if sidechains were found,
        and a dictionary with the residue IDs of the residues with missing sidechains and a list of the missing atoms.
    """
    ## initialise a set of backbone atom names and terminal oxygen name
    backboneAtoms: set = {"N", "C", "O", "CA", "OXT"}
    ## get amino acid names and dictionary containing heavy atom counts for each residue
    aminoAcidResNames = drListInitiator.get_amino_acid_residue_names()
    heavySideChainAtomCounts = drListInitiator.get_residue_heavy_atom_counts()

    problemPdbs: List = []
    for pdbDf, pdbName in zip(pdbDfs, pdbNames):
        ## get only the protein part of the pdb dataframe
        protDf = pdbDf[pdbDf["RES_NAME"].isin(aminoAcidResNames)]
        ## initialise an empty dict to store missing sidechains
        ## loop through chains and residues
        for chainId, chainDf in protDf.groupby(f"CHAIN_ID"):
            for resId, resDf in chainDf.groupby(f"RES_ID"):
                ## get residue name of this residue
                resName: str = resDf["RES_NAME"].iloc[0]
                correctHeavyAtomCount = heavySideChainAtomCounts.get(resName, None)
                if correctHeavyAtomCount is None:
                    continue
                ## get atom names of this residue
                resAtomNames: List[str] = resDf["ATOM_NAME"].tolist()
                ## exclude hydrogen atoms
                heavyAtomNames: List[str] = [atom for atom in resAtomNames if not atom.startswith(f"H")]
                ## exclude backbone atoms
                sideChainAtoms: List[str] = list(set([atom for atom in heavyAtomNames if atom not in backboneAtoms ]))
                ## check if number of heavy sidechain atoms matches expected value
                if  len(sideChainAtoms) != correctHeavyAtomCount:
                    problemPdbs.append(pdbName)
    return problemPdbs


#################################################################################################
def check_for_broken_chains(pdbDfs: List[pd.DataFrame], pdbNames: List[str]) -> List[str]:
    """
    Check for broken chains in the protein dataframe.
    
    Args:
        pdbDfs List[pd.DataFrame]: The protein dataframe.

    Returns:
        Tuple[bool, Optional[Dict[str, List[str]]]: A tuple containing a boolean indicating if broken chains were found,
        and a dictionary with chain IDs as keys and a list of broken residues or a specific message as values.
    """
    aminoAcidResNames = drListInitiator.get_amino_acid_residue_names()
    problemPdbs: List = []
    for pdbDf, pdbName in zip(pdbDfs, pdbNames):
        # Filter only the protein part of the dataframe
        protDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"].isin(aminoAcidResNames)]
        ## loop through chains
        for chainId, chainDf in protDf.groupby(f"CHAIN_ID"):
            ## get residue IDs
            resIds = chainDf["RES_ID"].unique().tolist()
            # Look for non-consecutive residue numbering
            isConsecutive, nonConsecutiveResidues = are_consecutive(resIds)
            ## update dict if needed
            if not isConsecutive:
                problemPdbs.append(pdbName)
    return problemPdbs
#################################################################################################
def check_for_termini_in_chain_middles(pdbDfs: List[pd.DataFrame], pdbNames: List[str]) -> List[str]:
    aminoAcidResNames = drListInitiator.get_amino_acid_residue_names()
    problemPdbs: List = []
    for pdbDf, pdbName in zip(pdbDfs, pdbNames):
        # Filter only the protein part of the dataframe
        protDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"].isin(aminoAcidResNames)]
        ## loop through chains
        for chainId, chainDf in protDf.groupby(f"CHAIN_ID"):
        # Look for termini in the middle of chains
            lastResIdInChain = chainDf["RES_ID"].to_list()[-1]
            residuesWithTerminalOxygenDf = chainDf[chainDf["ATOM_NAME"] == "OXT"]
            middleWithTerminalOxygenDf = residuesWithTerminalOxygenDf[residuesWithTerminalOxygenDf["RES_ID"] != lastResIdInChain]
            if len(middleWithTerminalOxygenDf) > 0:
                problemPdbs.append(pdbName)

    return problemPdbs
#################################################################################################

def are_consecutive(intList: List[int]) -> bool:
    """
    Check if the numbers in the input list are consecutive.
    
    Args:
        intList (List[int]): List of integers to check for consecutiveness.
    
    Returns:
        Tuple[bool, Optional[List[str]]: A tuple containing a boolean indicating if the numbers are consecutive,
        and a list of non-consecutive number pairs.
    """
    if not intList:
        return False  # An empty list is not considered to have consecutive numbers
    
    sortedList: List[int] = sorted(intList)
    nonConsecutives: List = []
    for i in range(len(sortedList) - 1):
        if sortedList[i + 1] - sortedList[i] != 1:
            nonConsecutives.append(f"{sortedList[i]} and {sortedList[i + 1]}")
            return False, nonConsecutives
    
    return True, None
#################################################################################################
def check_for_duplicate_atoms(pdbDfs: List[pd.DataFrame], pdbNames: List[str]) ->List[str]:
    """
    Checks for multiple conformers in a pdb dataframe.

    Args:
        pdbDf (pd.DataFrame): The pdb dataframe.

    Returns:
        Tuple[bool, Optional[Dict[str, List[str]]]]: A tuple containing a boolean indicating if multiple conformers were found,
        and a dictionary with the residue IDs of the residues with multiple conformers and a list of the duplicated atoms.
    """
    # Initialize lists
    aminoAcidResNames = drListInitiator.get_amino_acid_residue_names()
    
    problemPdbs: List = []
    for pdbDf, pdbName in zip(pdbDfs, pdbNames):
        # Loop through chains and residues in the pdb dataframe
        for chainId, chainDf in pdbDf.groupby(f"CHAIN_ID"):
            for resId, resDf in chainDf.groupby(f"RES_ID"):
                # Get residue name of this residue
                resName: str = resDf["RES_NAME"].tolist()[0]
                
                # Skip if amino acid residue
                if not resName in aminoAcidResNames:
                    continue
                
                # Get atom names of this residue
                resAtomNames: list = resDf["ATOM_NAME"].tolist()
                
                # Count the number of occurrences of each atom name
                counter: Counter = Counter(resAtomNames)
                
                # Find duplicated atoms
                duplicatedAtoms: List[str] = [atomName for atomName in resAtomNames if counter[atomName] > 1]
                
                # If there are duplicated atoms, add the residue ID and duplicated atoms to the dictionary
                if len(duplicatedAtoms) > 0:
                    problemPdbs.append(pdbName)
    return problemPdbs

#################################################################################
if __name__ == "__main__":
    ### for testing
    pdbDir = "/home/esp/scriptDevelopment/drMD/01_inputs/"
    dummyConfig = {"pathInfo": {"outputDir": "/home/esp/scriptDevelopment/drMD/03_outputs"}}
    pdb_triage(pdbDir, dummyConfig)

#################################################################################

