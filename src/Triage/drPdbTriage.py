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

    # Initialize a dictionary to store common PDB problems
    commonPdbProblems : Dict[bool] = {
        "isMultipleConformers": False,  # PDB contains multiple conformers
        "isBrokenChains": False,  # PDB contains broken chains
        "isMissingSidechains": False,  # PDB contains missing sidechains
        "isNonCanonicalAminoAcids": False,  # PDB contains non-canonical amino acids
        "isOrganimetallicLigands": False,  # PDB contains organimetallic ligand
        "isSharedChains": False  # Chains shared between ligand and protein
    }

    # Iterate through all PDB files in the directory
    for file in os.listdir(pdbDir):
        if not file.endswith(f".pdb"):
            continue
        pdbFile: FilePath = p.join(pdbDir, file)
        # Check for common problems in the PDB file
        problemsDict = pdb_triage_protocol(pdbFile, pdbDir, config)
        # Update the dictionary that looks at all PDB files
        commonPdbProblems = update_problem_dict(commonPdbProblems, problemsDict)
    # Print the results to the terminal
    report_problems(commonPdbProblems, pdbTriageLog)

    ## deactivate logging
    drLogger.close_logging()
#################################################################################################
def update_problem_dict(commonPdbProblems: Dict[str,bool], thisProblemDict: Dict[str,bool]) -> Dict[str,bool]:
    """
    Updates the commonPdbProblems dictionary with the values from thisProblemDict. 
    If a value in thisProblemDict is True, it means that the corresponding problem is present in the PDB file,
    so it is set to True in the commonPdbProblems dictionary as well.
    
    Args:
        commonPdbProblems (Dict[str, bool]): A dictionary containing the common problems found in all PDB files.
        thisProblemDict (Dict[str, bool]): A dictionary containing the problems found in a specific PDB file.
        
    Returns:
        Dict[str, bool]: The updated commonPdbProblems dictionary.
    """
    
    # Iterate through all keys in thisProblemDict
    for key, value  in thisProblemDict.items():
        # If the value in thisProblemDict is True, it means that the corresponding problem is present in the PDB file,
        # so set the value in commonPdbProblems to True as well
        if value:
            commonPdbProblems[key] = True
    
    # Return the updated commonPdbProblems dictionary
    return commonPdbProblems
#################################################################################################
def report_problems(commonPdbProblems: Dict[str, bool], pdbTriageLog: FilePath) -> None:
    """
    Prints out the common problems found in the PDB files based on the commonPdbProblems dictionary.
    
    Args:
        commonPdbProblems (Dict[str, bool]): A dictionary containing the common problems found in the PDB files.
    """
    logging.info(f"\n\n")
    if any(commonPdbProblems.values()):
        logging.info(f"The following common problems were found in the PDB files:")
        if commonPdbProblems["isMultipleConformers"]:
            logging.info(f"\n  * Multiple conformers found for sidechains of residues *")
            logging.info(f"\t> This often occurs in X-ray structures when electron density is found for multiple conformers")
            logging.info(f"\t> You can fix this in Pymol with the following command:")
            logging.info(f"\t> remove not (alt '' or alt A)")

        if commonPdbProblems["isBrokenChains"]:
            logging.info(f"\n  * Problems with gaps in the protein's backbone *")
            logging.info(f"\t> This often occurs with loops in X-ray structures")
            logging.info(f"\t> If you know the sequence of the protein, you can fix this in RF-Diffusion")

        if commonPdbProblems["isMissingSidechains"]:
            logging.info(f"\n  * Some Residues are missing sidechain atoms *")
            logging.info(f"\t> This often occurs in X-ray structures in areas of low electron density")
            logging.info(f"\t> If you know the sequence of the protein, you can fix this in Scwrl4")

        if commonPdbProblems["isNonCanonicalAminoAcids"]:
            logging.info(f"\n  * non-canonical amino acids have been identified *")
            logging.info(f"\t> This will cause parameterisation of your system to fail")
            logging.info(f"\t> You can create your own parameters for non-canonical amino acids")
            logging.info(f"\t> and supply them in the inputs directory")

        if commonPdbProblems["isOrganimetallicLigands"]: 
            logging.info(f"\n  * Organimetallic ligand have been identified *")
            logging.info(f"\t> This will cause parameterisation of your system to fail")
            logging.info(f"\t> You can create your own parameters for organometallic ligand")
            logging.info(f"\t> and supply them in the inputs directory")
        if commonPdbProblems["isSharedChains"]:
            logging.info(f"\n  * Chains with both protein and ligand residues detected *")
            logging.info(f"\t> This will cause parameterisation of your system to fail")
    else:
        logging.info(f"No common problems found in the PDB files.")

    if any(commonPdbProblems.values()):
        drSplash.print_pdb_error()
        drLogger.log_info(f"\n\n")
        drLogger.log_info(f"Problems with the PDB files will cause parameterisation to fail", True, True)
        drLogger.log_info(f"Consult the following log file for more details:", True, True)
        drLogger.log_info(f"\t{pdbTriageLog}", True, True)
        exit(1)

#################################################################################################

def pdb_triage_protocol(pdbFile: FilePath, inputDir: DirectoryPath, config: dict) -> Dict[str,bool]:
    """
    Runs before drMD main protocol
    Checks input pdb files for common problems that will cause drMD to crash
    These include:
    1. Broken chains (usually caused by unresolved loops in X-ray structure)
        - these will cause TLEAP to fail 
    2. Non-Canonical Amino Acids
        - these will cause TLEAP to fail
    3. Organimetallic Ligands
        - this will cause parmchk to fail
    4. Incorrect Chain IDs
        - if multiple different chains have the same chain ID, TLEAP will fail
    5. Residues with multiple conformers
        - If duplicate atoms exist in a residue, TLEAP will fail
    Args:   
        pdbFile (str): Path to the input PDB file.
    Returns:
        Dict[str,bool]: A dictionary containing the common problems found in the PDB file.
    """
    pdbName: str = p.basename(pdbFile)
    logging.info(f"\nChecking PDB file {pdbName} for common problems...")
    
    ## load pdb file as a dataframe
    pdbDf: pd.DataFrame = pdbUtils.pdb2df(pdbFile)

    ## check for non-canonical amino acids
    isNonCanonicalAminoAcids, nonCanonicalAminoAcids = check_for_non_canonical_amino_acids(pdbDf, inputDir, config)

    ## check for multiple conformers
    isMultipleConformers, multipleConformers = check_for_multiple_conformers(pdbDf)

    ## check for broken chains
    isBrokenChains, brokenChains = check_for_broken_chains(pdbDf)

    ## check for missing sidechains
    isMissingSidechains, missingSideChains = check_for_missing_sidechains(pdbDf)

    ## check for organimetallic ligand
    isOrganimetallicLigands, organimetallicLigands = check_for_organometallic_ligand(pdbDf)

    isSharedChains, sharedChains = check_for_shared_chains(pdbDf, config)

    ## report any problems found in pdb file
    if isMultipleConformers:
        logging.info(f"  * Multiple conformers found in {pdbName} for the following residues: *")
        for key, value in multipleConformers.items():
            logging.info(f"\t\t{key}: {value}")
    if isBrokenChains:
        logging.info(f"  * Broken chains found in {pdbName} in the following chains, between these residues: *")
        for key, value in brokenChains.items():
            logging.info(f"\t\t{key}: {value}")
    if isMissingSidechains:
        logging.info(f"  * Missing sidechain atoms found in {pdbName} for the following residues: *")
        for key, value in missingSideChains.items():
            logging.info(f"\t\t{key}: {value}")
    if isNonCanonicalAminoAcids:
        logging.info(f"  * Non-canonical amino acids found in {pdbName} for the following residues: *")
        for key, value in nonCanonicalAminoAcids.items():
            logging.info(f"\t\t{key}: {value}")
    if isOrganimetallicLigands:
        logging.info(f"  * Organimetallic ligand found in {pdbName} for the following residues: *")
        for key, value in organimetallicLigands.items():
            logging.info(f"\t\t{key}: {value}")
    if isSharedChains:
        logging.info(f"  * Shared chains found in {pdbName} for the following chains: *")
        for key, value in sharedChains.items():
            logging.info(f"\t\t{key}: {value}")


    problemsDict: Dict[str,bool] = {
        "isMultipleConformers": isMultipleConformers,
        "isBrokenChains": isBrokenChains,
        "isMissingSidechains": isMissingSidechains,
        "isNonCanonicalAminoAcids": isNonCanonicalAminoAcids,
        "isOrganimetallicLigands": isOrganimetallicLigands,
        "isSharedChains": isSharedChains
    }

    if not any(problemsDict.values()):
        logging.info(f"  * No common problems found in the PDB file *")

    return problemsDict

def check_for_shared_chains(pdbDf: pd.DataFrame, config: dict) -> tuple[bool, Optional[Dict[str, int]]]:
    """
    Check to see if ligand and proteins are in the same chain

    Args:   
        pdbDf (pd.DataFrame): The pdb dataframe.

    Returns:
        isSharedChains  (bool): A boolean indicating if ligand and proteins are in the same chain
        sharedChains (Optional[Dict[str, int]]): A dictionary with the residue IDs of the shared chains and the number of non-organic atoms in each
    """
    isSharedChains = False
    aminoAcidResNames = drListInitiator.get_amino_acid_residue_names()
    sharedChains: Dict = {}
    for chainId, chainDf in pdbDf.groupby(f"CHAIN_ID"):
        chainResidues = chainDf["RES_ID"].tolist()
        chainProteinResidues = chainDf[chainDf["RES_NAME"].isin(aminoAcidResNames)].RES_ID.tolist()
        chainLigandResidues = chainDf[~chainDf["RES_NAME"].isin(aminoAcidResNames)].RES_ID.tolist()
        if len(chainProteinResidues) > 0 and len(chainLigandResidues) > 0:
            isSharedChains = True
            sharedChains[chainId] = "Chain has both protein and ligand residues"
   
    return isSharedChains, sharedChains or None



#################################################################################################
def check_for_organometallic_ligand(pdbDf: pd.DataFrame, uaaInfo: Optional[Dict] = None) -> tuple[bool, Optional[Dict[str, int]]]:
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
    organoMetallicResidues: Dict = {}

    # Loop through chains and residues in the pdb dataframe
    for chainId, chainDf in pdbDf.groupby(f"CHAIN_ID"):
        for resId, resDf in chainDf.groupby(f"RES_ID"):
            resName: str = resDf["RES_NAME"].iloc[0]

            # Skip if amino acid residue
            if resName in aminoAcidResNames:
                continue
            try: 
                resElements: Set[str] = set(resDf["ELEMENT"])
            except:
                logging.info(f"No elements column found in pdb dataframe")
                return False,  None
    
            # If there are non-organic atoms, add the residue ID and atom count to the dictionary
            inorganicElements = [ele for ele in resElements if ele.upper() not in organicElements]


            
            if len(inorganicElements) > 0:
                organoMetallicResidues[f"{chainId}:{resName}:{str(resId)}"] = inorganicElements

    # Return boolean indicating if organometallic ligand were found and the dictionary
    return bool(organoMetallicResidues), organoMetallicResidues or None


    
#################################################################################################
def check_for_non_canonical_amino_acids(pdbDf: pd.DataFrame, inputDir: DirectoryPath, config: dict) -> tuple[bool, Optional[Dict[str, int]]]:
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

    # Dictionary to store residue IDs and messages
    nonCanonicalAminoAcids: Dict = {}
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
            # Look for missing frcmod and/or lib files
            resKey: str = f"{chainId}:{resName}:{str(resId)}"
            frcmodFile: FilePath = p.join(inputDir, f"{resName}.frcmod")
            libFile: FilePath = p.join(inputDir, f"{resName}.lib")

            if not p.isfile(frcmodFile) and not p.isfile(libFile):
                nonCanonicalAminoAcids[resKey] = "non-canonical amino acid missing frcmod and lib files"
            elif not p.isfile(frcmodFile):
                nonCanonicalAminoAcids[resKey] = "non-canonical amino acid missing frcmod file"
            elif not p.isfile(libFile):
                nonCanonicalAminoAcids[resKey] = "non-canonical amino acid missing lib file"

    # Return boolean indicating if non-canonical amino acids were found and the dictionary
    return bool(nonCanonicalAminoAcids), nonCanonicalAminoAcids or None
                
#################################################################################################
def check_for_missing_sidechains(pdbDf: pd.DataFrame) -> tuple[bool, Optional[Dict[str, List[str]]]]:
    """
    Checks if a pdb dataframe contains missing sidechains.
    
    Args:
        pdbDf (pd.DataFrame): The pdb dataframe.
        
    Returns:
        Tuple[bool, Optional[Dict[str, List[str]]]]: A tuple containing a boolean indicating if sidechains were found,
        and a dictionary with the residue IDs of the residues with missing sidechains and a list of the missing atoms.
    """
    ## initialise a set of backbone atom names and terminal oxygen name
    backboneAtoms: set = {"N", "C", "O", "CA", "OXT"}
    ## get amino acid names and dictionary containing heavy atom counts for each residue
    aminoAcidResNames = drListInitiator.get_amino_acid_residue_names()
    heavySideChainAtomCounts = drListInitiator.get_residue_heavy_atom_counts()
    ## get only the protein part of the pdb dataframe
    protDf = pdbDf[pdbDf["RES_NAME"].isin(aminoAcidResNames)]
    ## initialise an empty dict to store missing sidechains
    missingSidechains: Dict = {}
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
            ## ecxlude backbone atoms
            sideChainAtoms: List[str] = list(set([atom for atom in heavyAtomNames if atom not in backboneAtoms ]))
            ## check if number of heavy sidechain atoms matches expected value
            if  len(sideChainAtoms) != correctHeavyAtomCount:
                missingSidechains[f"{chainId}:{resName}:{str(resId)}"] = sideChainAtoms
    ## return boolean indicating if sidechains with missing atoms were found and the dictionary
    return bool(missingSidechains), missingSidechains or None

#################################################################################################
def check_for_broken_chains(pdbDf: pd.DataFrame) -> tuple[bool, Optional[Dict[str, List[str]]]]:
    """
    Check for broken chains in the protein dataframe.
    
    Args:
        pdbDf (pd.DataFrame): The protein dataframe.

    Returns:
        Tuple[bool, Optional[Dict[str, List[str]]]: A tuple containing a boolean indicating if broken chains were found,
        and a dictionary with chain IDs as keys and a list of broken residues or a specific message as values.
    """
    aminoAcidResNames = drListInitiator.get_amino_acid_residue_names()

    # Filter only the protein part of the dataframe
    protDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"].isin(aminoAcidResNames)]
    ## initialise an empty dict to store broken chains
    brokenChains: Dict = {}
    ## loop through chains
    for chainId, chainDf in protDf.groupby(f"CHAIN_ID"):
        ## get residue IDs
        resIds = chainDf["RES_ID"].unique().tolist()
        
        # Look for non-consecutive residue numbering
        isConsecutive, nonConsecutiveResidues = are_consecutive(resIds)
        ## update dict if needed
        if not isConsecutive:
            brokenChains[chainId] = nonConsecutiveResidues
        
        # Look for termini in the middle of chains

        resIdsWithOxt = chainDf[chainDf["ATOM_NAME"] == "OXT"]["RES_ID"].tolist()
        lastResidueId = resIds[-1]
        try:
            resIdsWithOxt.remove(lastResidueId)
        except:
            pass
        if len(resIdsWithOxt) == 0:
            continue
        
        brokenChains[chainId] = "OXT atom name found in residues " + ", ".join(map(str, resIdsWithOxt))

    return bool(brokenChains), brokenChains or None

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
def check_for_multiple_conformers(pdbDf: pd.DataFrame) -> tuple[bool, Optional[Dict[str, List[str]]]]:
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
    
    # Dictionary to store residue IDs and duplicated atoms
    multipleConformers: Dict = {}
    
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
                multipleConformers[f"{chainId}:{resName}:{str(resId)}"] = list(set(duplicatedAtoms))

    # Return boolean indicating if multiple conformers were found and the dictionary
    return bool(multipleConformers), multipleConformers or None





#################################################################################
if __name__ == "__main__":
    ### for testing
    pdbDir = "/home/esp/scriptDevelopment/drMD/01_inputs/"
    dummyConfig = {"pathInfo": {"outputDir": "/home/esp/scriptDevelopment/drMD/03_outputs"}}
    pdb_triage(pdbDir, dummyConfig)

#################################################################################

