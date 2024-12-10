## BASIC PYTHON LIBRARIES
import os
from os import path as p
import subprocess
from subprocess import run
import string
from shutil import copy
import logging
import pandas as pd
import numpy as np

## drMD MODULES
from UtilitiesCloset import drFixer, drSplash, drListInitiator
from ExaminationRoom import drLogger

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

## CLEAN CODE
from typing import Optional, Dict, List, Tuple, Union, Any
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

#####################################################################################

def prep_protocol(config: dict) -> Tuple[str, str, str]:
    """
    Prepares the protocol for the simulation.

    Args:
        config (dict): The configuration dictionary containing the necessary information.

    Returns:
        Tuple[str, str, str]: A tuple containing the merged PDB file, input coordinates file, and Amber parameter files.

    This function prepares the protocol for the simulation by performing the following steps:
    1. Retrieves the output directory from the configuration dictionary.
    2. Creates the preparation directory if it doesn't exist.
    3. Checks if the preparation is already complete by checking for the presence of the necessary files.
    4. If the preparation is already complete, returns the relevant files.
    5. If the preparation is not complete, proceeds with the main preparation protocol.
    6. Splits the input PDB file into protein and ligand files.
    7. Prepares the ligand parameters and outputs the ligand PDB files.
    8. Prepares the protein structure.
    9. Combines the protein and ligand PDB files.
    10. Merges the protein PDB files.
    11. Makes Amber parameter files using TLEAP.
    12. Returns the merged PDB file, input coordinates file, and Amber parameter files.

    If the configuration dictionary does not contain ligand information, the function proceeds with the following steps:
    1. Prepares the protein structure.
    2. Merges the protein PDB files.
    3. Makes Amber parameter files using TLEAP.
    4. Returns the merged PDB file, input coordinates file, and Amber parameter files.

    Note:
        The function assumes that the necessary directories and files are present in the configuration dictionary.
    """

    ## read config for path info
    outDir: DirectoryPath = config["pathInfo"]["outputDir"]
    protName: str = config["proteinInfo"]["proteinName"]
    ## create a prep dir if it doesn't exist
    prepDir: DirectoryPath = p.join(outDir,"00_prep")
    os.makedirs(prepDir,exist_ok=True)

    set_up_logging(outDir, protName)


    skipPrep, prepFiles = choose_to_skip_prep(config=config, prepDir=prepDir, protName=protName)
    if skipPrep:
        drLogger.log_info(f"Prep steps already complete for {protName}: Skipping ...",True)
        return prepFiles


    ######### MAIN PREP PROTOCOL #########
    if "ligandInfo" in config:
        solvatedPdb, inputCoords, amberParams = ligand_prep_protocol(config=config,
                                                                      protName=protName,
                                                                        prepDir=prepDir)
        

    
    else:
        solvatedPdb, inputCoords, amberParams = no_ligand_prep_protocol(config=config,
                                                                          protName=protName,
                                                                            prepDir=prepDir)
        

    
    drLogger.log_info(f"Prep steps complete for {protName}!",True)
    drLogger.close_logging()

    return solvatedPdb, inputCoords, amberParams
#####################################################################################
def no_ligand_prep_protocol(config: dict, protName: str, prepDir: DirectoryPath) -> Tuple[FilePath, FilePath, FilePath]:
    ## PREPARE PROTEIN STRUCTURE
    protPdb: FilePath = prepare_protein_structure(config=config, outDir = prepDir)  
    ## MERGE PROTEIN PDBS
    outName: str = config["pathInfo"]["outputName"]
    ## MAKE AMBER PARAMETER FILES WITH TLEAP
    drLogger.log_info(f"Solvating, Charge Balencing and Creating parameters for {protName}...")

    inputCoords, amberParams, solvatedPdb = make_amber_params(outDir = p.join(prepDir,"PROT"),
                                                    pdbFile= protPdb,
                                                outName= outName, 
                                                config = config)
    
    solvatedPdb = drFixer.reset_chains_residues(protPdb, solvatedPdb)

    return solvatedPdb, inputCoords, amberParams

#####################################################################################
def ligand_prep_protocol(config: dict, protName: str, prepDir: DirectoryPath) -> Tuple[str, str, str, str]:
        ## SPLIT INPUT PDB INTO PROT AND ONE FILE PER LIGAND
        inputPdb: FilePath = config["pathInfo"]["inputPdb"]
        split_input_pdb(inputPdb =inputPdb,
                        config = config,
                        outDir=prepDir)
        ## PREPARE LIGAND PARAMETERS, OUTPUT LIGAND PDBS
        ligandPdbs, ligandFileDict = prepare_ligand_parameters(config = config)
        ## PREPARE PROTEIN STRUCTURE
        protPdb = prepare_protein_structure(config=config, outDir = prepDir)
        ## RE-COMBINE PROTEIN AND LIGAND PDB FILES
        wholePrepDir: DirectoryPath = p.join(prepDir,"WHOLE")
        os.makedirs(wholePrepDir,exist_ok=True)
        allPdbs: List[FilePath] = [protPdb] + ligandPdbs
        if p.isdir(p.join(prepDir, "IONS")):
            allPdbs = allPdbs + [p.join(prepDir, "IONS", f"IONS.pdb")]
        outName: str = config["pathInfo"]["outputName"]
        mergedPdb: FilePath = p.join(wholePrepDir,f"{protName}.pdb")
        pdbUtils.mergePdbs(pdbList=allPdbs, outFile = mergedPdb)
        ## RESET ATOM NUMBERS AFTER MERGING
        mergedPdb: FilePath = drFixer.reset_atom_numbers(pdbFile = mergedPdb)
        ## MAKE AMBER PARAMETER FILES WITH TLEAP
        inputCoords, amberParams, solvatedPdb  = make_amber_params(outDir = wholePrepDir,
                            ligandFileDict=ligandFileDict,
                            pdbFile= mergedPdb,
                            config= config,
                            outName= outName)
        niceChainsPdb = drFixer.reset_chains_residues(mergedPdb, solvatedPdb)
        return niceChainsPdb, inputCoords, amberParams
#####################################################################################   
def choose_to_skip_prep(config: dict, prepDir: DirectoryPath, protName: str) -> Tuple[bool, Optional[Tuple[FilePath, FilePath, FilePath]]]:
    """
    Check if the preparation is already complete and return the relevant files if it is.

    Args:
        config (dict): The configuration dictionary containing the necessary information.
        prepDir (str): The path to the preparation directory.
        protName (str): The name of the protein.

    Returns:
        Tuple[FilePath, FilePath, FilePath]: A tuple containing the paths to the merged PDB file, input coordinates file, and Amber parameter files.
    """
    ## init some false booleans for check later
    amberParams = False
    inputCoords = False
    pdbFile = False
    ## check if prmtop, inpcrd, and pdb files exist in prep dir

    if "ligandInfo" in config:
        completePrepDir = p.join(prepDir,"WHOLE")
    else:
        completePrepDir = p.join(prepDir,"PROT")

    if p.isdir(p.join(completePrepDir)):
        for file in os.listdir(completePrepDir):
            if p.splitext(file)[1] == ".prmtop":
                amberParams: FilePath = p.join(completePrepDir,file)
            elif p.splitext(file)[1] == ".inpcrd":
                inputCoords: FilePath = p.join(completePrepDir,file)
            if p.splitext(file)[1] == ".pdb" and not file == f"{protName}.pdb":
                pdbFile: FilePath = p.join(completePrepDir, file)

    ## return to drOperator if files already have been made
    if p.isfile(amberParams) and p.isfile(inputCoords) and p.isfile(pdbFile):
        return True, (pdbFile, inputCoords, amberParams)
    else:
        return False, None
#####################################################################################
def set_up_logging(outDir, protName):
    ## set up logging
    logDir = p.join(p.dirname(outDir), "00_drMD_logs")
    os.makedirs(logDir, exist_ok=True)
    prepLog: FilePath = p.join(logDir,f"{protName}_prep.log")
    drLogger.setup_logging(prepLog)
    drLogger.log_info(f"Running Prep protocol for {protName}...",True)
#####################################################################################
def find_ligand_charge(ligDf: pd.DataFrame,
                        ligName: str,
                          outDir: DirectoryPath,
                            pH: float) -> int:
    """
    This function uses propka to identify charges on a ligand.
    
    Args:
        ligDf (pd.DataFrame): Ligand DataFrame.
        ligName (str): Ligand name.
        outDir (str): Output directory.
        pH (float): pH value.
        
    Returns:
        int: Total charge of the ligand.
    """
    logDir = p.join(p.dirname(outDir), "00_drMD_logs")
    drLogger.setup_logging(p.join(logDir, "02_ligand_pka_predictions.log"))
    drLogger.log_info(f"Running propka to predict ligand charges for {ligName} at pH {str(pH)}...", True)
    # Change working directory to the output directory
    ## propka3 needs this
    os.chdir(outDir)
    
    # Fix atom names in the ligand dataframe
    ligDf: pd.DataFrame = drFixer.fix_atom_names(ligDf)
    # Remove hydrogen atoms and rename ATOM to HETATM
    ligDf: pd.DataFrame = ligDf[ligDf["ELEMENT"] != "H"]
    ligDf["ATOM"] = "HETATM"
    # Create a temporary pdb file from the ligand dataframe
    tmpPdb: FilePath = p.join(outDir, f"{ligName}.pdb")
    pdbUtils.df2pdb(ligDf, tmpPdb, chain=False)
    
    # Run propka to predict charges on the ligand
    proPkaCommand: str = f"propka3 {tmpPdb}"

    run_with_log(proPkaCommand, "ligand PKA prediction with propka", None)
    
    # Read the propka output to extract charge at the specified pH
    proPkaFile: FilePath = f"{ligName}.pka"
    with open(proPkaFile, "r") as f:
        pkaPredictions: List[List[str]] = []
        extract = False
        for line in f:
            # Extract predictions from the propka output
            if line.startswith("SUMMARY OF THIS PREDICTION"):
                extract = True
            if extract and line.startswith("-"):
                break
            if extract and not line == "\n":
                pkaPredictions.append(line.split())
        pkaPredictions = pkaPredictions[2:]
        
        # Calculate total charge
        totalCharge = 0
        for pred in pkaPredictions:
            atomType: str = pred[-1][0]
            atomProb: float = float(pred[-2])
            if atomType == "O":
                if atomProb < pH:
                    totalCharge -= 1
            elif atomType == "N":
                if atomProb > pH:
                    totalCharge += 1

    # Clean up temporary files
    # [os.remove(p.join(outDir, file)) for file in os.listdir(outDir) if not p.splitext(file)[1] == ".yaml"]
    drLogger.close_logging()
    return totalCharge

#####################################################################################
def split_input_pdb(inputPdb: FilePath, config: Dict, outDir: DirectoryPath) -> None:
    """
    Split an input PDB file into separate PDB files for each ligand and a protein-only PDB file.

    Args:
        inputPdb (str): The path to the input PDB file.
        config (dict): The configuration dictionary containing information about the ligand.
        outDir (str): The output directory where the split PDB files will be saved.

    Returns:
        None
    """
    # Read whole pdb into a df
    pdbDf: pd.DataFrame = pdbUtils.pdb2df(inputPdb)
    # Write each ligand to a separate pdb file
    ligandInfo: List[dict] = config["ligandInfo"]
    for ligand in ligandInfo:
        ligandName: str = ligand["ligandName"]
        ligPrepDir: DirectoryPath = p.join(outDir, ligandName)
        os.makedirs(ligPrepDir, exist_ok=True)
        ligDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"] == ligandName]
        pdbUtils.df2pdb(ligDf, p.join(ligPrepDir, f"{ligandName}.pdb"), chain=False)
        pdbDf.drop(pdbDf[pdbDf["RES_NAME"] == ligandName].index, inplace=True)

    ## extract ions from pdb dataframe
    ionResidueNames = drListInitiator.get_ion_residue_names()
    ionDf = pdbDf[pdbDf["RES_NAME"].isin(ionResidueNames)]
    if len(ionDf) > 0:
        # Write ions to pdb file    
        ionDir = p.join(outDir, "IONS")
        os.makedirs(ionDir, exist_ok=True)
        pdbUtils.df2pdb(ionDf, p.join(ionDir, "IONS.pdb"))


    # Write protein only to pdb file
    protDf = pdbDf[~pdbDf["RES_NAME"].isin(ionResidueNames)]
    protPrepDir: DirectoryPath = p.join(outDir, "PROT")
    os.makedirs(protPrepDir, exist_ok=True)
    pdbUtils.df2pdb(protDf, p.join(protPrepDir, "PROT.pdb"))

#############################  PROTONATION & PDB CREATION ###############################


def ligand_protonation(
    ligand: Dict[str, Union[bool, int, str]],
    ligPrepDir: DirectoryPath,
    ligandName: str,
    ligandPdbs: List[FilePath],
    ligPdb: FilePath) -> Tuple[str, List[str]]:
    """
    Protonates a ligand if specified in the configuration, otherwise performs protonation using Open Babel and pdb4amber.

    Args:
        ligand (Dict[str, Union[bool, int, str]]): A dictionary containing information about the ligand.
                                                   The dictionary should contain the following keys:
                                                   - "protons": Whether to protonate the ligand (bool).
                                                   - "charge": The charge of the ligand (int).
                                                   - "ligandName": The name of the ligand (str).
        ligPrepDir (str): The directory where the ligand files are prepared (str).
        ligandName (str): The name of the ligand (str).
        ligandPdbs (List[str]): A list of ligand pdb files (List[str]).
        liPdb (FilePath): The path to the ligand pdb file (FilePath).
    Returns:
        Tuple[str, List[str]]: A tuple containing the protonated ligand pdb file (str) and the updated list of ligand pdb files (List[str]).
    """

    # If the ligand is already protonated, return the ligand pdb file
    if ligand["protons"]:
        ligandPdbs.append(ligPdb)
        return ligPdb, ligandPdbs

    # If the ligand is not protonated, perform protonation using Open Babel and pdb4amber
    else:
        # # find pdb ligand pdb file
        ligPdb_H: FilePath = p.join(ligPrepDir, f"{ligandName}_H.pdb")

        # Protonate the ligand using Open Babel
        obabelCommand: str = f"obabel {ligPdb} -O {ligPdb_H} -h"
        run_with_log(obabelCommand, "Ligand protonation with openbabel", ligPdb_H)

        ligandPdbs.append(ligPdb_H)
        return ligPdb_H, ligandPdbs

###############################  MOL2 CREATION #####################################
def ligand_mol2(
    ligand: Dict[str, Any],
    inputDir: DirectoryPath,
    ligandName: str,
    ligandParamDir: DirectoryPath,
    ligandPrepDir: DirectoryPath,
    ligPdb: FilePath,
    ligandFileDict: Dict[str, str]) -> Tuple[str, Dict[str, str]]:
    """
    Create a mol2 file for the ligand.

    This function looks for a mol2 file from the configuration, then in the ligParamDir.
    If not found, it creates a new mol2 file using antechamber.

    Args:
        ligand (Dict[str, Any]): Ligand information dictionary.
        inputDir (DirectoryPath): Path to the input directory.
        ligandName (str): Name of the ligand.
        ligandParamDir (DirectoryPath): Path to the ligand parameter directory.
        ligandPrepDir (DirectoryPath): Path to the ligand preparation directory.
        ligPdb (FilePath): Path to the ligand pdb file.
        ligandFileDict (Dict[str, str]): Dictionary containing ligand file information.
        prepLog (FilePath): Path to the preparation log file.

    Returns:
        Tuple[str, Dict[str, str]]: A tuple containing the path to the mol2 file and the updated lig_file_dict.
    """

    ligLib = False
    # Look for mol2 from config, then in ligParamDir, if not found, create new mol2
    if ligand.get("mol2"):  # Look in config
        ligMol2: FilePath = p.join(inputDir, f"{ligandName}.mol2")
        ligLib: FilePath = p.join(inputDir, f"{ligandName}.lib")

    elif p.isfile(p.join(ligandParamDir, f"{ligandName}.mol2")):  # Look in ligParamDir
        ligMol2: FilePath = p.join(ligandParamDir, f"{ligandName}.mol2")
    else:  # Convert to mol2 with antechamber
        charge: int = int(ligand["charge"])
        ligMol2: FilePath = p.join(ligandPrepDir, f"{ligandName}.mol2")
        # Create antechamber command
        antechamberCommand: str = (
            f"antechamber -i {ligPdb} -fi pdb -o {ligMol2} -fo mol2 -c bcc -s 2 -nc {charge}"
        )
        # Run antechamber and log the command
        run_with_log(antechamberCommand, "Ligand charge calculation with antechamber", ligMol2)
        # Copy to ligParamDir for future use
        copy(ligMol2, p.join(ligandParamDir, f"{ligandName}.mol2"))
    
    # Add mol2 path to lig_file_dict
    ligandFileDict.update({"mol2": ligMol2})
    if  ligLib:
        ligandFileDict.update({"lib": ligLib})

    return ligMol2, ligandFileDict
######################### TOPPAR CREATION ##########################################
def ligand_toppar(ligand: dict,
                   inputDir: DirectoryPath,
                     ligandName: str,
                       ligParamDir: DirectoryPath,
                         ligPrepDir: DirectoryPath,
                           ligMol2: FilePath,
                             ligFileDict: dict) -> dict:
    """
    Create a frcmod file for a given ligand.

    Args:
        ligand (dict): A dictionary containing information about the ligand.
        inputDir (str): The path to the input directory.
        ligandName (str): The name of the ligand.
        ligParamDir (str): The path to the directory where the ligand parameter files are located.
        ligPrepDir (str): The path to the directory where the ligand preparation files are located.
        ligMol2 (str): The path to the ligand mol2 file.
        ligFileDict (dict): A dictionary containing information about the ligand files.
        prepLog (str): The log file for preparation.

    Returns:
        dict: A dictionary containing information about the ligand files, including the path to the frcmod file.
    """

    # Look for frcmod from config, then in ligParamDir, if not found, create new frcmod
    if ligand.get("toppar"):  # Look in config
        # Use the frcmod file from the config
        ligFrcmod: FilePath = p.join(inputDir, f"{ligandName}.frcmod")

    elif p.isfile(p.join(ligParamDir, f"{ligandName}.frcmod")):  # Look in ligParamDir
        # Use the existing frcmod file in ligParamDir
        ligFrcmod: FilePath = p.join(ligParamDir, f"{ligandName}.frcmod")

    else:  # Create new frcmod using parmchk2
        # Create a new frcmod file using parmchk2
        ligFrcmod: FilePath = p.join(ligPrepDir, f"{ligandName}.frcmod")
        parmchk2Command: str = f"parmchk2 -i {ligMol2} -f mol2 -o {ligFrcmod}"
        run_with_log(parmchk2Command, "Ligand parameter generation with parmchk2", ligFrcmod)
        copy(ligFrcmod, p.join(ligParamDir, f"{ligandName}.frcmod"))

    # Add frcmod path to ligFileDict
    ligFileDict.update({"frcmod": ligFrcmod})

    return ligFileDict
#####################################################################################
def prepare_ligand_parameters(config: Dict) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """
    Prepares the ligand parameters for a given configuration file.

    Args:
        config (dict): The configuration file containing all information needed

    Returns:
        Tuple[List[str], Dict[str, Dict[str, str]]]: A tuple containing two elements:
            - ligandPdbs (List[str]): A list of ligand PDB files.
            - ligandFileDict (Dict[str, Dict[str, str]]): A dictionary containing the ligand file information.
    """
    # read inputs from config file
    outDir: DirectoryPath = config["pathInfo"]["outputDir"]
    ligandInfo: dict = config["ligandInfo"]
    inputDir: DirectoryPath = config["pathInfo"]["inputDir"]
    mainDir: DirectoryPath = p.dirname(config["pathInfo"]["outputDir"])
    # create a dir to save parameter files in (saves re-running on subsequent runs)
    ligParamDir: DirectoryPath = p.join(mainDir,"01_ligand_parameters")
    os.makedirs(ligParamDir,exist_ok=True)
    # initialise list to store pdb files and dict to store all infod
    ligandPdbs: List[FilePath] = []
    ligandFileDict: Dict = {}
    # for each ligand in config
    for ligand in ligandInfo:
        ## init an empty dict to store the ligand file paths
        ligFileDict: Dict = {}
        ## get ligand name
        ligandName: str = ligand["ligandName"]
        ## write to log
        drLogger.log_info(f"Preparing ligand {ligandName}...",True)
        # find files and directories
        ligPrepDir: DirectoryPath = p.join(outDir,"00_prep",ligandName)
        os.chdir(ligPrepDir)
        ## get ligand pdb location
        ligPdb: FilePath = p.join(ligPrepDir, f"{ligandName}.pdb")

        # Protonate the ligand
        drLogger.log_info(f"Protonating ligand {ligandName}...",True)
        ligPdb, ligandPdbs = ligand_protonation(ligand,ligPrepDir,ligandName,ligandPdbs, ligPdb)  

        # deal with atom names in ligand, make sure they are all unique
        ligPdb = ensure_ligand_atoms_are_unique(ligPdb)


        # Create mol2 file
        drLogger.log_info(f"Calculating partial charges for ligand {ligandName}...",True)
        ligMol2, ligFileDict = ligand_mol2(ligand,inputDir,ligandName,ligParamDir,
                                          ligPrepDir,ligPdb,ligFileDict)
        
        # Create frcmod file
        drLogger.log_info(f"Creating parameter files for ligand {ligandName}...",True)
        ligFileDict = ligand_toppar(ligand,inputDir,ligandName,ligParamDir,
                                    ligPrepDir,ligMol2,ligFileDict)

        ligandFileDict.update({ligandName:ligFileDict})
    return ligandPdbs, ligandFileDict

#####################################################################################
def ensure_ligand_atoms_are_unique(ligPdb: FilePath) -> FilePath:
    """
    Ensure that all atoms in the ligand have unique atom names.

    Args:
        ligPdb (str): Path to the input ligand PDB file.

    Returns:
        str: Path to the output PDB file with unique atom names.
    """
    # Read PDB file
    ligDf = pdbUtils.pdb2df(ligPdb)

    atomNames = ligDf["ATOM_NAME"].tolist()
    uniqueAtomName = set(atomNames)

    if len(atomNames) != len(uniqueAtomName):
        ligDf = rename_heteroatoms(ligDf)    

    pdbUtils.df2pdb(ligDf, ligPdb)
    return ligPdb
#####################################################################################
def rename_heteroatoms(pdbDf: pd.DataFrame) -> pd.DataFrame:
    # Generate new names for hydrogens
    letters: List[str] = list(string.ascii_uppercase)
    numbers: List[str] = [str(i) for i in range(1,10)]
    newNameSuffixes: List[str] = [number+letter for letter in letters for number in numbers]

    # Rename hydrogens in the DataFrame
    count: int = 0
    for index, _ in pdbDf.iterrows():
        ## Use element column if available
        if len(pdbDf.loc[index,"ELEMENT"]) > 0:
            pdbDf.loc[index,"ATOM_NAME"] = pdbDf.loc[index,"ELEMENT"]  + newNameSuffixes[count]
        ## Use first letter of atom name if element is not available (will break for 2 letter element names)
        else:
             pdbDf.loc[index,"ATOM_NAME"] = pdbDf.loc[index,"ATOM_NAME"][0]  + newNameSuffixes[count]
        count += 1

    return pdbDf

#####################################################################################
def rename_hydrogens(pdbFile: FilePath, outFile: DirectoryPath) -> None:
    """
    Renames all hydrogens in a PDB file to new names.

    Args:
        pdbFile (str): Path to the input PDB file.
        outFile (str): Path to the output PDB file.

    Returns:
        None
    """
    # Read PDB file and get the DataFrame of hydrogens
    pdbDf: pd.DataFrame = pdbUtils.pdb2df(pdbFile)
    hDf: pd.DataFrame = pdbDf[pdbDf["ELEMENT"]=="H"]

    # Generate new names for hydrogens
    letters: List[str] = list(string.ascii_uppercase)
    numbers: List[str] = [str(i) for i in range(1,10)]
    newNameHs: List[str] = ["H"+letter+number for letter in letters for number in numbers]

    # Rename hydrogens in the DataFrame
    count: int = 0
    for index, _ in hDf.iterrows():
        pdbDf.loc[index,"ATOM_NAME"] = newNameHs[count]
        count += 1

    # Write the modified DataFrame back to the PDB file
    pdbUtils.df2pdb(pdbDf, outFile, chain=False)
#####################################################################################
def prepare_protein_structure(config: Dict, outDir: DirectoryPath) -> FilePath:
    """
    Prepare the protein structure for simulations.

    This function prepares the protein structure by copying the input PDB file
    to the output directory if it doesn't exist.

    Args:
        config (Dict): The configuration dictionary.
        outDir (DirectoryPath): The output directory.
        prepLog (FilePath): The path to the preparation log file.

    Returns:
       Union[os.PathLike, str]: The PATH to a PDB file.
    """
    # Find files and directories
    protPrepDir: DirectoryPath = p.join(outDir, "PROT")  # Directory to prepare the protein
    os.makedirs(protPrepDir, exist_ok=True)
    os.chdir(protPrepDir)
    protPdb: FilePath = p.join(protPrepDir, "PROT.pdb")  # Path to the protein PDB file
    # Check for PROT.pdb in protPrepDir (won't be there if noLigand)
    if not p.isfile(protPdb):
        # Copy the input PDB file to the output directory
        copy(config["pathInfo"]["inputPdb"], protPdb)

    ## read per-protein config file and return protein PDB if already protonated
    proteinInfo: Dict = config.get("proteinInfo")
    isProteinProtonated: bool = proteinInfo.get("protons")
    if isProteinProtonated:
        newHisPdb = sort_out_histidine_names(protPdb)
        return newHisPdb

    ## use pdb2pqr to protonate the protein at a specific pH
    pH: int = str(float(config["miscInfo"]["pH"]))
    protPqr = p.join(protPrepDir, "PROT.pqr")
    pdb2pqrCommand: str = ["pdb2pqr",
                         "--ffout", "AMBER",
                         "--titration-state-method", "propka",
                          "--keep-chain",
                           "--with-ph", str(float(pH)),
                             protPdb, protPqr]
    
    run_with_log(pdb2pqrCommand, "Protein protonation with pdb2pqr and propka", protPqr)
    ## fix the betafactor and occupancy cols to zero
    protPqr = pdbUtils.pdb2df(protPqr)
    protPqr["OCCUPANCY"] = 0
    protPqr["BETAFACTOR"] = 0
    ## add the element column using the first letter of the ATOM_NAME 
    # (breaks for two-letter element symbols like "Cl", should be ok for proteins)
    protPqr["ELEMENT"] = protPqr["ATOM_NAME"].str[0]
    ## write to pdb file
    protonatedPdb: FilePath = p.join(protPrepDir, "PROT_protonated.pdb")

    pdbUtils.df2pdb(protPqr, protonatedPdb)

    return protonatedPdb

#################################################################################
def sort_out_histidine_names(protPdb: str) -> str:
    """
    Rename histidine residues in a PDB file, according to presence of hydrogen names.
    """
    pdbDf = pdbUtils.pdb2df(protPdb)
    histidineNames = {"HIS", "HID", "HIE", "HIP"}

    # Filter for histidine residues
    histDf = pdbDf[pdbDf["RES_NAME"].isin(histidineNames)]

    for chainId, chainDf in histDf.groupby("CHAIN_ID"):
        for resId, resDf in chainDf.groupby("RES_ID"):
            resAtomNames = resDf["ATOM_NAME"].tolist()
            if "HD1" in resAtomNames and "HE2" in resAtomNames:
                # Rename to HIP
                pdbDf.loc[(pdbDf["CHAIN_ID"] == chainId) & (pdbDf["RES_ID"] == resId), "RES_NAME"] = "HIP"
            elif "HD1" in resAtomNames:
                # Rename to HID
                pdbDf.loc[(pdbDf["CHAIN_ID"] == chainId) & (pdbDf["RES_ID"] == resId), "RES_NAME"] = "HID"
            elif "HE2" in resAtomNames:
                # Rename to HIE
                pdbDf.loc[(pdbDf["CHAIN_ID"] == chainId) & (pdbDf["RES_ID"] == resId), "RES_NAME"] = "HIE"

    # Save the modified DataFrame back to a PDB file
    outputPdb = protPdb.replace(".pdb", "_newHis.pdb")
    pdbUtils.df2pdb(pdbDf, outputPdb)

    return outputPdb
#################################################################################
def detect_disulphides(pdbFile: FilePath) -> List[str]:
    """
    Detect disulphides in a PDB file.
    """
    pdbDf: pd.DataFrame = pdbUtils.pdb2df(pdbFile)
    cyxDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"] == "CYX"]

    sgDf: pd.DataFrame = cyxDf[cyxDf["ATOM_NAME"] == "SG"]


    sgCoords: np.array = sgDf[["X", "Y", "Z"]].values


    distanceMatrixDf = pd.DataFrame(np.sqrt(np.sum((sgCoords[:, np.newaxis, :] - sgCoords[np.newaxis, :, :])**2, axis=-1)),
                                   index=sgDf["RES_ID"], columns=sgDf["RES_ID"])
    
    disulphidesDf = distanceMatrixDf[(distanceMatrixDf < 2.5) & (distanceMatrixDf != 0)]


    disulphidePairs = [sorted((index, col)) for index, row in disulphidesDf.iterrows() for col in disulphidesDf.columns if not pd.isna(row[col])]

    uniqueDisulphidePairs = {tuple(pair) for pair in disulphidePairs}

    return uniqueDisulphidePairs

    






#####################################################################################
def make_amber_params(
    outDir: DirectoryPath,
    pdbFile: FilePath,
    outName: str,
    config : Dict,
    ligandFileDict: Optional[Dict[str, Dict[str, str]]] = None
) -> Tuple[FilePath, FilePath, FilePath]:
    """
    Prepare the protein structure for simulations using Amber.

    Args:
        outDir (DirectoryPath): The output directory.
        pdbFile (FilePath): The path to the input PDB file.
        outName (str): The name of the prepared protein.
        prepLog (FilePath): The path to the preparation log file.
        ligandFileDict (Optional[Dict[str, Dict[str, str]]]): A dictionary containing the paths to
            the mol2 and frcmod files for each ligand. Defaults to None.

    Returns:
        Tuple[FilePath, FilePath, FilePath]: The paths to the input coordinate file and the Amber
            parameter files and a pdb file
    """

    drLogger.log_info(f"Preparing parameters for system: {outName}...",True)
    # Change the working directory to the output directory
    os.chdir(outDir)


    ## find dusulphides
    disulphidePairs: List[Tuple[int, int]] = detect_disulphides(pdbFile)

    boxGeometry: str = config["miscInfo"]["boxGeometry"]
    if boxGeometry == "cubic":
        solvateKeyword: str = "solvatebox"
    elif boxGeometry == "octahedral":
        solvateKeyword: str = "solvateoct"

    # Write the TLEAP input file
    tleapInput: str = p.join(outDir, "TLEAP.in")
    with open(tleapInput, "w") as f:
        # Load Amber force fields and TIP3P water model
        f.write("source leaprc.protein.ff19SB\n")
        f.write("source leaprc.gaff2\n")
        f.write("source leaprc.water.tip3p\n\n")

        # Load Amber force fields for ions
        f.write("loadamberparams frcmod.ions1lm_126_tip3p\n")
        f.write("loadamberparams frcmod.ions234lm_126_tip3p\n")

        # Load molecular mechanics parameters for ligand
        if ligandFileDict:
            for ligandName, ligandInfo in ligandFileDict.items():
                ligMol2: FilePath = ligandInfo["mol2"]
                ligFrcmod: FilePath = ligandInfo["frcmod"]
                ligLib: FilePath = ligandInfo.get("lib",False)
                if p.isfile(ligMol2):
                    f.write(f"{ligandName} = loadmol2 {ligMol2}\n")
                if p.isfile(ligFrcmod):
                    f.write(f"loadamberparams {ligFrcmod}\n")
                if p.isfile(ligLib):
                    f.write(f"loadoff {ligLib}\n")


        # Load the protein structure
        f.write(f"mol = loadpdb {pdbFile}\n")

        # Solvate the protein and add ions

        f.write(f"{solvateKeyword} mol TIP3PBOX 10.0\n")
        f.write("addions mol Na+ 0\n")
        f.write("addions mol Cl- 0\n")

        ## make disulphide bonds
        for disulphidePair in disulphidePairs:
            f.write(f"bond mol.{disulphidePair[0]}.SG mol.{disulphidePair[1]}.SG\n")

        # Save the solvated protein structure and parameter files
        solvatedPdb: str = f"{outName}_solvated.pdb"
        f.write(f"savepdb mol {solvatedPdb}\n")
        prmTop: FilePath = p.join(outDir, f"{outName}.prmtop")
        inputCoords: FilePath = p.join(outDir, f"{outName}.inpcrd")
        f.write(f"saveamberparm mol {prmTop} {inputCoords}\n")
        f.write("quit\n")

    # Execute TLEAP and log the output
    tleapOutput: FilePath = p.join(outDir, "TLEAP.out")
    amberParams: FilePath = p.join(outDir, f"{outName}.prmtop")
    tleapCommand: str = f"tleap -f {tleapInput} > {tleapOutput}"
    run_with_log(tleapCommand, "System Parameterisation with TLEAP", amberParams)

    # Reset chain and residue IDs in Amber PDB file
    solvatedPdb: FilePath = p.join(outDir, solvatedPdb)
    inputCoords: FilePath = p.join(outDir, f"{outName}.inpcrd")
    ## reset chain and residue IDs in amber PDB
    solvatedPdb: FilePath = p.join(outDir, solvatedPdb)
    return inputCoords, amberParams, solvatedPdb
#####################################################################################
def run_with_log(
    command: str,
    stepName: str,
    expectedOutput: Optional[str] = None
) -> None:
    """
    Execute a command and log its output.

    Args:
        command (str): The command to execute.
        expectedOutput (Optional[str], optional): The path to the expected output file. Defaults to None.

    Returns:
        None
    """
    ## log preparation step
    drLogger.log_info(f"Running: {stepName}")

    ## split command into list
    if isinstance(command, str):
        command = command.split()

    try: 
        # Execute the command and capture its output
        result: subprocess.CompletedProcess[str] = subprocess.run(
            command,
            capture_output=True,
            check=True,
            text=True, 
            env = os.environ
        )
    except Exception as errorMessage:
        drSplash.print_prep_failed(errorMessage, stepName)


    # Log the command output
    if result.stdout:
        drLogger.log_info(f"Command output:\n{result.stdout}")
    if result.stderr:
        logging.error(f"Command error output:\n{result.stderr}")

    
    # Check if the expected output file exists
    if expectedOutput is None or p.isfile(expectedOutput):
        drLogger.log_info(f"{stepName} Completed Successfully")
        return
    else:
        raise FileNotFoundError(f"Expected output:\n{' '*4}{expectedOutput}\n{' '*4}{' '*4} not found for {stepName}.")


#####################################################################################