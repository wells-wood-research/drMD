## BASIC PYTHON LIBRARIES
import os
from os import path as p
import pandas as pd

## drMD LIBRARIES
from UtilitiesCloset import drListInitiator

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

##  CLEAN CODE
from typing import Dict, Callable, List
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

##################################################################################################
def reset_atom_numbers(pdbFile: str) -> str:
    """
    Resets the atom numbers in a PDB file.

    Parameters:
        pdbFile (str): Path to the PDB file.

    Returns:
        str: Path to the modified PDB file.
    """

    pdbDf = pdbUtils.pdb2df(pdbFile)
    pdbDf["ATOM_ID"] = range(1, len(pdbDf) + 1)

    pdbUtils.df2pdb(pdbDf, pdbFile)

    return pdbFile



def reset_chains(refPdb: FilePath, inputPdb: FilePath, ligandNames = []) -> FilePath:

    refDf = pdbUtils.pdb2df(refPdb)
    inputDf = pdbUtils.pdb2df(inputPdb)

    refLigDf = refDf[refDf["RES_NAME"].isin(ligandNames)]

    inputLigDf = inputDf[inputDf["RES_NAME"].isin(ligandNames)]




    ligChainMap: dict = {}
    for chainId, chainDf in refLigDf.groupby("CHAIN_ID"):
        for resId, resDf in chainDf.groupby("RES_ID"):
            ligChainMap[resId] = chainId

    for chainId, chainDf in refLigDf.groupby("CHAIN_ID"):
        for resId, resDf in chainDf.groupby("RES_ID"):
            pass
##################################################################################################


def reset_chains_residues(templatePdb: FilePath, inputPdb: FilePath) -> FilePath:
    """
    New implementation of reset_chain_residues function
    resets the chain and resid columns of a PDB file using a template
    Crawls though CA atoms to find the chain and resid for proteins
    Uses "NOT A PROT or COUNTER-ION" logic to fnd ligand chain and resid

    Overwrites inputPdb with fixed residues

    Args:
        templatePdb (FilePath): Path to the template PDB file.
        inputPdb (FilePath): Path to the input PDB file.

    Returns:
        inputPdb (FilePath): Path to the modified PDB file.
    """



    ## load pdb files into dataframes
    templateDf: pd.DataFrame = pdbUtils.pdb2df(templatePdb)
    inputDf: pd.DataFrame = pdbUtils.pdb2df(inputPdb)

    ## reset chains and residues for protein residues
    protFixedDf: pd.DataFrame = reset_chain_residues_protein(templateDf, inputDf)

    # reset chains and residues for non-protein non-counter-ion residues
    ligFixedDf: pd.DataFrame = reset_chain_residues_ligands(templateDf, protFixedDf)

    waterIonsFixedDf: pd.DataFrame = reset_water_ions(ligFixedDf)

    ## overwrite original pdb file with fixed residues and chains
    pdbUtils.df2pdb(waterIonsFixedDf, inputPdb)
    return inputPdb
##################################################################################################
def reset_water_ions(inputDf: pd.DataFrame) -> pd.DataFrame:
    """
    Sets the CHAIN_ID for counter-ions and water to be different to that of the rest of the system
    """
    ## init sets for matching
    counterIons: set = {"Na+", "Cl-"}
    water: set = {"HOH", "WAT"}



    ## create dataframes for water and counter-ions
    nonWaterIonsDf = inputDf[(~inputDf["RES_NAME"].isin(counterIons)) &
                             (~inputDf["RES_NAME"].isin(water))]
    
    maxChain = nonWaterIonsDf["CHAIN_ID"].max()

    if maxChain in ["Y", "Z"] or len(maxChain) != 1:
        ionChain = " "
        waterChain = " "
    else:
        ionChain = chr(ord(maxChain) + 1)
        waterChain = chr(ord(maxChain) + 2)
    ## set CHAIN_ID for counter-ions and water to be different to that of the rest of the system
    outputDf = inputDf.copy()
    outputDf.loc[outputDf["RES_NAME"].isin(counterIons), "CHAIN_ID"] = ionChain
    outputDf.loc[outputDf["RES_NAME"].isin(water), "CHAIN_ID"] = waterChain

    return outputDf
##################################################################################################
def reset_chain_residues_ligands(templateDf: pd.DataFrame, inputDf: pd.DataFrame) -> pd.DataFrame:
    """
    Finds ligands in both template and input dataframes using "NOT A PROT or COUNTER-ION" logic
    Resets chain and resid columns for ligands

    Args:
        templateDf (pd.DataFrame): Dataframe of template PDB file.
        inputDf (pd.DataFrame): Dataframe of input PDB file.

    Returns:
       inputDf (pd.DataFrame): Updated dataframe with fixed chain and resid
    
    """

    ## init sets of amino acids and counter ions residue names
    aminoAcids = drListInitiator.get_amino_acid_residue_names()
    counterIonsAndWater = {"Na+", "Cl-", "HOH", "WAT"}


    ## create dataframes for ligands
    templateLigandsDf = templateDf[~templateDf["RES_NAME"].isin(aminoAcids) & 
                                   ~templateDf["RES_NAME"].isin(counterIonsAndWater)]
    
    inputLigandsDf = inputDf[~inputDf["RES_NAME"].isin(aminoAcids) &
                             ~inputDf["RES_NAME"].isin(counterIonsAndWater)]
    
    outputDf = inputDf.copy()
    ## loop over chains and residues for both target and template ligands
    for (inputChain, inputChainDf), (templateChain, templateChainDf) in zip(inputLigandsDf.groupby("CHAIN_ID"), templateLigandsDf.groupby("CHAIN_ID")):
        for (inputRes, inputResDf), (templateRes, templateResDf) in zip(inputChainDf.groupby("RES_ID"), templateChainDf.groupby("RES_ID")):
            ## set chain and resid for input dataframe
            outputDf.loc[inputDf["RES_ID"] == inputRes, "CHAIN_ID"] = templateChain
            outputDf.loc[inputDf["RES_ID"] == inputRes, "RES_ID"] = templateRes

    return outputDf
##################################################################################################
def reset_chain_residues_protein(templateDf: pd.DataFrame, inputDf: pd.DataFrame) -> pd.DataFrame:
    """
    Resets chain and resid columns for protein residues

    Args:
        templateDf (pd.DataFrame): Dataframe of template PDB file.
        inputDf (pd.DataFrame): Dataframe of input PDB file.

    Returns:
       inputDf (pd.DataFrame): Updated dataframe with fixed chain and resid
    
    """
    aminoAcids = drListInitiator.get_amino_acid_residue_names()
    ## create dataframes for CA atoms in both template and input dfs
    templateCaDf = templateDf[(templateDf["ATOM_NAME"] == "CA") &
                            (templateDf["RES_NAME"].isin(aminoAcids))]
    inputCaDf = inputDf[(inputDf["ATOM_NAME"] == "CA") &
                        (inputDf["RES_NAME"].isin(aminoAcids))]

    outputDf = inputDf.copy()
    ## loop over CA atoms for both template and input dfs
    for templateCa, inputCa in zip(templateCaDf.iterrows(), inputCaDf.iterrows()):
        ## extract chain and resid for both template and input dfs
        inputResidueId = inputCa[1]["RES_ID"]
        targetResidueId = templateCa[1]["RES_ID"]
        targetChainId = templateCa[1]["CHAIN_ID"]
        ## reset chain and resid in inputDf
        outputDf.loc[inputDf["RES_ID"] == inputResidueId, "CHAIN_ID"] = targetChainId
        outputDf.loc[inputDf["RES_ID"] == inputResidueId, "RES_ID"] = targetResidueId
    return outputDf

# ##################################################################################################
def fix_atom_names(df): 
    # deal with unwanted apostrophies (prime)
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].str.replace("'", "")
    # deal with numbers at the beginning of atom names
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].replace(r'^(\d+)(.+)$', r'\2\1', regex=True)
    # deal with "A" at the start of atom name
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].apply(lambda x: x.lstrip('A') if x.startswith('A') else x)

    ## ensure unique names
    count_series = df.groupby('ATOM_NAME').cumcount()
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'] + "_" +count_series.astype(str)
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].str.replace("_0", "")
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].str.replace("_", "")

    return df 

##################################################################################################

if __name__ == "__main__":
    goodPdb = "/home/esp/scriptDevelopment/drMD/03_test_outputs/A0A0D2XFD3_TPA_1/00_prep/WHOLE/A0A0D2XFD3_TPA_1.pdb"
    badPdb = "/home/esp/scriptDevelopment/drMD/03_test_outputs/A0A0D2XFD3_TPA_1/00_prep/WHOLE/A0A0D2XFD3_TPA_1_solvated.pdb"
    reset_chains_residues(goodPdb, badPdb)