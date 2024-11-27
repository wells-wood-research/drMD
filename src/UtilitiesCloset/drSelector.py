## BASIC PYTHON LIBRARIES
import pandas as pd
from typing import List
import numpy as np

## drMD LIBRARIES
from UtilitiesCloset import drListInitiator

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

##  CLEAN CODE
from typing import Dict, Callable, Optional, Tuple
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

#######################################################################
def get_atom_indexes(selection: Dict, pdbFile: FilePath) -> List[int]:
    """
    This function takes a selection dictionary and a PDB file as input,
    and returns a list of atom indexes that correspond to the selection.

    Args:
        selection (Dict): A dictionary containing the selection type and input.
        pdbFile (FilePath): The path to the PDB file.

    Returns:
        atomIndexes (List[int]): A list of atom indexes that correspond to the selection.
    """
    # Read PDB file into a DataFrame
    pdbDf: pd.DataFrame = pdbUtils.pdb2df(pdbFile)

    # Initialize lists of amino acid residue names, backbone atom names,
    # solvent residue names, and ion residue names
    aminoAcidResNames = drListInitiator.get_amino_acid_residue_names()
    backboneAtomNames = drListInitiator.get_backbone_atom_names()
    solventResNames = drListInitiator.get_solvent_residue_names()
    ionResidueNames = drListInitiator.get_ion_residue_names()

    # Initialize list of atom indexes
    atomIndexes: List[int] = []
    # Check selection type and find required atom indexes
    if selection["keyword"] == "all":
        # Find indexes for all atoms
        atomIndexes = pdbDf.index.tolist()
    elif selection["keyword"] == "backbone":
        # Find indexes for all backbone atoms
        backboneDf: pd.DataFrame = pdbDf[pdbDf["ATOM_NAME"].isin(backboneAtomNames) &
                                          pdbDf["RES_NAME"].isin(aminoAcidResNames)]
        atomIndexes: List[int] = backboneDf.index.tolist()
    elif selection["keyword"] == "protein":
        # Find indexes for all protein atoms
        proteinDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"].isin(aminoAcidResNames)]
        atomIndexes: List[int] = proteinDf.index.tolist()
    elif selection["keyword"] == "water":
        # Find indexes for water molecules
        waterDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"].isin(solventResNames) &
                                       ~pdbDf["RES_NAME"].isin([aminoAcidResNames])]
        atomIndexes: List[int] = waterDf.index.tolist()
    elif selection["keyword"] == "ions":
        # Find indexes for ions
        ionDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"].isin(ionResidueNames) &
                                    ~pdbDf["RES_NAME"].isin([aminoAcidResNames])]

        atomIndexes: List[int] = ionDf.index.tolist()
    elif selection["keyword"] == "ligand":
        # Find indexes for all ligand / organics / oddball molecules
        ligandDf: pd.DataFrame = pdbDf[~pdbDf["RES_NAME"].isin(aminoAcidResNames) &
                                        ~pdbDf["RES_NAME"].isin(solventResNames) & 
                                        ~pdbDf["RES_NAME"].isin(ionResidueNames)]
        atomIndexes: List[int] =  ligandDf.index.tolist()
    elif selection["keyword"] == "custom":
        # Find indexes for whole residue selections
        customSelection: List[Dict] = selection["customSelection"]
        ## check for "all" selections in selectionSytax
        for selection in customSelection:
            ## init empty list to store conditions
            selectionConditions: list  = []
            ## for each selection key...
            for selctionKey in ["CHAIN_ID", "RES_NAME", "RES_ID", "ATOM_NAME"]:
                ## "all" inputs create no condition
                if selection[selctionKey] == "_":
                    continue
                ## use list inputs to create a .isin() condition
                if isinstance(selection[selctionKey], list):
                    selectionConditions.append(pdbDf[selctionKey].isin(selection[selctionKey]))
                ## use single input to create a == condition
                elif isinstance(selection[selctionKey], (str,int)):
                    selectionConditions.append(pdbDf[selctionKey] == selection[selctionKey])
                    
            selectionConditionsCombined = np.logical_and.reduce(selectionConditions)
            selectionDf = pdbDf[selectionConditionsCombined]
            ## add to atomIndexes list
            atomIndexes += selectionDf.index.tolist()

    return atomIndexes

#######################################################################
def slice_pdb_file(trajectorySelection: Dict, inputPdb: FilePath, outputPdb: FilePath) -> FilePath:
    """
    This function takes a selection dictionary and a PDB file as input,
    And returns a sliced PDB file based on the selection
    """
    atomIndexes: list = []
    for selection in trajectorySelection:
        atomIndexes.extend(get_atom_indexes(selection["selection"], inputPdb))
    inputDf: pd.DataFrame = pdbUtils.pdb2df(inputPdb)

    outputDf: pd.DataFrame = inputDf.iloc[atomIndexes]

    pdbUtils.df2pdb(outputDf, outputPdb)

    return outputPdb



#######################################################################