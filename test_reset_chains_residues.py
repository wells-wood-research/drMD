import os
from os import path as p
from pdbUtils import pdbUtils

def inputs():
    templatePdb = "/home/esp/scriptDevelopment/drMD/00_MERT/outputs/P450_Aposim1_heme_pose1/00_prep/WHOLE/P450_Aposim1_heme_pose1.pdb"
    inputPdb = "/home/esp/scriptDevelopment/drMD/00_MERT/outputs/P450_Aposim1_heme_pose1/00_prep/WHOLE/P450_Aposim1_heme_pose1_solvated.pdb"

    return templatePdb, inputPdb

def main():
    templatePdb, inputPdb = inputs()

    templateDf = pdbUtils.pdb2df(templatePdb)
    inputDf = pdbUtils.pdb2df(inputPdb)


    protFixedDf = reset_chain_residues_protein(templateDf, inputDf)

    ligFixedDf = reset_chain_residues_ligands(templateDf, protFixedDf)


    pdbUtils.df2pdb(ligFixedDf, "./test.pdb")


def reset_chain_residues_ligands(templateDf, inputDf):
    aminoAcids = get_amino_acid_residue_names()
    counterIons = {"Na+", "Cl-"}

    templateLigandsDf = templateDf[~templateDf["RES_NAME"].isin(aminoAcids) & 
                                   ~templateDf["RES_NAME"].isin(counterIons)]
    
    inputLigandsDf = inputDf[~inputDf["RES_NAME"].isin(aminoAcids) &
                             ~inputDf["RES_NAME"].isin(counterIons)]

    for (inputChain, inputChainDf), (templateChain, templateChainDf) in zip(inputLigandsDf.groupby("CHAIN_ID"), templateLigandsDf.groupby("CHAIN_ID")):
        for (inputRes, inputResDf), (templateRes, templateResDf) in zip(inputChainDf.groupby("RES_ID"), templateChainDf.groupby("RES_ID")):
            inputDf.loc[inputDf["RES_ID"] == inputRes, "CHAIN_ID"] = templateChain
            inputDf.loc[inputDf["RES_ID"] == inputRes, "RES_ID"] = templateRes

    return inputDf

def reset_chain_residues_protein(templateDf, inputDf):
    templateCaDf = templateDf[templateDf["ATOM_NAME"] == "CA"]
    inputCaDf = inputDf[inputDf["ATOM_NAME"] == "CA"]

    for templateCa, inputCa in zip(templateCaDf.iterrows(), inputCaDf.iterrows()):
        inputResidueId = inputCa[1]["RES_ID"]

        targetResidueId = templateCa[1]["RES_ID"]
        targetChainId = templateCa[1]["CHAIN_ID"]
        inputDf.loc[inputDf["RES_ID"] == inputResidueId, "CHAIN_ID"] = targetChainId
        inputDf.loc[inputDf["RES_ID"] == inputResidueId, "RES_ID"] = targetResidueId
    return inputDf




def get_amino_acid_residue_names() -> set:
    """
    Returns a list of the amino acid names.
    """
    return  {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN',
            'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS',
            'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            'ASH', 'GLH', 'HIP', 'HIE', 'HID', 'CYX', 'CYM', 'LYN',   ## oddball protonations
            'ACE','NME','NHE'}              ## caps

if __name__ == "__main__":
    main()