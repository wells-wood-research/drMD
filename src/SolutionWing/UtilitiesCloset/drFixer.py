## BASIC PYTHON LIBRARIES
import os
from os import path as p
from pathlib import Path
import pandas as pd
import re

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


def reset_chains_residues(templatePdb: FilePath, inputPdb: FilePath, config: dict) -> FilePath:
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

    ## unpack config to get ncaa and ligand names (if any)
    ncaaNames = config["miscInfo"].get("nonCanonicalResidueNames", [])
    ligandInfo = config.get("ligandInfo", None)
    if not ligandInfo is None:
        ligandNames = [ligand["ligandName"] for ligand in ligandInfo]
    else:
        ligandNames = []

    ## load pdb files into dataframes
    templateDf: pd.DataFrame = pdbUtils.pdb2df(templatePdb)
    inputDf: pd.DataFrame = pdbUtils.pdb2df(inputPdb)

    ## reset chains and residues for protein residues
    protFixedDf: pd.DataFrame = reset_chain_residues_protein(templateDf, inputDf, ncaaNames)

    # reset chains and residues for non-protein non-counter-ion residues
    ligFixedDf: pd.DataFrame = reset_chain_residues_ligands(templateDf, protFixedDf, ligandNames)

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
def reset_chain_residues_ligands(templateDf: pd.DataFrame, inputDf: pd.DataFrame, ligandNames: list) -> pd.DataFrame:
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
    templateLigandsDf = templateDf[templateDf["RES_NAME"].isin(ligandNames)]
    
    inputLigandsDf = inputDf[inputDf["RES_NAME"].isin(ligandNames)]
    
    outputDf = inputDf.copy()
    ## loop over chains and residues for both target and template ligands
    for (inputChain, inputChainDf), (templateChain, templateChainDf) in zip(inputLigandsDf.groupby("CHAIN_ID"), templateLigandsDf.groupby("CHAIN_ID")):
        for (inputRes, inputResDf), (templateRes, templateResDf) in zip(inputChainDf.groupby("RES_ID"), templateChainDf.groupby("RES_ID")):
            ## set chain and resid for input dataframe
            outputDf.loc[inputDf["RES_ID"] == inputRes, "CHAIN_ID"] = templateChain
            outputDf.loc[inputDf["RES_ID"] == inputRes, "RES_ID"] = templateRes

    return outputDf

def _account_for_no_CA_ncaas(pdbDf: pd.DataFrame, ncaaNames: list[str]) -> pd.DataFrame:

    ncaaRefDfs = []
    ## account for non-canonical amino acids with no CA
    for chainId, chainDf in pdbDf.groupby("CHAIN_ID"):
        for resId, resDf in chainDf.groupby("RES_ID"):
            if resDf["RES_NAME"].iloc[0] in ncaaNames:
                resAtomNames = resDf["ATOM_NAME"].tolist()
                if "CA" not in resAtomNames:
                    referenceAtom = resAtomNames[0]
                    ncaaRefDfs.append(resDf[resDf["ATOM_NAME"] == referenceAtom])
    return pd.concat(ncaaRefDfs)
##################################################################################################
def reset_chain_residues_protein(templateDf: pd.DataFrame, inputDf: pd.DataFrame, ncaaNames: list[str]) -> pd.DataFrame:
    """
    Resets chain and resid columns for protein residues

    Args:
        templateDf (pd.DataFrame): Dataframe of template PDB file.
        inputDf (pd.DataFrame): Dataframe of input PDB file.

    Returns:
       inputDf (pd.DataFrame): Updated dataframe with fixed chain and resid
    
    """
    aminoAcids = drListInitiator.get_amino_acid_residue_names()
    ## include non-canonical amino acids
    aminoAcids = list(aminoAcids) + ncaaNames

    ## create dataframes for CA atoms in both template and input dfs
    templateCaDf = templateDf[(templateDf["ATOM_NAME"] == "CA") &
                            (templateDf["RES_NAME"].isin(aminoAcids))]
    inputCaDf = inputDf[(inputDf["ATOM_NAME"] == "CA") &
                        (inputDf["RES_NAME"].isin(aminoAcids))]
    
    if len(ncaaNames) > 0:
        ncaaTemplateRefDf = _account_for_no_CA_ncaas(templateDf, ncaaNames)
        ncaaInputRefDf = _account_for_no_CA_ncaas(inputDf, ncaaNames)
        templateCaDf = pd.concat([templateCaDf, ncaaTemplateRefDf])
        inputCaDf = pd.concat([inputCaDf, ncaaInputRefDf])

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
##################################################################################################



def _parse_pdb_atom_line(line: str) -> dict | None:
    """Parse an ATOM/HETATM line using fixed-width PDB columns, with token fallback."""
    if not line.startswith(("ATOM", "HETATM")):
        return None

    line = line.rstrip()

    try:
        atom_serial = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain_id = line[21:22].strip() or "A"
        res_id = line[22:26].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        occupancy = float(line[54:60].strip()) if line[54:60].strip() else 1.0
        temp_factor = float(line[60:66].strip()) if line[60:66].strip() else 0.0
        element = line[76:78].strip() or atom_name[0].upper()

        return {
            "record": line[:6].strip(),
            "atom_serial": atom_serial,
            "atom_name": atom_name,
            "res_name": res_name,
            "chain_id": chain_id,
            "res_id": res_id,
            "x": x,
            "y": y,
            "z": z,
            "occupancy": occupancy,
            "temp_factor": temp_factor,
            "element": element,
        }
    except (ValueError, TypeError):
        pass

    tokens = line.split()
    if len(tokens) < 12:
        return None

    try:
        atom_serial = int(tokens[1])
        atom_name = tokens[2]
        res_name = tokens[3]
        chain_id = tokens[4] if len(tokens[4]) == 1 and tokens[4].isalpha() else "A"
        res_id = tokens[5]
        x = float(tokens[6])
        y = float(tokens[7])
        z = float(tokens[8])
        occupancy = float(tokens[9])
        temp_factor = float(tokens[10])
        element = tokens[11] if len(tokens) > 11 else atom_name[0].upper()

        return {
            "record": tokens[0],
            "atom_serial": atom_serial,
            "atom_name": atom_name,
            "res_name": res_name,
            "chain_id": chain_id,
            "res_id": res_id,
            "x": x,
            "y": y,
            "z": z,
            "occupancy": occupancy,
            "temp_factor": temp_factor,
            "element": element,
        }
    except (ValueError, TypeError, IndexError):
        # Some malformed records merge the atom serial and atom name without a separator,
        # e.g. 'ATOM  500000EPW  WAT ...'. Split that into serial + name before renumbering.
        serial_match = re.match(r"^(ATOM|HETATM)\s*(\d+)([A-Za-z]+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*(\S+)?$", line)
        if serial_match is None:
            return None

        atom_serial = int(serial_match.group(2))
        atom_name = serial_match.group(3)
        res_name = serial_match.group(4)
        chain_id = serial_match.group(5)[0] if serial_match.group(5) and len(serial_match.group(5)) == 1 and serial_match.group(5).isalpha() else "A"
        res_id = serial_match.group(6)
        x = float(serial_match.group(7))
        y = float(serial_match.group(8))
        z = float(serial_match.group(9))
        occupancy = float(serial_match.group(10))
        temp_factor = float(serial_match.group(11))
        element = serial_match.group(12) if serial_match.group(12) else atom_name[0].upper()

        return {
            "record": serial_match.group(1),
            "atom_serial": atom_serial,
            "atom_name": atom_name,
            "res_name": res_name,
            "chain_id": chain_id,
            "res_id": res_id,
            "x": x,
            "y": y,
            "z": z,
            "occupancy": occupancy,
            "temp_factor": temp_factor,
            "element": element,
        }


def sanitize_pdb_lines(pdb_text: str) -> str:
    """Split merged TER/ATOM rows and ensure a newline at the end of every PDB record."""
    text = pdb_text.replace("\r\n", "\n").replace("\r", "\n")
    cleaned_lines: list[str] = []

    for raw_line in text.split("\n"):
        line = raw_line.rstrip()
        if not line:
            continue

        if line.startswith("TER") and "ATOM" in line and line.find("ATOM") > 0:
            prefix = line[: line.find("ATOM")].rstrip()
            suffix = line[line.find("ATOM") :].rstrip()
            if prefix:
                cleaned_lines.append(prefix)
            if suffix:
                cleaned_lines.append(suffix)
            continue

        cleaned_lines.append(line)

    sanitized_text = "\n".join(cleaned_lines)
    if sanitized_text and not sanitized_text.endswith("\n"):
        sanitized_text += "\n"
    return sanitized_text if sanitized_text.endswith("\n") else sanitized_text + "\n"


def _fit_pdb_field(value: object, width: int, *, align: str = "right") -> str:
    """Fit a value into a fixed-width PDB field without shifting later columns."""
    text = str(value)
    if len(text) > width:
        text = text[-width:] if align == "right" else text[:width]
    return text.rjust(width) if align == "right" else text.ljust(width)


def _normalize_residue_id(res_id: object) -> int:
    """Keep residue numbering valid in the 4-digit PDB residue field by wrapping at 9999."""
    try:
        value = int(str(res_id).strip())
    except (TypeError, ValueError):
        return 1

    wrapped = ((value - 1) % 9999) + 1
    return wrapped


def _wrap_pdb_counter(value: object, limit: int, *, default: int = 1) -> int:
    """Wrap serial IDs back within the valid PDB field width without mutating residue numbering."""
    try:
        number = int(str(value).strip())
    except (TypeError, ValueError):
        return default
    if number <= 0:
        return default
    if limit <= 0:
        return default
    return ((number - 1) % limit) + 1 if number > limit else number


def format_pdb_atom_record(
    *,
    record: str,
    atom_serial: int,
    atom_name: str,
    res_name: str,
    chain_id: str,
    res_id: int,
    x: float,
    y: float,
    z: float,
    occupancy: float,
    bfactor: float,
    element: str,
) -> str:
    """Write a fixed-width PDB atom line matching the canonical field layout."""
    record_name = str(record).strip().upper()
    if record_name not in {"ATOM", "HETATM"}:
        record_name = "ATOM"

    atom_serial = _wrap_pdb_counter(atom_serial, 100000)
    res_id = int(res_id) if str(res_id).strip() else 1

    line = [" "] * 80

    # 1-6: record type
    line[0:6] = list(f"{record_name:<6}")

    # 7-11: atom serial number
    line[6:11] = list(f"{int(atom_serial):>5d}")

    # 12: blank
    line[11] = " "

    # 13-16: atom name
    atom_name_field = str(atom_name).strip()[:4].ljust(4)
    line[12:16] = list(atom_name_field)

    # 17: alt loc
    line[16] = " "

    # 18-20: residue name
    res_name_field = str(res_name).strip()[:3].ljust(3)
    line[17:20] = list(res_name_field)

    # 21: blank
    line[20] = " "

    # 22: chain identifier
    line[21] = str(chain_id or "A").strip()[:1] or "A"

    # 23-26: residue sequence number
    line[22:26] = list(f"{int(res_id):>4d}")

    # 27: insertion code
    line[26] = " "

    # 28-30: blank
    for idx in range(27, 30):
        line[idx] = " "

    # 31-38: x
    line[30:38] = list(f"{float(x):>8.3f}")

    # 39-46: y
    line[38:46] = list(f"{float(y):>8.3f}")

    # 47-54: z
    line[46:54] = list(f"{float(z):>8.3f}")

    # 55-60: occupancy
    line[54:60] = list(f"{float(occupancy):>6.2f}")

    # 61-66: B-factor
    line[60:66] = list(f"{float(bfactor):>6.2f}")

    # 67-72: blank
    for idx in range(66, 72):
        line[idx] = " "

    # 73-76: segment identifier blank
    for idx in range(72, 76):
        line[idx] = " "

    # 77-78: element
    element_field = str(element).strip()[:2].rjust(2)
    line[76:78] = list(element_field)

    # 79-80: charge
    line[78:80] = list("  ")

    return "".join(line)


def write_pdb_atom(
    *,
    record: str,
    atom_serial: int,
    atom_name: str,
    res_name: str,
    chain: str,
    res_num: int,
    x: float,
    y: float,
    z: float,
    bfactor: float,
    occupancy: float,
    element: str,
) -> str:
    """Compatibility wrapper so the renumberer writes through the canonical formatter."""
    return format_pdb_atom_record(
        record=record,
        atom_serial=atom_serial,
        atom_name=atom_name,
        res_name=res_name,
        chain_id=chain,
        res_id=res_num,
        x=x,
        y=y,
        z=z,
        occupancy=occupancy,
        bfactor=bfactor,
        element=element,
    )


##################################################################################################

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


def reset_chains_residues(templatePdb: FilePath, inputPdb: FilePath, config: dict) -> FilePath:
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

    ## unpack config to get ncaa and ligand names (if any)
    ncaaNames = config["miscInfo"].get("nonCanonicalResidueNames", [])
    ligandInfo = config.get("ligandInfo", None)
    if not ligandInfo is None:
        ligandNames = [ligand["ligandName"] for ligand in ligandInfo]
    else:
        ligandNames = []

    ## load pdb files into dataframes
    templateDf: pd.DataFrame = pdbUtils.pdb2df(templatePdb)
    inputDf: pd.DataFrame = pdbUtils.pdb2df(inputPdb)

    ## reset chains and residues for protein residues
    protFixedDf: pd.DataFrame = reset_chain_residues_protein(templateDf, inputDf, ncaaNames)

    # reset chains and residues for non-protein non-counter-ion residues
    ligFixedDf: pd.DataFrame = reset_chain_residues_ligands(templateDf, protFixedDf, ligandNames)

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
def reset_chain_residues_ligands(templateDf: pd.DataFrame, inputDf: pd.DataFrame, ligandNames: list) -> pd.DataFrame:
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
    templateLigandsDf = templateDf[templateDf["RES_NAME"].isin(ligandNames)]
    
    inputLigandsDf = inputDf[inputDf["RES_NAME"].isin(ligandNames)]
    
    outputDf = inputDf.copy()
    ## loop over chains and residues for both target and template ligands
    for (inputChain, inputChainDf), (templateChain, templateChainDf) in zip(inputLigandsDf.groupby("CHAIN_ID"), templateLigandsDf.groupby("CHAIN_ID")):
        for (inputRes, inputResDf), (templateRes, templateResDf) in zip(inputChainDf.groupby("RES_ID"), templateChainDf.groupby("RES_ID")):
            ## set chain and resid for input dataframe
            outputDf.loc[inputDf["RES_ID"] == inputRes, "CHAIN_ID"] = templateChain
            outputDf.loc[inputDf["RES_ID"] == inputRes, "RES_ID"] = templateRes

    return outputDf

def _account_for_no_CA_ncaas(pdbDf: pd.DataFrame, ncaaNames: list[str]) -> pd.DataFrame:

    ncaaRefDfs = []
    ## account for non-canonical amino acids with no CA
    for chainId, chainDf in pdbDf.groupby("CHAIN_ID"):
        for resId, resDf in chainDf.groupby("RES_ID"):
            if resDf["RES_NAME"].iloc[0] in ncaaNames:
                resAtomNames = resDf["ATOM_NAME"].tolist()
                if "CA" not in resAtomNames:
                    referenceAtom = resAtomNames[0]
                    ncaaRefDfs.append(resDf[resDf["ATOM_NAME"] == referenceAtom])
    return pd.concat(ncaaRefDfs)
##################################################################################################
def reset_chain_residues_protein(templateDf: pd.DataFrame, inputDf: pd.DataFrame, ncaaNames: list[str]) -> pd.DataFrame:
    """
    Resets chain and resid columns for protein residues

    Args:
        templateDf (pd.DataFrame): Dataframe of template PDB file.
        inputDf (pd.DataFrame): Dataframe of input PDB file.

    Returns:
       inputDf (pd.DataFrame): Updated dataframe with fixed chain and resid
    
    """
    aminoAcids = drListInitiator.get_amino_acid_residue_names()
    ## include non-canonical amino acids
    aminoAcids = list(aminoAcids) + ncaaNames

    ## create dataframes for CA atoms in both template and input dfs
    templateCaDf = templateDf[(templateDf["ATOM_NAME"] == "CA") &
                            (templateDf["RES_NAME"].isin(aminoAcids))]
    inputCaDf = inputDf[(inputDf["ATOM_NAME"] == "CA") &
                        (inputDf["RES_NAME"].isin(aminoAcids))]
    
    if len(ncaaNames) > 0:
        ncaaTemplateRefDf = _account_for_no_CA_ncaas(templateDf, ncaaNames)
        ncaaInputRefDf = _account_for_no_CA_ncaas(inputDf, ncaaNames)
        templateCaDf = pd.concat([templateCaDf, ncaaTemplateRefDf])
        inputCaDf = pd.concat([inputCaDf, ncaaInputRefDf])

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
##################################################################################################



def _parse_pdb_atom_line(line: str) -> dict | None:
    """Parse an ATOM/HETATM line using fixed-width PDB columns, with token fallback."""
    if not line.startswith(("ATOM", "HETATM")):
        return None

    line = line.rstrip()

    try:
        atom_serial = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain_id = line[21:22].strip() or "A"
        res_id = line[22:26].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        occupancy = float(line[54:60].strip()) if line[54:60].strip() else 1.0
        temp_factor = float(line[60:66].strip()) if line[60:66].strip() else 0.0
        element = line[76:78].strip() or atom_name[0].upper()

        return {
            "record": line[:6].strip(),
            "atom_serial": atom_serial,
            "atom_name": atom_name,
            "res_name": res_name,
            "chain_id": chain_id,
            "res_id": res_id,
            "x": x,
            "y": y,
            "z": z,
            "occupancy": occupancy,
            "temp_factor": temp_factor,
            "element": element,
        }
    except (ValueError, TypeError):
        pass

    tokens = line.split()
    if len(tokens) < 9:
        return None

    try:
        atom_serial = int(tokens[1])
        atom_name = tokens[2]
        res_name = tokens[3]

        if len(tokens) >= 10 and tokens[4].isdigit():
            chain_id = "A"
            res_id = tokens[4]
            x = float(tokens[5])
            y = float(tokens[6])
            z = float(tokens[7])
            occupancy = float(tokens[8]) if len(tokens) > 8 else 1.0
            temp_factor = float(tokens[9]) if len(tokens) > 9 else 0.0
            element = tokens[10] if len(tokens) > 10 else atom_name[0].upper()
        else:
            chain_id = tokens[4] if len(tokens[4]) == 1 and tokens[4].isalpha() else "A"
            res_id = tokens[5] if len(tokens) > 5 else "1"
            x = float(tokens[6])
            y = float(tokens[7])
            z = float(tokens[8])
            occupancy = float(tokens[9]) if len(tokens) > 9 else 1.0
            temp_factor = float(tokens[10]) if len(tokens) > 10 else 0.0
            element = tokens[11] if len(tokens) > 11 else atom_name[0].upper()

        return {
            "record": tokens[0],
            "atom_serial": atom_serial,
            "atom_name": atom_name,
            "res_name": res_name,
            "chain_id": chain_id,
            "res_id": res_id,
            "x": x,
            "y": y,
            "z": z,
            "occupancy": occupancy,
            "temp_factor": temp_factor,
            "element": element,
        }
    except (ValueError, TypeError, IndexError):
        # Some malformed records merge the atom serial and atom name without a separator,
        # e.g. 'ATOM  500000EPW  WAT ...'. Split that into serial + name before renumbering.
        serial_match = re.match(r"^(ATOM|HETATM)\s*(\d+)([A-Za-z]+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*(\S+)?$", line)
        if serial_match is None:
            return None

        atom_serial = int(serial_match.group(2))
        atom_name = serial_match.group(3)
        res_name = serial_match.group(4)
        chain_id = serial_match.group(5)[0] if serial_match.group(5) and len(serial_match.group(5)) == 1 and serial_match.group(5).isalpha() else "A"
        res_id = serial_match.group(6)
        x = float(serial_match.group(7))
        y = float(serial_match.group(8))
        z = float(serial_match.group(9))
        occupancy = float(serial_match.group(10))
        temp_factor = float(serial_match.group(11))
        element = serial_match.group(12) if serial_match.group(12) else atom_name[0].upper()

        return {
            "record": serial_match.group(1),
            "atom_serial": atom_serial,
            "atom_name": atom_name,
            "res_name": res_name,
            "chain_id": chain_id,
            "res_id": res_id,
            "x": x,
            "y": y,
            "z": z,
            "occupancy": occupancy,
            "temp_factor": temp_factor,
            "element": element,
        }


def sanitize_pdb_lines(pdb_text: str) -> str:
    """Split merged TER/ATOM rows and ensure a newline at the end of every PDB record."""
    text = pdb_text.replace("\r\n", "\n").replace("\r", "\n")
    cleaned_lines: list[str] = []

    for raw_line in text.split("\n"):
        line = raw_line.rstrip()
        if not line:
            continue

        if line.startswith("TER") and "ATOM" in line and line.find("ATOM") > 0:
            prefix = line[: line.find("ATOM")].rstrip()
            suffix = line[line.find("ATOM") :].rstrip()
            if prefix:
                cleaned_lines.append(prefix)
            if suffix:
                cleaned_lines.append(suffix)
            continue

        cleaned_lines.append(line)

    sanitized_text = "\n".join(cleaned_lines)
    if sanitized_text and not sanitized_text.endswith("\n"):
        sanitized_text += "\n"
    return sanitized_text if sanitized_text.endswith("\n") else sanitized_text + "\n"


def _fit_pdb_field(value: object, width: int, *, align: str = "right") -> str:
    """Fit a value into a fixed-width PDB field without shifting later columns."""
    text = str(value)
    if len(text) > width:
        text = text[-width:] if align == "right" else text[:width]
    return text.rjust(width) if align == "right" else text.ljust(width)


def _normalize_residue_id(res_id: object) -> int:
    """Keep residue numbering valid in the 4-digit PDB residue field by wrapping at 9999."""
    try:
        value = int(str(res_id).strip())
    except (TypeError, ValueError):
        return 1

    wrapped = ((value - 1) % 9999) + 1
    return wrapped



def _coerce_pdb_integer(value: object, *, default: int = 1) -> int:
    """Extract the first integer-like token from malformed PDB fields."""
    if value is None:
        return default

    text = str(value).strip()
    if not text:
        return default

    match = re.search(r"[-+]?\d+", text)
    if match is None:
        return default

    try:
        return int(match.group(0))
    except ValueError:
        return default


def _coerce_pdb_float(value: object, *, default: float = 0.0) -> float:
    """Extract the first float-like token from malformed PDB coordinate fields."""
    if value is None:
        return default

    text = str(value).strip()
    if not text:
        return default

    match = re.search(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)", text)
    if match is None:
        return default

    try:
        return float(match.group(0))
    except ValueError:
        return default


def _repair_malformed_atom_line(line: str) -> str:
    """Repair a malformed ATOM/HETATM row by rebuilding it in canonical PDB columns."""
    if not line.startswith(("ATOM", "HETATM")):
        return line

    tokens = line.split()
    if len(tokens) < 9:
        return line

    try:
        record = tokens[0]
        atom_serial = _coerce_pdb_integer(tokens[1], default=1)
        atom_name = tokens[2]
        res_name = tokens[3]

        if len(tokens) >= 11 and tokens[4].isalpha() and len(tokens[4]) == 1:
            chain_id = tokens[4]
            res_id = tokens[5] if len(tokens) > 5 else "1"
            coord_tokens = [t for t in tokens[6:] if re.search(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)", t)]
            x = _coerce_pdb_float(coord_tokens[0], default=0.0) if len(coord_tokens) > 0 else 0.0
            y = _coerce_pdb_float(coord_tokens[1], default=0.0) if len(coord_tokens) > 1 else 0.0
            z = _coerce_pdb_float(coord_tokens[2], default=0.0) if len(coord_tokens) > 2 else 0.0
            occupancy = _coerce_pdb_float(tokens[9], default=1.0) if len(tokens) > 9 else 1.0
            bfactor = _coerce_pdb_float(tokens[10], default=0.0) if len(tokens) > 10 else 0.0
            element = tokens[11] if len(tokens) > 11 else atom_name[0].upper()
        else:
            chain_id = "A"
            res_id = tokens[4] if len(tokens) > 4 else "1"
            coord_tokens = [t for t in tokens[5:] if re.search(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)", t)]
            x = _coerce_pdb_float(coord_tokens[0], default=0.0) if len(coord_tokens) > 0 else 0.0
            y = _coerce_pdb_float(coord_tokens[1], default=0.0) if len(coord_tokens) > 1 else 0.0
            z = _coerce_pdb_float(coord_tokens[2], default=0.0) if len(coord_tokens) > 2 else 0.0
            occupancy = _coerce_pdb_float(tokens[8], default=1.0) if len(tokens) > 8 else 1.0
            bfactor = _coerce_pdb_float(tokens[9], default=0.0) if len(tokens) > 9 else 0.0
            element = tokens[10] if len(tokens) > 10 else atom_name[0].upper()

        return write_pdb_atom(
            record=record,
            atom_serial=atom_serial,
            atom_name=atom_name,
            res_name=res_name,
            chain=chain_id,
            res_num=_coerce_pdb_integer(res_id),
            x=x,
            y=y,
            z=z,
            bfactor=bfactor,
            occupancy=occupancy,
            element=element,
        )
    except (ValueError, TypeError, IndexError):
        return line


def repair_pdb_file(pdb_path: str) -> str:
    """Repair malformed ATOM/HETATM records in a PDB file before parsing it."""
    if not Path(pdb_path).exists():
        return pdb_path

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as pdb:
        original_text = pdb.read()

    fixed_lines: list[str] = []
    for raw_line in original_text.splitlines():
        line = raw_line.rstrip()
        if not line:
            continue
        if line.startswith(("ATOM", "HETATM")):
            fixed_lines.append(_repair_malformed_atom_line(line))
        else:
            fixed_lines.append(line)

    with open(pdb_path, "w", encoding="utf-8") as pdb:
        pdb.write("\n".join(fixed_lines) + "\n")

    return pdb_path


def renumber_pdb_atom_serials(pdb_path: str) -> None:
    """Rewrite malformed ATOM/HETATM rows into canonical 80-column PDB records."""
    repair_pdb_file(pdb_path)
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as pdb:
        original_text = pdb.read()

    sanitized_text = sanitize_pdb_lines(original_text)
    lines = sanitized_text.split("\n")

    sanitized_lines: list[str] = []
    next_atom_serial = 1
    residue_index_map: dict[tuple[str, str], int] = {}

    for raw_line in lines:
        if not raw_line.strip():
            continue

        line = raw_line.rstrip()
        if not line.startswith(("ATOM", "HETATM")):
            sanitized_lines.append(line + "\n")
            continue

        parsed = _parse_pdb_atom_line(line)
        if parsed is None:
            sanitized_lines.append(line + "\n")
            continue

        atom_serial = next_atom_serial
        chain_id = (parsed["chain_id"] or "A")[:1]
        res_key = (chain_id, str(parsed["res_id"]).strip() or "1")
        if res_key not in residue_index_map:
            residue_index_map[res_key] = len(residue_index_map) + 1
        res_id_int = residue_index_map[res_key]

        sanitized_lines.append(
            write_pdb_atom(
                record=parsed["record"],
                atom_serial=atom_serial,
                atom_name=parsed["atom_name"],
                res_name=parsed["res_name"],
                chain=chain_id,
                res_num=res_id_int,
                x=parsed["x"],
                y=parsed["y"],
                z=parsed["z"],
                bfactor=parsed["temp_factor"],
                occupancy=parsed["occupancy"],
                element=parsed["element"],
            ) + "\n"
        )
        next_atom_serial += 1

    with open(pdb_path, "w", encoding="utf-8") as pdb:
        pdb.writelines(sanitized_lines)
