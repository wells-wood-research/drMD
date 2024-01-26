## BASIC LIBS
import os
from os import path as p
from shutil import copy, rmtree
## INPUT LIBS
import yaml
import argpass
## CUSTOM LIBS
from pdbUtils import *
from prep_drMD import *
import cleanup_drMD as cleanup
## Multiprocessing
from concurrent.futures import ThreadPoolExecutor

CURRENT_DIR = p.dirname(p.abspath(__file__))
DRMD_SCRIPT_PATH = p.join(CURRENT_DIR, 'drMD.py')


#####################################################################################################
def read_inputs():
    ## create an argpass parser, read config file, snip off ".py" if on the end of file
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()

    batchConfig=args.config
    ## Read config.yaml into a dictionary
    with open(batchConfig,"r") as yamlFile:
        batchConfig = yaml.safe_load(yamlFile) 
    return batchConfig

#####################################################################################################
def process_pdb_file(pdbFile, pdbDir, outDir, yamlDir, simInfo, topDir, batchConfig):
    # Skip if not a PDB file
    fileData = p.splitext(pdbFile)
    if not fileData[1] == ".pdb":
        return

    # Extract basic info, make dirs
    protName = fileData[0]
    pdbPath = p.join(pdbDir, pdbFile)
    runDir = p.join(outDir, protName)
    os.makedirs(runDir, exist_ok=True)

    # Convert to DataFrame, extract rest of info
    pdbDf = pdb2df(pdbPath)
    proteinInfo, ligandInfo = extract_info(pdbDf, pdbDir, protName, yamlDir, batchConfig)

    # Get path info
    pathInfo = {
        "inputDir": pdbDir,
        "inputPdb": pdbPath,
        "outputDir": runDir,
        "outputName": protName
    }

    # Combine infos into one dict (ignore ligInfo if empty)
    if ligandInfo is None:
        config = {
            "pathInfo": pathInfo,
            "proteinInfo": proteinInfo,
            "simulationInfo": simInfo
        }
    else:
        config = {
            "pathInfo": pathInfo,
            "proteinInfo": proteinInfo,
            "ligandInfo": ligandInfo,
            "simulationInfo": simInfo
        }
    # Write config to YAML
    yamlFile = p.join(yamlDir, f"{protName}_config.yaml")
    with open(yamlFile, "w") as f:
        yaml.dump(config, f, default_flow_style=False)

    # Run drMD MD protocol
    os.chdir(topDir)
    drMdCommand = f"python {DRMD_SCRIPT_PATH} --config {yamlFile}"
    run(drMdCommand, shell=True)

######################################################################################################
def main():
    ## SORT OUT DIRECTORIES
    topDir = os.getcwd()
    batchConfig = read_inputs()
    outDir = batchConfig["pathInfo"]["outputDir"]
    yamlDir = p.join(outDir,"00_configs")
    pdbDir = batchConfig["pathInfo"]["inputDir"]
    simInfo = batchConfig["simulationInfo"]
    cpuCount = batchConfig["generalInfo"]["cpuCount"]
    os.makedirs(yamlDir,exist_ok=True)
    if cpuCount == 1:
        run_serial(batchConfig, pdbDir, outDir, yamlDir, simInfo, topDir)
    elif cpuCount > 1:
        run_paralell(cpuCount, batchConfig, pdbDir, outDir, yamlDir, simInfo, topDir)

###################################################################################################### 
def run_serial(batchConfig, pdbDir, outDir, yamlDir, simInfo, topDir):
    for pdbFile in os.listdir(pdbDir):
        fileData = p.splitext(pdbFile)
        if not fileData[1] == ".pdb":
            continue  
        process_pdb_file(pdbFile, pdbDir, outDir, yamlDir, simInfo, topDir, batchConfig)
    ## CLEAN UP
    clean_up_handler(batchConfig)
######################################################################################################
def run_paralell(cpuCount, batchConfig, pdbDir, outDir, yamlDir, simInfo, topDir):
    with ThreadPoolExecutor(max_workers=min(os.cpu_count(),cpuCount)) as executor:
        for pdbFile in os.listdir(pdbDir):
            fileData = p.splitext(pdbFile)
            if not fileData[1] == ".pdb":
                continue
            executor.submit(process_pdb_file, pdbFile, pdbDir, outDir, yamlDir, simInfo, topDir, batchConfig)
   ## CLEAN UP
    clean_up_handler(batchConfig)

######################################################################################################
def clean_up_handler(batchConfig):
    # READ FROM batchConfig
    simulationInfo = batchConfig["simulationInfo"]
    outDir = batchConfig["pathInfo"]["outputDir"]
    # RUN THROUGH OPTIONS IN cleanUpInfo
    if  not "cleanUpInfo" in batchConfig:
        return 
    cleanUpInfo = batchConfig["cleanUpInfo"]
    if "getEndpointPdbs" in cleanUpInfo:
        if cleanUpInfo["getEndpointPdbs"]:
            cleanup.get_endpoint_pdbs(simulationInfo, outDir,cleanUpInfo)
            if any(key in cleanUpInfo for key in ["removeWaters","removeIons"]):
                cleanup.remove_atoms_from_pdb(simulationInfo, cleanUpInfo, outDir)

######################################################################################################
def extract_info(pdbDf,pdbDir,protName,yamlDir,batchConfig): ## gets info from pdb file, writes a config file
    ## GET PROTEIN INFORMATION
    aminoAcids =   ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                    'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                    'SER', 'THR', 'TRP', 'TYR', 'VAL']
    protDf = pdbDf[pdbDf["RES_NAME"].isin(aminoAcids)]

    protH = False
    if (protDf["ELEMENT"] == "H").any():
        protH = True

    proteinInfo = {"nProteins":1,
                   "proteins":[{"proteinName":f"{protName}",
                                "protons":protH}]}   
    ## GET LIGAND INFORMATION 
    ligsDf = pdbDf[~pdbDf["RES_NAME"].isin(aminoAcids)]
    ligNames = ligsDf["RES_NAME"].unique().tolist()

    ## SKIP IF NOT LIGAND
    if len(ligNames) == 0:
        ligandInfo = None
        return proteinInfo, ligandInfo
    
    ## USE ligandInfo IF SUPPLIED IN BATCH CONFIG
    if "ligandInfo" in batchConfig:
        ligandInfo = batchConfig["ligandInfo"]
        return proteinInfo, ligandInfo
    ## CREATE ligandInfo AUTOMATICALLY (WORKS FOR SIMPLE LIGANDS)
    else:
        ligList = []
        for ligName in ligNames:
            ligDf = ligsDf[ligsDf["RES_NAME"] == ligName]
            # deal with protonation
            ligH = False
            if (ligDf["ELEMENT"] == "H").any():
                ligH = True

            # deal with mol2
            ligMol2 = p.join(pdbDir,f"{ligName}.mol2")
            mol2 = False
            if p.isfile(ligMol2):
                mol2 = True
            # deal with frcmod
            ligFrcmod = p.join(pdbDir,f"{ligName}.frcmod")
            frcmod = False
            if p.isfile(ligFrcmod):
                frcmod = True  
            # deal with charge
            charge = find_ligand_charge(ligDf,ligName,yamlDir,pH=7.4)
            # write to temporary dict, then to ligandInfo for configg
            tmpDict = {"ligandName":ligName,
                    "protons":   ligH,
                    "mol2":      mol2,
                    "toppar":    frcmod,
                    "charge":    charge}
            ligList.append(tmpDict)
        nLigands = len(ligList)
        ligandInfo = {"nLigands":nLigands,
                    "ligands":ligList}

        return proteinInfo, ligandInfo
######################################################################################################
main()