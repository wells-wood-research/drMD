## BASIC PYTHON LIBRARIES
import os
from os import path as p
import sys
import yaml
from glob import glob
import numpy as np
import inflect

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

## drMD MODULES
from ExaminationRoom import drLogger

##  CLEAN CODE
from typing import Dict, Callable, List, Tuple, Set, Union
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

##########################################################################################
def methods_writer_protocol(batchConfig: Dict, configDir: DirectoryPath, outDir: DirectoryPath) -> None:
    """
    Main protocol for automated methods section writer.
    Looks through both batch config and per-run config files for information

    Args:
        batchConfig (dict): The batch configuration dictionary.
        configDir (DirectoryPath): Path to the directory containing config files. 
        outDir (DirectoryPath): Path to the output directory.

    Returns:
        None
    """


    global inflecter
    inflecter = inflect.engine()

    ## find all config files
    configDicts: List[dict] = get_config_dicts(configDir)
    autoMethodsDir: DirectoryPath = p.join(outDir, "00_AutoMethods")
    ## make outDir if it doesn't exist
    os.makedirs(autoMethodsDir, exist_ok=True)
    ## create methods file
    methodsFile: FilePath = p.join(autoMethodsDir, "drMD_AutoMethods.md")

    ## log what we are doing
    drLogger.log_info(f"writing automatic methods section to {methodsFile}", True)

    ## write a header to the methods file
    with open(methodsFile, "w") as f:
        f.write("# Molecular Dynamics Protocol\n\n")
        f.write(f"**This methods section was automatically generated by drMD [Ref. {cite('drMD')}].**\n\n")
        f.write("This document contains all the information needed to recreate your simulations.\n\n")
        f.write("Feel free to use this as the basis for your methods section in your papers, thesis, etc.\n\n")


    ## write ligand-prep related methods
    write_ligand_parameterisation_methods(configDicts, methodsFile)
    ## write protein-prep related methods
    write_protein_preparation_methods(configDicts, methodsFile)
    ## write solvation related methods
    write_solvation_charge_balance_methods(batchConfig, configDicts, methodsFile)
    ## write forcefield related methods
    write_forecefields_methods(methodsFile)
    ## write simulation related methods
    simulationInfo: dict = configDicts[0]["simulationInfo"]
    write_simulation_methods(methodsFile, simulationInfo)
##########################################################################################

def get_config_dicts(configDir: DirectoryPath) -> List[Dict]:
    """
    Reads through config files in configDir and returns a list of config dicts.

    Args:
        configDir (DirectoryPath): Path to the directory containing config files.

    Returns:
        List[Dict]: List of config dicts.
    """
    ## list comp to create a list of config files
    configFiles: List[FilePath] = [p.join(configDir, file) for file in os.listdir(configDir) if p.splitext(file)[1] == ".yaml"]
    ## create a list of dicts from config files
    configDicts = []
    for configFile in configFiles:
        with open(configFile, 'r') as f:
            configDicts.append(yaml.safe_load(f))

    return configDicts
##########################################################################################
def write_ligand_parameterisation_methods(configDicts: List[dict], methodsFile: FilePath) -> None:
    """
    Looks through individual run config dicts and writes out the ligand parameterisation methods.

    Args:
        configDicts: (list)List of config dicts
        methodsFile: (FilePath) Path to methods file

    Returns:
        None
    """
    ## init some empty sets
    obabelProtonatedLigands: set = set()
    antechamberChargesLigands: set = set()
    parmchkParamsLigands: set = set()
    allLigandNames = set()
    ## loop through config dicts
    for config in configDicts:
        ## read ligand info | skip if not there
        ligandInfo = config.get("ligandInfo", False)
        if not ligandInfo:
            continue
        ## loop through ligands
        for ligand in ligandInfo:
            ## check if drMD has: Added protons | Calculated charges | created parameter files
            allLigandNames.add(ligand["ligandName"])
            if ligand["protons"]:
                obabelProtonatedLigands.add(ligand["ligandName"])
            if ligand["mol2"]:
                antechamberChargesLigands.add(ligand["ligandName"])
            if ligand["toppar"]:
                parmchkParamsLigands.add(ligand["ligandName"])

    ## if no ligand in ligandInfo, then return an empty string 
    if len(allLigandNames) == 0:
        return 
    
    ## open methods file and append to it
    with open(methodsFile, "a", encoding = "utf-8") as methods:
        ## when we have ligand, but none have been processed by drMD, write a warning to fill this in manually
        methods.write("## Ligand Preparation and Parameterisation\n\n")
        if len(obabelProtonatedLigands) > 0 and len(antechamberChargesLigands) > 0 and len(parmchkParamsLigands) > 0:
            methods.write(f"\n\n**WARNING** drMD did not run any automated procedures to parameterise your ligand(s).")
            methods.write("You will need to fill in this section manually.\n\n")
            return
        ## generic header
        methods.write(f"Ligand parameter generation was performed using the following procedure: ")
        
        ## if all ligands were fully processed by drMD
        if obabelProtonatedLigands == antechamberChargesLigands == parmchkParamsLigands:
            methods.write(f"All ligand were protonated using OpenBabel [Ref. {cite('obabel')}], ")
            methods.write(f"newly added hydrogen atoms were then renamed to ensure compatability with the AMBER forcefeild. ")
            methods.write(f"Partial charges of all ligand were calculated, and atom types were assigned using antechamber [Ref. {cite('antechamber')}], ")
            methods.write(f"and the parameters for the ligand were generated using parmchk [Ref. {cite('parmchk')}].")
        ## in odd cases where drMD hasn't fully processed ligands
        else:
            if len(obabelProtonatedLigands) > 0:
                ligand = format_list(list(obabelProtonatedLigands))
                methods.write(f"The ligand {ligand} were protonated using OpenBabel [Ref. {cite('obabel')}], novel hydrogen atoms were then renamed to ensure compatability with the AMBER forcefeild. ")
            if len(antechamberChargesLigands) > 0:
                ligand = format_list(list(antechamberChargesLigands))
                methods.write(f"For the ligand {ligand} partial charges were calculated and atom types were assigned using antechamber [Ref. {cite('antechamber')}]. ")
            if len(parmchkParamsLigands) > 0:
                ligand = format_list(list(parmchkParamsLigands))
                methods.write(f"The parameters for the ligand {ligand} were generated using parmchk [Ref. {cite('parmchk')}]. ")
        ## line break at the end
        methods.write("\n\n")
##########################################################################################


def write_protein_preparation_methods(configDicts: List[dict], methodsFile: FilePath) -> None:
    """
    Looks through per-run config dicts and writes out the protein preparation methods.

    Args:
        configDicts: (list) List of config dicts
        methodsFile: (FilePath) Path to methods file

    Returns:
        None
    """
    ## work out which proteins were protonated, collect in a dict
    proteinsProtonated = {}
    for config in configDicts:
        proteinInfo = config.get("proteinInfo", False)
        isProteinProtonated = proteinInfo.get("protons", False)
        inputPdbName = p.splitext(p.basename(config["pathInfo"]["inputPdb"]))[0]
        proteinsProtonated[inputPdbName] = isProteinProtonated


    ## get the pH of simulations
    pH = configDicts[0]["miscInfo"]["pH"]

    ## open methods file and append to it
    with open (methodsFile, "a") as methods:
        methods.write("## Protein Preparation\n\n")
        ## if all of the proteins were protonated before drMD, write a warning
        if  all(proteinsProtonated.values()):
            methods.write(f"\n\n**WARNING** drMD did not run any automated procedures to prepare your proteins. ")
            methods.write(f"You will need to fill in this section manually. ")
            methods.write("\n\n")
        ## if none of the proteins were protonated before drMD, write automated protonation procedure
        elif not any(proteinsProtonated.values()):
            methods.write(f"\nAll proteins were protonated using software pdb2pqr [Ref. {cite('pdb2pqr')}] ")
            methods.write(f"which uses ProPKA to calculate per-residue proton affinities [Ref. {cite('propka')}]. ")
            methods.write(f"Proteins were protonated using the following pH: {pH}. ")
            methods.write("This process also automatically creates disulfide bonds as appropriate. ")
        ## if some, but not all of the proteins were protonated before drMD, write automated protonation procedure
        else:
            protonatedProteinNames = [protName for protName in proteinsProtonated if proteinsProtonated[protName]]
            protonatedText = format_list(protonatedProteinNames)
            nonProtonatedProteinNames = [protName for protName in proteinsProtonated if not proteinsProtonated[protName]]
            nonProtonatedText = format_list(nonProtonatedProteinNames)

            methods.write(f"\nThe proteins {nonProtonatedText} were protonated using software pdb2pqr [Ref. {cite('pdb2pqr')}] ")
            methods.write(f"which uses ProPKA to calculate per-residue proton affinities [Ref. {cite('propka')}]. ")
            methods.write(f"Proteins were protonated using the following pH: {pH}. ")
            methods.write("This process also automatically creates disulfide bonds as appropriate. ")
            methods.write(f"\n\n**WARNING** drMD did not protonate proteins {protonatedText}. ")
            methods.write("You will need to fill in this section manually. ")
        

        ## line break at the end of this section
        methods.write("\n")
##########################################################################################
def count_waters(pdbFile: FilePath) -> int:
    """
    Helper function that counts water molecules in a PDB file.

    Args:
        pdbFile: (FilePath) Path to PDB file

    Returns:
        waterCount: (int) Number of water molecules

    """
    ## load pdb as a DataFrame
    pdbDf = pdbUtils.pdb2df(pdbFile)
    ## count water atoms
    waterDf = pdbDf[pdbDf["RES_NAME"] == "WAT"]
    ## return water molecule count
    return len(waterDf) / 3

##########################################################################################
def write_solvation_charge_balance_methods(batchConfig: Dict, configDicts: List[dict], methodsFile: FilePath) -> None:
    """
    Reads through per-run config dicts and writes out the solvation and charge balance methods.
    Looks through prep files to see how many water molecules and ions were added.
    
    Args:
        batchConfig (dict): The batch configuration dictionary.
        configDicts: (list) List of config dicts
        methodsFile: (FilePath) Path to methods file

    Returns:
        None
    """
    ## get box geometry
    boxGeometry = configDicts[0]["miscInfo"]["boxGeometry"]
    ## get an approximate count of water molecules added to the box
    approxWaterCount, counterIonCounts = get_solvation_atom_counts(configDicts, batchConfig)
    ## write solvation and charge balence methods
    with open(methodsFile, "a", encoding = "utf-8") as methods:
        methods.write("## Solvation and Charge Balencing\n\n")
        ## info on solvation box
        methods.write(f"All proteins were placed in {inflect.engine().an(boxGeometry)} {boxGeometry} solvation box with a 10 Å buffer between")
        methods.write(f" the protein and the nearest edge of the box. ")
        methods.write("The system was treated using periodic boundary conditions. ")
        ## average water count
        methods.write(f"Approximately {approxWaterCount} TIP3P water molecules were added to the solvation box. ")
        ## info on counter ions
        methods.write(f"Sodium and Chloride ions were added to the box to balance the charge of the system.\n")
        methods.write(f"A table showing the counts of counter ions is provided below:\n\n")
        methods.write("|\t Protein Name\t| \tSodium Ions\t| \tChloride Ions\t|\n")
        methods.write("|\t------------\t| \t------------\t| \t------------\t|\n")
        for protName in counterIonCounts:
            methods.write(f"|\t{protName}|{counterIonCounts[protName]['Sodium']}| {counterIonCounts[protName]['Chloride']}|\n")
        ## line break at the end of this section
        methods.write("\n\n")
################################################################################
def get_solvation_atom_counts(configDicts: List[dict], batchConfig: Dict) -> Tuple[int,Dict[str,int]]:
    """
    Counts counter-ions added to the solvation box.

    Args:   
        configDicts: (list) List of config dicts
        batchConfig (dict): The batch configuration dictionary.
    
    Returns:
        approxWaterCount: (int) Number of water molecules added to the solvation box
        counterIonCounts: (Dict[str,int]) Dictionary of counter-ions counts per protien
    """
    ## init empty list and dict to store data
    waterCounts = []
    counterIonCounts = {}
    ## find outDir
    outDir = batchConfig["pathInfo"]["outputDir"]
    ## for each config file...
    for config in configDicts:
        ## find a the solvated pdb file
        protName = config["proteinInfo"]["proteinName"]
        inputPdbName = p.splitext(p.basename(config["pathInfo"]["inputPdb"]))[0]
        prepDir = p.join(outDir, inputPdbName, "00_prep")

        ligandInfo = config.get("ligandInfo", False)
        if ligandInfo:
            solvationDir = p.join(prepDir,"WHOLE")
        else:
            solvationDir = p.join(prepDir,"PROT")
        solvatedPdb = [p.join(solvationDir, file) for file in os.listdir(solvationDir) if file.endswith("solvated.pdb")][0]

        ## count water molecules | append to list
        waterCount = count_waters(solvatedPdb)
        waterCounts.append(waterCount)

        ## count counter-ions | update dict
        nNa, nCl = count_ions(solvatedPdb)
        counterIonCounts[protName] = {"Sodium": nNa, "Chloride": nCl}

    ## average water count and round
    waterCounts = np.array(waterCounts)
    averageWaterCount = np.mean(waterCounts)
    approximateAverageWaterCount = int(round(averageWaterCount, 2 - len(str(int(averageWaterCount)))))

    return approximateAverageWaterCount, counterIonCounts
##########################################################################################

def count_ions(pdbFile: FilePath) -> Tuple[int, int]:
    """
    Helper function that counts counter-ions in a PDB file.

    Args:
        pdbFile: (FilePath) Path to PDB file

    Returns:
        naCount: (int) Number of sodium ions
        clCount: (int) Number of chloride ions
    """
    ## load pdb as a DataFrame
    pdbDf = pdbUtils.pdb2df(pdbFile)
    ## count counter-ions
    naDf = pdbDf[pdbDf["RES_NAME"] == "Na+"]
    clDf = pdbDf[pdbDf["RES_NAME"] == "Cl-"]

    return len(naDf) , len(clDf)

##########################################################################################
def format_list(inputList: List[str]) -> str:
    """
    Helper function that formats a list of strings into:
    "item1, item2, and item3"

    Args:
        inputList: (List[str]) List of strings
    
    Returns:
        formattedList: (str) Formatted string
    """
    ## if there is only one item, unpack into a str and return
    if len(inputList) == 1:
        return inputList[0]
    ## if there are 2 or more items, return a comma separated list
    else:
        return ", ".join(inputList[:-1]) + ", and " + inputList[-1]
##########################################################################################

def get_progression_word(stepIndex: int, maxSteps: int) -> str:
    """
    Helper function that returns the progression word for a simulation step.

    Args:
        stepIndex: (int) Index of the current simulation step
        maxSteps: (int) Total number of simulation steps
    
    Returns:
        word: (str) The progression word
    """
    if stepIndex == 0 and maxSteps == 1:
        return ""
    if stepIndex == 0:
        return "Initially, "
    elif stepIndex == maxSteps - 1:
        return "Finally, "
    else:
        return "Next, "
##########################################################################################
def get_simulation_type_text(sim: Dict, progressionWord: str) -> str:
    """
    Helper Function that converts simulation type to methods text.

    Args:
        sim: (dict) Simulation dictionary

    Returns:
        text: (str) Methods text
    """
    capitalise = False
    if len(progressionWord) == 0:
        capitalise = True

    ## read simulation type from simulation dictionary | choose appropriate methods text
    simulationType = sim["simulationType"]
    if simulationType.upper() == "NPT":
        article = "A" if capitalise else "a"
        return f"{progressionWord}{article} simulation was performed using the *isothermal-isobaric* (NpT) ensemble"
    elif simulationType.upper() == "NVT":
        article = "A" if capitalise else "a"
        return f"{progressionWord}{article} simulation was performed using the canonical (NVT) ensemble"
    elif simulationType.upper() == "EM":
        article = "An" if capitalise else "an"
        return f"{progressionWord}{article} energy minimisation step was performed using the steepest descent method"
    elif simulationType == "META":
        article = "A" if capitalise else "a"
        return f"{progressionWord}{article} metadynamics simulation was performed"
    
##########################################################################################
def get_restraints_methods_text(sim: Dict) -> str:
    """
    Generates the text for the restraints methods.

    Args:
        sim: (dict) Simulation dictionary

    Returns:
        text: (str) Methods section for restraints
    
    """
    ## if there are no restraints, return an empty string
    if not "restraintInfo" in sim:
        return ""
    ## if there are restraints, generate the text
    restraintInfo = sim["restraintInfo"]
    text = ""
    for restraint in restraintInfo:
        text += f"{inflecter.a(restraint['restraintType']).capitalize()} restraint"
        text += f" with a force constant of {restraint['parameters']['k']} {get_force_constant_units(restraint['restraintType'])} "
        text += f" {get_restraint_target(restraint)} "
        text += f"was applied to {selection_to_text(restraint['selection'])}. "
    return text

##########################################################################################
def get_force_constant_units(restraintType: str) -> str:
    """
    Gets appropriate force constant units for each restraint type.

    Args:
        restraintType: (str) The type of restraint

    Returns:
        forceConstantUnits: (str) The force constant units
    """

    if restraintType == "position":
        return "kJ mol<sup>-1</sup> nm<sup>-2</sup>"
    elif restraintType == "distance":
        return "kJ mol<sup>-1</sup> nm<sup>-2</sup>"
    elif restraintType == "angle":
        return "kJ mol<sup>-1</sup> rad<sup>-2</sup>"
    elif restraintType == "torsion":
        return "kJ mol<sup>-1</sup> rad<sup>-2</sup>"
##########################################################################################
def selection_to_text(selection: Dict) -> str:
    """
    Generates the text for an atom selection.

    Args:
        selection: (dict) Selection dictionary

    Returns:
        text: (str) Selection text
    """

    ## find the keyword for the selection
    keyword = selection["keyword"]
    ## for keyword selections, we can just return the keyword as part of the text
    if not keyword == "custom":
        if keyword == "all":
            return "all atoms in the system"
        else:
            return f"all {keyword} atoms in the system"  

    ## for custom selections, we need to generate the text
    text = "the following atoms: "
    ## init an empty list of selection texts
    selectionTexts = []
    ## get the selection list
    customSelections: List[dict] = selection["customSelection"]
    for customSelection in customSelections:
        selectionText = ""
        ## deal with atoms
        atomName = customSelection["ATOM_NAME"]
        if  atomName == "all":
            selectionText += "all atoms in"
        else:
            selectionText += f"atom{identifier_list_to_str(atomName)}"

        ## deal with resId and resName together
        residueId = customSelection["RES_ID"]
        residueName = customSelection["RES_NAME"]
        if not residueId == "all" and not residueName == "all":
            if isinstance(residueId, str) and isinstance(residueName, str):
                selectionText += f" in residue {residueName}{residueId}"
            else:
                selectionText += f" in residue {identifier_list_to_str(residueId)}{identifier_list_to_str(residueName)}"
        elif not residueId == "all":
            selectionText += f" in residue {identifier_list_to_str(residueId)}"
        elif not residueName == "all":
            selectionText += f" in residue {identifier_list_to_str(residueName)}"
        ## deal with chain
        chainId = customSelection["CHAIN_ID"]
        if not chainId == "all":
            selectionText += f" in chain{identifier_list_to_str(chainId)}"

        selectionTexts.append(selectionText)


    selectionTexts = format_list(selectionTexts)
    text += selectionTexts

    return text

##########################################################################################
def get_restraint_target(restraint: str) -> str:
    """
    Helper function that returns the target of a restraint.

    Args:
        restraint: (dict) Restraint dictionary

    Returns:
        text: (str) methods text for the restraint target
    """

    restraintType: str = restraint["restraintType"]

    if restraintType == "position":
        return " "
    elif restraintType == "distance":
        return f" and an equilibrium distance of {restraint['parameters']['r0']} Å"
    elif restraintType == "angle":
        return f" and an equilibrium angle of {restraint['parameters']['theta0']} degrees"    
    elif restraintType == "torsion":
        return f" and an equilibrium dihedral angle of {restraint['parameters']['phi0']} degrees"

##########################################################################################
def identifier_list_to_str(identifier: Union[str, list]) -> str:
    """
    Helper function that puts an 's' for a list or nothing for a string
    
    Args:
        identifier: (Union[str, list]) List or string

    Returns:
        text: (str) methods text for the identifier
    """
    if isinstance(identifier, str):
        return " " + identifier
    else:
        return f"s {format_list(identifier)}"

##########################################################################################
def write_per_step_simulation_methods(methodsFile: FilePath, sim: dict, stepIndex: int, maxSteps: int) -> None:
    """
    Writes the methods for a per step simulation

    Args:
        methodsFile: (FilePath) Path to methods file
        sim: (dict) Simulation dictionary
        stepIndex: (int) Index of the current step
        maxSteps: (int) Maximum number of steps

    Returns:
        None
    """

    with open(methodsFile, "a", encoding = "utf-8") as methods:
        ## progression word
        progressionWord = get_progression_word(stepIndex, maxSteps)
        ## simulation type
        methods.write(f"{get_simulation_type_text(sim, progressionWord)}.\n")
        ## deal with EM and maxIterations
        if sim["simulationType"] == "EM":
            if sim["maxIterations"] == -1:
                methods.write(f"This energy minimisation step was performed until it reached convergence. ")
            else:
                methods.write(f"This energy minimisation step was performed for {sim['maxIterations']} steps, ")
                methods.write(f"or until it reached convergence. ")
        ## deal with NPT, NVT, META
        else:
            methods.write(f"This simulation was performed for {sim['duration']} ") 
            if "temperature" in sim:
                methods.write(f" at {sim['temperature']} K. ") 
            elif "temperatureRange" in sim:
                methods.write(f"The temperature of this simulation was stepped through the range ")
                methods.write(f"{format_list([str(temp) + ' K' for temp in sim['temperatureRange']])} in even time increments. ")

            ## deal with heavy protons // timestep
            heavyProtons = sim.get("heavyProtons", False)
            if heavyProtons:
                methods.write(f"This simulations was performed using a mass of 4.03036 amu for hydrogen atoms. ")
                methods.write(f"The mass added to each hydrogen atom was subtracted from the mass of the heavy atom it was bonded to. ")
                methods.write(f"This, combined with the constraints placed upon bonds between heavy and hydrogen atoms, ")
                methods.write(f"allowed the simulation to be performed using a timestep of {sim['timestep']} [Ref. {cite('heavyProtons')}]. ")
            else:
                methods.write(f"This simulation was performed using a timestep of {sim['timestep']}. ")

        ## deal with metadynamics   
        if sim["simulationType"] == "META":
            write_metadynamics_simulation_methods(methodsFile, sim)

        ## deal with restraints
        methods.write(f"{get_restraints_methods_text(sim)}\n")
        ## line break at the end of the section
        methods.write("\n\n")
##########################################################################################
def write_metadynamics_simulation_methods(methodsFile: FilePath, sim: dict) -> None:
    """
    Writes the methods for a metadynamics simulation

    Args:
        methodsFile: (FilePath) Path to methods file
        sim: (dict) Simulation dictionary

    Returns:
        None
    """
    metaDynamicsInfo = sim["metaDynamicsInfo"]
    biases = metaDynamicsInfo["biases"]
    with open(methodsFile, "a", encoding = "utf-8") as methods:
        methods.write(f"This metadynamics simulation was performed using the well-tempered metadynamics method [Ref. {cite('metadynamics')}]. ")
        methods.write(f"This simulation was performed using a height parameter of {metaDynamicsInfo['height']}, ")
        methods.write(f"a biasFactor parameter of {metaDynamicsInfo['biasFactor']}, ")
        methods.write(f"a frequency parameter of {metaDynamicsInfo['frequency']}. ")

        if len(biases) > 1:
            methods.write(f"The following bias variables were used in this simulation: ")
        for bias in biases:
            if bias["biasVar"].upper() in ["RMSD", "DISTANCE"]:
                unit = "Å"
            elif bias["biasVar"].upper() in ["ANGLE", "TORSION"]:
                unit = "°"
  
            biasSelection = bias["selection"]
            methods.write(f"This simulation used a {bias['biasVar']} bias variable ")
            methods.write(f"with a maximum value of {bias['maxValue']} {unit } and ")
            methods.write(f"a minimum value of {bias['minValue']} {unit}. ")
            methods.write(f"Gaussians with a width parameter of {bias['biasWidth']} {unit} were used to perturb the bias variable. ")
            methods.write(f"This bias variable was applied to {selection_to_text(biasSelection)}. ")




##########################################################################################
def write_generic_simulation_methods(methodsFile: FilePath, simulationInfo: dict) -> None:
    """
    Writes methods section generic to all simulations run by drMD

    Args:
        methodsFile: (FilePath) Path to methods file

    Returns:
        None
    """

    with open(methodsFile, "a", encoding = "utf-8") as methods:
        methods.write(f"All simulations were performed using the OpenMM simulation toolkit [Ref. {cite('openmm')}]. ")
        ## Check for NVT and/or NPT
        simluationTypes = [step["simulationType"].upper() for step in simulationInfo]
        if "NVT" in simluationTypes or "NPT" in simluationTypes:
            ## LangevinMiddleIntegrator
            methods.write(f"All simulations were performed using the Langevin Middle Integrator [Ref. {cite('langevinMiddleIntegrator')}] ")
            methods.write("which was used to enforce constant temperature conditions in each simulation. ")
        if "NPT" in simluationTypes:
            ## MonteCarloBarostat
            methods.write("For simulations run under the *isothermal-isobaric* (NpT) ensemble, ")
            methods.write("the Monte-Carlo barostat was used to enforce a constant pressure of 1 atm. ")
        ## ParticleMeshEwald
        methods.write("In all simulations, long-range Coulombic interactions were modelled using the ")
        methods.write(f"Particle-Mesh Ewald (PME) method [Ref. {cite('pme')}], with a 10 Å cutoff distance. ")
        ## HBond constraints
        methods.write(f"In all simulations, constraints were applied to bonds between hydrogen atoms and heavy atoms using the SHAKE algorithm [Ref. {cite('shake')}]. ")
        ## water constraints
        methods.write(f"In all simulations, bonds lengths and angles of water molecules were constrained using the SETTLE algorithm [Ref. {cite('settle')}]. ")
        
##########################################################################################
def write_forecefields_methods(methodsFile):
    """
    Writes methods section talking about forcefields, generic to all simulations run by drMD
    
    Args:
        methodsFile: (FilePath) Path to methods file

    Returns:
        None    
    """
    with open(methodsFile, "a", encoding = "utf-8") as methods:
        methods.write("## Forcefield Information\n\n")
        methods.write(f"All protein residues were parameterised using the AMBER ff19SB and forcefield [Ref. {cite('ff19SB')}]. ")
        methods.write(f" These parameters were prepared using tleap from the Ambertools package [Ref. {cite('ambertools')}]. ")
        methods.write(f"Simulations were performed in explicit solvent. All water molecules parameterised using the TIP3P model [Ref. {cite('tip3pParams')}]. ")
        methods.write(f"Any ions in our system were treated using parameteres calculated to complement the TIP3P water model [Ref. {cite('ionParams')}]. \n\n")
##########################################################################################
def  write_simulation_methods(methodsFile: FilePath, simulationInfo: dict):
    """
    Protocol function for writing simulation methods to methods file

    """
    ## write a title for the section
    with open(methodsFile, "a", encoding = "utf-8") as methods:
        methods.write("## Simulation Protocols\n\n")
    ## write the per-step simulation methods
    for stepIndex, sim in enumerate(simulationInfo):
        write_per_step_simulation_methods(methodsFile, sim, stepIndex, len(simulationInfo))
    ## write the generic simulation methods
    write_generic_simulation_methods(methodsFile, simulationInfo)
##########################################################################################
def cite(key: str) -> str:
    """
    Used to add citations to methods text

    Args:
        key: (str) Key of the citation to be added

    Returns:
        (str) Formatted Citation
    """

    doiDict = {
        ## drMD citations
        "drMD": ["PLACEHOLDER"],
        ## openmm citation
        "openmm": ["10.1371/journal.pcbi.1005659"],
        ## prep step citations
        "pdb2pqr": ["10.1093/nar/gkm276", "10.1093/nar/gkh381"],
        "propka": ["10.1021/ct200133y", "10.1021/ct100578z"],
        "obabel": ["10.1186/1758-2946-3-33"],
        "antechamber": ["10.1002/jcc.20035", "10.1016/j.jmgm.2005.12.005"],
        "parmchk": ["10.1002/jcc.20035", "10.1016/j.jmgm.2005.12.005"],
        "ambertools" :["10.1021/acs.jcim.3c01153"],
        ## simulation step citations
        "langevinMiddleIntegrator": ["10.1021/acs.jpca.9b02771"],
        "heavyProtons": ["10.1002/(SICI)1096-987X(199906)20:8<786::AID-JCC5>3.0.CO;2-B"],
        "shake": ["10.1002/1096-987X(20010415)22:5<501::AID-JCC1021>3.0.CO;2-V"],
        "settle": ["10.1002/jcc.540130805"],
        "pme":["10.1063/1.464397"],
        "metadynamics": ["10.1103/PhysRevLett.100.020603"],
        ## parameter citations
        "ionParams": ["10.1021/ct500918t", "10.1021/ct400146w"],
        "tip3pParams": ["/10.1063/1.472061"],
        "ff19SB": ["10.1021/acs.jctc.9b00591"]
    }


    citation = doiDict.get(key, "CITATION NOT FOUND")

    if len(citation) > 1:
        citation = [f"[{key}({index+1})](https://doi.org/{ref})" for index, ref in enumerate(citation)]
    else:
        citation = [f"[{key}](https://doi.org/{citation[0]})"]

    return format_list(citation)
##########################################################################################

## FOR TESTING
if __name__ == "__main__":
    configDir = "/home/esp/scriptDevelopment/drMD/04_PET_proj_outputs/00_configs"
    outDir = "/home/esp/scriptDevelopment/drMD/04_PET_proj_outputs/00_methods"
    batchConfigYaml = "/home/esp/scriptDevelopment/drMD/prescriptions/PETase_MD_config.yaml"
    methods_writer_protocol(batchConfigYaml, configDir, outDir)
