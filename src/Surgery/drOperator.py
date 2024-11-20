## BASIC PYTHON LIBRARIES
import os
from os import path as p
from shutil import rmtree

## OPENMM LIBRARIES
import openmm.app as app
import openmm as openmm

## drMD LIBRARIES
from Surgery import drPrep, drSim, drMeta
from Triage import drConfigTriage
from ExaminationRoom import drLogger

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

##  CLEAN CODE
from typing import Dict, Callable
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath
#####################################################################################
def drMD_protocol(configYaml: FilePath) -> None:
    """
    Run the drMD protocol for a given configuration file.

    Args:
        configYaml (str): The path to the YAML configuration file.

    Returns:
        None

    This function reads the configuration file, prepares the protocol, and runs the simulation.
    """
    # Read the configuration file
    config: dict = drConfigTriage.read_config(configYaml)

    # Create the output directory if it doesn't exist
    outDir: str = config["pathInfo"]["outputDir"]
    os.makedirs(outDir, exist_ok=True)
    # Prepare the protocol
    solvatedPdb, inputCoords, amberParams = drPrep.prep_protocol(config)

    # Run the simulation
    run_simulation(config = config,
                    outDir = outDir,
                      inputCoords = inputCoords,
                        amberParams = amberParams,
                          pdbFile = solvatedPdb)
###########################################################################################
def run_simulation(config: dict, outDir: str, inputCoords: str, amberParams: str, pdbFile: str) -> None:
    """
    Run the simulation according to the given configuration.

    Args:
        config (dict): The configuration dictionary.
        outDir (str): The output directory.
        inputCoords (str): The path to the input coordinates file.
        amberParams (str): The path to the Amber parameters file.
        pdbFile (str): The path to the PDB file.

    Returns:
        None
    """
    ## set up logging
    logDir = p.join(p.dirname(outDir), "00_drMD_logs")
    protName = config["proteinInfo"]["proteinName"]
    drLogger.setup_logging(p.join(logDir,f"{protName}_simulations.log"))

    platform = choose_platform(config)

    # Load Amber files and create system
    prmtop: app.Topology = app.AmberPrmtopFile(amberParams)
    inpcrd: app.InpcrdFile = app.AmberInpcrdFile(inputCoords)

    # Loop over simulations
    simulations = config["simulationInfo"]
    for i in range(len(simulations)):
        sim: dict = simulations[i]
        simDir: str = p.join(outDir,sim["stepName"])

        saveFile = None
        # Decide whether to skip, resume, or start a new simulation
        skipResumeSim, foundSaveFile = skip_resume_or_simulate(simDir=simDir,
                                                           simulations = simulations,
                                                           i = i, 
                                                           outDir=outDir)

        if not foundSaveFile == None:
            saveFile = foundSaveFile

        # Skip or resume simulation
        if skipResumeSim == "skip":
            stepName: str = sim["stepName"]
            drLogger.log_info(f"Skipping {stepName} for run: {protName}", True)
            continue
        if skipResumeSim == "resume":
            drLogger.log_info(f"Resuming {stepName} from checkpoint file for run: {protName}", True)
            rename_out_files(simDir)    

        # Run simulation
        simulationFunction = choose_simulation_function(sim["simulationType"])
        try:
            saveFile = simulationFunction(prmtop = prmtop,
                                       inpcrd = inpcrd,
                                         sim = sim,
                                           saveFile = saveFile,
                                             outDir = outDir,
                                               platform = platform,
                                                 refPdb = pdbFile,
                                                   config = config)
        except Exception as e:
            drLogger.log_info(f"Error running simulation: {e}", True, True)
            continue


###########################################################################################
def choose_platform(config: Dict) -> openmm.Platform:
    # Set up platform
    usePlatform: str = config["hardwareInfo"]["platform"]
    if usePlatform == "CUDA":
        platform=openmm.Platform.getPlatformByName("CUDA")
    elif usePlatform == "OpenCL":
        platform=openmm.Platform.getPlatformByName("OpenCL")
    elif usePlatform == "CPU":
        platform=openmm.Platform.getPlatformByName("CPU")

    return platform
###########################################################################################
def rename_out_files(simDir: DirectoryPath) -> None:
    progressReportCsv = p.join(simDir, "progress_report.csv")
    vitalsReportCsv = p.join(simDir, "vitals_report.csv")
    trajectoryDcd = p.join(simDir, "trajectory.dcd")

    allTrajectoryDcds = [p.join(simDir, file) for file in os.listdir(simDir) if file.endswith(".dcd")]
    resumeNumber = len(allTrajectoryDcds)

    if p.isfile(progressReportCsv):
        os.rename(progressReportCsv, p.join(simDir, f"progress_report_partial_{str(resumeNumber)}.csv"))
    if p.isfile(vitalsReportCsv):
        os.rename(vitalsReportCsv, p.join(simDir, f"vitals_report_partial_{str(resumeNumber)}.csv"))
    if p.isfile(trajectoryDcd):
        os.rename(trajectoryDcd, p.join(simDir, f"trajectory_partial_{str(resumeNumber)}.dcd"))
    
###########################################################################################
def choose_simulation_function(simulationType: str) -> Callable:
    """
    Choose the appropriate simulation function based on the simulation type.

    Args:
        simulationType (str): The simulation type.

    Returns:
        Callable: The appropriate simulation function.

    This function chooses the appropriate simulation function based on the simulation type.
    """
    if simulationType.upper() == "EM":
        return drSim.run_energy_minimisation
    elif simulationType.upper() in ["NPT","NVT"]:
        return drSim.run_molecular_dynamics
    elif simulationType.upper() == "META":
        return drMeta.run_metadynamics
###########################################################################################
def skip_resume_or_simulate(simDir: str, simulations: list, i: int, outDir: str) -> tuple:
    """
    Check if the simulation directory exists and decide whether to skip, resume or start a new simulation.

    Args:
        simDir (str): The path to the simulation directory.
        simulations (list): List of simulation information dictionaries.
        i (int): Index of the current simulation in the list.
        outDir (str): The path to the output directory.

    Returns:
        tuple: A tuple containing the action to be taken (simulate, skip or resume) and the path to the save file.
    """


    ## if the simDir for this step doesn't exist, find the xml file for the previous step
    if not p.isdir(simDir):
        ## if this is the first simulation in the series and the simDir doesn't exist, run the step from scratch
        if i == 0:
            return "simulate", None

        previousSimName = simulations[i-1]["stepName"]
        previousSimDir = p.join(outDir, previousSimName) if i > 0 else False
        saveXml = p.join(previousSimDir, f"{previousSimName}.xml") if previousSimDir else False
        return "simulate", saveXml
    
    ## if the simDir already exists:
    ## 1. look for XML file in simDir, if found *skip* and return saveFile
    ## 2. look for CHK file in simDir, if found *resume* and return saveFile
    ## 3. if neither exists, delete the simDir -> find the last saveFile -> run the step from scratch
    else:
        thisSimName = simulations[i]["stepName"]
        saveXml = p.join(simDir, f"{thisSimName}.xml")  
        if p.isfile(saveXml):  
            return "skip", saveXml
        saveChk = p.join(simDir, "checkpoint.chk")
        if p.isfile(saveChk):
            return "resume", saveChk
        else:
            # rmtree(simDir)
            previousSimName = simulations[i-1]["stepName"]
            previousSimDir = p.join(outDir, previousSimName) if i > 0 else False
            saveXml = p.join(previousSimDir, f"{previousSimName}.xml") if previousSimDir else False
            return "simulate", saveXml


