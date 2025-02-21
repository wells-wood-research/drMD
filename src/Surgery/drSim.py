## BASIC PYTHON LIBRARIES
import os
from os import path as p
import numpy as np
import sys

## OPENMM LIBRARIES
import openmm.app as app
import openmm as openmm
import  openmm.unit  as unit
## MDTRAJ LIBRARIES // DISABLE WARNINGS
import mdtraj 
import warnings
from mdtraj.utils.validation import TypeCastPerformanceWarning
warnings.filterwarnings("ignore", category=TypeCastPerformanceWarning)

## drMD LIBRARIES
from Surgery import drRestraints, drFirstAid
from ExaminationRoom import drLogger, drCheckup
from UtilitiesCloset import drFixer, drSelector

## CLEAN CODE
from typing import Optional, Dict, List, Tuple, Union, Any
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

###########################################################################################
def initialise_simulation(prmtop: app.AmberPrmtopFile,
                           inpcrd: app.AmberInpcrdFile,
                             sim: Dict,
                               saveFile: FilePath,
                                 refPdb: FilePath,
                                   platform: str,
                                     hardwareInfo: Dict) -> openmm.System:
    """
    Create an openmm.System from a prmtop.

    Parameters
    ----------
    prmtop : app.Topology
        The topology of the system.

    Returns
    -------
    system : openmm.System
        The initialized system.

    """

    # Define the nonbonded method and cutoff.
    nonbondedMethod: openmm.NonbondedForce = app.PME
    nonbondedCutoff: unit.Quantity = 1 * unit.nanometer

    ## deal with heavy protons
    heavyProtons = sim.get("heavyProtons", False)
    if heavyProtons:
        protonMass = 1.00784 * 4 * unit.amu
        ## check to see if a large timestep has been used
        timestep = sim["timestep"].value_in_unit(unit.femtoseconds)
        if timestep <= 2:
            drLogger.log_info("WARNING: Timestep is too small for heavy protons", True, True)
            drLogger.log_info("WARNING: Use at least 2 fs to get a faster simulation", True, True)

    else:
        protonMass = 1.00784 * unit.amu
 

    # Create the system.
    system: openmm.System = prmtop.createSystem(nonbondedMethod=nonbondedMethod,
                                                nonbondedCutoff=nonbondedCutoff,
                                                constraints=app.HBonds,
                                                hydrogenMass=protonMass)
    
    ## use static temperature if specified, or use first value in temperatureRange to start with
    if "temperature" in sim:
        initailSimulationTemp = sim["temperature"]
    elif "temperatureRange" in sim:
        initailSimulationTemp = sim["temperatureRange"][0]

    ## deal with any restraints
    system: openmm.System = drRestraints.restraints_handler(system, prmtop, inpcrd, sim, saveFile, refPdb)
    # add constant pressure force to system (makes this an NpT simulation)    
    if sim["simulationType"].upper() in ["NPT","META"]:
        if "temperature" in sim:
            system.addForce(openmm.MonteCarloBarostat(1*unit.bar, initailSimulationTemp))
        elif "temperatureRange" in sim:
            system.addForce(openmm.MonteCarloBarostat(1*unit.bar, initailSimulationTemp))

    ## setup an intergrator
    if sim["simulationType"].upper() == "EM":
        integrator: openmm.Integrator = openmm.LangevinMiddleIntegrator(initailSimulationTemp,
                                                                  1/unit.picosecond,
                                                                  4*unit.femtosecond)
    else:
        integrator: openmm.Integrator = openmm.LangevinMiddleIntegrator(initailSimulationTemp,
                                                                         1/unit.picosecond,
                                                                           sim["timestep"])


    simulation: app.simulation.Simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    

    return simulation , integrator

###########################################################################################
def process_sim_data(sim: Dict) -> Dict:
    """
    Reads simulation variables from config file and processes them.

    Parameters
    ----------
    sim : Dict
        Dictionary containing simulation variables.

    Returns
    -------
    sim : Dict
        Updated dictionary with processed simulation variables.
    """

    if "processed" in sim:
        return sim

    # Set up unit translator
    timescale: Dict = {"fs":unit.femtoseconds,
                 "ps":unit.picoseconds,
                 "ns":unit.nanoseconds}
    
    # Read timestep data and process it
    timestepData = sim["timestep"].split()
    timestep: unit.Quantity = float(timestepData[0]) * timescale[timestepData[1]]
    sim["timestep"] = timestep

    # Read duration data and process it
    durationData = sim["duration"].split()
    duration: int = int(durationData[0]) * timescale[durationData[1]]
    sim["duration"] = duration

    # Read log interval data and process it
    logIntervalData = sim["logInterval"].split()
    logInterval: unit.Quantity = float(logIntervalData[0]) * timescale[logIntervalData[1]]
    logIntervalInSteps: int = int(round(logInterval / timestep))
    sim["logInterval"] = logIntervalInSteps

    # Calculate number of steps
    nSteps: int = int(duration / timestep)
    sim["nSteps"] = nSteps

    # Calculate number of log steps (every 500 steps)
    nLogSteps: int = round(nSteps / 500)
    sim["nLogSteps"] = nLogSteps

    # Process temperature
    if "teperature" in sim:
        temperature: unit.Quantity = sim["teperature"] * unit.kelvin
        sim["teperature"] = temperature
    elif "temperatureRange" in sim:
        tempRange = [int(val) * unit.kelvin for val in sim["temperatureRange"]]
        sim["temperatureRange"] = tempRange

    sim["processed"] = True
    return sim
###########################################################################################
def init_reporters(simDir: str,
                    nSteps: int,
                      reportInterval: int,
                        simulation: app.Simulation,
                          dcdAtomSelections: List[Dict],
                            refPdb: FilePath) -> app.Simulation:
    """
    Initializes and returns a dictionary of reporters for a simulation.
    Uses OpenMM StateDataReporter for temperature, volume and density and energies
    Uses mdtraj DCDReporter for trajectories (this enables atom selection to lower file sizes)

    Args:
        simDir (str): The directory where the simulation files will be saved.
        nSteps (int): The total number of simulation steps.
        reportInterval (int): The interval at which the reporters will report.

    Returns:
        dict: A dictionary containing the initialized reporters. The keys are the names of the reporters and the values are lists containing the reporter object and the corresponding file.

    The function initializes and returns a dictionary of reporters for a simulation. It takes in the simulation directory, the total number of simulation steps, and the report interval. The reporters are initialized with the appropriate file paths and report intervals. The dictionary contains the initialized reporters, with the keys being the names of the reporters and the values being lists containing the reporter object and the corresponding file.

    Example usage:
    ```
    reporters = init_reporters("/path/to/simulation", 1000, 100)
    ```
    """
    ## as well as temperature, volume and density
    vitalsCsv = p.join(simDir, "vitals_report.csv")
    vitalsStateReporter: app.StateDataReporter = app.StateDataReporter(file = vitalsCsv, reportInterval = reportInterval, step = True,
                                            time = True, potentialEnergy = True, kineticEnergy = True,
                                            totalEnergy = True, temperature = True, volume = True,
                                            density = True)
    ## progresCsv gives you information on how long the simulation has been running and
    ## how long is left
    progressCsv: str = p.join(simDir, "progress_report.csv")
    progressStateReporter: app.StateDataReporter = app.StateDataReporter(file = progressCsv, progress = True, remainingTime = True,
                                            speed = True, elapsedTime = True, totalSteps = nSteps,
                                            reportInterval = reportInterval)
    ## dcdFile is the trajectory of the simulation
    dcdFile: str = p.join(simDir, "trajectory.dcd")
    ## get an subset of atoms that will have their trajectory saved
    dcdAtomSelection = []
    for selection in dcdAtomSelections:
        dcdAtomSelection.extend(drSelector.get_atom_indexes(selection["selection"], refPdb))
    # convert to numpy array
    dcdAtomSelection = np.array(dcdAtomSelection)
    ## create a mdtraj trajectory reporter
    dcdTrajectoryReporter: mdtraj.reporters.DCDReporter = mdtraj.reporters.DCDReporter(dcdFile,
                                                              reportInterval,
                                                              dcdAtomSelection)
    ## chkFile works as a checkpoint so that the simulation can be resumed if something goes wrong
    chkFile: str = p.join(simDir, "checkpoint.chk")
    chkReporter: app.CheckpointReporter = app.CheckpointReporter(file = chkFile, reportInterval = reportInterval, writeState = False)
    reporters: dict = {"vitals":[vitalsStateReporter,vitalsCsv],
                        "progress": [progressStateReporter,progressCsv],
                        "trajectory": [dcdTrajectoryReporter, dcdFile],
                        "checkpoint": [chkReporter, chkFile]}
    for rep in reporters:
        simulation.reporters.append(reporters[rep][0])

    return  simulation
###########################################################################################
def load_simulation_state(simulation: app.Simulation, saveFile: FilePath) -> app.Simulation:
    """
    Load the simulation state from a checkpoint or XML file.

    Args:
        simulation (app.Simulation): The simulation object.
        saveFile (str): The path to the checkpoint or XML file.

    Returns:
        app.Simulation: The modified simulation object.

    This function loads the simulation state from a checkpoint or XML file.
    The function determines the file extension of the `saveFile` and calls
    the appropriate method to load the state.

    If the file extension is ".chk", the function calls the `loadCheckpoint`
    method of the simulation object.

    If the file extension is ".xml", the function calls the `loadState`
    method of the simulation object.
    """
    # Determine the file extension of the saveFile
    saveFileExt: str = p.splitext(saveFile)[1]

    # Load the simulation state from a checkpoint file
    if saveFileExt == ".chk":
        simulation.loadCheckpoint(saveFile)
    # Load the simulation state from an XML file
    elif saveFileExt == ".xml":
        simulation.loadState(saveFile)

    ## reset time and step count of simulation to zero
    simulation.context.setTime(0.0)
    simulation.context.setStepCount(0)

    # Return the modified simulation object
    return simulation
###########################################################################################
@drLogger.monitor_progress_decorator()
# @drFirstAid.firstAid_handler()
@drCheckup.check_up_handler()
def run_molecular_dynamics(prmtop: app.AmberPrmtopFile,
                           inpcrd: app.AmberInpcrdFile,
                             sim: Dict,
                               saveFile: FilePath,
                                 outDir: FilePath,
                                   platform: openmm.Platform,
                                     refPdb: FilePath,
                                       config: Dict) -> FilePath:
    """
    Sets up and runs a molecular dynamics simulation
    Using either the isothermal-isobaric (NpT) or cannonical (NVT) ensembles

    Args:
        prmtop (app.Topology): The topology of the system.  
        inpcrd (app.Topology): The coordinates of the system.
        sim (Dict): dictionary containing simulation parameters
        saveFile (str): The path to the checkpoint or XML file.
        outDir (str): The path to the output directory.
        platform (openmm.Platform): The simulation platform.
        refPdb (str): The path to the reference PDB file.
        config (Dict): Dictionary containing information for all simluations in this run
    """

    ##TODO: REMOVE
    stepName = sim["SMEG"]
    drLogger.log_info(f"Running {stepName} for {config['proteinInfo']['proteinName']}")
    protName = config["proteinInfo"]["proteinName"]
    sim = process_sim_data(sim)
    ## create simulation directory
    simDir: str = p.join(outDir, sim["stepName"])
    os.makedirs(simDir, exist_ok=True)

    ## initialise a new system from parameters
    hardwareInfo = config["hardwareInfo"]
    simulation, integrator = initialise_simulation(prmtop, inpcrd, sim, saveFile, refPdb, platform, hardwareInfo)
    # set up intergrator and system
    # load state from previous simulation (or continue from checkpoint)
    simulation: app.Simulation = load_simulation_state(simulation, saveFile)
    # set up reporters
    totalSteps: int = simulation.currentStep + sim["nSteps"]
    reportInterval: int = sim["logInterval"]

    simulation: app.Simulation = init_reporters(simDir=simDir,
                                nSteps=totalSteps,
                                reportInterval=reportInterval,
                                simulation=simulation,
                                dcdAtomSelections= config["miscInfo"]["trajectorySelections"],
                                refPdb=refPdb
                                )
    # run NVT / NPT simulation
    simulation: app.Simulation = step_simulation(simulation, integrator, sim)

    # find name to call outFiles
    protName: str = p.basename(p.dirname(simDir))
    # save result as pdb - reset chain and residue Ids
    state: openmm.State = simulation.context.getState(getPositions=True, getEnergy=True)
    endPointPdb: str = p.join(simDir, f"{protName}.pdb")
    with open(endPointPdb, 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    ## reset the chain and residue ids to original
    drFixer.reset_chains_residues(refPdb, endPointPdb)

    ## create a PDB file with the same atoms as the trajectory
    trajectoryPdb = p.join(simDir, "trajectory.pdb")
    drSelector.slice_pdb_file(config["miscInfo"]["trajectorySelections"], endPointPdb, trajectoryPdb)
    # save simulation as XML
    saveXml: str = p.join(simDir, f"{stepName}.xml")
    simulation.saveState(saveXml)

    return saveXml


##########################################################################################

def step_simulation(simulation: app.Simulation, integrator: openmm.Integrator, sim: Dict) ->  app.Simulation:
    ## for simulations with constant temperature
    if "temperature" in sim:
        try:
            simulation.step(sim["nSteps"])
            return simulation
        except Exception as e:
            raise e
    ## for simulations with stepped temperature
    nTempSteps = len(sim["temperatureRange"])
    nStepsPerTempStep = round(sim["nSteps"] / nTempSteps)

    try:
        for temperature in sim["temperatureRange"]:
            integrator.setTemperature(temperature)
            simulation.step(nStepsPerTempStep)
        return simulation
    except Exception as e:
        raise e

###########################################################################################
def run_energy_minimisation(prmtop: app.AmberPrmtopFile,
                           inpcrd: app.AmberInpcrdFile,
                             sim: Dict,
                               saveFile: FilePath,
                                 outDir: DirectoryPath,
                                   platform: openmm.Platform,
                                     refPdb: FilePath,
                                       config: Dict) -> FilePath:

    """
    Run energy minimisation on a system.

    Args:
        prmtop (str): The path to the topology file.
        inpcrd (str): The path to the coordinates file.
        sim (dict): The simulation parameters.
        simDir (str): The path to the simulation directory.
        platform (openmm.Platform): The simulation platform.
        refPdb (str): The path to the reference PDB file.

    Returns:
        str: The path to the XML file containing the final state of the simulation.

    This function runs energy minimisation on a system. The function initializes the
    system, sets up the integrator and simulation, runs the energy minimisation, saves
    the result as a PDB file, resets the chain and residue Ids, saves the simulation
    state as an XML file, and returns the path to the XML file.
    """
    stepName: str = sim["stepName"]
    protName: str = config["proteinInfo"]["proteinName"]

    drLogger.log_info(f"Running {stepName} Step for: {protName} {' '*10}", True)
    ## create simluation directory
    simDir: str = p.join(outDir, sim["stepName"])
    os.makedirs(simDir, exist_ok=True)

    ## initialise a new system from parameters
    hardwareInfo: Dict = config["hardwareInfo"]
    simulation, _ = initialise_simulation(prmtop, inpcrd, sim, saveFile, refPdb, platform, hardwareInfo)

    ## if needed, convert temperature to kelvin
    if isinstance(sim["temperature"], int):
        sim["temperature"] = sim["temperature"] * unit.kelvin

    ## set coordinates of simulation 
    simulation.context.setPositions(inpcrd.positions)

    # Set box vectors if provided
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    # Run energy minimisation
    simulation.minimizeEnergy(maxIterations=sim['maxIterations'])

    # save result as pdb - reset chain and residue Ids
    state: openmm.State = simulation.context.getState(getPositions=True, getEnergy=True)
    minimisedPdb: str = p.join(simDir, f"{protName}.pdb")
    with open(minimisedPdb, 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology,
                                     state.getPositions(),
                                     output)
    drFixer.reset_chains_residues(refPdb, minimisedPdb)


    # Save simulation as XML
    saveXml: str = p.join(simDir, f"{stepName}.xml")
    simulation.saveState(saveXml)
    return saveXml

