## BASIC PYTHON LIBRARIES
import os
from os import path as p
from functools import wraps

## OPENMM LIBRARIES
import openmm
from openmm import app
from openmm import OpenMMException
import  simtk.unit  as unit

## drMD LIBRARIES
from Surgery import drSim
from UtilitiesCloset import drSplash, drSplicer
from ExaminationRoom import drLogger

## CLEAN CODE
from typing import Tuple, Union, Dict, List, Any
from os import PathLike


#######################################################################
def firstAid_handler():
    """
    A decorator that handles firstAid retries.

    Implementation in drSim.py:
        @drTriage.firstAid_handler(drTriage.run_firstAid_simulation, maxRetries=5)
        def run_molecular_dynamics(prmtop: str, inpcrd: str, sim: dict, saveFile: str, outDir: str, platform: openmm.Platform, refPdb: str) -> str:
        
    How it works:

        1. Run the target simulation function
            - If the simulation succeeds first time round, return the saveFile

        2. If the simulation fails giving a ValueError or OpenMMException:

            a) Work out how long the simulation took before it crashed

            b) Rename the output files in the failed simulation directory

            c) Find the checkpoint file of the failed simulation        

            d) Run an energy minimisation step, starting with the geometry in the exploded checkpoint file

            e) Reset arguments for the target simulation function 

            f) adjust the target simulation number of steps parameter based on the output os step a)

            f) Run the target simulation function using the output geometry from step d)

        3. Step 2 will continue until either the target simulation succeeds or the max number of retries is reached

    Args:
        firstAid_function (callable): The function to be decorated.   
        maxRetries (int, optional): The maximum number of retries. Defaults to 10.

    Returns:
        None
    """
    ## set up decorator and wrapper for our simulation function
    def decorator(simulationFunction):
        @wraps(simulationFunction)
        def wrapper(*args, **kwargs):
            ## initialise a retry counter
            retries: int = 0
            ## keep running firstAid simulations until max retries is reached or the simulation succeeds
            maxRetries: int = kwargs["config"]["miscInfo"]["firstAidMaxRetries"]
            if maxRetries == 0:
                return
            while retries < maxRetries:
                ## try to run the simulation - if this runs without errors, the rest of this function is skipped
                try:
                    saveFile: Union[PathLike, str] = simulationFunction(*args, **kwargs)
                    ## in this case, some firstAid steps have been run
                    if retries > 0:
                        ## let user know whats going on
                        drLogger.log_info(f"Success after {retries} tries.", True)
                        ## merge partial reports and trajectories
                        runOutDir: Union[PathLike, str] = kwargs["outDir"]
                        simDir: Union[PathLike, str] = p.join(runOutDir, kwargs["sim"]["stepName"])
                        drSplicer.merge_partial_outputs(simDir=simDir,
                                                         pdbFile=kwargs["refPdb"],
                                                           simInfo=kwargs["sim"],
                                                            config = kwargs["config"])
                    ## return the saveFile (XML) of the simulation to be used by subsequent simulations
                    return saveFile
                ## if our simulation crashes due to a numerical error or an OpenMM exception
                ## run firstAid protocol to try and recover
                except (OpenMMException, ValueError) as errorOpenMM:
                    saveFile, retries = run_first_aid_protocol(retries, maxRetries, *args, **kwargs)
                except Exception as error:
                    drLogger.log_info(f"Unexpected Error for {kwargs['sim']['stepName']}:\n{error}", True, True)
                    raise error
            else:
                ## If we have got here, the firstAid has failed
                ## let user know and merge output reporters and trajectories
                drLogger.log_info(f"Max retries reached. Stopping.", True, True)
                drSplash.print_first_aid_failed("Particle coordinate NaN")
                raise OpenMMException
        return wrapper
    return decorator

#########################################################################################################################
def firstAid_handler():
    def decorator(simulationFunction):
        @wraps(simulationFunction)
        def wrapper(*args, **kwargs):
            retries: int = 0
            maxRetries: int = kwargs["config"]["miscInfo"].get("firstAidMaxRetries", 10)
            lastError = None
            while retries < maxRetries:
                try:
                    saveFile: Union[str, PathLike] = simulationFunction(*args, **kwargs)
                    if retries > 0:
                        drLogger.log_info(f"Success after {retries} tries.", True)
                        runOutDir: Union[str, PathLike] = kwargs["outDir"]
                        simDir: Union[str, PathLike] = p.join(runOutDir, kwargs["sim"]["stepName"])
                        drSplicer.merge_partial_outputs(simDir=simDir,
                                                        pdbFile=kwargs["refPdb"],
                                                        simInfo=kwargs["sim"],
                                                        config=kwargs["config"])
                    return saveFile
                except (OpenMMException, ValueError) as errorOpenMM:
                    lastError = errorOpenMM
                    saveFile, retries = run_first_aid_protocol(retries, maxRetries, *args, **kwargs)
                except Exception as error:
                    drLogger.log_info(f"Unexpected Error for {kwargs['sim']['stepName']}:\n{error}", True, True)
                    raise error
                retries += 1
            
            # If max retries are reached, raise the last caught OpenMMException or ValueError
            drLogger.log_info(f"Max retries reached. Stopping.", True, True)
            drSplash.print_first_aid_failed("Particle coordinate NaN")
            if lastError:
                raise lastError
            else:
                raise "Unknown Error"

        return wrapper
    return decorator


#########################################################################################################################

def run_first_aid_protocol(retries: int, maxRetries: int, *args, **kwargs):
    retries += 1
    ## let user know whats going on
    if retries == 1:
        drSplash.print_performing_first_aid()

    drLogger.log_info(f"Attempting simulation firstAid, try {retries} of {maxRetries}", True)
    ## find simulation output directory
    runOutDir: Union[PathLike, str] = kwargs["outDir"]
    simDir: Union[PathLike, str] = p.join(runOutDir, kwargs["sim"]["stepName"])
    ## look for checkpoint file of exploded simulation
    ## if none found, simulatoin has exploded on first step, try simulation again
    firstAidDir, explodedCheckpoint = pre_firstAid_processing(simDir, retries)
    ## if no checkpoint file is found, simulation has exploded on first step
    ## in this case, reload from the original saveFile from the preceding run
    if explodedCheckpoint is None:
        explodedCheckpoint = kwargs["saveFile"]

    ## get current timestep of simulation
    currentNsteps: int = get_nsteps_at_crash(explodedCheckpoint, kwargs["prmtop"])
    ## prepare arguments for firstAid simulation
    kwargs: Dict = prepare_arguments_for_firstAid(kwargs, firstAidDir, explodedCheckpoint, retries)
    ## run firstAid simulation using modified keyword arguments
    try:
        saveFile: Union[PathLike, str] = run_energy_minimisation_then_npt(*args, **kwargs)
    except Exception as e:
        drLogger.log_info(e, True,True)
    ## reset keyword arguments so simulation can be resumed
    reset_keyword_arguments(kwargs, runOutDir, currentNsteps)

    return saveFile, retries

#######################################################################
def get_nsteps_at_crash(explodedCheckpoint: Union[PathLike, str], prmtop: app.AmberPrmtopFile) -> int:
    """
    Gets the current number of steps that have been run in the simulation
    Later on we will use this to modify the number of steps we need when we restart the simulation

    Args:
        explodedCheckpoint (Union[PathLike, str]): The path to the exploded checkpoint file
        prmtop (app.AmberPrmtopFile): The prmtop file

    Returns:
        int: The number of steps that have been run in the simulation
    """
    try:
        ## create a system, integrator and simulation using prmtop
        system: openmm.System = prmtop.createSystem()
        integrator: openmm.Integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
        simulation: app.Simulation = app.Simulation(prmtop.topology, system, integrator)
        ## load the checkpoint file
        with open(explodedCheckpoint, 'rb') as f:
            simulation.context.loadCheckpoint(f.read())
        ## get the number of steps from the context
        nsteps: int = simulation.context.getStepCount()
    except:
        nsteps = 0
    return nsteps    

#######################################################################
def reset_keyword_arguments(kwargs: dict, outDir: Union[PathLike, str], currentNsteps: int):
    """
    Before resuming our simulation, we need to reset the keyword arguments
    so that the simulation can be resumed correctly

    Args:
        kwargs (dict): The keyword arguments
    
    """
    ## reset outDir to the original output directory
    kwargs["outDir"] = outDir
    ## remove firstAidTries argument
    del kwargs["firstAidTries"]
    ## reduce the number of steps in the resumed simulation
    ## by th enumber of steps that were performed prior to the crash
    kwargs["sim"]["nSteps"] = kwargs["sim"]["nSteps"] - currentNsteps

    return kwargs
#######################################################################
def prepare_arguments_for_firstAid(kwargs: dict,
                                  firstAidDir: Union[PathLike, str],
                                    explodedCheckpoint: Union[PathLike, str],
                                      retries: int) -> dict:
    """
    Before running the firstAid energy minimisation, we neeed to modify the keyword arguments 

    Args:
        kwargs (dict): The keyword arguments
        firstAidDir (Union[PathLike, str]): The path to the firstAid directory
        explodedCheckpoint (Union[PathLike, str]): The path to the exploded checkpoint file
        retries (int): The number of retries

    Returns:
        dict: The modified keyword arguments
    """
    ## set outDir to the firstAid directory, within the simulation directory
    kwargs["outDir"] = firstAidDir
    ## chage the saveFile to the checkpint of the crashed simulation 
    kwargs["saveFile"] = explodedCheckpoint
    ## pass firstAidTries argumet to the firstAid function
    kwargs["firstAidTries"] = retries

    return kwargs

#######################################################################
def pre_firstAid_processing(simDir: Union[PathLike, str],
                           retries: int) -> Tuple[Union[PathLike, str], Union[PathLike, str]]:
    """
    A simulation has just crashed, probably due to a numerical error!
    This function looks for the checkpoint file of the exploded simulation
    If no checkpoint file is found, it raises a FileNotFoundError
    This will result in the simulation being restarted
    This function also renames the output files in the exploded directory so they don't get overwritten
    This function also creates a firstAid directory for firstAid steps to be run in

    Args:
        simDir (Union[PathLike, str]): The path to the simulation directory
        retries (int): The number of retries

    Returns:
        firstAidDir (Union[PathLike, str]): The path to the firstAid directory
        explodedCheckpoint (Union[PathLike, str]): The path to the exploded checkpoint file
    """
    ## make a firstAid directory - this will be where each firstAid step is performed
    firstAidDir = p.join(simDir, "firstAid")
    os.makedirs(firstAidDir, exist_ok=True)

    ## rename output files in exploded directory
    rename_output_files(simDir, retries)

    ## find renamed checkpoint
    explodedCheckpoint = p.join(simDir, f"checkpoint_partial_{str(retries)}.chk")
    if not p.isfile(explodedCheckpoint):
        explodedCheckpoint = None
    return firstAidDir, explodedCheckpoint
#######################################################################

def rename_output_files(outDir: Union[PathLike, str], retries: int) -> None:
    """
    Renames reporter .csv and trajectory .dcd files in the output directory
    such that they have _partial_[retries].csv and _partial_[retries].dcd suffixes
    this ensures that they are not overwritten by the next simulation

    Args:
        outDir (Union[PathLike, str]): The path to the output directory
        retries (int): The number of retries

    Returns:
        None
    """
    ## look for reporter.csv and trajectory.dcd files
    extensions: List[str] = [".dcd", ".csv", ".chk"]
    for file in os.listdir(outDir):
        ## skip files with "partial" in the name
        if "partial" in file:
            continue

        ## rename files
        if p.splitext(file)[1] in extensions:
            os.rename(p.join(outDir, file),
                       p.join(outDir,f"{p.splitext(file)[0]}_partial_{str(retries)}{p.splitext(file)[1]}"))
#######################################################################
def run_firstAid_energy_minimisation(prmtop: str,
                                      inpcrd: str,
                                        sim: dict,
                                          saveFile: str,
                                            outDir: str,
                                              platform: openmm.Platform,
                                                refPdb: str,
                                                  firstAidTries: int,
                                                    config: dict) -> None:
    """
    Run a firstAid simulation

    Args:
        prmtop (str): The path to the topology file.
        inpcrd (str): The path to the coordinates file.
        sim (dict): The simulation parameters.
        saveFile (str): The path to the checkpoint or XML file.
        simDir (str): The path to the simulation directory.
        platform (openmm.Platform): The simulation platform.
        refPdb (str): The path to the reference PDB file.

    Returns:
        str: The path to the XML file containing the final state of the simulation.

    This function runs a firstAid simulation. 
    First, it makes a simulation dictionary containing parameters to be passed to drSim.run_energy_minimisation

    """
    firstAidSimInfo =  {
        "stepName" : f"firstAid_step_EM_{firstAidTries}",
        "simulationType": "NPT",
        "duration" : "10 ps",
        "timestep" : "0.5 ps",
        "temperature" : 300,
        "logInterval" : "2 ps",
        "maxIterations" : 1000
        }
    firstAidSimInfo = drSim.process_sim_data(firstAidSimInfo)

    ## run firstAid energy minimisation
    firstAidSaveFile = drSim.run_energy_minimisation(prmtop = prmtop,
                                                   inpcrd = inpcrd,
                                                   sim = firstAidSimInfo,
                                                    saveFile= saveFile,
                                                     outDir = outDir,
                                                     platform = platform,
                                                      refPdb= refPdb,
                                                      config= config)
    return firstAidSaveFile
#######################################################################
def run_firstAid_npt_simulation(prmtop: str, inpcrd: str, sim: dict, saveFile: str, outDir: str, platform: openmm.Platform, refPdb: str, firstAidTries: int, config: dict) -> None:

    firstAidSimInfo =  {
        "stepName" : f"firstAid_step_NpT_{firstAidTries}",
        "simulationType": "NPT",
        "duration" : "10 ps",
        "timestep" : "0.1 fs",
        "temperature" : 10,
        "logInterval" : "1 ps",
        }
    
    firstAidSimInfo = drSim.process_sim_data(firstAidSimInfo)

    firstAidSaveFile = drSim.run_molecular_dynamics(prmtop = prmtop,
                                                inpcrd = inpcrd,
                                                sim = firstAidSimInfo,
                                                saveFile= saveFile,
                                                    outDir = outDir,
                                                    platform = platform,
                                                    refPdb= refPdb,
                                                    config=config)


    return firstAidSaveFile

#######################################################################
def run_energy_minimisation_then_npt(prmtop: str, inpcrd: str, sim: dict, saveFile: str, outDir: str, platform: openmm.Platform, refPdb: str, firstAidTries: int, config: dict) -> str:


    emFailed: bool = False
    try:
        firstAidSaveFile_energyMinimisation = run_firstAid_energy_minimisation(prmtop = prmtop,
                                                                            inpcrd = inpcrd,
                                                                              sim = sim,
                                                                                outDir = outDir,
                                                                                  platform = platform,
                                                                                    refPdb = refPdb,
                                                                                      saveFile = saveFile,
                                                                                        firstAidTries=firstAidTries,
                                                                                          config=config)
    ## if em step fails, just run npt
    except Exception as e:
        emFailed = True
    if emFailed:
        firstAidSaveFile_energyMinimisation  = saveFile
    nptFailed = False
    try:
        firstAidSaveFile_nptSimulation = run_firstAid_npt_simulation(prmtop = prmtop,
                                                                  inpcrd = inpcrd,
                                                                    sim = sim,
                                                                      outDir = outDir,
                                                                        platform = platform,
                                                                          refPdb = refPdb,
                                                                            saveFile = firstAidSaveFile_energyMinimisation,
                                                                              firstAidTries=firstAidTries,
                                                                                config=config)
    ## if npt step fails, return em savefile
    except Exception as e:
        nptFailed = True
               
    ## depending on what step worked, return appropriate savefile
    if nptFailed and emFailed:
        return saveFile
    elif nptFailed:
        return firstAidSaveFile_energyMinimisation
    else:
        return firstAidSaveFile_nptSimulation
    
