## BASIC PYTHON LIBRARIES
import os
from os import path as p
import numpy as np
import yaml
import itertools
import argpass
import re
import shlex
import pytraj as pt
import subprocess
from subprocess import run
from StandardOperations import drManual, drESI
import textwrap
import io
import logging
import select
import sys

## ERROR HANDLING ##
import traceback
import inspect

## PARALLELISATION LIBRARIES
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from concurrent.futures.process import BrokenProcessPool



## CLEAN CODE
from typing import Optional, Dict, Tuple

class DirectoryPath:
    def __init__(self, path: str):
        if not os.path.isdir(path):
            raise ValueError(f"{path} is not a valid directory path")
        self.path = os.path.abspath(path)

    def __str__(self):
        return self.path

class FilePath:
    def __init__(self, path: str):
        if not os.path.isfile(path):
            raise ValueError(f"{path} is not a valid file path")
        self.path = os.path.abspath(path)

    def __str__(self):
        return self.path
    
def main(batchConfigYaml: Optional[FilePath] = None) -> None:
    '''
    Main function for drReferral
    processes input config file
    separated config file into separate files for each operation
    handles ESI desolvation simulations
    Runs drMD

    Args:
        Nothing
    Returns:
        Nothing
    '''
    if __name__ == "__main__":
        batchConfigYaml: FilePath = get_config_input_arg()
        ## read bacth config file into a dictionary
    try:
        batchConfig: dict = read_input_yaml(batchConfigYaml)
    except (FileNotFoundError, yaml.YAMLError, KeyError, TypeError, ValueError) as e:
        drSplash.print_config_error(e)
        ## unpack batchConfig into variables for this function

    ## Establish file directories and pdbFiles
    outDir: DirectoryPath = batchConfig["pathInfo"]["outputDir"]
    pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
    pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]
    drMD_src= __file__.removesuffix("drReferral.py")
    
    ## Checks to see if the config file has multiple "Operations" or is just a standard config file
    if batchConfig.get('Operations') != None:

        ## creates a new file for the new config files
        try:
            os.makedirs("Config_Files")
        except FileExistsError:
            print(f"Directory Config_Files already exists.")
        except PermissionError:
            print(f"Permission denied: Unable to create Config_Files.")
        except Exception as e:
            print(f"An error occurred: {e}")
        
        ## Process and run each "Operation"
        for Procedure in batchConfig['Operations']:
            
            ## finds all pdb files in the input file
            pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
            
            ## to run multiple drMD simulations one after the other, pdb files needed to be collated through aftercare
            ## This copys trajectory pdb files as well as the save pdb files.
            ## the trajectory files are deleted to prevent reducndancy
            if os.path.exists(f"{pdbDir}/trajectory.pdb"):
                os.remove(f"{pdbDir}/trajectory.pdb")
            pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]

            ##sets up a new outDir for each opperation
            newOutDir= outDir+ f"/{Procedure['OperationName']}"
            (batchConfig["pathInfo"]).update({"outputDir": newOutDir} )

            ## handles "condensed" config files that use standard operations configs found in "StandardOpperations"
            if Procedure.get('simulationInfo') == None:
                simulationInfo= drManual.StandardOperation(Procedure['OperationName'])
                Procedure.update( {'simulationInfo': simulationInfo['simulationInfo']})
                Type= drManual.GetSimulationType(Procedure['OperationName'])
                Procedure.update({'Type': Type})

            ## Handles the diffrent system types (solution phase, gas phase and ESI simulations)
            if Procedure['Type'] == "Solution":
                write_referral(batchConfig, Procedure)
                Procedure_Name= Procedure['OperationName']
                drMDCommand = f"python3 {drMD_src}SolutionWing/drMD.py --config Config_Files/{Procedure_Name}.yaml "
                Run_drMD(drMDCommand)
            elif Procedure['Type'] == "Vacuum":
                ## removes water from the PDB files
                for PDB in pdbFiles:
                    Water_Count= Count_Water(PDB)
                    if Water_Count !=0:
                        Remove_Solvent(PDB)
                write_referral(batchConfig, Procedure)
                Procedure_Name= Procedure['OperationName']
                drMDCommand = f"python3 {drMD_src}GasWing/drMD.py --config Config_Files/{Procedure_Name}.yaml "
                Run_drMD(drMDCommand)

            elif Procedure['Type'] == "ESI":
                ESI_Count= 0
                Procedure_Name= Procedure['OperationName']
                for PDB in pdbFiles:
                    Water_Count = Count_Water(PDB)

                ## while ther is still water in the pdb
                while Water_Count >=1:

                    ## sets up an incremented ESI simulation
                    ESI_Count +=1
                    pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
                    if os.path.exists(f"{pdbDir}/trajectory.pdb"):
                        os.remove(f"{pdbDir}/trajectory.pdb")
                    pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]
                    newOutDir= outDir+ f"/{Procedure_Name}/ESI_{ESI_Count}"
                    (batchConfig["pathInfo"]).update({"outputDir": newOutDir})

                    ## if no Cutoff is specified, use the defult of 10 angstroms
                    if  Procedure.get('Cutoff') == None:
                        Cutoff= 10
                    else:
                        Cutoff= Procedure['Cutoff']

                    ## Removes all waters outside the cutoff distance    
                    for PDBs in pdbFiles:
                        drESI.RemoveWater(Cutoff, PDBs, ESI_Count)
                    Procedure.update({'OperationName': f'ESI_{ESI_Count}'})

                    ## run the simulation with the new pdb file
                    write_referral(batchConfig, Procedure)
                    drMDCommand = f"python3 {drMD_src}ESIWing/drMD.py --config Config_Files/{Procedure_Name}_{ESI_Count}.yaml "
                    Run_drMD(drMDCommand)

                    ## input directory is the directory that contains the outputs from the previous simulation
                    newInDir= newOutDir+ "/00_collated_pdbs/"+ Procedure["simulationInfo"][-1]["stepName"]
                    (batchConfig["pathInfo"]).update({"inputDir": newInDir})
                    pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
                    pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]

                    ## counts how much water is left in the pdb structure
                    for PDB in pdbFiles:
                        Water_Count = Count_Water(PDB)

                        ## if a maximum iteration has been specified, stop the simulations when reached
                        if Procedure["maxIterations"] != None:
                            if ESI_Count >= Procedure["maxIterations"]:
                                Water_Count=0


            newInDir= newOutDir+ "/00_collated_pdbs/"+ Procedure["simulationInfo"][-1]["stepName"]
            (batchConfig["pathInfo"]).update({"inputDir": newInDir})
            pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
            pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]


    ## If the config file is the standard config file then run as normal
    elif batchConfig.get('Operations') == None and batchConfig.get('simulationInfo') != None:

        ## assume if boxGeomtry and size are not specified then vacuum simulations are to be ran
        if (batchConfig['miscInfo']).get('boxGeometry') == None and (batchConfig['miscInfo']).get('boxSize') == None:
            drMDCommand= f"python3 {drMD_src}GasWing/drMD.py --config {batchConfigYaml}"
            Run_drMD(drMDCommand)
            
        else:
            drMDCommand= f"python3 {drMD_src}SolutionWing/drMD.py --config {batchConfigYaml} "
            Run_drMD(drMDCommand)
            


######################################################################################################
def write_referral(Equipment: dict, Procedure: dict) -> None:
    '''
    Writes a new config file to be run in drMD for the operation


    Args:
        Equipment: this contains pathInfo and Hardware Info that should stay constant (anything before the operations)
        Procedure: this contains simulation info and restiants that should change between operations
    Returns:
        Nothing
    '''

    ## Checks if aftercareInfo is present
    ## if not sets up endpointInfo so that following simulations can run endpointPDBs
    if Procedure.get('aftercareInfo') == None:
        stepNames= [Procedure["simulationInfo"][-1]["stepName"]]
        endpointInfo= {"stepNames": stepNames}
        aftercareInfo= {"endPointInfo": endpointInfo}
    else:
        aftercareInfo = Procedure["aftercareInfo"]
    
    ##if miscInfo is constant
    if Procedure.get("miscInfo") == None:
        Referral= { 
        "pathInfo": Equipment["pathInfo"],
        "hardwareInfo": Equipment["hardwareInfo"],
        "miscInfo": Equipment["miscInfo"],
        "simulationInfo": Procedure["simulationInfo"],
        'aftercareInfo': aftercareInfo
        }
    
    ## If miscInfo changes between operations
    elif Procedure.get("miscInfo") != None:
        Referral= { 
        "pathInfo": Equipment["pathInfo"],
        "hardwareInfo": Equipment["hardwareInfo"],
        "miscInfo": Procedure["miscInfo"],
        "simulationInfo": Procedure["simulationInfo"],
        'aftercareInfo': aftercareInfo
        }
    
    ## if there are restaints to be added
    if Equipment.get("equilibriationRestraints") != None:
        Referral.update({"equilibriationRestraints": Procedure["equilibriationRestraints"]})
    Procedure_Name= Procedure["OperationName"]
    
    ## write the new config file
    with open(f'Config_Files/{Procedure_Name}.yaml', 'w') as Config:
        yaml.dump(Referral, Config)
    return 
######################################################################################################
def get_config_input_arg() -> FilePath:
    """
    Sets up argpass to read the config.yaml file from command line
    Reads a YAML file using the "--config" flag with argpass

    Returns:
    - configFile (FilePath)
    """
    # create an argpass parser, read config file,
    parser = argpass.ArgumentParser()
    parser.add_argument(f"--config")
    args = parser.parse_args()

    configFile: FilePath = args.config

    return configFile

######################################################################################################   
def read_input_yaml(configFile: FilePath) -> dict:
    """
    Reads YAML file into a dict

    Args:
    - configFile (str): Path to the YAML configuration file.

    Returns:
    - config (dict): Parsed YAML content as a dictionary.
    """
    yellow = "\033[33m"
    reset = "\033[0m"
    teal = "\033[38;5;37m"
    try:
        with open(configFile, "r") as yamlFile:
            config: dict = yaml.safe_load(yamlFile)
            return config
    except FileNotFoundError:
        print(f"-->{' '*4}Config file {configFile} not found.")
        exit(1)
    except yaml.YAMLError as exc:
        print(f"-->{' '*4}{yellow}Error while parsing YAML file:{reset}")
        if hasattr(exc, 'problem_mark'):
            mark = exc.problem_mark
            print(f"{' '*7}Problem found at line {mark.line + 1}, column {mark.column + 1}:")
            if exc.context:
                print(f"{' '*7}{exc.problem} {exc.context}")
            else:
                print(f"{' '*7}{exc.problem}")
            print(f"{' '*7}Please correct the data and retry.")
        else:
            print(f"{' '*7}Something went wrong while parsing the YAML file.")
            print(f"{' '*7}Please check the file for syntax errors.")
        print(f"\n{teal}TIP:{reset} Large language models (LLMs) like GPT-4 can be helpful for debugging YAML files.")
        print(f"{' '*5}If you get stuck with the formatting, ask a LLM for help!")
        exit(1)

######################################################################################################
def read_config(configYaml: str) -> dict:
    """
    Reads a config.yaml file and returns its contents as a dictionary.
    This is performed by drOperator on automatically generated configs,

    Args:
        configYaml (str): The path to the config.yaml file.

    Returns:
        dict: The contents of the config.yaml file as a dictionary.
    """
    try:
        # Open the config.yaml file and read its contents into a dictionary
        with open(configYaml, "r") as yamlFile:
            config: dict = yaml.safe_load(yamlFile)
    except FileNotFoundError:
        raise FileNotFoundError(f"config file {configYaml} not found")
    except yaml.YAMLError as exc:
        raise yaml.YAMLError(f"Error parsing YAML file:", exc)

    return config
######################################################################################################
def Count_Water(PDBfileName: str) -> int:
    """
    Reads a pdbfile to count the number of water molecules present

    Args:
        PDBfileName (str): file path to pdbfile 

    Returns:
        water_count(int): number of waters in the pdb file.
    """
    Water_Count=0
    with open(PDBfileName, 'r') as PDB:
        rows= PDB.readlines()
        for line in rows:
            if line.find('HOH') != -1 or line.find('WAT') != -1:
                Water_Count +=1
    return Water_Count
######################################################################################################
def Remove_Solvent(PDBfile: str) -> None:
    """
    removes all solvent and Ions in the pdb file for vacuum simulations using pytraj

    Args:
        PDBfile(str): file path to pdbfile 
    """
    
    traj= pt.load(PDBfile)
    stripped_traj = traj.strip(':HOH')
    stripped_traj = traj.strip('@H')
    ## if there are any Ions
    if (traj.top.select("@Cl-")).any != None:
        stripped_traj = traj.strip(':Cl-')
    elif (traj.top.select("@Na+")).any != None:
        stripped_traj = traj.strip(':Na+')
    pt.write_traj(f"{PDBfile}", stripped_traj, overwrite=True)

######################################################################################################
def Run_drMD(command, print_output= True,log_output=False, check=True, *args, **kwargs):
    """
    runs drMD as a command and prints the output 

    Args:
        command(str): command to run drMD
    """
    shell = isinstance(command.split(), str)
    logging.debug(f"Running command: {command}")
    process = subprocess.Popen(  # type: ignore
        command.split(),
        shell=shell,
        bufsize=1,  # Output is line buffered, required to print output in real time
        universal_newlines=True,  # Required for line buffering
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        #*args,
        #**kwargs,
    )
    stdoutbuf = io.StringIO()
    stderrbuf = io.StringIO()
    stdout_fileno = process.stdout.fileno()  # type: ignore
    stderr_fileno = process.stderr.fileno()  # type: ignore
    # This returns None until the process terminates
    while process.poll() is None:

        # select() waits until there is data to read (or an "exceptional case") on any of the streams
        readready, writeready, exceptionready = select.select(
            [process.stdout, process.stderr],
            [],
            [process.stdout, process.stderr],
            0.5,
        )

        # Check if what is ready is a stream, and if so, which stream.
        # Copy the stream to the buffer so we can use it,
        # and print it to stdout/stderr in real time if print_output is True.
        for stream in readready:
            if stream.fileno() == stdout_fileno:
                line = process.stdout.readline()  # type: ignore
                stdoutbuf.write(line)
                if print_output:
                    sys.stdout.write(line)
            elif stream.fileno() == stderr_fileno:
                line = process.stderr.readline()  # type: ignore
                stderrbuf.write(line)
                if print_output:
                    sys.stderr.write(line)
            else:
                raise Exception(
                    f"Unknown file descriptor in select result. Fileno: {stream.fileno()}"
                )

        # If what is ready is an exceptional situation, blow up I guess;
        # I haven't encountered this and this should probably do something more sophisticated.
        for stream in exceptionready:
            if stream.fileno() == stdout_fileno:
                raise Exception("Exception on stdout")
            elif stream.fileno() == stderr_fileno:
                raise Exception("Exception on stderr")
            else:
                raise Exception(
                    f"Unknown exception in select result. Fileno: {stream.fileno()}"
                )

    # Check for any remaining output after the process has exited.
    # Without this, the last line of output may not be printed,
    # if output is buffered (very normal)
    # and the process doesn't explictly flush upon exit
    # (also very normal, and will definitely happen if the process crashes or gets KILLed).
    for stream in [process.stdout, process.stderr]:
        for line in stream.readlines():
            if stream.fileno() == stdout_fileno:
                stdoutbuf.write(line)
                if print_output:
                    sys.stdout.write(line)
            elif stream.fileno() == stderr_fileno:
                stderrbuf.write(line)
                if print_output:
                    sys.stderr.write(line)

    # We'd like to just seek(0) on the stdout/stderr buffers, but "underlying stream is not seekable",
    # So we create new buffers above, write to them line by line, and replace the old ones with these.
    process.stdout.close()  # type: ignore
    stdoutbuf.seek(0)
    process.stdout = stdoutbuf
    process.stderr.close()  # type: ignore
    stderrbuf.seek(0)
    process.stderr = stderrbuf

    if check and process.returncode != 0:
        msg = f"Command failed with exit code {process.returncode}: {command}"
        logging.error(msg)
        logging.info(f"stdout: {process.stdout.getvalue()}")
        logging.info(f"stderr: {process.stderr.getvalue()}")
        raise Exception(msg)

    logging.info(f"Command completed with return code {process.returncode}: {command}")

    # The user may have already seen the output in std out/err,
    # but logging it here also logs it to syslog (if configured).
    if log_output:
        # Note that .getvalue() is not (always?) available on normal Popen stdout/stderr,
        # but it is available on our StringIO objects.
        # .getvalue() doesn't change the seek position.
        logging.info(f"stdout: {process.stdout.getvalue()}")
        logging.info(f"stderr: {process.stderr.getvalue()}")

    # Now that we've set stdout/err to StringIO objects,
    # we can return the Popen object as a MagicPopen object.
    
    return process

######################################################################################################
def ConfigChecker(Config: dict) -> dict:


    return Config
######################################################################################################

if __name__ == "__main__":
    main()