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
from StandardOperations import drManual
from ESIWing.SpecialistEquipment import  drESI
from drOperationReportIndex import build_operation_report_index
import textwrap
import MDAnalysis as mda
import io
import logging
import select
import sys
from shutil import copy, copy2, copyfile
from csv import DictReader, DictWriter

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
        if not p.isdir(path):
            raise ValueError(f"{path} is not a valid directory path")
        self.path = p.abspath(path)

    def __str__(self):
        return self.path

class FilePath:
    def __init__(self, path: str):
        if not p.isfile(path):
            raise ValueError(f"{path} is not a valid file path")
        self.path = p.abspath(path)

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
    except (FileNotFoundError, yaml.YAMLError, KeyError, typeError, ValueError) as e:
        drSplash.print_config_error(e)
        ## unpack batchConfig into variables for this function

    ## Establish file directories and pdbFiles
    outDir: DirectoryPath = batchConfig["pathInfo"]["outputDir"]
    pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
    pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]
    drMD_src = p.dirname(p.abspath(__file__))

    ## Checks to see if the config file has multiple "Operations" or is just a standard config file
    if batchConfig.get('Operations') != None:
        combinedReport= p.join(outDir, "Full_Methods_Report.md")
        try:
            os.remove(combinedReport)
        except FileNotFoundError:
            pass

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
            if p.exists(f"{pdbDir}/trajectory.pdb"):
                os.remove(f"{pdbDir}/trajectory.pdb")
            pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]

            ##sets up a new outDir for each opperation
            newOutDir= outDir+ f"/{Procedure['OperationName']}"
            (batchConfig["pathInfo"]).update({"outputDir": newOutDir} )

            ## handles "condensed" config files that use standard operations configs found in "StandardOpperations"
            if Procedure.get('simulationInfo') == None:
                simulationInfo= drManual.standardOperation(Procedure['OperationName'])
                if simulationInfo != None:
                    Procedure.update( {'simulationInfo': simulationInfo['simulationInfo']})
                type= drManual.getSimulationType(Procedure['OperationName'])
                Procedure.update({'type': type})

            ## Handles the diffrent system types (solution phase, gas phase and ESI simulations)
            if Procedure["type"] == "Solution":
                write_referral(batchConfig, Procedure)
                Procedure_Name= Procedure["OperationName"]
                drMDCommand = build_drmd_command("SolutionWing", f"Config_Files/{Procedure_Name}.yaml")
                Run_drMD(drMDCommand)
            elif Procedure["type"] == "Vacuum":
                ## removes water from the PDB files
                for PDB in pdbFiles:
                    Water_Count= Count_Water(PDB)
                    if Water_Count !=0:
                        Remove_Solvent(PDB)
                write_referral(batchConfig, Procedure)
                Procedure_Name= Procedure["OperationName"]
                drMDCommand = build_drmd_command("GasWing", f"Config_Files/{Procedure_Name}.yaml")
                Run_drMD(drMDCommand)

            elif Procedure["type"] == "ESI":
                if Procedure.get("mode") == None:
                    mode= "Positive"
                else:
                    mode= Procedure["mode"]
                ESI_Count= 0
                Procedure_Name= Procedure["OperationName"]
                for PDB in pdbFiles:
                    Water_Count = Count_Water(PDB)

                ## while ther is still water in the pdb
                while Water_Count >=1:

                    ## sets up an incremented ESI simulation
                    ESI_Count +=1
                    simulationInfo= drESI.ESIOperation(ESI_Count, Water_Count)
                    Procedure.update( {'simulationInfo': simulationInfo['simulationInfo']})
                    pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
                    if p.exists(f"{pdbDir}/trajectory.pdb"):
                        os.remove(f"{pdbDir}/trajectory.pdb")
                    pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]
                    newOutDir= outDir+ f"/{Procedure_Name}/ESI_{ESI_Count}"
                    (batchConfig["pathInfo"]).update({"outputDir": newOutDir})

                    ## if no cutoff is specified, use the defult of 10 angstroms
                    if  Procedure.get('cutoff') == None:
                        cutoff= 100
                        Procedure.update({'cutoff': cutoff})
                    else:
                        cutoff= Procedure['cutoff']

                    ## Removes all waters outside the cutoff distance 

                    Do_Prep= True
                    for PDBs in pdbFiles:
                        if PDBs.find(f"ESI_{ESI_Count}") != -1:
                            Do_Prep=False
                    if Do_Prep == False:
                        for PDBs in pdbFiles:
                            if PDBs.find(f"ESI_{ESI_Count}") == -1:
                                os.remove(PDBs)
                    else:
                        for PDBs in pdbFiles:
                            if PDBs.find(f"ESI_{ESI_Count}") == -1:
                                drESI.ESI_Handler(cutoff, PDBs, ESI_Count, mode)
                    Procedure.update({'OperationName': f'ESI_{ESI_Count}'})

                    ## run the simulation with the new pdb file
                    write_referral(batchConfig, Procedure)
                    drMDCommand = build_drmd_command("ESIWing", f"Config_Files/{Procedure_Name}_{ESI_Count}.yaml")
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
                        if Procedure.get("maxESIIter") != None:
                            if ESI_Count >= Procedure["maxESIIter"]:
                                Water_Count=0
                combine_simulation_ESI(outDir, batchConfig)
                
            elif Procedure["type"] == "VacuumCM":
                Procedure_Name= Procedure["OperationName"]
                totalDuration= Procedure["simulationInfo"][-1]["duration"] 
                if totalDuration.find("ns") != -1:
                    totalDuration= float(totalDuration.replace("ns", ""))*1000
                elif totalDuration.find("ps") != -1:
                    totalDuration= float(totalDuration.replace("ps", ""))
                totalSims= totalDuration/20
                for I in range(int(totalSims)):
                    Procedure["simulationInfo"][-1]["duration"]= "20 ps"
                    newOutDir= outDir+ f"/{Procedure_Name}/CM_{I+1}"
                    (batchConfig["pathInfo"]).update({"outputDir": newOutDir})
                    write_referral(batchConfig, Procedure)
                    drMDCommand = build_drmd_command("GasWing", f"Config_Files/{Procedure['OperationName']}.yaml")
                    Run_drMD(drMDCommand)
                    newInDir= newOutDir+ "/00_collated_pdbs/"+ Procedure["simulationInfo"][-1]["stepName"]
                    (batchConfig["pathInfo"]).update({"inputDir": newInDir})
                    pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
                    os.remove(f"{pdbDir}/trajectory.pdb")
                    pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]
                combineDir= outDir+ f"/{Procedure_Name}"
                combine_simulation(combineDir, batchConfig)
                run_combined_report(combineDir, f"Config_Files/{Procedure['OperationName']}.yaml")
                    



            combine_method_report(p.join(newOutDir, "methods.md"), outDir, Procedure)
            newInDir= newOutDir+ "/00_collated_pdbs/"+ Procedure["simulationInfo"][-1]["stepName"]
            (batchConfig["pathInfo"]).update({"inputDir": newInDir})
            pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
            pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]
            

        finalReportIndex = build_operation_report_index(str(outDir))
        if finalReportIndex is not None:
            print(f"\nCreated batch full report index: {finalReportIndex}")


    ## If the config file is the standard config file then run as normal
    elif batchConfig.get('Operations') == None and batchConfig.get('simulationInfo') != None:

        ## assume if boxGeomtry and size are not specified then vacuum simulations are to be ran
        if (batchConfig['miscInfo']).get('boxGeometry') == None and (batchConfig['miscInfo']).get('boxSize') == None:
            drMDCommand = build_drmd_command("GasWing", batchConfigYaml)
            Run_drMD(drMDCommand)

        else:
            drMDCommand = build_drmd_command("SolutionWing", batchConfigYaml)
            Run_drMD(drMDCommand)
            
    
    print("\n")          


######################################################################################################
def build_drmd_command(module_name: str, config_path: str) -> list[str]:
    """
    Build a subprocess command for launching drMD with the current Python interpreter.

    Args:
        module_name (str): Name of the drMD module directory, e.g. "GasWing".
        config_path (str): Path to the config file to pass to the child process.

    Returns:
        list[str]: Command tokens for subprocess execution.
    """
    script_path = p.join(p.dirname(p.abspath(__file__)), module_name, "drMD.py")
    return [sys.executable, script_path, "--config", config_path]

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
    elif Procedure['aftercareInfo'].get('endPointInfo') == None:
        aftercareInfo = Procedure["aftercareInfo"]
        aftercareInfo.update({"endPointInfo": {"stepNames": [Procedure["simulationInfo"][-1]["stepName"]]}})
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
    
    if Procedure.get("mode") != None and Procedure.get("cutoff") != None:
        ESIDict= {"mode": Procedure["mode"], "cutoff": Procedure["cutoff"]}
        Referral.update({"ESIInfo": ESIDict})

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
    if (traj.top.select(":Cl-")).any():
        stripped_traj = traj.strip(':Cl-')
    if (traj.top.select(":Na+")).any():
        stripped_traj = traj.strip(':Na+')
    pt.write_traj(f"{PDBfile}", stripped_traj, overwrite=True)

######################################################################################################
def Run_drMD(command, print_output= True,log_output=False, check=True, *args, **kwargs):
    """
    runs drMD as a command and prints the output 

    Args:
        command(str): command to run drMD
    """
    if isinstance(command, str):
        command_parts = shlex.split(command)
    else:
        command_parts = list(command)

    logging.debug(f"Running command: {command_parts}")
    process = subprocess.Popen(  # type: ignore
        command_parts,
        shell=False,
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


def combine_simulation(outDir: FilePath, batchConfig: dict = None) -> None:
    """
    Combines the results of multiple simulations into a single output directory.

    Args:
        outDir (FilePath): The base output directory containing simulation subDirectories.
    """
    # Get all subDirectories in the output directory
    subDirs = [d for d in os.listdir(outDir) if p.isdir(p.join(outDir, d)) and not d.startswith("00_")]
    pdbDir= batchConfig["pathInfo"]["inputDir"] 
    pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"] if pdbDir else [] 
    # Create a combined output directory
    for pdbFile in pdbFiles:
        
        protDir= p.join(outDir, p.basename(pdbFile).removesuffix(".pdb"), "00_combined_simulation")
        os.makedirs(protDir, exist_ok=True)
        trajectoryFiles = []
        with open (p.join(protDir, "vitals_report.csv"), "w") as vitals_report:
            fieldNames= ['#"Step"','Time (ps)',"Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)"]
            writer=DictWriter(vitals_report, fieldnames=fieldNames)
            writer.writeheader()
            stepOffset=0
            timeOffset=0
            for subDir in subDirs:
                subDirPath = p.join(outDir, subDir)
                for root, dirs, files in os.walk(subDirPath):
                    for file in files:
                        if root.find(p.basename(pdbFile).removesuffix(".pdb")) != -1 and root.find("CM_") != -1:
                            sourceFile = p.join(root, file)
                            relativePath = p.relpath(sourceFile, subDirPath)
                            targetFile = p.join(protDir, relativePath)
                            if file.endswith("trajectory.pdb"):
                                copyfile(sourceFile,   p.join(protDir, "trajectory.pdb"))
                            elif file.endswith(".dcd"):
                                trajectoryFiles.append(sourceFile)
                            elif file == "vitals_report.csv":
                                with open(sourceFile, "r", newline="") as input_csv:
                                    reader = DictReader(input_csv)

                                    lastStep = stepOffset
                                    lastTime = timeOffset
                                    for row in reader:
                                        row['#"Step"'] = str(
                                            int(row['#"Step"']) + stepOffset
                                        )
                                        row["Time (ps)"] = str(
                                            float(row["Time (ps)"]) + timeOffset
                                        )

                                        writer.writerow(row)
                                        lastStep = int(row['#"Step"'])
                                        lastTime = float(row["Time (ps)"])
                                    stepOffset = lastStep
                                    timeOffset = lastTime

        traj= pt.iterload(trajectoryFiles, top= p.join(protDir, "trajectory.pdb"))
        traj.save(p.join(protDir, "trajectory.dcd"), overwrite=True)

def combine_simulation_ESI(outDir: FilePath, batchConfig: dict = None) -> None:
    """
    Combines the results of multiple simulations into a single output directory.

    Args:
        outDir (FilePath): The base output directory containing simulation subDirectories.
    """
    # Get all subDirectories in the output directory
    subDirs = [d for d in os.listdir(outDir) if p.isdir(p.join(outDir, d)) and not d.startswith("00_")]
    pdbDir= batchConfig["pathInfo"]["inputDir"] 
    pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"] if pdbDir else [] 
    # Create a combined output directory
    for pdbFile in pdbFiles:
        
        protDir= p.join(outDir, p.basename(pdbFile).removesuffix(".pdb"), "00_combined_simulation")
        os.makedirs(protDir, exist_ok=True)
        trajectoryFiles = []
        topologyFiles = []
        with open (p.join(protDir, "vitals_report.csv"), "w") as vitals_report:
            fieldNames= ['#"Step"','Time (ps)',"Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)"]
            writer=DictWriter(vitals_report, fieldnames=fieldNames)
            writer.writeheader()
            stepOffset=0
            timeOffset=0
            for subDir in subDirs:
                subDirPath = p.join(outDir, subDir)
                for root, dirs, files in os.walk(subDirPath):
                    for file in files:
                        if root.find(p.basename(pdbFile).removesuffix(".pdb")) != -1 and root.find("ESI_") != -1 and root.find("00_collated_pdbs") == -1:
                            sourceFile = p.join(root, file)
                            relativePath = p.relpath(sourceFile, subDirPath)
                            targetFile = p.join(protDir, relativePath)
                            if file.endswith("trajectory.pdb"):
                                if root.find("ESI_1") != -1:
                                    copyfile(sourceFile,   p.join(protDir, "trajectory.pdb"))
                                else:
                                    topologyFiles.append(sourceFile)                   
                            elif file.endswith(".dcd"):
                                if root.find("ESI_1") != -1:
                                    initialTrajectoryFile= sourceFile
                                else:
                                    trajectoryFiles.append(sourceFile)
                            elif file == "vitals_report.csv":
                                with open(sourceFile, "r", newline="") as input_csv:
                                    reader = DictReader(input_csv)

                                    lastStep = stepOffset
                                    lastTime = timeOffset
                                    for row in reader:
                                        row['#"Step"'] = str(
                                            int(row['#"Step"']) + stepOffset
                                        )
                                        row["Time (ps)"] = str(
                                            float(row["Time (ps)"]) + timeOffset
                                        )

                                        writer.writerow(row)
                                        lastStep = int(row['#"Step"'])
                                        lastTime = float(row["Time (ps)"])
                                    stepOffset = lastStep
                                    timeOffset = lastTime
        initialUniverse= mda.Universe(p.join(protDir, "trajectory.pdb"), initialTrajectoryFile)
        masterUniverse= mda.Merge(initialUniverse.atoms)
        totalAtoms= initialUniverse.atoms.n_atoms
        with mda.Writer(p.join(protDir, "trajectory.dcd"), masterUniverse.atoms.n_atoms) as W:
            for frames in initialUniverse.trajectory:
                masterUniverse.atoms.positions= initialUniverse.atoms.positions
                W.write(masterUniverse.atoms)
            for topologyFile in topologyFiles:
                ESIcountTop= topologyFile.split("ESI_")[-1].split("/")[0]
                for trajectoryFile in trajectoryFiles:
                    ESIcountTraj= topologyFile.split("ESI_")[-1].split("/")[0]
                    if ESIcountTop == ESIcountTraj:
                        universe= mda.Universe(topologyFile, trajectoryFile)
                        smallAtomCount= universe.atoms.n_atoms
                        if smallAtomCount < totalAtoms:
                            for ts in universe.trajectory:
                                framePosition= np.zeros((totalAtoms, 3))
                                framePosition[:smallAtomCount, :] = universe.atoms.positions
                                framePosition[smallAtomCount:, :] = 999

                masterUniverse.atoms.positions= framePosition
                W.write(masterUniverse.atoms)

def combine_method_report(methodReport:str, outDir:str, operationConfig: dict) -> None:
    combinedReport= p.join(outDir, "Full_Methods_Report.md")
    operationName= operationConfig.get("OperationName", "Unknown Operation")
    if operationConfig.get("type") == "solution":
        version= "SolutionWing"
    elif operationConfig.get("type") == "vacuum" or operationConfig.get("type") == "vacuumCM":
        version= "GasWing"
    elif operationConfig.get("type") == "ESI":
        version= "ESIWing"
    else:
        version= "Unknown Version"
    with open(combinedReport, "a") as combinedFile:
        combinedFile.write(f"# Methods Report for {operationName}\n\n")
        combinedFile.write(f"The following operation was run using the {version} version of drMD.\n\n")
        with open(methodReport, "r") as methodFile:
            combinedFile.write(methodFile.read())
            combinedFile.write("\n\n")

        


    
######################################################################################################
def run_combined_report(combined_dir: str, batch_config: dict) -> None:
    report_script = p.join(
        p.dirname(p.abspath(__file__)),
        "GasWing",
        "drCombinedReport.py",
    )

    Run_drMD(
        [
            sys.executable,
            report_script,
            combined_dir,
            "--config",
            batch_config
        ]
    )
######################################################################################################

if __name__ == "__main__":
    main()