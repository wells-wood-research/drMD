## BASIC PYTHON LIBRARIES
import os
from os import path as p
import numpy as np
import yaml
import itertools

## ERROR HANDLING ##
import traceback
import inspect

## PARALLELISATION LIBRARIES
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from concurrent.futures.process import BrokenProcessPool

## CUSTOM DR MD MODULES
from Triage import drConfigTriage, drPdbTriage, drConfigWriter
from Surgery import drOperator
from ExaminationRoom import  drCleanup, drLogger
from UtilitiesCloset import drSplash, drMethodsWriter

## CLEAN CODE
from typing import Optional, Dict, Tuple
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

######################################################################################################
def main(batchConfigYaml: Optional[FilePath] = None) -> None:
    '''
    Main function for drMD
    Unpacks batchConfig dictionary 
    Based on desired CPU useage, decides to use multiprocessing or not
    Based on desired CPU useage, manages CPU useage per run

    Args:
        Nothing
    Returns:
        Nothing
    '''
    ## print drMD logo
    drSplash.print_drMD_logo()

    ## if run from command line, use argpass to get batch config file
    if __name__ == "__main__":
        batchConfigYaml: FilePath = drConfigTriage.get_config_input_arg()
    ## read bacth config file into a dictionary
    try:
        batchConfig: dict = drConfigTriage.read_input_yaml(batchConfigYaml)
    except (FileNotFoundError, yaml.YAMLError, KeyError, TypeError, ValueError) as e:
        drSplash.print_config_error(e)

    batchConfig  = drConfigTriage.validate_config(batchConfig)

    ## unpack batchConfig into variables for this function
    outDir: DirectoryPath = batchConfig["pathInfo"]["outputDir"]
    yamlDir: DirectoryPath = p.join(outDir,"00_configs")
    pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
    parallelCPU: int = batchConfig["hardwareInfo"]["parallelCPU"]
    subprocessCpus: int = batchConfig["hardwareInfo"]["subprocessCpus"]

    ## create logDir if it doesn't exist
    logDir: DirectoryPath = p.join(outDir, "00_drMD_logs")
    os.makedirs(logDir, exist_ok=True)

    skipPdbTriage = batchConfig["miscInfo"].get("skipPdbTriage", False)
    if not skipPdbTriage:
        ## run pdbTriage to detect commmon problems with pdb files
        pdbTriageLog = p.join(logDir,"pdb_triage.log")
        if not p.exists(pdbTriageLog):
            drPdbTriage.pdb_triage(pdbDir, batchConfig)

    ## set environment variables for OpenMP and OpenMM - this should limit their CPU useage
    manage_cpu_usage_for_subprocesses("ON",subprocessCpus)

    ## create yamlDir if it doesn't exist, this will be used to store per-run yaml files
    os.makedirs(yamlDir,exist_ok=True)
    ## run simulations in serial or paralell
    if parallelCPU == 1:
        run_serial(batchConfig)
    elif parallelCPU > 1:
        run_parallel(batchConfig)

    ## write a methods section if desired
    writeMyMethodsSection = batchConfig["miscInfo"].get("writeMyMethodsSection", False)

    ## set up logging for post simulation processes
    drLogger.setup_logging(p.join(batchConfig["pathInfo"]["outputDir"], "00_drMD_logs", "aftercare.log"))
    ## write methods section if desired
    if writeMyMethodsSection:
        try:
            drMethodsWriter.methods_writer_protocol(batchConfig, yamlDir, outDir)
        except Exception as e:
            drLogger.log_info(f"Error writing my methods section: {e}", True, True)
    ## perform post simulation operations
    drCleanup.clean_up_handler(batchConfig)

    drLogger.log_info("Simulations Complete!", True)
    ## close logging for post simulation processes
    drLogger.close_logging()

    ## unset envorment variables for OpenMP and OpenMM
    manage_cpu_usage_for_subprocesses("OFF")

######################################################################################################
def manage_cpu_usage_for_subprocesses(mode: str, subprocessCpus: Optional[int] = None) -> None:
    '''
    In ON mode, sets thread usage for OpenMP and OpenMM 
    In OFF mode, unsets thread usage

    Args:
        mode (string): "ON" or "OFF"
        subprocessCpus (int): will set thread usage if mode == "ON"
    Returns:
        Nothing
    '''
    if mode == "ON":
        if subprocessCpus is not None:
            os.environ['OMP_NUM_THREADS'] = str(subprocessCpus)
            os.environ['OPENMM_CPU_THREADS'] = str(subprocessCpus)
        else:
            raise ValueError("subprocessCpus must be provided when mode is 'ON'")
    elif mode == "OFF":
        os.environ.pop('OMP_NUM_THREADS', None)
        os.environ.pop('OPENMM_CPU_THREADS', None)
    else:
        raise ValueError("mode must be 'ON' or 'OFF'")



###################################################################################################### 
def run_serial(batchConfig: Dict) -> None:
    """
    Process each PDB file in the given directory serially.

    Args:
        batchConfig (dict): Batch configuration dictionary.
        pdbDir (str): Path to the directory containing PDB files.
        outDir (str): Path to the output directory.
        yamlDir (str): Path to the directory to write YAML configuration files.
        simInfo (dict): Simulation information dictionary.

    Returns:
        None
    """
    botchedSimulations = []
    ## unpack batchConfig to get pdbDir
    pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
    ## create a list of PDB files
    pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]
    # Iterate over each file in the PDB directory
    for pdbFile in pdbFiles:
        # Process the PDB file
        runConfigYaml: FilePath = drConfigWriter.make_per_protein_config(pdbFile, batchConfig)
        pdbName = p.splitext(p.basename(pdbFile))[0]
        try:
            drOperator.drMD_protocol(runConfigYaml)
        except Exception as e:
            errorData = handle_exceptions(e, pdbName)
            botchedSimulations.append(errorData)

    if len(botchedSimulations) > 0:
         drSplash.print_botched(botchedSimulations)
######################################################################################################
def handle_exceptions(e, pdbName):
    tb = traceback.extract_tb(e.__traceback__)
    if tb:
        tb.reverse()
        fullTraceBack = [f"{frame.filename}:{frame.lineno} in {frame.name}" for frame in tb]
        last_frame = tb[-1]
        functionName = last_frame.name
        lineNumber = last_frame.lineno
        lineOfCode = last_frame.line
        scriptName = last_frame.filename
    else:
        functionName = 'Unknown'
        lineNumber = 'Unknown'
        lineOfCode = 'Unknown'
    
    errorType = type(e).__name__
    errorData: dict = {
        "pdbName": pdbName,
        "errorType": errorType,
        "errorMessage": str(e),
        "functionName": functionName,
        "lineNumber": lineNumber,
        "lineOfCode": lineOfCode,
        "scriptName": scriptName,
        "fullTraceBack": fullTraceBack
    }
    return errorData


######################################################################################################
def run_parallel(batchConfig: Dict) -> None:
    """
    Process each PDB file in the given directory in parallel using multiple worker threads.

    Args:
        parallelCPU (int): Number of worker threads to use.
        batchConfig (dict): Batch configuration dictionary.
        pdbDir (str): Path to the directory containing PDB files.
        outDir (str): Path to the output directory.
        yamlDir (str): Path to the directory to write YAML configuration files.
        simInfo (dict): Simulation information dictionary.

    Returns:
        None
    """
    
    ## read input directory from batchConfig
    pdbDir = batchConfig["pathInfo"]["inputDir"]
    parallelCpus: int = batchConfig["hardwareInfo"]["parallelCPU"]

    # Get list of PDB files in the directory
    pdbFiles: list[str] = sorted([p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"])
    ## construct inputArgs for multiprocessing
    inputArgs: list[tuple] = [(pdbFile, batchConfig) for pdbFile in pdbFiles]
    ## create batched inputs
    batchedArgsWithPos = [(batch, pos) for pos, batch in enumerate(np.array_split(inputArgs, parallelCpus))]
    ## add a dummy batch to be used for printing logging
    batchedArgsWithPos = [(["dummy"], -1)] + batchedArgsWithPos

    try:
        ## run simulations in parallel
        botchedSimulations = process_map(per_core_worker, batchedArgsWithPos, 
                    max_workers=parallelCpus)
    except BrokenProcessPool:
        print("BrokenProcessPool: Terminating remaining processes")

    botchedSimulations = list(itertools.chain.from_iterable(botchedSimulations))

    if len(botchedSimulations) > 0:
         drSplash.print_botched(botchedSimulations)
######################################################################################################
def per_core_worker(batchedArgsWithPos: Tuple[Dict, int]) -> None:
    """
    Each core is passed a batch of arguments to process
    This function unpacks the arguments and handles the progress bar

    Args:
        batchedArgsWithPos List[(Tuple[Dict, int])]: 
            Alist of tuples containing input argunents for process_pdb_file and 
            the position of the core for the loading bar
    Returns:
        None
    """
    perWorkerBotchedSimulations: list[Dict] = []
    ## unpack batchedArgsWithPos into the batch of arguments for 
    batchedArgs, pos = batchedArgsWithPos

    ## create a list of colors for the loading bar
    cmap = plt.get_cmap('plasma', 32)
    colors = [mcolors.rgb2hex(cmap(i)) for i in range(32)]
    ## create a dummy progress bar to be used for printing logs
    if pos == -1:
        with tqdm(total=1, position=0, bar_format='{desc}', 
                  colour="#000000", leave=True) as dummy_progress:
            dummy_progress.set_description_str("Logs:")
            dummy_progress.refresh()
    ## run simulations in parallel
    else:
        with tqdm(desc=f"Core {str(pos)}", total=len(batchedArgs), 
                position=pos+1, colour=colors[pos % len(colors)], 
                leave=False) as progress:
            for args in batchedArgs:
                pdbFile, batchConfig = args
                pdbName = p.splitext(p.basename(pdbFile))[0]
                runConfigYaml: FilePath = drConfigWriter.make_per_protein_config(pdbFile, batchConfig)
                try:
                    drOperator.drMD_protocol(runConfigYaml)
        
                except Exception as e:
                    errorData = handle_exceptions(e, pdbName)
                    perWorkerBotchedSimulations.append(errorData)
                    continue
                progress.update(1)
            progress.close()  
    return perWorkerBotchedSimulations
######################################################################################################

if __name__ == "__main__":
    main()
