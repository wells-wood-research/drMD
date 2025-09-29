## BASIC PYTHON LIBRARIES
import argpass
import yaml
import os
from os import path as p
from pathlib import Path
## PARALLELISATION LIBRARIES
import multiprocessing as mp

## drMD LIBRARIES
from ExaminationRoom import drLogger
from UtilitiesCloset import drSplash

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

## CLEAN CODE
from typing import Tuple, Union, List, Dict
from os import PathLike
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath


#####################################################################################
def validate_config(config: dict) -> FilePath:
    """
    Main script for drConfigInspector:
    1. Accepts --config flag arg with argpass, reads yaml file into a dict
    2. Checks paths in pathInfo to see if they are real
    3. Checks inputs of cpuInfo to see if they are correct
    4. Checks each dockingOrder in dockingOrders to see if they are correct

    Returns:
    - config (dict)
    """

    # ## set up logging
    topDir: DirectoryPath = os.getcwd()
    # configTriageLog: FilePath = p.join(topDir,"config_triage.log")
    # drLogger.setup_logging(configTriageLog)
    # configTriageLogData.append(f"Checking config file...",True)


        ## DEFAULTS ###
    ## get some defaults for config
    configDefaults: dict = init_config_defaults(topDir)
    ## merge with config to apply defaults to empty fields
    config = {**configDefaults, **config}

    ## check each major section in config file
    ## throw errors with a nice splash screen if something goes wrong
    configDisorders = {}
    config, configDisorders["miscInfo"], miscInfoOk = check_miscInfo(config, configDefaults)

    config, configDisorders["pathInfo"], pathInfoOk = check_pathInfo(config, configDefaults)
    config, configDisorders["hardwareInfo"], hardwareInfoOk  = check_hardwareInfo(config, configDefaults)

    config, configDisorders["simulationInfo"], simulationInfoOk = check_simulationInfo(config)

    configDisorders["aftercareInfo"], aftercareInfoOk  = check_aftercareInfo(config)
  
    configDisorders["ligandInfo"], ligandInfoOk  = check_ligandInfo(config)
                

    allInfoOk = miscInfoOk * pathInfoOk * hardwareInfoOk * simulationInfoOk * aftercareInfoOk * ligandInfoOk

    ## TODO: write configDisorders to file in all cases
    if allInfoOk:
        return config
    else:
        drSplash.print_config_error(configDisorders)
       

#####################################################################################
def init_config_defaults(topDir: DirectoryPath) -> dict:
    configDefaults = {
        "pathInfo": {
            "inputDir": topDir,
            "outputDir": p.join(topDir, "outputs")
        },
        "hardwareInfo": {
            "parallelCPU": 1,
            "platform": "CPU",
            "subprocessCpus": 1
        },

        "miscInfo": {   
            "pH": 7,
            "firstAidMaxRetries": 10,
            "boxGeometry": "cubic",
            "boxSize": 10,
            "writeMyMethodsSection": True,
            "skipPdbTriage": False,
            "trajectorySelections": [{"selection": {"keyword": 'all'}}]
            },

        "simInfo": {
            "stepName": None,           ## THIS IS REQUIRED
            "simulationType": "NPT",
            "temperature": 300,
            "temperatureRange": None,   ## OPTIONAL, NO DEFAULT
            "maxIterations": -1,   ## ONLY TURN ON WHEN EM STEP  
            "duration": None,       ## THIS IS REQUIRED
            "timestep": "2 fs",           ## IF HEAVY PROTONS IS SET TO TRUE, USE "4 fs"  INSTEAD
            "logInterval": "10 ps",
            "heavyProtons": False,   ## 
        }
    }
    return configDefaults


#####################################################################################
def check_pathInfo(config: dict, configDefaults: dict) -> Tuple[dict, dict, bool]:
    """
    Checks for pathInfo entry in config
    Checks paths in pathInfo to see if they are real
    Don't check outputDir, this will be made automatically
    """
    pathInfoOk = True
    pathInfoDisorders = {}
    ## log this check

    ## check to see if pathInfo in config
    pathInfo = config.get("pathInfo", None)
    if pathInfo is None:
        config["pathInfo"] = configDefaults["pathInfo"]
        pathInfoDisorders["inputDir"] = "Automatic Default Used!"
        pathInfoDisorders["outputDir"] = "Automatic Default Used!"

        return config, pathInfoDisorders, False

    ## validate inputDir
    inputDir = pathInfo.get("inputDir", None)
    if inputDir is None:
        config["pathInfo"]["inputDir"] = configDefaults["pathInfo"]["inputDir"]
        pathInfoDisorders["inputDir"] = "Automatic Default Used!"
        pathInfoOk = False
    else:
        inputDirPathProblem = validate_path(f"inputDir", inputDir)
        if inputDirPathProblem is not None:
            pathInfoDisorders["inputDir"] = inputDirPathProblem
            pathInfoOk = False
        else:
            pathInfoDisorders["inputDir"] = None

    ## validate outputDir
    outputDir = pathInfo.get("outputDir", None)
    if outputDir is None:
        config["pathInfo"]["outputDir"] = configDefaults["pathInfo"]["outputDir"]
        pathInfoDisorders["outputDir"] = "Automatic Default Used!"
        pathInfoOk = False
    else:
        if not isinstance(Path(outputDir), PathLike):
            pathInfoDisorders["outputDir"] = "outputDir must be a PathLike string!"
            pathInfoOk = False
        else:
            pathInfoDisorders["outputDir"] = None

    return config, pathInfoDisorders, pathInfoOk




#####################################################################################
def check_hardwareInfo(config: dict, configDefaults: dict) -> Tuple[dict, dict, bool]:
    """
    Checks hardwareInfo in config
    Makes sure that CPU allocations are properly formatted
    Makes sure that the "Platform" specified is an allowed value 
    """
    haredwareInfoOk = True
    hardwareInfoDisorders = {}
    ## check if hardwareInfo in config
    hardwareInfo = config.get("hardwareInfo", None)
    if hardwareInfo is None:
        config["hardwareInfo"] = configDefaults["hardwareInfo"]
        for argName in ["parallelCPU", "platform", "subprocessCpus"]:
            config["hardwareInfo"][argName] = configDefaults["hardwareInfo"][argName]
            hardwareInfoDisorders[argName] = "Automatic Default Used!"
        return config, hardwareInfoDisorders, True

    ## validate parallelCPU
    parallelCPU = hardwareInfo.get("parallelCPU", None)
    if parallelCPU is None:
        ## use a default value
        config["hardwareInfo"]["parallelCPU"] = configDefaults["hardwareInfo"]["parallelCPU"]
        hardwareInfoDisorders["parallelCPU"] = "Automatic Default Used!"
    else:
        if not isinstance(parallelCPU, int):
            hardwareInfoDisorders["parallelCPU"] = "parallelCPU is not an int"
            haredwareInfoOk = False
        else:
            if parallelCPU < 1:
                hardwareInfoDisorders["parallelCPU"] = "parallelCPU is less than 1, this must be a positive integer"
                haredwareInfoOk = False
            else:
                hardwareInfoDisorders["parallelCPU"] = None


    ## validate subprocessCpus
    subprocessCpus = hardwareInfo.get("subprocessCpus", None)
    if subprocessCpus is None:
        # use a default value
        config["hardwareInfo"]["subprocessCpus"] = configDefaults["hardwareInfo"]["subprocessCpus"]
        hardwareInfoDisorders["subprocessCpus"] = "Automatic Default Used!"
    else:
        if not isinstance(subprocessCpus, int):
            hardwareInfoDisorders["subprocessCpus"] = "subprocessCpus is not an int"
            haredwareInfoOk = False
        else:
            if subprocessCpus < 1:
                hardwareInfoDisorders["subprocessCpus"] = "subprocessCpus is less than 1, this must be a positive integer"
                haredwareInfoOk = False
            else:
                hardwareInfoDisorders["subprocessCpus"] = None

    ## check to see if the number of CPU cores requested is less than the number of CPU cores available
    if isinstance(parallelCPU,int) and isinstance(subprocessCpus,int):
        systemCpus = mp.cpu_count()
        if parallelCPU * subprocessCpus > systemCpus:
            hardwareInfoDisorders["totalCpuUseage"] = "Number for CPU cores requested exceeds number of CPU cores available, change the values of parallelCPU and subprocessCpus"
            haredwareInfoOk = False
    ## validate platform
    platform = hardwareInfo.get("platform", None)
    if platform is None:
        ## use a default value
        config["hardwareInfo"]["platform"] = configDefaults["hardwareInfo"]["platform"]
        hardwareInfoDisorders["platform"] = "Automatic Default Used!"
    else:
        if platform not in ["CUDA", "OPENCL", "CPU"]:
            hardwareInfoDisorders["platform"] = "platform is not 'CUDA', 'OPENCL', or 'CPU'"
            haredwareInfoOk = False
        else:
            hardwareInfoDisorders["platform"] = None


    return config, hardwareInfoDisorders, haredwareInfoOk



#####################################################################################
def check_miscInfo(config:dict, configDefaults:dict) -> Tuple[dict,dict,bool]:
    miscInfoOk = True
    miscInfoDisorders = {}

    miscInfo = config.get("miscInfo", None)
    ## if no miscInfo in config, use defaults and report disorders
    if miscInfo is None:
        config["miscInfo"] = configDefaults["miscInfo"]
        for argName in ["pH",
                         "firstAidMaxRetries",
                           "boxGeometry",
                             "writeMyMethodsSection",
                               "skipPdbTriage",
                                 "trajectorySelections"]:
            miscInfoDisorders[argName] = "Automatic Default Used!"

        return config, miscInfoDisorders, True
    
    ## validate pH
    pH = miscInfo.get("pH", None)
    if pH is None:
        ## use a default value
        config["miscInfo"]["pH"] = configDefaults["miscInfo"]["pH"]
        miscInfoDisorders["pH"] = "No pH specified, using default of 7"
    else:
        if not isinstance(pH, (int, float)):
            miscInfoDisorders["pH"] = "pH must be an int or float between 0 and 14"
            miscInfoOk = False
        else:
            if pH < 0 or pH > 14:
                miscInfoDisorders["pH"] = "pH must be an int or float between 0 and 14"
                miscInfoOk = False
            else:
                miscInfoDisorders["pH"] = None

    ## validate firstAidMaxRetries
    firstAidMaxRetries = miscInfo.get("firstAidMaxRetries", None)
    if firstAidMaxRetries is None:
        ## use a default value
        config["miscInfo"]["firstAidMaxRetries"] = configDefaults["miscInfo"]["firstAidMaxRetries"]
        miscInfoDisorders["firstAidMaxRetries"] = "No firstAidMaxRetries specified, using default of 10"
    else:
        if not isinstance(firstAidMaxRetries, int):
            miscInfoDisorders["firstAidMaxRetries"] = "firstAidMaxRetries must be an int greater than 0"
            miscInfoOk = False
        else:
            if firstAidMaxRetries < 0:
                miscInfoDisorders["firstAidMaxRetries"] = "firstAidMaxRetries must be an int greater than or equal to 0"
                miscInfoOk = False
            else:
                miscInfoDisorders["firstAidMaxRetries"] = None

    ## validate boxGeometry
    boxGeometry = miscInfo.get("boxGeometry", None)
    if boxGeometry is None:
        ## use a default value
        config["miscInfo"]["boxGeometry"] = configDefaults["miscInfo"]["boxGeometry"]
        miscInfoDisorders["boxGeometry"] = "No boxGeometry specified, using default of 'cubic'"
    else:
        if boxGeometry not in ["cubic", "octahedral"]:
            miscInfoDisorders["boxGeometry"] = "boxGeometry must be either 'cubic' or 'octahedral'"
            miscInfoOk = False
        else:
            miscInfoDisorders["boxGeometry"] = None

    ## validate boxSize
    boxSize = miscInfo.get("boxSize", None)
    if boxSize is None:
        ## use a default value
        config["miscInfo"]["boxSize"] = configDefaults["miscInfo"]["boxSize"]
        miscInfoDisorders["boxSize"] = "No boxSize specified, using default of None"
    else:
        if not isinstance(boxSize, int):
            miscInfoDisorders["boxSize"] = "boxSize must be an int"
            miscInfoOk = False
        else:
            miscInfoDisorders["boxSize"] = None

    ## validate skipPdbTriage
    skipPdbTriage = miscInfo.get("skipPdbTriage", None)
    if skipPdbTriage is None:
        ## use a default value
        config["miscInfo"]["skipPdbTriage"] = configDefaults["miscInfo"]["skipPdbTriage"]
        miscInfoDisorders["skipPdbTriage"] = "No skipPdbTriage specified, using default of False"
    else:
        if not isinstance(skipPdbTriage, bool):
            miscInfoDisorders["boxGeometry"] = "skipPdbTriage must be either True or False"
            miscInfoOk = False
        else:
            miscInfoDisorders["skipPdbTriage"] = None
    
    ## validate writeMyMethodsSection
    writeMyMethodsSection = miscInfo.get("writeMyMethodsSection", None)
    if writeMyMethodsSection is None:
        ## use a default value
        config["miscInfo"]["writeMyMethodsSection"] = configDefaults["miscInfo"]["writeMyMethodsSection"]
        miscInfoDisorders["writeMyMethodsSection"] = "No writeMyMethodsSection specified, using default of True"
    else:
        if not isinstance(writeMyMethodsSection, bool):
            miscInfoDisorders["writeMyMethodsSection"] = "writeMyMethodsSection must be either True or False"
            miscInfoOk = False
        else:
            miscInfoDisorders["writeMyMethodsSection"] = None

    ## validate trajectorySelections
    trajectorySelections = miscInfo.get("trajectorySelections", None)
    if trajectorySelections is None:
        ## use a default value  
        config["miscInfo"]["trajectorySelections"] = configDefaults["miscInfo"]["trajectorySelections"]
        miscInfoDisorders["trajectorySelections"] = "No trajectorySelections specified, all atoms will be selected by default"
    else:
        if not isinstance(trajectorySelections, list):
            miscInfoDisorders["trajectorySelections"] = "trajectorySelections must be a list of selection dictionaries (see README for syntax!)"
            miscInfoOk = False
        elif len(trajectorySelections) == 0:
            miscInfoDisorders["trajectorySelections"] = "trajectorySelections must be a list of selection dictionaries (see README for syntax!)"
            miscInfoOk = False
        else:
            for selection in trajectorySelections:
                selectionDisorder = check_selection(selection)
                if len(selectionDisorder) > 0:
                    miscInfoDisorders["trajectorySelections"] = selectionDisorder
                    miscInfoOk = False
                else:
                    miscInfoDisorders["trajectorySelections"] = None

    return config, miscInfoDisorders, miscInfoOk

#####################################################################################
def check_simulationInfo(config: dict) -> Tuple[dict, bool]:
    """
    Checks for simulationInfo in config
    Depending on the type of simulation, checks your parameters
    """
    ## log this check
    simulationInfoOk = True
    simulationInfo = config.get("simulationInfo", None)


    ## check for problems with the simuationInfo dict
    if simulationInfo is None:
        return config,  "No simulationInfo found", False
    if not isinstance(simulationInfo, list):
        return config, "simulationInfo must be a list of dicts", False
    if len(simulationInfo) == 0:
        return config, "simulationInfo must have at least one entry", False
    

    simulationInfoDisorders = {}
    for counter, simulation in enumerate(simulationInfo):

        disorders = {}
        disorders, stepName, simulationType, sharedOptionsOk, simulationWithDefaults  = check_shared_simulation_options(simulation, disorders)
        simulationInfoOk *= sharedOptionsOk

        ## if we don't have stepName, further checks wont work.
        if stepName is None:
            disorders["stepName"] = "stepName must be specified"
            simulationInfoDisorders[f"Unnamed_step_{str(counter+1)}"] = "stepName must be specified"
            simulationInfo = False
            continue


        allStepsOk = True
        if simulationWithDefaults["simulationType"] in ["NVT", "NPT"]:
            disorders, mdOptionsOk, simulationWithDefaults  = check_nvt_npt_options(simulationWithDefaults, stepName, disorders)
            allStepsOk *= mdOptionsOk
        elif simulationWithDefaults["simulationType"] == "META":
            disorders, metaOptionsOk = check_metadynamics_options(simulationWithDefaults, disorders)
            allStepsOk *+ metaOptionsOk
        elif simulationWithDefaults["simulationType"] == "EM":
            disorders, emOptionsOk, simulationWithDefaults = check_em_options(simulationWithDefaults, disorders)
            allStepsOk *= emOptionsOk

        restraintsInfo = simulation.get("restraintInfo", None)
        if restraintsInfo:
            disorders, restraintInfoOk = check_restraintInfo(restraintsInfo, disorders)
            allStepsOk *= restraintInfoOk
        simulationInfoDisorders[stepName] = disorders

        simulationInfoOk *= allStepsOk

        ## update config dict to use defaults
        config["simulationInfo"][counter] = simulationWithDefaults

    return config, simulationInfoDisorders, simulationInfoOk

#################################################################################################
def check_em_options(simulation: dict,disorders: dict) -> Tuple[dict, bool]:
    emOptionsOk = True
    ## check duration
    maxIterations = simulation.get("maxIterations", None)
    if maxIterations is None:
        disorders["maxIterations"] = "maxIterations not set for EM step, using a default of -1"
        simulation["maxIterations"] = -1
    else:
        if not isinstance(maxIterations, int):
            disorders["maxIterations"] = "maxIterations must be an integer"
            emOptionsOk = False
        else:
            disorders["maxIterations"] = None

    return disorders, emOptionsOk, simulation


#################################################################################################
def check_restraintInfo(restraintInfo: dict, disorders: dict) -> Tuple[dict, bool]:

    restraintInfofOk = True
    if not isinstance(restraintInfo, list):
        disorders["restraintInfo"] = "restraintInfo must be a dictionary"
        return disorders, False
    if len(restraintInfo) == 0:
        disorders["restraintInfo"] = "restraintInfo must have at least one entry"
        return disorders, False
    
    ## check each entry in restraintInfo

    for restraintIndex, info in enumerate(restraintInfo):
        if not isinstance(info, dict):
            disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo must be a dictionary (see README for syntax!)"
            restraintInfofOk = False
        else:
            ## check restraintType
            restraintType = info.get("restraintType", None)
            if restraintType is None:
                disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo must have a 'restraintType' key"
                restraintInfofOk = False
            else:
                if not isinstance(restraintType, str):
                    disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo['restraintType'] must be a string"
                    restraintInfofOk = False
                if not restraintType in ["position", "distance", "angle", "torsion", "pulling", "centroid_bond", "torsion_flat_bottom"]:
                    disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo['restraintType'] must be one of 'position', 'distance', 'angle', 'torsion', 'pulling', 'centroid_bond', 'torsion_flat_bottom'"
                    restraintInfofOk = False
            ## check selection for restraint to act upon
            restrantSelection = info.get("selection", None)
            if  restrantSelection is None:
                disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo must have a 'selection' key"
                restraintInfofOk = False
            else:
                if not isinstance(restrantSelection, dict):
                    disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo['selection'] must be a dictionary (see README for syntax!)"
                    restraintInfofOk = False
                else:
                    selectionDisorder = check_selection({"selection": restrantSelection})
                    if len(selectionDisorder) > 0:
                        disorders["restraintInfo"][f"restraint_{restraintIndex}"] = selectionDisorder
                        restraintInfofOk = False
            ## check parameters for restraint
            restraintParamers = info.get("parameters", None)
            if restraintParamers is None:
                disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo must have a 'parameters' key"
                restraintInfofOk = False
            else:
                restraintParamProblems = check_restraint_parameters(restraintType, restraintParamers)
                if len(restraintParamProblems) > 0:
                    disorders["restraintInfo"][f"restraint_{restraintIndex}"] = restraintParamProblems
                    restraintInfofOk = False
                else:
                    disorders["restraintInfo"] = {}
                    disorders["restraintInfo"][f"restraint_{restraintIndex}"] = None

    return disorders, restraintInfofOk

    


#####################################################################################
def check_aftercareInfo(config: dict) -> Tuple[dict,bool]:
    """
    Checks for aftercareInfo in config
    """
    aftercareInfoOk = True
    aftercareInfoDisorders = {}
    ## check for aftercareInfo in config
    aftercareInfo = config.get("aftercareInfo", None)
    if aftercareInfo is None:
        return None, aftercareInfoOk
    
    endPointInfo = aftercareInfo.get("endPointInfo", None)
    if endPointInfo is not None:
        endPointInfoDisorders, endPointInfoOk = check_endPointInfo(endPointInfo)
        aftercareInfoOk *= endPointInfoOk 
        aftercareInfoDisorders["endPointInfo"] = endPointInfoDisorders

    clusterInfo = aftercareInfo.get("clusterInfo", None)
    if clusterInfo is not None:
        aftercareInfoDisorders["clusterInfo"], clusterInfoOk = check_clusterInfo(clusterInfo)
        aftercareInfoOk *= clusterInfoOk 


    collateVitalsReports = aftercareInfo.get("collateVitalsReports", None)
    if collateVitalsReports is not None:
        if not isinstance(collateVitalsReports, bool):
            aftercareInfoDisorders["collateVitalsReports"] = "collateVitalsReports must be a boolean"
            afterCareInfoOk = False

    return aftercareInfoDisorders, aftercareInfoOk



#####################################################################################
def check_endPointInfo(endPointInfo: dict) -> Tuple[dict, bool]:
    """
    Checks for endPointInfo in config
    """
    endPointInfoOk = True
    endPointDisorders = {}
    ## check if endPointInfo is a dictionary
    if not isinstance(endPointInfo, dict):
        return "endPointInfo must be a dictionary"
    ## check stepNames
    stepNames = endPointInfo.get("stepNames", None)
    if stepNames is None:
        endPointDisorders["stepNames"] = "endPointInfo must have a 'stepNames' entry"
        endPointInfoOk = False
    else:
        if not isinstance(stepNames, list):
            endPointDisorders["stepNames"] = "stepNames must be a list"
            endPointInfoOk = False
        ## ensure that stepNames is a list of strings
        if not all(isinstance(stepName, str) for stepName in stepNames):
            endPointDisorders["stepNames"] = "stepNames must be a list of strings"
            endPointInfoOk = False
        ## ensure that stepNames is not empty
        if  len(stepNames) == 0:
            endPointDisorders["stepNames"] = "stepNames must not be empty"
            endPointInfoOk = False
        else:
            endPointDisorders["stepNames"] = None
    ## check removeAtoms
    removeAtoms = endPointInfo.get("removeAtoms", None)
    if removeAtoms is not None:
        if not isinstance(removeAtoms, list):
            endPointDisorders["removeAtoms"] = "removeAtoms must be a list"
            endPointInfoOk = False
        if not all(isinstance(removeAtom, dict) for removeAtom in removeAtoms):
            endPointDisorders["removeAtoms"] = "removeAtoms must be a list of dictionaries"
            endPointInfoOk = False
        for selection in removeAtoms:
            removeAtomSelectionDisorders = check_selection(selection)
            if len(removeAtomSelectionDisorders) > 0:
                endPointDisorders["removeAtoms"] = removeAtomSelectionDisorders
                endPointInfoOk = False
            else:
                endPointDisorders["removeAtoms"] = None
        
    return endPointDisorders, endPointInfoOk
#####################################################################################
def check_clusterInfo(clusterInfo: dict) -> Tuple[dict,bool]:
    """
    Checks for clusterInfo in config
    """
    clusterInfoOk = True
    clusterInfoDisorders = {}
    ## check if clusterInfo is a dictionary
    if not isinstance(clusterInfo, dict):
        return "clusterInfo must be a dictionary", False
    ## check stepNames
    stepNames = clusterInfo.get("stepNames", None)
    if stepNames is None:
        clusterInfoDisorders["stepNames"] = "clusterInfo must have a 'stepNames' entry"
        clusterInfoOk = False
    else:
        if not isinstance(stepNames, list):
            clusterInfoDisorders["stepNames"] = "clusterInfo['stepNames'] must be a list"
            clusterInfoOk = False
        ## ensure that stepNames is a list of strings
        elif not all(isinstance(stepName, str) for stepName in stepNames):
            clusterInfoDisorders["stepNames"] = "clusterInfo['stepNames'] must be a list of strings"
            clusterInfoOk = False
        ## ensure that stepNames is not empty
        elif  len(stepNames) == 0:
            clusterInfoDisorders["stepNames"] = "clusterInfo['stepNames'] must not be empty"
            clusterInfoOk = False
        else:
            clusterInfoDisorders["stepNames"] = None

    ## check nClusters
    nClusters = clusterInfo.get("nClusters", None)
    if nClusters is None:
        clusterInfoDisorders["nClusters"] = "nClusters must be specified as a positive int"
        clusterInfoOk = False
    else:
        if not isinstance(nClusters, int):
            clusterInfoDisorders["nClusters"] = "nClusters must be an int"
            clusterInfoOk = False
        elif nClusters < 1:
            clusterInfoDisorders["nClusters"] = "nClusters must be specified as a positive int"
            clusterInfoOk = False
        else:
            clusterInfoDisorders["nClusters"] = None

    ## check clusterBy
    clusterBy = clusterInfo.get("clusterBy", None)
    if clusterBy is None:
        clusterInfoDisorders["clusterBy"] = "clusterInfo must have a 'clusterBy' entry"
        clusterInfoOk = False
    else:
        if not isinstance(clusterBy, list):
            clusterInfoDisorders["clusterBy"] = "clusterBy must be a list"
            clusterInfoOk = False
        elif not all(isinstance(clusterSelection, dict) for clusterSelection in clusterBy):
            clusterInfoDisorders["clusterBy"] = "clusterBy must be a list of dictionaries"
            clusterInfoOk = False
        elif len(clusterBy) == 0:
            clusterInfoDisorders["clusterBy"] = "clusterBy must not be empty"
            clusterInfoOk = False
        for clusterSelection in clusterBy:
            clusterSelectionDisorders = check_selection(clusterSelection)
            if len(clusterSelectionDisorders) > 0:
                clusterInfoDisorders["clusterBy"] = clusterSelectionDisorders
                clusterInfoOk = False
            else:
                clusterInfoDisorders["clusterBy"] = None
    ## check removeAtoms
    removeAtoms = clusterInfo.get("removeAtoms", None)
    if removeAtoms is not None:
        if not isinstance(removeAtoms, list):
            clusterInfoDisorders["removeAtoms"] = "endPointInfo['removeAtoms'] must be a list"
            clusterInfoOk = False
        elif not all(isinstance(removeAtom, dict) for removeAtom in removeAtoms):
            clusterInfoDisorders["removeAtoms"] = "endPointInfo['removeAtoms'] must be a list of dictionaries"
            clusterInfoOk = False
        for selection in removeAtoms:
            removeAtomSelectionDisorders = check_selection(selection)
            if len(removeAtomSelectionDisorders) > 0:
                clusterInfoDisorders["removeAtoms"] = removeAtomSelectionDisorders
                clusterInfoOk = False
            else:
                clusterInfoDisorders["removeAtoms"] = None

    return clusterInfoDisorders, clusterInfoOk


#########################################################################
def check_ligandInfo(config: dict) -> Tuple[dict, bool]:
    """
    Checks optional ligandInfo entry in config

    Args:
        config (dict): The main configuration dictionary

    Raises:
        TypeError: If ligandInfo is not a list of dictionaries, or if ligandName is not a string
        TypeError: If protons, charge, frcmod, or mol2 is not a boolean
        ValueError: If ligandInfo is an empty list
        ValueError: If ligandInfo does not have at least one entry
    """
    ligandInfoOk = True
    # Check if ligandInfo in config
    ligandInfo = config.get("ligandInfo", None)
    ## if there is no ligandInfo specified, not a problem, return None
    if ligandInfo is None:
        return None, ligandInfoOk
    elif len(ligandInfo) == 0:
        ligandInfoOk = False
        return "ligandInfo must be a contain at least one ligand dictionary (see README for more info)", ligandInfoOk
    
    ligandInfoDisorders = {}
    # Check each entry in ligandInfo
    for ligandIndex, ligand in enumerate(ligandInfo):
        # Check if ligand is a dictionary
        if not isinstance(ligand, dict):
            ligandInfoDisorders[f"ligand_{ligandIndex}"] = "ligand entry must be a dictionary"
            ligandInfoOk = False 
            continue
        ## check ligandName
        ligandName = ligand.get("ligandName", None)
        if ligandName is None:
            ligandInfoDisorders[f"ligand_{ligandIndex}"]["ligandName"] = "each ligand entry must have a ligandName entry as a unique string"
            ligandInfoOk = False
            ## set a temporary ligand name for reporting disorders
            ligandName = f"ligand_{ligandIndex}"
        else:
            ligandInfoDisorders[ligandName] = {}
        if not isinstance(ligandName, str):
            ligandInfoDisorders[ligandName]["ligandName"] = "each ligand entry must have a ligandName entry as a unique string"
            ligandInfoOk = False
        else:
            ligandInfoDisorders[ligandName]["ligandName"] = None

        ## check boolean flags
        for argName in ["protons", "frcmod", "mol2"]:
            argValue = ligand.get(argName, None)
            if argValue is None:
                ligandInfoDisorders[ligandName][argName] = f"each ligand entry must have a {argName} entry as a bool"
                ligandInfoOk = False
            elif not isinstance(argValue, bool):
                ligandInfoDisorders[ligandName][argName] = f"{argName} must be a bool"
                ligandInfoOk = False
            else:
                ligandInfoDisorders[ligandName][argName] = None
        ## check charge
        charge = ligand.get("charge", None)
        if charge is None:
            ligandInfoDisorders[ligandName]["charge"] = "each ligand entry must have a charge entry as an int"
            ligandInfoOk = False
        elif not isinstance(charge, int):
            ligandInfoDisorders[ligandName]["charge"] = "charge must be an integer"
            ligandInfoOk = False
        else:
            ligandInfoDisorders[ligandName]["charge"] = None
 
    return ligandInfoDisorders, ligandInfoOk






#########################################################################
def check_restraint_parameters(restraintType: str, parameters: dict) -> None:

    parameterProblems = []
    # Check for force constants
    forceConstant = parameters.get("k", None)
    if not forceConstant:
        parameterProblems.append("No force constat k provided for restraints")
    else:
        if not isinstance(forceConstant, (int, float)):
            parameterProblems.append("force constants must be a positive number")
        elif forceConstant <= 0:
            parameterProblems.append("force constants must be a positive number")

    ## deal with phi0 in torsions
    if restraintType.upper() == "TORSION":
        phi0 = parameters.get("phi0", None)
        if  phi0 is None:
            parameterProblems.append("phi0 parameter must be provided for torsion restraints")
        else:
            if not isinstance(phi0, (int, float)):
                parameterProblems.append("phi0 parameter must be a number for torsion restraints")
            if phi0 < -180 or phi0 > 180:
                parameterProblems.append("phi0 parameter must be between -180 and 180 for torsion restraints")
        
    elif restraintType.upper() == "DISTANCE":
        r0 = parameters.get("r0", None)
        if  r0 is None:
            parameterProblems.append("r0 parameter must be provided for distance restraints")
        else:
            if not isinstance(r0, (int, float)):
                parameterProblems.append("r0 parameter must be a number for distance restraints")
            if r0 < 0:
                parameterProblems.append("r0 parameter must be positive for distance restraints")

    elif restraintType.upper() == "ANGLE":
        theta0 = parameters.get("theta0", None)
        if  theta0 is None:
            parameterProblems.append("theta0 parameter must be provided for angle restraints")
        else:
            if not isinstance(theta0, (int, float)):
                parameterProblems.append("theta0 parameter must be a number for angle restraints")
            if theta0 < 0 or theta0 > 360:
                parameterProblems.append("theta0 parameter must be between 0 and 360 for angle restraints")


    return parameterProblems

#########################################################################
def check_selection(selection: dict) -> list:
    selectionDisorders = []
    ## check keyword
    subSelection = selection.get("selection", None)
    if subSelection is None:
        return ["No selection specified in selection"]

    keyword = subSelection.get("keyword", None)
    if keyword is None:
        return ["No keyword specified in selection"]

    if not isinstance(keyword, str):
        return [f"keyword must be a string, not {type(keyword)}"]
    
    if not keyword in ["all", "protein", "ligand", "water", "ions", "custom", "backbone", "custom_centroid"]:
        return [f"selection keywords incorrect see README.md for more details"]
    
    ## check custom selection syntax
    if keyword == "custom":
        customSelection = subSelection.get("customSelection", None)
        if customSelection is None:
            selectionDisorders.append("No customSelction specified in selection")
            return selectionDisorders
        if not isinstance(customSelection, list):
            selectionDisorders.append("customSelection must be a list of selection dictionaries (see README for more details)")
            return selectionDisorders
        ## check each selection
        for customSelctionDict in customSelection:
            ## check chainId
            chainId = customSelctionDict.get("CHAIN_ID", None)
            if chainId is None:
                selectionDisorders.append("No CHAIN_ID specified in selection")
            ## check resName
            resName = customSelctionDict.get("RES_NAME", None)
            if resName is None:
                selectionDisorders.append("No RES_NAME specified in selection")
            elif resName == "all":
                pass
            elif not isinstance(resName, str):
                selectionDisorders.append("RES_NAME must be a three-letter string")
            elif len(resName) > 3:
                selectionDisorders.append("RES_NAME must be an up-to three-letter string")
            ## check resId
            resId = customSelctionDict.get("RES_ID", None)
            if resId is None:
                selectionDisorders.append("No RES_ID specified in selection")
            elif resId == "_":
                pass
            elif not isinstance(resId, int):
                selectionDisorders.append("RES_ID must be an integer")
            ## check atomName
            atomName = customSelctionDict.get("ATOM_NAME", None)
            if atomName is None:
                selectionDisorders.append("No ATOM_NAME specified in selection")
            elif not isinstance(atomName, str):
                selectionDisorders.append("ATOM_NAME must be a string")
            elif len(atomName) > 4:
                selectionDisorders.append("ATOM_NAME must be a string less than 4 characters")

        ## todo: check custom_centroid selection syntax

    return selectionDisorders


#########################################################################
def check_metadynamics_options(simulation: dict, disorders: dict) -> Tuple[dict,bool]:
    metaOptionsOk = True
    ## check metaDynamicsInfo
    metaDynamicsInfo = simulation.get("metaDynamicsInfo", None)
    if metaDynamicsInfo is None:
        disorders["metaDynamicsInfo"] = "No metadynamicsInfo found in step with simluationType of META"
        return disorders, False
    if not isinstance(metaDynamicsInfo, dict):
        disorders["metaDynamicsInfo"] = "metaDynamicsInfo must be a dictionary"
        return disorders, False
    
    ## check height parameter
    disorders["metaDynamicsInfo"] = {}
    height = metaDynamicsInfo.get("height", None)
    if height is  None:
        disorders["metaDynamicsInfo"]["height"] = "No height specified in metaDynamicsInfo"
        metaOptionsOk = False
    else:
        if not isinstance(height, (int, float)):
            disorders["metaDynamicsInfo"]["height"] = "height must be a number"
            metaOptionsOk = False
        elif height <= 0:
            disorders["metaDynamicsInfo"]["height"] = "height must be positive"
            metaOptionsOk = False
        else:
            disorders["metaDynamicsInfo"]["height"] = None

    
    ## check biasFactor parameter
    biasFactor = metaDynamicsInfo.get("biasFactor", None)
    if biasFactor is None:
        disorders["metaDynamicsInfo"]["biasFactor"] = "No biasFactor specified in metaDynamicsInfo"
        metaOptionsOk = False
    else:
        if not isinstance(biasFactor, (int, float)):
            disorders["metaDynamicsInfo"]["biasFactor"] = "biasFactor must be a number"
            metaOptionsOk = False
        elif biasFactor <= 0:
            disorders["metaDynamicsInfo"]["biasFactor"] = "biasFactor must be positive"
            metaOptionsOk = False
        else:
            disorders["metaDynamicsInfo"]["biasFactor"] = None
    
    ## check biases parameter
    biases = metaDynamicsInfo.get("biases", None)
    ## make sure biases is present in metaDynamicsInfo
    if biases is None:
        disorders["metaDynamicsInfo"]["biases"] = "No biases specified in metaDynamicsInfo"
        metaOptionsOk = False
    else:
        ## make sure biases is a list
        if not isinstance(biases, list):
            disorders["metaDynamicsInfo"]["biases"] = "biases must be a list of biases (check README for more details)"
            metaOptionsOk = False
        ## make sure biases is not empty
        elif len(biases) == 0:
            disorders["metaDynamicsInfo"]["biases"] = "biases must contain at least one bias variable"
            metaOptionsOk = False
        else:
            ## check through each bias
            disorders["metaDynamicsInfo"]["biases"] = {}
            for biasCount, bias in enumerate(biases):
                disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"] = []
                ## make sure bias is a dictionary
                if not isinstance(bias, dict):
                    disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].append("biases must be a list of biases (check README for more details)")
                    metaOptionsOk = False
                ## check for biasVar entry in bias
                biasVar = bias.get("biasVar", None)
                if biasVar is None:
                    disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].append("No biasVar specified in bias")
                    metaOptionsOk = False
                ## make sure biasVar is a string with an allowed value
                else:
                    if not isinstance(biasVar, str):
                        disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].append("biasVar must be 'rmsd', 'torsion', 'distance', or 'angle'")
                        metaOptionsOk = False
                    elif not biasVar.lower() in ["rmsd", "torsion", "distance", "angle"]:
                        disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].append("biasVar must be 'rmsd', 'torsion', 'distance', or 'angle'")
                        metaOptionsOk = False
                ## check for minValue in bias
                minValue = bias.get("minValue", None)
                if minValue is None:
                    disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].append("No minValue specified in bias")
                    metaOptionsOk = False
                else:
                    if not isinstance(minValue, (int, float)):
                        disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].append("minValue must be a number")
                        metaOptionsOk = False
                ## check for maxValue in bias
                maxValue = bias.get("maxValue", None)
                if maxValue is None:
                    disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].append("No maxValue specified in bias")
                    metaOptionsOk = False
                else:
                    if not isinstance(maxValue, (int, float)):
                        disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].append("maxValue must be a number")
                        metaOptionsOk = False
                ## check for bias in bias
                biasWidth = bias.get("biasWidth", None)
                if biasWidth is None:
                    disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].append("No biasWidth specified in bias")
                    metaOptionsOk = False
                else:
                    if not isinstance(biasWidth, (int, float)):
                        disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].append("biasWidth must be a number")
                        metaOptionsOk = False
                ## check fir biasSelection in bias
                biasSelection = bias.get("selection", None)
                if biasSelection is None:
                    disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].append("No selection specified in bias")
                    metaOptionsOk = False
                else:
                    ## make sure selection is a dictionary
                    if not isinstance(biasSelection, dict):
                        disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].append("selection must be a dictionary")
                        metaOptionsOk = False
                    else:
                        ## check selection is correctly formatted
                        selectionDisorders = check_selection({"selection": biasSelection})
                        if len(selectionDisorders) > 0:
                            disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"].extend(selectionDisorders)
                            metaOptionsOk = False
                if len(disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"]) == 0:
                    disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"] = None

    return disorders, metaOptionsOk
#########################################################################
def check_nvt_npt_options(simulation: dict, stepName: str, disorders: dict) -> Tuple[dict,bool,dict]:
    ## check for required args for a nvt or npt simulation
    mdOptionsOk = True
    ## check duration
    duration = simulation.get("duration", None)
    if duration is None:
        disorders["duration"] = "Duration must be specified. example `duration = '1 ns'`"
        mdOptionsOk = False
    else:
        timeCheckProblems = check_time_input(duration, "duration", stepName)
        if timeCheckProblems is not None:
            disorders["duration"] = timeCheckProblems
            mdOptionsOk = False
        else:
            disorders["duration"] = None

    # check logInterval
    logInterval = simulation.get("logInterval", None)
    if logInterval is None:
        simulation["logInterval"] = "10 ps"
        disorders["logInterval"] = "No logInterval specified in simulation, using default of 10 ps"
        
    else:
        timeCheckProblems = check_time_input(logInterval, "logInterval", stepName)
        if timeCheckProblems is not None:
            disorders["logInterval"] = timeCheckProblems
            mdOptionsOk = False
        else:
            disorders["logInterval"] = None

    ## check heavyProtons
    heavyProtons = simulation.get("heavyProtons", None)
    if heavyProtons is None:
        disorders["heavyProtons"] = "No heavyProtons specified in simulation, using default of False"
        simulation["heavyProtons"] = False
    else:
        if not isinstance(heavyProtons, bool):
            disorders["heavyProtons"] = "heavyProtons must be a boolean"
            mdOptionsOk = False
        else:
            disorders["heavyProtons"] = None


    ## check timestep
    timestep = simulation.get("timestep", None)
    if timestep is None:
        if simulation["heavyProtons"] is True:
            simulation["timestep"] = "4 fs"
            disorders["timestep"] = "No timestep specified in simulation, heavyProtons set to True so using default of 4 fs"
        else:
            simulation["timestep"] = "2 fs"
            disorders["timestep"] = "No timestep specified in simulation, using default of 2 fs"
    else:
        timeCheckProblems = check_time_input(timestep, "timestep", stepName)
        if timeCheckProblems is not None:
            disorders["timestep"] = timeCheckProblems
            mdOptionsOk = False
        else:
            disorders["timestep"] = None
        
    return disorders, mdOptionsOk, simulation

        

#########################################################################
def check_time_input(timeInputValue: str, timeInputName: str, stepName: str) -> None:
    problemText = f"{timeInputName} in {stepName} must be in the format \"10 ns\""

    if not isinstance(timeInputValue, str):
        return problemText
    timeInputData = timeInputValue.split()
    if len(timeInputData) != 2:
        return problemText
    ## check if time is a number
    try:
        _ = float(timeInputData[0])
    except ValueError:
        return problemText
    
    if not timeInputData[1] in ["fs","ps","ns","ms"]:
        return problemText
    
    return None

#########################################################################
def check_shared_simulation_options(simulation: dict, disorders: dict) -> Tuple[dict, str, str, bool, dict]:

    sharedOptionsOk = True
    ## check simulation step name
    stepName = simulation.get("stepName", None)
    if stepName is None: 
        disorders["stepName"] = "No stepName specified in simulation"
        sharedOptionsOk = False
    elif not isinstance(stepName, str):
        disorders["stepName"] = "stepName must be a string"
        sharedOptionsOk = False
    elif " " in stepName:
        disorders["stepName"] = "No whitespace allowed in stepName"
        sharedOptionsOk = False
    else:
        disorders["stepName"] = None


    ## check simulationType
    simulationType = simulation.get("simulationType", None)
    if simulationType is None:
        simulation["simulationType"] = "NPT"
        disorders["simulationType"] = "No simulationType specified in simulation, using NPT as default"
    elif not simulationType.upper() in ["EM", "NVT", "NPT", "META"]:
        disorders["simulationType"] = "simulationType in simulation must be one of the following: 'EM', 'NVT', 'NPT', 'META'"
        sharedOptionsOk = False
    else:
        disorders["simulationType"] = None
    ## check for either temperature or temparatureRange in simulation
    temperature = simulation.get("temperature", None)
    tempRange = simulation.get("temperatureRange", None)
    ## check for both temperature and tempRange (this is not allowed!)
    if temperature and tempRange:
        disorders["temperature"] = "Cannot specify both temperature and temperatureRange in simulation"
        sharedOptionsOk = False
    ## check for neither temperature or tempRange
    ## Use a default value of 300 K
    elif not temperature and not tempRange:
        simulation["temperature"] = 300
        disorders["temperature"] = "No temperature or temperatureRange specified in simulation, using 300 K as default"
    ## check temperature
    elif temperature:
        if not isinstance(temperature, int):
            disorders["temperature"] = "Temperature in simulation must be an int"
            sharedOptionsOk = False
        if temperature < 0:
            disorders["temperature"] = "Temperature in simulation must be a positive int"
            sharedOptionsOk = False
        else:
            disorders["temperature"] = None
    ## check tempRange
    elif tempRange:
        if not isinstance(tempRange, list):
            disorders["tempRange"] = "TemperatureRange in simulation must be a list of ints"
            sharedOptionsOk = False
        if len(tempRange) == 0:
            disorders["tempRange"] = "TemperatureRange in simulation must be a list of at least one int"
            sharedOptionsOk = False
        for temp in tempRange:
            if not isinstance(temp, int):
                disorders["tempRange"] = "TemperatureRange in simulation must be a list of ints"
                sharedOptionsOk = False
            if temp < 0:
                disorders["tempRange"] = "TemperatureRange in simulation must be a list of positive ints"
                sharedOptionsOk = False

    return disorders, stepName, simulationType, sharedOptionsOk, simulation

#########################################################################
def validate_path(argName: str, argPath: Union[FilePath, DirectoryPath]) -> str:
    """
    Check to see if a path variable is indeed the correct type
    Check to see if the path exists
    """
    if  not isinstance(Path(argPath), (PathLike, str)) :
        return f"The config argument {argName} = {argPath} is not a PathLike."
    # Check if the path exists
    if not p.exists(argPath):
        return f"The config argument {argName} = {argPath} does not exist."
    return None
#####################################################################################
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
#####################################################################################

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


#####################################################################################
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

#####################################################################################
