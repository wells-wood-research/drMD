## BASIC PYTHON LIBRARIES
import os
from os import path as p
import yaml
from subprocess import run



## ERROR HANDLING ##
import traceback
import inspect



    

def standardOperation(OperationName: str) -> dict:
    """
    adds an operation config file parameters to the batch config in drReferral if the operation name is the same as a standard operation.
    If a user wants to add custom standard operations, add the following code and have the yaml file in the dir

    elif OperationName == "Operation name":
        with open(f"Operation name config file", "r") as yamlFile:
            Procedure: dict = yaml.safe_load(yamlFile)
            os.chdir(os.path.expanduser('~'))
            os.chdir(CurrentWorkingDir)
        return Procedure
    
    Args:
        OperationName: file path to pdbfile 
    """

    ## sets up file paths
    CurrentWorkingDir= os.getcwd()
    FilePath= os.path.dirname(os.path.abspath(__file__))
    os.chdir(os.path.expanduser('~'))
    os.chdir(FilePath)

    ##Checks operation name for a matching operation yaml file
    if OperationName == "Solvent_Equilibriation":
        with open(f"Solvent_Equilibriation.yaml", "r") as yamlFile:
            Procedure: dict = yaml.safe_load(yamlFile)
            os.chdir(os.path.expanduser('~'))
            os.chdir(CurrentWorkingDir)
        return Procedure
    elif OperationName == "Vacuum_Equilibriation":
        with open(f"Vacuum_Equilibriation.yaml", "r") as yamlFile:
            Procedure: dict = yaml.safe_load(yamlFile)
            os.chdir(os.path.expanduser('~'))
            os.chdir(CurrentWorkingDir)
            return Procedure
    elif OperationName == "Vacuum_Metadynamics":
        with open(f"Vacuum_Metadynamics.yaml", "r") as yamlFile:
                Procedure: dict = yaml.safe_load(yamlFile)
        os.chdir(os.path.expanduser('~'))
        os.chdir(CurrentWorkingDir)
        return Procedure
    elif OperationName == "Solvent_Metadynamics":
        with open(f"Solvent_Metadynamics.yaml", "r") as yamlFile:
                Procedure: dict = yaml.safe_load(yamlFile)
        os.chdir(os.path.expanduser('~'))
        os.chdir(CurrentWorkingDir)
        return Procedure
    elif OperationName.find("ESI") != -1:
        os.chdir(os.path.expanduser('~'))
        os.chdir(CurrentWorkingDir)
        return 
    else:
        print("Error: Standard Procedure not found")
        os.chdir(os.path.expanduser('~'))
        os.chdir(CurrentWorkingDir)
        exit


def getSimulationType(Name: str) -> str:
    """
    adds simulation Type based on simulation name

    
    Args:
        Name: Name of the standard operation
    """
    if Name.find("Vacuum") != -1:
        Type = "Vacuum"
    elif Name.find("Solvent") != -1:
        Type = "Solution"
    elif Name.find("ESI") != -1:
         Type = "ESI"
    return Type