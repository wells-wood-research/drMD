## BASIC PYTHON LIBRARIES
import os
from os import path as p
import numpy as np
import yaml
import itertools
import argpass
import re
import subprocess
from subprocess import run
from ESIWing.SpecialistEquipment import drDroplet



## ERROR HANDLING ##
import traceback
import inspect

## PARALLELISATION LIBRARIES
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from concurrent.futures.process import BrokenProcessPool


def ESIOperation(ESI_Count, Water_Count) -> dict:
    """
    adds simulation Type and config if the simulation is ESI

    
    Args:
        None
    """
    CurrentWorkingDir= os.getcwd()
    FilePath= os.path.dirname(os.path.abspath(__file__))
    os.chdir(os.path.expanduser('~'))
    os.chdir(FilePath)
    if ESI_Count <= 1 :
        with open(f"ESI_Inital.yaml", "r") as yamlFile:
            Procedure: dict = yaml.safe_load(yamlFile)
            os.chdir(os.path.expanduser('~'))
            os.chdir(CurrentWorkingDir)
    else:
        if Water_Count >= 50:
            with open(f"ESI.yaml", "r") as yamlFile:
                Procedure: dict = yaml.safe_load(yamlFile)
                os.chdir(os.path.expanduser('~'))
                os.chdir(CurrentWorkingDir)
        else:
            with open(f"ESI_Final.yaml", "r") as yamlFile:
                Procedure: dict = yaml.safe_load(yamlFile)
                os.chdir(os.path.expanduser('~'))
                os.chdir(CurrentWorkingDir)
    return Procedure



def ESI_Handler(Cutoff: int, PDBfile: str, ESI_Count, Mode):
    '''
    Removes all atoms from a pdb file further than the cutoff distance, then cleans the pdb file for further simulations using cpptraj

    Args:
        Cutoff(int): Cutoff distance for ESI simulations
        PDBfile(str): pdbfile path
        ESI_Count(int): increment of ESI simulation
    Returns:
        Nothing
    '''
    



    ## handles the save name of the output pdb file 
    saveName= PDBfile.removesuffix(".pdb")+"_strip.pdb"
    

    if ESI_Count >= 2:
        Strip(PDBfile, saveName, Cutoff)
    else:
        Prep(PDBfile, saveName)
    try:
      os.remove(f"{PDBfile}")
    except:
      print(f"no file called {PDBfile} could be found")
    if ESI_Count <=1:
        drDroplet.DropletFormation(saveName, Mode)
    os.rename(saveName, PDBfile)
    return 

def Strip(PDBfile, saveName, Cutoff):
        ## finds the number of residues in the protein
    with open(PDBfile, 'r') as pdb:
      Residue_Number=0
      lines= pdb.readlines()
      for rows in lines:
        if rows.find("OXT") != -1:
             Values=rows.split()
             Residue_Number= Values[5]
    
    ##writes the input file to be run in cpptraj 
    cpptrajInputs: str = p.join( "Strip.in")
    with open(cpptrajInputs, "w") as f:
        f.write(f"parm {PDBfile}\n")
        f.write(f"trajin {PDBfile}\n")
        f.write(f"reference  {PDBfile}\n")
        f.write(f"center :1-{Residue_Number} mass\n")
        f.write(f"strip !(:1-{Residue_Number}<:{Cutoff})\n")
        f.write(f"fixatomorder\n")
        f.write(f'trajout {saveName} nobox\n')

    CpptrajOutput: str= p.join("Cpptraj.out")

    ## command to run cpptraj
    CpptrajCommand: str = f"cpptraj -i Strip.in"
    try: 
       # Execute the command and capture its output
        result: subprocess.CompletedProcess[str] = subprocess.run(
            CpptrajCommand.split(),
            capture_output=True,
            check=True,
            text=True, 
            env = os.environ
            )
    except Exception as errorMessage:
            print(errorMessage)
    os.remove(f"{PDBfile}")
    return

def Prep(PDBfile, saveName):
            ## finds the number of residues in the protein
    with open(PDBfile, 'r') as pdb:
      Residue_Number=0
      lines= pdb.readlines()
      for rows in lines:
        if rows.find("OXT") != -1:
             Values=rows.split()
             Residue_Number= Values[5]
    
    ##writes the input file to be run in cpptraj 
    cpptrajInputs: str = p.join( "Strip.in")
    with open(cpptrajInputs, "w") as f:
        f.write(f"parm {PDBfile}\n")
        f.write(f"trajin {PDBfile}\n")
        f.write(f"reference  {PDBfile}\n")
        f.write(f"center :1-{Residue_Number} mass\n")
        f.write(f"strip !(:1-{Residue_Number})\n")
        f.write(f"fixatomorder\n")
        f.write(f'trajout {saveName} nobox\n')

    CpptrajOutput: str= p.join("Cpptraj.out")

    ## command to run cpptraj
    CpptrajCommand: str = f"cpptraj -i Strip.in"
    try: 
       # Execute the command and capture its output
        result: subprocess.CompletedProcess[str] = subprocess.run(
            CpptrajCommand.split(),
            capture_output=True,
            check=True,
            text=True, 
            env = os.environ
            )
    except Exception as errorMessage:
            print(errorMessage)
    
    return