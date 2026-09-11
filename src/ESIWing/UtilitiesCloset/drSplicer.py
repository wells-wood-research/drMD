## BASIC PYTHON LIBRARIES
import os
from os import path as p
from typing import Union, Dict
import pandas as pd

## OPENMM LIBRARIES
from openmm import unit

## MDTRAJ LIBRARIES
import mdtraj as md

## drMD LIBRARIES
from ExaminationRoom import drLogger
from UtilitiesCloset import drSelector

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

##  CLEAN CODE
from typing import List, Dict
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath


#####################################################################################
def merge_partial_outputs(simDir: DirectoryPath, pdbFile: FilePath, simInfo: Dict, config: Dict) -> None:
    """
    After firstAid protocols have been run, we need to merge the partial reports and trajectories

    Args:
        simDir (Union[PathLike, str]): The path to the simulation directory
        pdbFile (Union[PathLike, str]): The path to the pdb file

    Returns:
        None
    """
    drLogger.log_info(f"-->{' '*4}Merging partial outputs...", True)
    # merge vitals reports
    vitalsDf = merge_partial_reports(simDir, "vitals_report", removePartials=True)
    vitalsDf = fix_merged_vitals(vitalsDf, simInfo)
    vitalsDf.to_csv(p.join(simDir, "vitals_report.csv"))
    ## merge progress reports
    merge_partial_reports(simDir, "progress_report", removePartials=True)
    ## merge trajectories
    merge_partial_trajectories(config=config,
                               simDir = simDir,
                               pdbFile = pdbFile,
                                   removePartials=True)
#####################################################################################
def make_trajectory_pdb(trajectorySelections: List[Dict], pdbFile: FilePath, outDir: DirectoryPath) -> None:
    """
    Creates a trajectory PDB file based on the provided trajectory selections and PDB file.

    Args:
        trajectorySelections (List[Dict]): A list of dictionaries containing the trajectory selections.
        pdbFile (FilePath): The path to the PDB file.
        outDir (DirectoryPath): The output directory for the trajectory PDB file.

    Returns:
        None
    """
    dcdAtomSelection: List = []
    for selection in trajectorySelections:
        dcdAtomSelection.extend(drSelector.get_atom_indexes(selection["selection"], pdbFile))


    pdbDf = pdbUtils.pdb2df(pdbFile)
    dcdDf = pdbDf.iloc[dcdAtomSelection]

    trajectoryPdb = p.join(outDir, "trajectory.pdb")
    pdbUtils.df2pdb(dcdDf, trajectoryPdb)

    return trajectoryPdb
#####################################################################################
def merge_partial_reports(simDir: DirectoryPath, matchString: str, removePartials: bool = False) -> None:
    """
    Merges partial reports into a single report

    Args:
        simDir (Union[PathLike, str]): The path to the simulation directory
        matchString (str): The string to match
        removePartials (bool, optional): Whether to remove partial reports. Defaults to False.  
    """
    ## rename the last report to be made
    lastReport = p.join(simDir, f"{matchString}.csv")
    if p.isfile(lastReport):
        os.rename(lastReport, p.join(simDir, f"{matchString}_partial_999.csv"))

    ## collect all partial reports into a list exept for last one written
    reports = []
    for file in os.listdir(simDir):
        if ("partial" in file) and (file.startswith(matchString)) and (p.splitext(file)[1] == ".csv") and (p.getsize(p.join(simDir, file)) > 0):
            reports.append(p.join(simDir, file))

    ## sort the list so that reports are in chronological order
    reports = sorted(reports)
    ## merge reports
    dfsToConcat = []
    for report in reports:
        dfsToConcat.append(pd.read_csv(report))

    # concatonate | write back to csv
    # remove partial reports to tidy up 
    [os.remove(report) for report in reports if removePartials]
    # concat dataframes and write to file
    df = pd.concat(dfsToConcat, ignore_index=True)

    df.to_csv(p.join(simDir, f"{matchString}.csv"))

    return df
#######################################################################
def merge_dcd_files(dcdFiles: list[FilePath],
                    pdbFile: FilePath,
                    outputDcd: FilePath) -> FilePath:
    
    traj = md.load_dcd(dcdFiles[0], top = pdbFile)
    for file in dcdFiles[1:]:
        newTraj = md.load_dcd(file, top = pdbFile)
        traj = traj + newTraj
    traj.save_dcd(outputDcd)
#######################################################################
def merge_partial_trajectories(config: Dict,
                               simDir: DirectoryPath,
                                pdbFile: FilePath,
                                  removePartials: bool = False) -> None:
    """
    Merges partial trajectories into a single trajectory

    Args:
        simDir (Union[PathLike, str]): The path to the simulation directory
        pdbFile (Union[PathLike, str]): The path to the pdb file
        removePartials (bool, optional): Whether to remove partial trajectories. Defaults to False.  
    """
    trajectorySelections = config["miscInfo"]["trajectorySelections"]
    trajectoryPdb = make_trajectory_pdb(trajectorySelections, pdbFile, simDir)
    ## rename last trajectory to be made
    lastTrajectory = p.join(simDir, "trajectory.dcd")
    if p.isfile(lastTrajectory):
        os.rename(lastTrajectory, p.join(simDir, "trajectory_partial_99.dcd"))

    ## collect all partial trajectories into a list exept for last one written
    ## collect all partial reports into a list exept for last one written
    trajectories = []
    for file in os.listdir(simDir):
        if ("partial" in file) and (file.startswith("trajectory")) and (p.splitext(file)[1] == ".dcd"):
            trajectories.append(p.join(simDir, file))

    ## sort the list so that trajectories are in chronological order
    trajectories = sorted(trajectories)

    ## merge trajectories
    merge_dcd_files(trajectories, trajectoryPdb,  lastTrajectory)

    ## delete partial trajectories
    if removePartials:
        for trajectory in trajectories:
            if p.basename(trajectory) == "trajectory.dcd":
                continue
            os.remove(trajectory)

#######################################################################
def fix_merged_vitals(vitalsDf: pd.DataFrame, simInfo: Dict) -> pd.DataFrame:
    ## read stuff from simInfo
    logInterval: int = simInfo["logInterval"]
    timeStep: openmm.Quantity = simInfo["timestep"]
    duration: openmm.Quantity = simInfo["duration"]
    ## convert to ints 
    logInterval_ps = int(logInterval * timeStep.value_in_unit(unit.picoseconds))
    duration_ps: int = round(duration.value_in_unit(unit.picoseconds))

    ## construct time range
    timeRange_ps = range(logInterval_ps, duration_ps + logInterval_ps, logInterval_ps)

    ## construct step range
    stepsRange = [val * logInterval for val in timeRange_ps]    

    vitalsDf['Time (ps)'] = timeRange_ps
    vitalsDf['#"Step"'] = stepsRange
    return vitalsDf
#######################################################################

if __name__ == "__main__":
    simDir = "/home/esp/scriptDevelopment/drMD/03_outputs/SC_bpy_with_Pd/05_production"
    pdbFile = "/home/esp/scriptDevelopment/drMD/03_outputs/SC_bpy_with_Pd/00_prep/WHOLE/SC_bpy_with_Pd_solvated.pdb"

    simInfo = {
        "logInterval": 1250,
        "timestep": 4 * unit.femtoseconds,
        "duration": 500 * unit.picoseconds
    }
    config = {
        "loggingInfo": {
            "trajectorySelections": [{"selection": {"keyword" : "protein"}}, 
                                     {"selection": {"keyword" : "custom",
                                                  "customSelection": [{"CHAIN_ID": "A", "RES_NAME": "C8X", "RES_ID": "all", "ATOM_NAME": "all"},
                                    {"CHAIN_ID": "A", "RES_NAME": "PD", "RES_ID": "all", "ATOM_NAME": "all"}]  }}]
             
             
            
        }
    }

    merge_partial_outputs(simDir, pdbFile, simInfo, config)