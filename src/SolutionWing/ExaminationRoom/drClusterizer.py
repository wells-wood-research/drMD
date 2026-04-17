import mdtraj as md
## BASIC PYTHON LIBRARIES
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import numpy as np
import os
from os import path as p

## drMD LIBRARIES
from ExaminationRoom import drLogger
from UtilitiesCloset import drSelector, drListInitiator

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

## CLEAN CODE
from typing import Dict, Union, Any, List
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath
#######################################################################
def clustering_manager(pathInfo: Dict, clusterInfo: Dict) -> List[FilePath]: 
    """
    Identifies the directories to cluster and performs clustering on them.

    Args:
        pathInfo (Dict): The path information dictionary.
        clusterInfo (Dict): The clustering information dictionary.

    Returns:
        clusterPdbs (List[FilePath]): The list of cluster PDB files.
    """
    ## unpack pathInfo to get outDir
    outDir: DirectoryPath = pathInfo["outputDir"]

    ## define and creare a cluster directory
    clusterDir: DirectoryPath = p.join(outDir,"00_clustered_pdbs")
    if p.exists(clusterDir):
        drLogger.log_info(f"Cluster directory already exists: {clusterDir}. Skipping...", True)
        return None
    os.makedirs(clusterDir, exist_ok=True)

    ## list of dirs created by drMD that we don't want to cluster
    notRunDirs = drListInitiator.get_not_a_run_dir()
    ## create list of dirs to cluster
    runDirs: List[DirectoryPath] = [p.join(outDir, dir) for dir in os.listdir(outDir) if not dir in notRunDirs]
    dirsToCluster: List[DirectoryPath] = [p.join(runDir,stepDir) for stepDir in clusterInfo["stepNames"] for runDir in runDirs]

    ## init a list of to store pdb files that have been created by clustering
    allClusterPdbs: list[FilePath] = []
    ## Run clustering on each specified directory
    for dirToCluster in dirsToCluster:
        drLogger.log_info(f"Clustering trajectory for system:{p.basename(p.dirname(dirToCluster))} and step:  {p.basename(dirToCluster)}", True)
        clusterPdbs: List[FilePath] = rmsd_clustering_protocol(dirToCluster, clusterInfo, clusterDir)
        ## add cluster pdbs to list 
        allClusterPdbs.extend(clusterPdbs)
    
    return allClusterPdbs

#######################################################################
def rmsd_clustering_protocol(inDir: DirectoryPath, clusterInfo: Dict[str, Union[int, str]], clusterDir: DirectoryPath) -> List[FilePath]:
    """
    Clusters a trajectory based on the RMSD values of a subset of atoms.

    Args:
        inDir (str): The input directory containing the trajectory files.
        clusterInfo (Dict[str, Union[int, str]]): A dictionary containing the number of clusters (int)
            and the selection of atoms to be used for clustustering (str).

    Returns:
        clusterPdbs (List[FilePath]): The list of cluster PDB files.
    """


    # Get protName to name output files
    stepName: str = p.basename(inDir)
    protName: str = p.basename(p.dirname(inDir))

    drLogger.log_info(f"Clustering trajectory for system: {protName} and step: {stepName}", True)

    thisClusterDir: DirectoryPath = p.join(clusterDir, protName, stepName)
    os.makedirs(thisClusterDir, exist_ok=True)

    ## unpack clusterInfos
    nClusters: int = clusterInfo["nClusters"]
    clusterBy: str = clusterInfo["clusterBy"]

    ## find trajectory file and matching pdb file
    dcdFile: FilePath = p.join(inDir, "trajectory.dcd")
    pdbFile: FilePath = p.join(inDir, "trajectory.pdb")


    # Check if the trajectory.dcd and output.pdb files exist, return an empty list if not
    if not p.isfile(dcdFile) or not p.isfile(pdbFile):
        drLogger.log_info(f"  * No PDB file or DCD file found in {inDir} *", True, True)
        return []

    # Get the atom indexes for the selected atoms
    clusterSelectionAtomIndexes: List[int] = []
    for clusterBySelection in clusterBy:
        clusterSelection = clusterBySelection["selection"]
        clusterSelectionAtomIndexes.extend(drSelector.get_atom_indexes(clusterSelection, pdbFile))

    # Load trajectory
    traj: md.Trajectory = md.load(dcdFile, top=pdbFile)
    # Optionally, superimpose all frames to the first to remove translational and rotational motions
    traj.superpose(traj[0])

    # Convert trajectory to RMSD matrix
    rmsdMatrix: np.ndarray = convert_traj_to_rmsdMatrix(traj,clusterSelectionAtomIndexes)

    # Perform clustering and save the clusters to PDB files
    clusterPdbs: List[FilePath] = kmeans_clusters_to_pdb(rmsdMatrix, nClusters, thisClusterDir, traj, protName)
    
    return clusterPdbs
#######################################################################
def convert_traj_to_rmsdMatrix(traj: md.Trajectory, atomIndexes: List[int]) -> np.ndarray:
    """
    Converts a trajectory to an RMSD matrix.

    Args:
        traj (md.Trajectory): The input trajectory.
        atomIndexes (List[int]): The atom indexes to be included in the subtrajectory.

    Returns:
        rmsdMatrix (np.ndarray): The RMSD matrix.
    """
    # Create a subtrajectory only containing the selected atoms
    sectionTraj: md.Trajectory = traj.atom_slice(atomIndexes)

    # Compute the pairwise RMSD matrix for all frames
    nFrames: int = sectionTraj.n_frames  # Get the number of frames
    rmsdMatrix = np.empty((nFrames, nFrames))  # Create an empty RMSD matrix
    for i in range(nFrames):
        rmsdMatrix[i] = md.rmsd(sectionTraj, sectionTraj, frame=i)  # Compute the RMSD for each pair of frames

    return rmsdMatrix


#######################################################################
def kmeans_clusters_to_pdb(rmsdMatrix: np.ndarray,
                            bestK: int,
                              outDir: DirectoryPath,
                                traj: md.Trajectory,
                                protName: str) -> List[FilePath]:
    """
    Saves the cluster centers of a k-means clustering to PDB files.

    Args:
        rmsdMatrix (np.ndarray): The RMSD matrix.
        bestK (int): The number of clusters.
        outDir (str): The output directory.
        traj (md.Trajectory): The input trajectory.

    Returns:
        clusterPdbs (List[FilePath]): The list of cluster PDB files.
    """
    ## If the user has specified bestK to be -1
    ## We will use silhouette scores to find the best number of clusters
    ## NOTE that this takes a while and produces less clusters than you may want!
    if bestK == -1:
        bestK = find_best_k_with_silhouette(rmsdMatrix)

    ## create a KMeans object and fit it to the RMSD matrix
    kmeansModel = KMeans(n_clusters=bestK)  # Create a KMeans object
    _ = kmeansModel.fit_predict(rmsdMatrix)  # Fit the KMeans object to the RMSD matrix

    ## create an empty list to store pdb file locations
    clusterPdbs: List[FilePath] = []
    for i in range(bestK):
        # Find the frame closest to each cluster center
        clusterCentroidIndex: int = np.argmin(np.linalg.norm(rmsdMatrix - kmeansModel.cluster_centers_[i], axis=1))
        # Save the representative frame as a PDB file
        clusterPdb: Union[PathLike, str] = p.join(outDir, f"{protName}_cluster_{i+1}.pdb")
        traj[clusterCentroidIndex].save_pdb(clusterPdb)
        ## add the PDB file location to the list
        clusterPdbs.append(clusterPdb)

    ## fix element names (wont work for 2-letter elements)
    for clusterPdb in clusterPdbs:
        clusterDf = pdbUtils.pdb2df(clusterPdb)
        clusterDf["ELEMENT"] = clusterDf["ATOM_NAME"].apply(lambda x: x[0])
        pdbUtils.df2pdb(clusterDf, clusterPdb)

    return clusterPdbs

#######################################################################

def find_best_k_with_silhouette(rmsdMatrix: np.ndarray) -> int:
    """
    Finds the best number of clusters using the silhouette score.

    Args:
        rmsdMatrix (np.ndarray): The RMSD matrix.

    Returns:
        bestK (int): The best number of clusters.
    """


    # Perform silhouette score analysis to find the best number of clusters
    silhouetteScores: List[float] = []
    nClusterRange: range = range(2, 25)  
    for nClusters in nClusterRange:
        kmeans = KMeans(n_clusters=nClusters)
        cluster_labels: np.ndarray = kmeans.fit_predict(rmsdMatrix)
        silhouette_avg: float = silhouette_score(rmsdMatrix, cluster_labels)
        silhouetteScores.append(silhouette_avg)
    bestK: int = nClusterRange[np.argmax(silhouetteScores)]

    return bestK    

#######################################################################