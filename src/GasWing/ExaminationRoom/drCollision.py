#!/usr/bin/env python3
"""Compute a protein convex hull and solvent-accessible surface area (SASA) from a PDB file."""

import argparse
import warnings
import subprocess
from pathlib import Path
import shutil
import tempfile
import math
import os
import yaml
from os import path as p
import numpy as np
from pdbUtils import pdbUtils
import pandas as pd
import time
import mdtraj as md
import MDAnalysis as mda
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath
## plotting libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from ExaminationRoom.drCollisionsReport import create_ccs_report_html
from PIL import Image, ImageFilter, ImageEnhance, ImageDraw, ImageFont, ImageSequence
from PIL.Image import Resampling

## BASIC PYTHON LIBRARIES
import os
from os import path as p 
from shutil import copy, rmtree
import MDAnalysis as mda
import mdtraj as md
import time

## drMD LIBRARIES
from ExaminationRoom import drLogger, drClusterizer
from UtilitiesCloset import drSelector, drListInitiator, drMethodsWriter

## PARALLELISATION LIBRARIES
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from concurrent.futures.process import BrokenProcessPool


## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils



## CLEAN CODE
from typing import List, Dict, Union, Any, Optional
from os import PathLike
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

## get ShadowScreen Executable
THISDIR = p.dirname(__file__)
SHADOWSCREEN_EXE = p.join(THISDIR, "shadowScreen")








def visualise_shadow(pdbFile, probe: str) -> float:
    import glob
    from PIL import Image, ImageOps
    outDir= p.dirname(pdbFile)
    pdbName= p.basename(pdbFile).removesuffix(".pdb")
    output_image_path = p.join(outDir, f"{pdbName}_shadow.gif")
    shadowScreenCommand: list = [SHADOWSCREEN_EXE+"AAImage", pdbFile, "250", "0.3", "--gas", probe, "--output-dir", outDir, "--debug", "--smooth"]
    shadowDir= outDir+"/debug_shadows"
    try: 
        # Execute the command and capture its output
        result: subprocess.CompletedProcess[str] = subprocess.run(
            shadowScreenCommand,
            capture_output=True,
            check=True,
            text=True, 
            env = os.environ
        )
    except Exception as errorMessage:
        print(errorMessage)

    ppmFiles = [p.join(shadowDir, frameFile) for frameFile in os.listdir(shadowDir) if p.splitext(frameFile)[1] == ".ppm"]
    ppmFiles.sort(key=lambda x: int(x.split("_")[-1].split(".")[0]))
    
    invertedFrames = []
    for image_path in ppmFiles:
        # Open inside a context manager to ensure the file handle closes properly
        with Image.open(image_path) as img:
            img.load() # <-- CRITICAL: Forces Pillow to read the raw pixel data into RAM instantly
            rgbImg = img.convert("RGB")
            invertedRgb = ImageOps.invert(rgbImg)
            
            finalFrame = invertedRgb.convert("P", palette=Image.ADAPTIVE)
            invertedFrames.append(finalFrame) # Safely stored entirely in memory

    if invertedFrames:
        invertedFrames[0].save(
            output_image_path,
            format="GIF",
            save_all=True,
            append_images=invertedFrames[1:],
            duration=50,
            loop=0
        )

    normalFrames= []
    clusterNumber= pdbName.removesuffix(".pdb").split("_")[-1]
    output_image_path = p.join(outDir, f"cluster_{clusterNumber}.gif")
    for image_path in ppmFiles:
        # Open inside a context manager to ensure the file handle closes properly
        with Image.open(image_path) as img:
            img.load() # <-- CRITICAL: Forces Pillow to read the raw pixel data into RAM instantly
            normalFrames.append(img) # Safely stored entirely in memory
    if normalFrames:
        normalFrames[0].save(
            output_image_path,
            format="GIF",
            save_all=True,
            append_images=normalFrames[1:],
            duration=50,
            loop=0
        )


    # Now it is 100% safe to clear the directory, as Python no longer relies on the disk files
    shutil.rmtree(shadowDir, ignore_errors=True)
    pdbGif= pdbName.removesuffix(".pdb")+".gif"
    shutil.rmtree(outDir + "/debug_gifs", ignore_errors=True)
    return output_image_path
def configure_scripts() -> None:
    permissionCommand: list = ["chmod", "-R", "777", THISDIR]
    try: 
        # Execute the command and capture its output
        result: subprocess.CompletedProcess[str] = subprocess.run(
                permissionCommand,
                capture_output=True,
                check=True,
                text=True, 
                env = os.environ
            )
    except Exception as errorMessage:
            print(errorMessage)
def main(pdbFile, representation, probe, rotations, resolution) -> None:

    totalTime= time.time()
    configure_scripts()
    ccsTime= time.time()
    shadowScreenCommand: list = [f"{SHADOWSCREEN_EXE}{representation}", pdbFile, str(rotations), str(resolution), "--gas", probe, "--ccs"]
    try: 
        # Execute the command and capture its output
        result: subprocess.CompletedProcess[str] = subprocess.run(
                shadowScreenCommand,
                capture_output=True,
                check=True,
                text=True, 
                env = os.environ
            )
    except Exception as errorMessage:
            print(errorMessage)

    ccsRaw= float(result.stdout)
    if probe == "N2":
        ccsCorrected= 0.843*ccsRaw**1.0503
    elif probe == "He":
        ccsCorrected= 0.841*ccsRaw**1.0565


    get_frame_shadow(pdbFile, representation, probe)
    return ccsCorrected


def read_image(filename):
    with open(filename, "rb") as f:

        def next_token():
            while True:
                token = b""
                c = f.read(1)
                while c and c.isspace():
                    c = f.read(1)
                if not c:
                    return None
                if c == b"#":
                    f.readline()
                    continue
                while c and not c.isspace():
                    token += c
                    c = f.read(1)
                return token

        magic = next_token().decode()

        width = int(next_token())
        height = int(next_token())
        maxval = int(next_token())

        if maxval > 255:
            raise ValueError("Only 8-bit images supported")

        if magic == "P2":
            pixels = bytes(int(next_token()) for _ in range(width * height))

        elif magic == "P5":
            pixels = f.read(width * height)

        elif magic == "P3":
            rgb = bytes(int(next_token()) for _ in range(width * height * 3))
            pixels = bytearray()
            for i in range(0, len(rgb), 3):
                r, g, b = rgb[i:i+3]
                pixels.append((299*r + 587*g + 114*b) // 1000)

        elif magic == "P6":
            rgb = f.read(width * height * 3)
            pixels = bytearray()
            for i in range(0, len(rgb), 3):
                r, g, b = rgb[i:i+3]
                pixels.append((299*r + 587*g + 114*b) // 1000)

        else:
            raise ValueError(f"Unsupported image format: {magic}")

        return width, height, pixels

def ppm_to_bitmap(width, height, pixels, threshold=128):
    """
    Convert a grayscale (P2/P5) image into a packed 1-bit bitmap.

    Parameters
    ----------
    width : int
        Image width.
    height : int
        Image height.
    pixels : iterable of int
        One grayscale value (0-255) per pixel.
    threshold : int
        Pixels >= threshold become 1, otherwise 0.

    Returns
    -------
    bytes
        Packed bitmap (8 pixels per byte).
    """

    threshold = int(threshold)

    out = bytearray()
    current_byte = 0
    bit_count = 0

    for value in pixels:
        # Ensure we have an integer
        value = int(value)

        # Threshold to binary
        bit = 1 if value >= threshold else 0

        current_byte = (current_byte << 1) | bit
        bit_count += 1

        if bit_count == 8:
            out.append(current_byte)
            current_byte = 0
            bit_count = 0

    # Pad the final byte if necessary
    if bit_count:
        current_byte <<= (8 - bit_count)
        out.append(current_byte)

    return bytes(out)

def get_frame_shadow(pdbFile, representation, probe):
    outDir= p.dirname(pdbFile)
    pdbBaseName = p.basename(pdbFile).removesuffix(".pdb")
    # Use a per-call temp directory so parallel workers do not race on debug_shadows.
    shadowOutDir = tempfile.mkdtemp(prefix=f"shadow_{pdbBaseName}_", dir=outDir)
    shadowScreenCommand: list = [SHADOWSCREEN_EXE+"Frame", pdbFile, "100", "0.25", "--gas", probe, "--debug", "--output-dir", shadowOutDir]
    try: 
        # Execute the command and capture its output
        result: subprocess.CompletedProcess[str] = subprocess.run(
                shadowScreenCommand,
                capture_output=True,
                check=True,
                text=True, 
                env = os.environ
            )
    except Exception as errorMessage:
            print(errorMessage)

    pdbName= pdbBaseName.removeprefix(outDir)
    framePpm= p.join(shadowOutDir, "debug_shadows", pdbName +"_max_ccs.ppm")
    
    width, height, pixels = read_image(framePpm)

    bitmap = ppm_to_bitmap(width, height, pixels, 128)

    frameDir= outDir+"/frames/frame_"+ pdbName +".bin"
    os.makedirs(outDir+"/frames", exist_ok=True)
    with open(frameDir, "wb") as f:
        f.write(width.to_bytes(4, "little"))
        f.write(height.to_bytes(4, "little"))
        f.write(bitmap)
    shutil.rmtree(shadowOutDir, ignore_errors=True)
if __name__ == "__main__":
    main()


######################################################################################################
"""
New code for CCS reporting and clustering.

Following code contains the following functions:
- plot_ccs_shadow_data: Plots the CCS data from a CSV file generated by the ShadowScreen
- cluster_ccs_shadow_data: Clusters the CCS data from a CSV file generated by the ShadowScreen and generates a KMeans distribution plot
- plot_kmeans_distribution: Plots the KMeans cluster distribution for the CCS array
- _process_ccs_frame_chunk: Processes a chunk of trajectory frames to compute CCS values
- calculate_ccs_for_protein: Calculates CCS values for all requested simulation steps for a single protein
- calculate_ccs_shadow: Uses the Shadow Screen to calculate CCS for all proteins in the batch configuration
"""
######################################################################################################
def plot_ccs_shadow_data(csvFile, gas: str, logTimeUnit: str) -> None:
    """
    Plots the CCS data from a CSV file generated by the ShadowScreen

    Args:
        csvFile (FilePath): The path to the CSV file containing CCS data
        gas (str): The gas used for the simulation
        logTimeUnit (str): The time unit for the log interval
    """


    # Read the CSV file into a DataFrame
    df = pd.read_csv(csvFile)

    # Plotting
    cmap = plt.get_cmap("viridis")
    x = pd.to_numeric(df.iloc[:, 0], errors="coerce").astype(float)
    y = pd.to_numeric(df.iloc[:, 1], errors="coerce").astype(float)
    shapeFactor = (
        pd.to_numeric(df.iloc[:, 3], errors="coerce").astype(float)
        / pd.to_numeric(df.iloc[:, 2], errors="coerce").astype(float)
    )

    fig, ax = plt.subplots(figsize=(10, 6))
    scatter = ax.scatter(x, y, c=shapeFactor, cmap=cmap, s=60, edgecolor='k', linewidth=0.3)

    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(shapeFactor)
    x = x[valid]
    y = y[valid]
    shapeFactor = shapeFactor[valid]

    if len(x) > 1:
        segments = [
            [(x[i], y[i]), (x[i + 1], y[i + 1])]
            for i in range(len(x) - 1)
        ]
        color_values = np.linspace(0, 1, len(segments))
        if len(shapeFactor) > 1:
            color_values = (shapeFactor[:-1] - np.nanmin(shapeFactor)) / (
                np.nanmax(shapeFactor) - np.nanmin(shapeFactor)
            )
        line_collection = plt.matplotlib.collections.LineCollection(
            segments,
            colors=cmap(color_values),
            linewidths=1.2,
            alpha=0.8,
        )
        ax.add_collection(line_collection)
    # ax.plot(df.iloc[:, 0], df.iloc[:, 2], marker='x', label='Hull Area (Å²)')
    # ax.plot(df.iloc[:, 0], df.iloc[:, 3], marker='s', label='SASA Area (Å²)')

    ax.set_title('CCS Analysis Over Time')
    ax.set_xlabel(f'Time Step ({logTimeUnit})')
    ax.set_ylabel(f'$^{{TH}}CCS_{{{gas}}} (Å^2)$')
    ax.legend()
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label('Shape Factor ( Solvent Accessable Surface Area/ Convex Hull Area)')
    ax.grid()
    fig.tight_layout()
    output_path = "CCSAnalysisOverTime".join(csvFile.split(".")[:-1]) + ".svg"
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path
######################################################################################################
def cluster_ccs_shadow_data(
    csvFile,
    trajectory: md.Trajectory,
    gas: str,
    proteinName: str,
    visualiseRep: bool,
    savePlots: bool = False,
) -> None:
    csvDir = p.dirname(p.abspath(csvFile))
    ccsDf = pd.read_csv(csvFile)
    ccsDataArray = ccsDf.to_numpy()
    ccsArray = ccsDataArray[:, 1:3]
    clusterDir = p.abspath(p.join(csvDir, "00_clustered_pdbs"))

    if len(ccsArray) < 2:
        drLogger.log_info(
            f"Not enough CCS points to cluster for {proteinName.removeprefix('/')}. Need at least 2 points.",
            True,
        )
        return None

    max_allowed_k = min(24, len(ccsArray) - 1)
    if max_allowed_k < 2:
        bestK = 1
    else:
        bestK = drClusterizer.find_best_k_with_silhouette(ccsArray)
        bestK = max(2, min(bestK, max_allowed_k))

    kmeans = KMeans(n_clusters=bestK, random_state=42)
    cluster_labels = kmeans.fit_predict(ccsArray)
    kmeans_centers = kmeans.cluster_centers_

    os.makedirs(clusterDir, exist_ok=True)
    existingClusterPdbs = [
        p.join(clusterDir, pdbFile)
        for pdbFile in os.listdir(clusterDir)
        if p.splitext(pdbFile)[1] == ".pdb"
    ]
    for existingPdb in existingClusterPdbs:
        os.remove(existingPdb)

    drClusterizer.kmeans_clusters_to_pdb(
        ccsArray,
        bestK=bestK,
        outDir=clusterDir,
        traj=trajectory,
        protName=str(proteinName.removeprefix("/")),
    )
    if savePlots is True:
        plot_kmeans_distribution(
            ccsArray,
            cluster_labels,
            kmeans_centers,
            clusterDir,
            proteinName,
            gas,
        )

    if visualiseRep is True:
        clusterPdbs = [p.join(clusterDir, pdbFile) for pdbFile in os.listdir(clusterDir) if p.splitext(pdbFile)[1] == ".pdb"]
        for clusterPdb in clusterPdbs:
            ## calculates the accruate CCS and each frame of the shadow GIF for the representative structure
            _ = visualise_shadow(clusterPdb, gas)
    
    return kmeans
######################################################################################
def plot_kmeans_distribution(
    data: np.ndarray,
    labels: np.ndarray,
    centers: np.ndarray,
    outDir,
    proteinName: str,
    gas: str,
    ) -> None:
    """
    Plot the KMeans cluster distribution for the CCS array.

    Args:
        data (np.ndarray): The 2D data array used for clustering.
        labels (np.ndarray): Cluster labels for each point.
        centers (np.ndarray): KMeans cluster centers.
        outDir (DirectoryPath): Directory to save the plot.
        proteinName (str): Protein name used for plot title.
        gas (str): Gas type used for the CCS calculation.
    """
    if data.shape[1] < 2:
        drLogger.log_info("KMeans distribution plot requires at least two dimensions.", True)
        return

    fig, ax = plt.subplots(figsize=(10, 7))
    scatter = ax.scatter(
        data[:, 0],
        data[:, 1],
        c=labels,
        cmap='tab20',
        s=50,
        edgecolor="k",
        linewidth=0.3,
        alpha=0.85,
    )


    ax.scatter(
        centers[:, 0],
        centers[:, 1],
        c=range(len(centers)),
        cmap='tab20',
        marker="X",
        s=150,
        edgecolor="Black",
        linewidth=1.0,
        label="Cluster centers",
    )

    ax.set_title(f"KMeans CCS Distribution for {proteinName.removeprefix('/')} ({gas})")
    ax.set_xlabel("CCS (Å²)")
    ax.set_ylabel("Radius of Gyration (Å)")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.4)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("Cluster label")

    output_path = p.join(outDir, f"{proteinName.removeprefix('/')}__kmeans_ccs_distribution.svg")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


######################################################################################################
def process_ccs_frame_chunk(frameArgs: tuple):
    """
    Process a chunk of trajectory frames to compute CCS values.
    """
    frameIndices, ccsDir, trajectoryPdb, trajectoryDcd, model, gas, logTimeStep = frameArgs
    results: List[tuple] = []
    universe: mda.Universe = mda.Universe(trajectoryPdb, trajectoryDcd)

    for timeStep in frameIndices:
        universe.trajectory[timeStep]
        trajPdb = p.join(ccsDir, f"{timeStep}.pdb")
        Rg = universe.atoms.radius_of_gyration()
        universe.atoms.write(trajPdb)
        ccs= main(trajPdb, model, gas, "1000", "0.75")
        results.append((timeStep * logTimeStep, ccs, Rg))

    return results


######################################################################################################
def calculate_ccs_for_protein(proteinArgs: tuple) -> None:
    """
    Calculate CCS values for all requested simulation steps for a single protein.
    """
    (
        pdbName,
        batchConfig,
        outDir,
        model,
        gas,
        trajectoryStride,
        simulationSteps,
        visualiseRep
    ) = proteinArgs

    subprocessCpus = max(1, int(batchConfig["hardwareInfo"].get("subprocessCpus", 1)))
    simInfo = batchConfig["simulationInfo"]
    savePlots = batchConfig.get("aftercareInfo", {}).get("ccsReporterInfo", {}).get("savePlots", False)

    for simulationStep in simulationSteps:
        timeForStep = time.time()
        logInterval = next((step["logInterval"] for step in simInfo if step["stepName"] == simulationStep), None)
        if not logInterval:
            continue

        logTimeStep = float(logInterval.split(" ")[0])
        logTimeUnit = logInterval.split(" ")[1]
        simDir = p.join(outDir, pdbName, simulationStep)
        ccsDir = p.join(outDir, "00_ccs_analysis", pdbName, simulationStep)

        try:
            os.makedirs(ccsDir, exist_ok=True)
        except FileExistsError:
            pass
        except PermissionError:
            print(f"Permission denied: Unable to create {ccsDir}.")
            continue
        except Exception as exc:
            print(f"An error occurred: {exc}")
            continue

        if not p.isdir(simDir):
            print(f"Directory {simDir} does not exist. Skipping CCS calculation for {pdbName} at step {simulationStep}.")
            continue

        trajectoryDcd = p.join(simDir, "trajectory.dcd")
        trajectoryPdb = p.join(simDir, "trajectory.pdb")
        traj: md.Trajectory = md.load(trajectoryDcd, top=trajectoryPdb)
        universe: mda.Universe = mda.Universe(trajectoryPdb, trajectoryDcd)
        trajectoryLength = len(universe.trajectory)
        timeSteps = list(range(0, trajectoryLength, trajectoryStride))

        if not timeSteps:
            continue
        if p.exists(p.join(ccsDir, "ccsReport.csv")):
            drLogger.log_info(f"Skipping CCS Analysis for {pdbName}", True)
        else:
            frameChunks = np.array_split(timeSteps, max(1, min(subprocessCpus, len(timeSteps))))
            frameChunks = [list(chunk) for chunk in frameChunks if len(chunk) > 0]
            chunkArgs = [
                (chunk, ccsDir, trajectoryPdb, trajectoryDcd, model, gas, logTimeStep)
                for chunk in frameChunks
            ]

            chunkResults = process_map(process_ccs_frame_chunk, chunkArgs, max_workers=max(1, len(chunkArgs)))
            allResults = [row for chunkResult in chunkResults for row in chunkResult]

            with open(p.join(ccsDir, "ccsReport.csv"), "w", encoding="utf-8") as f:
                f.write(f"Time ({logTimeUnit}), CCS (Å²), Radius of Gyration (Å²)\n")
                for timeValue, ccs, Rg in sorted(allResults, key=lambda item: item[0]):
                    f.write(f"{timeValue}, {ccs}, {Rg}\n")

            ccsPdbs = [p.join(ccsDir, pdbFile) for pdbFile in os.listdir(ccsDir) if p.splitext(pdbFile)[1] == ".pdb"]
            for ccsPdb in ccsPdbs:
                os.remove(ccsPdb)
        timeForStep = time.time() - timeForStep
        #print(f"CCS calculation for {pdbName} at step {simulationStep} completed in {timeForStep:.2f} seconds.")
        if  savePlots is True:
            ccsplot = plot_ccs_shadow_data(p.join(ccsDir, "ccsReport.csv"), gas, logTimeUnit)
        kmeansData = cluster_ccs_shadow_data(
            p.join(ccsDir, "ccsReport.csv"),
            traj,
            gas,
            pdbName,
            visualiseRep,
            savePlots,
        )
        create_ccs_report_html(
            ccsDir=ccsDir,
            rawCsv=p.join(ccsDir, "ccsReport.csv"),
            kmeansData=kmeansData,
            gas=gas,
            logTimeUnit=logTimeUnit,
            stride= trajectoryStride
        )


######################################################################################################
def calculate_ccs_shadow(batchConfig: dict):
    """
    Uses the Shadow Screen to calculate CCS.

    Proteins are distributed across the available parallelCPU workers, and each protein's
    frames are split into chunks for the available subprocessCpus workers.
    """



    parallelCPU = max(1, int(batchConfig["hardwareInfo"].get("parallelCPU", 1)))

    aftercareInfo = batchConfig.get("aftercareInfo", {})
    ccsReporterInfo = aftercareInfo.get("ccsReporterInfo", False)
    if not ccsReporterInfo:
        return

    outDir: DirectoryPath = batchConfig["pathInfo"]["outputDir"]
    pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
    simInfo = batchConfig["simulationInfo"]
    
    pdbFilesInt = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]
    pdbNames = []
    for file in pdbFilesInt:
        pdbNames.append(file.removesuffix(".pdb").removeprefix(pdbDir+"/"))

    trajectoryStride = ccsReporterInfo["stride"]
    model = ccsReporterInfo["model"]
    gas = ccsReporterInfo["gas"]
    simulationSteps = ccsReporterInfo["steps"]
    visualiseRep= ccsReporterInfo.get("visualiseRepresentation")

    for simulationStep in simulationSteps:
        logInterval = next((step["logInterval"] for step in simInfo if step["stepName"] == simulationStep), None)
        logIntervalValue = int(logInterval.split(" ")[0])*trajectoryStride
        interval= str(logIntervalValue) + logInterval.split(" ")[1]


    proteinArgs = [
        (pdbName, batchConfig, outDir, model, gas, trajectoryStride, simulationSteps, visualiseRep)
        for pdbName in pdbNames
    ]

    process_map(
        calculate_ccs_for_protein,
        proteinArgs,
        max_workers=max(1, min(parallelCPU, len(proteinArgs))),
    )
    drMethodsWriter.add_analysis_step_to_log("CCS")
    drMethodsWriter.add_parameter_to_analysis_log("CCS", "stepNames", simulationSteps)
    drMethodsWriter.add_parameter_to_analysis_log("CCS", "gas", gas)   
    drMethodsWriter.add_parameter_to_analysis_log("CCS", "stride", trajectoryStride)
    drMethodsWriter.add_parameter_to_analysis_log("CCS", "model", model)
    return
######################################################################################################