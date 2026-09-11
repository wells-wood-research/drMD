## BASIC PYTHON LIBRARIES
import os
from os import path as p
import pandas as pd
import numpy as np
from functools import wraps
import textwrap
from textwrap import fill
import sys
from shutil import move
import shutil
import warnings
from scipy.stats import linregress
import re

## PLOTTING LIBRARIES
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
from matplotlib.ticker import MaxNLocator
import plotly.graph_objects as go
from plotly.subplots import make_subplots

## PDF LIBS
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML, CSS

## MDTRAJ AND MDANALYSIS LIBRARIES
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis import rms

## drMD LIBRARIES
from ExaminationRoom import drLogger
from UtilitiesCloset import drSplicer, drSelector, drMethodsWriter
from ExaminationRoom.drRMSD import main as drRMSD


## CLEAN CODE
from typing import Dict, Tuple, List
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

## DISABLE WARNINGS
import logging
logging.getLogger('weasyprint').setLevel(logging.ERROR)
warnings.filterwarnings('ignore')

REPORT_PLOT_HEIGHT = 440
REPORT_PLOT_WIDTH = "100%"


def _wrap_plotly_html_fragment(plot_html: str) -> str:
    """Wrap a Plotly div/script fragment in a marginless full-size HTML page."""
    return (
        "<!DOCTYPE html>\n"
        "<html><head><meta charset=\"utf-8\">"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
        "<style>"
        f"html,body{{margin:0;padding:0;width:{REPORT_PLOT_WIDTH};height:{REPORT_PLOT_HEIGHT}px;min-height:{REPORT_PLOT_HEIGHT}px;overflow:hidden;background:#021108;}}"
        f"#plot-root{{width:{REPORT_PLOT_WIDTH};height:{REPORT_PLOT_HEIGHT}px;min-height:{REPORT_PLOT_HEIGHT}px;}}"
        f"#plot-root .plotly-graph-div{{width:{REPORT_PLOT_WIDTH}!important;height:{REPORT_PLOT_HEIGHT}px!important;min-height:{REPORT_PLOT_HEIGHT}px!important;}}"
        "</style></head><body><div id=\"plot-root\">"
        f"{plot_html}"
        "</div></body></html>"
    )


def _copy_shared_css(output_dir: str, css_files: list[str]) -> None:
    """Copy shared CSS assets next to generated reports so they can be loaded by URL."""
    base_dir = p.dirname(__file__)
    output_path = p.abspath(output_dir)
    os.makedirs(output_path, exist_ok=True)
    for css_name in css_files:
        source_path = p.join(base_dir, css_name)
        target_path = p.join(output_path, css_name)
        if p.isfile(source_path):
            shutil.copy2(source_path, target_path)

######################################################################
def check_vitals(simDir: DirectoryPath,
                  vitalsFiles: Dict[str, FilePath]) -> None:
    """
    Looks at the following properties of the simulation:
        - RMSD
        - Energies (Potential, Kinetic, Total)
        - Properties (Temperature, Volume, Density)
    Checks that each of these properties have converged
    
    
    Also gathers time data for the simulation
    Plots these data and creates a report PDF

    Args:
        simDir (DirectoryPath): The directory of the simulation
        vitalsFiles (Dict[str, FilePath]): A dictionary of the vitals files
    """
    drLogger.log_info(f"Checking vitals...", True)
    ## read OpenMM reporters into dataframes
    vitalsDf: pd.DataFrame = pd.read_csv(vitalsFiles["vitals"])
    progressDf: pd.DataFrame = pd.read_csv(vitalsFiles["progress"])

    ## make energy values relative to start of simulation
    for colName in vitalsDf.columns:
        if "Energy" in colName:
            vitalsDf[colName] = vitalsDf[colName] - vitalsDf[colName][0]


    ## use mdtraj to calculate RMSD for non water and ions, get that data into a dataframe
    rmsdDf: pd.DataFrame = calculate_rmsd_mda(trajectoryDcd = vitalsFiles["trajectory"],
                                                trajectoryPdb = vitalsFiles["pdb"])
    
    ## join the dataframes 
    vitalsDf: pd.DataFrame = pd.concat([rmsdDf, vitalsDf],axis=1)

    smoothedDf: pd.DataFrame = smooth_data(vitalsDf)

    ## chunk dataframes into 1 ns blocks
    chunkedDfs: pd.DataFrame = chunk_dataframe_by_timestep(smoothedDf, 1000)
    ## check for convergance for each chunk
    propertiesConvergedInfo, plottingData = check_convergance_chunks(chunkedDfs)
    ## decide whether a property has converged
    converganceDiagnosis: Dict[str, str] = diagnose_convergance(propertiesConvergedInfo)
    ## plot all traces with line-of-best-fit for each property
    create_vitals_html(
        simDir=simDir,
        vitalsDf=vitalsDf,
        progressDf=progressDf,
        smoothedDf=smoothedDf,
        plottingData=plottingData,
        converganceDiagnosis=converganceDiagnosis,
    )

    drRMSD(simDir)
    ## tidy up reporters to avoid clutter
    tidy_up(simDir)
    drMethodsWriter.add_analysis_step_to_log("vitals")
    drMethodsWriter.add_parameter_to_analysis_log("vitals", "stepName", "all")
    drMethodsWriter.add_parameter_to_analysis_log("vitals", "convergedProperties", vitalsDf.columns.tolist())

###############################################################################################
###############################################################################################
def create_vitals_html(simDir, vitalsDf, progressDf, smoothedDf, plottingData, converganceDiagnosis):
    """
    Uses Jinja2 and Plotly to create an interactive vitals report HTML.

    Args:   
        simDir (DirectoryPath): The directory of the simulation
    """
    ## get the instruments directory path
    instrumentsDir: DirectoryPath = p.dirname(__file__)
    env: Environment = Environment(loader=FileSystemLoader(instrumentsDir), autoescape=True)
    template = env.get_template("vitals_template.html")

    stepName: str = p.basename(simDir)
    systemName: str = p.basename(p.dirname(simDir))

    metricOrder = list(plottingData.keys())
    plotDir = p.join(simDir, "00_plotly_plots")
    os.makedirs(plotDir, exist_ok=True)
    _copy_shared_css(simDir, ["report_theme.css"])
    metricPlots = []
    for metric_name in metricOrder:
        metricFigure = build_vitals_metric_figure(
            vitalsDf=vitalsDf,
            smoothedDf=smoothedDf,
            plottingData=plottingData,
            metricName=metric_name,
            diagnosis=converganceDiagnosis.get(metric_name, "N/A"),
        )
        safeMetricName = re.sub(r"[^A-Za-z0-9_.-]+", "_", metric_name).strip("_") or "metric"
        plotPath = p.join(plotDir, f"vitals_{safeMetricName}_plot.html")
        plotHtml = metricFigure.to_html(
            full_html=False,
            include_plotlyjs="cdn",
            config={
                "responsive": True,
                "displaylogo": False,
                "scrollZoom": True,
            },
        )
        plotHtml = _wrap_plotly_html_fragment(plotHtml)
        with open(plotPath, "w", encoding="utf-8") as plotFile:
            plotFile.write(plotHtml)

        metricPlots.append({
            "label": metric_name,
            "diagnosis": converganceDiagnosis.get(metric_name, "N/A"),
            "plotPath": p.relpath(plotPath, simDir),
        })

    timeDataDf = extract_time_data(vitalsDf, progressDf)
    timeSummaryRows = [
        {
            "label": row["index"],
            "value": row["Value"],
        }
        for _, row in timeDataDf.iterrows()
    ]

    vitalsCsvName = "vitals_report.csv"
    vitalsCsvPath = p.join(simDir, vitalsCsvName)

    context = {
        "vitals_data": {
            "systemName": systemName,
            "stepName": stepName,
            "nFrames": len(vitalsDf),
            "reportPath": vitalsCsvName,
            "reportName": vitalsCsvName,
            "reportCsvLink": vitalsCsvName if p.isfile(vitalsCsvPath) else None,
            "metricPlots": metricPlots,
            "timeSummaryRows": timeSummaryRows,
        }
    }
    rendered_html = template.render(context)
    
    outHtml = p.join(simDir, "vitals_report.html")
    with open(outHtml, "w", encoding="utf-8") as reportFile:
        reportFile.write(rendered_html)

    drLogger.log_info(f"Created vitals report HTML: {outHtml}", True)


def build_vitals_trace_figure(vitalsDf: pd.DataFrame,
                              smoothedDf: pd.DataFrame,
                              plottingData: dict,
                              converganceDiagnosis: dict):
    """
    Builds the interactive Plotly figure for all vitals traces.
    """
    metricOrder = list(plottingData.keys())

    if not metricOrder:
        fig = go.Figure()
        fig.update_layout(
            template="plotly_dark",
            title="Vitals Traces",
            height=400,
        )
        fig.add_annotation(
            text="No vitals trace data available.",
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            showarrow=False,
        )
        return fig

    traceColors: dict = {
        "Potential Energy (kJ/mole)": "#00FF00",
        "Kinetic Energy (kJ/mole)": "#00FFFF",
        "Total Energy (kJ/mole)": "#FFFF00",
        "Temperature (K)": "#FFA500",
        "Box Volume (nm^3)": "#FF00FF",
        "Density (g/mL)": "#FF7F50",
        "Backbone RMSD (Angstrom)": "#FFC0CB",
    }

    subplotTitles = [f"{columnName}<br>{converganceDiagnosis.get(columnName, 'N/A')}" for columnName in metricOrder]
    fig = make_subplots(
        rows=len(metricOrder),
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.02,
        subplot_titles=subplotTitles,
    )

    timeValues = vitalsDf["Time (ps)"]

    for rowIndex, columnName in enumerate(metricOrder, start=1):
        traceColor = traceColors.get(columnName, "#00FF00")

        fig.add_trace(
            go.Scatter(
                x=timeValues,
                y=vitalsDf[columnName],
                mode="lines",
                line=dict(color=traceColor, width=1),
                opacity=0.35,
                name=f"{columnName} raw",
                hovertemplate=f"{columnName} raw<br>Time: %{{x:.2f}} ps<br>Value: %{{y:.2f}}<extra></extra>",
                showlegend=False,
            ),
            row=rowIndex,
            col=1,
        )

        fig.add_trace(
            go.Scatter(
                x=timeValues,
                y=smoothedDf[columnName],
                mode="lines",
                line=dict(color=traceColor, width=2.5),
                name=f"{columnName} smoothed",
                hovertemplate=f"{columnName} smoothed<br>Time: %{{x:.2f}} ps<br>Value: %{{y:.2f}}<extra></extra>",
                showlegend=False,
            ),
            row=rowIndex,
            col=1,
        )

        for bestFitLine in plottingData[columnName].values():
            bestFitTrace = bestFitLine["dy/dx"] * bestFitLine["Time (ps)"] + bestFitLine["intercept"]
            fig.add_trace(
                go.Scatter(
                    x=bestFitLine["Time (ps)"],
                    y=bestFitTrace,
                    mode="lines",
                    line=dict(color="#F5F5F5", width=1, dash="dash"),
                    opacity=0.55,
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=rowIndex,
                col=1,
            )

        fig.update_yaxes(title_text=columnName, row=rowIndex, col=1, automargin=True)

    fig.update_xaxes(title_text="Time (ps)", row=len(metricOrder), col=1)
    fig.update_layout(
        template="plotly_dark",
        title="Vitals Traces",
        height=max(1800, 240 * len(metricOrder)),
        margin=dict(l=70, r=30, t=90, b=60),
        hovermode="x unified",
        showlegend=False,
    )

    return fig


def build_vitals_metric_figure(vitalsDf: pd.DataFrame,
                              smoothedDf: pd.DataFrame,
                              plottingData: dict,
                              metricName: str,
                              diagnosis: str):
    """Build a separate Plotly figure for one vitals metric."""
    traceColors: dict = {
        "Potential Energy (kJ/mole)": "#00FF00",
        "Kinetic Energy (kJ/mole)": "#00FFFF",
        "Total Energy (kJ/mole)": "#FFFF00",
        "Temperature (K)": "#FFA500",
        "Box Volume (nm^3)": "#FF00FF",
        "Density (g/mL)": "#FF7F50",
        "Backbone RMSD (Angstrom)": "#FFC0CB",
    }

    fig = go.Figure()
    timeValues = vitalsDf["Time (ps)"]
    traceColor = traceColors.get(metricName, "#00FF00")

    fig.add_trace(
        go.Scatter(
            x=timeValues,
            y=vitalsDf[metricName],
            mode="lines",
            line=dict(color=traceColor, width=1),
            opacity=0.35,
            name=f"{metricName} raw",
            hovertemplate=f"{metricName} raw<br>Time: %{{x:.2f}} ps<br>Value: %{{y:.2f}}<extra></extra>",
            showlegend=False,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=timeValues,
            y=smoothedDf[metricName],
            mode="lines",
            line=dict(color=traceColor, width=2.5),
            name=f"{metricName} smoothed",
            hovertemplate=f"{metricName} smoothed<br>Time: %{{x:.2f}} ps<br>Value: %{{y:.2f}}<extra></extra>",
            showlegend=False,
        )
    )

    for bestFitLine in plottingData.get(metricName, {}).values():
        bestFitTrace = bestFitLine["dy/dx"] * bestFitLine["Time (ps)"] + bestFitLine["intercept"]
        fig.add_trace(
            go.Scatter(
                x=bestFitLine["Time (ps)"],
                y=bestFitTrace,
                mode="lines",
                line=dict(color="#F5F5F5", width=1, dash="dash"),
                opacity=0.55,
                hoverinfo="skip",
                showlegend=False,
            )
        )

    fig.update_layout(
        template="plotly_dark",
        title=f"{metricName}<br>{diagnosis}",
        xaxis_title="Time (ps)",
        yaxis_title=metricName,
        height=REPORT_PLOT_HEIGHT,
        margin=dict(l=60, r=30, t=60, b=60),
        hovermode="x unified",
        showlegend=False,
    )
    return fig


def build_time_summary_figure(vitalsDf: pd.DataFrame, progressDf: pd.DataFrame):
    """
    Builds a Plotly table summarizing the timing information for the run.
    """
    timeDataDf: pd.DataFrame = extract_time_data(vitalsDf, progressDf)

    fig = go.Figure(
        data=[
            go.Table(
                header=dict(
                    values=["Metric", "Value"],
                    fill_color="#444444",
                    font=dict(color="#00FF00", size=14),
                    align="left",
                ),
                cells=dict(
                    values=[timeDataDf["index"].tolist(), timeDataDf["Value"].tolist()],
                    fill_color="#2a2a2a",
                    font=dict(color="#f0f0f0", size=13),
                    align="left",
                    height=34,
                ),
            )
        ]
    )

    fig.update_layout(
        template="plotly_dark",
        title="Run Time Summary",
        margin=dict(l=20, r=20, t=70, b=20),
        height=REPORT_PLOT_HEIGHT,
    )

    return fig

###############################################################################################

def smooth_data(vitalsDf: pd.DataFrame, windowSize = 1000) -> pd.DataFrame:
    """
    Smooths the data in the dataframe

    Args:
        vitalsDf (pd.DataFrame): The dataframe to smooth
        windowSize (int, optional): The window size for the rolling average. Defaults to 10.

    Returns:
        pd.DataFrame: The smoothed dataframe
    """
    ## initialise data column names
    dataColumnNames: list = ["Total Energy (kJ/mole)", "Potential Energy (kJ/mole)", "Kinetic Energy (kJ/mole)",
                    "Temperature (K)", "Box Volume (nm^3)", "Density (g/mL)", "Backbone RMSD (Angstrom)"]
    ## copy the dataframe
    smoothedDf = vitalsDf.copy()
    ## smooth the data
    for columnName in dataColumnNames:
        smoothedDf[columnName] = smoothedDf[columnName].rolling(window=windowSize, min_periods=1).median()

    return smoothedDf
###############################################################################################
def plot_system_info(systemName: str, stepName: str, outDir):
    """
    Plots system information into a table
    
    Args:
        systemName (str): The name of the system
        stepName (str): The name of the step
        outDir (str): The output directory
    """

    # Set up plot size and convert cm to inches
    widthInch: float = 25 / 2.54
    heightInch: float = 3 / 2.54  # Increased height for wrapping

    # Prepare data for a three-column table
    data = [
        ["VITALS", "SYSTEM NAME", "STEP NAME"],
        ["REPORT", systemName, stepName]
    ]

    # Set up plot
    fig, ax = plt.subplots(figsize=(widthInch, heightInch))
    fig.patch.set_facecolor('#1a1a1a')  # Much darker grey
    ax.axis('off')
    
    # Create the table
    table = ax.table(cellText=data,
                      cellLoc='left', loc='center',
                      colWidths=[0.2, 0.4, 0.4])
    table.auto_set_font_size(False)

    # Adjust font size and row height
    for (i, j), cell in table.get_celld().items():
        if i == 0:  # Header row
            if j == 0:  # "VITALS"
                cell.get_text().set_fontsize(26)
            else:
                cell.get_text().set_fontsize(10)

            cell.set_height(0.25)
        else:  # Data row
            if j == 0:  # "REPORT"
                cell.get_text().set_fontsize(26)
            else:
                cell.get_text().set_fontsize(16)
                # Wrap text for systemName and stepName
                cell.get_text().set_text(fill(data[i][j],30))
            cell.set_height(0.4)


        cell.set_edgecolor('none')  # Disable edges
        cell.set_facecolor('#1a1a1a')
        cell.get_text().set_color("yellow")

    # Adjust layout to fill space
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # Save the plot as a PNG image with reduced border
    savePng = p.join(outDir, "systemInfo.png")
    plt.savefig(savePng, bbox_inches="tight", pad_inches=0, facecolor='#1a1a1a', dpi=300)
    plt.close()

    return savePng
###############################################################################################
def diagnose_convergance(propertiesConvergedInfo) -> Dict[str, str]:
    """
    Decides whether properties of the simulation have converged
    Creates a dict of ticks, crosses and "converged in last N ns"

    Args:
        propertiesConvergedInfo (dict): A dictionary of the properties that have converged

    Returns:
        converganceDiagnosis (dict): A dictionary of the properties that have converged
    """

    ## init an empty dictionary
    converganceDiagnosis: dict = {}

    ## check if all properties converged
    for columnName in propertiesConvergedInfo:
        ## check for all converged - give this a tick ✔
        if all(converged for converged in propertiesConvergedInfo[columnName].values()):
            converganceDiagnosis[columnName] = "✔"
        ## check for all not converged - give this a cross ✘
        elif all(not converged for converged in propertiesConvergedInfo[columnName].values()):
            converganceDiagnosis[columnName] = "✘"
        ## for a property that has not converged, check if it has converged in the last N nanoseconds
        else:
            n = len(propertiesConvergedInfo[columnName].values())
            while True:
                n -= 1
                if n == 0:
                    converganceDiagnosis[columnName] = "✘"
                    break
                if all(converged for converged in list(propertiesConvergedInfo[columnName].values())[-n:]):
                    converganceDiagnosis[columnName] = f"Converged\n in last \n{str(n)} ns"
                    break

    return converganceDiagnosis

###############################################################################################
def check_convergance_chunks(chunkedDfs: list) -> Tuple[dict,dict]:
    """
    Loops through each dataframe in chunkedDfs and checks if the properties have converged

    Args:
        chunkedDfs (list): A list of dataframes

    Returns:
        converganceDiagnosis (dict): A dictionary of the properties that have converged
        plottingData (dict): A dictionary of the plotting data
    """


    ## initiate a dict of gradiant tolerances for simulation properties
    ##TODO: maybe tweak these...
    columnToleranceInfo: dict = {
                "Potential Energy (kJ/mole)": {"gradientTol": 1000, "unit": "kJ/mole"},
                "Kinetic Energy (kJ/mole)": {"gradientTol": 1000, "unit": "kJ/mole"},
                "Total Energy (kJ/mole)": {"gradientTol": 1000, "unit": "kJ/mole"},
                "Temperature (K)": {"gradientTol": 1, "unit": "K"},
                "Box Volume (nm^3)": {"gradientTol": 1, "unit": "nm^3"},
                "Density (g/mL)": {"gradientTol": 0.1, "unit": "g/mL"},
                "Backbone RMSD (Angstrom)": {"gradientTol": 1, "unit": "A"},
                  }

    ## init some empty dicts to store data
    propertiesConvergedInfo: dict = {}
    plottingData: dict = {}
    ## look through each property
    for columnName in columnToleranceInfo:
        ## make a new entries to store data
        propertiesConvergedInfo[columnName] = {}
        plottingData[columnName] = {}
        ## loop through chunked dataframes
        for chunkIndex, chunkedDf in enumerate(chunkedDfs):
            ## calculate the line of best fit
            slope, intercept, r_squared = calculate_line_of_best_fit(chunkedDf["Time (ps)"], chunkedDf[columnName])
            ## add results to the plotting data dictionary
            plottingData[columnName][chunkIndex] = {"dy/dx": slope,
                                                     "intercept": intercept,
                                                       "r_squared": r_squared,
                                                       "Time (ps)": chunkedDf["Time (ps)"]}
            ## get slope per nanosecond
            slopePerNs = slope * 1000
            ## check if the slope is within the tolerance
            if slopePerNs > columnToleranceInfo[columnName]["gradientTol"]:
                chunkConverged = False
            else:
                chunkConverged = True
            ## add result to the propertiesConvergedInfo dictionary
            propertiesConvergedInfo[columnName][chunkIndex] = chunkConverged

    return propertiesConvergedInfo, plottingData
###############################################################################################
def chunk_dataframe_by_timestep(vitalsDf: pd.DataFrame, timeChunkSize: int = 1000) -> list[pd.DataFrame]:
    """
    Chops a dataframe into chunks based on time columns (use 1 ns chunks as default)

    Args:
        vitalsDf (pd.DataFrame): The vitals dataframe
        timeChunkSize (int, optional): The size of the time chunk. Defaults to 1000.

    Returns:
        chunkedDfs (list[pd.DataFrame]): A list of dataframes
    
    """


    ## init an empty list to store the chunked dataframes
    chunkedDfs = []
    ## copy timeChunkSize to timeChunkIncrement
    timeChunkIncrement = timeChunkSize
    ## init a counter to keep track chunk location
    lastTimeChunk = 0
    ## chop dataframe into chunks
    while True:
        chunkDf = vitalsDf[(vitalsDf["Time (ps)"] <= timeChunkSize) & (vitalsDf["Time (ps)"] > lastTimeChunk)]
        if len(chunkDf) == 0:
            break
        chunkedDfs.append(chunkDf)
        lastTimeChunk = timeChunkSize
        timeChunkSize += timeChunkIncrement
    return chunkedDfs

###############################################################################################
def calculate_line_of_best_fit(x, y):
    """
    Calculates the line of best fit for a set of x and y values

    Args:
        x (list): The x values
        y (list): The y values

    Returns:
        slope (float): The slope of the line of best fit
        intercept (float): The intercept of the line of best fit
        r_squared (float): The r squared value
    """
    ## calculate the line of best fit
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    ## calculate r squared
    r_squared = r_value**2
    return slope, intercept, r_squared



######################################################################
def tidy_up(simDir: DirectoryPath):
    """
    Tidies up the reporters in the simulation directory

    Args:
        simDir (DirectoryPath): The directory of the simulation
    """
    ## make a new directory to tidy up reporters and png files
    tidyDir: DirectoryPath = p.join(simDir, "00_reporters_and_plots")
    os.makedirs(tidyDir, exist_ok=True)
    ## move all reporters and png files to the new directory
    for file in os.listdir(simDir):
        if p.splitext(file)[1] in [".csv", ".png"]:
            move(p.join(simDir, file), p.join(tidyDir, file))

    plt.close("all")
######################################################################
def check_up_handler():
    """
    Decorator for the check_up function

    Returns:
        decorator (function): The decorated function
    """
    def decorator(simulationFunction):
        @wraps(simulationFunction)
        def wrapper(*args, **kwargs):
            saveFile: FilePath = simulationFunction(*args, **kwargs)
            try:
                vitalsFiles, simDir = find_vitals_files(simInfo=kwargs["sim"],
                                                     outDir=kwargs["outDir"], 
                                                     pdbFile=kwargs["refPdb"],
                                                     config=kwargs["config"])

                check_vitals(simDir = simDir,
                            vitalsFiles = vitalsFiles)
            except FileNotFoundError as e:
                drLogger.log_info(f"Error running checkup: File not found: {e}", True, True)
                raise e
            except Exception as e:
                drLogger.log_info(f"Error running checkup: {e}", True,True)
                raise e
            return saveFile
        return wrapper
    return decorator
######################################################################
def find_vitals_files(simInfo: Dict,
                       outDir: DirectoryPath,
                       pdbFile: FilePath,
                       config: Dict) -> Tuple[Dict[str, FilePath], DirectoryPath]:
    """
    Looks through the simulation directory to find the vitals reporter files

    Args:
        simInfo (Dict): The simulation info
        outDir (DirectoryPath): The output directory
        pdbFile (FilePath): The reference pdb file
        config (Dict): The config dictionary

    Returns:
        vitalsFiles (Dict[str, FilePath]): A dictionary of the vitals files
        simDir (DirectoryPath): The simulation directory        
    
    """
    
    ## get the simulation directory
    simDir: DirectoryPath= p.join(outDir, simInfo["stepName"])

    ## check to see if multiple partial trajectories exist
    trajectoryDcds: list = [p.join(simDir, file) for file in os.listdir(simDir) if file.endswith(".dcd")]
    ## if more than one trajectory is found, merge partial outputs before running the health check
    if len(trajectoryDcds) > 1:
        drSplicer.merge_partial_outputs(simDir = simDir,
                                            pdbFile=pdbFile,
                                            config=config,
                                            simInfo = simInfo)

    ## find vitals reporter file
    vitalsReport: FilePath = p.join(simDir, "vitals_report.csv")
    if not p.isfile(vitalsReport):
        raise FileNotFoundError(f"->\tReporter file not found at {vitalsReport}")
    ## find progress reporter file
    progressReport: FilePath = p.join(simDir, "progress_report.csv")
    if not p.isfile(progressReport):
        raise FileNotFoundError(f"->\tReporter file not found at {progressReport}")
    ## find trajectory file
    trajectoryDcd: FilePath = p.join(simDir, "trajectory.dcd")
    if not p.isfile(trajectoryDcd):
        raise FileNotFoundError(f"->\Trajectory file not found at {trajectoryDcd}")
    
    ## find the pdb file
    pdbFile = p.join(simDir, "trajectory.pdb")
    if not p.isfile(pdbFile):
        raise FileNotFoundError(f"->\tPDB file not found at {simDir}")
    
    ## collect files into a dictionary
    vitalsFiles = {"vitals": vitalsReport, "progress": progressReport, "trajectory": trajectoryDcd, "pdb": pdbFile}
    
    return vitalsFiles, simDir
######################################################################
def convert_seconds(seconds: int) -> str:
    """
    Converts seconds as an int to HH:MM:SS format

    Args:
        seconds (int): The number of seconds
    """
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:02}:{:02}:{:02}".format(int(hours), int(minutes), int(seconds))
######################################################################
def extract_time_data(vitalsDf: pd.DataFrame, progressDf: pd.DataFrame) -> pd.DataFrame:
    """
    Reads vitals and progress data to extract key time-related data
    
    Args:
        vitalsDf (pd.DataFrame): The vitals dataframe
        progressDf (pd.DataFrame): The progress dataframe

    Returns:
        timeDf (pd.DataFrame): The time dataframe
    """

    # get elapsed time
    timeElapsed: int = progressDf["Elapsed Time (s)"].tail(1).values[0]
    timeElapsed: str = convert_seconds(timeElapsed)
    # get average speed
    averageSpeed = str(round(progressDf["Speed (ns/day)"].mean(),0))
    # get nSteps
    nSteps = str(round(vitalsDf["#\"Step\""].tail(1).values[0], 0))
    # get simulation duration
    duration = str(round(vitalsDf["Time (ps)"].tail(1).values[0], 0))
    # create a DataFrame
    timeDf = pd.DataFrame({'Elapsed Time (HH:MM:SS)': [timeElapsed],
                            'Average Speed (ns/day)': [averageSpeed],
                            'Total Steps': [nSteps],
                            'Simulation Duration (ps)': [duration]})

    timeDf: pd.DataFrame = timeDf.transpose()
    timeDf.columns = ['Value']  # rename the column
    timeDf.reset_index(level=0, inplace=True)  # reset the index

    return timeDf

######################################################################
def plot_time_data(timeDf: pd.DataFrame, outDir: str):
    """
    Plots time data into a table
    
    Args:
        timeDf (pd.DataFrame): The time dataframe
        outDir (str): The output directory
    """

    # Set up plot size and convert cm to inches
    widthInch: float = 7.75 / 2.54
    heightInch: float = 18.7 / 2.54

    # Prepare data for one-column table with adjustable gaps
    data = []
    for idx, row in timeDf.iterrows():
        data.append([row['index'].upper()])
        data.append([row['Value']])
        data.append([''])  # Add an empty row for the gap
    data = data[:-1]

    # Set up plot
    fig, ax = plt.subplots(figsize=(widthInch, heightInch))
    fig.patch.set_facecolor('#1a1a1a')  # Much darker grey
    brightGreen: str = '#00FF00'  # Brighter green
    ax.axis('off')
    
    # Create the table without column labels
    table = ax.table(cellText=data, cellLoc='center', loc='center', colLabels=None)
    table.auto_set_font_size(False)

    # Adjust font size and row height for index, value, and gap rows
    for i, cell in enumerate(table.get_celld().values()):
        if i % 3 == 0:  # Index rows
            cell.get_text().set_fontsize(12)
            cell.set_height(0.05)  # Adjust height for index
        elif i % 3 == 1:  # Value rows
            cell.get_text().set_fontsize(30)
            cell.set_height(0.1)  # Adjust height for value
        else:  # Gap rows
            cell.get_text().set_text('')
            cell.set_height(0.1)  # Adjust height for gap

        cell.set_edgecolor('none')  # Disable edges
        cell.set_facecolor('#1a1a1a')
        cell.get_text().set_color("yellow")

    # Adjust layout to fill space
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # Save the plot as a PNG image with reduced border
    savePng = p.join(outDir, "time_info.png")
    plt.savefig(savePng, bbox_inches="tight", pad_inches=0, facecolor='#1a1a1a', dpi=300)
    plt.close()

    return savePng
#########################################################################################################
def calculate_rmsd_mda(trajectoryDcd: FilePath, trajectoryPdb: FilePath) -> pd.DataFrame:
    """
    Uses MDAnalysis to calculate backbone RMSD

    Args:
        trajectoryDcd (FilePath): The trajectory DCD file
        trajectoryPdb (FilePath): The trajectory PDB file

    Returns:
        rmsdDf (pd.DataFrame): The RMSD dataframe
    """


    universe: mda.Universe = mda.Universe(trajectoryPdb, trajectoryDcd)

    rmsdCalculation = mda.analysis.rms.RMSD(universe, select="backbone")

    rmsdCalculation.run()


    rmsdDf = pd.DataFrame(rmsdCalculation.rmsd)
    rmsdDf.columns = ["Frame Index", "Timestep (ps)", "Backbone RMSD (Angstrom)"]

    rmsdDf.drop(columns=["Frame Index"], inplace=True)
    rmsdDf['Backbone RMSD (Angstrom)'] = rmsdDf['Backbone RMSD (Angstrom)'].round(2)

    return rmsdDf
#########################################################################################################
def calculate_ccs_shadow(trajectoryDcd: FilePath, trajectoryPdb: FilePath) -> pd.DataFrame:
    """
    Uses the Shadow Screen to calculate CCS

    Args:
        trajectoryDcd (FilePath): The trajectory DCD file
        trajectoryPdb (FilePath): The trajectory PDB file

    Returns:
        ccsDf (pd.DataFrame): The CCS dataframe
    """

    universe: mda.Universe = mda.Universe(trajectoryPdb, trajectoryDcd)
    trajectoryLength= len(universe.trajectory)
#########################################################################################################
def plot_traces(vitalsDf: pd.DataFrame, plottingData: dict, convergedDiagnosis: dict, outDir: DirectoryPath) -> FilePath:
    """
    Creates the traces plot for the vitals report

    Args:
        vitalsDf (pd.DataFrame): The vitals dataframe
        plottingData (dict): The plotting data
        convergedDiagnosis (dict): The converged diagnosis
        outDir (DirectoryPath): The output directory

    Returns:
        savePng (FilePath): The path to the saved plot
    
    """
    ## init some colors to be used
    darkGrey: str = '#1a1a1a'
    white :str = '#FFFFFF'

    traceColors: dict = {
        "brightRed": '#FF0000',
        "brightGreen": '#00FF00',
        "brightYellow": '#FFFF00',
        "brightCyan": '#00FFFF',
        "brightMagenta": '#FF00FF',
        "brightPink": '#FFC0CB',
        "brightOrange": '#FFA500',
    }

    # Set up plot size and convert cm to inches (26 x 26 cm)
    plt.figure(figsize=(26 / 2.54, 26 / 2.54))
    ## create the plot with 7 rows and 2 columns
    fig, axes = plt.subplots(7, 2, gridspec_kw={'width_ratios': [3, 1]})
    ## set the background color
    fig.patch.set_facecolor(darkGrey)

    # Adjust spacing between subplots
    plt.subplots_adjust(hspace=0.7, wspace=-0.2)

    # Plot the traces
    for i, (ax, columnName, traceColor) in enumerate(zip(axes[:, 0], plottingData, traceColors.values())):
        ## set the background color for each subplot
        ax.set_facecolor(darkGrey)
        ## make a glowing line around the trace
        for n in range(1, 6):
            ax.plot(vitalsDf["Time (ps)"], vitalsDf[columnName], linestyle='-', color=traceColor, linewidth=1+n, alpha=0.1)

        ## plot the trace
        ax.plot(vitalsDf["Time (ps)"], vitalsDf[columnName], linestyle='-', color=traceColor, linewidth=1)
        ## set the title, put it in the correct place
        ax.set_title(columnName, fontsize=8, color=traceColor, loc="left", x=0.035, y=0.9)
        ## deal with the x axis number of ticks
        ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=3))
        ## set the x axis range 
        x_min = vitalsDf["Time (ps)"].min()
        x_max = vitalsDf["Time (ps)"].max()
        x_range = x_max - x_min
        padding = x_range * 0.05  
        ax.set_xlim(x_min - padding, x_max + padding)

        ## plot x axis for last plot
        if i == 6:
            ax.spines['bottom'].set_color(white)
            ax.spines['bottom'].set_position(('outward', 10))  
            ax.set_facecolor(darkGrey)
            ax.set_xlabel("Time (ps)", fontsize=8, color=white)
            ax.tick_params(axis='x', colors=white)
            # Set x-ticks to include the first and last data points
            ax.set_xticks(np.linspace(x_min, x_max))  
            # Set x-ticks to integer
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))

            for label in ax.get_xticklabels():
                label.set_fontsize(8)
                label.set_color(white)
        ## remove all x-ticks and x-axis for other plots
        else:
            ax.spines['bottom'].set_visible(False)
            ax.set_xticks([])
            ax.set_xticklabels([])

        ## hide top spines
        ax.spines['top'].set_visible(False)
        ## set left and right spines and y-ticks to traceColor 
        ax.spines['left'].set_color(traceColor)
        ax.tick_params(axis='y', colors=traceColor)
        ax.spines['right'].set_color(traceColor)

        ## set labels color and fontsize
        for label in ax.get_yticklabels():
            label.set_fontsize(8)
            label.set_color(traceColor)

        ## plot best fit lines on top of the main trace
        bestFitLines = plottingData[columnName]
        for _, bestFitLine in bestFitLines.items():
            ## use slope and intercept to plot the best fit line
            bestFitTrace = bestFitLine["dy/dx"] * bestFitLine["Time (ps)"] + bestFitLine["intercept"]
            ax.plot(bestFitLine["Time (ps)"], bestFitTrace, linestyle="--", color=white, linewidth=0.5)

        ## put a table in displaying convergance diagnosis
        ax_table = axes[i, 1]
        ax_table.axis('off')  # Turn off the axis   
        value = convergedDiagnosis[columnName]
        if len(value) == 1:
            fontSize = 12
        else:
            fontSize = 8
        ax_table.text(0.45, 0.6, str(value), fontsize=fontSize, color=traceColor, ha='left', va='center', wrap=True)

    ## add title for property converged table
    fig.text(0.77, 0.92, 'Property\nConverged?', fontsize=8, color='white', ha='left', va= "bottom")
    line = plt.Line2D([0.77, 0.87], [0.91, 0.91], color='white', linewidth=0.8, linestyle='-')
    fig.add_artist(line)

    ## disable extra axis stuff and spines
    for ax_table in axes[:, 1]:
        ax_table.axis('off')  # Turn off the axis
        ax_table.set_xticks([])  # Remove x-ticks
        ax_table.set_xticklabels([])  # Remove x-tick labels
        ax_table.spines['top'].set_visible(False)  # Hide top spine
        ax_table.spines['bottom'].set_visible(False)  # Hide bottom spine
        ax_table.spines['left'].set_visible(False)  # Hide left spine
        ax_table.spines['right'].set_visible(False)  # Hide right spine



    outPng = p.join(outDir, "convergance_check_tests.png")
    plt.savefig(outPng, bbox_inches=None, dpi=300)
    plt.close()

    return outPng

######################################################################
if __name__ == "__main__":

    ## get current working directory and import it
    cwd = os.getcwd()
    sys.path.append(cwd)


    apoDir = "/home/esp/scriptDevelopment/drAnalysis/PETase_MD_for_drMD_paper/04_PET_proj_outputs/6eqe_PET-tetramer_1"
    holoDirs = "/home/esp/scriptDevelopment/drAnalysis/PETase_MD_for_drMD_paper/04_PET_proj_outputs/6eqe_1"

    apoRunDirs = sorted([p.join(apoDir, dir) for dir in os.listdir(apoDir) if not "energy" in dir and not "prep" in dir])

    holoDirs = sorted([p.join(holoDirs, dir) for dir in os.listdir(holoDirs) if not "energy" in dir and not "prep" in dir])

    allRunDirs = apoRunDirs + holoDirs

    for simDir in allRunDirs:

        ## find reporter csv files
        tidyDir = p.join(simDir, "00_reporters_and_plots")
        if p.isdir(tidyDir):
            vitalsCsv = p.join(tidyDir, "vitals_report.csv")
            progressCsv = p.join(tidyDir, "progress_report.csv")
        else:
            vitalsCsv = p.join(simDir,"vitals_report.csv")
            progressCsv = p.join(simDir,"progress_report.csv")
        ## find structure files
        trajectoryDcd = p.join(simDir,"trajectory.dcd")
        pdbFile = p.join(simDir,"trajectory.pdb")

        vitals = {"vitals": vitalsCsv, "progress": progressCsv, "trajectory": trajectoryDcd, "pdb": pdbFile}

        check_vitals(simDir,
                    vitals)
        


