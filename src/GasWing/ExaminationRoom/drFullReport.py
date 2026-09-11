## BASIC PYTHON LIBRARIES
import os
from os import path as p
from shutil import copy2

## drMD LIBRARIES
from ExaminationRoom import drLogger
from UtilitiesCloset import drListInitiator

## CLEAN CODE
from typing import Any, List, Dict
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath
from jinja2 import Environment, FileSystemLoader


######################################################################################################
def full_report_handler(batchConfig: dict) -> None:
    """
    Handler for generating a full report of the batch simulation.

    Args:
        batchConfig (dict): The batch configuration dictionary.
    Returns:
        None
    """

    aftercareInfo: Dict = batchConfig.get("aftercareInfo", {})
    fullReport: bool = aftercareInfo.get("fullReport", False)
    if not fullReport:
        return

    pathInfo: Dict = batchConfig["pathInfo"]
    outDir: DirectoryPath = pathInfo["outputDir"]
    notRunDirs = set(drListInitiator.get_not_a_run_dir())

    proteinDirs: List[DirectoryPath] = []
    for entryName in os.listdir(outDir):
        entryPath = p.join(outDir, entryName)
        if entryName in notRunDirs or not p.isdir(entryPath) or entryName.startswith("00_"):
            continue
        proteinDirs.append(entryPath)

    if not proteinDirs:
        drLogger.log_info("No protein simulation directories found for batch report generation.", True)
        return

    fullReportsRoot: DirectoryPath = p.join(outDir, "00_full_reports")
    os.makedirs(fullReportsRoot, exist_ok=True)

    templateDir: DirectoryPath = p.dirname(__file__)
    env = Environment(loader=FileSystemLoader(templateDir), autoescape=True)
    reportTemplate = env.get_template("full_report_template.html")
    indexTemplate = env.get_template("full_report_index_template.html")

    _copy_css_assets(fullReportsRoot, ["full_report_template.css"])

    indexEntries: List[Dict[str, str]] = []

    for proteinDir in proteinDirs:
        proteinName: str = p.basename(proteinDir)
        htmlReports = _discover_html_reports(outDir, proteinName)
        if not htmlReports:
            continue

        proteinReportDir: DirectoryPath = p.join(fullReportsRoot, proteinName)
        stepEntries = _group_html_reports_by_step(outDir, proteinName, proteinReportDir, htmlReports)

        os.makedirs(proteinReportDir, exist_ok=True)
        _copy_css_assets(proteinReportDir, ["full_report_template.css"])
        reportPath = p.join(proteinReportDir, "full_report.html")
        with open(reportPath, "w", encoding="utf-8") as reportFile:
            reportFile.write(
                reportTemplate.render(
                    proteinName=proteinName,
                    stepEntries=stepEntries,
                    reportCount=sum(stepEntry["reportCount"] for stepEntry in stepEntries),
                    batchDir=p.relpath(outDir, p.dirname(outDir)),
                )
            )

        indexEntries.append(
            {
                "proteinName": proteinName,
                "reportRelPath": p.relpath(reportPath, fullReportsRoot),
                "reportCount": str(sum(stepEntry["reportCount"] for stepEntry in stepEntries)),
            }
        )

        drLogger.log_info(f"Created full report HTML: {reportPath}", True)

    if not indexEntries:
        drLogger.log_info("No HTML simulation reports were found to collate into a full batch report.", True)
        return

    indexPath = p.join(fullReportsRoot, "full_report_index.html")
    with open(indexPath, "w", encoding="utf-8") as indexFile:
        indexFile.write(indexTemplate.render(indexEntries=indexEntries))

    drLogger.log_info(f"Created batch full report HTML: {indexPath}", True)


def _copy_css_assets(target_dir: DirectoryPath, css_files: List[str]) -> None:
    """Copy external CSS files next to generated pages so relative stylesheet links continue to work."""
    os.makedirs(target_dir, exist_ok=True)
    template_dir = p.dirname(__file__)
    for css_name in css_files:
        source_path = p.join(template_dir, css_name)
        target_path = p.join(target_dir, css_name)
        if p.isfile(source_path):
            copy2(source_path, target_path)


def _discover_html_reports(outDir: DirectoryPath, proteinName: str) -> List[FilePath]:
    """Find all HTML report files associated with a single protein."""

    discoveredReports: List[FilePath] = []
    proteinRunRoot = p.join(outDir, proteinName)
    if p.isdir(proteinRunRoot):
        for root, dirs, files in os.walk(proteinRunRoot):
            dirs[:] = [d for d in dirs if not d.startswith("00_")]
            if "00_full_reports" in p.normpath(root).split(os.sep):
                continue

            for fileName in files:
                if fileName.endswith(".html"):
                    discoveredReports.append(p.join(root, fileName))

    ccsProteinRoot = p.join(outDir, "00_ccs_analysis", proteinName)
    if p.isdir(ccsProteinRoot):
        for root, dirs, files in os.walk(ccsProteinRoot):
            dirs[:] = [d for d in dirs if not d.startswith("00_")]
            for fileName in files:
                if fileName.endswith(".html"):
                    discoveredReports.append(p.join(root, fileName))

    return sorted(set(discoveredReports))


def _group_html_reports_by_step(
    outDir: DirectoryPath,
    proteinName: str,
    proteinReportDir: DirectoryPath,
    htmlReports: List[FilePath],
) -> List[Dict[str, Any]]:
    """Group report files by their corresponding simulation step."""

    stepGroups: Dict[str, Dict[str, Any]] = {}

    for reportPath in htmlReports:
        sourceRelPath = p.relpath(reportPath, outDir)
        stepName = _extract_step_name(sourceRelPath, proteinName)
        if not stepName:
            continue

        stepReportDir: DirectoryPath = p.join(proteinReportDir, stepName)
        bundleName = p.splitext(p.basename(reportPath))[0]
        bundleRoot: DirectoryPath = p.join(stepReportDir, "00_plots_and_data", bundleName)
        os.makedirs(bundleRoot, exist_ok=True)
        bundleReportPath: DirectoryPath = p.join(bundleRoot, p.basename(reportPath))
        _copy_report_tree(reportPath, bundleReportPath)

        copiedReportPath: DirectoryPath = p.join(stepReportDir, p.basename(reportPath))
        _write_report_wrapper(
            copiedReportPath,
            p.relpath(bundleReportPath, stepReportDir),
        )

        reportEntry = {
            "title": _humanize_report_title(p.basename(reportPath)),
            "sourceRelPath": sourceRelPath,
            "originalRelPath": p.relpath(reportPath, proteinReportDir),
            "copyRelPath": p.relpath(copiedReportPath, proteinReportDir),
        }

        if stepName not in stepGroups:
            stepGroups[stepName] = {
                "stepName": stepName,
                "reports": [],
            }

        stepGroups[stepName]["reports"].append(reportEntry)

    stepEntries: List[Dict[str, Any]] = []
    for stepName in sorted(stepGroups.keys(), key=_step_sort_key):
        reports = sorted(
            stepGroups[stepName]["reports"],
            key=lambda entry: (entry["title"], entry["sourceRelPath"]),
        )
        stepEntries.append(
            {
                "stepName": stepName,
                "reportCount": len(reports),
                "reports": [
                    {**reportEntry, "isActive": reportIndex == 0}
                    for reportIndex, reportEntry in enumerate(reports)
                ],
                "isActive": len(stepEntries) == 0,
            }
        )

    return stepEntries


def _copy_report_tree(reportPath: str, copiedReportPath: str) -> None:
    """Copy a report and any sibling or nested files it references."""

    sourceDir = p.dirname(reportPath)
    destinationDir = p.dirname(copiedReportPath)
    os.makedirs(destinationDir, exist_ok=True)

    if p.isdir(sourceDir):
        for root, dirs, files in os.walk(sourceDir):
            relRoot = p.relpath(root, sourceDir)
            targetRoot = destinationDir if relRoot == "." else p.join(destinationDir, relRoot)
            os.makedirs(targetRoot, exist_ok=True)
            for fileName in files:
                if p.splitext(fileName)[1].lower() not in {".html", ".csv", ".bin", ".gif", ".css"}:
                    continue
                sourceFile = p.join(root, fileName)
                destinationFile = p.join(targetRoot, fileName)
                copy2(sourceFile, destinationFile)
        return

    copy2(reportPath, copiedReportPath)


def _write_report_wrapper(wrapperPath: str, targetRelPath: str) -> None:
    """Create a lightweight HTML wrapper that forwards to the copied report bundle."""

    os.makedirs(p.dirname(wrapperPath), exist_ok=True)
    wrapper_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta http-equiv="refresh" content="0; url={targetRelPath}">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Report</title>
</head>
<body>
    <p>Opening report...</p>
</body>
</html>
"""
    with open(wrapperPath, "w", encoding="utf-8") as wrapperFile:
        wrapperFile.write(wrapper_html)

def _extract_step_name(sourceRelPath: str, proteinName: str) -> str:
    """Extract the step directory that owns a report file."""

    normalizedParts = p.normpath(sourceRelPath).split(os.sep)

    if len(normalizedParts) >= 4 and normalizedParts[0] == "00_ccs_analysis" and normalizedParts[1] == proteinName:
        return normalizedParts[2]

    if normalizedParts and normalizedParts[0] == proteinName and len(normalizedParts) >= 2:
        return normalizedParts[1]

    if proteinName in normalizedParts:
        proteinIndex = normalizedParts.index(proteinName)
        if proteinIndex + 1 < len(normalizedParts):
            return normalizedParts[proteinIndex + 1]

    return ""


def _humanize_report_title(fileName: str) -> str:
    """Convert a report filename into a tab label."""

    lowerName = fileName.lower()
    if lowerName == "ccs_report.html":
        return "CCS Report"
    if lowerName == "vitals_report.html":
        return "Vitals Report"
    if lowerName == "rmsd_report.html":
        return "RMSD Report"

    return p.splitext(fileName)[0].replace("_", " ").strip().title()


def _step_sort_key(stepName: str) -> tuple[int, str]:
    """Sort step directories by their numeric prefix when present."""

    prefix = stepName.split("_", 1)[0]
    if prefix.isdigit():
        return (int(prefix), stepName)
    return (10**9, stepName)