## BASIC PYTHON LIBRARIES
import os
from os import path as p
from shutil import copy2

## CLEAN CODE
from typing import Dict, List
from jinja2 import Environment, FileSystemLoader


def build_operation_report_index(batch_output_dir: str) -> str | None:
    """Create a portable HTML index for the collated protein full reports."""

    if not p.isdir(batch_output_dir):
        raise NotADirectoryError(f"Batch output directory not found: {batch_output_dir}")

    template_dir = p.dirname(__file__)
    env = Environment(loader=FileSystemLoader(template_dir), autoescape=True)
    template = env.get_template("operation_report_index_template.html")

    portable_reports_dir = p.join(batch_output_dir, "00_portable_reports")
    os.makedirs(portable_reports_dir, exist_ok=True)
    _copy_css_assets(portable_reports_dir, ["operation_report_index.css"])

    protein_entries, operation_entries = _collect_report_entries(batch_output_dir)

    if not operation_entries and not protein_entries:
        return None

    output_path = p.join(portable_reports_dir, "full_report_index.html")
    with open(output_path, "w", encoding="utf-8") as index_file:
        index_file.write(
            template.render(
                batchDir=p.basename(batch_output_dir),
                operationEntries=operation_entries,
                proteinEntries=protein_entries,
            )
        )

    return output_path


def _copy_css_assets(target_dir: str, css_files: list[str]) -> None:
    """Copy standalone CSS files next to the generated HTML pages so relative links remain valid."""
    os.makedirs(target_dir, exist_ok=True)
    template_dir = p.dirname(__file__)
    for css_name in css_files:
        source_path = p.join(template_dir, css_name)
        target_path = p.join(target_dir, css_name)
        if p.isfile(source_path):
            copy2(source_path, target_path)


def _collect_report_entries(batch_output_dir: str) -> tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    """Collect per-protein report links from the batch output tree.

    The project stores each operation in its own directory, and each one contains a nested
    00_full_reports folder with full_report_index.html entries for the protein reports inside it.
    """

    protein_entries: List[Dict[str, str]] = []
    operation_entries: List[Dict[str, str]] = []
    portable_reports_dir = p.join(batch_output_dir, "00_portable_reports")

    for entry_name in sorted(os.listdir(batch_output_dir)):
        operation_dir = p.join(batch_output_dir, entry_name)
        if not p.isdir(operation_dir) or entry_name.startswith("00_"):
            continue

        full_reports_root = p.join(operation_dir, "00_full_reports")
        if not p.isdir(full_reports_root):
            continue

        operation_index = p.join(full_reports_root, "full_report_index.html")
        if not p.isfile(operation_index):
            continue

        report_rel_path = p.relpath(operation_index, portable_reports_dir)
        operation_entries.append(
            {
                "operationName": entry_name,
                "reportCount": str(_count_reports_in_operation(full_reports_root)),
                "reportRelPath": report_rel_path,
            }
        )

        for protein_name in sorted(os.listdir(full_reports_root)):
            protein_dir = p.join(full_reports_root, protein_name)
            if not p.isdir(protein_dir):
                continue

            report_path = p.join(protein_dir, "full_report.html")
            if not p.isfile(report_path):
                continue

            protein_entries.append(
                {
                    "operationName": entry_name,
                    "proteinName": protein_name,
                    "reportRelPath": p.relpath(report_path, portable_reports_dir),
                }
            )

    return protein_entries, operation_entries


def _count_reports_in_operation(full_reports_root: str) -> int:
    """Count HTML full-report documents stored under an operation's 00_full_reports tree."""

    report_count = 0
    for root, dirs, files in os.walk(full_reports_root):
        dirs[:] = [d for d in dirs if not d.startswith("00_")]
        for file_name in files:
            if file_name.lower().endswith(".html"):
                report_count += 1
    return report_count


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build a batch index of operation full reports.")
    parser.add_argument("batch_output_dir", help="Batch output directory containing operation folders")
    args = parser.parse_args()

    result = build_operation_report_index(args.batch_output_dir)
    if result:
        print(f"Created batch full report index: {result}")
    else:
        print("No operation full reports were found to index.")