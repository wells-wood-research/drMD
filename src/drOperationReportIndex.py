## BASIC PYTHON LIBRARIES
import os
import re
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

    Each operation is stored in its own directory, but its 00_full_reports folder may be nested
    (e.g. ESI/ESI_1/00_full_reports or VacuumCM/CM_1/00_full_reports), so every 00_full_reports
    folder below an operation directory is collected. A per-operation full_report_index.html is
    not required; any <protein>/full_report.html found is indexed.
    """

    protein_entries: List[Dict[str, str]] = []
    operation_entries: List[Dict[str, str]] = []
    portable_reports_dir = p.join(batch_output_dir, "00_portable_reports")

    for entry_name in sorted(os.listdir(batch_output_dir)):
        operation_dir = p.join(batch_output_dir, entry_name)
        if not p.isdir(operation_dir) or entry_name.startswith("00_"):
            continue

        for full_reports_root in _find_full_report_roots(operation_dir):
            run_rel_path = p.relpath(p.dirname(full_reports_root), batch_output_dir)
            operation_name = run_rel_path.replace(os.sep, " / ")

            operation_protein_entries: List[Dict[str, str]] = []
            for protein_name in sorted(os.listdir(full_reports_root)):
                report_path = p.join(full_reports_root, protein_name, "full_report.html")
                if not p.isfile(report_path):
                    continue

                operation_protein_entries.append(
                    {
                        "operationName": operation_name,
                        "proteinName": protein_name,
                        "reportRelPath": p.relpath(report_path, portable_reports_dir),
                    }
                )

            if not operation_protein_entries:
                continue

            operation_index = p.join(full_reports_root, "full_report_index.html")
            if p.isfile(operation_index):
                operation_rel_path = p.relpath(operation_index, portable_reports_dir)
            else:
                operation_rel_path = operation_protein_entries[0]["reportRelPath"]

            operation_entries.append(
                {
                    "operationName": operation_name,
                    "reportCount": str(_count_reports_in_operation(full_reports_root)),
                    "reportRelPath": operation_rel_path,
                }
            )
            protein_entries.extend(operation_protein_entries)

    return protein_entries, operation_entries


def _find_full_report_roots(operation_dir: str) -> List[str]:
    """Return every 00_full_reports directory beneath an operation directory, in sorted order."""

    full_report_roots: List[str] = []
    for root, dirs, _ in os.walk(operation_dir):
        if "00_full_reports" in dirs:
            full_report_roots.append(p.join(root, "00_full_reports"))
        # don't descend into other 00_* folders (configs, logs, collated pdbs, the reports themselves)
        dirs[:] = sorted(d for d in dirs if not d.startswith("00_"))
    return sorted(full_report_roots, key=_natural_sort_key)


def _natural_sort_key(path: str) -> list:
    """Sort ESI_2 before ESI_10."""
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", path)]


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