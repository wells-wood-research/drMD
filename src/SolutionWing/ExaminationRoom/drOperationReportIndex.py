## BASIC PYTHON LIBRARIES
import os
from os import path as p
from shutil import copy2

## CLEAN CODE
from typing import Dict, List
from jinja2 import Environment, FileSystemLoader


def build_operation_report_index(batch_output_dir: str) -> str | None:
    """Create an HTML index that links each operation's full report index."""

    if not p.isdir(batch_output_dir):
        raise NotADirectoryError(f"Batch output directory not found: {batch_output_dir}")

    template_dir = p.dirname(__file__)
    env = Environment(loader=FileSystemLoader(template_dir), autoescape=True)
    template = env.get_template("operation_report_index_template.html")

    _copy_css_assets(batch_output_dir, ["operation_report_index.css"])
    index_entries: List[Dict[str, str]] = []

    for entry_name in sorted(os.listdir(batch_output_dir)):
        if not _is_operation_dir(entry_name):
            continue

        operation_dir = p.join(batch_output_dir, entry_name)
        operation_index = p.join(operation_dir, "00_full_reports", "full_report_index.html")
        if not p.isfile(operation_index):
            continue

        index_entries.append(
            {
                "operationName": entry_name,
                "reportCount": _count_operation_reports(operation_dir),
                "reportRelPath": p.relpath(operation_index, batch_output_dir),
            }
        )

    if not index_entries:
        return None

    output_path = p.join(batch_output_dir, "full_report_index.html")
    with open(output_path, "w", encoding="utf-8") as index_file:
        index_file.write(
            template.render(
                batchDir=p.basename(batch_output_dir),
                indexEntries=index_entries,
            )
        )

    return output_path


def _copy_css_assets(target_dir: str, css_files: list[str]) -> None:
    """Copy standalone CSS files next to generated HTML output."""
    os.makedirs(target_dir, exist_ok=True)
    template_dir = p.dirname(__file__)
    for css_name in css_files:
        source_path = p.join(template_dir, css_name)
        target_path = p.join(target_dir, css_name)
        if p.isfile(source_path):
            copy2(source_path, target_path)


def _is_operation_dir(entry_name: str) -> bool:
    if entry_name.startswith("00_"):
        return False
    return True


def _count_operation_reports(operation_dir: str) -> int:
    full_reports_root = p.join(operation_dir, "00_full_reports")
    if not p.isdir(full_reports_root):
        return 0

    report_count = 0
    for root, dirs, files in os.walk(full_reports_root):
        dirs[:] = [d for d in dirs if not d.startswith("00_")]
        for file_name in files:
            if file_name.endswith(".html"):
                report_count += 1
    return report_count


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build a batch index of operation full reports.")
    parser.add_argument("batch_output_dir", help="Batch output directory containing operation folders")
    args = parser.parse_args()

    result = build_operation_report_index(args.batch_output_dir)
    if result:
        print(f"Created operation full report index: {result}")
    else:
        print("No operation full reports were found to index.")