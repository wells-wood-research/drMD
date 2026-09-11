from __future__ import annotations

import argparse
from pathlib import Path
from html import escape
import yaml
from ExaminationRoom.drRMSD import main as run_rmsd_report
from ExaminationRoom.drCollision import calculate_ccs_shadow as run_collision_report

def create_combined_report(operationDir: Path, batchConfig: dict) -> Path:
    
    inputDir= batchConfig["pathInfo"]["inputDir"]
    pdbNames = [p.stem for p in Path(inputDir).glob("*.pdb")]
    for pdbName in pdbNames:
        simulation_dir = Path(operationDir) / pdbName / "00_combined_simulation"
        print(simulation_dir)
        required_files = (
            simulation_dir / "trajectory.pdb",
            simulation_dir / "trajectory.dcd",
            simulation_dir / "vitals_report.csv",
        )


        # Reuse the existing RMSD/RMSF analysis and report generator.
        run_rmsd_report(str(simulation_dir))

        rmsd_report = simulation_dir / "rmsd_report.html"
        vitals_report = simulation_dir / "vitals_report.csv"



    batchConfig['simulationInfo'][-1]['stepName'] = "00_combined_simulation"
    batchConfig['aftercareInfo']['ccsReporterInfo']['steps'] = ["00_combined_simulation"]
    batchConfig['pathInfo']['outputDir'] = str(operationDir)
    run_collision_report(batchConfig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a report for combined simulation output."
    )
    parser.add_argument("simulation_dir", type=Path)
    parser.add_argument("--config", type=Path)
    args = parser.parse_args()

    batchConfig = None
    if args.config:
        with open(args.config, "r") as f:
            batchConfig: dict = yaml.safe_load(f)

    create_combined_report(args.simulation_dir, batchConfig)


if __name__ == "__main__":
    main()