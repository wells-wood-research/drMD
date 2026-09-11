from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Tuple
from os import path as p
import MDAnalysis as mda
import numpy as np
import plotly.graph_objects as go
from jinja2 import Environment, FileSystemLoader
from MDAnalysis.analysis import align, rms
from UtilitiesCloset import drMethodsWriter
from UtilitiesCloset.drCustomClasses import FilePath


DEFAULT_RMSD_SELECTION = "name CA"
DEFAULT_RMSF_SELECTION = "protein and name CA"
DEFAULT_STRIDE = 1
REPORT_PLOT_HEIGHT = 440
REPORT_PLOT_WIDTH = "100%"


def _read_text_if_file(file_path: Path) -> str | None:
    if not file_path.is_file():
        return None
    return file_path.read_text(encoding="utf-8", errors="replace")


def _wrap_plotly_html_fragment(plot_html: str) -> str:
    """Wrap a Plotly div/script fragment in a marginless full-size HTML page."""
    return (
        "<!DOCTYPE html>\n"
        "<html><head><meta charset=\"utf-8\">"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
        "<style>"
        f"html,body{{margin:0;padding:0;width:{REPORT_PLOT_WIDTH};height:{REPORT_PLOT_HEIGHT}px;min-height:{REPORT_PLOT_HEIGHT}px;overflow:hidden;background:#000000;}}"
        f"#plot-root{{width:{REPORT_PLOT_WIDTH};height:{REPORT_PLOT_HEIGHT}px;min-height:{REPORT_PLOT_HEIGHT}px;}}"
        f"#plot-root .plotly-graph-div{{width:{REPORT_PLOT_WIDTH}!important;height:{REPORT_PLOT_HEIGHT}px!important;min-height:{REPORT_PLOT_HEIGHT}px!important;}}"
        "</style></head><body><div id=\"plot-root\">"
        f"{plot_html}"
        "</div></body></html>"
    )


def _copy_shared_css(output_dir: Path, css_files: list[str]) -> None:
    """Copy shared CSS assets next to generated reports so they can be loaded by URL."""
    output_dir.mkdir(parents=True, exist_ok=True)
    base_dir = Path(__file__).resolve().parent
    for css_name in css_files:
        source_path = base_dir / css_name
        target_path = output_dir / css_name
        if source_path.is_file():
            shutil.copy2(source_path, target_path)


def calculate_rmsd(
    universe: mda.Universe,
    reference: mda.Universe,
    selection: str = DEFAULT_RMSD_SELECTION,
    ref_frame: int = 0,
) -> np.ndarray:
    rmsd_runner = rms.RMSD(universe, reference, select=selection, ref_frame=ref_frame)
    rmsd_runner.run()

    results = np.asarray(rmsd_runner.results.rmsd)
    if results.ndim == 1:
        results = results[None, :]

    if results.shape[1] >= 3:
        frame_numbers = results[:, 0]
        times = results[:, 1]
        rmsd_values = results[:, 2]
    else:
        frame_numbers = np.arange(results.shape[0], dtype=int)
        times = np.zeros(results.shape[0], dtype=float)
        rmsd_values = results[:, 0]

    return np.column_stack([frame_numbers, times, rmsd_values])


def calculate_rmsf(
    pdb_file: Path,
    universe: mda.Universe,
    selection: str = DEFAULT_RMSF_SELECTION,
    stride: int = DEFAULT_STRIDE,
) -> Tuple[np.ndarray, np.ndarray]:
    sampled_positions = np.array([ts.positions.copy() for ts in universe.trajectory[:: max(1, stride)]])
    if sampled_positions.size == 0:
        raise ValueError(f"No frames were available for RMSF calculation in {pdb_file}")

    sampled_universe = mda.Universe(str(pdb_file))
    sampled_universe.load_new(sampled_positions)

    average = align.AverageStructure(sampled_universe, sampled_universe, select=selection, ref_frame=0).run()
    reference = average.results.universe
    align.AlignTraj(sampled_universe, reference, select=selection, in_memory=True).run()

    selected_atoms = sampled_universe.select_atoms(selection)
    if len(selected_atoms) == 0:
        raise ValueError(f"Selection '{selection}' returned no atoms for {pdb_file}")

    rmsf_runner = rms.RMSF(selected_atoms).run()
    return selected_atoms.resids.copy(), np.asarray(rmsf_runner.results.rmsf)


def write_output_files(
    rmsd_path: Path,
    rmsf_path: Path,
    rmsd_initial_data: np.ndarray,
    rmsd_average_data: np.ndarray,
    rmsf_resids: np.ndarray,
    rmsf_values: np.ndarray,
) -> None:
    rmsd_row_count = max(rmsd_initial_data.shape[0], rmsd_average_data.shape[0])
    rmsd_combined = np.full((rmsd_row_count, 4), np.nan, dtype=float)

    if rmsd_initial_data.size:
        initial_rows = rmsd_initial_data.shape[0]
        rmsd_combined[:initial_rows, 0] = rmsd_initial_data[:, 0]
        rmsd_combined[:initial_rows, 1] = rmsd_initial_data[:, 1]
        rmsd_combined[:initial_rows, 2] = rmsd_initial_data[:, 2]

    if rmsd_average_data.size:
        average_rows = rmsd_average_data.shape[0]
        rmsd_combined[:average_rows, 3] = rmsd_average_data[:, 2]

    np.savetxt(
        rmsd_path.with_suffix(".csv"),
        rmsd_combined,
        header="Frame,Time(ps),RMSD initial(A),RMSD average(A)",
        fmt=["%.0f", "%.3f", "%.4f", "%.4f"],
        delimiter=",",
    )

    np.savetxt(
        rmsf_path.with_suffix(".csv"),
        np.column_stack([rmsf_resids, rmsf_values]),
        header="ResidueID, RMSF(A)",
        fmt=["%d", "%.6f"],
        delimiter=",",
    )


def build_rmsd_figure(rmsd_data: np.ndarray, title: str, line_color: str) -> go.Figure:
    fig = go.Figure()
    if rmsd_data.size:
        fig.add_trace(
            go.Scatter(
                x=rmsd_data[:, 1],
                y=rmsd_data[:, 2],
                mode="lines",
                name="RMSD",
                line=dict(color=line_color, width=2),
                hovertemplate="Time: %{x:.3f} ps<br>RMSD: %{y:.4f} Å<extra></extra>",
            )
        )
    fig.update_layout(
        template="plotly_dark",
        title=title,
        xaxis_title="Time (ps)",
        yaxis_title="RMSD (Å)",
        height=REPORT_PLOT_HEIGHT,
        margin=dict(l=60, r=30, t=60, b=60),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    if not rmsd_data.size:
        fig.add_annotation(
            text="No RMSD data available.",
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            showarrow=False,
        )
    return fig


def build_rmsf_figure(rmsf_resids: np.ndarray, rmsf_values: np.ndarray) -> go.Figure:
    fig = go.Figure()
    if rmsf_values.size:
        fig.add_trace(
            go.Scatter(
                x=rmsf_resids,
                y=rmsf_values,
                mode="lines+markers",
                name="RMSF",
                line=dict(color="#ffe45e", width=2),
                marker=dict(size=6),
                hovertemplate="Residue %{x}<br>RMSF: %{y:.4f} Å<extra></extra>",
            )
        )
    fig.update_layout(
        template="plotly_dark",
        title="Residue RMSF",
        xaxis_title="Residue ID",
        yaxis_title="RMSF (Å)",
        height=REPORT_PLOT_HEIGHT,
        margin=dict(l=60, r=30, t=60, b=60),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    if not rmsf_values.size:
        fig.add_annotation(
            text="No RMSF data available.",
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            showarrow=False,
        )
    return fig


def render_rmsd_report(
    report_path: Path,
    sim_dir: Path,
    topology_path: Path,
    trajectory_path: Path,
    rmsd_selection: str,
    rmsf_selection: str,
    stride: int,
    rmsd_initial_data: np.ndarray,
    rmsd_average_data: np.ndarray,
    rmsf_resids: np.ndarray,
    rmsf_values: np.ndarray,
    rmsd_path: Path,
    rmsf_path: Path,
) -> None:
    template_dir = Path(__file__).resolve().parent
    env = Environment(loader=FileSystemLoader(str(template_dir)), autoescape=True)
    template = env.get_template("rmsd_template.html")

    rmsd_initial_plot = build_rmsd_figure(
        rmsd_initial_data,
        "RMSD to Initial Structure",
        "#00ff87",
    ).to_html(
        full_html=False,
        include_plotlyjs="cdn",
        config={"responsive": True, "displaylogo": False, "scrollZoom": True},
    )
    rmsd_average_plot = build_rmsd_figure(
        rmsd_average_data,
        "RMSD to Average Structure",
        "#7ee3ff",
    ).to_html(
        full_html=False,
        include_plotlyjs="cdn",
        config={"responsive": True, "displaylogo": False, "scrollZoom": True},
    )
    rmsf_plot = build_rmsf_figure(rmsf_resids, rmsf_values).to_html(
        full_html=False,
        include_plotlyjs="cdn",
        config={"responsive": True, "displaylogo": False, "scrollZoom": True},
    )

    rmsd_initial_plot = _wrap_plotly_html_fragment(rmsd_initial_plot)
    rmsd_average_plot = _wrap_plotly_html_fragment(rmsd_average_plot)
    rmsf_plot = _wrap_plotly_html_fragment(rmsf_plot)

    rmsd_csv = rmsd_path.with_suffix(".csv")
    rmsf_csv = rmsf_path.with_suffix(".csv")

    plot_dir = report_path.parent / "00_plotly_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    _copy_shared_css(report_path.parent, ["report_theme.css"])
    rmsd_initial_plot_path = plot_dir / "rmsd_initial_plot.html"
    rmsd_average_plot_path = plot_dir / "rmsd_average_plot.html"
    rmsf_plot_path = plot_dir / "rmsf_plot.html"

    rmsd_initial_plot_path.write_text(rmsd_initial_plot, encoding="utf-8")
    rmsd_average_plot_path.write_text(rmsd_average_plot, encoding="utf-8")
    rmsf_plot_path.write_text(rmsf_plot, encoding="utf-8")

    embedded_csv_files = []
    for csv_label, csv_path in (
        ("RMSD", rmsd_csv),
        ("RMSF", rmsf_csv),
    ):
        if csv_path.is_file():
            embedded_csv_files.append(
                {
                    "label": csv_label,
                    "fileName": csv_path.name,
                    "href": p.relpath(csv_path, report_path.parent),
                }
            )

    context = {
        "rmsd_data": {
            "systemName": sim_dir.parent.name,
            "stepName": sim_dir.name,
            "nFrames": int(rmsd_initial_data.shape[0]),
            "reportPath": report_path.name,
            "description": "Interactive RMSD and RMSF analysis for a single trajectory.",
            "inputRoot": str(sim_dir),
            "analysisStep": "RMSD / RMSF",
            "reportTitle": f"RMSD Report - {sim_dir.name}",
            "topologyPdb": str(topology_path),
            "trajectoryDcd": str(trajectory_path),
            "outputDir": str(report_path.parent),
            "rmsdSelection": rmsd_selection,
            "rmsfSelection": rmsf_selection,
            "stride": stride,
            "frameCount": int(rmsd_initial_data.shape[0]),
            "rmsdInitialMean": f"{float(np.nanmean(rmsd_initial_data[:, 2])):.3f}" if rmsd_initial_data.size else "nan",
            "rmsdAverageMean": f"{float(np.nanmean(rmsd_average_data[:, 2])):.3f}" if rmsd_average_data.size else "nan",
            "rmsfMax": f"{float(np.nanmax(rmsf_values)):.3f}" if rmsf_values.size else "nan",
            "rmsdReport": rmsd_path.name,
            "rmsfReport": rmsf_path.name,
            "rmsdInitialPlot": p.relpath(rmsd_initial_plot_path, report_path.parent),
            "rmsdAveragePlot": p.relpath(rmsd_average_plot_path, report_path.parent),
            "rmsfPlot": p.relpath(rmsf_plot_path, report_path.parent),
            "embeddedCsvFiles": embedded_csv_files,
            "trajectoryRows": [
                {
                    "run": sim_dir.name,
                    "frames": int(rmsd_initial_data.shape[0]),
                    "rmsdInitialMean": f"{float(np.nanmean(rmsd_initial_data[:, 2])):.3f}" if rmsd_initial_data.size else "nan",
                    "rmsdAverageMean": f"{float(np.nanmean(rmsd_average_data[:, 2])):.3f}" if rmsd_average_data.size else "nan",
                    "rmsfPoints": int(rmsf_values.size),
                    "rmsdFile": rmsd_path.name,
                    "rmsfFile": rmsf_path.name,
                }
            ],
        }
    }

    report_path.write_text(template.render(**context), encoding="utf-8")


def main(simDir: FilePath) -> None:
    output_dir = simDir
    stride = DEFAULT_STRIDE
    rmsd_selection = DEFAULT_RMSD_SELECTION
    rmsf_selection = DEFAULT_RMSF_SELECTION
    drMethodsWriter.add_analysis_step_to_log("RMSD")
    drMethodsWriter.add_parameter_to_analysis_log("RMSD", "stepName", "all")
    drMethodsWriter.add_parameter_to_analysis_log("RMSD", "selectionRMSD", rmsd_selection)
    drMethodsWriter.add_parameter_to_analysis_log("RMSD", "selectionRMSF", rmsf_selection)
    drMethodsWriter.add_parameter_to_analysis_log("RMSD", "stride", stride)

    sim_dir = Path(simDir).expanduser().resolve()
    if not sim_dir.is_dir():
        raise NotADirectoryError(f"Simulation directory does not exist: {sim_dir}")

    trajectory_path = sim_dir / "trajectory.dcd"
    topology_path = sim_dir / "trajectory.pdb"
    if not trajectory_path.is_file():
        raise FileNotFoundError(f"Trajectory file does not exist: {trajectory_path}")
    if not topology_path.is_file():
        raise FileNotFoundError(f"Topology file does not exist: {topology_path}")

    output_path = Path(output_dir).expanduser().resolve() if output_dir else sim_dir
    output_path.mkdir(parents=True, exist_ok=True)

    rmsd_path = output_path / "rmsd.dat"
    rmsf_path = output_path / "rmsf.dat"
    report_path = output_path / "rmsd_report.html"

    universe = mda.Universe(str(topology_path), str(trajectory_path))

    reference_initial = mda.Universe(str(topology_path), str(trajectory_path))
    rmsd_initial_data = calculate_rmsd(universe, reference_initial, rmsd_selection, ref_frame=0)

    average_structure = align.AverageStructure(universe, universe, select=rmsd_selection, ref_frame=0).run()
    reference_average = average_structure.results.universe
    rmsd_average_data = calculate_rmsd(universe, reference_average, rmsd_selection, ref_frame=0)

    rmsf_resids, rmsf_values = calculate_rmsf(topology_path, universe, rmsf_selection, stride)

    write_output_files(
        rmsd_path,
        rmsf_path,
        rmsd_initial_data,
        rmsd_average_data,
        rmsf_resids,
        rmsf_values,
    )
    render_rmsd_report(
        report_path=report_path,
        sim_dir=sim_dir,
        topology_path=topology_path,
        trajectory_path=trajectory_path,
        rmsd_selection=rmsd_selection,
        rmsf_selection=rmsf_selection,
        stride=stride,
        rmsd_initial_data=rmsd_initial_data,
        rmsd_average_data=rmsd_average_data,
        rmsf_resids=rmsf_resids,
        rmsf_values=rmsf_values,
        rmsd_path=rmsd_path,
        rmsf_path=rmsf_path,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate RMSD and RMSF report for a simulation directory.")
    parser.add_argument("simDir", help="Path to the simulation directory")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(simDir=args.simDir)
