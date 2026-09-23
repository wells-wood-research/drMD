"""CCS report rendering helpers for ExaminationRoom."""

import base64
import mimetypes
import os
import shutil
from os import path as p
from typing import Optional

import pandas as pd
from jinja2 import Environment, FileSystemLoader
import numpy as np
from ExaminationRoom import drLogger
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath
import plotly.express as px
import plotly.figure_factory as ff
from sklearn.cluster import KMeans

REPORT_PLOT_HEIGHT = 780
REPORT_PLOT_WIDTH = "100%"


def _copy_shared_css(output_dir: str, css_files: list[str]) -> None:
    """Copy shared CSS assets next to generated reports so they can be loaded by URL."""
    output_path = p.abspath(output_dir)
    os.makedirs(output_path, exist_ok=True)
    base_dir = p.dirname(__file__)
    for css_name in css_files:
        source_path = p.join(base_dir, css_name)
        target_path = p.join(output_path, css_name)
        if p.isfile(source_path):
            shutil.copy2(source_path, target_path)


def _wrap_plotly_html_fragment(plot_html: str) -> str:
    """Wrap a Plotly fragment and bridge hover events to the parent report."""
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
        "</div><script>"
        "(function(){"
        "function flattenCustomData(value){"
        "var values=[];"
        "function visit(item){"
        "if(item===null||typeof item==='undefined'){return;}"
        "if(Array.isArray(item)){item.forEach(visit);return;}"
        "if(typeof item==='object'){Object.keys(item).forEach(function(key){visit(item[key]);});return;}"
        "values.push(String(item));"
        "}"
        "visit(value);"
        "return values;"
        "}"
        "function frameTokenFromPoint(point){"
        "if(!point){return null;}"
        "var customData=point.customdata;"
        "var candidates=flattenCustomData(customData);"
        "if(typeof point.pointNumber!=='undefined'){candidates.push(String(point.pointNumber));}"
        "if(typeof point.id!=='undefined'){candidates.push(String(point.id));}"
        "for(var i=0;i<candidates.length;i++){var candidate=candidates[i];if(!candidate){continue;}if(candidate.indexOf('cluster_')===0||candidate.indexOf('frame_')===0||/^\\d+$/.test(candidate)){return candidate;}}"
        "return null;"
        "}"
        "function postFrameToken(point){"
        "var frameToken=frameTokenFromPoint(point);"
        "if(!frameToken||!window.parent||window.parent===window){return;}"
        "window.parent.postMessage({type:'ccs-frame-hover',frameToken:frameToken},'*');"
        "}"
        "function attachPlotEvents(){"
        "var plotDivs=document.querySelectorAll('.plotly-graph-div');"
        "if(!plotDivs.length){return;}"
        "plotDivs.forEach(function(plotDiv){"
        "if(typeof plotDiv.on!=='function'){return;}"
        "plotDiv.on('plotly_hover',function(eventData){"
        "if(eventData&&eventData.points&&eventData.points[0]){postFrameToken(eventData.points[0]);}"
        "});"
        "plotDiv.on('plotly_click',function(eventData){"
        "if(eventData&&eventData.points&&eventData.points[0]){postFrameToken(eventData.points[0]);}"
        "});"
        "});"
        "}"
        "if(document.readyState==='loading'){document.addEventListener('DOMContentLoaded',attachPlotEvents);}"
        "else{attachPlotEvents();}"
        "})();"
        "</script></body></html>"
    )


def interactive_ccs_over_time(df, logTimeUnit, stride):

    # Remove any leading/trailing whitespace from column names
    df = df.copy()
    df = df.reset_index(drop=True)
    df.columns = df.columns.str.strip()

    df["Frame ID"] = df.index
    df["Frame Token"] = df["Frame ID"].astype(str)
    # Human-friendly filename for each frame's binary pixel map. Must match
    # the actual on-disk filenames built in create_ccs_report_html, which are
    # named by original (pre-stride) frame number, i.e. row index * stride --
    # not the raw row index -- since only every `stride`-th frame gets saved.
    df["Frame File"] = df["Frame ID"].apply(lambda i: f"frame_{int(i) * (stride or 1)}.bin")
    # Calculate SASA / Convex Hull ratio



    # Row-position identifier, used by the report's JS to look up each
    # point's binary pixel map image. Must stay in sync with the "Frame ID"
    # keys used to build ccsPixelMapMap in create_ccs_report_html.

    # Interactive scatter plot
    fig = px.scatter(
    df,
    x=f"Time ({logTimeUnit})",
    y="CCS (Å²)",
    title="CCS over Time",
    custom_data=["Frame Token"],
    labels={
        f"Time ({logTimeUnit})": f"Time ({logTimeUnit})",
        "CCS (Å²)": "CCS (Å²)",
    },
    hover_data={
            # Placed first so it lands at customdata[0] for the pixel-map viewer.
            "Frame File": True,
            "Frame ID": True,
            f"Time ({logTimeUnit})": ":.2f",
            "CCS (Å²)": ":.2f",
        },
    )   

    # Connect points in time
    fig.update_traces(
    mode="lines+markers",
    marker=dict(size=6),
    customdata=df["Frame Token"].to_numpy(),
    )

    fig.update_layout(
    template="plotly_dark",
    hovermode="x unified",
    height=REPORT_PLOT_HEIGHT,
    )
    return fig


def create_kmeans_plot(df, kmeans, logTimeUnit, stride):

    df = df.copy()
    df = df.reset_index(drop=True)
    df.columns = df.columns.str.strip()

    df["Frame ID"] = df.index
    df["Frame Token"] = df["Frame ID"].astype(str)
    # Human-friendly filename for each frame's binary pixel map
    df["Frame File"] = df["Frame ID"].apply(lambda i: f"frame_{int(i)}.bin")
    # Cluster as a string so Plotly treats it as a discrete/qualitative
    # variable (one color swatch per cluster) rather than a continuous one.
    # Shift cluster labels so they start at 1 instead of 0 for user-facing
    # labels and filenames (cluster_1.gif, cluster_2.gif, ...).
    df["Cluster"] = (kmeans.labels_.astype(int) + 1).astype(str)

    # Cluster labels are 0..k-1, matching the row order of kmeans.cluster_centers_
    clusterOrder = sorted(df["Cluster"].unique(), key=int)
    palette = px.colors.qualitative.G10
    colorMap = {cluster: palette[i % len(palette)] for i, cluster in enumerate(clusterOrder)}

    fig = px.scatter(
        df,
        x="CCS (Å²)",
        y="Radius of Gyration (Å²)",
        color="Cluster",
        category_orders={"Cluster": clusterOrder},
        color_discrete_map=colorMap,
        title="KMeans Clustering",
        custom_data=["Frame Token"],
        labels={
            "CCS (Å²)": "CCS (Å²)",
            "Radius of Gyration (Å²)": "Radius of Gyration (Å²)",
            "Cluster": "Cluster",
        },
        hover_data={
            # Placed first so it lands at customdata[0] for the pixel-map viewer.
            "Frame File": True,
            "Frame ID": True,
            f"Time ({logTimeUnit})": True,
            "CCS (Å²)": ":.2f",
        },
    )

    fig.update_traces(marker=dict(size=7), customdata=df["Frame Token"].to_numpy())

    # Plot each cluster centroid in the same color as its cluster's points.
    # Assumes kmeans was fit on ["CCS (Å²)", "Convex Hull Area (Å²)"], in
    # that order, so cluster_centers_ columns line up with the x/y axes above.
    centers = kmeans.cluster_centers_
    for clusterIdx, (centerX, centerY) in enumerate(centers):
        # `clusterIdx` is 0-based here; user-visible cluster labels are
        # 1-based, so increment when constructing labels/filenames.
        clusterLabel = str(clusterIdx + 1)
        if clusterLabel not in colorMap:
            continue
        fig.add_scatter(
            x=[centerX],
            y=[centerY],
            mode="markers",
            marker=dict(
                symbol="x",
                size=16,
                color=colorMap[clusterLabel],
                line=dict(width=2, color="white"),
            ),
            name=f"Centroid {clusterLabel}",
            showlegend=False,
            hovertemplate=(
                f"Cluster {clusterLabel} centroid<br>"
                "CCS: %{x:.2f}<br>Radius of Gyration: %{y:.2f}<extra></extra>"
            ),
            # customdata uses the mapping key used in `ccsPixelMapMap`.
            customdata=[[f"cluster_{clusterIdx+1}" ]],
        )

    # Dummy trace so the legend explains what the "x" markers mean, since the
    # per-cluster centroid traces above are hidden from the legend to avoid
    # duplicating each cluster's entry.
    fig.add_scatter(
        x=[None],
        y=[None],
        mode="markers",
        marker=dict(symbol="x", size=12, color="white", line=dict(width=2, color="white")),
        name="Cluster Centroid",
        showlegend=True,
    )

    fig.update_layout(
        template="plotly_dark",
        height=REPORT_PLOT_HEIGHT,
        legend=dict(
            title="Cluster",
            orientation="h",
            yanchor="top",
            y=-0.25,
            xanchor="center",
            x=0.5,
        ),
        margin=dict(b=110),
    )

    return fig


def create_kde_plot(df, logTimeUnit):

    df = df.copy()
    df.columns = df.columns.str.strip()

    ccsValues = df["CCS (Å²)"].dropna()

    # create_distplot needs at least a couple of distinct points to fit a KDE
    if len(ccsValues) < 2 or ccsValues.nunique() < 2:
        fig = px.scatter(title="CCS Density (insufficient data for KDE)")
        fig.update_layout(template="plotly_dark", height=REPORT_PLOT_HEIGHT)
        return fig

    fig = ff.create_distplot(
        [ccsValues.tolist()],
        group_labels=["CCS (Å²)"],
        colors=["#00cc96"],
        show_hist=False,
        show_rug=True,
    )

    fig.update_layout(
        template="plotly_dark",
        title="CCS Kernel Density Estimate",
        xaxis_title="CCS (Å²)",
        yaxis_title="Density",
        height=REPORT_PLOT_HEIGHT,
        showlegend=False,
    )

    return fig


def create_ccs_report_html(
    ccsDir: DirectoryPath,
    rawCsv: FilePath,
    kmeansData,
    gas: Optional[str] =None,
    logTimeUnit: Optional[str] = None,
    stride: Optional[int] = None,
) -> Optional[FilePath]:
    """Render the CCS HTML report from the CCS report template."""
    if not p.exists(rawCsv):
        drLogger.log_info(f"CCS report CSV missing: {rawCsv}. Skipping HTML report generation.", True)
        return None
    ccsDf = pd.read_csv(rawCsv)
    plotDir = p.join(ccsDir, "00_plotly_plots")
    os.makedirs(plotDir, exist_ok=True)
    ccsPlot = interactive_ccs_over_time(ccsDf, logTimeUnit, stride)
    ccsPlotPath = p.join(plotDir, "ccs_over_time_plot.html")
    ccsPlotHtml = ccsPlot.to_html(
        full_html=False,
        include_plotlyjs="cdn",
        config={
            "responsive": True,
            "displaylogo": False,
            "scrollZoom": True,
        },
    )
    ccsPlotHtml = _wrap_plotly_html_fragment(ccsPlotHtml)
    with open(ccsPlotPath, "w", encoding="utf-8") as plot_file:
        plot_file.write(ccsPlotHtml)

    kmeansPlot = create_kmeans_plot(ccsDf, kmeansData, logTimeUnit, stride)
    kmeansPlotPath = p.join(plotDir, "kmeans_plot.html")
    kmeansPlotHtml = kmeansPlot.to_html(
        full_html=False,
        include_plotlyjs="cdn",
        config={
            "responsive": True,
            "displaylogo": False,
            "scrollZoom": True,
        },
    )
    kmeansPlotHtml = _wrap_plotly_html_fragment(kmeansPlotHtml)
    with open(kmeansPlotPath, "w", encoding="utf-8") as plot_file:
        plot_file.write(kmeansPlotHtml)

    kdePlot = create_kde_plot(ccsDf, logTimeUnit)
    kdePlotPath = p.join(plotDir, "kde_plot.html")
    kdePlotHtml = kdePlot.to_html(
        full_html=False,
        include_plotlyjs="cdn",
        config={
            "responsive": True,
            "displaylogo": False,
            "scrollZoom": True,
        },
    )
    kdePlotHtml = _wrap_plotly_html_fragment(kdePlotHtml)
    with open(kdePlotPath, "w", encoding="utf-8") as plot_file:
        plot_file.write(kdePlotHtml)
    ccsHeaders = [str(col) for col in ccsDf.columns.tolist()]
    ccsData = [row.tolist() for _, row in ccsDf.iterrows()]

    # Per-frame binary pixel map images for the report's hover-updating side
    # panel. ASSUMPTION: these live in a "frames" subfolder next to the raw
    # CCS CSV, named "frame_<row index>.bin" (an 8-byte header of two
    # little-endian uint32s — width, then height — followed by a packed
    # 1-bit-per-pixel bitmap; see the reader in ccs_report_template.html) —
    # adjust pixelMapDir / the filename pattern below if your pipeline
    # writes them somewhere else.
    #
    # Files are base64-embedded directly into the rendered HTML (rather than
    # linked by relative path) so the report is a single self-contained file
    # that still works if it's moved or shared outside ccsDir. This does
    # inflate the HTML file size by roughly 4/3 the total size of all the
    # embedded .bin files, so it's worth keeping an eye on nFrames for very
    # long trajectories.
    pixelMapDir = p.join(p.dirname(rawCsv), "frames")
    ccsPixelMapMap = {}

    # Build a map: frame index (string) -> {file: filename, href: relative path, data: embedded base64}
    # The report may be opened directly from disk (file://), where fetch() to local .bin files is blocked by the browser,
    # so we embed the actual binary payload in the HTML as base64 for the JS viewer to decode immediately.
    for frameIdx in range(len(ccsDf)):
        fname = f"frame_{int(frameIdx) * (stride or 1)}.bin"
        candidate = p.join(pixelMapDir, fname)
        if p.exists(candidate) and p.isfile(candidate):
            with open(candidate, "rb") as frame_file:
                payload = frame_file.read()
            ccsPixelMapMap[str(frameIdx)] = {
                "file": fname,
                "href": p.relpath(candidate, ccsDir),
                "data": base64.b64encode(payload).decode("ascii"),
            }

    # Also include any per-cluster GIFs named cluster_<i>.gif so centroid
    # markers in the KMeans plot can display/download them. Keys use
    # "cluster_<i>" so they don't clash with numeric frame keys.
    try:
        n_clusters = None
        if hasattr(kmeansData, "n_clusters") and kmeansData.n_clusters is not None:
            n_clusters = int(kmeansData.n_clusters)
        elif hasattr(kmeansData, "labels_") and kmeansData.labels_ is not None:
            n_clusters = len(set(kmeansData.labels_))
        if n_clusters is not None:
            for clusterIdx in range(n_clusters):
                # Files are named cluster_1.gif ... cluster_n.gif
                one_based = clusterIdx + 1
                fname = f"00_clustered_pdbs/cluster_{int(one_based)}.gif"
                # Try the clustered pdbs directory next to the CSV first
                candidate = p.join(p.dirname(rawCsv), fname)
                # Fallback: also accept it under the frames/ subfolder
                if not (p.exists(candidate) and p.isfile(candidate)):
                    candidate = p.join(pixelMapDir, fname)
                b64 = None
                if p.exists(candidate) and p.isfile(candidate):
                    ccsPixelMapMap[f"cluster_{one_based}"] = {
                        "file": fname,
                        "href": p.relpath(candidate, ccsDir),
                    }
    except Exception:
        # Non-fatal: if kmeansData isn't present or is unexpected, continue.
        pass

    if not ccsPixelMapMap:
        drLogger.log_info(
        f"No CCS binary pixel map files found under {pixelMapDir}; "
        "the report's pixel map panel will stay empty until they're generated there.",
        True,
        )
    
    templateDir = p.dirname(__file__)
    env = Environment(loader=FileSystemLoader(templateDir), autoescape=True)
    template = env.get_template("ccs_report_template.html")
    _copy_shared_css(ccsDir, ["report_theme.css"])
    protein = ccsDir.split("/")[-2]
    stepName = ccsDir.split("/")[-1]
    context = {
        "ccs_data": {
            "protein": protein,
            "stepName": stepName,
            "gas": gas,
            "ccsMethod": "ShadowScreen",
            "timeUnit": logTimeUnit,
            "nFrames": len(ccsDf),
            "rawCsv": p.relpath(rawCsv, ccsDir),
            "rawCsvName": p.basename(rawCsv),
            "rawCsvContent": None,
            "ccsPlot": p.relpath(ccsPlotPath, ccsDir),
            "kmeansPlot": p.relpath(kmeansPlotPath, ccsDir),
            "kdePlot": p.relpath(kdePlotPath, ccsDir),
            "ccsPixelMapMap": ccsPixelMapMap,
            "colorBar": None,
            "ccsHeaders": ccsHeaders,
            "ccsData": ccsData,
        }
    }

    rendered_html = template.render(context)
    report_path = p.join(ccsDir, "ccs_report.html")
    with open(report_path, "w", encoding="utf-8") as report_file:
        report_file.write(rendered_html)

    drLogger.log_info(f"Created CCS report HTML: {report_path}", True)
    return report_path