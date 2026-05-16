"""
plotting.py — Visualisation Module for TrajMedic.

Generates publication-style Matplotlib figures and interactive
Plotly charts for RMSD, Rg, and COM drift data.

All Matplotlib figures are saved to the outputs/ directory and also
returned as Figure objects so Streamlit can display them inline.
"""

from __future__ import annotations

import os
import logging
from pathlib import Path
from typing import List, Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")   # headless backend — safe for Streamlit
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import plotly.graph_objects as go

logger = logging.getLogger(__name__)

# ── Global style ──────────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi":         150,
    "font.family":        "sans-serif",
    "font.size":          11,
    "axes.titlesize":     13,
    "axes.labelsize":     11,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "grid.alpha":         0.3,
    "lines.linewidth":    1.6,
    "legend.framealpha":  0.7,
    "savefig.bbox":       "tight",
    "savefig.facecolor":  "white",
})

PALETTE = {
    "blue":    "#2980b9",
    "green":   "#27ae60",
    "orange":  "#e67e22",
    "red":     "#c0392b",
    "purple":  "#8e44ad",
    "grey":    "#95a5a6",
}

OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)


# ──────────────────────────────────────────────────────────────────────────────
# RMSD Plot
# ──────────────────────────────────────────────────────────────────────────────

def plot_rmsd(
    times_ns: List[float],
    rmsd_values: List[float],
    jump_frame_indices: Optional[List[int]] = None,
    equilibration_time_ns: float = -1.0,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Backbone RMSD vs simulation time.

    Parameters
    ----------
    times_ns : list[float]
        Simulation time in ns for each data point.
    rmsd_values : list[float]
        RMSD in Ångströms for each frame.
    jump_frame_indices : list[int], optional
        Indices of frames flagged as sudden jumps.
    equilibration_time_ns : float
        If > 0, a vertical dashed line marks equilibration onset.
    save_path : str, optional
        Overrides the default output path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, ax = plt.subplots(figsize=(8, 3.8))

    t = np.array(times_ns)
    r = np.array(rmsd_values)

    ax.plot(t, r, color=PALETTE["blue"], alpha=0.85, label="Backbone RMSD")

    # Running mean
    if len(r) >= 20:
        window = max(5, len(r) // 20)
        rm = np.convolve(r, np.ones(window) / window, mode="valid")
        ax.plot(t[window - 1:], rm, color=PALETTE["red"],
                linewidth=1.0, linestyle="--", alpha=0.9, label=f"Running mean ({window}-frame)")

    # Flag jump events
    if jump_frame_indices:
        jump_mask = np.zeros(len(t), dtype=bool)
        for ji in jump_frame_indices:
            # find nearest index in our sampled times array
            nearest = np.argmin(np.abs(np.arange(len(t)) - ji))
            if nearest < len(t):
                jump_mask[nearest] = True
        ax.scatter(t[jump_mask], r[jump_mask],
                   color=PALETTE["red"], s=40, zorder=5,
                   marker="x", label="Detected jump")

    # Equilibration marker
    if equilibration_time_ns > 0:
        ax.axvline(equilibration_time_ns, color=PALETTE["green"],
                   linestyle=":", linewidth=1.5,
                   label=f"Equilibration ≈ {equilibration_time_ns:.1f} ns")

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("Backbone RMSD vs Time")
    ax.legend(fontsize=9, loc="upper right")

    _finalize(fig, save_path or str(OUTPUT_DIR / "rmsd.png"))
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Radius of Gyration Plot
# ──────────────────────────────────────────────────────────────────────────────

def plot_rg(
    times_ns: List[float],
    rg_values: List[float],
    spike_frame_indices: Optional[List[int]] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Radius of gyration vs simulation time.

    Parameters
    ----------
    times_ns : list[float]
    rg_values : list[float]
        Rg in Ångströms.
    spike_frame_indices : list[int], optional
        Frames flagged as Rg outliers.
    save_path : str, optional

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, ax = plt.subplots(figsize=(8, 3.8))

    t  = np.array(times_ns)
    rg = np.array(rg_values)

    ax.plot(t, rg, color=PALETTE["purple"], alpha=0.85, label="Rg")
    ax.axhline(np.mean(rg), color=PALETTE["grey"], linestyle="--",
               linewidth=1.2, label=f"Mean = {np.mean(rg):.2f} Å")

    # Standard deviation band
    mu, sd = np.mean(rg), np.std(rg)
    ax.fill_between(t, mu - sd, mu + sd, color=PALETTE["purple"],
                    alpha=0.10, label=f"±1σ = {sd:.2f} Å")

    # Spike markers
    if spike_frame_indices:
        spike_mask = np.zeros(len(t), dtype=bool)
        for si in spike_frame_indices:
            nearest = np.argmin(np.abs(np.arange(len(t)) - si))
            if nearest < len(t):
                spike_mask[nearest] = True
        ax.scatter(t[spike_mask], rg[spike_mask],
                   color=PALETTE["red"], s=45, zorder=5,
                   marker="^", label="Spike detected")

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Radius of Gyration (Å)")
    ax.set_title("Protein Compactness (Radius of Gyration) vs Time")
    ax.legend(fontsize=9, loc="upper right")

    _finalize(fig, save_path or str(OUTPUT_DIR / "rg.png"))
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# COM Drift Plot
# ──────────────────────────────────────────────────────────────────────────────

def plot_com_drift(
    times_ns: List[float],
    displacements: List[float],
    warning_threshold: float = 20.0,
    critical_threshold: float = 50.0,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Center-of-mass displacement vs simulation time.

    Parameters
    ----------
    times_ns : list[float]
    displacements : list[float]
        COM displacement in Ångströms from frame-0.
    warning_threshold : float
        Horizontal line marking the warning level.
    critical_threshold : float
        Horizontal line marking the critical level.
    save_path : str, optional

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, ax = plt.subplots(figsize=(8, 3.8))

    t  = np.array(times_ns)
    d  = np.array(displacements)

    ax.plot(t, d, color=PALETTE["orange"], alpha=0.85, label="COM displacement")

    ax.axhline(warning_threshold, color=PALETTE["orange"], linestyle="--",
               linewidth=1.0, alpha=0.7, label=f"Warning ({warning_threshold:.0f} Å)")
    ax.axhline(critical_threshold, color=PALETTE["red"], linestyle="--",
               linewidth=1.0, alpha=0.7, label=f"Critical ({critical_threshold:.0f} Å)")

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("COM Displacement (Å)")
    ax.set_title("Center-of-Mass Drift vs Time")
    ax.legend(fontsize=9, loc="upper right")

    _finalize(fig, save_path or str(OUTPUT_DIR / "com_drift.png"))
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Health Score Gauge (Plotly)
# ──────────────────────────────────────────────────────────────────────────────

def plot_health_gauge(score: float) -> go.Figure:
    """
    Render an interactive Plotly gauge for the trajectory health score.

    Parameters
    ----------
    score : float  0–100

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if score >= 75:
        bar_color = "#2ecc71"
    elif score >= 50:
        bar_color = "#f39c12"
    else:
        bar_color = "#e74c3c"

    fig = go.Figure(go.Indicator(
        mode="gauge+number+delta",
        value=score,
        delta={"reference": 75, "suffix": " pts", "font": {"size": 14}},
        title={"text": "Trajectory Health Score", "font": {"size": 16}},
        gauge={
            "axis": {"range": [0, 100], "tickwidth": 1},
            "bar":  {"color": bar_color},
            "bgcolor": "white",
            "steps": [
                {"range": [0,  50], "color": "#fde8e8"},
                {"range": [50, 75], "color": "#fef6e4"},
                {"range": [75, 100],"color": "#e8f8ee"},
            ],
            "threshold": {
                "line": {"color": "black", "width": 2},
                "thickness": 0.75,
                "value": score,
            },
        },
        number={"suffix": "/100", "font": {"size": 28}},
    ))

    fig.update_layout(
        margin=dict(l=20, r=20, t=60, b=20),
        paper_bgcolor="rgba(0,0,0,0)",
        height=250,
    )
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Interactive Plotly versions (for Streamlit)
# ──────────────────────────────────────────────────────────────────────────────

def plotly_rmsd(
    times_ns: List[float],
    rmsd_values: List[float],
) -> go.Figure:
    """Interactive Plotly RMSD chart."""
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=times_ns, y=rmsd_values,
        mode="lines", name="RMSD",
        line=dict(color="#2980b9", width=1.5),
    ))
    fig.update_layout(
        xaxis_title="Time (ns)", yaxis_title="RMSD (Å)",
        title="Backbone RMSD vs Time",
        template="simple_white", height=350,
    )
    return fig


def plotly_rg(
    times_ns: List[float],
    rg_values: List[float],
) -> go.Figure:
    """Interactive Plotly Rg chart."""
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=times_ns, y=rg_values,
        mode="lines", name="Rg",
        line=dict(color="#8e44ad", width=1.5),
    ))
    mu = float(np.mean(rg_values))
    fig.add_hline(y=mu, line_dash="dash", line_color="grey",
                  annotation_text=f"Mean {mu:.2f} Å")
    fig.update_layout(
        xaxis_title="Time (ns)", yaxis_title="Rg (Å)",
        title="Radius of Gyration vs Time",
        template="simple_white", height=350,
    )
    return fig


def plotly_com_drift(
    times_ns: List[float],
    displacements: List[float],
) -> go.Figure:
    """Interactive Plotly COM drift chart."""
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=times_ns, y=displacements,
        mode="lines", name="COM Displacement",
        line=dict(color="#e67e22", width=1.5),
    ))
    fig.update_layout(
        xaxis_title="Time (ns)", yaxis_title="Displacement (Å)",
        title="Center-of-Mass Drift vs Time",
        template="simple_white", height=350,
    )
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

def _finalize(fig: plt.Figure, path: str) -> None:
    """Save figure and close it to free memory."""
    try:
        fig.savefig(path)
        logger.debug(f"Saved figure: {path}")
    except Exception as exc:
        logger.warning(f"Could not save figure to {path}: {exc}")
    # Do NOT close — Streamlit needs the Figure object returned to caller.
