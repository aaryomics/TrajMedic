"""
helpers.py — Shared utility functions for TrajMedic.

Provides trajectory loading, validation, and metadata extraction
using MDAnalysis. This is the entry point for all file I/O.
"""

import os
import logging
from typing import Optional, Tuple, Dict, Any

import numpy as np
import MDAnalysis as mda
from MDAnalysis.exceptions import NoDataError

logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# Trajectory Loading
# ──────────────────────────────────────────────────────────────────────────────

def load_universe(
    topology_path: str,
    trajectory_path: str
) -> Optional[mda.Universe]:
    """
    Load an MDAnalysis Universe from topology and trajectory files.

    Parameters
    ----------
    topology_path : str
        Path to the GROMACS .tpr topology file.
    trajectory_path : str
        Path to the .xtc trajectory file.

    Returns
    -------
    mda.Universe or None
        Loaded Universe object, or None on failure.
    """
    if not os.path.isfile(topology_path):
        logger.error(f"Topology file not found: {topology_path}")
        return None
    if not os.path.isfile(trajectory_path):
        logger.error(f"Trajectory file not found: {trajectory_path}")
        return None

    try:
        u = mda.Universe(topology_path, trajectory_path)
        logger.info(
            f"Universe loaded: {u.atoms.n_atoms} atoms, "
            f"{u.trajectory.n_frames} frames"
        )
        return u
    except Exception as exc:
        logger.error(f"Failed to load Universe: {exc}")
        return None


# ──────────────────────────────────────────────────────────────────────────────
# Metadata Extraction
# ──────────────────────────────────────────────────────────────────────────────

def extract_metadata(u: mda.Universe) -> Dict[str, Any]:
    """
    Extract basic trajectory and system metadata.

    Parameters
    ----------
    u : mda.Universe
        A loaded MDAnalysis Universe.

    Returns
    -------
    dict
        Dictionary containing atom count, frame count, duration,
        timestep, box dimensions, and protein residue count.
    """
    traj = u.trajectory

    # Attempt to select protein residues
    try:
        protein = u.select_atoms("protein")
        n_residues = protein.residues.n_residues
    except Exception:
        n_residues = 0

    # Box dimensions from the first frame
    traj[0]
    box = traj.ts.dimensions  # [lx, ly, lz, alpha, beta, gamma] or None

    # Estimate duration from dt and n_frames
    try:
        dt_ps = traj.dt          # picoseconds per frame
        duration_ns = (traj.n_frames * dt_ps) / 1000.0
    except (AttributeError, NoDataError):
        dt_ps = None
        duration_ns = None

    metadata = {
        "n_atoms":      u.atoms.n_atoms,
        "n_frames":     traj.n_frames,
        "n_residues":   n_residues,
        "dt_ps":        dt_ps,
        "duration_ns":  duration_ns,
        "box_dims":     box.tolist() if box is not None else None,
    }

    logger.debug(f"Metadata: {metadata}")
    return metadata


# ──────────────────────────────────────────────────────────────────────────────
# Frame Sampling Utilities
# ──────────────────────────────────────────────────────────────────────────────

def get_time_array(u: mda.Universe) -> np.ndarray:
    """
    Return an array of simulation times in nanoseconds for all frames.

    Parameters
    ----------
    u : mda.Universe

    Returns
    -------
    np.ndarray
        Array of shape (n_frames,) with time in ns.
    """
    times = []
    for ts in u.trajectory:
        times.append(ts.time)   # time in picoseconds
    return np.array(times) / 1000.0   # convert to ns


def sample_frames(u: mda.Universe, max_frames: int = 500) -> np.ndarray:
    """
    Return evenly spaced frame indices to cap expensive calculations.

    Parameters
    ----------
    u : mda.Universe
    max_frames : int
        Maximum number of frames to analyse. If the trajectory has
        fewer frames than max_frames, all frames are returned.

    Returns
    -------
    np.ndarray
        Array of integer frame indices.
    """
    n = u.trajectory.n_frames
    if n <= max_frames:
        return np.arange(n)
    return np.linspace(0, n - 1, max_frames, dtype=int)


# ──────────────────────────────────────────────────────────────────────────────
# Selection Helpers
# ──────────────────────────────────────────────────────────────────────────────

def get_protein_selection(u: mda.Universe) -> Optional[mda.AtomGroup]:
    """
    Return the protein atom group, or None if no protein atoms exist.

    Parameters
    ----------
    u : mda.Universe

    Returns
    -------
    mda.AtomGroup or None
    """
    ag = u.select_atoms("protein")
    if len(ag) == 0:
        logger.warning("No 'protein' atoms found. Trying 'all'.")
        ag = u.select_atoms("all")
    return ag if len(ag) > 0 else None


def get_calpha_selection(u: mda.Universe) -> Optional[mda.AtomGroup]:
    """
    Return C-alpha atoms, falling back to all backbone if unavailable.

    Parameters
    ----------
    u : mda.Universe

    Returns
    -------
    mda.AtomGroup or None
    """
    ca = u.select_atoms("protein and name CA")
    if len(ca) == 0:
        logger.warning("No CA atoms found — trying backbone.")
        ca = u.select_atoms("backbone")
    return ca if len(ca) > 0 else None


# ──────────────────────────────────────────────────────────────────────────────
# Numeric Utilities
# ──────────────────────────────────────────────────────────────────────────────

def detect_jumps(
    series: np.ndarray,
    threshold_sigma: float = 4.0
) -> np.ndarray:
    """
    Identify indices in a 1-D time series where the value deviates
    more than `threshold_sigma` standard deviations from the mean.

    Parameters
    ----------
    series : np.ndarray
        1-D numerical array.
    threshold_sigma : float
        Number of standard deviations above which a point is flagged.

    Returns
    -------
    np.ndarray
        Array of integer indices flagged as outliers.
    """
    if len(series) < 3:
        return np.array([], dtype=int)
    mu, sigma = np.mean(series), np.std(series)
    if sigma == 0:
        return np.array([], dtype=int)
    z_scores = np.abs((series - mu) / sigma)
    return np.where(z_scores > threshold_sigma)[0]


def score_from_penalty(total_penalty: float) -> float:
    """
    Convert total penalty points to a 0–100 health score.

    Parameters
    ----------
    total_penalty : float
        Sum of all penalty points (must be >= 0).

    Returns
    -------
    float
        Health score clamped to [0, 100].
    """
    return float(np.clip(100.0 - total_penalty, 0, 100))


def health_label(score: float) -> Tuple[str, str]:
    """
    Map a health score to a human-readable label and a hex color.

    Parameters
    ----------
    score : float

    Returns
    -------
    (label, hex_color)
    """
    if score >= 75:
        return "Good", "#2ecc71"
    elif score >= 50:
        return "Moderate", "#f39c12"
    else:
        return "Poor", "#e74c3c"
