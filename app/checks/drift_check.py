"""
drift_check.py — Center-of-Mass (COM) Drift Analyser.

Tracks the displacement of the protein's center of mass from its initial
position over the course of the trajectory.

Significant COM drift typically indicates:
  - Missing centering step in preprocessing
  - Improper application of PBC corrections
  - Forgotten -center flag in gmx trjconv

Fix: gmx trjconv -center -pbc mol
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List

import numpy as np
import MDAnalysis as mda

from app.utils.helpers import get_protein_selection, sample_frames

logger = logging.getLogger(__name__)

# Threshold for "significant" drift in Ångströms.
# Box sizes are typically 60–120 Å; >20 Å drift is noteworthy.
_DRIFT_WARNING_THRESHOLD_ANG  = 20.0
_DRIFT_CRITICAL_THRESHOLD_ANG = 50.0


@dataclass
class DriftCheckResult:
    """Results from the COM drift analysis."""

    passed: bool = True
    penalty: float = 0.0
    severity: str = "None"       # None | Low | Moderate | High
    max_drift_ang: float = 0.0
    mean_drift_ang: float = 0.0
    drift_trend: str = "Stable"  # Stable | Increasing | Oscillating
    times_ns: List[float] = field(default_factory=list)
    displacements_ang: List[float] = field(default_factory=list)
    details: str = ""
    fix_commands: List[str] = field(default_factory=list)


def run_drift_check(
    u: mda.Universe,
    warning_threshold: float = _DRIFT_WARNING_THRESHOLD_ANG,
    critical_threshold: float = _DRIFT_CRITICAL_THRESHOLD_ANG,
    max_frames: int = 500,
) -> DriftCheckResult:
    """
    Calculate protein COM displacement over time.

    The first sampled frame is used as the reference position.
    Displacement is the Euclidean distance from that reference COM.

    Parameters
    ----------
    u : mda.Universe
        Loaded MDAnalysis Universe.
    warning_threshold : float
        Drift in Å above which a Low/Moderate warning is issued.
    critical_threshold : float
        Drift in Å above which a High warning is issued.
    max_frames : int
        Maximum frames to sample (evenly spaced).

    Returns
    -------
    DriftCheckResult
    """
    result = DriftCheckResult()

    protein = get_protein_selection(u)
    if protein is None:
        result.details = "No protein atoms found; COM drift check skipped."
        logger.warning(result.details)
        return result

    frame_indices = sample_frames(u, max_frames)
    coms: List[np.ndarray] = []
    times: List[float] = []

    for fi in frame_indices:
        u.trajectory[fi]
        coms.append(protein.center_of_mass().copy())
        times.append(u.trajectory.time / 1000.0)   # ps → ns

    coms_arr = np.array(coms)          # (n_frames, 3)
    ref_com  = coms_arr[0]             # reference = first sampled frame

    displacements = np.linalg.norm(coms_arr - ref_com, axis=1)

    result.times_ns           = times
    result.displacements_ang  = displacements.tolist()
    result.max_drift_ang      = float(np.max(displacements))
    result.mean_drift_ang     = float(np.mean(displacements))
    result.drift_trend        = _classify_trend(displacements)

    # ── Severity classification ───────────────────────────────────────────────
    if result.max_drift_ang < warning_threshold:
        result.passed   = True
        result.severity = "None"
        result.penalty  = 0.0
        result.details  = (
            f"COM drift is within acceptable limits. "
            f"Max displacement = {result.max_drift_ang:.2f} Å."
        )

    elif result.max_drift_ang < critical_threshold:
        result.passed   = False
        result.severity = "Moderate"
        result.penalty  = 10.0
        result.details  = (
            f"Moderate COM drift detected "
            f"(max = {result.max_drift_ang:.2f} Å, "
            f"mean = {result.mean_drift_ang:.2f} Å). "
            f"Trend: {result.drift_trend}. "
            "Consider applying -center in gmx trjconv."
        )

    else:
        result.passed   = False
        result.severity = "High"
        result.penalty  = 20.0
        result.details  = (
            f"Large COM drift detected "
            f"(max = {result.max_drift_ang:.2f} Å, "
            f"mean = {result.mean_drift_ang:.2f} Å). "
            f"Trend: {result.drift_trend}. "
            "The protein has drifted significantly from its starting position — "
            "centering must be applied before analysis."
        )

    # ── Fix recommendations ───────────────────────────────────────────────────
    if not result.passed:
        result.fix_commands = _build_fix_commands()

    logger.info(
        f"Drift check: severity={result.severity}, "
        f"max_drift={result.max_drift_ang:.2f} Å, "
        f"trend={result.drift_trend}"
    )
    return result


def _classify_trend(displacements: np.ndarray) -> str:
    """
    Classify the drift trend as Stable, Increasing, or Oscillating.

    Uses a simple linear regression slope relative to the mean
    and the ratio of std-dev to mean to separate monotonic from
    oscillating trajectories.

    Parameters
    ----------
    displacements : np.ndarray
        1-D array of COM displacements in Å.

    Returns
    -------
    str
        'Stable' | 'Increasing' | 'Oscillating'
    """
    if len(displacements) < 4:
        return "Stable"

    n = len(displacements)
    x = np.arange(n, dtype=float)
    slope = np.polyfit(x, displacements, 1)[0]

    # Normalise slope by mean displacement to get relative trend
    mean_disp = np.mean(displacements) + 1e-9   # avoid div-by-zero
    rel_slope  = slope * n / mean_disp

    # Oscillation: high coefficient of variation and non-monotonic
    cv = np.std(displacements) / mean_disp
    if cv > 0.4 and rel_slope < 2.0:
        return "Oscillating"
    elif rel_slope > 1.5:
        return "Increasing"
    else:
        return "Stable"


def _build_fix_commands() -> List[str]:
    """Return the standard fix commands for COM drift."""
    return [
        "# Centre the protein in the box and remove COM translation\n"
        "gmx trjconv -s topol.tpr -f traj.xtc -o traj_center.xtc \\\n"
        "            -center -pbc mol",

        "# Optionally fit out rotation too\n"
        "gmx trjconv -s topol.tpr -f traj_center.xtc \\\n"
        "            -o traj_fit.xtc \\\n"
        "            -fit rot+trans",
    ]
