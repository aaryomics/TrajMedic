"""
compactness_check.py — Protein Compactness Analyser (Radius of Gyration).

The radius of gyration (Rg) measures how "spread out" protein atoms are
from the protein's center of mass.

Sudden spikes or a monotonically increasing Rg trajectory can indicate:
  - Partial or full protein unfolding
  - PBC fragmentation artefacts inflating Rg
  - Simulation instability (e.g., bad force-field parameters)

A well-folded, stable protein typically shows a roughly constant Rg
with small fluctuations around the experimental value.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms

from app.utils.helpers import get_protein_selection, sample_frames, detect_jumps

logger = logging.getLogger(__name__)

# Rule-of-thumb thresholds for Rg variance (fraction of mean).
# > 5 % variance is a yellow flag; > 12 % is a red flag.
_VAR_WARN_FRAC    = 0.05
_VAR_CRITICAL_FRAC = 0.12

# Spike detection: flag frames where Rg deviates > 3 σ from the mean.
_SPIKE_SIGMA = 3.0


@dataclass
class CompactnessCheckResult:
    """Results from the Rg-based compactness analysis."""

    passed: bool = True
    penalty: float = 0.0
    severity: str = "None"        # None | Low | Moderate | High
    mean_rg_ang: float = 0.0
    std_rg_ang: float = 0.0
    max_rg_ang: float = 0.0
    min_rg_ang: float = 0.0
    n_spike_frames: int = 0
    spike_frame_indices: List[int] = field(default_factory=list)
    times_ns: List[float] = field(default_factory=list)
    rg_values_ang: List[float] = field(default_factory=list)
    details: str = ""
    fix_commands: List[str] = field(default_factory=list)


def run_compactness_check(
    u: mda.Universe,
    max_frames: int = 500,
) -> CompactnessCheckResult:
    """
    Compute the radius of gyration for the protein over time.

    Parameters
    ----------
    u : mda.Universe
        Loaded MDAnalysis Universe.
    max_frames : int
        Maximum number of frames to sample.

    Returns
    -------
    CompactnessCheckResult
    """
    result = CompactnessCheckResult()

    protein = get_protein_selection(u)
    if protein is None:
        result.details = "No protein atoms found; compactness check skipped."
        logger.warning(result.details)
        return result

    frame_indices = sample_frames(u, max_frames)
    rg_values: List[float] = []
    times: List[float] = []

    for fi in frame_indices:
        u.trajectory[fi]
        rg = protein.radius_of_gyration()   # Å
        rg_values.append(float(rg))
        times.append(u.trajectory.time / 1000.0)   # ps → ns

    rg_arr = np.array(rg_values)

    result.times_ns      = times
    result.rg_values_ang = rg_values
    result.mean_rg_ang   = float(np.mean(rg_arr))
    result.std_rg_ang    = float(np.std(rg_arr))
    result.max_rg_ang    = float(np.max(rg_arr))
    result.min_rg_ang    = float(np.min(rg_arr))

    # ── Spike detection ───────────────────────────────────────────────────────
    spike_indices = detect_jumps(rg_arr, threshold_sigma=_SPIKE_SIGMA)
    result.n_spike_frames     = len(spike_indices)
    result.spike_frame_indices = [frame_indices[i] for i in spike_indices]

    # Coefficient of variation (std / mean) for severity classification
    cv = result.std_rg_ang / max(result.mean_rg_ang, 1e-9)

    # ── Severity classification ───────────────────────────────────────────────
    if cv < _VAR_WARN_FRAC and result.n_spike_frames == 0:
        result.passed   = True
        result.severity = "None"
        result.penalty  = 0.0
        result.details  = (
            f"Protein compactness appears stable. "
            f"Rg = {result.mean_rg_ang:.2f} ± {result.std_rg_ang:.2f} Å "
            f"(CV = {cv*100:.1f}%)."
        )

    elif cv < _VAR_CRITICAL_FRAC or result.n_spike_frames < 5:
        result.passed   = False
        result.severity = "Low" if result.n_spike_frames == 0 else "Moderate"
        result.penalty  = 8.0
        result.details  = (
            f"Moderate Rg variability detected. "
            f"Rg = {result.mean_rg_ang:.2f} ± {result.std_rg_ang:.2f} Å "
            f"(CV = {cv*100:.1f}%). "
            f"Spike frames: {result.n_spike_frames}. "
            "May indicate partial unfolding or PBC artefacts."
        )

    else:
        result.passed   = False
        result.severity = "High"
        result.penalty  = 15.0
        result.details  = (
            f"Significant Rg instability detected. "
            f"Rg = {result.mean_rg_ang:.2f} ± {result.std_rg_ang:.2f} Å "
            f"(CV = {cv*100:.1f}%). "
            f"Spike frames: {result.n_spike_frames}. "
            "Possible unfolding event or severe PBC artefact."
        )

    # ── Fix recommendations ───────────────────────────────────────────────────
    if not result.passed:
        result.fix_commands = _build_fix_commands(result.severity)

    logger.info(
        f"Compactness check: Rg_mean={result.mean_rg_ang:.2f} Å, "
        f"CV={cv*100:.1f}%, severity={result.severity}"
    )
    return result


def _build_fix_commands(severity: str) -> List[str]:
    """Return recommended commands for compactness problems."""
    cmds = []

    cmds.append(
        "# Re-wrap to compact representation to fix visual fragmentation\n"
        "gmx trjconv -s topol.tpr -f traj.xtc -o traj_compact.xtc \\\n"
        "            -pbc mol -ur compact"
    )

    if severity == "High":
        cmds.append(
            "# Visualise Rg with GROMACS native tool for cross-check\n"
            "gmx gyrate -s topol.tpr -f traj.xtc -o gyrate.xvg"
        )
        cmds.append(
            "# If unfolding is suspected, verify with secondary structure\n"
            "gmx do_dssp -s topol.tpr -f traj.xtc -o dssp.xpm"
        )

    return cmds
