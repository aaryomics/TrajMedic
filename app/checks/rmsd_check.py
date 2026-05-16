"""
rmsd_check.py — Backbone RMSD Stability Analyser.

Root-Mean-Square Deviation (RMSD) of the protein backbone measures how
much the structure has changed relative to a reference (typically the
first frame or the crystal structure).

Interpretation guidelines:
  - A flat plateau after an initial rise indicates equilibration.
  - Continuous, unbounded growth suggests the protein is still evolving
    (simulation too short, or instability).
  - Sudden discrete jumps are PBC or preprocessing artefacts.
  - Very low RMSD (<0.5 Å) may indicate frozen/constrained atoms.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

from app.utils.helpers import get_calpha_selection, sample_frames, detect_jumps

logger = logging.getLogger(__name__)

_RMSD_LOW_THRESHOLD    = 3.0    # Å  — trajectory probably equilibrated
_RMSD_MEDIUM_THRESHOLD = 6.0    # Å  — notable structural change
_RMSD_HIGH_THRESHOLD   = 10.0   # Å  — severe instability / likely artefact
_JUMP_SIGMA            = 4.0    # σ  — for inter-frame jump detection


@dataclass
class RMSDCheckResult:
    """Results from the backbone RMSD analysis."""

    passed: bool = True
    penalty: float = 0.0
    severity: str = "None"
    mean_rmsd_ang: float = 0.0
    max_rmsd_ang: float = 0.0
    equilibration_frame: int = -1     # -1 = not detected
    equilibration_time_ns: float = -1.0
    n_jump_frames: int = 0
    jump_frame_indices: List[int] = field(default_factory=list)
    times_ns: List[float] = field(default_factory=list)
    rmsd_values_ang: List[float] = field(default_factory=list)
    interpretation: str = ""
    details: str = ""
    fix_commands: List[str] = field(default_factory=list)


def run_rmsd_check(
    u: mda.Universe,
    max_frames: int = 500,
) -> RMSDCheckResult:
    """
    Compute backbone RMSD after superimposing each frame on the first.

    Steps:
    1.  Select backbone / CA atoms.
    2.  Align every frame to frame-0 (fit rot+trans).
    3.  Compute RMSD from the aligned positions.
    4.  Classify severity and generate interpretation.

    Parameters
    ----------
    u : mda.Universe
        Loaded MDAnalysis Universe.
    max_frames : int
        Maximum frames to sample.

    Returns
    -------
    RMSDCheckResult
    """
    result = RMSDCheckResult()

    ca = get_calpha_selection(u)
    if ca is None or len(ca) < 5:
        result.details = "Too few CA / backbone atoms for RMSD analysis."
        logger.warning(result.details)
        return result

    frame_indices = sample_frames(u, max_frames)

    # ── Compute RMSD ──────────────────────────────────────────────────────────
    # Store the reference (frame-0) positions
    u.trajectory[frame_indices[0]]
    ref_positions = ca.positions.copy()

    rmsd_values: List[float] = []
    times: List[float] = []

    for fi in frame_indices:
        u.trajectory[fi]
        mobile_positions = ca.positions.copy()

        # Superimpose mobile onto reference (Kabsch / least-squares fit)
        rmsd_val = _calc_rmsd(ref_positions, mobile_positions)
        rmsd_values.append(rmsd_val)
        times.append(u.trajectory.time / 1000.0)

    rmsd_arr = np.array(rmsd_values)

    result.times_ns       = times
    result.rmsd_values_ang = rmsd_values
    result.mean_rmsd_ang   = float(np.mean(rmsd_arr))
    result.max_rmsd_ang    = float(np.max(rmsd_arr))

    # ── Jump detection ────────────────────────────────────────────────────────
    # Compute frame-to-frame RMSD differences to find sudden jumps
    inter_frame_delta = np.abs(np.diff(rmsd_arr))
    jump_indices_delta = detect_jumps(inter_frame_delta, threshold_sigma=_JUMP_SIGMA)
    # Map back to the sampled frame indices
    result.jump_frame_indices = [frame_indices[i] for i in jump_indices_delta]
    result.n_jump_frames      = len(result.jump_frame_indices)

    # ── Equilibration detection ───────────────────────────────────────────────
    eq_idx = _find_equilibration_point(rmsd_arr)
    if eq_idx >= 0:
        result.equilibration_frame    = int(frame_indices[eq_idx])
        result.equilibration_time_ns  = float(times[eq_idx])

    # ── Severity classification ───────────────────────────────────────────────
    _classify_rmsd(result, rmsd_arr)

    # ── Human-readable interpretation ─────────────────────────────────────────
    result.interpretation = _interpret(result)

    logger.info(
        f"RMSD check: mean={result.mean_rmsd_ang:.2f} Å, "
        f"max={result.max_rmsd_ang:.2f} Å, "
        f"jumps={result.n_jump_frames}, severity={result.severity}"
    )
    return result


# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

def _calc_rmsd(ref: np.ndarray, mobile: np.ndarray) -> float:
    """
    Calculate RMSD between two equal-sized sets of 3-D coordinates
    after optimal superimposition (Kabsch algorithm).

    Parameters
    ----------
    ref : np.ndarray  shape (N, 3)
    mobile : np.ndarray  shape (N, 3)

    Returns
    -------
    float  RMSD in Ångströms.
    """
    # Centre both sets
    ref_c    = ref    - ref.mean(axis=0)
    mobile_c = mobile - mobile.mean(axis=0)

    # Rotation matrix via SVD (Kabsch)
    H   = mobile_c.T @ ref_c
    U, S, Vt = np.linalg.svd(H)
    d   = np.linalg.det(Vt.T @ U.T)
    D   = np.diag([1, 1, d])
    R   = Vt.T @ D @ U.T

    mobile_rot = mobile_c @ R.T
    diff       = ref_c - mobile_rot
    rmsd_val   = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
    return rmsd_val


def _find_equilibration_point(rmsd: np.ndarray, window: int = 20) -> int:
    """
    Estimate the frame index at which the trajectory equilibrates.

    Uses a rolling-window variance approach: equilibration is assumed
    at the first window whose variance drops below 10 % of the overall
    variance and whose mean has stabilised.

    Returns -1 if equilibration cannot be determined.

    Parameters
    ----------
    rmsd : np.ndarray  1-D RMSD array.
    window : int       Rolling window size (in frames).

    Returns
    -------
    int  Frame index of equilibration onset, or -1.
    """
    if len(rmsd) < window * 2:
        return -1

    global_var = np.var(rmsd)
    if global_var == 0:
        return -1

    for i in range(window, len(rmsd) - window):
        local_var = np.var(rmsd[i: i + window])
        if local_var < 0.10 * global_var:
            return i

    return -1


def _classify_rmsd(result: RMSDCheckResult, rmsd_arr: np.ndarray) -> None:
    """
    Mutate *result* with severity, penalty, and details based on RMSD values.
    """
    max_r  = result.max_rmsd_ang
    jumps  = result.n_jump_frames

    if max_r <= _RMSD_LOW_THRESHOLD and jumps == 0:
        result.passed   = True
        result.severity = "None"
        result.penalty  = 0.0
        result.details  = (
            f"RMSD is low and stable (max = {max_r:.2f} Å). "
            "Trajectory appears well-equilibrated."
        )

    elif max_r <= _RMSD_MEDIUM_THRESHOLD and jumps <= 3:
        result.passed   = True
        result.severity = "Low"
        result.penalty  = 5.0
        result.details  = (
            f"RMSD within acceptable range (max = {max_r:.2f} Å). "
            f"Minor fluctuations or {jumps} small jump(s) detected."
        )

    elif max_r <= _RMSD_HIGH_THRESHOLD or jumps <= 10:
        result.passed   = False
        result.severity = "Moderate"
        result.penalty  = 15.0
        result.details  = (
            f"Elevated RMSD (max = {max_r:.2f} Å). "
            f"{jumps} inter-frame jump(s) detected. "
            "Possible instability or incomplete equilibration."
        )

    else:
        result.passed   = False
        result.severity = "High"
        result.penalty  = 25.0
        result.details  = (
            f"Severe RMSD instability (max = {max_r:.2f} Å). "
            f"{jumps} large jump(s) detected. "
            "Trajectory may not be suitable for analysis without remediation."
        )

        result.fix_commands = [
            "# Fit each frame to a reference to remove rigid-body motions\n"
            "gmx trjconv -s topol.tpr -f traj.xtc -o traj_fit.xtc \\\n"
            "            -fit rot+trans",

            "# Verify equilibration by re-running gmx rms\n"
            "gmx rms -s topol.tpr -f traj_fit.xtc -o rmsd.xvg \\\n"
            "        -tu ns",
        ]


def _interpret(result: RMSDCheckResult) -> str:
    """Generate a plain-English interpretation sentence."""
    lines = []

    if result.severity == "None":
        lines.append(
            f"Trajectory is stable with a mean backbone RMSD of "
            f"{result.mean_rmsd_ang:.2f} Å."
        )
    elif result.severity == "Low":
        lines.append(
            f"Backbone RMSD shows minor variability "
            f"(mean {result.mean_rmsd_ang:.2f} Å, max {result.max_rmsd_ang:.2f} Å)."
        )
    else:
        lines.append(
            f"Large RMSD fluctuations detected "
            f"(mean {result.mean_rmsd_ang:.2f} Å, max {result.max_rmsd_ang:.2f} Å)."
        )

    if result.equilibration_time_ns > 0:
        lines.append(
            f"Trajectory appears to equilibrate around "
            f"{result.equilibration_time_ns:.1f} ns."
        )
    else:
        lines.append(
            "Equilibration point could not be determined — "
            "the simulation may need more sampling."
        )

    if result.n_jump_frames > 0:
        lines.append(
            f"{result.n_jump_frames} sudden inter-frame RMSD jump(s) detected; "
            "these may be preprocessing artefacts."
        )

    return "  ".join(lines)
