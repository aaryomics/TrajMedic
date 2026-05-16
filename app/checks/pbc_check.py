"""
pbc_check.py — Periodic Boundary Condition Fragmentation Detector.

Detects protein fragmentation artifacts caused by improper PBC handling.
These appear as unrealistically large inter-residue distances between
consecutive C-alpha atoms, often spanning the length of the simulation box.

Common cause: missing `gmx trjconv -pbc nojump` or `-pbc mol` step.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
import MDAnalysis as mda

from app.utils.helpers import get_calpha_selection, sample_frames

logger = logging.getLogger(__name__)

# Threshold in Ångströms for a "suspicious" CA–CA sequential distance.
# The expected CA–CA distance along the backbone is ~3.8 Å.
# We flag anything above 15 Å as a potential PBC artifact.
_CA_JUMP_THRESHOLD_ANG = 15.0


@dataclass
class PBCCheckResult:
    """Results from the PBC fragmentation analysis."""

    passed: bool = True
    penalty: float = 0.0
    severity: str = "None"                # None | Low | Moderate | High
    confidence: float = 0.0              # 0–1 probability of true fragmentation
    n_affected_frames: int = 0
    affected_frames: List[int] = field(default_factory=list)
    max_jump_ang: float = 0.0
    mean_jump_ang: float = 0.0
    details: str = ""
    fix_commands: List[str] = field(default_factory=list)


def run_pbc_check(
    u: mda.Universe,
    jump_threshold: float = _CA_JUMP_THRESHOLD_ANG,
    max_frames: int = 500,
) -> PBCCheckResult:
    """
    Analyse inter-residue C-alpha distances to detect PBC fragmentation.

    For each sampled frame, compute the distances between consecutive
    C-alpha atoms along the protein chain.  A large jump (> jump_threshold)
    implies that the protein has been split across the periodic boundary.

    Parameters
    ----------
    u : mda.Universe
        Loaded MDAnalysis Universe.
    jump_threshold : float
        Distance in Ångströms above which a CA–CA pair is flagged.
    max_frames : int
        Maximum number of frames to evaluate (for speed).

    Returns
    -------
    PBCCheckResult
        Populated result dataclass.
    """
    result = PBCCheckResult()

    ca = get_calpha_selection(u)
    if ca is None or len(ca) < 2:
        result.details = "Insufficient CA atoms for PBC check."
        logger.warning(result.details)
        return result

    frame_indices = sample_frames(u, max_frames)
    affected_frames: List[int] = []
    all_max_jumps: List[float] = []

    for fi in frame_indices:
        u.trajectory[fi]
        positions = ca.positions          # shape (n_ca, 3)

        # Compute distances between consecutive CA atoms
        diffs = np.diff(positions, axis=0)            # (n_ca-1, 3)
        dists = np.linalg.norm(diffs, axis=1)         # (n_ca-1,)

        max_dist = float(np.max(dists))
        all_max_jumps.append(max_dist)

        if max_dist > jump_threshold:
            affected_frames.append(int(fi))

    all_max_jumps_arr = np.array(all_max_jumps)

    result.affected_frames   = affected_frames
    result.n_affected_frames = len(affected_frames)
    result.max_jump_ang      = float(np.max(all_max_jumps_arr))
    result.mean_jump_ang     = float(np.mean(all_max_jumps_arr))

    frac_affected = result.n_affected_frames / max(len(frame_indices), 1)

    # ── Severity classification ───────────────────────────────────────────────
    if frac_affected == 0:
        result.passed     = True
        result.severity   = "None"
        result.confidence = 0.0
        result.penalty    = 0.0
        result.details    = (
            f"No PBC fragmentation detected. "
            f"Max CA–CA jump = {result.max_jump_ang:.2f} Å "
            f"(threshold {jump_threshold:.1f} Å)."
        )

    elif frac_affected < 0.10:
        result.passed     = True      # minor — pass with note
        result.severity   = "Low"
        result.confidence = 0.4
        result.penalty    = 10.0
        result.details    = (
            f"Occasional CA–CA jumps detected in "
            f"{result.n_affected_frames} frames "
            f"({frac_affected*100:.1f}% of sampled frames). "
            f"Max jump = {result.max_jump_ang:.2f} Å."
        )

    elif frac_affected < 0.40:
        result.passed     = False
        result.severity   = "Moderate"
        result.confidence = 0.75
        result.penalty    = 20.0
        result.details    = (
            f"Significant PBC fragmentation in "
            f"{result.n_affected_frames} frames "
            f"({frac_affected*100:.1f}%). "
            f"Max jump = {result.max_jump_ang:.2f} Å. "
            "The protein is likely split across the periodic boundary."
        )

    else:
        result.passed     = False
        result.severity   = "High"
        result.confidence = 0.95
        result.penalty    = 30.0
        result.details    = (
            f"Severe PBC fragmentation detected in "
            f"{result.n_affected_frames} frames "
            f"({frac_affected*100:.1f}%). "
            f"Max CA–CA jump = {result.max_jump_ang:.2f} Å. "
            "Trajectory must be re-wrapped before any analysis."
        )

    # ── Fix recommendations ───────────────────────────────────────────────────
    if not result.passed or result.severity != "None":
        result.fix_commands = _build_fix_commands(result.severity)

    logger.info(
        f"PBC check: severity={result.severity}, "
        f"penalty={result.penalty}, "
        f"affected={result.n_affected_frames}/{len(frame_indices)}"
    )
    return result


def _build_fix_commands(severity: str) -> List[str]:
    """
    Return a list of recommended GROMACS commands based on severity.

    Parameters
    ----------
    severity : str
        One of: None, Low, Moderate, High.

    Returns
    -------
    list[str]
    """
    commands = []

    if severity in ("Low", "Moderate", "High"):
        commands.append(
            "# Step 1 — Remove PBC jumps (keeps protein whole across boundaries)\n"
            "gmx trjconv -s topol.tpr -f traj.xtc -o traj_nojump.xtc \\\n"
            "            -pbc nojump"
        )

    if severity in ("Moderate", "High"):
        commands.append(
            "# Step 2 — Centre the protein and re-wrap the trajectory\n"
            "gmx trjconv -s topol.tpr -f traj_nojump.xtc \\\n"
            "            -o traj_center.xtc \\\n"
            "            -center -pbc mol -ur compact"
        )

    if severity == "High":
        commands.append(
            "# Step 3 — Fit rotational/translational degrees of freedom\n"
            "gmx trjconv -s topol.tpr -f traj_center.xtc \\\n"
            "            -o traj_fit.xtc \\\n"
            "            -fit rot+trans"
        )

    return commands
