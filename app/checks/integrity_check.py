"""
integrity_check.py — Structural Integrity and Anomaly Detector.

Flags extreme per-atom displacement events, "exploded" simulation frames,
and anomalous residue separations that can result from:
  - Numerical instabilities (e.g., bad timestep or force-field parameters)
  - Corrupted trajectory files
  - Missed energy minimisation before production run
  - Poorly restrained systems

Checks performed:
  1. Per-frame maximum single-atom displacement (inter-frame).
  2. Overall simulation box dimension consistency.
  3. Detection of "explosion" events (atoms flying beyond 2× box diagonal).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
import MDAnalysis as mda

from app.utils.helpers import get_protein_selection, sample_frames, detect_jumps

logger = logging.getLogger(__name__)

# Per-frame maximum atom displacement threshold.
# A single atom moving > 10 Å between consecutive frames is suspicious.
_ATOM_DISP_THRESHOLD_ANG = 10.0

# If any atom is located further than (2 × box_diagonal) from the COM,
# the frame is considered "exploded".
_EXPLOSION_FACTOR = 2.0


@dataclass
class IntegrityCheckResult:
    """Results from the structural integrity analysis."""

    passed: bool = True
    penalty: float = 0.0
    severity: str = "None"         # None | Low | Moderate | High
    n_anomalous_frames: int = 0
    anomalous_frames: List[int] = field(default_factory=list)
    n_explosion_frames: int = 0
    explosion_frames: List[int] = field(default_factory=list)
    max_single_atom_disp: float = 0.0
    box_dimension_variance: float = 0.0   # CV of box diagonal over time
    details: str = ""
    fix_commands: List[str] = field(default_factory=list)


def run_integrity_check(
    u: mda.Universe,
    disp_threshold: float = _ATOM_DISP_THRESHOLD_ANG,
    max_frames: int = 500,
) -> IntegrityCheckResult:
    """
    Detect structural anomalies and simulation explosion events.

    Parameters
    ----------
    u : mda.Universe
        Loaded MDAnalysis Universe.
    disp_threshold : float
        Maximum per-atom inter-frame displacement (Å) before flagging.
    max_frames : int
        Maximum frames to sample.

    Returns
    -------
    IntegrityCheckResult
    """
    result = IntegrityCheckResult()

    protein = get_protein_selection(u)
    if protein is None:
        result.details = "No protein atoms; integrity check skipped."
        logger.warning(result.details)
        return result

    frame_indices = sample_frames(u, max_frames)

    anomalous_frames: List[int]  = []
    explosion_frames: List[int]  = []
    max_disps: List[float]       = []
    box_diagonals: List[float]   = []
    prev_positions: Optional[np.ndarray] = None

    for i, fi in enumerate(frame_indices):
        u.trajectory[fi]
        curr_positions = protein.positions.copy()   # (n_atoms, 3)

        # ── Box diagonal ──────────────────────────────────────────────────────
        dims = u.trajectory.ts.dimensions
        if dims is not None:
            diag = np.sqrt(np.sum(dims[:3] ** 2))
            box_diagonals.append(diag)
        else:
            diag = None

        # ── Explosion check ───────────────────────────────────────────────────
        if diag is not None:
            com = curr_positions.mean(axis=0)
            atom_dists_from_com = np.linalg.norm(curr_positions - com, axis=1)
            if np.any(atom_dists_from_com > _EXPLOSION_FACTOR * diag):
                explosion_frames.append(int(fi))

        # ── Inter-frame displacement check ────────────────────────────────────
        if prev_positions is not None:
            # Only valid if frame step is small (consecutive sampled frames)
            per_atom_disp = np.linalg.norm(
                curr_positions - prev_positions, axis=1
            )
            max_disp = float(np.max(per_atom_disp))
            max_disps.append(max_disp)

            if max_disp > disp_threshold:
                anomalous_frames.append(int(fi))

        prev_positions = curr_positions

    # ── Aggregate stats ───────────────────────────────────────────────────────
    result.anomalous_frames    = anomalous_frames
    result.n_anomalous_frames  = len(anomalous_frames)
    result.explosion_frames    = explosion_frames
    result.n_explosion_frames  = len(explosion_frames)
    result.max_single_atom_disp = float(np.max(max_disps)) if max_disps else 0.0

    if len(box_diagonals) > 1:
        bd = np.array(box_diagonals)
        result.box_dimension_variance = float(np.std(bd) / np.mean(bd))
    else:
        result.box_dimension_variance = 0.0

    # ── Classify severity ─────────────────────────────────────────────────────
    total_issues = result.n_anomalous_frames + result.n_explosion_frames * 3
    n_sampled    = max(len(frame_indices) - 1, 1)
    issue_frac   = total_issues / n_sampled

    if issue_frac == 0:
        result.passed   = True
        result.severity = "None"
        result.penalty  = 0.0
        result.details  = "No structural integrity issues detected."

    elif issue_frac < 0.05:
        result.passed   = True
        result.severity = "Low"
        result.penalty  = 5.0
        result.details  = (
            f"Minor anomalies in {result.n_anomalous_frames} frame(s). "
            f"Max single-atom inter-frame displacement: "
            f"{result.max_single_atom_disp:.2f} Å."
        )

    elif issue_frac < 0.20 or result.n_explosion_frames < 3:
        result.passed   = False
        result.severity = "Moderate"
        result.penalty  = 10.0
        result.details  = (
            f"Structural anomalies in {result.n_anomalous_frames} frame(s), "
            f"{result.n_explosion_frames} potential explosion event(s). "
            f"Max displacement: {result.max_single_atom_disp:.2f} Å."
        )

    else:
        result.passed   = False
        result.severity = "High"
        result.penalty  = 20.0
        result.details  = (
            f"Severe structural instability: "
            f"{result.n_anomalous_frames} anomalous frames, "
            f"{result.n_explosion_frames} explosion event(s). "
            f"Max displacement: {result.max_single_atom_disp:.2f} Å. "
            "Trajectory is likely corrupted or un-minimised."
        )

    if not result.passed:
        result.fix_commands = _build_fix_commands(result.severity)

    logger.info(
        f"Integrity check: severity={result.severity}, "
        f"anomalous={result.n_anomalous_frames}, "
        f"explosions={result.n_explosion_frames}"
    )
    return result


def _build_fix_commands(severity: str) -> List[str]:
    """Return fix recommendations for structural integrity failures."""
    cmds = [
        "# Verify energy minimisation was performed before production\n"
        "gmx grompp -f minim.mdp -c confout.gro -p topol.top -o em.tpr\n"
        "gmx mdrun  -v -deffnm em",
    ]

    if severity == "High":
        cmds += [
            "# Check for bad contacts using gmx check\n"
            "gmx check -f traj.xtc",

            "# Reduce timestep and re-run simulation\n"
            "# In your .mdp file: set dt = 0.001  (1 fs instead of 2 fs)",
        ]

    return cmds
