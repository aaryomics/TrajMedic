"""
main.py — TrajMedic Command-Line Entry Point.

Orchestrates the full trajectory diagnostic pipeline:

  1. Load trajectory and topology
  2. Extract system metadata
  3. Run all checks (PBC, drift, compactness, RMSD, integrity)
  4. Aggregate health score
  5. Generate plots
  6. Write reports

Usage
-----
  python -m app.main --tpr examples/sample.tpr --xtc examples/traj.xtc

The Streamlit UI (app/ui/dashboard.py) calls run_full_analysis() directly.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from app.utils.helpers import (
    load_universe,
    extract_metadata,
    score_from_penalty,
    health_label,
)
from app.checks import (
    run_pbc_check,
    run_drift_check,
    run_compactness_check,
    run_rmsd_check,
    run_integrity_check,
)
from app.utils.plotting import (
    plot_rmsd,
    plot_rg,
    plot_com_drift,
)
from app.utils.report import generate_txt_report, generate_html_report

# ── Logging setup ─────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(name)s — %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("trajmedic")


# ──────────────────────────────────────────────────────────────────────────────
# Core pipeline
# ──────────────────────────────────────────────────────────────────────────────

def run_full_analysis(
    tpr_path: str,
    xtc_path: str,
    max_frames: int = 500,
    progress_callback=None,
) -> Dict[str, Any]:
    """
    Run the complete TrajMedic diagnostic pipeline.

    This function is called by both the CLI and the Streamlit dashboard.

    Parameters
    ----------
    tpr_path : str
        Path to the GROMACS .tpr topology file.
    xtc_path : str
        Path to the .xtc trajectory file.
    max_frames : int
        Maximum frames to sample per check (for speed).
    progress_callback : callable, optional
        A callable(step: int, total: int, message: str) for progress updates.
        Used by the Streamlit UI progress bar.

    Returns
    -------
    dict containing:
        metadata, results, health_score, health_label,
        figures (matplotlib), plotly_figs, txt_report, html_report
    """
    total_steps = 8

    def _progress(step: int, msg: str) -> None:
        if progress_callback:
            progress_callback(step, total_steps, msg)
        logger.info(f"[{step}/{total_steps}] {msg}")

    # ── Step 1 — Load universe ─────────────────────────────────────────────
    _progress(1, "Loading trajectory and topology…")
    u = load_universe(tpr_path, xtc_path)
    if u is None:
        raise RuntimeError(
            f"Failed to load universe from:\n  TPR: {tpr_path}\n  XTC: {xtc_path}"
        )

    # ── Step 2 — Extract metadata ──────────────────────────────────────────
    _progress(2, "Extracting system metadata…")
    metadata = extract_metadata(u)

    # ── Step 3 — PBC check ─────────────────────────────────────────────────
    _progress(3, "Checking for PBC fragmentation…")
    pbc_res = run_pbc_check(u, max_frames=max_frames)

    # ── Step 4 — COM drift ─────────────────────────────────────────────────
    _progress(4, "Analysing center-of-mass drift…")
    drift_res = run_drift_check(u, max_frames=max_frames)

    # ── Step 5 — Compactness ───────────────────────────────────────────────
    _progress(5, "Computing radius of gyration…")
    comp_res = run_compactness_check(u, max_frames=max_frames)

    # ── Step 6 — RMSD ──────────────────────────────────────────────────────
    _progress(6, "Calculating backbone RMSD…")
    rmsd_res = run_rmsd_check(u, max_frames=max_frames)

    # ── Step 7 — Structural integrity ──────────────────────────────────────
    _progress(7, "Running structural integrity check…")
    integ_res = run_integrity_check(u, max_frames=max_frames)

    # ── Aggregate health score ─────────────────────────────────────────────
    total_penalty = (
        pbc_res.penalty
        + drift_res.penalty
        + comp_res.penalty
        + rmsd_res.penalty
        + integ_res.penalty
    )
    score  = score_from_penalty(total_penalty)
    label, color = health_label(score)

    results = {
        "pbc":         pbc_res,
        "drift":       drift_res,
        "compactness": comp_res,
        "rmsd":        rmsd_res,
        "integrity":   integ_res,
    }

    # ── Step 8 — Plots + Reports ───────────────────────────────────────────
    _progress(8, "Generating plots and reports…")

    figs = {}     # matplotlib figures
    pfigs = {}    # plotly figures

    if rmsd_res.times_ns:
        from app.utils.plotting import plot_rmsd, plotly_rmsd
        figs["rmsd"]  = plot_rmsd(
            rmsd_res.times_ns,
            rmsd_res.rmsd_values_ang,
            jump_frame_indices=rmsd_res.jump_frame_indices,
            equilibration_time_ns=rmsd_res.equilibration_time_ns,
        )
        pfigs["rmsd"] = plotly_rmsd(rmsd_res.times_ns, rmsd_res.rmsd_values_ang)

    if comp_res.times_ns:
        from app.utils.plotting import plot_rg, plotly_rg
        figs["rg"]  = plot_rg(
            comp_res.times_ns,
            comp_res.rg_values_ang,
            spike_frame_indices=comp_res.spike_frame_indices,
        )
        pfigs["rg"] = plotly_rg(comp_res.times_ns, comp_res.rg_values_ang)

    if drift_res.times_ns:
        from app.utils.plotting import plot_com_drift, plotly_com_drift
        figs["com_drift"]  = plot_com_drift(
            drift_res.times_ns,
            drift_res.displacements_ang,
        )
        pfigs["com_drift"] = plotly_com_drift(
            drift_res.times_ns,
            drift_res.displacements_ang,
        )

    # Plot paths for HTML embedding
    plot_paths = {
        "rmsd":      "outputs/rmsd.png",
        "rg":        "outputs/rg.png",
        "com_drift": "outputs/com_drift.png",
    }

    txt_report  = generate_txt_report(metadata, results, score, label)
    html_report = generate_html_report(
        metadata, results, score, label, plot_paths=plot_paths
    )

    # ── Explicitly close the Universe so Windows releases the file handle ──
    # On Windows, MDAnalysis keeps the .xtc/.tpr file open until the Universe
    # is garbage-collected or explicitly closed. Without this, deleting or
    # overwriting the session files later causes a PermissionError.
    try:
        u.trajectory.close()
    except Exception:
        pass

    return {
        "metadata":     metadata,
        "results":      results,
        "health_score": score,
        "health_label": label,
        "health_color": color,
        "figures":      figs,
        "plotly_figs":  pfigs,
        "txt_report":   txt_report,
        "html_report":  html_report,
    }


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="trajmedic",
        description="TrajMedic — GROMACS Trajectory Diagnostic Tool",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpr",        required=True, help="Path to .tpr topology file")
    p.add_argument("--xtc",        required=True, help="Path to .xtc trajectory file")
    p.add_argument("--max-frames", type=int, default=500,
                   help="Max frames to sample per check")
    p.add_argument("--no-html",    action="store_true",
                   help="Skip HTML report generation")
    return p


def main() -> None:
    parser = _build_parser()
    args   = parser.parse_args()

    logger.info("=" * 60)
    logger.info("  TrajMedic — Trajectory Diagnostic Tool")
    logger.info("=" * 60)
    logger.info(f"  TPR : {args.tpr}")
    logger.info(f"  XTC : {args.xtc}")

    try:
        output = run_full_analysis(
            tpr_path=args.tpr,
            xtc_path=args.xtc,
            max_frames=args.max_frames,
        )
    except RuntimeError as exc:
        logger.error(str(exc))
        sys.exit(1)

    # Print summary to stdout
    score = output["health_score"]
    label = output["health_label"]
    print("\n" + "=" * 60)
    print(f"  Trajectory Health Score: {score:.1f}/100  [{label.upper()}]")
    print("=" * 60)

    for key, res in output["results"].items():
        status = "PASS" if res.passed else "FAIL"
        print(f"  [{status}] {key:<12}  severity={res.severity}  penalty=-{res.penalty:.0f}")

    print("\n  Reports saved to outputs/")
    print("  Plots  saved to outputs/*.png")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
