"""
report.py — Diagnostic Report Generator for TrajMedic.

Produces:
  1. A plain-text (.txt) report suitable for archiving or sharing.
  2. An HTML report with embedded plots and styled sections.

The report is generated from the aggregated results of all checks.
"""

from __future__ import annotations

import datetime
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
# Public entry points
# ──────────────────────────────────────────────────────────────────────────────

def generate_txt_report(
    metadata: Dict[str, Any],
    results: Dict[str, Any],
    health_score: float,
    health_label: str,
    output_path: Optional[str] = None,
) -> str:
    """
    Build and save a plain-text diagnostic report.

    Parameters
    ----------
    metadata : dict
        Output of helpers.extract_metadata().
    results : dict
        Mapping of check name → result dataclass.
    health_score : float
        Overall 0–100 health score.
    health_label : str
        'Good' | 'Moderate' | 'Poor'
    output_path : str, optional
        Override default output path.

    Returns
    -------
    str
        Full content of the generated report.
    """
    lines: List[str] = []

    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines += [
        "=" * 70,
        "  TrajMedic — Trajectory Diagnostic Report",
        f"  Generated: {ts}",
        "=" * 70,
        "",
    ]

    # ── System metadata ───────────────────────────────────────────────────────
    lines += [
        "SYSTEM SUMMARY",
        "-" * 40,
        f"  Atoms           : {metadata.get('n_atoms', 'N/A')}",
        f"  Frames          : {metadata.get('n_frames', 'N/A')}",
        f"  Protein residues: {metadata.get('n_residues', 'N/A')}",
        f"  Timestep (ps)   : {metadata.get('dt_ps', 'N/A')}",
        f"  Duration (ns)   : {_fmt(metadata.get('duration_ns'))}",
    ]
    box = metadata.get("box_dims")
    if box:
        lines.append(f"  Box (Å / deg)   : {[round(v, 2) for v in box]}")
    lines.append("")

    # ── Health score ──────────────────────────────────────────────────────────
    lines += [
        "TRAJECTORY HEALTH SCORE",
        "-" * 40,
        f"  Score : {health_score:.1f} / 100   [{health_label.upper()}]",
        "",
    ]

    # ── Per-check results ─────────────────────────────────────────────────────
    check_labels = {
        "pbc":         "PBC / Fragmentation Check",
        "drift":       "Center-of-Mass Drift",
        "compactness": "Protein Compactness (Rg)",
        "rmsd":        "Backbone RMSD Stability",
        "integrity":   "Structural Integrity",
    }

    lines += ["DIAGNOSTIC DETAILS", "-" * 40]
    for key, label in check_labels.items():
        res = results.get(key)
        if res is None:
            continue
        status = "✓ PASS" if res.passed else "✗ FAIL"
        lines += [
            f"  [{status}] {label}",
            f"     Severity : {res.severity}",
            f"     Penalty  : -{res.penalty:.0f} pts",
            f"     {res.details}",
            "",
        ]

    # ── Recommendations ───────────────────────────────────────────────────────
    lines += ["FIX RECOMMENDATIONS", "-" * 40]
    has_any = False
    for key, res in results.items():
        if res and getattr(res, "fix_commands", []):
            lines.append(f"  [{check_labels.get(key, key)}]")
            for cmd in res.fix_commands:
                lines.append(_indent(cmd, 4))
                lines.append("")
            has_any = True
    if not has_any:
        lines.append("  No fixes required — trajectory appears healthy.")
    lines.append("")

    lines += [
        "=" * 70,
        "  End of Report — TrajMedic v1.0",
        "=" * 70,
    ]

    content = "\n".join(lines)
    path = output_path or str(OUTPUT_DIR / "trajmedic_report.txt")
    _write(content, path)
    return content


def generate_html_report(
    metadata: Dict[str, Any],
    results: Dict[str, Any],
    health_score: float,
    health_label: str,
    plot_paths: Optional[Dict[str, str]] = None,
    output_path: Optional[str] = None,
) -> str:
    """
    Build and save a styled HTML diagnostic report with embedded plots.

    Parameters
    ----------
    metadata, results, health_score, health_label : (see generate_txt_report)
    plot_paths : dict, optional
        Mapping of 'rmsd' | 'rg' | 'com_drift' → absolute file path.
    output_path : str, optional

    Returns
    -------
    str
        Full HTML content.
    """
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    score_color = (
        "#2ecc71" if health_score >= 75 else
        "#f39c12" if health_score >= 50 else
        "#e74c3c"
    )

    check_labels = {
        "pbc":         "PBC / Fragmentation",
        "drift":       "COM Drift",
        "compactness": "Protein Compactness (Rg)",
        "rmsd":        "Backbone RMSD",
        "integrity":   "Structural Integrity",
    }

    # ── Build check rows ──────────────────────────────────────────────────────
    rows_html = ""
    for key, label in check_labels.items():
        res = results.get(key)
        if res is None:
            continue
        icon   = "✅" if res.passed else "❌"
        sev_cls = res.severity.lower()
        rows_html += f"""
        <tr class="sev-{sev_cls}">
          <td>{icon} {label}</td>
          <td>{res.severity}</td>
          <td>-{res.penalty:.0f}</td>
          <td>{res.details}</td>
        </tr>"""

    # ── Build recommendations ─────────────────────────────────────────────────
    recs_html = ""
    for key, res in results.items():
        if res and getattr(res, "fix_commands", []):
            recs_html += f"<h4>{check_labels.get(key, key)}</h4>"
            for cmd in res.fix_commands:
                recs_html += f"<pre><code>{cmd}</code></pre>"

    if not recs_html:
        recs_html = "<p>✅ No fixes required.</p>"

    # ── Embed plots ───────────────────────────────────────────────────────────
    plots_html = ""
    if plot_paths:
        for pname, ppath in plot_paths.items():
            if ppath and os.path.isfile(ppath):
                import base64
                with open(ppath, "rb") as f:
                    b64 = base64.b64encode(f.read()).decode()
                plots_html += (
                    f'<img src="data:image/png;base64,{b64}" '
                    f'style="max-width:100%;margin:10px 0;" /><br/>'
                )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <title>TrajMedic Report</title>
  <style>
    body {{ font-family: 'Segoe UI', Arial, sans-serif; background:#f5f7fa;
            color:#2c3e50; margin:0; padding:0; }}
    .container {{ max-width:900px; margin:40px auto; background:#fff;
                  border-radius:10px; box-shadow:0 2px 12px rgba(0,0,0,.1);
                  padding:40px; }}
    h1 {{ color:#1a252f; margin-bottom:4px; }}
    h2 {{ color:#2c3e50; border-bottom:2px solid #ecf0f1;
          padding-bottom:6px; margin-top:32px; }}
    h3 {{ color:#34495e; }}
    .badge {{
      display:inline-block; padding:10px 24px; border-radius:6px;
      font-size:2em; font-weight:700; color:#fff;
      background:{score_color}; margin:12px 0;
    }}
    table {{ width:100%; border-collapse:collapse; margin-top:12px; }}
    th {{ background:#ecf0f1; padding:8px 12px; text-align:left;
          font-size:0.85em; text-transform:uppercase; letter-spacing:.5px; }}
    td {{ padding:8px 12px; vertical-align:top; font-size:0.92em;
          border-bottom:1px solid #ecf0f1; }}
    tr.sev-none {{ background:#f0fff4; }}
    tr.sev-low  {{ background:#fffbea; }}
    tr.sev-moderate {{ background:#fff3e0; }}
    tr.sev-high {{ background:#fdecea; }}
    pre {{ background:#2c3e50; color:#ecf0f1; padding:14px; border-radius:6px;
           overflow-x:auto; font-size:0.85em; }}
    .meta-grid {{ display:grid; grid-template-columns:1fr 1fr;
                  gap:8px 24px; }}
    .meta-item {{ font-size:0.93em; }}
    .meta-label {{ color:#7f8c8d; font-size:0.8em; }}
    footer {{ text-align:center; color:#aaa; font-size:0.8em;
              margin-top:32px; }}
  </style>
</head>
<body>
<div class="container">
  <h1>🧬 TrajMedic Diagnostic Report</h1>
  <p style="color:#7f8c8d;">Generated: {ts}</p>

  <h2>System Summary</h2>
  <div class="meta-grid">
    <div class="meta-item">
      <div class="meta-label">Total Atoms</div>
      <strong>{metadata.get('n_atoms', 'N/A')}</strong>
    </div>
    <div class="meta-item">
      <div class="meta-label">Frames</div>
      <strong>{metadata.get('n_frames', 'N/A')}</strong>
    </div>
    <div class="meta-item">
      <div class="meta-label">Protein Residues</div>
      <strong>{metadata.get('n_residues', 'N/A')}</strong>
    </div>
    <div class="meta-item">
      <div class="meta-label">Duration</div>
      <strong>{_fmt(metadata.get('duration_ns'))} ns</strong>
    </div>
  </div>

  <h2>Trajectory Health Score</h2>
  <div class="badge">{health_score:.1f} / 100 — {health_label}</div>

  <h2>Diagnostic Results</h2>
  <table>
    <thead>
      <tr><th>Check</th><th>Severity</th><th>Penalty</th><th>Details</th></tr>
    </thead>
    <tbody>{rows_html}</tbody>
  </table>

  <h2>Fix Recommendations</h2>
  {recs_html}

  <h2>Visualisations</h2>
  {plots_html if plots_html else '<p>No plots available.</p>'}

  <footer>TrajMedic v1.0 — AI-Assisted Trajectory Diagnostics</footer>
</div>
</body>
</html>"""

    path = output_path or str(OUTPUT_DIR / "trajmedic_report.html")
    _write(html, path)
    return html


# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

def _fmt(value: Optional[float], decimals: int = 2) -> str:
    if value is None:
        return "N/A"
    return f"{value:.{decimals}f}"


def _indent(text: str, n: int) -> str:
    prefix = " " * n
    return "\n".join(prefix + line for line in text.splitlines())


def _write(content: str, path: str) -> None:
    try:
        with open(path, "w", encoding="utf-8") as f:
            f.write(content)
        logger.info(f"Report saved to {path}")
    except Exception as exc:
        logger.error(f"Failed to write report: {exc}")
