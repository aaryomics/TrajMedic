"""
dashboard.py — TrajMedic Streamlit Dashboard.

A modern, research-oriented UI that lets users upload trajectory files,
run the full diagnostic pipeline, and explore results interactively.

Launch with:
    streamlit run app/ui/dashboard.py
"""

from __future__ import annotations

import os
import sys
import logging
from pathlib import Path
from typing import Optional

# ── Ensure project root is on the path ───────────────────────────────────────
_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_ROOT))

import streamlit as st

logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# Page config (must be first Streamlit call)
# ──────────────────────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="TrajMedic",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ──────────────────────────────────────────────────────────────────────────────
# Custom CSS
# ──────────────────────────────────────────────────────────────────────────────

st.markdown("""
<style>
/* ── Global ── */
html, body, [class*="css"] { font-family: 'Inter', 'Segoe UI', sans-serif; }

/* ── Sidebar ── */
section[data-testid="stSidebar"] {
    background: #1a252f;
    color: #ecf0f1;
}
section[data-testid="stSidebar"] * { color: #ecf0f1 !important; }
section[data-testid="stSidebar"] h1,
section[data-testid="stSidebar"] h2,
section[data-testid="stSidebar"] h3 { color: #1abc9c !important; }

/* ── Health score card ── */
.health-card {
    border-radius: 12px;
    padding: 24px 32px;
    text-align: center;
    margin-bottom: 16px;
}
.health-card .score { font-size: 3.5em; font-weight: 800; }
.health-card .label { font-size: 1.2em; font-weight: 600; margin-top: 4px; }

/* ── Warning boxes ── */
.warn-high     { background:#fdecea; border-left:5px solid #e74c3c; padding:10px 14px; border-radius:4px; margin:6px 0; color:#1a1a1a !important; }
.warn-moderate { background:#fff3e0; border-left:5px solid #f39c12; padding:10px 14px; border-radius:4px; margin:6px 0; color:#1a1a1a !important; }
.warn-low      { background:#fffbea; border-left:5px solid #f1c40f; padding:10px 14px; border-radius:4px; margin:6px 0; color:#1a1a1a !important; }
.warn-none     { background:#f0fff4; border-left:5px solid #2ecc71; padding:10px 14px; border-radius:4px; margin:6px 0; color:#1a1a1a !important; }
.warn-high *, .warn-moderate *, .warn-low *, .warn-none * { color:#1a1a1a !important; }

/* ── Code blocks ── */
pre { background:#1e2a35 !important; color:#ecf0f1 !important;
      border-radius:6px; font-size:0.85em; }

/* ── Metric cards ── */
.metric-card {
    background:#f8f9fa; border-radius:8px; padding:14px 18px;
    border:1px solid #e0e0e0; text-align:center;
}
.metric-card .val { font-size:1.6em; font-weight:700; color:#2c3e50; }
.metric-card .lbl { font-size:0.78em; color:#7f8c8d; text-transform:uppercase;
                    letter-spacing:.5px; }
</style>
""", unsafe_allow_html=True)


# ──────────────────────────────────────────────────────────────────────────────
# Sidebar
# ──────────────────────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown("# 🧬 TrajMedic")
    st.markdown("**AI-Assisted Trajectory Diagnostics**")
    st.markdown("---")

    st.markdown("### 📂 Upload Files")
    # File size limits (enforced in code; global ceiling set in .streamlit/config.toml)
    _TPR_MAX_MB = 5       # TPR files are tiny — flag anything suspicious
    _XTC_MAX_GB = 4       # Practical ceiling for single-session analysis

    tpr_file = st.file_uploader(
        "Topology (.tpr)",
        type=["tpr"],
        help=f"GROMACS binary run input file (.tpr) — max {_TPR_MAX_MB} MB",
    )
    xtc_file = st.file_uploader(
        "Trajectory (.xtc)",
        type=["xtc"],
        help=f"GROMACS compressed trajectory file (.xtc) — max {_XTC_MAX_GB} GB",
    )

    # ── Per-file size validation ──────────────────────────────────────────
    _upload_ok = True

    if tpr_file is not None:
        tpr_mb = tpr_file.size / 1024 / 1024
        if tpr_mb > _TPR_MAX_MB:
            st.error(
                f"⚠️ TPR file is {tpr_mb:.1f} MB — expected < {_TPR_MAX_MB} MB. "
                "A standard .tpr topology is usually under 5 MB. "
                "Make sure you're not uploading a trajectory by mistake.",
                icon="🚨",
            )
            _upload_ok = False
        else:
            st.caption(f"✅ TPR: {tpr_mb:.2f} MB")

    if xtc_file is not None:
        xtc_mb  = xtc_file.size / 1024 / 1024
        xtc_gb  = xtc_mb / 1024
        if xtc_gb > _XTC_MAX_GB:
            st.error(
                f"⚠️ XTC file is {xtc_gb:.2f} GB — exceeds the {_XTC_MAX_GB} GB limit. "
                "Consider sub-sampling with: "
                "`gmx trjconv -skip N` to reduce file size.",
                icon="🚨",
            )
            _upload_ok = False
        elif xtc_mb < 1:
            st.caption(f"✅ XTC: {xtc_mb:.1f} MB")
        else:
            st.caption(f"✅ XTC: {xtc_gb:.2f} GB" if xtc_gb >= 1 else f"✅ XTC: {xtc_mb:.1f} MB")

    st.markdown("### ⚙️ Settings")
    max_frames = st.slider(
        "Max frames to analyse",
        min_value=50,
        max_value=2000,
        value=500,
        step=50,
        help=(
            "Fewer frames = faster analysis. "
            "Frames are sampled evenly across the trajectory."
        ),
    )

    st.markdown("---")
    run_btn = st.button("▶ Run Diagnostics", use_container_width=True, type="primary",
                        disabled=not _upload_ok)

    st.markdown("---")
    st.markdown(
        "<small>TrajMedic v1.0 · [GitHub](https://github.com) · MIT License</small>",
        unsafe_allow_html=True,
    )


# ──────────────────────────────────────────────────────────────────────────────
# Main content area
# ──────────────────────────────────────────────────────────────────────────────

st.title("🧬 TrajMedic — Trajectory Diagnostic Dashboard")
st.caption(
    "Upload your GROMACS trajectory files and run a complete health check "
    "to detect preprocessing issues, instabilities, and PBC artefacts."
)

# Welcome state
if not (tpr_file and xtc_file and _upload_ok):
    st.info(
        "👈 Upload a **.tpr** topology and **.xtc** trajectory in the sidebar, "
        "then click **Run Diagnostics**.",
        icon="ℹ️",
    )

    st.markdown("### What TrajMedic checks")
    cols = st.columns(3)
    checks_info = [
        ("🔗 PBC Fragmentation",
         "Detects protein split across periodic boundary conditions."),
        ("🧭 COM Drift",
         "Tracks center-of-mass displacement to find uncentered trajectories."),
        ("🔵 Compactness (Rg)",
         "Radius of gyration to detect unfolding or inflation artefacts."),
        ("📏 RMSD Stability",
         "Backbone RMSD to assess equilibration and structural jumps."),
        ("🏗️ Structural Integrity",
         "Identifies atom explosions and anomalous displacements."),
        ("💊 Health Score",
         "Single 0–100 score summarising all detected issues."),
    ]
    for i, (title, desc) in enumerate(checks_info):
        with cols[i % 3]:
            st.markdown(
                f'<div class="metric-card">'
                f'<div style="font-size:1.5em">{title.split()[0]}</div>'
                f'<div style="font-weight:600;margin:4px 0">{" ".join(title.split()[1:])}</div>'
                f'<div style="font-size:0.85em;color:#555">{desc}</div>'
                f'</div>',
                unsafe_allow_html=True,
            )
            st.markdown("")
    st.stop()

# ──────────────────────────────────────────────────────────────────────────────
# Run analysis when button clicked
# ──────────────────────────────────────────────────────────────────────────────

if run_btn:
    # ── Write uploaded files to a fixed session folder inside outputs/ ─────
    # We deliberately avoid TemporaryDirectory here because on Windows,
    # MDAnalysis keeps a file handle open on the .xtc until the Universe is
    # garbage-collected, which causes a PermissionError when Python tries to
    # delete a temp folder that still has an open handle.
    session_dir = Path(_ROOT) / "outputs" / "session"
    session_dir.mkdir(parents=True, exist_ok=True)

    tpr_path = str(session_dir / "input.tpr")
    xtc_path = str(session_dir / "input.xtc")

    with open(tpr_path, "wb") as f:
        f.write(tpr_file.read())
    with open(xtc_path, "wb") as f:
        f.write(xtc_file.read())

    # ── Progress bar ───────────────────────────────────────────────────────
    progress_bar = st.progress(0)
    status_text  = st.empty()

    def update_progress(step: int, total: int, message: str) -> None:
        progress_bar.progress(step / total)
        status_text.markdown(f"⏳ *{message}*")

    try:
        from app.main import run_full_analysis
        output = run_full_analysis(
            tpr_path=tpr_path,
            xtc_path=xtc_path,
            max_frames=max_frames,
            progress_callback=update_progress,
        )
    except Exception as exc:
        st.error(f"❌ Analysis failed: {exc}", icon="🚨")
        st.stop()

    progress_bar.progress(1.0)
    status_text.markdown("✅ *Analysis complete!*")

    # Cache output in session state so it persists across reruns
    st.session_state["output"] = output


# ──────────────────────────────────────────────────────────────────────────────
# Display results from session state
# ──────────────────────────────────────────────────────────────────────────────

if "output" not in st.session_state:
    st.stop()

output       = st.session_state["output"]
metadata     = output["metadata"]
results      = output["results"]
health_score = output["health_score"]
label        = output["health_label"]
color        = output["health_color"]
pfigs        = output["plotly_figs"]
txt_report   = output["txt_report"]
html_report  = output["html_report"]

st.markdown("---")

# ── Section 1: Health Score ────────────────────────────────────────────────
st.markdown("## 💊 Trajectory Health Score")

hs_col, meta_col = st.columns([1, 2])

with hs_col:
    from app.utils.plotting import plot_health_gauge
    st.plotly_chart(plot_health_gauge(health_score), use_container_width=True)

with meta_col:
    st.markdown("### System Metadata")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Atoms",   f"{metadata.get('n_atoms', 'N/A'):,}")
    m2.metric("Frames",  f"{metadata.get('n_frames', 'N/A'):,}")
    m3.metric("Residues",f"{metadata.get('n_residues', 'N/A')}")
    dur = metadata.get("duration_ns")
    m4.metric("Duration", f"{dur:.2f} ns" if dur else "N/A")

    box = metadata.get("box_dims")
    if box:
        st.caption(
            f"Box: {box[0]:.1f} × {box[1]:.1f} × {box[2]:.1f} Å  "
            f"({box[3]:.1f}° {box[4]:.1f}° {box[5]:.1f}°)"
        )

# ── Section 2: Warnings ───────────────────────────────────────────────────
st.markdown("---")
st.markdown("## ⚠️ Detected Issues")

check_labels = {
    "pbc":         ("🔗 PBC Fragmentation", "pbc_check"),
    "drift":       ("🧭 Center-of-Mass Drift", "drift_check"),
    "compactness": ("🔵 Protein Compactness (Rg)", "compactness_check"),
    "rmsd":        ("📏 Backbone RMSD Stability", "rmsd_check"),
    "integrity":   ("🏗️ Structural Integrity", "integrity_check"),
}

any_issues = False
for key, (label_txt, _) in check_labels.items():
    res = results.get(key)
    if res is None:
        continue

    sev  = res.severity.lower()
    icon = {"none": "✅", "low": "🟡", "moderate": "🟠", "high": "🔴"}.get(sev, "❓")
    css  = f"warn-{sev}"

    st.markdown(
        f'<div class="{css}"><strong>{icon} {label_txt}</strong> — '
        f'Severity: <em>{res.severity}</em> | Penalty: -{res.penalty:.0f} pts<br/>'
        f'<span style="font-size:0.9em">{res.details}</span></div>',
        unsafe_allow_html=True,
    )
    if not res.passed:
        any_issues = True

if not any_issues:
    st.success(
        "🎉 All checks passed! Your trajectory looks healthy.",
        icon="✅",
    )

# ── Section 3: Interactive Plots ──────────────────────────────────────────
st.markdown("---")
st.markdown("## 📊 Interactive Visualisations")

tab_rmsd, tab_rg, tab_drift = st.tabs(
    ["📏 RMSD", "🔵 Radius of Gyration", "🧭 COM Drift"]
)

with tab_rmsd:
    if "rmsd" in pfigs:
        rmsd_res = results["rmsd"]
        st.plotly_chart(pfigs["rmsd"], use_container_width=True)
        st.info(rmsd_res.interpretation)
        c1, c2, c3 = st.columns(3)
        c1.metric("Mean RMSD", f"{rmsd_res.mean_rmsd_ang:.2f} Å")
        c2.metric("Max RMSD",  f"{rmsd_res.max_rmsd_ang:.2f} Å")
        eq = rmsd_res.equilibration_time_ns
        c3.metric("Equilibration", f"{eq:.1f} ns" if eq > 0 else "—")
    else:
        st.warning("No RMSD data available.")

with tab_rg:
    if "rg" in pfigs:
        comp_res = results["compactness"]
        st.plotly_chart(pfigs["rg"], use_container_width=True)
        c1, c2, c3 = st.columns(3)
        c1.metric("Mean Rg",   f"{comp_res.mean_rg_ang:.2f} Å")
        c2.metric("Std Rg",    f"{comp_res.std_rg_ang:.2f} Å")
        c3.metric("Rg Spikes", str(comp_res.n_spike_frames))
    else:
        st.warning("No Rg data available.")

with tab_drift:
    if "com_drift" in pfigs:
        drift_res = results["drift"]
        st.plotly_chart(pfigs["com_drift"], use_container_width=True)
        c1, c2, c3 = st.columns(3)
        c1.metric("Max Drift",  f"{drift_res.max_drift_ang:.2f} Å")
        c2.metric("Mean Drift", f"{drift_res.mean_drift_ang:.2f} Å")
        c3.metric("Trend",      drift_res.drift_trend)
    else:
        st.warning("No COM drift data available.")

# ── Section 4: Recommendations ────────────────────────────────────────────
st.markdown("---")
st.markdown("## 🛠️ GROMACS Fix Recommendations")

has_any_commands = False
for key, (label_txt, _) in check_labels.items():
    res = results.get(key)
    if res and getattr(res, "fix_commands", []):
        with st.expander(f"{label_txt}", expanded=False):
            for cmd in res.fix_commands:
                st.code(cmd, language="bash")
        has_any_commands = True

if not has_any_commands:
    st.success(
        "No fixes required — trajectory is within acceptable parameters.",
        icon="✅",
    )

# ── Section 5: Per-check detail expanders ─────────────────────────────────
st.markdown("---")
st.markdown("## 🔬 Detailed Check Results")

for key, (label_txt, _) in check_labels.items():
    res = results.get(key)
    if res is None:
        continue

    icon = "✅" if res.passed else "❌"
    with st.expander(f"{icon} {label_txt} — Severity: {res.severity}", expanded=False):
        col_a, col_b = st.columns(2)
        col_a.markdown(f"**Status:** {'Pass' if res.passed else 'Fail'}")
        col_b.markdown(f"**Penalty:** -{res.penalty:.0f} pts")
        st.markdown(f"**Details:** {res.details}")

        # Show any extra numeric fields
        skip = {"passed", "penalty", "severity", "details",
                 "fix_commands", "times_ns", "rmsd_values_ang",
                 "rg_values_ang", "displacements_ang", "affected_frames",
                 "explosion_frames", "anomalous_frames", "jump_frame_indices",
                 "spike_frame_indices", "interpretation"}
        extras = {
            k: v for k, v in vars(res).items()
            if k not in skip and not isinstance(v, list)
        }
        if extras:
            st.json(extras)

# ── Section 6: Download reports ───────────────────────────────────────────
st.markdown("---")
st.markdown("## 📥 Download Reports")

dl1, dl2 = st.columns(2)
with dl1:
    st.download_button(
        label="📄 Download TXT Report",
        data=txt_report,
        file_name="trajmedic_report.txt",
        mime="text/plain",
        use_container_width=True,
    )
with dl2:
    st.download_button(
        label="🌐 Download HTML Report",
        data=html_report,
        file_name="trajmedic_report.html",
        mime="text/html",
        use_container_width=True,
    )
