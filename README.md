# 🧬 TrajMedic

> **AI-Assisted Trajectory Diagnostics for GROMACS Molecular Dynamics Simulations**

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org)
[![MDAnalysis](https://img.shields.io/badge/MDAnalysis-2.6+-green.svg)](https://www.mdanalysis.org)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.28+-red.svg)](https://streamlit.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

TrajMedic is a modular, research-grade Python tool for diagnosing common
preprocessing problems and trajectory integrity issues in GROMACS molecular
dynamics simulations.  It produces a **Trajectory Health Score**, actionable
**GROMACS fix commands**, publication-quality **plots**, and a downloadable
**diagnostic report** — all accessible through a polished **Streamlit dashboard**.

---

## ✨ Features

| Module | What it detects |
|--------|----------------|
| 🔗 **PBC Check** | Protein fragmentation due to improper periodic boundary handling |
| 🧭 **COM Drift** | Center-of-mass translation indicating missing `-center` step |
| 🔵 **Compactness (Rg)** | Unfolding events and artificial Rg inflation |
| 📏 **RMSD Stability** | Backbone RMSD jumps, equilibration detection |
| 🏗️ **Integrity Check** | Atom explosions, numerical instabilities, corrupted frames |
| 💊 **Health Score** | Aggregated 0–100 trajectory quality metric |
| 🛠️ **Recommendations** | Exact `gmx trjconv` commands to fix each detected issue |
| 📊 **Interactive Plots** | Plotly + Matplotlib charts (RMSD, Rg, COM drift) |
| 📄 **Reports** | Downloadable TXT and styled HTML reports |

---

## 📸 Screenshots

> _Add dashboard screenshots here once running_

```
[Dashboard overview]          [Health score gauge]
[RMSD interactive plot]       [Recommendation panel]
```

---

## 🧪 Scientific Motivation

Even experienced researchers routinely submit trajectories to analysis
pipelines without proper preprocessing.  Common mistakes include:

- Forgetting `gmx trjconv -pbc nojump` → protein split across box boundary
- Omitting `-center` → COM drifts across the box
- Skipping `-fit rot+trans` → artificial RMSD inflation
- Inadequate energy minimisation → explosive numerical instabilities

TrajMedic acts as a **first-pass sanity check** that catches these issues
before they invalidate downstream analyses (free-energy calculations, RMSF,
principal component analysis, etc.).

---

## ⚡ Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/your-username/trajmedic.git
cd trajmedic
```

### 2. Create and activate a virtual environment

```bash
python -m venv .venv
source .venv/bin/activate          # Linux / macOS
# .venv\Scripts\activate           # Windows
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

---

## 🚀 Usage

### Option A — Streamlit Dashboard (recommended)

```bash
streamlit run app/ui/dashboard.py
```

Open your browser at `http://localhost:8501`, upload your `.tpr` and `.xtc`
files in the sidebar, and click **Run Diagnostics**.

---

### Option B — Command-Line Interface

```bash
python -m app.main --tpr examples/sample.tpr --xtc examples/traj.xtc
```

**Options:**

```
--tpr PATH          Path to .tpr topology file   (required)
--xtc PATH          Path to .xtc trajectory file (required)
--max-frames INT    Max frames to sample [default: 500]
--no-html           Skip HTML report generation
```

**Example output:**

```
============================================================
  Trajectory Health Score: 68.0/100  [MODERATE]
============================================================
  [FAIL] pbc          severity=High    penalty=-30
  [PASS] drift        severity=None    penalty=0
  [FAIL] compactness  severity=Low     penalty=-8
  [PASS] rmsd         severity=None    penalty=0
  [PASS] integrity    severity=None    penalty=0

  Reports saved to outputs/
  Plots  saved to outputs/*.png
============================================================
```

---

### Option C — Python API

```python
from app.main import run_full_analysis

output = run_full_analysis(
    tpr_path="examples/sample.tpr",
    xtc_path="examples/traj.xtc",
    max_frames=500,
)

print(f"Health score: {output['health_score']:.1f}/100")
print(f"Label:        {output['health_label']}")

# Access individual check results
pbc_result = output["results"]["pbc"]
print(f"PBC severity: {pbc_result.severity}")
print(f"Details:      {pbc_result.details}")
```

---

## 🗂️ Project Structure

```
trajmedic/
│
├── app/
│   ├── main.py                  # Pipeline orchestrator & CLI entry point
│   │
│   ├── checks/
│   │   ├── pbc_check.py         # PBC fragmentation detector
│   │   ├── drift_check.py       # Center-of-mass drift analyser
│   │   ├── rmsd_check.py        # Backbone RMSD stability
│   │   ├── compactness_check.py # Radius of gyration / unfolding
│   │   └── integrity_check.py   # Explosion & anomaly detector
│   │
│   ├── utils/
│   │   ├── helpers.py           # Universe loading, sampling, scoring
│   │   ├── plotting.py          # Matplotlib + Plotly figures
│   │   └── report.py            # TXT and HTML report generator
│   │
│   └── ui/
│       └── dashboard.py         # Streamlit dashboard
│
├── examples/                    # Sample .tpr / .xtc files
├── outputs/                     # Generated plots and reports
├── requirements.txt
├── README.md
└── LICENSE
```

---

## 🏥 Health Score System

| Score | Label    | Meaning |
|-------|----------|---------|
| 75–100 | 🟢 Good     | Trajectory is suitable for analysis |
| 50–74  | 🟡 Moderate | Issues present; fix before key analyses |
| 0–49   | 🔴 Poor     | Significant problems; preprocessing required |

**Penalty table:**

| Issue | Penalty |
|-------|---------|
| PBC fragmentation (High) | -30 pts |
| COM drift (High)         | -20 pts |
| RMSD instability (High)  | -25 pts |
| Rg anomaly (High)        | -15 pts |
| Integrity failure (High) | -20 pts |

---

## 🔬 Using MDAnalysis Test Files

If you don't have GROMACS data handy, use the bundled MDAnalysis test files:

```python
import MDAnalysis.tests.datafiles as mda_data
from app.main import run_full_analysis

output = run_full_analysis(
    tpr_path=mda_data.TPR,
    xtc_path=mda_data.XTC,
)
```

---

## 🗺️ Roadmap

- [ ] **Automatic trajectory fixing** — apply `gmx trjconv` commands directly from the dashboard
- [ ] **Batch processing** — analyse entire directories of trajectories
- [ ] **Multi-chain / antibody support** — handle Fv fragments and multi-domain systems
- [ ] **AI-generated interpretations** — LLM-powered natural language explanations of results
- [ ] **AMBER / NAMD / OpenMM support** — via MDAnalysis universal topology readers
- [ ] **RMSF and per-residue plots** — residue-level flexibility visualisation
- [ ] **Secondary structure tracking** — DSSP-based unfolding detection
- [ ] **Comparative analysis** — overlay multiple trajectories on one dashboard
- [ ] **CI integration** — run TrajMedic as part of a simulation quality-control pipeline

---

## 🤝 Contributing

Contributions are welcome!  Please open an issue to discuss proposed changes
before submitting a pull request.

```bash
# Run linting
pip install ruff
ruff check app/

# Run tests (once test suite is added)
pytest tests/
```

---

## 📚 Dependencies

| Package | Purpose |
|---------|---------|
| [MDAnalysis](https://www.mdanalysis.org) | Trajectory parsing and atom selections |
| [NumPy](https://numpy.org) | Numerical computations |
| [Pandas](https://pandas.pydata.org) | Data manipulation |
| [Matplotlib](https://matplotlib.org) | Publication-quality static plots |
| [Plotly](https://plotly.com) | Interactive dashboard charts |
| [Streamlit](https://streamlit.io) | Web-based dashboard UI |
| [SciPy](https://scipy.org) | Statistical utilities |

---

## 📄 License

This project is licensed under the **MIT License** — see [LICENSE](LICENSE) for details.

---

## 🙏 Acknowledgements

- [MDAnalysis](https://www.mdanalysis.org) development team
- [GROMACS](https://www.gromacs.org) community
- The structural biology community for documenting these preprocessing pitfalls

---

*TrajMedic is a research utility, not a substitute for expert judgment.  
Always validate results against domain knowledge and experimental data.*
