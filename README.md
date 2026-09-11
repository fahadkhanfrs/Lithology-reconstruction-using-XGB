# SMALT: Subsurface Stratigraphic Modeling & Active Learning Toolkit

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-21%20passed%20%7C%20100%25-brightgreen.svg)](tests/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**SMALT** is an open, verifiable, and mathematically grounded Python toolkit for 1D/2D/3D stratigraphic forward modeling, vertical facies succession analysis, and contextual active learning. The framework is calibrated against fluvio-deltaic outcrop analogs from the Upper Cretaceous **Blackhawk Formation** (and Ferron Sandstone) to bridge classical sedimentological principles with probabilistic subsurface characterization.

---

## Key Features

- **Standardized Stratigraphic Data Pipeline (Phase 0):** Automated parsing, quality control (QC), depth-monotonicity validation, zero-thickness layer filtering, and synthetic Gamma Ray imputation ($30-35\text{ API}$ sand, $70-75\text{ API}$ silt, $105\text{ API}$ mud).
- **1D Vertical Stratigraphic Markov Engine (Phase 1):** Rigorous discrete-step Markov modeling supporting both:
  - *Regular Chains ($P_{\text{reg}}$):* Fixed-step transitions preserving lithologic layer thickness and unit persistence.
  - *Embedded Chains ($P_{\text{emb}}$):* Pure boundary-crossing transitions ($P_{ii} = 0$) isolating genetic depositional cycles.
- **Sedimentological Invariance & Net-to-Gross Diagnostics:** Invariant left-eigensolver for stationary facies occupancy ($\boldsymbol{\pi}\mathbf{P} = \boldsymbol{\pi}$) matching empirical field Net-to-Gross within $4.63\%$.
- **Cyclicity & Predictability Metrics:** Directional asymmetry tensors ($\mathbf{A} = \mathbf{P} - \mathbf{P}^T$) and normalized state transition entropy ($H_i$) capturing upward-fining channel fills vs. coarsening-upward splays.
- **Publication-Ready Visualization:** Annotated 300-DPI heatmaps, transition succession networks, and empirical validation plots.
- **Context-Aware Active Learning & XGBoost Reconstruction:** Adaptive uncertainty sampling and hybrid sampling strategies for cost-effective core/well-log data acquisition.

---

## Project Structure

```text
Lithology-reconstruction-using-XGB/
├── data/
│   ├── raw_lithologs/             # Raw field outcrop litholog CSVs
│   ├── processed/                 # Unified parquet/csv datasets
│   └── loader.py                  # Phase 0 Data standardization & QC engine
├── docs/
│   ├── figures/                   # 300-DPI publication figures & heatmaps
│   ├── notes/                     # In-depth technical derivations (Phase 1)
│   └── FACULTY_DISCUSSION_BRIEF.md# Executive briefing for faculty / stakeholders
├── results/
│   └── phase1_markov_summary.json # Complete numerical transition tensors & metrics
├── scripts/
│   └── run_phase1_markov.py       # Markov pipeline runner and figure generator
├── smalt/
│   ├── __init__.py
│   ├── data/                      # Data subpackage (aliased loader)
│   ├── geostat/
│   │   └── markov.py              # 1D Stratigraphic Markov Chain core engine
│   └── viz/
│       └── markov_viz.py          # Heatmaps, network graphs, and diagnostics
├── src/                           # Machine learning & active learning modules
│   ├── litholog_model.py          # Active learning & XGBoost facies classifier
│   └── ...
├── tests/
│   ├── test_phase0_loader.py      # Unit & QC pipeline tests (10/10 passing)
│   └── test_phase1_markov.py      # Mathematical invariant & Markov tests (11/11 passing)
├── HANDOFF.md                     # Phase status & verification registry
├── main.py                        # ML experiment driver script
└── README.md
```

---

## Standardized Facies Schema (5-State SMALT)

Field lithologies are aggregated into a standardized 5-state discrete space:

| Code | Facies Label | Depositional Sub-Environment | Typical Base GR |
|:---:|:---|:---|:---:|
| `0` | **Coal** | Waterlogged peat mire / swamp | $20\text{ API}$ |
| `1` | **Channel Sandstone** | High-energy channel core & dunes | $35\text{ API}$ |
| `2` | **Fine Sandstone / Splay** | Crevasse splay & carbonaceous heterolithics | $130\text{ API}$ (carbon mud) |
| `3` | **Siltstone** | Waning flow, levee & channel margin | $70\text{ API}$ |
| `4` | **Overbank Mudstone** | Floodplain shale & suspension mud | $105\text{ API}$ |

---

## Quick Start

### 1. Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/fahadkhanfrs/Lithology-reconstruction-using-XGB.git
cd Lithology-reconstruction-using-XGB
pip install -r requirements.txt
```

*(Core dependencies: `numpy`, `pandas`, `scipy`, `matplotlib`, `pyarrow`, `scikit-learn`, `xgboost`, `pytest`)*

### 2. Run the Phase 0 & Phase 1 Markov Pipeline

Process the raw lithologs, fit vertical transition matrices, compute stationary distributions, and generate diagnostic figures:

```bash
python scripts/run_phase1_markov.py
```

### 3. Run Active Learning & ML Experiments

Execute baseline reconstruction, uniform sampling, and hybrid uncertainty-based sampling:

```bash
python main.py
```

### 4. Run Automated Test Suite

Execute the full suite of 21 unit and property-based verification tests:

```bash
python -m pytest tests/ -v
```

---

## Summary of Verified Results

### 1D Stratigraphic Markov Succession (Blackhawk Outcrops)

- **Total Validated Observations:** 247 regularized depth points across 3 wells (`litholog1`, `litholog9`, `litholog11`).
- **Channel Sandstone Transition Asymmetry:** $100\%$ of boundary transitions from Channel Sandstone ($P_{\text{emb}}$) proceed upward into finer-grained facies ($40\%$ Fine Sand/Splay, $30\%$ Coal, $20\%$ Siltstone, $10\%$ Overbank Mudstone), quantitatively validating the upward-fining channel-fill model.
- **Macroscopic Net-to-Gross Agreement:**
  - Stationary Theoretical ($Sand + Splay$): **$52.51\%$**
  - Empirical Observed ($Sand + Splay$): **$55.06\%$**
  - Relative Difference: **$4.63\%$** ($< 5\%$ tolerance).

---

## Documentation

- **Executive Brief:** [`docs/FACULTY_DISCUSSION_BRIEF.md`](docs/FACULTY_DISCUSSION_BRIEF.md) - High-level walk-through for meetings and reviews.
- **Technical Note:** [`docs/notes/phase1_markov.md`](docs/notes/phase1_markov.md) - Full mathematical derivations and sedimentological references.
- **Handoff Registry:** [`HANDOFF.md`](HANDOFF.md) - Exact numerical matrices, condition numbers, and test acceptance criteria.

---

## References

1. **Krumbein, W. C., & Dacey, M. F. (1969).** Markov chains and embedded Markov chains in geology. *Mathematical Geology*, 1(1), 79-96.
2. **Doveton, J. H. (1971).** An application of Markov chain analysis to the Ayrshire Coal Measures succession. *Scottish Journal of Geology*, 7(1), 11-27.
3. **Sahoo, H., Gani, M. R., & Gani, N. D. (2016).** 3D facies architecture and sequence stratigraphy of a fluvio-deltaic succession, Cretaceous Ferron Sandstone, Utah. *Sedimentology*, 63(6), 1403-1437.
