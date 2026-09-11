# Faculty Discussion Brief: SMALT Stratigraphic Modeling Toolkit

**To:** Prof. Hiranya Sahoo  
**From:** Research Student  
**Date:** August 2026  
**Subject:** Completion of Phase 0 (Data Standardization & QC) and Phase 1 (1D Stratigraphic Markov Engine)  
**Target Reading Time:** 5 - 10 minutes  

---

## 1. High-Level Overview

Prof. Sahoo, the core objective of the **SMALT** (*Subsurface Stratigraphic Modeling & Active Learning Toolkit*) project is to build an open, verifiable, and mathematically grounded 2D/3D stratigraphic forward modeling and active learning framework. 

Our modeling architecture is directly anchored in field observations from the fluvio-deltaic Upper Cretaceous **Blackhawk Formation** (and Ferron Sandstone analogs). Rather than treating subsurface interpolation as a black-box machine learning problem, SMALT combines sedimentological rules (e.g., channel geometry, upward fining, avulsion) with probabilistic geostatistics.

We have successfully completed **Phase 0** and **Phase 1**, establishing a clean data pipeline and a verified 1D vertical Markov engine. Below is a quick walk-through of the code structure, the physical principles implemented, and our validation results.

---

## 2. Phase 0 - Data Standardization & Quality Control

**Core Module:** [`smalt/data/loader.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/data/loader.py) (aliased with [`data/loader.py`](file:///d:/Lithology-reconstruction-using-XGB/data/loader.py))  
**Class:** [`LithologLoader`](file:///d:/Lithology-reconstruction-using-XGB/data/loader.py#L63-L85)

### Unified Stratigraphic Schema
Raw field lithologs often come with inconsistent formats (some as continuous depth points, others as interval beds with top/bottom depths). I standardized all incoming data into a strict 6-column schema stored in parquet and CSV format at [`data/processed/lithologs_unified.parquet`](file:///d:/Lithology-reconstruction-using-XGB/data/processed/lithologs_unified.parquet):

| Column Name | Type | Physical Meaning / Description |
|---|---|---|
| `litholog_id` | String | Unique well / outcrop section identifier (e.g., `litholog1`, `litholog9`, `litholog11`) |
| `strike_pos_m` | Float | Spatial coordinate along strike in meters ($x$-position) |
| `depth_m` | Float | Vertical depth in meters ($z$-coordinate, 1-meter regularized interval) |
| `facies_code` | Integer | Standardized discrete state index ($0$ to $4$) |
| `facies_name` | String | Standardized facies label |
| `gamma_ray` | Float | Synthetic or measured Gamma Ray log response (API units) |

### 5-State Facies Aggregation
To ensure robust statistical counts across our outcrop sections, I grouped the field facies into 5 primary depositional states:
- **State 0 (Coal):** Mire / waterlogged peat swamp
- **State 1 (Channel Sandstone):** High-energy channel core and bar deposits (trough cross-stratified and planar sandstones)
- **State 2 (Fine Sandstone / Splay):** Crevasse splay and carbonaceous channel margin heterolithics
- **State 3 (Siltstone):** Waning flow and levee deposits
- **State 4 (Overbank Mudstone):** Floodplain suspension muds and shales

### Automated Quality Control (QC) Pipeline
Before any data reaches the mathematical engine, the [`LithologLoader.run_quality_control()`](file:///d:/Lithology-reconstruction-using-XGB/data/loader.py#L289-L330) method enforces three automated gatekeeper checks:
1. **Strict Depth Monotonicity:** Ensures depth values are strictly increasing without reversals or duplicate depth measurements ($z_{t+1} > z_t$).
2. **Zero-Thickness & Negative Interval Filtering:** Detects and flags erroneous intervals where $\text{Bottom} \le \text{Top}$ or depth $< 0$.
3. **Synthetic Gamma Ray (GR) Imputation:** Where real downhole GR logs are missing, the loader synthesizes realistic API responses based on physical shale content:
   - **Sandstone:** $\sim 30 - 35\text{ API}$ ($\text{Base } 35.0 + \mathcal{N}(0, 5)$)
   - **Siltstone:** $\sim 70 - 75\text{ API}$ ($\text{Base } 70.0 + \mathcal{N}(0, 5)$)
   - **Overbank Mudstone:** $\sim 105\text{ API}$ ($\text{Base } 105.0 + \mathcal{N}(0, 5)$)
   - **Coal:** $\sim 20\text{ API}$ (low natural radioactivity)

---

## 3. Phase 1 - 1D Stratigraphic Markov Engine

**Core Modules:** [`smalt/geostat/markov.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/geostat/markov.py) and [`smalt/viz/markov_viz.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/viz/markov_viz.py)  
**Class:** [`StratigraphicMarkovChain`](file:///d:/Lithology-reconstruction-using-XGB/smalt/geostat/markov.py#L35-L104)

Phase 1 captures the vertical facies transition tendencies from the digitized lithologs. Here is what each key method does in plain language:

### 1. `fit()` - Tallying Vertical Upward Successions
- **Sedimentological Logic:** Geological strata are deposited sequentially from bottom to top. The [`fit()`](file:///d:/Lithology-reconstruction-using-XGB/smalt/geostat/markov.py#L106-L240) method sorts each well from base to top (decreasing depth $z$) and tallies every step from the underlying bed $S_{t-1}$ to the overlying bed $S_t$.
- **Two Modeling Modes:**
  - *Regular Chain ($P_{\text{reg}}$):* 1-meter fixed steps, preserving self-transitions ($P_{ii} > 0$) to capture bed thicknesses.
  - *Embedded Chain ($P_{\text{emb}}$):* Isolates true lithologic boundary crossings ($P_{ii} = 0$) to analyze purely genetic facies transitions.
- **Physical Result:** In our outcrop data, when a Channel Sandstone transitions across a bed boundary, it transitions $100\%$ of the time into finer facies ($40\%$ Fine Sand/Splay, $30\%$ Coal/Marsh, $20\%$ Siltstone, $10\%$ Mudstone), successfully quantifying the classic **fining-upward fluvial cycle**.

### 2. `compute_stationary_distribution()` - Long-Run Facies Proportions
- **What it does:** Solves the invariant left-eigenvector equation ($\boldsymbol{\pi} \mathbf{P} = \boldsymbol{\pi}$, $\sum \pi_i = 1.0$) to determine the theoretical long-run facies volume fractions if deposition were to continue indefinitely under these transition probabilities.
- **Net-to-Gross Proxy:** This serves as a vital diagnostic cross-check against our empirical field measurements:
  - *Theoretical Stationary Net-to-Gross (Sand + Splay):* **52.51%**
  - *Observed Empirical Net-to-Gross:* **55.06%**
  - *Relative Difference:* Only **4.63%**, confirming that the Markov transition matrix accurately preserves the macroscopic sand budget.

### 3. `compute_asymmetry_matrix()` & `compute_state_entropy()` - Directionality & Predictability
- **Asymmetry Matrix ($\mathbf{A} = \mathbf{P} - \mathbf{P}^T$):** Quantifies whether a transition is directionally favored ($A_{ij} > 0$) versus its reverse downward counterpart. Positive values for Channel Sand $\to$ Splay/Mud isolate fining-upward channel abandonment cycles from coarsening-upward splay progradation.
- **State Transition Entropy ($H_i$):** Measures the geological disorder and predictability of transitions out of facies $i$. Channel Sandstone shows low normalized entropy ($H_1^{\text{norm}} = 0.28$ in regular chain), reflecting its predictable fining-upward behavior, whereas mudstone/floodplain shows higher entropy due to variable avulsion and splay events.

### 4. `plot_transition_matrix()` - Publication-Quality Visualization
- Uses [`smalt/viz/markov_viz.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/viz/markov_viz.py) to generate annotated 300-DPI heatmaps showing transition probabilities and transition observation counts ($n$) per cell.
- Diagnostic figures are saved in [`docs/figures/`](file:///d:/Lithology-reconstruction-using-XGB/docs/figures):
  - Transition Heatmaps: `phase1_transition_matrix_regular.png` & `phase1_transition_matrix_embedded.png`
  - Validation Proportions: `phase1_stationary_vs_empirical.png`
  - Directional Network: `phase1_facies_succession_network.png`

---

## 4. Test Suite & Verification

**Test Directory:** [`tests/`](file:///d:/Lithology-reconstruction-using-XGB/tests)  
**Execution:** `python -m pytest tests/ -v`

We have implemented an automated unit and property-based test suite with **21 passing tests (100% pass rate)**:

1. **Phase 0 Tests (10/10 passed - [`tests/test_phase0_loader.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_phase0_loader.py)):**
   - Verified strict column schema integrity and zero null values.
   - Tested automatic catching of negative depths, zero thickness, and inverted top/bottom layers.
   - Verified synthetic GR distribution matches specified Gaussian noise models ($\mu \pm 5\text{ API}$).
   - Validated end-to-end multi-well export to Parquet and CSV.

2. **Phase 1 Tests (11/11 passed - [`tests/test_phase1_markov.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_phase1_markov.py)):**
   - **Row-Stochasticity Invariant:** Every matrix row strictly satisfies $\sum_{j} P_{ij} = 1.0$ within numerical tolerance of $10^{-7}$.
   - **Non-Negativity:** All transition probabilities satisfy $0.0 \le P_{ij} \le 1.0$.
   - **Embedded Chain Invariant:** Diagonal elements are strictly zero ($\text{diag}(P_{\text{emb}}) = 0.0$).
   - **Synthetic Sequence Recovery:** Confirmed that feeding a synthetic Sand $\to$ Silt $\to$ Mud sequence accurately yields directional probabilities of $1.0$ in upward order.
   - **Eigensolver Convergence:** Stationary distribution satisfies $\|\boldsymbol{\pi}\mathbf{P} - \boldsymbol{\pi}\|_{\infty} < 10^{-6}$.

---

## 5. Next Immediate Step: Transition to Phase 2

With our vertical 1D succession statistics mathematically verified, we are now ready to commence **Phase 2 (2D Cross-Sectional Geostatistical Modeling)**. 

The vertical transition probability matrix will directly seed the vertical stacking rules of a 2D object-based fluvial simulator, where channel bodies are populated using realistic Blackhawk channel geometry parameters (width-to-thickness aspect ratios $W/T \approx 35$) conditioned on our outcrop lithologs.

---

### Quick Reference Links
- Data Loader: [`smalt/data/loader.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/data/loader.py)
- Markov Engine: [`smalt/geostat/markov.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/geostat/markov.py)
- Visualizations: [`smalt/viz/markov_viz.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/viz/markov_viz.py)
- Test Suite: [`tests/`](file:///d:/Lithology-reconstruction-using-XGB/tests/)
- Comprehensive Technical Note: [`docs/notes/phase1_markov.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/notes/phase1_markov.md)
- Complete Numerical Summary: [`results/phase1_markov_summary.json`](file:///d:/Lithology-reconstruction-using-XGB/results/phase1_markov_summary.json)
