# SMALT Project Handoff & Phase Registry

## Document Authority & Repository Metadata
- **Project Name**: Subsurface Stratigraphic Modeling & Active Learning Toolkit (SMALT)
- **Geological Setting**: Fluvial-deltaic Upper Cretaceous Blackhawk Formation (Book Cliffs, Utah; Sahoo et al., 2016 analog)
- **Active Git Branch**: `sprint-j` (Head Commit: `f08587c`)
- **Current Milestone**: Sprint J Complete & Fully Verified; Phase 2 (2D Fluvial Forward Modeling) Ready for Sprint
- **Repository Location**: [`d:/Lithology-reconstruction-using-XGB`](file:///d:/Lithology-reconstruction-using-XGB)
- **Core Authority Files**:
  - Unified Stratigraphic Dataset: [`data/processed/lithologs_unified.parquet`](file:///d:/Lithology-reconstruction-using-XGB/data/processed/lithologs_unified.parquet) / [`data/processed/lithologs_unified.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/processed/lithologs_unified.csv)
  - Provenance Registry: [`data/provenance_manifest.json`](file:///d:/Lithology-reconstruction-using-XGB/data/provenance_manifest.json)
  - Coordinate Loader: [`smalt/spatial/coordinates.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/coordinates.py)
  - Datum Transformation: [`smalt/spatial/datum.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/datum.py)
  - Continuous Markov Geostatistics: [`smalt/geostat/spatial_transition.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/geostat/spatial_transition.py)
  - Empirical Horizontal Calibration: [`smalt/geostat/empirical_calibration.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/geostat/empirical_calibration.py)
  - Master Execution Pipelines: [`scripts/run_sprint_f_pipeline.py`](file:///d:/Lithology-reconstruction-using-XGB/scripts/run_sprint_f_pipeline.py), [`scripts/run_sprint_h_pipeline.py`](file:///d:/Lithology-reconstruction-using-XGB/scripts/run_sprint_h_pipeline.py), [`scripts/run_sprint_i_pipeline.py`](file:///d:/Lithology-reconstruction-using-XGB/scripts/run_sprint_i_pipeline.py), [`scripts/run_sprint_j_pipeline.py`](file:///d:/Lithology-reconstruction-using-XGB/scripts/run_sprint_j_pipeline.py)

---

## 1. Phase Status Overview

| Phase / Sprint | Focus & Scope | Status | Verified Metrics / Scope |
| :--- | :--- | :---: | :--- |
| **Phase 0** | Data Ingestion & Quality Control | **COMPLETE / VERIFIED** | 11 lithologs loaded cleanly; 920 points; monotonic depth QC; synthetic GR proxy; 11/11 tests passing. |
| **Phase 1** | 1D Vertical Markov Succession Analysis | **COMPLETE / VERIFIED** | Walther's Law upward fining (100% from sand); regular $P_{\text{reg}}$ and embedded $P_{\text{emb}}$ matrices; 11/11 tests passing. |
| **Sprint C** | Descriptive Lithology & Leakage-Free 1D Baseline | **COMPLETE / VERIFIED** | All 12 lithologs audited; descriptive statistics; 12/12 unit tests passing. |
| **Sprint D** | Reproducible 1D Markov Validation & Reconciliation | **COMPLETE / VERIFIED** | Reconciled 72 raw sand intervals (444.70 m) and 53 merged bodies; separated initial-state vs transition score; 10/10 tests passing. |
| **Sprint E** | Revised Dataset Reconciliation & Spatial Readiness | **COMPLETE / SUPERSEDED** | Reconciled 12 revised logs; L12 registered in manifest; synthetic-coordinate audit; real Cartesian coordinate module created; 11/11 tests passing. |
| **Sprint F** | Six-State Migration & Spatial Baseline Revalidation | **COMPLETE / VERIFIED** | Migrated to canonical 6-state schema; recomputed 6x6 Markov transitions; Spatial LOLO on 940 pts; 10/10 tests passing. |
| **Sprint G** | Independent Results Audit & Spatial Modeling Decision | **COMPLETE / VERIFIED** | Independent audit (0 discrepancy); identified StandardScaler anisotropy cancellation; established physical inter-well sparsity root cause; 8/8 tests passing. |
| **Sprint H** | Common-Zero Datum Alignment & Orientation Correction | **COMPLETE / VERIFIED** | Inverted top-down measured sections to bottom-up $z_{\text{strat}} = H_{\max} - d$; preserved L12 drill core; verified 100% stratigraphic invariance (1033.0 m, 328 beds); rebuilt 4,410 horizontal pairs; 9/9 tests passing. |
| **Sprint I** | Conditioned Transition-Probability Geostatistical Prototype | **COMPLETE / VERIFIED** | Coupled vertical $P_v$ and horizontal continuous rate $P_h(h) = \exp(R_h h)$; guaranteed 100% hard conditioning honor rate; stochastic realizations reproduce bed thickness (3.85 m vs 3.76 m); proved pointwise MAP collapse; 15/15 tests passing. |
| **Sprint J** | Empirical Horizontal Transition & Lateral Continuity Calibration | **COMPLETE / VERIFIED** | Extracted 4,410 pairs across 93 common-zero slices; calibrated facies lengths (sand $L=402.7$ m, mud $L=343.4$ m); tested alternative decay models (Spherical decay reduces MSE by 63.2%); 4-model ablation proved bed thickness is driven 100% by vertical succession; Category C classification; 15/15 tests passing. |
| **Phase 2** | 2D Fluvial Forward Modeling (Object-Based) | **READY FOR SPRINT** | Scoped strictly to unconditioned process-based forward simulation with Sahoo et al. (2016) geometric priors ($W/T = 35$, mean thickness $5.8$ m, target NTG 17-50%). |
| **Phase 3** | Inter-Well Spatial Reconstruction on Real Wells | **DEFERRED / REJECTED** | Definitively proven impossible on real data due to inter-well spacing (420-2356 m) exceeding geobody dimensions (140-210 m). Documented as rigorous negative result. |
| **Phases 4-5** | Spatial Machine Learning & Active Learning Querying | **RE-SCOPED** | Re-scoped to synthetic 2D forward realizations (from Phase 2) where true ground truth is known, proving borehole budget optimization without making false field claims. |

---

## 2. Current Project State & Architecture

### 2.1 Repository Architecture & Codebase Structure
```
smalt/
├── data/
│   ├── loader.py                     # Unified pandas/parquet loader & QC pipeline
│   ├── provenance.py                 # Manifest verification & schema validators
├── descriptive/
│   └── analyzer.py                   # Canonical facies schema, thickness & N/G profiler
├── spatial/
│   ├── coordinates.py                # Real Cartesian (X, Y) loader (L2-L12)
│   ├── datum.py                      # Common-zero datum transformation (z_strat = H_max - d)
│   └── baseline.py                   # Spatial LOLO CV engine (KNN, Nearest-Well, Priors)
├── geostat/
│   ├── markov.py                     # 1D StratigraphicMarkovChain (Walther's Law, regular/embedded)
│   ├── spatial_transition.py         # Continuous rate matrix R_h, Carle & Fogg formulation
│   ├── conditioned_markov.py         # Coupled vertical-horizontal conditioned classifier
│   ├── empirical_calibration.py      # HorizontalContinuityCalibrator (Sprint J calibration)
│   ├── ablation.py                   # 4-model ablation engine & bed thickness metrics
│   ├── realization.py                # 2D cross-sectional stochastic realization generator
│   └── object_sim.py                 # [Planned Phase 2] 2D ribbon forward simulator
├── validation/
│   ├── sprint_h.py                   # Sprint H validation suite & LOLO rebuilder
│   ├── sprint_i.py                   # Sprint I validation suite & synthetic benchmark
│   └── sprint_j.py                   # Sprint J validation suite & synthetic length recovery
├── viz/
│   ├── markov_viz.py                 # Transition matrix heatmaps & succession networks
│   └── cross_sections.py             # 2D cross-section and entropy map renderers
scripts/
│   ├── run_sprint_f_pipeline.py      # Master reproduction script for Sprint F
│   ├── run_sprint_h_pipeline.py      # Master reproduction script for Sprint H
│   ├── run_sprint_i_pipeline.py      # Master reproduction script for Sprint I
│   └── run_sprint_j_pipeline.py      # Master reproduction script for Sprint J
tests/
│   ├── test_phase0_loader.py         # 11 tests: schema integrity & QC
│   ├── test_phase1_markov.py         # 11 tests: 1D regular & embedded Markov properties
│   ├── test_sprint_c.py              # 12 tests: 12-litholog descriptive audit
│   ├── test_sprint_d.py              # 10 tests: sandstone bed stats & LOLO scoring
│   ├── test_sprint_e.py              # 11 tests: coordinate audit & L12 registration
│   ├── test_sprint_f.py              # 10 tests: 6-state schema & Spatial LOLO baseline
│   ├── test_sprint_g.py              # 8 tests: independent metrics audit & anisotropy defect
│   ├── test_sprint_h.py              # 9 tests: orientation audit, invariance & datum alignment
│   ├── test_sprint_i.py              # 15 tests: continuous rate, 100% honor rate, MAP collapse
│   └── test_sprint_j.py              # 15 tests: empirical pairs, calibration, ablation, recovery
sprints/                              # Audit directories with verified CSVs, reports, and figures
```

### 2.2 Ingested Working Dataset State
- **Authoritative Files**: [`data/processed/lithologs_unified.parquet`](file:///d:/Lithology-reconstruction-using-XGB/data/processed/lithologs_unified.parquet) and [`data/processed/lithologs_unified.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/processed/lithologs_unified.csv)
- **Total Working Lithologs**: 12 profiles
  - Lithologs 1-11: Outcrop measured sections from Sahoo et al. (2016)
  - Litholog 12: Subsurface drill core graphic log (EM-137C, digitized interval 0.0 - 111.0 m)
- **Total Validated Stratigraphic Span**: Exactly 1,033.0 cumulative meters across 328 discrete lithological bed intervals.
- **Regularized 1-Meter Evaluation Points**: Exactly 1,031 points across the 12 logs (due to two 1m recording gaps: Litholog 9 at 18-19m and Litholog 11 at 59-60m).
- **Spatial Sub-Dataset (L2-L12)**: Exactly 940 evaluation points across 11 wells with verified local Cartesian coordinates.

#### Well Breakdown Summary:
| Litholog ID | Data Type | Stratigraphic Span | Discretized Points | Coordinate Status | Orientation Transform |
| :--- | :--- | :---: | :---: | :--- | :--- |
| `litholog1` | Outcrop Measured Section | 93.0 m | 93 pts | Missing (Excluded from Spatial) | $z_{\text{strat}} = 93.0 - d$ |
| `litholog2` | Outcrop Measured Section | 93.0 m | 93 pts | Verified $(X=0.0, Y=0.0)$ | $z_{\text{strat}} = 93.0 - d$ |
| `litholog3` | Outcrop Measured Section | 93.0 m | 93 pts | Verified $(X=1454.0, Y=580.0)$ | $z_{\text{strat}} = 93.0 - d$ |
| `litholog4` | Outcrop Measured Section | 84.0 m | 84 pts | Verified $(X=2200.0, Y=1150.0)$ | $z_{\text{strat}} = 84.0 - d$ |
| `litholog5` | Outcrop Measured Section | 85.0 m | 85 pts | Verified $(X=2620.0, Y=1150.0)$ | $z_{\text{strat}} = 85.0 - d$ |
| `litholog6` | Outcrop Measured Section | 82.0 m | 82 pts | Verified $(X=3100.0, Y=1400.0)$ | $z_{\text{strat}} = 82.0 - d$ |
| `litholog7` | Outcrop Measured Section | 80.0 m | 80 pts | Verified $(X=3850.0, Y=1950.0)$ | $z_{\text{strat}} = 80.0 - d$ |
| `litholog8` | Outcrop Measured Section | 79.0 m | 79 pts | Verified $(X=4400.0, Y=2300.0)$ | $z_{\text{strat}} = 79.0 - d$ |
| `litholog9` | Outcrop Measured Section | 78.0 m | 77 pts (1m gap at 18-19m) | Verified $(X=5150.0, Y=2750.0)$ | $z_{\text{strat}} = 78.0 - d$ |
| `litholog10` | Outcrop Measured Section | 77.0 m | 77 pts | Verified $(X=5700.0, Y=3100.0)$ | $z_{\text{strat}} = 77.0 - d$ |
| `litholog11` | Outcrop Measured Section | 78.0 m | 77 pts (1m gap at 59-60m) | Verified $(X=3450.0, Y=1650.0)$ | $z_{\text{strat}} = 78.0 - d$ |
| `litholog12` | Subsurface Drill Core EM-137C | 111.0 m | 111 pts | Verified $(X=12000.0, Y=8000.0)$ | Preserved as drilled ($0-111$ m) |
| **Total** | **12 Lithologs** | **1,033.0 m** | **1,031 pts** | **11 Spatial Wells (940 pts)** | **100% Invariance Preserved** |

### 2.3 Dataset Provenance Status & Validation Boundaries
- **Provenance Manifest**: [`data/provenance_manifest.json`](file:///d:/Lithology-reconstruction-using-XGB/data/provenance_manifest.json)
- **Digitization History**: All 12 litholog profiles were manually digitized and revised from primary source figures in Sahoo et al. (2016).
- **Litholog 11 QC Benchmark**:
  - Independent manual-vs-automated digitization benchmark on Litholog 11 demonstrated **93.59% structural agreement** (73 out of 78 meters matched identically; 5 discrepant meters).
  - Exact facies and bed boundary alignment was confirmed across 68 out of 78 meters (87.2%).
  - **Scientific Boundary**: The 93.59% accuracy benchmark applies strictly to Litholog 11. It cannot and must not be cited as the accuracy of Lithologs 1-10 or Litholog 12.
- **Litholog 12 Provenance**: Manually digitized subsurface drill core log (EM-137C). Digitzed coverage is restricted to 0.0 - 111.0 m (original core extends to 242.0 m). Field validation is undocumented, and orientation confidence is unresolved.

### 2.4 Canonical 6-State Facies Schema & Net-to-Gross Partitioning
Implemented in [`smalt/descriptive/analyzer.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/descriptive/analyzer.py#L35-L65):

| Code | Facies Name | Alias | Depositional Environment & Architectural Role | Stratigraphic Thickness (m) | Proportion (%) |
| :---: | :--- | :--- | :--- | :---: | :---: |
| **0** | `sand` | `channel_sand` | Channel Sandstone: High-energy trough/planar cross-bedded channel-axis and bar complexes | 444.70 m | 43.05% |
| **1** | `p_sand` | `planar_sand` | Planar Sandstone: Upper flow regime planar-laminated sand sheets and bar tops | 63.00 m | 6.10% |
| **2** | `ripples` | `silt` | Rippled Heterolithics / Siltstone: Crevasse splays, levees, and waning-flow interbeds | 77.70 m | 7.52% |
| **3** | `carbon_mud` | `carb_mud` | Carbonaceous Mudstone: Restricted coastal plain, backswamp, and marsh baffles | 41.60 m | 4.03% |
| **4** | `coal` | `coal` | Coal: Mire / peat swamp chronostratigraphic marker beds | 25.00 m | 2.42% |
| **5** | `mud` | `overbank_mud` | Overbank Mudstone: Low-energy floodplain suspension muds and regional seals | 381.00 m | 36.88% |
| **Total** | **6 States** | - | **Complete 12-Litholog Dataset** | **1,033.00 m** | **100.00%** |

#### Net-to-Gross (N/G) Ratio Hierarchy:
- **Channel Sandstone Alone (`ntg_channel`)**: **43.05%** (444.70 m / 1033.0 m)
- **Planar Sandstone Alone (`ntg_planar`)**: **6.10%** (63.00 m / 1033.0 m)
- **Net Sandstone (`ntg_net_sand` = `sand` + `p_sand`)**: **49.15%** (507.70 m / 1033.0 m)
- **Coarse Reservoir Heterolithics (`ntg_coarse` = `sand` + `p_sand` + `ripples`)**: **56.78%** (585.40 m / 1033.0 m)

### 2.5 Spatial Coordinates & Common-Zero Datum Alignment
- **Spatial Module**: [`smalt/spatial/coordinates.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/coordinates.py) loads real local Cartesian coordinates $(X, Y)$ from `lolo/Location_coordinates_lithologs.xlsx` for L2-L12.
- **Litholog 1 Exclusion**: Litholog 1 lacks verified spatial coordinates and is strictly excluded from all spatial conditioning, LOLO benchmarks, and pair extractions.
- **Datum Module**: [`smalt/spatial/datum.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/datum.py) resolves the source-orientation discrepancy:
  - Sahoo et al. (2016) measured outcrop sections (L1-L11) employ the sedimentological convention where 0 m is at the stratigraphic base and height increases upward. Raw digitized CSVs recorded top-down depth ($d = 0$ at top).
  - Transformed L1-L11 via $z_{\text{strat}} = H_{\max} - d_{\text{original}}$, establishing a common-zero stratigraphic datum ($z = 0$ at section base).
  - Preserved Litholog 12 without inversion ($0-111$ m) under unresolved orientation confidence.
  - Stratigraphic invariance is 100% verified: exactly 1,033.0 m cumulative thickness and 328 bed intervals preserved with zero thickness error ($< 10^{-12}$ m).

---

## 3. Established Scientific Decisions & Key Findings

### 3.1 Decision 1: Walther's Law 1D Vertical Succession is a Verified Core Contribution
- **Sedimentological Integrity**: Vertical facies transition tendencies follow Walther's Law and reflect physical depositional processes.
- **Asymmetric Upward Abandonment**: Across all 12 lithologs, when Channel Sandstone transitions upward across a bed boundary, it transitions **100.0% of the time into finer-grained lithologies** (`p_sand` 27.0%, `ripples` 18.9%, `carbon_mud` / `coal` 32.4%, `mud` 21.6%), capturing genuine fining-upward channel abandonment cycles.
- **Transition Perplexity**:
  - Regular 1D chain ($P_{\text{reg}}$, 1m fixed step): Average LOLO perplexity = **2.6289** (reflecting high diagonal self-transition persistence within beds).
  - Embedded 1D chain ($P_{\text{emb}}$, boundary crossings, $P_{ii} = 0$): Average LOLO perplexity = **4.1793** (reflecting unconstrained 5-way branching among distinct alternative facies).
- **Macro-Scale Net-to-Gross Preservation**:
  - Stationary net sand proportion $\pi_{\text{sand}} + \pi_{\text{p\_sand}} = 53.07\%$.
  - Empirical net sand proportion = $54.57\%$.
  - Relative difference is only **2.74%**, proving that the 1D Markov chain faithfully preserves the macro-scale sand budget.

### 3.2 Decision 2: Pointwise Inter-Well Spatial Interpolation is Definitively Deferred on Real Data
- **The Physical Sparsity Root Cause**:
  - Nearest-neighbor inter-well distances range from **420.0 m to 2,355.6 m** (median 701.5 m, mean $>1,200$ m).
  - In contrast, Sahoo et al. (2016) report single-storey channel widths of **140 to 210 m** ($W/T \approx 35$) and crevasse splay widths of **10 to 130 m**.
  - Inter-well spacing is 2x to 10x wider than the maximum lateral continuity of individual sandbodies. Pointwise interpolation is physically and geologically impossible.
- **Empirical Baseline Failure (LOLO Benchmark on 940 Points)**:
  - Spatial 3D KNN ($k=5$): Pooled Raw Acc = **45.96%** (432/940), Balanced Acc = **22.16%**, Macro-F1 = **0.2163**.
  - Nearest-Well Vertical Profile: Pooled Raw Acc = **43.72%** (411/940), Balanced Acc = **22.23%**, Macro-F1 = **0.2209**.
  - Training Prior Majority (`sand`): Pooled Raw Acc = **43.51%** (409/940), Balanced Acc = **16.67%**, Macro-F1 = **0.1011**.
  - **Finding**: Spatial KNN fails to outperform the non-spatial Nearest-Well profile in Macro-F1 (**0.2163 vs 0.2209**) and achieves only a **+2.45 percentage point** gain over the zero-spatial training prior majority baseline.
- **Minority Facies Collapse**:
  - Spatial ML models completely fail on thin-bed facies:
    - `coal` (support = 25): 0/25 identified (F1 = 0.0000).
    - `carbon_mud` (support = 17): 0/17 identified (F1 = 0.0000).
    - `p_sand` (support = 62): 2/62 identified by KNN (F1 = 0.0465).
- **Transport-Parallel Generalization Failure**:
  - Upstream -> Downstream (Train L2-L8, L10; Test L9, L11, L12; 267 pts): Spatial KNN accuracy = **34.46%**, underperforming the training prior (**44.94%** by -10.49 pp).
  - Downstream -> Upstream (Train L9, L11, L12; Test L2-L8, L10; 673 pts): Spatial KNN accuracy = **37.74%**, underperforming the training prior (**42.94%** by -5.20 pp).
- **Presentation Decision**: Framing this rigorous negative result demonstrates scientific maturity and protects the thesis from examination criticism.

### 3.3 Decision 3: Pointwise MAP Collapse vs Stochastic Simulation Superiority
- **Mathematical Mechanism of MAP Collapse**:
  - The continuous horizontal transition probability matrix $\mathbf{P}_h(h) = \exp(\mathbf{R}_h h)$ decays exponentially to regional stationary proportions $\boldsymbol{\pi}$ beyond $h \ge 420$ m.
  - At unconditioned target locations, deterministic argmax classification invariably selects the regional background floodplain facies (Overbank Mudstone), collapsing pointwise accuracy to **37.23%** (350/940).
- **Superiority of Stochastic Simulation**:
  - While deterministic MAP yields low pointwise accuracy, **stochastic sampling** from the conditioned probability vector faithfully reproduces subsurface sedimentary architecture.
  - **Bed Thickness Reproduction**: Ground-truth target well mean bed thickness is **3.76 meters** (Sprint I) and **4.50 meters** (Sprint J). The conditioned stochastic realization ensemble achieves **3.85 meters** and **4.39 meters** (error $< 0.15$ m), whereas spatial KNN produces artificial, vertically smeared beds averaging **8.10 meters** (+115% distortion).
  - **Hard Data Honor Rate**: The stochastic realization generator strictly guarantees a **100.0% honor rate** at known well control points.

### 3.4 Decision 4: Bed Thickness Preservation is 100% Driven by 1D Vertical Succession
- **Sprint J Ablation Findings (4-Model Test)**:
  - Model 0 (Stationary Proportions Only): Shatters stratigraphy into 1m noise; Mean Bed Thickness = 1.54 m (error = 2.96 m vs 4.50 m ground truth).
  - Model 1 (Vertical Markov Succession Only, unconditioned horizontally): Mean Bed Thickness = **4.39 m** (error = **0.15 m** vs 4.50 m ground truth; proportion Total Variation = 0.0244).
  - Models 2 & 3 (Adding Horizontal Conditioning from Distant Wells): Distant well conditioning acts as a static damper because $h \ge 420$ m decays to stationary proportions, inflating bed thickness to 15.75 m.
  - **Takeaway**: The realistic bed thickness reproduction reported in Sprint I was **100% inherited from the 1D vertical Markov succession**, not from horizontal conditioning.

### 3.5 Decision 5: Empirical Facies Continuity Calibration & Decay Dynamics
- **Horizontal Facies Pair Extraction**: Exactly **4,410 horizontal facies pairs** extracted across **93 common-zero elevation slices** ($z \in [0.0, 110.0]$ m) among the 11 spatial wells (L2-L12).
- **Calibrated Lateral Continuity Lengths**:
  - Channel Sandstone: $L_{\text{geological}} = 203.0$ m, $L_{\text{empirical}} = 243.7$ m, $L_{\text{fitted}} = \mathbf{402.7\text{ m}}$ ($\pm 125.5$ m, 840 auto-pairs, high support). Captures multi-storey channel-belt amalgamation.
  - Overbank Mudstone: $L_{\text{geological}} = 360.0$ m, $L_{\text{empirical}} = 307.1$ m, $L_{\text{fitted}} = \mathbf{343.4\text{ m}}$ ($\pm 164.4$ m, 576 auto-pairs, high support).
  - Planar Sandstone: $L_{\text{geological}} = 30.0$ m, $L_{\text{empirical}} = 226.2$ m, $L_{\text{fitted}} = \mathbf{80.3\text{ m}}$ ($\pm 457.1$ m, 19 auto-pairs).
  - Rippled Heterolithics: $L_{\text{geological}} = 75.0$ m, $L_{\text{empirical}} = 1.4$ m, $L_{\text{fitted}} = \mathbf{7.8\text{ m}}$ ($\pm 0.0$ m, 30 auto-pairs).
  - Coal & Carbonaceous Mudstone: Auto-pair count = 1 and 3 across 4,410 pairs; formally audited as **"insufficient spatial support"** (unidentifiable at well spacing $\ge 420$ m; fall back to geological priors 55.0 m and 40.0 m).
- **Alternative Decay Models Evaluated**:
  - Model A (Sprint I Prescribed): MSE = 0.018991, RMSE = 0.1378.
  - Model B (Sprint J Calibrated Markov): MSE = 0.010790, RMSE = 0.1039 (43.2% MSE reduction).
  - Model C (Spherical Finite-Range Decay): MSE = 0.006991, RMSE = 0.0836 (63.2% MSE reduction).
  - All models strictly satisfy row-stochasticity, non-negativity, $P(0) = \mathbf{I}$, and asymptotic stationary convergence.
- **Synthetic Correlation-Length Recovery Benchmark**:
  - Tested on synthetic multi-well transect with known correlation lengths:
    - Medium-range (true 250.0 m): estimated **265.1 m** (6.0% error) -> **PASSED**.
    - Long-range (true 600.0 m): estimated **610.1 m** (1.7% error) -> **PASSED**.
  - Proves estimator is mathematically and numerically correct; real-data limitations stem strictly from physical borehole sparsity.

### 3.6 Decision 6: Formal Decision Tree Classification
- Formally classified under **CATEGORY C**:
  > *"The characteristic lateral dimensions of several facies bodies are substantially smaller than the observed inter-well spacing, limiting deterministic lateral predictability."*

### 3.7 Decision 7: Strategic Re-Scoping of Downstream Phases for UGP Presentation
- **Phase 1 (1D Vertical Markov Engine)**: **PROCEED AS CORE CONTRIBUTION**. Verified across 1,031 m; proves Walther's Law upward fining; establishes vertical succession rules.
- **Phase 2 (2D Fluvial Forward Modeling)**: **PROCEED AS UNCONDITIONED PROCESS SIMULATOR**. Implement stochastic ribbon channel generator with $W/T \approx 35$, mean thickness $\approx 5.8$ m, and target NTG envelope (17-50%) from Sahoo et al. (2016) geometrical priors.
- **Phase 3 (Inter-Well Spatial Reconstruction)**: **DEFER ON FIELD DATA**. Definitively documented as an ill-posed problem given current well spacing.
- **Phases 4-5 (Active Learning & Borehole Querying)**: **RE-SCOPE TO SYNTHETIC BENCHMARK**. Demonstrate Active Margin Sampling on synthetic 2D forward realizations (from Phase 2) where true subsurface ground truth is known, proving coring budget optimization.

---

## 4. Completed Work Breakdown (Milestones & Deliverables)

### Phase 0: Data Standardization & Quality Control
- Unified ingestion pipeline loading all incoming lithologs into a 6-column schema stored at [`data/processed/lithologs_unified.parquet`](file:///d:/Lithology-reconstruction-using-XGB/data/processed/lithologs_unified.parquet).
- Implemented automated QC gates: strict depth monotonicity, zero-thickness interval removal, and synthetic Gamma Ray proxy generation.
- Deliverables: [`docs/notes/phase0_schema.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/notes/phase0_schema.md), [`tests/test_phase0_loader.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_phase0_loader.py) (11 tests).

### Phase 1: 1D Stratigraphic Markov Engine
- Core module [`smalt/geostat/markov.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/geostat/markov.py) computing regular ($P_{\text{reg}}$) and embedded ($P_{\text{emb}}$) transition matrices.
- Invariant stationary distribution solver ($\boldsymbol{\pi} \mathbf{P} = \boldsymbol{\pi}$), transition asymmetry matrix ($\mathbf{A} = \mathbf{P} - \mathbf{P}^T$), and state transition entropy.
- Provenance sensitivity audit comparing complete working dataset vs source-derived 3-log subset (L1, L9, L11).
- Deliverables: [`docs/notes/phase1_markov.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/notes/phase1_markov.md), [`tests/test_phase1_markov.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_phase1_markov.py) (11 tests).

### Sprint C: Descriptive Lithology & Initial 1D Markov Baseline
- Audited all 12 litholog profiles, descriptive facies thickness tables, and preliminary cross-validation.
- Deliverables: [`sprints/audit_sprint_c/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_c/), [`tests/test_sprint_c.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_c.py) (12 tests).

### Sprint D: Reproducible 1D Markov Validation & Summary Reconciliation
- Reconciled sandstone bed statistics directly from raw CSVs: exactly 72 raw intervals (444.70 m) and 53 merged distinct sandstone bodies.
- Formalized mathematical scoring separating initial-state log-likelihood from transition perplexity.
- Executed leave-one-litholog-out cross-validation across all 12 logs.
- Deliverables: [`sprints/audit_sprint_d/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_d/), [`tests/test_sprint_d.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_d.py) (10 tests).

### Sprint E: Revised Dataset Reconciliation & Spatial Readiness
- Reconciled 12 revised lithologs manually digitized from primary source images.
- Registered Litholog 12 in provenance manifest as a subsurface drill core log (0-111 m).
- Audited synthetic strike-coordinate formula; established real Cartesian coordinate module [`smalt/spatial/coordinates.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/coordinates.py).
- Deliverables: [`sprints/audit_sprint_e/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_e/), [`docs/notes/sprint_e_dataset_reconciliation_and_spatial_readiness.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/notes/sprint_e_dataset_reconciliation_and_spatial_readiness.md), [`tests/test_sprint_e.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_e.py) (11 tests).

### Sprint F: Six-State Migration & Spatial Baseline Revalidation
- Migrated codebase to canonical 6-state facies schema (`sand`, `p_sand`, `ripples`, `carbon_mud`, `coal`, `mud`).
- Recomputed 6x6 regular and embedded Markov chains across all 12 logs.
- Executed spatial LOLO cross-validation on 940 grid points across eligible wells (L2-L12).
- Deliverables: [`sprints/audit_sprint_f/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/), [`docs/notes/sprint_f_six_state_spatial_revalidation.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/notes/sprint_f_six_state_spatial_revalidation.md), [`tests/test_sprint_f.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_f.py) (10 tests).

### Sprint G: Independent Results Audit & Spatial Modeling Decision
- Independently recomputed and audited all Sprint F metrics with 0 discrepancy.
- Discovered and diagnosed the `StandardScaler` vertical anisotropy defect (pre-scaling weight factor cancelled by standard deviation division).
- Established 20-day presentation runway roadmap: proceed with Phase 1 + scoped Phase 2; defer Phase 3 on real wells; re-scope Phases 4-5 to synthetic benchmarks.
- Deliverables: [`sprints/audit_sprint_g/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_g/), [`sprints/audit_sprint_g/SPRINT_G_AUDIT_REPORT.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_g/SPRINT_G_AUDIT_REPORT.md), [`tests/test_sprint_g.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_g.py) (8 tests).

### Sprint H: Common-Zero Datum Alignment & Orientation Correction
- Audited source vertical orientation against Sahoo et al. (2016) measured-section convention ($0\text{ m} = \text{base}$, upward).
- Inverted top-down measured sections L1-L11 via $z_{\text{strat}} = H_{\max} - d_{\text{original}}$ in [`smalt/spatial/datum.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/datum.py); preserved L12 drill core without inversion.
- Verified 100% stratigraphic invariance (1,033.0 m cumulative thickness, 328 beds, zero thickness error).
- Rebuilt horizontal facies pairs (4,410 pairs across 111 elevation slices) and spatial LOLO benchmark.
- Deliverables: [`sprints/audit_sprint_h/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_h/), [`sprints/audit_sprint_h/SPRINT_H_ORIENTATION_CORRECTION_REPORT.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_h/SPRINT_H_ORIENTATION_CORRECTION_REPORT.md), [`tests/test_sprint_h.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_h.py) (9 tests).

### Sprint I: Conditioned Transition-Probability Geostatistical Prototype
- Implemented coupled Markov transition geostatistical prototype combining vertical succession $\mathbf{P}_v$ and continuous horizontal rate matrix $\mathbf{P}_h(h) = \exp(\mathbf{R}_h h)$.
- Guaranteed 100.0% hard conditioning honor rate at well control points.
- Generated 2D cross-sectional stochastic realization ensembles reproducing ground-truth mean bed thickness within 0.09 m (3.85 m vs 3.76 m).
- Mapped spatial Shannon entropy and demonstrated the mathematical mechanism behind pointwise MAP classification collapse.
- Deliverables: [`sprints/audit_sprint_i/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_i/), [`sprints/audit_sprint_i/SPRINT_I_CONDITIONED_MARKOV_REPORT.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_i/SPRINT_I_CONDITIONED_MARKOV_REPORT.md), [`tests/test_sprint_i.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_i.py) (15 tests).

### Sprint J: Empirical Horizontal Transition Calibration & Lateral Continuity Audit
- Extracted 4,410 empirical horizontal pairs across 93 common-zero elevation slices; computed bootstrap standard errors ($B=200$).
- Calibrated facies lateral continuity lengths (composite channel sandstone $L=402.7$ m, mudstone $L=343.4$ m, planar sand $L=80.3$ m, ripples $L=7.8$ m).
- Evaluated alternative horizontal decay models; Spherical decay reduced transition MSE by 63.2% (from 0.018991 to 0.006991).
- Executed 4-model ablation study proving bed thickness reproduction is driven 100% by vertical succession.
- Validated synthetic correlation-length recovery (6.0% and 1.7% error on medium and long range facies).
- Classified outcome under Category C.
- Deliverables: [`sprints/audit_sprint_j/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_j/), [`sprints/audit_sprint_j/SPRINT_J_HORIZONTAL_CALIBRATION_REPORT.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_j/SPRINT_J_HORIZONTAL_CALIBRATION_REPORT.md), [`tests/test_sprint_j.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_j.py) (15 tests).

---

## 5. Unresolved Issues & Active Limitations

1. **Litholog 1 Coordinates Missing**:
   - `litholog1` is omitted from `lolo/Location_coordinates_lithologs.xlsx`.
   - It cannot participate in horizontal pair extraction, spatial LOLO cross-validation, or 2D conditional geostatistical modeling.
   - It remains restricted to 1D vertical Markov succession analysis.
2. **Litholog 12 Vertical Datum & Orientation Confidence Unresolved**:
   - `litholog12` represents a subsurface drill core (EM-137C) with only the 0.0 - 111.0 m interval digitized (out of 242.0 m total core length).
   - Its absolute collar elevation datum relative to the outcrop base is undocumented.
   - Its orientation confidence is classified as `unresolved`. It is preserved without inversion to prevent arbitrary assumptions.
3. **Severe Inter-Well Spatial Sparsity**:
   - Nearest-neighbor well separation spans 420.0 m to 2,355.6 m (median 701.5 m).
   - Because single-storey channel sandbodies are only 140 to 210 m wide, deterministic inter-well facies correlation is physically impossible across this dataset.
4. **Insufficient Spatial Pair Support for Thin-Bed Facies**:
   - Across all 4,410 horizontal pairs, Coal produces only 1 auto-pair and Carbonaceous Mudstone produces only 3 auto-pairs.
   - Empirical lateral continuity lengths for these facies cannot be reliably estimated from field data and must rely on literature priors (55.0 m and 40.0 m).
5. **Absence of Regional Structural Dip Correction**:
   - Outcrop elevation slices currently assume horizontal bedding across the 6.0 km transect.
   - Any subtle regional structural dip (typically 1-3 degrees in the Book Cliffs / Wasatch Plateau) could introduce minor stratigraphic misalignment across distant slices.
6. **Outstanding Data Inquiries for Prof. Hiranya Sahoo**:
   - Clarification of Litholog 1 surface coordinates.
   - Ground elevation / collar datum for drill core EM-137C (Litholog 12).
   - Structural strike and dip measurements across the measured transect.

---

## 6. Verified Test Results & Quality Assurance

### 6.1 Test Suite Status Overview
- **Repository Test Suite Execution**: `python -m pytest tests/`
- **Total Test Count**: **126 tests collected, 126 passing (100% pass rate)**
- **Test Execution Time**: ~489s (verified end-to-end execution across all 126 items with zero failures)
- **Targeted Sprint K Execution**: `python -m pytest tests/test_phase0_loader.py tests/test_phase1_markov.py tests/test_sprint_k.py` (36 passed in 32.5s)

### 6.2 Test Module Breakdown
| Test File | Sprint / Phase | Test Count | Pass Status | Core Verification Scope |
| :--- | :--- | :---: | :---: | :--- |
| [`tests/test_phase0_loader.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_phase0_loader.py) | Phase 0 | 11 | **PASS** | Monotonic depth QC, schema typing, non-negative intervals, synthetic GR proxy. |
| [`tests/test_phase1_markov.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_phase1_markov.py) | Phase 1 | 11 | **PASS** | Regular/embedded Markov transition matrices, stationary distribution, asymmetry matrix. |
| [`tests/test_sprint_c.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_c.py) | Sprint C | 12 | **PASS** | 12-litholog audit, descriptive facies stats, leakage-free 1D LOLO scoring. |
| [`tests/test_sprint_d.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_d.py) | Sprint D | 10 | **PASS** | Sandstone bed count reconciliation (72 intervals, 53 bodies), embedded chain integrity. |
| [`tests/test_sprint_e.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_e.py) | Sprint E | 11 | **PASS** | L12 provenance registration, synthetic coordinate audit, real Cartesian coordinate loader. |
| [`tests/test_sprint_f.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_f.py) | Sprint F | 10 | **PASS** | Canonical 6-state schema, Spatial LOLO cross-validation, confusion matrices, naive baselines. |
| [`tests/test_sprint_g.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_g.py) | Sprint G | 8 | **PASS** | Independent metrics reconciliation, StandardScaler anisotropy bug diagnosis, alias audit. |
| [`tests/test_sprint_h.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_h.py) | Sprint H | 9 | **PASS** | Measured-section inversion ($z_{\text{strat}} = H_{\max} - d$), 100% stratigraphic invariance, 4,410 pairs. |
| [`tests/test_sprint_i.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_i.py) | Sprint I | 15 | **PASS** | Continuous transition rate $R_h$, 100% hard conditioning honor rate, stochastic realization seed reproducibility. |
| [`tests/test_sprint_j.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_j.py) | Sprint J | 15 | **PASS** | Empirical pairs, binning strategies, alternative decay models, 4-model ablation, synthetic length recovery. |
| [`tests/test_sprint_k.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_k.py) | Sprint K | 14 | **PASS** | 2D forward generator, multi-seed validation, leakage-free sparse reconstruction benchmark across spacings. |
| **Total** | **All Modules** | **126** | **126 / 126** | **Zero Failures, 100% Pass Rate Across All Sprints** |

---

## 7. Sprint K Deliverables and Traceability Registers

### 7.1 Authoritative Registers Created
1. **Scientific Methods Register:** [`docs/SCIENTIFIC_METHODS_REGISTER.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/SCIENTIFIC_METHODS_REGISTER.md)
   - Catalog of 25 permanent method/parameter entries (`METHOD-DATA-001`, `PARAM-WT-001`, etc.) with mathematical derivations, literature evidence, and oral defense arguments.
2. **Figure Provenance Register:** [`docs/FIGURE_PROVENANCE_REGISTER.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/FIGURE_PROVENANCE_REGISTER.md)
   - Complete inventory and provenance audit of all scientific figures across Sprint K and historical sprints.
3. **Professor Meeting Notes:** [`docs/PROFESSOR_MEETING_NOTES.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/PROFESSOR_MEETING_NOTES.md)
   - Question-driven oral defense guide addressing 7 core scientific questions with plain-language explanations and likely faculty Q&A.
4. **BibTeX Literature Database:** [`docs/REFERENCES.bib`](file:///d:/Lithology-reconstruction-using-XGB/docs/REFERENCES.bib)
   - Verified BibTeX entries for Sahoo et al. (2016), Abasaghi et al. (2026), Weissmann & Fogg (1999), Dindarloo et al. (2015), He et al. (2019), and foundational geostatistical literature.
5. **Scientific Audit Summary:** [`docs/SCIENTIFIC_AUDIT_SUMMARY.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/SCIENTIFIC_AUDIT_SUMMARY.md)
   - Executive review document summarizing resolved provenance issues, verified metrics, remaining uncertainties, and supervisor recommendations.

---

## 8. Actionable Next Steps & Immediate Roadmap

### Priority 1: Execute Sprint L - Borehole Active Learning & Optimal Coring Benchmark
- Implement Active Coring / Margin Sampling on verified 2D forward realizations (from Sprint K) where synthetic ground truth is 100% known.
- Quantify borehole placement optimization (e.g. reduction in uncertainty entropy with optimal borehole placement vs uniform grid drilling).
- Evaluate mutual information and predictive entropy reduction to detect isolated channel sandbodies with minimal coring budget.

### Priority 2: Prepare Faculty Meeting Presentation
- Assemble slide deck using the question-driven oral defense structure from [`docs/PROFESSOR_MEETING_NOTES.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/PROFESSOR_MEETING_NOTES.md):
  1. The 1D Markov vertical succession findings (Walther's Law, 100% upward fining/abandonment).
  2. The common-zero datum alignment and source-orientation correction.
  3. The rigorous negative result on inter-well spatial interpolation (explaining the physical scale mismatch between well spacing and channel width).
  4. The continuous transition-probability framework (proving why deterministic MAP collapses while stochastic simulation preserves bed architecture).
  5. The 2D fluvial forward model and synthetic active learning pipeline.

