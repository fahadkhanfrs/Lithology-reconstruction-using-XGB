# SMALT Sprint G - Independent Results Audit & Spatial Modeling Decision Report

**Repository:** `fahadkhanfrs/Lithology-reconstruction-using-XGB`  
**Branch:** `sprint-g` (Starting commit: `f6998cd`)  
**Role:** Principal Geostatistical Machine Learning Scientist and Technical Supervisor for SMALT  
**Project:** Undergraduate Project (UGP) under Prof. Hiranya Sahoo, Department of Earth Sciences, IIT Kanpur  
**Timeline:** Approximately 20 days to final UGP presentation  
**Guiding Principle:** Scientific rigor, mathematical correctness, and reproducible geology over artificial model complexity.

---

## 1. Executive Summary

This audit independently investigates the six-state spatial reconstruction results from Sprint F, evaluates the root causes of the low classification metrics (Spatial 3D KNN pooled Macro-F1 = 0.2232, raw accuracy = 47.77% vs naive majority prior = 43.51%), and establishes a scientifically defensible decision framework for the remaining 20-day runway before the UGP presentation.

### Key Audit Conclusions:
1. **Metrics & Math Reconciled**: All reported Sprint F pooled metrics (raw accuracy 0.4777, balanced accuracy 0.2281, Macro-F1 0.2232 across 940 points) were independently recomputed from the serialized confusion matrices and confirmed with zero discrepancy. Unweighted mean fold balanced accuracy (24.89%) was distinguished from pooled balanced accuracy (22.81%).
2. **Weak Spatial Performance is Geologically Inevitable**: Spatial 3D KNN's inability to outperform the 1D Nearest-Well baseline (Macro-F1 0.2232 vs 0.2276) is not a hyperparameter tuning failure. It is driven by insurmountable geological sparsity: inter-well nearest-neighbor distances (minimum 420.0 m, median 701.5 m) are 2x to 10x wider than individual single-storey channel widths (140-210 m) and splay bodies (10-130 m) documented by Sahoo et al. (2016).
3. **Severe Class Imbalance Distorts Raw Accuracy**: Channel Sandstone (43.5%) and Overbank Mudstone (37.2%) constitute 80.7% of the dataset. Point-wise models achieve 44-48% raw accuracy by predicting majority facies, while completely failing on thin-bed minority facies (`coal` F1 = 0.0000, 0/25 identified; `carbon_mud` F1 = 0.0000, 0/17 identified; `p_sand` F1 = 0.0465, 2/62 identified).
4. **Anisotropy Pre-Scaling Cancellation**: A technical defect was identified in `smalt/spatial/baseline.py`: multiplying the vertical coordinate by an anisotropy factor of 10.0 prior to `StandardScaler` is numerically cancelled by standard deviation division, leaving standardized features unaffected. Corrective testing across actual post-scaling vertical weights (0.1 to 50.0) proved that Macro-F1 remains fundamentally locked between 0.21 and 0.23, demonstrating that spatial KNN is physically unsuited for this problem.
5. **Clear SMALT Phase Roadmap**:
   - **Phase 1 (1D Vertical Markov Analysis)**: Fully verified, scientifically robust, and ready for presentation. Waltherian vertical succession operates independently of horizontal coordinates.
   - **Phase 2 (2D Fluvial Forward Modeling)**: Proceed as an unconditioned stochastic forward model using Sahoo et al. (2016) geometric priors ($W/T \approx 35$).
   - **Phase 3 (Conditional Spatial Reconstruction)**: Definitively rejected/deferred for the real outcrop dataset due to severe inter-well sparsity and unanchored datums. Presenting this negative result is scientifically honest and protects against examination criticism.
   - **Phases 4-5 (Active Learning)**: Scoped down to synthetic benchmark demonstration only, illustrating core-budget optimization on forward models where ground truth is known.

---

## 2. Audit Methodology and Files Inspected

### Files Inspected:
- Raw digitized lithologs: [`data/raw_lithologs/`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/) (`litholog1.csv` through `litholog12.csv`).
- Legacy raw lithologs: `data/litholog1.csv`, `data/litholog9.csv`, `data/litholog11.csv`, `data/log.csv`.
- Spatial coordinates: [`lolo/Location_coordinates_lithologs.xlsx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/Location_coordinates_lithologs.xlsx).
- Sprint F deliverables: [`sprints/audit_sprint_f/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/) (confusion matrices, per-class metrics, LOLO results, datum sensitivity).
- Preprocessing and model code: [`data/loader.py`](file:///d:/Lithology-reconstruction-using-XGB/data/loader.py), [`smalt/spatial/baseline.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/baseline.py), [`smalt/spatial/coordinates.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/coordinates.py).
- Reference literature: [`docs/Sahoo et al 2016.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/Sahoo et al 2016.md) and [`docs/SMALT_Study_and_Build_Plan.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/SMALT_Study_and_Build_Plan.md).

### Verification Protocol:
1. Recomputed all 6x6 confusion matrix sums, diagonal traces, per-class precisions, recalls, F1 scores, balanced accuracies, and macro-F1 values using NumPy and scikit-learn.
2. Verified fold-level cross-validation results and distinguished pooled metrics from unweighted fold averages.
3. Inspected raw CSV strings across all 328 stratigraphic intervals to verify alias propagation and detect ambiguous label handling.
4. Calculated inter-well Euclidean distance matrices and compared them with architectural dimensions from Sahoo et al. (2016).
5. Executed targeted mathematical diagnostics on feature scaling and coordinate anisotropy.

---

## 3. Metrics Reconciliation Table

*Generated by [`scripts/run_sprint_g_audit.py`](file:///d:/Lithology-reconstruction-using-XGB/scripts/run_sprint_g_audit.py); archived in [`sprints/audit_sprint_g/metrics_reconciliation_audit.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_g/metrics_reconciliation_audit.csv)*.

| Metric Scope | Model / Evaluation Target | Reported Sprint F | Independently Recomputed | Discrepancy | Audit Status | Audit Commentary |
| :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| **Total Evaluation Points** | All Spatial Models (L2-L12) | 940 | 940 | 0 | **PASS** | Exact match. L1 (91 pts) strictly excluded due to missing coordinates. |
| **Spatial 3D KNN Raw Accuracy** | Pooled (k=5, Distance-Weighted) | 0.4777 | 0.4777 (449/940) | 0.0000 | **PASS** | Exact match. Pooled over all 940 blind test intervals. |
| **Spatial 3D KNN Balanced Accuracy**| Pooled (Macro-Recall) | 0.2281 | 0.2281 | 0.0000 | **PASS** | Exact match. Mean of recalls across all 6 classes. |
| **Spatial 3D KNN Macro-F1** | Pooled (Macro-F1) | 0.2232 | 0.2232 | 0.0000 | **PASS** | Exact match. Mean of F1 scores across all 6 classes. |
| **Spatial 3D KNN Fold Mean Accuracy**| Unweighted Average of 11 Folds | N/A | 0.4717 | Distinguishable | **PASS** | Fold mean (47.17%) differs from pooled (47.77%) due to unequal well lengths (77-111 m). |
| **Spatial 3D KNN Fold Mean Bal Acc**| Unweighted Average of 11 Folds | N/A | 0.2489 | Distinguishable | **PASS** | Fold mean (24.89%) exceeds pooled (22.81%) because missing classes in test wells reduce denominator. |
| **Nearest-Well Raw Accuracy** | Pooled Profile Baseline | 0.4426 | 0.4426 (416/940) | 0.0000 | **PASS** | Exact match. |
| **Nearest-Well Balanced Accuracy** | Pooled Profile Baseline | 0.2291 | 0.2291 | 0.0000 | **PASS** | Slightly outperforms 3D KNN (0.2291 vs 0.2281). |
| **Nearest-Well Macro-F1** | Pooled Profile Baseline | 0.2276 | 0.2276 | 0.0000 | **PASS** | Outperforms 3D KNN (0.2276 vs 0.2232). |
| **Training Prior Raw Accuracy** | Pooled Majority Class (`sand`) | 0.4351 | 0.4351 (409/940) | 0.0000 | **PASS** | Exact match. Zero-spatial baseline achieves 43.51%. |
| **Training Prior Balanced Accuracy** | Pooled Majority Class (`sand`) | 0.1667 | 0.1667 (1/6) | 0.0000 | **PASS** | Exactly 1.0 / 6 = 0.1667. |
| **Training Prior Macro-F1** | Pooled Majority Class (`sand`) | 0.1011 | 0.1011 | 0.0000 | **PASS** | Exactly 0.6064 / 6 = 0.1011. |
| **Minority Class Support: `coal`** | Ground Truth (940 pts) | 25 | 25 | 0 | **PASS** | Exact match (2.66% of dataset). |
| **Minority Class Recall: `coal`** | Spatial KNN & Nearest-Well | 0.0000 | 0.0000 (0/25) | 0.0000 | **PASS** | Complete collapse: 0/25 identified. |
| **Minority Class Support: `carbon_mud`**| Ground Truth (940 pts) | 17 | 17 | 0 | **PASS** | Exact match (1.81% of dataset). |
| **Minority Class Recall: `carbon_mud`** | Spatial KNN & Nearest-Well | 0.0000 | 0.0000 (0/17) | 0.0000 | **PASS** | Complete collapse: 0/17 identified. |
| **Sample-Level Prediction Persistence**| Storage Audit | Serialized | Serialized CM only | Storage Gap | **NOTED** | 6x6 confusion matrices are serialized in JSON, but point-by-point $(x,y,z,y,\hat{y})$ CSV was omitted. |

---

## 4. Facies Mapping & Alias Audit

*Full analysis in [`sprints/audit_sprint_g/facies_mapping_and_aliases_audit.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_g/facies_mapping_and_aliases_audit.csv)*.

### Authoritative Dataset Hygiene:
Inspection of all 12 revised raw CSVs in [`data/raw_lithologs/`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/) confirmed that every one of the 328 recorded intervals uses strictly the canonical tokens:
- `mud`: 109 intervals
- `sand`: 92 intervals
- `p_sand`: 43 intervals
- `ripples`: 31 intervals
- `coal`: 29 intervals
- `carbon_mud`: 24 intervals

### Alias Dictionary Vulnerabilities in `data/loader.py`:
1. **Generic `silt` and `siltstone` $\to$ `ripples` (HIGH AMBIGUITY)**:
   - In `data/loader.py`, lines 66-67 map `"silt": "ripples"` and `"siltstone": "ripples"`.
   - *Sedimentological Concern*: In Sahoo et al. (2016), Facies 3 is defined as *"Thinly interbedded mudstones, siltstones, and rippled sandstones"*, whereas Facies 4 is defined as *"Mudstones and siltstones"*.
   - A generic field label of "siltstone" without ripple structures is ambiguous between Facies 3 (`ripples`) and Facies 6 (`mud`).
   - *Audit Finding*: The current revised files in `data/raw_lithologs/` contain zero instances of `silt` or `siltstone`. However, legacy files (`data/litholog9.csv`, `data/litholog11.csv`) used `silt`.
   - *Recommendation*: The ingestion pipeline should flag ambiguous generic labels like `silt` or `siltstone` for manual sedimentological confirmation rather than silently coercing them into Facies 3.
2. **`splay` and `fine sandstone` $\to$ `carbon_mud` (CRITICAL DEFECT IN LEGACY ALIAS MAP)**:
   - In `data/loader.py`, lines 76-79 map `"fine sandstone / splay": "carbon_mud"`, `"fine_sandstone": "carbon_mud"`, and `"splay": "carbon_mud"`.
   - *Sedimentological Concern*: Crevasse splay sandstones are coarse, high-energy clastic overbank sheets (related to Facies 1/2), NOT waterlogged organic-rich swamp mudstones (Facies 4/`carbon_mud`)!
   - *Audit Finding*: Fortunately, there are zero occurrences of `splay` or `fine sandstone` in `data/raw_lithologs/`. However, this legacy alias mapping is a severe trap for external data ingestion.
   - *Recommendation*: Remove `splay` $\to$ `carbon_mud` immediately from `data/loader.py` to prevent corrupting future well additions.

---

## 5. Validation & Leakage Audit

### 1. Leave-One-Litholog-Out (LOLO) Leak-Free Integrity:
- Spatial LOLO in [`smalt/spatial/baseline.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/baseline.py) was verified to be strictly leak-free.
- In each fold, the test litholog is completely excluded from:
  - Feature normalization: `StandardScaler` is fitted strictly on $X_{\text{train}}$ and applied to $X_{\text{test}}$.
  - Class frequency priors: Majority prior is computed strictly from $y_{\text{train}}$.
  - Model fitting: KNN is trained strictly on points from the remaining 10 wells.
- Zero sequential leakage: unlike the legacy Phase 4/XGBoost script audited in Sprint B, no true previous-facies labels ($S_{t-1}$) or synthetic Gamma Ray proxies are used as features.

### 2. The StandardScaler Anisotropy Cancellation Defect:
- In `smalt/spatial/baseline.py`, line 148 executes:
  ```python
  X_train_scaled = X_train_raw.copy()
  X_train_scaled[:, 2] *= self.vertical_weight  # default = 10.0
  scaler = StandardScaler()
  X_train_norm = scaler.fit_transform(X_train_scaled)
  ```
- *Mathematical Fact*: `StandardScaler` standardizes by $(x - \mu)/\sigma$. If column 2 is scaled by constant $c = 10.0$, then $\mu' = 10\mu$ and $\sigma' = 10\sigma$, meaning:
  $$\frac{10 Z - 10 \mu_Z}{10 \sigma_Z} = \frac{Z - \mu_Z}{\sigma_Z}$$
- Multiplying before `StandardScaler` is **identically cancelled** (difference $< 10^{-15}$).
- Instead, `StandardScaler` gave equal unit variance to horizontal coordinates ($\sigma_X \approx 4000$ m) and vertical coordinate ($\sigma_Z \approx 28$ m), implicitly creating an anisotropic ratio of $\sigma_X / \sigma_Z \approx 140$.
- *Experimental Follow-Up*: To determine whether fixing this defect improves spatial KNN, we tested explicit post-scaling weights ($w_z \in \{0.1, 1.0, 5.0, 10.0, 20.0, 50.0\}$):
  - Weight 0.1: Macro-F1 = 0.2232, Raw Acc = 47.77%
  - Weight 1.0: Macro-F1 = 0.2236, Raw Acc = 47.87%
  - Weight 10.0: Macro-F1 = 0.2257, Raw Acc = 47.34%
  - Weight 50.0: Macro-F1 = 0.2142, Raw Acc = 44.68%
- *Audit Finding*: Adjusting the vertical weight does not rescue KNN. Macro-F1 is bounded between 0.21 and 0.23 across all weights.

---

## 6. Root-Cause Evidence Table

*Systematic audit of 7 hypotheses; archived in [`sprints/audit_sprint_g/root_cause_evidence_matrix.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_g/root_cause_evidence_matrix.csv)*.

| Hypothesis ID & Name | Description | Classification | Observed Facts | Scientific Interpretation |
| :--- | :--- | :---: | :--- | :--- |
| **H1: Sparse Sampling** | 11 coordinate-bearing wells provide insufficient spatial density for 6 detailed facies. | **SUPPORTED** | Min inter-well distance = 420.0 m (L4-L5); median = 701.5 m; max = 2,355.6 m (L12-L9). Sahoo et al. (2016) report single-storey channel widths of 140-210 m and splay widths of 10-130 m. | Inter-well distance is 2x to 10x wider than individual channel bodies. Point-wise spatial interpolation cannot bridge this geological gap without dense conditioning. |
| **H2: Stratigraphic Misalignment** | Treating measured depth or relative-to-base as equivalent geological elevation causes severe misalignment. | **SUPPORTED** | Prof. Sahoo's stratigraphic datum remains unconfirmed. Outcrops were measured along differing cliff topography without a tied marker. Perturbing datum by $\pm 20$ m yields flat Macro-F1 (0.2164 to 0.2270). | Equal raw depth across wells 420-5000 m apart corresponds to uncorrelated time-stratigraphic intervals. Flat sensitivity proves models are predicting background proportions rather than correlative stratigraphy. |
| **H3: Coordinate Uncertainty** | Local Cartesian coordinates have unverified projection, origin, units, and lack Litholog 1. | **SUPPORTED** | `Location_coordinates_lithologs.xlsx` contains raw X, Y values without CRS header, projection name, or geodetic datum. Row 0 (L1) is completely NaN. | While relative local distances are plausible (~420-5000 m), absolute spatial orientation and regional tie cannot be validated without metadata from Prof. Sahoo. |
| **H4: Spatial Separation & Core Impact** | Subsurface core L12 is located 2.4-3.8 km from outcrops; inter-cluster separation (~14 km) degrades regional transfer. | **SUPPORTED** | On L12, KNN achieves 45.95% accuracy (below prior 50.45%), and Nearest-Well drops to 31.53%. Cross-cluster transfer: Upstream $\to$ Downstream achieves 39.33% (prior 44.94%); Downstream $\to$ Upstream achieves 38.93% (prior 42.94%). | Spatial separation completely destroys predictive skill across the regional divide. Models perform significantly worse than non-spatial majority baselines. |
| **H5: Class Imbalance** | Channel Sandstone and Overbank Mudstone dominate 80.7% of points, causing complete collapse on minority facies. | **SUPPORTED** | Ground truth support: `sand`=409 (43.5%), `mud`=350 (37.2%), `ripples`=77 (8.2%), `p_sand`=62 (6.6%), `coal`=25 (2.7%), `carbon_mud`=17 (1.8%). For `coal`, both models achieved 0/25 correct (F1=0.0000). For `carbon_mud`, both achieved 0/17 correct (F1=0.0000). | Raw accuracy (47.77%) is entirely an illusion created by majority facies prediction. Models have zero capability to detect economically or sedimentologically critical thin beds (coal, carbonaceous mud). |
| **H6: Model Unsuitability** | Euclidean 3D KNN is physically and sedimentologically unsuited for facies architecture reconstruction. | **SUPPORTED** | Testing actual post-scaling vertical weights (0.1 to 50.0) confirmed Macro-F1 remains trapped at 0.21-0.23. KNN fails to outperform 1D Nearest-Well (0.2232 vs 0.2276 Macro-F1). | Point-wise KNN cannot enforce facies continuity, channel geometry ($W/T \approx 35$), or Markovian bed successions. Hyperparameter tuning cannot fix an unphysical modeling paradigm. |
| **H7: Validation Design** | Leave-One-Litholog-Out (LOLO) tests inter-well blind prediction, exposing the true limitations of sparse spatial ML. | **SUPPORTED** | LOLO completely withholds each test well from fitting, scaling, and training. All 940 test evaluations are strictly blind. | LOLO is the scientifically rigorous protocol. High accuracy in legacy studies was an artifact of random train/test splitting within the same wells and leaking sequential features (`prev_facies`). |

---

## 7. Baseline Comparison

*Summary archived in [`sprints/audit_sprint_g/baseline_comparison_summary.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_g/baseline_comparison_summary.csv)*.

| Model / Baseline | Model Type | Conditioning Inputs | Evaluated Population | Pooled Raw Acc | Pooled Bal Acc | Pooled Macro-F1 | Minority F1 (`coal`) | Net Sand Error | Scientific Verdict |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Spatial 3D KNN (k=5)** | Spatial 3D Classifier | Local $(X, Y)$ + Relative Elevation ($Z_{\text{rel}}$) | 940 pts (L2-L12, LOLO) | **0.4777** | **0.2281** | **0.2232** | 0.0000 (0/25) | 8.65% | Fails to outperform 1D Nearest-Well. Adds only +4.26% raw accuracy over zero-spatial prior. Fails completely on minority facies. |
| **Nearest-Well Vertical Profile** | 1D Transferred Baseline | Nearest Well ID + Relative Elevation ($Z_{\text{rel}}$) | 940 pts (L2-L12, LOLO) | **0.4426** | **0.2291** | **0.2276** | 0.0000 (0/25) | 7.91% | Matches or slightly exceeds 3D KNN in Macro-F1 (0.2276 vs 0.2232). Demonstrates that KNN adds zero spatial interpolation value beyond copying the nearest well. |
| **Training Prior Facies** | Zero-Spatial Naive Baseline | Training majority class (`sand`) | 940 pts (L2-L12, LOLO) | **0.4351** | **0.1667** | **0.1011** | 0.0000 (0/25) | 15.65% | Achieves 43.51% raw accuracy without any spatial features, proving that 47.77% KNN accuracy is an artifact of class imbalance. |
| **1D Vertical Markov Succession** | Vertical Stochastic Succession | Vertical step transition $P(S_t \mid S_{t-1})$ | 1031 pts (all 12 logs, LOLO) | *N/A (1D sequence)* | *N/A* | *N/A* | *Preserved (72% to mud)* | Preserves 49.15% N/G | **ROBUST & VERIFIED**. Operates vertically under Walther's Law; does not require horizontal coordinates or unconfirmed datums. Core scientific contribution. |

---

## 8. Coordinate & Datum Limitations

1. **Stratigraphic Datum Unanchored**:
   - Outcrop sections were measured on disparate canyon walls without a correlated marker datum (e.g. Bear Canyon coal seam or Star Point top).
   - Relative-to-base ($Z_{\text{rel}} = \max(\text{depth}) - \text{depth}$) is an arbitrary geometric alignment, not a time-stratigraphic correlation.
   - Sensitivity testing proved that perturbing datum offsets between -20 m and +20 m produces flat performance (Macro-F1 0.2164 to 0.2270), demonstrating that models are fitting background proportions rather than correlative stratigraphy.
2. **Missing Litholog 1 Coordinates**:
   - Litholog 1 has no coordinates in `Location_coordinates_lithologs.xlsx` and must remain strictly excluded from spatial modeling.
3. **Unverified CRS and Units**:
   - The coordinates in `Location_coordinates_lithologs.xlsx` lack EPSG codes, projection parameters, and unit metadata. While distances are consistent with meters, absolute georeferencing is impossible without confirmation from Prof. Sahoo.
4. **Litholog 12 Digitization Coverage**:
   - L12 represents only the upper 111.0 m of an original 242.0 m core. More than 130 m of lower core stratigraphy is unrecorded.

---

## 9. Recommended Phase-by-Phase Decision for SMALT

*Decision matrix archived in [`sprints/audit_sprint_g/phase_decision_matrix.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_g/phase_decision_matrix.csv)*.

With approximately 20 days remaining before the UGP presentation, the following phase-by-phase roadmap is established:

```
[Phase 1: 1D Markov Analysis] ───────────────> PROCEED TO PRESENTATION (Core Empirical Anchor)
[Phase 2: 2D Fluvial Forward Modeling] ──────> PROCEED (Scoped Down to Unconditioned Process Model)
[Phase 3: Inter-Well Spatial Reconstruction] ─> DEFER / SCIENTIFICALLY REJECT with Current Data
[Phases 4-5: ML & Active Learning] ──────────> RE-SCOPE to Synthetic Forward Benchmarks Only
```

### Detailed Justification:

### Phase 1: 1D Vertical Markov Succession Analysis & Descriptive Sedimentology
- **Decision: PROCEED TO FINAL PRESENTATION (CORE CONTRIBUTION)**
- **Justification**:
  - Operates vertically under Walther's Law within each measured section; requires zero horizontal coordinates or cross-canyon datums.
  - Reconciled across 1031 m of validated, manually digitized section in 12 lithologs under the Sahoo et al. (2016) six-state schema.
  - Regular transition perplexity (2.6289) and embedded bed-boundary perplexity (4.1793) provide rigorous, publication-grade quantitative benchmarks.
  - Asymmetric transitions (coal transitioning upward exclusively into mud or carbonaceous mud, never directly into channel sand) provide authentic geological insight.
- **Action for next 20 days**: Freeze the 1D Markov codebase. Prepare high-resolution presentation slides from [`sprints/audit_sprint_f/figures/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/figures/).

### Phase 2: 2D Fluvial Object-Based Forward Modeling
- **Decision: PROCEED AS AN UNCONDITIONED STOCHASTIC PROCESS PROTOTYPE**
- **Justification**:
  - Sahoo et al. (2016) provides explicit architectural priors: channel aspect ratio $W/T \approx 35$, mean thickness $\sim 5.8$ m, splay width ranges (10-130 m), and target Net-to-Gross envelope (17% to 50%).
  - An unconditioned forward model generating synthetic 2D cross-sections (e.g. 1-2 km wide) proves understanding of sedimentary architecture.
  - *Critical Constraint*: Must be presented strictly as a *process-based forward simulator* or *synthetic training benchmark*, NOT as a reconstruction of the specific subsurface between the 12 wells.
- **Action for next 20 days**: Implement a clean 2D ribbon channel generator conditioned on $W/T = 35$ and target N/G.

### Phase 3: Inter-Well Conditional Spatial Reconstruction
- **Decision: DEFER / SCIENTIFICALLY REJECT WITH CURRENT DATA**
- **Justification**:
  - Our audit has proven that point-wise spatial interpolation (KNN or tabular ML) across wells 420 m to 5,000 m apart fails mathematically and geologically (Macro-F1 0.2232; cross-cluster transfer raw accuracy 38.9%; minority classes F1 0.0000).
  - Claiming successful spatial reconstruction on this sparse dataset would be scientifically indefensible and easily dismantled by geological examiners.
  - *The Strong Defense*: Presenting this negative result with rigorous diagnostic proof (inter-well spacing vs channel width, datum invariance, majority class inflation) demonstrates deep scientific maturity.
- **Action for next 20 days**: Do NOT attempt to tune KNN or fit spatial ML to real wells. Document the exact data blockers (marker datum, CRS, inter-well density) as formal requests for future work.

### Phases 4 & 5: Spatial Machine Learning & Active Learning
- **Decision: RE-SCOPE TO SYNTHETIC BENCHMARK DEMONSTRATION**
- **Justification**:
  - Active Margin Sampling cannot be validated on the real 11 wells because spatial ML collapses on real data.
  - However, Active Learning can be demonstrated cleanly on *synthetic 2D cross-sections generated in Phase 2*, where true facies geometries are completely known.
  - This demonstrates the algorithm's capability to identify channel pinch-outs and optimize borehole drilling budgets without making false claims on real field data.
- **Action for next 20 days**: If time permits, run a simple Active Learning query demonstration on a synthetic 2D grid. If time is tight, focus entirely on Phase 1 + Phase 2 + Spatial Audit.

---

## 10. Tests Executed & Reproducibility Instructions

### Tests Executed:
- Dedicated Sprint G test suite: [`tests/test_sprint_g.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_g.py) (8/8 tests passing).
- Full repository test suite: `pytest tests/` (**73 passed, 0 failed in 21.21s**).
  - `tests/test_phase0_loader.py`: 11 passed
  - `tests/test_phase1_markov.py`: 11 passed
  - `tests/test_sprint_c.py`: 12 passed
  - `tests/test_sprint_d.py`: 10 passed
  - `tests/test_sprint_e.py`: 11 passed
  - `tests/test_sprint_f.py`: 10 passed
  - `tests/test_sprint_g.py`: 8 passed

### Reproducibility Instructions:
To reproduce all audit metrics and verify the test suite from scratch:
```bash
# 1. Activate project environment and branch
git checkout sprint-g

# 2. Run Sprint G audit generation script
python scripts/run_sprint_g_audit.py

# 3. Execute Sprint G test suite
python -m pytest tests/test_sprint_g.py -v

# 4. Execute full repository test suite
python -m pytest tests/
```
