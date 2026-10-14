# SMALT Sprint K - Stage 4 Final Scientific Reasoning and Reproducibility Audit Report

**Project:** Subsurface Stratigraphic Modeling & Active Learning Toolkit (SMALT)  
**Document ID:** `SMALT-AUDIT-SPRINT-K-STAGE4-001`  
**Date:** October 2026  
**Auditor:** Principal Geostatistical Scientist, Sedimentologist, and Research-Software Auditor  
**Repository Branch:** `sprint-k`  
**Baseline Commit:** `f08587ca107b518c8baac560339381c53af9ed1f`  
**Cross-References:**
- Scientific Claims Register: [`sprints/audit_sprint_k/scientific_claims_register.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/scientific_claims_register.csv)
- Parameter Provenance Register: [`sprints/audit_sprint_k/parameter_provenance.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/parameter_provenance.csv)
- Metric Definitions and Semantics: [`sprints/audit_sprint_k/metric_definitions.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/metric_definitions.md)
- Reproducibility Manifest: [`sprints/audit_sprint_k/reproducibility_manifest.json`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/reproducibility_manifest.json)
- Primary Literature BibTeX: [`docs/REFERENCES.bib`](file:///d:/Lithology-reconstruction-using-XGB/docs/REFERENCES.bib)
- Scientific Methods Register: [`docs/SCIENTIFIC_METHODS_REGISTER.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/SCIENTIFIC_METHODS_REGISTER.md)

---

## 1. Executive Summary

Stage 4 conducted an exhaustive scientific reasoning, provenance, and reproducibility audit of SMALT Sprint K Stages 1-3. The audit evaluated every substantive scientific claim, parameter value, metric definition, empirical calculation, simulation campaign, benchmark protocol, and software safeguard against primary field observations (Sahoo et al., 2016), geostatistical literature (Carle & Fogg, 1996; Weissmann & Fogg, 1999; Dindarloo et al., 2015; Abasaghi et al., 2026; He et al., 2019), and active repository code.

### Key Audit Findings:
1. **Mathematical Reproducibility Confirmed:** All empirical tables, forward simulation metrics, sensitivity curves, and benchmark scores reproduce exactly from repository source files and code without data invention or cherry-picking. Automated regression testing confirms 100% test passage (19/19 Sprint K tests; 41/41 targeted tests in 38.5s).
2. **Net-to-Gross Target Discrepancy Mechanistically Resolved:** The apparent discrepancy where a 45.0% target yielded $38.09\%$ realized channel sand is proven to be an algorithmic stopping rule artifact. The generator checks single-facies `sand` (Facies 0) against the tolerance threshold ($0.45 - 0.08 = 0.37$). Because channels simultaneously place basal planar sand (`p_sand`, Facies 1; configured at 20%), total pure channel-belt sand ($N/G_{\text{pure}} = \text{sand} + \text{p\_sand}$) realizes at **$47.55\% \pm 3.13\%$**, matching the empirical Upstream target ($50.52\%$) within $2.97\%$ error.
3. **Log-Normal Goodness-of-Fit Rigorously Corrected:** The ordinary Kolmogorov-Smirnov test result ($D = 0.1059, p = 0.8023$) was confirmed to be conservatively biased because parameters $\mu = 1.7870, \sigma = 0.8823$ were estimated on the test data. A newly implemented parameter-estimation-aware Lilliefors bootstrap test yields $p = \mathbf{0.4145}$ (and $p = \mathbf{0.6940}$ when accounting for 1.0 m integer discretization). While non-rejection supports log-normality as a defensible engineering parameterization, cluster-bootstrapping across the 8 originating lithologs establishes wide parameter confidence intervals ($\mu \in [1.5118, 2.1619]$, corresponding to median bed thickness between $4.53\text{ m}$ and $8.69\text{ m}$).
4. **Physical Scale Mismatch (Category C) Formulated as a Major Negative Result:** The benchmark proves that when inter-well spacing ($400 - 1000\text{ m}$) exceeds channel belt width ($W \approx 287\text{ m}$), deterministic interpolation collapses toward background mudstone proportions ($75.9\% \to 66.5\% \to 57.7\%$). Nearest Well outperforms 2D Spatial KNN because it preserves vertical Waltherian succession and bed thicknesses, whereas 2D KNN horizontally blurs thin beds.
5. **Strict Data Integrity and Safeguards Enforced:** Confirmed that all 12 lithologs were manually digitized from published figures (no AI generation); Litholog 1 is strictly excluded from spatial distance calculations due to missing coordinates; Core EM-137C is strictly bounded to 0-111 m across 63 beds; synthetic coordinates are quarantined; and synthetic gamma ray is prohibited from being used as a predictive machine-learning feature.

---

## 2. Classification Standards of Scientific Evidence

To maintain scientific defensibility, every statement and parameter in SMALT is classified under one of the following evidence standards:

- **Verified Implementation Fact:** An invariant guaranteed by software logic, deterministic formulas, or automated regression tests.
- **Empirical Observation / Calibrated Parameter:** Directly measured or numerically estimated from the 12 digitized field lithologs (1033.0 m, 296 beds).
- **Literature-Derived Prior:** Extracted from peer-reviewed publications (e.g. channel aspect ratio $W/T \approx 35 \pm 5$ from Sahoo et al., 2016).
- **Modeling Assumption:** Introduced where physical data do not uniquely constrain subsurface geometry (e.g. parabolic channel scour cross-section; overbank mud background).
- **Unresolved Scientific Question:** A geological phenomenon that cannot be constrained from the available 1D data (e.g. true 3D channel sinuosity and paleoflow variance).
- **Contradicted or Downgraded Claim:** A claim that was previously overstated or mislabeled (e.g. claiming ordinary KS $p=0.80$ proves log-normality; claiming L2-L10 were AI-reconstructed).

---

## 3. Scientific Claims and Evidence Audit

A formal claim register comprising 20 substantive claims was generated in [`sprints/audit_sprint_k/scientific_claims_register.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/scientific_claims_register.csv). Below is the audit summary across key domains:

### 3.1 Field Data Provenance & Digitization
- **`CLAIM-DATA-001` (Manual Digitization):** Verified. All 12 lithologs represent manual human digitization from published cliff plates and core photographs in Sahoo et al. (2016). Outcrop log 11 was benchmarked against automated optical extraction, achieving $93.59\%$ agreement ($73/78\text{ m}$). Mislabeled historical references to "AI-reconstructed" logs were corrected across all metadata and code.
- **`CLAIM-DATA-002` (Canonical Six-State Schema):** Verified. Preserves all six sedimentological facies from Sahoo et al. (2016, Table 1). Facies codes $\{0, 1, 2, 3, 4, 5\}$ are immutable. Planar sand and carbonaceous mud are never aggregated into generic buckets.
- **`CLAIM-DATA-003` (Litholog 1 Spatial Quarantine):** Verified. Field coordinates for Litholog 1 were omitted from `lolo/Location_coordinates_lithologs.xlsx`. L1 is included in 1D vertical Markov tallies but strictly excluded from 2D spatial models.
- **`CLAIM-DATA-004` (Core EM-137C 0-111 m Coverage):** Verified. EM-137C drilled 242.0 m, but only the upper 111.0 m across 63 beds was digitized. The undigitized lower 131.0 m is strictly unobserved; no observations are fabricated.
- **`CLAIM-DATUM-001` (Common-Zero Coal Datum):** Verified. All sections are vertically flattened to the base of the regional Axel Anderson coal zone ($Z = 0.0\text{ m}$), aligning logs along a chronostratigraphic paleohorizontal surface.

### 3.2 Channel Thickness and Distribution Modeling
- **`CLAIM-THICK-001` (Log-Normal MLE Parameterization):** Verified as a candidate parameterization. MLE on 34 Upstream channel sandstone beds yielded $\mu = 1.7870, \sigma = 0.8823$.
- **`CLAIM-THICK-002` (Goodness-of-Fit and Lilliefors Bootstrap):** Verified and corrected. The ordinary KS test ($D = 0.1059, p = 0.8023$) was conservatively biased. The parameter-estimation-aware Lilliefors bootstrap test yields $p = \mathbf{0.4145}$ (and $p = \mathbf{0.6940}$ under 1m discretization). Non-rejection supports log-normality as a defensible model, but does not prove a universal physical law. Cluster-bootstrapping across 8 lithologs establishes that the true effective sample size is approximately 8, with $95\%\text{ CI}$ for $\mu$ spanning $[1.5118, 2.1619]$ (median thickness $4.53 - 8.69\text{ m}$).
- **`CLAIM-THICK-003` (Measurement Discretization):** Verified. All 34 channel sandstone beds in the Upstream outcrop logs have exact integer-meter thicknesses ($1, 2, 3, \dots, 32\text{ m}$), reflecting 1m drafting resolution on published outcrop plates.

### 3.3 Lateral Architecture and Net-to-Gross
- **`CLAIM-WT-001` (Aspect Ratio Prior $W/T \approx 35 \pm 5$):** Verified as a literature prior. Sahoo et al. (2016, page 52) calculated $W/T = 35$ from fully preserved channel sandbodies and validated it against LiDAR measurements (Rittersbacher et al., 2013). Sensitivity experiments across $W/T \in \{20, 35, 50\}$ confirm its first-order control on lateral sand amalgamation.
- **`CLAIM-SPLAY-001` (Crevasse Splay Widths $20-100\text{ m}$):** Verified. Honors outcrop measurements where $91\%$ ($39/43$) of splays were under 100 m wide (Sahoo et al., 2016, page 61, Figure 10D).
- **`CLAIM-COAL-001` (Coal Thickness $1.0\text{ m}$):** Verified. Matches empirical median ($1.00\text{ m}$) and mode ($1.00\text{ m}$) in Upstream outcrop logs.
- **`CLAIM-NTG-001` (Net-to-Gross Stopping Rule):** Verified. The configured 45% target terminates when single-facies `sand` reaches $37.0\%$, realizing $38.09\% \pm 2.28\%$ `sand` and $47.55\% \pm 3.13\%$ pure sand (`sand` + `p_sand`).

### 3.4 Benchmark Findings and Spatial Negative Result
- **`CLAIM-BENCH-001` (Nearest Well Superiority):** Verified on tested synthetic fluvial stratigraphy. Nearest Well preserves vertical Waltherian succession and bed thicknesses intact, whereas 2D KNN horizontally blurs thin beds.
- **`CLAIM-BENCH-002` (Monotonic Accuracy Decay with Spacing):** Verified. Nearest Well raw accuracy decays from $75.91\%$ at 400 m to $66.50\%$ at 600 m and $57.74\%$ at 1000 m.
- **`CLAIM-SPATIAL-001` (Category C Spatial Sparsity Finding):** Verified. When borehole spacing exceeds channel belt width ($420 - 1000\text{ m}$ vs $287\text{ m}$), deterministic inter-well interpolation is physically unfeasible without geometric object constraints, collapsing to background mudstone proportions.
- **`CLAIM-HET-001` (Downstream Outcrop Ripple Discrepancy):** Verified. The $16.1\%$ ripple deficit in the downstream preset is caused by Litholog 9's anomalous 23.0 m continuous ripple unit, reflecting small-sample variance ($N=2$ logs).

---

## 4. Parameter Provenance and Calibration Audit

Every parameter in `ForwardSimulatorConfig` and its child dataclasses was audited and recorded in [`sprints/audit_sprint_k/parameter_provenance.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/parameter_provenance.csv):

| Component | Parameter | Configured Value | Units | Provenance Category | Audit Verdict |
| :--- | :--- | :---: | :---: | :--- | :--- |
| **GridConfig** | `x_extent_m` | 2000.0 | m | Modeling Assumption | Modeling Assumption |
| | `z_extent_m` | 50.0 | m | Modeling Assumption | Modeling Assumption |
| | `dx_m` | 10.0 | m | Computational Setting | Computational Setting |
| | `dz_m` | 0.5 | m | Computational Setting | Computational Setting |
| **ChannelBeltConfig** | `thickness_mu` | 1.7870 | - | Estimated from Data | Calibrated Parameter |
| | `thickness_sigma` | 0.8823 | - | Estimated from Data | Calibrated Parameter |
| | `thickness_min_m` | 1.5 | m | Modeling Assumption | Modeling Assumption |
| | `thickness_max_m` | 25.0 | m | Modeling Assumption | Modeling Assumption |
| | `width_to_thickness_ratio` | 35.0 | - | Literature-Derived | Literature Prior |
| | `width_to_thickness_std` | 5.0 | - | Literature-Derived | Literature Prior |
| | `core_facies_code` | 0 (`sand`) | code | Measured from Lithologs | Verified Empirical |
| | `basal_facies_code` | 1 (`p_sand`) | code | Measured from Lithologs | Verified Empirical |
| | `basal_fraction` | 0.20 | - | Estimated from Data | Calibrated Parameter |
| | `include_basal_planar_sand` | True | bool | Modeling Assumption | Modeling Assumption |
| **CrevasseSplayConfig** | `enabled` | True | bool | Modeling Assumption | Modeling Assumption |
| | `width_min_m` | 20.0 | m | Literature-Derived | Literature Prior |
| | `width_max_m` | 100.0 | m | Literature-Derived | Literature Prior |
| | `thickness_min_m` | 0.5 | m | Estimated from Data | Calibrated Parameter |
| | `thickness_max_m` | 2.5 | m | Estimated from Data | Calibrated Parameter |
| | `lobe_facies_code` | 2 (`ripples`) | code | Measured from Lithologs | Verified Empirical |
| | `margin_facies_code` | 3 (`carbon_mud`) | code | Measured from Lithologs | Verified Empirical |
| | `margin_fraction` | 0.25 | - | Modeling Assumption | Modeling Assumption |
| **CoalHorizonConfig** | `enabled` | True | bool | Modeling Assumption | Modeling Assumption |
| | `num_horizons` | 2 | count | Estimated from Data | Calibrated Parameter |
| | `thickness_m` | 1.0 | m | Estimated from Data | Calibrated Parameter |
| | `facies_code` | 4 (`coal`) | code | Measured from Lithologs | Verified Empirical |
| | `relative_elevations` | `[0.25, 0.75]` | - | Modeling Assumption | Modeling Assumption |
| **BackgroundConfig** | `primary_facies_code` | 5 (`mud`) | code | Measured from Lithologs | Verified Empirical |
| **ForwardSimulatorConfig** | `target_net_to_gross` | 0.45 | - | Literature-Derived | Calibrated Parameter |
| | `tolerance_net_to_gross` | 0.08 | - | Computational Setting | Computational Setting |
| | `max_channel_attempts` | 50 | count | Computational Setting | Computational Setting |

---

## 5. Metric Semantics and Resolution of Net-to-Gross Target vs Realized Discrepancy

Detailed metric definitions were codified in [`sprints/audit_sprint_k/metric_definitions.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/metric_definitions.md).

### 5.1 Net-to-Gross Hierarchy
1. **$N/G_{\text{channel}}$ (Channel Core Sandstone):**
   $$N/G_{\text{channel}} = \frac{N_{\text{cells}}(\text{Facies } 0)}{N_{\text{total\_cells}}} = \frac{\sum \text{Thickness}(\text{sand})}{\text{Total Column}}$$
2. **$N/G_{\text{pure}}$ (Pure Channel-Belt Sandstone):**
   $$N/G_{\text{pure}} = \frac{N_{\text{cells}}(\text{Facies } 0) + N_{\text{cells}}(\text{Facies } 1)}{N_{\text{total\_cells}}} = \frac{\sum \text{Thickness}(\text{sand} + \text{p\_sand})}{\text{Total Column}}$$
3. **$N/G_{\text{coarse}}$ (Total Coarse Siliciclastics):**
   $$N/G_{\text{coarse}} = \frac{N_{\text{cells}}(\text{Facies } 0) + N_{\text{cells}}(\text{Facies } 1) + N_{\text{cells}}(\text{Facies } 2)}{N_{\text{total\_cells}}}$$

### 5.2 Resolution of the 45% Target Discrepancy
- The stopping rule checks:
  $$\text{current\_sand} \ge \text{target\_ntg} - \text{tolerance\_ntg} = 0.45 - 0.08 = \mathbf{0.37}$$
- `current_sand` counts strictly Facies 0 (`sand`), excluding `p_sand` (Facies 1).
- Loop termination occurs as soon as `sand` reaches $37.0\%$, realizing a mean of **$38.09\% \pm 2.28\%$**.
- Basal planar sand adds an additional **$9.46\% \pm 1.20\%$** to the channel belt.
- Realized pure channel-belt sand is **$47.55\% \pm 3.13\%$**, matching the empirical field target ($50.52\%$) within $2.97\%$ absolute error.
- **Decision:** Do not silently modify labels or force artificial fits. Expose `ntg_channel` ($38\%$) and `ntg_pure` ($48\%$) transparently.

---

## 6. Empirical Reference Data and Statistical Resampling Audit

### 6.1 Empirical Group Summaries

Data verified against [`sprints/audit_sprint_k/litholog_group_facies_proportions.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/litholog_group_facies_proportions.csv):

| Litholog Group | Wells | Beds | Total Meters | $N/G_{\text{channel}}$ | $N/G_{\text{pure}}$ | `ripples` | `coal` | `mud` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Upstream** (L2-L8, L10) | 8 | 166 | 673.0 m | $42.94\%$ | $50.52\%$ | $4.90\%$ | $1.93\%$ | $41.75\%$ |
| **Downstream Outcrop** (L9, L11) | 2 | 44 | 156.0 m | $40.38\%$ | $48.08\%$ | $21.15\%$ | $1.92\%$ | $24.36\%$ |
| **Downstream Composite** (L9, L11, L12) | 3 | 107 | 267.0 m | $45.21\%$ | $49.70\%$ | $16.03\%$ | $4.04\%$ | $26.25\%$ |
| **L12 Core Alone** (EM-137C) | 1 | 63 | 111.0 m | $51.98\%$ | $51.98\%$ | $8.83\%$ | $7.03\%$ | $28.92\%$ |
| **Pooled All** (L1-L12) | 12 | 296 | 1033.0 m | $43.05\%$ | $49.15\%$ | $7.63\%$ | $2.98\%$ | $38.34\%$ |
| **Pooled Spatial** (L2-L12) | 11 | 273 | 940.0 m | $43.59\%$ | $50.29\%$ | $8.06\%$ | $2.53\%$ | $37.35\%$ |

### 6.2 Channel Bed-Thickness Goodness-of-Fit and Clustering Audit
- **Empirical Sample:** $N = 34$ channel sandstone beds originating from 8 Upstream lithologs (Litholog 2: 5 beds; L3: 3 beds; L4: 3 beds; L5: 5 beds; L6: 5 beds; L7: 6 beds; L8: 3 beds; L10: 4 beds).
- **MLE Parameterization:** $\hat{\mu} = 1.7870, \hat{\sigma} = 0.8823$.
- **Ordinary Kolmogorov-Smirnov Test:** $D = 0.1059, p = 0.8023$ (inflated due to parameter estimation from test sample).
- **Lilliefors Parametric Bootstrap Test:** $p_{\text{Lilliefors}} = \mathbf{0.4145}$ ($B=2000$ iterations refitting $\mu, \sigma$ for each synthetic sample).
- **Discretization-Aware Lilliefors Bootstrap Test:** $p_{\text{discrete}} = \mathbf{0.6940}$ (accounting for 1.0 m integer measurement rounding).
- **Cluster Bootstrap across 8 Lithologs:**
  - $\mu$: Mean $= 1.7986$, $95\%\text{ CI} = [1.5118, 2.1619]$ (median thickness $4.53 - 8.69\text{ m}$).
  - $\sigma$: Mean $= 0.8645$, $95\%\text{ CI} = [0.7293, 0.9922]$.
- **Audit Conclusion:** The log-normal distribution cannot be rejected ($p > 0.05$), supporting its use as a defensible engineering model. However, because beds are clustered within only 8 lithologs, parameter uncertainty is substantial and must be reported transparently.

---

## 7. Geological Validation Campaign Audit

Audit of [`smalt/forward/validation.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/forward/validation.py) and output CSV tables confirmed:
- **Seed Reproducibility:** Realizations are 100% bitwise deterministic given the same seed.
- **Campaign Execution:** Exactly 30 successful realizations were generated per preset (90 runs total), without failure or silent omissions.
- **Reference Discrepancy Tracking:** Each preset is evaluated against its corresponding empirical group reference. The Upstream preset matches empirical targets closely (discrepancy $-4.85\%$ on `sand`, $+1.88\%$ on `p_sand`, $-0.20\%$ on `coal`).
- **Downstream Preset Deviation:** The large ripple discrepancy ($-16.13\%$) is correctly traced to Litholog 9's 23.0 m ripple unit, not simulator failure.
- **Sensitivity Experiments:** 30 runs across $W/T \in \{20, 35, 50\}$ and 30 runs across $\Delta z \in \{0.25, 0.5, 1.0\}\text{ m}$ were executed and saved in [`sensitivity_wt_and_resolution.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/sensitivity_wt_and_resolution.csv).

---

## 8. Reconstruction Benchmark Validity and Leakage Safeguards Audit

Audit of [`smalt/forward/benchmark.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/forward/benchmark.py) confirmed:
- **Information Quarantine:** Predictive models receive strictly the `SyntheticBoreholeDataset` containing extracted 1D borehole columns. The full 2D grid and geobody bounding boxes remain completely hidden.
- **Observed Column Withholding:** All columns containing sampled boreholes are strictly masked out ($M[:, c] = \text{False}$) during metric computation. Scored cells count equals $N_{\text{total}} - (N_{\text{boreholes}} \times N_z)$.
- **Domain Stratification:** Metrics are computed separately for the interior interpolation domain ($\min(x_b) \le x \le \max(x_b)$) and exterior extrapolation margins.
- **Multi-Metric Evaluation:** Raw accuracy is supplemented by balanced accuracy (macro recall), Macro-F1, and per-facies precision/recall/F1 with explicit handling of absent classes (`zero_division = 0`).
- **Realization-Level Aggregation:** Metrics are computed per realization before computing across-realization means and standard deviations, preventing pseudoreplication.

---

## 9. Coordinate, Data Provenance, and Anti-Leakage Safeguards Audit

The repository enforces six strict scientific safeguards:
1. **Litholog 1 Spatial Quarantine:** L1 has no physical coordinates in `Location_coordinates_lithologs.xlsx` and is excluded from all 2D/3D spatial interpolation scripts. Verified by `test_15`.
2. **Core EM-137C Digitized Boundary:** Bounded strictly to 0-111 m across 63 beds; undigitized 111-242 m section is never fabricated. Verified by `test_15`.
3. **Synthetic Coordinate Quarantine:** All synthetic borehole coordinates are prefixed with `synthetic_` / `synth_` in code and exported datasets, preventing accidental mixing with field wells. Verified by `test_16`.
4. **Synthetic Gamma Ray Leakage Prohibition:** Synthetic gamma ray was deterministically synthesized from facies labels (e.g. sand = 35 API, mud = 105 API). Predictive models are prohibited from using synthetic GR as an input feature (100% label leakage). Verified by `test_17`.
5. **No Deep Learning / Active Learning Without Real Logs:** In the absence of authentic multi-curve downhole physical wireline logs, neural networks and complex classifiers are prohibited.
6. **Separation of Synthetic and Real Findings:** Synthetic benchmark findings evaluate algorithm behavior under known geological geometries; they are never pooled with historical Sprint E-J real-data results.

---

## 10. Automated Test Suite Execution and Verification

All 19 tests in [`tests/test_sprint_k.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_k.py) pass with zero failures:
```powershell
python -m pytest tests/test_sprint_k.py
```
**Outcome:** **19 passed in 38.71s (100% pass rate)**.

Targeted regression suite:
```powershell
python -m pytest tests/test_phase0_loader.py tests/test_phase1_markov.py tests/test_sprint_k.py
```
**Outcome:** **41 passed in 38.52s (100% pass rate)**:
- `test_phase0_loader.py`: 11 passed (depth monotonicity, provenance schema, zero thickness).
- `test_phase1_markov.py`: 11 passed (transition stochasticity, stationary eigenvector, asymmetry).
- `test_sprint_k.py`: 19 passed (facies schema, grid shape, orientation, reproducibility, geometry, run lengths, sampling, ground truth isolation, NTG hierarchy, configuration validation, regression, multi-seed validation, baseline models, spacing decay, leakage prevention, L1/L12 boundary invariants, coordinate quarantine, GR leakage prohibition, Lilliefors goodness-of-fit).

---

## 11. Prioritized Scientific Issues, Discrepancies, and Actionable Resolutions

| Priority | Issue / Discrepancy | Evidence & Root Cause | Actionable Resolution |
| :---: | :--- | :--- | :--- |
| **High** | **Apparent 7% Net-to-Gross deficit in simulator** | Stopping rule checks single-facies `sand` at $37\%$ threshold; basal `p_sand` (20%) is excluded from check. | Document single-facies `ntg_channel` ($38\%$) vs complete channel-belt `ntg_pure` ($48\%$) in all reports and tables. |
| **High** | **Ordinary KS p-value ($p=0.8023$) overstated** | Parameters $\mu, \sigma$ were estimated on the same sample; 34 beds originate from only 8 lithologs. | Replace ordinary KS p-value with Lilliefors bootstrap p-value ($p=0.4145$) and report cluster bootstrap confidence intervals. |
| **Medium** | **Downstream outcrop ripple discrepancy ($5.02\%$ vs $21.15\%$)** | Small-sample variance of $N=2$ outcrop logs; Litholog 9 contains an anomalous 23m ripple unit. | Prioritize Upstream preset ($N=8$ logs) for primary benchmarking; document L9 anomaly. |
| **Medium** | **Measurement discretization in outcrop logs** | All 34 channel beds are exact integer meters due to 1m drafting tick marks in published figures. | Document discretization effects; report that $\Delta z = 0.5\text{ m}$ grid captures outcrop resolution accurately. |
| **Low** | **Channel aspect ratio $W/T \approx 35$ identifiability** | Literature-derived prior from regional analogs; cannot be uniquely identified from 1D logs. | Maintain $W/T$ as an explicit configurable uncertainty parameter; report sensitivity across $[20, 35, 50]$. |

---

## 12. Final Supervisory Checkpoint and Credibility Verdict

### Stage 4 Audit Verdict: **SCIENTIFICALLY DEFENSIBLE AND AUDITABLE**

The forward simulator, multi-seed validation campaign, and synthetic reconstruction benchmark satisfy all criteria of scientific rigor:
1. Every claim is cross-referenced to its source data, literature citation, or code location.
2. Parameter provenance distinguishes empirical observations from literature priors and computational assumptions.
3. Metric definitions and stopping rules are mathematically explicit.
4. Data integrity, coordinate quarantine, and leakage prevention are verified by automated tests.
5. Negative results (Category C spatial sparsity) are formulated with mathematical precision.

Per instructions, **all changes remain uncommitted and execution is stopped at this checkpoint for supervisory review**.
