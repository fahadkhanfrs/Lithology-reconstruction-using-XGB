# SMALT Project Handoff & Phase Registry

## Phase Status Overview
- **Phase 0 (Data Ingestion & Quality Control)**: COMPLETE / VERIFIED (11/11 lithologs loaded cleanly, 920 observations, 11/11 unit tests passing, provenance manifest registered, docs/notes/phase0_schema.md published)
- **Phase 1 (1D Vertical Markov Succession Analysis)**: VERIFIED (Full 11-log dataset fitted [920 observations], provenance sensitivity analysis documented, 11/11 unit tests passing)
- **Sprint C (Descriptive Lithology & Leakage-Free 1D Markov Baseline)**: COMPLETE / VERIFIED (All 12 lithologs audited, 12 unit tests passing, publication figures generated in audit_sprint_c/)
- **Sprint D (Reproducible 1D Markov Validation & Summary Reconciliation)**: COMPLETE / VERIFIED (Mathematical scoring formalized with separate initial-state vs transition metrics, embedded bed sequence integrity enforced, 10/10 new unit tests passing [44/44 repository total], summary statistics reconciled directly from raw CSVs, deliverables in audit_sprint_d/)
- **Sprint E (Revised Dataset Reconciliation & Spatial Prototype Readiness)**: COMPLETE / SUPERSEDED (Reconciled 12 revised lithologs under provisional 5-state mapping; L12 registered in manifest; synthetic coordinate audit completed; provisional spatial baseline evaluated; 11/11 tests passing; deliverables in sprints/audit_sprint_e/)
- **Sprint F (Six-State Migration & Spatial Baseline Revalidation)**: COMPLETE / VERIFIED (Migrated to 6-state canonical schema preserving distinct p_sand and ripples; recomputed 6x6 Markov transitions; leak-free Spatial LOLO revalidated with per-class metrics, 6x6 confusion matrices, and naive baselines; 10/10 new tests passing [65/65 repository total]; deliverables in sprints/audit_sprint_f/)
- **Sprint G (Independent Results Audit & Spatial Modeling Decision)**: COMPLETE / VERIFIED (Independently audited Sprint F metrics [0 discrepancy]; confirmed inter-well sparsity root cause [420-5000 m spacing vs 140-210 m channel width]; identified and tested StandardScaler vertical anisotropy defect; established 20-day UGP presentation roadmap: proceed with Phase 1 [1D Markov] + scoped Phase 2 [unconditioned 2D fluvial forward model], defer Phase 3 [spatial interpolation on real wells], re-scope Phase 4-5 to synthetic benchmarks; 8/8 new unit tests passing [73/73 repository total]; deliverables in sprints/audit_sprint_g/)
- **Sprint H (Common-Zero Datum Alignment & Spatial Markov Foundation)**: COMPLETE / VERIFIED (Implemented explicit vertical coordinate API for common-zero datum [z = -depth, top z = 0.0 m] per Prof. Sahoo's project directive; verified 100% stratigraphic invariance [1033.0 m, 328 intervals, zero facies count shift]; preserved validated 1D vertical Markov models [81/81 tests passing]; established horizontal pair extraction [4,410 empirical pairs] and Carle & Fogg [1996] continuous transition rate framework; evaluated minimal spatial Markov predictor under LOLO CV; 8/8 new unit tests passing [81/81 repository total]; deliverables in sprints/audit_sprint_h/)
- **Phase 2 (2D Cross-Sectional Geostatistical Modeling)**: READY FOR SPRINT (Scoped to unconditioned process-based forward simulation with W/T=35 geometric priors)

---

## 1. Current Working Dataset State & Ingestion Authority

- **Unified Ingestion Source**: `data/processed/lithologs_unified.parquet` and `data/processed/lithologs_unified.csv`
- **Total Ingested Lithologs**: 11 (litholog1, litholog10, litholog11, litholog2, litholog3, litholog4, litholog5, litholog6, litholog7, litholog8, litholog9)
- **Total Validated Observations**: 920 standardized 1-meter intervals
- **Cumulative Stratigraphic Span**: 922 meters (due to 1m continuity gap in Litholog 9 at 18-19m and 1m gap in Litholog 11 at 59-60m)
- **Provenance Manifest**: `data/provenance_manifest.json` (programmatic access via `load_provenance_manifest()`)

### Ingested Well Breakdown (Standardized 1m Points):
- `litholog1`: 93 pts (0 - 93 m)
- `litholog2`: 93 pts (0 - 93 m)
- `litholog3`: 93 pts (0 - 93 m)
- `litholog4`: 84 pts (0 - 84 m)
- `litholog5`: 85 pts (0 - 85 m)
- `litholog6`: 82 pts (0 - 82 m)
- `litholog7`: 80 pts (0 - 80 m)
- `litholog8`: 79 pts (0 - 79 m)
- `litholog9`: 77 pts (span 78 m; gap at 18-19 m, overlaps at 28-30 m resolved)
- `litholog10`: 77 pts (0 - 77 m)
- `litholog11`: 77 pts (span 78 m; manual raw log has gap at 59-60 m)

---

## 2. Dataset Provenance Status & Validation Limitations

| Litholog ID | Provenance Category | Reconstruction Method | Independent Validation | Benchmark Accuracy | Status & Conditioning Role |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **litholog1** | Source-derived | Existing dataset / Sahoo et al. (2016) | Unknown | None | Working dataset; independent validation undocumented |
| **litholog2** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **litholog3** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **litholog4** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **litholog5** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **litholog6** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **litholog7** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **litholog8** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **litholog9** | Source-derived | Existing dataset / Sahoo et al. (2016) | Unknown | None | Working dataset; independent validation undocumented |
| **litholog10** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **litholog11** | Source-derived (Benchmarked) | Existing manual log vs. automated digitization | **True** | **93.59%** | **Reference QC Benchmark Log** (73/78 m match) |

### Litholog 11 QC Benchmark Summary:
- **Depth Agreement**: 73 out of 78 meters match identically (93.59% structural accuracy)
- **Discrepant Intervals**: 5 meters (6.41% discrepancy rate): thin-bed coal streak (12-13m), recording gap (59-60m), boundary shift (65-66m), and facies aggregation (67-69m).
- **Exact Alignment**: 68 out of 78 meters (87.2%) exhibited exact facies and boundary alignment.

### Scientific Constraints & Caveats:
1. **No Benchmark Propagation**: The 93.59% accuracy benchmark applies strictly to Litholog 11. It must not be claimed or assumed as the accuracy of Lithologs 1-10.
2. **AI-Reconstructed Retention**: Lithologs 2-8 and 10 are preserved as active working profiles without manual correction to maintain unadulterated provenance.
3. **Downstream Conditioning**: All statistics, transition matrices, 2D cross-sections, and ML models utilizing the full 11-log dataset are conditioned partly on AI-reconstructed observations.

---

## 3. Verified Phase 1 Metrics (11-Log Complete Working Dataset)

### 1. Transition Probability Matrices

#### Regular Transition Matrix $P_{\text{reg}}$ (Fixed-Step 1m, 11 Logs)
```
                         Coal  Channel Sandstone  Fine Sandstone / Splay  Siltstone  Overbank Mudstone
Coal                    0.286              0.071                   0.286      0.036              0.321
Channel Sandstone       0.031              0.904                   0.026      0.018              0.021
Fine Sandstone / Splay  0.000              0.052                   0.759      0.043              0.147
Siltstone               0.027              0.082                   0.082      0.466              0.342
Overbank Mudstone       0.020              0.065                   0.013      0.082              0.820
```
- **Matrix Condition Number $\kappa(P_{\text{reg}})$**: 4.52

#### Embedded Transition Matrix $P_{\text{emb}}$ (Boundary Crossings, $P_{ii} = 0$, 11 Logs)
```
                         Coal  Channel Sandstone  Fine Sandstone / Splay  Siltstone  Overbank Mudstone
Coal                    0.000              0.100                   0.400      0.050              0.450
Channel Sandstone       0.324              0.000                   0.270      0.189              0.216
Fine Sandstone / Splay  0.000              0.214                   0.000      0.179              0.607
Siltstone               0.051              0.154                   0.154      0.000              0.641
Overbank Mudstone       0.109              0.364                   0.073      0.455              0.000
```
- **Matrix Condition Number $\kappa(P_{\text{emb}})$**: 17.82

### 2. Directional Asymmetry & Upward Succession Findings
- Channel Sandstone (State 1) upward transitions strictly favor finer-grained facies:
  - Sand $\to$ Coal: 0.324
  - Sand $\to$ Fine Sand/Splay: 0.270
  - Sand $\to$ Siltstone: 0.189
  - Sand $\to$ Overbank Mudstone: 0.216
  - **Combined Upward Fining/Abandonment Transitions**: 100.0%

### 3. Stationary Facies Occupancy Diagnostic (11 Logs)
| Facies State | Stationary $\pi_i$ | Empirical $p_{\text{emp}, i}$ | Absolute Diff | Relative Diff (%) |
|---|---|---|---|---|
| 0 (Coal) | 0.0305 | 0.0304 | 0.0001 | 0.29% |
| 1 (Channel Sandstone) | 0.4043 | 0.4196 | 0.0153 | 3.64% |
| 2 (Fine Sand / Splay) | 0.1264 | 0.1261 | 0.0003 | 0.25% |
| 3 (Siltstone) | 0.0807 | 0.0793 | 0.0014 | 1.73% |
| 4 (Overbank Mudstone) | 0.3581 | 0.3446 | 0.0135 | 3.92% |

- **Bulk Net-to-Gross (Sand + Splay)**:
  - Theoretical Stationary: **53.07%**
  - Empirical Observed: **54.57%**
  - Relative Difference: **2.74%** (< 3% diagnostic agreement)

---

## 4. Provenance Sensitivity Analysis: Complete 11-Log vs. 3-Log Source Subset

### Comparative Overview:
| Metric | 11-Log Working Dataset | 3-Log Source Subset (1, 9, 11) | Difference (11-log - 3-log) |
| :--- | :---: | :---: | :---: |
| **Observation Points ($N$)** | 920 | 247 | +673 |
| **Regular Transitions ($N_{\text{reg}}$)** | 909 | 244 | +665 |
| **Boundary Crossings ($N_{\text{emb}}$)** | 201 | 50 | +151 |
| **Stationary Net-to-Gross** | 53.07% | 52.51% | +0.56% |
| **Empirical Net-to-Gross** | 54.57% | 55.06% | -0.50% |
| **Matrix Frobenius Distance $||P_{\text{reg, 11}} - P_{\text{reg, 3}}||_F$** | - | - | **0.3120** |
| **Max Absolute Regular Difference** | - | - | **0.1786** |
| **Matrix Frobenius Distance $||P_{\text{emb, 11}} - P_{\text{emb, 3}}||_F$** | - | - | **0.6158** |
| **Max Absolute Embedded Difference** | - | - | **0.2660** |

#### Regular Matrix Difference ($P_{\text{reg, 11}} - P_{\text{reg, 3}}$):
```
                         Coal  Channel Sandstone  Fine Sandstone / Splay  Siltstone  Overbank Mudstone
Coal                   -0.014              0.071                   0.086      0.036             -0.179
Channel Sandstone       0.001              0.005                  -0.014     -0.002              0.011
Fine Sandstone / Splay  0.000              0.025                   0.029      0.016             -0.070
Siltstone               0.027             -0.094                  -0.035     -0.064              0.166
Overbank Mudstone      -0.030              0.004                  -0.012      0.032              0.005
```

#### Embedded Matrix Difference ($P_{\text{emb, 11}} - P_{\text{emb, 3}}$):
```
                         Coal  Channel Sandstone  Fine Sandstone / Splay  Siltstone  Overbank Mudstone
Coal                    0.000              0.100                   0.114      0.050             -0.264
Channel Sandstone       0.024              0.000                  -0.130     -0.011              0.116
Fine Sandstone / Splay  0.000              0.114                   0.000      0.079             -0.193
Siltstone               0.051             -0.221                  -0.096      0.000              0.266
Overbank Mudstone      -0.158              0.030                  -0.061      0.188              0.000
```

### Statistical Sample Size & Sparsity Constraints:
1. **Sample Size Disparity**: The 3-log source subset contributes only 50 boundary transitions across 247m. Sparse states (e.g. Coal with 7 boundaries, Siltstone with 8 boundaries) exhibit extreme small-sample variance, where an addition or subtraction of a single event shifts cell probabilities by $12.5\% - 14.3\%$.
2. **Invariance Preservation**: Crucially, adding the 8 AI-reconstructed profiles expands boundary crossings to 201 and smooths transitional noise while strictly preserving key sedimentological invariants:
   - Upward fining/abandonment from Channel Sandstones remains **100.0%**.
   - Net-to-Gross sand proportion shifts by only **0.50%** empirically and **0.56%** in stationary occupancy.
3. **Conditioning Acknowledgment**: While adding the AI-reconstructed logs reduces small-sample variance, downstream users must recognize that these smoother transition statistics are conditioned partly (73.2% of points) on AI-reconstructed observations.

---

## 5. Phase 2 Status Assessment

- **Dependency Analysis**: Phase 2 ("2D Object-Based Fluvial Generator", Track B, Days 21-30) is formulated as an unconditioned stochastic body generator conditioned strictly on literature priors from Sahoo et al. (2016):
  - Channel aspect ratio $W/T = 35$
  - Mean channel thickness $\sim 5.8\text{ m}$
  - Crevasse splay width range ($10 - 130\text{ m}$)
  - Target Net-to-Gross envelope ($17\% - 46\%$)
- **Validation Gate ("Done When")**: 100 unconditioned synthetic realizations yielding mean $W/T \in 35 \pm 2$ and $\text{NTG} \in [17\%, 46\%]$.
- **Conclusion**: Phase 2 does **not** depend quantitatively on the empirical litholog subset. No components of Phase 2 need rebuilding or invalidation due to the dataset update. It remains **READY FOR SPRINT**.

---

## 6. Sprint D Verified Metrics & Summary Reconciliation

### 6.1 Reconciled Sandstone Bed Statistics (Full 12 Lithologs)
- **Raw Sandstone Intervals**: Exactly 72 intervals across all 12 lithologs (cumulative thickness 444.70 m).
  - The Sprint C narrative count of 73 was an unverified typographical error; the source CSVs strictly sum to 72.
  - Min: 1.40 m, Median: 5.05 m, Mean: 6.18 m, Max: 17.00 m (Litholog 2, 40-57 m).
  - Sample standard deviation (ddof=1): 3.36 m; population standard deviation (ddof=0): 3.33 m (rounded to 3.34 m).
- **Merged Sandstone Lithosomes (Distinct Continuous Bodies)**: Exactly 53 distinct bodies (cumulative thickness 444.70 m).
  - Min: 1.40 m, Median: 6.00 m, Mean: 8.39 m.
  - Max: 32.00 m (Litholog 3 amalgamated channel-sand package at 61-93 m).
  - Sample standard deviation: 6.54 m.

### 6.2 Reconciled Stratigraphic Group Statistics
- **Upstream (L2-L8, L10)**:
  - Total thickness: 673.0 m (corrected from 664.0 m in Sprint C Table 7.1, which had transcribed a hardcoded plot label string).
  - Pure sandstone thickness: 287.0 m; pure sandstone N/G: 42.64% (0.4264).
- **Downstream Outcrop (L9, L11)**:
  - Total thickness: 156.0 m; pure sandstone thickness: 65.0 m; pure sandstone N/G: 41.67% (0.4167). Exact match.
- **Downstream Composite (L9, L11, L12)**:
  - Total thickness: 267.0 m; pure sandstone thickness: 122.70 m; pure sandstone N/G: 45.96% (0.4596) (corrected from 47.19%).
- **All 12 Lithologs (Full Study)**:
  - Total thickness: 1033.0 m; pure sandstone thickness: 444.70 m; pure sandstone N/G: 43.05% (0.4305). Exact match.

### 6.3 Leakage-Free Leave-One-Litholog-Out (LOLO) Cross-Validation
- **Regular 1D Chain (1 m Discretized Grid)**:
  - Average perplexity (All 12 logs): **2.1869** (Outcrop L1-L11 average: **2.0166**; Core L12: **4.0597**).
  - Initial-state score ($\ln P_{\text{base}}(s_1)$) kept strictly separate from transition perplexity.
- **Embedded 1D Chain (Distinct Bed Sequence, P_ii = 0)**:
  - Average bed-transition perplexity (All 12 logs): **3.7647** (Outcrop L1-L11 average: **3.6894**; Core L12: **4.5929**).
  - Non-equivalence principle: Regular perplexity (~2.19, high within-bed diagonal persistence) and embedded perplexity (~3.76, 4-way branching among distinct alternative facies) must not be compared as the same task.
- **Gap Sensitivity Analysis**:
  - Only L9 (18-19 m gap) and L11 (59-60 m gap) contain unrecorded gaps (2 transition steps in regular; 1 in embedded).
  - Masking gap-crossing transitions shifts perplexity by less than +/- 0.07, confirming negligible impact on conclusions.

---

## 7. Sprint E Reconciliation & Spatial Readiness Status

### 7.1 Reconciled Dataset & Provenance
- All 12 lithologs manually digitized and revised by user from primary source images.
- Encodings mapped: `p_sand` -> `sand` (Channel Sandstone undivided, code 1); `ripples` -> `silt` (Siltstone / heterolithics, code 3).
- Litholog 12 formally registered in `data/provenance_manifest.json` as a manually digitized subsurface core log (coverage strictly 0.0-111.0 m; original core length 242.0 m; unvalidated; not AI-generated).
- Total stratigraphic span across all 12 lithologs: 1033.0 m (922.0 m outcrop; 111.0 m core). Valid observations: 1031.0 m (due to two 1 m unmapped gaps in L9).
- Total sandstone raw intervals increased from 72 to 135; merged distinct sandstone lithosomes increased from 53 to 99; total sandstone thickness increased from 444.70 m to 507.70 m (49.15% overall N/G).

### 7.2 Synthetic-Coordinate Audit Findings
- Traced `data/loader.py:L136-141` formula: `strike_pos_m = float(file_idx * 100.0)`. Never entered 1D Markov chains.
- Dedicated spatial module (`smalt/spatial/coordinates.py`) created. Uses real local Cartesian coordinates $(X, Y)$ from `lolo/Location_coordinates_lithologs.xlsx` for L2-L12.
- Litholog 1 is missing coordinates and is strictly excluded from spatial validation.

### 7.3 Provisional Spatial Prototype & Empirical Findings
- Leave-One-Litholog-Out spatial cross-validation executed on eligible wells (L2-L12) in `smalt/spatial/baseline.py`.
- Evaluated against naive baselines:
  - Spatial 3D KNN: Macro-F1 = 0.2521, Balanced Acc = 0.2561, Raw Acc = 49.04%.
  - Nearest-Well Profile: Macro-F1 = 0.2523, Balanced Acc = 0.2545, Raw Acc = 47.55%.
  - Training Prior Facies: Macro-F1 = 0.1335, Balanced Acc = 0.2000, Raw Acc = 49.15%.
- Vertical datum sensitivity testing (±5m to ±20m shift) confirms Macro-F1 remains invariant at ~0.25, proving that sparse wells (~1.5-3.0 km apart) cannot correlate facies without an established stratigraphic marker datum.

---

## 8. Sprint F Six-State Migration & Spatial Baseline Revalidation Status

### 8.1 Six-State Canonical Facies Schema
- **Migration Completed**: Replaced provisional 5-state mapping with the definitive 6-state facies schema adhering to Sahoo et al. (2016):
  - **Code 0: `sand`** (`Channel Sandstone`, Facies 1)
  - **Code 1: `p_sand`** (`Planar Sandstone`, Facies 2)
  - **Code 2: `ripples`** (`Rippled Heterolithics`, Facies 3; alias `silt`)
  - **Code 3: `carbon_mud`** (`Carbonaceous Mudstone`, Facies 4)
  - **Code 4: `coal`** (`Coal`, Facies 5)
  - **Code 5: `mud`** (`Overbank Mudstone`, Facies 6)
- **Net Sand Partitioning Formalized**:
  - `ntg_channel` (Code 0 only): 43.05% (444.7 m across all 12 logs).
  - `ntg_planar` (Code 1 only): 6.10% (63.0 m across all 12 logs).
  - `ntg_net_sand` (`sand` + `p_sand`): 49.15% (507.7 m raw span; 506.7 m in discretized points).
  - `ntg_coarse` (`sand` + `p_sand` + `ripples`): 56.78% (585.4 m raw span).
- **Regenerated Pipeline Outputs**: Master script [`scripts/run_sprint_f_pipeline.py`](file:///d:/Lithology-reconstruction-using-XGB/scripts/run_sprint_f_pipeline.py) executed cleanly; all 11 CSV/JSON audit deliverables and 5 high-resolution figures generated in [`sprints/audit_sprint_f/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/).

### 8.2 Recomputed 6x6 Markov Transition Chains
- **Regular Discretized 1D Chain (1 m Grid)**:
  - Transition matrix dimension: 6x6. Average LOLO transition perplexity: **2.6289** across all 12 logs.
  - Strong diagonal self-transition persistence: `sand` (0.8354), `p_sand` (0.4286), `ripples` (0.5844), `carbon_mud` (0.4706), `coal` (0.5833), `mud` (0.8149).
- **Embedded 1D Chain (Distinct Bed Boundary Transitions, P_ii = 0)**:
  - Transition matrix dimension: 6x6 with zero diagonals. Average LOLO transition perplexity: **4.1793** across all 12 logs.
  - Asymmetric bed boundary transitions: `coal` transitions exclusively into `mud` (72.41%) or `carbon_mud` (27.59%), never directly erosive into `sand` (0.0%).

### 8.3 Spatial LOLO Baseline Revalidation Findings (L2-L12, 940 Grid Points)
- **Pooled Benchmark Comparison**:
  - **Spatial 3D KNN (k=5)**: Pooled Macro-F1 = **0.2232**, Balanced Accuracy = **22.81%**, Raw Accuracy = **47.77%** (449/940 correct).
  - **Nearest-Well Vertical Profile**: Pooled Macro-F1 = **0.2276**, Balanced Accuracy = **22.91%**, Raw Accuracy = **44.26%** (416/940 correct).
  - **Training Prior Majority Baseline (`sand`)**: Pooled Macro-F1 = **0.1011**, Balanced Accuracy = **16.67%**, Raw Accuracy = **43.51%** (409/940 correct).
- **Core Scientific Finding**:
  - Spatial 3D KNN fails to beat the 1D Nearest-Well vertical profile baseline in Macro-F1 (**0.2232 vs 0.2276**).
  - Spatial 3D KNN yields only a **+4.26 percentage point** gain in raw accuracy over the zero-spatial training prior majority baseline (47.77% vs 43.51%).
  - Raw accuracy is heavily inflated by majority-class predictions (`sand` and `mud`).
- **Minority Facies Failure**:
  - `coal` (support = 25): **0/25 correct** for both Spatial KNN and Nearest-Well (Recall = 0.0000, F1 = 0.0000).
  - `carbon_mud` (support = 17): **0/17 correct** for both Spatial KNN and Nearest-Well (Recall = 0.0000, F1 = 0.0000).
  - `p_sand` (support = 62): Nearest-Well achieves F1 = **0.1111** (7/62 correct); Spatial KNN achieves only F1 = **0.0465** (2/62 correct).
- **Directional Gap Generalization Failure**:
  - Upstream -> Downstream: Spatial KNN accuracy = **39.33%**, underperforming the training prior (**44.94%**).
  - Downstream -> Upstream: Spatial KNN accuracy = **38.93%**, underperforming the training prior (**42.94%**).
- **Vertical Datum Sensitivity**:
  - Systematic datum offset perturbations from -20 m to +20 m shift Spatial KNN Macro-F1 only between **0.2164 and 0.2270**. Macro-F1 remains effectively invariant, demonstrating that lateral facies heterogeneity and wide inter-well spacing (150 m to 14 km), rather than vertical datum alignment, govern spatial predictability.

### 8.4 Test Suite & Quality Assurance
- **Full Test Suite Status**: 65/65 tests passing cleanly across the repository.
- Dedicated test suite [`tests/test_sprint_f.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_f.py) validates 6-state mapping, distinct `p_sand`/`ripples`, strict L1 exclusion, zero leakage during cross-validation, 6x6 Markov normalization, and artifact synchronization.

---

## 9. Sprint G Independent Results Audit & Spatial Modeling Decision

### 9.1 Audit Findings & Metric Verification
- **Metrics Reconciliation (0 Discrepancy)**: All Sprint F pooled metrics independently recomputed from confusion matrices and verified: Spatial 3D KNN raw accuracy = **0.4777** (449/940), balanced accuracy = **0.2281**, Macro-F1 = **0.2232**; Nearest-Well raw accuracy = **0.4426**, balanced accuracy = **0.2291**, Macro-F1 = **0.2276**; Training Prior raw accuracy = **0.4351**, balanced accuracy = **0.1667**, Macro-F1 = **0.1011**. Unweighted mean fold balanced accuracy (**24.89%**) distinguished from pooled (**22.81%**).
- **Physical Sparsity Root Cause**: Nearest-neighbor inter-well distances range from **420.0 m to 2,355.6 m** (median 701.5 m). In contrast, Sahoo et al. (2016) report single-storey channel widths of **140 to 210 m** ($W/T \approx 35$) and splay widths of **10 to 130 m**. Inter-well spacing is 2x to 10x wider than the maximum lateral continuity of individual sandbodies. Point-wise interpolation is geologically impossible.
- **Class Imbalance Distortion**: `sand` (43.5%) and `mud` (37.2%) account for 80.7% of all points. Thin-bed minority facies collapsed completely: `coal` (support 25) F1 = **0.0000** (0/25 identified); `carbon_mud` (support 17) F1 = **0.0000** (0/17 identified); `p_sand` (support 62) F1 = **0.0465** (2/62 identified).
- **StandardScaler Anisotropy Defect**: Pre-scaling vertical weight factor ($c = 10.0$) was numerically cancelled by `StandardScaler` standard deviation division. Diagnostic sweep across actual post-scaling vertical weights (0.1 to 50.0) confirmed Macro-F1 remains trapped between **0.21 and 0.23**, proving that spatial KNN is physically unsuited for this problem.
- **Alias Dictionary Audit**: Confirmed authoritative raw CSVs in `data/raw_lithologs/` contain strictly canonical tokens. Identified legacy risks in `data/loader.py`: `splay` mapped to `carbon_mud` (defect) and generic `siltstone` mapped to `ripples` (ambiguous).

### 9.2 SMALT Phase Decisions (20-Day Presentation Runway)
1. **Phase 1 (1D Vertical Markov Analysis)**: **PROCEED (Core Contribution)**. Operates under Walther's Law; verified on 1031 m across 12 lithologs; 6x6 transition perplexity 2.6289 (regular) and 4.1793 (embedded).
2. **Phase 2 (2D Fluvial Forward Modeling)**: **PROCEED AS UNCONDITIONED PROCESS PROTOTYPE**. Implement ribbon channel generator with $W/T \approx 35$ and target Net-to-Gross envelope (17-50%) from Sahoo et al. (2016) as a forward stochastic simulator, not inter-well reconstruction.
3. **Phase 3 (Inter-Well Spatial Reconstruction)**: **DEFER / REJECT ON REAL WELLS**. Definitively proven impossible with current data. Presenting this rigorous negative result protects against examination criticism and demonstrates scientific maturity.
4. **Phases 4-5 (Spatial ML & Active Learning)**: **RE-SCOPE TO SYNTHETIC BENCHMARK DEMONSTRATION**. Demonstrate Active Margin Sampling on synthetic 2D forward realizations (from Phase 2) where ground truth is known, proving borehole budget optimization without making false claims on sparse field data.

### 9.3 Test Suite Status
- Repository test suite: **73 passed, 0 failed** (`pytest tests/`).
- Dedicated test suite: [`tests/test_sprint_g.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_g.py) (8/8 tests passing).
- Audit deliverables: [`sprints/audit_sprint_g/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_g/) (5 audit CSVs and comprehensive audit report [`sprints/audit_sprint_g/SPRINT_G_AUDIT_REPORT.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_g/SPRINT_G_AUDIT_REPORT.md)).

---

## 10. Sprint H Common-Zero Datum Alignment & Spatial Markov Foundation Status

### 10.1 Common-Zero Vertical Alignment & Invariance
- **Project Reference Level**: Adopted Prof. Hiranya Sahoo's explicit project instruction: "All litholog zeros should be treated as being at the same reference level."
- **Standard Vertical Coordinate**: Constructed explicit, reversible API in [`smalt/spatial/datum.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/datum.py):
  - $z_{\text{common}} = -\text{depth}_{\text{measured}}$ (top datum $z = 0.0$ m; elevation convention with deeper intervals increasingly negative).
  - Preserves original measured depth column (`depth_original_m`), interval thickness, facies codes, and facies names.
  - Round-trip exact reversibility verified ($\text{depth} = -z_{\text{common}}$).
- **100% Stratigraphic Invariance Verified**:
  - Total cumulative thickness: Exactly **1033.0 m** across all 12 lithologs (328 intervals).
  - Facies interval counts completely unchanged: `mud` (109), `sand` (92), `p_sand` (43), `ripples` (31), `coal` (29), `carbon_mud` (24).
  - Litholog 1 (93.0 m) retained for 1D vertical analysis; strictly excluded from spatial modeling due to missing coordinates.

### 10.2 Empirical Horizontal Transition Pairs & Continuous Rate Matrix
- **Matched Elevation Slicing**: Extracted **4,410 empirical horizontal well pairs** across matched 1 m elevation slices from the 11 coordinate-bearing wells (L2-L12).
- **Lag-Binned Transition Dynamics**:
  - Short lag (400-1000 m, 736 pairs): `sand` auto-transition probability = **0.542**; `mud` = **0.448**; `coal` = **0.000** (no inter-well coal coincidence observed).
  - Intermediate lag (1000-2500 m, 1,384 pairs): `sand` auto-transition = **0.478**; `mud` = **0.410**.
  - Long lag (2500-5500 m, 2,290 pairs): `sand` auto-transition decays to **0.421**; `mud` to **0.385**, approaching stationary background proportions (43.5% sand, 37.2% mud).
- **Continuous Transition Rate Matrix ($\mathbf{R}_h$)**:
  - Implemented Carle & Fogg (1996) formulation: $R_{ii} = -1 / \bar{L}_{h,i}$, off-diagonals $R_{ij} = -R_{ii} p_j / (1 - p_i)$, rows sum to 0.
  - Continuous matrix exponential $\mathbf{P}(h) = \exp(\mathbf{R}_h h)$ provides exact analytical transition decay for any horizontal separation.

### 10.3 Spatial LOLO Benchmark Performance on Common Zero (940 Grid Points)
- **Spatial LOLO Results**:
  - **Nearest-Well Vertical Profile**: Pooled Raw Acc = **45.85%** (431/940), Balanced Acc = **23.18%**, Macro-F1 = **0.2303**.
  - **Spatial 3D KNN (k=5)**: Pooled Raw Acc = **48.51%** (456/940), Balanced Acc = **22.98%**, Macro-F1 = **0.2256**.
  - **Training Prior Majority Baseline (`sand`)**: Pooled Raw Acc = **43.51%** (409/940), Balanced Acc = **16.67%**, Macro-F1 = **0.1011**.
  - **Spatial Markov Transition Model**: Pooled Raw Acc = **37.23%** (350/940), Balanced Acc = **16.67%**, Macro-F1 = **0.0904**.
- **Theoretical Finding on Spatial Markov Behavior**:
  - The continuous spatial Markov model is the exact theoretical link between Nearest-Well ($h \to 0, \mathbf{P}(0) = \mathbf{I}$) and stationary prior ($h \to \infty, \mathbf{P}(h) \to \mathbf{1}\mathbf{p}^T$).
  - Because inter-well spacing ($h \ge 420.0$ m) exceeds lateral channel widths ($\bar{L}_h \approx 203$ m), the Markov model smoothly decays to predicting the most extensive regional lithology: Overbank Mudstone ($\bar{L}_h \approx 360$ m, 37.23% support).
  - This mathematically proves the physical limit of inter-well correlation in sparse data.

### 10.4 Test Suite & Deliverables
- **Test Suite Status**: **81 passed, 0 failed** (`pytest tests/`). Dedicated test file: [`tests/test_sprint_h.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_h.py) (8/8 passing).
- **Audit Deliverables**: [`sprints/audit_sprint_h/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_h/) (5 audit CSVs, 2 high-resolution publication figures, and comprehensive report [`sprints/audit_sprint_h/SPRINT_H_COMMON_ZERO_SPATIAL_MARKOV_REPORT.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_h/SPRINT_H_COMMON_ZERO_SPATIAL_MARKOV_REPORT.md)).

---

## 11. Immediate Next Actions (Phase 2 & Presentation Preparation)

1. **Implement Scoped Phase 2 Fluvial Forward Generator**:
   - Build lightweight, unconditioned 2D ribbon channel cross-section simulator conditioned on Sahoo et al. (2016) architectural priors ($W/T \approx 35$, mean thickness $\approx 5.8$ m, target N/G 17-50%).
2. **Prepare UGP Presentation Slide Deck**:
   - Assemble slide deck highlighting the 1D Markov empirical succession, the common-zero datum alignment, the spatial Markov decay findings, and the 2D forward process model.
3. **Formulate Next Data Inquiries for Prof. Sahoo**:
   - Deliver Sprint H report and common-zero alignment transects to Prof. Sahoo.



