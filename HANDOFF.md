# SMALT Project Handoff & Phase Registry

## Phase Status Overview
- **Phase 0 (Data Ingestion & Quality Control)**: COMPLETE / VERIFIED (11/11 lithologs loaded cleanly, 920 observations, 11/11 unit tests passing, provenance manifest registered, docs/notes/phase0_schema.md published)
- **Phase 1 (1D Vertical Markov Succession Analysis)**: VERIFIED (Full 11-log dataset fitted [920 observations], provenance sensitivity analysis documented, 11/11 unit tests passing)
- **Phase 2 (2D Cross-Sectional Geostatistical Modeling)**: READY FOR SPRINT (Validated solely against Sahoo et al. [2016] geological priors [W/T = 35, mean thickness 5.8m, NTG 17%-46%] and synthetic realizations; independent of litholog subset)

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

## 6. Recommended Next Action

1. **Sprint Phase 2**: Implement `smalt/geostat/object_sim.py` according to Sahoo et al. (2016) geometrical priors.
