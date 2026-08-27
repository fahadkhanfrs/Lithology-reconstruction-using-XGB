# SMALT Project Handoff & Phase Registry

## Phase Status Overview
- **Phase 0 (Data Ingestion & Quality Control)**: VERIFIED (10/10 tests passing)
- **Phase 1 (1D Vertical Markov Succession Analysis)**: VERIFIED (11/11 tests passing)
- **Phase 2 (2D Cross-Sectional Geostatistical Modeling)**: READY FOR SPRINT

---

## Verified Phase 1 Metrics & Acceptance Summary

### 1. Dataset & Input Authority
- **Source**: `data/processed/lithologs_unified.parquet`
- **Total Validated Observations**: 247
- **Wells Fitted**: 3 (litholog1, litholog11, litholog9)
- **Total Regular Transitions**: 244
- **Total Boundary Crossing Transitions**: 50

### 2. State Space & Facies Encoding (5-State SMALT Aggregation)
- `0`: Coal
- `1`: Channel Sandstone
- `2`: Fine Sandstone / Splay
- `3`: Siltstone
- `4`: Overbank Mudstone

### 3. Transition Probability Matrices

#### Regular Transition Matrix $P_{\text{reg}}$ (Fixed-Step 1m)
```
                         Coal  Channel Sandstone  Fine Sandstone / Splay  Siltstone  Overbank Mudstone
Coal                    0.300              0.000                   0.200      0.000              0.500
Channel Sandstone       0.030              0.899                   0.040      0.020              0.010
Fine Sandstone / Splay  0.000              0.027                   0.730      0.027              0.216
Siltstone               0.000              0.176                   0.118      0.529              0.176
Overbank Mudstone       0.049              0.062                   0.025      0.049              0.815
```
- **Matrix Condition Number $\kappa(P_{\text{reg}})$**: 4.98

#### Embedded Transition Matrix $P_{\text{emb}}$ (Boundary Crossings, $P_{ii} = 0$)
```
                         Coal  Channel Sandstone  Fine Sandstone / Splay  Siltstone  Overbank Mudstone
Coal                    0.000              0.000                   0.286      0.000              0.714
Channel Sandstone       0.300              0.000                   0.400      0.200              0.100
Fine Sandstone / Splay  0.000              0.100                   0.000      0.100              0.800
Siltstone               0.000              0.375                   0.250      0.000              0.375
Overbank Mudstone       0.267              0.333                   0.133      0.267              0.000
```
- **Matrix Condition Number $\kappa(P_{\text{emb}})$**: 417.25

### 4. Asymmetry & Upward Succession Findings
- Channel Sandstone (State 1) upward transitions strictly favor finer-grained facies:
  - Sand $\to$ Coal: 0.300
  - Sand $\to$ Fine Sand/Splay: 0.400
  - Sand $\to$ Siltstone: 0.200
  - Sand $\to$ Overbank Mudstone: 0.100
  - **Combined Upward Fining/Abandonment Transitions**: 100.0%

### 5. Stationary Facies Occupancy Diagnostic
| Facies State | Stationary $\pi_i$ | Empirical $p_{\text{emp}, i}$ | Absolute Diff | Relative Diff (%) |
|---|---|---|---|---|
| 0 (Coal) | 0.0423 | 0.0405 | 0.0018 | 4.55% |
| 1 (Channel Sandstone) | 0.3762 | 0.4008 | 0.0246 | 6.13% |
| 2 (Fine Sand / Splay) | 0.1489 | 0.1498 | 0.0009 | 0.60% |
| 3 (Siltstone) | 0.0634 | 0.0688 | 0.0054 | 7.83% |
| 4 (Overbank Mudstone) | 0.3691 | 0.3401 | 0.0290 | 8.53% |

- **Bulk Net-to-Gross (Sand + Splay)**:
  - Theoretical Stationary: **52.51%**
  - Empirical Observed: **55.06%**
  - Relative Difference: **4.63%** (< 5% diagnostic agreement)

### 6. Phase 1 Artifact Registry
- Technical Note: `docs/notes/phase1_markov.md`
- Engine Module: `smalt/geostat/markov.py`
- Visualization Module: `smalt/viz/markov_viz.py`
- Pytest Suite: `tests/test_phase1_markov.py`
- Summary JSON: `results/phase1_markov_summary.json`
- Figures:
  - `docs/figures/phase1_transition_matrix_regular.png`
  - `docs/figures/phase1_transition_matrix_embedded.png`
  - `docs/figures/phase1_stationary_vs_empirical.png`
  - `docs/figures/phase1_facies_succession_network.png`
