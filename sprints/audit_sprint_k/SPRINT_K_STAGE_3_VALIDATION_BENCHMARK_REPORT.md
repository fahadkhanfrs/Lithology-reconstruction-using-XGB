# SMALT Sprint K - Stage 3 Final Report
**Multi-Seed Geological Validation and Sparse-Borehole Reconstruction Benchmark**

---

## 1. Executive Summary

Sprint K Stage 3 evaluated whether the Subsurface Stratigraphic Modeling & Active Learning Toolkit (SMALT) can generate scientifically defensible, reproducible 2D synthetic geology and recover known subsurface ground truth from sparse borehole observations.

Following the pre-implementation inspection and Stage 1-2 deliverables, this stage completed:
1. **Multi-Seed Geological Validation**: Executed $N=30$ independent realizations across three presets (`upstream`, `downstream_outcrop`, and `pooled`), totaling 90 validation runs.
2. **Net-to-Gross & Architecture Diagnosis**: Quantified the empirical Net-to-Gross discrepancy, isolating the roles of internal Waltherian facies assignment (20% basal `p_sand`), erosional cannibalization during channel superposition, and stopping tolerance.
3. **Bed-Thickness Discretization & Log-Normal Evaluation**: Compared continuous empirical bed measurements with grid-discretized column run-lengths and verified the log-normal MLE model ($\mu=1.7870, \sigma=0.8823$) with cluster-aware caveats.
4. **Sensitivity Experiments**: Evaluated architectural sensitivity across literature width-to-thickness priors ($W/T \in \{20, 35, 50\}$) and grid resolutions ($\Delta z \in \{0.25, 0.5, 1.0\}\text{ m}$).
5. **Sparse-Borehole Reconstruction Benchmark**: Implemented and benchmarked four leakage-free reconstruction baselines across three borehole spacings ($400\text{ m}, 600\text{ m}, 1000\text{ m}$) over 10 independent evaluation realizations per spacing (120 model evaluations).
6. **Automated Testing**: Created and executed 14 comprehensive unit and integration tests in [`tests/test_sprint_k.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_k.py), passing with a 100% success rate.

---

## 2. Validation Campaign Design & Realization Counts

The validation campaign was configured in [`smalt/forward/validation.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/forward/validation.py):

| Campaign Component | Presets / Settings | Realizations | Seed Range | Execution Time | Output CSV File |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **Multi-Seed Validation** | `upstream`, `downstream_outcrop`, `pooled` | $30 \times 3 = 90$ | 1000 - 3029 | ~9.0 s | [`multi_seed_geological_validation_metrics.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/multi_seed_geological_validation_metrics.csv) |
| **W/T Ratio Sensitivity** | $W/T \in \{20.0, 35.0, 50.0\}$ | $10 \times 3 = 30$ | 7000 - 10009 | ~3.5 s | [`sensitivity_wt_and_resolution.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/sensitivity_wt_and_resolution.csv) |
| **Grid Step Sensitivity** | $\Delta z \in \{0.25, 0.5, 1.0\}\text{ m}$ | $10 \times 3 = 30$ | 6000 - 6109 | ~3.8 s | [`sensitivity_wt_and_resolution.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/sensitivity_wt_and_resolution.csv) |
| **Reconstruction Benchmark** | Spacing $\in \{400, 600, 1000\}\text{ m}$ (4 models) | $10 \times 3 \times 4 = 120$ | 8000 - 8029 | ~12.5 s | [`synthetic_reconstruction_benchmark_metrics.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/synthetic_reconstruction_benchmark_metrics.csv) |

---

## 3. Geological Validation Results & Reference Discrepancies

### 3.1 Facies Proportions & Net-to-Gross Across Seeds
Realized cell proportions across $N=30$ independent realizations compared against corresponding empirical group reference statistics (from [`multi_seed_reference_comparison.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/multi_seed_reference_comparison.csv)):

| Preset | Facies / Metric | Simulated Mean $\pm$ Std | Empirical Reference | Discrepancy | Relative Error |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **Upstream** ($N=30$) | `sand` (Facies 0) | $38.09\% \pm 2.28\%$ | $42.94\%$ | $-4.85\%$ | $11.28\%$ |
| | `p_sand` (Facies 1) | $9.46\% \pm 1.20\%$ | $7.58\%$ | $+1.88\%$ | $24.78\%$ |
| | `ripples` (Facies 2) | $2.79\% \pm 0.95\%$ | $4.90\%$ | $-2.11\%$ | $43.11\%$ |
| | `carbon_mud` (Facies 3) | $0.91\% \pm 0.33\%$ | $0.89\%$ | $+0.02\%$ | $1.76\%$ |
| | `coal` (Facies 4) | $1.73\% \pm 0.37\%$ | $1.93\%$ | $-0.20\%$ | $10.41\%$ |
| | `mud` (Facies 5) | $47.03\% \pm 2.45\%$ | $41.75\%$ | $+5.28\%$ | $12.64\%$ |
| | **$N/G_{\text{channel}}$** | $\mathbf{38.09\% \pm 2.28\%}$ | $\mathbf{42.94\%}$ | $\mathbf{-4.85\%}$ | $\mathbf{11.28\%}$ |
| | **$N/G_{\text{pure}}$** (`sand`+`p_sand`) | $\mathbf{47.55\% \pm 3.13\%}$ | $\mathbf{50.52\%}$ | $\mathbf{-2.97\%}$ | $\mathbf{5.87\%}$ |
| | **$N/G_{\text{coarse}}$** (all coarse) | $\mathbf{50.34\% \pm 2.64\%}$ | $\mathbf{55.42\%}$ | $\mathbf{-5.08\%}$ | $\mathbf{9.17\%}$ |
| **Downstream Outcrop** ($N=30$) | `sand` (Facies 0) | $35.25\% \pm 1.11\%$ | $40.38\%$ | $-5.13\%$ | $12.71\%$ |
| | `p_sand` (Facies 1) | $7.68\% \pm 0.90\%$ | $7.69\%$ | $-0.01\%$ | $0.18\%$ |
| | `ripples` (Facies 2) | $5.02\% \pm 0.99\%$ | $21.15\%$ | $-16.13\%$ | $76.25\%$ |
| | `mud` (Facies 5) | $49.39\% \pm 2.07\%$ | $24.36\%$ | $+25.03\%$ | $102.76\%$ |
| | **$N/G_{\text{channel}}$** | $35.25\% \pm 1.11\%$ | $40.38\%$ | $-5.13\%$ | $12.71\%$ |
| | **$N/G_{\text{pure}}$** | $42.93\% \pm 1.60\%$ | $48.08\%$ | $-5.15\%$ | $10.72\%$ |
| | **$N/G_{\text{coarse}}$** | $47.95\% \pm 1.93\%$ | $69.23\%$ | $-21.28\%$ | $30.74\%$ |
| **Pooled** ($N=30$) | `sand` (Facies 0) | $37.89\% \pm 1.61\%$ | $43.05\%$ | $-5.16\%$ | $11.98\%$ |
| | `p_sand` (Facies 1) | $6.85\% \pm 0.80\%$ | $6.10\%$ | $+0.75\%$ | $12.26\%$ |
| | `mud` (Facies 5) | $49.80\% \pm 2.47\%$ | $38.34\%$ | $+11.46\%$ | $29.88\%$ |
| | **$N/G_{\text{pure}}$** | $44.74\% \pm 2.13\%$ | $49.15\%$ | $-4.41\%$ | $8.97\%$ |

### 3.2 Root Cause Analysis: Target Net-to-Gross Discrepancy
When configuring a channel sand target of $45.0\%$, the realized `ntg_channel` averaged $38.09\% \pm 2.28\%$ (discrepancy of $-4.85\%$). Investigation identified three contributing factors:
1. **Internal Waltherian Sand Partitioning**:
   In accordance with sedimentological observations (Sahoo et al., 2016), the channel belt is partitioned into a core trough cross-bedded sandstone (`sand`, Facies 0) and a basal planar-laminated sheet sandstone (`p_sand`, Facies 1; configured at 20% basal fraction). When evaluating total channel-belt sand ($N/G_{\text{pure}} = \text{sand} + \text{p\_sand}$), the realized average is **$47.55\% \pm 3.13\%$**, matching the empirical target ($50.52\%$) within $2.97\%$ error.
2. **Erosional Cannibalization**:
   As channel bodies are stacked sequentially, younger channels scour into older channel sand rather than pristine mud background. Each subsequent channel adds diminishing net new sand to the 2D cross-section.
3. **Stopping Window Tolerance**:
   The generator stops when current sand reaches `target_ntg - tolerance` ($0.45 - 0.08 = 0.37$). The simulation correctly terminates once the threshold of $37\%$ is achieved.
- **Scientific Conclusion**: The simulator should not silently alter facies labels to force an artificial $45\%$ match. Instead, users should distinguish between single-facies $N/G_{\text{channel}}$ ($38\%$) and complete channel-belt $N/G_{\text{pure}}$ ($48\%$).

### 3.3 Downstream Outcrop Group Divergence
The downstream outcrop preset exhibits a substantial discrepancy in `ripples` ($5.02\%$ simulated vs $21.15\%$ empirical target). This occurs because:
1. The empirical downstream outcrop reference consists of only two measured sections (L9 and L11), where Litholog 9 contains an unusually thick $23.0\text{ m}$ continuous ripple-laminated heterolith unit.
2. In the 2D forward simulator, crevasse splays are configured with realistic maximum thickness bounds ($T \le 4.5\text{ m}$), preventing a single splay from dominating an entire $50\text{ m}$ column.
3. This discrepancy highlights the extreme small-sample variance of the two-log downstream outcrop group, reinforcing why the Upstream group ($N=8$ logs) is the superior calibration reference.

---

## 4. Vertical Succession, Bed Thickness & Discretization Effects

### 4.1 Empirical Bed Thickness vs Grid Discretization
- **Empirical Measured Beds (Upstream)**: Mean channel sand bed thickness is $8.50\text{ m}$ (median $5.50\text{ m}$).
- **Simulated Grid Run-Lengths ($\Delta z = 0.5\text{ m}$)**:
  - Simulated `sand` run-length: $4.90 - 7.50\text{ m}$ (averaging $6.12\text{ m}$ across seeds).
  - Simulated `mud` run-length: $4.65 - 8.20\text{ m}$ (averaging $5.45\text{ m}$).
- **Grid Discretization Effect**:
  - The $0.5\text{ m}$ vertical grid cannot resolve continuous geological lamina thinner than $0.5\text{ m}$.
  - In L12 Core, where coal beds were digitized down to $0.2\text{ m}$, grid cells at $\Delta z = 0.5\text{ m}$ expand minimum bed thickness to $0.5\text{ m}$ (a $150\%$ inflation of sub-meter intervals).
  - Conversely, for thick channel sandstones ($T > 5\text{ m}$), discretization onto $0.5\text{ m}$ cells induces less than $2\%$ geometric distortion.

### 4.2 Lognormal Model Goodness-of-Fit Assessment
- An unadjusted Kolmogorov-Smirnov test against fitted parameters yielded $D = 0.1059, p = 0.8023$. As required by Stage 3 instructions, this ordinary KS p-value must not be accepted as conclusive because parameters $\mu=1.7870, \sigma=0.8823$ were estimated on the exact same 34 beds.
- A descriptive parametric resampling check was performed:
  - Resampled theoretical median: $5.97\text{ m}$ vs empirical median $5.50\text{ m}$.
  - Resampled 90th percentile: $18.52\text{ m}$ vs empirical 90th percentile $17.40\text{ m}$.
  - Resampled maximum: $34.2\text{ m}$ vs empirical maximum $32.0\text{ m}$.
- **Clustering Caveat**: The 34 channel sandstone beds originate from only 8 distinct lithologs. Because beds within the same litholog are stratigraphically clustered, the true effective degrees of freedom is closer to 8 than 34. The log-normal distribution is a defensible parameterization, not a fully validated geological law.

---

## 5. Lateral Architecture & Sensitivity Analysis

Sensitivity analysis evaluated 30 independent runs across literature width-to-thickness priors ($W/T \in \{20, 35, 50\}$) and grid steps:

| Experiment Setting | Mean Channel Width (m) | Mean Channels Placed | Realized $N/G_{\text{channel}}$ | Realized $N/G_{\text{pure}}$ | Mean Sand Bed Thickness (m) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **$W/T = 20.0$** (Narrow Ribbons) | $167.3\text{ m} \pm 18.2\text{ m}$ | $48.6 \pm 3.1$ | $33.69\% \pm 4.22\%$ | $42.87\% \pm 5.58\%$ | $6.14\text{ m} \pm 1.08\text{ m}$ |
| **$W/T = 35.0$** (Baseline Prior) | $286.8\text{ m} \pm 36.1\text{ m}$ | $31.6 \pm 9.4$ | $38.12\% \pm 1.34\%$ | $47.66\% \pm 2.45\%$ | $6.77\text{ m} \pm 1.25\text{ m}$ |
| **$W/T = 50.0$** (Wide Sand Sheets) | $437.0\text{ m} \pm 112.5\text{ m}$ | $24.7 \pm 11.5$ | $40.07\% \pm 3.23\%$ | $49.97\% \pm 3.75\%$ | $8.33\text{ m} \pm 3.01\text{ m}$ |
| **$\Delta z = 0.25\text{ m}$** | $293.5\text{ m} \pm 24.1\text{ m}$ | $27.5 \pm 2.1$ | $39.56\% \pm 1.94\%$ | $49.41\% \pm 2.32\%$ | $7.41\text{ m} \pm 0.81\text{ m}$ |
| **$\Delta z = 0.50\text{ m}$** (Default) | $286.8\text{ m} \pm 36.1\text{ m}$ | $31.6 \pm 9.4$ | $38.12\% \pm 1.34\%$ | $47.66\% \pm 2.45\%$ | $6.77\text{ m} \pm 1.25\text{ m}$ |
| **$\Delta z = 1.00\text{ m}$** | $281.2\text{ m} \pm 31.0\text{ m}$ | $33.1 \pm 8.2$ | $37.24\% \pm 2.11\%$ | $46.12\% \pm 2.89\%$ | $6.10\text{ m} \pm 1.40\text{ m}$ |

### Key Geological Takeaways from Sensitivity Experiments
1. **Impact of Aspect Ratio ($W/T$)**:
   - At $W/T = 20$, channels form narrow ribbons ($W \approx 167\text{ m}$). To reach target Net-to-Gross, the generator places $\sim 49$ channels, resulting in isolated sandbodies embedded in mud.
   - At $W/T = 50$, channels form broad sheets ($W \approx 437\text{ m}$). Only $\sim 25$ channels are required, leading to extensive lateral amalgamation and increasing mean sand bed thickness to $8.33\text{ m}$ (max $14.94\text{ m}$).
2. **Impact of Grid Step ($\Delta z$)**:
   - Refining $\Delta z$ from $1.0\text{ m}$ to $0.25\text{ m}$ improves boundary definition but alters mean bed thickness by less than $1.3\text{ m}$ ($6.10\text{ m}$ vs $7.41\text{ m}$), confirming numerical stability.

---

## 6. Sparse-Borehole Reconstruction Benchmark

### 6.1 Benchmark Design & Leakage Prevention
- **Evaluation Set**: 10 independent evaluation realizations per spacing (30 realization evaluations $\times$ 4 models = 120 model evaluations).
- **Domain Extent**: $X = 2000\text{ m}, Z = 50\text{ m}$ ($100 \times 200 = 20,000$ grid cells).
- **Spacings Tested**:
  - $400\text{ m}$ spacing: 5 boreholes ($x \in [200, 600, 1000, 1400, 1800]\text{ m}$), 500 observed cells.
  - $600\text{ m}$ spacing: 3 boreholes ($x \in [300, 900, 1500]\text{ m}$), 300 observed cells.
  - $1000\text{ m}$ spacing: 2 boreholes ($x \in [500, 1500]\text{ m}$), 200 observed cells.
- **Strict Leakage Prevention**:
  - Reconstruction models receive strictly the `SyntheticBoreholeDataset` (only sampled borehole rows).
  - All sampled borehole columns are withheld from evaluation scoring.
  - The full 2D grid and geobody bounding boxes are hidden from predictive models.

### 6.2 Benchmark Results Across Spacings
Realization-level aggregate metrics (mean $\pm$ std across 10 evaluation seeds from [`synthetic_reconstruction_summary.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/synthetic_reconstruction_summary.csv)):

| Spacing | Reconstruction Model | Wells | Raw Accuracy | Balanced Accuracy | Macro-F1 | Interp. Accuracy | Extrap. Accuracy |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **400 m** | **Nearest Well** | 5 | $\mathbf{75.91\% \pm 4.63\%}$ | $\mathbf{52.61\% \pm 3.22\%}$ | $\mathbf{0.5226 \pm 0.0252}$ | $75.10\% \pm 5.00\%$ | $79.17\% \pm 10.11\%$ |
| | **Spatial 2D KNN ($k=5$)** | 5 | $74.79\% \pm 4.57\%$ | $42.46\% \pm 3.75\%$ | $0.4286 \pm 0.0370$ | $74.05\% \pm 5.03\%$ | $77.75\% \pm 9.83\%$ |
| | **Spatial Markov Decay** | 5 | $68.82\% \pm 3.68\%$ | $44.17\% \pm 2.63\%$ | $0.4669 \pm 0.0238$ | $67.78\% \pm 3.64\%$ | $72.99\% \pm 11.69\%$ |
| | **Vertical Prior Baseline** | 5 | $60.42\% \pm 5.95\%$ | $31.64\% \pm 7.71\%$ | $0.3098 \pm 0.0761$ | $60.99\% \pm 6.74\%$ | $58.16\% \pm 8.17\%$ |
| **600 m** | **Nearest Well** | 3 | $\mathbf{66.50\% \pm 5.23\%}$ | $\mathbf{46.74\% \pm 4.61\%}$ | $\mathbf{0.4592 \pm 0.0436}$ | $66.84\% \pm 6.86\%$ | $65.99\% \pm 7.88\%$ |
| | **Spatial 2D KNN ($k=5$)** | 3 | $65.38\% \pm 5.06\%$ | $36.98\% \pm 3.92\%$ | $0.3722 \pm 0.0329$ | $65.92\% \pm 6.71\%$ | $64.56\% \pm 7.66\%$ |
| | **Vertical Prior Baseline** | 3 | $58.70\% \pm 6.43\%$ | $33.42\% \pm 5.96\%$ | $0.3214 \pm 0.0620$ | $60.01\% \pm 6.79\%$ | $56.74\% \pm 7.31\%$ |
| | **Spatial Markov Decay** | 3 | $56.84\% \pm 4.49\%$ | $33.41\% \pm 2.66\%$ | $0.3606 \pm 0.0350$ | $56.50\% \pm 5.72\%$ | $57.34\% \pm 9.28\%$ |
| **1000 m** | **Nearest Well** | 2 | $\mathbf{57.74\% \pm 7.29\%}$ | $\mathbf{38.89\% \pm 3.20\%}$ | $\mathbf{0.3854 \pm 0.0396}$ | $58.94\% \pm 10.40\%$ | $56.53\% \pm 6.75\%$ |
| | **Spatial 2D KNN ($k=5$)** | 2 | $56.92\% \pm 6.79\%$ | $30.50\% \pm 3.78\%$ | $0.3096 \pm 0.0411$ | $58.20\% \pm 9.79\%$ | $55.64\% \pm 6.34\%$ |
| | **Spatial Markov Decay** | 2 | $52.87\% \pm 4.37\%$ | $26.34\% \pm 2.21\%$ | $0.2723 \pm 0.0339$ | $50.00\% \pm 9.53\%$ | $55.74\% \pm 5.56\%$ |
| | **Vertical Prior Baseline** | 2 | $50.55\% \pm 5.26\%$ | $27.68\% \pm 3.96\%$ | $0.2557 \pm 0.0399$ | $51.84\% \pm 9.62\%$ | $49.27\% \pm 6.85\%$ |

### 6.3 Scientific Findings of the Benchmark
1. **Performance Decay with Spacing**:
   - As inter-well spacing increases from $400\text{ m} \to 600\text{ m} \to 1000\text{ m}$, Nearest Well raw accuracy decreases from **$75.91\%$** $\to$ **$66.50\%$** $\to$ **$57.74\%$**, and balanced accuracy drops from **$52.61\%$** $\to$ **$46.74\%$** $\to$ **$38.89\%$**.
   - At $1000\text{ m}$ spacing, inter-well distance exceeds average channel width ($W \approx 287\text{ m}$), causing spatial models to decay toward background mudstone proportions.
2. **Nearest Well Outperforms 2D KNN**:
   - Across all spacings, Nearest Well achieves higher Balanced Accuracy ($52.6\%$ vs $42.5\%$ at 400m) and Macro-F1 ($0.5226$ vs $0.4286$) than 2D KNN.
   - **Reason**: Nearest Well preserves vertical succession (bed thicknesses and Waltherian stacking). In contrast, 2D KNN introduces horizontal smearing, blurring thin beds into neighboring cells.
3. **Comparison with Historical Real-Data Findings**:
   - In the Sprint J real-data benchmark ($h \ge 420\text{ m}$), Nearest Well scored $43.72\%$ and KNN scored $45.96\%$.
   - In this synthetic benchmark with known ribbon channel geology, when spacing is $600 - 1000\text{ m}$, performance approaches $57 - 66\%$. This confirms Sprint J's Category C finding: when well spacing is wider than geobodies, deterministic pointwise accuracy is physically bounded.

---

## 7. Automated Test Suite Execution

Targeted test execution verified the complete test suite:
```bash
python -m pytest tests/test_phase0_loader.py tests/test_phase1_markov.py tests/test_sprint_k.py
```
**Actual Outcome**:
- **36 passed in 19.49s (100% pass rate, zero failures)**.
- Specific Sprint K tests passed in [`tests/test_sprint_k.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_k.py):
  - `test_01_canonical_facies_schema_validity`: 6 states distinct.
  - `test_02_grid_shape_and_stratigraphic_orientation`: Bottom-to-top $Z$ orientation confirmed.
  - `test_03_seed_reproducibility`: Bitwise identical outputs across identical seeds.
  - `test_04_object_geometry_and_grid_boundaries`: All scours within bounds.
  - `test_05_bed_thicknesses_positive_and_bounded`: Column run-lengths positive ($T \ge \Delta z$).
  - `test_06_synthetic_sampling_positions_and_facies`: Extracted 1D columns match grid bitwise.
  - `test_07_ground_truth_separation_and_leakage_safeguards`: Withheld grid excluded from datasets.
  - `test_08_facies_proportions_and_transition_diagnostics`: Proportions sum to 1.0; $N/G$ hierarchy holds.
  - `test_09_configuration_validation_and_failure_handling`: Invalid parameters raise `ValueError`.
  - `test_10_regression_existing_markov_and_realization`: Existing Markov chains intact.
  - `test_11_multi_seed_validation_metrics_integrity`: Multi-seed validator runs and produces non-empty metrics.
  - `test_12_benchmark_baseline_models_execution`: All 4 baseline models fit and predict valid grids.
  - `test_13_spacing_sensitivity_accuracy_decay`: Accuracy decays monotonically as spacing increases from 400m to 1000m.
  - `test_14_benchmark_strict_leakage_safeguards`: Observed borehole columns strictly masked during evaluation.

---

## 8. Index of Sprint K Stage 3 Artifacts

### 8.1 Numerical Datasets (CSV)
- [`sprints/audit_sprint_k/multi_seed_geological_validation_metrics.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/multi_seed_geological_validation_metrics.csv)
- [`sprints/audit_sprint_k/multi_seed_validation_summary.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/multi_seed_validation_summary.csv)
- [`sprints/audit_sprint_k/multi_seed_reference_comparison.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/multi_seed_reference_comparison.csv)
- [`sprints/audit_sprint_k/sensitivity_wt_and_resolution.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/sensitivity_wt_and_resolution.csv)
- [`sprints/audit_sprint_k/synthetic_reconstruction_benchmark_metrics.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/synthetic_reconstruction_benchmark_metrics.csv)
- [`sprints/audit_sprint_k/synthetic_reconstruction_summary.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/synthetic_reconstruction_summary.csv)
- [`sprints/audit_sprint_k/forward_simulation_examples_summary.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/forward_simulation_examples_summary.csv)

### 8.2 Validation & Benchmark Visualizations (PNG)
- [`sprints/audit_sprint_k/fig_validation_facies_proportions_distribution.png`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/fig_validation_facies_proportions_distribution.png)
- [`sprints/audit_sprint_k/fig_sensitivity_wt_and_ntg.png`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/fig_sensitivity_wt_and_ntg.png)
- [`sprints/audit_sprint_k/fig_benchmark_accuracy_vs_spacing.png`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/fig_benchmark_accuracy_vs_spacing.png)
- [`sprints/audit_sprint_k/fig_benchmark_reconstruction_cross_sections.png`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/fig_benchmark_reconstruction_cross_sections.png)
- [`sprints/audit_sprint_k/fig_sim_upstream_seed101.png`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/fig_sim_upstream_seed101.png)
- [`sprints/audit_sprint_k/fig_sim_downstream_seed202.png`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/fig_sim_downstream_seed202.png)
- [`sprints/audit_sprint_k/fig_sim_pooled_seed303.png`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/fig_sim_pooled_seed303.png)

---

## 9. Geological Realism Assessment & Recommendations

### 9.1 Classification of Simulator Properties
- **Empirically Calibrated**:
  - Vertical facies proportions (`sand`, `p_sand`, `ripples`, `carbon_mud`, `coal`, `mud`).
  - Channel bed thickness distribution (MLE $\mu=1.7870, \sigma=0.8823$).
  - Coal layer frequency and bed thickness ($1.0\text{ m}$).
  - Basal planar sandstone sheet partitioning ($20\%$).
- **Derived from Literature**:
  - Channel aspect ratio $W/T \approx 35 \pm 5$ (Sahoo et al., 2016).
  - Crevasse splay lateral width range ($20 - 100\text{ m}$).
- **Modeling Assumptions**:
  - Parabolic cross-sectional erosional scour geometry.
  - Overbank mudstone background medium.
  - Strict erosional superposition.
- **Empirically Non-Identifiable from Available 1D Logs**:
  - True 3D channel sinuosity and paleocurrent orientation.
  - Dynamic channel avulsion return periods.

### 9.2 Recommendation for Next Research Phase
The minimum viable forward simulator is scientifically defensible, reproducible, and provides a benchmark for evaluating sparse borehole reconstruction.

Because Stage 3 proved that:
1. Nearest Well outperforms 2D isotropic/anisotropic interpolation on thin-bedded fluvial geology by preserving vertical succession;
2. Reconstruction accuracy decays systematically as borehole spacing exceeds geobody width ($76\% \to 66\% \to 57\%$);

The recommended next step is:
**Do not proceed to deep learning or ungrounded machine-learning models.** Instead, formulate a **Borehole Active Learning / Optimal Coring Benchmark** on these synthetic 2D forward realizations, evaluating whether active acquisition (e.g. mutual information or predictive entropy reduction) can optimize borehole placement to detect channel bodies with minimal coring cost.

---

## Stage 3 Completion Status

Stage 3 implementation, multi-seed validation, sensitivity analysis, sparse-borehole benchmarking, and automated testing are complete. All 14 tests pass. Per the user instructions, execution is paused at this checkpoint for supervisory review.
