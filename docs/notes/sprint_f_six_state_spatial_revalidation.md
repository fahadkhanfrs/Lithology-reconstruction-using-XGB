# SMALT Sprint F Technical Note: Six-State Migration and Spatial Baseline Revalidation

**Project:** SMALT (Subsurface Stratigraphic Modeling & Active Learning Toolkit)  
**Repository:** `fahadkhanfrs/Lithology-reconstruction-using-XGB`  
**Branch:** `sprint-f`  
**Date:** October 2026  
**Auditor / Senior Geostatistical ML Engineer:** Antigravity Pair  

---

## 1. Executive Summary

Sprint F completes the migration from the legacy five-state schema to the full **six-state stratigraphic facies schema** matching the field classification of Sahoo et al. (2016). Planar-laminated sandstone (`p_sand`, Facies 2) and heterolithic rippled sandstone and siltstone (`ripples`, Facies 3) are now preserved as distinct states and are no longer collapsed into channel sandstone or generic siltstone.

All affected pipeline components across data loading, descriptive analysis, 1D Markov chains, and spatial Leave-One-Litholog-Out (LOLO) cross-validation have been re-executed from the revised raw CSVs.

Key findings:
1. **Six-State Facies Architecture:** The 12 lithologs represent 1031.0 m of valid continuous observations:
   - Channel Sandstone (`sand`, Facies 1): 443.0 m (42.97%)
   - Planar Sandstone (`p_sand`, Facies 2): 63.7 m (6.18%)
   - Rippled Heterolithics (`ripples`, Facies 3): 78.7 m (7.63%)
   - Carbonaceous Mudstone (`carbon_mud`, Facies 4): 19.3 m (1.87%)
   - Coal Seams (`coal`, Facies 5): 28.5 m (2.76%)
   - Overbank Mudstone (`mud`, Facies 6): 397.8 m (38.58%)
   - Total Net Sand (`sand` + `p_sand`): 506.7 m (49.15%)
   - Total Coarse Sediments (`sand` + `p_sand` + `ripples`): 585.4 m (56.78%)
2. **Spatial Evaluation Reality Check:** Spatial 3D KNN (k=5) achieves a pooled Macro-F1 of **0.2232** (raw accuracy 47.77%, balanced accuracy 22.81%). It does not outperform the naive Nearest-Well vertical profile baseline (pooled Macro-F1 **0.2276**, accuracy 44.26%, balanced accuracy 22.91%).
3. **Severe Minority Class Failure:** Both Spatial 3D KNN and Nearest-Well baselines completely fail to predict thin chronostratigraphic marker beds (`coal`, support = 25, F1 = 0.0000; `carbon_mud`, support = 17, F1 = 0.0000). On planar sandstones (`p_sand`, support = 62), the nearest-well baseline actually outperforms Spatial KNN (F1 0.1111 vs 0.0465).
4. **Scientific Invariant Upheld:** Raw accuracy (47.77%) is inflated by predicting majority classes (Channel Sandstone and Overbank Mudstone). The training prior majority baseline alone yields 43.51% raw accuracy. Spatial 3D KNN adds only +4.26 percentage points of accuracy over a trivial single-state predictor.
5. **Coordinate & Datum Boundaries:** Provisional coordinates from [`lolo/Location_coordinates_lithologs.xlsx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/Location_coordinates_lithologs.xlsx) allow inter-well distances to be calculated, but CRS, units, and geographic origin remain unverified. Litholog 1 remains strictly excluded from spatial modeling due to absent coordinates. The vertical datum remains unanchored.

---

## 2. Six-State Facies Schema & Geological Definitions

The schema establishes a 1:1 mapping with the six sedimentary facies identified by Sahoo et al. (2016) in the Cretaceous Ferron Sandstone:

| Code | Facies # | Technical ID | Canonical Name | Color Hex | Sahoo et al. (2016) Geological Description & Architecture |
| :---: | :---: | :---: | :---: | :---: | :--- |
| **0** | Facies 1 | `sand` | Channel Sandstone | `#F4D03F` | Multi-storey channel-belt sandstone. Trough cross-bedded medium to fine sandstone. Represents high-energy fluvial channel-axis deposits. |
| **1** | Facies 2 | `p_sand` | Planar Sandstone | `#EB984E` | Planar and parallel laminated sandstone sheets. Upper flow-regime plane bed deposition or channel bar-top sheets. |
| **2** | Facies 3 | `ripples` | Rippled Heterolithics | `#5DADE2` | Interbedded rippled sandstone and siltstone. Current and wave ripple laminations, representing channel margins, levee flanks, and subaqueous splay transitions. |
| **3** | Facies 4 | `carbon_mud` | Carbonaceous Mudstone | `#6C3483` | Organic-rich muddy sediment deposited in poorly drained backswamps, oxbow lakes, and mire margins. |
| **4** | Facies 5 | `coal` | Coal | `#1C2833` | Autochthonous mire peat / coal seam. Represents prolonged organic accumulation during low clastic influx. Key chronostratigraphic markers. |
| **5** | Facies 6 | `mud` | Overbank Mudstone | `#95A5A6` | Structureless to faintly laminated fine siliciclastic mudstone. Floodplain and interdistributary bay deposits. |

### Normalization and Rejection Policy
- Label normalization is strict: case-insensitive, whitespace-trimmed, and hyphen-to-underscore normalized.
- Aliases map legacy variants safely:
  - `silt` and `siltstone` map to `ripples` (Facies 3).
  - `planar_sand` and `parallel_sand` map to `p_sand` (Facies 2).
  - `channel_sand` maps to `sand` (Facies 1).
- Unrecognized strings (e.g. `limestone`, `granite`, null) raise `ValueError` immediately. No silent fallback or unmapped recoding occurs.

---

## 3. Regenerated Files and Deliverables

All deliverables in [`sprints/audit_sprint_f/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/) have been regenerated by [`scripts/run_sprint_f_pipeline.py`](file:///d:/Lithology-reconstruction-using-XGB/scripts/run_sprint_f_pipeline.py):

| Deliverable Artifact | Description & Scope |
| :--- | :--- |
| [`data/processed/lithologs_unified.parquet`](file:///d:/Lithology-reconstruction-using-XGB/data/processed/lithologs_unified.parquet) | Unified 1031-row tabular dataset with 6-state integer codes (0..5). |
| [`data/processed/lithologs_unified.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/processed/lithologs_unified.csv) | Plain-text export of unified 6-state stratigraphy. |
| [`sprints/audit_sprint_f/facies_schema_six_state.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/facies_schema_six_state.csv) | Full 6-state technical reference table. |
| [`sprints/audit_sprint_f/per_litholog_six_state_statistics.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/per_litholog_six_state_statistics.csv) | Continuous thickness, proportions, bed statistics, and 1m discretized deltas for all 12 logs. |
| [`sprints/audit_sprint_f/group_six_state_statistics.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/group_six_state_statistics.csv) | Comparative table for All Logs, Upstream (8 logs), Downstream Outcrop (2 logs), and Downstream Composite (3 logs). |
| [`sprints/audit_sprint_f/markov_transition_matrices_six_state.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/markov_transition_matrices_six_state.csv) | 6x6 regular and embedded count, probability, and stationary vectors across all scopes. |
| [`sprints/audit_sprint_f/markov_lolo_six_state_results.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/markov_lolo_six_state_results.csv) | Leakage-free 1D Markov cross-validation log-likelihoods and perplexities for all 12 logs. |
| [`sprints/audit_sprint_f/provisional_coordinates_audit.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/provisional_coordinates_audit.csv) | Complete coordinate provenance audit table (L1 missing, L2-L12 local coordinates). |
| [`sprints/audit_sprint_f/spatial_lolo_six_state_results.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/spatial_lolo_six_state_results.csv) | Per-fold spatial evaluation table comparing KNN, Nearest-Well, and Prior across L2-L12. |
| [`sprints/audit_sprint_f/spatial_per_class_metrics.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/spatial_per_class_metrics.csv) | Per-class precision, recall, F1-score, and support for all three models across all 6 classes. |
| [`sprints/audit_sprint_f/spatial_confusion_matrices.json`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/spatial_confusion_matrices.json) | Full 6x6 confusion matrices for KNN, Nearest-Well, and Training Prior. |
| [`sprints/audit_sprint_f/directional_spatial_six_state_results.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/directional_spatial_six_state_results.csv) | Cross-group spatial prediction results (Upstream <-> Downstream). |
| [`sprints/audit_sprint_f/datum_sensitivity_six_state.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/datum_sensitivity_six_state.csv) | Stratigraphic datum offset sensitivity results (-20m to +20m) for relative-to-base and relative-to-top. |

### Publication Figures Regenerated (6-State Schema)
1. [`sprints/audit_sprint_f/figures/transect_vertical_successions_six_state.png`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/figures/transect_vertical_successions_six_state.png): Side-by-side stratigraphic columns L1-L12 with 6-color palette, provenance badges, and unanchored datum caveat.
2. [`sprints/audit_sprint_f/figures/facies_proportions_and_ntg_six_state.png`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/figures/facies_proportions_and_ntg_six_state.png): Stacked continuous facies proportions and annotated pure Net Sand across all 12 logs.
3. [`sprints/audit_sprint_f/figures/upstream_vs_downstream_facies_ntg_six_state.png`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/figures/upstream_vs_downstream_facies_ntg_six_state.png): 6-bar proximal vs distal comparative distribution and Net-to-Gross comparison.
4. [`sprints/audit_sprint_f/figures/sandstone_thickness_distributions_six_state.png`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/figures/sandstone_thickness_distributions_six_state.png): Fluvial sandstone lithosome thickness distributions distinguishing Channel and Planar beds.
5. [`sprints/audit_sprint_f/figures/markov_transition_heatmaps_six_state.png`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_f/figures/markov_transition_heatmaps_six_state.png): 6x6 heatmaps for Regular (1m grid) and Embedded (boundary-crossing beds) transition probability matrices.

---

## 4. Coordinate Provenance and Remaining Limitations

### Coordinate Traceability
- Coordinates originate from [`lolo/Location_coordinates_lithologs.xlsx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/Location_coordinates_lithologs.xlsx).
- Columns ingested: Well identifier (`L2`..`L12`), Local Cartesian Easting $X$ (meters), Local Cartesian Northing $Y$ (meters).
- **Litholog 1 is missing** from the spreadsheet and has no coordinates. It is strictly excluded from all spatial pipelines.
- Legacy synthetic coordinates (`strike_pos_m = file_idx * 100`) have been audited and verified to never enter spatial model features. Features are strictly $[X, Y, Z_{rel}]$.

### Spatial Inter-Well Distances
- Upstream wells (L2-L8, L10) form a localized cluster spanning ~2,100 m East-West and ~1,100 m North-South.
  - Closest well pair: L7 to L8 (155.6 m).
  - L2 to L3: 400.3 m.
  - L4 to L5: 350.0 m.
- Downstream outcrop wells (L9 and L11) are located ~7.5 km to the Northwest.
  - Distance between L9 and L11: 457.7 m.
- Downstream drill-core EM-137C (L12) is located ~14.4 km to the North-Northwest.
  - Distance between L12 and L9: 6,861.9 m.
  - Distance between L12 and nearest upstream well (L10): 14,037.6 m.

### Unverified Geodetic Status
- **CRS is unverified:** The projection (e.g. UTM Zone 12N NAD83 or local mine grid) is undocumented in the spreadsheet.
- **Units are assumed meters:** Values range from -1,120 to +1,600 in $X$ and 0 to 14,407 in $Y$.
- **Geographic Origin is unverified:** The coordinate $(0, 0)$ is an arbitrary project-local datum.
- **Geological Consequence:** Because inter-well lateral distances (150 m to 14 km) exceed typical channel-belt widths (often 50-300 m) by 1 to 2 orders of magnitude, spatial interpolation without correlated chronostratigraphic markers is severely underconstrained.

---

## 5. Stratigraphic Datum Assumptions & Sensitivity Results

### Current Vertical Assumption
- Modern cliff tops and drill-core collars are **not** an isochronous stratigraphic datum.
- Equal numerical depths across different logs do **not** represent geological time-equivalence.
- In the absence of a confirmed chronostratigraphic marker (such as a correlated coal seam or flooding surface), the pipeline adopts a provisional **relative-to-base** reference:
  $$Z_{rel} = \max(Depth) - Depth + \Delta_{datum}$$
- An anisotropic vertical weight factor of 10.0 is applied ($Z_{scaled} = 10 \cdot Z_{rel}$) to reflect that stratigraphic correlation scales are 10-100x finer vertically than horizontally.

### Datum Offset Sensitivity Analysis (6-State Target)
We evaluated sensitivity to datum shifts by introducing systematic vertical offsets ($\Delta_{datum} \in \{-20, -10, -5, 0, +5, +10, +20\}$ meters) between upstream and downstream logs:

| Vertical Reference | Datum Offset ($\Delta$) | KNN Macro-F1 | KNN Balanced Acc | KNN Accuracy | Nearest-Well Macro-F1 | Nearest-Well Accuracy | Prior Macro-F1 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `relative_to_base` | -20 m | 0.2219 | 22.72% | 48.09% | 0.2276 | 44.26% | 0.1011 |
| `relative_to_base` | -10 m | 0.2172 | 22.22% | 47.23% | 0.2276 | 44.26% | 0.1011 |
| `relative_to_base` | -5 m | 0.2197 | 22.46% | 47.34% | 0.2276 | 44.26% | 0.1011 |
| **`relative_to_base`** | **0 m (Default)** | **0.2232** | **22.81%** | **47.77%** | **0.2276** | **44.26%** | **0.1011** |
| `relative_to_base` | +5 m | 0.2198 | 22.49% | 47.77% | 0.2276 | 44.26% | 0.1011 |
| `relative_to_base` | +10 m | 0.2270 | 23.19% | 48.19% | 0.2276 | 44.26% | 0.1011 |
| `relative_to_base` | +20 m | 0.2164 | 22.14% | 46.81% | 0.2276 | 44.26% | 0.1011 |
| `relative_to_top` | -20 m | 0.2278 | 23.20% | 49.04% | 0.2303 | 45.85% | 0.1011 |
| `relative_to_top` | 0 m | 0.2256 | 22.98% | 48.51% | 0.2303 | 45.85% | 0.1011 |
| `relative_to_top` | +20 m | 0.2274 | 23.15% | 48.83% | 0.2303 | 45.85% | 0.1011 |

**Datum Sensitivity Conclusion:** Across a 40-meter vertical offset perturbation, Spatial 3D KNN Macro-F1 remains locked in a narrow band between 0.2164 and 0.2270. This demonstrates that performance is not constrained by fine datum calibration, but by lateral facies heterogeneity and wide inter-well spacing.

---

## 6. Spatial LOLO Cross-Validation Results

### Per-Litholog Performance (11 Held-Out Wells, 940 Total Grid Points)

| Held-Out Well | Nearest Well | Distance (m) | Points | KNN Acc | KNN BalAcc | KNN Macro-F1 | Near-Well Acc | Near-Well Macro-F1 | Prior Macro-F1 | Net Sand True | KNN Sand Err |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **L2** | L3 | 400.3 | 93 | 69.89% | 36.02% | 0.3642 | 69.89% | 0.3675 | 0.1837 | 0.6989 | 0.0000 |
| **L3** | L2 | 400.3 | 93 | 74.19% | 41.71% | 0.3703 | 69.89% | 0.3675 | 0.2112 | 0.7312 | 0.0430 |
| **L4** | L5 | 350.0 | 84 | 39.29% | 15.32% | 0.1442 | 36.90% | 0.1727 | 0.0899 | 0.3929 | 0.0476 |
| **L5** | L4 | 350.0 | 84 | 38.10% | 14.11% | 0.1375 | 36.90% | 0.1727 | 0.0940 | 0.3929 | 0.0357 |
| **L6** | L7 | 180.3 | 84 | 39.29% | 17.63% | 0.1770 | 39.29% | 0.1786 | 0.1128 | 0.3929 | 0.0238 |
| **L7** | L8 | 155.6 | 79 | 44.30% | 24.82% | 0.1989 | 39.24% | 0.1898 | 0.1274 | 0.4430 | 0.0253 |
| **L8** | L7 | 155.6 | 79 | 43.04% | 21.79% | 0.2025 | 39.24% | 0.1898 | 0.0871 | 0.4684 | 0.0380 |
| **L9** | L11 | 457.7 | 78 | 44.87% | 32.91% | 0.2762 | 41.03% | 0.2275 | 0.1084 | 0.4487 | 0.0769 |
| **L10** | L8 | 806.2 | 77 | 33.77% | 16.83% | 0.1486 | 38.96% | 0.2042 | 0.0714 | 0.5065 | 0.0779 |
| **L11** | L9 | 457.7 | 78 | 46.15% | 28.86% | 0.2813 | 41.03% | 0.2275 | 0.1032 | 0.5128 | 0.0513 |
| **L12** | L9 | 6,861.9 | 111 | 45.95% | 23.81% | 0.2326 | 31.53% | 0.1706 | 0.1341 | 0.5045 | 0.0991 |

### Overall Aggregate Summary Metrics

| Evaluation Metric | Spatial 3D KNN (k=5) | Nearest-Well Vertical Profile | Training Prior (Majority Class) |
| :--- | :---: | :---: | :---: |
| **Pooled Macro-F1** | **0.2232** | **0.2276** | 0.1011 |
| **Unweighted Mean Macro-F1** | 0.2303 | 0.2244 | 0.1203 |
| **Balanced Accuracy** | 22.81% | 22.91% | 16.67% |
| **Raw Accuracy** | 47.77% (449/940) | 44.26% (416/940) | 43.51% (409/940) |
| **Mean Net Sand Error** | 4.71% | 5.38% | 8.35% |

---

## 7. Per-Class Precision, Recall, and F1 Breakdown

The table below shows why Macro-F1 is low across all models:

| Facies State | Support | Spatial 3D KNN Prec / Rec / F1 | Nearest-Well Prec / Rec / F1 | Training Prior Prec / Rec / F1 |
| :--- | :---: | :---: | :---: | :---: |
| **Channel Sandstone (`sand`)** | 409 | 0.5760 / 0.6112 / **0.5931** | 0.5756 / 0.5306 / **0.5522** | 0.4351 / 1.0000 / **0.6064** |
| **Planar Sandstone (`p_sand`)** | 62 | 0.0833 / 0.0323 / **0.0465** | 0.1094 / 0.1129 / **0.1111** | 0.0000 / 0.0000 / **0.0000** |
| **Rippled Heterolithics (`ripples`)** | 77 | 0.2254 / 0.2078 / **0.2162** | 0.2250 / 0.2338 / **0.2293** | 0.0000 / 0.0000 / **0.0000** |
| **Carbonaceous Mudstone (`carbon_mud`)** | 17 | 0.0000 / 0.0000 / **0.0000** | 0.0000 / 0.0000 / **0.0000** | 0.0000 / 0.0000 / **0.0000** |
| **Coal (`coal`)** | 25 | 0.0000 / 0.0000 / **0.0000** | 0.0000 / 0.0000 / **0.0000** | 0.0000 / 0.0000 / **0.0000** |
| **Overbank Mudstone (`mud`)** | 350 | 0.4536 / 0.5171 / **0.4833** | 0.4508 / 0.4971 / **0.4728** | 0.0000 / 0.0000 / **0.0000** |

### Confusion Matrix (Spatial 3D KNN, Pooled 940 Points)

```
                       Predicted Facies
                  Sand  p_Sand Ripples CarbMud  Coal   Mud   Total
True Sand          250       2      20       6     2   129     409
True p_Sand         11       2       3       1     0    45      62
True Ripples        35       1      16       0     0    25      77
True CarbMud         4       2       2       0     0     9      17
True Coal           14       0       1       0     0    10      25
True Mud           120      17      29       3     0   181     350
Total Predicted    434      24      71      10     2   399     940
```

### Critical Geostatistical Diagnostics
1. **Total Collapse on Mires and Lakes:** Zero coal points (0/25) and zero carbonaceous mudstone points (0/17) were correctly predicted by either spatial model. Coal is falsely predicted as channel sandstone (14 points) or overbank mudstone (10 points).
2. **Nearest-Well Superiority on Bars:** On planar sandstones (`p_sand`), the Nearest-Well vertical profile baseline achieves F1 = 0.1111 (7 correct), outperforming Spatial 3D KNN which achieves F1 = 0.0465 (only 2 correct).
3. **Net Sand Estimation:** Although pointwise multi-class facies prediction fails on thin beds, overall Net Sand prediction error is relatively small (~4.7% for KNN, ~5.4% for Nearest-Well). This confirms that aggregate reservoir volume / sand proportion can be estimated far more reliably than exact 3D spatial facies geometry.

---

## 8. Directional Cross-Group Experiments

We evaluated cross-group transferability between proximal (Upstream, 8 logs, 673 m) and distal (Downstream, 3 logs, 267 m):

| Experiment | Training Scope | Target Scope | Test Points | KNN Macro-F1 | KNN BalAcc | KNN Accuracy | Prior Accuracy |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Upstream -> Downstream** | L2-L8, L10 (8 wells) | L9, L11, L12 (3 wells) | 267 | 0.2050 | 24.19% | 39.33% | 44.94% |
| **Downstream -> Upstream** | L9, L11, L12 (3 wells) | L2-L8, L10 (8 wells) | 673 | 0.1999 | 20.67% | 38.93% | 42.94% |

**Directional Conclusion:** In both directions, Spatial 3D KNN achieves raw accuracy (38.9%-39.3%) that is **lower** than the training prior baseline (42.9%-44.9%). The spatial model fails to transfer across the 7.5-14 km inter-group gap.

---

## 9. Verification & Test Suite Execution

All 65 automated tests across the test suite pass cleanly:
- [`tests/test_phase0_loader.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_phase0_loader.py): 11 tests passed.
- [`tests/test_phase1_markov.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_phase1_markov.py): 11 tests passed.
- [`tests/test_sprint_c.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_c.py): 12 tests passed.
- [`tests/test_sprint_d.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_d.py): 10 tests passed.
- [`tests/test_sprint_e.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_e.py): 11 tests passed.
- [`tests/test_sprint_f.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_f.py): 10 tests passed.

---

## 10. Known Limitations & Next Data Required from Prof. Hiranya Sahoo

1. **Stratigraphic Datum:** The current model assumes modern well bases are aligned. This is a severe geological simplification. To anchor stratigraphy, we need Prof. Sahoo to confirm a correlated chronostratigraphic datum (e.g. Subsequence boundary, maximum flooding surface, or correlated coal seam).
2. **Litholog 1 Coordinates:** Litholog 1 cannot be modeled spatially until its coordinates or field outcrop location are supplied.
3. **Geodetic CRS & Geographic Origin:** The local Cartesian coordinates in `Location_coordinates_lithologs.xlsx` need confirmation of their datum/projection (e.g. UTM Zone 12N) and whether units are exact meters.
4. **Subsurface Depth of Core EM-137C (L12):** The 0-111 m digitized section represents the uppermost portion of a 242 m borehole. Absolute depth below ground level and top elevation are needed for true 3D spatial positioning.
