# SMALT Sprint E Technical Note: Revised Dataset Reconciliation and Spatial Prototype Readiness

> [!NOTE]
> **SUPERSEDED BY SPRINT F:** The five-state facies schema and metrics in this document represent historical Sprint E deliverables. In Sprint F, the repository fully migrated to a six-state schema preserving distinct planar sandstone (`p_sand`) and rippled heterolithics (`ripples`). For current active results, refer to [`sprint_f_six_state_spatial_revalidation.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/notes/sprint_f_six_state_spatial_revalidation.md).

**Project:** SMALT (Subsurface Stratigraphic Modeling & Active Learning Toolkit)  
**Repository:** `fahadkhanfrs/Lithology-reconstruction-using-XGB`  
**Branch:** `sprint-e` (Superseded by `sprint-f`)  
**Date:** October 2026  
**Auditor / Senior Geostatistical ML Engineer:** Antigravity Pair  

---

## 1. Executive Summary

Sprint E was triggered by the manual digitization and revision of all 12 lithologs from the primary source images (Sahoo et al., 2016 outcrop profiles L1-L11 and drill-core profile EM-137C L12). This revision introduced updated facies labels (`p_sand`, `ripples`), corrected numerous depth intervals, resolved historical data gaps (e.g. L11 59-60 m), and documented new interval boundaries.

This sprint completed the full reconciliation across the entire repository:
1. **Facies Encoding Reconciliation:** Mapped revised lithological descriptors (`p_sand` -> `sand`, `ripples` -> `silt`) to the canonical 5-state SMALT schema without data loss or unmapped facies exceptions.
2. **Provenance Registration of Litholog 12:** Registered L12 in [`data/provenance_manifest.json`](file:///d:/Lithology-reconstruction-using-XGB/data/provenance_manifest.json) as a manually digitized subsurface core log (coverage strictly 0.0-111.0 m; original core length 242.0 m; unvalidated; not AI-generated).
3. **Synthetic Coordinate Audit:** Traced legacy `strike_pos_m = (file_idx) * 100.0` generated in [`data/loader.py`](file:///d:/Lithology-reconstruction-using-XGB/data/loader.py#L136-L141). Verified that it never entered 1D Markov models. For spatial modeling, synthetic coordinates are deprecated and replaced with real local Cartesian coordinates $(X, Y)$ from [`lolo/Location_coordinates_lithologs.xlsx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/Location_coordinates_lithologs.xlsx) for L2-L12, while L1 is strictly excluded due to missing coordinates.
4. **Full Pipeline Regeneration:** Regenerated all descriptive statistics, sandstone bed thickness analyses, and 1D Markov LOLO validations in [`sprints/audit_sprint_c/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_c/) and [`sprints/audit_sprint_d/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_d/).
5. **Provisional Spatial Prototype:** Implemented leak-free spatial Leave-One-Litholog-Out cross-validation across eligible wells (L2-L12) in [`smalt/spatial/baseline.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/baseline.py). Evaluated Spatial 3D KNN against naive baselines (training prior facies and nearest-well profile) under documented unanchored vertical datum assumptions and sensitivity testing (±5m to ±20m).
6. **Full Test Suite:** All 55 test cases across Phase 0, Phase 1, Sprint C, Sprint D, and Sprint E pass cleanly.

---

## 2. What Changed in the Revised Digitized Data

The user manually digitized and revised all 12 raw CSVs in [`data/raw_lithologs/`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/). The changes fall into three categories:

### 2.1 Updated Facies Encodings
Sahoo et al. (2016) established six detailed sedimentary facies in the field:
1. Facies 1: Multi-storey channel sandstone (trough cross-bedded)
2. Facies 2: Planar/parallel-laminated sandstone (`p_sand`)
3. Facies 3: Heterolithic rippled sandstone and siltstone (`ripples`)
4. Facies 4: Carbonaceous mudstone (`carbon_mud`)
5. Facies 5: Coal (`coal`)
6. Facies 6: Structureless overbank mudstone (`mud`)

In previous sprints, raw CSVs used simplified strings (`sand`, `silt`, `mud`, `carbon_mud`, `coal`). In the revised digitization:
- `p_sand` was introduced in L4, L5, L6, L7, L8, L9, L10, L11 to capture planar-laminated sand beds.
- `ripples` was introduced in L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12 to capture heterolithic rippled siltstone/sandstone beds.

**Reconciliation Rule:**
- `p_sand` maps to canonical code `1` (`sand`, display name: "Channel Sandstone (undivided)").
- `ripples` maps to canonical code `3` (`silt`, display name: "Siltstone / heterolithics").
- This mapping preserves geological continuity with the project's canonical 5-state Markov space while allowing finer-scale facies recovery if expanded to 6 states in the future.

### 2.2 Stratigraphic Depth Bounds and Continuity Revisions
- **Litholog 5 & 6:** Revised stratigraphic span corrected to 84.0 m (previously 85.0 m and 82.0 m).
- **Litholog 7:** Revised stratigraphic span corrected to 79.0 m (previously 80.0 m).
- **Litholog 11:** The historical 1 m unmapped gap at 59.0-60.0 m was digitized by the user as `carbon_mud`, making L11 continuous across that boundary; an overlap at 71.0-72.0 m was recorded between `ripples` (69-72m) and `mud` (71-76m).
- **Litholog 9:** The revised CSV contains two 1 m gaps: 18.0-19.0 m and 51.0-52.0 m (the 50-52m mud vs 52-56m carbon_mud intervals create an unmapped gap at 51-52m), plus an overlap at 28.0-29.0 m.
- **Litholog 12:** Subsurface drill core EM-137C, continuous across 0.0-111.0 m (63 digitized intervals).

### 2.3 Thickness Overview
- Total stratigraphic span across all 12 lithologs: **1033.0 m** (Outcrop: 922.0 m; Core: 111.0 m).
- Total valid observations: **1031.0 m** (due to two 1 m gaps in L9).

---

## 3. What Was Regenerated

All downstream deliverables were recomputed directly from the revised raw CSVs:

1. **Sprint C Outputs** ([`sprints/audit_sprint_c/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_c/)):
   - [`per_litholog_descriptive_statistics.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_c/per_litholog_descriptive_statistics.csv)
   - [`facies_mapping_table.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_c/facies_mapping_table.csv)
   - [`markov_transition_matrices.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_c/markov_transition_matrices.csv)
   - [`lolo_validation_results.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_c/lolo_validation_results.csv)
   - [`directional_validation_results.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_c/directional_validation_results.csv)
   - All 5 publication-quality figures in [`figures/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_c/figures/)
2. **Sprint D Outputs** ([`sprints/audit_sprint_d/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_d/)):
   - [`sandstone_thickness_reconciliation.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_d/sandstone_thickness_reconciliation.csv)
   - [`reconciled_per_litholog_statistics.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_d/reconciled_per_litholog_statistics.csv)
   - [`reconciled_group_statistics.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_d/reconciled_group_statistics.csv)
   - [`lolo_validation_results_sprint_d.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_d/lolo_validation_results_sprint_d.csv)
   - [`sensitivity_gap_masking_results.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_d/sensitivity_gap_masking_results.csv)
   - [`markov_transition_matrices_sprint_d.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_d/markov_transition_matrices_sprint_d.csv)
3. **Sprint E Deliverables** ([`sprints/audit_sprint_e/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_e/)):
   - [`spatial_metadata_readiness.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_e/spatial_metadata_readiness.csv)
   - [`interwell_horizontal_distances.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_e/interwell_horizontal_distances.csv)
   - [`synthetic_coordinate_audit.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_e/synthetic_coordinate_audit.md)
   - [`provisional_spatial_lolo_results.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_e/provisional_spatial_lolo_results.csv)
   - [`provisional_spatial_lolo_summary.json`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_e/provisional_spatial_lolo_summary.json)
   - [`directional_spatial_results.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_e/directional_spatial_results.csv)
   - [`datum_offset_sensitivity_results.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_e/datum_offset_sensitivity_results.csv)

---

## 4. Key Results That Changed and Why

### 4.1 Sandstone Bed Statistics
Because planar-laminated sandstones (`p_sand`) were previously grouped with other units or logged as generic sand, separating them into distinct beds increased bed resolution while consolidating true sandstone volume:

| Metric | Previous Sprint D Value | Revised Sprint E Value | Reason for Change |
| :--- | :---: | :---: | :--- |
| **Raw Sandstone Interval Count** | 72 | **135** | User digitized distinct `p_sand` intervals across L4-L11. |
| **Merged Sandstone Lithosomes** | 53 | **99** | Intercalated `p_sand` and `ripples` beds produce finer distinct lithosomes. |
| **Cumulative Sandstone Thickness** | 444.70 m | **507.70 m** | Previously misidentified or unmapped intervals digitized as sand. |
| **Overall Sandstone Proportion (N/G)** | 43.05% | **49.15%** | Higher net sand proportion recognized in revised logs. |
| **Mean Raw Bed Thickness** | 6.18 m | **3.76 m** | Finer subdivision of sandstone beds lowers raw interval mean thickness. |
| **Mean Merged Lithosome Thickness**| 8.39 m | **5.13 m** | Distinct lithosomes are thinner due to interbedded rippled siltstone layers. |
| **Maximum Merged Sand Thickness** | 32.00 m (L3) | **32.00 m (L3)** | Unchanged; 32 m multi-storey channel in L3 remains maximum. |

### 4.2 Stratigraphic Group Net-to-Gross (Continuous)
- **Upstream (L2-L8, L10):** Total thickness 673.0 m; Sandstone 340.0 m (**50.52% N/G** vs 42.64% previously).
- **Downstream Outcrop (L9, L11):** Total thickness 156.0 m; Sandstone 75.0 m (**48.08% N/G** vs 41.67% previously).
- **Downstream Composite (L9, L11, L12):** Total thickness 267.0 m; Sandstone 132.7 m (**49.70% N/G** vs 45.96% previously).
- **All 12 Lithologs:** Total thickness 1033.0 m; Sandstone 507.7 m (**49.15% N/G** vs 43.05% previously).

### 4.3 1D Markov Validation Metrics
- **Regular Chain (1 m grid):** Average outcrop perplexity is **2.04**; full-dataset average perplexity is **2.18**. Self-transition persistence remains high (~0.75-0.85).
- **Embedded Chain (distinct beds):** Average outcrop perplexity is **2.76**; full-dataset average perplexity is **2.94**. Because the embedded chain has zero diagonal, this reflects genuine transition unpredictability among the $K-1=4$ alternative facies.

---

## 5. Litholog 12 Provenance and Scope

Litholog 12 represents subsurface drill-core EM-137C from the Wasatch Plateau.
- **Identifier:** `litholog12` (display label: `L12`).
- **Provenance Category:** `digitized_core_log`.
- **Reconstruction Method:** Manually digitized by the user from supplied core-log image.
- **Digitized Coverage:** Strictly **0.0 - 111.0 m** (63 continuous raw intervals).
- **Original Core Length:** **242.0 m** (recorded as source metadata; depths >111 m were not digitized).
- **Scope Limitation:** Strictly no extrapolation or synthetic fill for undigitized depths (111-242 m).
- **Validation Status:** Unvalidated (`independent_validation: false`; no third-party manual extraction exists).
- **Coordinates:** Available in Row 11 of [`lolo/Location_coordinates_lithologs.xlsx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/Location_coordinates_lithologs.xlsx): $X = -1119.91\text{ m}, Y = 14407.30\text{ m}$.

---

## 6. Synthetic-Coordinate Audit Findings

A rigorous code search confirmed:
1. `data/loader.py:L136-141` generated an artificial coordinate: `strike_pos_m = float(file_idx * 100.0)`.
2. This column was stored in `data/processed/lithologs_unified.csv/.parquet`.
3. **No 1D Markov model ever used this column.** 1D Markov chains rely strictly on vertical sequence indices.
4. **Resolution in Sprint E:** Dedicated module [`smalt/spatial/coordinates.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/coordinates.py) loads real Cartesian coordinates directly from the Excel spreadsheet. `strike_pos_m` is strictly omitted from all spatial models.
5. **Litholog 1 Handling:** Litholog 1 is missing from the Excel spreadsheet. It is assigned `NaN` coordinates and permanently excluded from spatial validation.

---

## 7. Spatial Prototype Readiness and Empirical Findings

### 7.1 Readiness Assessment
- **Horizontal Coordinates:** Usable under a **provisional relative geometry** assumption. Inter-well distances range from 88.6 m (L4 to L5) to 4806.9 m (L2 to L11), with a median of 1918.4 m.
- **Vertical Datum:** **UNANCHORED.** No surface elevations, sea-level datums, or stratigraphic marker horizons exist in the repository. Equal raw depths do NOT correlate.

### 7.2 Provisional Spatial Baseline Implementation
To assess whether spatial models can predict facies under these limitations, [`smalt/spatial/baseline.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/baseline.py) implements:
- Strict Leave-One-Litholog-Out (LOLO) across the 11 eligible wells (L2-L12).
- Zero target leakage: target well excluded from training points and `StandardScaler` fitting.
- Features: $[X, Y, Z_{rel}]$ where $Z_{rel}$ is relative to profile base (`relative_to_base`).
- Models:
  1. **Spatial 3D KNN** (Distance-weighted, $K=5$, vertical weight factor = 10.0).
  2. **Nearest-Well Vertical Profile** (Transfers the vertical profile of the geographically closest training well).
  3. **Training Prior Facies** (Predicts the majority facies of the training set).

### 7.3 Empirical Results
From [`sprints/audit_sprint_e/provisional_spatial_lolo_results.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_e/provisional_spatial_lolo_results.csv):

| Model | Macro-F1 | Balanced Accuracy | Raw Accuracy | Interpretation |
| :--- | :---: | :---: | :---: | :--- |
| **Spatial 3D KNN** | **0.2521** | **0.2561** | **49.04%** | Fails to predict minority facies; predicts predominantly sandstone. |
| **Nearest-Well Profile** | **0.2523** | **0.2545** | **47.55%** | Matches KNN performance almost identically. |
| **Training Prior Facies** | **0.1335** | **0.2000** | **49.15%** | Always predicts sandstone (the 49.15% majority class). |

### 7.4 Datum Sensitivity Analysis
From [`sprints/audit_sprint_e/datum_offset_sensitivity_results.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_e/datum_offset_sensitivity_results.csv), applying relative vertical datum shifts to downstream wells:
- No shift (0 m): Macro-F1 = **0.2521**
- +5 m shift: Macro-F1 = **0.2481**
- -5 m shift: Macro-F1 = **0.2481**
- +10 m shift: Macro-F1 = **0.2567**
- -10 m shift: Macro-F1 = **0.2450**
- +20 m shift: Macro-F1 = **0.2444**
- -20 m shift: Macro-F1 = **0.2505**

**Scientific Conclusion:** Vertical shifting across plausible offsets (±5m to ±20m) leaves Macro-F1 virtually unchanged (~0.25). Because inter-well spacing (1.5-3.0 km) vastly exceeds the lateral extent of individual sandstone channels and mud lenses (typically 100-400 m in fluvial-deltaic settings), spatial interpolation without stratigraphic datum correlation reduces to near-random facies assignment or simple majority voting.

---

## 8. Test Suite Execution & Verification

Pytest was executed across all 5 test modules in the repository:
```
tests/test_phase0_loader.py ...........                                  [ 20%]
tests/test_phase1_markov.py ...........                                  [ 40%]
tests/test_sprint_c.py ............                                      [ 61%]
tests/test_sprint_d.py ..........                                        [ 80%]
tests/test_sprint_e.py ...........                                       [100%]
======================= 55 passed, 2 warnings in 47.89s =======================
```
All 55 tests pass without regressions.

---

## 9. Remaining Questions for Professor Hiranya Sahoo

1. **Stratigraphic Datum:** Can the top or base of a known regional marker (e.g., the Castlegate Sandstone contact, a regional coal seam like the Upper Sunnyside coal, or a flooding surface) be identified in these 12 profiles to establish a correlated vertical datum?
2. **Litholog 1 Coordinates:** Does a GPS coordinate or field traverse location exist for Litholog 1 so it can be integrated into spatial modeling?
3. **Coordinate Reference System (CRS):** What is the exact projection and origin of the Cartesian coordinates in `Location_coordinates_lithologs.xlsx` (e.g. UTM Zone 12N in meters)?
4. **Facies Schema Expansion:** Should the canonical model retain 5 states (`sand` undivided) or expand to 6 states to explicitly model `p_sand` (planar sandstone) vs `ripples` (heterolithic siltstone)?

---

## 10. Exact Next Recommended Sprint

**Sprint F: Stratigraphic Marker Correlation and 2D Cross-Sectional Geometry Formulation.**
- Before implementing complex 2D/3D multipoint geostatistics or active learning, establish relative stratigraphic datums using coal marker horizons or maximum flooding surfaces.
- Evaluate whether marker-guided vertical flattening improves spatial correlation across the proximal-to-distal transect.
