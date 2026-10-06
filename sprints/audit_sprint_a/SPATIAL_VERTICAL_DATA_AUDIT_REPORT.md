# SMALT Sprint A - Spatial and Vertical Data Audit Report

**Toolkit:** Subsurface Stratigraphic Modeling & Active Learning Toolkit (SMALT)  
**Analog System:** Cretaceous Blackhawk Formation (Wasatch Plateau, Utah; Sahoo et al., 2016)  
**Audit Authority:** Geostatistical Data Engineering (Pre-Simulation Quality Gate)  
**Execution Mode:** Read-Only Audit (Zero modification of existing codebase or raw data)  
**Audit Output Directory:** [`audit/`](file:///d:/Lithology-reconstruction-using-XGB/audit/)  

---

## 1. Repository Inventory

A comprehensive inspection of the repository was performed across all directories, scripts, data files, and documentation. No file paths or column schemas were assumed in advance.

### 1.1 Litholog and Well-Log Data Files

The repository contains three categories of subsurface and outcrop data:

| File Path | Format | Record Count | Column Names | Description & Provenance |
| :--- | :---: | :---: | :--- | :--- |
| [`data/raw_lithologs/litholog1.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog1.csv) | CSV | 20 rows | `Top`, `Bottom`, `Facies` | Source-derived outcrop section (Sahoo et al., 2016). Span: 0 to 93 m. |
| [`data/raw_lithologs/litholog2.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog2.csv) | CSV | 16 rows | `Top`, `Bottom`, `Facies` | AI-reconstructed section from published figures. Span: 0 to 93 m. |
| [`data/raw_lithologs/litholog3.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog3.csv) | CSV | 15 rows | `Top`, `Bottom`, `Facies` | AI-reconstructed section from published figures. Span: 0 to 93 m. |
| [`data/raw_lithologs/litholog4.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog4.csv) | CSV | 20 rows | `Top`, `Bottom`, `Facies` | AI-reconstructed section from published figures. Span: 0 to 84 m. |
| [`data/raw_lithologs/litholog5.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog5.csv) | CSV | 20 rows | `Top`, `Bottom`, `Facies` | AI-reconstructed section from published figures. Span: 0 to 85 m. |
| [`data/raw_lithologs/litholog6.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog6.csv) | CSV | 17 rows | `Top`, `Bottom`, `Facies` | AI-reconstructed section from published figures. Span: 0 to 82 m. |
| [`data/raw_lithologs/litholog7.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog7.csv) | CSV | 16 rows | `Top`, `Bottom`, `Facies` | AI-reconstructed section from published figures. Span: 0 to 80 m. |
| [`data/raw_lithologs/litholog8.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog8.csv) | CSV | 28 rows | `Top`, `Bottom`, `Facies` | AI-reconstructed section from published figures. Span: 0 to 79 m. |
| [`data/raw_lithologs/litholog9.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog9.csv) | CSV | 20 rows | `Top`, `Bottom`, `Facies` | Source-derived section (Sahoo et al., 2016). Span: 0 to 78 m (discontinuities resolved: 18-19m taken as `c_sand`, 28-30m overlap taken as `p_sand`). |
| [`data/raw_lithologs/litholog10.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog10.csv) | CSV | 19 rows | `Top`, `Bottom`, `Facies` | AI-reconstructed section from published figures. Span: 0 to 77 m. |
| [`data/raw_lithologs/litholog11.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog11.csv) | CSV | 18 rows | `Top`, `Bottom`, `Facies` | Source-derived benchmark section. Span: 0 to 78 m (unmapped gap at 59-60m resolved as `carbon_mud`; 93.59% accuracy benchmark). |
| [`data/raw_lithologs/litholog12.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog12.csv) | CSV | 63 rows | `Top`, `Bottom`, `Facies` | Digitized subsurface drill core (EM-137C core; Sahoo et al., 2016). Span: 0 to 111 m (Blackhawk Formation interval). |
| [`data/litholog1.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/litholog1.csv) | CSV | 20 rows | `Top`, `Bottom`, `Facies` | Exact duplicate of `data/raw_lithologs/litholog1.csv`. |
| [`data/litholog9.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/litholog9.csv) | CSV | 20 rows | `Top`, `Bottom`, `Facies` | Exact duplicate of `data/raw_lithologs/litholog9.csv`. |
| [`data/litholog11.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/litholog11.csv) | CSV | 18 rows | `Top`, `Bottom`, `Facies` | Exact duplicate of `data/raw_lithologs/litholog11.csv`. |
| [`data/log.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/log.csv) | CSV | 401 rows | `Depth`, `RxoRt`, `RLL3`, `SP`, `RILD`, `ILM`, `RILM`, `DPHI`, `PE`, `NPHI`, `RHOC`, `DPOR`, `CNLS`, `GR` | Legacy continuous wireline log (195.0 to 395.0 m) from initial XGBoost repository; not an outcrop litholog. |

### 1.2 Processed and Standardized Stratigraphic Datasets

| File Path | Format | Record Count | Column Names | Purpose & Lineage |
| :--- | :---: | :---: | :--- | :--- |
| [`data/processed/lithologs_unified.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/processed/lithologs_unified.csv) | CSV | 920 rows | `litholog_id`, `strike_pos_m`, `depth_m`, `facies_code`, `facies_name`, `gamma_ray` | Standardized 1-meter discretized dataset across all 11 logs generated by Phase 0 `LithologLoader`. |
| [`data/processed/lithologs_unified.parquet`](file:///d:/Lithology-reconstruction-using-XGB/data/processed/lithologs_unified.parquet) | Parquet | 920 rows | Same as above | Binary column-oriented version of the unified 11-log dataset. |
| [`data/provenance_manifest.json`](file:///d:/Lithology-reconstruction-using-XGB/data/provenance_manifest.json) | JSON | 156 lines | Schema tracking dictionary | Provenance registry classifying source-derived logs (L1, L9, L11) vs. AI-reconstructed logs (L2-L8, L10). |

### 1.3 Professor's Coordinate Spreadsheet and Supplementary Deliverables

Found in directory [`lolo/`](file:///d:/Lithology-reconstruction-using-XGB/lolo/):

| File Path | Format | Dimensions / Extent | Content Description | Author / Creator Metadata |
| :--- | :---: | :---: | :--- | :--- |
| [`lolo/Location_coordinates_lithologs.xlsx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/Location_coordinates_lithologs.xlsx) | XLSX | 13 rows, 2 sheets | Sheet `Location Data of 12 logs` contains columns: `[None, 'X-coordionate', 'Y-coordinate', None]`. Contains coordinates for L2-L12. Sheet `Sheet3` is empty. | Creator: `hsahoo1` (Prof. Hiranya Sahoo), Created: 2013-03-16, Modified: 2026-10-03. |
| [`lolo/litholog_locations.pptx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/litholog_locations.pptx) | PPTX | 1 slide | Annotated outcrop map based on Sahoo et al. (2016) Figure 3B/C. Contains color-coded ovals (Red = Upstream, Green = Downstream), orientation labels ("Left", "Right", "Upstream", "Downstream"), and group lists. | Author: `Hiranya Sahoo`, Created: 2026-02-18, Modified: 2026-10-03. |
| [`lolo/litholog12.pdf`](file:///d:/Lithology-reconstruction-using-XGB/lolo/litholog12.pdf) | PDF | 1 page | High-resolution sedimentological graphic core log titled "Figure DR6. Detailed sedimentological description of the EM-137C core in the study area", measuring 0 to 242 meters. | Sahoo et al. (2016) Data Repository item DR6. |

### 1.4 Existing Scripts Loading, Aligning, Resampling, or Modeling Lithologs

| File Path | Language | Key Classes / Functions | Operational Functionality & Assumptions |
| :--- | :---: | :--- | :--- |
| [`data/loader.py`](file:///d:/Lithology-reconstruction-using-XGB/data/loader.py) | Python | `LithologLoader`, `standardize_dataframe()`, `check_layer_anomalies()` | Ingests raw CSVs, enforces positive thickness, discretizes to 1m regular steps, resolves overlaps with `keep='last'`, generates synthetic GR. **Critical finding:** `_extract_strike_pos()` assigns an artificial synthetic coordinate `(num - 1) * 100.0` from litholog numbers. |
| [`smalt/data/loader.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/data/loader.py) | Python | `LithologLoader` alias | Package re-export of `data/loader.py`. |
| [`scripts/run_phase1_markov.py`](file:///d:/Lithology-reconstruction-using-XGB/scripts/run_phase1_markov.py) | Python | Execution script | Ingests `data/processed/lithologs_unified.parquet`, executes 1D Markov chain fitting, computes stationary distributions, exports JSON summary. |
| [`smalt/geostat/markov.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/geostat/markov.py) | Python | `StratigraphicMarkovChain`, `fit()`, `compute_stationary_distribution()` | Calculates 1D vertical transition matrices. Reverses sequence per well (decreasing depth from base to top) to enforce upward depositional fining. Does not model spatial or inter-well relationships. |
| [`smalt/viz/markov_viz.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/viz/markov_viz.py) | Python | `plot_transition_matrix()`, `plot_stationary_vs_empirical()` | Generates annotated 300-DPI heatmaps and network succession diagrams. |
| [`src/load_data.py`](file:///d:/Lithology-reconstruction-using-XGB/src/load_data.py) | Python | `load_data()` | Legacy script loading `data/log.csv` for initial XGBoost experiment. |
| [`tests/test_phase0_loader.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_phase0_loader.py) | Python | Pytest suite (10 tests) | Validates Phase 0 schema compliance, non-negative depth assertions, and Parquet round-trip integrity. |
| [`tests/test_phase1_markov.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_phase1_markov.py) | Python | Pytest suite (11 tests) | Validates Markov matrix row-stochasticity, eigensolver convergence, and embedded chain zero diagonals. |

### 1.5 Existing Documentation Describing Datums, Coordinates, Groups, or Grids

| File Path | Relevant Section / Content | Stratigraphic / Spatial Findings |
| :--- | :--- | :--- |
| [`docs/Stratigraphic_Reservoir_Characterization_Report.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/Stratigraphic_Reservoir_Characterization_Report.md) | Sections 3, 4, 5 | Identifies a prominent "Chronostratigraphic Coal Datum (~37 - 44 m)" across L1-L7 and L11; secondary coal at 4-18m; proposes hypothetical 3D grid 6.0 km × 1.0 km × 93 m. |
| [`docs/FACULTY_DISCUSSION_BRIEF.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/FACULTY_DISCUSSION_BRIEF.md) | Sections 2, 3, 5 | Summarizes Phase 0 schema and Phase 1 1D Markov engine; previews Phase 2 2D object-based fluvial modeling with channel aspect ratio $W/T \approx 35$. |
| [`docs/notes/phase0_schema.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/notes/phase0_schema.md) | Sections 2, 3, 5 | Documents 5-facies state space, synthetic GR model, 1m discretization, and L9 overlap resolution. |
| [`docs/notes/phase1_markov.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/notes/phase1_markov.md) | Sections 1 to 5 | Mathematical formulations for regular ($P_{\text{reg}}$) and embedded ($P_{\text{emb}}$) Markov transition probability matrices. |
| [`docs/Sahoo et al 2016.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/Sahoo%20et%20al%202016.md) | Sections 2, 3, 4, Figs 3, 6-9 | Foundational paper. Notes 8 contiguous cliff faces in Cottonwood Creek, Wasatch Plateau; mean paleocurrent $47^\circ$; "Star Point top" is the basal regional stratigraphic contact; Axel Anderson coal zone. |
| [`docs/SMALT_Study_and_Build_Plan.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/SMALT_Study_and_Build_Plan.md) | Tracks A-D (Phases 0-8) | Repository roadmap outlining 2D cross-sectional simulation, inter-well conditioning, and active margin sampling. |
| [`lolo/litholog_locations.pptx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/litholog_locations.pptx) | Slide 0 | Documents explicit upstream group (L2, L3, L4, L5, L6, L7, L8, L10) and downstream group (L12, L9, L11). Notes red bold numbers 1-8 are canyon numbers. |

---

## 2. Coordinate Audit

Coordinates were extracted directly from the professor's spreadsheet ([`lolo/Location_coordinates_lithologs.xlsx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/Location_coordinates_lithologs.xlsx)) and verified against slide annotations in [`lolo/litholog_locations.pptx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/litholog_locations.pptx).

### 2.1 Proposed Spatial Metadata Table

The audited coordinates are structured into the proposed schema and saved as a new artifact at [`audit/spatial_metadata_table.csv`](file:///d:/Lithology-reconstruction-using-XGB/audit/spatial_metadata_table.csv):| litholog_id | x | y | coordinate_reference_system | coordinate_units | upstream_downstream_group | coordinate_source | coordinate_available |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **L1** | `NaN` | `NaN` | Undocumented | Assumed Meters (Excluded from spatial) | Permanently Excluded | None (absent from coordinate spreadsheet; excluded from spatial) | **False** |
| **L2** | -4057.6840 | 13755.6090 | Undocumented | Assumed Meters | Upstream | `Location_coordinates_lithologs.xlsx` | **True** |
| **L3** | -3577.0463 | 13317.0003 | Undocumented | Assumed Meters | Upstream | `Location_coordinates_lithologs.xlsx` | **True** |
| **L4** | -3121.5293 | 12356.8640 | Undocumented | Assumed Meters | Upstream | `Location_coordinates_lithologs.xlsx` | **True** |
| **L5** | -2709.2719 | 12276.4754 | Undocumented | Assumed Meters | Upstream | `Location_coordinates_lithologs.xlsx` | **True** |
| **L6** | -2737.7434 | 11361.2388 | Undocumented | Assumed Meters | Upstream | `Location_coordinates_lithologs.xlsx` | **True** |
| **L7** | -2163.3639 | 10733.2563 | Undocumented | Assumed Meters | Upstream | `Location_coordinates_lithologs.xlsx` | **True** |
| **L8** | -1476.1128 | 10873.8634 | Undocumented | Assumed Meters | Upstream | `Location_coordinates_lithologs.xlsx` | **True** |
| **L9** | -137.3543 | 12266.4076 | Undocumented | Assumed Meters | Downstream | `Location_coordinates_lithologs.xlsx` | **True** |
| **L10** | -469.4270 | 11082.7203 | Undocumented | Assumed Meters | Upstream | `Location_coordinates_lithologs.xlsx` | **True** |
| **L11** | 704.0729 | 12169.9234 | Undocumented | Assumed Meters | Downstream | `Location_coordinates_lithologs.xlsx` | **True** |
| **L12** | -1119.9100 | 14407.3000 | Undocumented | Assumed Meters | Downstream | `Location_coordinates_lithologs.xlsx` | **True** |

### 2.2 Integrity Checks and Observations

1. **Duplicate Identifiers**:
   - Zero duplicate IDs were detected in the coordinate spreadsheet. Each entry from L2 to L12 is unique.
2. **Missing Coordinates (Litholog 1 - Excluded from Spatial Modeling)**:
   - **Litholog 1 is absent** from `Location_coordinates_lithologs.xlsx`.
   - In accordance with project decisions, Litholog 1 is permanently excluded from 2D/3D lateral spatial modeling (no spatial coordinates are required or requested), while remaining fully validated for 1D vertical descriptive and succession analysis.
3. **Duplicate Coordinate Pairs**:
   - Zero duplicate $(X, Y)$ coordinate pairs exist. All 11 recorded positions are distinct spatial locations.
4. **Data Types and Non-Numeric Values**:
   - All 11 coordinate rows contain valid floating-point numeric values.
5. **Coordinate Reference System (CRS) & Projection**:
   - **Undocumented**: Neither `Location_coordinates_lithologs.xlsx` nor `litholog_locations.pptx` specifies a CRS, EPSG code, geodetic datum (e.g., NAD27, NAD83, WGS84), or map projection (e.g., UTM Zone 12N, State Plane Utah Central).
   - Numerical values ($X \in [-4057.7, 704.1]$, $Y \in [10733.3, 14407.3]$) do not match standard UTM northing/easting (where Wasatch Plateau Easting is $\sim 480,000\text{ m}$ and Northing is $\sim 4,360,000\text{ m}$). They appear to be a local project origin, local mine grid, or relative offset coordinates.
6. **Coordinate Units**:
   - Header text: `X-coordionate` [sic] and `Y-coordinate`.
   - The Euclidean span from L2 to L11 is $\Delta X = 4761.8\text{ units}$, $\Delta Y = -1585.7\text{ units}$, yielding a transect length of $\sqrt{4761.8^2 + 1585.7^2} = 5018.8\text{ units}$. Consistent with Sahoo et al. (2016) reporting a "$\sim 6\text{ km strike-transect}$" (with L1 extending further northwest), the physical units are officially assumed to be meters for all modeling tracks. Foot units would yield only $1.5\text{ km}$, directly contradicting the publication geometry.
7. **Spatial Ordering vs. Litholog Numbering**:
   - **Litholog numerical indices do not follow spatial ordering**:
     - L6 ($X = -2737.74$) is located west of L5 ($X = -2709.27$).
     - L10 ($X = -469.43, Y = 11082.72$) is located south-southwest of L9 ($X = -137.35, Y = 12266.41$).
     - L12 ($X = -1119.91, Y = 14407.30$) is located far to the north.
     - In `data/loader.py`, line 140 previously calculated `strike_pos_m = float((num - 1) * 100.0)`. This audit confirms that this synthetic assumption was an unverified heuristic that violated the true spatial configuration.
8. **Nature and Identity of Litholog 12**:
   - L12 corresponds to the **EM-137C core** (Sahoo et al., 2016, Supplementary Data Repository Figure DR6; [`lolo/litholog12.pdf`](file:///d:/Lithology-reconstruction-using-XGB/lolo/litholog12.pdf)), a subsurface drill core located $\sim 2.3\text{ km}$ north-northeast of the outcrop cliff transect.
   - The Blackhawk Formation interval (0.0 to 111.0 m) has been digitized and verified as tabular data in [`data/raw_lithologs/litholog12.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog12.csv) (63 contiguous beds).t.

---

## 3. Vertical Reference and Log-Length Audit

Each of the 11 litholog CSV files in [`data/raw_lithologs/`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/) was audited individually for vertical formatting, continuity, interval thickness, and stratigraphic datuming.

### 3.1 Litholog Vertical Attribute Table

| Litholog ID | Row Count | Available Columns | Min Depth (m) | Max Depth (m) | Stratigraphic Span (m) | Depth Direction | Min Facies Interval (m) | Max Facies Interval (m) | Sampling Structure | Detected Gaps (Top > Prev Bottom) | Detected Overlaps (Top < Prev Bottom) | Missing Values | Aligned to Common Datum? |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **litholog1** | 20 | `Top`, `Bottom`, `Facies` | 0 | 93 | 93 | Increasing downward | 1 | 17 | Irregular bed intervals | None | None | 0 | **Unresolved** |
| **litholog2** | 16 | `Top`, `Bottom`, `Facies` | 0 | 93 | 93 | Increasing downward | 1 | 18 | Irregular bed intervals | None | None | 0 | **Unresolved** |
| **litholog3** | 15 | `Top`, `Bottom`, `Facies` | 0 | 93 | 93 | Increasing downward | 1 | 12 | Irregular bed intervals | None | None | 0 | **Unresolved** |
| **litholog4** | 20 | `Top`, `Bottom`, `Facies` | 0 | 84 | 84 | Increasing downward | 1 | 11 | Irregular bed intervals | None | None | 0 | **Unresolved** |
| **litholog5** | 20 | `Top`, `Bottom`, `Facies` | 0 | 85 | 85 | Increasing downward | 1 | 12 | Irregular bed intervals | None | None | 0 | **Unresolved** |
| **litholog6** | 17 | `Top`, `Bottom`, `Facies` | 0 | 82 | 82 | Increasing downward | 1 | 13 | Irregular bed intervals | None | None | 0 | **Unresolved** |
| **litholog7** | 16 | `Top`, `Bottom`, `Facies` | 0 | 80 | 80 | Increasing downward | 2 | 11 | Irregular bed intervals | None | None | 0 | **Unresolved** |
| **litholog8** | 28 | `Top`, `Bottom`, `Facies` | 0 | 79 | 79 | Increasing downward | 1 | 7 | Irregular bed intervals | None | None | 0 | **Unresolved** |
| **litholog9** | 20 | `Top`, `Bottom`, `Facies` | 0 | 78 | 78 | Increasing downward | 1 | 11 | Irregular bed intervals | `(18, 19)` -> `c_sand` | `(29, 28)`, `(30, 29)` -> `p_sand` | 0 | **Unresolved** |
| **litholog10** | 19 | `Top`, `Bottom`, `Facies` | 0 | 77 | 77 | Increasing downward | 1 | 10 | Irregular bed intervals | None | None | 0 | **Unresolved** |
| **litholog11** | 18 | `Top`, `Bottom`, `Facies` | 0 | 78 | 78 | Increasing downward | 1 | 12 | Irregular bed intervals | `(59, 60)` -> `carbon_mud` | None | 0 | **Unresolved** |
| **litholog12** | 63 | `Top`, `Bottom`, `Facies` | 0.0 | 111.0 | 111.0 | Increasing downward | 0.2 | 8.5 | Irregular bed intervals | None | None | 0 | **Unresolved** |

*Note on Litholog 12 ([`data/raw_lithologs/litholog12.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog12.csv)):* Digitized tabular CSV is now available covering the 0.0 to 111.0 m Blackhawk Formation core interval (63 contiguous beds). The source vector core graphic ([`lolo/litholog12.pdf`](file:///d:/Lithology-reconstruction-using-XGB/lolo/litholog12.pdf)) extends to 242 m, with the deeper 111-242 m interval representing underlying regional strata.

### 3.2 Vertical Reference & Datum Analysis

1. **Depth Convention & Direction**:
   - In all 11 lithologs, `Top` and `Bottom` values are strictly positive integers. Depths increase downward, with `Top = 0` assigned to the uppermost bed of each section.
2. **Top-of-Log Geological Inconsistency**:
   - In all 11 logs, the top layer (`0 m` to $4-11\text{ m}$) is recorded as `mud`. However, in outcrop exposures along the Wasatch Plateau escarpment, the physical top of a cliff represents the modern erosional landscape surface, not a single chronostratigraphic depositional surface.
3. **Bottom-of-Log Stratigraphic Asymmetry**:
   - The base of each litholog terminates at disparate vertical depths and in completely different lithofacies:
     - Lithologs 1, 2, and 3 terminate at depth $93\text{ m}$ in massive channel sandstone bodies ($11\text{ m}$, $10\text{ m}$, and $11\text{ m}$ thick, respectively).
     - Lithologs 4, 5, 6, 7, 8, 9, and 10 terminate at depths ranging from $77\text{ m}$ to $85\text{ m}$ in overbank mudstones.
     - Litholog 11 terminates at depth $78\text{ m}$ in siltstone.
4. **Coal Marker Distribution**:
   - While Sahoo et al. (2016) and [`docs/Stratigraphic_Reservoir_Characterization_Report.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/Stratigraphic_Reservoir_Characterization_Report.md) highlight a prominent coal marker bed at $\sim 37-44\text{ m}$ (observed in L2, L3, L4, L5, L6, L7, and L11), this marker:
     - Is **completely absent** in Lithologs 8, 9, and 10.
     - Is not present at $\sim 37-44\text{ m}$ in Litholog 1 (which records coals at 5-6m, 9-10m, 12-14m, 17-19m, 58-59m, and 81-82m).
5. **Stratigraphic Datum Status**:
   - In Sahoo et al. (2016) Figures 6, 7, 8, and 9, outcrop cliff panels are explicitly correlated and flattened against the **"Star Point top"** (the regional contact between the underlying Star Point Sandstone and the basal Blackhawk Formation).
   - In the raw CSVs, however, **no elevation above sea level, collar elevation, or distance to the Star Point top is recorded**.
   - **Conclusion**: The dataset **does not appear aligned to a common stratigraphic datum**. Treating equal numerical depths in different lithologs ($z_A = z_B$) as identical stratigraphic horizons is scientifically invalid and remains **UNRESOLVED**.

---

## 4. Spatial Map

A provisional spatial map of all documented litholog locations was generated using the professor's coordinates and grouping metadata. The plot honors equal axis scaling, clearly annotates litholog identifiers, and highlights upstream and downstream domains without interpolating or fabricating coordinates.

![Provisional Spatial Map of Litholog Locations](file:///d:/Lithology-reconstruction-using-XGB/audit/provisional_coordinate_map.png)

*Figure 1: Provisional coordinate map of known litholog locations (L2-L12) based on [`lolo/Location_coordinates_lithologs.xlsx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/Location_coordinates_lithologs.xlsx) and [`lolo/litholog_locations.pptx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/litholog_locations.pptx). Axes represent recorded coordinates (presumed meters; CRS undocumented). Axis aspect ratio is strictly 1:1.*

### 4.1 Observations from the Spatial Configuration

1. **Upstream vs. Downstream Alignment**:
   - Upstream lithologs (Red: L2, L3, L4, L5, L6, L7, L8, L10) form an elongated arcuate escarpment path extending from the northwest ($X \approx -4058, Y \approx 13756$) toward the south-southeast ($X \approx -2163, Y \approx 10733$), then turning east-northeast toward L8 and L10.
   - Downstream lithologs (Green: L9, L11, L12) lie distinctly down-paleoflow. L9 and L11 reside on the eastern canyon flank, while L12 (EM-137C core) is located $2.3\text{ km}$ north-northeast of the main transect.
2. **Consistency with Published Paleocurrent**:
   - In Sahoo et al. (2016) Figure 3C, paleocurrent measurements exhibit a unimodal distribution with a vector mean of **$47^\circ$ (Northeast)**.
   - Sediment transport flowed from the southwest toward the northeast. The professor's classification of L2-L8 and L10 as "Upstream" and L9, L11, and L12 as "Downstream" perfectly mirrors this physical paleoflow geometry.
3. **Provisional Status**:
   - Because the CRS and geodetic datum remain undocumented, this map must remain designated as **provisional** until the coordinate projection is verified by Prof. Sahoo.

---

## 5. Shortest-Log Assessment

In exploratory workflows, it is sometimes suggested to truncate all logs to the length of the shortest observed log to force a uniform rectangular data matrix. Here, we evaluate this proposal strictly without truncating any data.

### 5.1 Observed Depth Range Disparity

The 11 outcrop lithologs exhibit substantial depth range variation:
- **Longest logs**: Lithologs 1, 2, and 3 span **$93.0\text{ m}$**.
- **Intermediate logs**: Litholog 5 ($85.0\text{ m}$), Litholog 4 ($84.0\text{ m}$), Litholog 6 ($82.0\text{ m}$), Litholog 7 ($80.0\text{ m}$), Litholog 8 ($79.0\text{ m}$).
- **Shortest logs**: Litholog 9 ($78.0\text{ m}$ span, $77\text{ m}$ valid), Litholog 11 ($78.0\text{ m}$ span, $77\text{ m}$ valid), and Litholog 10 (**$77.0\text{ m}$** span).

### 5.2 Methodological Comparison

#### (a) Retaining Each Complete Litholog for Model Estimation
- **Advantages**:
  - Preserves 100% of collected geological observations (all 920 validated 1-meter intervals).
  - Essential for 1D vertical Markov chain estimation (Phase 1), which models upward succession transition probabilities independently within each vertical section.
  - Retains the massive basal channel sandstones in Lithologs 1, 2, and 3 (depths $82-93\text{ m}$), which constitute the primary high-net-to-gross reservoir sweet spot of the entire system.
- **Limitations**:
  - Cannot be directly placed onto an unaligned 2D/3D regular simulation grid without stratigraphic datum correction.

#### (b) Restricting Evaluation to Shortest Log ($77\text{ m}$) vs. Verified Common Stratigraphic Interval
- **Arbitrary Truncation to Shortest Log ($77\text{ m}$)**:
  - If truncated from the top ($z = 0$ to $77\text{ m}$), this would **discard the bottom $16\text{ meters}$ of Lithologs 1, 2, and 3**.
  - In L1, L2, and L3, the interval from $82\text{ m}$ to $93\text{ m}$ consists of $10-11\text{ meters}$ of continuous, porous channel sandstone. Truncating at $77\text{ m}$ would systematically delete the thickest reservoir fairways from the western sector, severely distorting the Net-to-Gross sand budget (reducing L3 N/G from $73.1\%$ to an artificially lower value).
  - Because numerical depth $z$ does not represent a chronostratigraphic horizon, cutting at $77\text{ m}$ has no geological meaning.
- **Restricting to a Verified Common Stratigraphic Interval**:
  - Restricting modeling to a genetically bounded interval (e.g., between the regional Star Point Sandstone top and a continuous chronostratigraphic coal seam) is standard sequence stratigraphic practice.
  - However, doing so requires knowing the exact depth or elevation of those bounding surfaces in every well.

### 5.3 Required Alignment Information Before Cross-Well Comparisons

Before any cross-well correlation, horizontal variogram estimation, or 2D/3D spatial interpolation can be considered scientifically valid, the following alignment information is strictly required:
1. **Stratigraphic Anchor Horizon**: A confirmed common geological datum (such as the top of the Star Point Sandstone or a continuous regional coal seam) identified by depth or relative offset in every well.
2. **Collar Elevation / Topography**: Absolute ground elevation (meters above sea level) or relative cliff-top elevation for each litholog location.
3. **Structural Dip and Azimuth**: Regional tectonic dip (known to be gently eastward/northeastward in the Wasatch Plateau) must be removed to structurally flatten the succession before calculating lateral facies continuity.
4. **Wheeler Domain Transform / Stratigraphic Normalization**: Proportional slicing or chronostratigraphic flattening between bounding datums to prevent artificial cross-stratal smearing during geostatistical simulation.

---

## 6. Actionable Status and Unresolved Questions

To enable progress from 1D geostatistics (Phase 1) to 2D/3D spatial simulation (Phase 2), the status of key technical questions is summarized below, distinguishing between resolved project decisions and remaining items requiring Professor Hiranya Sahoo's input:

### 6.1 Resolved Technical Decisions

1. **Coordinate Units**:
   - **Resolved**: Coordinate units are confirmed and assumed as **international meters** across all datasets. This is fully consistent with the $\sim 5\text{ km}$ transect length between L2 and L11 documented by Sahoo et al. (2016).
2. **Litholog 1 Spatial Coordinates**:
   - **Resolved**: Litholog 1 has no spatial coordinates and is **permanently excluded from 2D/3D lateral spatial modeling**. It is retained strictly for 1D vertical descriptive and Markov succession modeling. No additional spatial coordinates are required or requested.
3. **Availability of Litholog 12 CSV**:
   - **Resolved**: The tabular CSV for Litholog 12 has been digitized and verified at [`data/raw_lithologs/litholog12.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog12.csv), covering the 0.0 to 111.0 m Blackhawk Formation interval (63 contiguous beds). Tabular data is already available.
4. **Discontinuities in Raw Data (Lithologs 9 and 11)**:
   - **Resolved**:
     - **Litholog 9 (18-19 m gap)**: Taken as `c_sand` (channel sandstone).
     - **Litholog 9 (28-30 m overlap)**: Taken as `p_sand` (planar / splay sandstone).
     - **Litholog 11 (59-60 m unmapped gap)**: Taken as `carbon_mud` (carbonaceous mudstone).
     - *(Note: Downstream encodings and manual updates to raw data files will be reconciled accordingly).*

### 6.2 Remaining Unresolved Questions for Professor Hiranya Sahoo

1. **Coordinate Reference System (CRS) & Projection**:
   - What specific coordinate reference system or projection governs the numeric coordinates in `Location_coordinates_lithologs.xlsx` (e.g., UTM Zone 12N with a local false origin, State Plane Utah Central, or a local coal mine survey grid)?
2. **Role of Litholog 12 (EM-137C Core)**:
   - Is Litholog 12 intended to be included as a hard conditioning well in the spatial simulation grid, or is it reserved for regional down-dip blind testing?
   - Does the lower interval of the EM-137C core (111-242 m in `litholog12.pdf`) represent underlying non-Blackhawk stratigraphy (Star Point Sandstone / Mancos Shale), or is it pending digitization?
3. **Vertical Reference Datum & Elevation**:
   - What physical surface corresponds to `Depth = 0 m` for each litholog (e.g., local erosional cliff top, base of an overlying formation, or modern ground surface)?
   - Are absolute elevations (meters above sea level) or relative heights above the Star Point Sandstone contact available for each litholog?
4. **Coal Marker Horizon Correlation**:
   - Does the coal seam observed at $\sim 37-44\text{ m}$ in Lithologs 2-7 and 11 represent a single continuous chronostratigraphic marker bed (e.g., Axel Anderson coal zone)? If so, what accounts for its absence in Lithologs 8, 9, and 10?

---

## 7. Readiness Decision

### 7.1 Decision Summary

| Modeling Track | Readiness Status | Justification & Preconditions |
| :--- | :---: | :--- |
| **Phase 1: 1D Vertical Markov Modeling** | **READY / VERIFIED** | Validated on independent vertical successions; unaffected by lateral spatial coordinates or inter-well datuming. All unit tests passing. |
| **Phase 2/3: 2D/3D Spatial & Conditional Simulation** | **NOT READY / BLOCKED** | **BLOCKED** by undocumented coordinate reference system (CRS/projection) and lack of a verified vertical stratigraphic datum. (Coordinate units assumed meters; L1 excluded from spatial; L12 CSV available). |

### 7.2 Detailed Readiness Rationale

The dataset in its current state is **NOT READY for 2D or 3D spatial simulation, inter-well spatial interpolation, or object-based conditioning**.

While the 1D succession statistics are mathematically verified, attempting to define a spatial simulation grid ($X, Y, Z$) under present conditions would require making unverified assumptions:
1. Treating $z = 0$ as a flat plane would distort geological bodies across the $5\text{ km}$ transect, causing channels deposited at different geological times to artificially intersect or terminate.
2. Ingesting $X$ and $Y$ without a verified CRS prevents re-projection, spatial integration with regional GIS/seismic data, or valid distance calculations against true north.
3. Litholog 1 is permanently excluded from spatial modeling due to missing coordinates, focusing 2D/3D spatial modeling on the remaining transect (L2-L12).

### 7.3 Exact Additional Metadata Required Before Defining Simulation Grid

Before any spatial simulation grid can be initialized, the following items must be provided and verified:
1. **Verified CRS and Geodetic Datum** (EPSG code or local coordinate grid origin definition).
2. **Stratigraphic datum elevation or structural contact depth** (e.g., depth to Star Point Sandstone top or a regional marker coal) for each section.

*(Note: Coordinate units are resolved as assumed meters, Litholog 1 is resolved as excluded from spatial modeling, and Litholog 12 CSV is resolved and available at [`data/raw_lithologs/litholog12.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog12.csv)).*

---

## Repository Audit Manifest & Created Files

In accordance with strict read-only audit constraints, **zero existing repository files were modified, deleted, renamed, or overwritten**.

The following new audit deliverables were produced and stored in the newly created audit directory:
1. [`audit/spatial_metadata_table.csv`](file:///d:/Lithology-reconstruction-using-XGB/audit/spatial_metadata_table.csv) - Proposed spatial metadata table containing audited coordinates, CRS status, units, upstream/downstream groups, and availability flags for L1 through L12.
2. [`audit/provisional_coordinate_map.png`](file:///d:/Lithology-reconstruction-using-XGB/audit/provisional_coordinate_map.png) - High-resolution, publication-quality provisional coordinate map with equal aspect scaling, color-coded upstream/downstream groups, and clear provenance annotations.
3. [`audit/SPATIAL_VERTICAL_DATA_AUDIT_REPORT.md`](file:///d:/Lithology-reconstruction-using-XGB/audit/SPATIAL_VERTICAL_DATA_AUDIT_REPORT.md) - This comprehensive audit report.

*Audit completed on: 2026-10-05 | SMALT Geostatistical Data Engineering*
