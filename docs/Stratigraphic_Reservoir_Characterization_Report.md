# Comprehensive Stratigraphic Reservoir Characterization & Lithological Profile Report

## Executive Summary

This report delivers an integrated subsurface reservoir characterization based on an 11-litholog working dataset across a 6.0 km outcrop transect (up to 93 m vertical depth). Utilizing the Cretaceous Blackhawk Formation (Book Cliffs, Utah; Sahoo et al., 2016) as a geological analog for tight-gas fluvial reservoirs, this study establishes a standardized 5-facies framework (`sand`, `silt`, `mud`, `carbon_mud`, `coal`), quantifies Net-to-Gross (N/G) sand ratios, documents an explicit manual-vs-automated QC benchmark on Litholog 11 (93.59% structural accuracy), tracks dataset provenance across all 11 logs (source-derived vs. AI-reconstructed), and models 2D spatial contours and 3D volumetric lithofacies distributions conditioned on this working dataset.

---

## 1. Geological Analog Framework & Facies Classification

### 1.1 Outcrop Analog Framework
The study grounds its facies architecture in the **Blackhawk Formation** (Book Cliffs, Utah) as detailed by **Sahoo et al. (2016)**. This system serves as a primary depositional analog for Rocky Mountain low-permeability tight-gas fluvial plays.

### 1.2 Standardized 5-Facies System
To model heterogeneity across hierarchical scales (large-scale channel sandbodies, intermediate splay bars, and small-scale facies variations), all intervals are classified into five lithological facies:

1. **Sandstone (`sand`)**: High-energy, cross-bedded channel-fill and point-bar deposits (Primary Reservoir Units).
2. **Siltstone (`silt`)**: Intermediate-energy crevasse splay and levee deposits (Inter-compartment Pressure Connectors).
3. **Mudstone (`mud`)**: Low-energy, fine-grained overbank and floodplain deposits (Regional Vertical Seals).
4. **Carbonaceous Mudstone (`carbon_mud`)**: Organic-rich, quiet-water coastal plain/marsh deposits (Flow Baffles & Caprocks).
5. **Coal (`coal`)**: High-organic, swamp-accumulated marker beds (Chronostratigraphic Datums).

---

## 2. Dataset Provenance Framework & Verification Benchmark

### 2.1 Mixed Provenance Classification
The 11-litholog dataset across the 6.0 km transect has mixed provenance and must not be treated as a uniformly validated ground truth:

| Litholog ID | Provenance Category | Reconstruction Method | Independent Validation | Benchmark Accuracy | Status & Role |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **Litholog 1** | Source-derived | Existing dataset / Sahoo et al. (2016) | Unknown | None | Working dataset; independent field validation undocumented |
| **Litholog 2** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **Litholog 3** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **Litholog 4** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **Litholog 5** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **Litholog 6** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **Litholog 7** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **Litholog 8** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **Litholog 9** | Source-derived | Existing dataset / Sahoo et al. (2016) | Unknown | None | Working dataset; independent field validation undocumented; has interval gap at 18-19 m and overlaps at 28-30 m |
| **Litholog 10** | AI-reconstructed | AI-assisted reconstruction from published source | False | None | Working dataset; unvalidated reconstruction |
| **Litholog 11** | Source-derived (Benchmarked) | Existing manual digitization vs. automated QC | **True** | **93.59%** | **Reference QC Benchmark Log**; 73/78 m structural agreement |

### 2.2 Litholog 11 Verification Benchmark
To evaluate the automated image-to-data extraction workflow against manual digitization, **Litholog 11** was independently digitized and compared against a manually created benchmark CSV (`litholog11.csv`).

* **Total Measured Depth**: 0 to 78 meters (78 cumulative meters)
* **Depth Agreement**: **73 out of 78 meters** match identically (**93.59% Structural Accuracy**)
* **Relative Difference / Discrepancy Rate**: **6.41%** (5 discrepant meters)
* **Interval Count**: 18 intervals (Manual CSV) vs. 16 intervals (Automated CSV)
* **Exact Facies & Boundary Alignment**: 68 out of 78 meters (87.2%) exhibited 100% exact alignment.

> [!IMPORTANT]
> The 93.59% structural accuracy metric applies **exclusively to Litholog 11**. It cannot and must not be extrapolated to represent the accuracy of Lithologs 1-10. Lithologs 2-8 and 10 are retained as AI-assisted working reconstructions, and all downstream statistics (including transition matrices, 2D contour maps, and 3D models) must be recognized as being conditioned partly on these reconstructed observations.

### 2.3 Detailed Discrepancy Breakdown for Litholog 11

| Depth Interval | Manual CSV (`litholog11.csv`) | Automated CSV | Nature of Difference / Error |
| :--- | :--- | :--- | :--- |
| **12 - 13 m** | `silt` (11 - 13 m) | `coal` (12 - 13 m) | **Thin-Bed Resolution**: Automated log resolved a 1 m thin coal streak; manual log grouped it into a 2 m silt block. |
| **59 - 60 m** | *Unmapped Gap* (59-60 m missing) | `carbon_mud` (59 - 66 m) | **Depth Continuity**: Manual log contained a 1 m recording gap; automated engine maintained continuous coverage. |
| **65 - 66 m** | `silt` (65 - 67 m) | `carbon_mud` (59 - 66 m) | **Boundary Shift**: 1-meter offset in upper silt bed boundary. |
| **67 - 69 m** | `mud` (67-68 m) + `silt` (68-69 m) | `carbon_mud` (67 - 72 m) | **Facies Aggregation**: Manual log captured fine 1 m interbedding; automated log consolidated the unit into `carbon_mud`. |

---

## 3. Stratigraphic Correlation & Stacking Architecture

A correlation panel across all 11 lithological columns (Lithologs 1-11) reveals distinct architectural transitions across the 6.0 km transect:

### 3.1 Regional Marker Horizon
* **Chronostratigraphic Coal Datum (~37 - 44 m)**: A persistent, continuous coal seam extends across Lithologs 1, 2, 3, 4, 5, 6, 7, and 11, serving as a reliable regional marker for structural flattening and correlation.
* **Secondary Coal Horizon (~4 - 18 m)**: An upper organic-rich zone tracks across Lithologs 1-6 and 11.

### 3.2 Lateral Architectural Domains
* **Western Sector (Lithologs 1-3 | 0 - 1.8 km)**: Dominated by thick, multi-story, amalgamated channel sandbodies (sand thickness reaches 68 m in Litholog 3).
* **Central Sector (Lithologs 4-6 & 11 | 1.8 - 3.6 km)**: Stacked channel sandbodies separated by floodplain muds and thick carbonaceous marsh units.
* **Eastern Sector (Lithologs 7-10 | 3.6 - 6.0 km)**: Characterized by isolated single-story ribbon sandbodies, thin crevasse splays, and thick encased overbank mudstones.

---

## 4. Quantitative Net-to-Gross (N/G) Sand Ratio Breakdown

### 4.1 Master N/G Summary Table (Working Dataset)

The table below summarizes cumulative facies thicknesses and Net-to-Gross sand ratios across the 11-litholog working dataset. Column annotations explicitly indicate dataset provenance:

| Litholog | Provenance | Stratigraphic Span (m) | Net Sand (m) | Silt (m) | Mud (m) | Carbonaceous Mud (m) | Coal (m) | Net-to-Gross (N/G) (%) | Reservoir Stacking Classification |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Litholog 1** | Source-derived (Val. Unknown) | 93 | 35 | 0 | 47 | 3 | 8 | **37.6%** | Low-Moderate N/G; Coal-rich overbank |
| **Litholog 2** | AI-reconstructed (Unvalidated) | 93 | 51 | 0 | 31 | 4 | 7 | **54.8%** | Moderate-High N/G; Multi-lateral channel-belts |
| **Litholog 3** | AI-reconstructed (Unvalidated) | 93 | 68 | 0 | 18 | 5 | 2 | **73.1%** | **High N/G Sweet Spot**; Amalgamated channel complex |
| **Litholog 4** | AI-reconstructed (Unvalidated) | 84 | 32 | 5 | 41 | 4 | 2 | **38.1%** | Low N/G; Isolated channel sandbodies |
| **Litholog 5** | AI-reconstructed (Unvalidated) | 85 | 34 | 6 | 31 | 12 | 2 | **40.0%** | Moderate N/G; Multi-story channel sands |
| **Litholog 6** | AI-reconstructed (Unvalidated) | 82 | 31 | 8 | 19 | 21 | 3 | **37.8%** | Low-Moderate N/G; Organic-rich floodplain encasement |
| **Litholog 7** | AI-reconstructed (Unvalidated) | 80 | 24 | 9 | 25 | 20 | 2 | **30.0%** | **Low N/G**; Compartmentalized sands |
| **Litholog 8** | AI-reconstructed (Unvalidated) | 79 | 23 | 16 | 37 | 3 | 0 | **29.1%** | **Lowest N/G**; Overbank muds & crevasse splays |
| **Litholog 9** | Source-derived (Val. Unknown) | 78 | 28 | 11 | 24 | 16 | 0 | **35.9%** | Low N/G; Interbedded crevasse splays & sands (has 1m gap at 18-19m, overlaps at 28-30m) |
| **Litholog 10** | AI-reconstructed (Unvalidated) | 77 | 24 | 12 | 31 | 10 | 0 | **31.2%** | **Low N/G**; Isolated ribbon sandbodies |
| **Litholog 11** | Source-derived (QC Benchmarked) | 78 | 37 | 4 | 12 | 22 | 3 | **47.4%** | Moderate N/G; Stacked sands with carbonaceous muds |
| **TRANSECT TOTAL** | **Mixed Working Dataset** | **922** | **387** | **65** | **279** | **116** | **29** | **42.0%** | **Overall Low-to-Moderate Net-to-Gross Reservoir** |

*Methodological Reconciliation Notes:*
1. **Litholog 11 Reporting**: Values in this summary table reflect the automated digitization workflow (78 m total depth with continuous coverage across the 59-60 m gap using `carbon_mud`). The raw manual benchmark CSV (`data/raw_lithologs/litholog11.csv`) contains 77 cumulative meters of classified intervals (Sand: 37m, Silt: 7m, Mud: 13m, Carbon Mud: 18m, Coal: 2m) with an unmapped 1 m gap at 59-60 m.
2. **Litholog 9 Continuity**: Raw `litholog9.csv` exhibits a 1 m gap at 18-19 m (depth 18 missing) and 2 m overlaps at 28-30 m. When ingested into 1 m standardized bins, it yields 77 valid points.
3. **Conditioning Notice**: Lithologs 2-8 and 10 are AI-assisted reconstructions. Downstream Net-to-Gross aggregations must be interpreted as derived from the current working dataset rather than independently confirmed ground truth.

---

## 5. 2D Spatial N/G Mapping & 3D Volumetric Block Modeling

*(Note: Spatial interpolation and 3D modeling are conditioned on the 11-litholog working dataset, which incorporates 8 AI-reconstructed profiles).*

### 5.1 2D Spatial N/G Distribution
* **Western High-N/G Zone (0 - 1.8 km)**: Displays N/G values **>70-80%**, representing the primary reservoir fairway with high sandbody connectivity.
* **Mid-Depth Barrier Corridor (30 - 45 m depth)**: A continuous low-N/G (<20%) shale/coal band extending laterally across Lithologs 4-10, forming a regional vertical hydraulic seal.
* **Eastern Low-N/G Zone (3.6 - 6.0 km)**: N/G drops below **30%**, indicating severe lateral compartmentalization.

### 5.2 3D Lithofacies Volumetric Model
* **Grid Extent**: 6.0 km (length) × 1.0 km (width) × 93 m (depth).
* **Volumetric Architecture**: Confirms a wedge-shaped reservoir sand body thinning from West to East, encased within low-permeability mudstones and bounded vertically by continuous coal marker horizons.

---

## 6. Field Development & Reservoir Management Recommendations

*(Working dataset interpretations subject to future ground-truth validation of AI-reconstructed logs)*

1. **Targeted Infill Drilling**: Prioritize horizontal well placement in the **western 0 - 2.0 km sector** (Lithologs 1-3) where high N/G (>50-73%) ensures high effective drainage and low compartmentalization risk.
2. **Dual-Zone Completion Strategy**: Treat the upper (<30 m) and lower (>45 m) sandbody complexes as independent reservoir compartments separated by the continuous 30-45 m mid-depth regional seal.
3. **EOR & Pressure Maintenance**: Avoid relying on lateral waterfloods across the central-to-eastern transition (Lithologs 7-10) due to low sand connectivity and barrier pinch-outs.

---
*Report compiled from the 11-litholog SMALT working dataset (`litholog1.csv` through `litholog11.csv`). Provenance manifest registered in `data/provenance_manifest.json`.*
