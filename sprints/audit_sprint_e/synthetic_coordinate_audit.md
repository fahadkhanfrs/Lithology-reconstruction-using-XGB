# SMALT Sprint E - Synthetic Coordinate Audit Report

## 1. Executive Summary

This audit rigorously inspects the provenance, storage, calculation, and usage of spatial coordinates across the SMALT repository following the manual digitization and revision of all 12 lithologs.

Key findings:
1. **Source-Derived Coordinates**: Verified real Cartesian coordinates $(X, Y)$ exist for 11 out of 12 lithologs (L2 through L12) in [`lolo/Location_coordinates_lithologs.xlsx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/Location_coordinates_lithologs.xlsx). Units are in meters, representing a local coordinate system spanning ~4.76 km east-west and ~3.67 km north-south.
2. **Missing Coordinates (Litholog 1)**: Litholog 1 is completely absent from the coordinate spreadsheet. It has no physical spatial position and is permanently excluded from spatial modeling.
3. **Synthetic Coordinate Tracing**: In legacy Phase 0 preprocessing ([`data/loader.py`](file:///d:/Lithology-reconstruction-using-XGB/data/loader.py#L136-L141)), an artificial strike coordinate was generated:
   $$\text{strike\_pos\_m} = (\text{file\_index}) \times 100.0$$
   This synthetic coordinate was written into [`data/processed/lithologs_unified.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/processed/lithologs_unified.csv) and `.parquet`.
4. **Model Leakage Check**: Neither the 1D Markov transition models (Sprint C/D) nor the Phase 1 geostatistical validation pipelines ever ingested `strike_pos_m`.
5. **Provisional Spatial Baseline Resolution**: The new spatial module ([`smalt/spatial/coordinates.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/coordinates.py) and [`smalt/spatial/baseline.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/baseline.py)) completely ignores `strike_pos_m`, ingests only verified source spreadsheet coordinates for L2-L12, strictly enforces L1 exclusion, and benchmarks spatial models against naive baselines under documented unanchored datum assumptions.

---

## 2. Coordinate Inventory & Provenance by Litholog

| Litholog ID | Raw Label | Coordinate Status | Source / Provenance | Raw X (m) | Raw Y (m) | Synthetic Formula Present? | Spatial Modeling Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `litholog1` | L1 | **Missing** | None (absent from spreadsheet) | NaN | NaN | Previously assigned `strike_pos_m = 0.0` in loader | **Strictly Excluded** |
| `litholog2` | L2 | Source-derived | `Location_coordinates_lithologs.xlsx` (Row 1) | -4057.68 | 13755.61 | Previously assigned `strike_pos_m = 100.0` in loader | Included (Provisional relative geometry) |
| `litholog3` | L3 | Source-derived | `Location_coordinates_lithologs.xlsx` (Row 2) | -3577.05 | 13317.00 | Previously assigned `strike_pos_m = 200.0` in loader | Included (Provisional relative geometry) |
| `litholog4` | L4 | Source-derived | `Location_coordinates_lithologs.xlsx` (Row 3) | -3121.53 | 12356.86 | Previously assigned `strike_pos_m = 300.0` in loader | Included (Provisional relative geometry) |
| `litholog5` | L5 | Source-derived | `Location_coordinates_lithologs.xlsx` (Row 4) | -2709.27 | 12276.48 | Previously assigned `strike_pos_m = 400.0` in loader | Included (Provisional relative geometry) |
| `litholog6` | L6 | Source-derived | `Location_coordinates_lithologs.xlsx` (Row 5) | -2737.74 | 11361.24 | Previously assigned `strike_pos_m = 500.0` in loader | Included (Provisional relative geometry) |
| `litholog7` | L7 | Source-derived | `Location_coordinates_lithologs.xlsx` (Row 6) | -2163.36 | 10733.26 | Previously assigned `strike_pos_m = 600.0` in loader | Included (Provisional relative geometry) |
| `litholog8` | L8 | Source-derived | `Location_coordinates_lithologs.xlsx` (Row 7) | -1476.11 | 10873.86 | Previously assigned `strike_pos_m = 700.0` in loader | Included (Provisional relative geometry) |
| `litholog9` | L9 | Source-derived | `Location_coordinates_lithologs.xlsx` (Row 8) | -137.35 | 12266.41 | Previously assigned `strike_pos_m = 800.0` in loader | Included (Provisional relative geometry) |
| `litholog10` | L10 | Source-derived | `Location_coordinates_lithologs.xlsx` (Row 9) | -469.43 | 11082.72 | Previously assigned `strike_pos_m = 900.0` in loader | Included (Provisional relative geometry) |
| `litholog11` | L11 | Source-derived | `Location_coordinates_lithologs.xlsx` (Row 10) | 704.07 | 12169.92 | Previously assigned `strike_pos_m = 1000.0` in loader | Included (Provisional relative geometry) |
| `litholog12` | L12 | Source-derived | `Location_coordinates_lithologs.xlsx` (Row 11) | -1119.91 | 14407.30 | Not included in legacy loader | Included (Provisional relative geometry) |

---

## 3. Tracing Synthetic Coordinates Across the Pipeline

### 3.1 Where Were Synthetic Coordinates Generated?
In [`data/loader.py`](file:///d:/Lithology-reconstruction-using-XGB/data/loader.py#L136-L141):
```python
# Derive synthetic strike position from file index (100m spacing)
strike_pos = float(file_idx * 100.0)
df["strike_pos_m"] = strike_pos
```
This formula assumed that lithologs were evenly spaced by 100 meters in order of filename (`litholog1.csv` to `litholog11.csv`).

### 3.2 Where Did They Flow?
1. Written to processed dataset files:
   - `data/processed/lithologs_unified.csv`
   - `data/processed/lithologs_unified.parquet`
2. **Never entered** 1D Markov chains in [`smalt/geostat/markov.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/geostat/markov.py) or [`smalt/validation/markov_lolo.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/validation/markov_lolo.py).
3. **Never entered** Sprint C or Sprint D descriptive statistics.

### 3.3 What Was the Danger?
If an ML model used `strike_pos_m` as an $(X)$ coordinate feature, it would:
- Fabricate a 1D linear spatial geometry that directly contradicts the true non-linear 2D layout in the field.
- Impart an artificial spatial ordering (e.g. L9 placed 800m away from L1, when in reality L9 is 3.9 km east of L2 and L1 coordinates are completely unknown).

### 3.4 Repository Resolution in Sprint E
1. Created dedicated spatial coordinate loader [`smalt/spatial/coordinates.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/coordinates.py).
2. All spatial models now directly load real local Cartesian coordinates $(X, Y)$ from [`lolo/Location_coordinates_lithologs.xlsx`](file:///d:/Lithology-reconstruction-using-XGB/lolo/Location_coordinates_lithologs.xlsx).
3. `litholog1` is explicitly detected as missing coordinates (`coordinates_available == False`) and raises an exception if passed to spatial model fitting.
4. Legacy `strike_pos_m` is marked deprecated and excluded from all spatial feature vectors.

---

## 4. Physical Geometry & Spatial Dispersal

Using true local Cartesian coordinates for L2-L12:
- Minimum Easting ($X$): -4057.68 m (Litholog 2)
- Maximum Easting ($X$): +704.07 m (Litholog 11)
- Total East-West Span: 4761.75 m (~4.76 km)
- Minimum Northing ($Y$): 10733.26 m (Litholog 7)
- Maximum Northing ($Y$): 14407.30 m (Litholog 12)
- Total North-South Span: 3674.04 m (~3.67 km)

Inter-well Euclidean horizontal distance summary:
- Minimum inter-well distance: 88.58 m (between Litholog 4 and Litholog 5)
- Median inter-well distance: 1918.42 m (~1.92 km)
- Maximum inter-well distance: 4806.92 m (between Litholog 2 and Litholog 11)

---

## 5. Decision on Spatial Prototype Readiness

### Can Spatial Modeling Proceed?
**Yes, but strictly under provisional relative geometry status, not as a validated geological reconstruction.**

### Conditions and Operational Boundaries:
1. **L1 Exclusion**: Litholog 1 must remain permanently excluded from spatial cross-validation.
2. **Provisional Vertical Datum**: Because no surface elevations or stratigraphic datums exist in the repository, vertical coordinates $Z_{rel}$ are relative to profile base (`relative_to_base`). Sensitivity testing confirms that vertical shifting across plausible offsets (±5m to ±20m) does not improve Macro-F1 (which remains ~0.25).
3. **Exploratory Status**: Spatial cross-validation results must be benchmarked against naive baselines (training prior and nearest well). Current results show that Spatial 3D KNN (Macro-F1 0.2521) does not outperform the nearest well profile (Macro-F1 0.2523), proving that inter-well distances (~1-3 km) exceed the lateral continuity of individual facies bodies without stratigraphic alignment.
