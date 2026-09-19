# Phase 0 Technical Note: Standardized Stratigraphic Schema & Facies Encoding

**Subsurface Stratigraphic Modeling & Active Learning Toolkit (SMALT)**  
**Author:** Mohd Fahad • IIT Kanpur  
**Document Authority:** Phase 0 Data Ingestion & Quality Control ("Read & Prove It")

---

## 1. Overview & Objective

Characterizing subsurface stratigraphic architecture requires converting heterogeneous, mixed-provenance outcrop lithologs into a mathematically rigorous, discrete state space. Phase 0 implements an automated, typed ingestion, quality control, standardization, and synthetic Gamma Ray imputation pipeline for all 11 measured sections of the Cretaceous Blackhawk Formation (Book Cliffs and Wasatch Plateau, Utah; Sahoo et al., 2016).

---

## 2. Standardized 5-Facies State Space

Field lithologies and published descriptions are aggregated into a deterministic 5-state discrete space:

| Code | Facies Identifier (`facies_name`) | Depositional Sub-Environment | Hydraulic & Flow Role | Typical Base GR Proxy |
|:---:|:---|:---|:---|:---:|
| `0` | **Coal** (`coal`) | Waterlogged peat mire / coastal plain swamp | Chronostratigraphic datum & regional seal | $20.0\text{ API}$ |
| `1` | **Channel Sandstone** (`sand`) | High-energy fluvial channel core & dunes | Primary porous reservoir fairway | $35.0\text{ API}$ |
| `2` | **Fine Sandstone / Splay** (`carbon_mud`) | Crevasse splay, levee & carbonaceous heterolithics | Intermediate connector / baffle | $130.0\text{ API}$ |
| `3` | **Siltstone** (`silt`) | Waning flow overbank & channel margin | Transition connector / mild baffle | $70.0\text{ API}$ |
| `4` | **Overbank Mudstone** (`mud`) | Low-energy floodplain & suspension settling | Regional vertical hydraulic seal | $105.0\text{ API}$ |

### 2.1 Canonical State Equivalence: Phase 0 vs. Phase 1
Phase 0 uses short lowercase technical identifiers in tabular data schemas, while Phase 1 uses descriptive sedimentological labels for Markov transition matrices and publication figures. Both represent the exact same five discrete states:

| State Code | Phase 0 Technical Identifier (`facies_name`) | Phase 1 Markov Display Name (`DEFAULT_FACIES_MAP`) | Sahoo et al. (2016) Source Analogy | Equivalence Status |
|:---:|:---|:---|:---|:---:|
| `0` | `coal` | `Coal` | Facies 6 (Coal) | **Identical** |
| `1` | `sand` | `Channel Sandstone` | Facies 1 & 2 (Trough cross-stratified & parallel-laminated sand) | **Identical** |
| `2` | `carbon_mud` | `Fine Sandstone / Splay` | Facies 3 & 5 (Splay sand & carbonaceous mud) | **Identical** |
| `3` | `silt` | `Siltstone` | Facies 4 (Coarse silt fraction) | **Identical** |
| `4` | `mud` | `Overbank Mudstone` | Facies 4 (Fine mud / shale fraction) | **Identical** |

Downstream geostatistical code in `smalt/geostat/markov.py` operates strictly on the integer `facies_code` (`0, 1, 2, 3, 4`), guaranteeing 100% mathematical and data-flow consistency between Phase 0 and Phase 1.

### 2.1 Facies Alias Resolution Dictionary
To guarantee deterministic ingestion from variant text descriptions in raw field CSVs, the loader implements an explicit normalization dictionary (`FACIES_ALIASES` in `data/loader.py`):

```python
FACIES_ALIASES: Dict[str, str] = {
    "coal": "coal",
    "sand": "sand",
    "sandstone": "sand",
    "channel sandstone": "sand",
    "channel_sandstone": "sand",
    "carbon_mud": "carbon_mud",
    "carbonaceous_mud": "carbon_mud",
    "carbonaceous_mudstone": "carbon_mud",
    "carbon mud": "carbon_mud",
    "fine sandstone / splay": "carbon_mud",
    "fine_sandstone": "carbon_mud",
    "fine sandstone": "carbon_mud",
    "splay": "carbon_mud",
    "silt": "silt",
    "siltstone": "silt",
    "mud": "mud",
    "mudstone": "mud",
    "shale": "mud",
    "overbank mudstone": "mud",
    "overbank_mudstone": "mud",
}
```

---

## 3. Standardized Output Schema

All ingested logs are standardized to the following 6-column schema (`CRITICAL_COLUMNS`):

| Column Name | Data Type | Physical Units | Description | Null Allowed? |
|:---|:---|:---|:---|:---:|
| `litholog_id` | `string` | - | Unique section identifier (e.g., `litholog1`) | No |
| `strike_pos_m` | `float64` | Meters (m) | Coordinate along depositional strike transect | No |
| `depth_m` | `float64` | Meters (m) | Stratigraphic depth from section top (1m resolution) | No |
| `facies_code` | `int64` | Integer $[0, 4]$ | Standardized 5-state integer class | No |
| `facies_name` | `string` | - | Standardized facies label | No |
| `gamma_ray` | `float64` | API Units | Measured or synthetically imputed Gamma Ray proxy | No |

---

## 4. Synthetic Gamma Ray Imputation Model

Where raw logs lack continuous wireline Gamma Ray logs, SMALT generates a physically constrained synthetic Gamma Ray proxy conforming to:

$$\text{GR}(z) = \text{base\_GR}(\text{facies}(z)) + \epsilon, \quad \epsilon \sim \mathcal{N}(0, \sigma^2=25)$$

Baseline facies parameters ($\mu \pm \sigma$):
- **Coal (`coal`)**: $20.0 \pm 5.0\text{ API}$
- **Sandstone (`sand`)**: $35.0 \pm 5.0\text{ API}$
- **Siltstone (`silt`)**: $70.0 \pm 5.0\text{ API}$
- **Mudstone (`mud`)**: $105.0 \pm 5.0\text{ API}$
- **Carbonaceous Mudstone (`carbon_mud`)**: $130.0 \pm 5.0\text{ API}$

### 4.1 Gamma Ray Scoping & Purpose
- **Synthetic Proxy Status**: The `gamma_ray` column is strictly an artificially generated petrophysical proxy. It is not, and must never be represented as, measured downhole wireline log data from the field.
- **Required Role in SMALT**: This proxy is explicitly required by the SMALT Study Plan (`docs/SMALT_Study_and_Build_Plan.md`, Track C, Phase 4) as a feature input for the baseline tabular spatial classifier `ml/xgboost_model.py` ($[X, Z, \text{GR}, \text{prev\_facies}]$) and margin-sampling active learning loops.
- **No Scientific Intrusion**: Gamma Ray values are not used in Phase 0 quality gates or Phase 1 Markov transition probability matrices, which operate strictly on lithofacies classifications.

---

## 5. Quality Control & Validation Rules

The `LithologLoader` enforces strict quality control gates prior to export:

1. **Non-Negative Depths**: Asserts $\text{Top} \ge 0$ and $\text{Bottom} \ge 0$. Any negative depth raises `ValueError`.
2. **Strict Layer Thickness**: Asserts $\text{Bottom} - \text{Top} > 0$. Zero-thickness or inverted intervals raise `ValueError`.
3. **Monotonic Depth Ordering**: Asserts $\frac{d(\text{depth\_m})}{di} > 0$ strictly for each individual litholog.
4. **Interval Discretization & Overlap Resolution**: Raw $(Top, Bottom)$ intervals are discretized into 1-meter integer steps ($d \in [\lfloor Top \rfloor, \lceil Bottom \rceil)$). Duplicate depth entries are deduplicated deterministically with `keep='last'`.
5. **Round-Trip Parquet/CSV Integrity**: Exported datasets must survive read/write round-trips with zero null values and identical column schemas.

### 5.1 Litholog 9 Overlap Resolution Audit
In `data/raw_lithologs/litholog9.csv`, three consecutive raw intervals exhibit 1-meter boundary overlaps:
- Raw row 4: `25, 29, sand` (produces depths 25, 26, 27, 28)
- Raw row 5: `28, 30, silt` (produces depths 28, 29)
- Raw row 6: `29, 33, mud` (produces depths 29, 30, 31, 32)

Deduplication behavior under `keep='last'`:
- **Depth 28 m**: Appears in row 4 (`sand`) and row 5 (`silt`). The later stratigraphic unit (`silt`) is retained.
- **Depth 29 m**: Appears in row 5 (`silt`) and row 6 (`mud`). The later stratigraphic unit (`mud`) is retained.
- **Resulting Vertical Profile**: $27\text{m}=\text{sand}$, $28\text{m}=\text{silt}$, $29\text{m}=\text{mud}$, $30\text{m}=\text{mud}$.
- **Sedimentological Consistency**: This preserves the transitional siltstone bed between underlying mudstone and overlying channel sand in upward stratigraphic sequence (decreasing depth: $30\text{m mud} \to 29\text{m mud} \to 28\text{m silt} \to 27\text{m sand}$).
- **Raw Data Preservation**: The raw CSV file remains completely untouched. The deduplication is executed purely in-memory during standardized DataFrame generation.

---

## 6. Dataset Provenance & Benchmark Scope

The working dataset incorporates mixed provenance:
- **Litholog 11**: Independently validated QC benchmark section (73/78 m match, 93.59% structural agreement against manual digitization).
- **Lithologs 1 & 9**: Source-derived from published figures; independent primary validation undocumented / unknown.
- **Lithologs 2-8 & 10**: AI-assisted reconstructions from Sahoo et al. (2016) outcrop photomosaics; unvalidated.

> [!IMPORTANT]
> The 93.59% accuracy benchmark applies strictly to Litholog 11. It must not be extrapolated to Lithologs 1-10. All 11 logs are preserved as the active working dataset with provenance tracked in `data/provenance_manifest.json`.
