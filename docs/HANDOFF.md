# SMALT Project Handoff & Phase Registry

## Phase Status Overview
- **Phase 0 (Data Ingestion & Quality Control)**: VERIFIED (11/11 lithologs loaded cleanly, 920 observations, 10/10 unit tests passing, provenance manifest registered)
- **Phase 1 (1D Vertical Markov Succession Analysis)**: STALE - RECOMPUTATION REQUIRED (Previously verified on 3-well subset [247 observations]; must be recomputed on complete 11-well working dataset [920 observations])
- **Phase 2 (2D Cross-Sectional Geostatistical Modeling)**: READY FOR SPRINT (Validated solely against Sahoo et al. [2016] geological priors [W/T = 35, mean thickness 5.8m, NTG 17%-46%] and synthetic realizations; independent of litholog subset)

---

## 1. Current Working Dataset State & Ingestion Authority

- **Unified Ingestion Source**: `data/processed/lithologs_unified.parquet` and `data/processed/lithologs_unified.csv`
- **Total Ingested Lithologs**: 11 (`litholog1.csv` through `litholog11.csv`)
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

## 3. Stale Phase 1 Empirical Outputs (Old 3-Well Subset)

The following empirical outputs in the repository were generated using only 3 wells (`litholog1`, `litholog9`, `litholog11`; 247 observations) and are currently **STALE**:

1. **`results/phase1_markov_summary.json`**:
   - `metadata.n_wells = 3`, `wells_used = ["litholog1", "litholog11", "litholog9"]`, `n_observations = 247`.
   - Regular and embedded transition matrices ($P_{\text{reg}}, P_{\text{emb}}$) and count matrices ($N_{\text{reg}}, N_{\text{emb}}$).
2. **Diagnostic Visualizations (`docs/figures/`)**:
   - `phase1_transition_matrix_regular.png`
   - `phase1_transition_matrix_embedded.png`
   - `phase1_stationary_vs_empirical.png`
   - `phase1_facies_succession_network.png`
3. **Documentation Metrics (`docs/notes/phase1_markov.md`)**:
   - Empirical counts and stationary distributions derived from the 247-observation subset.

---

## 4. Phase 2 Status Assessment

- **Dependency Analysis**: Phase 2 ("2D Object-Based Fluvial Generator", Track B, Days 21-30) is formulated as an unconditioned stochastic body generator conditioned strictly on literature priors from Sahoo et al. (2016):
  - Channel aspect ratio $W/T = 35$
  - Mean channel thickness $\sim 5.8\text{ m}$
  - Crevasse splay width range ($10 - 130\text{ m}$)
  - Target Net-to-Gross envelope ($17\% - 46\%$)
- **Validation Gate ("Done When")**: 100 unconditioned synthetic realizations yielding mean $W/T \in 35 \pm 2$ and $\text{NTG} \in [17\%, 46\%]$.
- **Conclusion**: Phase 2 does **not** depend quantitatively on the empirical litholog subset. No components of Phase 2 need rebuilding or invalidation due to the dataset update. It remains **READY FOR SPRINT**.

---

## 5. Recommended Next Action

1. **Recompute Phase 1 Pipeline**: Execute `python scripts/run_phase1_markov.py` against the unified 11-litholog dataset (`data/processed/lithologs_unified.parquet`, 920 rows) to regenerate fresh transition matrices, stationary diagnostics, publication figures, and JSON summaries.
2. **Sprint Phase 2**: Implement `smalt/geostat/object_sim.py` according to Sahoo et al. (2016) geometrical priors.
