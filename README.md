# SMALT: Subsurface Stratigraphic Modeling & Active Learning Toolkit

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-126%20passed%20%7C%20100%25-brightgreen.svg)](tests/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**SMALT** is an open, verifiable, and mathematically grounded Python research toolkit for 1D vertical stratigraphic Markov analysis, continuous transition-probability geostatistics, 2D process-based fluvial forward modeling, and sparse-borehole reconstruction benchmarking. The framework is calibrated against fluvio-deltaic outcrop analogs from the Upper Cretaceous **Blackhawk Formation** (Wasatch Plateau, Utah) to bridge classical sedimentological principles with probabilistic subsurface characterization.

---

## Authoritative Documentation & Registers

- **Scientific Methods Register:** [`docs/SCIENTIFIC_METHODS_REGISTER.md`](docs/SCIENTIFIC_METHODS_REGISTER.md) - Permanent identifiers, mathematical derivations, parameter provenance, and oral defense arguments for every substantive method.
- **Figure Provenance Register:** [`docs/FIGURE_PROVENANCE_REGISTER.md`](docs/FIGURE_PROVENANCE_REGISTER.md) - Complete inventory and provenance audit of all scientific visualizations.
- **Professor Meeting Notes:** [`docs/PROFESSOR_MEETING_NOTES.md`](docs/PROFESSOR_MEETING_NOTES.md) - Question-driven oral defense notes for student-supervisor meetings.
- **Scientific Audit Summary:** [`docs/SCIENTIFIC_AUDIT_SUMMARY.md`](docs/SCIENTIFIC_AUDIT_SUMMARY.md) - Executive audit summary, resolved discrepancies, and remaining uncertainties.
- **BibTeX References Database:** [`docs/REFERENCES.bib`](docs/REFERENCES.bib) - Complete, verified primary literature citations.
- **Sprint K Stage 3 Report:** [`sprints/audit_sprint_k/SPRINT_K_STAGE_3_VALIDATION_BENCHMARK_REPORT.md`](sprints/audit_sprint_k/SPRINT_K_STAGE_3_VALIDATION_BENCHMARK_REPORT.md) - Multi-seed geological validation and sparse-borehole benchmark report.
- **Handoff Registry:** [`HANDOFF.md`](HANDOFF.md) - Full sprint chronology, numerical condition numbers, and test acceptance criteria.

---

## Key Scientific Modules

1. **Standardized Stratigraphic Data Pipeline (Phase 0 / Sprint E):**
   - Automated quality control, depth monotonicity enforcement ($z_{t+1} > z_t$), and provenance validation across 12 manually digitized lithologs (1033.0 m, 296 beds).
   - Provenance manifest strictly registers L12 core depth restriction (0-111 m) and Litholog 1 spatial isolation ([`data/provenance_manifest.json`](data/provenance_manifest.json)).
2. **Canonical Six-State Facies Schema (Sprint F):**
   - Retains all six primary facies defined in Sahoo et al. (2016) Table 1: `sand` (0), `p_sand` (1), `ripples` (2), `carbon_mud` (3), `coal` (4), `mud` (5).
3. **Common-Zero Datum Alignment (Sprint H):**
   - Normalizes measured outcrop and core sections to the base of the regional Axel Anderson coal zone as a chronostratigraphic datum line.
4. **1D Stratigraphic Markov Engine (Phase 1):**
   - Evaluates vertical upward transitions in depositional order; supports both regular (1m discretized) and embedded (boundary-crossing, $P_{ii}=0$) chains.
   - Computes invariant stationary distributions ($\boldsymbol{\pi}\mathbf{P} = \boldsymbol{\pi}$), transition asymmetry matrices ($\mathbf{A} = \mathbf{P} - \mathbf{P}^T$), and directional entropy (Hattori, 1976).
5. **Decoupled Hybrid Markov-Monte Carlo Architecture:**
   - Decouples transition sequencing (embedded Markov chain) from continuous thickness sampling (LogNormal MLE: $\mu=1.7870, \sigma=0.8823$), eliminating the "sticky state" trapping pathology (Dindarloo et al., 2015).
6. **2D Fluvial Forward Simulator (Sprint K):**
   - Unconditioned 2D process-based forward simulation engine ($X \times Z = 2000\text{ m} \times 50\text{ m}$) generating reproducible stratigraphy with known ground truth everywhere:
     - Channel belt ribbons with parabolic scours and $W/T \approx 35 \pm 5$ aspect ratios (Sahoo et al., 2016).
     - Internal 20% basal planar sandstone partition (`p_sand`).
     - Flanking crevasse splays ($W \in [20, 100]\text{ m}$).
     - Regional continuous coal markers ($1.0\text{ m}$ thickness).
     - Deterministic erosional superposition and 100% bitwise seed reproducibility.
7. **Sparse-Borehole Reconstruction Benchmark:**
   - Evaluates reconstruction baselines across borehole spacings ($400\text{ m}, 600\text{ m}, 1000\text{ m}$) with strict leakage safeguards (borehole columns strictly masked out during evaluation).

---

## Canonical Six-State Facies Schema

| Code | Facies Name | Description | Color Hex | Color Swatch |
|:---:|:---|:---|:---:|:---:|
| `0` | **sand** | Channel Sandstone (Trough cross-stratified medium-coarse sand) | `#E69F00` | Gold |
| `1` | **p_sand** | Planar Sandstone (Parallel-laminated basal/sheet sand) | `#F0E442` | Yellow |
| `2` | **ripples** | Rippled Heterolithics (Interbedded rippled sand/mud, splay lobes) | `#56B4E9` | Light Blue |
| `3` | **carbon_mud** | Carbonaceous Mudstone (Organic-rich floodplain margins) | `#009E73` | Green |
| `4` | **coal** | Coal (Waterlogged peat mire / swamp) | `#000000` | Black |
| `5` | **mud** | Overbank Mudstone (Floodplain suspension mud and shale) | `#999999` | Grey |

---

## Quick Start & Verification

### 1. Execute Automated Test Suite

Execute the targeted Sprint K verification suite:
```bash
python -m pytest tests/test_phase0_loader.py tests/test_phase1_markov.py tests/test_sprint_k.py -v
```
*(Outcome: 36 passed in ~32.5s, 100% pass rate).*

Execute the complete repository test suite:
```bash
python -m pytest tests/ -q
```
*(Outcome: 126 passed in ~8m, 100% pass rate across 11 test modules).*

### 2. Run Forward Simulation Validation Campaign

Execute the multi-seed validation campaign and sensitivity analysis:
```bash
python -m smalt.forward.validation
```

### 3. Run Sparse-Borehole Reconstruction Benchmark

Run the 120-evaluation benchmark across borehole spacings:
```bash
python -m smalt.forward.benchmark
```

---

## Primary Literature References

1. **Sahoo, H., Gani, M. R., Hampson, G. J., Gani, N. D., & Ranson, A. (2016).** Facies- to sandbody-scale heterogeneity in a tight-gas fluvial reservoir analog: Blackhawk Formation, Wasatch Plateau, Utah, USA. *Marine and Petroleum Geology*, 78, 48-69. DOI: 10.1016/j.marpetgeo.2016.02.005.
2. **Abasaghi, F., Omidpour, A., & Khan, Z. A. (2026).** Integrating Markov Chain Analysis, Cyclostratigraphy and Sequence Stratigraphy to Decipher the Oligo-Miocene Asmari Formation's Evolution in the Zagros Basin. *Basin Research*, 38(1), e70139. DOI: 10.1111/bre.70139.
3. **Weissmann, G. S., & Fogg, G. E. (1999).** Multi-scale alluvial fan heterogeneity modeled with transition probability geostatistics in a sequence stratigraphic framework. *Journal of Hydrology*, 226(1-2), 48-65. DOI: 10.1016/S0022-1694(99)00160-2.
4. **Dindarloo, S. R., Bagherieh, A., Hower, J. C., Calder, J. H., & Wagner, N. J. (2015).** Coal modeling using Markov Chain and Monte Carlo simulation: Analysis of microlithotype and lithotype succession. *Sedimentary Geology*, 329, 1-11. DOI: 10.1016/j.sedgeo.2015.08.005.
5. **He, J., La Croix, A. D., Wang, J., Ding, W., & Underschultz, J. R. (2019).** Using neural networks and the Markov Chain approach for facies analysis and prediction from well logs in the Precipice Sandstone and Evergreen Formation, Surat Basin, Australia. *Marine and Petroleum Geology*, 101, 410-427. DOI: 10.1016/j.marpetgeo.2018.12.022.
6. **Carle, S. F., & Fogg, G. E. (1996).** Transition probability-based indicator geostatistics. *Mathematical Geology*, 28(4), 453-476.
7. **Carle, S. F., & Fogg, G. E. (1997).** Modeling spatial variability with one- and multidimensional continuous-lag Markov chains. *Mathematical Geology*, 29(7), 891-918.
8. **Powers, M. C., & Easterling, R. G. (1982).** Improved methodology for using embedded Markov chains to describe cyclicity in sedimentary successions. *Journal of Sedimentary Petrology*, 52(3), 913-923.
9. **Hattori, I. (1976).** Entropy in Markov chains and discrimination of cyclic patterns in lithologic successions. *Mathematical Geology*, 8(4), 477-497.
