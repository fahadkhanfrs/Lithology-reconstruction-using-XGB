# SMALT Sprint K - Supervisory Review Gate Report

**Document ID:** `SMALT-GATE-SPRINT-K-001`  
**Date:** October 2026  
**Reviewer:** Principal Geostatistical Scientist, Sedimentologist, and Research-Software Auditor  
**Repository Branch:** `sprint-k`  
**Baseline Commit:** `f08587ca107b518c8baac560339381c53af9ed1f`  
**Gate Status:** Review-Only Checkpoint Completed; Discrepancies Registered; Ready for Supervisory Review  
**Discrepancy Register:** [`sprints/audit_sprint_k/stage4_supervisory_discrepancies.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/stage4_supervisory_discrepancies.csv)  

---

## 1. Repository State and Preservation Confirmation

### 1.1 Branch and Baseline Status
- **Current Branch:** `sprint-k`
- **Current HEAD Commit:** `f08587ca107b518c8baac560339381c53af9ed1f` (Sprint J baseline: *Sprint J: Empirical horizontal transition calibration and lateral continuity audit*).
- **Baseline Confirmation:** `f08587ca107b518c8baac560339381c53af9ed1f` remains the verified, intended baseline for Sprint K.
- **Commit Status:** Zero commits have been made on `sprint-k`. All Stage 1, 2, 3, and 4 changes remain uncommitted in the working tree. No changes have been pushed, reset, cleaned, or discarded.

### 1.2 Full Working-Tree Status
The working tree contains user modifications and newly generated Sprint K audit packages:

1. **Tracked Modified Files (Preserved):**
   - [`HANDOFF.md`](file:///d:/Lithology-reconstruction-using-XGB/HANDOFF.md)
   - [`README.md`](file:///d:/Lithology-reconstruction-using-XGB/README.md)
   - [`data/provenance_manifest.json`](file:///d:/Lithology-reconstruction-using-XGB/data/provenance_manifest.json) (corrected provenance: all 12 logs manually digitized; L12 coverage 0-111 m)
   - [`docs/HANDOFF.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/HANDOFF.md)
   - [`docs/SMALT_Study_and_Build_Plan.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/SMALT_Study_and_Build_Plan.md)
   - [`smalt/descriptive/analyzer.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/descriptive/analyzer.py) (metadata descriptions)

2. **Untracked Modules and Directories (Preserved):**
   - [`docs/FIGURE_PROVENANCE_REGISTER.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/FIGURE_PROVENANCE_REGISTER.md)
   - [`docs/PROFESSOR_MEETING_NOTES.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/PROFESSOR_MEETING_NOTES.md)
   - [`docs/REFERENCES.bib`](file:///d:/Lithology-reconstruction-using-XGB/docs/REFERENCES.bib)
   - [`docs/SCIENTIFIC_AUDIT_SUMMARY.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/SCIENTIFIC_AUDIT_SUMMARY.md)
   - [`docs/SCIENTIFIC_METHODS_REGISTER.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/SCIENTIFIC_METHODS_REGISTER.md)
   - [`smalt/descriptive/group_analysis.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/descriptive/group_analysis.py)
   - [`smalt/forward/`](file:///d:/Lithology-reconstruction-using-XGB/smalt/forward/) (`__init__.py`, `config.py`, `simulator.py`, `sampling.py`, `validation.py`, `benchmark.py`)
   - [`sprints/audit_sprint_k/`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/) (all 26 deliverables: CSVs, figures, manifests, reports)
   - [`tests/test_sprint_k.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_k.py)

### 1.3 Test Verification and Fresh Execution Log
To distinguish fresh runs from historical report copies, a fresh test execution was performed during this review pass:

- **Command Executed:**
  ```powershell
  python -m pytest tests/test_sprint_k.py
  ```
- **Execution Timestamp:** 2026-10-11T07:18:06Z
- **Actual Fresh Result:** **19 passed in 41.26s (100% pass rate, zero failures)**.
- **Coverage:** Canonical 6-facies schema, grid shape and bottom-to-top orientation, seed reproducibility, object geometry, run lengths, borehole sampling, ground truth isolation, Net-to-Gross hierarchy, configuration validation, regression, multi-seed validation, baseline models, spacing decay, leakage prevention, L1/L12 boundary invariants, coordinate quarantine, GR leakage prohibition, and Lilliefors goodness-of-fit.

Targeted regression suite:
- **Command Executed:**
  ```powershell
  python -m pytest tests/test_phase0_loader.py tests/test_phase1_markov.py tests/test_sprint_k.py
  ```
- **Actual Fresh Result:** **41 passed in 38.52s (100% pass rate, zero failures)**.

Full repository suite:
- **Command:** `python -m pytest tests/`
- **Historical Result:** 112 passed, 4 warnings in 683.04s.

---

## 2. Reconciled Discrepancy Findings and Evidence

The discrepancy register is recorded in [`sprints/audit_sprint_k/stage4_supervisory_discrepancies.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/stage4_supervisory_discrepancies.csv). Below is the comprehensive evidence and analysis for each issue.

---

### Issue A: Common Stratigraphic Datum (Highest Priority)

#### 1. Report Claim
In the Stage 4 report, Claims Register (`CLAIM-DATUM-001`), and Methods Register (`METHOD-DATUM-001`):
> *"All measured outcrop lithologs are vertically shifted so that the base of the continuous, regional Axel Anderson coal zone (Sanchez & Brown, 1986; Sahoo et al., 2016) serves as the common zero stratigraphic datum ($Z = 0.0\text{ m}$). Formula: $Z_{\text{strat}} = z_{\text{measured}} - z_{\text{base\_coal}}$."*

#### 2. Repository Code & Field Evidence
Tracing this claim through [`smalt/spatial/datum.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/spatial/datum.py) and primary field data reveals three facts:
1. **Code Behavior:** In `datum.py` (lines 209-211 and lines 257-262), the function `align_to_common_datum()` does **not** search for a coal seam in each litholog, nor does it subtract `z_base_coal`. It applies:
   $$z_{\text{strat}} = H_{\text{max}} - d$$
   where $H_{\text{max}}$ is the total height of the measured section and $d$ is depth measured from the top. This transformation places $z=0$ at the **base of the measured section**, directly implementing Prof. Hiranya Sahoo's provisional instruction: *"All litholog zeros should be treated as being at the same reference level."* The formula $Z_{\text{strat}} = z_{\text{measured}} - z_{\text{base\_coal}}$ does not exist anywhere in the code.
2. **Field Geological Evidence:** Sahoo et al. (2016, Section 3, Section 4.3.1, Figure 6) document that the coal seam at the base of Cliff Face 1 (L2-L8) is the Axel Anderson coal zone. However, across Cliff Faces 7-8 (L9, L11) and subsurface drill core EM-137C, an continuous basal Axel Anderson coal is either absent, unverified, or at a different stratigraphic interval (Sprint A audit, line 269). Treating all sections as aligned to a single continuous chronostratigraphic coal base is an **inferred geological correlation assumption**, not an independently surveyed elevation datum.
3. **Forward Simulator Independence:** The Sprint K 2D forward simulator is an **unconditioned process-based forward generator**. It creates a synthetic domain of $X=2000\text{ m}, Z=50\text{ m}$ from scratch and does **not** ingest transformed field depths or coordinates as predictive inputs. Field data were used solely to compute 1D summary statistics (facies proportions, bed thickness MLE).

#### 3. Audit Verdict
**Status: Contradicted / Overstated Report Claim.**  
The Stage 4 report conflated:
- Prof. Sahoo's provisional section-base alignment rule ($z=0$ at source base);
- The local geological observation that the Axel Anderson coal lies at the base of Cliff Face 1; and
- A non-existent formula ($z - z_{\text{base\_coal}}$).

#### 4. Recommended Resolution
Downgrade the claim in reports and registers from *"verified common chronostratigraphic coal datum"* to:
> *"Provisional section-base datum assumption instructed by faculty supervisor, supported in Cliff Face 1 by the basal Axel Anderson coal zone."*

Maintain strict semantic separation between:
- **Measured Depth ($d$):** Downward distance from wellhead or cliff top.
- **Stratigraphic Height ($z_{\text{strat}}$):** Upward distance from measured section base ($H_{\text{max}} - d$).
- **Synthetic Model Grid Elevation ($Z$):** Height coordinate in the synthetic forward simulator ($0$ to $50\text{ m}$).

---

### Issue B: Lognormal Goodness-of-Fit and Lilliefors Verification

#### 1. Report Claim
Stage 4 reported that the log-normal channel sandstone thickness distribution ($\mu = 1.7870, \sigma = 0.8823$) was evaluated with an ordinary KS test ($p = 0.8023$) and a parameter-estimation-aware Lilliefors bootstrap test yielding $p = \mathbf{0.4145}$ (and $p = \mathbf{0.6940}$ with 1m discretization).

#### 2. Independent Reproduction Results
Executing 5,000 Monte Carlo bootstrap replicates on the 34 Upstream channel sandstone beds reproduced the statistics:
- **Sample Size:** $N = 34$ channel sandstone beds from 8 Upstream lithologs.
- **MLE Parameters:** $\hat{\mu} = 1.7870, \hat{\sigma} = 0.8823$.
- **Empirical KS Statistic:** $D_{\text{obs}} = 0.1059$.
- **Ordinary KS Test:** $p = 0.8023$ (inflated due to parameter estimation on the same data).
- **Lilliefors Parametric Bootstrap ($B=5000$, unclipped):** $p = \mathbf{0.4196 \pm 0.0070}$.
- **Discretization-Aware Lilliefors Bootstrap ($B=5000$, 1.0 m integer rounding):** $p = \mathbf{0.6930 \pm 0.0065}$.
- **Clipped $[1.5, 25.0]\text{ m}$ Lilliefors Bootstrap ($B=5000$):** $p = \mathbf{0.4430 \pm 0.0070}$.
- **Clipped $[1.5, 25.0]\text{ m}$ + 1.0 m Discretized ($B=5000$):** $p = \mathbf{0.8476 \pm 0.0051}$.

#### 3. Audit Verdict
**Status: Partially Verified with Necessary Scientific Caveats.**
- The term "Lilliefors parametric bootstrap" is mathematically justified: parameters $\mu^*$ and $\sigma^*$ were refitted on every simulated replicate, exactly following the Lilliefors (1967) procedure for continuous distributions with unknown parameters.
- **Critical Caution:** Non-rejection ($p \approx 0.42 - 0.69 > 0.05$) indicates that the empirical sample of 34 beds is **statistically consistent with** a log-normal distribution, but **does not prove** that the true depositional process is log-normal. Many right-skewed heavy-tailed distributions (e.g. Gamma, Weibull, or truncated Pareto) would also fail to be rejected on a sample of $N=34$ beds.
- Truncation at $[1.5, 25.0]\text{ m}$ is an engineering configuration in `ChannelBeltConfig` to prevent unphysical bed thicknesses in a 50 m column, whereas empirical beds range from $1.0\text{ m}$ to $32.0\text{ m}$.

---

### Issue C: Bed-Thickness Parameter Uncertainty and Clustered Sampling

#### 1. Report Estimates Verified
- Point MLE: $\hat{\mu} = 1.7870$, $\hat{\sigma} = 0.8823$ (natural logarithm of bed thickness in meters).
- Implied Theoretical Median: $\exp(\hat{\mu}) = \exp(1.7870) = \mathbf{5.97\text{ m}}$.
- Empirical Sample Median: $\mathbf{5.50\text{ m}}$.
- Empirical Sample Mean: $\mathbf{8.50\text{ m}}$.
- Empirical Sample Range: $1.0\text{ m}$ to $32.0\text{ m}$.

#### 2. Cluster Bootstrap across Litholog Clusters
Because the 34 beds originate from only 8 distinct lithologs (Litholog 2: 5 beds; L3: 3; L4: 3; L5: 5; L6: 5; L7: 6; L8: 3; L10: 4), beds within the same section are stratigraphically clustered. Resampling the 8 litholog clusters with replacement ($B=5000$):
- **Cluster Bootstrap $\mu$:** Mean $= 1.7984$, $95\%\text{ CI} = [\mathbf{1.4975}, \mathbf{2.1900}]$.
- **Cluster Bootstrap $\sigma$:** Mean $= 0.8621$, $95\%\text{ CI} = [\mathbf{0.7253}, \mathbf{0.9958}]$.
- **Implied Median Bed Thickness ($e^\mu$):** $95\%\text{ CI} = [\mathbf{4.47\text{ m}}, \mathbf{8.94\text{ m}}]$.

#### 3. Audit Verdict
**Status: Verified as a Candidate Engineering Model.**  
The true effective sample size is approximately 8, not 34. The median channel sandstone bed thickness can plausibly range anywhere between $4.5\text{ m}$ and $8.9\text{ m}$ across the basin. The log-normal distribution must remain classified strictly as a *calibrated candidate engineering parameterization*, never as a proven universal geological law.

---

### Issue D: Net-to-Gross Target Semantics and Stopping Rule

#### 1. Code Tracing and Target Provenance
- In [`smalt/forward/config.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/forward/config.py): `target_net_to_gross = 0.45`, `tolerance_net_to_gross = 0.08`.
- In Sahoo et al. (2016, page 61), net-to-gross across pseudowells ranges from $17\%$ to $46\%$, with upstream cliff faces averaging around $40-46\%$.
- In SMALT's empirical Upstream outcrop logs (L2-L8, L10):
  - Single-facies `sand` (Facies 0) = $289.0\text{ m} / 673.0\text{ m} = \mathbf{42.94\%}$.
  - Basal `p_sand` (Facies 1) = $51.0\text{ m} / 673.0\text{ m} = \mathbf{7.58\%}$.
  - Pure channel-belt sand (`sand` + `p_sand`) = $340.0\text{ m} / 673.0\text{ m} = \mathbf{50.52\%}$.
  - Total coarse siliciclastics (`sand` + `p_sand` + `ripples`) = $373.0\text{ m} / 673.0\text{ m} = \mathbf{55.42\%}$.
- The value $0.45$ ($45.0\%$) was selected as a **calibrated midpoint prior** between single-facies sand ($42.9\%$) and pure channel-belt sand ($50.5\%$), falling within Sahoo et al.'s reported 17-46% range.

#### 2. Stopping Rule Mechanics
In [`smalt/forward/simulator.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/forward/simulator.py#L280-L286):
```python
current_sand_cells = np.sum(grid == channel_cfg.core_facies_code)  # Facies 0 only
current_ntg = current_sand_cells / (nx * nz)
if placed_channels >= 3 and current_ntg >= (target_ntg - self.config.tolerance_net_to_gross):
    break
```
- The stopping rule checks `current_ntg >= 0.45 - 0.08 = 0.37` on Facies 0 (`sand`) alone.
- Placement stops as soon as `sand` reaches $37.0\%$, realizing **$38.09\% \pm 2.28\%$**.
- Basal planar sand simultaneously accumulates to **$9.46\% \pm 1.20\%$**.
- Total pure sand achieves **$47.55\% \pm 3.13\%$**, matching the empirical field target ($50.52\%$) within $2.97\%$ error.

#### 3. Audit Verdict
**Status: Partially Verified / Semantic Clarification Required.**  
There is an explicit semantic ambiguity: the configuration parameter is named `target_net_to_gross = 0.45`, but the code evaluates it against single-facies `sand`.
- If the target was intended for single-facies channel sandstone ($42.9\%$), the tolerance window ($0.37$) terminates channel placement prematurely, causing the $38\%$ realization.
- If the target was intended for total channel-belt sand ($50.5\%$), evaluating `grid == 0` alone undercounts net sand during the stopping check, but the realized pure sand ($47.55\%$) honors the empirical target.
- **Decision:** Do not change code during this review gate. Clarify in documentation that `target_net_to_gross` in the current simulator acts as a lower-bound threshold on single-facies sand.

---

### Issue E: Downstream Ripple Discrepancy and False 23m Bed Claim

#### 1. Report Claim
Stage 3 (line 73) and Stage 4 (`CLAIM-HET-001`) stated:
> *"The downstream outcrop preset exhibits a ripple deficit (5.02% simulated vs 21.15% empirical) because Litholog 9 contains an unusually thick 23.0 m continuous ripple-laminated heterolith unit."*

#### 2. Repository Code & Field Evidence
Filtering the raw bed intervals in [`sprints/audit_sprint_k/litholog_group_raw_beds.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/litholog_group_raw_beds.csv) and [`data/raw_lithologs/litholog9.csv`](file:///d:/Lithology-reconstruction-using-XGB/data/raw_lithologs/litholog9.csv):
- **Litholog 9:** Total thickness = $77.0\text{ m}$. Contains 4 ripple beds with thicknesses **6.0 m, 3.0 m, 2.0 m, and 4.0 m**. Total ripple thickness = **$15.0\text{ m}$**.
- **Litholog 11:** Total thickness = $79.0\text{ m}$. Contains 4 ripple beds with thicknesses **3.0 m, 5.0 m, 7.0 m, and 3.0 m**. Total ripple thickness = **$18.0\text{ m}$**.
- **Group Total:** $15.0\text{ m} + 18.0\text{ m} = 33.0\text{ m}$ out of $156.0\text{ m} = \mathbf{21.15\%}$.
- **Finding:** **There is NO 23.0 m continuous ripple bed anywhere in Litholog 9 or Litholog 11.** The maximum ripple bed in L9 is $6.0\text{ m}$, and in L11 is $7.0\text{ m}$.

#### 3. Root Cause of Simulator Deficit
The simulator's realized ripple proportion ($5.02\% \pm 0.99\%$) is constrained by **structural limitations of the forward simulator's crevasse splay model**:
1. Splays are placed strictly as lateral wings flanking channel margins.
2. In a 50 m column, only ~25-30 channels are placed.
3. Splays have configured thickness bounds of $1.0 - 4.5\text{ m}$ and width bounds of $40 - 130\text{ m}$.
4. Many flanking splays are erosively truncated by subsequent channel placements.
5. Consequently, splays can never physically occupy $21\%$ of the 2D cross-section under the current architecture.

#### 4. Audit Verdict
**Status: Contradicted Report Assertion.**  
The claim that Litholog 9 contains a 23.0 m continuous ripple unit is **false** and must be retracted. The downstream ripple mismatch is a genuine forward-modeling limitation reflecting both small-sample variance ($N=2$ outcrop logs) and structural constraints of the channel-attached splay generator.

---

### Issue F: Parameter Provenance Consistency Across Documents

Auditing [`parameter_provenance.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/parameter_provenance.csv) against [`scientific_claims_register.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/scientific_claims_register.csv) and [`docs/SCIENTIFIC_METHODS_REGISTER.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/SCIENTIFIC_METHODS_REGISTER.md) identified three provenance classifications requiring refinement:

1. **Basal Planar Sand Fraction (`basal_fraction = 0.20`):**
   - Empirical Upstream data show `p_sand` represents $15.0\%$ of total channel-belt sand ($51.0\text{ m} / 340.0\text{ m}$).
   - Setting `basal_fraction = 0.20` is a rounded **calibrated engineering choice**, not a directly measured physical invariant.
2. **Coal Horizons (`num_horizons = 2`, `thickness_m = 1.0`):**
   - Empirical median and mode of coal beds are exactly $1.00\text{ m}$.
   - Upstream logs have 12 coal beds across 8 logs (1.5 beds per 84m log).
   - Setting 2 continuous sheet horizons placed at relative elevations $[0.25, 0.75]$ is a **modeling assumption** designed to match total coal volume (~1.7% simulated vs 1.9% field).
3. **Channel Aspect Ratio Prior ($W/T = 35 \pm 5$):**
   - Correctly derived from Sahoo et al. (2016, p. 52).
   - This is a **regional literature prior**, not an empirically identifiable parameter from SMALT's 1D lithologs.

---

## 3. Prioritized Required Corrections

The required corrections are ranked by priority. All proposed modifications are documented below and **left unapplied** until supervisory review:

| Priority | Issue ID | Required Correction | Target Files |
| :---: | :--- | :--- | :--- |
| **P1** | `DISC-DATUM-001` | Correct the Axel Anderson coal datum claim. State that sections are aligned at the base ($z=0$) per Prof. Sahoo's provisional rule, not via $Z = z - z_{\text{base\_coal}}$. Clarify that the synthetic simulator is unconditioned and does not consume field depths. | [`sprints/audit_sprint_k/SPRINT_K_STAGE_4_SCIENTIFIC_AUDIT_REPORT.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/SPRINT_K_STAGE_4_SCIENTIFIC_AUDIT_REPORT.md)<br>[`sprints/audit_sprint_k/scientific_claims_register.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/scientific_claims_register.csv)<br>[`docs/SCIENTIFIC_METHODS_REGISTER.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/SCIENTIFIC_METHODS_REGISTER.md) |
| **P2** | `DISC-RIPPLE-001` | Retract the false claim of a 23.0 m continuous ripple unit in Litholog 9. Report true empirical bed thicknesses (max 6.0m in L9, max 7.0m in L11). Acknowledge that the downstream preset ripple deficit arises from structural forward simulator constraints (channel-attached splays) and small-sample variance. | [`sprints/audit_sprint_k/SPRINT_K_STAGE_4_SCIENTIFIC_AUDIT_REPORT.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/SPRINT_K_STAGE_4_SCIENTIFIC_AUDIT_REPORT.md)<br>[`sprints/audit_sprint_k/scientific_claims_register.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/scientific_claims_register.csv)<br>[`sprints/audit_sprint_k/SPRINT_K_STAGE_3_VALIDATION_BENCHMARK_REPORT.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/SPRINT_K_STAGE_3_VALIDATION_BENCHMARK_REPORT.md) |
| **P3** | `DISC-THICK-001` | Reclassify log-normal bed thickness from "validated physical law" to "candidate engineering parameterization that cannot be rejected". Report cluster-bootstrap parameter uncertainty ($\mu \in [1.50, 2.19]$, median $e^\mu \in [4.47, 8.94]\text{ m}$). | [`sprints/audit_sprint_k/SPRINT_K_STAGE_4_SCIENTIFIC_AUDIT_REPORT.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/SPRINT_K_STAGE_4_SCIENTIFIC_AUDIT_REPORT.md)<br>[`sprints/audit_sprint_k/scientific_claims_register.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/scientific_claims_register.csv) |
| **P4** | `DISC-NTG-001` | Document the semantic ambiguity of `target_net_to_gross = 0.45`. Clarify that the stopping rule evaluates single-facies `sand` at $37\%$ threshold, while realized pure sand ($47.55\%$) honors total channel-belt sand ($50.52\%$). | [`sprints/audit_sprint_k/metric_definitions.md`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/metric_definitions.md)<br>[`sprints/audit_sprint_k/parameter_provenance.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/parameter_provenance.csv) |
| **P5** | `DISC-PARAM-001` | Refine parameter provenance descriptions: distinguish empirical measurements (coal median 1.0m) from calibrated choices (20% basal partition, 2 coal horizons) and literature priors ($W/T \approx 35$). | [`sprints/audit_sprint_k/parameter_provenance.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/parameter_provenance.csv) |

---

## 4. Key Questions Requiring Professor Sahoo's Geological Judgment

Before presenting to the faculty committee or proceeding to Sprint L, the following four questions require Prof. Hiranya Sahoo's specific guidance:

1. **Question 1 (Stratigraphic Datum):**
   *"In Sprint H, you instructed that 'All litholog zeros should be treated as being at the same reference level.' In our code, this was implemented by aligning the base of each measured section ($z=0$ at cliff base). Sahoo et al. (2016) note that along Cliff Face 1, the Axel Anderson coal zone marks the base of the studied interval. Should we formally define our vertical datum as the base of Cliff Face 1 (the Axel Anderson coal zone), or maintain the provisional section-base alignment rule across all 11 spatial wells?"*

2. **Question 2 (Net-to-Gross Target Definition):**
   *"In our forward simulator, does your target Net-to-Gross (e.g. 45%) represent single-facies trough cross-bedded channel sandstone (`sand` alone, which is empirically 42.9% in upstream outcrop), or does it represent the complete channel-belt sandbody (`sand` plus basal planar-laminated `p_sand`, which is empirically 50.5%)? Currently, our generator checks single-facies sand and stops at 38%, but pure channel-belt sand reaches 47.6%."*

3. **Question 3 (Downstream Crevasse Splay Architecture):**
   *"In our downstream outcrop data (Lithologs 9 and 11), rippled heterolithics occupy 21.15% of the column, whereas our forward simulator realizes only 5.0% ripples because splays are currently modeled strictly as lateral wings attached to channel margins. In coastal plain / delta plain settings, are crevasse splays predominantly channel-attached, or should our simulator introduce unattached overbank splay lobe events to model widespread floodplain heteroliths?"*

4. **Question 4 (Channel Bed-Thickness Distribution Model):**
   *"Our Maximum Likelihood fit on the 34 upstream channel sandstone beds gives a Log-Normal distribution ($\mu=1.7870, \sigma=0.8823$, median 5.97 m). A parameter-estimation-aware Lilliefors bootstrap test yields $p=0.42$ (cannot be rejected). However, because these 34 beds originate from only 8 lithologs, cluster-bootstrapping shows that median bed thickness has a 95% confidence interval of 4.5 m to 8.9 m. Do you recommend presenting this log-normal model as a calibrated candidate engineering parameterization, or should we also evaluate alternative distributions (e.g. Gamma or truncated Pareto)?"*

---

## 5. Exact Files Requiring Modification After Approval

Once the faculty supervisor reviews and approves the resolutions, the following files will be updated with the exact corrections detailed below:

### 1. `sprints/audit_sprint_k/SPRINT_K_STAGE_4_SCIENTIFIC_AUDIT_REPORT.md`
- **Section 3.1 (`CLAIM-DATUM-001`):** Replace Axel Anderson coal datum claim with the provisional section-base datum alignment instructed by Prof. Sahoo, noting Cliff Face 1 geological support.
- **Section 3.3 (`CLAIM-HET-001`):** Replace reference to "Litholog 9's 23.0 m continuous ripple unit" with true bed thicknesses (max 6.0 m in L9, max 7.0 m in L11; total ripples 15.0 m and 18.0 m).
- **Section 6.2:** Update Lilliefors p-values ($p = 0.4196 \pm 0.0070$ continuous, $p = 0.6930 \pm 0.0065$ discrete) and cluster-bootstrap median interval ($[4.47, 8.94]\text{ m}$).

### 2. `sprints/audit_sprint_k/scientific_claims_register.csv`
- **Row 6 (`CLAIM-DATUM-001`):** Update `exact_claim_in_report` and `estimation_or_calculation_procedure` to reflect section-base alignment under Prof. Sahoo's instruction.
- **Row 7 (`CLAIM-THICK-002`):** Update uncertainty to reflect cluster-bootstrap interval $[4.47, 8.94]\text{ m}$.
- **Row 17 (`CLAIM-HET-001`):** Correct claim text to state that downstream ripple deficit reflects channel-attached splay architecture and small-sample variance ($N=2$ logs), removing the 23m unit assertion.

### 3. `sprints/audit_sprint_k/SPRINT_K_STAGE_3_VALIDATION_BENCHMARK_REPORT.md`
- **Section 3.3 (line 73):** Remove sentence asserting that Litholog 9 contains a 23.0 m continuous ripple unit. State that Lithologs 9 and 11 contain 15.0 m and 18.0 m of ripples across 4 beds each (averaging 21.15% of the 156 m total section).

### 4. `docs/SCIENTIFIC_METHODS_REGISTER.md`
- **`METHOD-DATUM-001`:** Clarify that `z_common = 0` aligns the source base of each measured section ($z = H_{\text{max}} - d$) per Prof. Sahoo's instruction, while the Axel Anderson coal is observed at the base of Cliff Face 1.

---

## 6. Review Gate Conclusion

The Sprint K Stage 4 audit is verified, all discrepancies are registered with complete empirical evidence, and actionable questions for supervisory review are formulated.

**Per instructions, no implementation of Sprint L has been performed, zero files were committed or pushed, and execution is paused at this checkpoint.**
