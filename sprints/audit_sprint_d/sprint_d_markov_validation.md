# SMALT Sprint D Note: Reproducible 1D Markov Validation & Summary Reconciliation

**Project:** SMALT (Subsurface Stratigraphic Modeling & Active Learning Toolkit)  
**Repository:** `fahadkhanfrs/Lithology-reconstruction-using-XGB`  
**Branch:** `sprint-d`  
**Date:** October 2026  
**Auditor / Engineer:** Geostatistical Software Engineering Pair  

---

## 1. Executive Summary & Sprint Objectives

Sprint D builds directly on Sprint C to establish a mathematically rigorous, fully reproducible, and leakage-free baseline for vertical 1D Markov chain validation across all 12 available lithologs in the Blackhawk Formation dataset.

Sprint D resolves three key technical objectives:
1. **Scoring Formalization:** Strictly separates initial-state log scores ($\ln P_{\text{base}}(s_1)$) from mean transition log scores ($\frac{1}{N_{\text{trans}}} \sum \ln P_{ij}$) and transition perplexity ($\exp(-\text{mean})$) for both regular and embedded chains.
2. **Embedded Bed Integrity:** Formulates the embedded Markov chain over genuine sequences of distinct lithological beds (merging consecutive raw intervals of identical facies), ensuring $P_{ii}^{\text{emb}} \equiv 0$ by construction without conflating self-transition suppression with distinct-bed transitions.
3. **Summary Statistics Reconciliation:** Traces and reconciles every reported descriptive metric directly from the 12 raw litholog CSVs, resolving all prior reporting discrepancies (including the 73 vs 72 sandstone bed count typo, the 664.0 m vs 673.0 m Upstream thickness discrepancy, and Downstream N/G calculations).

All raw interval CSVs (`data/raw_lithologs/`) and previous Sprint C outputs (`audit_sprint_c/`) have been strictly preserved. Versioned deliverables are published in `audit_sprint_d/`.

---

## 2. Mathematical Scoring Formulations

### 2.1 Regular 1D Markov Chain (1 m Discretized Grid)

The regular chain operates on the 1 m discretized grid sequence of each well, sorted upward stratigraphically (from deepest basal grid cell $t=1$ to shallowest top cell $t=T$). Because adjacent 1 m cells frequently belong to the same lithological body, **self-transitions ($s_t = s_{t+1}$) are explicitly retained**.

#### Training Estimation:
For each Leave-One-Litholog-Out (LOLO) fold $f \in \{1, \dots, 12\}$, the training transition matrix $P^{(f)}$ is estimated exclusively from the remaining 11 training lithologs using Laplace (additive) smoothing:

$$
P_{ij}^{(f)} = \frac{N_{ij}^{(f)} + \alpha}{\sum_{k=1}^K N_{ik}^{(f)} + K\alpha}, \qquad \alpha = 0.1, \quad K = 5
$$

where $N_{ij}^{(f)}$ is the tally of transitions from facies $i$ to facies $j$ across all training profiles.

#### Evaluation Metrics:
Let $s = (s_1, s_2, \dots, s_T)$ be the target sequence of length $T$, with $N_{\text{trans}} = T - 1$ observed step transitions.

1. **Empirical Initial State (Basal Cell) Log Score:**
   $$
   \ln P_{\text{base}}^{(f)}(s_1), \quad \text{where} \quad P_{\text{base}}^{(f)}(i) = \frac{N_{\text{base}, i}^{(f)} + \alpha}{N_{\text{wells}}^{(f)} + K\alpha}
   $$
   The initial state is evaluated against the empirical distribution of basal facies observed across training wells, not conflated with the long-term stationary distribution $\pi^{(f)}$.

2. **Total Transition Log Score:**
   $$
   \mathcal{S}_{\text{trans}}(s) = \sum_{t=1}^{T-1} \ln P_{s_t, s_{t+1}}^{(f)}
   $$

3. **Mean Transition Log Score:**
   $$
   \bar{\ell}_{\text{trans}}(s) = \frac{1}{T-1} \sum_{t=1}^{T-1} \ln P_{s_t, s_{t+1}}^{(f)}
   $$

4. **Transition Perplexity:**
   $$
   \mathcal{PP}_{\text{trans}}(s) = \exp\left( -\bar{\ell}_{\text{trans}}(s) \right) = \left( \prod_{t=1}^{T-1} P_{s_t, s_{t+1}}^{(f)} \right)^{-\frac{1}{T-1}}
   $$
   Perplexity represents the effective branching factor of the Markov model per 1 m vertical step.

5. **Total Sequence Log Score:**
   $$
   \mathcal{S}_{\text{seq}}(s) = \ln P_{\text{base}}^{(f)}(s_1) + \mathcal{S}_{\text{trans}}(s)
   $$
   The initial-state score is reported separately and never divided into the mean transition likelihood.

---

### 2.2 Embedded 1D Markov Chain (Distinct Bed Sequence)

The embedded Markov chain describes the succession of distinct lithological beds, stepping upward stratigraphically. It does not model within-bed persistence (thickness); it models **which facies succeeds the current bed once a lithological boundary is crossed**.

#### Bed Sequence Construction & Gap Handling:
Consecutive raw intervals sharing the exact same facies are merged into a single distinct bed. 
- **Thickness Accumulation:** When merging intervals, rock thicknesses are directly accumulated ($\text{thickness} = \sum (z_{\text{bottom}} - z_{\text{top}})$). Unrecorded gaps (e.g., the 18-19 m gap in Litholog 9) are never absorbed into bed thickness.
- **Self-Transition Suppression:** By construction, $s_t \ne s_{t+1}$ for all $t$ in the distinct bed sequence. Consequently, the embedded transition matrix has a **strictly zero diagonal**:
  $$
  P_{ii}^{\text{emb}} \equiv 0 \quad \forall i \in \{0, \dots, K-1\}
  $$

#### Training Estimation:
Additive smoothing is applied strictly across the $K-1$ alternative off-diagonal destinations:

$$
P_{ij}^{\text{emb}, (f)} = \begin{cases} 
0, & j = i \\
\frac{N_{ij}^{\text{emb}, (f)} + \alpha}{\sum_{k \ne i} N_{ik}^{\text{emb}, (f)} + (K - 1)\alpha}, & j \ne i
\end{cases}
$$

#### Evaluation Metrics:
For an embedded sequence of $B$ distinct beds ($N_{\text{trans}} = B - 1$ genuine boundary transitions):
- **Basal Bed Initial Score:** $\ln P_{\text{base}}^{\text{emb}, (f)}(s_1)$, evaluated against the empirical basal bed distribution of the training wells.
- **Mean Bed-Transition Log Score:** $\bar{\ell}_{\text{trans}}^{\text{emb}}(s) = \frac{1}{B-1} \sum_{t=1}^{B-1} \ln P_{s_t, s_{t+1}}^{\text{emb}, (f)}$.
- **Bed-Transition Perplexity:** $\mathcal{PP}_{\text{trans}}^{\text{emb}}(s) = \exp(-\bar{\ell}_{\text{trans}}^{\text{emb}}(s))$.

#### Critical Non-Equivalence Principle:
Regular-chain perplexity ($\sim 1.5 - 2.5$) and embedded-chain perplexity ($\sim 2.8 - 6.5$) **must never be directly compared**. 
- Regular perplexity evaluates step-by-step 1 m vertical persistence, where self-transition diagonal probabilities are high ($\sim 0.70 - 0.90$).
- Embedded perplexity evaluates boundary-crossing choice among $K-1 = 4$ distinct lithological alternatives. A perplexity of 3.5 in an embedded chain indicates an effective branching factor of 3.5 out of 4 possible facies transitions.

---

## 3. Sandstone Bed Statistics Reconciliation

A primary mandate of Sprint D was to trace and reconcile all sandstone bed metrics directly from source CSVs and resolve discrepancies in prior reporting.

### 3.1 Reconciliation Summary Table

| Metric | Raw Intervals (Sprint E Reconciled) | Merged Lithosomes (Sprint E Reconciled) | Sprint C Reported Baseline | Status | Audit Explanation |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Total Sandstone Count** | **135** | **99** | 73 (narrative) / 72 (sum) | **Reconciled in Sprint E** | Due to user digitization of planar sandstones (`p_sand`) across L4-L11, raw intervals contain 135 sandstone intervals. Merging consecutive sand intervals produces 99 distinct lithosomes. |
| **Cumulative Thickness (m)** | **507.70** | **507.70** | 444.70 | **Reconciled in Sprint E** | Total sandstone thickness across all 12 lithologs is 507.70 m in both raw intervals and merged lithosomes. |
| **Minimum Bed Thickness (m)** | **1.00** | **1.00** | 1.40 | **Reconciled in Sprint E** | Observed in L4, L5, L6, L7, L8, L9, L10, L11 as 1.0 m thin sandstone beds. |
| **Median Bed Thickness (m)** | **3.00** | **2.00** | 5.05 | **Reconciled in Sprint E** | Raw interval median is 3.00 m. For merged distinct lithosomes, median thickness is 2.00 m due to thin interbedded units. |
| **Mean Bed Thickness (m)** | **3.76** | **5.13** | 6.18 | **Reconciled in Sprint E** | Raw interval mean is 3.76 m ($507.70 / 135$). Merged distinct lithosomes mean is 5.13 m ($507.70 / 99$). |
| **Maximum Bed Thickness (m)** | **17.00** | **32.00** | 17.00 | **Reconciled** | Raw interval maximum is 17.00 m (Litholog 2, 40 - 57 m). When consecutive sand intervals are merged, Litholog 3 contains a 32.00 m amalgamated channel-sand package (61 - 93 m). |
| **Sample Std Dev (ddof=1) (m)** | **2.96** | **5.90** | 3.36 | **Reconciled in Sprint E** | Sample standard deviation for raw intervals is 2.96 m. For merged lithosomes, sample standard deviation is 5.90 m. |
| **Population Std Dev (ddof=0) (m)**| **2.95** | **5.87** | 3.34 | **Reconciled in Sprint E** | Population standard deviation is 2.95 m for raw intervals; 5.87 m for merged lithosomes. |

*Artifact Reference:* [`sprints/audit_sprint_d/sandstone_thickness_reconciliation.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_d/sandstone_thickness_reconciliation.csv)

---

## 4. Stratigraphic Group Reconciliation

Sprint C Table 7.1 reported group summary statistics with minor discrepancies caused by hardcoded plot labels and manual entries. In Sprint E, all group metrics were recomputed directly from the revised source intervals.

### 4.1 Group Reconciliation Table

| Group Name | Wells Included | Total Thickness (m) | Sand Thickness (m) | Pure Sand N/G | Historical Reported N/G | Discrepancy & Root Cause |
| :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| **Upstream** | L2, L3, L4, L5, L6, L7, L8, L10 (8 logs) | **673.0** | **340.0** | **50.52%** | 43.22% | In Sprint E revision, distinct `p_sand` intervals increased pure sand thickness from 287.0 m to 340.0 m across the 673.0 m Upstream transect. |
| **Downstream Outcrop** | L9, L11 (2 logs) | **156.0** | **75.0** | **48.08%** | 41.67% | Total sand increased to 75.0 m ($35.0 + 40.0\text{ m}$) due to revised `p_sand` intervals. |
| **Downstream Composite** | L9, L11, L12 (3 logs) | **267.0** | **132.7** | **49.70%** | 45.96% / 47.19% | Total thickness is $156.0 + 111.0 = 267.0\text{ m}$. Total sand is $75.0 + 57.7 = 132.7\text{ m}$, yielding pure N/G of 49.70%. |
| **All Outcrop** | L1 - L11 (11 logs) | **922.0** | **450.0** | **48.81%** | 41.97% | Cumulative outcrop thickness 922.0 m, pure sandstone 450.0 m. |
| **All 12 Lithologs** | L1 - L12 (12 logs) | **1033.0** | **507.7** | **49.15%** | 43.05% | Cumulative study thickness 1033.0 m, pure sandstone 507.70 m. |

*Artifact Reference:* [`sprints/audit_sprint_d/reconciled_group_statistics.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_d/reconciled_group_statistics.csv)

### 4.2 Group Facies Proportions (Continuous Stratigraphy - Sprint E Reconciled)

| Facies | Upstream (673.0 m) | Downstream Outcrop (156.0 m) | Downstream Composite (267.0 m) | All 12 Lithologs (1033.0 m) |
| :--- | :---: | :---: | :---: | :---: |
| **Coal** | 1.93% (13.0 m) | 1.92% (3.0 m) | 4.04% (10.8 m) | 2.98% (30.8 m) |
| **Channel Sandstone** | 50.52% (340.0 m) | 48.08% (75.0 m) | 49.70% (132.7 m) | 49.15% (507.7 m) |
| **Carbonaceous Mudstone** | 0.89% (6.0 m) | 4.49% (7.0 m) | 3.97% (10.6 m) | 1.90% (19.6 m) |
| **Siltstone** | 4.90% (33.0 m) | 21.15% (33.0 m) | 16.03% (42.8 m) | 7.63% (78.8 m) |
| **Overbank Mudstone** | 41.75% (281.0 m) | 24.36% (38.0 m) | 26.25% (70.1 m) | 38.34% (396.1 m) |

---

## 5. Leakage-Free Leave-One-Litholog-Out (LOLO) Results

Validation was conducted strictly out-of-fold. For each target litholog $k$:
- $k$ was completely excluded from transition tallies, basal distributions, and smoothing prior fits.
- Training transition matrices $P^{(f)}$ and basal distributions $p_{\text{base}}^{(f)}$ were fitted strictly on $\{L_j : j \ne k\}$.
- Synthetic Gamma Ray and ground truth `prev_facies` were excluded.

### 5.1 Regular 1D Markov LOLO Results (1 m Discretized Grid)

| Target Log | Sequence Length | Transitions | Basal Facies | Basal Log Score | Total Trans Log Score | Mean Trans Log Score | Trans Perplexity | Total Seq Score | Gap Transitions | Gap-Masked Perplexity | Stationary TV Dist |
| :--- | :---: | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **litholog1** | 93 | 92 | sand | -1.7004 | -65.2631 | -0.7094 | **2.0327** | -66.9635 | 0 | 2.0327 | 0.2354 |
| **litholog2** | 93 | 92 | sand | -1.7004 | -57.7919 | -0.6282 | **1.8742** | -59.4923 | 0 | 1.8742 | 0.1780 |
| **litholog3** | 93 | 92 | sand | -1.7004 | -35.7899 | -0.3890 | **1.4755** | -37.4903 | 0 | 1.4755 | 0.3363 |
| **litholog4** | 84 | 83 | mud | -0.6341 | -53.6765 | -0.6467 | **1.9092** | -54.3106 | 0 | 1.9092 | 0.1557 |
| **litholog5** | 85 | 84 | mud | -0.6341 | -56.2831 | -0.6700 | **1.9543** | -56.9172 | 0 | 1.9543 | 0.0477 |
| **litholog6** | 82 | 81 | mud | -0.6341 | -60.4687 | -0.7465 | **2.1097** | -61.1028 | 0 | 2.1097 | 0.1699 |
| **litholog7** | 80 | 79 | mud | -0.6341 | -56.2669 | -0.7122 | **2.0386** | -56.9010 | 0 | 2.0386 | 0.1789 |
| **litholog8** | 79 | 78 | mud | -0.6341 | -70.7429 | -0.9070 | **2.4768** | -71.3770 | 0 | 2.4768 | 0.2651 |
| **litholog9** | 78 | 77 | mud | -0.6341 | -58.1324 | -0.7550 | **2.1275** | -58.7664 | 2 | 2.1635 | 0.1470 |
| **litholog10** | 77 | 76 | mud | -0.6341 | -59.5534 | -0.7836 | **2.1893** | -60.1874 | 0 | 2.1893 | 0.1564 |
| **litholog11** | 78 | 77 | silt | -2.3470 | -53.1969 | -0.6909 | **1.9954** | -55.5439 | 2 | 1.9544 | 0.2033 |
| **litholog12** | 111 | 110 | silt | -2.3470 | -154.1215 | -1.4011 | **4.0597** | -156.4686 | 0 | 4.0597 | 0.1667 |

*Average Outcrop Regular Perplexity (L1 - L11):* **2.0166**  
*Full Dataset Average Regular Perplexity (L1 - L12):* **2.1869**

---

### 5.2 Embedded 1D Markov LOLO Results (Distinct Bed Sequence)

| Target Log | Distinct Beds | Transitions | Basal Facies | Basal Log Score | Total Trans Log Score | Mean Trans Log Score | Bed Trans Perplexity | Total Seq Score | Gap Transitions | Gap-Masked Perplexity | Stationary TV Dist |
| :--- | :---: | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **litholog1** | 20 | 19 | sand | -1.7004 | -23.8454 | -1.2550 | **3.5079** | -25.5458 | 0 | 3.5079 | 0.2739 |
| **litholog2** | 15 | 14 | sand | -1.7004 | -26.1240 | -1.8660 | **6.4624** | -27.8244 | 0 | 6.4624 | 0.2834 |
| **litholog3** | 10 | 9 | sand | -1.7004 | -9.6734 | -1.0748 | **2.9295** | -11.3739 | 0 | 2.9295 | 0.2520 |
| **litholog4** | 18 | 17 | mud | -0.6341 | -19.4568 | -1.1445 | **3.1410** | -20.0909 | 0 | 3.1410 | 0.1581 |
| **litholog5** | 18 | 17 | mud | -0.6341 | -21.1887 | -1.2464 | **3.4779** | -21.8228 | 0 | 3.4779 | 0.0898 |
| **litholog6** | 16 | 15 | mud | -0.6341 | -21.6800 | -1.4453 | **4.2431** | -22.3141 | 0 | 4.2431 | 0.1983 |
| **litholog7** | 15 | 14 | mud | -0.6341 | -18.8604 | -1.3472 | **3.8466** | -19.4945 | 0 | 3.8466 | 0.1952 |
| **litholog8** | 26 | 25 | mud | -0.6341 | -26.2847 | -1.0514 | **2.8617** | -26.9188 | 0 | 2.8617 | 0.2223 |
| **litholog9** | 17 | 16 | mud | -0.6341 | -20.3440 | -1.2715 | **3.5661** | -20.9781 | 1 | 3.4977 | 0.1584 |
| **litholog10** | 19 | 18 | mud | -0.6341 | -19.1481 | -1.0638 | **2.8975** | -19.7822 | 0 | 2.8975 | 0.1654 |
| **litholog11** | 16 | 15 | silt | -2.3470 | -19.4208 | -1.2947 | **3.6501** | -21.7679 | 1 | 3.5923 | 0.1718 |
| **litholog12** | 63 | 62 | silt | -2.3470 | -94.5218 | -1.5245 | **4.5929** | -96.8689 | 0 | 4.5929 | 0.1601 |

*Average Outcrop Embedded Perplexity (L1 - L11):* **3.6894**  
*Full Dataset Average Embedded Perplexity (L1 - L12):* **3.7647**

*Artifact Reference:* [`audit_sprint_d/lolo_validation_results_sprint_d.csv`](file:///d:/Lithology-reconstruction-using-XGB/audit_sprint_d/lolo_validation_results_sprint_d.csv)

---

## 6. Sensitivity Analysis: Gap Masking

Only two lithologs in the dataset possess unrecorded physical measurement gaps:
1. **Litholog 9:** 1 m gap between 18.0 m and 19.0 m (unrecorded interval between raw intervals).
2. **Litholog 11:** 1 m gap between 59.0 m and 60.0 m (recording gap in manual field log).

All other 10 lithologs have zero gap-crossing transitions ($N_{\text{gap}} = 0$).

### 6.1 Sensitivity Findings

| Litholog ID | Model Form | Evaluated Transitions | Gap-Crossing Transitions | Standard Mean Log Score | Standard Perplexity | Gap-Masked Mean Log Score | Gap-Masked Perplexity | Delta Perplexity | Audit Interpretation |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **litholog9** | Regular 1m | 77 | 2 | -0.7550 | 2.1275 | -0.7717 | 2.1635 | **+0.0360** | Masking the 2 gap-crossing cells slightly increases perplexity (+0.036), confirming negligible impact. |
| **litholog9** | Embedded Beds | 16 | 1 | -1.2715 | 3.5661 | -1.2521 | 3.4977 | **-0.0684** | Masking the single gap-crossing transition slightly improves perplexity (-0.068), well within noise. |
| **litholog11** | Regular 1m | 77 | 2 | -0.6909 | 1.9954 | -0.6701 | 1.9544 | **-0.0410** | Masking the 2 gap-crossing cells slightly improves perplexity (-0.041). |
| **litholog11** | Embedded Beds | 15 | 1 | -1.2947 | 3.6501 | -1.2788 | 3.5923 | **-0.0578** | Masking the single gap-crossing transition slightly improves perplexity (-0.058). |

**Conclusion:** Across all models, gap-masking alters perplexity by less than $\pm 0.07$ and mean log scores by less than $\pm 0.02$. The default continuous treatment is scientifically robust and introduces no qualitative distortion.

*Artifact Reference:* [`audit_sprint_d/sensitivity_gap_masking_results.csv`](file:///d:/Lithology-reconstruction-using-XGB/audit_sprint_d/sensitivity_gap_masking_results.csv)

---

## 7. Status of Special Lithologs (L1 and L12)

### 7.1 Litholog 1 (Unassigned Coordinates)
- **Vertical Eligibility:** Fully eligible and included in all vertical descriptive summaries, Markov transition tallies, and vertical LOLO validation folds.
- **Spatial Restriction:** Because no physical coordinates or elevation datum exist in the repository or Sahoo et al. (2016) for L1, **it must remain strictly excluded from any future 2D/3D spatial modeling, spatial conditioning, or spatial cross-validation**.

### 7.2 Litholog 12 (Digitized Core Log)
- **Stratigraphic Scope:** Strictly limited to the digitized $0 - 111\text{ m}$ interval present in `data/raw_lithologs/litholog12.csv`.
- **Core Span Limitation:** The original physical core spans 242 m, but depths $> 111\text{ m}$ were never digitized into the repository. No unobserved depths $> 111\text{ m}$ may be assumed or invented.
- **Downstream Distality:** L12 exhibits higher regular perplexity (4.0597) and embedded perplexity (4.5929) relative to the outcrop logs due to higher thin-bed interbedding (63 distinct beds across 111 m, average bed thickness 1.76 m).

---

## 8. Test Verification Summary

A dedicated pytest suite [`tests/test_sprint_d.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_sprint_d.py) was implemented covering all 10 requirements of Sprint D Section 6.

All 44 tests across the entire repository now pass cleanly:
1. `tests/test_phase0_loader.py`: 11 passed
2. `tests/test_phase1_markov.py`: 11 passed
3. `tests/test_sprint_c.py`: 12 passed
4. `tests/test_sprint_d.py`: 10 passed
   - Regular transition scoring retains self-transitions: **PASSED**
   - Embedded evaluation uses distinct-bed transitions only: **PASSED**
   - Embedded transition matrices have strictly zero diagonal within $10^{-12}$: **PASSED**
   - All transition-matrix rows sum to 1 within $10^{-7}$: **PASSED**
   - Target litholog is 100% excluded from fitted training counts: **PASSED**
   - Gap-crossing transitions handled per documented policy: **PASSED**
   - Fixed toy sequence gives analytically verified log scores down to $10^{-12}$: **PASSED**
   - Continuous bed thickness and 1m discretized proportions reported separately: **PASSED**
   - Deterministic training and scoring across runs: **PASSED**
   - Deliverables exist in `audit_sprint_d/`: **PASSED**

---

## 9. Non-Claim Prohibitions & Scientific Boundaries

To preserve strict scientific integrity, the following boundaries are formally recorded:
1. **Transition Likelihood is Not Reconstruction Accuracy:** Out-of-fold transition log-likelihood, perplexity, and stationary-distribution divergence are 1D vertical succession diagnostics. They **must not be claimed as classification accuracy, predictive success, or evidence of inter-well spatial reconstruction**.
2. **Stationary Distribution is Descriptive Only:** Stationary vector $\pi$ describes asymptotic proportion under Markovian assumptions; it is not a machine learning classification metric.
3. **No Spatial Inferences:** 1D vertical Markov models do not establish lateral continuity, channel belt width, or 2D connectivity. Spatial modeling remains blocked pending coordinate harmonization and stratigraphic datum definition.
