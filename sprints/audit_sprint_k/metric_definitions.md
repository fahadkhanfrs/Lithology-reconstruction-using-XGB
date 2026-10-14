# SMALT Metric Definitions and Semantics Register

**Document ID:** `SMALT-METRIC-DEF-001`  
**Version:** 1.0.0 (Audited Sprint K Baseline)  
**Academic Context:** Undergraduate Earth Sciences Research Project (SMALT)  
**Scope:** Mathematical formulas, denominators, facies memberships, and algorithmic stopping rules  

---

## 1. Executive Summary

This document provides unambiguous mathematical definitions, units, denominators, and facies memberships for every quantitative metric used in SMALT descriptive analysis, geological forward simulation, validation campaigns, and synthetic reconstruction benchmarking.

Ambiguities in reservoir sandstone volume reporting (such as "Net-to-Gross") are resolved by establishing a strict tripartite hierarchy that distinguishes single-facies channel sandstone from total channel-belt sandstone and total coarse siliciclastics.

---

## 2. Canonical Facies Membership

All metrics operate on the closed, authoritative 6-state canonical schema defined in Sahoo et al. (2016, Table 1):

| Facies Code | Facies Name | Sedimentological Description | Grain Size / Structure | Color Hex Code |
| :---: | :--- | :--- | :--- | :---: |
| `0` | `sand` | Channel Sandstone (core channel belt) | Medium-coarse, trough cross-stratified | `#f5d442` (Gold) |
| `1` | `p_sand` | Planar Sandstone (basal channel sheet) | Medium-fine, horizontal planar laminated | `#e6a817` (Dark Amber) |
| `2` | `ripples` | Rippled Heterolithics (crevasse splay lobe) | Very fine sand / silt, ripple cross-laminated | `#48bfe3` (Cyan/Sky) |
| `3` | `carbon_mud` | Carbonaceous Mudstone (distal splay margin) | Organic-rich mudstone / siltstone | `#7b2cbf` (Deep Violet) |
| `4` | `coal` | Coal (mire / swamp horizon) | Vitrinitic / clarain organic peat | `#1b1b1b` (Near Black) |
| `5` | `mud` | Overbank Mudstone (floodplain suspension) | Siltstone and clay shale | `#6c757d` (Neutral Slate) |

---

## 3. Tripartite Net-to-Gross Hierarchy

### 3.1 `ntg_channel` (Channel Core Sandstone Fraction)
- **Mathematical Formula:**
  $$\text{ntg\_channel} = \frac{N_{\text{cells}}(\text{Facies } 0)}{N_{\text{total\_cells}}} = \frac{\sum \text{Thickness}(\text{Facies } 0)}{\text{Total Gross Thickness}}$$
- **Facies Membership:** Strictly `{0: sand}` (Trough cross-bedded channel sandstone).
- **Denominator:** Total gross stratigraphic column height (in 1D logs) or total cross-sectional grid cells $N_x \times N_z$ (in 2D realizations).
- **Units:** Dimensionless volume fraction in $[0, 1]$ (reported as percentage in text).
- **Sedimentological Meaning:** Isolates high-energy, high-permeability trough cross-bedded sand deposited in the central channel thalweg. Excludes basal planar sand and rippled heterolithics.
- **Empirical Reference (Upstream Group):** $42.94\%$ ($289.0\text{ m} / 673.0\text{ m}$).
- **Forward Simulator Realization (Upstream Preset, $N=30$):** $38.09\% \pm 2.28\%$.

### 3.2 `ntg_pure` (Pure Channel-Belt Sandstone Fraction)
- **Mathematical Formula:**
  $$\text{ntg\_pure} = \frac{N_{\text{cells}}(\text{Facies } 0) + N_{\text{cells}}(\text{Facies } 1)}{N_{\text{total\_cells}}} = \frac{\sum \text{Thickness}(\text{Facies } 0) + \sum \text{Thickness}(\text{Facies } 1)}{\text{Total Gross Thickness}}$$
- **Facies Membership:** `{0: sand, 1: p_sand}` (Channel sandstone plus basal planar sandstone).
- **Denominator:** Total gross column height or total grid cells $N_x \times N_z$.
- **Units:** Dimensionless volume fraction in $[0, 1]$.
- **Sedimentological Meaning:** Measures the total sand volume enclosed within the macroform channel belt (the complete channel fill). This corresponds directly to the "channelized sandbodies" net-to-gross evaluated in pseudowell transects by Sahoo et al. (2016, Section 5, Figure 13).
- **Empirical Reference (Upstream Group):** $50.52\%$ ($(289.0 + 51.0)\text{ m} / 673.0\text{ m}$).
- **Forward Simulator Realization (Upstream Preset, $N=30$):** $47.55\% \pm 3.13\%$.

### 3.3 `ntg_coarse` (Total Coarse Siliciclastic Fraction)
- **Mathematical Formula:**
  $$\text{ntg\_coarse} = \frac{N_{\text{cells}}(\text{Facies } 0) + N_{\text{cells}}(\text{Facies } 1) + N_{\text{cells}}(\text{Facies } 2)}{N_{\text{total\_cells}}}$$
- **Facies Membership:** `{0: sand, 1: p_sand, 2: ripples}` (Channel sand, basal planar sand, and crevasse splay ripples).
- **Denominator:** Total gross column height or total grid cells $N_x \times N_z$.
- **Units:** Dimensionless volume fraction in $[0, 1]$.
- **Sedimentological Meaning:** Encompasses all permeable sandstone and coarse heterolithics, representing potential reservoir plus thief-zone intervals.
- **Empirical Reference (Upstream Group):** $55.42\%$ ($(289.0 + 51.0 + 33.0)\text{ m} / 673.0\text{ m}$).
- **Forward Simulator Realization (Upstream Preset, $N=30$):** $50.34\% \pm 2.64\%$.

### 3.4 Invariant Mathematical Hierarchy
Across all empirical litholog groups and all synthetic realizations, the following strict mathematical ordering is guaranteed:
$$\text{ntg\_channel} \le \text{ntg\_pure} \le \text{ntg\_coarse}$$
Equality holds if and only if intermediate facies counts (`p_sand` or `ripples`) are exactly zero.

---

## 4. Algorithmic Stopping Rule Mechanics and the Apparent Target Discrepancy

### 4.1 Configured Target vs Stopping Condition
In `ForwardSimulatorConfig`, the target parameter is:
```python
target_net_to_gross: float = 0.45      # Configured Target (45.0%)
tolerance_net_to_gross: float = 0.08   # Stopping Window (8.0%)
max_channel_attempts: int = 50         # Loop Safety Bound
```

The stochastic channel placement loop in `smalt/forward/simulator.py` (lines 280-286) executes as follows:
```python
# Pre-compute current sand proportion
current_sand_cells = np.sum(grid == channel_cfg.core_facies_code)
current_ntg = current_sand_cells / (nx * nz)

# Stopping condition: if target reached after at least 3 channels, stop
if placed_channels >= 3 and current_ntg >= (target_ntg - self.config.tolerance_net_to_gross):
    break
```

### 4.2 Exact Mathematical Breakdown of Stopping Behavior
1. **Quantity Evaluated by Stopping Rule:**
   `current_ntg` evaluates strictly `grid == channel_cfg.core_facies_code` where `core_facies_code = 0` (`sand`).
   It **excludes** `p_sand` (Facies 1).
2. **Effective Stopping Threshold:**
   $$\text{Threshold} = \text{target\_net\_to\_gross} - \text{tolerance\_net\_to\_gross} = 0.45 - 0.08 = \mathbf{0.370} \quad (37.0\%)$$
3. **Loop Termination:**
   As soon as the placed channels bring Facies 0 (`sand` alone) to $\ge 37.0\%$ (provided at least 3 channels have been placed), the generator terminates. Across $N=30$ random seeds, the simulation terminates at a realized `sand` proportion of:
   $$\text{Realized } \text{ntg\_channel} = \mathbf{38.09\% \pm 2.28\%}$$
4. **Contribution of Basal Planar Sand:**
   Each channel scour partitions the lowest $20\%$ of its vertical column to `p_sand` (`basal_fraction = 0.20`).
   Consequently, while `sand` accumulates to $38.09\%$, `p_sand` simultaneously accumulates to **$9.46\% \pm 1.20\%$**.
5. **Realized Pure Sandstone:**
   $$\text{Realized } \text{ntg\_pure} = \text{ntg\_channel} + \text{p\_sand} = 38.09\% + 9.46\% = \mathbf{47.55\% \pm 3.13\%}$$
   This closely matches the empirical field target for pure channel-belt sand ($50.52\%$) with only $2.97\%$ absolute discrepancy ($5.87\%$ relative error).

### 4.3 Erosional Cannibalization and Denominator Invariance
- **Denominator Invariance:** The denominator $N_{\text{total\_cells}} = N_x \times N_z$ (e.g. $200 \times 100 = 20,000$ cells) remains strictly constant throughout the entire simulation.
- **Cannibalization Effect:** When younger channels scour across pre-existing channel sand, older sand cells are replaced by younger sand cells. This superposition adds zero net sand to the domain. Thus, as channel packing density increases, each additional channel yields diminishing net additions to `current_ntg`.
- **Scientific Audit Verdict:** The simulator is functioning as designed. The apparent discrepancy between the 45% target and the 38% realized channel sand is completely explained by:
  - The single-facies definition of `current_ntg` in the stopping rule;
  - The $0.37$ termination threshold created by the tolerance subtraction;
  - The exclusion of the $9.46\%$ basal planar sandstone from the stopping check.

---

## 5. Bed Thickness and Geometric Dimension Metrics

### 5.1 Continuous Empirical Bed Thickness (1D Logs)
- **Definition:** The physical vertical distance between conformable or erosional stratigraphic contacts of different facies in digitized lithologs:
  $$T_{\text{bed}} = z_{\text{top}} - z_{\text{bottom}} \quad (\text{meters})$$
- **Filtering Invariant:** Adjacent intervals with identical canonical facies codes are merged into a single continuous bed prior to thickness analysis (unless crossing an unobserved gap). Zero-thickness intervals are filtered out.
- **Measurement Resolution:**
  - Outcrop logs (L1-L11): Discrete 1.0 m intervals (all 34 channel sandstone beds are exact integers $\ge 1.0\text{ m}$).
  - Core log (L12): Discrete 0.2 m intervals.

### 5.2 Discretized Column Run-Length (2D Forward Grid)
- **Definition:** Continuous vertical sequences of cells sharing the same facies within a single vertical grid column $j$:
  $$T_{\text{run}} = N_{\text{consecutive\_cells}} \times \Delta z \quad (\text{meters})$$
- **Discretization Boundary:**
  - Minimum possible run length is $\Delta z = 0.5\text{ m}$.
  - Any continuous bed thinner than $\Delta z$ cannot be represented on the grid.
  - Sub-meter beds (e.g. 0.2 m coal seams in L12 core) are inflated to $0.5\text{ m}$ (a $150\%$ inflation).
  - Thick channel sandstones ($T > 5.0\text{ m}$) experience less than $2\%$ geometric distortion.

### 5.3 Channel Geometry Metrics
- **Center Location ($x_c$):** Uniform random coordinate $x_c \sim \mathcal{U}(0, X_{\text{extent}})$ in meters.
- **Base Elevation ($z_{\text{base}}$):** Uniform random coordinate $z_{\text{base}} \sim \mathcal{U}(0, Z_{\text{extent}} - T)$ in meters.
- **Channel Scour Thickness ($T$):** Sampled from empirical LogNormal distribution:
  $$T = \text{clip}\left(\exp(\mu + \sigma \cdot Z), T_{\text{min}}, T_{\text{max}}\right) \quad (\text{meters})$$
  with Upstream parameters $\mu = 1.7870, \sigma = 0.8823, T_{\text{min}} = 1.5\text{ m}, T_{\text{max}} = 25.0\text{ m}$.
- **Width-to-Thickness Ratio ($W/T$):** Sampled from Gaussian literature prior:
  $$W/T = \text{clip}\left(\mathcal{N}(\mu_{W/T}, \sigma_{W/T}), 15.0, 70.0\right) \quad (\text{dimensionless})$$
  with baseline $\mu_{W/T} = 35.0, \sigma_{W/T} = 5.0$.
- **Channel Width ($W$):** Derived deterministically as:
  $$W = T \times (W/T) \quad (\text{meters})$$
- **Parabolic Scour Envelope:**
  $$z_{\text{top}}(x) = z_{\text{base}} + T \cdot \max\left(0, 1 - \frac{4(x - x_c)^2}{W^2}\right) \quad (\text{meters})$$

---

## 6. Reconstruction Benchmark Evaluation Metrics

### 6.1 Evaluation Masks and Information Barriers
- **Total Cells ($N_{\text{total}}$):** $N_x \times N_z = 200 \times 100 = 20,000$ cells.
- **Observed Borehole Mask ($M_{\text{obs}}$):** Set of all grid columns $c$ containing a sampled borehole:
  $$c \in \left\{ \text{round}(x_{\text{borehole}} / \Delta x) \right\}$$
- **Evaluation Mask ($M_{\text{eval}}$):** Withholds all observed borehole columns from scoring:
  $$M_{\text{eval}}[:, c] = \text{False} \quad \forall c \in M_{\text{obs}}$$
  Models are evaluated strictly on the $N_{\text{eval}} = N_{\text{total}} - (N_{\text{boreholes}} \times N_z)$ unobserved cells.
- **Interpolation Mask ($M_{\text{interp}}$):** Subsets $M_{\text{eval}}$ strictly within the interior well network:
  $$\min(x_{\text{borehole}}) \le x \le \max(x_{\text{borehole}})$$
- **Extrapolation Mask ($M_{\text{extrap}}$):** Subsets $M_{\text{eval}}$ on the exterior domain margins:
  $$x < \min(x_{\text{borehole}}) \quad \text{or} \quad x > \max(x_{\text{borehole}})$$

### 6.2 Performance Formulas
Let $y_i$ be true withheld facies and $\hat{y}_i$ be predicted facies at cell $i \in M_{\text{eval}}$, with $K=6$ canonical classes:

1. **Raw Accuracy:**
   $$\text{Accuracy}_{\text{raw}} = \frac{1}{N_{\text{eval}}} \sum_{i=1}^{N_{\text{eval}}} \mathbb{I}(y_i = \hat{y}_i)$$
   *Unit:* Fraction in $[0, 1]$ (or percentage). Subject to majority-class mudstone bias.

2. **Balanced Accuracy (Macro-Averaged Recall):**
   $$\text{Accuracy}_{\text{bal}} = \frac{1}{K_{\text{present}}} \sum_{k=1}^K \frac{\sum_{i: y_i = k} \mathbb{I}(\hat{y}_i = k)}{N_k}$$
   where $N_k = \sum_i \mathbb{I}(y_i = k)$ is the true cell count of facies $k$. Prevents dominant background mud from masking poor channel/coal recovery.

3. **Macro-F1 Score:**
   $$\text{Macro-F1} = \frac{1}{K} \sum_{k=0}^{K-1} \frac{2 \cdot P_k \cdot R_k}{P_k + R_k + \epsilon}$$
   where $P_k$ and $R_k$ are per-facies precision and recall, with `zero_division = 0`.

4. **Per-Facies Precision and Recall:**
   $$P_k = \frac{\sum_i \mathbb{I}(y_i = k \wedge \hat{y}_i = k)}{\sum_i \mathbb{I}(\hat{y}_i = k)}, \quad R_k = \frac{\sum_i \mathbb{I}(y_i = k \wedge \hat{y}_i = k)}{\sum_i \mathbb{I}(y_i = k)}$$

---

## 7. Metric Semantics Summary Table

| Metric Name | Numerator | Denominator | Facies Membership | Implementation Code |
| :--- | :--- | :--- | :--- | :--- |
| `ntg_channel` | $\sum \mathbb{I}(\text{grid} == 0)$ | $N_{\text{total\_cells}}$ | `{0: sand}` | `smalt/forward/simulator.py:372` |
| `ntg_pure` | $\sum \mathbb{I}(\text{grid} \in \{0, 1\})$ | $N_{\text{total\_cells}}$ | `{0: sand, 1: p_sand}` | `smalt/forward/simulator.py:373` |
| `ntg_coarse` | $\sum \mathbb{I}(\text{grid} \in \{0, 1, 2\})$ | $N_{\text{total\_cells}}$ | `{0: sand, 1: p_sand, 2: ripples}` | `smalt/forward/simulator.py:374` |
| `sand_prop` | $\sum \mathbb{I}(\text{grid} == 0)$ | $N_{\text{total\_cells}}$ | `{0: sand}` | `smalt/forward/simulator.py:365` |
| `p_sand_prop` | $\sum \mathbb{I}(\text{grid} == 1)$ | $N_{\text{total\_cells}}$ | `{1: p_sand}` | `smalt/forward/simulator.py:365` |
| `ripples_prop` | $\sum \mathbb{I}(\text{grid} == 2)$ | $N_{\text{total\_cells}}$ | `{2: ripples}` | `smalt/forward/simulator.py:365` |
| `carbon_mud_prop` | $\sum \mathbb{I}(\text{grid} == 3)$ | $N_{\text{total\_cells}}$ | `{3: carbon_mud}` | `smalt/forward/simulator.py:365` |
| `coal_prop` | $\sum \mathbb{I}(\text{grid} == 4)$ | $N_{\text{total\_cells}}$ | `{4: coal}` | `smalt/forward/simulator.py:365` |
| `mud_prop` | $\sum \mathbb{I}(\text{grid} == 5)$ | $N_{\text{total\_cells}}$ | `{5: mud}` | `smalt/forward/simulator.py:365` |
| `raw_accuracy` | Correct predictions on withheld cells | $N_{\text{eval\_cells}}$ | All 6 classes | `smalt/forward/benchmark.py:277` |
| `balanced_accuracy` | Unweighted sum of class recalls | Number of present classes | All present classes | `smalt/forward/benchmark.py:278` |
| `macro_f1` | Unweighted average of class F1 scores | $K=6$ classes | All 6 classes | `smalt/forward/benchmark.py:279` |
