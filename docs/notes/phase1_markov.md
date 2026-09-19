# Phase 1 Technical Note: 1D Vertical Markov Transition Matrices and Facies Succession Analysis

## 1. Sedimentological Justification & Geological Context

### 1.1 Stratigraphic Succession as a Discrete-Step Markov Process
Vertical stratigraphic columns preserve the chronological sequence of depositional environments according to the principle of stratigraphic superposition. When a continuous core or outcrop litholog is digitized at regular vertical increments $\Delta z$ (or at discrete lithological bed boundaries), the succession of lithofacies states can be formulated as a discrete-space, discrete-step Markov process (Krumbein & Dacey, 1969; Doveton, 1971).

Let $S_t \in \{0, 1, \dots, K-1\}$ denote the lithofacies state observed at stratigraphic level $t$. Under the first-order Markovian assumption, the conditional probability of encountering facies $S_t$ given the complete depositional history depends strictly on the immediately preceding underlying facies $S_{t-1}$:

$$P(S_t = s_t \mid S_{t-1} = s_{t-1}, S_{t-2} = s_{t-2}, \dots, S_0 = s_0) = P(S_t = s_t \mid S_{t-1} = s_{t-1})$$

In stratigraphic data processing, depth is traversed from larger depth to smaller depth ($z \to z - \Delta z$) so that each transition represents a lower-to-upper stratigraphic step, preserving the true chronological order of sedimentation.

### 1.2 Asymmetric Fluvial-Deltaic Cycles vs. Random Transitions
In fluvial and deltaic depositional systems (such as the Upper Cretaceous Blackhawk Formation and Ferron Sandstone analogs), vertical transitions between lithologies are non-random and directionally asymmetric:
- **Fining-Upward Channel Fills**: High-energy channel complexes typically initiate with an erosional basal scour, followed by coarse-to-medium trough cross-stratified channel sandstones (high hydrodynamic energy). As the channel migrates laterally or fills during waning flow and abandonment, sedimentation transitions upward into fine-grained ripple-laminated sandstones, siltstones, and floodplain/overbank mudstones and mire coals. This generates a pronounced asymmetric transition pathway:
  $$\text{Channel Sandstone} \longrightarrow \text{Fine Sandstone/Splay} \longrightarrow \text{Siltstone} \longrightarrow \text{Overbank Mudstone / Coal}$$
- **Coarsening-Upward Progradational Packages**: Conversely, crevasse splays, delta front progradation, and bay-fill packages exhibit coarsening-upward successions (mudstone $\to$ siltstone $\to$ splay sandstone) as sediment lobes prograde into interdistributary bays or floodplain lakes.

Because natural sedimentation involves both fining-upward channel avulsion and coarsening-upward splay/delta progradation, Markov transition analysis provides an objective empirical measure of the dominant vertical succession statistics rather than imposing an artificial assumption of universal symmetry.

### 1.3 SMALT 5-State Facies Aggregation
Sahoo et al. (2016) classified the fluvio-deltaic stratigraphy into six distinct lithofacies:
1. Trough cross-stratified sandstone (high-energy channel core)
2. Parallel-laminated sandstone (upper-stage plane bed / channel bar top)
3. Heterolithic mudstone/siltstone/rippled sandstone (channel margin, levee, crevasse splay)
4. Mudstone/siltstone (floodplain overbank)
5. Carbonaceous mudstone (poorly drained marsh / mire fringe)
6. Coal (waterlogged peat mire)

The Sahoo lithologs constitute the geological observational dataset. SMALT's 5-state encoding represents a project-defined aggregation designed for robust stochastic simulation and active learning:
- **State 0 (Coal)**: Sahoo Facies 6 (Coal)
- **State 1 (Channel Sandstone)**: Sahoo Facies 1 & Facies 2 (Trough cross-stratified and parallel-laminated sandstones)
- **State 2 (Fine Sandstone / Splay)**: Sahoo Facies 3 & Facies 5 (Heterolithic splay sandstone and carbonaceous mudstone)
- **State 3 (Siltstone)**: Sahoo Facies 4 (Coarser siltstone fraction)
- **State 4 (Overbank Mudstone)**: Sahoo Facies 4 (Fine mudstone / floodplain shale)

---

## 2. Mathematical Derivation

### 2.1 Transition Count Matrix
Given $W$ digitized lithologs, let $s_{w, t}$ denote the facies state at vertical step $t$ in well $w \in \{1, \dots, W\}$, ordered from bottom to top (decreasing depth $z$). The total transition count from facies $i$ to facies $j$ is tallied into the $K \times K$ transition count matrix $\mathbf{N} = [n_{ij}]$:

$$n_{ij} = \sum_{w=1}^W \sum_{t=1}^{T_w - 1} \mathbb{I}(s_{w, t} = i \;\land\; s_{w, t+1} = j)$$

where $T_w$ is the number of vertical observations in well $w$, and $\mathbb{I}(\cdot)$ is the indicator function.

### 2.2 Maximum Likelihood Estimator: Regular (Fixed-Step) Chain
In a fixed-step Markov chain ($\Delta z = \text{constant}$, e.g. 1 m), self-transitions ($i \to i$) are retained. The diagonal elements $n_{ii}$ capture the vertical thickness and persistence of individual lithologic units. With Laplace smoothing parameter $\alpha \ge 0$, the maximum likelihood transition probability estimator is:

$$P_{ij} = \frac{n_{ij} + \alpha}{\sum_{k=0}^{K-1} (n_{ik} + K\alpha)}, \qquad \forall i, j \in \{0, \dots, K-1\}$$

Row-stochasticity requires that each row sums strictly to unity:

$$\sum_{j=0}^{K-1} P_{ij} = 1.0, \qquad \forall i \in \{0, \dots, K-1\}$$

### 2.3 Embedded (Boundary-Crossing) Markov Chain
An embedded Markov chain isolates pure lithologic succession by removing thickness-induced self-transitions ($P_{ii} = 0$) before normalization and smoothing (Krumbein & Dacey, 1969). The diagonal counts are strictly zeroed ($n_{ii} = 0$), and smoothing is applied exclusively across the remaining $K-1$ states:

$$P_{ij} = \begin{cases} 
\dfrac{n_{ij} + \alpha}{\sum_{k \neq i} (n_{ik} + (K-1)\alpha)}, & \text{if } i \neq j \\[10pt]
0, & \text{if } i = j 
\end{cases}$$

### 2.4 Unobserved / Absorbing State Fallback
If facies state $i$ has zero observed outgoing transitions in the dataset ($\sum_k n_{ik} = 0$) and $\alpha = 0$, standard normalization encounters a division-by-zero anomaly. To prevent numerical breakdown during downstream stochastic simulation:
- **Regular Chain**: Assign uniform probability across all $K$ states:
  $$P_{ij} = \frac{1}{K}, \qquad \forall j \in \{0, \dots, K-1\}$$
- **Embedded Chain**: Assign uniform probability across the $K-1$ off-diagonal states:
  $$P_{ij} = \begin{cases} \dfrac{1}{K-1}, & i \neq j \\[6pt] 0, & i = j \end{cases}$$

### 2.5 Stationary Distribution ($\boldsymbol{\pi}$)
The stationary distribution $\boldsymbol{\pi} = [\pi_0, \pi_1, \dots, \pi_{K-1}]$ represents the invariant left eigenvector of the transition matrix $\mathbf{P}$ corresponding to eigenvalue $\lambda = 1$:

$$\boldsymbol{\pi} \mathbf{P} = \boldsymbol{\pi} \quad \iff \quad \mathbf{P}^T \boldsymbol{\pi}^T = \boldsymbol{\pi}^T$$

subject to the probability simplex constraints:

$$\sum_{i=0}^{K-1} \pi_i = 1.0 \quad \text{and} \quad \pi_i \ge 0, \quad \forall i$$

#### Physical Interpretation
For an ergodic regular chain with step size $\Delta z$, $\boldsymbol{\pi}$ represents the theoretical asymptotic facies occupancy distribution (long-run volumetric proportions) in the depositional sequence. Comparing $\boldsymbol{\pi}$ against the empirical facies proportions $\mathbf{p}_{\text{emp}}$ serves as a diagnostic test of chain ergodicity and finite-sample boundary effects:

$$\epsilon_{\text{rel}}(i) = \frac{|\pi_i - p_{\text{emp}, i}|}{p_{\text{emp}, i}}$$

### 2.6 Row/State Transition Entropy
To quantify the predictability and disorder of transitions originating from each facies state, we define the row/state transition entropy $H_i$:

$$H_i = - \sum_{j=0, P_{ij} > 0}^{K-1} P_{ij} \ln P_{ij}$$

The theoretical maximum entropy occurs when all transitions from state $i$ are equiprobable ($P_{ij} = 1/K$ for regular, $1/(K-1)$ for embedded):

$$H_i^{\max} = \begin{cases} \ln K, & \text{regular chain} \\ \ln(K-1), & \text{embedded chain} \end{cases}$$

The normalized state transition entropy is:

$$H_i^{\text{norm}} = \frac{H_i}{H_i^{\max}} \in [0, 1]$$

A low $H_i^{\text{norm}}$ indicates high deterministic predictability (e.g., strong selective transition to a specific succeeding facies), whereas $H_i^{\text{norm}} \approx 1$ indicates maximum geological uncertainty.

### 2.7 Directional Asymmetry Matrix
To isolate directional cyclicity and determine whether vertical transitions favor fining-upward over coarsening-upward pathways, we compute the directional asymmetry matrix $\mathbf{A}$:

$$\mathbf{A} = \mathbf{P} - \mathbf{P}^T, \qquad A_{ij} = P_{ij} - P_{ji}$$

- $A_{ij} > 0$: The upward transition $i \to j$ is favored over the downward reverse $j \to i$.
- $A_{ij} < 0$: The reverse transition is favored.
- $A_{ij} = 0$: Reversible / symmetric transition.

---

## 3. Formal Literature References

1. **Krumbein, W. C., & Dacey, M. F. (1969).** Markov chains and embedded Markov chains in geology. *Mathematical Geology*, 1(1), 79-96. https://doi.org/10.1007/BF02047072
2. **Doveton, J. H. (1971).** An application of Markov chain analysis to the Ayrshire Coal Measures succession. *Scottish Journal of Geology*, 7(1), 11-27. https://doi.org/10.1144/sjg07010011
3. **Sahoo, H., Gani, M. R., & Gani, N. D. (2016).** 3D facies architecture and sequence stratigraphy of a fluvio-deltaic succession, Cretaceous Ferron Sandstone, Utah. *Sedimentology*, 63(6), 1403-1437. https://doi.org/10.1111/sed.12267
4. **Elfeki, A., & Dekking, M. (2001).** A Markov chain model for subsurface characterization: theory and applications. *Mathematical Geology*, 33(5), 569-589.
5. **Carle, S. F., & Fogg, G. E. (1996).** Transition probability-based geostatistics. *Mathematical Geology*, 28(4), 453-476.

---

## 4. Empirical Recomputed Results (Complete 11-Litholog Dataset)

> [!NOTE]
> Previous Phase 1 results generated from the incomplete 3-well subset (lithologs 1, 9, 11; 247 observations) are **OBSOLETE**. The authoritative empirical statistics below reflect the complete 11-litholog working dataset (920 observations, 909 regular transitions, 201 boundary crossings).

### 4.1 Regular Transition Matrix $P_{\text{reg}}$ (Fixed-Step 1m, 11 Logs)
```
                         Coal  Channel Sandstone  Fine Sandstone / Splay  Siltstone  Overbank Mudstone
Coal                    0.286              0.071                   0.286      0.036              0.321
Channel Sandstone       0.031              0.904                   0.026      0.018              0.021
Fine Sandstone / Splay  0.000              0.052                   0.759      0.043              0.147
Siltstone               0.027              0.082                   0.082      0.466              0.342
Overbank Mudstone       0.020              0.065                   0.013      0.082              0.820
```
- **Matrix Condition Number $\kappa(P_{\text{reg}})$**: 4.52
- **Row Stochasticity**: Strictly verified $|\sum_j P_{ij} - 1| \le 1 \times 10^{-7}$

### 4.2 Embedded Transition Matrix $P_{\text{emb}}$ (Boundary Crossings, $P_{ii} = 0$, 11 Logs)
```
                         Coal  Channel Sandstone  Fine Sandstone / Splay  Siltstone  Overbank Mudstone
Coal                    0.000              0.100                   0.400      0.050              0.450
Channel Sandstone       0.324              0.000                   0.270      0.189              0.216
Fine Sandstone / Splay  0.000              0.214                   0.000      0.179              0.607
Siltstone               0.051              0.154                   0.154      0.000              0.641
Overbank Mudstone       0.109              0.364                   0.073      0.455              0.000
```
- **Matrix Condition Number $\kappa(P_{\text{emb}})$**: 17.82

### 4.3 Directional Upward Fining Invariant
Channel Sandstone (State 1) upward transitions strictly favor finer-grained facies:
- Sand $\to$ Coal: 0.324
- Sand $\to$ Fine Sandstone / Splay: 0.270
- Sand $\to$ Siltstone: 0.189
- Sand $\to$ Overbank Mudstone: 0.216
- **Combined Upward Fining / Abandonment Transitions**: **100.0%** (1.000)

### 4.4 Stationary Facies Occupancy Diagnostic
| Facies State | Stationary $\pi_i$ | Empirical $p_{\text{emp}, i}$ | Absolute Diff | Relative Diff (%) |
|---|---|---|---|---|
| 0 (Coal) | 0.0305 | 0.0304 | 0.0001 | 0.29% |
| 1 (Channel Sandstone) | 0.4043 | 0.4196 | 0.0153 | 3.64% |
| 2 (Fine Sand / Splay) | 0.1264 | 0.1261 | 0.0003 | 0.25% |
| 3 (Siltstone) | 0.0807 | 0.0793 | 0.0014 | 1.73% |
| 4 (Overbank Mudstone) | 0.3581 | 0.3446 | 0.0135 | 3.92% |

- **Bulk Net-to-Gross (Sand + Splay)**:
  - Theoretical Stationary: **53.07%**
  - Empirical Observed: **54.57%**
  - Relative Difference: **2.74%** (< 3% diagnostic agreement)

### 4.5 Comparison Against Obsolete 3-Log Results (Sensitivity Analysis)
- **Stationary Net-to-Gross**: 53.07% (11-log) vs. 52.51% (3-log obsolete) [Difference: **+0.56%**]
- **Empirical Net-to-Gross**: 54.57% (11-log) vs. 55.06% (3-log obsolete) [Difference: **-0.50%**]
- **Regular Matrix Frobenius Difference $\|P_{11} - P_3\|_F$**: **0.3120**
- **Embedded Matrix Frobenius Difference $\|P_{11} - P_3\|_F$**: **0.6158**
- **Sparsity Constraints**: The obsolete 3-log subset contained only 50 boundary transitions (with sparse rows such as 7 coal and 8 siltstone boundaries, where a single count shifted transition probabilities by 12.5% to 14.3%). Expanding to the 11-log dataset (201 boundary crossings) regularizes the transition matrix while preserving the 100% upward fining invariant.
