# SMALT Professor Meeting Notes: Scientific Defensibility & Oral Defense Guide

**Project:** Subsurface Stratigraphic Modeling & Active Learning Toolkit (SMALT)  
**Document ID:** `SMALT-NOTES-PROF-001`  
**Target Audience:** Undergraduate Earth Sciences Student preparing for faculty meetings  
**Design Intent:** Plain-language oral explanations with full mathematical and geological rigor  
**BibTeX Reference:** [`docs/REFERENCES.bib`](file:///d:/Lithology-reconstruction-using-XGB/docs/REFERENCES.bib)  
**Methods Cross-Reference:** [`docs/SCIENTIFIC_METHODS_REGISTER.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/SCIENTIFIC_METHODS_REGISTER.md)  
**Figure Cross-Reference:** [`docs/FIGURE_PROVENANCE_REGISTER.md`](file:///d:/Lithology-reconstruction-using-XGB/docs/FIGURE_PROVENANCE_REGISTER.md)

---

## Meeting Strategy: The "Golden Rules" of Academic Defense

1. **Never bluff or invent data.** If a parameter is an assumption, state it proudly: *"That is an explicit modeling assumption. We tested its sensitivity, but our current lithologs do not uniquely constrain it."*
2. **Distinguish empirical observations from literature priors.** What we measured in the 12 logs is empirical; $W/T \approx 35$ is a literature prior from Sahoo et al. (2016).
3. **Own our negative results.** Our discovery that deterministic spatial interpolation collapses because well spacing exceeds channel width is a major scientific finding, not a failure.
4. **Never mention code syntax in a geological discussion.** Explain the physical process, the boundary conditions, and the geological implications.

---

## Question 1: Data Provenance & Standardization
### *"Where did these 12 lithologs come from, and why are certain sections excluded from spatial modeling?"*

#### 1. What problem it addresses
Field geological data come from different sources, have inconsistent formatting, and often contain subtle digitizing or spatial errors that corrupt mathematical models if not standardized.

#### 2. Why this method was selected
We established a strict 6-column unified schema and an automated Quality Control (QC) pipeline to enforce depth monotonicity ($z_{t+1} > z_t$), positive interval thicknesses, and provenance tracking.

#### 3. Basic intuition and essential principles
- All 12 lithologs were manually digitized from published figures and core photographs.
- Outcrop log 11 was benchmarked against automated optical extraction, matching 73 of 78 meters (93.59% agreement).
- The 12 logs span 1033.0 meters of stratigraphy across 296 distinct beds.

#### 4. Inputs and outputs
- *Inputs:* Raw digitized interval tables from outcrop cliff plates (L1-L11) and drill core EM-137C (L12).
- *Outputs:* Standardized CSV/Parquet files and the project provenance manifest ([`data/provenance_manifest.json`](file:///d:/Lithology-reconstruction-using-XGB/data/provenance_manifest.json)).

#### 5. Supporting literature
Sahoo et al. (2016), Marine and Petroleum Geology; Sanchez & Brown (1986).

#### 6. Differences from cited studies
Sahoo et al. analyzed photomosaics and measured sections visually; SMALT digitizes them into discrete mathematical intervals for algorithmic modeling.

#### 7. Evidence supporting current implementation
All 12 logs pass automated monotonicity and column validation tests in [`tests/test_phase0_loader.py`](file:///d:/Lithology-reconstruction-using-XGB/tests/test_phase0_loader.py).

#### 8. What remains uncertain
- Lithologs 2 through 10 have not undergone independent optical benchmarking.
- Core EM-137C total penetration was 242 meters, but only the upper 111 meters was digitized.

#### 9. Alternative methods considered
- *Alternative:* Merging thin beds into coarse facies. *Rejected:* Destroys thin coal markers and splay heterolithics.
- *Alternative:* Synthesizing missing depths. *Rejected:* Violates data integrity rules.

#### 10. Likely Professor Questions and Defensible Oral Answers
- **Q1: "Did you use an AI tool to hallucinate these lithologs?"**  
  *Your Answer:* "No, Professor. An earlier sprint had a confusing label called 'ai_reconstructed', but our Sprint E audit confirmed that all 12 lithologs were manually digitized by human measurement from the published cliff plates in Sahoo et al. (2016) and the EM-137C core photo. They are real digitized observations."
- **Q2: "Why did you exclude Litholog 1 from your spatial maps?"**  
  *Your Answer:* "Litholog 1 is a valid 93-meter measured section, but its physical $(X, Y)$ coordinates were omitted from the location register. Rather than inventing arbitrary coordinates, we quarantined L1: it is used for 1D vertical Markov tallies, but strictly excluded from spatial distance calculations."
- **Q3: "Does your database include the complete EM-137C core?"**  
  *Your Answer:* "The original drill core penetrated 242 meters, but our digitized database covers strictly the top 111 meters across 63 digitized beds. We report this boundary explicitly in our manifest and never fabricate data for the lower 131 meters."

---

## Question 2: Vertical Succession & Walther's Law
### *"How do you mathematically prove and model vertical fining-upward cycles?"*

#### 1. What problem it addresses
Visual litholog inspection suggests fining-upward tendencies, but subjective interpretations cannot be fed into numerical simulation engines. We need a mathematical model of vertical succession.

#### 2. Why this method was selected
A first-order Markov chain tallied in upward stratigraphic order quantifies the transition probability from an underlying bed to an overlying bed, testing Walther's Law of facies succession.

#### 3. Basic intuition and essential formula
- Deposition proceeds upward through time.
- Transition probability:
  $$P_{ij} = \frac{n_{ij}}{\sum_k n_{ik}}$$
- To isolate true facies boundaries from bed thickness, we use an **Embedded Markov Chain** where self-transitions on the diagonal are set to zero ($P_{ii} = 0$).
- Asymmetry matrix:
  $$\mathbf{A} = \mathbf{P} - \mathbf{P}^T$$
  Positive values ($A_{ij} > 0$) reveal preferred upward transitions.

#### 4. Inputs and outputs
- *Inputs:* Discrete vertical facies sequences sorted from base to top.
- *Outputs:* Transition count matrices, transition probability matrices, stationary distributions, and asymmetry heatmaps.

#### 5. Supporting literature
- Abasaghi et al. (2026), *Basin Research* (Quasi-independence model, Chi-square testing, directional entropy).
- Powers & Easterling (1982), *Journal of Sedimentary Petrology* (Embedded Markov chain counting).
- Hattori (1976), *Mathematical Geology* (Pre- and post-depositional entropy).

#### 6. Differences from cited studies
- Abasaghi et al. studied a carbonate ramp forced by Milankovitch orbital cycles. SMALT studies a siliciclastic coastal plain where cyclicity is driven by channel avulsion and compaction.

#### 7. Evidence supporting current implementation
In our Upstream outcrop data, when a Channel Sandstone transitions across a bed boundary, it transitions **100% of the time into finer facies** ($40\%$ splay ripples, $30\%$ coal, $20\%$ siltstone, $10\%$ overbank mud). It never transitions directly into another coarse channel sand without an intervening contact.

#### 8. What remains uncertain
A first-order Markov chain assumes memoryless transitions (one step back). It does not capture multi-cycle sequence-stratigraphic packages.

#### 9. Alternative methods considered
- *Alternative:* Symmetric variograms. *Rejected:* Variograms are non-directional; they cannot distinguish fining upward from coarsening upward.

#### 10. Likely Professor Questions and Defensible Oral Answers
- **Q1: "Why do you tally transitions from bottom to top rather than top to bottom?"**  
  *Your Answer:* "Because sediments accumulate from bottom to top over time. Tallying upward follows the physical arrow of geological time, allowing us to quantify upward-fining channel abandonment directly."
- **Q2: "What does the stationary distribution prove?"**  
  *Your Answer:* "The stationary distribution solves the eigenvector equation $\boldsymbol{\pi}\mathbf{P} = \boldsymbol{\pi}$. It gives the theoretical long-term facies proportions. In our upstream data, the theoretical stationary pure sand proportion is 50.5%, matching our observed field measurements within 0.1%, which proves our transition matrix preserves the macroscopic sand budget."
- **Q3: "Can you use this 1D Markov chain to predict facies laterally between wells?"**  
  *Your Answer:* "No, Professor, and we explicitly avoid that mistake. A 1D vertical Markov chain reflects vertical Waltherian stacking. It cannot be flipped horizontally to predict lateral facies distributions because lateral channel dimensions are governed by geomorphic width-to-thickness ratios, not vertical bed transitions."

---

## Question 3: Bed Thickness Modeling & The "Sticky State" Trap
### *"How do you model bed thicknesses without causing simulations to trap in thick states?"*

#### 1. What problem it addresses
In computer simulation, if you model stratigraphy using fixed vertical grid steps (e.g. 0.5 m cells), thick beds produce huge self-transition probabilities ($P_{ii} > 0.90$). A stochastic simulator entering that state will "stick" there, generating unphysically thick megabeds.

#### 2. Why this method was selected
Following Dindarloo et al. (2015), we decouple facies transition sequencing from bed thickness modeling using a **hybrid Markov-Monte Carlo approach**: the embedded Markov chain picks the next facies, and an independent continuous distribution samples the bed thickness.

#### 3. Basic intuition and essential formula
- Channel sandstone bed thicknesses follow a Log-Normal distribution:
  $$T \sim \text{LogNormal}(\mu, \sigma), \quad \hat{\mu} = 1.7870, \quad \hat{\sigma} = 0.8823$$
- Estimated via Maximum Likelihood Estimation (MLE) on 34 Upstream channel sandstone beds.
- Theoretical median: $e^{1.7870} = 5.97\text{ m}$ (Empirical median = $5.50\text{ m}$).
- Theoretical mean: $8.81\text{ m}$ (Empirical mean = $8.50\text{ m}$).

#### 4. Inputs and outputs
- *Inputs:* 34 continuous measured bed thicknesses from the Upstream litholog group.
- *Outputs:* Calibrated parametric thickness generator in [`smalt/forward/config.py`](file:///d:/Lithology-reconstruction-using-XGB/smalt/forward/config.py).

#### 5. Supporting literature
- Dindarloo et al. (2015), *Sedimentary Geology* (Proved that regular chains cause pathological over-persistence; established hybrid Markov-Monte Carlo framework).
- Sahoo et al. (2016), Figure 10A (Empirical right-skewed channel thickness histogram).

#### 6. Differences from cited studies
Dindarloo et al. studied coal lithotypes in a 2.64 m core interval. SMALT applies the hybrid principle to fluvial channel sandstones and crevasse splays across a 50-meter stratigraphic column.

#### 7. Evidence supporting current implementation
Resampling checks confirm that the fitted log-normal model reproduces empirical medians (5.97m vs 5.50m) and 90th percentiles (18.5m vs 17.4m).

#### 8. What remains uncertain (The Clustered Sampling Caveat)
An unadjusted Kolmogorov-Smirnov test gave $D=0.1059, p=0.8023$. However, because $\mu$ and $\sigma$ were fitted on the same 34 beds, and because those 34 beds are clustered across only 8 lithologs, observations are not strictly independent. The log-normal model is a *defensible engineering parameterization*, not an established universal law.

#### 9. Alternative methods considered
- *Alternative:* Exponential distribution. *Rejected:* Underestimates the frequency of thick 5-15 m amalgamated channel complexes.
- *Alternative:* Non-parametric empirical sampling. *Rejected:* Cannot simulate realistic bed thicknesses outside the exact observed numbers.

#### 10. Likely Professor Questions and Defensible Oral Answers
- **Q1: "Why not just use the diagonal probabilities in a regular Markov chain to model thickness?"**  
  *Your Answer:* "Because as Dindarloo et al. (2015) demonstrated in Sedimentary Geology, regular chains suffer from diagonal dominance. When the diagonal probability exceeds 90%, stochastic simulation creates runaway 'sticky states' that produce artificially massive beds. Decoupling the transition sequence from continuous thickness sampling eliminates that pathology."
- **Q2: "Your KS test p-value is 0.80. Does that prove bed thicknesses are log-normal?"**  
  *Your Answer:* "No, Professor, and we do not claim that. That KS test is unadjusted: the parameters were estimated from the same 34 beds, which inflates the p-value. Furthermore, because those 34 beds come from eight wells, they are stratigraphically clustered. The log-normal model is a reasonable, well-calibrated approximation, but it is not a proven physical law."

---

## Question 4: Spatial Correlation vs Sparsity (The Negative Result)
### *"Can we interpolate facies between these wells using spatial geostatistics or machine learning?"*

#### 1. What problem it addresses
Evaluating whether the existing 12 lithologs contain sufficient spatial correlation to enable direct inter-well subsurface interpolation.

#### 2. Why this method was selected
We applied continuous transition probability geostatistics (Carle & Fogg, 1996; Weissmann & Fogg, 1999) and horizontal transition decay calibration to 4,410 empirical borehole pairs across common-zero elevation slices.

#### 3. Basic intuition and essential formula
- Continuous spatial transition probability:
  $$\mathbf{T}(h) = \exp(\mathbf{R} h)$$
  where diagonal rate $r_{jj} = -1 / L_j$ is inverse mean facies length.
- Spatial correlation length: Channel width $W \approx 287\text{ m}$.
- Field well spacing: $h \ge 420.0\text{ m}$ (mean $> 800\text{ m}$).

#### 4. Inputs and outputs
- *Inputs:* 11 spatial lithologs aligned to the Axel Anderson coal datum.
- *Outputs:* Empirical horizontal transition curves and Leave-One-Log-Out (LOLO) cross-validation scores.

#### 5. Supporting literature
- Weissmann & Fogg (1999), *Journal of Hydrology* (T-PROGS transition probability framework).
- Carle & Fogg (1996, 1997), *Mathematical Geology*.

#### 6. Differences from cited studies
Weissmann and Fogg had 60 closely spaced test holes and regional soil maps, allowing them to capture continuous lateral transitions. In SMALT, our wells are spaced farther apart than the channel bodies themselves.

#### 7. Evidence supporting current implementation (Category C Proof)
- Because $h_{\text{min}} = 420\text{ m} > W = 287\text{ m}$, pairs of wells almost never penetrate the same channel body.
- Horizontal transition probabilities decay to background proportions by 420 meters.
- Under Maximum A Posteriori (MAP) estimation, conditioned geostatistical models collapse to the background floodplain mudstone mode everywhere between wells ($>90\%$ mudstone), achieving ~44% accuracy.

#### 8. What remains uncertain
True 3D channel sinuosity and exact lateral body terminations cannot be resolved from 1D borehole columns alone.

#### 9. Alternative methods considered
- *Alternative:* Training XGBoost or Neural Networks on spatial coordinates. *Rejected:* With only 11 wells spaced 800 m apart, ML models overfit or memorize training points, producing unphysical spatial artifacts.

#### 10. Likely Professor Questions and Defensible Oral Answers
- **Q1: "Why did your spatial interpolation models fail to predict channel sands between wells?"**  
  *Your Answer:* "It was not a mathematical failure; it was a physical scale limitation. Sahoo et al. (2016) showed that Blackhawk channel belts average 287 meters wide. Our closest wells are 420 meters apart. When well spacing is wider than the geobodies, there is zero correlation between adjacent wells. Any deterministic model mathematically collapses to background mudstone."
- **Q2: "What is the scientific value of this negative result?"**  
  *Your Answer:* "It establishes an essential scientific boundary: it proves that machine learning interpolation cannot 'connect the dots' when well spacing exceeds facies dimensions. This redirected SMALT to process-based forward simulation and active learning to determine where to drill next."

---

## Question 5: Forward Geological Modeling Architecture
### *"How does the 2D forward simulator work, and what literature priors constrain it?"*

#### 1. What problem it addresses
Because sparse field wells cannot resolve inter-well architecture, we need a process-based forward model to create realistic synthetic geology with 100% known ground truth to benchmark coring strategies.

#### 2. Why this method was selected
A process-like 2D generator honors sedimentological rules: channel aspect ratios, erosional scours, flanking splays, and coal horizons.

#### 3. Basic intuition and essential formula
- Domain: $X = 2000\text{ m}, Z = 50\text{ m}, \Delta x = 10\text{ m}, \Delta z = 0.5\text{ m}$.
- Parabolic channel scour:
  $$z_{\text{top}}(x) = z_{\text{base}} + T \cdot \left(1 - \frac{4 (x - x_c)^2}{W^2}\right)$$
- Channel width scaled via literature prior:
  $$W = T \cdot (W/T), \quad W/T \sim \mathcal{N}(35.0, 5.0)$$
- Lower 20% partitioned into basal planar sandstone (`p_sand`).
- Younger channels deterministically erode older deposits (superposition).

#### 4. Inputs and outputs
- *Inputs:* Calibrated preset configurations (`upstream`, `downstream_outcrop`, `pooled`).
- *Outputs:* Complete 2D facies grids (20,000 cells) and cross-section visualizations.

#### 5. Supporting literature
- Sahoo et al. (2016): $W/T \approx 35$ (page 52); crevasse splay width 20-100 m (Figure 10D).
- Gibling (2006): Global compilation of fluvial channel aspect ratios.

#### 6. Differences from cited studies
Sahoo et al. observed static outcrop cliff faces. SMALT converts these observational priors into an algorithmic stochastic generator.

#### 7. Evidence supporting current implementation
Across 30 seeds, Upstream realizations yield realized pure sand $N/G = 47.55\% \pm 3.13\%$, mean channel thickness $6.12\text{ m}$, and mean channel width $286.8\text{ m}$, matching Sahoo et al.'s field observations.

#### 8. What remains uncertain
- 3D channel sinuosity and avulsion return periods are not modeled in 2D.
- Basal planar sand fraction (20%) is a calibrated modeling choice.

#### 9. Alternative methods considered
- *Alternative:* Rectangular channel boxes. *Rejected:* Unphysical flat bases; violates concave erosional scour geometry.

#### 10. Likely Professor Questions and Defensible Oral Answers
- **Q1: "Where did your width-to-thickness ratio of 35 come from?"**  
  *Your Answer:* "It comes directly from page 52 of Sahoo et al. (2016). The authors measured fully preserved channel sandbodies in this exact Wasatch Plateau outcrop and found $W/T \approx 35$, which agreed with LiDAR measurements by Rittersbacher et al. (2013) and Gibling's (2006) coastal plain compilation."
- **Q2: "Why did your simulator realize 38% channel sand when configured for 45%?"**  
  *Your Answer:* "Because of the stopping rule: the simulator stops placing channels when single-facies `sand` crosses $45\% - 8\% = 37\%$. When you add the 20% basal planar sand, the total pure channel sandstone is 47.55%, perfectly matching the geological target."

---

## Question 6: Sparse Reconstruction & Active Learning Benchmark
### *"How do you evaluate subsurface facies recovery without ground-truth leakage, and what does the benchmark prove?"*

#### 1. What problem it addresses
Evaluating whether baseline algorithms can recover unknown subsurface facies from sparse boreholes, establishing a benchmark for future Active Learning coring optimization.

#### 2. Why this method was selected
Controlled synthetic benchmarking with strict leakage safeguards: models receive only sparse 1D borehole observations; all borehole columns are masked out during scoring.

#### 3. Basic intuition and essential formula
- Tested spacings: 400 m (5 wells), 600 m (3 wells), 1000 m (2 wells).
- Evaluated on withheld inter-well cells:
  $$\text{Balanced Accuracy} = \frac{1}{K} \sum_{k=1}^K \text{Recall}_k$$
- Evaluated four baseline models: Nearest Well, Vertical Prior, Spatial 2D KNN ($k=5$), and Spatial Markov Decay.

#### 4. Inputs and outputs
- *Inputs:* `SyntheticBoreholeDataset` (only sampled borehole rows).
- *Outputs:* [`synthetic_reconstruction_benchmark_metrics.csv`](file:///d:/Lithology-reconstruction-using-XGB/sprints/audit_sprint_k/synthetic_reconstruction_benchmark_metrics.csv) (120 model evaluations).

#### 5. Supporting literature
He et al. (2019) (Ablation testing of well data); Weissmann & Fogg (1999).

#### 6. Differences from cited studies
He et al. had continuous downhole wireline petrophysical logs; SMALT benchmarks reconstruction using only discrete lithofacies from sparse lithologs.

#### 7. Evidence supporting current implementation
Across 120 evaluations:
- At 400m spacing: Nearest Well Balanced Acc = $52.61\%$, KNN = $42.46\%$.
- At 600m spacing: Nearest Well Balanced Acc = $46.74\%$, KNN = $36.98\%$.
- At 1000m spacing: Nearest Well Balanced Acc = $38.89\%$, KNN = $30.50\%$.

#### 8. What remains uncertain
Optimal non-uniform borehole placement (active acquisition) has not yet been executed; it is the scheduled priority for Sprint L.

#### 9. Alternative methods considered
- *Alternative:* Scoring models including borehole columns. *Rejected:* Gives artificial 100% scores at well locations, masking true inter-well error.

#### 10. Likely Professor Questions and Defensible Oral Answers
- **Q1: "Why does Nearest Well outperform 2D KNN?"**  
  *Your Answer:* "Because Nearest Well preserves vertical succession and bed thickness integrity. In contrast, 2D KNN introduces horizontal smearing, blurring sharp 0.5-meter coal and splay beds into neighboring cells and destroying thin-bed geometry."
- **Q2: "How did you prevent ground-truth leakage in the benchmark?"**  
  *Your Answer:* "First, predictive models receive only the extracted borehole columns; the full grid is completely hidden. Second, during scoring, all columns containing boreholes are masked out so the model is evaluated strictly on unobserved inter-well cells."
- **Q3: "What is the next step for SMALT?"**  
  *Your Answer:* "Now that we have a verified 2D forward simulator and benchmark harness, the next step is Active Learning: using predictive entropy reduction to optimize where to drill the next borehole to locate channel sandbodies with minimum drilling cost."

---

## Question 7: Audit of the Four Professor-Provided Papers
### *"How do the four papers relate to SMALT, and where do their methods transfer or stop?"*

| Paper | Primary Research Question & Setting | Data Used | Relevant Method Motivations | Non-Transferable Scope Limits |
| :--- | :--- | :--- | :--- | :--- |
| **Abasaghi et al. (2026)** | Cyclostratigraphy and 1D Markov cyclicity in Oligo-Miocene Asmari carbonate ramp, Zagros Basin. | Spectral gamma-ray logs and thin-sections from 6 deep wells. | Quasi-independence null model, Chi-square testing, $Z_{ij}$ transition paths, Hattori entropy. | Analyzes 1D vertical wellbores on a carbonate ramp; does NOT construct 2D/3D lateral spatial reconstructions. |
| **Weissmann & Fogg (1999)** | Multi-scale hydrofacies modeling in Pleistocene Kings River Alluvial Fan, California. | 1000m core across 60 shallow holes, 7 well logs, soil survey maps. | Continuous transition probability $\mathbf{T}(h) = \exp(\mathbf{R}h)$, sequence-stratigraphic non-stationarity. | Exceptional data density (60 core holes); alluvial fan setting. Cannot calibrate horizontal rates from 12 sparse lithologs. |
| **Dindarloo et al. (2015)** | Microlithotype and lithotype modeling in Pennsylvanian Kentucky coal beds. | 2.64m megascopic coal profile and 256 polished block intervals. | Decouples transition sequencing (embedded Markov chain) from continuous thickness sampling (Monte Carlo PDF). | Studied 1D vertical peat swamp profiles; contains no lateral sandbody dimensions, $W/T$ ratios, or spatial geostatistics. |
| **He et al. (2019)** | Neural network facies prediction from wireline logs in Jurassic Surat Basin, Australia. | 2000m core across 8 wells; 6 physical wireline logs across 31 wells at 0.15m resolution. | Markov chain analysis constrains vertical facies associations; ablated feature benchmarks. | Had authentic 6-log petrophysical suite. SMALT has NO physical wireline logs; using synthetic GR would be circular leakage. |

#### Likely Professor Questions on Literature Audit
- **Q1: "Why can't you use He et al.'s neural network approach on SMALT's dataset?"**  
  *Your Answer:* "He et al. (2019) had an authentic petrophysical suite of six continuous wireline logs (density, sonic, resistivity, neutron, photoelectric, gamma ray) across 31 wells. In SMALT, we have manually digitized lithologs with zero physical wireline logs. Our gamma-ray column is synthetic, imputed directly from facies labels. Training a neural network on synthetic gamma-ray to predict facies would be 100% circular ground-truth leakage."
- **Q2: "Does Abasaghi et al. (2026) validate using 1D Markov chains to reconstruct cross-sections?"**  
  *Your Answer:* "No, Professor. Abasaghi et al. validate 1D vertical Markov succession and entropy analysis in vertical wellbores. They do not simulate 2D or 3D lateral stratigraphy. Citing them as validation for lateral inter-well prediction would be a scientific overreach."
- **Q3: "How did Dindarloo et al. (2015) influence your forward model?"**  
  *Your Answer:* "Dindarloo et al. proved that fixed-grid Markov chains suffer from diagonal dominance and produce unrealistically thick 'sticky' beds. They introduced a hybrid solution: using an embedded chain for transition sequence and an independent continuous distribution for bed thickness. We adopted that exact hybrid principle in our simulator."
