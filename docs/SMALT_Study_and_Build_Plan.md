# **Subsurface Stratigraphic Modeling & Active Learning Toolkit (SMALT)** 

### _A Study-and-Build Plan_ 

**Prepared for:** Mohd Fahad • IIT Kanpur 

**Target Output:** Verified Python Repository, Active Learning Core-Optimization Engine, and Publication-Ready Manuscript Draft 

## **The Goal, in One Paragraph** 

Characterizing subsurface sedimentary architecture from sparse borehole logs is a fundamental challenge in petroleum geology and carbon storage. While deterministic models over-smooth rock boundaries, Active Learning paired with Stochastic Geostatistical Simulation can identify high-uncertainty stratigraphic transitions (such as channel pinch-outs and caprock seals) with minimal coring cost. Before deploying complex machine learning pipelines, we require a trustworthy, field-validated Python toolkit grounded in the fluvial-deltaic outcrop data of the Blackhawk Formation (Sahoo et al., 2016). This plan builds that toolkit phase-by-phase, converting raw lithologs into a 2D/3D conditional Markov simulator, an active margin sampling engine, and a deep-learning comparison framework. 

## **The Three Gates (Every 10-Day Phase)** 

- Read & Prove It: A written deliverable or mathematical derivation in docs/notes/ proving the underlying sedimentological or machine learning paper was thoroughly digested. 

- Build It: Concrete, documented Python modules with paper equations cited directly in docstrings. 

- Done When: Objective verification — a test passes, a synthetic benchmark is recovered, or a 2D cross-section matches known outcrop ground truth within specified tolerances. 

## **The Repository Architecture** 



<!-- Start of picture text -->
smalt/<br>├── data/<br>│   ├── raw_lithologs/       # Digitized 11 Wasatch Plateau outcrop sections (.csv)<br>│   ├── processed/           # Standardized facies, depth-normalized tables<br>│   └── synthetic_grids/     # 2D/3D geostatistically generated realizations<br>├── geostat/<br>│   ├── markov.py            # 1D/2D vertical & horizontal transition probability matrices<br>│   ├── object_sim.py        # Fluvial geometry generator (W/T=35, thickness bounds)<br>│   └── conditional_grid.py  # Inter-well 2D conditional spatial field generator<br>├── ml/<br>│   ├── active_sampling.py   # Margin sampling uncertainty & querying engine<br>│   ├── xgboost_model.py     # Fast tabular spatial classifier<br>│   └── transformer_bench.py # KD-Tree + Transformer baseline (Hang et al., 2025)<br>├── viz/<br>│   ├── cross_sections.py    # 2D facies panels & channel geometry plots<br>│   └── entropy_maps.py      # Normalized Information Entropy uncertainty maps<br>├── tests/<br>│   ├── test_markov.py       # Transition matrix row-stochasticity tests<br>│   └── test_active_loop.py  # Synthetic core budget efficiency benchmarks<br>├── docs/<br>│   ├── notes/               # Paper derivations and literature summaries<br>│   └── HANDOFF.md           # Running 1-page progress & block status<br><!-- End of picture text -->

```
└── run_smalt_pipeline.py    # Master end-to-end execution script
```

## **Cross-Cutting Discipline** 

- 10-Day Sprint Cycle: Each phase is executed over 10 days. 

- Docstring Citations: Every scientific function must cite the paper name and equation (e.g., Sahoo et al., 2016, Eq. 3 or Hang et al., 2025, Eq. 14). 

- Ground-Truth Verification: No module is merged until verified against synthetic benchmarks or the 11 real outcrop lithologs. 

- Living Handoff: HANDOFF.md is updated every sprint listing completed gates, verified metrics, and active blockers. 

## **Sprint Roadmap & Phase Breakdown** 

```
TRACK A: Data & Geostatistical Foundations
```

```
  Phase 0 (Days 1–10)   ──> Data Schema & 11-Litholog Standardization
  Phase 1 (Days 11–20)  ──> Markov Chain Transition Matrices & Vertical Succession
```

```
TRACK B: 2D/3D Stochastic Simulator Engine
  Phase 2 (Days 21–30)  ──> 2D Object-Based Fluvial Generator (W/T = 35)
  Phase 3 (Days 31–40)  ──> Inter-Well Conditional Grid Simulation & Synthetic Logs
```

```
TRACK C: Active Learning & Classification Engine
  Phase 4 (Days 41–50)  ──> XGBoost Baseline & Margin Sampling Uncertainty
  Phase 5 (Days 51–60)  ──> Active Core-Budget Optimization Pipeline
  Phase 6 (Days 61–70)  ──> Deep Learning Benchmark (Transformer vs. XGBoost)
```

```
TRACK D: Validation, Entropy Mapping & Manuscript
  Phase 7 (Days 71–80)  ──> Normalized Information Entropy & Blind Outcrop Validation
  Phase 8 (Days 81–90)  ──> End-to-End Capstone Execution & Paper Draft
```

## **Track A: Data & Geostatistical Foundations** 

**<mark>Phase 0 | Data Schema & 11-Litholog Standardization</mark>** _<mark>(Days 1–10)</mark>_ 

##### **Read & Prove It** 

Summarize data preparation and magnitude/facies encoding standards in docs/notes/phase0_schema.md. 

##### **Build It** 

- data/loader.py: Typed Pandas pipeline loading all 11 digitized lithologs into a uniform schema (litholog_id, depth_m, facies_code, facies_name, strike_pos_m). 

- data/clean.py: Depth normalization, zero-thickness cleaning, and synthetic Gamma Ray proxy generation (~30 API for Sandstone, ~100 API for Mudstone, ~15 API for Coal). 

##### **Done When** 

All 11 lithologs load cleanly into the standardized schema, pass unit tests for monotonic depth order, and survive CSV/Parquet round-trips. 

#### **<mark>Phase 1 | Markov Chain Transition Matrices & Vertical Succession</mark>** _<mark>(Days 11–20)</mark>_ 

##### **Read & Prove It** 

Derive the 1D embedded Markov transition probability matrix P(Sₜ | Sₜ₋₁) in docs/notes/phase1_markov.md, citing Krumbein & Dacey (1969) and Doveton (1971). 

##### **Build It** 

- geostat/markov.py: Compute vertical facies transition matrices P, tally upward fining/coarsening counts, and calculate stationary facies distributions across the 11 lithologs. 

- viz/markov_viz.py: Heatmap renderer for transition matrices and vertical succession state diagrams. 

##### **Done When** 

All rows in P sum strictly to 1.0000 (row-stochasticity unit test passes), and fining-upward sequences (Channel Sandstone → Siltstone → Mudstone) correctly dominate transition probabilities. 

## **Track B: 2D/3D Stochastic Simulator Engine** 

#### **<mark>Phase 2 | 2D Object-Based Fluvial Generator (W/T = 35)</mark>** _<mark>(Days 21–30)</mark>_ 

##### **Read & Prove It** 

Document fluvial geometry modeling rules in docs/notes/phase2_fluvial_priors.md using Sahoo et al. (2016) parameters: channel sandstone aspect ratio W/T = 35, mean thickness ~5.8 m, and splay width ranges (10–130 m). 

##### **Build It** 

- geostat/object_sim.py: Object-based body generator creating channel sandstones, crevasse splays, and overbank mudstones conditioned on aspect ratios and target Net-to-Gross (17%–46%). 

**Done When** 

100 unconditioned synthetic channel realizations produce mean width-to-thickness ratios within 35 ± 2 and Net-toGross values strictly within the 17%–46% envelope. 

#### **<mark>Phase 3 | Inter-Well Conditional Grid Simulation</mark>** _<mark>(Days 31–40)</mark>_ 

##### **Read & Prove It** 

Write the 2D conditional Markov simulation framework in docs/notes/phase3_conditional_sim.md, citing Elfeki & Dekking (2001) Coupled Markov Chains (CMC) and Carle & Fogg (1996) T-PROGS. 

##### **Build It** 

- geostat/conditional_grid.py: Generate 2D cross-sectional grids (500 m × 100 m at 1 m resolution) anchored to known vertical lithologs as hard constraints, filling intermediate space via distance-weighted spatial transition probabilities. 

##### **Done When** 

Simulated 2D grid cells matching anchor well coordinates are 100% identical to the input lithologs, while inter-well channels pinch out according to W/T geometric constraints. 

## **Track C: Active Learning & Classification Engine** 

**<mark>Phase 4 | XGBoost Baseline & Margin Sampling Uncertainty</mark>** _<mark>(Days 41–50)</mark>_ 

##### **Read & Prove It** 

Derive the Margin Sampling uncertainty metric in docs/notes/phase4_margin_sampling.md: Margin = P(ŷ₁ | x) − P(ŷ₂ | x), where ŷ₁ and ŷ₂ are the top two predicted facies classes. 

##### **Build It** 

- ml/xgboost_model.py: Train a baseline XGBoost classifier on sparse 2D grid spatial coordinates (X, Z, Gamma Ray, prev_facies). 

- ml/uncertainty.py: Compute probability margin heatmaps across the entire 2D cross-section. 

##### **Done When** 

Baseline XGBoost trains on a 10% spatial sample in under 5 seconds and outputs valid prediction probability distributions across all facies classes. 

#### **<mark>Phase 5 | Active Core-Budget Optimization Pipeline</mark>** _<mark>(Days 51–60)</mark>_ 

##### **Read & Prove It** 

Outline the Active Learning core-location selection algorithm in docs/notes/phase5_active_learning.md. 

##### **Build It** 

- ml/active_sampling.py: Iteratively query the top-N spatial coordinates with the smallest prediction margin (highest ambiguity), append queried true labels to the training set, and retrain XGBoost. 

##### **Done When** 

Active Margin Sampling achieves >90% classification accuracy using 40%–50% fewer labeled cells compared to standard random spatial sampling on synthetic 2D grid benchmarks. 

#### **<mark>Phase 6 | Deep Learning Benchmark: Transformer vs. XGBoost</mark>** _<mark>(Days 61–70)</mark>_ 

##### **Read & Prove It** 

Summarize the KD-Tree + Transformer architecture in docs/notes/phase6_transformer.md, citing Hang et al. (2025). 

##### **Build It** 

- ml/transformer_bench.py: Implement a PyTorch KD-Tree + Transformer Encoder benchmark as described by Hang et al. (2025). 

- Compare XGBoost + Active Learning against the Transformer model under sparse (3–5 well) vs. dense (10+ well) regimes. 

##### **Done When** 

Benchmarks quantitatively confirm that XGBoost + Active Learning outperforms Transformers under ultra-sparse seed conditions (<5 wells), while Transformers excel when pre-trained on 100+ synthetic simulator outputs. 

## **Track D: Validation, Entropy Mapping & Publication** 

**<mark>Phase 7 | Normalized Information Entropy & Blind Outcrop Validation</mark>** _<mark>(Days 71–80)</mark>_ 

##### **Read & Prove It** 

Derive Normalized Information Entropy H(X) in docs/notes/phase7_entropy.md, citing Hang et al. (2025, Eq. 14): H(X) = − [Σ p(xᵢ) ln(p(xᵢ))] / ln(S), summed over i = 1 to S. 

##### **Build It** 

- viz/entropy_maps.py: Generate 2D/3D spatial Information Entropy heatmaps to highlight caprock seal boundaries and channel margins. 

- Run blind cross-validation: Train on 7 lithologs, withhold 4 real outcrop lithologs as an un-touched test set, and evaluate inter-well prediction accuracy. 

##### **Done When** 

Information Entropy correctly spikes (H(X) > 0.7) along channel boundaries and mudstone seal margins, and blind outcrop test accuracy exceeds 85%. 

#### **<mark>Phase 8 | Capstone Pipeline Integration & Paper Draft</mark>** _<mark>(Days 81–90)</mark>_ 

##### **Build It** 

- run_smalt_pipeline.py: Single command that loads lithologs → builds Markov matrices → runs 2D conditional simulation → executes Active Learning coring loops → evaluates blind outcrop accuracy → renders publication-ready figures. 

##### **Done When** 

The master script runs end-to-end without errors, producing high-resolution figures for: 2D Conditional Facies Cross-Section; Active Learning Query Progression vs. Coring Budget Saved; Information Entropy Uncertainty Map; Blind Outcrop Validation Confusion Matrix. 

