"""
SMALT Sprint K Test Suite: Geological Forward Simulator & Synthetic Borehole Sampling.

Verifies:
1. Canonical 6-state facies schema validity and non-merging.
2. Correct grid shape and bottom-to-top stratigraphic orientation.
3. 100% deterministic seed reproducibility.
4. Object geometry constraints, boundary safeguards, and truncation rules.
5. Positive, bounded bed thicknesses and run lengths.
6. Valid borehole sampling positions, coordinate alignment, and extracted facies.
7. Strict ground-truth separation, hidden-grid protection, and leakage prevention.
8. Correct facies proportions, Net-to-Gross definitions, and transition diagnostics.
9. Configuration validation, preset integrity, and clear exception handling.
10. Regression of existing Markov chain and conditioned realization modules.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest

from smalt.descriptive.analyzer import CANONICAL_FACIES_SCHEMA, normalize_facies_label
from smalt.geostat.markov import StratigraphicMarkovChain
from smalt.forward import (
    GridConfig,
    ChannelBeltConfig,
    CrevasseSplayConfig,
    CoalHorizonConfig,
    BackgroundConfig,
    ForwardSimulatorConfig,
    get_preset_config,
    SyntheticRealization2D,
    FluvialForwardSimulator2D,
    CANONICAL_CODE_TO_NAME,
    SyntheticBoreholeDataset,
    SyntheticBoreholeSampler,
)


# -----------------------------------------------------------------------------
# 1. Canonical Facies-Code Validity
# -----------------------------------------------------------------------------
def test_01_canonical_facies_schema_validity():
    """Test 1: All six canonical facies codes (0 to 5) are distinct and validated."""
    expected_mapping = {
        0: "sand",
        1: "p_sand",
        2: "ripples",
        3: "carbon_mud",
        4: "coal",
        5: "mud",
    }
    assert len(CANONICAL_CODE_TO_NAME) == 6
    for code, name in expected_mapping.items():
        assert CANONICAL_CODE_TO_NAME[code] == name
        assert CANONICAL_FACIES_SCHEMA[name]["code"] == code

    # Test realization produces only valid canonical facies codes
    sim = FluvialForwardSimulator2D()
    realization = sim.generate(seed=101)
    unique_codes = np.unique(realization.grid)
    for c in unique_codes:
        assert c in expected_mapping, f"Encountered non-canonical facies code {c}"


# -----------------------------------------------------------------------------
# 2. Correct Grid Shape and Stratigraphic Orientation
# -----------------------------------------------------------------------------
def test_02_grid_shape_and_stratigraphic_orientation():
    """Test 2: Grid dimensions and vertical stratigraphic orientation (bottom-to-top)."""
    grid_cfg = GridConfig(x_extent_m=1000.0, z_extent_m=40.0, dx_m=20.0, dz_m=1.0)
    cfg = ForwardSimulatorConfig(grid=grid_cfg)
    sim = FluvialForwardSimulator2D(cfg)
    realization = sim.generate(seed=42)

    assert realization.grid.shape == (40, 50)
    assert realization.nx == 50
    assert realization.nz == 40
    assert len(realization.x_coords_m) == 50
    assert len(realization.z_coords_m) == 40

    # Stratigraphic orientation: row 0 is base (near z=0), row nz-1 is top (near z=z_extent)
    assert realization.z_coords_m[0] < realization.z_coords_m[-1]
    assert np.all(np.diff(realization.z_coords_m) > 0), "Stratigraphic z coordinates must strictly increase upward"
    assert realization.z_coords_m[0] == pytest.approx(0.5, abs=1e-3)
    assert realization.z_coords_m[-1] == pytest.approx(39.5, abs=1e-3)


# -----------------------------------------------------------------------------
# 3. Seed Reproducibility
# -----------------------------------------------------------------------------
def test_03_seed_reproducibility():
    """Test 3: Identical seed produces bitwise identical grid; different seeds differ."""
    sim = FluvialForwardSimulator2D(get_preset_config("upstream"))

    res1 = sim.generate(seed=777)
    res2 = sim.generate(seed=777)
    res3 = sim.generate(seed=888)

    # Identical seed -> 100% bitwise equality
    np.testing.assert_array_equal(res1.grid, res2.grid)
    assert res1.facies_proportions == res2.facies_proportions
    assert res1.realized_net_to_gross == res2.realized_net_to_gross

    # Different seed -> distinct realization
    assert not np.array_equal(res1.grid, res3.grid), "Different seeds should produce different realizations"


# -----------------------------------------------------------------------------
# 4. Valid Object Geometry and Grid Boundaries
# -----------------------------------------------------------------------------
def test_04_object_geometry_and_grid_boundaries():
    """Test 4: Placed objects respect domain boundaries and geometric rules."""
    sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
    res = sim.generate(seed=42)

    assert len(res.objects_metadata) > 0
    x_max = res.config.grid.x_extent_m
    z_max = res.config.grid.z_extent_m

    for obj in res.objects_metadata:
        if obj["object_type"] == "channel_belt":
            assert obj["thickness_m"] > 0.0
            assert obj["width_m"] > 0.0
            assert obj["wt_ratio"] > 0.0
            assert 0.0 <= obj["xc_m"] <= x_max
            assert 0.0 <= obj["z_base_m"] <= z_max
            assert obj["x_bounds"][0] >= 0.0
            assert obj["x_bounds"][1] <= x_max
            assert obj["z_bounds"][0] >= 0.0
            assert obj["z_bounds"][1] <= z_max
        elif obj["object_type"] == "coal_horizon":
            assert 0.0 <= obj["z_center_m"] <= z_max
            assert obj["thickness_m"] > 0.0


# -----------------------------------------------------------------------------
# 5. Positive, Valid Bed Thicknesses
# -----------------------------------------------------------------------------
def test_05_bed_thicknesses_positive_and_bounded():
    """Test 5: Sampled bed thicknesses in vertical columns are strictly positive and bounded."""
    sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
    res = sim.generate(seed=42)
    diag = res.vertical_diagnostics

    assert "mean_bed_thickness_m" in diag
    means = diag["mean_bed_thickness_m"]

    # Every observed facies has positive mean bed thickness
    for f_name, mean_t in means.items():
        if res.facies_proportions[f_name] > 0.0:
            assert mean_t >= res.config.grid.dz_m, f"Facies {f_name} mean bed thickness {mean_t} < dz"
            assert mean_t <= res.config.grid.z_extent_m, f"Facies {f_name} mean bed thickness {mean_t} > z_extent"


# -----------------------------------------------------------------------------
# 6. Valid Sampling Positions and Extracted Facies
# -----------------------------------------------------------------------------
def test_06_synthetic_sampling_positions_and_facies():
    """Test 6: Borehole sampler extracts exact grid columns at specified coordinates."""
    sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
    res = sim.generate(seed=42)
    sampler = SyntheticBoreholeSampler(res)

    dataset = sampler.sample_regular(spacing_m=500.0, offset_m=100.0)
    assert dataset.num_boreholes == 4  # 100, 600, 1100, 1600
    assert dataset.x_locations_m == [100.0, 600.0, 1100.0, 1600.0]
    assert dataset.num_observations == 4 * res.nz

    # Verify extracted facies matches grid array exactly
    for b_id, x_val in zip(dataset.borehole_ids, dataset.x_locations_m):
        col_df = dataset.get_borehole(b_id)
        assert len(col_df) == res.nz
        # Downward depth + strat_height = z_extent
        np.testing.assert_allclose(
            col_df["depth_m"] + col_df["strat_height_m"],
            res.config.grid.z_extent_m,
            atol=1e-5,
        )

        col_idx = int(np.clip(np.round(x_val / res.config.grid.dx_m), 0, res.nx - 1))
        expected_facies = res.grid[:, col_idx]
        actual_facies = col_df["facies_code"].to_numpy()
        np.testing.assert_array_equal(actual_facies, expected_facies)


# -----------------------------------------------------------------------------
# 7. Separation Between Observed Boreholes and Withheld Ground Truth
# -----------------------------------------------------------------------------
def test_07_ground_truth_separation_and_leakage_safeguards():
    """Test 7: Borehole dataset does not expose withheld grid, geobodies, or full stats."""
    sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
    res = sim.generate(seed=42)
    sampler = SyntheticBoreholeSampler(res)
    dataset = sampler.sample_regular(spacing_m=400.0)

    # Dataset container does NOT hold the 2D grid
    assert not hasattr(dataset, "grid")
    assert not hasattr(dataset, "objects_metadata")
    assert not hasattr(dataset, "realized_net_to_gross")

    # Dataframe contains only observation columns
    permitted_cols = {
        "litholog_id",
        "synthetic_x_m",
        "depth_m",
        "strat_height_m",
        "facies_code",
        "facies_name",
        "grid_col_idx",
        "grid_row_idx",
    }
    assert set(dataset.data.columns) == permitted_cols

    # Synthetic coordinates explicitly tagged and quarantined
    assert dataset.metadata["is_synthetic"] is True
    assert dataset.metadata["coordinate_system"] == "project_local_synthetic_quarantined"
    for b_id in dataset.borehole_ids:
        assert b_id.startswith("synth_well_")


# -----------------------------------------------------------------------------
# 8. Correct Facies Proportions and Transition Diagnostics
# -----------------------------------------------------------------------------
def test_08_facies_proportions_and_transition_diagnostics():
    """Test 8: Proportions sum to 1.0 and Net-to-Gross metrics obey strict hierarchy."""
    sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
    res = sim.generate(seed=42)

    total_prop = sum(res.facies_proportions.values())
    assert total_prop == pytest.approx(1.0, abs=1e-4)

    ntg = res.realized_net_to_gross
    assert 0.0 <= ntg["ntg_channel"] <= 1.0
    assert 0.0 <= ntg["ntg_pure"] <= 1.0
    assert 0.0 <= ntg["ntg_coarse"] <= 1.0
    # Net-to-Gross nesting: channel <= pure sand <= total coarse
    assert ntg["ntg_channel"] <= ntg["ntg_pure"] + 1e-9
    assert ntg["ntg_pure"] <= ntg["ntg_coarse"] + 1e-9

    # Transition count matrix dimensions: 6 x 6
    N_trans = np.array(res.vertical_diagnostics["transition_count_matrix"])
    assert N_trans.shape == (6, 6)
    assert np.all(N_trans >= 0)
    # Total vertical transitions across nx columns
    expected_transitions = res.nx * (res.nz - 1)
    assert np.sum(N_trans) == expected_transitions


# -----------------------------------------------------------------------------
# 9. Configuration Validation and Failure Handling
# -----------------------------------------------------------------------------
def test_09_configuration_validation_and_failure_handling():
    """Test 9: Invalid simulator configurations raise clear ValueErrors."""
    # Invalid grid extent
    with pytest.raises(ValueError, match="Grid extents must be positive"):
        GridConfig(x_extent_m=-100.0, z_extent_m=50.0).validate()

    # Invalid thickness bounds
    with pytest.raises(ValueError, match="Invalid thickness bounds"):
        ChannelBeltConfig(thickness_min_m=10.0, thickness_max_m=5.0).validate()

    # Invalid target Net-to-Gross
    with pytest.raises(ValueError, match="target_net_to_gross must be in"):
        ForwardSimulatorConfig(target_net_to_gross=1.5).validate()

    # Unknown preset
    with pytest.raises(ValueError, match="Unknown preset"):
        get_preset_config("nonexistent_preset")


# -----------------------------------------------------------------------------
# 10. Regression of Existing Markov and Conditioned Realization Modules
# -----------------------------------------------------------------------------
def test_10_regression_existing_markov_and_realization():
    """Test 10: Existing StratigraphicMarkovChain and provenance metadata remain intact."""
    # Existing StratigraphicMarkovChain
    dummy_df = pd.DataFrame({
        "litholog_id": ["well_1"] * 10,
        "depth_m": list(range(10, 0, -1)),
        "facies_code": [0, 0, 1, 1, 2, 3, 4, 5, 5, 5],
    })
    mc = StratigraphicMarkovChain(num_classes=6)
    mc.fit(dummy_df)
    assert mc.transition_matrix_.shape == (6, 6)
    np.testing.assert_allclose(mc.transition_matrix_.sum(axis=1), np.ones(6), atol=1e-6)

    # Provenance manifest retains manual digitization
    manifest_path = Path("data/provenance_manifest.json")
    assert manifest_path.exists()
    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    assert "All 12 lithologs were manually digitized" in manifest["metadata"]["provenance_policy"]
    l12 = manifest["lithologs"]["litholog12"]
    assert l12["digitized_intervals_count"] == 63
    assert l12["original_core_length_m"] == 242.0
    assert l12["independent_third_party_validation"] == "none (unvalidated)"


# -----------------------------------------------------------------------------
# 11. Multi-Seed Validation Pipeline & Metrics Integrity
# -----------------------------------------------------------------------------
def test_11_multi_seed_validation_metrics_integrity():
    """Test 11: Multi-seed validator runs and generates valid non-empty metrics."""
    from smalt.forward.validation import MultiSeedGeologicalValidator

    validator = MultiSeedGeologicalValidator()
    # Run 3 seeds across 2 presets for fast automated testing
    df_metrics = validator.run_multi_seed_campaign(presets=["upstream", "pooled"], num_seeds=3, base_seed=9500)
    assert len(df_metrics) == 6
    assert set(df_metrics["preset"].unique()) == {"upstream", "pooled"}

    # Proportions are non-negative and sum to approximately 1.0
    for _, row in df_metrics.iterrows():
        p_sum = sum(row[f"{f}_prop"] for f in ["sand", "p_sand", "ripples", "carbon_mud", "coal", "mud"])
        assert p_sum == pytest.approx(1.0, abs=1e-3)
        assert row["ntg_channel"] <= row["ntg_pure"] + 1e-6
        assert row["ntg_pure"] <= row["ntg_coarse"] + 1e-6


# -----------------------------------------------------------------------------
# 12. Synthetic Benchmark Baseline Models Execution & Output Grid Shape
# -----------------------------------------------------------------------------
def test_12_benchmark_baseline_models_execution():
    """Test 12: All 4 benchmark baseline models fit and predict valid grids."""
    from smalt.forward.benchmark import (
        NearestWellModel,
        VerticalPriorModel,
        SpatialKNNModel,
        SpatialMarkovModel,
    )

    sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
    realization = sim.generate(seed=42)
    sampler = SyntheticBoreholeSampler(realization)
    dataset = sampler.sample_regular(spacing_m=400.0, offset_m=200.0)

    models = {
        "Nearest_Well": NearestWellModel(),
        "Vertical_Prior": VerticalPriorModel(),
        "Spatial_KNN": SpatialKNNModel(n_neighbors=5, vertical_anisotropy=20.0),
        "Spatial_Markov": SpatialMarkovModel(decay_length_m=200.0),
    }

    for name, model in models.items():
        model.fit(dataset)
        pred_grid = model.predict(realization)
        assert pred_grid.shape == realization.shape, f"Model {name} output shape {pred_grid.shape} != {realization.shape}"
        assert np.all(pred_grid >= 0) and np.all(pred_grid <= 5), f"Model {name} predicted invalid facies codes"


# -----------------------------------------------------------------------------
# 13. Spacing Sensitivity Monotonicity
# -----------------------------------------------------------------------------
def test_13_spacing_sensitivity_accuracy_decay():
    """Test 13: Reconstruction accuracy decays as borehole spacing increases from 400m to 1000m."""
    from smalt.forward.benchmark import SyntheticReconstructionBenchmark

    bench = SyntheticReconstructionBenchmark()
    sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
    realization = sim.generate(seed=7000)

    # Evaluate at 400m vs 1000m on the same realization
    rec_400 = bench.evaluate_realization(realization, spacing_m=400.0)
    rec_1000 = bench.evaluate_realization(realization, spacing_m=1000.0)

    nw_400 = next(r for r in rec_400 if r["model"] == "Nearest_Well")
    nw_1000 = next(r for r in rec_1000 if r["model"] == "Nearest_Well")

    assert nw_400["raw_accuracy"] >= nw_1000["raw_accuracy"], (
        f"Expected 400m accuracy ({nw_400['raw_accuracy']}) >= 1000m accuracy ({nw_1000['raw_accuracy']})"
    )


# -----------------------------------------------------------------------------
# 14. Strict Leakage Prevention in Benchmark Pipeline
# -----------------------------------------------------------------------------
def test_14_benchmark_strict_leakage_safeguards():
    """Test 14: Benchmark harness masks sampled borehole columns during evaluation."""
    from smalt.forward.benchmark import SyntheticReconstructionBenchmark

    bench = SyntheticReconstructionBenchmark()
    sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
    realization = sim.generate(seed=8001)

    rec = bench.evaluate_realization(realization, spacing_m=400.0)
    for r in rec:
        # Evaluation cells count must be strictly less than total grid cells
        total_cells = realization.nx * realization.nz
        assert r["num_eval_cells"] < total_cells
        assert r["num_eval_cells"] == total_cells - (r["num_boreholes"] * realization.nz)


# -----------------------------------------------------------------------------
# 15. L1 Spatial Quarantine and L12 Depth Boundary Invariants
# -----------------------------------------------------------------------------
def test_15_l1_spatial_exclusion_and_l12_depth_boundary_invariants():
    """Test 15: L1 excluded from spatial models; L12 depth strictly <= 111.0m."""
    import pandas as pd
    import json

    with open("data/provenance_manifest.json") as f:
        manifest = json.load(f)

    # L1 must be marked as coordinates_available = False
    assert manifest["lithologs"]["litholog1"]["coordinates_available"] is False, (
        "Litholog 1 must not be marked as having coordinates"
    )

    # Coordinates Excel must not contain L1
    df_coords = pd.read_excel("lolo/Location_coordinates_lithologs.xlsx")
    well_names = df_coords.iloc[:, 0].dropna().astype(str).tolist()
    assert "L1" not in well_names, "L1 must not appear in physical location register"

    # L12 (EM-137C) coverage strictly <= 111.0 m across 63 intervals
    df12 = pd.read_csv("data/raw_lithologs/litholog12.csv")
    assert len(df12) == 63, f"Expected 63 intervals in L12, got {len(df12)}"
    assert df12["Bottom"].max() <= 111.0 + 1e-6, f"L12 max depth {df12['Bottom'].max()} exceeds 111.0m"
    assert manifest["lithologs"]["litholog12"]["digitized_interval_m"] == "0-111m"
    assert manifest["lithologs"]["litholog12"]["digitized_intervals_count"] == 63


# -----------------------------------------------------------------------------
# 16. Synthetic Coordinate Quarantine Invariant
# -----------------------------------------------------------------------------
def test_16_synthetic_coordinate_quarantine_invariant():
    """Test 16: Synthetic borehole dataset enforces synthetic_ coordinate prefix."""
    sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
    realization = sim.generate(seed=123)
    sampler = SyntheticBoreholeSampler(realization)
    dataset = sampler.sample_regular(spacing_m=400.0, offset_m=200.0)

    # Coordinate columns must be prefixed with synthetic_ or synth_
    for col in dataset.data.columns:
        if "x" in col and "idx" not in col:
            assert col.startswith("synthetic_") or col.startswith("synth_"), (
                f"Synthetic coordinate column {col} lacks quarantine prefix"
            )

    # Borehole IDs must carry synthetic prefix
    for bid in dataset.borehole_ids:
        assert bid.startswith("synthetic_") or bid.startswith("synth_"), (
            f"Borehole ID {bid} lacks quarantine prefix"
        )


# -----------------------------------------------------------------------------
# 17. Synthetic Gamma Ray Leakage Prohibition Invariant
# -----------------------------------------------------------------------------
def test_17_synthetic_gamma_ray_leakage_prohibition():
    """Test 17: Benchmark models predict solely from coordinates, never requiring synthetic GR."""
    from smalt.forward.benchmark import (
        NearestWellModel,
        VerticalPriorModel,
        SpatialKNNModel,
        SpatialMarkovModel,
    )

    sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
    realization = sim.generate(seed=456)
    sampler = SyntheticBoreholeSampler(realization)
    dataset = sampler.sample_regular(spacing_m=400.0, offset_m=200.0)

    # Verify that dataset columns do NOT contain gamma_ray
    assert "gamma_ray" not in dataset.data.columns, (
        "Synthetic borehole dataset should not expose gamma_ray as a feature to benchmark models"
    )

    models = [NearestWellModel(), VerticalPriorModel(), SpatialKNNModel(), SpatialMarkovModel()]
    for model in models:
        model.fit(dataset)
        pred = model.predict(realization)
        assert pred.shape == realization.shape


# -----------------------------------------------------------------------------
# 18. Lognormal Parameter-Estimation-Aware Goodness-of-Fit Test
# -----------------------------------------------------------------------------
def test_18_lognormal_parameter_estimation_aware_goodness_of_fit():
    """Test 18: Parameter-estimation-aware Lilliefors bootstrap test on Upstream channel beds."""
    from scipy import stats

    df_beds = pd.read_csv("sprints/audit_sprint_k/litholog_group_raw_beds.csv")
    upstream_sand = df_beds[(df_beds["group"] == "upstream") & (df_beds["facies_code"] == 0)]
    thicks = upstream_sand["thickness_m"].to_numpy()
    assert len(thicks) == 34

    log_thicks = np.log(thicks)
    mu_mle = float(np.mean(log_thicks))
    sigma_mle = float(np.std(log_thicks, ddof=1))
    assert mu_mle == pytest.approx(1.7870, abs=1e-3)
    assert sigma_mle == pytest.approx(0.8823, abs=1e-3)

    # Fast Lilliefors parametric bootstrap (B=300 for fast automated test)
    res_ks = stats.kstest(thicks, "lognorm", args=(sigma_mle, 0, np.exp(mu_mle)))
    np.random.seed(42)
    B = 300
    N = len(thicks)
    D_boot = np.empty(B)
    for b in range(B):
        sim = np.random.lognormal(mean=mu_mle, sigma=sigma_mle, size=N)
        sim_log = np.log(sim)
        m_b = np.mean(sim_log)
        s_b = np.std(sim_log, ddof=1)
        sim_sorted = np.sort(sim)
        cdf_vals = stats.lognorm.cdf(sim_sorted, s_b, scale=np.exp(m_b))
        d_plus = np.max((np.arange(1, N + 1) / N) - cdf_vals)
        d_minus = np.max(cdf_vals - (np.arange(0, N) / N))
        D_boot[b] = max(d_plus, d_minus)

    p_lilliefors = float(np.mean(D_boot >= res_ks.statistic))
    # Lilliefors p-value does not reject log-normal (p > 0.05)
    assert p_lilliefors > 0.05, f"Lilliefors bootstrap p-value {p_lilliefors} <= 0.05"


# -----------------------------------------------------------------------------
# 19. Metric Definitions Mathematical Consistency Invariant
# -----------------------------------------------------------------------------
def test_19_metric_definitions_mathematical_consistency():
    """Test 19: ntg_channel <= ntg_pure <= ntg_coarse invariant holds strictly."""
    sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
    for seed in [111, 222, 333]:
        res = sim.generate(seed=seed)
        ntg_chan = res.realized_net_to_gross["ntg_channel"]
        ntg_pure = res.realized_net_to_gross["ntg_pure"]
        ntg_coarse = res.realized_net_to_gross["ntg_coarse"]

        # Exact definitions
        sand_prop = res.facies_proportions["sand"]
        p_sand_prop = res.facies_proportions["p_sand"]
        ripples_prop = res.facies_proportions["ripples"]

        assert ntg_chan == pytest.approx(sand_prop, abs=1e-5)
        assert ntg_pure == pytest.approx(sand_prop + p_sand_prop, abs=1e-5)
        assert ntg_coarse == pytest.approx(sand_prop + p_sand_prop + ripples_prop, abs=1e-5)
        assert ntg_chan <= ntg_pure + 1e-6
        assert ntg_pure <= ntg_coarse + 1e-6


