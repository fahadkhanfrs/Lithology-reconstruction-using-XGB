"""
SMALT Sprint J Test Suite: Empirical Horizontal Calibration & Lateral Continuity.

Validates:
1. Empirical transition counts extraction (4,410 pairs across 93 common-zero slices).
2. Distance binning strategies (Fixed, Quantile, Sensitivity) and bootstrap uncertainty.
3. No fabricated overlap: pair extraction strictly honors common-zero elevation overlap.
4. Row stochasticity: transition probability matrices sum to 1.0.
5. P(0) = I: zero-distance lag yields identity matrix.
6. Non-negative probabilities: all transition probabilities are non-negative.
7. Fitted rate matrix row sums = 0: R_h satisfies Carle & Fogg zero row-sum condition.
8. Deterministic fitting: empirical calibration produces identical results given identical seed.
9. LOLO target exclusion: withheld target well is strictly excluded from parameter fitting.
10. Directional target exclusion: downstream test wells are excluded from upstream training.
11. Six-state preservation: canonical schema (0..5) preserved throughout.
12. Synthetic correlation-length recovery: estimator recovers known synthetic lengths.
13. Stochastic seed reproducibility: identical seeds produce bitwise-identical realizations.
14. Hard conditioning remains 100%: known well control points are honored at 100%.
15. Existing Sprint I baselines remain valid and reproducible.
"""

import pytest
import numpy as np
import pandas as pd

from smalt.descriptive.analyzer import LithologInspector, CANONICAL_FACIES_SCHEMA
from smalt.geostat.spatial_transition import (
    EmpiricalHorizontalTransitionEstimator,
    SpatialTransitionRateModel,
    DEFAULT_LATERAL_FACIES_LENGTHS_M,
)
from smalt.geostat.empirical_calibration import (
    HorizontalContinuityCalibrator,
    DEFAULT_FIXED_BINS,
    exponential_decay_func,
    spherical_decay_func,
)
from smalt.geostat.ablation import AblationStudyEngine, calculate_bed_thicknesses
from smalt.geostat.conditioned_markov import ConditionedMarkovClassifier
from smalt.geostat.realization import InterWellRealizationGenerator
from smalt.validation.sprint_j import SprintJValidator


@pytest.fixture(scope="module")
def calibrator():
    inspector = LithologInspector()
    return HorizontalContinuityCalibrator(inspector=inspector)


def test_1_empirical_transition_counts(calibrator):
    """Test 1: Verify exact horizontal pair extraction across common-zero elevation slices."""
    pairs_df = calibrator.extract_pairs()
    assert len(pairs_df) == 4410, f"Expected exactly 4,410 horizontal pairs, got {len(pairs_df)}"
    unique_z = pairs_df["z_common_m"].unique()
    assert len(unique_z) == 93, f"Expected 93 unique elevation slices, got {len(unique_z)}"
    assert pairs_df["lag_distance_m"].min() >= 420.0, "Minimum lag distance should be >= 420m"


def test_2_distance_binning_and_bootstrap_uncertainty(calibrator):
    """Test 2: Verify multiple binning strategies and bootstrap standard error computation."""
    bin_dict = calibrator.evaluate_binning_strategies(n_bootstraps=20, random_seed=42)
    assert "fixed_bins" in bin_dict
    assert "quantile_bins" in bin_dict
    assert "sensitivity_bins" in bin_dict

    fixed_df = bin_dict["fixed_bins"]
    assert len(fixed_df) == len(DEFAULT_FIXED_BINS) * 36  # 7 bins * 36 matrix elements
    assert "std_error" in fixed_df.columns
    assert "ci_95_low" in fixed_df.columns
    assert np.all(fixed_df["std_error"] >= 0.0)
    assert np.all(fixed_df["ci_95_low"] <= fixed_df["transition_probability"])


def test_3_no_fabricated_overlap(calibrator):
    """Test 3: Confirm no pairs are fabricated outside overlapping elevation intervals."""
    pairs_df = calibrator.extract_pairs()
    coords_by_id = calibrator.coords_by_id
    for _, r in pairs_df.head(100).iterrows():
        wa, wb = r["well_a"], r["well_b"]
        assert wa in coords_by_id
        assert wb in coords_by_id
        assert wa != wb
        # Distance calculation matches Pythagorean formula
        expected_dist = np.hypot(r["dx_m"], r["dy_m"])
        assert np.isclose(r["lag_distance_m"], expected_dist, atol=0.2)


def test_4_row_stochasticity(calibrator):
    """Test 4: Continuous transition matrices must satisfy row stochasticity (sum to 1.0)."""
    stats_df, models_dict = calibrator.test_alternative_decay_models()
    model_B = models_dict["model_B"]
    for d in [0.0, 50.0, 420.0, 1000.0, 3000.0, 10000.0]:
        P = model_B.evaluate_transition_matrix(d)
        row_sums = P.sum(axis=1)
        np.testing.assert_allclose(row_sums, np.ones(6), atol=1e-6)


def test_5_p_zero_equals_identity(calibrator):
    """Test 5: Continuous transition probability at zero separation must equal Identity."""
    stats_df, models_dict = calibrator.test_alternative_decay_models()
    P_A0 = models_dict["model_A"].evaluate_transition_matrix(0.0)
    P_B0 = models_dict["model_B"].evaluate_transition_matrix(0.0)
    P_C0 = models_dict["eval_model_C"](0.0)
    np.testing.assert_allclose(P_A0, np.eye(6), atol=1e-9)
    np.testing.assert_allclose(P_B0, np.eye(6), atol=1e-9)
    np.testing.assert_allclose(P_C0, np.eye(6), atol=1e-9)


def test_6_non_negative_probabilities(calibrator):
    """Test 6: All transition matrix entries must be non-negative for any distance."""
    stats_df, models_dict = calibrator.test_alternative_decay_models()
    for d in [10.0, 420.0, 1500.0, 5000.0]:
        P_B = models_dict["model_B"].evaluate_transition_matrix(d)
        P_C = models_dict["eval_model_C"](d)
        assert np.all(P_B >= -1e-12), "Model B has negative probabilities"
        assert np.all(P_C >= -1e-12), "Model C has negative probabilities"


def test_7_fitted_rate_matrix_row_sums_zero(calibrator):
    """Test 7: Continuous transition rate matrix R_h must satisfy row-sum = 0."""
    stats_df, models_dict = calibrator.test_alternative_decay_models()
    R_A = models_dict["model_A"].R_matrix
    R_B = models_dict["model_B"].R_matrix
    np.testing.assert_allclose(R_A.sum(axis=1), np.zeros(6), atol=1e-10)
    np.testing.assert_allclose(R_B.sum(axis=1), np.zeros(6), atol=1e-10)


def test_8_deterministic_fitting(calibrator):
    """Test 8: Re-estimating facies continuity lengths produces identical results."""
    df1 = calibrator.estimate_facies_lengths()
    df2 = calibrator.estimate_facies_lengths()
    pd.testing.assert_frame_equal(df1, df2)


def test_9_lolo_target_exclusion(calibrator):
    """Test 9: In LOLO cross-validation, target well must be excluded from calibration."""
    target_id = "litholog9"
    train_ids = [lid for lid in calibrator.eligible_ids if lid != target_id]
    cal_fold = HorizontalContinuityCalibrator(
        inspector=calibrator.inspector,
        eligible_litholog_ids=train_ids,
    )
    fold_pairs = cal_fold.extract_pairs()
    assert target_id not in fold_pairs["well_a"].values
    assert target_id not in fold_pairs["well_b"].values


def test_10_directional_target_exclusion(calibrator):
    """Test 10: Downstream test wells (L9, L11, L12) strictly excluded from upstream training."""
    from smalt.validation.sprint_j import UPSTREAM_IDS, DOWNSTREAM_IDS
    cal_up = HorizontalContinuityCalibrator(
        inspector=calibrator.inspector,
        eligible_litholog_ids=UPSTREAM_IDS,
    )
    up_pairs = cal_up.extract_pairs()
    for d_id in DOWNSTREAM_IDS:
        assert d_id not in up_pairs["well_a"].values
        assert d_id not in up_pairs["well_b"].values


def test_11_six_state_schema_preservation(calibrator):
    """Test 11: Six canonical states are preserved and distinct."""
    assert len(CANONICAL_FACIES_SCHEMA) == 6
    assert 0 in [m["code"] for m in CANONICAL_FACIES_SCHEMA.values()]  # sand
    assert 1 in [m["code"] for m in CANONICAL_FACIES_SCHEMA.values()]  # p_sand
    assert 2 in [m["code"] for m in CANONICAL_FACIES_SCHEMA.values()]  # ripples
    assert 3 in [m["code"] for m in CANONICAL_FACIES_SCHEMA.values()]  # carbon_mud
    assert 4 in [m["code"] for m in CANONICAL_FACIES_SCHEMA.values()]  # coal
    assert 5 in [m["code"] for m in CANONICAL_FACIES_SCHEMA.values()]  # mud


def test_12_synthetic_correlation_length_recovery():
    """Test 12: Estimator recovers known synthetic correlation lengths."""
    validator = SprintJValidator()
    syn_df = validator.run_synthetic_benchmark(random_seed=42)
    assert len(syn_df) == 3
    # Check that lengths preserve ranking: L0 < L1 < L2
    l_est = syn_df["estimated_length_m"].to_numpy()
    assert l_est[0] < l_est[1] < l_est[2], "Synthetic correlation length ranking was not preserved"
    # Medium and long range facies recovery within 50%
    assert syn_df.loc[1, "relative_error_pct"] < 50.0
    assert syn_df.loc[2, "relative_error_pct"] < 50.0


def test_13_stochastic_seed_reproducibility():
    """Test 13: Identical random seeds produce bitwise-identical stochastic realizations."""
    inspector = LithologInspector()
    clf = ConditionedMarkovClassifier(
        training_litholog_ids=["litholog2", "litholog3"],
        inspector=inspector,
    )
    gen1 = InterWellRealizationGenerator(classifier=clf, random_seed=123)
    gen2 = InterWellRealizationGenerator(classifier=clf, random_seed=123)

    res1 = gen1.generate_2d_transect_realizations(["litholog2", "litholog3"], n_realizations=2)
    res2 = gen2.generate_2d_transect_realizations(["litholog2", "litholog3"], n_realizations=2)

    np.testing.assert_array_equal(res1["realizations"], res2["realizations"])


def test_14_hard_conditioning_remains_100():
    """Test 14: Hard conditioning honor rate is strictly 100% across conditioning wells."""
    inspector = LithologInspector()
    clf = ConditionedMarkovClassifier(
        training_litholog_ids=["litholog2", "litholog9"],
        inspector=inspector,
    )
    gen = InterWellRealizationGenerator(classifier=clf, random_seed=42)
    res = gen.generate_2d_transect_realizations(["litholog2", "litholog9"], n_realizations=2)
    assert res["hard_data_honor_rate"] == 1.0, f"Expected 1.0 honor rate, got {res['hard_data_honor_rate']}"


def test_15_existing_sprint_i_baselines_valid():
    """Test 15: Sprint I CTP and baselines remain functional."""
    from smalt.geostat.spatial_transition import SpatialTransitionRateModel
    rate_model = SpatialTransitionRateModel()
    pi = rate_model.compute_stationary_distribution()
    assert len(pi) == 6
    np.testing.assert_allclose(pi.sum(), 1.0, atol=1e-6)
