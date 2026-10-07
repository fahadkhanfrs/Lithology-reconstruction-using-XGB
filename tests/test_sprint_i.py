"""
Sprint I Comprehensive Test Suite.

Covers all 15 required verification constraints:
1. common-zero coordinate handling
2. empirical horizontal transition counting
3. no observations outside overlap being fabricated
4. transition matrix row sums = 1
5. P(0) = I
6. non-negative probabilities
7. long-distance convergence behavior where applicable
8. hard-data conditioning
9. hard-data honor rate (100%)
10. probability vector sums to 1
11. deterministic random seed
12. target litholog excluded from model fitting
13. L1 excluded from spatial validation
14. six-state schema preserved
15. synthetic benchmark sanity test
"""

import pytest
import numpy as np
import pandas as pd

from smalt.descriptive.analyzer import LithologInspector, CANONICAL_FACIES_SCHEMA
from smalt.spatial.datum import align_to_common_datum
from smalt.geostat.spatial_transition import (
    EmpiricalHorizontalTransitionEstimator,
    SpatialTransitionRateModel,
    DEFAULT_LATERAL_FACIES_LENGTHS_M,
)
from smalt.geostat.conditioned_markov import ConditionedMarkovClassifier
from smalt.geostat.realization import InterWellRealizationGenerator
from smalt.validation.sprint_i import SprintIValidator


@pytest.fixture
def inspector():
    return LithologInspector()


@pytest.fixture
def empirical_estimator(inspector):
    return EmpiricalHorizontalTransitionEstimator(inspector=inspector)


@pytest.fixture
def rate_model():
    return SpatialTransitionRateModel()


def test_1_common_zero_coordinate_handling(inspector):
    """1. Verifies common-zero coordinate alignment for all measured sections (z >= 0, base at 0)."""
    for lid in ["litholog2", "litholog9", "litholog11"]:
        disc_df = inspector.discretize_litholog_1m(lid)
        aligned = align_to_common_datum(disc_df, litholog_id=lid, apply_source_orientation=True)
        assert "z_common_m" in aligned.columns
        assert aligned["z_common_m"].min() >= 0.0
        assert "depth_original_m" in aligned.columns
        # Original thickness preserved
        assert len(aligned) == len(disc_df)


def test_2_empirical_horizontal_transition_counting(empirical_estimator):
    """2. Verifies empirical horizontal transition extraction and pair count."""
    df_pairs = empirical_estimator.extract_horizontal_facies_pairs()
    assert len(df_pairs) > 0
    assert "lag_distance_m" in df_pairs.columns
    assert "facies_code_a" in df_pairs.columns
    assert "facies_code_b" in df_pairs.columns
    # Check that all facies codes are within 0..5
    assert set(df_pairs["facies_code_a"].unique()).issubset(set(range(6)))
    assert set(df_pairs["facies_code_b"].unique()).issubset(set(range(6)))


def test_3_no_observations_outside_overlap_fabricated(empirical_estimator):
    """3. Verifies that only true overlapping elevation levels are paired, no synthetic extrapolation."""
    df_pairs = empirical_estimator.extract_horizontal_facies_pairs()
    # Unique z values must come strictly from common elevation slices
    grid_df = empirical_estimator.build_common_elevation_grid()
    grid_z = set(grid_df["z_common_m"].unique())
    pair_z = set(df_pairs["z_common_m"].unique())
    assert pair_z.issubset(grid_z)


def test_4_transition_matrix_row_sums_equal_one(rate_model):
    """4. Verifies that continuous transition matrices P(h) have row sums equal to 1.0."""
    for h in [0.0, 10.0, 100.0, 500.0, 2000.0, 10000.0]:
        P_h = rate_model.evaluate_transition_matrix(h)
        row_sums = P_h.sum(axis=1)
        np.testing.assert_allclose(row_sums, np.ones(6), atol=1e-10)


def test_5_p0_equals_identity(rate_model):
    """5. Verifies that at zero lag distance h = 0, P(0) equals the Identity matrix."""
    P_0 = rate_model.evaluate_transition_matrix(0.0)
    np.testing.assert_allclose(P_0, np.eye(6), atol=1e-12)


def test_6_non_negative_probabilities(rate_model):
    """6. Verifies that all transition probabilities are non-negative."""
    for h in [0.1, 50.0, 250.0, 1000.0, 5000.0]:
        P_h = rate_model.evaluate_transition_matrix(h)
        assert np.all(P_h >= 0.0)
        assert np.all(P_h <= 1.0)


def test_7_long_distance_convergence_behavior(rate_model):
    """7. Verifies that as h -> inf, rows of P(h) converge toward the invariant stationary distribution."""
    P_inf = rate_model.evaluate_transition_matrix(1e6)
    pi_stat = rate_model.compute_stationary_distribution()
    for i in range(6):
        np.testing.assert_allclose(P_inf[i, :], pi_stat, atol=1e-3)


def test_8_hard_data_conditioning(inspector):
    """8. Verifies that at conditioning well coordinates, conditioning facies is returned with P=1.0."""
    train_ids = ["litholog2", "litholog3", "litholog4"]
    classifier = ConditionedMarkovClassifier(training_litholog_ids=train_ids, inspector=inspector)
    coords = classifier.coords_by_id["litholog2"]

    # Pick a point from litholog2
    disc_df = inspector.discretize_litholog_1m("litholog2")
    aligned = align_to_common_datum(disc_df, litholog_id="litholog2", apply_source_orientation=True)
    sample_row = aligned.iloc[0]
    z_val = float(sample_row["z_common_m"])
    true_code = int(sample_row["facies_code"])

    res = classifier.evaluate_conditional_probability_at_point(
        x_m=coords[0],
        y_m=coords[1],
        z_common_m=z_val,
        hard_tolerance_m=1.0,
    )
    assert res["is_hard_conditioned"] is True
    assert res["map_facies_code"] == true_code
    assert res["probabilities"][true_code] == 1.0
    assert res["predictive_entropy"] == 0.0


def test_9_hard_data_honor_rate_equals_100_percent(inspector):
    """9. Verifies that the classifier's hard data honor rate across all training points is 100%."""
    train_ids = ["litholog2", "litholog3"]
    classifier = ConditionedMarkovClassifier(training_litholog_ids=train_ids, inspector=inspector)
    honor_rate = classifier.compute_hard_data_honor_rate()
    assert honor_rate == 1.0


def test_10_probability_vector_sums_to_one(inspector):
    """10. Verifies that conditional probability vectors sum to 1.0 at any arbitrary query location."""
    train_ids = ["litholog2", "litholog3", "litholog4"]
    classifier = ConditionedMarkovClassifier(training_litholog_ids=train_ids, inspector=inspector)
    # Query at arbitrary inter-well location
    res = classifier.evaluate_conditional_probability_at_point(
        x_m=500000.0,
        y_m=4300000.0,
        z_common_m=45.0,
        hard_tolerance_m=0.0,
    )
    p_vec = res["probabilities"]
    assert len(p_vec) == 6
    np.testing.assert_allclose(p_vec.sum(), 1.0, atol=1e-10)
    assert np.all(p_vec >= 0.0)


def test_11_deterministic_random_seed(inspector):
    """11. Verifies that InterWellRealizationGenerator produces identical realizations with identical seeds."""
    train_ids = ["litholog2", "litholog3"]
    classifier = ConditionedMarkovClassifier(training_litholog_ids=train_ids, inspector=inspector)

    gen1 = InterWellRealizationGenerator(classifier=classifier, random_seed=42)
    res1 = gen1.generate_2d_transect_realizations(
        anchor_well_ids=["litholog2", "litholog3"],
        n_realizations=3,
        x_resolution_m=50.0,
    )

    gen2 = InterWellRealizationGenerator(classifier=classifier, random_seed=42)
    res2 = gen2.generate_2d_transect_realizations(
        anchor_well_ids=["litholog2", "litholog3"],
        n_realizations=3,
        x_resolution_m=50.0,
    )

    np.testing.assert_array_equal(res1["realizations"], res2["realizations"])


def test_12_target_litholog_excluded_from_model_fitting(inspector):
    """12. Verifies that target litholog in LOLO validation is strictly absent from training parameters."""
    target_id = "litholog9"
    train_ids = [f"litholog{i}" for i in range(2, 9)]  # L2-L8
    classifier = ConditionedMarkovClassifier(training_litholog_ids=train_ids, inspector=inspector)

    assert target_id not in classifier.training_ids
    assert target_id not in classifier.coords_by_id
    assert target_id not in classifier.training_obs_df["litholog_id"].values


def test_13_litholog1_excluded_from_spatial_validation(inspector):
    """13. Verifies that Litholog 1 (missing coordinates) is strictly excluded from spatial models."""
    estimator = EmpiricalHorizontalTransitionEstimator(inspector=inspector)
    assert "litholog1" not in estimator.eligible_ids

    validator = SprintIValidator(inspector=inspector)
    assert "litholog1" not in validator.eligible_ids


def test_14_six_state_schema_preserved(inspector, rate_model):
    """14. Verifies that the six-state canonical schema is strictly preserved without facies collapse."""
    assert len(CANONICAL_FACIES_SCHEMA) == 6
    assert rate_model.num_classes == 6
    # p_sand and ripples must be distinct states
    assert CANONICAL_FACIES_SCHEMA["p_sand"]["code"] == 1
    assert CANONICAL_FACIES_SCHEMA["ripples"]["code"] == 2
    assert CANONICAL_FACIES_SCHEMA["p_sand"]["canonical_name"] == "Planar Sandstone"
    assert CANONICAL_FACIES_SCHEMA["ripples"]["canonical_name"] == "Rippled Heterolithics"


def test_15_synthetic_benchmark_sanity_test(inspector):
    """15. Executes the synthetic sanity benchmark and verifies all test items pass."""
    validator = SprintIValidator(inspector=inspector)
    synth_res = validator.run_synthetic_benchmark()
    df_bench = synth_res["benchmark_table"]
    assert len(df_bench) >= 5
    assert all(df_bench["status"] == "PASSED")
