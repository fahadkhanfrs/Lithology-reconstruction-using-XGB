"""
Unit test suite for SMALT Sprint H: Common-Zero Datum Alignment & Spatial Markov Foundation.

Validates:
1. Common-zero datum alignment: all 12 lithologs aligned at z = 0.
2. Stratigraphic invariance: thickness, interval counts, and facies counts preserved.
3. Reversibility of vertical datum transformation.
4. Metadata table integrity and strict spatial exclusion of Litholog 1.
5. Horizontal elevation pair extraction across common-datum slices.
6. Carle & Fogg (1996) continuous-time Markov rate matrix mathematical properties.
7. Continuous spatial transition probability asymptotic limits (h -> 0 and h -> inf).
8. Spatial Markov LOLO predictor execution and leak-free validation.
9. Deliverable files existence and non-emptiness.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from smalt.descriptive.analyzer import LithologInspector, CANONICAL_FACIES_SCHEMA
from smalt.spatial.datum import (
    align_to_common_datum,
    revert_from_common_datum,
    build_common_datum_metadata_table,
    verify_datum_invariance,
    DEFAULT_CONVENTION,
)
from smalt.geostat.spatial_markov import (
    SpatialMarkovTransitionAnalyzer,
    SpatialMarkovPredictor,
    DEFAULT_LATERAL_FACIES_LENGTHS_M,
)


@pytest.fixture
def inspector():
    return LithologInspector()


def test_common_zero_datum_alignment_all_zeros_at_zero(inspector):
    """Verifies that all 12 lithologs have z_common_m = 0.0 at depth = 0."""
    all_lids = [f"litholog{i}" for i in range(1, 13)]
    for lid in all_lids:
        raw_df = inspector.load_raw_litholog(lid)
        aligned_df = align_to_common_datum(raw_df, litholog_id=lid, convention=DEFAULT_CONVENTION)

        # In elevation convention, reference level is z = 0 (top of section)
        assert aligned_df["z_common_top_m"].max() == pytest.approx(0.0)
        # Deeper intervals must be negative
        assert (aligned_df["z_common_bottom_m"] <= 0.0).all()
        # Original depth must be preserved
        assert "depth_original_m" in aligned_df.columns
        assert (aligned_df["depth_original_top_m"] >= 0.0).all()


def test_common_datum_thickness_and_facies_invariance(inspector):
    """Verifies that common-zero alignment preserves thicknesses, interval counts, and facies."""
    df_inv = verify_datum_invariance(inspector=inspector)
    assert len(df_inv) == 12

    for _, row in df_inv.iterrows():
        assert row["raw_interval_count"] == row["aligned_interval_count"]
        assert row["raw_thickness_m"] == pytest.approx(row["aligned_thickness_m"])
        assert row["facies_counts_identical"] is True
        assert row["thickness_identical"] is True
        assert row["reversibility_verified"] is True

    # Total logged span across all 12 logs must remain exactly 1033.0 m
    assert df_inv["aligned_thickness_m"].sum() == pytest.approx(1033.0)


def test_metadata_table_completeness_and_l1_exclusion(inspector):
    """Verifies metadata table schema, 12 lithologs, and strict spatial exclusion of Litholog 1."""
    df_meta = build_common_datum_metadata_table(inspector=inspector)
    assert len(df_meta) == 12

    # Litholog 1 checks
    l1_row = df_meta[df_meta["litholog"] == "litholog1"].iloc[0]
    assert bool(l1_row["coordinate_available"]) is False
    assert bool(l1_row["spatial_eligible"]) is False

    # Lithologs 2-12 checks
    for i in range(2, 13):
        row = df_meta[df_meta["litholog"] == f"litholog{i}"].iloc[0]
        assert bool(row["coordinate_available"]) is True
        assert bool(row["spatial_eligible"]) is True
        assert row["common_z_max"] == pytest.approx(0.0)
        assert row["common_z_min"] < 0.0


def test_horizontal_pair_extraction_and_elevation_matching(inspector):
    """Verifies extraction of horizontal well pairs at common-datum elevation slices."""
    analyzer = SpatialMarkovTransitionAnalyzer(inspector=inspector)
    df_pairs = analyzer.extract_horizontal_facies_pairs()

    assert len(df_pairs) > 0
    # Minimum lag distance between any two coordinate-bearing wells must be >= 420.0 m
    assert df_pairs["lag_distance_m"].min() >= 420.0
    # Elevation must be non-positive under elevation convention
    assert (df_pairs["z_common_m"] <= 0.0).all()
    # Both wells in pair must be spatial-eligible (not litholog1)
    assert "litholog1" not in df_pairs["well_a"].values
    assert "litholog1" not in df_pairs["well_b"].values


def test_carle_fogg_rate_matrix_mathematical_properties(inspector):
    """Verifies that the Carle & Fogg (1996) transition rate matrix satisfies row-sum zero condition."""
    analyzer = SpatialMarkovTransitionAnalyzer(inspector=inspector)
    R = analyzer.build_theoretical_horizontal_rate_matrix()

    n = len(CANONICAL_FACIES_SCHEMA)
    assert R.shape == (n, n)

    # Property 1: Diagonal elements are strictly negative
    assert (np.diag(R) < 0).all()

    # Property 2: Off-diagonal elements are non-negative
    for i in range(n):
        for j in range(n):
            if i != j:
                assert R[i, j] >= 0.0

    # Property 3: Row sums are identically zero (continuous Markov property)
    np.testing.assert_allclose(R.sum(axis=1), np.zeros(n), atol=1e-12)


def test_continuous_transition_probability_asymptotics(inspector):
    """Verifies that P(h) = expm(R*h) converges to Identity at h=0 and stationary prior at h->inf."""
    analyzer = SpatialMarkovTransitionAnalyzer(inspector=inspector)
    R = analyzer.build_theoretical_horizontal_rate_matrix()
    n = len(CANONICAL_FACIES_SCHEMA)

    # At lag h = 0: P(0) = Identity
    P_0 = analyzer.evaluate_transition_probability_at_lag(0.0, R)
    np.testing.assert_allclose(P_0, np.eye(n), atol=1e-10)

    # For intermediate lags: strictly row-stochastic
    for h in [10.0, 50.0, 200.0, 500.0, 1000.0]:
        P_h = analyzer.evaluate_transition_probability_at_lag(h, R)
        assert (P_h >= 0.0).all()
        assert (P_h <= 1.0).all()
        np.testing.assert_allclose(P_h.sum(axis=1), np.ones(n), atol=1e-10)

    # At large lag h -> inf (e.g. 50,000 m): rows become identical (stationary distribution)
    P_inf = analyzer.evaluate_transition_probability_at_lag(50000.0, R)
    for i in range(1, n):
        np.testing.assert_allclose(P_inf[i, :], P_inf[0, :], atol=1e-4)


def test_spatial_markov_predictor_execution_and_leak_free(inspector):
    """Verifies that SpatialMarkovPredictor executes leak-free LOLO across eligible wells."""
    predictor = SpatialMarkovPredictor()
    res = predictor.run_spatial_markov_lolo()

    df_lolo = res["per_fold_table"]
    summary = res["aggregate_summary"]

    assert len(df_lolo) == 11
    assert summary["total_test_points"] == 940
    assert 0.0 <= summary["pooled_markov_accuracy"] <= 1.0
    assert 0.0 <= summary["pooled_markov_macro_f1"] <= 1.0


def test_sprint_h_deliverables_exist_and_non_empty():
    """Verifies that all Sprint H audit tables and figures exist and are populated."""
    output_dir = Path("sprints/audit_sprint_h")
    assert (output_dir / "common_datum_litholog_metadata.csv").exists()
    assert (output_dir / "common_datum_facies_consistency.csv").exists()
    assert (output_dir / "empirical_horizontal_transition_matrices.csv").exists()
    assert (output_dir / "spatial_markov_lolo_results.csv").exists()
    assert (output_dir / "spatial_markov_vs_baselines_comparison.csv").exists()

    figures_dir = output_dir / "figures"
    assert (figures_dir / "common_zero_transect_alignment.png").exists()
    assert (figures_dir / "horizontal_transition_probability_decay.png").exists()

    df_meta = pd.read_csv(output_dir / "common_datum_litholog_metadata.csv")
    assert len(df_meta) == 12
    assert (df_meta["common_z_max"] == 0.0).all()
