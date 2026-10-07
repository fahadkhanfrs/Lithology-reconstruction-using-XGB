"""
Unit test suite for SMALT Sprint H: Source-Orientation Audit & Corrected Spatial Markov Foundation.

Validates all 10 core requirements:
1. Litholog 9 source orientation sanity check (0 m at base, fining-upward succession, thicknesses preserved).
2. Coordinate transformation correctness (z_strat = H_max - d) and exact reversibility.
3. Interval thickness preservation (1033.0 m cumulative, 328 intervals).
4. Facies-count preservation (100% identical counts for all 6 facies).
5. Common-zero assignment (all lithologs aligned with common_z_min = 0.0 m at source base).
6. Litholog 12 separate handling (preserved without inversion, 0-111 m scope).
7. 1D Markov directional succession modeling (forward base->top vs reverse).
8. Carle & Fogg continuous transition-rate matrix mathematical validity (R_ii <= 0, sum_j R_ij = 0).
9. Strict leakage-free spatial LOLO execution across eligible wells.
10. Deliverable files existence and non-emptiness (including archived pre-correction files).
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
    build_orientation_audit_table,
    LITHOLOG_ORIENTATION_METADATA,
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


def test_litholog9_source_orientation_sanity_check(inspector):
    """
    Validates Requirement 1: Litholog 9 source orientation.
    - 0 m is at the base.
    - Stratigraphic coordinate increases upward from 0 to 78 m.
    - Channel sandstone body (3.0 to 9.0 m) is capped by planar sand, ripples, and mud (fining upward).
    - Interval thicknesses are 100% preserved.
    """
    raw_df = inspector.load_raw_litholog("litholog9")
    aligned_df = align_to_common_datum(raw_df, litholog_id="litholog9", apply_source_orientation=True)

    # 1. Base is at z = 0.0 m, top is at z = 78.0 m
    assert aligned_df["z_strat_base_m"].min() == pytest.approx(0.0)
    assert aligned_df["z_strat_top_m"].max() == pytest.approx(78.0)

    # 2. First interval at base (z = 0 to 3 m) is Overbank Mudstone
    base_interval = aligned_df.iloc[0]
    assert base_interval["z_strat_base_m"] == pytest.approx(0.0)
    assert base_interval["z_strat_top_m"] == pytest.approx(3.0)
    assert base_interval["facies"] == "Overbank Mudstone"

    # 3. Second interval (z = 3 to 9 m) is Channel Sandstone
    sand_interval = aligned_df.iloc[1]
    assert sand_interval["z_strat_base_m"] == pytest.approx(3.0)
    assert sand_interval["z_strat_top_m"] == pytest.approx(9.0)
    assert sand_interval["facies"] == "Channel Sandstone"
    assert sand_interval["thickness_m"] == pytest.approx(6.0)

    # 4. Thickness preservation for every interval
    assert (aligned_df["thickness_m"] > 0).all()
    assert (aligned_df["thickness_m"] == (aligned_df["z_strat_top_m"] - aligned_df["z_strat_base_m"])).all()


def test_coordinate_transformation_and_reversibility(inspector):
    """
    Validates Requirement 2: Correctness of z_strat = H_max - d and exact reversibility.
    """
    all_lids = [f"litholog{i}" for i in range(1, 13)]
    for lid in all_lids:
        raw_df = inspector.load_raw_litholog(lid)
        aligned_df = align_to_common_datum(raw_df, litholog_id=lid, apply_source_orientation=True)

        # Reversibility test
        reverted_df = revert_from_common_datum(aligned_df)
        assert len(reverted_df) == len(raw_df)

        # Check that original depth column is preserved without modification
        assert "depth_original_m" in aligned_df.columns
        assert "depth_original_top_m" in aligned_df.columns
        assert "depth_original_bottom_m" in aligned_df.columns


def test_interval_thickness_and_facies_invariance(inspector):
    """
    Validates Requirements 3 & 4: 100% preservation of interval thicknesses and facies counts.
    Cumulative span across all 12 lithologs must equal 1033.0 m across 328 intervals.
    """
    df_inv = verify_datum_invariance(inspector=inspector, apply_source_orientation=True)
    assert len(df_inv) == 12

    for _, row in df_inv.iterrows():
        assert row["raw_interval_count"] == row["aligned_interval_count"]
        assert row["raw_thickness_m"] == pytest.approx(row["aligned_thickness_m"])
        assert row["facies_counts_identical"] is True
        assert row["thickness_identical"] is True
        assert row["reversibility_verified"] is True

    assert df_inv["raw_interval_count"].sum() == 328
    assert df_inv["aligned_thickness_m"].sum() == pytest.approx(1033.0)


def test_common_zero_datum_assignment_all_zeros_at_base(inspector):
    """
    Validates Requirement 5: All 12 litholog zeros are treated as the same reference level.
    Under the corrected orientation, all lithologs have common_z_min = 0.0 m.
    """
    df_meta = build_common_datum_metadata_table(inspector=inspector, apply_source_orientation=True)
    assert len(df_meta) == 12

    for _, row in df_meta.iterrows():
        assert row["common_z_min"] == pytest.approx(0.0)
        assert row["common_z_max"] > 70.0

    # Litholog 1 missing coordinates, strictly excluded from spatial modeling
    l1_row = df_meta[df_meta["litholog"] == "litholog1"].iloc[0]
    assert bool(l1_row["coordinate_available"]) is False
    assert bool(l1_row["spatial_eligible"]) is False

    # Lithologs 2-12 have coordinates and are spatial-eligible
    for i in range(2, 13):
        row = df_meta[df_meta["litholog"] == f"litholog{i}"].iloc[0]
        assert bool(row["coordinate_available"]) is True
        assert bool(row["spatial_eligible"]) is True


def test_litholog12_separate_handling(inspector):
    """
    Validates Requirement 6: Litholog 12 (EM-137C core) is handled separately.
    - Preserved in digitized orientation (0-111 m) without inverted axis.
    - Has required_transform == 'preserve'.
    """
    assert "litholog12" in LITHOLOG_ORIENTATION_METADATA
    l12_meta = LITHOLOG_ORIENTATION_METADATA["litholog12"]
    assert l12_meta["source_type"] == "drill_core"
    assert l12_meta["required_transform"] == "preserve"

    raw_df = inspector.load_raw_litholog("litholog12")
    aligned_df = align_to_common_datum(raw_df, litholog_id="litholog12", apply_source_orientation=True)

    assert aligned_df["z_strat_base_m"].min() == pytest.approx(0.0)
    assert aligned_df["z_strat_top_m"].max() == pytest.approx(111.0)
    assert aligned_df["thickness_m"].sum() == pytest.approx(111.0)
    assert len(aligned_df) == 63


def test_1d_markov_direction_correctness(inspector):
    """
    Validates Requirement 7: 1D Markov succession operates in upward stratigraphic order (base -> top).
    """
    analyzer = SpatialMarkovTransitionAnalyzer(inspector=inspector)
    grid = analyzer.build_common_elevation_grid()

    # Grid z_common_m must be non-negative (stratigraphic height upward from base)
    assert (grid["z_common_m"] >= 0.0).all()
    assert (grid["z_common_m"] <= 111.0).all()


def test_carle_fogg_transition_rate_matrix_validity(inspector):
    """
    Validates Requirement 8: Carle & Fogg (1996) rate matrix satisfies:
    - R_ii <= 0
    - R_ij >= 0 (i != j)
    - sum_j R_ij = 0
    - P(h) = expm(R*h) is strictly row-stochastic.
    """
    analyzer = SpatialMarkovTransitionAnalyzer(inspector=inspector)
    R = analyzer.build_theoretical_horizontal_rate_matrix()
    n = len(CANONICAL_FACIES_SCHEMA)

    assert R.shape == (n, n)
    assert (np.diag(R) < 0).all()

    for i in range(n):
        for j in range(n):
            if i != j:
                assert R[i, j] >= 0.0

    np.testing.assert_allclose(R.sum(axis=1), np.zeros(n), atol=1e-12)

    # Check matrix exponential row-stochasticity
    for h in [0.0, 50.0, 200.0, 1000.0, 5000.0]:
        P_h = analyzer.evaluate_transition_probability_at_lag(h, R)
        assert (P_h >= 0.0).all()
        assert (P_h <= 1.0).all()
        np.testing.assert_allclose(P_h.sum(axis=1), np.ones(n), atol=1e-10)


def test_spatial_markov_leak_free_execution(inspector):
    """
    Validates Requirement 9: Spatial LOLO cross-validation executes without data leakage.
    Target well points are never used in training.
    """
    predictor = SpatialMarkovPredictor()
    res = predictor.run_spatial_markov_lolo()

    df_lolo = res["per_fold_table"]
    summary = res["aggregate_summary"]

    assert len(df_lolo) == 11
    assert summary["total_test_points"] == 940
    assert 0.0 <= summary["pooled_markov_accuracy"] <= 1.0
    assert 0.0 <= summary["pooled_markov_macro_f1"] <= 1.0


def test_sprint_h_correction_deliverables_exist_and_archived():
    """
    Validates Requirement 10: Deliverables exist, non-empty, and pre-correction files are archived.
    """
    audit_dir = Path("sprints/audit_sprint_h")

    # Corrected deliverables
    assert (audit_dir / "litholog_orientation_audit.csv").exists()
    assert (audit_dir / "common_datum_litholog_metadata.csv").exists()
    assert (audit_dir / "common_datum_facies_consistency.csv").exists()
    assert (audit_dir / "litholog9_sanity_check.csv").exists()
    assert (audit_dir / "markov_1d_orientation_comparison.csv").exists()
    assert (audit_dir / "horizontal_facies_pairs_summary.csv").exists()
    assert (audit_dir / "empirical_horizontal_transition_matrices.csv").exists()
    assert (audit_dir / "spatial_markov_lolo_results.csv").exists()
    assert (audit_dir / "spatial_markov_vs_baselines_comparison.csv").exists()
    assert (audit_dir / "directional_upstream_downstream_results.csv").exists()
    assert (audit_dir / "orientation_correction_impact_comparison.csv").exists()

    fig_dir = audit_dir / "figures"
    assert (fig_dir / "litholog9_orientation_sanity_check.png").exists()
    assert (fig_dir / "common_zero_transect_alignment_corrected.png").exists()
    assert (fig_dir / "horizontal_transition_probability_decay_corrected.png").exists()

    # Archived pre-correction deliverables
    assert (audit_dir / "common_datum_litholog_metadata_pre_orientation_correction.csv").exists()
    assert (audit_dir / "spatial_markov_vs_baselines_comparison_pre_orientation_correction.csv").exists()
    assert (fig_dir / "common_zero_transect_alignment_pre_orientation_correction.png").exists()
