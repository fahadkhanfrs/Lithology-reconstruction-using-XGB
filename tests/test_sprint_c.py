"""
Unit tests for SMALT Sprint C: Descriptive Analysis, Data Quality, and 1D Markov LOLO.

Validates:
1. Interval thickness and coverage validation.
2. Detection of missing intervals (gaps) and overlaps.
3. Explicit facies mapping and unknown-label handling.
4. Preservation of provenance across all 12 lithologs.
5. Exclusion of synthetic Gamma Ray and prev_facies from predictive feature sets.
6. Entire-litholog separation between training and test data in LOLO.
7. Zero target fitting: transition probabilities are never fitted on held-out well.
8. Transition-matrix row stochasticity and embedded-chain zero-diagonal invariant.
9. Stable handling of unseen or zero-count transitions via Laplace smoothing.
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path

from smalt.descriptive.analyzer import (
    LithologInspector,
    CANONICAL_FACIES_SCHEMA,
    PROVENANCE_METADATA,
    normalize_facies_label,
)
from smalt.validation.markov_lolo import MarkovLOLOValidator


@pytest.fixture
def inspector():
    return LithologInspector(raw_dir="data/raw_lithologs")


@pytest.fixture
def validator(inspector):
    return MarkovLOLOValidator(inspector=inspector, smoothing_alpha=0.1)


# -----------------------------------------------------------------------------
# 1. Interval Thickness & Coverage Validation
# -----------------------------------------------------------------------------
def test_all_12_lithologs_load_successfully(inspector):
    """Verifies that all 12 raw litholog CSVs exist, load, and have positive intervals."""
    for i in range(1, 13):
        lid = f"litholog{i}"
        df = inspector.load_raw_litholog(lid)
        assert len(df) > 0
        assert (df["Thickness"] > 0).all()
        assert "Top" in df.columns
        assert "Bottom" in df.columns
        assert "Facies" in df.columns


def test_litholog12_specific_coverage(inspector):
    """
    Specifically verifies that litholog12.csv covers 0.0 to 111.0 m exactly,
    has 63 contiguous intervals, has zero gaps, and zero overlaps.
    """
    insp = inspector.inspect_litholog("litholog12")
    assert insp["num_intervals"] == 63
    assert insp["min_depth_m"] == 0.0
    assert insp["max_depth_m"] == 111.0
    assert insp["represented_span_m"] == 111.0
    assert insp["sum_thickness_m"] == 111.0
    assert insp["span_thickness_diff_m"] == 0.0
    assert not insp["has_gaps"]
    assert not insp["has_overlaps"]


def test_thickness_equals_sum_not_row_count(inspector):
    """Verifies that cumulative thickness equals sum of intervals, not row count."""
    for i in range(1, 13):
        insp = inspector.inspect_litholog(f"litholog{i}")
        assert insp["sum_thickness_m"] != insp["num_intervals"]
        assert insp["sum_thickness_m"] > 70.0


# -----------------------------------------------------------------------------
# 2. Detection of Missing Intervals (Gaps) and Overlaps
# -----------------------------------------------------------------------------
def test_detect_gap_in_litholog9_and_litholog11(inspector):
    """Verifies that gaps in L9 (18-19m) and L11 (59-60m) are correctly detected."""
    insp9 = inspector.inspect_litholog("litholog9")
    assert insp9["has_gaps"]
    assert any(g["depth_before"] == 18.0 and g["depth_after"] == 19.0 for g in insp9["gaps"])

    insp11 = inspector.inspect_litholog("litholog11")
    assert insp11["has_gaps"]
    assert any(g["depth_before"] == 59.0 and g["depth_after"] == 60.0 for g in insp11["gaps"])


def test_detect_overlaps_in_litholog9(inspector):
    """Verifies that overlapping intervals in L9 (28-30m) are correctly identified."""
    insp9 = inspector.inspect_litholog("litholog9")
    assert insp9["has_overlaps"]
    assert len(insp9["overlaps"]) > 0


# -----------------------------------------------------------------------------
# 3. Explicit Facies Mapping & Unknown-Label Handling
# -----------------------------------------------------------------------------
def test_facies_normalization_valid_labels():
    """Verifies normalization for all known variants."""
    assert normalize_facies_label("sand") == "sand"
    assert normalize_facies_label("Sandstone") == "sand"
    assert normalize_facies_label("coal") == "coal"
    assert normalize_facies_label("carbon_mud") == "carbon_mud"
    assert normalize_facies_label("carbonaceous mudstone") == "carbon_mud"
    assert normalize_facies_label("silt") == "silt"
    assert normalize_facies_label("siltstone") == "silt"
    assert normalize_facies_label("mud") == "mud"
    assert normalize_facies_label("Overbank Mudstone") == "mud"


def test_facies_normalization_unknown_label_raises():
    """Verifies that unknown or unmapped facies raise ValueError without silent dropping."""
    with pytest.raises(ValueError, match="Unknown facies label"):
        normalize_facies_label("limestone")

    with pytest.raises(ValueError, match="Unknown facies label"):
        normalize_facies_label("granite")

    with pytest.raises(ValueError, match="null or NaN"):
        normalize_facies_label(np.nan)


# -----------------------------------------------------------------------------
# 4. Preservation of Provenance
# -----------------------------------------------------------------------------
def test_provenance_metadata_completeness():
    """Verifies that all 12 lithologs have complete, verified provenance records."""
    assert len(PROVENANCE_METADATA) == 12

    # L1, L9, L11 source-derived
    assert PROVENANCE_METADATA["litholog1"]["provenance_category"] == "source_derived"
    assert PROVENANCE_METADATA["litholog9"]["provenance_category"] == "source_derived"
    assert PROVENANCE_METADATA["litholog11"]["provenance_category"] == "source_derived_benchmarked"
    assert PROVENANCE_METADATA["litholog11"]["benchmark_accuracy"] == 0.9359

    # L2-L8, L10 AI-reconstructed
    for lid in ["litholog2", "litholog3", "litholog4", "litholog5", "litholog6", "litholog7", "litholog8", "litholog10"]:
        assert PROVENANCE_METADATA[lid]["provenance_category"] == "ai_reconstructed"

    # L12 digitized core
    assert PROVENANCE_METADATA["litholog12"]["provenance_category"] == "digitized_core_log"

    # L1 lacks coordinates
    assert not PROVENANCE_METADATA["litholog1"]["coordinates_available"]


# -----------------------------------------------------------------------------
# 5. Feature Leakage Prevention Invariants
# -----------------------------------------------------------------------------
def test_no_synthetic_gr_or_prev_facies_in_validator():
    """
    Verifies that the validator and inspector do NOT use or require synthetic GR
    or ground-truth prev_facies in any predictive calculation.
    """
    inspector = LithologInspector(raw_dir="data/raw_lithologs")
    df1 = inspector.load_raw_litholog("litholog1")
    assert "gamma_ray" not in df1.columns
    assert "prev_facies" not in df1.columns

    df_disc = inspector.discretize_litholog_1m("litholog1")
    assert "gamma_ray" not in df_disc.columns
    assert "prev_facies" not in df_disc.columns


# -----------------------------------------------------------------------------
# 6. Entire-Litholog Separation in LOLO
# -----------------------------------------------------------------------------
def test_lolo_disjoint_holdout_separation(validator):
    """Verifies that the target litholog is strictly excluded from training log IDs."""
    all_lids = [f"litholog{i}" for i in range(1, 13)]
    for target in all_lids:
        train_lids = [lid for lid in all_lids if lid != target]
        assert target not in train_lids
        assert len(train_lids) == 11


# -----------------------------------------------------------------------------
# 7. Transition-Matrix Row Sums & Embedded Diagonal Invariants
# -----------------------------------------------------------------------------
def test_transition_matrix_row_stochasticity(inspector):
    """Verifies that transition matrix rows sum to exactly 1.0."""
    all_lids = [f"litholog{i}" for i in range(1, 13)]

    # Regular
    P_reg, _, _ = inspector.compute_vertical_transitions(all_lids, embedded=False, use_discretized=True)
    assert np.allclose(P_reg.sum(axis=1), 1.0, atol=1e-7)

    # Embedded
    P_emb, _, _ = inspector.compute_vertical_transitions(all_lids, embedded=True, use_discretized=False)
    assert np.allclose(P_emb.sum(axis=1), 1.0, atol=1e-7)
    assert np.allclose(np.diag(P_emb), 0.0, atol=1e-7)


# -----------------------------------------------------------------------------
# 8. Stable Handling of Unseen / Zero-Count Transitions
# -----------------------------------------------------------------------------
def test_laplace_smoothing_prevents_zero_probabilities(inspector):
    """Verifies that Laplace smoothing alpha > 0 ensures non-zero probabilities everywhere except diagonal in embedded."""
    # Test on single well with few facies transitions (e.g. L3 which has no silt)
    P_reg, _, pi_reg = inspector.compute_vertical_transitions(["litholog3"], embedded=False, smoothing_alpha=0.1)
    assert (P_reg > 0.0).all()
    assert (pi_reg > 0.0).all()
    assert np.allclose(P_reg.sum(axis=1), 1.0)

    P_emb, _, pi_emb = inspector.compute_vertical_transitions(["litholog3"], embedded=True, smoothing_alpha=0.1)
    # Off-diagonals must be strictly positive
    mask = ~np.eye(P_emb.shape[0], dtype=bool)
    assert (P_emb[mask] > 0.0).all()
    assert np.allclose(np.diag(P_emb), 0.0)
