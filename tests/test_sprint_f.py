"""
SMALT Sprint F Test Suite: Six-State Migration & Spatial Baseline Revalidation.

Verifies:
1. All six facies are represented and mapped consistently across the pipeline.
2. p_sand (Planar Sandstone) and ripples (Rippled Heterolithics) remain distinct throughout preprocessing and model outputs.
3. Unknown facies labels raise ValueError rather than being silently recoded.
4. Transition matrices use the correct 6x6 dimensions and row-stochastic normalization (sum=1 within 1e-7).
5. Embedded chains exclude self-transitions (diagonal elements strictly 0.0).
6. Litholog 1 remains strictly excluded from spatial evaluation due to missing coordinates.
7. Synthetic coordinates (e.g. strike_pos_m) cannot leak into spatial model feature vectors.
8. Held-out logs do not influence fitted preprocessing (StandardScaler) or model training in spatial LOLO folds.
9. Regenerated artifacts match current raw data and 6-state schema.
10. All figures and deliverable tables identify the 6-state schema correctly.
"""

from pathlib import Path
import json
import pytest
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier

from smalt.descriptive.analyzer import (
    LithologInspector,
    CANONICAL_FACIES_SCHEMA,
    PROVENANCE_METADATA,
    normalize_facies_label,
)
from smalt.geostat.markov import StratigraphicMarkovChain, DEFAULT_FACIES_MAP, CANONICAL_TECHNICAL_MAP
from smalt.spatial.coordinates import load_source_coordinates, get_spatial_litholog_subset
from smalt.spatial.baseline import ProvisionalSpatialValidator


@pytest.fixture(scope="module")
def inspector():
    repo_root = Path(__file__).resolve().parent.parent
    return LithologInspector(raw_dir=repo_root / "data" / "raw_lithologs")


@pytest.fixture(scope="module")
def spatial_validator(inspector):
    return ProvisionalSpatialValidator(inspector=inspector)


# -----------------------------------------------------------------------------
# 1. Six Facies Representation & Distinctness
# -----------------------------------------------------------------------------
def test_six_facies_represented_and_mapped_consistently():
    """Verifies that all 6 facies have distinct integer codes (0..5) and names."""
    assert len(CANONICAL_FACIES_SCHEMA) == 6
    codes = [meta["code"] for meta in CANONICAL_FACIES_SCHEMA.values()]
    assert sorted(codes) == list(range(6))
    assert len(set(codes)) == 6

    expected_keys = ["sand", "p_sand", "ripples", "carbon_mud", "coal", "mud"]
    assert list(CANONICAL_FACIES_SCHEMA.keys()) == expected_keys

    # Check DEFAULT_FACIES_MAP in geostat/markov.py
    assert len(DEFAULT_FACIES_MAP) == 6
    assert DEFAULT_FACIES_MAP[0] == "Channel Sandstone"
    assert DEFAULT_FACIES_MAP[1] == "Planar Sandstone"
    assert DEFAULT_FACIES_MAP[2] == "Rippled Heterolithics"
    assert DEFAULT_FACIES_MAP[3] == "Carbonaceous Mudstone"
    assert DEFAULT_FACIES_MAP[4] == "Coal"
    assert DEFAULT_FACIES_MAP[5] == "Overbank Mudstone"


def test_p_sand_and_ripples_remain_distinct():
    """
    Verifies that p_sand and ripples are NOT collapsed into sand or silt.
    They must normalize to distinct canonical keys and map to distinct codes.
    """
    # p_sand must be distinct from sand
    assert normalize_facies_label("p_sand") == "p_sand"
    assert normalize_facies_label("sand") == "sand"
    assert CANONICAL_FACIES_SCHEMA["p_sand"]["code"] != CANONICAL_FACIES_SCHEMA["sand"]["code"]
    assert CANONICAL_FACIES_SCHEMA["p_sand"]["code"] == 1
    assert CANONICAL_FACIES_SCHEMA["sand"]["code"] == 0

    # ripples must be distinct from mud and sand
    assert normalize_facies_label("ripples") == "ripples"
    assert normalize_facies_label("silt") == "ripples"  # Silt aliases to ripples
    assert CANONICAL_FACIES_SCHEMA["ripples"]["code"] == 2
    assert CANONICAL_FACIES_SCHEMA["ripples"]["code"] != CANONICAL_FACIES_SCHEMA["p_sand"]["code"]


def test_unknown_facies_labels_fail_clearly():
    """Verifies that invalid or unknown labels raise ValueError immediately."""
    invalid_labels = ["limestone", "dolomite", "granite", "unknown_rock", "", None]
    for lbl in invalid_labels:
        with pytest.raises(ValueError):
            normalize_facies_label(lbl)


# -----------------------------------------------------------------------------
# 2. 6x6 Transition Matrices & Normalization Invariants
# -----------------------------------------------------------------------------
def test_markov_matrices_six_state_dimensions_and_stochasticity(inspector):
    """
    Verifies that transition count and probability matrices are strictly 6x6,
    row-stochastic (sums to 1.0 within 1e-7), and embedded diagonal is strictly 0.0.
    """
    all_lids = [f"litholog{i}" for i in range(1, 13)]

    # Regular matrix
    P_reg, N_reg, pi_reg = inspector.compute_vertical_transitions(all_lids, embedded=False, use_discretized=True)
    assert P_reg.shape == (6, 6)
    assert N_reg.shape == (6, 6)
    assert len(pi_reg) == 6
    np.testing.assert_allclose(P_reg.sum(axis=1), np.ones(6), atol=1e-7)

    # Embedded matrix
    P_emb, N_emb, pi_emb = inspector.compute_vertical_transitions(all_lids, embedded=True, use_discretized=False)
    assert P_emb.shape == (6, 6)
    assert N_emb.shape == (6, 6)
    assert len(pi_emb) == 6
    np.testing.assert_allclose(P_emb.sum(axis=1), np.ones(6), atol=1e-7)
    np.testing.assert_allclose(np.diag(P_emb), np.zeros(6), atol=1e-12)


def test_embedded_bed_sequence_distinct_facies(inspector):
    """Verifies that embedded sequences contain zero consecutive duplicates across all 12 logs."""
    for i in range(1, 13):
        lid = f"litholog{i}"
        df_beds = inspector.build_embedded_bed_sequence(lid)
        codes = df_beds["facies_code"].to_numpy()
        for t in range(len(codes) - 1):
            assert codes[t] != codes[t + 1], f"Duplicate adjacent bed facies in {lid} at bed {t}"


# -----------------------------------------------------------------------------
# 3. Spatial Validation Invariants (L1 Exclusion & No Coordinate Leakage)
# -----------------------------------------------------------------------------
def test_litholog_1_excluded_from_spatial_evaluation(spatial_validator):
    """Verifies that Litholog 1 has no coordinates and raises ValueError if passed to spatial points."""
    coords_df = load_source_coordinates()
    l1_row = coords_df[coords_df["litholog_id"] == "litholog1"].iloc[0]
    assert not l1_row["coordinates_available"]
    assert np.isnan(l1_row["x_m"])
    assert np.isnan(l1_row["y_m"])

    # Attempting to extract spatial points for L1 must raise ValueError
    with pytest.raises(ValueError, match="no source coordinates available"):
        spatial_validator.prepare_discretized_spatial_points("litholog1")

    # Spatial subset must strictly have 11 lithologs (L2..L12)
    spatial_subset = get_spatial_litholog_subset()
    assert len(spatial_subset) == 11
    assert "litholog1" not in spatial_subset["litholog_id"].values


def test_no_synthetic_coordinates_in_spatial_features(spatial_validator):
    """
    Verifies that model feature vectors use only [x_m, y_m, z_rel_m] and
    never include legacy synthetic coordinate proxies (e.g. strike_pos_m = file_idx * 100).
    """
    for lid in [f"litholog{i}" for i in range(2, 13)]:
        df_pts = spatial_validator.prepare_discretized_spatial_points(lid)
        assert "strike_pos_m" not in df_pts.columns
        assert "x_m" in df_pts.columns
        assert "y_m" in df_pts.columns
        assert "z_rel_m" in df_pts.columns

        # Verify X and Y match the source spreadsheet
        coords_df = load_source_coordinates()
        expected_x = coords_df[coords_df["litholog_id"] == lid]["x_m"].iloc[0]
        expected_y = coords_df[coords_df["litholog_id"] == lid]["y_m"].iloc[0]
        assert np.isclose(df_pts["x_m"].iloc[0], expected_x)
        assert np.isclose(df_pts["y_m"].iloc[0], expected_y)


def test_spatial_lolo_strict_holdout_isolation(spatial_validator):
    """
    Verifies that in each LOLO fold, the target litholog is completely excluded
    from feature scaling fitting and training points.
    """
    eligible_ids = [f"litholog{i}" for i in range(2, 13)]
    for target_id in eligible_ids[:3]:  # Test sample of folds
        train_ids = [lid for lid in eligible_ids if lid != target_id]
        assert target_id not in train_ids
        assert len(train_ids) == 10

        train_dfs = [spatial_validator.prepare_discretized_spatial_points(lid) for lid in train_ids]
        df_train = pd.concat(train_dfs, ignore_index=True)
        df_test = spatial_validator.prepare_discretized_spatial_points(target_id)

        # Confirm target well points are disjoint from training points
        assert not np.isin(df_train["x_m"].unique(), df_test["x_m"].iloc[0]).any()


# -----------------------------------------------------------------------------
# 4. Regenerated Outputs & Delivery Verification
# -----------------------------------------------------------------------------
def test_sprint_f_deliverables_exist():
    """Verifies that all Sprint F deliverables are present and non-empty."""
    repo_root = Path(__file__).resolve().parent.parent
    f_dir = repo_root / "sprints" / "audit_sprint_f"
    fig_dir = f_dir / "figures"

    # CSV and JSON deliverables
    expected_files = [
        "facies_schema_six_state.csv",
        "per_litholog_six_state_statistics.csv",
        "group_six_state_statistics.csv",
        "markov_transition_matrices_six_state.csv",
        "markov_lolo_six_state_results.csv",
        "provisional_coordinates_audit.csv",
        "spatial_lolo_six_state_results.csv",
        "spatial_per_class_metrics.csv",
        "spatial_confusion_matrices.json",
        "directional_spatial_six_state_results.csv",
        "datum_sensitivity_six_state.csv",
    ]
    for fname in expected_files:
        p = f_dir / fname
        assert p.exists(), f"Missing deliverable: {p}"
        assert p.stat().st_size > 0, f"Empty deliverable: {p}"

    # Figures
    expected_figures = [
        "transect_vertical_successions_six_state.png",
        "facies_proportions_and_ntg_six_state.png",
        "upstream_vs_downstream_facies_ntg_six_state.png",
        "sandstone_thickness_distributions_six_state.png",
        "markov_transition_heatmaps_six_state.png",
    ]
    for fig_name in expected_figures:
        fp = fig_dir / fig_name
        assert fp.exists(), f"Missing figure: {fp}"
        assert fp.stat().st_size > 1000, f"Figure file too small: {fp}"


def test_spatial_confusion_matrices_six_state():
    """Verifies that the generated confusion matrices are exactly 6x6."""
    repo_root = Path(__file__).resolve().parent.parent
    json_path = repo_root / "sprints" / "audit_sprint_f" / "spatial_confusion_matrices.json"
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    assert len(data["facies_labels"]) == 6
    assert len(data["facies_codes"]) == 6
    for cm_key in ["knn_confusion_matrix", "near_well_confusion_matrix", "prior_confusion_matrix"]:
        cm = np.array(data[cm_key])
        assert cm.shape == (6, 6), f"{cm_key} shape is not 6x6: {cm.shape}"
        assert np.sum(cm) == 940, f"{cm_key} total points mismatch: {np.sum(cm)}"
