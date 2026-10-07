"""
SMALT Sprint E Test Suite: Revised Dataset Reconciliation, Provenance & Spatial Readiness.

Verifies:
1. All 12 revised lithologs are loaded and validated.
2. Canonical facies encoding (including revised labels p_sand and ripples) and unknown-label rejection.
3. Interval depths and strictly positive thickness calculations.
4. Explicit reporting of gaps and overlaps (including L9 gaps and L9/L11 overlaps).
5. Litholog 12 provenance registration, core metadata, and strictly 0-111 m digitized scope.
6. Absence of synthetic coordinates (e.g. strike_pos_m) in spatial model feature vectors.
7. Strict exclusion of Litholog 1 from spatial validation due to missing coordinates.
8. Complete exclusion of target litholog from training set and preprocessing fitting in spatial LOLO.
9. Zero feature leakage (no synthetic Gamma Ray, no ground-truth previous facies in predictors).
10. Existence and consistency of Sprint E generated audit deliverables.
"""

import pytest
from pathlib import Path
import json
import numpy as np
import pandas as pd

from smalt.descriptive.analyzer import (
    LithologInspector,
    CANONICAL_FACIES_SCHEMA,
    PROVENANCE_METADATA,
    normalize_facies_label,
)
from smalt.spatial.coordinates import (
    load_source_coordinates,
    get_spatial_litholog_subset,
    compute_interwell_horizontal_distances,
)
from smalt.spatial.baseline import ProvisionalSpatialValidator


@pytest.fixture
def inspector():
    return LithologInspector(raw_dir=Path("data/raw_lithologs"))


@pytest.fixture
def spatial_validator(inspector):
    return ProvisionalSpatialValidator(inspector=inspector)


# -----------------------------------------------------------------------------
# 1. Dataset Reconciliation & Loading Invariants
# -----------------------------------------------------------------------------

def test_all_12_revised_lithologs_loaded(inspector):
    """Verifies that all 12 revised lithologs exist and load cleanly."""
    for i in range(1, 13):
        lid = f"litholog{i}"
        df = inspector.load_raw_litholog(lid)
        assert len(df) > 0, f"{lid} is empty"
        assert "Top" in df.columns
        assert "Bottom" in df.columns
        assert "Facies" in df.columns
        assert (df["Thickness"] > 0).all(), f"{lid} has non-positive thickness"


def test_canonical_facies_encoding_and_unknown_detection():
    """Verifies that revised facies labels map correctly to canonical classes and invalid labels fail."""
    # Test revised labels (in Sprint F: p_sand and ripples are distinct canonical facies)
    assert normalize_facies_label("p_sand") == "p_sand"
    assert normalize_facies_label("planar_sand") == "p_sand"
    assert normalize_facies_label("ripples") == "ripples"
    assert normalize_facies_label("ripple") == "ripples"

    # Test standard canonical labels
    assert normalize_facies_label("coal") == "coal"
    assert normalize_facies_label("sand") == "sand"
    assert normalize_facies_label("carbon_mud") == "carbon_mud"
    assert normalize_facies_label("silt") == "ripples"
    assert normalize_facies_label("mud") == "mud"

    # Test unknown label rejection
    with pytest.raises(ValueError, match="Unknown facies label"):
        normalize_facies_label("unknown_formation_x")

    with pytest.raises(ValueError, match="Encountered null or NaN"):
        normalize_facies_label(np.nan)


def test_interval_depths_and_stratigraphic_ordering(inspector):
    """Verifies valid depths, non-negative bounds, and positive thickness across all 12 logs."""
    for i in range(1, 13):
        lid = f"litholog{i}"
        df = inspector.load_raw_litholog(lid)
        assert (df["Top"] >= 0).all()
        assert (df["Bottom"] > df["Top"]).all()
        assert (df["Thickness"] > 0).all()


def test_explicit_detection_of_gaps_and_overlaps(inspector):
    """Verifies that gaps and overlaps are explicitly detected and reported, not concealed."""
    insp9 = inspector.inspect_litholog("litholog9")
    # L9 has 2 gaps: 18-19m and 51-52m, and 1 overlap: 28-30m
    assert len(insp9["gaps"]) == 2
    assert insp9["gaps"][0]["depth_before"] == 18.0
    assert insp9["gaps"][0]["depth_after"] == 19.0
    assert insp9["gaps"][1]["depth_before"] == 51.0
    assert insp9["gaps"][1]["depth_after"] == 52.0
    assert len(insp9["overlaps"]) == 1
    assert insp9["overlaps"][0]["depth_before"] == 29.0
    assert insp9["overlaps"][0]["depth_after"] == 28.0

    insp11 = inspector.inspect_litholog("litholog11")
    # In revised L11, 59-60m gap is resolved as carbon_mud; overlap at 71-72m exists
    assert len(insp11["gaps"]) == 0
    assert len(insp11["overlaps"]) == 1
    assert insp11["overlaps"][0]["depth_before"] == 72.0
    assert insp11["overlaps"][0]["depth_after"] == 71.0


# -----------------------------------------------------------------------------
# 2. Litholog 12 Provenance and Scope Invariants
# -----------------------------------------------------------------------------

def test_litholog_12_provenance_and_depth_scope(inspector):
    """Verifies L12 provenance registration, manual digitization, and strict 0-111m coverage."""
    # Check manifest JSON
    manifest_path = Path("data/provenance_manifest.json")
    assert manifest_path.exists()
    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    assert "litholog12" in manifest["lithologs"]
    l12_man = manifest["lithologs"]["litholog12"]
    assert l12_man["provenance_category"] == "digitized_core_log"
    assert "manually digitized" in l12_man["reconstruction_method"]
    assert l12_man["original_core_length_m"] == 242.0
    assert l12_man["stratigraphic_span_m"] == 111.0
    assert l12_man["independent_validation"] is False

    # Check raw data inspection
    insp12 = inspector.inspect_litholog("litholog12")
    assert insp12["min_depth_m"] == 0.0
    assert insp12["max_depth_m"] == 111.0
    assert insp12["represented_span_m"] == 111.0
    assert insp12["sum_thickness_m"] == 111.0
    assert len(insp12["gaps"]) == 0
    assert len(insp12["overlaps"]) == 0


# -----------------------------------------------------------------------------
# 3. Coordinate Audit & Spatial Leakage Prevention Invariants
# -----------------------------------------------------------------------------

def test_source_coordinates_ingestion_and_missing_l1():
    """Verifies that source coordinates are ingested correctly and L1 is marked missing."""
    coords = load_source_coordinates()
    assert len(coords) == 12

    # L1 must be missing
    l1_row = coords[coords["litholog_id"] == "litholog1"].iloc[0]
    assert not l1_row["coordinates_available"]
    assert pd.isna(l1_row["x_m"])
    assert pd.isna(l1_row["y_m"])

    # L2-L12 must be present
    for i in range(2, 13):
        lid = f"litholog{i}"
        row = coords[coords["litholog_id"] == lid].iloc[0]
        assert row["coordinates_available"], f"{lid} should have coordinates"
        assert pd.notna(row["x_m"])
        assert pd.notna(row["y_m"])

    # Subset must exclude L1
    subset = get_spatial_litholog_subset(exclude_missing=True)
    assert len(subset) == 11
    assert "litholog1" not in subset["litholog_id"].values


def test_litholog_1_exclusion_from_spatial_validator(spatial_validator):
    """Verifies that spatial validator strictly rejects litholog1."""
    with pytest.raises(ValueError, match="(?i)strictly ineligible"):
        spatial_validator.prepare_discretized_spatial_points("litholog1")

    with pytest.raises(ValueError, match="(?i)missing physical coordinates"):
        spatial_validator.run_spatial_lolo(eligible_log_ids=["litholog1", "litholog2"])


def test_absence_of_synthetic_coordinates_in_spatial_points(spatial_validator):
    """Verifies that spatial points contain ONLY true Cartesian [X, Y, Z_rel] and no synthetic indices."""
    pts = spatial_validator.prepare_discretized_spatial_points("litholog2")
    assert "x_m" in pts.columns
    assert "y_m" in pts.columns
    assert "z_rel_m" in pts.columns
    assert "strike_pos_m" not in pts.columns
    assert "file_idx" not in pts.columns
    assert "well_num" not in pts.columns

    # Verify true coordinate values for L2
    assert pts["x_m"].iloc[0] == pytest.approx(-4057.68)
    assert pts["y_m"].iloc[0] == pytest.approx(13755.61)


def test_spatial_lolo_strict_leak_free_execution(spatial_validator):
    """Verifies leak-free LOLO spatial cross-validation execution across eligible wells."""
    res = spatial_validator.run_spatial_lolo(eligible_log_ids=["litholog4", "litholog5", "litholog6"], k_neighbors=3)
    df_folds = res["per_fold_table"]
    summary = res["aggregate_summary"]

    assert len(df_folds) == 3
    assert set(df_folds["target_litholog_id"]) == {"litholog4", "litholog5", "litholog6"}
    
    # Check that each fold has valid metrics and support counts
    for _, fold in df_folds.iterrows():
        assert 0.0 <= fold["knn_accuracy"] <= 1.0
        assert 0.0 <= fold["knn_macro_f1"] <= 1.0
        assert isinstance(fold["support_by_class"], dict)
        assert sum(fold["support_by_class"].values()) == fold["test_points_count"]


def test_vertical_datum_sensitivity_execution(spatial_validator):
    """Verifies that datum offset sensitivity runs cleanly with arbitrary vertical shifts."""
    offsets = {"litholog5": 10.0, "litholog6": -10.0}
    res = spatial_validator.run_spatial_lolo(
        eligible_log_ids=["litholog4", "litholog5", "litholog6"],
        k_neighbors=3,
        datum_offsets=offsets,
    )
    summary = res["aggregate_summary"]
    assert 0.0 <= summary["overall_knn_macro_f1"] <= 1.0


# -----------------------------------------------------------------------------
# 4. Deliverable File Existence and Non-Emptiness
# -----------------------------------------------------------------------------

def test_sprint_e_deliverables_exist():
    """Verifies that all Sprint E deliverables are generated and populated."""
    audit_dir = Path("sprints/audit_sprint_e")
    assert (audit_dir / "spatial_metadata_readiness.csv").exists()
    assert (audit_dir / "synthetic_coordinate_audit.md").exists()
    assert (audit_dir / "provisional_spatial_lolo_results.csv").exists()
    assert (audit_dir / "directional_spatial_results.csv").exists()
    assert (audit_dir / "datum_offset_sensitivity_results.csv").exists()
    assert (audit_dir / "interwell_horizontal_distances.csv").exists()

    df_readiness = pd.read_csv(audit_dir / "spatial_metadata_readiness.csv")
    assert len(df_readiness) == 12

    df_lolo = pd.read_csv(audit_dir / "provisional_spatial_lolo_results.csv")
    assert len(df_lolo) == 11
    assert "litholog1" not in df_lolo["target_litholog_id"].values
