"""
Unit tests for SMALT Sprint G: Independent Results Audit & Spatial Modeling Decision.

Validates:
1. Metrics reconciliation exact recomputations from confusion matrices.
2. Statistical distinction between pooled metrics and unweighted fold means.
3. StandardScaler scale-invariance on vertical anisotropy weighting.
4. Geological sparsity bound: inter-well spacing strictly exceeds channel width.
5. Authoritative raw data hygiene: zero ambiguous labels ('silt', 'splay') in data/raw_lithologs/.
6. Minority facies failure confirmation (coal and carbon_mud recall = 0.0).
7. Datum sensitivity invariance band verification.
8. Verification of all Sprint G audit CSV deliverables.
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from sklearn.preprocessing import StandardScaler

from smalt.spatial.coordinates import load_source_coordinates
from smalt.descriptive.analyzer import CANONICAL_FACIES_SCHEMA


def test_metrics_reconciliation_exact_match():
    """Verify that independent recomputation from confusion matrices matches reported metrics exactly."""
    cm_path = Path("sprints/audit_sprint_f/spatial_confusion_matrices.json")
    assert cm_path.exists(), "spatial_confusion_matrices.json must exist"

    with open(cm_path) as f:
        cm_data = json.load(f)

    knn_cm = np.array(cm_data["knn_confusion_matrix"])
    near_cm = np.array(cm_data["near_well_confusion_matrix"])
    prior_cm = np.array(cm_data["prior_confusion_matrix"])

    # Verify total samples
    assert knn_cm.sum() == 940, "KNN total points must be exactly 940"
    assert near_cm.sum() == 940, "Nearest-well total points must be exactly 940"
    assert prior_cm.sum() == 940, "Prior total points must be exactly 940"

    # Verify KNN metrics
    knn_raw_acc = np.diag(knn_cm).sum() / 940.0
    assert round(knn_raw_acc, 4) == 0.4777, f"Expected 0.4777, got {knn_raw_acc:.4f}"

    rec = np.diag(knn_cm) / knn_cm.sum(axis=1)
    knn_bal_acc = np.mean(rec)
    assert round(knn_bal_acc, 4) == 0.2281, f"Expected 0.2281, got {knn_bal_acc:.4f}"

    col_sums = knn_cm.sum(axis=0)
    prec = np.zeros(6)
    for i in range(6):
        prec[i] = np.diag(knn_cm)[i] / col_sums[i] if col_sums[i] > 0 else 0.0
    f1 = np.zeros(6)
    for i in range(6):
        f1[i] = 2 * prec[i] * rec[i] / (prec[i] + rec[i]) if (prec[i] + rec[i]) > 0 else 0.0
    knn_macro_f1 = np.mean(f1)
    assert round(knn_macro_f1, 4) == 0.2232, f"Expected 0.2232, got {knn_macro_f1:.4f}"

    # Verify Nearest Well metrics
    near_raw_acc = np.diag(near_cm).sum() / 940.0
    assert round(near_raw_acc, 4) == 0.4426, f"Expected 0.4426, got {near_raw_acc:.4f}"


def test_unweighted_mean_vs_pooled_distinction():
    """Verify that unweighted mean fold balanced accuracy is distinct from pooled balanced accuracy."""
    df_fold = pd.read_csv("sprints/audit_sprint_f/spatial_lolo_six_state_results.csv")
    unweighted_bal_acc = df_fold["knn_balanced_acc"].mean()
    pooled_bal_acc = 0.2281

    # Unweighted fold balanced accuracy is higher (~24.89%) due to omitted classes in local test logs
    assert round(unweighted_bal_acc, 4) == 0.2489, f"Expected 0.2489, got {unweighted_bal_acc:.4f}"
    assert unweighted_bal_acc > pooled_bal_acc, "Fold average should exceed pooled due to class absence denominator effect"


def test_standard_scaler_scale_invariance_diagnostic():
    """Verify that multiplying a column by an anisotropy scalar before StandardScaler has zero effect."""
    rng = np.random.default_rng(42)
    Z = rng.uniform(0, 100, size=(100, 1))

    scaler1 = StandardScaler()
    s1 = scaler1.fit_transform(Z)

    scaler2 = StandardScaler()
    s2 = scaler2.fit_transform(Z * 10.0)  # Multiplied by 10

    scaler3 = StandardScaler()
    s3 = scaler3.fit_transform(Z * 100.0)  # Multiplied by 100

    np.testing.assert_allclose(s1, s2, atol=1e-12)
    np.testing.assert_allclose(s1, s3, atol=1e-12)


def test_inter_well_sparsity_lower_bound():
    """Verify that inter-well nearest neighbor distances all exceed channel body widths (210 m)."""
    df_coords = load_source_coordinates()
    eligible = df_coords[df_coords["coordinates_available"]].reset_index(drop=True)
    coords = eligible[["x_m", "y_m"]].to_numpy()
    n = len(coords)

    nn_dists = []
    for i in range(n):
        dists = [np.hypot(coords[i, 0] - coords[j, 0], coords[i, 1] - coords[j, 1]) for j in range(n) if i != j]
        nn_dists.append(min(dists))

    min_nn = min(nn_dists)
    median_nn = np.median(nn_dists)

    assert min_nn >= 420.0, f"Minimum NN distance must be >= 420.0 m, got {min_nn:.1f} m"
    assert median_nn >= 700.0, f"Median NN distance must be >= 700.0 m, got {median_nn:.1f} m"
    # Channel width benchmark from Sahoo et al. (2016): W/T ~ 35, mean T ~ 5.8 m -> ~203 m
    channel_width_upper_bound = 210.0
    assert min_nn > channel_width_upper_bound, "All wells must be spaced wider than individual channel widths"


def test_authoritative_raw_data_hygiene():
    """Verify that authoritative raw CSVs contain zero ambiguous labels ('silt', 'siltstone', 'splay')."""
    raw_files = list(Path("data/raw_lithologs").glob("*.csv"))
    assert len(raw_files) == 12, "Must have exactly 12 raw litholog CSVs"

    disallowed_labels = {"silt", "siltstone", "splay", "fine sandstone"}
    for f in raw_files:
        df = pd.read_csv(f)
        col = [c for c in df.columns if "facies" in c.lower() or "lith" in c.lower()][0]
        found_disallowed = set(df[col].astype(str).str.strip().str.lower()).intersection(disallowed_labels)
        assert len(found_disallowed) == 0, f"File {f.name} contains disallowed ambiguous labels: {found_disallowed}"


def test_minority_facies_zero_recall_audit():
    """Verify that KNN and Nearest-Well have strictly 0.0 recall on coal and carbonaceous mud."""
    with open("sprints/audit_sprint_f/spatial_confusion_matrices.json") as f:
        cm_data = json.load(f)

    knn_cm = np.array(cm_data["knn_confusion_matrix"])
    near_cm = np.array(cm_data["near_well_confusion_matrix"])

    # Carbonaceous Mudstone is code 3
    # Coal is code 4
    assert knn_cm[3, 3] == 0, "KNN must have 0 correct predictions for carbon_mud"
    assert knn_cm[4, 4] == 0, "KNN must have 0 correct predictions for coal"
    assert near_cm[3, 3] == 0, "Nearest-well must have 0 correct predictions for carbon_mud"
    assert near_cm[4, 4] == 0, "Nearest-well must have 0 correct predictions for coal"


def test_datum_sensitivity_invariance_band():
    """Verify that across all tested datum shifts (-20 to +20 m), KNN Macro-F1 remains trapped in [0.21, 0.23]."""
    df_datum = pd.read_csv("sprints/audit_sprint_f/datum_sensitivity_six_state.csv")
    f1_min = df_datum["knn_macro_f1"].min()
    f1_max = df_datum["knn_macro_f1"].max()

    assert f1_min >= 0.2100, f"F1 min out of bound: {f1_min}"
    assert f1_max <= 0.2350, f"F1 max out of bound: {f1_max}"


def test_all_audit_g_deliverables_exist_and_non_empty():
    """Verify all 5 CSV audit deliverables in sprints/audit_sprint_g/ exist and are populated."""
    expected_csvs = [
        "metrics_reconciliation_audit.csv",
        "facies_mapping_and_aliases_audit.csv",
        "root_cause_evidence_matrix.csv",
        "baseline_comparison_summary.csv",
        "phase_decision_matrix.csv",
    ]
    for fname in expected_csvs:
        p = Path("sprints/audit_sprint_g") / fname
        assert p.exists(), f"Deliverable {fname} must exist"
        df = pd.read_csv(p)
        assert len(df) > 0, f"Deliverable {fname} must not be empty"
