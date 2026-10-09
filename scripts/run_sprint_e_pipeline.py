"""
SMALT Sprint E Pipeline Runner: Revised Dataset Reconciliation & Spatial Prototype Readiness.

Executes:
1. Spatial metadata readiness audit table generation across all 12 lithologs.
2. Synthetic coordinate audit verification.
3. Provisional spatial baseline execution (Leave-One-Litholog-Out across eligible wells L2-L12).
4. Naive baseline comparisons (Training Prior Facies and Nearest-Well Vertical Profile).
5. Directional spatial experiments (Upstream <-> Downstream).
6. Stratigraphic vertical datum offset sensitivity analysis (+-5m, +-10m, +-20m).
7. Deliverable export to sprints/audit_sprint_e/.
"""

import sys
from pathlib import Path
import json
import numpy as np
import pandas as pd

# Add repo root to sys.path
repo_root = Path(__file__).resolve().parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from smalt.descriptive.analyzer import (
    LithologInspector,
    CANONICAL_FACIES_SCHEMA,
    PROVENANCE_METADATA,
)
from smalt.spatial.coordinates import (
    load_source_coordinates,
    get_spatial_litholog_subset,
    compute_interwell_horizontal_distances,
)
from smalt.spatial.baseline import ProvisionalSpatialValidator


def main():
    print("=" * 80)
    print("STARTING SMALT SPRINT E DATASET RECONCILIATION & SPATIAL READINESS PIPELINE")
    print("=" * 80)

    output_dir = repo_root / "sprints" / "audit_sprint_e"
    output_dir.mkdir(parents=True, exist_ok=True)

    inspector = LithologInspector(raw_dir=repo_root / "data" / "raw_lithologs")
    coords_df = load_source_coordinates()
    all_lids = [f"litholog{i}" for i in range(1, 13)]

    # =========================================================================
    # 1. Spatial Metadata Readiness Table
    # =========================================================================
    print("\n[1/5] Compiling spatial metadata readiness table across all 12 lithologs...")
    readiness_records = []

    for lid in all_lids:
        insp = inspector.inspect_litholog(lid)
        prov = PROVENANCE_METADATA[lid]
        coord_row = coords_df[coords_df["litholog_id"] == lid].iloc[0]

        has_coords = bool(coord_row["coordinates_available"])
        x_val = coord_row["x_m"] if has_coords else np.nan
        y_val = coord_row["y_m"] if has_coords else np.nan

        num = int(lid.replace("litholog", ""))
        raw_label = f"L{num}"

        # Quality flag
        gaps_cnt = len(insp["gaps"])
        overlaps_cnt = len(insp["overlaps"])
        if gaps_cnt > 0 and overlaps_cnt > 0:
            qual_flag = f"Gaps ({gaps_cnt}); Overlaps ({overlaps_cnt})"
        elif gaps_cnt > 0:
            qual_flag = f"Gaps ({gaps_cnt})"
        elif overlaps_cnt > 0:
            qual_flag = f"Overlaps ({overlaps_cnt})"
        else:
            qual_flag = "Continuous"

        # Spatial eligibility
        spatial_eligible = has_coords
        if not spatial_eligible:
            inclusion_reason = "EXCLUDED: No physical coordinates in source spreadsheet."
        else:
            inclusion_reason = "INCLUDED: Source coordinates available in Location_coordinates_lithologs.xlsx (provisional relative geometry)."

        record = {
            "litholog_id": lid,
            "raw_label": raw_label,
            "provenance_category": prov["provenance_category"],
            "observation_type": "core" if lid == "litholog12" else "outcrop",
            "digitized_interval_m": f"{insp['min_depth_m']:.1f}-{insp['max_depth_m']:.1f}m",
            "stratigraphic_span_m": insp["represented_span_m"],
            "valid_observations_m": insp["sum_thickness_m"],
            "quality_flag": qual_flag,
            "coordinates_available": has_coords,
            "x_m": round(x_val, 2) if pd.notna(x_val) else np.nan,
            "y_m": round(y_val, 2) if pd.notna(y_val) else np.nan,
            "coord_units_crs": "Meters; Local Cartesian (project-local CRS; unprojected geographic origin unverified)" if has_coords else "None",
            "spatial_position_semantics": "Provisional relative inter-well geometry; absolute georeferencing unverified" if has_coords else "None",
            "vertical_datum_status": "Unanchored; no common stratigraphic datum confirmed; equal raw depths do not correlate",
            "vertical_validation_eligible": True,
            "spatial_validation_eligible": spatial_eligible,
            "inclusion_exclusion_status": inclusion_reason,
        }
        readiness_records.append(record)

    df_readiness = pd.DataFrame(readiness_records)
    readiness_path = output_dir / "spatial_metadata_readiness.csv"
    df_readiness.to_csv(readiness_path, index=False)
    print(f"Saved spatial metadata readiness table to: {readiness_path}")

    # =========================================================================
    # 2. Inter-well Distances Table
    # =========================================================================
    print("\n[2/5] Computing inter-well horizontal distances...")
    dist_df = compute_interwell_horizontal_distances()
    dist_path = output_dir / "interwell_horizontal_distances.csv"
    dist_df.to_csv(dist_path, index=False)
    print(f"Saved inter-well distances to: {dist_path}")

    # =========================================================================
    # 3. Provisional Spatial LOLO Baseline (L2 - L12)
    # =========================================================================
    print("\n[3/5] Running leak-free provisional spatial LOLO baseline (L2 - L12)...")
    validator = ProvisionalSpatialValidator(
        inspector=inspector,
        vertical_reference="relative_to_base",
        vertical_weight=10.0,
    )

    eligible_lids = [lid for lid in all_lids if lid != "litholog1"]
    lolo_res = validator.run_spatial_lolo(eligible_log_ids=eligible_lids, k_neighbors=5)
    df_lolo = lolo_res["per_fold_table"]
    lolo_summary = lolo_res["aggregate_summary"]

    lolo_path = output_dir / "provisional_spatial_lolo_results.csv"
    df_lolo.to_csv(lolo_path, index=False)
    print(f"Saved provisional spatial LOLO results to: {lolo_path}")
    print(f"  Overall KNN Macro-F1: {lolo_summary['overall_knn_macro_f1']:.4f}")
    print(f"  Overall KNN Balanced Acc: {lolo_summary['overall_knn_balanced_acc']:.4f}")
    print(f"  Overall Nearest-Well Macro-F1: {lolo_summary['overall_near_well_macro_f1']:.4f}")
    print(f"  Overall Prior Facies Macro-F1: {lolo_summary['overall_prior_macro_f1']:.4f}")

    # Export aggregate summary JSON
    summary_path = output_dir / "provisional_spatial_lolo_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(lolo_summary, f, indent=2)

    # =========================================================================
    # 4. Directional Spatial Validation Experiments
    # =========================================================================
    print("\n[4/5] Running directional spatial validation experiments...")
    upstream_ids = ["litholog2", "litholog3", "litholog4", "litholog5", "litholog6", "litholog7", "litholog8", "litholog10"]
    downstream_ids = ["litholog9", "litholog11", "litholog12"]

    dir_up_to_down = validator.run_directional_spatial_experiment(
        training_group_ids=upstream_ids,
        target_group_ids=downstream_ids,
        experiment_name="Upstream_to_Downstream",
        k_neighbors=5,
    )

    dir_down_to_up = validator.run_directional_spatial_experiment(
        training_group_ids=downstream_ids,
        target_group_ids=upstream_ids,
        experiment_name="Downstream_to_Upstream",
        k_neighbors=5,
    )

    dir_records = [
        {
            "experiment": dir_up_to_down["experiment_name"],
            "training_wells": ",".join(upstream_ids),
            "target_wells": ",".join(downstream_ids),
            "total_test_points": dir_up_to_down["total_test_points"],
            "knn_macro_f1": dir_up_to_down["knn_macro_f1"],
            "knn_balanced_acc": dir_up_to_down["knn_balanced_accuracy"],
            "knn_accuracy": dir_up_to_down["knn_accuracy"],
            "prior_macro_f1": dir_up_to_down["prior_macro_f1"],
            "prior_balanced_acc": dir_up_to_down["prior_balanced_accuracy"],
            "prior_accuracy": dir_up_to_down["prior_accuracy"],
        },
        {
            "experiment": dir_down_to_up["experiment_name"],
            "training_wells": ",".join(downstream_ids),
            "target_wells": ",".join(upstream_ids),
            "total_test_points": dir_down_to_up["total_test_points"],
            "knn_macro_f1": dir_down_to_up["knn_macro_f1"],
            "knn_balanced_acc": dir_down_to_up["knn_balanced_accuracy"],
            "knn_accuracy": dir_down_to_up["knn_accuracy"],
            "prior_macro_f1": dir_down_to_up["prior_macro_f1"],
            "prior_balanced_acc": dir_down_to_up["prior_balanced_accuracy"],
            "prior_accuracy": dir_down_to_up["prior_accuracy"],
        },
    ]

    df_directional = pd.DataFrame(dir_records)
    dir_path = output_dir / "directional_spatial_results.csv"
    df_directional.to_csv(dir_path, index=False)
    print(f"Saved directional spatial results to: {dir_path}")

    # =========================================================================
    # 5. Datum Offset Sensitivity Analysis
    # =========================================================================
    print("\n[5/5] Running stratigraphic vertical datum offset sensitivity analysis...")
    shift_scenarios = [
        ("no_shift", 0.0),
        ("downstream_plus_5m", 5.0),
        ("downstream_minus_5m", -5.0),
        ("downstream_plus_10m", 10.0),
        ("downstream_minus_10m", -10.0),
        ("downstream_plus_20m", 20.0),
        ("downstream_minus_20m", -20.0),
    ]

    sensitivity_records = []
    for scen_name, shift_val in shift_scenarios:
        # Apply shift to downstream wells relative to upstream
        offsets = {lid: shift_val for lid in downstream_ids}
        res = validator.run_spatial_lolo(
            eligible_log_ids=eligible_lids,
            k_neighbors=5,
            datum_offsets=offsets,
        )
        agg = res["aggregate_summary"]
        sensitivity_records.append({
            "scenario": scen_name,
            "downstream_datum_offset_m": shift_val,
            "knn_macro_f1": agg["overall_knn_macro_f1"],
            "knn_balanced_acc": agg["overall_knn_balanced_acc"],
            "knn_accuracy": agg["overall_knn_accuracy"],
            "near_well_macro_f1": agg["overall_near_well_macro_f1"],
            "near_well_balanced_acc": agg["overall_near_well_balanced_acc"],
            "near_well_accuracy": agg["overall_near_well_accuracy"],
            "prior_macro_f1": agg["overall_prior_macro_f1"],
        })

    df_sensitivity = pd.DataFrame(sensitivity_records)
    sens_path = output_dir / "datum_offset_sensitivity_results.csv"
    df_sensitivity.to_csv(sens_path, index=False)
    print(f"Saved datum offset sensitivity results to: {sens_path}")

    print("\n" + "=" * 80)
    print("SPRINT E PIPELINE EXECUTION COMPLETE")
    print(f"All deliverables successfully exported to: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
