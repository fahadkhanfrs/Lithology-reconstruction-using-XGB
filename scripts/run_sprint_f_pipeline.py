"""
SMALT Sprint F Master Pipeline Runner: Six-State Migration & Spatial Baseline Revalidation.

Executes:
1. Serialized processed dataset export with 6-state canonical schema.
2. 6-state facies schema specification artifact export.
3. Descriptive statistics, continuous vs discretized comparison, and bed statistics across all 12 logs.
4. Upstream vs Downstream group comparative analysis.
5. 6x6 Regular and Embedded vertical Markov transition matrices & stationary vectors.
6. 1D Markov Leave-One-Litholog-Out (LOLO) cross-validation and gap sensitivity.
7. Provisional coordinate audit tracing (L1 excluded, CRS unverified caveats).
8. Spatial LOLO cross-validation (L2 - L12) comparing Spatial 3D KNN, Nearest-Well, and Prior baselines.
9. Per-class precision, recall, F1, support, and 6x6 confusion matrices.
10. Directional cross-group spatial experiments (Upstream <-> Downstream).
11. Stratigraphic vertical datum sensitivity analysis.
12. Publication-quality figure generation.
13. Complete deliverable artifact export to sprints/audit_sprint_f/.
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

from data.loader import LithologLoader
from smalt.descriptive.analyzer import (
    LithologInspector,
    CANONICAL_FACIES_SCHEMA,
    PROVENANCE_METADATA,
)
from smalt.descriptive.plotting import (
    plot_vertical_successions,
    plot_facies_proportions,
    plot_upstream_vs_downstream,
    plot_sandstone_thickness_distributions,
    plot_transition_matrix_heatmaps,
)
from smalt.validation.markov_lolo import MarkovLOLOValidator
from smalt.spatial.coordinates import (
    load_source_coordinates,
    get_spatial_litholog_subset,
    compute_interwell_horizontal_distances,
)
from smalt.spatial.baseline import ProvisionalSpatialValidator


def main():
    print("=" * 80)
    print("STARTING SMALT SPRINT F SIX-STATE MIGRATION & SPATIAL BASELINE REVALIDATION")
    print("=" * 80)

    output_dir = repo_root / "sprints" / "audit_sprint_f"
    figures_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    inspector = LithologInspector(raw_dir=repo_root / "data" / "raw_lithologs")
    validator = MarkovLOLOValidator(inspector=inspector, smoothing_alpha=0.1)
    coords_df = load_source_coordinates()
    all_lids = [f"litholog{i}" for i in range(1, 13)]

    # =========================================================================
    # 1. Regenerate Unified Processed Dataset
    # =========================================================================
    print("\n[1/10] Exporting unified processed datasets with 6-state schema...")
    loader = LithologLoader(raw_dir=repo_root / "data" / "raw_lithologs", output_dir=repo_root / "data" / "processed")
    parquet_path, csv_path = loader.export()
    print(f"  -> Exported to {parquet_path} and {csv_path}")

    # =========================================================================
    # 2. Export 6-State Facies Schema Reference
    # =========================================================================
    print("\n[2/10] Exporting 6-state facies schema reference...")
    schema_records = []
    for fkey, meta in CANONICAL_FACIES_SCHEMA.items():
        schema_records.append({
            "code": meta["code"],
            "facies_number": meta["facies_number"],
            "technical_identifier": fkey,
            "canonical_name": meta["canonical_name"],
            "lithology_type": meta["lithology_type"],
            "color_hex": meta["color_hex"],
            "sahoo_2016_analogy": meta["geological_caveat"],
        })
    df_schema = pd.DataFrame(schema_records).sort_values("code")
    schema_csv_path = output_dir / "facies_schema_six_state.csv"
    df_schema.to_csv(schema_csv_path, index=False)
    print(f"  -> Exported to {schema_csv_path}")

    # =========================================================================
    # 3. Descriptive Lithology & Continuous vs Discretized Analysis
    # =========================================================================
    print("\n[3/10] Inspecting raw intervals and computing 6-state descriptive statistics...")
    per_litholog_records = []
    for lid in all_lids:
        insp = inspector.inspect_litholog(lid)
        comp = inspector.compare_continuous_vs_discretized(lid)

        gaps_str = "; ".join([f"{g['depth_before']}-{g['depth_after']}m ({g['gap_size_m']}m)" for g in insp["gaps"]]) if insp["gaps"] else "None"
        overlaps_str = "; ".join([f"{o['depth_before']}-{o['depth_after']}m ({o['overlap_size_m']}m)" for o in insp["overlaps"]]) if insp["overlaps"] else "None"

        record = {
            "litholog_id": lid,
            "provenance_category": insp["provenance_category"],
            "provenance_description": insp["provenance_description"],
            "independent_validation": insp["independent_validation"],
            "group": insp["group"],
            "coordinates_available": insp["coordinates_available"],
            "coordinates_status": insp["coordinates_status"],
            "num_raw_intervals": insp["num_intervals"],
            "min_depth_m": insp["min_depth_m"],
            "max_depth_m": insp["max_depth_m"],
            "represented_span_m": insp["represented_span_m"],
            "continuous_thickness_m": insp["sum_thickness_m"],
            "span_thickness_diff_m": insp["span_thickness_diff_m"],
            "has_gaps": insp["has_gaps"],
            "has_overlaps": insp["has_overlaps"],
            "gaps_summary": gaps_str,
            "overlaps_summary": overlaps_str,
            # 6-state continuous thicknesses
            "sand_thickness_m": insp["facies_thickness_m"]["sand"],
            "p_sand_thickness_m": insp["facies_thickness_m"]["p_sand"],
            "ripples_thickness_m": insp["facies_thickness_m"]["ripples"],
            "carbon_mud_thickness_m": insp["facies_thickness_m"]["carbon_mud"],
            "coal_thickness_m": insp["facies_thickness_m"]["coal"],
            "mud_thickness_m": insp["facies_thickness_m"]["mud"],
            # 6-state continuous proportions
            "sand_prop_cont": insp["facies_proportions"]["sand"],
            "p_sand_prop_cont": insp["facies_proportions"]["p_sand"],
            "ripples_prop_cont": insp["facies_proportions"]["ripples"],
            "carbon_mud_prop_cont": insp["facies_proportions"]["carbon_mud"],
            "coal_prop_cont": insp["facies_proportions"]["coal"],
            "mud_prop_cont": insp["facies_proportions"]["mud"],
            # Explicit Net-to-Gross definitions
            "ntg_channel": insp["ntg_channel"],
            "ntg_planar": insp["ntg_planar"],
            "ntg_net_sand": insp["ntg_pure"],
            "ntg_coarse": insp["ntg_coarse"],
            # Sandstone bed metrics
            "channel_sand_beds_count": insp["sand_beds_count"],
            "channel_sand_mean_bed_thick_m": insp["sand_mean_bed_thickness_m"],
            "channel_sand_max_bed_thick_m": insp["sand_max_bed_thickness_m"],
            "planar_sand_beds_count": insp["p_sand_beds_count"],
            "planar_sand_mean_bed_thick_m": insp["p_sand_mean_bed_thickness_m"],
            "planar_sand_max_bed_thick_m": insp["p_sand_max_bed_thickness_m"],
            "all_sand_beds_count": insp["all_sand_beds_count"],
            "all_sand_mean_bed_thick_m": insp["all_sand_mean_bed_thickness_m"],
            "all_sand_max_bed_thick_m": insp["all_sand_max_bed_thickness_m"],
            # Coal bed metrics
            "coal_beds_count": insp["coal_beds_count"],
            "coal_mean_bed_thick_m": insp["coal_mean_bed_thickness_m"],
            "coal_max_bed_thick_m": insp["coal_max_bed_thickness_m"],
            # Discretized 1m comparison
            "discretized_samples_1m": comp["discretized_samples_1m"],
            "discretized_ntg_net_sand": comp["discretized_ntg_pure"],
            "delta_ntg_net_sand": comp["delta_ntg_pure"],
        }
        per_litholog_records.append(record)

    df_per_litholog = pd.DataFrame(per_litholog_records)
    per_litholog_csv = output_dir / "per_litholog_six_state_statistics.csv"
    df_per_litholog.to_csv(per_litholog_csv, index=False)
    print(f"  -> Exported to {per_litholog_csv}")

    # =========================================================================
    # 4. Group Aggregation (Upstream vs Downstream)
    # =========================================================================
    print("\n[4/10] Computing group-level 6-state aggregations...")
    upstream_ids = [f"litholog{i}" for i in [2, 3, 4, 5, 6, 7, 8, 10]]
    downstream_outcrop_ids = ["litholog9", "litholog11"]
    downstream_all_ids = ["litholog9", "litholog11", "litholog12"]
    all_12_ids = all_lids

    group_defs = [
        ("All 12 Lithologs", all_12_ids),
        ("Upstream (L2-L8, L10)", upstream_ids),
        ("Downstream Outcrop (L9, L11)", downstream_outcrop_ids),
        ("Downstream Composite (L9, L11, L12)", downstream_all_ids),
    ]

    group_records = []
    for gname, lids in group_defs:
        tot_thick = 0.0
        fac_thick = {k: 0.0 for k in CANONICAL_FACIES_SCHEMA.keys()}
        tot_intervals = 0
        tot_ch_beds = 0
        tot_pl_beds = 0
        tot_coal_beds = 0

        for lid in lids:
            insp = inspector.inspect_litholog(lid)
            tot_thick += insp["sum_thickness_m"]
            tot_intervals += insp["num_intervals"]
            tot_ch_beds += insp["sand_beds_count"]
            tot_pl_beds += insp["p_sand_beds_count"]
            tot_coal_beds += insp["coal_beds_count"]
            for k in CANONICAL_FACIES_SCHEMA.keys():
                fac_thick[k] += insp["facies_thickness_m"][k]

        grec = {
            "group_name": gname,
            "num_lithologs": len(lids),
            "litholog_ids": ", ".join(lids),
            "total_intervals": tot_intervals,
            "total_thickness_m": round(tot_thick, 3),
            "sand_thickness_m": round(fac_thick["sand"], 3),
            "p_sand_thickness_m": round(fac_thick["p_sand"], 3),
            "ripples_thickness_m": round(fac_thick["ripples"], 3),
            "carbon_mud_thickness_m": round(fac_thick["carbon_mud"], 3),
            "coal_thickness_m": round(fac_thick["coal"], 3),
            "mud_thickness_m": round(fac_thick["mud"], 3),
            "sand_prop": round(fac_thick["sand"] / tot_thick, 4),
            "p_sand_prop": round(fac_thick["p_sand"] / tot_thick, 4),
            "ripples_prop": round(fac_thick["ripples"] / tot_thick, 4),
            "carbon_mud_prop": round(fac_thick["carbon_mud"] / tot_thick, 4),
            "coal_prop": round(fac_thick["coal"] / tot_thick, 4),
            "mud_prop": round(fac_thick["mud"] / tot_thick, 4),
            "ntg_channel": round(fac_thick["sand"] / tot_thick, 4),
            "ntg_planar": round(fac_thick["p_sand"] / tot_thick, 4),
            "ntg_net_sand": round((fac_thick["sand"] + fac_thick["p_sand"]) / tot_thick, 4),
            "ntg_coarse": round((fac_thick["sand"] + fac_thick["p_sand"] + fac_thick["ripples"]) / tot_thick, 4),
            "channel_sand_beds_total": tot_ch_beds,
            "planar_sand_beds_total": tot_pl_beds,
            "coal_beds_total": tot_coal_beds,
        }
        group_records.append(grec)

    df_groups = pd.DataFrame(group_records)
    group_csv = output_dir / "group_six_state_statistics.csv"
    df_groups.to_csv(group_csv, index=False)
    print(f"  -> Exported to {group_csv}")

    # =========================================================================
    # 5. 6x6 Vertical Stratigraphic Markov Transition Matrices
    # =========================================================================
    print("\n[5/10] Estimating 6x6 Regular & Embedded Markov transition matrices...")
    matrix_records = []
    names = [meta["canonical_name"] for meta in CANONICAL_FACIES_SCHEMA.values()]

    for gname, lids in [("All_12_Lithologs", all_12_ids), ("Upstream", upstream_ids), ("Downstream_Composite", downstream_all_ids)]:
        for embedded in [False, True]:
            use_disc = not embedded
            P, N, pi = inspector.compute_vertical_transitions(
                litholog_ids=lids,
                embedded=embedded,
                use_discretized=use_disc,
                smoothing_alpha=0.1,
            )
            for i, from_name in enumerate(names):
                for j, to_name in enumerate(names):
                    matrix_records.append({
                        "dataset_scope": gname,
                        "model_type": "Embedded (Continuous Beds)" if embedded else "Regular (1m Discretized)",
                        "embedded": embedded,
                        "from_state_code": i,
                        "from_state_name": from_name,
                        "to_state_code": j,
                        "to_state_name": to_name,
                        "transition_count_N": float(N[i, j]),
                        "transition_prob_P": round(float(P[i, j]), 5),
                        "from_stationary_pi": round(float(pi[i]), 5),
                        "to_stationary_pi": round(float(pi[j]), 5),
                    })

    df_matrices = pd.DataFrame(matrix_records)
    matrices_csv = output_dir / "markov_transition_matrices_six_state.csv"
    df_matrices.to_csv(matrices_csv, index=False)
    print(f"  -> Exported to {matrices_csv}")

    # =========================================================================
    # 6. 1D Markov Leave-One-Litholog-Out (LOLO) Cross-Validation
    # =========================================================================
    print("\n[6/10] Running 1D Markov LOLO cross-validation across all 12 lithologs...")
    lolo_markov_records = []
    for embedded in [False, True]:
        use_disc = not embedded
        df_lolo_m = validator.run_lolo_cross_validation(
            eligible_log_ids=all_12_ids,
            embedded=embedded,
            use_discretized=use_disc,
        )
        lolo_markov_records.append(df_lolo_m)

    df_lolo_markov = pd.concat(lolo_markov_records, ignore_index=True)
    lolo_markov_csv = output_dir / "markov_lolo_six_state_results.csv"
    df_lolo_markov.to_csv(lolo_markov_csv, index=False)
    print(f"  -> Exported to {lolo_markov_csv}")

    # =========================================================================
    # 7. Spatial Coordinates Audit & Ingestion
    # =========================================================================
    print("\n[7/10] Compiling provisional coordinate audit table...")
    coord_audit_records = []
    for lid in all_lids:
        insp = inspector.inspect_litholog(lid)
        prov = PROVENANCE_METADATA[lid]
        crow = coords_df[coords_df["litholog_id"] == lid].iloc[0]
        has_c = bool(crow["coordinates_available"])

        coord_audit_records.append({
            "litholog_id": lid,
            "raw_label": crow["raw_label"],
            "provenance_category": prov["provenance_category"],
            "group": prov["group"],
            "digitized_span_m": insp["represented_span_m"],
            "coordinates_available": has_c,
            "x_m": round(crow["x_m"], 2) if has_c else np.nan,
            "y_m": round(crow["y_m"], 2) if has_c else np.nan,
            "crs_status": "Project-local Cartesian (unverified geographic origin/units)" if has_c else "None",
            "spatial_modeling_eligibility": "ELIGIBLE" if has_c else "EXCLUDED (no coordinates)",
            "notes": crow["notes"],
        })
    df_coord_audit = pd.DataFrame(coord_audit_records)
    coord_audit_csv = output_dir / "provisional_coordinates_audit.csv"
    df_coord_audit.to_csv(coord_audit_csv, index=False)
    print(f"  -> Exported to {coord_audit_csv}")

    # =========================================================================
    # 8. Spatial LOLO Cross-Validation (6-State Schema, L2-L12)
    # =========================================================================
    print("\n[8/10] Running Spatial LOLO cross-validation (KNN, Nearest-Well, Prior baselines)...")
    spatial_val = ProvisionalSpatialValidator(inspector=inspector, vertical_reference="relative_to_base", vertical_weight=10.0)
    spatial_results = spatial_val.run_spatial_lolo(k_neighbors=5)

    df_spatial_folds = spatial_results["per_fold_table"]
    spatial_summary = spatial_results["aggregate_summary"]

    spatial_folds_csv = output_dir / "spatial_lolo_six_state_results.csv"
    df_spatial_folds.to_csv(spatial_folds_csv, index=False)
    print(f"  -> Exported per-fold table to {spatial_folds_csv}")

    # Save confusion matrices to JSON
    cm_dict = {
        "facies_labels": [meta["canonical_name"] for meta in CANONICAL_FACIES_SCHEMA.values()],
        "facies_codes": [meta["code"] for meta in CANONICAL_FACIES_SCHEMA.values()],
        "knn_confusion_matrix": spatial_summary["knn_confusion_matrix"],
        "near_well_confusion_matrix": spatial_summary["near_well_confusion_matrix"],
        "prior_confusion_matrix": spatial_summary["prior_confusion_matrix"],
    }
    cm_json_path = output_dir / "spatial_confusion_matrices.json"
    with open(cm_json_path, "w", encoding="utf-8") as f:
        json.dump(cm_dict, f, indent=2)
    print(f"  -> Exported confusion matrices to {cm_json_path}")

    # Compile Per-Class Metrics Table
    per_class_records = []
    for model_name, key in [("Spatial 3D KNN", "knn_per_class"), ("Nearest-Well Profile", "near_well_per_class"), ("Training Prior Facies", "prior_per_class")]:
        pdict = spatial_summary[key]
        for fname, metrics in pdict.items():
            per_class_records.append({
                "model_name": model_name,
                "facies_name": fname,
                "facies_code": metrics["code"],
                "precision": metrics["precision"],
                "recall": metrics["recall"],
                "f1_score": metrics["f1"],
                "test_support": metrics["support"],
            })
    df_per_class = pd.DataFrame(per_class_records)
    per_class_csv = output_dir / "spatial_per_class_metrics.csv"
    df_per_class.to_csv(per_class_csv, index=False)
    print(f"  -> Exported per-class metrics to {per_class_csv}")

    # =========================================================================
    # 9. Directional Spatial Prediction Experiments & Datum Sensitivity
    # =========================================================================
    print("\n[9/10] Running directional cross-group experiments & stratigraphic datum sensitivity...")
    # Directional Experiments
    dir_records = []
    # Upstream -> Downstream (eligible wells only: L1 excluded)
    exp1 = spatial_val.run_directional_spatial_experiment(
        training_group_ids=upstream_ids,
        target_group_ids=downstream_all_ids,
        experiment_name="Upstream_to_Downstream",
        k_neighbors=5,
    )
    # Downstream -> Upstream
    exp2 = spatial_val.run_directional_spatial_experiment(
        training_group_ids=downstream_all_ids,
        target_group_ids=upstream_ids,
        experiment_name="Downstream_to_Upstream",
        k_neighbors=5,
    )
    for exp in [exp1, exp2]:
        dir_records.append({
            "experiment_name": exp["experiment_name"],
            "training_wells_count": exp["training_wells_count"],
            "target_wells_count": exp["target_wells_count"],
            "total_test_points": exp["total_test_points"],
            "knn_macro_f1": exp["knn_macro_f1"],
            "knn_balanced_accuracy": exp["knn_balanced_accuracy"],
            "knn_accuracy": exp["knn_accuracy"],
            "prior_macro_f1": exp["prior_macro_f1"],
            "prior_balanced_accuracy": exp["prior_balanced_accuracy"],
            "prior_accuracy": exp["prior_accuracy"],
        })
    df_directional = pd.DataFrame(dir_records)
    directional_csv = output_dir / "directional_spatial_six_state_results.csv"
    df_directional.to_csv(directional_csv, index=False)
    print(f"  -> Exported directional results to {directional_csv}")

    # Stratigraphic Datum Sensitivity
    offsets = [-20.0, -10.0, -5.0, 0.0, 5.0, 10.0, 20.0]
    sens_records = []
    for ref in ["relative_to_base", "relative_to_top"]:
        for d_off in offsets:
            validator_sens = ProvisionalSpatialValidator(
                inspector=inspector,
                vertical_reference=ref,
                vertical_weight=10.0,
            )
            dummy_offsets = {lid: (d_off if "downstream" in PROVENANCE_METADATA[lid]["group"] else 0.0) for lid in coords_df[coords_df["coordinates_available"]]["litholog_id"]}
            res_sens = validator_sens.run_spatial_lolo(datum_offsets=dummy_offsets)
            s_summ = res_sens["aggregate_summary"]
            sens_records.append({
                "vertical_reference": ref,
                "downstream_datum_offset_m": d_off,
                "knn_macro_f1": s_summ["overall_knn_macro_f1"],
                "knn_balanced_accuracy": s_summ["overall_knn_balanced_acc"],
                "knn_accuracy": s_summ["overall_knn_accuracy"],
                "near_well_macro_f1": s_summ["overall_near_well_macro_f1"],
                "near_well_balanced_accuracy": s_summ["overall_near_well_balanced_acc"],
                "near_well_accuracy": s_summ["overall_near_well_accuracy"],
                "prior_macro_f1": s_summ["overall_prior_macro_f1"],
                "prior_accuracy": s_summ["overall_prior_accuracy"],
            })
    df_sens = pd.DataFrame(sens_records)
    sens_csv = output_dir / "datum_sensitivity_six_state.csv"
    df_sens.to_csv(sens_csv, index=False)
    print(f"  -> Exported datum sensitivity to {sens_csv}")

    # =========================================================================
    # 10. Publication-Quality Figure Generation
    # =========================================================================
    print("\n[10/10] Regenerating all 5 publication-quality figures with 6-state schema...")
    plot_vertical_successions(
        inspector=inspector,
        output_path=figures_dir / "transect_vertical_successions_six_state.png",
        litholog_ids=all_lids,
    )
    print("  -> Fig 1: Transect vertical successions generated.")

    plot_facies_proportions(
        inspector=inspector,
        output_path=figures_dir / "facies_proportions_and_ntg_six_state.png",
        litholog_ids=all_lids,
    )
    print("  -> Fig 2: Facies proportions and net-to-gross generated.")

    plot_upstream_vs_downstream(
        inspector=inspector,
        output_path=figures_dir / "upstream_vs_downstream_facies_ntg_six_state.png",
    )
    print("  -> Fig 3: Upstream vs Downstream comparison generated.")

    plot_sandstone_thickness_distributions(
        inspector=inspector,
        output_path=figures_dir / "sandstone_thickness_distributions_six_state.png",
    )
    print("  -> Fig 4: Sandstone thickness distributions generated.")

    # Generate heatmaps using All 12 logs 6x6 matrices
    P_reg_all, _, pi_reg_all = inspector.compute_vertical_transitions(all_12_ids, embedded=False, use_discretized=True, smoothing_alpha=0.1)
    P_emb_all, _, pi_emb_all = inspector.compute_vertical_transitions(all_12_ids, embedded=True, use_discretized=False, smoothing_alpha=0.1)
    plot_transition_matrix_heatmaps(
        P_reg=P_reg_all,
        P_emb=P_emb_all,
        pi_reg=pi_reg_all,
        pi_emb=pi_emb_all,
        output_path=figures_dir / "markov_transition_heatmaps_six_state.png",
    )
    print("  -> Fig 5: 6x6 Markov transition heatmaps generated.")

    print("\n" + "=" * 80)
    print("SMALT SPRINT F PIPELINE COMPLETED SUCCESSFULLY!")
    print(f"All deliverables saved to: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
