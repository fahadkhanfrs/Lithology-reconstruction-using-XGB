"""
SMALT Sprint C Pipeline Runner: Descriptive Lithology & Leakage-Free 1D Markov Validation.

Executes:
1. Complete raw interval inspection & validation for Lithologs 1 - 12.
2. Continuous vs 1m discretized facies proportions and Net-to-Gross calculations.
3. Upstream (proximal) vs Downstream (distal) comparative analysis.
4. Regular & Embedded vertical Markov transition matrix estimation.
5. Leakage-free Leave-One-Litholog-Out (LOLO) cross-validation.
6. Directional (Upstream <-> Downstream) validation experiments.
7. Publication-quality figure generation.
8. Deliverable CSV export to audit_sprint_c/.
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
from smalt.descriptive.plotting import (
    plot_vertical_successions,
    plot_facies_proportions,
    plot_upstream_vs_downstream,
    plot_sandstone_thickness_distributions,
    plot_transition_matrix_heatmaps,
)
from smalt.validation.markov_lolo import MarkovLOLOValidator


def main():
    print("=" * 80)
    print("STARTING SMALT SPRINT C DESCRIPTIVE & MARKOV VALIDATION PIPELINE")
    print("=" * 80)

    output_dir = repo_root / "sprints" / "audit_sprint_c"
    figures_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    inspector = LithologInspector(raw_dir=repo_root / "data" / "raw_lithologs")
    validator = MarkovLOLOValidator(inspector=inspector, smoothing_alpha=0.1)

    all_lids = [f"litholog{i}" for i in range(1, 13)]

    # -------------------------------------------------------------------------
    # 1. Per-Litholog Descriptive Statistics & Inspection
    # -------------------------------------------------------------------------
    print("\n[1/6] Running raw interval inspection & continuous vs discretized comparisons...")
    stats_records = []
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
            "sand_thickness_m": insp["facies_thickness_m"]["sand"],
            "sand_prop_cont": insp["facies_proportions"]["sand"],
            "mud_thickness_m": insp["facies_thickness_m"]["mud"],
            "mud_prop_cont": insp["facies_proportions"]["mud"],
            "carbon_mud_thickness_m": insp["facies_thickness_m"]["carbon_mud"],
            "carbon_mud_prop_cont": insp["facies_proportions"]["carbon_mud"],
            "silt_thickness_m": insp["facies_thickness_m"]["silt"],
            "silt_prop_cont": insp["facies_proportions"]["silt"],
            "coal_thickness_m": insp["facies_thickness_m"]["coal"],
            "coal_prop_cont": insp["facies_proportions"]["coal"],
            "ntg_pure_cont": insp["ntg_pure"],
            "ntg_coarse_cont": insp["ntg_coarse"],
            "coal_beds_count": insp["coal_beds_count"],
            "coal_mean_bed_thickness_m": insp["coal_mean_bed_thickness_m"],
            "coal_max_bed_thickness_m": insp["coal_max_bed_thickness_m"],
            "sand_beds_count": insp["sand_beds_count"],
            "sand_min_bed_thickness_m": insp["sand_min_bed_thickness_m"],
            "sand_median_bed_thickness_m": insp["sand_median_bed_thickness_m"],
            "sand_mean_bed_thickness_m": insp["sand_mean_bed_thickness_m"],
            "sand_max_bed_thickness_m": insp["sand_max_bed_thickness_m"],
            "sand_std_bed_thickness_m": insp["sand_std_bed_thickness_m"],
            "discretized_samples_1m": comp["discretized_samples_1m"],
            "sand_prop_disc": comp["discretized_proportions"]["sand"],
            "mud_prop_disc": comp["discretized_proportions"]["mud"],
            "carbon_mud_prop_disc": comp["discretized_proportions"]["carbon_mud"],
            "silt_prop_disc": comp["discretized_proportions"]["silt"],
            "coal_prop_disc": comp["discretized_proportions"]["coal"],
            "ntg_pure_disc": comp["discretized_ntg_pure"],
            "delta_ntg_pure_disc_minus_cont": comp["delta_ntg_pure"],
        }
        stats_records.append(record)

    df_stats = pd.DataFrame(stats_records)
    stats_csv = output_dir / "per_litholog_descriptive_statistics.csv"
    df_stats.to_csv(stats_csv, index=False)
    print(f"Saved descriptive statistics to: {stats_csv}")

    # -------------------------------------------------------------------------
    # 2. Facies Mapping Table
    # -------------------------------------------------------------------------
    print("\n[2/6] Generating explicit facies mapping table...")
    mapping_rows = [
        {
            "raw_facies_label": "coal",
            "canonical_code": 0,
            "canonical_facies_name": "Coal",
            "lithological_class": "Biogenic organic deposit / mire facies",
            "color_hex": "#1C2833",
            "base_gr_api_phase0": 20.0,
            "harmonization_status": "Unambiguous",
            "audit_justification": (
                "Direct 1:1 correspondence. Coal intervals form laterally extensive markers across the Blackhawk "
                "Formation. Readily distinguished in both outcrop and core."
            ),
        },
        {
            "raw_facies_label": "sand",
            "canonical_code": 1,
            "canonical_facies_name": "Channel Sandstone (undivided)",
            "lithological_class": "Coarse siliciclastic / channel-belt sandstone",
            "color_hex": "#F4D03F",
            "base_gr_api_phase0": 35.0,
            "harmonization_status": "Documented Caveat",
            "audit_justification": (
                "Mapped as 'Channel Sandstone' in Phase 0 schema. Raw logs do not distinguish channel-axis, "
                "channel-margin, crevasse splay, or sheet sandstones. All undivided sandstone is grouped here."
            ),
        },
        {
            "raw_facies_label": "carbon_mud",
            "canonical_code": 2,
            "canonical_facies_name": "Carbonaceous Mudstone",
            "lithological_class": "Organic-rich muddy sediment / poorly drained swamp",
            "color_hex": "#6C3483",
            "base_gr_api_phase0": 130.0,
            "harmonization_status": "FLAGGED CONTRADICTION",
            "audit_justification": (
                "CRITICAL CONTRADICTION IN PHASE 0: CANONICAL_FACIES_NAMES[2] labeled this 'Fine Sandstone / Splay', "
                "while CANONICAL_TECHNICAL_NAMES[2] used 'carbon_mud' with base GR 130 API. Geologically, carbonaceous "
                "mudstone is organic-rich mud, NOT fine sandstone. Sprint C restores the proper name Carbonaceous Mudstone."
            ),
        },
        {
            "raw_facies_label": "silt",
            "canonical_code": 3,
            "canonical_facies_name": "Siltstone",
            "lithological_class": "Fine-grained siliciclastic / floodplain-levee transition",
            "color_hex": "#73C6B6",
            "base_gr_api_phase0": 70.0,
            "harmonization_status": "Documented Caveat",
            "audit_justification": (
                "Unambiguous fine siliciclastic lithology. Omitted from legacy src/litholog_model.py encode_facies "
                "dictionary, causing silent dropna() row deletion. Fully retained in Sprint C."
            ),
        },
        {
            "raw_facies_label": "mud",
            "canonical_code": 4,
            "canonical_facies_name": "Overbank Mudstone",
            "lithological_class": "Fine siliciclastic / floodplain-overbank mudstone",
            "color_hex": "#95A5A6",
            "base_gr_api_phase0": 105.0,
            "harmonization_status": "Documented Caveat",
            "audit_justification": (
                "Dominant floodplain overbank mudstone. May encompass minor lacustrine, prodelta, or abandoned "
                "channel-fill clay plugs."
            ),
        },
    ]
    df_map = pd.DataFrame(mapping_rows)
    map_csv = output_dir / "facies_mapping_table.csv"
    df_map.to_csv(map_csv, index=False)
    print(f"Saved facies mapping table to: {map_csv}")

    # -------------------------------------------------------------------------
    # 3. Regular & Embedded Transition Matrices Across Full Dataset
    # -------------------------------------------------------------------------
    print("\n[3/6] Estimating full-dataset Markov transition matrices...")
    # 1. Embedded on continuous intervals (true bed boundary changes)
    P_emb_cont, N_emb_cont, pi_emb_cont = inspector.compute_vertical_transitions(
        all_lids, embedded=True, use_discretized=False, smoothing_alpha=0.1
    )
    # 2. Regular on 1m discretized grid
    P_reg_disc, N_reg_disc, pi_reg_disc = inspector.compute_vertical_transitions(
        all_lids, embedded=False, use_discretized=True, smoothing_alpha=0.1
    )
    # 3. Embedded on 1m discretized grid
    P_emb_disc, N_emb_disc, pi_emb_disc = inspector.compute_vertical_transitions(
        all_lids, embedded=True, use_discretized=True, smoothing_alpha=0.1
    )

    facies_names = [CANONICAL_FACIES_SCHEMA[f]["canonical_name"] for f in CANONICAL_FACIES_SCHEMA.keys()]

    # Format matrices into tidy CSV
    mat_records = []
    matrix_types = [
        ("Regular_Discretized_1m", P_reg_disc, N_reg_disc, pi_reg_disc),
        ("Embedded_Continuous_Intervals", P_emb_cont, N_emb_cont, pi_emb_cont),
        ("Embedded_Discretized_1m", P_emb_disc, N_emb_disc, pi_emb_disc),
    ]

    for mtype, P, N, pi in matrix_types:
        for i, from_name in enumerate(facies_names):
            row = {
                "matrix_type": mtype,
                "from_facies_code": i,
                "from_facies_name": from_name,
                "stationary_probability": round(float(pi[i]), 5),
                "total_transitions_from_state": float(N[i, :].sum()),
            }
            for j, to_name in enumerate(facies_names):
                row[f"prob_to_{to_name}"] = round(float(P[i, j]), 5)
                row[f"count_to_{to_name}"] = float(N[i, j])
            mat_records.append(row)

    df_mat = pd.DataFrame(mat_records)
    mat_csv = output_dir / "markov_transition_matrices.csv"
    df_mat.to_csv(mat_csv, index=False)
    print(f"Saved transition matrices to: {mat_csv}")

    # -------------------------------------------------------------------------
    # 4. Leakage-Free Leave-One-Litholog-Out (LOLO) Validation
    # -------------------------------------------------------------------------
    print("\n[4/6] Running leakage-free Leave-One-Litholog-Out (LOLO) cross-validation...")
    # Experiment A: Vertical LOLO on all 12 lithologs (Regular, 1m discretized)
    df_lolo_reg = validator.run_lolo_cross_validation(
        eligible_log_ids=all_lids, embedded=False, use_discretized=True
    )
    df_lolo_reg["experiment_type"] = "LOLO_vertical_all12_regular_1m"

    # Experiment B: Vertical LOLO on all 12 lithologs (Embedded, continuous beds)
    df_lolo_emb = validator.run_lolo_cross_validation(
        eligible_log_ids=all_lids, embedded=True, use_discretized=False
    )
    df_lolo_emb["experiment_type"] = "LOLO_vertical_all12_embedded_continuous"

    # Combine LOLO records
    df_lolo_all = pd.concat([df_lolo_reg, df_lolo_emb], ignore_index=True)
    lolo_csv = output_dir / "lolo_validation_results.csv"
    df_lolo_all.to_csv(lolo_csv, index=False)
    print(f"Saved LOLO validation results to: {lolo_csv}")

    # -------------------------------------------------------------------------
    # 5. Directional Validation Experiments (Upstream <-> Downstream)
    # -------------------------------------------------------------------------
    print("\n[5/6] Running directional (Upstream <-> Downstream) experiments...")
    upstream_ids = [f"litholog{i}" for i in [2, 3, 4, 5, 6, 7, 8, 10]]
    downstream_ids = ["litholog9", "litholog11", "litholog12"]

    dir_up_down = validator.run_directional_experiment(
        training_group_ids=upstream_ids,
        target_group_ids=downstream_ids,
        experiment_name="Upstream_to_Downstream_Regular_1m",
        embedded=False,
        use_discretized=True,
    )

    dir_down_up = validator.run_directional_experiment(
        training_group_ids=downstream_ids,
        target_group_ids=upstream_ids,
        experiment_name="Downstream_to_Upstream_Regular_1m",
        embedded=False,
        use_discretized=True,
    )

    directional_records = []
    # Up -> Down individual results
    for t_res in dir_up_down["per_target_results"]:
        t_res["experiment_name"] = dir_up_down["experiment_name"]
        t_res["training_group_size"] = dir_up_down["n_training_logs"]
        t_res["target_group_size"] = dir_up_down["n_target_logs"]
        t_res["composite_group_frobenius"] = dir_up_down["composite_frobenius_divergence"]
        t_res["composite_group_tv_distance"] = dir_up_down["composite_stationary_tv_distance"]
        directional_records.append(t_res)

    # Down -> Up individual results
    for t_res in dir_down_up["per_target_results"]:
        t_res["experiment_name"] = dir_down_up["experiment_name"]
        t_res["training_group_size"] = dir_down_up["n_training_logs"]
        t_res["target_group_size"] = dir_down_up["n_target_logs"]
        t_res["composite_group_frobenius"] = dir_down_up["composite_frobenius_divergence"]
        t_res["composite_group_tv_distance"] = dir_down_up["composite_stationary_tv_distance"]
        directional_records.append(t_res)

    df_directional = pd.DataFrame(directional_records)
    dir_csv = output_dir / "directional_validation_results.csv"
    df_directional.to_csv(dir_csv, index=False)
    print(f"Saved directional validation results to: {dir_csv}")

    # -------------------------------------------------------------------------
    # 6. Generate Publication-Quality Figures
    # -------------------------------------------------------------------------
    print("\n[6/6] Generating figures for undergraduate research presentation...")
    plot_vertical_successions(inspector, figures_dir / "vertical_facies_successions.png")
    plot_facies_proportions(inspector, figures_dir / "facies_proportions_comparison.png")
    plot_upstream_vs_downstream(inspector, figures_dir / "upstream_vs_downstream_facies.png")
    plot_sandstone_thickness_distributions(inspector, figures_dir / "sandstone_thickness_distributions.png")
    plot_transition_matrix_heatmaps(P_reg_disc, P_emb_cont, pi_reg_disc, pi_emb_cont, figures_dir / "markov_transition_matrices_heatmaps.png")

    print(f"Figures saved in: {figures_dir}")
    print("=" * 80)
    print("SPRINT C PIPELINE EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
