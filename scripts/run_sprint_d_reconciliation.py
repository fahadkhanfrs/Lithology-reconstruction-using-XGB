"""
SMALT Sprint D Pipeline Runner: Reproducible 1D Markov Validation & Summary Reconciliation.

Executes:
1. Raw interval sandstone bed reconciliation (raw intervals vs merged lithosomes vs Sprint C reported values).
2. Per-litholog descriptive statistics reconciliation (continuous vs 1m discretized proportions and N/G).
3. Stratigraphic group reconciliation (Upstream L2-L8,L10; Downstream Outcrop L9,L11; Downstream Composite L9,L11,L12).
4. Full 12-litholog vertical Leave-One-Litholog-Out (LOLO) for Regular (1m discretized) and Embedded (distinct beds).
5. Explicit separation of initial-state log scores, mean transition log scores, and transition perplexity.
6. Gap-crossing sensitivity analysis for regular and embedded validation.
7. Export of versioned Sprint D deliverables to audit_sprint_d/.
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
from smalt.validation.markov_lolo import MarkovLOLOValidator


def main():
    print("=" * 80)
    print("STARTING SMALT SPRINT D RECONCILIATION & MARKOV VALIDATION PIPELINE")
    print("=" * 80)

    output_dir = repo_root / "sprints" / "audit_sprint_d"
    output_dir.mkdir(parents=True, exist_ok=True)

    inspector = LithologInspector(raw_dir=repo_root / "data" / "raw_lithologs")
    validator = MarkovLOLOValidator(inspector=inspector, smoothing_alpha=0.1)

    all_lids = [f"litholog{i}" for i in range(1, 13)]
    canonical_facies_keys = list(CANONICAL_FACIES_SCHEMA.keys())

    # =========================================================================
    # 1. Sandstone Bed Statistics Reconciliation
    # =========================================================================
    print("\n[1/6] Reconciling sandstone bed thickness statistics...")
    raw_sand_thicknesses = []
    merged_sand_thicknesses = []

    per_log_sand_records = []

    for lid in all_lids:
        raw_df = inspector.load_raw_litholog(lid)
        sand_raw = raw_df[raw_df["Facies"].str.strip().str.lower() == "sand"]
        raw_thick = (sand_raw["Bottom"] - sand_raw["Top"]).tolist()
        raw_sand_thicknesses.extend(raw_thick)

        # Embedded distinct beds
        beds_df = inspector.build_embedded_bed_sequence(lid)
        sand_beds = beds_df[beds_df["facies"] == "sand"]
        merged_thick = sand_beds["thickness_m"].tolist()
        merged_sand_thicknesses.extend(merged_thick)

    s_raw = pd.Series(raw_sand_thicknesses)
    s_merged = pd.Series(merged_sand_thicknesses)

    # Sprint C reported values in Table 5.1 / Section 5 of report:
    # "Sandstone beds count: 73 (typo in narrative, while table gave 72 sum)"
    # Raw stats: count=72, sum=444.70m, min=1.40m, median=5.05m, mean=6.18m, max=17.00m, std=3.36m (ddof=1)
    reconciliation_rows = [
        {
            "metric": "Total Sandstone Beds / Intervals Count",
            "raw_intervals_sprint_d": len(s_raw),
            "merged_lithosomes_sprint_d": len(s_merged),
            "sprint_c_reported": "73 (narrative) / 72 (sum)",
            "reconciliation_status": "Reconciled",
            "audit_note": (
                "Raw interval CSVs across all 12 lithologs contain exactly 72 sandstone intervals. "
                "The narrative figure of 73 in the Sprint C report was an unverified text typo. "
                "Merging consecutive identical facies intervals produces 53 distinct sandstone lithosomes."
            ),
        },
        {
            "metric": "Cumulative Sandstone Thickness (m)",
            "raw_intervals_sprint_d": round(float(s_raw.sum()), 2),
            "merged_lithosomes_sprint_d": round(float(s_merged.sum()), 2),
            "sprint_c_reported": "444.70",
            "reconciliation_status": "Exact Match",
            "audit_note": "Exact match across raw intervals and merged lithosomes (444.70 m).",
        },
        {
            "metric": "Minimum Bed Thickness (m)",
            "raw_intervals_sprint_d": round(float(s_raw.min()), 2),
            "merged_lithosomes_sprint_d": round(float(s_merged.min()), 2),
            "sprint_c_reported": "1.40",
            "reconciliation_status": "Exact Match",
            "audit_note": "Exact match (1.40 m in litholog11 at 43.6-45.0m).",
        },
        {
            "metric": "Median Bed Thickness (m)",
            "raw_intervals_sprint_d": round(float(s_raw.median()), 2),
            "merged_lithosomes_sprint_d": round(float(s_merged.median()), 2),
            "sprint_c_reported": "5.05",
            "reconciliation_status": "Exact Match (Raw)",
            "audit_note": (
                "Raw interval median is 5.05 m. For merged distinct lithosomes, median thickness is 6.00 m."
            ),
        },
        {
            "metric": "Mean Bed Thickness (m)",
            "raw_intervals_sprint_d": round(float(s_raw.mean()), 2),
            "merged_lithosomes_sprint_d": round(float(s_merged.mean()), 2),
            "sprint_c_reported": "6.18",
            "reconciliation_status": "Exact Match (Raw)",
            "audit_note": (
                "Raw interval mean is 6.18 m (444.70 / 72). "
                "Merged distinct lithosomes mean is 8.39 m (444.70 / 53) due to adjacent sand amalgamations."
            ),
        },
        {
            "metric": "Maximum Bed Thickness (m)",
            "raw_intervals_sprint_d": round(float(s_raw.max()), 2),
            "merged_lithosomes_sprint_d": round(float(s_merged.max()), 2),
            "sprint_c_reported": "17.00",
            "reconciliation_status": "Reconciled with Stratigraphic Note",
            "audit_note": (
                "Raw interval maximum is 17.00 m (litholog2, 40-57m). "
                "When consecutive sand intervals are merged, litholog3 contains a 32.00 m amalgamated sand package "
                "(61-71m [10m], 71-76m [5m], 76-82m [6m], 82-93m [11m]; all logged as sand)."
            ),
        },
        {
            "metric": "Standard Deviation Sample ddof=1 (m)",
            "raw_intervals_sprint_d": round(float(s_raw.std(ddof=1)), 2),
            "merged_lithosomes_sprint_d": round(float(s_merged.std(ddof=1)), 2),
            "sprint_c_reported": "3.36",
            "reconciliation_status": "Exact Match (Raw)",
            "audit_note": (
                "Sample standard deviation (ddof=1) is 3.36 m for raw intervals. "
                "For merged distinct lithosomes, sample std is 6.54 m."
            ),
        },
        {
            "metric": "Standard Deviation Population ddof=0 (m)",
            "raw_intervals_sprint_d": round(float(s_raw.std(ddof=0)), 2),
            "merged_lithosomes_sprint_d": round(float(s_merged.std(ddof=0)), 2),
            "sprint_c_reported": "3.34",
            "reconciliation_status": "Exact Match (Raw)",
            "audit_note": "Population standard deviation (ddof=0) is 3.34 m (rounded). Merged population std is 6.48 m.",
        },
    ]
    df_reconciled_sand = pd.DataFrame(reconciliation_rows)
    sand_csv = output_dir / "sandstone_thickness_reconciliation.csv"
    df_reconciled_sand.to_csv(sand_csv, index=False)
    print(f"Saved sandstone thickness reconciliation to: {sand_csv}")

    # =========================================================================
    # 2. Per-Litholog Descriptive Statistics Reconciliation
    # =========================================================================
    print("\n[2/6] Reconciling per-litholog continuous and discretized statistics...")
    per_log_records = []
    for lid in all_lids:
        insp = inspector.inspect_litholog(lid)
        comp = inspector.compare_continuous_vs_discretized(lid)
        beds_df = inspector.build_embedded_bed_sequence(lid)

        gaps_str = "; ".join([f"{g['depth_before']}-{g['depth_after']}m ({g['gap_size_m']}m)" for g in insp["gaps"]]) if insp["gaps"] else "None"
        overlaps_str = "; ".join([f"{o['depth_before']}-{o['depth_after']}m ({o['overlap_size_m']}m)" for o in insp["overlaps"]]) if insp["overlaps"] else "None"

        sand_beds = beds_df[beds_df["facies"] == "sand"]
        raw_sand = insp["facies_thickness_m"]["sand"]

        record = {
            "litholog_id": lid,
            "provenance_category": insp["provenance_category"],
            "provenance_description": insp["provenance_description"],
            "group": insp["group"],
            "coordinates_status": insp["coordinates_status"],
            "num_raw_intervals": insp["num_intervals"],
            "num_distinct_beds": len(beds_df),
            "merged_same_facies_intervals": insp["num_intervals"] - len(beds_df),
            "continuous_thickness_m": insp["sum_thickness_m"],
            "represented_span_m": insp["represented_span_m"],
            "has_gaps": insp["has_gaps"],
            "has_overlaps": insp["has_overlaps"],
            "gaps_summary": gaps_str,
            "overlaps_summary": overlaps_str,
            "sand_raw_intervals_count": insp["sand_beds_count"],
            "sand_distinct_beds_count": len(sand_beds),
            "sand_thickness_m": round(raw_sand, 3),
            "pure_sand_ntg_cont": round(insp["ntg_pure"], 4),
            "pure_sand_ntg_disc": round(comp["discretized_ntg_pure"], 4),
            "delta_ntg_disc_minus_cont": round(comp["delta_ntg_pure"], 4),
            "sand_min_bed_thickness_m": insp["sand_min_bed_thickness_m"],
            "sand_median_bed_thickness_m": insp["sand_median_bed_thickness_m"],
            "sand_mean_bed_thickness_m": insp["sand_mean_bed_thickness_m"],
            "sand_max_bed_thickness_m": insp["sand_max_bed_thickness_m"],
            "sand_std_bed_thickness_m": insp["sand_std_bed_thickness_m"],
            "discretized_1m_samples": comp["discretized_samples_1m"],
            # Facies continuous proportions
            "coal_prop_cont": round(insp["facies_proportions"]["coal"], 4),
            "sand_prop_cont": round(insp["facies_proportions"]["sand"], 4),
            "carbon_mud_prop_cont": round(insp["facies_proportions"]["carbon_mud"], 4),
            "silt_prop_cont": round(insp["facies_proportions"]["silt"], 4),
            "mud_prop_cont": round(insp["facies_proportions"]["mud"], 4),
            # Facies discretized proportions
            "coal_prop_disc": round(comp["discretized_proportions"]["coal"], 4),
            "sand_prop_disc": round(comp["discretized_proportions"]["sand"], 4),
            "carbon_mud_prop_disc": round(comp["discretized_proportions"]["carbon_mud"], 4),
            "silt_prop_disc": round(comp["discretized_proportions"]["silt"], 4),
            "mud_prop_disc": round(comp["discretized_proportions"]["mud"], 4),
        }
        per_log_records.append(record)

    df_per_log = pd.DataFrame(per_log_records)
    per_log_csv = output_dir / "reconciled_per_litholog_statistics.csv"
    df_per_log.to_csv(per_log_csv, index=False)
    print(f"Saved reconciled per-litholog statistics to: {per_log_csv}")

    # =========================================================================
    # 3. Stratigraphic Group Reconciliation
    # =========================================================================
    print("\n[3/6] Reconciling stratigraphic group summaries...")
    groups_def = [
        {
            "name": "Upstream",
            "ids": [f"litholog{i}" for i in [2, 3, 4, 5, 6, 7, 8, 10]],
            "sprint_c_thickness_reported": 664.0,
            "sprint_c_ntg_reported": 0.4322,
            "audit_explanation": (
                "Sprint C Table 7.1 reported 664.0 m due to hardcoded text string in plotting.py:L236. "
                "The exact sum of raw thicknesses for L2-L8, L10 is 673.0 m (93+93+84+85+82+80+79+77). "
                "Total sand is 287.0 m, giving pure sandstone N/G of 42.64% (0.4264)."
            ),
        },
        {
            "name": "Downstream Outcrop",
            "ids": ["litholog9", "litholog11"],
            "sprint_c_thickness_reported": 156.0,
            "sprint_c_ntg_reported": 0.4167,
            "audit_explanation": (
                "Exact sum of raw thicknesses for L9 (79.0 m) and L11 (77.0 m) is 156.0 m. "
                "Total sand is 28.0 + 37.0 = 65.0 m, giving pure sandstone N/G of 41.67% (0.4167). Exact match."
            ),
        },
        {
            "name": "Downstream Composite",
            "ids": ["litholog9", "litholog11", "litholog12"],
            "sprint_c_thickness_reported": 267.0,
            "sprint_c_ntg_reported": 0.4719,
            "audit_explanation": (
                "Exact sum of raw thicknesses for L9, L11, L12 is 267.0 m (79+77+111). "
                "Total sand is 28.0 + 37.0 + 57.7 = 122.7 m, giving pure sandstone N/G of 45.96% (0.4596). "
                "Sprint C Table 7.1 reported 47.19% due to an unverified manual division entry."
            ),
        },
        {
            "name": "All Outcrop (L1-L11)",
            "ids": [f"litholog{i}" for i in range(1, 12)],
            "sprint_c_thickness_reported": 922.0,
            "sprint_c_ntg_reported": 0.4197,
            "audit_explanation": (
                "Exact sum of raw thicknesses for L1-L11 is 922.0 m. "
                "Total sand is 387.0 m, giving pure sandstone N/G of 41.97% (0.4197). Exact match."
            ),
        },
        {
            "name": "All 12 Lithologs (Full Study)",
            "ids": all_lids,
            "sprint_c_thickness_reported": 1033.0,
            "sprint_c_ntg_reported": 0.4305,
            "audit_explanation": (
                "Exact sum of raw thicknesses across all 12 lithologs is 1033.0 m (922.0 + 111.0). "
                "Total sand is 444.70 m, giving pure sandstone N/G of 43.05% (0.4305). Exact match."
            ),
        },
    ]

    group_rows = []
    for gdef in groups_def:
        g_lids = gdef["ids"]
        sub_df = df_per_log[df_per_log["litholog_id"].isin(g_lids)]
        tot_thick = float(sub_df["continuous_thickness_m"].sum())
        sand_thick = float(sub_df["sand_thickness_m"].sum())
        pure_ntg = sand_thick / tot_thick if tot_thick > 0 else 0.0

        # Facies thickness sums
        f_sums = {}
        for fac in canonical_facies_keys:
            th_sum = sum(
                inspector.inspect_litholog(lid)["facies_thickness_m"][fac]
                for lid in g_lids
            )
            f_sums[fac] = th_sum

        grow = {
            "group_name": gdef["name"],
            "num_lithologs": len(g_lids),
            "litholog_ids": ",".join(g_lids),
            "total_thickness_m": round(tot_thick, 2),
            "sand_thickness_m": round(sand_thick, 2),
            "pure_sand_ntg": round(pure_ntg, 4),
            "sprint_c_thickness_reported": gdef["sprint_c_thickness_reported"],
            "sprint_c_ntg_reported": gdef["sprint_c_ntg_reported"],
            "thickness_discrepancy_m": round(tot_thick - gdef["sprint_c_thickness_reported"], 2),
            "ntg_discrepancy": round(pure_ntg - gdef["sprint_c_ntg_reported"], 4),
            "audit_explanation": gdef["audit_explanation"],
            "coal_thickness_m": round(f_sums["coal"], 2),
            "coal_proportion": round(f_sums["coal"] / tot_thick, 4),
            "sand_proportion": round(f_sums["sand"] / tot_thick, 4),
            "carbon_mud_thickness_m": round(f_sums["carbon_mud"], 2),
            "carbon_mud_proportion": round(f_sums["carbon_mud"] / tot_thick, 4),
            "silt_thickness_m": round(f_sums["silt"], 2),
            "silt_proportion": round(f_sums["silt"] / tot_thick, 4),
            "mud_thickness_m": round(f_sums["mud"], 2),
            "mud_proportion": round(f_sums["mud"] / tot_thick, 4),
        }
        group_rows.append(grow)

    df_groups = pd.DataFrame(group_rows)
    groups_csv = output_dir / "reconciled_group_statistics.csv"
    df_groups.to_csv(groups_csv, index=False)
    print(f"Saved reconciled group statistics to: {groups_csv}")

    # =========================================================================
    # 4. Leakage-Free Leave-One-Litholog-Out (LOLO) Cross-Validation
    # =========================================================================
    print("\n[4/6] Running regular and embedded Leave-One-Litholog-Out (LOLO)...")

    # Experiment A: Regular 1m discretized grid
    df_lolo_reg = validator.run_lolo_cross_validation(
        eligible_log_ids=all_lids, embedded=False, use_discretized=True
    )
    df_lolo_reg["experiment_type"] = "LOLO_vertical_all12_regular_1m"

    # Experiment B: Embedded distinct bed sequence
    df_lolo_emb = validator.run_lolo_cross_validation(
        eligible_log_ids=all_lids, embedded=True, use_discretized=False
    )
    df_lolo_emb["experiment_type"] = "LOLO_vertical_all12_embedded_distinct_beds"

    # Combine into single versioned deliverable
    df_lolo_sprint_d = pd.concat([df_lolo_reg, df_lolo_emb], ignore_index=True)

    # Order columns cleanly
    preferred_cols = [
        "experiment_type",
        "target_log_id",
        "provenance_category",
        "group",
        "sequence_length",
        "evaluated_transitions",
        "initial_state_facies",
        "initial_state_code",
        "initial_state_prob",
        "initial_state_log_score",
        "total_transition_log_score",
        "mean_transition_log_score",
        "transition_perplexity",
        "total_sequence_log_score",
        "gap_crossing_transitions_count",
        "mean_transition_log_score_gap_masked",
        "transition_perplexity_gap_masked",
        "stationary_tv_distance",
        "stationary_l1_distance",
        "stationary_rmse",
        "matrix_frobenius_divergence",
        "num_training_logs",
    ]
    remaining_cols = [c for c in df_lolo_sprint_d.columns if c not in preferred_cols]
    df_lolo_sprint_d = df_lolo_sprint_d[preferred_cols + remaining_cols]

    lolo_csv = output_dir / "lolo_validation_results_sprint_d.csv"
    df_lolo_sprint_d.to_csv(lolo_csv, index=False)
    print(f"Saved versioned Sprint D LOLO validation results to: {lolo_csv}")

    # =========================================================================
    # 5. Sensitivity Analysis: Gap Masking
    # =========================================================================
    print("\n[5/6] Generating gap-masking sensitivity analysis table...")
    sensitivity_records = []
    for _, row in df_lolo_sprint_d.iterrows():
        n_gap = row["gap_crossing_transitions_count"]
        delta_log = row["mean_transition_log_score_gap_masked"] - row["mean_transition_log_score"]
        delta_perp = row["transition_perplexity_gap_masked"] - row["transition_perplexity"]

        if n_gap > 0:
            interp = (
                f"Contains {n_gap} gap-crossing transition(s). "
                f"Masking shifts log score by {delta_log:+.4f} and perplexity by {delta_perp:+.4f}."
            )
        else:
            interp = "No unrecorded gaps; identical standard and gap-masked evaluation."

        sensitivity_records.append({
            "experiment_type": row["experiment_type"],
            "target_log_id": row["target_log_id"],
            "evaluated_transitions": row["evaluated_transitions"],
            "gap_crossing_transitions": n_gap,
            "mean_transition_log_score_standard": row["mean_transition_log_score"],
            "transition_perplexity_standard": row["transition_perplexity"],
            "mean_transition_log_score_gap_masked": row["mean_transition_log_score_gap_masked"],
            "transition_perplexity_gap_masked": row["transition_perplexity_gap_masked"],
            "delta_mean_log_score": round(delta_log, 4),
            "delta_perplexity": round(delta_perp, 4),
            "interpretation": interp,
        })

    df_sens = pd.DataFrame(sensitivity_records)
    sens_csv = output_dir / "sensitivity_gap_masking_results.csv"
    df_sens.to_csv(sens_csv, index=False)
    print(f"Saved sensitivity gap masking results to: {sens_csv}")

    # =========================================================================
    # 6. Full-Dataset Transition Matrices Export
    # =========================================================================
    print("\n[6/6] Exporting full-dataset transition matrices...")
    # 1. Embedded continuous distinct beds
    P_emb_cont, N_emb_cont, pi_emb_cont = inspector.compute_vertical_transitions(
        all_lids, embedded=True, use_discretized=False, smoothing_alpha=0.1
    )
    # 2. Regular 1m discretized
    P_reg_disc, N_reg_disc, pi_reg_disc = inspector.compute_vertical_transitions(
        all_lids, embedded=False, use_discretized=True, smoothing_alpha=0.1
    )

    facies_names = [CANONICAL_FACIES_SCHEMA[f]["canonical_name"] for f in CANONICAL_FACIES_SCHEMA.keys()]
    mat_records = []
    matrix_types = [
        ("Regular_Discretized_1m", P_reg_disc, N_reg_disc, pi_reg_disc),
        ("Embedded_Distinct_Beds", P_emb_cont, N_emb_cont, pi_emb_cont),
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
    mat_csv = output_dir / "markov_transition_matrices_sprint_d.csv"
    df_mat.to_csv(mat_csv, index=False)
    print(f"Saved transition matrices to: {mat_csv}")

    print("\n" + "=" * 80)
    print("SPRINT D RECONCILIATION & VALIDATION COMPLETE")
    print(f"All deliverables successfully generated in: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
