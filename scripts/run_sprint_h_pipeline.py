"""
SMALT Sprint H Master Pipeline: Common-Zero Datum Alignment & Spatial Markov Foundation.

Executes:
1. Standardizes all 12 lithologs under the common-zero datum reference (z = 0 at depth = 0).
2. Validates datum alignment invariance (interval counts, thickness, facies counts).
3. Produces publication-quality diagnostic plot of all 12 lithologs aligned at z = 0.
4. Constructs empirical horizontal transition data across matched elevation slices.
5. Computes lag-binned empirical horizontal transition probability matrices.
6. Evaluates continuous spatial Markov transition model against baselines via LOLO CV.
7. Produces transition decay diagnostic figures.
"""

import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from smalt.descriptive.analyzer import LithologInspector, CANONICAL_FACIES_SCHEMA
from smalt.spatial.datum import (
    align_to_common_datum,
    build_common_datum_metadata_table,
    verify_datum_invariance,
    DEFAULT_CONVENTION,
)
from smalt.spatial.baseline import ProvisionalSpatialValidator
from smalt.geostat.spatial_markov import (
    SpatialMarkovTransitionAnalyzer,
    SpatialMarkovPredictor,
    DEFAULT_LATERAL_FACIES_LENGTHS_M,
)

OUTPUT_DIR = Path("sprints/audit_sprint_h")
FIGURES_DIR = OUTPUT_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Canonical colors and names for the 6 facies
FACIES_PALETTE = {
    0: ("#DAA520", "Channel Sandstone"),
    1: ("#FF8C00", "Planar Sandstone"),
    2: ("#4682B4", "Rippled Heterolithics"),
    3: ("#2F4F4F", "Carbonaceous Mudstone"),
    4: ("#1C1C1C", "Coal"),
    5: ("#556B2F", "Overbank Mudstone"),
}


def step1_datum_alignment_and_validation(inspector: LithologInspector):
    print("--- Step 1: Common-Zero Datum Alignment & Invariance Audit ---")
    df_meta = build_common_datum_metadata_table(inspector=inspector)
    meta_path = OUTPUT_DIR / "common_datum_litholog_metadata.csv"
    df_meta.to_csv(meta_path, index=False)
    print(f"Saved metadata table: {meta_path} ({len(df_meta)} rows)")

    df_inv = verify_datum_invariance(inspector=inspector)
    inv_path = OUTPUT_DIR / "common_datum_facies_consistency.csv"
    df_inv.to_csv(inv_path, index=False)
    print(f"Saved consistency table: {inv_path} ({len(df_inv)} rows)")
    return df_meta, df_inv


def step2_generate_diagnostic_alignment_plot(inspector: LithologInspector, df_meta: pd.DataFrame):
    print("--- Step 2: Diagnostic Common-Zero Transect Plot ---")
    all_lids = [f"litholog{i}" for i in range(1, 13)]

    fig, ax = plt.subplots(figsize=(16, 9), dpi=300)

    # Plot common datum horizontal line
    ax.axhline(0.0, color="#B22222", linestyle="--", linewidth=2.0, zorder=5, label="Common-Zero Reference Level (z = 0 m)")

    col_width = 0.55
    x_positions = np.arange(len(all_lids))

    for idx, lid in enumerate(all_lids):
        raw_df = inspector.load_raw_litholog(lid)
        aligned_df = align_to_common_datum(raw_df, litholog_id=lid, convention=DEFAULT_CONVENTION)
        x_pos = x_positions[idx]

        for _, row in aligned_df.iterrows():
            z_top = row["z_common_top_m"]
            z_bot = row["z_common_bottom_m"]
            thick = abs(z_bot - z_top)
            code = int(row["facies_code"])
            color, _ = FACIES_PALETTE[code]

            rect = plt.Rectangle(
                (x_pos - col_width / 2.0, z_bot),
                col_width,
                thick,
                facecolor=color,
                edgecolor="#333333",
                linewidth=0.5,
                zorder=3,
            )
            ax.add_patch(rect)

        # Base label
        meta_row = df_meta[df_meta["litholog"] == lid].iloc[0]
        base_z = meta_row["common_z_min"]
        ax.text(
            x_pos,
            base_z - 2.5,
            f"{meta_row['total_thickness']:.1f} m",
            ha="center",
            va="top",
            fontsize=8.5,
            fontweight="bold",
            color="#222222",
        )

    # Styling and formatting
    ax.set_xticks(x_positions)
    labels = []
    for lid in all_lids:
        if lid == "litholog1":
            labels.append("Litholog 1\n(Vertical Only;\nNo Coords)")
        else:
            labels.append(f"{lid.replace('litholog', 'Litholog ')}\n(Spatial Eligible)")
    ax.set_xticklabels(labels, fontsize=9.5, fontweight="medium")

    ax.set_ylabel("Common-Datum Elevation z_common (m)\n[z = -depth; Reference Level z = 0.0 m]", fontsize=11, fontweight="bold")
    ax.set_title(
        "SMALT Common-Zero Vertical Datum Alignment Across All 12 Lithologs\n"
        "Provisional Reference Standard Supplied by Prof. Hiranya Sahoo: All Log Zeros at Common Level",
        fontsize=13,
        fontweight="bold",
        pad=15,
    )
    ax.set_ylim(-125.0, 10.0)
    ax.grid(axis="y", linestyle=":", alpha=0.5, color="#888888")

    # Legend patches for 6 facies + reference datum line
    legend_elements = [
        plt.Line2D([0], [0], color="#B22222", linestyle="--", linewidth=2.0, label="Common-Zero Level (z = 0 m)")
    ]
    for code, (color, name) in FACIES_PALETTE.items():
        legend_elements.append(mpatches.Patch(facecolor=color, edgecolor="#333333", label=f"Facies {code+1}: {name}"))

    ax.legend(
        handles=legend_elements,
        loc="upper right",
        bbox_to_anchor=(0.99, 0.98),
        framealpha=0.95,
        fontsize=9.0,
        title="Facies Schema & Datum",
        title_fontsize=10.0,
    )

    plt.tight_layout()
    fig_path = FIGURES_DIR / "common_zero_transect_alignment.png"
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"Saved diagnostic alignment plot: {fig_path}")


def step3_spatial_transition_data_construction(analyzer: SpatialMarkovTransitionAnalyzer):
    print("--- Step 3: Spatial Transition Pair Extraction & Lag Matrices ---")
    df_pairs = analyzer.extract_horizontal_facies_pairs()
    print(f"Extracted {len(df_pairs)} horizontal well pairs across all common elevation slices.")

    df_matrices = analyzer.compute_lag_binned_horizontal_matrices()
    mat_path = OUTPUT_DIR / "empirical_horizontal_transition_matrices.csv"
    df_matrices.to_csv(mat_path, index=False)
    print(f"Saved lag-binned horizontal transition matrices: {mat_path} ({len(df_matrices)} rows)")

    # Produce transition probability decay curves vs distance
    lags = np.linspace(0.0, 5000.0, 100)
    R = analyzer.build_theoretical_horizontal_rate_matrix()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

    # Plot theoretical auto-transition probability decay P_ii(h)
    for code in range(len(CANONICAL_FACIES_SCHEMA)):
        color, name = FACIES_PALETTE[code]
        p_ii = [analyzer.evaluate_transition_probability_at_lag(h, R)[code, code] for h in lags]
        ax1.plot(lags, p_ii, color=color, linewidth=2.2, label=f"{name} (L={DEFAULT_LATERAL_FACIES_LENGTHS_M[code]:.0f}m)")

    ax1.set_xlabel("Horizontal Lag Distance h (meters)", fontsize=10.5, fontweight="bold")
    ax1.set_ylabel("Auto-Transition Probability P(i, i; h)", fontsize=10.5, fontweight="bold")
    ax1.set_title("Theoretical Horizontal Auto-Transition Decay P(i, i; h)\nCarle & Fogg (1996) Rate Formulation: P(h) = expm(R*h)", fontsize=11.5, fontweight="bold")
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="upper right", fontsize=8.5)
    ax1.set_ylim(-0.02, 1.05)

    # Plot empirical pair matches vs distance
    # Group empirical pairs into distance bins of 250m
    bin_edges = np.arange(400.0, 5250.0, 350.0)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    empirical_same_facies_prop = []
    pair_counts = []

    for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
        sub = df_pairs[(df_pairs["lag_distance_m"] >= lo) & (df_pairs["lag_distance_m"] < hi)]
        if len(sub) > 0:
            empirical_same_facies_prop.append(float(sub["is_same_facies"].mean()))
            pair_counts.append(len(sub))
        else:
            empirical_same_facies_prop.append(np.nan)
            pair_counts.append(0)

    ax2.plot(bin_centers, empirical_same_facies_prop, marker="o", color="#B22222", linewidth=2.0, label="Empirical Inter-Well Pairs (Common Datum)")
    ax2.axhline(0.35, color="#555555", linestyle="--", label="Expected Stationary Coincidence (~35%)")
    ax2.axvline(420.0, color="#2E8B57", linestyle=":", label="Min Inter-Well Spacing (420 m)")

    ax2.set_xlabel("Inter-Well Distance h (meters)", fontsize=10.5, fontweight="bold")
    ax2.set_ylabel("Proportion of Identical Facies Pairs", fontsize=10.5, fontweight="bold")
    ax2.set_title("Empirical Spatial Facies Coincidence vs Inter-Well Distance\nObservations Across Matched Common-Zero Elevation Slices", fontsize=11.5, fontweight="bold")
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right", fontsize=8.5)
    ax2.set_ylim(0.0, 0.7)

    plt.tight_layout()
    decay_fig_path = FIGURES_DIR / "horizontal_transition_probability_decay.png"
    plt.savefig(decay_fig_path, dpi=300)
    plt.close()
    print(f"Saved transition probability decay plot: {decay_fig_path}")


def step4_spatial_markov_lolo_evaluation(predictor: SpatialMarkovPredictor):
    print("--- Step 4: Spatial Markov LOLO Evaluation on Common Zero ---")
    res_markov = predictor.run_spatial_markov_lolo()
    df_lolo = res_markov["per_fold_table"]
    lolo_path = OUTPUT_DIR / "spatial_markov_lolo_results.csv"
    df_lolo.to_csv(lolo_path, index=False)
    print(f"Saved Spatial Markov LOLO results: {lolo_path} ({len(df_lolo)} rows)")

    # Compare with Spatial 3D KNN under common-zero datum
    validator = ProvisionalSpatialValidator(vertical_reference="common_zero")
    res_knn = validator.run_spatial_lolo(k_neighbors=5)
    summary_knn = res_knn["aggregate_summary"]

    summary_markov = res_markov["aggregate_summary"]

    comparison_records = [
        {
            "model_name": "Spatial Markov Transition Model",
            "model_type": "Continuous Geostatistical Transition Model (P(h) = expm(R*h))",
            "vertical_datum": "common_zero_datum (Prof. Sahoo standard)",
            "pooled_accuracy": summary_markov["pooled_markov_accuracy"],
            "pooled_balanced_acc": summary_markov["pooled_markov_balanced_acc"],
            "pooled_macro_f1": summary_markov["pooled_markov_macro_f1"],
            "conditioning_principle": "Continuous lateral transition decay conditioned on Sahoo et al. (2016) W/T aspect ratios",
        },
        {
            "model_name": "Spatial 3D KNN (k=5)",
            "model_type": "Spatial 3D Classifier",
            "vertical_datum": "common_zero_datum (Prof. Sahoo standard)",
            "pooled_accuracy": summary_knn["overall_knn_accuracy"],
            "pooled_balanced_acc": summary_knn["overall_knn_balanced_acc"],
            "pooled_macro_f1": summary_knn["overall_knn_macro_f1"],
            "conditioning_principle": "Inverse-distance weighted point-wise interpolation in (X, Y, Z_common)",
        },
        {
            "model_name": "Nearest-Well Vertical Profile",
            "model_type": "1D Transferred Profile Baseline (Zero Horizontal Decay)",
            "vertical_datum": "common_zero_datum (Prof. Sahoo standard)",
            "pooled_accuracy": summary_markov["pooled_near_well_accuracy"],
            "pooled_balanced_acc": summary_markov["pooled_near_well_balanced_acc"],
            "pooled_macro_f1": summary_markov["pooled_near_well_macro_f1"],
            "conditioning_principle": "Direct copy of nearest well facies at matching common elevation slice",
        },
        {
            "model_name": "Training Prior Majority Facies",
            "model_type": "Zero-Spatial Naive Baseline (Infinite Horizontal Decay)",
            "vertical_datum": "common_zero_datum (Prof. Sahoo standard)",
            "pooled_accuracy": summary_markov["pooled_prior_accuracy"],
            "pooled_balanced_acc": summary_markov["pooled_prior_balanced_acc"],
            "pooled_macro_f1": summary_markov["pooled_prior_macro_f1"],
            "conditioning_principle": "Always predicts training majority facies (Channel Sandstone)",
        },
    ]

    df_comp = pd.DataFrame(comparison_records)
    comp_path = OUTPUT_DIR / "spatial_markov_vs_baselines_comparison.csv"
    df_comp.to_csv(comp_path, index=False)
    print(f"Saved baseline comparison table: {comp_path} ({len(df_comp)} rows)")
    print("\nSummary Results:")
    print(df_comp[["model_name", "pooled_accuracy", "pooled_balanced_acc", "pooled_macro_f1"]].to_string(index=False))


if __name__ == "__main__":
    print("===============================================================")
    print("Executing SMALT Sprint H Master Pipeline...")
    print("Common-Zero Datum Alignment & Spatial Markov Foundation")
    print("===============================================================")

    inspector = LithologInspector()
    df_meta, df_inv = step1_datum_alignment_and_validation(inspector)
    step2_generate_diagnostic_alignment_plot(inspector, df_meta)

    analyzer = SpatialMarkovTransitionAnalyzer(inspector=inspector)
    step3_spatial_transition_data_construction(analyzer)

    predictor = SpatialMarkovPredictor(analyzer=analyzer)
    step4_spatial_markov_lolo_evaluation(predictor)

    print("\n===============================================================")
    print("Sprint H Master Pipeline Completed Successfully.")
    print("===============================================================")
