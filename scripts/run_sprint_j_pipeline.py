"""
SMALT Sprint J Master Pipeline: Empirical Horizontal Calibration & Continuity Audit.

Executes all 28 requirements for Sprint J:
1. Recomputes empirical horizontal transitions (4,410 pairs across 93 common-zero slices).
2. Evaluates multiple distance binning strategies (Fixed, Quantile, Sensitivity) with bootstrap standard errors.
3. Estimates facies-specific lateral continuity lengths (L_geological, L_empirical, L_fitted).
4. Tests alternative decay models (Model A: Sprint I, Model B: Calibrated Markov, Model C: Spherical Decay).
5. Evaluates spatial decay at critical distances (420m, 700m, 1000m, 2000m, 5000m) vs stationary proportions.
6. Conducts mandatory 4-model ablation study (Model 0, 1, 2, 3).
7. Executes leak-free LOLO cross-validation across L2-L12 (940 points).
8. Runs directional Upstream <-> Downstream validation.
9. Runs controlled synthetic length recovery benchmark.
10. Generates all 9 required CSV deliverables in sprints/audit_sprint_j/.
11. Generates all 11 publication-quality figures in sprints/audit_sprint_j/figures/.
"""

import sys
import os
from pathlib import Path

# Ensure repository root is on sys.path
repo_root = Path(__file__).resolve().parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

from smalt.descriptive.analyzer import LithologInspector, CANONICAL_FACIES_SCHEMA
from smalt.geostat.empirical_calibration import (
    HorizontalContinuityCalibrator,
    CODE_TO_NAME,
    DEFAULT_FIXED_BINS,
)
from smalt.geostat.ablation import AblationStudyEngine
from smalt.validation.sprint_j import SprintJValidator

# Output directory structure
OUTPUT_DIR = Path("sprints/audit_sprint_j")
FIG_DIR = OUTPUT_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Facies palette
FACIES_COLORS = [meta["color_hex"] for meta in CANONICAL_FACIES_SCHEMA.values()]
FACIES_CMAP = mcolors.ListedColormap(FACIES_COLORS)


def run_sprint_j_master():
    print("=" * 80)
    print("SMALT SPRINT J: EMPIRICAL HORIZONTAL TRANSITION & CONTINUITY CALIBRATION")
    print("=" * 80)

    inspector = LithologInspector()
    calibrator = HorizontalContinuityCalibrator(inspector=inspector)

    # -------------------------------------------------------------------------
    # Step 6 & 7: Empirical Horizontal Transitions & Multiple Binning Strategies
    # -------------------------------------------------------------------------
    print("\n[Step 1/11] Computing empirical horizontal pairs and evaluating binning strategies...")
    df_pairs = calibrator.extract_pairs()
    print(f"  Extracted {len(df_pairs)} matched pairs across {len(df_pairs['z_common_m'].unique())} elevation slices.")
    df_pairs.to_csv(OUTPUT_DIR / "empirical_horizontal_transitions.csv", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'empirical_horizontal_transitions.csv'}")

    bin_dict = calibrator.evaluate_binning_strategies(n_bootstraps=200, random_seed=42)
    fixed_bins_df = bin_dict["fixed_bins"]

    # -------------------------------------------------------------------------
    # Step 8 & 9: Facies-Specific Lateral Continuity Estimation (3 Lengths)
    # -------------------------------------------------------------------------
    print("\n[Step 2/11] Estimating facies-specific lateral continuity lengths...")
    lengths_df = calibrator.estimate_facies_lengths(binned_df=fixed_bins_df)
    lengths_df.to_csv(OUTPUT_DIR / "facies_horizontal_length_comparison.csv", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'facies_horizontal_length_comparison.csv'}")
    print(lengths_df[["facies_name", "L_geological_m", "L_empirical_m", "L_fitted_m", "auto_pair_count", "identifiability_status"]].to_string())

    # -------------------------------------------------------------------------
    # Step 10: Alternative Decay Models & Fit Statistics
    # -------------------------------------------------------------------------
    print("\n[Step 3/11] Testing alternative horizontal decay models (A, B, C)...")
    stats_df, models_dict = calibrator.test_alternative_decay_models()
    stats_df.to_csv(OUTPUT_DIR / "transition_model_fit_statistics.csv", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'transition_model_fit_statistics.csv'}")
    print(stats_df.to_string())

    # -------------------------------------------------------------------------
    # Step 11 & 12: Distance Decay & Model Comparison at Observed Distances
    # -------------------------------------------------------------------------
    print("\n[Step 4/11] Evaluating model decay at observed key well distances...")
    decay_df = calibrator.evaluate_decay_at_key_distances(models_dict)
    decay_df.to_csv(OUTPUT_DIR / "current_vs_empirical_transition.csv", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'current_vs_empirical_transition.csv'}")

    # -------------------------------------------------------------------------
    # Step 13: Mandatory Ablation Study (Models 0, 1, 2, 3)
    # -------------------------------------------------------------------------
    print("\n[Step 5/11] Running 4-model ablation study (Models 0, 1, 2, 3)...")
    ablation_engine = AblationStudyEngine(
        anchor_well_ids=["litholog2", "litholog9"],
        inspector=inspector,
        random_seed=42,
    )
    ablation_res = ablation_engine.run_ablation_realizations(n_realizations=5)
    ablation_df = ablation_res["ablation_summary"]
    ablation_df.to_csv(OUTPUT_DIR / "ablation_results.csv", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'ablation_results.csv'}")
    print(ablation_df[["model_key", "proportion_tv_distance_mean", "mean_bed_thickness_m", "bed_thickness_error_m", "channel_sand_connectivity_m", "hard_data_honor_rate"]].to_string())

    # -------------------------------------------------------------------------
    # Step 14 & 15: Leave-One-Litholog-Out (LOLO) Cross-Validation
    # -------------------------------------------------------------------------
    print("\n[Step 6/11] Executing leak-free LOLO cross-validation across L2-L12...")
    validator = SprintJValidator(inspector=inspector)
    lolo_res = validator.run_lolo_cross_validation()
    lolo_res["summary_comparison"].to_csv(OUTPUT_DIR / "lolo_results.csv", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'lolo_results.csv'}")
    print(lolo_res["summary_comparison"].to_string())

    # -------------------------------------------------------------------------
    # Step 14B & 14C: Directional Upstream <-> Downstream Validation
    # -------------------------------------------------------------------------
    print("\n[Step 7/11] Running transport-parallel directional validation...")
    directional_df = validator.run_directional_validation()
    up_down_df = directional_df[directional_df["experiment_direction"] == "Upstream -> Downstream"].copy()
    down_up_df = directional_df[directional_df["experiment_direction"] == "Downstream -> Upstream"].copy()

    up_down_df.to_csv(OUTPUT_DIR / "upstream_downstream_results.csv", index=False)
    down_up_df.to_csv(OUTPUT_DIR / "downstream_upstream_results.csv", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'upstream_downstream_results.csv'}")
    print(f"  Saved: {OUTPUT_DIR / 'downstream_upstream_results.csv'}")
    print(directional_df.to_string())

    # -------------------------------------------------------------------------
    # Step 21: Synthetic Length Recovery Benchmark
    # -------------------------------------------------------------------------
    print("\n[Step 8/11] Running synthetic length recovery benchmark...")
    syn_df = validator.run_synthetic_benchmark(random_seed=42)
    syn_df.to_csv(OUTPUT_DIR / "synthetic_length_recovery.csv", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'synthetic_length_recovery.csv'}")
    print(syn_df.to_string())

    # -------------------------------------------------------------------------
    # Step 19: Spatial Sampling Diagnostic
    # -------------------------------------------------------------------------
    print("\n[Step 9/11] Running spatial sampling diagnostic (inter-well spacing vs facies lengths)...")
    diag_res = calibrator.analyze_spatial_sampling_density()
    spacing_stats = diag_res["inter_well_spacing_m"]
    print(f"  Inter-well spacing: Min={spacing_stats['min']:.1f}m, Median={spacing_stats['median']:.1f}m, Max={spacing_stats['max']:.1f}m")

    # -------------------------------------------------------------------------
    # Step 22: Generating 11 Publication Figures
    # -------------------------------------------------------------------------
    print("\n[Step 10/11] Generating all 11 publication-quality scientific figures...")
    generate_all_figures(
        calibrator=calibrator,
        fixed_bins_df=fixed_bins_df,
        lengths_df=lengths_df,
        models_dict=models_dict,
        decay_df=decay_df,
        ablation_res=ablation_res,
        lolo_res=lolo_res,
        directional_df=directional_df,
        spacing_stats=spacing_stats,
    )

    print("\n[Step 11/11] All 9 CSV deliverables and 11 figures successfully generated!")
    print("=" * 80)


def generate_all_figures(
    calibrator: HorizontalContinuityCalibrator,
    fixed_bins_df: pd.DataFrame,
    lengths_df: pd.DataFrame,
    models_dict: Dict[str, Any],
    decay_df: pd.DataFrame,
    ablation_res: Dict[str, Any],
    lolo_res: Dict[str, Any],
    directional_df: pd.DataFrame,
    spacing_stats: Dict[str, float],
):
    plt.rcParams["font.sans-serif"] = "DejaVu Sans"
    plt.rcParams["font.size"] = 10

    # -------------------------------------------------------------------------
    # Figure 1: Empirical Horizontal Transition Probability vs Distance
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    for c in range(6):
        c_bins = fixed_bins_df[
            (fixed_bins_df["from_facies_code"] == c) &
            (fixed_bins_df["to_facies_code"] == c) &
            (fixed_bins_df["pair_count"] > 0)
        ].sort_values("midpoint_m")
        if len(c_bins) > 0:
            h_m = c_bins["midpoint_m"].to_numpy()
            p_val = c_bins["transition_probability"].to_numpy()
            se = c_bins["std_error"].to_numpy()
            ax.errorbar(
                h_m, p_val, yerr=se,
                label=CODE_TO_NAME[c],
                color=FACIES_COLORS[c],
                marker="o", linewidth=2, capsize=4,
            )
    ax.set_title("Figure 1: Empirical Horizontal Auto-Transition Probability vs Distance (with Bootstrap SE)", fontweight="bold")
    ax.set_xlabel("Horizontal Separation Distance (m)")
    ax.set_ylabel("Auto-Transition Probability P_ii(h)")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="upper right", frameon=True)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "empirical_horizontal_transition_vs_distance.png", dpi=300)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Figure 2: Current Sprint I Decay vs Empirical Observations
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    model_A = models_dict["model_A"]
    h_smooth = np.linspace(0.0, 5000.0, 200)
    for c in [0, 5]:  # Focus on major facies: Sand and Mud
        c_bins = fixed_bins_df[
            (fixed_bins_df["from_facies_code"] == c) &
            (fixed_bins_df["to_facies_code"] == c) &
            (fixed_bins_df["pair_count"] > 0)
        ].sort_values("midpoint_m")
        ax.plot(
            c_bins["midpoint_m"], c_bins["transition_probability"],
            "o", color=FACIES_COLORS[c], label=f"Empirical {CODE_TO_NAME[c]}", markersize=7
        )
        p_model_a = [model_A.evaluate_transition_matrix(h)[c, c] for h in h_smooth]
        ax.plot(
            h_smooth, p_model_a,
            "--", color=FACIES_COLORS[c], linewidth=2.5,
            label=f"Sprint I Model A ({CODE_TO_NAME[c]}, L={lengths_df.loc[c, 'L_geological_m']}m)"
        )
    ax.axvline(x=420.0, color="gray", linestyle=":", label="Min Inter-Well Spacing (420m)")
    ax.set_title("Figure 2: Sprint I Prescribed Decay vs Empirical Observations (Premature Decay)", fontweight="bold")
    ax.set_xlabel("Horizontal Separation Distance (m)")
    ax.set_ylabel("Auto-Transition Probability P_ii(h)")
    ax.set_ylim(0.2, 1.05)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="upper right", frameon=True)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "current_sprint_i_decay_vs_empirical.png", dpi=300)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Figure 3: Sprint J Fitted Decay vs Empirical Observations
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    model_B = models_dict["model_B"]
    for c in [0, 5]:
        c_bins = fixed_bins_df[
            (fixed_bins_df["from_facies_code"] == c) &
            (fixed_bins_df["to_facies_code"] == c) &
            (fixed_bins_df["pair_count"] > 0)
        ].sort_values("midpoint_m")
        ax.errorbar(
            c_bins["midpoint_m"], c_bins["transition_probability"], yerr=c_bins["std_error"],
            fmt="o", color=FACIES_COLORS[c], label=f"Empirical {CODE_TO_NAME[c]} (+/- SE)", markersize=7, capsize=4
        )
        p_model_b = [model_B.evaluate_transition_matrix(h)[c, c] for h in h_smooth]
        ax.plot(
            h_smooth, p_model_b,
            "-", color=FACIES_COLORS[c], linewidth=2.5,
            label=f"Sprint J Calibrated Model B ({CODE_TO_NAME[c]}, L={lengths_df.loc[c, 'L_fitted_m']}m)"
        )
    ax.axvline(x=420.0, color="gray", linestyle=":", label="Min Inter-Well Spacing (420m)")
    ax.set_title("Figure 3: Sprint J Calibrated Horizontal Decay vs Empirical Observations", fontweight="bold")
    ax.set_xlabel("Horizontal Separation Distance (m)")
    ax.set_ylabel("Auto-Transition Probability P_ii(h)")
    ax.set_ylim(0.2, 1.05)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="upper right", frameon=True)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "sprint_j_fitted_decay_vs_empirical.png", dpi=300)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Figure 4: Facies-Specific Correlation-Length Comparison
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(11, 6))
    x_idx = np.arange(6)
    width = 0.25
    ax.bar(x_idx - width, lengths_df["L_geological_m"], width, label="L_geological (Geometric Prior)", color="#5DADE2")
    ax.bar(x_idx, lengths_df["L_empirical_m"], width, label="L_empirical (e-folding)", color="#F4D03F")
    ax.bar(x_idx + width, lengths_df["L_fitted_m"], width, label="L_fitted (Continuous Rate Fit)", color="#E74C3C")
    ax.set_xticks(x_idx)
    ax.set_xticklabels([CODE_TO_NAME[c] for c in range(6)], rotation=15, ha="right")
    ax.set_ylabel("Characteristic Lateral Length (m)")
    ax.set_title("Figure 4: Facies-Specific Lateral Continuity Length Comparison", fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.5, axis="y")
    ax.legend(loc="upper right", frameon=True)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "facies_specific_correlation_length_comparison.png", dpi=300)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Figure 5: Well Spacing vs Estimated Facies Continuity
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    y_pos = np.arange(6)
    ax.barh(y_pos, lengths_df["L_geological_m"], color=FACIES_COLORS, alpha=0.8, edgecolor="black", label="Facies Lateral Dimension")
    ax.axvline(x=spacing_stats["min"], color="red", linestyle="--", linewidth=2.5, label=f"Min Well Spacing ({spacing_stats['min']:.0f}m)")
    ax.axvline(x=spacing_stats["median"], color="black", linestyle="-", linewidth=2.5, label=f"Median Well Spacing ({spacing_stats['median']:.0f}m)")
    ax.set_yticks(y_pos)
    ax.set_yticklabels([CODE_TO_NAME[c] for c in range(6)])
    ax.set_xlabel("Horizontal Distance (m)")
    ax.set_title("Figure 5: Inter-Well Spacing vs Facies Dimensions (Sub-Grid Sparsity Diagnostic)", fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.5, axis="x")
    ax.legend(loc="lower right", frameon=True)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "well_spacing_vs_facies_length.png", dpi=300)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Figure 6: Transition Matrices at Selected Distances
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    dists_plot = [420.0, 1000.0, 2500.0, 5000.0]
    for idx, (dist, ax_cur) in enumerate(zip(dists_plot, axes.flatten())):
        P_mat = models_dict["model_B"].evaluate_transition_matrix(dist)
        im = ax_cur.imshow(P_mat, cmap="YlGnBu", vmin=0.0, vmax=1.0)
        ax_cur.set_title(f"Lag h = {dist:.0f} m", fontweight="bold")
        ax_cur.set_xticks(range(6))
        ax_cur.set_yticks(range(6))
        labels = [meta["canonical_name"].split()[0] for meta in CANONICAL_FACIES_SCHEMA.values()]
        ax_cur.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax_cur.set_yticklabels(labels, fontsize=8)
        for i in range(6):
            for j in range(6):
                ax_cur.text(j, i, f"{P_mat[i, j]:.2f}", ha="center", va="center", color="black" if P_mat[i, j] < 0.6 else "white", fontsize=8)
    fig.suptitle("Figure 6: Continuous Horizontal Transition Matrices at Selected Well Distances", fontweight="bold", fontsize=14)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "transition_matrices_at_selected_distances.png", dpi=300)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Figure 7: Ablation-Model Comparison
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    ab_df = ablation_res["ablation_summary"]
    labels_ab = ["M0: Stationary", "M1: Vert Markov", "M2: Vert + Sp I", "M3: Vert + Sp J"]
    ax1.bar(labels_ab, ab_df["bed_thickness_error_m"], color=["#BDC3C7", "#2ECC71", "#E67E22", "#3498DB"], edgecolor="black")
    ax1.set_title("Mean Bed Thickness Error vs Ground Truth", fontweight="bold")
    ax1.set_ylabel("Absolute Error (m)")
    ax1.grid(True, linestyle="--", alpha=0.5, axis="y")

    ax2.bar(labels_ab, ab_df["proportion_tv_distance_mean"], color=["#BDC3C7", "#2ECC71", "#E67E22", "#3498DB"], edgecolor="black")
    ax2.set_title("Facies Proportion Total Variation Distance", fontweight="bold")
    ax2.set_ylabel("TV Distance")
    ax2.grid(True, linestyle="--", alpha=0.5, axis="y")
    plt.tight_layout()
    fig.savefig(FIG_DIR / "ablation_model_comparison.png", dpi=300)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Figure 8: Stochastic Realization Comparison
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(4, 1, figsize=(14, 12), sharex=True)
    all_reals = ablation_res["models_realizations"]
    m_keys = [
        "Model_0_Stationary_Proportions",
        "Model_1_Vertical_Markov_Only",
        "Model_2_Vertical_Plus_Sprint_I",
        "Model_3_Vertical_Plus_Sprint_J",
    ]
    grid_meta = ablation_res["grid_metadata"]
    x_coords = grid_meta["x_coords_m"]
    z_coords = grid_meta["z_coords_m"]

    for idx, (m_k, ax_cur) in enumerate(zip(m_keys, axes)):
        grid_data = all_reals[m_k]["realizations"][0]
        im = ax_cur.imshow(
            grid_data, cmap=FACIES_CMAP, origin="lower",
            extent=[x_coords[0], x_coords[-1], z_coords[0], z_coords[-1]],
            aspect="auto", vmin=0, vmax=5,
        )
        ax_cur.set_title(all_reals[m_k]["label"], fontweight="bold", fontsize=11)
        ax_cur.set_ylabel("Common-Zero Elev (m)")
        for pos in grid_meta["well_positions_x_m"]:
            ax_cur.axvline(x=pos, color="black", linestyle="--", alpha=0.8)

    axes[-1].set_xlabel("Distance Along Transect (m)")
    plt.tight_layout()
    fig.savefig(FIG_DIR / "stochastic_realization_comparison.png", dpi=300)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Figure 9: Uncertainty Comparison (Entropy and Ensemble Variation)
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(labels_ab, ab_df["shannon_entropy_mean"], color=["#34495E", "#1ABC9C", "#9B59B6", "#2980B9"], edgecolor="black")
    ax.set_title("Figure 9: Spatial Shannon Entropy Across Ablation Architectures", fontweight="bold")
    ax.set_ylabel("Mean Shannon Entropy (nats)")
    ax.grid(True, linestyle="--", alpha=0.5, axis="y")
    plt.tight_layout()
    fig.savefig(FIG_DIR / "uncertainty_comparison.png", dpi=300)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Figure 10: Upstream -> Downstream Results
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    up_df = directional_df[directional_df["experiment_direction"] == "Upstream -> Downstream"]
    bars = ax.bar(up_df["model_name"], up_df["accuracy"] * 100.0, color=["#95A5A6", "#3498DB", "#2ECC71", "#E74C3C"], edgecolor="black")
    ax.axhline(y=up_df.iloc[0]["accuracy"] * 100.0, color="gray", linestyle="--", label="Training Prior")
    ax.set_title("Figure 10: Upstream -> Downstream Generalization Accuracy (%)", fontweight="bold")
    ax.set_ylabel("Test Accuracy (%)")
    ax.set_ylim(0, 60)
    ax.grid(True, linestyle="--", alpha=0.5, axis="y")
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2.0, yval + 1.0, f"{yval:.1f}%", ha="center", va="bottom", fontweight="bold")
    plt.tight_layout()
    fig.savefig(FIG_DIR / "upstream_to_downstream_results.png", dpi=300)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Figure 11: Downstream -> Upstream Results
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    down_df = directional_df[directional_df["experiment_direction"] == "Downstream -> Upstream"]
    bars = ax.bar(down_df["model_name"], down_df["accuracy"] * 100.0, color=["#95A5A6", "#3498DB", "#2ECC71", "#E74C3C"], edgecolor="black")
    ax.axhline(y=down_df.iloc[0]["accuracy"] * 100.0, color="gray", linestyle="--", label="Training Prior")
    ax.set_title("Figure 11: Downstream -> Upstream Generalization Accuracy (%)", fontweight="bold")
    ax.set_ylabel("Test Accuracy (%)")
    ax.set_ylim(0, 60)
    ax.grid(True, linestyle="--", alpha=0.5, axis="y")
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2.0, yval + 1.0, f"{yval:.1f}%", ha="center", va="bottom", fontweight="bold")
    plt.tight_layout()
    fig.savefig(FIG_DIR / "downstream_to_upstream_results.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    run_sprint_j_master()
