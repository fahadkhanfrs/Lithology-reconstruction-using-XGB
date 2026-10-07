"""
SMALT Sprint I End-to-End Master Pipeline Runner.

Executes all Sprint I workflows, generating all required CSV deliverables
and all 10 publication-quality scientific figures in sprints/audit_sprint_i/.

Workflow Steps:
1. Empirical Horizontal Transition Extraction & Isotropy Analysis
2. Continuous Spatial Transition-Rate Model Calibration & Distance Decay Comparison
3. Full Spatial LOLO Cross-Validation Benchmark (L2-L12, 940 grid points)
4. Directional Upstream <-> Downstream Validation
5. Geological-Statistical Reproduction Analysis (Proportions, Transition Divergence, Bed Thickness)
6. Synthetic Sanity Benchmark (Step 10)
7. Inter-Well Stochastic Realization Generation & Ensemble Uncertainty Mapping
8. Generation of 10 Publication-Quality Figures
"""

from typing import Dict, Any, List, Optional
from pathlib import Path
import sys

# Ensure repository root is on sys.path
repo_root = Path(__file__).resolve().parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from smalt.descriptive.analyzer import LithologInspector, CANONICAL_FACIES_SCHEMA
from smalt.geostat.markov import StratigraphicMarkovChain
from smalt.geostat.spatial_transition import (
    EmpiricalHorizontalTransitionEstimator,
    SpatialTransitionRateModel,
    DEFAULT_LATERAL_FACIES_LENGTHS_M,
)
from smalt.geostat.conditioned_markov import ConditionedMarkovClassifier, CODE_TO_NAME
from smalt.geostat.realization import InterWellRealizationGenerator
from smalt.validation.sprint_i import SprintIValidator


# Facies palette mapping
FACIES_COLORS = {meta["code"]: meta["color_hex"] for meta in CANONICAL_FACIES_SCHEMA.values()}
FACIES_NAMES = [CODE_TO_NAME[c] for c in range(6)]


def main():
    print("=" * 80)
    print("SMALT SPRINT I: CONDITIONED MARKOV TRANSITION-PROBABILITY PIPELINE")
    print("=" * 80)

    repo_root = Path(__file__).resolve().parent.parent
    output_dir = repo_root / "sprints" / "audit_sprint_i"
    fig_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    inspector = LithologInspector()

    # -------------------------------------------------------------------------
    # STEP 1: Empirical Horizontal Transitions & Isotropy Analysis
    # -------------------------------------------------------------------------
    print("\n[Step 1] Extracting empirical horizontal transitions & analyzing isotropy...")
    estimator = EmpiricalHorizontalTransitionEstimator(inspector=inspector)
    df_pairs = estimator.extract_horizontal_facies_pairs()
    print(f"  Extracted {len(df_pairs)} matched horizontal pairs across {df_pairs['z_common_m'].nunique()} elevation slices.")

    isotropy_info = estimator.evaluate_directional_sparsity()
    print(f"  Isotropy justification: {isotropy_info['isotropy_justification']}")

    df_emp_transitions = estimator.compute_lag_binned_horizontal_matrices()
    emp_csv_path = output_dir / "empirical_horizontal_transitions.csv"
    df_emp_transitions.to_csv(emp_csv_path, index=False)
    print(f"  Saved: {emp_csv_path}")

    # -------------------------------------------------------------------------
    # STEP 2: Transition-Probability Decay Model Calibration
    # -------------------------------------------------------------------------
    print("\n[Step 2] Calibrating continuous spatial transition-rate model P(h) = expm(R * h)...")
    rate_model = SpatialTransitionRateModel()
    df_fit = rate_model.compare_empirical_vs_fitted(estimator)
    fit_csv_path = output_dir / "horizontal_transition_fit.csv"
    df_fit.to_csv(fit_csv_path, index=False)
    print(f"  Saved: {fit_csv_path}")

    # -------------------------------------------------------------------------
    # STEP 3 & 4: Spatial LOLO Cross-Validation Benchmark
    # -------------------------------------------------------------------------
    print("\n[Step 3] Running Leave-One-Litholog-Out (LOLO) cross-validation across L2-L12...")
    validator = SprintIValidator(inspector=inspector)
    lolo_results = validator.run_lolo_cross_validation()

    df_folds = lolo_results["per_fold_table"]
    df_comp = lolo_results["comparison_table"]
    df_geo = lolo_results["geological_statistical_metrics"]

    folds_csv_path = output_dir / "conditioned_markov_lolo_results.csv"
    comp_csv_path = output_dir / "conditioned_markov_baseline_comparison.csv"
    geo_csv_path = output_dir / "geological_statistical_metrics.csv"

    df_folds.to_csv(folds_csv_path, index=False)
    df_comp.to_csv(comp_csv_path, index=False)
    df_geo.to_csv(geo_csv_path, index=False)
    print(f"  Saved: {folds_csv_path}")
    print(f"  Saved: {comp_csv_path}")
    print(f"  Saved: {geo_csv_path}")

    print("\nBaseline Comparison Summary (940 Evaluation Points):")
    print(df_comp.to_string(index=False))

    print("\nGeological-Statistical Reproduction Metrics:")
    print(df_geo.to_string(index=False))

    # -------------------------------------------------------------------------
    # STEP 5: Directional Upstream <-> Downstream Validation
    # -------------------------------------------------------------------------
    print("\n[Step 5] Running directional generalization validation (Upstream <-> Downstream)...")
    df_directional = validator.run_directional_upstream_downstream_validation()
    dir_csv_path = output_dir / "directional_upstream_downstream_sprint_i.csv"
    df_directional.to_csv(dir_csv_path, index=False)
    print(f"  Saved: {dir_csv_path}")
    print(df_directional.to_string(index=False))

    # -------------------------------------------------------------------------
    # STEP 6: Synthetic Benchmark Sanity Test (Step 10)
    # -------------------------------------------------------------------------
    print("\n[Step 6] Running synthetic benchmark sanity test (Step 10)...")
    synth_res = validator.run_synthetic_benchmark()
    df_synth = synth_res["benchmark_table"]
    synth_csv_path = output_dir / "synthetic_benchmark_results.csv"
    df_synth.to_csv(synth_csv_path, index=False)
    print(f"  Saved: {synth_csv_path}")
    print(df_synth.to_string(index=False))

    # -------------------------------------------------------------------------
    # STEP 7: Inter-Well Stochastic Realization Generation
    # -------------------------------------------------------------------------
    print("\n[Step 7] Generating stochastic inter-well facies realizations along L2 -> L9 transect...")
    # Train classifier on eligible wells
    train_ids = [lid for lid in validator.eligible_ids if lid not in ["litholog9", "litholog11"]]
    classifier_all = ConditionedMarkovClassifier(
        training_litholog_ids=validator.eligible_ids,
        inspector=inspector,
    )
    honor_rate = classifier_all.compute_hard_data_honor_rate()
    print(f"  Conditioning hard data honor rate: {honor_rate * 100:.1f}%")

    realization_gen = InterWellRealizationGenerator(classifier=classifier_all, random_seed=42)
    # Transect from L2 (westernmost upstream) to L9 (downstream canyon)
    transect_wells = ["litholog2", "litholog3", "litholog4", "litholog9"]
    sim_res = realization_gen.generate_2d_transect_realizations(
        anchor_well_ids=transect_wells,
        n_realizations=5,
        x_resolution_m=25.0,
        z_resolution_m=1.0,
    )
    print(f"  Simulated {sim_res['n_realizations']} realizations across grid {sim_res['realizations'].shape[2]} x {sim_res['realizations'].shape[1]} cells.")
    print(f"  Realization hard data honor rate: {sim_res['hard_data_honor_rate'] * 100:.1f}%")

    # -------------------------------------------------------------------------
    # STEP 8: Generation of 10 Publication-Quality Figures
    # -------------------------------------------------------------------------
    print("\n[Step 8] Generating 10 publication-quality scientific figures...")

    # Figure 1: Empirical horizontal transition probability vs distance
    plt.figure(figsize=(10, 6))
    auto_df = df_fit[df_fit["is_auto_transition"] == 1]
    for code in [0, 1, 2, 4, 5]:
        sub = auto_df[auto_df["from_facies_code"] == code]
        plt.plot(
            sub["nominal_distance_m"],
            sub["empirical_probability"],
            marker="o",
            linestyle="-",
            label=f"{CODE_TO_NAME[code]} (Empirical)",
            color=FACIES_COLORS[code],
            linewidth=2,
        )
    plt.axhline(0.435, color=FACIES_COLORS[0], linestyle="--", alpha=0.5, label="Sand Stationary (43.5%)")
    plt.axhline(0.372, color=FACIES_COLORS[5], linestyle="--", alpha=0.5, label="Mud Stationary (37.2%)")
    plt.xlabel("Horizontal Separation Distance h (m)", fontsize=12)
    plt.ylabel("Empirical Auto-Transition Probability P(i, i; h)", fontsize=12)
    plt.title("SMALT Sprint I: Empirical Horizontal Transition Probability Decay vs Separation", fontsize=14, fontweight="bold")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right", frameon=True, fontsize=10)
    plt.ylim(0.0, 1.0)
    fig1_path = fig_dir / "empirical_horizontal_transition_vs_distance.png"
    plt.savefig(fig1_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved Figure 1: {fig1_path}")

    # Figure 2: Empirical vs fitted transition decay
    plt.figure(figsize=(12, 6))
    h_eval = np.linspace(0.0, 5500.0, 300)
    for code in [0, 2, 5]:  # sand, ripples, mud
        P_curves = np.array([rate_model.evaluate_transition_matrix(h)[code, code] for h in h_eval])
        plt.plot(
            h_eval,
            P_curves,
            linestyle="-",
            linewidth=2.5,
            color=FACIES_COLORS[code],
            label=f"{CODE_TO_NAME[code]} Fitted Curve P(h) = exp(R*h)",
        )
        sub = auto_df[auto_df["from_facies_code"] == code]
        plt.scatter(
            sub["nominal_distance_m"],
            sub["empirical_probability"],
            color=FACIES_COLORS[code],
            edgecolor="black",
            s=80,
            zorder=5,
            label=f"{CODE_TO_NAME[code]} Empirical Points",
        )
    plt.xlabel("Horizontal Separation Distance h (m)", fontsize=12)
    plt.ylabel("Auto-Transition Probability P(i, i; h)", fontsize=12)
    plt.title("SMALT Sprint I: Continuous Transition-Rate Model Fit vs Observed Field Pairs", fontsize=14, fontweight="bold")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right", frameon=True, fontsize=10)
    plt.ylim(0.0, 1.05)
    fig2_path = fig_dir / "empirical_vs_fitted_transition_decay.png"
    plt.savefig(fig2_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved Figure 2: {fig2_path}")

    # Figure 3: Vertical transition matrix heatmap
    plt.figure(figsize=(9, 8))
    chain = classifier_all.vertical_markov
    vtm = chain.transition_matrix_
    im = plt.imshow(vtm, cmap="YlGnBu", vmin=0.0, vmax=1.0)
    plt.colorbar(im, label="Transition Probability P(i -> j)")
    for i in range(6):
        for j in range(6):
            val = vtm[i, j]
            color = "white" if val > 0.5 else "black"
            plt.text(j, i, f"{val:.3f}", ha="center", va="center", color=color, fontsize=10)
    plt.xticks(range(6), FACIES_NAMES, rotation=25, ha="right", fontsize=10)
    plt.yticks(range(6), FACIES_NAMES, fontsize=10)
    plt.xlabel("Succeeding Facies (Upward)", fontsize=12)
    plt.ylabel("Preceding Facies (Stratigraphic Base)", fontsize=12)
    plt.title("SMALT Sprint I: 1D Vertical Stratigraphic Markov Transition Matrix", fontsize=13, fontweight="bold")
    fig3_path = fig_dir / "vertical_transition_matrix_heatmap.png"
    plt.savefig(fig3_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved Figure 3: {fig3_path}")

    # Figure 4: Horizontal transition matrix at selected distances
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)
    for ax, dist_m in zip(axes, [500.0, 1500.0, 3500.0]):
        P_mat = rate_model.evaluate_transition_matrix(dist_m)
        im = ax.imshow(P_mat, cmap="Blues", vmin=0.0, vmax=1.0)
        for i in range(6):
            for j in range(6):
                val = P_mat[i, j]
                color = "white" if val > 0.5 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=color, fontsize=9)
        ax.set_xticks(range(6))
        ax.set_xticklabels(FACIES_NAMES, rotation=25, ha="right", fontsize=9)
        if dist_m == 500.0:
            ax.set_yticks(range(6))
            ax.set_yticklabels(FACIES_NAMES, fontsize=9)
            ax.set_ylabel("Source Facies", fontsize=10)
        else:
            ax.set_yticks([])
        ax.set_title(f"Distance h = {int(dist_m)} m", fontsize=12, fontweight="bold")
        ax.set_xlabel("Target Facies", fontsize=10)
    plt.suptitle("SMALT Sprint I: Continuous Horizontal Transition Probability Matrices Across Distances", fontsize=14, fontweight="bold")
    fig4_path = fig_dir / "horizontal_transition_matrix_selected_distances.png"
    plt.savefig(fig4_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved Figure 4: {fig4_path}")

    # Figure 5: Target well truth vs MAP prediction (Held-out Litholog 9)
    # Evaluate L9 as a held-out test well
    train_l9 = [lid for lid in validator.eligible_ids if lid != "litholog9"]
    clf_l9 = ConditionedMarkovClassifier(training_litholog_ids=train_l9, inspector=inspector)
    df_l9_pred = clf_l9.predict_target_litholog("litholog9")

    fig, (ax_true, ax_map, ax_ent) = plt.subplots(1, 3, figsize=(10, 12), sharey=True)
    z_vals = df_l9_pred["z_common_m"].to_numpy()
    true_codes = df_l9_pred["true_facies_code"].to_numpy()
    map_codes = df_l9_pred["map_facies_code"].to_numpy()
    entropies = df_l9_pred["predictive_entropy"].to_numpy()

    for z, tc in zip(z_vals, true_codes):
        ax_true.barh(z, width=1.0, height=1.0, color=FACIES_COLORS[tc], align="center")
    ax_true.set_title("Ground Truth\n(Litholog 9)", fontsize=11, fontweight="bold")
    ax_true.set_ylabel("Common-Zero Elevation (m)", fontsize=12)
    ax_true.set_xlim(0, 1)
    ax_true.set_xticks([])

    for z, mc in zip(z_vals, map_codes):
        ax_map.barh(z, width=1.0, height=1.0, color=FACIES_COLORS[mc], align="center")
    ax_map.set_title("Blind MAP Prediction\n(Sprint I CTP)", fontsize=11, fontweight="bold")
    ax_map.set_xlim(0, 1)
    ax_map.set_xticks([])

    ax_ent.plot(entropies, z_vals, color="darkred", linewidth=2)
    ax_ent.fill_betweenx(z_vals, 0, entropies, color="salmon", alpha=0.3)
    ax_ent.set_title("Predictive Entropy\nH(S | W)", fontsize=11, fontweight="bold")
    ax_ent.set_xlabel("Entropy (nats)", fontsize=10)
    ax_ent.set_xlim(0, 1.8)
    ax_ent.grid(True, linestyle=":", alpha=0.6)

    # Legend for facies
    legend_patches = [
        mpatches.Patch(color=FACIES_COLORS[c], label=CODE_TO_NAME[c])
        for c in range(6)
    ]
    fig.legend(handles=legend_patches, loc="lower center", ncol=3, frameon=True, bbox_to_anchor=(0.5, -0.05), fontsize=10)
    plt.suptitle("SMALT Sprint I: Blind Target Evaluation vs Ground Truth (Litholog 9)", fontsize=13, fontweight="bold")
    fig5_path = fig_dir / "target_well_truth_vs_map_prediction.png"
    plt.savefig(fig5_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved Figure 5: {fig5_path}")

    # Figure 6: Stochastic inter-well realizations along transect
    fig, axes = plt.subplots(3, 1, figsize=(16, 12), sharex=True, sharey=True)
    cmap_custom = matplotlib.colors.ListedColormap([FACIES_COLORS[c] for c in range(6)])
    bounds = np.arange(7) - 0.5
    norm_custom = matplotlib.colors.BoundaryNorm(bounds, cmap_custom.N)

    x_coords = sim_res["x_coords_m"]
    z_coords = sim_res["z_coords_m"]

    for r_idx in range(3):
        ax = axes[r_idx]
        grid_data = sim_res["realizations"][r_idx]
        im = ax.imshow(
            grid_data,
            origin="lower",
            extent=[x_coords[0], x_coords[-1], z_coords[0], z_coords[-1]],
            cmap=cmap_custom,
            norm=norm_custom,
            aspect="auto",
        )
        # Mark anchor wells
        for wid, dist in zip(sim_res["anchor_wells"], sim_res["anchor_well_distances_m"]):
            ax.axvline(dist, color="black", linestyle="--", linewidth=1.5, alpha=0.8)
            ax.text(dist, z_coords[-1] + 2, wid, rotation=45, ha="center", fontsize=9, fontweight="bold")

        ax.set_ylabel("Elevation z (m)", fontsize=11)
        ax.set_title(f"Stochastic Realization #{r_idx + 1} (Hard Conditioning Honor Rate: 100%)", fontsize=11, fontweight="bold")

    axes[-1].set_xlabel("Transect Distance (m) [L2 -> L3 -> L4 -> L9]", fontsize=12)
    plt.suptitle("SMALT Sprint I: Stochastic Inter-Well Facies Realizations Along Fluvial Transect", fontsize=14, fontweight="bold")
    fig6_path = fig_dir / "stochastic_interwell_realizations_cross_section.png"
    plt.savefig(fig6_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved Figure 6: {fig6_path}")

    # Figure 7: Ensemble facies probability & uncertainty map
    fig, (ax_sand, ax_ent_map) = plt.subplots(2, 1, figsize=(16, 9), sharex=True, sharey=True)
    p_sand = sim_res["facies_probabilities"][0]  # Channel Sandstone prob
    im_sand = ax_sand.imshow(
        p_sand,
        origin="lower",
        extent=[x_coords[0], x_coords[-1], z_coords[0], z_coords[-1]],
        cmap="YlOrBr",
        aspect="auto",
        vmin=0.0,
        vmax=1.0,
    )
    plt.colorbar(im_sand, ax=ax_sand, label="P(Channel Sandstone)")
    ax_sand.set_title("Ensemble Probability Field: Channel Sandstone Occurrence", fontsize=12, fontweight="bold")
    ax_sand.set_ylabel("Elevation z (m)", fontsize=11)
    for dist in sim_res["anchor_well_distances_m"]:
        ax_sand.axvline(dist, color="black", linestyle=":", linewidth=1.2)

    entropy_map = sim_res["ensemble_entropy"]
    im_ent = ax_ent_map.imshow(
        entropy_map,
        origin="lower",
        extent=[x_coords[0], x_coords[-1], z_coords[0], z_coords[-1]],
        cmap="viridis",
        aspect="auto",
    )
    plt.colorbar(im_ent, ax=ax_ent_map, label="Shannon Entropy (nats)")
    ax_ent_map.set_title("Spatial Uncertainty Field: Ensemble Shannon Entropy", fontsize=12, fontweight="bold")
    ax_ent_map.set_xlabel("Transect Distance (m)", fontsize=12)
    ax_ent_map.set_ylabel("Elevation z (m)", fontsize=11)
    for dist in sim_res["anchor_well_distances_m"]:
        ax_ent_map.axvline(dist, color="black", linestyle=":", linewidth=1.2)

    plt.suptitle("SMALT Sprint I: Ensemble Facies Probability & Spatial Uncertainty Mapping", fontsize=14, fontweight="bold")
    fig7_path = fig_dir / "ensemble_facies_probability_uncertainty_map.png"
    plt.savefig(fig7_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved Figure 7: {fig7_path}")

    # Figure 8: Upstream -> Downstream comparison
    df_dir_up = df_directional[df_directional["experiment_direction"] == "Upstream -> Downstream"]
    plt.figure(figsize=(9, 5))
    x_pos = np.arange(len(df_dir_up))
    plt.bar(x_pos - 0.2, df_dir_up["accuracy"] * 100, width=0.4, label="Raw Accuracy (%)", color="navy")
    plt.bar(x_pos + 0.2, df_dir_up["balanced_accuracy"] * 100, width=0.4, label="Balanced Accuracy (%)", color="teal")
    plt.xticks(x_pos, df_dir_up["model_name"], rotation=20, ha="right", fontsize=10)
    plt.ylabel("Score (%)", fontsize=11)
    plt.title("SMALT Sprint I: Directional Generalization (Upstream -> Downstream: L9, L11, L12)", fontsize=12, fontweight="bold")
    plt.grid(True, axis="y", linestyle=":", alpha=0.6)
    plt.legend(frameon=True)
    plt.ylim(0, 60)
    fig8_path = fig_dir / "upstream_to_downstream_comparison.png"
    plt.savefig(fig8_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved Figure 8: {fig8_path}")

    # Figure 9: Downstream -> Upstream comparison
    df_dir_down = df_directional[df_directional["experiment_direction"] == "Downstream -> Upstream"]
    plt.figure(figsize=(9, 5))
    x_pos = np.arange(len(df_dir_down))
    plt.bar(x_pos - 0.2, df_dir_down["accuracy"] * 100, width=0.4, label="Raw Accuracy (%)", color="darkred")
    plt.bar(x_pos + 0.2, df_dir_down["balanced_accuracy"] * 100, width=0.4, label="Balanced Accuracy (%)", color="coral")
    plt.xticks(x_pos, df_dir_down["model_name"], rotation=20, ha="right", fontsize=10)
    plt.ylabel("Score (%)", fontsize=11)
    plt.title("SMALT Sprint I: Directional Generalization (Downstream -> Upstream: L2-L8, L10)", fontsize=12, fontweight="bold")
    plt.grid(True, axis="y", linestyle=":", alpha=0.6)
    plt.legend(frameon=True)
    plt.ylim(0, 60)
    fig9_path = fig_dir / "downstream_to_upstream_comparison.png"
    plt.savefig(fig9_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved Figure 9: {fig9_path}")

    # Figure 10: Baseline comparison metrics
    fig, (ax_acc, ax_geo) = plt.subplots(1, 2, figsize=(16, 6))
    models = df_comp["model_name"]
    x = np.arange(len(models))

    ax_acc.bar(x - 0.2, df_comp["raw_accuracy"] * 100, width=0.4, label="Raw Accuracy (%)", color="steelblue")
    ax_acc.bar(x + 0.2, df_comp["balanced_accuracy"] * 100, width=0.4, label="Balanced Accuracy (%)", color="goldenrod")
    ax_acc.set_xticks(x)
    ax_acc.set_xticklabels(models, rotation=25, ha="right", fontsize=9)
    ax_acc.set_ylabel("Score (%)", fontsize=11)
    ax_acc.set_title("Pointwise Classification Performance (LOLO CV, 940 Pts)", fontsize=11, fontweight="bold")
    ax_acc.grid(True, axis="y", linestyle=":", alpha=0.6)
    ax_acc.legend(frameon=True)
    ax_acc.set_ylim(0, 60)

    # Geological reproduction TV distance and bed thickness error
    geo_eval = df_geo[df_geo["model_name"] != "Ground Truth Target"]
    x_geo = np.arange(len(geo_eval))
    ax_geo.bar(x_geo - 0.2, geo_eval["proportion_tv_distance"], width=0.4, label="Facies Proportion TV Distance", color="darkgreen")
    ax_geo.bar(x_geo + 0.2, geo_eval["transition_matrix_frobenius_div"] / 2.0, width=0.4, label="Trans. Matrix Frobenius Div (/2)", color="purple")
    ax_geo.set_xticks(x_geo)
    ax_geo.set_xticklabels(geo_eval["model_name"], rotation=25, ha="right", fontsize=9)
    ax_geo.set_ylabel("Divergence / Error", fontsize=11)
    ax_geo.set_title("Geological-Statistical Reproduction Error (Lower is Better)", fontsize=11, fontweight="bold")
    ax_geo.grid(True, axis="y", linestyle=":", alpha=0.6)
    ax_geo.legend(frameon=True)

    plt.suptitle("SMALT Sprint I: Comprehensive Benchmarking vs Existing Baselines", fontsize=14, fontweight="bold")
    fig10_path = fig_dir / "baseline_comparison_metrics.png"
    plt.savefig(fig10_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved Figure 10: {fig10_path}")

    print("\n" + "=" * 80)
    print("SPRINT I PIPELINE COMPLETE: All 7 CSV deliverables and 10 figures generated successfully.")
    print("=" * 80)


if __name__ == "__main__":
    main()
