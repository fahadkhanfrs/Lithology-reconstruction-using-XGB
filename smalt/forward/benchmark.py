"""
SMALT Synthetic Borehole Reconstruction Benchmark (Sprint K Stage 3).

Evaluates baseline reconstruction algorithms on sparse synthetic boreholes:
1. Nearest-Well Profile Baseline.
2. Vertical Stratigraphic Prior Baseline.
3. Spatial 2D Anisotropic KNN Baseline.
4. Spatial Markov Decay Baseline.

Guarantees:
- Strict leakage prevention: models receive only sparse borehole data.
- Withheld subsurface ground truth is evaluated through a dedicated harness.
- Realization-level metric aggregation (Raw Acc, Balanced Acc, Macro-F1, per-facies metrics).
- Spacing sensitivity evaluated across 400 m, 600 m, and 1000 m.
- Evaluation split between interpolation domain and extrapolation margins.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    precision_recall_fscore_support,
    confusion_matrix,
)
from sklearn.neighbors import KNeighborsClassifier
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from smalt.descriptive.analyzer import CANONICAL_FACIES_SCHEMA
from smalt.forward.config import get_preset_config
from smalt.forward.simulator import (
    FluvialForwardSimulator2D,
    SyntheticRealization2D,
    CANONICAL_CODE_TO_NAME,
)
from smalt.forward.sampling import (
    SyntheticBoreholeSampler,
    SyntheticBoreholeDataset,
)

FACIES_NAMES = ["sand", "p_sand", "ripples", "carbon_mud", "coal", "mud"]


# =============================================================================
# Baseline Reconstruction Models
# =============================================================================
class NearestWellModel:
    """
    Nearest-Well Profile Baseline: assigns the facies profile of the horizontally
    nearest sampled borehole at each elevation z.
    """

    def __init__(self) -> None:
        self.dataset: Optional[SyntheticBoreholeDataset] = None

    def fit(self, dataset: SyntheticBoreholeDataset) -> "NearestWellModel":
        self.dataset = dataset
        return self

    def predict(self, realization: SyntheticRealization2D) -> np.ndarray:
        if self.dataset is None:
            raise RuntimeError("Model must be fitted before predict.")

        nz, nx = realization.shape
        pred_grid = np.zeros((nz, nx), dtype=np.int32)
        x_wells = np.array(self.dataset.x_locations_m)

        for j, x_val in enumerate(realization.x_coords_m):
            # Find nearest borehole
            nearest_idx = int(np.argmin(np.abs(x_wells - x_val)))
            nearest_id = self.dataset.borehole_ids[nearest_idx]
            bh_df = self.dataset.get_borehole(nearest_id)
            # In ascending z order
            pred_grid[:, j] = bh_df["facies_code"].to_numpy()

        return pred_grid


class VerticalPriorModel:
    """
    Vertical Prior Baseline: computes the empirical modal facies at each height z
    across all sampled boreholes (modeling vertical Waltherian stratigraphy without
    lateral continuity).
    """

    def __init__(self) -> None:
        self.vertical_modes: Optional[np.ndarray] = None
        self.global_prior_code: int = 5

    def fit(self, dataset: SyntheticBoreholeDataset) -> "VerticalPriorModel":
        df = dataset.data
        # Group by grid row index (z elevation)
        row_modes = []
        max_row = df["grid_row_idx"].max()

        # Global mode fallback
        self.global_prior_code = int(df["facies_code"].mode().iloc[0])

        for r in range(max_row + 1):
            sub = df[df["grid_row_idx"] == r]
            if not sub.empty:
                mode_code = int(sub["facies_code"].mode().iloc[0])
                row_modes.append(mode_code)
            else:
                row_modes.append(self.global_prior_code)

        self.vertical_modes = np.array(row_modes, dtype=np.int32)
        return self

    def predict(self, realization: SyntheticRealization2D) -> np.ndarray:
        if self.vertical_modes is None:
            raise RuntimeError("Model must be fitted before predict.")

        nz, nx = realization.shape
        pred_grid = np.zeros((nz, nx), dtype=np.int32)
        for j in range(nx):
            pred_grid[:, j] = self.vertical_modes[:nz]
        return pred_grid


class SpatialKNNModel:
    """
    Spatial 2D Anisotropic KNN Baseline: k=5 nearest neighbors on (x, z)
    with geological vertical exaggeration (z scaled by aspect ratio).
    """

    def __init__(self, n_neighbors: int = 5, vertical_anisotropy: float = 20.0) -> None:
        self.n_neighbors = n_neighbors
        self.anisotropy = vertical_anisotropy
        self.clf = KNeighborsClassifier(n_neighbors=n_neighbors, weights="distance")

    def fit(self, dataset: SyntheticBoreholeDataset) -> "SpatialKNNModel":
        df = dataset.data
        X_train = np.column_stack([
            df["synthetic_x_m"].to_numpy(),
            df["strat_height_m"].to_numpy() * self.anisotropy,
        ])
        y_train = df["facies_code"].to_numpy()
        self.clf.fit(X_train, y_train)
        return self

    def predict(self, realization: SyntheticRealization2D) -> np.ndarray:
        nz, nx = realization.shape
        # Create full coordinate grid
        xx, zz = np.meshgrid(realization.x_coords_m, realization.z_coords_m)
        X_eval = np.column_stack([
            xx.ravel(),
            zz.ravel() * self.anisotropy,
        ])
        y_pred = self.clf.predict(X_eval)
        return y_pred.reshape((nz, nx)).astype(np.int32)


class SpatialMarkovModel:
    """
    Spatial Markov Decay Baseline: combines nearest borehole vertical column
    with horizontal Markov decay toward sampled stationary prior proportions.
    """

    def __init__(self, decay_length_m: float = 200.0) -> None:
        self.decay_length_m = decay_length_m
        self.dataset: Optional[SyntheticBoreholeDataset] = None
        self.prior_probs: np.ndarray = np.zeros(6)

    def fit(self, dataset: SyntheticBoreholeDataset) -> "SpatialMarkovModel":
        self.dataset = dataset
        # Compute sampled prior proportions across 6 classes
        counts = np.bincount(dataset.data["facies_code"].to_numpy(), minlength=6)
        self.prior_probs = counts / np.sum(counts)
        return self

    def predict(self, realization: SyntheticRealization2D) -> np.ndarray:
        if self.dataset is None:
            raise RuntimeError("Model must be fitted before predict.")

        nz, nx = realization.shape
        pred_grid = np.zeros((nz, nx), dtype=np.int32)
        x_wells = np.array(self.dataset.x_locations_m)
        L = self.decay_length_m

        for j, x_val in enumerate(realization.x_coords_m):
            nearest_idx = int(np.argmin(np.abs(x_wells - x_val)))
            dist = float(np.abs(x_wells[nearest_idx] - x_val))
            nearest_id = self.dataset.borehole_ids[nearest_idx]
            col_obs = self.dataset.get_borehole(nearest_id)["facies_code"].to_numpy()

            # Decay factor rho(h) = exp(-h / L)
            rho = np.exp(-dist / L)

            # If rho is high (close to well), follow nearest well.
            # If rho decays towards 0, predict training modal prior.
            if rho >= 0.5:
                pred_grid[:, j] = col_obs
            else:
                modal_prior = int(np.argmax(self.prior_probs))
                pred_grid[:, j] = modal_prior

        return pred_grid


# =============================================================================
# Benchmark Evaluation Pipeline
# =============================================================================
class SyntheticReconstructionBenchmark:
    """
    Orchestrates the sparse-borehole reconstruction benchmark across baselines,
    spacings, and independent evaluation realizations.
    """

    def __init__(self, output_dir: str = "sprints/audit_sprint_k") -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def evaluate_realization(
        self,
        realization: SyntheticRealization2D,
        spacing_m: float = 400.0,
        models: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """
        Evaluates models on a single realization with sparse borehole sampling.
        """
        sampler = SyntheticBoreholeSampler(realization)
        # Sample with lateral offset of half spacing to ensure boreholes are interior
        dataset = sampler.sample_regular(spacing_m=spacing_m, offset_m=spacing_m / 2.0)

        eval_models = models or {
            "Nearest_Well": NearestWellModel(),
            "Vertical_Prior": VerticalPriorModel(),
            "Spatial_KNN": SpatialKNNModel(n_neighbors=5, vertical_anisotropy=20.0),
            "Spatial_Markov": SpatialMarkovModel(decay_length_m=200.0),
        }

        # Define evaluation mask: strictly withhold columns containing sampled boreholes
        nz, nx = realization.shape
        dx = realization.config.grid.dx_m
        sampled_col_indices = set(int(np.clip(np.round(x / dx), 0, nx - 1)) for x in dataset.x_locations_m)

        # Domain masks
        all_mask = np.ones((nz, nx), dtype=bool)
        for c in sampled_col_indices:
            all_mask[:, c] = False  # Exclude observed boreholes from test evaluation

        # Interpolation domain: between min and max borehole locations
        min_bx = min(dataset.x_locations_m)
        max_bx = max(dataset.x_locations_m)
        interp_mask = all_mask.copy()
        for j, x_val in enumerate(realization.x_coords_m):
            if not (min_bx <= x_val <= max_bx):
                interp_mask[:, j] = False

        # Extrapolation domain: margins beyond outer boreholes
        extrap_mask = all_mask.copy()
        for j, x_val in enumerate(realization.x_coords_m):
            if min_bx <= x_val <= max_bx:
                extrap_mask[:, j] = False

        y_true_all = realization.grid[all_mask]
        y_true_interp = realization.grid[interp_mask] if np.any(interp_mask) else np.array([])
        y_true_extrap = realization.grid[extrap_mask] if np.any(extrap_mask) else np.array([])

        records = []

        for name, model in eval_models.items():
            model.fit(dataset)
            pred_grid = model.predict(realization)

            # Predictions on withheld cells
            y_pred_all = pred_grid[all_mask]

            # Overall metrics
            raw_acc = float(accuracy_score(y_true_all, y_pred_all))
            bal_acc = float(balanced_accuracy_score(y_true_all, y_pred_all))
            macro_f1 = float(f1_score(y_true_all, y_pred_all, average="macro", zero_division=0))

            # Domain-specific accuracy
            if len(y_true_interp) > 0:
                y_pred_interp = pred_grid[interp_mask]
                interp_acc = float(accuracy_score(y_true_interp, y_pred_interp))
            else:
                interp_acc = raw_acc

            if len(y_true_extrap) > 0:
                y_pred_extrap = pred_grid[extrap_mask]
                extrap_acc = float(accuracy_score(y_true_extrap, y_pred_extrap))
            else:
                extrap_acc = raw_acc

            # Per-facies metrics
            p, r, f1, _ = precision_recall_fscore_support(
                y_true_all, y_pred_all, labels=range(6), zero_division=0
            )

            rec = {
                "seed": realization.seed,
                "spacing_m": spacing_m,
                "num_boreholes": dataset.num_boreholes,
                "model": name,
                "raw_accuracy": round(raw_acc, 4),
                "balanced_accuracy": round(bal_acc, 4),
                "macro_f1": round(macro_f1, 4),
                "interp_accuracy": round(interp_acc, 4),
                "extrap_accuracy": round(extrap_acc, 4),
                "num_eval_cells": len(y_true_all),
            }

            for code, fname in CANONICAL_CODE_TO_NAME.items():
                rec[f"{fname}_precision"] = round(float(p[code]), 4)
                rec[f"{fname}_recall"] = round(float(r[code]), 4)
                rec[f"{fname}_f1"] = round(float(f1[code]), 4)

            records.append(rec)

        return records

    def run_benchmark_campaign(
        self,
        spacings: List[float] = [400.0, 600.0, 1000.0],
        num_eval_realizations: int = 10,
        base_seed: int = 8000,
    ) -> pd.DataFrame:
        """
        Executes benchmark across all spacings and evaluation realizations.
        """
        sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
        all_records = []

        for spacing in spacings:
            for s_idx in range(num_eval_realizations):
                seed = base_seed + int(spacing) * 10 + s_idx
                realization = sim.generate(seed=seed)
                run_records = self.evaluate_realization(realization, spacing_m=spacing)
                all_records.extend(run_records)

        df_bench = pd.DataFrame(all_records)
        csv_path = self.output_dir / "synthetic_reconstruction_benchmark_metrics.csv"
        df_bench.to_csv(csv_path, index=False)
        return df_bench

    def summarize_benchmark(self, df_bench: pd.DataFrame) -> pd.DataFrame:
        """
        Computes mean and std across independent realizations for each model and spacing.
        """
        summary_rows = []
        metrics = ["raw_accuracy", "balanced_accuracy", "macro_f1", "interp_accuracy", "extrap_accuracy"]

        for (spacing, model), group in df_bench.groupby(["spacing_m", "model"]):
            row = {
                "spacing_m": spacing,
                "model": model,
                "num_realizations": len(group),
                "num_boreholes": group["num_boreholes"].iloc[0],
            }
            for m in metrics:
                vals = group[m].to_numpy()
                row[f"{m}_mean"] = round(float(np.mean(vals)), 4)
                row[f"{m}_std"] = round(float(np.std(vals, ddof=1)), 4)
            summary_rows.append(row)

        df_summary = pd.DataFrame(summary_rows)
        csv_path = self.output_dir / "synthetic_reconstruction_summary.csv"
        df_summary.to_csv(csv_path, index=False)
        return df_summary

    def plot_benchmark_figures(self, df_summary: pd.DataFrame) -> None:
        """
        Plots accuracy vs spacing comparison and example reconstruction cross-sections.
        """
        # Figure 1: Accuracy vs Borehole Spacing
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        models = sorted(df_summary["model"].unique())
        spacings = sorted(df_summary["spacing_m"].unique())

        # Raw Accuracy
        for m in models:
            sub = df_summary[df_summary["model"] == m].sort_values("spacing_m")
            axes[0].errorbar(
                sub["spacing_m"],
                sub["raw_accuracy_mean"] * 100,
                yerr=sub["raw_accuracy_std"] * 100,
                marker="o",
                capsize=4,
                label=m,
                linewidth=1.8,
            )

        axes[0].set_title("Raw Accuracy vs Borehole Spacing (Withheld Cells)", fontsize=11, fontweight="bold")
        axes[0].set_xlabel("Borehole Spacing (m)", fontsize=10)
        axes[0].set_ylabel("Raw Accuracy (%)", fontsize=10)
        axes[0].grid(True, linestyle="--", alpha=0.5)
        axes[0].legend(loc="upper right", fontsize=9)
        axes[0].set_xticks(spacings)

        # Balanced Accuracy
        for m in models:
            sub = df_summary[df_summary["model"] == m].sort_values("spacing_m")
            axes[1].errorbar(
                sub["spacing_m"],
                sub["balanced_accuracy_mean"] * 100,
                yerr=sub["balanced_accuracy_std"] * 100,
                marker="s",
                capsize=4,
                label=m,
                linewidth=1.8,
            )

        axes[1].set_title("Balanced Accuracy vs Borehole Spacing", fontsize=11, fontweight="bold")
        axes[1].set_xlabel("Borehole Spacing (m)", fontsize=10)
        axes[1].set_ylabel("Balanced Accuracy (%)", fontsize=10)
        axes[1].grid(True, linestyle="--", alpha=0.5)
        axes[1].legend(loc="upper right", fontsize=9)
        axes[1].set_xticks(spacings)

        plt.suptitle("SMALT Sprint K Synthetic Reconstruction Benchmark", fontsize=13, fontweight="bold")
        plt.tight_layout()
        plt.savefig(self.output_dir / "fig_benchmark_accuracy_vs_spacing.png", dpi=200)
        plt.close(fig)

        # Figure 2: Example Cross-Sections (Ground Truth vs Models)
        sim = FluvialForwardSimulator2D(get_preset_config("upstream"))
        eval_real = sim.generate(seed=8000)
        sampler = SyntheticBoreholeSampler(eval_real)
        ds = sampler.sample_regular(spacing_m=400.0, offset_m=200.0)

        nw_model = NearestWellModel().fit(ds)
        knn_model = SpatialKNNModel().fit(ds)
        sm_model = SpatialMarkovModel().fit(ds)

        grid_truth = eval_real.grid
        grid_nw = nw_model.predict(eval_real)
        grid_knn = knn_model.predict(eval_real)
        grid_sm = sm_model.predict(eval_real)

        fig, axes = plt.subplots(4, 1, figsize=(14, 12), sharex=True, sharey=True)
        colors = [CANONICAL_FACIES_SCHEMA[CANONICAL_CODE_TO_NAME[c]]["color_hex"] for c in range(6)]
        cmap = matplotlib.colors.ListedColormap(colors)
        bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        extent = [0.0, eval_real.config.grid.x_extent_m, 0.0, eval_real.config.grid.z_extent_m]

        plots = [
            (grid_truth, "True Synthetic Geology (Ground Truth) with 400m Boreholes"),
            (grid_nw, "Nearest-Well Profile Reconstruction"),
            (grid_knn, "Spatial 2D Anisotropic KNN Reconstruction (k=5)"),
            (grid_sm, "Spatial Markov Continuous Decay Reconstruction"),
        ]

        for ax, (grid, title) in zip(axes, plots):
            ax.imshow(grid, origin="lower", extent=extent, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
            for x_loc in ds.x_locations_m:
                ax.axvline(x=x_loc, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
            ax.set_ylabel("Height Z (m)", fontsize=9)
            ax.set_title(title, fontsize=10, fontweight="bold")

        axes[-1].set_xlabel("Lateral Distance X (m)", fontsize=10)
        plt.suptitle("Synthetic Facies Reconstruction Cross-Sections Comparison (400m Spacing)", fontsize=12, fontweight="bold")
        plt.tight_layout()
        plt.savefig(self.output_dir / "fig_benchmark_reconstruction_cross_sections.png", dpi=200)
        plt.close(fig)


if __name__ == "__main__":
    bench = SyntheticReconstructionBenchmark()
    print("Running Synthetic Reconstruction Benchmark Campaign...")
    df_bench = bench.run_benchmark_campaign(spacings=[400.0, 600.0, 1000.0], num_eval_realizations=10)
    print("Benchmark campaign completed. Summarizing metrics...")
    df_summary = bench.summarize_benchmark(df_bench)
    print("Rendering benchmark figures...")
    bench.plot_benchmark_figures(df_summary)
    print("Stage 3 Benchmark Deliverables Generated Successfully.")
    print(df_summary.to_string())
