"""
SMALT Spatial Markov Foundation & Transition Probability Modeling.

Extends the verified 1D vertical Markov succession framework to 2D/3D space
using the common-zero vertical datum alignment and geostatistical transition
probability theory (Carle & Fogg, 1996; Elfeki & Dekking, 2001).

Mathematical Framework:
1. Vertical Transition Probability:
   P_z(i, j; h_z) = P[S(x, y, z + h_z) = j | S(x, y, z) = i]
2. Horizontal Transition Probability:
   P_h(i, j; h) = P[S(x + h_x, y + h_y, z) = j | S(x, y, z) = i, ||h|| = h]
3. Transition Rate Matrix Formulation (Carle & Fogg, 1996):
   P(h) = expm(R * h)
   where R_ii = -1 / L_i (L_i = mean lateral length of facies i),
   R_ij = r_ij / L_i for j != i (r_ij = embedded transition probability),
   and sum_j R_ij = 0.
4. Asymptotic Limits:
   - h -> 0: P(h) -> I (Identity matrix; perfect self-correlation at well location)
   - h -> inf: P(h) -> 1 * p^T (Decays to stationary background proportions p)
5. Conditional Spatial Prediction:
   Given nearest conditioning well at distance h and elevation z with facies i,
   p_target(z) = e_i^T * P(h)
"""

from typing import Dict, Any, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import expm
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score

from smalt.descriptive.analyzer import LithologInspector, CANONICAL_FACIES_SCHEMA
from smalt.spatial.coordinates import load_source_coordinates
from smalt.spatial.datum import align_to_common_datum, DEFAULT_CONVENTION


# Default lateral lengths from Sahoo et al. (2016) aspect ratios (W/T ~ 35 for channel sand)
DEFAULT_LATERAL_FACIES_LENGTHS_M = {
    0: 203.0,  # Channel Sandstone (W/T ~ 35 * mean thickness 5.8 m)
    1: 30.0,   # Planar Sandstone (upper flow regime sheet sand, ~20-50 m)
    2: 75.0,   # Rippled Heterolithics (crevasse splay / levee margin, 10-130 m)
    3: 40.0,   # Carbonaceous Mudstone (waterlogged swamp margin, 20-60 m)
    4: 55.0,   # Coal (peat mire seam, 30-100 m)
    5: 360.0,  # Overbank Mudstone (widespread floodplain fines, 200-500 m)
}


class SpatialMarkovTransitionAnalyzer:
    """
    Extracts horizontal facies pairs across common-zero datum elevation slices,
    tallies empirical horizontal transition counts, and estimates continuous
    transition probability matrices P(h) via transition rate matrices.
    """

    def __init__(
        self,
        inspector: Optional[LithologInspector] = None,
        convention: str = DEFAULT_CONVENTION,
    ):
        self.inspector = inspector or LithologInspector()
        self.convention = convention
        self.coords_df = load_source_coordinates()
        self.coords_by_id = {
            r["litholog_id"]: (float(r["x_m"]), float(r["y_m"]))
            for _, r in self.coords_df.iterrows()
            if r["coordinates_available"] and r["litholog_id"] != "litholog1"
        }
        self.eligible_ids = sorted(list(self.coords_by_id.keys()))
        self._elevation_grid_cache: Optional[pd.DataFrame] = None

    def build_common_elevation_grid(self) -> pd.DataFrame:
        """
        Builds a unified 1m elevation slice grid across all eligible wells (L2-L12)
        aligned under the common-zero datum.

        Returns DataFrame indexed by (z_common_m, litholog_id) with columns:
        - litholog_id
        - depth_original_m
        - z_common_m
        - facies
        - facies_code
        - x_m
        - y_m
        """
        if self._elevation_grid_cache is not None:
            return self._elevation_grid_cache.copy()

        dfs = []
        for lid in self.eligible_ids:
            disc_df = self.inspector.discretize_litholog_1m(lid)
            aligned_df = align_to_common_datum(disc_df, litholog_id=lid, convention=self.convention)
            x_m, y_m = self.coords_by_id[lid]
            aligned_df["x_m"] = x_m
            aligned_df["y_m"] = y_m
            dfs.append(aligned_df)

        df_grid = pd.concat(dfs, ignore_index=True)
        self._elevation_grid_cache = df_grid
        return df_grid.copy()

    def extract_horizontal_facies_pairs(
        self,
        max_elevation_diff_m: float = 0.5,
    ) -> pd.DataFrame:
        """
        Extracts all horizontal well pairs at matching common-zero elevation slices.

        For each unique elevation z, finds all pairs of wells (A, B) active at z,
        recording facies_A, facies_B, horizontal lag distance h, and vector (dx, dy).
        """
        df_grid = self.build_common_elevation_grid()
        unique_z = np.sort(df_grid["z_common_m"].unique())

        pairs = []
        for z_val in unique_z:
            slice_df = df_grid[np.abs(df_grid["z_common_m"] - z_val) <= max_elevation_diff_m].reset_index(drop=True)
            n_wells = len(slice_df)
            if n_wells < 2:
                continue

            for i in range(n_wells):
                row_a = slice_df.iloc[i]
                for j in range(i + 1, n_wells):
                    row_b = slice_df.iloc[j]
                    dx = row_b["x_m"] - row_a["x_m"]
                    dy = row_b["y_m"] - row_a["y_m"]
                    dist = float(np.hypot(dx, dy))
                    theta_deg = float(np.degrees(np.arctan2(dy, dx))) % 360.0

                    # Record pair (symmetrized for omnidirectional tallying)
                    pairs.append({
                        "z_common_m": z_val,
                        "well_a": row_a["litholog_id"],
                        "well_b": row_b["litholog_id"],
                        "facies_code_a": int(row_a["facies_code"]),
                        "facies_code_b": int(row_b["facies_code"]),
                        "facies_name_a": row_a["facies"],
                        "facies_name_b": row_b["facies"],
                        "lag_distance_m": round(dist, 1),
                        "dx_m": round(dx, 1),
                        "dy_m": round(dy, 1),
                        "azimuth_deg": round(theta_deg, 1),
                        "is_same_facies": int(row_a["facies_code"]) == int(row_b["facies_code"]),
                    })

        df_pairs = pd.DataFrame(pairs)
        return df_pairs

    def tally_horizontal_transition_counts(
        self,
        min_lag_m: float = 0.0,
        max_lag_m: float = 1e6,
        alpha: float = 0.01,
    ) -> Dict[str, Any]:
        """
        Tallies empirical horizontal transition counts N_ij within a lag distance window [min_lag_m, max_lag_m).
        Computes row-stochastic transition probability matrix T_ij.
        """
        df_pairs = self.extract_horizontal_facies_pairs()
        subset = df_pairs[
            (df_pairs["lag_distance_m"] >= min_lag_m) & (df_pairs["lag_distance_m"] < max_lag_m)
        ]

        n_classes = len(CANONICAL_FACIES_SCHEMA)
        counts = np.zeros((n_classes, n_classes), dtype=float)

        # Symmetrized tally: pair (A, B) contributes A -> B and B -> A
        for _, r in subset.iterrows():
            ca = int(r["facies_code_a"])
            cb = int(r["facies_code_b"])
            counts[ca, cb] += 1.0
            counts[cb, ca] += 1.0

        # Apply Laplace smoothing and compute probabilities
        smoothed_counts = counts + alpha
        row_sums = smoothed_counts.sum(axis=1, keepdims=True)
        prob_matrix = smoothed_counts / row_sums

        # Marginal stationary facies proportions from pair counts
        total_counts = counts.sum()
        marginals = counts.sum(axis=1) / total_counts if total_counts > 0 else np.full(n_classes, 1.0 / n_classes)

        return {
            "min_lag_m": min_lag_m,
            "max_lag_m": max_lag_m,
            "pair_count": len(subset),
            "transition_count_matrix": counts,
            "transition_prob_matrix": np.round(prob_matrix, 4),
            "marginal_proportions": np.round(marginals, 4),
            "auto_transition_proportions": np.round(np.diag(prob_matrix), 4),
        }

    def compute_lag_binned_horizontal_matrices(
        self,
        bins: Optional[List[Tuple[float, float, str]]] = None,
    ) -> pd.DataFrame:
        """
        Computes empirical horizontal transition matrices across standard lag bins.
        """
        lag_bins = bins or [
            (400.0, 1000.0, "Short Lag (400 - 1000 m, Local Pairs)"),
            (1000.0, 2500.0, "Intermediate Lag (1000 - 2500 m)"),
            (2500.0, 5500.0, "Long Lag (2500 - 5500 m, Regional Span)"),
            (400.0, 5500.0, "Omnidirectional (All Pairs >= 400 m)"),
        ]

        code_to_meta = {meta["code"]: meta for meta in CANONICAL_FACIES_SCHEMA.values()}
        records = []
        for min_d, max_d, label in lag_bins:
            res = self.tally_horizontal_transition_counts(min_d, max_d)
            p_mat = res["transition_prob_matrix"]
            for i in range(len(CANONICAL_FACIES_SCHEMA)):
                for j in range(len(CANONICAL_FACIES_SCHEMA)):
                    records.append({
                        "lag_bin_label": label,
                        "min_lag_m": min_d,
                        "max_lag_m": max_d,
                        "pair_observations": res["pair_count"],
                        "from_facies_code": i,
                        "from_facies_name": code_to_meta[i]["canonical_name"],
                        "to_facies_code": j,
                        "to_facies_name": code_to_meta[j]["canonical_name"],
                        "transition_count": res["transition_count_matrix"][i, j],
                        "transition_probability": p_mat[i, j],
                    })

        return pd.DataFrame(records)

    def build_theoretical_horizontal_rate_matrix(
        self,
        lateral_lengths: Optional[Dict[int, float]] = None,
        embedded_proportions: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Constructs the Carle & Fogg (1996) continuous horizontal transition rate matrix R_h.

        R_h(i, i) = -1 / L_i
        R_h(i, j) = r_ij / L_i (for j != i)
        where sum_j R_h(i, j) = 0.
        """
        n_classes = len(CANONICAL_FACIES_SCHEMA)
        lengths = lateral_lengths or DEFAULT_LATERAL_FACIES_LENGTHS_M

        # If embedded proportions r_ij not provided, use equal off-diagonal branching
        if embedded_proportions is None:
            r = np.ones((n_classes, n_classes)) / (n_classes - 1)
            np.fill_diagonal(r, 0.0)
        else:
            r = embedded_proportions.copy()
            np.fill_diagonal(r, 0.0)
            row_sums = r.sum(axis=1, keepdims=True)
            r = np.where(row_sums > 0, r / row_sums, 1.0 / (n_classes - 1))
            np.fill_diagonal(r, 0.0)

        R = np.zeros((n_classes, n_classes), dtype=float)
        for i in range(n_classes):
            L_i = float(lengths.get(i, 100.0))
            R[i, i] = -1.0 / L_i
            for j in range(n_classes):
                if i != j:
                    R[i, j] = r[i, j] / L_i

        # Assert mathematical row-sum zero condition
        np.testing.assert_allclose(R.sum(axis=1), np.zeros(n_classes), atol=1e-12)
        return R

    def evaluate_transition_probability_at_lag(
        self,
        h_m: float,
        R_matrix: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Computes continuous horizontal transition probability matrix P(h) = expm(R * h).
        Guarantees row-stochasticity.
        """
        R = R_matrix if R_matrix is not None else self.build_theoretical_horizontal_rate_matrix()
        P_h = expm(R * h_m)
        # Ensure numerical non-negativity and row-stochasticity
        P_h = np.clip(P_h, 0.0, 1.0)
        row_sums = P_h.sum(axis=1, keepdims=True)
        P_h = P_h / row_sums
        return P_h


# -----------------------------------------------------------------------------
# Minimal Conditional Spatial Markov Classifier Prototype
# -----------------------------------------------------------------------------

class SpatialMarkovPredictor:
    """
    Minimal conditional spatial Markov classifier.
    Predicts facies probability vector at an unobserved target location (x, y, z)
    by conditioning on the nearest observed well at distance h using the
    continuous horizontal transition probability matrix P(h) = expm(R_h * h).
    """

    def __init__(
        self,
        analyzer: Optional[SpatialMarkovTransitionAnalyzer] = None,
        lateral_lengths: Optional[Dict[int, float]] = None,
        convention: str = DEFAULT_CONVENTION,
    ):
        self.analyzer = analyzer or SpatialMarkovTransitionAnalyzer(convention=convention)
        self.convention = convention
        self.lateral_lengths = lateral_lengths or DEFAULT_LATERAL_FACIES_LENGTHS_M
        self.R_matrix = self.analyzer.build_theoretical_horizontal_rate_matrix(
            lateral_lengths=self.lateral_lengths
        )

    def run_spatial_markov_lolo(
        self,
        eligible_log_ids: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Executes Leave-One-Litholog-Out (LOLO) cross-validation evaluating the
        minimal spatial Markov transition model against naive baselines on the
        common-zero datum.
        """
        eligible_ids = eligible_log_ids or self.analyzer.eligible_ids
        df_grid = self.analyzer.build_common_elevation_grid()

        fold_records = []
        all_y_true = []
        all_y_pred_markov = []
        all_y_pred_near_well = []
        all_y_pred_prior = []

        coords = self.analyzer.coords_by_id

        for target_id in eligible_ids:
            train_ids = [lid for lid in eligible_ids if lid != target_id]
            df_target = df_grid[df_grid["litholog_id"] == target_id].copy().reset_index(drop=True)
            df_train = df_grid[df_grid["litholog_id"].isin(train_ids)].copy().reset_index(drop=True)

            target_x, target_y = coords[target_id]

            # Find nearest training well
            dists = [
                (float(np.hypot(coords[lid][0] - target_x, coords[lid][1] - target_y)), lid)
                for lid in train_ids
            ]
            min_dist, nearest_well_id = min(dists, key=lambda x: x[0])
            df_near = df_train[df_train["litholog_id"] == nearest_well_id].copy()

            # Precompute transition probability matrix P(h) at inter-well distance min_dist
            P_h = self.analyzer.evaluate_transition_probability_at_lag(min_dist, self.R_matrix)

            # Training set majority class prior
            classes, counts = np.unique(df_train["facies_code"].to_numpy(), return_counts=True)
            prior_majority = classes[np.argmax(counts)]

            y_true_fold = df_target["facies_code"].to_numpy()
            y_pred_markov = []
            y_pred_near = []
            y_pred_prior = []

            near_z = df_near["z_common_m"].to_numpy()
            near_codes = df_near["facies_code"].to_numpy()

            for z_val in df_target["z_common_m"].to_numpy():
                # Find matching or closest elevation in nearest well
                idx_closest = int(np.argmin(np.abs(near_z - z_val)))
                near_code = int(near_codes[idx_closest])

                # Spatial Markov prediction: e_i^T * P(h)
                prob_vec = P_h[near_code, :]
                markov_pred = int(np.argmax(prob_vec))

                y_pred_markov.append(markov_pred)
                y_pred_near.append(near_code)
                y_pred_prior.append(prior_majority)

            y_pred_markov = np.array(y_pred_markov)
            y_pred_near = np.array(y_pred_near)
            y_pred_prior = np.array(y_pred_prior)

            # Metrics
            acc_markov = float(accuracy_score(y_true_fold, y_pred_markov))
            bal_markov = float(balanced_accuracy_score(y_true_fold, y_pred_markov))
            f1_markov = float(f1_score(y_true_fold, y_pred_markov, average="macro", zero_division=0))

            acc_near = float(accuracy_score(y_true_fold, y_pred_near))
            bal_near = float(balanced_accuracy_score(y_true_fold, y_pred_near))
            f1_near = float(f1_score(y_true_fold, y_pred_near, average="macro", zero_division=0))

            acc_prior = float(accuracy_score(y_true_fold, y_pred_prior))
            bal_prior = float(balanced_accuracy_score(y_true_fold, y_pred_prior))
            f1_prior = float(f1_score(y_true_fold, y_pred_prior, average="macro", zero_division=0))

            fold_records.append({
                "target_litholog_id": target_id,
                "nearest_well_id": nearest_well_id,
                "nearest_distance_m": round(min_dist, 1),
                "test_points_count": len(y_true_fold),
                "spatial_markov_accuracy": round(acc_markov, 4),
                "spatial_markov_balanced_acc": round(bal_markov, 4),
                "spatial_markov_macro_f1": round(f1_markov, 4),
                "near_well_accuracy": round(acc_near, 4),
                "near_well_balanced_acc": round(bal_near, 4),
                "near_well_macro_f1": round(f1_near, 4),
                "prior_accuracy": round(acc_prior, 4),
                "prior_balanced_acc": round(bal_prior, 4),
                "prior_macro_f1": round(f1_prior, 4),
            })

            all_y_true.extend(y_true_fold)
            all_y_pred_markov.extend(y_pred_markov)
            all_y_pred_near_well.extend(y_pred_near)
            all_y_pred_prior.extend(y_pred_prior)

        yt = np.array(all_y_true)
        yp_m = np.array(all_y_pred_markov)
        yp_n = np.array(all_y_pred_near_well)
        yp_p = np.array(all_y_pred_prior)

        summary = {
            "total_test_points": len(yt),
            "pooled_markov_accuracy": round(float(accuracy_score(yt, yp_m)), 4),
            "pooled_markov_balanced_acc": round(float(balanced_accuracy_score(yt, yp_m)), 4),
            "pooled_markov_macro_f1": round(float(f1_score(yt, yp_m, average="macro", zero_division=0)), 4),
            "pooled_near_well_accuracy": round(float(accuracy_score(yt, yp_n)), 4),
            "pooled_near_well_balanced_acc": round(float(balanced_accuracy_score(yt, yp_n)), 4),
            "pooled_near_well_macro_f1": round(float(f1_score(yt, yp_n, average="macro", zero_division=0)), 4),
            "pooled_prior_accuracy": round(float(accuracy_score(yt, yp_p)), 4),
            "pooled_prior_balanced_acc": round(float(balanced_accuracy_score(yt, yp_p)), 4),
            "pooled_prior_macro_f1": round(float(f1_score(yt, yp_p, average="macro", zero_division=0)), 4),
        }

        return {
            "per_fold_table": pd.DataFrame(fold_records),
            "aggregate_summary": summary,
        }
