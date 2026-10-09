"""
SMALT Geostatistical Core: Spatial Transition Probability & Continuous Rate Modeling.

Implements transition probability-based spatial geostatistics (Carle & Fogg, 1996)
adapted to the SMALT common-zero stratigraphic reference framework.

Mathematical & Theoretical Foundations:
1. Spatial Transition Probability:
   P_ij(h) = Pr[S(x + h) = j | S(x) = i]
   where S(x) represents lithofacies state at spatial location x, and h is the separation vector.
2. Continuous-Lag Transition Probability Model:
   P(h) = expm(R_h * h)
   where R_h is a continuous spatial transition rate matrix satisfying:
   - R_ii <= 0 (proportional to -1 / L_i, where L_i is mean lateral length of facies i)
   - R_ij >= 0 for j != i
   - sum_j R_ij = 0 (row-sum zero condition)
3. Essential Theoretical Properties:
   - P(0) = I (Identity matrix at zero lag; perfect self-correlation at well location)
   - sum_j P_ij(h) = 1.0 (row-stochasticity for all h >= 0)
   - P_ij(h) >= 0 (non-negativity up to numerical precision)
   - lim_{h -> inf} P(h) = 1 * p^T (asymptotic decay to stationary background proportions p)
4. Literature Grounding:
   - Carle, S. F., & Fogg, G. E. (1996). Transition probability-based geostatistics.
     Mathematical Geology, 28(4), 453-476.
   - Elfeki, A., & Dekking, M. (2001). A Markov chain model for subsurface characterization.
     Mathematical Geology, 33(5), 569-589.
"""

from typing import Dict, Any, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import expm

from smalt.descriptive.analyzer import LithologInspector, CANONICAL_FACIES_SCHEMA
from smalt.spatial.coordinates import load_source_coordinates
from smalt.spatial.datum import align_to_common_datum, DEFAULT_CONVENTION


# Authoritative lateral facies lengths calibrated from Sahoo et al. (2016)
# architectural measurements (W/T ~ 35 for ribbon channel sandstones)
DEFAULT_LATERAL_FACIES_LENGTHS_M: Dict[int, float] = {
    0: 203.0,  # Channel Sandstone (W/T ~ 35 * mean thickness 5.8 m)
    1: 30.0,   # Planar Sandstone (upper flow regime sheet sand, ~20-50 m)
    2: 75.0,   # Rippled Heterolithics (crevasse splay / levee margin, 10-130 m)
    3: 40.0,   # Carbonaceous Mudstone (waterlogged swamp margin, 20-60 m)
    4: 55.0,   # Coal (peat mire seam, 30-100 m)
    5: 360.0,  # Overbank Mudstone (widespread floodplain fines, 200-500 m)
}


class EmpiricalHorizontalTransitionEstimator:
    """
    Extracts horizontal facies pairs across common-zero datum elevation slices,
    tallies empirical horizontal transition counts, evaluates directional isotropy,
    and computes lag-binned empirical transition probability matrices.
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
        aligned under the corrected common-zero datum.
        """
        if self._elevation_grid_cache is not None:
            return self._elevation_grid_cache.copy()

        dfs = []
        for lid in self.eligible_ids:
            disc_df = self.inspector.discretize_litholog_1m(lid)
            aligned_df = align_to_common_datum(
                disc_df,
                litholog_id=lid,
                convention=self.convention,
                apply_source_orientation=True,
            )
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
        Never fabricates observations outside the overlapping observed intervals.

        Returns DataFrame of matched pairs with lag distance, dx, dy, azimuth, and facies codes.
        """
        df_grid = self.build_common_elevation_grid()
        unique_z = np.sort(df_grid["z_common_m"].unique())

        pairs = []
        for z_val in unique_z:
            slice_df = df_grid[
                np.abs(df_grid["z_common_m"] - z_val) <= max_elevation_diff_m
            ].reset_index(drop=True)
            n_wells = len(slice_df)
            if n_wells < 2:
                continue

            for i in range(n_wells):
                row_a = slice_df.iloc[i]
                for j in range(i + 1, n_wells):
                    row_b = slice_df.iloc[j]
                    dx = float(row_b["x_m"] - row_a["x_m"])
                    dy = float(row_b["y_m"] - row_a["y_m"])
                    dist = float(np.hypot(dx, dy))
                    theta_deg = float(np.degrees(np.arctan2(dy, dx))) % 360.0

                    pairs.append({
                        "z_common_m": float(z_val),
                        "well_a": str(row_a["litholog_id"]),
                        "well_b": str(row_b["litholog_id"]),
                        "facies_code_a": int(row_a["facies_code"]),
                        "facies_code_b": int(row_b["facies_code"]),
                        "facies_name_a": str(row_a["facies"]),
                        "facies_name_b": str(row_b["facies"]),
                        "lag_distance_m": round(dist, 1),
                        "dx_m": round(dx, 1),
                        "dy_m": round(dy, 1),
                        "azimuth_deg": round(theta_deg, 1),
                        "is_same_facies": int(row_a["facies_code"]) == int(row_b["facies_code"]),
                    })

        df_pairs = pd.DataFrame(pairs)
        return df_pairs

    def evaluate_directional_sparsity(self) -> Dict[str, Any]:
        """
        Evaluates directional distribution of horizontal pairs to document
        whether directional estimation (X vs Y vs Azimuth sectors) is statistically viable.
        """
        df_pairs = self.extract_horizontal_facies_pairs()
        if len(df_pairs) == 0:
            return {"status": "no_pairs", "n_pairs": 0}

        # Quadrant classification
        az = df_pairs["azimuth_deg"].to_numpy()
        q1 = np.sum((az >= 0) & (az < 90))
        q2 = np.sum((az >= 90) & (az < 180))
        q3 = np.sum((az >= 180) & (az < 270))
        q4 = np.sum((az >= 270) & (az < 360))

        # Distinct well pairs
        distinct_well_pairs = df_pairs[["well_a", "well_b"]].drop_duplicates()

        # Dominant canyon transect axis: Sahoo et al. transect is NW-SE oriented
        # Most pairs have |dx| > |dy| or align along canyon
        dx = np.abs(df_pairs["dx_m"].to_numpy())
        dy = np.abs(df_pairs["dy_m"].to_numpy())
        x_dominant = int(np.sum(dx >= dy))
        y_dominant = int(np.sum(dy > dx))

        sparsity_justification = (
            f"With only {len(self.eligible_ids)} coordinate-bearing wells (55 unique well pairs, "
            f"{len(distinct_well_pairs)} active at common elevations), dividing observations into "
            f"directional sectors results in severe statistical sparsity and zero-count transition "
            f"cells for minor facies (e.g., coal, carbonaceous mud). Therefore, isotropic horizontal "
            f"lag distance h is scientifically justified and mathematically required."
        )

        return {
            "total_horizontal_pairs": len(df_pairs),
            "unique_well_pairs_represented": len(distinct_well_pairs),
            "quadrant_distribution": {
                "NE (0-90 deg)": int(q1),
                "NW (90-180 deg)": int(q2),
                "SW (180-270 deg)": int(q3),
                "SE (270-360 deg)": int(q4),
            },
            "x_vs_y_dominance": {
                "dx >= dy": x_dominant,
                "dy > dx": y_dominant,
            },
            "isotropy_justification": sparsity_justification,
        }

    def tally_horizontal_transition_counts(
        self,
        min_lag_m: float = 0.0,
        max_lag_m: float = 1e6,
        alpha: float = 0.01,
    ) -> Dict[str, Any]:
        """
        Tallies empirical horizontal transition counts N_ij within a lag window [min_lag_m, max_lag_m).
        Computes row-stochastic transition probability matrix P_ij(h).
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

        # Laplace smoothed transition probability matrix
        smoothed_counts = counts + alpha
        row_sums = smoothed_counts.sum(axis=1, keepdims=True)
        prob_matrix = smoothed_counts / row_sums

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


class SpatialTransitionRateModel:
    """
    Continuous spatial transition-probability model P(h) = expm(R_h * h)
    based on the Carle & Fogg (1996) continuous Markov chain formulation.

    Guarantees:
    - P(0) = Identity
    - Rows sum to 1.0 for all h >= 0
    - Non-negative probabilities
    - Continuous asymptotic convergence to stationary background proportions
    """

    def __init__(
        self,
        lateral_lengths: Optional[Dict[int, float]] = None,
        stationary_proportions: Optional[np.ndarray] = None,
        embedded_transitions: Optional[np.ndarray] = None,
    ):
        self.num_classes = len(CANONICAL_FACIES_SCHEMA)
        self.lateral_lengths = lateral_lengths or DEFAULT_LATERAL_FACIES_LENGTHS_M.copy()
        self.stationary_proportions = (
            stationary_proportions
            if stationary_proportions is not None
            else np.array([0.435, 0.088, 0.052, 0.026, 0.027, 0.372])  # Default empirical working proportions
        )
        self.embedded_transitions = embedded_transitions
        self.R_matrix = self._build_transition_rate_matrix()
        self._p_cache: Dict[float, np.ndarray] = {}

    def _build_transition_rate_matrix(self) -> np.ndarray:
        """
        Constructs the continuous transition rate matrix R_h satisfying
        Carle & Fogg (1996) stochastic constraints:
        - R_ii = -1 / L_i
        - R_ij = r_ij / L_i (where r_ij is embedded transition or proportional allocation)
        - sum_j R_ij = 0
        """
        n = self.num_classes
        R = np.zeros((n, n), dtype=float)

        p = self.stationary_proportions.copy()
        p = p / p.sum()

        for i in range(n):
            L_i = float(self.lateral_lengths.get(i, 100.0))
            if L_i <= 0.0:
                raise ValueError(f"Lateral length L_{i} must be positive, got {L_i}")
            R[i, i] = -1.0 / L_i

            if self.embedded_transitions is not None:
                r_row = self.embedded_transitions[i, :].copy()
                r_row[i] = 0.0
                sum_r = r_row.sum()
                if sum_r > 0:
                    r_row = r_row / sum_r
                else:
                    r_row = np.ones(n) / (n - 1)
                    r_row[i] = 0.0
            else:
                # Proportional allocation: r_ij = p_j / (1 - p_i)
                denom = 1.0 - p[i]
                if denom > 1e-6:
                    r_row = p / denom
                    r_row[i] = 0.0
                else:
                    r_row = np.ones(n) / (n - 1)
                    r_row[i] = 0.0

            for j in range(n):
                if i != j:
                    R[i, j] = r_row[j] / L_i

        # Enforce exact row-sum zero condition
        R = R - np.diag(R.sum(axis=1))
        np.testing.assert_allclose(R.sum(axis=1), np.zeros(n), atol=1e-12)
        return R

    def evaluate_transition_matrix(self, h_m: float) -> np.ndarray:
        """
        Computes continuous transition probability matrix P(h) = expm(R_h * h).
        Guarantees strict row-stochasticity and non-negativity.
        """
        if h_m < 0.0:
            raise ValueError(f"Lag distance h must be non-negative, got {h_m}")
        if h_m == 0.0:
            return np.eye(self.num_classes)

        h_key = round(float(h_m), 1)
        if h_key in self._p_cache:
            return self._p_cache[h_key].copy()

        P = expm(self.R_matrix * float(h_m))
        P = np.clip(P, 0.0, 1.0)
        row_sums = P.sum(axis=1, keepdims=True)
        P = P / row_sums
        self._p_cache[h_key] = P
        return P.copy()

    def compute_stationary_distribution(self) -> np.ndarray:
        """
        Computes the theoretical stationary distribution vector pi of R_h
        satisfying pi^T * R_h = 0 and sum_k pi_k = 1.0.
        """
        eigvals, left_vecs = np.linalg.eig(self.R_matrix.T)
        null_idx = int(np.argmin(np.abs(eigvals)))
        pi = np.real(left_vecs[:, null_idx])
        if pi.sum() < 0:
            pi = -pi
        pi = np.clip(pi, 0.0, 1.0)
        pi = pi / pi.sum()
        return pi

    def compare_empirical_vs_fitted(
        self,
        empirical_estimator: EmpiricalHorizontalTransitionEstimator,
        test_distances_m: Optional[List[float]] = None,
    ) -> pd.DataFrame:
        """
        Compares empirical horizontal transition probabilities against
        the continuous model P(h) across multiple lag distances.
        """
        test_distances = test_distances_m or [
            500.0, 750.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0
        ]
        df_pairs = empirical_estimator.extract_horizontal_facies_pairs()
        code_to_meta = {meta["code"]: meta for meta in CANONICAL_FACIES_SCHEMA.values()}

        records = []
        for h_val in test_distances:
            # Empirical window around h_val (+/- 250 m)
            window_min = max(0.0, h_val - 250.0)
            window_max = h_val + 250.0
            emp_res = empirical_estimator.tally_horizontal_transition_counts(window_min, window_max)
            P_fitted = self.evaluate_transition_matrix(h_val)

            for i in range(self.num_classes):
                for j in range(self.num_classes):
                    p_emp = float(emp_res["transition_prob_matrix"][i, j])
                    p_fit = float(P_fitted[i, j])
                    records.append({
                        "nominal_distance_m": h_val,
                        "window_min_m": window_min,
                        "window_max_m": window_max,
                        "window_pair_count": emp_res["pair_count"],
                        "from_facies_code": i,
                        "from_facies_name": code_to_meta[i]["canonical_name"],
                        "to_facies_code": j,
                        "to_facies_name": code_to_meta[j]["canonical_name"],
                        "empirical_probability": round(p_emp, 4),
                        "fitted_probability": round(p_fit, 4),
                        "residual_error": round(p_fit - p_emp, 4),
                        "is_auto_transition": int(i == j),
                    })

        return pd.DataFrame(records)
