"""
SMALT Geostatistical Core: 1D Stratigraphic Markov Chain Engine.

Mathematically verified, row-stochastic transition probability modeling,
embedded boundary-crossing transitions, stationary facies occupancy distributions,
directional asymmetry analysis, and state transition entropy.

References:
- Krumbein, W. C., & Dacey, M. F. (1969). Markov chains and embedded Markov chains in geology.
  Mathematical Geology, 1(1), 79-96.
- Doveton, J. H. (1971). An application of Markov chain analysis to the Ayrshire Coal Measures succession.
  Scottish Journal of Geology, 7(1), 11-27.
- Sahoo, H., Gani, M. R., & Gani, N. D. (2016). 3D facies architecture and sequence stratigraphy
  of a fluvio-deltaic succession, Cretaceous Ferron Sandstone, Utah. Sedimentology, 63(6), 1403-1437.
"""

from pathlib import Path
from typing import Dict, Optional, Tuple, Union
import json
import logging
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

DEFAULT_FACIES_MAP: Dict[int, str] = {
    0: "Coal",
    1: "Channel Sandstone",
    2: "Fine Sandstone / Splay",
    3: "Siltstone",
    4: "Overbank Mudstone",
}


class StratigraphicMarkovChain:
    """
    1D Vertical Stratigraphic Markov Chain for facies succession modeling.

    Models the vertical sequence of lithofacies states S_t as a discrete-space,
    discrete-step Markov process where transitions are tallied strictly in upward
    stratigraphic order (decreasing depth_m / increasing stratigraphic height).

    Supports both:
    1. Regular (fixed-step) Markov chains (retaining self-transitions P_ii > 0)
    2. Embedded Markov chains (suppressing diagonal self-transitions P_ii = 0)

    Attributes:
        num_classes (int): Number of discrete lithofacies states K.
        facies_map (Dict[int, str]): Metadata dictionary mapping integer codes to names.
        embedded (bool): Whether self-transitions (P_ii) are suppressed.
        smoothing_alpha (float): Additive Laplace smoothing prior parameter.
        count_matrix_ (np.ndarray): Tallied K x K transition count matrix N.
        transition_matrix_ (np.ndarray): Row-stochastic K x K transition probability matrix P.
        stationary_dist_ (np.ndarray): Invariant stationary distribution vector pi (1 x K).
        n_wells_ (int): Total number of distinct lithologs/wells fitted.
        n_transitions_ (int): Total number of valid vertical transitions tallied.
        n_observations_ (int): Total number of raw depth observations processed.
    """

    def __init__(
        self,
        num_classes: int = 5,
        facies_map: Optional[Dict[int, str]] = None,
        embedded: bool = False,
        smoothing_alpha: float = 0.0,
    ) -> None:
        """
        Initializes the StratigraphicMarkovChain model.

        Args:
            num_classes: Number of discrete states K >= 2.
            facies_map: Optional dictionary mapping integer facies codes {0, ..., K-1}
                to human-readable string names. Defaults to the 5-state SMALT schema.
            embedded: If True, diagonal transitions are excluded prior to normalization (P_ii = 0).
            smoothing_alpha: Non-negative Laplace smoothing parameter alpha >= 0.0.
        """
        if num_classes < 2 and facies_map is None:
            raise ValueError(f"num_classes must be >= 2, got {num_classes}")

        if facies_map is not None:
            self.facies_map = {int(k): str(v) for k, v in facies_map.items()}
            max_code = max(self.facies_map.keys())
            self.num_classes = max(num_classes, max_code + 1, len(self.facies_map))
        else:
            if num_classes == 5:
                self.facies_map = DEFAULT_FACIES_MAP.copy()
            else:
                self.facies_map = {i: f"Facies_{i}" for i in range(num_classes)}
            self.num_classes = num_classes

        self.embedded = bool(embedded)
        self.smoothing_alpha = float(smoothing_alpha)

        if self.smoothing_alpha < 0.0:
            raise ValueError(f"smoothing_alpha must be non-negative, got {smoothing_alpha}")

        # Model state variables
        self.count_matrix_: Optional[np.ndarray] = None
        self.transition_matrix_: Optional[np.ndarray] = None
        self.stationary_dist_: Optional[np.ndarray] = None
        self.n_wells_: int = 0
        self.n_transitions_: int = 0
        self.n_observations_: int = 0

    def fit(
        self,
        df: pd.DataFrame,
        depth_col: str = "depth_m",
        facies_col: str = "facies_code",
        well_col: str = "litholog_id",
    ) -> "StratigraphicMarkovChain":
        """
        Fits the Markov transition probability matrix from a digitized litholog DataFrame.

        Transitions are strictly tallied in upward stratigraphic order (decreasing depth_m).
        Observations and wells are not silently discarded.

        Args:
            df: Validated DataFrame containing digitized litholog data.
            depth_col: Name of the depth column in meters (default: 'depth_m').
            facies_col: Name of the integer facies code column (default: 'facies_code').
            well_col: Name of the litholog/well identifier column (default: 'litholog_id').

        Returns:
            self: Fitted StratigraphicMarkovChain instance.
        """
        if df.empty:
            raise ValueError("Input DataFrame is empty.")

        if depth_col not in df.columns:
            raise ValueError(f"Required depth column '{depth_col}' not found in DataFrame.")
        if facies_col not in df.columns:
            raise ValueError(f"Required facies column '{facies_col}' not found in DataFrame.")

        self.n_observations_ = len(df)
        K = self.num_classes
        N = np.zeros((K, K), dtype=np.float64)

        if well_col in df.columns:
            groups = [group for _, group in df.groupby(well_col, sort=False)]
            self.n_wells_ = len(groups)
        else:
            groups = [df]
            self.n_wells_ = 1

        total_transitions = 0
        for well_df in groups:
            # Sort descending by depth so iteration proceeds upward stratigraphically:
            # lower bed (z) -> upper bed (z - delta_z)
            sorted_well = well_df.sort_values(by=depth_col, ascending=False)
            facies_seq = sorted_well[facies_col].to_numpy()

            if len(facies_seq) < 2:
                continue

            for from_state, to_state in zip(facies_seq[:-1], facies_seq[1:]):
                from_idx = int(from_state)
                to_idx = int(to_state)

                if not (0 <= from_idx < K):
                    raise ValueError(
                        f"Facies code {from_idx} in '{facies_col}' is outside valid range [0, {K-1}]."
                    )
                if not (0 <= to_idx < K):
                    raise ValueError(
                        f"Facies code {to_idx} in '{facies_col}' is outside valid range [0, {K-1}]."
                    )

                N[from_idx, to_idx] += 1.0
                total_transitions += 1

        # Embedded chain: suppress diagonal counts before normalization & smoothing
        if self.embedded:
            np.fill_diagonal(N, 0.0)
            self.count_matrix_ = N.copy()
            self.n_transitions_ = int(np.sum(N))
        else:
            self.count_matrix_ = N.copy()
            self.n_transitions_ = total_transitions

        # Maximum Likelihood Estimation with Laplace Prior
        P = np.zeros((K, K), dtype=np.float64)
        alpha = self.smoothing_alpha

        for i in range(K):
            if self.embedded:
                row_counts = N[i, :].copy()
                row_counts[i] = 0.0
                row_sum_raw = np.sum(row_counts)

                if row_sum_raw == 0.0 and alpha == 0.0:
                    # Unobserved/absorbing state in embedded chain: uniform fallback over K-1 states
                    logger.warning(
                        f"State {i} ({self.facies_map.get(i, i)}) has zero observed outgoing transitions. "
                        f"Assigning uniform fallback across other {K - 1} states."
                    )
                    P[i, :] = 1.0 / (K - 1)
                    P[i, i] = 0.0
                else:
                    denom = row_sum_raw + (K - 1) * alpha
                    P[i, :] = (row_counts + alpha) / denom
                    P[i, i] = 0.0
            else:
                row_counts = N[i, :]
                row_sum_raw = np.sum(row_counts)

                if row_sum_raw == 0.0 and alpha == 0.0:
                    # Unobserved/absorbing state in regular chain: uniform fallback over all K states
                    logger.warning(
                        f"State {i} ({self.facies_map.get(i, i)}) has zero observed outgoing transitions. "
                        f"Assigning uniform fallback across all {K} states."
                    )
                    P[i, :] = 1.0 / K
                else:
                    denom = row_sum_raw + K * alpha
                    P[i, :] = (row_counts + alpha) / denom

        # Guard against numerical drift: renormalize strictly to 1.0
        row_sums = P.sum(axis=1, keepdims=True)
        P = P / row_sums

        # Numerical Invariant Verification
        if not np.allclose(P.sum(axis=1), 1.0, atol=1e-7):
            raise AssertionError(
                f"Row-stochasticity invariant violated: max row sum error = "
                f"{np.max(np.abs(P.sum(axis=1) - 1.0)):.2e} > 1e-7"
            )

        if not np.all((P >= 0.0) & (P <= 1.0 + 1e-7)):
            raise AssertionError("Probability bounds violation: elements outside [0, 1].")

        P = np.clip(P, 0.0, 1.0)
        P = P / P.sum(axis=1, keepdims=True)

        if self.embedded and not np.all(np.diag(P) == 0.0):
            raise AssertionError("Embedded chain diagonal invariant violated: non-zero diagonal present.")

        self.transition_matrix_ = P
        self.compute_stationary_distribution()
        return self

    def compute_stationary_distribution(self) -> np.ndarray:
        """
        Computes the stationary distribution vector pi = [pi_0, ..., pi_{K-1}].

        Solves the invariant left-eigenvector equation:
            pi * P = pi  <=>  P^T * pi^T = pi^T
        subject to:
            sum(pi_i) = 1.0 and pi_i >= 0.

        Cross-checks result against stationary invariant ||pi * P - pi||_inf < 1e-6.

        Returns:
            pi: Stationary distribution array of shape (K,).
        """
        if self.transition_matrix_ is None:
            raise ValueError("Model must be fitted before computing stationary distribution.")

        P = self.transition_matrix_
        K = self.num_classes

        # Solve left eigensystem via transpose P^T
        eigvals, eigvecs = np.linalg.eig(P.T)
        idx = np.argmin(np.abs(eigvals - 1.0))
        pi = np.real(eigvecs[:, idx])

        # Clean numerical zero noise before evaluating sign
        pi[np.abs(pi) < 1e-12] = 0.0

        # Enforce correct global sign (dominant non-zero components must be positive)
        if np.sum(pi) < 0.0:
            pi = -pi
        elif np.sum(pi) == 0.0 and np.any(pi != 0):
            first_nonzero = pi[np.flatnonzero(pi)[0]]
            if first_nonzero < 0:
                pi = -pi

        pi = np.clip(pi, 0.0, None)

        sum_pi = np.sum(pi)
        if sum_pi == 0.0:
            pi = np.ones(K, dtype=np.float64) / K
        else:
            pi = pi / sum_pi

        # Numerical invariance check: ||pi * P - pi||_inf < 1e-6
        invariant_diff = float(np.max(np.abs(pi @ P - pi)))
        if invariant_diff >= 1e-6:
            # Fallback to least-squares solver for singular/degenerate systems:
            # Solve (P^T - I) pi = 0 with sum(pi) = 1
            A_eq = np.vstack([P.T - np.eye(K), np.ones((1, K))])
            b_eq = np.zeros(K + 1)
            b_eq[-1] = 1.0
            pi_lstsq, _, _, _ = np.linalg.lstsq(A_eq, b_eq, rcond=None)
            pi_lstsq = np.clip(pi_lstsq, 0.0, None)
            if np.sum(pi_lstsq) > 0:
                pi = pi_lstsq / np.sum(pi_lstsq)
                invariant_diff = float(np.max(np.abs(pi @ P - pi)))

            if invariant_diff >= 1e-6:
                raise AssertionError(
                    f"Stationary distribution invariant violated: ||pi * P - pi||_inf = {invariant_diff:.2e} >= 1e-6"
                )

        self.stationary_dist_ = pi
        return pi

    def compute_asymmetry_matrix(self) -> np.ndarray:
        """
        Computes the directional asymmetry matrix A = P - P^T.

        A_ij = P_ij - P_ji quantifies the net directional preference of the upward
        stratigraphic transition i -> j relative to the downward reverse j -> i.

        Returns:
            A: Asymmetry matrix of shape (K, K).
        """
        if self.transition_matrix_ is None:
            raise ValueError("Model must be fitted before computing asymmetry matrix.")
        return self.transition_matrix_ - self.transition_matrix_.T

    def compute_state_entropy(self, normalized: bool = False) -> np.ndarray:
        """
        Computes the row/state transition entropy H_i for each facies state.

        H_i = - sum_{j, P_ij > 0} P_ij * ln(P_ij)

        Args:
            normalized: If True, returns H_i / H_max in [0, 1], where
                H_max = ln(K) for regular chains and ln(K-1) for embedded chains.

        Returns:
            H: Array of entropy values of shape (K,).
        """
        if self.transition_matrix_ is None:
            raise ValueError("Model must be fitted before computing state entropy.")

        P = self.transition_matrix_
        K = self.num_classes
        H = np.zeros(K, dtype=np.float64)

        for i in range(K):
            row = P[i, :]
            pos_mask = row > 0.0
            if np.any(pos_mask):
                H[i] = -np.sum(row[pos_mask] * np.log(row[pos_mask]))

        if normalized:
            h_max = np.log(K - 1) if self.embedded and K > 1 else (np.log(K) if K > 1 else 1.0)
            if h_max > 0.0:
                H = H / h_max

        return H

    def sample_next_state(
        self, current_state: int, rng: Optional[np.random.Generator] = None
    ) -> int:
        """
        Draws a single categorical state realization from transition row P[current_state, :].

        Args:
            current_state: Integer index of the current facies state.
            rng: Optional numpy Generator for reproducible stochastic sampling.

        Returns:
            next_state: Integer index of the realized succeeding facies state.
        """
        if self.transition_matrix_ is None:
            raise ValueError("Model must be fitted before sampling transitions.")

        if not (0 <= current_state < self.num_classes):
            raise ValueError(
                f"current_state {current_state} is out of bounds [0, {self.num_classes - 1}]."
            )

        gen = rng if rng is not None else np.random.default_rng()
        probs = self.transition_matrix_[int(current_state), :]
        return int(gen.choice(self.num_classes, p=probs))

    def to_dataframe(self) -> pd.DataFrame:
        """
        Returns the transition probability matrix P as a labeled pandas DataFrame.

        Rows: 'From' facies S_{t-1} (Lower bed)
        Columns: 'To' facies S_t (Upper bed)
        """
        if self.transition_matrix_ is None:
            raise ValueError("Model must be fitted before exporting DataFrame.")

        labels = [self.facies_map.get(i, f"State_{i}") for i in range(self.num_classes)]
        return pd.DataFrame(self.transition_matrix_, index=labels, columns=labels)

    def export_json(self, filepath: Union[str, Path]) -> Path:
        """
        Exports the fitted Markov chain parameters, counts, and metadata to JSON.

        Args:
            filepath: Destination file path.

        Returns:
            path: Resolved Path object of the exported file.
        """
        if self.transition_matrix_ is None or self.count_matrix_ is None:
            raise ValueError("Model must be fitted before exporting JSON.")

        path = Path(filepath)
        path.parent.mkdir(parents=True, exist_ok=True)

        asymmetry = self.compute_asymmetry_matrix()
        entropy_raw = self.compute_state_entropy(normalized=False)
        entropy_norm = self.compute_state_entropy(normalized=True)

        payload = {
            "metadata": {
                "num_classes": self.num_classes,
                "embedded": self.embedded,
                "smoothing_alpha": self.smoothing_alpha,
                "n_wells": self.n_wells_,
                "n_transitions": self.n_transitions_,
                "n_observations": self.n_observations_,
                "facies_map": {str(k): v for k, v in self.facies_map.items()},
            },
            "count_matrix_N": self.count_matrix_.tolist(),
            "transition_matrix_P": self.transition_matrix_.tolist(),
            "stationary_distribution_pi": (
                self.stationary_dist_.tolist() if self.stationary_dist_ is not None else []
            ),
            "asymmetry_matrix_A": asymmetry.tolist(),
            "state_entropy_H": entropy_raw.tolist(),
            "state_entropy_H_norm": entropy_norm.tolist(),
        }

        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)

        return path
