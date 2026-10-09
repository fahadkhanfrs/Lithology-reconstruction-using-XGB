"""
SMALT Geostatistical Core: Stochastic Inter-Well Facies Realization Generator.

Implements sequential / coupled Markov conditional simulation for inter-well
stratigraphic domains, generating ensembles of equiprobable stochastic realizations
honoring observed lithologs as hard data (Elfeki & Dekking, 2001, 2005; Carle & Fogg, 1996).

Key Capabilities:
1. Hard Data Honoring: 100% exact facies reproduction at conditioning well columns.
2. Coupled Markov Transition Sampling: Simulates inter-well cells by sampling from
   conditional categorical distributions combining horizontal decay and vertical succession.
3. Multi-Realization Ensemble: Produces M stochastic realizations under a fixed random seed.
4. Ensemble Spatial Uncertainty: Quantifies per-cell modal facies, facies occurrence probabilities,
   and spatial Shannon entropy across the ensemble.
"""

from typing import Dict, Any, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import entropy as scipy_entropy

from smalt.descriptive.analyzer import LithologInspector, CANONICAL_FACIES_SCHEMA
from smalt.geostat.conditioned_markov import ConditionedMarkovClassifier
from smalt.spatial.datum import align_to_common_datum, DEFAULT_CONVENTION


class InterWellRealizationGenerator:
    """
    Generates stochastic 2D/3D inter-well facies realizations conditioned on observed lithologs.
    """

    def __init__(
        self,
        classifier: ConditionedMarkovClassifier,
        random_seed: int = 42,
    ):
        """
        Args:
            classifier: ConditionedMarkovClassifier instance with fitted models and conditioning wells.
            random_seed: Base integer seed for deterministic reproducibility.
        """
        self.classifier = classifier
        self.random_seed = int(random_seed)
        self.num_classes = len(CANONICAL_FACIES_SCHEMA)

    def generate_2d_transect_realizations(
        self,
        anchor_well_ids: List[str],
        n_realizations: int = 5,
        x_resolution_m: float = 25.0,
        z_resolution_m: float = 1.0,
        z_min_m: Optional[float] = None,
        z_max_m: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        Generates stochastic 2D cross-sectional realizations along a transect connecting anchor wells.

        Args:
            anchor_well_ids: Ordered list of well IDs defining the transect path (e.g. ['litholog2', 'litholog9']).
            n_realizations: Number of stochastic realizations to generate (M >= 3).
            x_resolution_m: Lateral cell spacing along the transect in meters.
            z_resolution_m: Vertical cell spacing in meters.
            z_min_m: Lower vertical datum limit (default: min across anchor wells).
            z_max_m: Upper vertical datum limit (default: max across anchor wells).

        Returns:
            Dictionary containing:
            - realizations: Array of shape (n_realizations, n_z, n_x) of facies codes.
            - x_coords_m: Lateral distance coordinates along transect (1D).
            - z_coords_m: Vertical common-zero elevation coordinates (1D).
            - well_positions_x_m: Projected transect locations for each anchor well.
            - ensemble_mode: Modal facies grid (n_z, n_x).
            - ensemble_entropy: Shannon entropy grid (n_z, n_x).
            - facies_probabilities: Probability grids (n_classes, n_z, n_x).
            - hard_data_honor_rate: Verification metric (strictly 1.0).
        """
        rng = np.random.RandomState(self.random_seed)
        coords_by_id = self.classifier.coords_by_id

        # Verify all anchor wells have spatial coordinates
        for wid in anchor_well_ids:
            if wid not in coords_by_id:
                raise ValueError(f"Anchor well '{wid}' is not in conditioning coordinates.")

        # Compute cumulative distance along transect
        transect_pts = [coords_by_id[wid] for wid in anchor_well_ids]
        cum_dists = [0.0]
        for i in range(1, len(transect_pts)):
            d = float(np.hypot(
                transect_pts[i][0] - transect_pts[i-1][0],
                transect_pts[i][1] - transect_pts[i-1][1],
            ))
            cum_dists.append(cum_dists[-1] + d)

        total_length = cum_dists[-1]
        x_grid = np.arange(0.0, total_length + 0.1, x_resolution_m)
        nx = len(x_grid)

        # Map each anchor well to nearest transect grid index
        well_grid_indices = {}
        for wid, d in zip(anchor_well_ids, cum_dists):
            idx = int(np.argmin(np.abs(x_grid - d)))
            well_grid_indices[wid] = idx

        # Determine vertical grid extent
        df_obs = self.classifier.training_obs_df
        anchor_obs = df_obs[df_obs["litholog_id"].isin(anchor_well_ids)]
        z_min = z_min_m if z_min_m is not None else float(anchor_obs["z_common_m"].min())
        z_max = z_max_m if z_max_m is not None else float(anchor_obs["z_common_m"].max())
        z_grid = np.arange(z_min, z_max + 0.1, z_resolution_m)
        nz = len(z_grid)

        # Pre-populate hard conditioning mask and values
        hard_mask = np.zeros((nz, nx), dtype=bool)
        hard_values = np.full((nz, nx), -1, dtype=int)

        for wid, x_idx in well_grid_indices.items():
            w_obs = anchor_obs[anchor_obs["litholog_id"] == wid]
            for _, r in w_obs.iterrows():
                z_val = float(r["z_common_m"])
                z_idx = int(np.argmin(np.abs(z_grid - z_val)))
                if 0 <= z_idx < nz:
                    hard_mask[z_idx, x_idx] = True
                    hard_values[z_idx, x_idx] = int(r["facies_code"])

        # Interpolate spatial (X, Y) coordinates along transect
        interp_xy = []
        for x_val in x_grid:
            # Find segment
            seg_idx = 0
            for s in range(len(cum_dists) - 1):
                if cum_dists[s] <= x_val <= cum_dists[s+1]:
                    seg_idx = s
                    break
            seg_len = max(1e-6, cum_dists[seg_idx+1] - cum_dists[seg_idx])
            t = (x_val - cum_dists[seg_idx]) / seg_len
            x_real = transect_pts[seg_idx][0] + t * (transect_pts[seg_idx+1][0] - transect_pts[seg_idx][0])
            y_real = transect_pts[seg_idx][1] + t * (transect_pts[seg_idx+1][1] - transect_pts[seg_idx][1])
            interp_xy.append((x_real, y_real))

        # Generate M realizations
        realizations = np.zeros((n_realizations, nz, nx), dtype=int)

        for m in range(n_realizations):
            grid_m = np.copy(hard_values)

            # Sequential simulation: column-by-column, bottom-to-top (stratigraphic order)
            for j in range(nx):
                rx, ry = interp_xy[j]
                prev_facies = None

                for i in range(nz):
                    if hard_mask[i, j]:
                        # Exact hard conditioning
                        grid_m[i, j] = hard_values[i, j]
                        prev_facies = hard_values[i, j]
                        continue

                    z_val = z_grid[i]

                    # Evaluate conditioned probability distribution
                    # Coupling vertical neighbor below (prev_facies) and horizontal wells
                    cond_res = self.classifier.evaluate_conditional_probability_at_point(
                        x_m=rx,
                        y_m=ry,
                        z_common_m=z_val,
                        prev_vertical_facies=prev_facies,
                        hard_tolerance_m=0.0,
                    )
                    prob_vec = cond_res["probabilities"]

                    # Draw from categorical distribution
                    drawn_code = int(rng.choice(self.num_classes, p=prob_vec))
                    grid_m[i, j] = drawn_code
                    prev_facies = drawn_code

            realizations[m] = grid_m

        # Verify hard data honor rate
        n_hard_pts = int(np.sum(hard_mask))
        if n_hard_pts > 0:
            honored_count = 0
            for m in range(n_realizations):
                honored_count += int(np.sum(realizations[m][hard_mask] == hard_values[hard_mask]))
            total_evals = n_realizations * n_hard_pts
            hard_honor_rate = float(honored_count / total_evals)
        else:
            hard_honor_rate = 1.0

        # Compute ensemble summary grids
        # 1. Class probability maps: P(S = k)
        prob_maps = np.zeros((self.num_classes, nz, nx), dtype=float)
        for k in range(self.num_classes):
            prob_maps[k] = np.mean(realizations == k, axis=0)

        # 2. Ensemble mode (MAP facies)
        ensemble_mode = np.argmax(prob_maps, axis=0)

        # 3. Ensemble Shannon entropy
        ensemble_entropy = np.zeros((nz, nx), dtype=float)
        for i in range(nz):
            for j in range(nx):
                p_vec = prob_maps[:, i, j]
                ensemble_entropy[i, j] = float(scipy_entropy(p_vec + 1e-15))

        return {
            "realizations": realizations,
            "x_coords_m": x_grid,
            "z_coords_m": z_grid,
            "anchor_wells": anchor_well_ids,
            "anchor_well_distances_m": cum_dists,
            "well_grid_indices": well_grid_indices,
            "hard_data_mask": hard_mask,
            "ensemble_mode": ensemble_mode,
            "ensemble_entropy": np.round(ensemble_entropy, 4),
            "facies_probabilities": np.round(prob_maps, 4),
            "hard_data_honor_rate": round(hard_honor_rate, 4),
            "n_realizations": n_realizations,
        }
