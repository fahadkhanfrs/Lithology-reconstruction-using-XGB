"""
SMALT Geostatistical Core: Conditioned Transition-Probability / Coupled-Markov Prototype.

Implements the "SMALT Conditioned Transition-Probability Prototype" (SMALT-CTP),
combining 1D vertical Markov succession (Krumbein & Dacey, 1969; Doveton, 1971),
continuous horizontal transition rate modeling (Carle & Fogg, 1996), and
coupled Markov conditioning (Elfeki & Dekking, 2001, 2005) for sparse multi-well networks.

Theoretical Distinctions & Adaptations:
1. Literature-Derived Principles:
   - Carle & Fogg (1996): Continuous spatial transition probability P_h(h) = expm(R_h * h)
     as a function of separation distance h.
   - Elfeki & Dekking (2001): Coupling vertical transition P_v and horizontal transition P_h
     via normalized product / Bayesian likelihood pooling.
2. SMALT Adaptation for Sparse Wells:
   - Rather than assuming two parallel bounding wells on a 2D regular grid (standard CMC),
     SMALT conditions on an irregular 3D network of sparse observed lithologs (L2-L12).
   - Horizontal conditioning from multiple neighboring wells is aggregated via inverse-distance
     weighting of their respective transition probability vectors.
   - Hard data conditioning is enforced with 100% honor rate at known well locations.
   - In cross-validation (LOLO), the target litholog is strictly withheld from parameter
     estimation, transition rate fitting, and conditioning sets.
"""

from typing import Dict, Any, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import entropy as scipy_entropy

from smalt.descriptive.analyzer import LithologInspector, CANONICAL_FACIES_SCHEMA
from smalt.geostat.markov import StratigraphicMarkovChain
from smalt.geostat.spatial_transition import (
    SpatialTransitionRateModel,
    DEFAULT_LATERAL_FACIES_LENGTHS_M,
)
from smalt.spatial.coordinates import load_source_coordinates
from smalt.spatial.datum import align_to_common_datum, DEFAULT_CONVENTION

CODE_TO_NAME: Dict[int, str] = {
    meta["code"]: meta["canonical_name"] for meta in CANONICAL_FACIES_SCHEMA.values()
}


class ConditionedMarkovClassifier:
    """
    Conditioned Transition-Probability / Coupled-Markov spatial predictor.

    Predicts the conditional facies probability vector P(S(x, y, z) = k | W_obs)
    for all six canonical facies at an unobserved location (x, y, z), honoring
    conditioning wells as hard data and coupling horizontal and vertical transitions.
    """

    def __init__(
        self,
        training_litholog_ids: List[str],
        inspector: Optional[LithologInspector] = None,
        convention: str = DEFAULT_CONVENTION,
        lateral_lengths: Optional[Dict[int, float]] = None,
        distance_power: float = 1.0,
        vertical_coupling_weight: float = 1.0,
    ):
        """
        Args:
            training_litholog_ids: List of litholog IDs available for training/conditioning.
                                   Target validation wells MUST NOT be included.
            inspector: LithologInspector instance.
            convention: Vertical datum alignment convention (default: common_zero).
            lateral_lengths: Dictionary of mean lateral lengths per facies.
            distance_power: Exponent for inverse-distance weighting of conditioning wells.
            vertical_coupling_weight: Power weighting for vertical Markov prior vs horizontal likelihood.
        """
        self.training_ids = sorted(list(training_litholog_ids))
        self.inspector = inspector or LithologInspector()
        self.convention = convention
        self.lateral_lengths = lateral_lengths or DEFAULT_LATERAL_FACIES_LENGTHS_M.copy()
        self.distance_power = distance_power
        self.vertical_coupling_weight = vertical_coupling_weight

        self.coords_df = load_source_coordinates()
        self.coords_by_id = {
            r["litholog_id"]: (float(r["x_m"]), float(r["y_m"]))
            for _, r in self.coords_df.iterrows()
            if r["coordinates_available"] and r["litholog_id"] in self.training_ids
        }

        # Fit vertical Markov chain strictly on training lithologs
        self.vertical_markov = self._fit_vertical_markov()
        self.stationary_proportions = self.vertical_markov.stationary_dist_.flatten()

        # Build continuous spatial transition rate model
        self.spatial_rate_model = SpatialTransitionRateModel(
            lateral_lengths=self.lateral_lengths,
            stationary_proportions=self.stationary_proportions,
        )

        # Build cached training observations table for fast conditioning lookup
        self.training_obs_df = self._build_training_observations()

        # Build fast O(1) elevation slice and hard-data lookup indices
        self._well_z_to_code: Dict[Tuple[str, int], int] = {}
        self._obs_by_z: Dict[int, List[Dict[str, Any]]] = {}
        if len(self.training_obs_df) > 0:
            for _, r in self.training_obs_df.iterrows():
                z_int = int(round(float(r["z_common_m"])))
                lid = str(r["litholog_id"])
                code = int(r["facies_code"])
                self._well_z_to_code[(lid, z_int)] = code
                if z_int not in self._obs_by_z:
                    self._obs_by_z[z_int] = []
                self._obs_by_z[z_int].append({
                    "x_m": float(r["x_m"]),
                    "y_m": float(r["y_m"]),
                    "facies_code": code,
                    "litholog_id": lid,
                })

    def _fit_vertical_markov(self) -> StratigraphicMarkovChain:
        """
        Fits 1D vertical StratigraphicMarkovChain exclusively on training lithologs.
        """
        chain = StratigraphicMarkovChain(num_classes=len(CANONICAL_FACIES_SCHEMA), smoothing_alpha=0.01)
        train_dfs = []
        for lid in self.training_ids:
            disc_df = self.inspector.discretize_litholog_1m(lid)
            train_dfs.append(disc_df)
        combined = pd.concat(train_dfs, ignore_index=True)
        chain.fit(combined)
        return chain

    def _build_training_observations(self) -> pd.DataFrame:
        """
        Extracts 1m discretized, common-zero aligned observations for all training wells.
        """
        dfs = []
        for lid in self.training_ids:
            if lid not in self.coords_by_id:
                continue
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

        if len(dfs) == 0:
            return pd.DataFrame()
        return pd.concat(dfs, ignore_index=True)

    def evaluate_conditional_probability_at_point(
        self,
        x_m: float,
        y_m: float,
        z_common_m: float,
        prev_vertical_facies: Optional[int] = None,
        hard_tolerance_m: float = 1.0,
    ) -> Dict[str, Any]:
        """
        Computes conditional facies probability distribution at (x, y, z):
        P(S(x, y, z) = k | conditioning lithologs)

        Returns dictionary with probabilities for each facies, MAP facies, and entropy.
        """
        n_classes = len(CANONICAL_FACIES_SCHEMA)
        z_int = int(round(z_common_m))

        # 1. Hard Conditioning Check
        if hard_tolerance_m > 0.0:
            for lid, (wx, wy) in self.coords_by_id.items():
                dist_to_well = float(np.hypot(x_m - wx, y_m - wy))
                if dist_to_well <= hard_tolerance_m:
                    true_code = self._well_z_to_code.get((lid, z_int))
                    if true_code is not None:
                        prob_vec = np.zeros(n_classes, dtype=float)
                        prob_vec[true_code] = 1.0
                        return {
                            "x_m": round(x_m, 1),
                            "y_m": round(y_m, 1),
                            "z_common_m": round(z_common_m, 1),
                            "probabilities": prob_vec,
                            "map_facies_code": true_code,
                            "map_facies_name": CODE_TO_NAME[true_code],
                            "predictive_entropy": 0.0,
                            "is_hard_conditioned": True,
                            "nearest_well_id": lid,
                            "nearest_distance_m": round(dist_to_well, 1),
                        }

        # 2. Extract active conditioning observations at matching elevation slice
        slice_obs = self._obs_by_z.get(z_int)
        if not slice_obs:
            if len(self._obs_by_z) > 0:
                nearest_z = min(self._obs_by_z.keys(), key=lambda k: abs(k - z_int))
                slice_obs = self._obs_by_z[nearest_z]
            else:
                slice_obs = []

        # 3. Horizontal Transition Likelihood Aggregation
        well_dists = []
        well_probs = []
        for r in slice_obs:
            wx = r["x_m"]
            wy = r["y_m"]
            w_code = r["facies_code"]
            dist = max(1.0, float(np.hypot(x_m - wx, y_m - wy)))
            P_h = self.spatial_rate_model.evaluate_transition_matrix(dist)
            well_dists.append(dist)
            well_probs.append(P_h[w_code, :])

        if len(well_dists) > 0:
            inv_dists = np.array([1.0 / (d ** self.distance_power) for d in well_dists])
            weights = inv_dists / inv_dists.sum()
            horizontal_prob = np.sum(np.array(well_probs) * weights[:, np.newaxis], axis=0)
            nearest_dist = min(well_dists)
            nearest_idx = np.argmin(well_dists)
            nearest_well = slice_obs[nearest_idx]["litholog_id"]
        else:
            horizontal_prob = self.stationary_proportions.copy()
            nearest_dist = 9999.0
            nearest_well = "none"

        # 4. Vertical Succession Prior
        if prev_vertical_facies is not None and 0 <= prev_vertical_facies < n_classes:
            vertical_prior = self.vertical_markov.transition_matrix_[prev_vertical_facies, :]
        else:
            vertical_prior = self.stationary_proportions.copy()

        # 5. Coupled Conditioning (Elfeki & Dekking 2001 coupling)
        # P_coupled(k) proportional to P_v(k)^beta * P_h(k)
        coupled_unnorm = (vertical_prior ** self.vertical_coupling_weight) * horizontal_prob
        coupled_sum = coupled_unnorm.sum()
        if coupled_sum > 0:
            final_prob = coupled_unnorm / coupled_sum
        else:
            final_prob = np.ones(n_classes) / n_classes

        # Final cleanup for numerical stability
        final_prob = np.clip(final_prob, 0.0, 1.0)
        final_prob = final_prob / final_prob.sum()

        map_code = int(np.argmax(final_prob))
        ent = float(scipy_entropy(final_prob + 1e-15))

        return {
            "x_m": round(x_m, 1),
            "y_m": round(y_m, 1),
            "z_common_m": round(z_common_m, 1),
            "probabilities": final_prob,
            "map_facies_code": map_code,
            "map_facies_name": CODE_TO_NAME[map_code],
            "predictive_entropy": round(ent, 4),
            "is_hard_conditioned": False,
            "nearest_well_id": nearest_well,
            "nearest_distance_m": round(nearest_dist, 1),
        }

    def compute_hard_data_honor_rate(self) -> float:
        """
        Evaluates predictions at all training well observation points to verify
        that hard conditioning is honored with 100% precision.
        """
        df_obs = self.training_obs_df
        if len(df_obs) == 0:
            return 1.0

        n_total = len(df_obs)
        n_honored = 0
        for _, r in df_obs.iterrows():
            res = self.evaluate_conditional_probability_at_point(
                x_m=float(r["x_m"]),
                y_m=float(r["y_m"]),
                z_common_m=float(r["z_common_m"]),
                hard_tolerance_m=1.0,
            )
            if res["map_facies_code"] == int(r["facies_code"]) and res["is_hard_conditioned"]:
                n_honored += 1

        return float(n_honored / n_total)

    def predict_target_litholog(
        self,
        target_litholog_id: str,
    ) -> pd.DataFrame:
        """
        Predicts facies probability field along a held-out target litholog.
        Strictly leak-free: target well is evaluated as an unobserved location.
        """
        coords = load_source_coordinates()
        coord_row = coords[coords["litholog_id"] == target_litholog_id]
        if len(coord_row) == 0 or not bool(coord_row.iloc[0]["coordinates_available"]):
            raise ValueError(f"Target litholog {target_litholog_id} has no coordinates or is unknown.")

        tx = float(coord_row.iloc[0]["x_m"])
        ty = float(coord_row.iloc[0]["y_m"])

        disc_df = self.inspector.discretize_litholog_1m(target_litholog_id)
        aligned_df = align_to_common_datum(
            disc_df,
            litholog_id=target_litholog_id,
            convention=self.convention,
            apply_source_orientation=True,
        )

        records = []
        prev_pred = None
        for _, r in aligned_df.iterrows():
            z_val = float(r["z_common_m"])
            true_code = int(r["facies_code"])

            # Evaluate conditioned probability without hard-conditioning to target itself
            res = self.evaluate_conditional_probability_at_point(
                x_m=tx,
                y_m=ty,
                z_common_m=z_val,
                prev_vertical_facies=prev_pred,
                hard_tolerance_m=0.0,  # Zero tolerance so target is treated as unobserved
            )
            prev_pred = res["map_facies_code"]

            probs = res["probabilities"]
            records.append({
                "litholog_id": target_litholog_id,
                "x_m": tx,
                "y_m": ty,
                "z_common_m": z_val,
                "true_facies_code": true_code,
                "true_facies_name": CODE_TO_NAME[true_code],
                "map_facies_code": res["map_facies_code"],
                "map_facies_name": res["map_facies_name"],
                "P_sand": round(float(probs[0]), 4),
                "P_p_sand": round(float(probs[1]), 4),
                "P_ripples": round(float(probs[2]), 4),
                "P_carbon_mud": round(float(probs[3]), 4),
                "P_coal": round(float(probs[4]), 4),
                "P_mud": round(float(probs[5]), 4),
                "predictive_entropy": res["predictive_entropy"],
                "nearest_well_id": res["nearest_well_id"],
                "nearest_distance_m": res["nearest_distance_m"],
            })

        return pd.DataFrame(records)
