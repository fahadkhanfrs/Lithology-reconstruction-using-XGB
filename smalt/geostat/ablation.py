"""
SMALT Sprint J: Ablation Study & Comparative Realization Engine.

Systematically compares four architectural configurations:
- Model 0: Stationary Facies Proportions Only (independent categorical sampling from pi)
- Model 1: Vertical Markov Succession Only (1D stratigraphic Markov chain, unconditioned laterally)
- Model 2: Vertical Markov + Sprint I Horizontal Model (prescribed continuous rate matrix R_A)
- Model 3: Vertical Markov + Sprint J Empirically Calibrated Horizontal Model (calibrated rate matrix R_B)

Evaluates:
1. Facies proportions (Total Variation distance vs true target distribution)
2. Vertical transition matrix divergence (Frobenius norm)
3. Mean bed thickness and bed thickness error
4. Bed thickness distribution distance (Wasserstein-1 distance)
5. Hard-data honor rate (strictly 100% at conditioning well locations)
6. Spatial connectivity (mean lateral run length of channel sand bodies)
7. Spatial Shannon entropy
8. Multi-realization ensemble variability
"""

from typing import Dict, Any, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance, entropy as scipy_entropy

from smalt.descriptive.analyzer import LithologInspector, CANONICAL_FACIES_SCHEMA
from smalt.geostat.markov import StratigraphicMarkovChain
from smalt.geostat.spatial_transition import (
    SpatialTransitionRateModel,
    DEFAULT_LATERAL_FACIES_LENGTHS_M,
)
from smalt.geostat.conditioned_markov import ConditionedMarkovClassifier
from smalt.geostat.realization import InterWellRealizationGenerator
from smalt.geostat.empirical_calibration import HorizontalContinuityCalibrator
from smalt.spatial.coordinates import load_source_coordinates
from smalt.spatial.datum import align_to_common_datum, DEFAULT_CONVENTION

CODE_TO_NAME: Dict[int, str] = {
    meta["code"]: meta["canonical_name"] for meta in CANONICAL_FACIES_SCHEMA.values()
}


def calculate_bed_thicknesses(facies_sequence: np.ndarray, cell_height_m: float = 1.0) -> List[float]:
    """
    Computes individual bed thicknesses in meters from a 1D vertical facies profile.
    A bed is defined as a contiguous run of the same facies state.
    """
    if len(facies_sequence) == 0:
        return []
    runs = []
    current_len = 1
    for k in range(1, len(facies_sequence)):
        if facies_sequence[k] == facies_sequence[k - 1]:
            current_len += 1
        else:
            runs.append(float(current_len * cell_height_m))
            current_len = 1
    runs.append(float(current_len * cell_height_m))
    return runs


def compute_vertical_transition_matrix(
    facies_grid: np.ndarray,
    num_classes: int = 6,
    alpha: float = 0.01,
) -> np.ndarray:
    """
    Computes empirical vertical transition matrix (upward transitions) across
    all vertical columns of a 2D realization grid (shape nz x nx).
    """
    nz, nx = facies_grid.shape
    counts = np.zeros((num_classes, num_classes), dtype=float)
    for col in range(nx):
        seq = facies_grid[:, col]
        for z in range(nz - 1):
            s_from = int(seq[z])
            s_to = int(seq[z + 1])
            counts[s_from, s_to] += 1.0
    smoothed = counts + alpha
    return smoothed / smoothed.sum(axis=1, keepdims=True)


def compute_lateral_connectivity_length(
    facies_grid: np.ndarray,
    target_facies_code: int = 0,
    cell_width_m: float = 25.0,
) -> float:
    """
    Computes mean lateral run-length (horizontal connectivity in meters) of target facies.
    """
    nz, nx = facies_grid.shape
    run_lengths = []
    for z in range(nz):
        row = facies_grid[z, :]
        curr = 0
        for x in range(nx):
            if row[x] == target_facies_code:
                curr += 1
            else:
                if curr > 0:
                    run_lengths.append(curr * cell_width_m)
                    curr = 0
        if curr > 0:
            run_lengths.append(curr * cell_width_m)

    return float(np.mean(run_lengths)) if len(run_lengths) > 0 else 0.0


class AblationStudyEngine:
    """
    Runs the mandatory ablation experiment comparing Model 0, Model 1, Model 2, and Model 3.
    """

    def __init__(
        self,
        anchor_well_ids: List[str] = ["litholog2", "litholog9"],
        inspector: Optional[LithologInspector] = None,
        convention: str = DEFAULT_CONVENTION,
        random_seed: int = 42,
    ):
        self.anchor_ids = anchor_well_ids
        self.inspector = inspector or LithologInspector()
        self.convention = convention
        self.random_seed = int(random_seed)
        self.num_classes = len(CANONICAL_FACIES_SCHEMA)

        # Build calibrator to obtain empirical lengths and rates
        self.calibrator = HorizontalContinuityCalibrator(
            inspector=self.inspector,
            convention=self.convention,
            eligible_litholog_ids=self.anchor_ids,
        )
        self.pi = self.calibrator.get_stationary_proportions()
        self.lengths_df = self.calibrator.estimate_facies_lengths()

        # Extract true target bed thicknesses from anchor wells for ground truth benchmark
        self.true_beds = []
        self.true_proportions = np.zeros(self.num_classes, dtype=float)
        total_pts = 0
        for wid in self.anchor_ids:
            disc = self.inspector.discretize_litholog_1m(wid)
            aligned = align_to_common_datum(disc, litholog_id=wid, convention=self.convention, apply_source_orientation=True)
            seq = aligned["facies_code"].to_numpy()
            self.true_beds.extend(calculate_bed_thicknesses(seq))
            for c in range(self.num_classes):
                self.true_proportions[c] += np.sum(seq == c)
            total_pts += len(seq)
        if total_pts > 0:
            self.true_proportions = self.true_proportions / total_pts

        self.true_mean_bed_thickness = float(np.mean(self.true_beds)) if len(self.true_beds) > 0 else 3.76

        # Target true vertical transition matrix
        t_counts = np.zeros((self.num_classes, self.num_classes), dtype=float)
        for wid in self.anchor_ids:
            disc = self.inspector.discretize_litholog_1m(wid)
            aligned = align_to_common_datum(disc, litholog_id=wid, convention=self.convention, apply_source_orientation=True)
            seq = aligned["facies_code"].to_numpy()
            for k in range(len(seq) - 1):
                t_counts[int(seq[k]), int(seq[k+1])] += 1.0
        t_smoothed = t_counts + 0.01
        self.true_transition_matrix = t_smoothed / t_smoothed.sum(axis=1, keepdims=True)

    def run_ablation_realizations(
        self,
        n_realizations: int = 5,
        x_resolution_m: float = 25.0,
        z_resolution_m: float = 1.0,
    ) -> Dict[str, Any]:
        """
        Executes stochastic realizations for Models 0, 1, 2, 3 under identical grid and seeds.
        """
        # Set up classifiers for Model 2 and Model 3
        # Model 2: Prescribed Sprint I lengths
        clf_model2 = ConditionedMarkovClassifier(
            training_litholog_ids=self.anchor_ids,
            inspector=self.inspector,
            convention=self.convention,
            lateral_lengths=DEFAULT_LATERAL_FACIES_LENGTHS_M,
        )
        gen_model2 = InterWellRealizationGenerator(classifier=clf_model2, random_seed=self.random_seed)
        res_m2 = gen_model2.generate_2d_transect_realizations(
            anchor_well_ids=self.anchor_ids,
            n_realizations=n_realizations,
            x_resolution_m=x_resolution_m,
            z_resolution_m=z_resolution_m,
        )

        # Model 3: Empirically fitted Sprint J lengths
        lengths_B = {
            int(r["facies_code"]): float(r["L_fitted_m"])
            for _, r in self.lengths_df.iterrows()
        }
        clf_model3 = ConditionedMarkovClassifier(
            training_litholog_ids=self.anchor_ids,
            inspector=self.inspector,
            convention=self.convention,
            lateral_lengths=lengths_B,
        )
        gen_model3 = InterWellRealizationGenerator(classifier=clf_model3, random_seed=self.random_seed)
        res_m3 = gen_model3.generate_2d_transect_realizations(
            anchor_well_ids=self.anchor_ids,
            n_realizations=n_realizations,
            x_resolution_m=x_resolution_m,
            z_resolution_m=z_resolution_m,
        )

        nz, nx = res_m2["realizations"][0].shape
        x_grid = res_m2["x_coords_m"]
        z_grid = res_m2["z_coords_m"]
        well_cols = list(res_m2["well_grid_indices"].values())
        well_pos = res_m2["anchor_well_distances_m"]

        hard_mask = res_m2["hard_data_mask"]
        hard_ref = res_m2["realizations"][0]

        # Model 0: Stationary Facies Proportions Only
        rng0 = np.random.RandomState(self.random_seed)
        m0_realizations = np.zeros((n_realizations, nz, nx), dtype=int)
        for m in range(n_realizations):
            for z in range(nz):
                for x in range(nx):
                    if hard_mask[z, x]:
                        m0_realizations[m, z, x] = hard_ref[z, x]
                    else:
                        m0_realizations[m, z, x] = int(rng0.choice(self.num_classes, p=self.pi))

        # Model 1: Vertical Markov Succession Only (Unconditioned laterally)
        P_v = clf_model2.vertical_markov.transition_matrix_
        rng1 = np.random.RandomState(self.random_seed)
        m1_realizations = np.zeros((n_realizations, nz, nx), dtype=int)
        for m in range(n_realizations):
            for x in range(nx):
                s_prev = int(rng1.choice(self.num_classes, p=self.pi))
                m1_realizations[m, 0, x] = s_prev
                for z in range(1, nz):
                    p_row = P_v[s_prev, :]
                    s_next = int(rng1.choice(self.num_classes, p=p_row))
                    m1_realizations[m, z, x] = s_next
                    s_prev = s_next
                for z in range(nz):
                    if hard_mask[z, x]:
                        m1_realizations[m, z, x] = hard_ref[z, x]

        all_models = {
            "Model_0_Stationary_Proportions": {
                "realizations": m0_realizations,
                "label": "Model 0: Stationary Proportions",
                "is_hard_conditioned": True,
            },
            "Model_1_Vertical_Markov_Only": {
                "realizations": m1_realizations,
                "label": "Model 1: Vertical Markov Only",
                "is_hard_conditioned": True,
            },
            "Model_2_Vertical_Plus_Sprint_I": {
                "realizations": res_m2["realizations"],
                "label": "Model 2: Vertical + Sprint I Horizontal",
                "is_hard_conditioned": True,
            },
            "Model_3_Vertical_Plus_Sprint_J": {
                "realizations": res_m3["realizations"],
                "label": "Model 3: Vertical + Sprint J Calibrated",
                "is_hard_conditioned": True,
            },
        }

        # Compute evaluation metrics for each model across all realizations
        records = []
        for m_key, m_info in all_models.items():
            reals = m_info["realizations"]
            m_props_list = []
            m_tv_list = []
            m_trans_div_list = []
            m_bed_thick_list = []
            m_bed_err_list = []
            m_wass_list = []
            m_conn_list = []
            m_entropy_list = []
            m_honor_list = []

            for m in range(n_realizations):
                grid = reals[m]

                # Proportions
                p_sim = np.zeros(self.num_classes, dtype=float)
                for c in range(self.num_classes):
                    p_sim[c] = np.mean(grid == c)
                tv = 0.5 * np.sum(np.abs(p_sim - self.true_proportions))
                m_tv_list.append(tv)

                # Vertical transition matrix
                P_sim = compute_vertical_transition_matrix(grid, num_classes=self.num_classes)
                frob_div = float(np.linalg.norm(P_sim - self.true_transition_matrix))
                m_trans_div_list.append(frob_div)

                # Bed thicknesses
                sim_beds = []
                for x in range(nx):
                    sim_beds.extend(calculate_bed_thicknesses(grid[:, x], cell_height_m=z_resolution_m))
                mean_t = float(np.mean(sim_beds)) if len(sim_beds) > 0 else 0.0
                t_err = abs(mean_t - self.true_mean_bed_thickness)
                m_bed_thick_list.append(mean_t)
                m_bed_err_list.append(t_err)

                # Wasserstein distance between bed thickness distributions
                wass = float(wasserstein_distance(sim_beds, self.true_beds)) if len(sim_beds) > 0 else 0.0
                m_wass_list.append(wass)

                # Channel sand lateral connectivity
                conn = compute_lateral_connectivity_length(grid, target_facies_code=0, cell_width_m=x_resolution_m)
                m_conn_list.append(conn)

                # Spatial Shannon entropy
                probs_grid = np.zeros((self.num_classes, nz, nx), dtype=float)
                for c in range(self.num_classes):
                    probs_grid[c] = (grid == c).astype(float)
                ent_val = float(np.mean(scipy_entropy(p_sim)))
                m_entropy_list.append(ent_val)

                # Hard data honor rate
                if np.sum(hard_mask) > 0:
                    match_count = int(np.sum(grid[hard_mask] == hard_ref[hard_mask]))
                    m_honor_list.append(match_count / float(np.sum(hard_mask)))
                else:
                    m_honor_list.append(1.0)

            records.append({
                "model_key": m_key,
                "model_label": m_info["label"],
                "proportion_tv_distance_mean": round(float(np.mean(m_tv_list)), 4),
                "proportion_tv_distance_std": round(float(np.std(m_tv_list)), 4),
                "transition_matrix_div_mean": round(float(np.mean(m_trans_div_list)), 4),
                "transition_matrix_div_std": round(float(np.std(m_trans_div_list)), 4),
                "mean_bed_thickness_m": round(float(np.mean(m_bed_thick_list)), 2),
                "bed_thickness_error_m": round(float(np.mean(m_bed_err_list)), 2),
                "bed_thickness_error_std": round(float(np.std(m_bed_err_list)), 2),
                "bed_distribution_wasserstein_m": round(float(np.mean(m_wass_list)), 2),
                "channel_sand_connectivity_m": round(float(np.mean(m_conn_list)), 1),
                "shannon_entropy_mean": round(float(np.mean(m_entropy_list)), 4),
                "hard_data_honor_rate": round(float(np.mean(m_honor_list)), 4),
            })

        ablation_df = pd.DataFrame(records)

        # Central Test: Improvement over Vertical Markov Only (Model 1)
        m1_row = ablation_df[ablation_df["model_key"] == "Model_1_Vertical_Markov_Only"].iloc[0]
        improvements = []
        for _, r in ablation_df.iterrows():
            delta_tv = float(r["proportion_tv_distance_mean"] - m1_row["proportion_tv_distance_mean"])
            delta_bed_err = float(r["bed_thickness_error_m"] - m1_row["bed_thickness_error_m"])
            delta_div = float(r["transition_matrix_div_mean"] - m1_row["transition_matrix_div_mean"])
            delta_conn = float(r["channel_sand_connectivity_m"] - m1_row["channel_sand_connectivity_m"])

            improvements.append({
                "model_key": r["model_key"],
                "model_label": r["model_label"],
                "delta_tv_distance_vs_m1": round(delta_tv, 4),
                "delta_bed_thickness_error_m": round(delta_bed_err, 2),
                "delta_transition_divergence": round(delta_div, 4),
                "delta_connectivity_m": round(delta_conn, 1),
                "meaningful_improvement": bool(delta_bed_err < -0.2 or delta_tv < -0.02),
            })

        imp_df = pd.DataFrame(improvements)

        return {
            "ablation_summary": ablation_df,
            "improvement_vs_vertical_markov": imp_df,
            "grid_metadata": {
                "x_coords_m": x_grid,
                "z_coords_m": z_grid,
                "well_positions_x_m": well_pos,
                "well_cols": well_cols,
                "true_mean_bed_thickness_m": self.true_mean_bed_thickness,
            },
            "models_realizations": all_models,
        }
