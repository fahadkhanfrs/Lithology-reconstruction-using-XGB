"""
SMALT Provisional Spatial Prototype & Benchmarking Baseline.

Implements:
1. Strict Leave-One-Litholog-Out (LOLO) spatial validation on eligible wells (L2 - L12).
2. Explicit exclusion of Litholog 1 (missing coordinates).
3. Strictly leak-free spatial features: [X, Y, Z_rel] with configurable vertical datum reference.
4. Naive baselines: Training Prior Facies & Nearest-Well Vertical Profile.
5. Spatial Distance-Weighted K-Nearest Neighbors (Spatial 3D KNN) classifier.
6. Per-fold metrics: Macro-F1, Balanced Accuracy, Accuracy, Net-to-Gross error, Confusion Matrix.
7. Directional cross-group validation (Upstream <-> Downstream).
8. Datum offset sensitivity analysis.
"""

from typing import Dict, Any, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import (
    f1_score,
    balanced_accuracy_score,
    accuracy_score,
    confusion_matrix,
    precision_recall_fscore_support,
)

from smalt.descriptive.analyzer import LithologInspector, CANONICAL_FACIES_SCHEMA
from smalt.spatial.coordinates import load_source_coordinates, get_spatial_litholog_subset


class ProvisionalSpatialValidator:
    """
    Executes leak-free provisional spatial classification baselines on
    subsurface lithologs using verified source coordinates.
    """

    def __init__(
        self,
        inspector: Optional[LithologInspector] = None,
        vertical_reference: str = "common_zero",
        vertical_weight: float = 10.0,
    ):
        """
        Args:
            inspector: LithologInspector instance.
            vertical_reference: 'common_zero' (standard project datum: z = -depth),
                                'relative_to_top' (synonym for common_zero), or
                                'relative_to_base' (legacy local datum: z = max(depth) - depth).
            vertical_weight: Scaling factor for vertical distance vs horizontal distance
                             (reflects anisotropic geological correlation where vertical
                             scale ~10-100x finer than lateral scale).
        """
        self.inspector = inspector or LithologInspector()
        self.vertical_reference = vertical_reference
        self.vertical_weight = vertical_weight
        self.coordinates_df = load_source_coordinates()
        self.coords_by_id = {
            r["litholog_id"]: (r["x_m"], r["y_m"])
            for _, r in self.coordinates_df.iterrows()
            if r["coordinates_available"]
        }

    def prepare_discretized_spatial_points(
        self,
        litholog_id: str,
        datum_offset_m: float = 0.0,
    ) -> pd.DataFrame:
        """
        Extracts 1m discretized points for a litholog and assigns spatial [X, Y, Z_rel].

        Raises ValueError if litholog_id lacks source coordinates (e.g. litholog1).
        """
        if litholog_id not in self.coords_by_id:
            raise ValueError(
                f"Litholog '{litholog_id}' has no source coordinates available. "
                "Strictly ineligible for spatial modeling."
            )

        x_m, y_m = self.coords_by_id[litholog_id]
        disc_df = self.inspector.discretize_litholog_1m(litholog_id)
        
        # Calculate vertical coordinate relative to chosen datum
        if self.vertical_reference in ("common_zero", "common_zero_datum", "stratigraphic_height"):
            from smalt.spatial.datum import align_to_common_datum
            aligned = align_to_common_datum(disc_df, litholog_id=litholog_id, apply_source_orientation=True)
            z_rel = aligned["z_common_m"] + datum_offset_m
        elif self.vertical_reference == "legacy_inverted_elevation":
            z_rel = (-disc_df["depth_m"]) + datum_offset_m
        elif self.vertical_reference == "relative_to_top":
            z_rel = (-disc_df["depth_m"]) + datum_offset_m
        elif self.vertical_reference == "relative_to_base":
            max_d = disc_df["depth_m"].max()
            z_rel = (max_d - disc_df["depth_m"]) + datum_offset_m
        else:
            raise ValueError(f"Unknown vertical_reference: {self.vertical_reference}")

        df_pts = disc_df.copy()
        df_pts["x_m"] = x_m
        df_pts["y_m"] = y_m
        df_pts["z_rel_m"] = z_rel

        return df_pts

    def run_spatial_lolo(
        self,
        eligible_log_ids: Optional[List[str]] = None,
        k_neighbors: int = 5,
        datum_offsets: Optional[Dict[str, float]] = None,
    ) -> Dict[str, Any]:
        """
        Executes Leave-One-Litholog-Out (LOLO) spatial cross-validation.

        For each fold:
          - Train: All points from other eligible wells
          - Test: All points from held-out target well
          - Preprocessing (StandardScaler) is strictly fit on training points only.
        """
        if eligible_log_ids is None:
            # Default: all eligible wells with coordinates (L2 - L12)
            eligible_log_ids = [lid for lid in self.coords_by_id.keys()]

        # Check for L1 exclusion
        if "litholog1" in eligible_log_ids:
            raise ValueError(
                "Litholog 1 is missing physical coordinates and must be excluded from spatial folds."
            )

        offsets = datum_offsets or {}

        # Pre-extract data points for all candidate wells
        well_data = {}
        for lid in eligible_log_ids:
            well_data[lid] = self.prepare_discretized_spatial_points(
                lid, datum_offset_m=offsets.get(lid, 0.0)
            )

        fold_results = []
        all_y_true = []
        all_y_pred_knn = []
        all_y_pred_prior = []
        all_y_pred_near_well = []

        for target_id in eligible_log_ids:
            train_ids = [lid for lid in eligible_log_ids if lid != target_id]
            
            # Combine training points
            train_dfs = [well_data[lid] for lid in train_ids]
            df_train = pd.concat(train_dfs, ignore_index=True)
            df_test = well_data[target_id].copy()

            # Coordinates features (anisotropic weighting: horizontal vs vertical)
            X_train_raw = df_train[["x_m", "y_m", "z_rel_m"]].to_numpy()
            X_train_scaled = X_train_raw.copy()
            X_train_scaled[:, 2] *= self.vertical_weight

            X_test_raw = df_test[["x_m", "y_m", "z_rel_m"]].to_numpy()
            X_test_scaled = X_test_raw.copy()
            X_test_scaled[:, 2] *= self.vertical_weight

            scaler = StandardScaler()
            X_train_norm = scaler.fit_transform(X_train_scaled)
            X_test_norm = scaler.transform(X_test_scaled)

            y_train = df_train["facies_code"].to_numpy()
            y_test = df_test["facies_code"].to_numpy()

            # -------------------------------------------------------------
            # Baseline 1: Training Set Prior (Majority Class)
            # -------------------------------------------------------------
            classes, counts = np.unique(y_train, return_counts=True)
            majority_class = classes[np.argmax(counts)]
            y_pred_prior = np.full_like(y_test, fill_value=majority_class)

            # -------------------------------------------------------------
            # Baseline 2: Nearest Spatial Well Profile
            # -------------------------------------------------------------
            target_x, target_y = self.coords_by_id[target_id]
            dists = [
                (np.sqrt((self.coords_by_id[lid][0] - target_x)**2 + (self.coords_by_id[lid][1] - target_y)**2), lid)
                for lid in train_ids
            ]
            nearest_well_id = min(dists, key=lambda x: x[0])[1]
            df_near = well_data[nearest_well_id]
            
            # For each test z_rel, find closest z_rel in nearest well
            near_z = df_near["z_rel_m"].to_numpy()
            near_y = df_near["facies_code"].to_numpy()
            y_pred_near_well = []
            for z_val in df_test["z_rel_m"].to_numpy():
                idx_closest = int(np.argmin(np.abs(near_z - z_val)))
                y_pred_near_well.append(near_y[idx_closest])
            y_pred_near_well = np.array(y_pred_near_well)

            # -------------------------------------------------------------
            # Model 3: Spatial 3D KNN (Inverse-Distance Weighted)
            # -------------------------------------------------------------
            knn = KNeighborsClassifier(n_neighbors=k_neighbors, weights="distance")
            knn.fit(X_train_norm, y_train)
            y_pred_knn = knn.predict(X_test_norm)

            # Calculate metrics for this fold
            k_labels = list(range(len(CANONICAL_FACIES_SCHEMA)))
            
            # Class support in test fold
            support_dict = {
                CANONICAL_FACIES_SCHEMA[k]["canonical_name"]: int(np.sum(y_test == meta["code"]))
                for k, meta in CANONICAL_FACIES_SCHEMA.items()
            }

            # Explicit Net-to-Gross definitions (6-state schema)
            # Pure / Net sand: Channel Sandstone (0) + Planar Sandstone (1)
            true_net_sand_ntg = float(np.mean(np.isin(y_test, [0, 1])))
            pred_net_sand_ntg_knn = float(np.mean(np.isin(y_pred_knn, [0, 1])))
            pred_net_sand_ntg_near = float(np.mean(np.isin(y_pred_near_well, [0, 1])))
            pred_net_sand_ntg_prior = float(np.mean(np.isin(y_pred_prior, [0, 1])))

            # Channel sandstone NTG
            true_channel_ntg = float(np.mean(y_test == 0))
            pred_channel_ntg_knn = float(np.mean(y_pred_knn == 0))

            fold_record = {
                "target_litholog_id": target_id,
                "nearest_well_id": nearest_well_id,
                "test_points_count": len(y_test),
                "true_sand_ntg": round(true_net_sand_ntg, 4),
                "true_channel_ntg": round(true_channel_ntg, 4),
                # KNN Metrics
                "knn_accuracy": round(accuracy_score(y_test, y_pred_knn), 4),
                "knn_balanced_acc": round(balanced_accuracy_score(y_test, y_pred_knn), 4),
                "knn_macro_f1": round(f1_score(y_test, y_pred_knn, average="macro", zero_division=0), 4),
                "knn_sand_ntg_error": round(abs(pred_net_sand_ntg_knn - true_net_sand_ntg), 4),
                "knn_channel_ntg_error": round(abs(pred_channel_ntg_knn - true_channel_ntg), 4),
                # Nearest Well Metrics
                "near_well_accuracy": round(accuracy_score(y_test, y_pred_near_well), 4),
                "near_well_balanced_acc": round(balanced_accuracy_score(y_test, y_pred_near_well), 4),
                "near_well_macro_f1": round(f1_score(y_test, y_pred_near_well, average="macro", zero_division=0), 4),
                "near_well_sand_ntg_error": round(abs(pred_net_sand_ntg_near - true_net_sand_ntg), 4),
                # Prior Facies Metrics
                "prior_accuracy": round(accuracy_score(y_test, y_pred_prior), 4),
                "prior_balanced_acc": round(balanced_accuracy_score(y_test, y_pred_prior), 4),
                "prior_macro_f1": round(f1_score(y_test, y_pred_prior, average="macro", zero_division=0), 4),
                "prior_sand_ntg_error": round(abs(pred_net_sand_ntg_prior - true_net_sand_ntg), 4),
                "support_by_class": support_dict,
            }
            fold_results.append(fold_record)

            all_y_true.extend(y_test)
            all_y_pred_knn.extend(y_pred_knn)
            all_y_pred_near_well.extend(y_pred_near_well)
            all_y_pred_prior.extend(y_pred_prior)

        all_y_true = np.array(all_y_true)
        all_y_pred_knn = np.array(all_y_pred_knn)
        all_y_pred_near_well = np.array(all_y_pred_near_well)
        all_y_pred_prior = np.array(all_y_pred_prior)

        k_labels = list(range(len(CANONICAL_FACIES_SCHEMA)))

        def _calc_per_class(y_t, y_p):
            p, r, f, s = precision_recall_fscore_support(y_t, y_p, labels=k_labels, zero_division=0)
            res = {}
            for meta in CANONICAL_FACIES_SCHEMA.values():
                c = meta["code"]
                res[meta["canonical_name"]] = {
                    "code": c,
                    "precision": round(float(p[c]), 4),
                    "recall": round(float(r[c]), 4),
                    "f1": round(float(f[c]), 4),
                    "support": int(s[c]),
                }
            return res

        # Overall aggregate summary
        summary = {
            "num_evaluated_wells": len(eligible_log_ids),
            "total_evaluated_points": len(all_y_true),
            "vertical_reference": self.vertical_reference,
            "vertical_weight": self.vertical_weight,
            "overall_knn_macro_f1": round(f1_score(all_y_true, all_y_pred_knn, average="macro", zero_division=0), 4),
            "overall_knn_balanced_acc": round(balanced_accuracy_score(all_y_true, all_y_pred_knn), 4),
            "overall_knn_accuracy": round(accuracy_score(all_y_true, all_y_pred_knn), 4),
            "overall_near_well_macro_f1": round(f1_score(all_y_true, all_y_pred_near_well, average="macro", zero_division=0), 4),
            "overall_near_well_balanced_acc": round(balanced_accuracy_score(all_y_true, all_y_pred_near_well), 4),
            "overall_near_well_accuracy": round(accuracy_score(all_y_true, all_y_pred_near_well), 4),
            "overall_prior_macro_f1": round(f1_score(all_y_true, all_y_pred_prior, average="macro", zero_division=0), 4),
            "overall_prior_balanced_acc": round(balanced_accuracy_score(all_y_true, all_y_pred_prior), 4),
            "overall_prior_accuracy": round(accuracy_score(all_y_true, all_y_pred_prior), 4),
            "knn_per_class": _calc_per_class(all_y_true, all_y_pred_knn),
            "near_well_per_class": _calc_per_class(all_y_true, all_y_pred_near_well),
            "prior_per_class": _calc_per_class(all_y_true, all_y_pred_prior),
            "knn_confusion_matrix": confusion_matrix(all_y_true, all_y_pred_knn, labels=k_labels).tolist(),
            "near_well_confusion_matrix": confusion_matrix(all_y_true, all_y_pred_near_well, labels=k_labels).tolist(),
            "prior_confusion_matrix": confusion_matrix(all_y_true, all_y_pred_prior, labels=k_labels).tolist(),
        }

        df_folds = pd.DataFrame(fold_results)
        return {
            "per_fold_table": df_folds,
            "aggregate_summary": summary,
        }

    def run_directional_spatial_experiment(
        self,
        training_group_ids: List[str],
        target_group_ids: List[str],
        experiment_name: str,
        k_neighbors: int = 5,
    ) -> Dict[str, Any]:
        """
        Executes a cross-group spatial prediction experiment (e.g. Upstream -> Downstream).
        """
        for lid in training_group_ids + target_group_ids:
            if lid not in self.coords_by_id:
                raise ValueError(f"Litholog '{lid}' lacks coordinates.")

        train_dfs = [self.prepare_discretized_spatial_points(lid) for lid in training_group_ids]
        df_train = pd.concat(train_dfs, ignore_index=True)

        test_dfs = [self.prepare_discretized_spatial_points(lid) for lid in target_group_ids]
        df_test = pd.concat(test_dfs, ignore_index=True)

        X_train_scaled = df_train[["x_m", "y_m", "z_rel_m"]].to_numpy().copy()
        X_train_scaled[:, 2] *= self.vertical_weight

        X_test_scaled = df_test[["x_m", "y_m", "z_rel_m"]].to_numpy().copy()
        X_test_scaled[:, 2] *= self.vertical_weight

        scaler = StandardScaler()
        X_train_norm = scaler.fit_transform(X_train_scaled)
        X_test_norm = scaler.transform(X_test_scaled)

        y_train = df_train["facies_code"].to_numpy()
        y_test = df_test["facies_code"].to_numpy()

        knn = KNeighborsClassifier(n_neighbors=k_neighbors, weights="distance")
        knn.fit(X_train_norm, y_train)
        y_pred = knn.predict(X_test_norm)

        # Baseline prior
        classes, counts = np.unique(y_train, return_counts=True)
        majority_class = classes[np.argmax(counts)]
        y_pred_prior = np.full_like(y_test, fill_value=majority_class)

        k_labels = list(range(len(CANONICAL_FACIES_SCHEMA)))

        return {
            "experiment_name": experiment_name,
            "training_wells_count": len(training_group_ids),
            "target_wells_count": len(target_group_ids),
            "total_test_points": len(y_test),
            "knn_macro_f1": round(f1_score(y_test, y_pred, average="macro", zero_division=0), 4),
            "knn_balanced_accuracy": round(balanced_accuracy_score(y_test, y_pred), 4),
            "knn_accuracy": round(accuracy_score(y_test, y_pred), 4),
            "prior_macro_f1": round(f1_score(y_test, y_pred_prior, average="macro", zero_division=0), 4),
            "prior_balanced_accuracy": round(balanced_accuracy_score(y_test, y_pred_prior), 4),
            "prior_accuracy": round(accuracy_score(y_test, y_pred_prior), 4),
            "confusion_matrix": confusion_matrix(y_test, y_pred, labels=k_labels).tolist(),
        }
