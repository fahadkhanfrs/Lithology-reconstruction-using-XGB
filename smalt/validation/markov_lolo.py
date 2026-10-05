"""
SMALT Leakage-Free 1D Stratigraphic Markov Validation Engine (Sprint C).

Implements strict Leave-One-Litholog-Out (LOLO) cross-validation and directional
(upstream <-> downstream) vertical succession benchmarking.

Mathematical and Scientific Invariants:
1. Complete holdout: facies observations from the target litholog are strictly
   excluded from all stages of model fitting, transition counting, smoothing,
   and stationary distribution estimation.
2. Directionality: transitions are tallied strictly in upward stratigraphic order
   (decreasing depth_m), matching depositional succession.
3. Metric safety: zero probabilities are regularized via Laplace smoothing prior (alpha).
4. Diagnostic clarity: stationary proportion differences are reported as succession
   diagnostics, NOT spatial predictions. Ordinary classification accuracy (which
   requires a point-wise predictor) is not claimed.
"""

from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import pandas as pd

from smalt.descriptive.analyzer import (
    LithologInspector,
    CANONICAL_FACIES_SCHEMA,
    PROVENANCE_METADATA,
)


class MarkovLOLOValidator:
    """
    Executes leakage-free Leave-One-Litholog-Out (LOLO) and directional
    Markov chain validation across digitized stratigraphic sequences.
    """

    def __init__(
        self,
        inspector: Optional[LithologInspector] = None,
        smoothing_alpha: float = 0.1,
    ):
        self.inspector = inspector or LithologInspector()
        self.smoothing_alpha = float(smoothing_alpha)
        self.K = len(CANONICAL_FACIES_SCHEMA)

    def fit_training_markov(
        self,
        training_log_ids: List[str],
        embedded: bool = False,
        use_discretized: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Fits transition probabilities and stationary distribution strictly on training logs.

        Returns:
            P_train: K x K row-stochastic transition matrix.
            N_train: K x K raw transition count matrix.
            pi_train: 1 x K invariant stationary distribution vector.
        """
        return self.inspector.compute_vertical_transitions(
            litholog_ids=training_log_ids,
            embedded=embedded,
            use_discretized=use_discretized,
            smoothing_alpha=self.smoothing_alpha,
        )

    def evaluate_target_sequence(
        self,
        target_log_id: str,
        P_train: np.ndarray,
        pi_train: np.ndarray,
        embedded: bool = False,
        use_discretized: bool = True,
    ) -> Dict[str, Any]:
        """
        Evaluates a held-out target litholog sequence under the training Markov model.

        The target facies sequence is NEVER used to adjust P_train or pi_train.
        """
        if use_discretized:
            df = self.inspector.discretize_litholog_1m(target_log_id)
            df_sorted = df.sort_values(by="depth_m", ascending=False)
            seq = df_sorted["facies_code"].to_numpy()
        else:
            df = self.inspector.load_raw_litholog(target_log_id)
            df_sorted = df.sort_values(by="Bottom", ascending=False)
            seq = df_sorted["Facies"].map(lambda f: CANONICAL_FACIES_SCHEMA[f]["code"]).to_numpy()

        T = len(seq)
        if T < 2:
            raise ValueError(f"Target sequence for {target_log_id} has fewer than 2 observations.")

        # Filter out self-transitions if embedded
        transitions = []
        for u, v in zip(seq[:-1], seq[1:]):
            if embedded and u == v:
                continue
            transitions.append((int(u), int(v)))

        n_trans = len(transitions)

        # 1. Sequence Log-Likelihood
        # Initial state likelihood + product of transition probabilities
        initial_prob = max(float(pi_train[seq[0]]), 1e-12)
        log_initial = np.log(initial_prob)

        trans_log_probs = []
        for u, v in transitions:
            p_uv = max(float(P_train[u, v]), 1e-12)
            trans_log_probs.append(np.log(p_uv))

        total_log_likelihood = float(log_initial + sum(trans_log_probs))
        mean_trans_log_likelihood = float(np.mean(trans_log_probs)) if n_trans > 0 else 0.0
        perplexity = float(np.exp(-mean_trans_log_likelihood)) if n_trans > 0 else 1.0

        # 2. Target Facies Proportions
        target_counts = np.zeros(self.K, dtype=np.float64)
        for c in seq:
            target_counts[c] += 1.0
        target_proportions = target_counts / target_counts.sum()

        # 3. Stationary Divergence
        tv_distance = 0.5 * float(np.sum(np.abs(pi_train - target_proportions)))
        l1_distance = float(np.sum(np.abs(pi_train - target_proportions)))
        rmse_stationary = float(np.sqrt(np.mean((pi_train - target_proportions) ** 2)))

        # 4. Target Empirical Transition Matrix (for divergence comparison only)
        P_target, _, _ = self.inspector.compute_vertical_transitions(
            litholog_ids=[target_log_id],
            embedded=embedded,
            use_discretized=use_discretized,
            smoothing_alpha=self.smoothing_alpha,
        )
        frobenius_divergence = float(np.linalg.norm(P_train - P_target, "fro"))

        prov = PROVENANCE_METADATA.get(target_log_id, {})

        return {
            "target_log_id": target_log_id,
            "provenance_category": prov.get("provenance_category", "unknown"),
            "coordinates_available": prov.get("coordinates_available", False),
            "group": prov.get("group", "unassigned"),
            "embedded": embedded,
            "use_discretized": use_discretized,
            "sequence_length": T,
            "evaluated_transitions": n_trans,
            "sequence_log_likelihood": round(total_log_likelihood, 4),
            "mean_transition_log_likelihood": round(mean_trans_log_likelihood, 4),
            "sequence_perplexity": round(perplexity, 4),
            "stationary_tv_distance": round(tv_distance, 4),
            "stationary_l1_distance": round(l1_distance, 4),
            "stationary_rmse": round(rmse_stationary, 4),
            "matrix_frobenius_divergence": round(frobenius_divergence, 4),
            "target_proportions": {
                name: round(float(target_proportions[meta["code"]]), 4)
                for name, meta in CANONICAL_FACIES_SCHEMA.items()
            },
            "training_stationary_proportions": {
                name: round(float(pi_train[meta["code"]]), 4)
                for name, meta in CANONICAL_FACIES_SCHEMA.items()
            },
        }

    def run_lolo_cross_validation(
        self,
        eligible_log_ids: List[str],
        embedded: bool = False,
        use_discretized: bool = True,
    ) -> pd.DataFrame:
        """
        Executes Leave-One-Litholog-Out (LOLO) across all provided eligible lithologs.

        For each fold k in 1..N:
          - Train: {L_j : j != k}
          - Test: {L_k}
        """
        results = []
        for target_id in eligible_log_ids:
            training_ids = [lid for lid in eligible_log_ids if lid != target_id]
            P_train, _, pi_train = self.fit_training_markov(
                training_log_ids=training_ids,
                embedded=embedded,
                use_discretized=use_discretized,
            )
            eval_metrics = self.evaluate_target_sequence(
                target_log_id=target_id,
                P_train=P_train,
                pi_train=pi_train,
                embedded=embedded,
                use_discretized=use_discretized,
            )
            eval_metrics["num_training_logs"] = len(training_ids)
            eval_metrics["training_log_ids"] = ",".join(training_ids)
            results.append(eval_metrics)

        return pd.DataFrame(results)

    def run_directional_experiment(
        self,
        training_group_ids: List[str],
        target_group_ids: List[str],
        experiment_name: str,
        embedded: bool = False,
        use_discretized: bool = True,
    ) -> Dict[str, Any]:
        """
        Executes a directional validation experiment (e.g. Upstream -> Downstream).
        """
        # Fit strictly on training group
        P_train, N_train, pi_train = self.fit_training_markov(
            training_log_ids=training_group_ids,
            embedded=embedded,
            use_discretized=use_discretized,
        )

        # Evaluate each target well individually
        target_evals = []
        for tid in target_group_ids:
            res = self.evaluate_target_sequence(
                target_log_id=tid,
                P_train=P_train,
                pi_train=pi_train,
                embedded=embedded,
                use_discretized=use_discretized,
            )
            target_evals.append(res)

        # Fit target group transition matrix as composite empirical reference
        P_target_group, N_target_group, pi_target_group = self.inspector.compute_vertical_transitions(
            litholog_ids=target_group_ids,
            embedded=embedded,
            use_discretized=use_discretized,
            smoothing_alpha=self.smoothing_alpha,
        )

        group_frobenius = float(np.linalg.norm(P_train - P_target_group, "fro"))
        group_tv_distance = 0.5 * float(np.sum(np.abs(pi_train - pi_target_group)))

        return {
            "experiment_name": experiment_name,
            "embedded": embedded,
            "use_discretized": use_discretized,
            "training_logs": training_group_ids,
            "target_logs": target_group_ids,
            "n_training_logs": len(training_group_ids),
            "n_target_logs": len(target_group_ids),
            "composite_frobenius_divergence": round(group_frobenius, 4),
            "composite_stationary_tv_distance": round(group_tv_distance, 4),
            "training_stationary": {
                name: round(float(pi_train[meta["code"]]), 4)
                for name, meta in CANONICAL_FACIES_SCHEMA.items()
            },
            "target_composite_stationary": {
                name: round(float(pi_target_group[meta["code"]]), 4)
                for name, meta in CANONICAL_FACIES_SCHEMA.items()
            },
            "per_target_results": target_evals,
        }
