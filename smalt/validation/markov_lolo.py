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
        p_init_train: Optional[np.ndarray] = None,
        embedded: bool = False,
        use_discretized: bool = True,
    ) -> Dict[str, Any]:
        """
        Evaluates a held-out target litholog sequence under the training Markov model.

        The target facies sequence is NEVER used to adjust P_train, pi_train, or p_init_train.

        Scores:
        - initial_state_log_score: ln P_init(s_1)
        - total_transition_log_score: sum_{t=1}^{M} ln P(s_t, s_{t+1})
        - mean_transition_log_score: (1/M) sum_{t=1}^{M} ln P(s_t, s_{t+1})
        - transition_perplexity: exp(-mean_transition_log_score)
        - total_sequence_log_score: initial_state_log_score + total_transition_log_score

        Gap sensitivity:
        - gap_crossing_transitions_count: number of transitions touching unrecorded gaps
        - mean_transition_log_score_gap_masked: transition log score excluding gap transitions
        - transition_perplexity_gap_masked: perplexity excluding gap transitions
        """
        if embedded and not use_discretized:
            # Build distinct bed sequence from continuous intervals
            beds_df = self.inspector.build_embedded_bed_sequence(target_log_id)
            seq = beds_df["facies_code"].to_numpy()
            T = len(seq)
            if T < 2:
                raise ValueError(f"Target bed sequence for {target_log_id} has fewer than 2 distinct beds.")

            transitions = list(zip(seq[:-1], seq[1:]))
            n_trans = len(transitions)
            crosses_gap_flags = beds_df["crosses_gap"].iloc[1:].to_numpy() if T > 1 else np.array([])
            gap_crossing_count = int(np.sum(crosses_gap_flags)) if len(crosses_gap_flags) > 0 else 0

            # Initial state (basal bed) distribution
            initial_code = int(seq[0])
            initial_facies = str(beds_df.iloc[0]["facies"])
            if p_init_train is not None:
                initial_prob = max(float(p_init_train[initial_code]), 1e-12)
            else:
                initial_prob = max(float(pi_train[initial_code]), 1e-12)
            initial_log_score = float(np.log(initial_prob))

            # Transitions
            trans_log_probs = [float(np.log(max(float(P_train[u, v]), 1e-12))) for u, v in transitions]
            total_trans_log_score = float(np.sum(trans_log_probs)) if n_trans > 0 else 0.0
            mean_trans_log_score = float(np.mean(trans_log_probs)) if n_trans > 0 else 0.0
            trans_perplexity = float(np.exp(-mean_trans_log_score)) if n_trans > 0 else 1.0

            # Gap-masked sensitivity
            if gap_crossing_count > 0 and (n_trans - gap_crossing_count) > 0:
                masked_probs = [lp for lp, gf in zip(trans_log_probs, crosses_gap_flags) if not gf]
                mean_trans_gap_masked = float(np.mean(masked_probs))
                perp_gap_masked = float(np.exp(-mean_trans_gap_masked))
            else:
                mean_trans_gap_masked = mean_trans_log_score
                perp_gap_masked = trans_perplexity

            total_seq_score = float(initial_log_score + total_trans_log_score)

        else:
            # Regular discretized 1m grid sequence (retains self-transitions)
            df = self.inspector.discretize_litholog_1m(target_log_id)
            df_sorted = df.sort_values(by="depth_m", ascending=False).reset_index(drop=True)
            seq = df_sorted["facies_code"].to_numpy()
            is_gap_filled = df_sorted["is_gap_filled"].to_numpy()
            T = len(seq)
            if T < 2:
                raise ValueError(f"Target discretized sequence for {target_log_id} has fewer than 2 observations.")

            transitions = list(zip(seq[:-1], seq[1:]))
            n_trans = len(transitions)

            gap_flags = [bool(is_gap_filled[t] or is_gap_filled[t + 1]) for t in range(T - 1)]
            gap_crossing_count = int(sum(gap_flags))

            # Initial state (basal grid cell)
            initial_code = int(seq[0])
            initial_facies = str(df_sorted.iloc[0]["facies"])
            if p_init_train is not None:
                initial_prob = max(float(p_init_train[initial_code]), 1e-12)
            else:
                initial_prob = max(float(pi_train[initial_code]), 1e-12)
            initial_log_score = float(np.log(initial_prob))

            # Transitions
            trans_log_probs = [float(np.log(max(float(P_train[u, v]), 1e-12))) for u, v in transitions]
            total_trans_log_score = float(np.sum(trans_log_probs)) if n_trans > 0 else 0.0
            mean_trans_log_score = float(np.mean(trans_log_probs)) if n_trans > 0 else 0.0
            trans_perplexity = float(np.exp(-mean_trans_log_score)) if n_trans > 0 else 1.0

            # Gap-masked sensitivity
            if gap_crossing_count > 0 and (n_trans - gap_crossing_count) > 0:
                masked_probs = [lp for lp, gf in zip(trans_log_probs, gap_flags) if not gf]
                mean_trans_gap_masked = float(np.mean(masked_probs))
                perp_gap_masked = float(np.exp(-mean_trans_gap_masked))
            else:
                mean_trans_gap_masked = mean_trans_log_score
                perp_gap_masked = trans_perplexity

            total_seq_score = float(initial_log_score + total_trans_log_score)

        # Target Facies Proportions
        target_counts = np.zeros(self.K, dtype=np.float64)
        for c in seq:
            target_counts[c] += 1.0
        target_proportions = target_counts / target_counts.sum()

        # Stationary Divergence (Diagnostic only, not classification score)
        tv_distance = 0.5 * float(np.sum(np.abs(pi_train - target_proportions)))
        l1_distance = float(np.sum(np.abs(pi_train - target_proportions)))
        rmse_stationary = float(np.sqrt(np.mean((pi_train - target_proportions) ** 2)))

        # Target Empirical Transition Matrix (strictly for divergence comparison)
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
            "initial_state_facies": initial_facies,
            "initial_state_code": initial_code,
            "initial_state_prob": round(initial_prob, 5),
            "initial_state_log_score": round(initial_log_score, 4),
            "total_transition_log_score": round(total_trans_log_score, 4),
            "mean_transition_log_score": round(mean_trans_log_score, 4),
            "transition_perplexity": round(trans_perplexity, 4),
            "total_sequence_log_score": round(total_seq_score, 4),
            "gap_crossing_transitions_count": gap_crossing_count,
            "mean_transition_log_score_gap_masked": round(mean_trans_gap_masked, 4),
            "transition_perplexity_gap_masked": round(perp_gap_masked, 4),
            "sequence_log_likelihood": round(total_seq_score, 4),  # Backward compatibility
            "mean_transition_log_likelihood": round(mean_trans_log_score, 4),  # Backward compatibility
            "sequence_perplexity": round(trans_perplexity, 4),  # Backward compatibility
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
            p_init_train = self.inspector.compute_initial_state_distribution(
                litholog_ids=training_ids,
                embedded=embedded,
                use_discretized=use_discretized,
                smoothing_alpha=self.smoothing_alpha,
            )
            eval_metrics = self.evaluate_target_sequence(
                target_log_id=target_id,
                P_train=P_train,
                pi_train=pi_train,
                p_init_train=p_init_train,
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
        p_init_train = self.inspector.compute_initial_state_distribution(
            litholog_ids=training_group_ids,
            embedded=embedded,
            use_discretized=use_discretized,
            smoothing_alpha=self.smoothing_alpha,
        )

        # Evaluate each target well individually
        target_evals = []
        for tid in target_group_ids:
            res = self.evaluate_target_sequence(
                target_log_id=tid,
                P_train=P_train,
                pi_train=pi_train,
                p_init_train=p_init_train,
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
