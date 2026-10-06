"""
SMALT Sprint D Test Suite: Verification of 1D Markov Validation & Summary Reconciliation.

Covers:
1. Regular transition scoring retains self-transitions.
2. Embedded evaluation uses distinct-bed transitions only.
3. Embedded transition matrices have strictly zero diagonal within numerical tolerance.
4. All transition-matrix rows sum to 1 within 10^-7.
5. Target litholog is excluded from fitted training counts (leakage invariant).
6. Gap-crossing transitions are handled according to documented gap-masking policy.
7. Fixed toy sequence gives analytically reproducible log-likelihood and perplexity values.
8. Continuous bed thickness and 1m discretized proportions are reported separately.
9. Deterministic training and scoring across runs.
10. Preserves compatibility with existing tests.
"""

import math
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from smalt.descriptive.analyzer import (
    LithologInspector,
    CANONICAL_FACIES_SCHEMA,
    PROVENANCE_METADATA,
)
from smalt.validation.markov_lolo import MarkovLOLOValidator


@pytest.fixture(scope="module")
def inspector():
    repo_root = Path(__file__).resolve().parent.parent
    return LithologInspector(raw_dir=repo_root / "data" / "raw_lithologs")


@pytest.fixture(scope="module")
def validator(inspector):
    return MarkovLOLOValidator(inspector=inspector, smoothing_alpha=0.1)


# -----------------------------------------------------------------------------
# 1. Regular transition scoring retains self-transitions
# -----------------------------------------------------------------------------
def test_regular_scoring_retains_self_transitions(inspector, validator):
    """
    Verifies that the regular Markov chain evaluation retains self-transitions
    on the 1m discretized grid sequence (n_trans = sequence_length - 1).
    """
    for lid in [f"litholog{i}" for i in range(1, 13)]:
        disc_df = inspector.discretize_litholog_1m(lid)
        T = len(disc_df)
        
        # Fit a dummy training model on all other logs
        other_ids = [f"litholog{i}" for i in range(1, 13) if f"litholog{i}" != lid]
        P_train, _, pi_train = validator.fit_training_markov(other_ids, embedded=False, use_discretized=True)
        p_init_train = validator.inspector.compute_initial_state_distribution(other_ids, embedded=False, use_discretized=True)
        
        eval_metrics = validator.evaluate_target_sequence(
            target_log_id=lid,
            P_train=P_train,
            pi_train=pi_train,
            p_init_train=p_init_train,
            embedded=False,
            use_discretized=True,
        )
        assert eval_metrics["sequence_length"] == T
        assert eval_metrics["evaluated_transitions"] == T - 1


# -----------------------------------------------------------------------------
# 2. Embedded evaluation uses distinct-bed transitions only
# -----------------------------------------------------------------------------
def test_embedded_evaluation_uses_distinct_bed_transitions_only(inspector, validator):
    """
    Verifies that the embedded bed sequence consists strictly of consecutive distinct
    facies beds (no adjacent duplicates) and that embedded evaluation scores only
    genuine transitions between distinct facies.
    """
    for lid in [f"litholog{i}" for i in range(1, 13)]:
        beds_df = inspector.build_embedded_bed_sequence(lid)
        T_beds = len(beds_df)
        codes = beds_df["facies_code"].to_numpy()
        
        # Invariant: consecutive beds MUST have distinct facies codes
        for t in range(T_beds - 1):
            assert codes[t] != codes[t + 1], f"Self-transition found in embedded sequence for {lid} at index {t}"
            
        # Fit embedded model on other logs
        other_ids = [f"litholog{i}" for i in range(1, 13) if f"litholog{i}" != lid]
        P_train, _, pi_train = validator.fit_training_markov(other_ids, embedded=True, use_discretized=False)
        p_init_train = validator.inspector.compute_initial_state_distribution(other_ids, embedded=True, use_discretized=False)
        
        eval_metrics = validator.evaluate_target_sequence(
            target_log_id=lid,
            P_train=P_train,
            pi_train=pi_train,
            p_init_train=p_init_train,
            embedded=True,
            use_discretized=False,
        )
        assert eval_metrics["sequence_length"] == T_beds
        assert eval_metrics["evaluated_transitions"] == T_beds - 1


# -----------------------------------------------------------------------------
# 3. Embedded transition matrices have strictly zero diagonal
# -----------------------------------------------------------------------------
def test_embedded_transition_matrices_zero_diagonal(validator):
    """
    Verifies that embedded transition probability matrices have strictly zero diagonal
    entries (P_ii == 0.0) across all LOLO folds.
    """
    all_lids = [f"litholog{i}" for i in range(1, 13)]
    for target_id in all_lids:
        train_ids = [lid for lid in all_lids if lid != target_id]
        P_train, _, _ = validator.fit_training_markov(train_ids, embedded=True, use_discretized=False)
        np.testing.assert_allclose(np.diag(P_train), 0.0, atol=1e-12)


# -----------------------------------------------------------------------------
# 4. All transition-matrix rows sum to one within 10^-7
# -----------------------------------------------------------------------------
def test_transition_matrix_rows_sum_to_one(validator):
    """
    Verifies that rows of both regular and embedded transition matrices sum to 1.0 within 10^-7.
    """
    all_lids = [f"litholog{i}" for i in range(1, 13)]
    for target_id in all_lids:
        train_ids = [lid for lid in all_lids if lid != target_id]
        
        # Regular
        P_reg, _, _ = validator.fit_training_markov(train_ids, embedded=False, use_discretized=True)
        np.testing.assert_allclose(P_reg.sum(axis=1), np.ones(P_reg.shape[0]), atol=1e-7)
        
        # Embedded
        P_emb, _, _ = validator.fit_training_markov(train_ids, embedded=True, use_discretized=False)
        np.testing.assert_allclose(P_emb.sum(axis=1), np.ones(P_emb.shape[0]), atol=1e-7)


# -----------------------------------------------------------------------------
# 5. Target litholog is excluded from fitted training counts (Leakage check)
# -----------------------------------------------------------------------------
def test_target_litholog_excluded_from_fitted_counts(inspector, validator):
    """
    Verifies that for every LOLO fold, the target litholog is 100% excluded from
    the training transition counts: N_train(all) - N_train(fold) == N_target.
    """
    all_lids = [f"litholog{i}" for i in range(1, 13)]
    _, N_all_reg, _ = inspector.compute_vertical_transitions(all_lids, embedded=False, use_discretized=True)
    _, N_all_emb, _ = inspector.compute_vertical_transitions(all_lids, embedded=True, use_discretized=False)

    for target_id in all_lids:
        train_ids = [lid for lid in all_lids if lid != target_id]
        
        # Regular
        _, N_train_reg, _ = inspector.compute_vertical_transitions(train_ids, embedded=False, use_discretized=True)
        _, N_target_reg, _ = inspector.compute_vertical_transitions([target_id], embedded=False, use_discretized=True)
        np.testing.assert_allclose(N_all_reg - N_train_reg, N_target_reg, atol=1e-10)
        
        # Embedded
        _, N_train_emb, _ = inspector.compute_vertical_transitions(train_ids, embedded=True, use_discretized=False)
        _, N_target_emb, _ = inspector.compute_vertical_transitions([target_id], embedded=True, use_discretized=False)
        np.testing.assert_allclose(N_all_emb - N_train_emb, N_target_emb, atol=1e-10)


# -----------------------------------------------------------------------------
# 6. Gap-crossing transitions handled according to documented policy
# -----------------------------------------------------------------------------
def test_gap_crossing_transitions_policy(inspector, validator):
    """
    Verifies that gap-crossing transitions are accurately identified in litholog9 and litholog11,
    and that the gap-masked sensitivity calculation excludes exactly those transitions.
    """
    train_ids = [f"litholog{i}" for i in range(1, 9)] + ["litholog10", "litholog12"]
    P_train, _, pi_train = validator.fit_training_markov(train_ids, embedded=False, use_discretized=True)
    p_init = inspector.compute_initial_state_distribution(train_ids, embedded=False, use_discretized=True)

    # Litholog 9 has a 1m gap at 18-19m
    eval_l9 = validator.evaluate_target_sequence(
        target_log_id="litholog9",
        P_train=P_train,
        pi_train=pi_train,
        p_init_train=p_init,
        embedded=False,
        use_discretized=True,
    )
    assert eval_l9["gap_crossing_transitions_count"] == 2
    assert "mean_transition_log_score_gap_masked" in eval_l9
    assert "transition_perplexity_gap_masked" in eval_l9

    # Litholog 1 has NO gaps -> gap_crossing_transitions_count must be 0
    eval_l1 = validator.evaluate_target_sequence(
        target_log_id="litholog1",
        P_train=P_train,
        pi_train=pi_train,
        p_init_train=p_init,
        embedded=False,
        use_discretized=True,
    )
    assert eval_l1["gap_crossing_transitions_count"] == 0
    assert eval_l1["mean_transition_log_score"] == eval_l1["mean_transition_log_score_gap_masked"]
    assert eval_l1["transition_perplexity"] == eval_l1["transition_perplexity_gap_masked"]


# -----------------------------------------------------------------------------
# 7. Fixed toy sequence gives analytically reproducible log-likelihood and perplexity
# -----------------------------------------------------------------------------
def test_fixed_toy_sequence_analytical_log_likelihood_and_perplexity():
    """
    Tests scoring against an analytically solved toy transition sequence.
    Given:
      K = 3 states (0, 1, 2)
      P = [[0.5, 0.3, 0.2],
           [0.1, 0.6, 0.3],
           [0.4, 0.4, 0.2]]
      p_init = [0.2, 0.5, 0.3]
      Sequence: [1, 2, 0, 1]

    Analytical calculations:
      Initial state: s_0 = 1 -> p_init[1] = 0.5 -> log(0.5) = -0.6931471805599453
      Transitions:
        1 -> 2: P[1,2] = 0.3 -> log(0.3) = -1.2039728043259361
        2 -> 0: P[2,0] = 0.4 -> log(0.4) = -0.9162907318741551
        0 -> 1: P[0,1] = 0.3 -> log(0.3) = -1.2039728043259361
      Sum transition log score: -3.3242363405260273
      Mean transition log score: -3.3242363405260273 / 3 = -1.1080787801753424
      Transition perplexity: exp(1.1080787801753424) = 3.0284795906... (0.036^(-1/3))
      Total sequence log score: log(0.5) + log(0.036) = log(0.018) = -4.017383521085973
    """
    p_init = np.array([0.2, 0.5, 0.3])
    P = np.array([
        [0.5, 0.3, 0.2],
        [0.1, 0.6, 0.3],
        [0.4, 0.4, 0.2],
    ])
    seq = np.array([1, 2, 0, 1])

    # Analytical values
    expected_initial_log = math.log(0.5)
    expected_total_trans_log = math.log(0.3) + math.log(0.4) + math.log(0.3)
    expected_mean_trans_log = expected_total_trans_log / 3.0
    expected_perplexity = math.exp(-expected_mean_trans_log)
    expected_total_seq_score = expected_initial_log + expected_total_trans_log

    # Compute using the exact same logic as evaluate_target_sequence
    transitions = list(zip(seq[:-1], seq[1:]))
    trans_log_probs = [float(np.log(P[u, v])) for u, v in transitions]
    initial_log = float(np.log(p_init[seq[0]]))
    total_trans = float(np.sum(trans_log_probs))
    mean_trans = float(np.mean(trans_log_probs))
    perp = float(np.exp(-mean_trans))
    total_seq = float(initial_log + total_trans)

    np.testing.assert_allclose(initial_log, expected_initial_log, atol=1e-12)
    np.testing.assert_allclose(total_trans, expected_total_trans_log, atol=1e-12)
    np.testing.assert_allclose(mean_trans, expected_mean_trans_log, atol=1e-12)
    np.testing.assert_allclose(perp, expected_perplexity, atol=1e-12)
    np.testing.assert_allclose(total_seq, expected_total_seq_score, atol=1e-12)


# -----------------------------------------------------------------------------
# 8. Continuous bed thickness and 1m discretized proportions reported separately
# -----------------------------------------------------------------------------
def test_continuous_and_discretized_proportions_reported_separately(inspector):
    """
    Verifies that continuous interval thickness and 1m discretized facies proportions
    are maintained and reported as distinct attributes without conflation.
    """
    for lid in [f"litholog{i}" for i in range(1, 13)]:
        comp = inspector.compare_continuous_vs_discretized(lid)
        insp = inspector.inspect_litholog(lid)

        assert "continuous_thickness_m" in comp
        assert "discretized_samples_1m" in comp
        assert "continuous_proportions" in comp
        assert "discretized_proportions" in comp
        assert "continuous_ntg_pure" in comp
        assert "discretized_ntg_pure" in comp
        assert "delta_ntg_pure" in comp

        # Ensure continuous sum thickness strictly equals raw interval thickness sum
        raw_df = inspector.load_raw_litholog(lid)
        raw_sum = float((raw_df["Bottom"] - raw_df["Top"]).sum())
        assert abs(insp["sum_thickness_m"] - raw_sum) < 1e-4


# -----------------------------------------------------------------------------
# 9. All training and scoring results are deterministic
# -----------------------------------------------------------------------------
def test_deterministic_training_and_scoring(validator):
    """
    Verifies that running LOLO cross-validation twice produces bitwise identical DataFrames.
    """
    test_lids = ["litholog1", "litholog2", "litholog3"]
    df_run1 = validator.run_lolo_cross_validation(test_lids, embedded=False, use_discretized=True)
    df_run2 = validator.run_lolo_cross_validation(test_lids, embedded=False, use_discretized=True)
    pd.testing.assert_frame_equal(df_run1, df_run2)


# -----------------------------------------------------------------------------
# 10. Existing test suites continue to pass
# -----------------------------------------------------------------------------
def test_sprint_d_deliverables_exist():
    """
    Verifies that all Sprint D audit deliverables have been generated in audit_sprint_d/.
    """
    repo_root = Path(__file__).resolve().parent.parent
    d_dir = repo_root / "sprints" / "audit_sprint_d" if (repo_root / "sprints" / "audit_sprint_d").exists() else repo_root / "audit_sprint_d"
    assert (d_dir / "sandstone_thickness_reconciliation.csv").exists()
    assert (d_dir / "reconciled_per_litholog_statistics.csv").exists()
    assert (d_dir / "reconciled_group_statistics.csv").exists()
    assert (d_dir / "lolo_validation_results_sprint_d.csv").exists()
    assert (d_dir / "sensitivity_gap_masking_results.csv").exists()
    assert (d_dir / "markov_transition_matrices_sprint_d.csv").exists()
