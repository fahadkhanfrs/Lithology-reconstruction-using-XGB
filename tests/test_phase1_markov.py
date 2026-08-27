"""
Unit and property-based test suite for Phase 1: StratigraphicMarkovChain.

Tests mathematical invariants, row-stochasticity, embedded zero-diagonal property,
upward directionality, unobserved state fallback, stationary eigensolver,
and reproducible sampling.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest

from smalt.geostat.markov import StratigraphicMarkovChain, DEFAULT_FACIES_MAP
from smalt.viz.markov_viz import (
    plot_transition_matrix,
    plot_facies_succession_network,
    plot_stationary_vs_empirical,
)


@pytest.fixture
def synthetic_fining_upward_df():
    """
    Creates a synthetic 10-well dataset with guaranteed fining-upward cycles:
    Lower (deep) -> Upper (shallow): Sand (1) -> Silt (3) -> Mud (4).
    """
    rows = []
    for well_id in range(1, 11):
        # 5 cycles per well; depth goes from 150m down to 0m (shallow)
        depth = 150.0
        for _ in range(5):
            # Base of cycle: Channel Sand (1) (e.g. 5m thick)
            for _ in range(5):
                rows.append({"litholog_id": f"well_{well_id}", "depth_m": depth, "facies_code": 1})
                depth -= 1.0
            # Middle: Siltstone (3) (e.g. 3m thick)
            for _ in range(3):
                rows.append({"litholog_id": f"well_{well_id}", "depth_m": depth, "facies_code": 3})
                depth -= 1.0
            # Top: Overbank Mudstone (4) (e.g. 2m thick)
            for _ in range(2):
                rows.append({"litholog_id": f"well_{well_id}", "depth_m": depth, "facies_code": 4})
                depth -= 1.0
    return pd.DataFrame(rows)


@pytest.fixture
def real_litholog_df():
    """Loads unified litholog dataset if available."""
    path = Path("data/processed/lithologs_unified.parquet")
    if path.exists():
        return pd.read_parquet(path)
    # Fallback minimal 3-well test frame
    rows = []
    for w in [1, 9, 11]:
        for d, f in zip(range(0, 30), [4, 4, 3, 1, 1, 1, 0, 2, 2, 4] * 3):
            rows.append({"litholog_id": f"litholog_{w}", "depth_m": float(d), "facies_code": f})
    return pd.DataFrame(rows)


def test_row_stochasticity_strict(synthetic_fining_upward_df, real_litholog_df):
    """
    Asserts that every row of matrix P sums to 1.0 within 1e-7 tolerance
    for both regular and embedded chains across synthetic and real data.
    """
    for df in [synthetic_fining_upward_df, real_litholog_df]:
        for embedded in [False, True]:
            for alpha in [0.0, 0.5]:
                model = StratigraphicMarkovChain(num_classes=5, embedded=embedded, smoothing_alpha=alpha)
                model.fit(df)
                P = model.transition_matrix_
                assert P is not None
                row_sums = P.sum(axis=1)
                assert np.allclose(row_sums, 1.0, atol=1e-7), (
                    f"Row stochasticity failed (embedded={embedded}, alpha={alpha}): "
                    f"max error = {np.max(np.abs(row_sums - 1.0)):.2e}"
                )


def test_probability_bounds(synthetic_fining_upward_df, real_litholog_df):
    """Asserts that all elements of P satisfy 0.0 <= P_ij <= 1.0."""
    for df in [synthetic_fining_upward_df, real_litholog_df]:
        for embedded in [False, True]:
            model = StratigraphicMarkovChain(num_classes=5, embedded=embedded)
            model.fit(df)
            P = model.transition_matrix_
            assert np.all((P >= 0.0) & (P <= 1.0)), "Probability bounds violated (P < 0 or P > 1)."


def test_upward_directionality_synthetic(synthetic_fining_upward_df):
    """
    Tests that upward stratigraphic ordering (decreasing depth) recovers the prescribed
    fining-upward transitions: Sand (1) -> Silt (3) -> Mud (4).
    """
    model = StratigraphicMarkovChain(num_classes=5, embedded=True, smoothing_alpha=0.0)
    model.fit(synthetic_fining_upward_df)
    P = model.transition_matrix_

    # Upward transitions: Sand (1) -> Silt (3) and Silt (3) -> Mud (4)
    assert P[1, 3] > 0.85, f"Expected P[1, 3] > 0.85, got {P[1, 3]:.3f}"
    assert P[3, 4] > 0.85, f"Expected P[3, 4] > 0.85, got {P[3, 4]:.3f}"

    # Reverse transitions: Mud (4) -> Silt (3) and Silt (3) -> Sand (1)
    assert P[4, 3] < 0.10, f"Expected reverse P[4, 3] < 0.10, got {P[4, 3]:.3f}"
    assert P[3, 1] < 0.10, f"Expected reverse P[3, 1] < 0.10, got {P[3, 1]:.3f}"


def test_embedded_markov_zero_diagonal(synthetic_fining_upward_df):
    """Asserts that embedded Markov chains have exactly zero on the diagonal (P_ii == 0.0)."""
    model = StratigraphicMarkovChain(num_classes=5, embedded=True, smoothing_alpha=0.0)
    model.fit(synthetic_fining_upward_df)
    P = model.transition_matrix_
    assert np.all(np.diag(P) == 0.0), f"Embedded diagonal non-zero: {np.diag(P)}"


def test_stationary_distribution_recovery():
    """
    Verifies stationary distribution computation:
    1. On a known doubly stochastic / uniform matrix, pi == [1/K, ..., 1/K].
    2. Verifies the invariance condition ||pi * P - pi||_inf < 1e-6.
    """
    K = 5
    # Create uniform cyclic transition matrix
    P_uniform = np.zeros((K, K))
    for i in range(K):
        P_uniform[i, (i + 1) % K] = 1.0

    model = StratigraphicMarkovChain(num_classes=K)
    model.transition_matrix_ = P_uniform
    pi = model.compute_stationary_distribution()

    expected_pi = np.ones(K) / K
    assert np.allclose(pi, expected_pi, atol=1e-6), (
        f"Stationary distribution mismatch on uniform matrix: {pi}"
    )
    assert np.allclose(pi @ P_uniform, pi, atol=1e-6), "pi * P == pi invariant failed."


def test_unobserved_absorbing_state_fallback():
    """
    Verifies that unobserved states (zero outgoing counts in data) are handled
    gracefully with uniform fallback without throwing NaNs or violating row stochasticity.
    """
    # Create dataset with only classes 1, 2, 3, 4 (Class 0 is completely absent)
    df_missing_0 = pd.DataFrame(
        {
            "litholog_id": ["well_1"] * 20,
            "depth_m": list(range(20, 0, -1)),
            "facies_code": [1, 2, 3, 4] * 5,
        }
    )

    # 1. Regular chain fallback: Row 0 should be uniform 1/K across all 5 states (0.2)
    model_reg = StratigraphicMarkovChain(num_classes=5, embedded=False, smoothing_alpha=0.0)
    model_reg.fit(df_missing_0)
    P_reg = model_reg.transition_matrix_

    assert not np.isnan(P_reg).any(), "P_reg contains NaN values."
    assert np.allclose(P_reg[0, :], 0.2, atol=1e-7), f"Expected row 0 to be 0.2, got {P_reg[0, :]}"
    assert np.allclose(P_reg.sum(axis=1), 1.0, atol=1e-7)

    # 2. Embedded chain fallback: Row 0 should have P_00 = 0.0 and P_0j = 1/(K-1) = 0.25 for j != 0
    model_emb = StratigraphicMarkovChain(num_classes=5, embedded=True, smoothing_alpha=0.0)
    model_emb.fit(df_missing_0)
    P_emb = model_emb.transition_matrix_

    assert not np.isnan(P_emb).any(), "P_emb contains NaN values."
    assert P_emb[0, 0] == 0.0, "Embedded P_00 must be 0.0"
    assert np.allclose(P_emb[0, 1:], 0.25, atol=1e-7), f"Expected off-diagonal row 0 to be 0.25, got {P_emb[0, :]}"
    assert np.allclose(P_emb.sum(axis=1), 1.0, atol=1e-7)


def test_reproducible_sampling(synthetic_fining_upward_df):
    """Verifies that sample_next_state with a fixed RNG generator produces deterministic outputs."""
    model = StratigraphicMarkovChain(num_classes=5, embedded=False)
    model.fit(synthetic_fining_upward_df)

    rng1 = np.random.default_rng(seed=42)
    samples1 = [model.sample_next_state(1, rng=rng1) for _ in range(50)]

    rng2 = np.random.default_rng(seed=42)
    samples2 = [model.sample_next_state(1, rng=rng2) for _ in range(50)]

    assert samples1 == samples2, "Sampling with identical RNG seed produced non-deterministic outputs."


def test_dynamic_class_dimensions():
    """Verifies that dynamic K (e.g. K=3 or K=6) is supported without hardcoding."""
    for K in [3, 6]:
        custom_map = {i: f"CustomFacies_{i}" for i in range(K)}
        df = pd.DataFrame(
            {
                "litholog_id": ["w1"] * 30,
                "depth_m": list(range(30, 0, -1)),
                "facies_code": [i % K for i in range(30)],
            }
        )
        model = StratigraphicMarkovChain(num_classes=K, facies_map=custom_map, embedded=True)
        model.fit(df)
        assert model.transition_matrix_.shape == (K, K)
        assert model.count_matrix_.shape == (K, K)
        assert len(model.compute_stationary_distribution()) == K
        assert np.allclose(model.transition_matrix_.sum(axis=1), 1.0, atol=1e-7)


def test_state_entropy_and_asymmetry(synthetic_fining_upward_df):
    """Verifies calculation of state entropy and directional asymmetry matrix."""
    model = StratigraphicMarkovChain(num_classes=5, embedded=False)
    model.fit(synthetic_fining_upward_df)

    H_raw = model.compute_state_entropy(normalized=False)
    H_norm = model.compute_state_entropy(normalized=True)
    A = model.compute_asymmetry_matrix()

    assert len(H_raw) == 5
    assert len(H_norm) == 5
    assert np.all(H_raw >= 0.0)
    assert np.all((H_norm >= 0.0) & (H_norm <= 1.0 + 1e-7))

    # Asymmetry matrix properties: A_ij = -A_ji and A_ii = 0
    assert np.all(np.diag(A) == 0.0)
    assert np.allclose(A, -A.T, atol=1e-7)


def test_export_json(synthetic_fining_upward_df, tmp_path):
    """Verifies exporting model parameters, counts, and metrics to JSON."""
    model = StratigraphicMarkovChain(num_classes=5, embedded=True, smoothing_alpha=0.1)
    model.fit(synthetic_fining_upward_df)

    json_path = tmp_path / "markov_summary.json"
    exported_path = model.export_json(json_path)

    assert exported_path.exists()
    with open(exported_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    assert "count_matrix_N" in data
    assert "transition_matrix_P" in data
    assert "stationary_distribution_pi" in data
    assert "asymmetry_matrix_A" in data
    assert "state_entropy_H" in data
    assert data["metadata"]["n_wells"] == 10
    assert data["metadata"]["embedded"] is True


def test_visualization_smoke(synthetic_fining_upward_df, tmp_path):
    """Smoke test verifying visualization routines render cleanly without error."""
    model = StratigraphicMarkovChain(num_classes=5, embedded=False)
    model.fit(synthetic_fining_upward_df)

    fig1_path = tmp_path / "tm.png"
    fig2_path = tmp_path / "net.png"
    fig3_path = tmp_path / "stat.png"

    plot_transition_matrix(model, title="Test Transition Matrix", save_path=fig1_path)
    plot_facies_succession_network(model, threshold=0.05, save_path=fig2_path)
    plot_stationary_vs_empirical(model, synthetic_fining_upward_df, save_path=fig3_path)

    assert fig1_path.exists()
    assert fig2_path.exists()
    assert fig3_path.exists()
