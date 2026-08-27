"""
Unit tests for LithologLoader (Phase 0 data ingestion and standardization).
"""

from pathlib import Path
import tempfile
import numpy as np
import pandas as pd
import pytest

from data.loader import (
    LithologLoader,
    DEFAULT_FACIES_MAPPING,
    DEFAULT_BASE_GR,
    CRITICAL_COLUMNS,
)


@pytest.fixture
def sample_raw_dir(tmp_path):
    """Creates a temporary raw lithologs directory with valid sample CSV files."""
    raw_dir = tmp_path / "raw_lithologs"
    raw_dir.mkdir(parents=True, exist_ok=True)

    # Log 1: Layer interval format
    log1_content = """Top,Bottom,Facies
0,5,mud
5,10,sand
10,12,coal
12,18,carbon_mud
18,22,silt
"""
    (raw_dir / "litholog1.csv").write_text(log1_content, encoding="utf-8")

    # Log 2: Layer interval with explicit GR and strike_pos_m
    log2_content = """Top,Bottom,Facies,strike_pos_m,gamma_ray
0,4,sand,150.0,38.2
4,8,mud,150.0,102.5
8,10,silt,150.0,71.0
"""
    (raw_dir / "litholog2.csv").write_text(log2_content, encoding="utf-8")

    return raw_dir


@pytest.fixture
def default_loader(sample_raw_dir, tmp_path):
    """Fixture providing a LithologLoader pointing to the temporary sample directory."""
    output_dir = tmp_path / "processed"
    return LithologLoader(raw_dir=sample_raw_dir, output_dir=output_dir, random_state=42)


def test_standardize_column_names(default_loader):
    """Verifies that standardized columns match the exact required list."""
    df = default_loader.process_all()
    assert list(df.columns) == CRITICAL_COLUMNS


def test_monotonic_depth_ordering(default_loader):
    """Verifies depth_m is strictly monotonically increasing for each litholog_id."""
    df = default_loader.process_all()
    for litholog_id, group in df.groupby("litholog_id"):
        diffs = group["depth_m"].diff().dropna()
        assert len(diffs) > 0, f"Litholog {litholog_id} has insufficient rows."
        assert (diffs > 0).all(), f"Litholog {litholog_id} depth is not strictly monotonically increasing."


def test_no_nan_values_in_critical_columns(default_loader):
    """Verifies that no critical columns contain NaN or null values."""
    df = default_loader.process_all()
    for col in CRITICAL_COLUMNS:
        assert not df[col].isna().any(), f"Column '{col}' contains NaN values."
        assert not df[col].isnull().any(), f"Column '{col}' contains null values."


def test_facies_encoding_consistency(default_loader):
    """Verifies that facies_code and facies_name match the facies encoding dictionary."""
    df = default_loader.process_all()
    for _, row in df.iterrows():
        fname = row["facies_name"]
        fcode = row["facies_code"]
        assert fname in DEFAULT_FACIES_MAPPING, f"Unexpected facies name: {fname}"
        assert DEFAULT_FACIES_MAPPING[fname] == fcode, f"Facies code {fcode} mismatch for {fname}"


def test_synthetic_gamma_ray_generation(sample_raw_dir, tmp_path):
    """Verifies synthetic GR generation conforms to base_GR + N(0, 5)."""
    # Create large log to verify statistical distribution
    large_log_dir = tmp_path / "large_raw"
    large_log_dir.mkdir(parents=True, exist_ok=True)

    layers = []
    current_depth = 0
    facies_cycle = ["sand", "mud", "coal", "carbon_mud", "silt"]
    for i in range(100):
        f = facies_cycle[i % len(facies_cycle)]
        layers.append(f"{current_depth},{current_depth + 10},{f}")
        current_depth += 10

    (large_log_dir / "litholog_large.csv").write_text(
        "Top,Bottom,Facies\n" + "\n".join(layers), encoding="utf-8"
    )

    loader = LithologLoader(raw_dir=large_log_dir, output_dir=tmp_path / "out", random_state=123)
    df = loader.process_all()

    for facies_name, base_val in DEFAULT_BASE_GR.items():
        subset = df[df["facies_name"] == facies_name]["gamma_ray"]
        assert len(subset) > 50, f"Insufficient samples for {facies_name}"

        mean_gr = subset.mean()
        std_gr = subset.std()

        # Mean within 3 standard errors of the mean
        assert abs(mean_gr - base_val) < 1.5, (
            f"Synthetic GR mean {mean_gr:.2f} deviated from base {base_val} for {facies_name}"
        )
        # Std dev close to 5.0 (between 3.5 and 6.5)
        assert 3.5 <= std_gr <= 6.5, (
            f"Synthetic GR std {std_gr:.2f} deviated from expected 5.0 for {facies_name}"
        )


def test_qc_negative_depth_anomaly(tmp_path):
    """Verifies that negative depths in raw data raise a ValueError."""
    bad_dir = tmp_path / "bad_neg_depth"
    bad_dir.mkdir()
    (bad_dir / "log_neg.csv").write_text("Top,Bottom,Facies\n-5,10,sand\n", encoding="utf-8")

    loader = LithologLoader(raw_dir=bad_dir, output_dir=tmp_path / "out")
    with pytest.raises(ValueError, match="Negative depth anomaly detected"):
        loader.process_all()


def test_qc_zero_thickness_anomaly(tmp_path):
    """Verifies that zero-thickness or inverted layers raise a ValueError."""
    bad_dir = tmp_path / "bad_zero_thick"
    bad_dir.mkdir()
    (bad_dir / "log_zero.csv").write_text("Top,Bottom,Facies\n10,10,sand\n", encoding="utf-8")

    loader = LithologLoader(raw_dir=bad_dir, output_dir=tmp_path / "out")
    with pytest.raises(ValueError, match="Zero-thickness or negative depth anomaly detected"):
        loader.process_all()


def test_qc_negative_thickness_anomaly(tmp_path):
    """Verifies that Bottom < Top raises a ValueError."""
    bad_dir = tmp_path / "bad_inv_thick"
    bad_dir.mkdir()
    (bad_dir / "log_inv.csv").write_text("Top,Bottom,Facies\n20,15,mud\n", encoding="utf-8")

    loader = LithologLoader(raw_dir=bad_dir, output_dir=tmp_path / "out")
    with pytest.raises(ValueError, match="Zero-thickness or negative depth anomaly detected"):
        loader.process_all()


def test_export_parquet_and_csv(default_loader, tmp_path):
    """Verifies exporting cleaned data to .parquet and .csv files."""
    df_cleaned = default_loader.process_all()
    parquet_path, csv_path = default_loader.export(df=df_cleaned)

    assert parquet_path.exists(), f"Parquet file {parquet_path} was not created."
    assert csv_path.exists(), f"CSV file {csv_path} was not created."

    df_p = pd.read_parquet(parquet_path)
    df_c = pd.read_csv(csv_path)

    assert len(df_p) == len(df_cleaned)
    assert len(df_c) == len(df_cleaned)
    assert list(df_p.columns) == CRITICAL_COLUMNS
    assert list(df_c.columns) == CRITICAL_COLUMNS


def test_real_workspace_raw_lithologs():
    """Runs loader against real workspace raw_lithologs files if present."""
    workspace_raw_dir = Path("data/raw_lithologs")
    if workspace_raw_dir.exists() and list(workspace_raw_dir.glob("*.csv")):
        loader = LithologLoader(raw_dir=workspace_raw_dir, output_dir="data/processed")
        df = loader.process_all()
        p_path, c_path = loader.export(df)

        assert p_path.exists()
        assert c_path.exists()
        assert len(df) > 0
        for litholog_id, group in df.groupby("litholog_id"):
            assert (group["depth_m"].diff().dropna() > 0).all()
        for col in CRITICAL_COLUMNS:
            assert not df[col].isna().any()
