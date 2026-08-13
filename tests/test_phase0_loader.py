import os
from pathlib import Path
import pytest
import pandas as pd
import numpy as np

from data.loader import LithologLoader, FACIES_DICT, REVERSE_FACIES_DICT, BASE_GR_FACIES, STANDARDIZED_COLUMNS


@pytest.fixture
def sample_raw_dir(tmp_path):
    """Creates temporary raw litholog files (both interval and point depth format)."""
    raw_dir = tmp_path / "raw_lithologs"
    raw_dir.mkdir()

    # File 1: Interval format
    df1 = pd.DataFrame({
        "Top": [0, 5, 10],
        "Bottom": [5, 10, 15],
        "Facies": ["sand", "mud", "coal"]
    })
    df1.to_csv(raw_dir / "log1.csv", index=False)

    # File 2: Point depth format
    df2 = pd.DataFrame({
        "litholog_id": ["log2", "log2", "log2"],
        "strike_pos_m": [10.0, 10.0, 10.0],
        "depth_m": [0.0, 1.0, 2.0],
        "facies_name": ["mud", "carbon_mud", "sand"],
        "gamma_ray": [112.5, 96.0, 34.2]
    })
    df2.to_csv(raw_dir / "log2.csv", index=False)

    return raw_dir


def test_litholog_loader_load_and_process(sample_raw_dir, tmp_path):
    """Verifies that all loaded logs pass standardized schema, monotonic depth, NaN check, and facies encoding match."""
    processed_dir = tmp_path / "processed"
    loader = LithologLoader(raw_dir=sample_raw_dir, processed_dir=processed_dir, seed=42)
    df_clean = loader.load_and_process()

    # 1. Standardized columns check
    assert list(df_clean.columns) == STANDARDIZED_COLUMNS

    # 2. Monotonic depth ordering check for each litholog_id
    for litholog_id, group in df_clean.groupby("litholog_id"):
        depths = group["depth_m"].values
        assert len(depths) > 0
        if len(depths) > 1:
            diffs = np.diff(depths)
            assert np.all(diffs > 0), f"Depths for {litholog_id} are not strictly monotonic increasing."

    # 3. No NaN values in critical columns
    critical_cols = ["litholog_id", "depth_m", "facies_code", "facies_name", "gamma_ray"]
    for col in critical_cols:
        assert df_clean[col].isna().sum() == 0, f"Column '{col}' contains NaN values."

    # 4. Facies encoding dictionary match
    for _, row in df_clean.iterrows():
        f_name = row["facies_name"]
        f_code = row["facies_code"]
        assert f_name in FACIES_DICT
        assert FACIES_DICT[f_name] == f_code
        assert REVERSE_FACIES_DICT[f_code] == f_name


def test_non_monotonic_depth_raises(tmp_path):
    """Verifies that non-monotonic depth ordering triggers a ValueError."""
    df_invalid = pd.DataFrame({
        "litholog_id": ["log_bad", "log_bad", "log_bad"],
        "strike_pos_m": [0.0, 0.0, 0.0],
        "depth_m": [10.0, 5.0, 15.0],  # Out of order
        "facies_code": [0, 1, 2],
        "facies_name": ["sand", "mud", "coal"],
        "gamma_ray": [35.0, 110.0, 20.0]
    })
    loader = LithologLoader(raw_dir=tmp_path / "raw", processed_dir=tmp_path / "proc")
    with pytest.raises(ValueError, match="not strictly monotonically increasing"):
        loader.run_quality_control(df_invalid)


def test_negative_depth_anomaly_raises(tmp_path):
    """Verifies that negative depth values trigger a ValueError."""
    df_invalid = pd.DataFrame({
        "litholog_id": ["log_neg", "log_neg"],
        "strike_pos_m": [0.0, 0.0],
        "depth_m": [-5.0, 0.0],  # Negative depth anomaly
        "facies_code": [0, 1],
        "facies_name": ["sand", "mud"],
        "gamma_ray": [35.0, 110.0]
    })
    loader = LithologLoader(raw_dir=tmp_path / "raw", processed_dir=tmp_path / "proc")
    with pytest.raises(ValueError, match="Negative depth anomaly detected"):
        loader.run_quality_control(df_invalid)


def test_zero_thickness_anomaly_raises(tmp_path):
    """Verifies that zero or negative interval thickness in raw CSV raises ValueError."""
    raw_dir = tmp_path / "raw_zero_thick"
    raw_dir.mkdir()
    df = pd.DataFrame({
        "Top": [0, 10],
        "Bottom": [10, 10],  # Zero thickness interval (10 to 10)
        "Facies": ["sand", "mud"]
    })
    df.to_csv(raw_dir / "zero_thick.csv", index=False)

    loader = LithologLoader(raw_dir=raw_dir, processed_dir=tmp_path / "proc")
    with pytest.raises(ValueError, match="Zero-thickness or negative thickness anomaly"):
        loader.load_raw_files()


def test_synthetic_gamma_ray_generation(tmp_path):
    """Verifies synthetic Gamma Ray values are generated when raw GR is missing."""
    df_missing_gr = pd.DataFrame({
        "litholog_id": ["log_syn", "log_syn", "log_syn"],
        "strike_pos_m": [0.0, 0.0, 0.0],
        "depth_m": [1.0, 2.0, 3.0],
        "facies_code": [0, 1, 2],
        "facies_name": ["sand", "mud", "coal"],
        "gamma_ray": [np.nan, np.nan, np.nan]  # All GR missing
    })
    loader = LithologLoader(raw_dir=tmp_path / "raw", processed_dir=tmp_path / "proc", seed=123)
    df_qc = loader.run_quality_control(df_missing_gr)

    assert df_qc["gamma_ray"].isna().sum() == 0

    # Sand base GR is 35, Mud is 110, Coal is 20 (with N(0, 5) noise)
    gr_vals = df_qc["gamma_ray"].tolist()
    assert 20.0 < gr_vals[0] < 50.0   # ~35
    assert 90.0 < gr_vals[1] < 130.0  # ~110
    assert 5.0 < gr_vals[2] < 35.0    # ~20


def test_export_data_parquet_and_csv(sample_raw_dir, tmp_path):
    """Verifies exports to .parquet and .csv files and verifies data integrity upon reloading."""
    processed_dir = tmp_path / "processed"
    loader = LithologLoader(raw_dir=sample_raw_dir, processed_dir=processed_dir, seed=42)
    df_clean = loader.load_and_process()

    parquet_file = processed_dir / "lithologs_unified.parquet"
    csv_file = processed_dir / "lithologs_unified.csv"

    assert parquet_file.exists()
    assert csv_file.exists()

    df_parquet = pd.read_parquet(parquet_file)
    df_csv = pd.read_csv(csv_file)

    assert len(df_parquet) == len(df_clean)
    assert len(df_csv) == len(df_clean)
    assert list(df_parquet.columns) == STANDARDIZED_COLUMNS
    assert list(df_csv.columns) == STANDARDIZED_COLUMNS
