import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import pandas as pd
import numpy as np

# Facies encoding dictionary and base GR values for synthetic generation
FACIES_DICT = {
    "sand": 0,
    "mud": 1,
    "coal": 2,
    "carbon_mud": 3,
    "silt": 4,
    "siltstone": 4,
}

REVERSE_FACIES_DICT = {
    0: "sand",
    1: "mud",
    2: "coal",
    3: "carbon_mud",
    4: "silt",
}

BASE_GR_FACIES = {
    "sand": 35.0,
    "mud": 110.0,
    "coal": 20.0,
    "carbon_mud": 95.0,
    "silt": 75.0,
    "siltstone": 75.0,
    0: 35.0,
    1: 110.0,
    2: 20.0,
    3: 95.0,
    4: 75.0,
}

STANDARDIZED_COLUMNS = [
    'litholog_id',
    'strike_pos_m',
    'depth_m',
    'facies_code',
    'facies_name',
    'gamma_ray'
]


class LithologLoader:
    """
    Loader and Quality Control pipeline for raw litholog data.
    
    Loads CSV files from raw directory, standardizes columns, applies QC checks
    (monotonic depth, anomaly detection, missing Gamma Ray synthesis), and exports
    unified datasets to parquet and csv formats.
    """

    def __init__(
        self,
        raw_dir: Union[str, Path] = "data/raw_lithologs",
        processed_dir: Union[str, Path] = "data/processed",
        seed: Optional[int] = 42,
    ):
        self.raw_dir = Path(raw_dir)
        self.processed_dir = Path(processed_dir)
        self.seed = seed

    def _ensure_raw_dir(self) -> List[Path]:
        """Ensure raw_dir exists. If empty, auto-populate from data/ if raw CSV files exist there."""
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        csv_files = list(self.raw_dir.glob("*.csv"))

        if not csv_files:
            # Check root data/ for sample litholog CSV files
            data_dir = self.raw_dir.parent
            if data_dir.exists():
                for sample_file in data_dir.glob("litholog*.csv"):
                    dest = self.raw_dir / sample_file.name
                    if not dest.exists():
                        shutil.copy(sample_file, dest)
            csv_files = list(self.raw_dir.glob("*.csv"))

        return sorted(csv_files)

    def load_raw_files(self) -> pd.DataFrame:
        """
        Loads all CSV files from raw_dir and standardizes column names.
        Supports both interval layer format (Top, Bottom, Facies) and point depth format.
        """
        csv_files = self._ensure_raw_dir()
        if not csv_files:
            raise FileNotFoundError(f"No CSV files found in raw directory: {self.raw_dir}")

        processed_dfs = []

        for csv_path in csv_files:
            df_raw = pd.read_csv(csv_path)
            litholog_id = csv_path.stem

            # Normalize column names (case-insensitive strip)
            col_map = {col: col.strip().lower() for col in df_raw.columns}
            df = df_raw.rename(columns=col_map)

            # Check if interval format (top, bottom, facies)
            if 'top' in df.columns and 'bottom' in df.columns:
                rows = []
                # Sort by top, bottom to maintain depth order
                df = df.sort_values(by=['top', 'bottom']).reset_index(drop=True)
                last_depth = None

                for _, row in df.iterrows():
                    top = row['top']
                    bottom = row['bottom']
                    
                    # QC check for zero-thickness or negative thickness intervals in raw data
                    if bottom <= top:
                        raise ValueError(
                            f"Zero-thickness or negative thickness anomaly detected in {csv_path.name}: "
                            f"Top={top}, Bottom={bottom}"
                        )

                    # Handle overlapping intervals if any
                    if last_depth is not None and top < last_depth:
                        top = last_depth
                        if bottom <= top:
                            continue
                    
                    facies_val = row.get('facies', row.get('facies_name', row.get('facies_code', None)))
                    strike_pos = row.get('strike_pos_m', row.get('strike_pos', row.get('strike', 0.0)))
                    gr_val = row.get('gamma_ray', row.get('gr', row.get('gamma', np.nan)))

                    for d in range(int(top), int(bottom)):
                        rows.append({
                            'litholog_id': litholog_id,
                            'strike_pos_m': float(strike_pos),
                            'depth_m': float(d),
                            'facies_input': facies_val,
                            'gamma_ray': float(gr_val) if pd.notna(gr_val) else np.nan
                        })
                        last_depth = d + 1
                df_points = pd.DataFrame(rows)

            else:
                # Point depth format
                df_points = pd.DataFrame()
                
                # Litholog ID
                if 'litholog_id' in df.columns:
                    df_points['litholog_id'] = df['litholog_id'].astype(str)
                elif 'id' in df.columns:
                    df_points['litholog_id'] = df['id'].astype(str)
                else:
                    df_points['litholog_id'] = litholog_id

                # Strike position
                if 'strike_pos_m' in df.columns:
                    df_points['strike_pos_m'] = df['strike_pos_m'].astype(float)
                elif 'strike_pos' in df.columns:
                    df_points['strike_pos_m'] = df['strike_pos'].astype(float)
                else:
                    df_points['strike_pos_m'] = 0.0

                # Depth
                if 'depth_m' in df.columns:
                    df_points['depth_m'] = df['depth_m'].astype(float)
                elif 'depth' in df.columns:
                    df_points['depth_m'] = df['depth'].astype(float)
                else:
                    raise KeyError(f"Missing depth column in {csv_path.name}")

                # Facies input
                if 'facies_name' in df.columns:
                    df_points['facies_input'] = df['facies_name']
                elif 'facies' in df.columns:
                    df_points['facies_input'] = df['facies']
                elif 'facies_code' in df.columns:
                    df_points['facies_input'] = df['facies_code']
                else:
                    raise KeyError(f"Missing facies column in {csv_path.name}")

                # Gamma Ray
                if 'gamma_ray' in df.columns:
                    df_points['gamma_ray'] = df['gamma_ray'].astype(float)
                elif 'gr' in df.columns:
                    df_points['gamma_ray'] = df['gr'].astype(float)
                elif 'gamma' in df.columns:
                    df_points['gamma_ray'] = df['gamma'].astype(float)
                else:
                    df_points['gamma_ray'] = np.nan

            # Parse and map facies code & name
            df_points['facies_code'], df_points['facies_name'] = self._standardize_facies(df_points['facies_input'])
            df_points = df_points.drop(columns=['facies_input'])

            processed_dfs.append(df_points[STANDARDIZED_COLUMNS])

        unified_df = pd.concat(processed_dfs, ignore_index=True)
        return unified_df

    def _standardize_facies(self, facies_series: pd.Series) -> Tuple[pd.Series, pd.Series]:
        """Standardize facies to matching code and name based on FACIES_DICT."""
        codes = []
        names = []
        for val in facies_series:
            if isinstance(val, (int, np.integer)) or (isinstance(val, float) and val.is_integer() and not pd.isna(val)):
                code = int(val)
                if code not in REVERSE_FACIES_DICT:
                    raise ValueError(f"Unknown facies code: {code}")
                name = REVERSE_FACIES_DICT[code]
            else:
                name_str = str(val).strip().lower()
                if name_str not in FACIES_DICT:
                    raise ValueError(f"Unknown facies name: {val}")
                code = FACIES_DICT[name_str]
                name = name_str
            codes.append(code)
            names.append(name)
        return pd.Series(codes, index=facies_series.index), pd.Series(names, index=facies_series.index)

    def run_quality_control(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Quality control checks:
        1. Negative depth check
        2. Monotonic increasing depth check per litholog_id
        3. Synthetic Gamma Ray generation if missing (GR = base_GR + N(0, 5))
        """
        df_qc = df.copy()

        # Check negative depth anomalies
        if (df_qc['depth_m'] < 0).any():
            invalid_rows = df_qc[df_qc['depth_m'] < 0]
            raise ValueError(f"Negative depth anomaly detected in litholog(s): {invalid_rows['litholog_id'].unique().tolist()}")

        # Check monotonic depth ordering per litholog_id
        for litholog_id, group in df_qc.groupby('litholog_id'):
            depths = group['depth_m'].values
            if len(depths) > 1:
                diffs = np.diff(depths)
                if np.any(diffs <= 0):
                    raise ValueError(
                        f"Depth is not strictly monotonically increasing for litholog_id '{litholog_id}'."
                    )

        # Synthetic Gamma Ray generation if missing
        missing_gr_mask = df_qc['gamma_ray'].isna()
        if missing_gr_mask.any():
            if self.seed is not None:
                np.random.seed(self.seed)
            
            synthetic_gr = df_qc.loc[missing_gr_mask, 'facies_name'].map(
                lambda f_name: BASE_GR_FACIES[f_name] + np.random.normal(loc=0.0, scale=5.0)
            )
            df_qc.loc[missing_gr_mask, 'gamma_ray'] = synthetic_gr

        return df_qc

    def export_data(
        self,
        df: pd.DataFrame,
        parquet_path: Optional[Union[str, Path]] = None,
        csv_path: Optional[Union[str, Path]] = None,
    ) -> Tuple[Path, Path]:
        """Exports processed DataFrame to parquet and csv formats."""
        self.processed_dir.mkdir(parents=True, exist_ok=True)

        out_parquet = Path(parquet_path) if parquet_path else self.processed_dir / "lithologs_unified.parquet"
        out_csv = Path(csv_path) if csv_path else self.processed_dir / "lithologs_unified.csv"

        df.to_parquet(out_parquet, index=False)
        df.to_csv(out_csv, index=False)

        return out_parquet, out_csv

    def load_and_process(
        self,
        parquet_path: Optional[Union[str, Path]] = None,
        csv_path: Optional[Union[str, Path]] = None,
    ) -> pd.DataFrame:
        """Runs complete ingestion, quality control, and export pipeline."""
        df_raw = self.load_raw_files()
        df_clean = self.run_quality_control(df_raw)
        self.export_data(df_clean, parquet_path=parquet_path, csv_path=csv_path)
        return df_clean
