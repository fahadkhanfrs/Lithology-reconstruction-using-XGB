"""
LithologLoader module for Phase 0 lithology data ingestion,
standardization, quality control, synthetic GR generation, and export.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import re
import numpy as np
import pandas as pd


DEFAULT_FACIES_MAPPING: Dict[str, int] = {
    "coal": 0,
    "sand": 1,
    "carbon_mud": 2,
    "silt": 3,
    "mud": 4,
}

# Aliases to map variant naming to standard facies names
FACIES_ALIASES: Dict[str, str] = {
    "coal": "coal",
    "sand": "sand",
    "sandstone": "sand",
    "channel sandstone": "sand",
    "channel_sandstone": "sand",
    "carbon_mud": "carbon_mud",
    "carbonaceous_mud": "carbon_mud",
    "carbonaceous_mudstone": "carbon_mud",
    "carbon mud": "carbon_mud",
    "fine sandstone / splay": "carbon_mud",
    "fine_sandstone": "carbon_mud",
    "fine sandstone": "carbon_mud",
    "splay": "carbon_mud",
    "silt": "silt",
    "siltstone": "silt",
    "mud": "mud",
    "mudstone": "mud",
    "shale": "mud",
    "overbank mudstone": "mud",
    "overbank_mudstone": "mud",
}

DEFAULT_BASE_GR: Dict[str, float] = {
    "coal": 20.0,
    "sand": 35.0,
    "carbon_mud": 130.0,
    "silt": 70.0,
    "mud": 105.0,
}

CRITICAL_COLUMNS: List[str] = [
    "litholog_id",
    "strike_pos_m",
    "depth_m",
    "facies_code",
    "facies_name",
    "gamma_ray",
]


class LithologLoader:
    """
    Loads, cleans, validates, and standardizes raw litholog CSV data.
    """

    def __init__(
        self,
        raw_dir: Union[str, Path] = "data/raw_lithologs",
        output_dir: Union[str, Path] = "data/processed",
        facies_mapping: Optional[Dict[str, int]] = None,
        base_gr: Optional[Dict[str, float]] = None,
        strike_positions: Optional[Dict[str, float]] = None,
        random_state: int = 42,
    ):
        self.raw_dir = Path(raw_dir)
        self.output_dir = Path(output_dir)
        self.facies_mapping = facies_mapping or DEFAULT_FACIES_MAPPING.copy()
        self.code_to_facies = {v: k for k, v in self.facies_mapping.items()}
        self.base_gr = base_gr or DEFAULT_BASE_GR.copy()
        self.strike_positions = strike_positions or {}
        self.random_state = random_state
        self.rng = np.random.default_rng(random_state)

    def _normalize_facies_name(self, name: Union[str, int, float]) -> str:
        """Normalizes facies name strings and maps common aliases."""
        if pd.isna(name):
            raise ValueError("Facies name cannot be null or NaN.")
        clean_str = str(name).strip().lower().replace(" ", "_")
        if clean_str in FACIES_ALIASES:
            return FACIES_ALIASES[clean_str]
        if clean_str in self.facies_mapping:
            return clean_str
        raise ValueError(
            f"Unknown facies '{name}'. Expected one of: {list(self.facies_mapping.keys())}"
        )

    def _extract_litholog_id(self, filepath: Path, df: pd.DataFrame) -> str:
        """Extracts or derives a clean litholog identifier."""
        for col in df.columns:
            if col.lower().strip() in ["litholog_id", "litholog", "well_id", "well"]:
                return str(df[col].iloc[0])
        return filepath.stem

    def _extract_strike_pos(self, litholog_id: str, df: pd.DataFrame, file_idx: int) -> float:
        """Extracts or assigns strike position in meters."""
        for col in df.columns:
            if col.lower().strip() in ["strike_pos_m", "strike_pos", "strike", "x_pos"]:
                val = df[col].iloc[0]
                if not pd.isna(val):
                    return float(val)

        # Check explicit user mapping
        if litholog_id in self.strike_positions:
            return float(self.strike_positions[litholog_id])

        # Infer from litholog number if present (e.g. litholog1 -> 0, litholog9 -> 800, etc.)
        match = re.search(r"(\d+)", litholog_id)
        if match:
            num = int(match.group(1))
            return float((num - 1) * 100.0)

        return float(file_idx * 100.0)

    def check_layer_anomalies(self, df_raw: pd.DataFrame, source_name: str = "") -> None:
        """
        Validates raw layer intervals for zero-thickness, negative thickness,
        or negative depth anomalies.
        """
        cols_lower = {str(c).strip().lower(): c for c in df_raw.columns}

        top_col = cols_lower.get("top") or cols_lower.get("top_m") or cols_lower.get("from")
        bottom_col = cols_lower.get("bottom") or cols_lower.get("bottom_m") or cols_lower.get("to")

        if top_col and bottom_col:
            tops = pd.to_numeric(df_raw[top_col], errors="coerce")
            bottoms = pd.to_numeric(df_raw[bottom_col], errors="coerce")

            if tops.isna().any() or bottoms.isna().any():
                raise ValueError(f"[{source_name}] Raw layer depths contain non-numeric or NaN values.")

            # Check negative depths
            if (tops < 0).any() or (bottoms < 0).any():
                neg_tops = tops[tops < 0].tolist()
                neg_bottoms = bottoms[bottoms < 0].tolist()
                raise ValueError(
                    f"[{source_name}] Negative depth anomaly detected: tops={neg_tops}, bottoms={neg_bottoms}"
                )

            # Check zero or negative thickness
            thickness = bottoms - tops
            if (thickness <= 0).any():
                invalid = df_raw[thickness <= 0][[top_col, bottom_col]].to_dict(orient="records")
                raise ValueError(
                    f"[{source_name}] Zero-thickness or negative depth anomaly detected in layers: {invalid}"
                )

        # Point depth check if depth column exists
        depth_col = cols_lower.get("depth") or cols_lower.get("depth_m")
        if depth_col:
            depths = pd.to_numeric(df_raw[depth_col], errors="coerce")
            if depths.isna().any():
                raise ValueError(f"[{source_name}] Raw depth column contains NaN values.")
            if (depths < 0).any():
                raise ValueError(f"[{source_name}] Negative depth anomaly detected: depths < 0")

    def standardize_dataframe(
        self, df_raw: pd.DataFrame, filepath: Optional[Union[str, Path]] = None, file_idx: int = 0
    ) -> pd.DataFrame:
        """
        Standardizes raw DataFrame to target schema:
        ['litholog_id', 'strike_pos_m', 'depth_m', 'facies_code', 'facies_name', 'gamma_ray']
        """
        source_name = Path(filepath).name if filepath else f"df_{file_idx}"
        self.check_layer_anomalies(df_raw, source_name=source_name)

        path_obj = Path(filepath) if filepath else Path(f"litholog_{file_idx}")
        litholog_id = self._extract_litholog_id(path_obj, df_raw)
        strike_pos_m = self._extract_strike_pos(litholog_id, df_raw, file_idx)

        cols_lower = {str(c).strip().lower(): c for c in df_raw.columns}
        top_col = cols_lower.get("top") or cols_lower.get("top_m") or cols_lower.get("from")
        bottom_col = cols_lower.get("bottom") or cols_lower.get("bottom_m") or cols_lower.get("to")
        facies_col = cols_lower.get("facies") or cols_lower.get("facies_name") or cols_lower.get("lithology")
        gr_col = (
            cols_lower.get("gamma_ray")
            or cols_lower.get("gamma")
            or cols_lower.get("gr")
            or cols_lower.get("gr_api")
        )

        rows = []

        # Case 1: Layer-interval format (Top, Bottom, Facies)
        if top_col and bottom_col and facies_col:
            for _, r in df_raw.iterrows():
                top = int(np.floor(float(r[top_col])))
                bottom = int(np.ceil(float(r[bottom_col])))
                facies_raw = r[facies_col]
                facies_name = self._normalize_facies_name(facies_raw)
                facies_code = self.facies_mapping[facies_name]
                raw_gr = float(r[gr_col]) if gr_col and not pd.isna(r[gr_col]) else np.nan

                for d in range(top, bottom):
                    rows.append(
                        {
                            "litholog_id": str(litholog_id),
                            "strike_pos_m": float(strike_pos_m),
                            "depth_m": float(d),
                            "facies_code": int(facies_code),
                            "facies_name": str(facies_name),
                            "gamma_ray": raw_gr,
                        }
                    )
            df_std = pd.DataFrame(rows)

        # Case 2: Discrete depth point format
        elif "depth_m" in cols_lower or "depth" in cols_lower:
            d_col = cols_lower.get("depth_m") or cols_lower.get("depth")
            df_work = df_raw.copy()

            df_std = pd.DataFrame()
            df_std["litholog_id"] = (
                df_work[cols_lower["litholog_id"]].astype(str)
                if "litholog_id" in cols_lower
                else str(litholog_id)
            )
            df_std["strike_pos_m"] = (
                df_work[cols_lower["strike_pos_m"]].astype(float)
                if "strike_pos_m" in cols_lower
                else float(strike_pos_m)
            )
            df_std["depth_m"] = df_work[d_col].astype(float)

            # Determine facies
            if facies_col:
                df_std["facies_name"] = df_work[facies_col].apply(self._normalize_facies_name)
                df_std["facies_code"] = df_std["facies_name"].map(self.facies_mapping)
            elif "facies_code" in cols_lower:
                df_std["facies_code"] = df_work[cols_lower["facies_code"]].astype(int)
                df_std["facies_name"] = df_std["facies_code"].map(self.code_to_facies)
                if df_std["facies_name"].isna().any():
                    raise ValueError("Unknown facies_code encountered in point dataset.")
            else:
                raise ValueError("Neither facies name nor facies code found in raw data.")

            if gr_col:
                df_std["gamma_ray"] = pd.to_numeric(df_work[gr_col], errors="coerce")
            else:
                df_std["gamma_ray"] = np.nan
        else:
            raise ValueError(
                f"Unsupported raw format in {source_name}. Expected either (Top, Bottom, Facies) or Depth columns."
            )

        # Deduplicate depth points within this litholog (keeping last/first resolved) and sort
        df_std = df_std.drop_duplicates(subset=["depth_m"], keep="last")
        df_std = df_std.sort_values(by="depth_m").reset_index(drop=True)

        return df_std[CRITICAL_COLUMNS]

    def generate_synthetic_gr(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Generates synthetic Gamma Ray values (GR = base_GR_facies + N(0, 5))
        for missing (NaN) GR values.
        """
        df = df.copy()
        nan_mask = df["gamma_ray"].isna()
        if not nan_mask.any():
            return df

        for facies_name, base_val in self.base_gr.items():
            mask = nan_mask & (df["facies_name"] == facies_name)
            count = mask.sum()
            if count > 0:
                noise = self.rng.normal(loc=0.0, scale=5.0, size=count)
                df.loc[mask, "gamma_ray"] = base_val + noise

        # Fallback for any other facies not in base_gr
        rem_mask = df["gamma_ray"].isna()
        if rem_mask.any():
            count = rem_mask.sum()
            noise = self.rng.normal(loc=0.0, scale=5.0, size=count)
            df.loc[rem_mask, "gamma_ray"] = 50.0 + noise

        return df

    def run_quality_control(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Executes strict quality control checks:
        1. Checks for negative depth anomalies.
        2. Asserts depth_m is strictly monotonically increasing for each litholog_id.
        3. Generates synthetic Gamma Ray values where missing.
        4. Asserts no NaN values in critical columns.
        5. Asserts facies encoding matches dictionary.
        """
        df = df.copy()

        # Check negative depth
        if (df["depth_m"] < 0).any():
            neg_depths = df[df["depth_m"] < 0]["depth_m"].tolist()
            raise ValueError(f"Negative depth anomaly found: {neg_depths}")

        # Check strict monotonic depth increase per litholog_id
        for litholog_id, group in df.groupby("litholog_id"):
            diffs = group["depth_m"].diff().dropna()
            if len(diffs) > 0 and not (diffs > 0).all():
                non_mono_idx = diffs[diffs <= 0].index
                raise AssertionError(
                    f"Litholog '{litholog_id}' violates strict monotonic depth ordering at indices {non_mono_idx.tolist()}."
                )

        # Impute missing Gamma Ray
        df = self.generate_synthetic_gr(df)

        # Asserts no NaNs in critical columns
        for col in CRITICAL_COLUMNS:
            if df[col].isna().any():
                nan_count = df[col].isna().sum()
                raise AssertionError(f"Critical column '{col}' contains {nan_count} NaN values after QC.")

        # Asserts facies code matches facies name according to dictionary
        for facies_name, code in self.facies_mapping.items():
            mismatches = df[df["facies_name"] == facies_name]["facies_code"] != code
            if mismatches.any():
                raise AssertionError(f"Facies code mismatch for facies '{facies_name}'. Expected {code}.")

        return df

    def load_raw_files(self) -> List[Tuple[Path, pd.DataFrame]]:
        """Loads all CSV files from raw_dir."""
        if not self.raw_dir.exists():
            raise FileNotFoundError(f"Raw directory '{self.raw_dir}' does not exist.")

        csv_files = sorted(list(self.raw_dir.glob("*.csv")))
        if not csv_files:
            raise FileNotFoundError(f"No CSV files found in '{self.raw_dir}'.")

        dataframes = []
        for p in csv_files:
            df = pd.read_csv(p)
            dataframes.append((p, df))
        return dataframes

    def process_all(self) -> pd.DataFrame:
        """
        Loads all raw litholog CSVs, standardizes them, runs QC checks,
        and returns a unified DataFrame.
        """
        raw_files = self.load_raw_files()
        standardized_dfs = []

        for idx, (filepath, df_raw) in enumerate(raw_files):
            std_df = self.standardize_dataframe(df_raw, filepath=filepath, file_idx=idx)
            standardized_dfs.append(std_df)

        unified_df = pd.concat(standardized_dfs, ignore_index=True)
        cleaned_df = self.run_quality_control(unified_df)
        return cleaned_df

    def export(
        self,
        df: Optional[pd.DataFrame] = None,
        parquet_path: Optional[Union[str, Path]] = None,
        csv_path: Optional[Union[str, Path]] = None,
    ) -> Tuple[Path, Path]:
        """
        Exports cleaned data to parquet and csv.
        Defaults to data/processed/lithologs_unified.parquet and .csv.
        """
        if df is None:
            df = self.process_all()

        self.output_dir.mkdir(parents=True, exist_ok=True)

        p_path = Path(parquet_path) if parquet_path else self.output_dir / "lithologs_unified.parquet"
        c_path = Path(csv_path) if csv_path else self.output_dir / "lithologs_unified.csv"

        p_path.parent.mkdir(parents=True, exist_ok=True)
        c_path.parent.mkdir(parents=True, exist_ok=True)

        df.to_parquet(p_path, index=False)
        df.to_csv(c_path, index=False)

        return p_path, c_path
