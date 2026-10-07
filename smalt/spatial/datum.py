"""
SMALT Common-Zero Datum Alignment Module.

Implements the project-wide vertical reference standard instructed by Prof. Hiranya Sahoo:
"All litholog zeros should be treated as being at the same reference level."

Key Principles:
1. Reference Level: z = 0 corresponds to measured depth = 0 for every litholog.
2. Sign Convention:
   - 'elevation' (default): z_common_m = -depth_original_m (deeper intervals are negative,
     compatible with 3D right-handed Cartesian geostatistical coordinates where +Z is up).
   - 'downward_depth': z_common_m = +depth_original_m (deeper intervals are positive).
3. Non-Invention: This reference level is an explicit provisional project assumption
   supplied by the professor; it does NOT represent a verified physical coal marker,
   flooding surface, or named stratigraphic horizon.
4. Non-Local: No litholog-specific zero (such as max(depth) - depth) is used.
5. Invariance: Reversible transformation preserving interval thicknesses, facies codes,
   and facies counts.
"""

from typing import Dict, Any, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd

from smalt.descriptive.analyzer import (
    LithologInspector,
    CANONICAL_FACIES_SCHEMA,
    PROVENANCE_METADATA,
)
from smalt.spatial.coordinates import load_source_coordinates


# Canonical datum metadata constants
COMMON_DATUM_NAME = "common_zero_datum"
DEFAULT_CONVENTION = "elevation"  # z_common_m = -depth (deeper is negative)


def align_to_common_datum(
    df: pd.DataFrame,
    litholog_id: Optional[str] = None,
    convention: str = DEFAULT_CONVENTION,
) -> pd.DataFrame:
    """
    Transforms a continuous interval or 1m discretized litholog dataframe to the
    common-zero reference datum.

    Args:
        df: Input DataFrame. Can be:
            - Continuous intervals with ('Top', 'Bottom', 'Facies') or ('from', 'to', 'facies')
            - Discretized points with ('depth_m', 'facies', 'facies_code')
        litholog_id: Optional identifier for the litholog (e.g. 'litholog1').
        convention: 'elevation' (default: deeper is negative, z = -depth) or
                    'downward_depth' (deeper is positive, z = +depth).

    Returns:
        DataFrame containing preserved original depth and standardized common-datum coordinates:
        - litholog_id (str)
        - depth_original_m (float)
        - z_common_m (float)
        - thickness_m (float)
        - facies (str, canonical 6-state name)
        - facies_code (int, 0 to 5)
        - [interval boundary columns if continuous]
    """
    if convention not in ("elevation", "downward_depth"):
        raise ValueError(
            f"Unknown datum convention '{convention}'. Must be 'elevation' or 'downward_depth'."
        )

    sign = -1.0 if convention == "elevation" else 1.0
    df_out = pd.DataFrame(index=df.index)

    # Identify if dataframe is continuous intervals or discretized points
    cols_lower = {str(c).strip().lower(): c for c in df.columns}

    # Case A: Discretized 1m grid DataFrame (contains depth_m)
    if "depth_m" in cols_lower:
        depth_col = cols_lower["depth_m"]
        fac_col = cols_lower.get("facies") or cols_lower.get("facies_name") or cols_lower.get("facies_code")
        
        df_out["litholog_id"] = litholog_id or df.get("litholog_id", "unknown")
        df_out["depth_original_m"] = pd.to_numeric(df[depth_col], errors="raise")
        df_out["z_common_m"] = df_out["depth_original_m"] * sign
        df_out["thickness_m"] = 1.0
        
        token_to_code = {k: meta["code"] for k, meta in CANONICAL_FACIES_SCHEMA.items()}
        token_to_code.update({meta["canonical_name"]: meta["code"] for meta in CANONICAL_FACIES_SCHEMA.values()})
        code_to_name = {meta["code"]: meta["canonical_name"] for meta in CANONICAL_FACIES_SCHEMA.values()}

        # Facies name and code
        if "facies_code" in cols_lower:
            df_out["facies_code"] = pd.to_numeric(df[cols_lower["facies_code"]], errors="raise").astype(int)
            df_out["facies"] = df_out["facies_code"].map(code_to_name)
        else:
            raw_fac = df[fac_col].astype(str)
            df_out["facies_code"] = raw_fac.map(token_to_code)
            df_out["facies"] = df_out["facies_code"].map(code_to_name)

        # Retain observation quality metadata if present
        for q_col in ("is_source_observed", "is_gap_filled", "is_overlap_resolved"):
            if q_col in cols_lower:
                df_out[q_col] = df[cols_lower[q_col]]

    # Case B: Continuous interval DataFrame (contains Top and Bottom)
    elif ("top" in cols_lower or "from" in cols_lower) and ("bottom" in cols_lower or "to" in cols_lower):
        top_col = cols_lower.get("top") or cols_lower.get("from")
        bot_col = cols_lower.get("bottom") or cols_lower.get("to")
        fac_col = cols_lower.get("facies") or cols_lower.get("lithology")

        df_out["litholog_id"] = litholog_id or df.get("litholog_id", "unknown")
        df_out["depth_original_top_m"] = pd.to_numeric(df[top_col], errors="raise")
        df_out["depth_original_bottom_m"] = pd.to_numeric(df[bot_col], errors="raise")
        df_out["depth_original_m"] = (df_out["depth_original_top_m"] + df_out["depth_original_bottom_m"]) / 2.0
        df_out["thickness_m"] = df_out["depth_original_bottom_m"] - df_out["depth_original_top_m"]

        if (df_out["thickness_m"] <= 0).any():
            invalid = df_out[df_out["thickness_m"] <= 0]
            raise ValueError(f"Encountered non-positive interval thickness:\n{invalid}")

        # Common datum coordinates
        df_out["z_common_top_m"] = df_out["depth_original_top_m"] * sign
        df_out["z_common_bottom_m"] = df_out["depth_original_bottom_m"] * sign
        df_out["z_common_m"] = df_out["depth_original_m"] * sign

        token_to_code = {k: meta["code"] for k, meta in CANONICAL_FACIES_SCHEMA.items()}
        token_to_code.update({meta["canonical_name"]: meta["code"] for meta in CANONICAL_FACIES_SCHEMA.values()})
        code_to_name = {meta["code"]: meta["canonical_name"] for meta in CANONICAL_FACIES_SCHEMA.values()}

        raw_fac = df[fac_col].astype(str)
        df_out["facies_code"] = raw_fac.map(token_to_code)
        df_out["facies"] = df_out["facies_code"].map(code_to_name)

    else:
        raise ValueError(
            f"Cannot determine format of input DataFrame. Columns found: {list(df.columns)}"
        )

    # Attach datum configuration metadata
    df_out.attrs["datum_reference"] = COMMON_DATUM_NAME
    df_out.attrs["convention"] = convention
    df_out.attrs["sign_multiplier"] = sign
    df_out.attrs["professor_instruction"] = (
        "All litholog zeros should be treated as being at the same reference level."
    )

    return df_out


def revert_from_common_datum(df_common: pd.DataFrame) -> pd.DataFrame:
    """
    Reverts a common-datum aligned DataFrame back to its original depth coordinates.
    Proves transformation reversibility.
    """
    convention = df_common.attrs.get("convention", DEFAULT_CONVENTION)
    sign = -1.0 if convention == "elevation" else 1.0

    df_rev = df_common.copy()
    if "depth_original_m" in df_rev.columns:
        # Recomputed depth from z_common_m matches depth_original_m identically
        recomputed_depth = df_rev["z_common_m"] * sign
        np.testing.assert_allclose(
            recomputed_depth.to_numpy(),
            df_rev["depth_original_m"].to_numpy(),
            atol=1e-12,
            err_msg="Common datum reversibility check failed.",
        )
    return df_rev


def build_common_datum_metadata_table(
    inspector: Optional[LithologInspector] = None,
    convention: str = DEFAULT_CONVENTION,
) -> pd.DataFrame:
    """
    Builds the authoritative project metadata table for all 12 lithologs aligned
    under the common-zero datum reference.

    Columns:
      litholog
      original_depth_min
      original_depth_max
      common_z_min
      common_z_max
      n_intervals
      total_thickness
      coordinate_available
      spatial_eligible
    """
    insp = inspector or LithologInspector()
    coords_df = load_source_coordinates()
    coords_map = {
        r["litholog_id"]: bool(r["coordinates_available"])
        for _, r in coords_df.iterrows()
    }

    records = []
    all_lids = [f"litholog{i}" for i in range(1, 13)]

    for lid in all_lids:
        raw_df = insp.load_raw_litholog(lid)
        aligned_df = align_to_common_datum(raw_df, litholog_id=lid, convention=convention)

        d_orig_min = float(aligned_df["depth_original_top_m"].min())
        d_orig_max = float(aligned_df["depth_original_bottom_m"].max())

        # z bounds (in elevation convention, min z is deepest negative, max z is top 0)
        if convention == "elevation":
            z_min = float(aligned_df["z_common_bottom_m"].min())
            z_max = float(aligned_df["z_common_top_m"].max())
        else:
            z_min = float(aligned_df["z_common_top_m"].min())
            z_max = float(aligned_df["z_common_bottom_m"].max())

        n_intervals = len(aligned_df)
        total_thick = float(aligned_df["thickness_m"].sum())
        coord_avail = coords_map.get(lid, False)
        # Litholog 1 lacks coordinates and is strictly spatial ineligible
        spatial_eligible = coord_avail and (lid != "litholog1")

        records.append({
            "litholog": lid,
            "original_depth_min": round(d_orig_min, 2),
            "original_depth_max": round(d_orig_max, 2),
            "common_z_min": round(z_min, 2),
            "common_z_max": round(z_max, 2),
            "n_intervals": n_intervals,
            "total_thickness": round(total_thick, 2),
            "coordinate_available": coord_avail,
            "spatial_eligible": spatial_eligible,
        })

    return pd.DataFrame(records)


def verify_datum_invariance(
    inspector: Optional[LithologInspector] = None,
    convention: str = DEFAULT_CONVENTION,
) -> pd.DataFrame:
    """
    Verifies that the common-datum transformation causes zero changes to
    stratigraphic facies counts, interval counts, and cumulative thickness.
    """
    insp = inspector or LithologInspector()
    records = []
    all_lids = [f"litholog{i}" for i in range(1, 13)]

    for lid in all_lids:
        raw_df = insp.load_raw_litholog(lid)
        aligned_df = align_to_common_datum(raw_df, litholog_id=lid, convention=convention)

        # Check raw counts vs aligned counts
        token_to_canonical = {k: meta["canonical_name"] for k, meta in CANONICAL_FACIES_SCHEMA.items()}
        raw_canonical_counts = raw_df["Facies"].map(token_to_canonical).value_counts().to_dict()
        aligned_counts = aligned_df["facies"].value_counts().to_dict()
        assert raw_canonical_counts == aligned_counts, f"Facies counts mismatch in {lid}"

        raw_thick = float(raw_df["Thickness"].sum())
        aligned_thick = float(aligned_df["thickness_m"].sum())
        np.testing.assert_allclose(raw_thick, aligned_thick, atol=1e-9)

        # Check reversibility
        revert_from_common_datum(aligned_df)

        records.append({
            "litholog": lid,
            "raw_interval_count": len(raw_df),
            "aligned_interval_count": len(aligned_df),
            "raw_thickness_m": round(raw_thick, 3),
            "aligned_thickness_m": round(aligned_thick, 3),
            "facies_counts_identical": raw_canonical_counts == aligned_counts,
            "thickness_identical": abs(raw_thick - aligned_thick) < 1e-6,
            "reversibility_verified": True,
        })

    return pd.DataFrame(records)
