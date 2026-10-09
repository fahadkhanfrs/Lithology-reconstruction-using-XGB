"""
SMALT Common-Zero Datum Alignment & Source Orientation Audit Module.

Implements the project-wide vertical reference standard instructed by Prof. Hiranya Sahoo:
"All litholog zeros should be treated as being at the same reference level."

Sprint H Orientation Correction:
1. Source Orientation Distinction:
   - Measured outcrop sections (L1-L11): Sahoo et al. (2016) source convention is 0 m at BASE,
     stratigraphic height increasing upward. Raw CSVs were digitized with 0 at the top and depth
     increasing downward, requiring axis reversal: z_strat = H_max - d.
   - Drill core log (L12 / EM-137C): Preserved in its digitized orientation (0-111 m) without
     inversion, isolated pending confirmation from Prof. Sahoo.
2. Common Reference Level:
   - Project reference level z_common = 0 corresponds to the source 0 m mark of each litholog.
   - For measured sections, z_common = 0 is at the source base; height increases upward.
3. Invariance:
   - Reversible transformation preserving interval thicknesses, facies codes, and counts.
   - Original measured depths are preserved in depth_original_m without destructive overwrite.
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


# Authoritative source orientation registry
LITHOLOG_ORIENTATION_METADATA: Dict[str, Dict[str, Any]] = {
    "litholog1": {
        "source_type": "measured_section",
        "source_zero_location": "base",
        "source_direction": "upward",
        "current_csv_direction": "downward",
        "required_transform": "reverse_stratigraphic_axis",
        "orientation_confidence": "high",
        "evidence_source": "Sahoo et al. (2016) measured section; 0 m at base, upward stratigraphic succession; CSV digitized from top down.",
    },
    "litholog2": {
        "source_type": "measured_section",
        "source_zero_location": "base",
        "source_direction": "upward",
        "current_csv_direction": "downward",
        "required_transform": "reverse_stratigraphic_axis",
        "orientation_confidence": "high",
        "evidence_source": "Sahoo et al. (2016) measured section; 0 m at base, upward stratigraphic succession; CSV digitized from top down.",
    },
    "litholog3": {
        "source_type": "measured_section",
        "source_zero_location": "base",
        "source_direction": "upward",
        "current_csv_direction": "downward",
        "required_transform": "reverse_stratigraphic_axis",
        "orientation_confidence": "high",
        "evidence_source": "Sahoo et al. (2016) measured section; 0 m at base, upward stratigraphic succession; CSV digitized from top down.",
    },
    "litholog4": {
        "source_type": "measured_section",
        "source_zero_location": "base",
        "source_direction": "upward",
        "current_csv_direction": "downward",
        "required_transform": "reverse_stratigraphic_axis",
        "orientation_confidence": "high",
        "evidence_source": "Sahoo et al. (2016) measured section; 0 m at base, upward stratigraphic succession; CSV digitized from top down.",
    },
    "litholog5": {
        "source_type": "measured_section",
        "source_zero_location": "base",
        "source_direction": "upward",
        "current_csv_direction": "downward",
        "required_transform": "reverse_stratigraphic_axis",
        "orientation_confidence": "high",
        "evidence_source": "Sahoo et al. (2016) measured section; 0 m at base, upward stratigraphic succession; CSV digitized from top down.",
    },
    "litholog6": {
        "source_type": "measured_section",
        "source_zero_location": "base",
        "source_direction": "upward",
        "current_csv_direction": "downward",
        "required_transform": "reverse_stratigraphic_axis",
        "orientation_confidence": "high",
        "evidence_source": "Sahoo et al. (2016) measured section; 0 m at base, upward stratigraphic succession; CSV digitized from top down.",
    },
    "litholog7": {
        "source_type": "measured_section",
        "source_zero_location": "base",
        "source_direction": "upward",
        "current_csv_direction": "downward",
        "required_transform": "reverse_stratigraphic_axis",
        "orientation_confidence": "high",
        "evidence_source": "Sahoo et al. (2016) measured section; 0 m at base, upward stratigraphic succession; CSV digitized from top down.",
    },
    "litholog8": {
        "source_type": "measured_section",
        "source_zero_location": "base",
        "source_direction": "upward",
        "current_csv_direction": "downward",
        "required_transform": "reverse_stratigraphic_axis",
        "orientation_confidence": "high",
        "evidence_source": "Sahoo et al. (2016) measured section; 0 m at base, upward stratigraphic succession; CSV digitized from top down.",
    },
    "litholog9": {
        "source_type": "measured_section",
        "source_zero_location": "base",
        "source_direction": "upward",
        "current_csv_direction": "downward",
        "required_transform": "reverse_stratigraphic_axis",
        "orientation_confidence": "high",
        "evidence_source": "Authoritative Sahoo et al. (2016) source figure confirms 0 m at base, fining-upward succession; CSV digitized from top down.",
    },
    "litholog10": {
        "source_type": "measured_section",
        "source_zero_location": "base",
        "source_direction": "upward",
        "current_csv_direction": "downward",
        "required_transform": "reverse_stratigraphic_axis",
        "orientation_confidence": "high",
        "evidence_source": "Sahoo et al. (2016) measured section; 0 m at base, upward stratigraphic succession; CSV digitized from top down.",
    },
    "litholog11": {
        "source_type": "measured_section",
        "source_zero_location": "base",
        "source_direction": "upward",
        "current_csv_direction": "downward",
        "required_transform": "reverse_stratigraphic_axis",
        "orientation_confidence": "high",
        "evidence_source": "Sahoo et al. (2016) measured section; 0 m at base, upward stratigraphic succession; CSV digitized from top down.",
    },
    "litholog12": {
        "source_type": "drill_core",
        "source_zero_location": "base",
        "source_direction": "upward",
        "current_csv_direction": "upward",
        "required_transform": "preserve",
        "orientation_confidence": "unresolved",
        "evidence_source": "Supplementary Figure DR6 (lolo/litholog12.pdf) EM-137C core; graphic column drawn 0-242 m with 0 at bottom; digitized 0-111 m preserved as measured from 0 mark without inversion.",
    },
}

COMMON_DATUM_NAME = "common_zero_datum"
DEFAULT_CONVENTION = "stratigraphic_height"


def build_orientation_audit_table() -> pd.DataFrame:
    """Returns the structured orientation audit table across all 12 lithologs."""
    records = []
    for lid, meta in LITHOLOG_ORIENTATION_METADATA.items():
        records.append({
            "litholog_id": lid,
            "source_type": meta["source_type"],
            "source_zero_location": meta["source_zero_location"],
            "source_direction": meta["source_direction"],
            "current_csv_direction": meta["current_csv_direction"],
            "required_transform": meta["required_transform"],
            "orientation_confidence": meta["orientation_confidence"],
            "evidence_source": meta["evidence_source"],
        })
    return pd.DataFrame(records)


def align_to_common_datum(
    df: pd.DataFrame,
    litholog_id: Optional[str] = None,
    apply_source_orientation: bool = True,
    convention: str = DEFAULT_CONVENTION,
) -> pd.DataFrame:
    """
    Transforms a continuous interval or 1m discretized litholog dataframe to the
    common-zero reference datum, applying source-orientation corrections where needed.

    Args:
        df: Input DataFrame (continuous intervals with Top/Bottom/Facies or discretized with depth_m).
        litholog_id: Identifier for the litholog (e.g. 'litholog9').
        apply_source_orientation: If True, applies orientation audit transforms (reversing
            stratigraphic axis for L1-L11 measured sections, preserving L12).
            If False, falls back to legacy uncorrected convention.
        convention: 'stratigraphic_height' (default: z_common >= 0, upward from base)
                    or legacy 'elevation' / 'downward_depth'.

    Returns:
        DataFrame containing preserved original depth and standardized common-datum coordinates.
    """
    lid = litholog_id or str(df.get("litholog_id", ["unknown"])[0] if "litholog_id" in df else "unknown")
    cols_lower = {str(c).strip().lower(): c for c in df.columns}
    orient_meta = LITHOLOG_ORIENTATION_METADATA.get(lid, {
        "required_transform": "reverse_stratigraphic_axis" if lid != "litholog12" else "preserve",
    })
    req_transform = orient_meta["required_transform"] if apply_source_orientation else "preserve"

    token_to_code = {k: meta["code"] for k, meta in CANONICAL_FACIES_SCHEMA.items()}
    token_to_code.update({meta["canonical_name"]: meta["code"] for meta in CANONICAL_FACIES_SCHEMA.values()})
    code_to_name = {meta["code"]: meta["canonical_name"] for meta in CANONICAL_FACIES_SCHEMA.values()}

    # Case A: Discretized 1m grid DataFrame (contains depth_m)
    if "depth_m" in cols_lower:
        depth_col = cols_lower["depth_m"]
        fac_col = cols_lower.get("facies") or cols_lower.get("facies_name") or cols_lower.get("facies_code")
        
        d_orig = pd.to_numeric(df[depth_col], errors="raise").to_numpy(dtype=float)
        H_max = float(np.ceil(np.max(d_orig) + 1.0))

        if apply_source_orientation and req_transform == "reverse_stratigraphic_axis":
            # Invert: base is at z = 0, top is at z = H_max - 1
            z_strat = (H_max - 1.0) - d_orig
        else:
            if convention == "elevation" and not apply_source_orientation:
                z_strat = -d_orig
            else:
                z_strat = d_orig

        df_out = pd.DataFrame(index=range(len(d_orig)))
        df_out["litholog_id"] = lid
        df_out["depth_original_m"] = d_orig
        df_out["z_strat_m"] = z_strat
        df_out["z_common_m"] = z_strat
        df_out["thickness_m"] = 1.0

        if "facies_code" in cols_lower:
            df_out["facies_code"] = pd.to_numeric(df[cols_lower["facies_code"]], errors="raise").astype(int)
            df_out["facies"] = df_out["facies_code"].map(code_to_name)
        else:
            raw_fac = df[fac_col].astype(str)
            df_out["facies_code"] = raw_fac.map(token_to_code)
            df_out["facies"] = df_out["facies_code"].map(code_to_name)

        for q_col in ("is_source_observed", "is_gap_filled", "is_overlap_resolved"):
            if q_col in cols_lower:
                df_out[q_col] = df[cols_lower[q_col]].to_numpy()

        # Sort in upward stratigraphic order (increasing z_common_m)
        df_out = df_out.sort_values(by="z_common_m", ascending=True).reset_index(drop=True)

    # Case B: Continuous interval DataFrame (contains Top and Bottom)
    elif ("top" in cols_lower or "from" in cols_lower) and ("bottom" in cols_lower or "to" in cols_lower):
        top_col = cols_lower.get("top") or cols_lower.get("from")
        bot_col = cols_lower.get("bottom") or cols_lower.get("to")
        fac_col = cols_lower.get("facies") or cols_lower.get("lithology")

        d_top = pd.to_numeric(df[top_col], errors="raise").to_numpy(dtype=float)
        d_bot = pd.to_numeric(df[bot_col], errors="raise").to_numpy(dtype=float)
        thick = d_bot - d_top
        d_mid = (d_top + d_bot) / 2.0

        if (thick <= 0).any():
            invalid_idx = np.where(thick <= 0)[0]
            raise ValueError(f"Encountered non-positive interval thickness at indices {invalid_idx}")

        H_max = float(np.max(d_bot))

        if apply_source_orientation and req_transform == "reverse_stratigraphic_axis":
            # For measured sections: source 0 m is at the base
            # z_strat = H_max - d
            z_base = H_max - d_bot
            z_top = H_max - d_top
            z_mid = H_max - d_mid
        else:
            if convention == "elevation" and not apply_source_orientation:
                z_base = -d_bot
                z_top = -d_top
                z_mid = -d_mid
            else:
                z_base = d_top
                z_top = d_bot
                z_mid = d_mid

        raw_fac = df[fac_col].astype(str)
        f_codes = raw_fac.map(token_to_code).to_numpy()
        f_names = [code_to_name.get(c, "unknown") for c in f_codes]

        df_out = pd.DataFrame({
            "litholog_id": lid,
            "depth_original_top_m": d_top,
            "depth_original_bottom_m": d_bot,
            "depth_original_m": d_mid,
            "z_strat_base_m": z_base,
            "z_strat_top_m": z_top,
            "z_strat_mid_m": z_mid,
            "z_common_base_m": z_base,
            "z_common_top_m": z_top,
            "z_common_m": z_mid,
            "thickness_m": thick,
            "facies": f_names,
            "facies_code": f_codes,
        })

        # Sort in upward stratigraphic order (increasing z_common_base_m)
        df_out = df_out.sort_values(by="z_common_base_m", ascending=True).reset_index(drop=True)

    else:
        raise ValueError(f"Cannot determine format of input DataFrame. Columns found: {list(df.columns)}")

    df_out.attrs["datum_reference"] = COMMON_DATUM_NAME
    df_out.attrs["convention"] = convention
    df_out.attrs["apply_source_orientation"] = apply_source_orientation
    df_out.attrs["required_transform"] = req_transform
    df_out.attrs["H_max"] = H_max
    df_out.attrs["professor_instruction"] = "All litholog zeros should be treated as being at the same reference level."

    return df_out


def revert_from_common_datum(df_common: pd.DataFrame) -> pd.DataFrame:
    """
    Reverts a common-datum aligned DataFrame back to its original depth coordinates.
    Validates exact mathematical reversibility.
    """
    req_transform = df_common.attrs.get("required_transform", "preserve")
    apply_source_orientation = df_common.attrs.get("apply_source_orientation", True)
    H_max = df_common.attrs.get("H_max")

    df_rev = df_common.copy()
    if "z_strat_mid_m" in df_rev.columns:
        if apply_source_orientation and req_transform == "reverse_stratigraphic_axis":
            recomputed_mid = H_max - df_rev["z_strat_mid_m"]
        else:
            recomputed_mid = df_rev["z_strat_mid_m"]
        np.testing.assert_allclose(
            recomputed_mid.to_numpy(),
            df_rev["depth_original_m"].to_numpy(),
            atol=1e-12,
            err_msg="Continuous common datum reversibility check failed.",
        )
    elif "z_common_m" in df_rev.columns and "depth_original_m" in df_rev.columns:
        if apply_source_orientation and req_transform == "reverse_stratigraphic_axis":
            recomputed_d = (H_max - 1.0) - df_rev["z_common_m"]
        else:
            recomputed_d = df_rev["z_common_m"]
        np.testing.assert_allclose(
            recomputed_d.to_numpy(),
            df_rev["depth_original_m"].to_numpy(),
            atol=1e-12,
            err_msg="Discretized common datum reversibility check failed.",
        )
    return df_rev


def build_common_datum_metadata_table(
    inspector: Optional[LithologInspector] = None,
    apply_source_orientation: bool = True,
) -> pd.DataFrame:
    """
    Builds the authoritative project metadata table for all 12 lithologs aligned
    under the common-zero datum reference.
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
        aligned_df = align_to_common_datum(raw_df, litholog_id=lid, apply_source_orientation=apply_source_orientation)

        d_orig_min = float(aligned_df["depth_original_top_m"].min())
        d_orig_max = float(aligned_df["depth_original_bottom_m"].max())
        z_min = float(aligned_df["z_common_base_m"].min())
        z_max = float(aligned_df["z_common_top_m"].max())

        n_intervals = len(aligned_df)
        total_thick = float(aligned_df["thickness_m"].sum())
        coord_avail = coords_map.get(lid, False)
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
    apply_source_orientation: bool = True,
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
        aligned_df = align_to_common_datum(raw_df, litholog_id=lid, apply_source_orientation=apply_source_orientation)

        token_to_canonical = {k: meta["canonical_name"] for k, meta in CANONICAL_FACIES_SCHEMA.items()}
        raw_canonical_counts = raw_df["Facies"].map(token_to_canonical).value_counts().to_dict()
        aligned_counts = aligned_df["facies"].value_counts().to_dict()
        assert raw_canonical_counts == aligned_counts, f"Facies counts mismatch in {lid}"

        raw_thick = float(raw_df["Thickness"].sum())
        aligned_thick = float(aligned_df["thickness_m"].sum())
        np.testing.assert_allclose(raw_thick, aligned_thick, atol=1e-9)

        # Verify reversibility
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
