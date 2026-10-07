"""
SMALT Descriptive Lithology & Interval Validation Engine (Sprint C).

Provides mathematically rigorous, geologically grounded inspection of raw
stratigraphic intervals, continuous bed thicknesses, facies proportions,
net-to-gross ratios, and vertical facies succession modeling.

Strictly preserves provenance, documents facies mapping caveats, and avoids
conflating continuous bed observations with 1-meter discretized grid cells.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
import numpy as np
import pandas as pd


# Documented provenance manifest for all 12 lithologs
PROVENANCE_METADATA: Dict[str, Dict[str, Any]] = {
    "litholog1": {
        "litholog_id": "litholog1",
        "provenance_category": "source_derived",
        "provenance_description": "Original field section from Sahoo et al. (2016). Source-derived outcrop log.",
        "independent_validation": "unknown",
        "benchmark_accuracy": None,
        "coordinates_available": False,
        "coordinates_status": "No coordinates documented; permanently excluded from spatial modeling.",
        "group": "unassigned_no_coords",
    },
    "litholog2": {
        "litholog_id": "litholog2",
        "provenance_category": "ai_reconstructed",
        "provenance_description": "AI-assisted reconstruction from published photomosaics/figures.",
        "independent_validation": "unvalidated",
        "benchmark_accuracy": None,
        "coordinates_available": True,
        "coordinates_status": "Documented in Location_coordinates_lithologs.xlsx (Row 1). CRS/units unverified.",
        "group": "upstream",
    },
    "litholog3": {
        "litholog_id": "litholog3",
        "provenance_category": "ai_reconstructed",
        "provenance_description": "AI-assisted reconstruction from published photomosaics/figures.",
        "independent_validation": "unvalidated",
        "benchmark_accuracy": None,
        "coordinates_available": True,
        "coordinates_status": "Documented in Location_coordinates_lithologs.xlsx (Row 2). CRS/units unverified.",
        "group": "upstream",
    },
    "litholog4": {
        "litholog_id": "litholog4",
        "provenance_category": "ai_reconstructed",
        "provenance_description": "AI-assisted reconstruction from published photomosaics/figures.",
        "independent_validation": "unvalidated",
        "benchmark_accuracy": None,
        "coordinates_available": True,
        "coordinates_status": "Documented in Location_coordinates_lithologs.xlsx (Row 3). CRS/units unverified.",
        "group": "upstream",
    },
    "litholog5": {
        "litholog_id": "litholog5",
        "provenance_category": "ai_reconstructed",
        "provenance_description": "AI-assisted reconstruction from published photomosaics/figures.",
        "independent_validation": "unvalidated",
        "benchmark_accuracy": None,
        "coordinates_available": True,
        "coordinates_status": "Documented in Location_coordinates_lithologs.xlsx (Row 4). CRS/units unverified.",
        "group": "upstream",
    },
    "litholog6": {
        "litholog_id": "litholog6",
        "provenance_category": "ai_reconstructed",
        "provenance_description": "AI-assisted reconstruction from published photomosaics/figures.",
        "independent_validation": "unvalidated",
        "benchmark_accuracy": None,
        "coordinates_available": True,
        "coordinates_status": "Documented in Location_coordinates_lithologs.xlsx (Row 5). CRS/units unverified.",
        "group": "upstream",
    },
    "litholog7": {
        "litholog_id": "litholog7",
        "provenance_category": "ai_reconstructed",
        "provenance_description": "AI-assisted reconstruction from published photomosaics/figures.",
        "independent_validation": "unvalidated",
        "benchmark_accuracy": None,
        "coordinates_available": True,
        "coordinates_status": "Documented in Location_coordinates_lithologs.xlsx (Row 6). CRS/units unverified.",
        "group": "upstream",
    },
    "litholog8": {
        "litholog_id": "litholog8",
        "provenance_category": "ai_reconstructed",
        "provenance_description": "AI-assisted reconstruction from published photomosaics/figures.",
        "independent_validation": "unvalidated",
        "benchmark_accuracy": None,
        "coordinates_available": True,
        "coordinates_status": "Documented in Location_coordinates_lithologs.xlsx (Row 7). CRS/units unverified.",
        "group": "upstream",
    },
    "litholog9": {
        "litholog_id": "litholog9",
        "provenance_category": "source_derived",
        "provenance_description": "Original field section from Sahoo et al. (2016). Source-derived outcrop log.",
        "independent_validation": "unknown",
        "benchmark_accuracy": None,
        "coordinates_available": True,
        "coordinates_status": "Documented in Location_coordinates_lithologs.xlsx (Row 8). CRS/units unverified.",
        "group": "downstream",
    },
    "litholog10": {
        "litholog_id": "litholog10",
        "provenance_category": "ai_reconstructed",
        "provenance_description": "AI-assisted reconstruction from published photomosaics/figures.",
        "independent_validation": "unvalidated",
        "benchmark_accuracy": None,
        "coordinates_available": True,
        "coordinates_status": "Documented in Location_coordinates_lithologs.xlsx (Row 9). CRS/units unverified.",
        "group": "upstream",
    },
    "litholog11": {
        "litholog_id": "litholog11",
        "provenance_category": "source_derived_benchmarked",
        "provenance_description": "Benchmarked section from Sahoo et al. (2016). 93.59% match against manual extraction.",
        "independent_validation": "validated",
        "benchmark_accuracy": 0.9359,
        "coordinates_available": True,
        "coordinates_status": "Documented in Location_coordinates_lithologs.xlsx (Row 10). CRS/units unverified.",
        "group": "downstream",
    },
    "litholog12": {
        "litholog_id": "litholog12",
        "provenance_category": "digitized_core_log",
        "provenance_description": "Subsurface core log EM-137C. Manually digitized by user from core-log image. Covers 0-111 m (original length 242 m).",
        "independent_validation": "unvalidated",
        "benchmark_accuracy": None,
        "coordinates_available": True,
        "coordinates_status": "Documented in Location_coordinates_lithologs.xlsx (Row 11: X=-1119.91, Y=14407.30).",
        "group": "downstream",
    },
}

# Canonical 5-state facies schema
CANONICAL_FACIES_SCHEMA: Dict[str, Dict[str, Any]] = {
    "coal": {
        "code": 0,
        "canonical_name": "Coal",
        "lithology_type": "Biogenic organic deposit / mire facies",
        "color_hex": "#1C2833",
        "geological_caveat": "Unambiguous biogenic coal bed. Forms critical chronostratigraphic markers.",
    },
    "sand": {
        "code": 1,
        "canonical_name": "Channel Sandstone (undivided)",
        "lithology_type": "Coarse siliciclastic / channel-belt sandstone",
        "color_hex": "#F4D03F",
        "geological_caveat": (
            "Mapped as 'Channel Sandstone' in Phase 0 schema. Represents all undivided sandstone intervals "
            "in raw logs; crevasse splay and levee sandstones cannot be distinguished without grain-size/structure logs."
        ),
    },
    "carbon_mud": {
        "code": 2,
        "canonical_name": "Carbonaceous Mudstone",
        "lithology_type": "Organic-rich muddy sediment / poorly-drained swamp",
        "color_hex": "#6C3483",
        "geological_caveat": (
            "CRITICAL CONTRADICTION: Phase 0 CANONICAL_FACIES_NAMES labeled code 2 as 'Fine Sandstone / Splay', "
            "yet assigned Base GR = 130.0 API and technical ID 'carbon_mud'. Geologically, this is carbonaceous "
            "mudstone (organic-rich mud), NOT a fine sandstone or crevasse splay."
        ),
    },
    "silt": {
        "code": 3,
        "canonical_name": "Siltstone",
        "lithology_type": "Fine-grained siliciclastic / floodplain-levee transition",
        "color_hex": "#73C6B6",
        "geological_caveat": "Unambiguous siltstone. Dropped by legacy XGBoost script due to dictionary omission.",
    },
    "mud": {
        "code": 4,
        "canonical_name": "Overbank Mudstone",
        "lithology_type": "Fine siliciclastic / floodplain-overbank mudstone",
        "color_hex": "#95A5A6",
        "geological_caveat": "Overbank floodplain mudstone; may include minor abandoned channel-plug muds.",
    },
}

RAW_FACIES_ALIASES: Dict[str, str] = {
    "coal": "coal",
    "sand": "sand",
    "sandstone": "sand",
    "channel_sandstone": "sand",
    "channel sandstone": "sand",
    "p_sand": "sand",
    "planar_sand": "sand",
    "planar_sandstone": "sand",
    "c_sand": "sand",
    "channel_sand": "sand",
    "ripples": "silt",
    "rippled_sandstone": "silt",
    "ripple": "silt",
    "carbon_mud": "carbon_mud",
    "carbonaceous_mud": "carbon_mud",
    "carbonaceous mud": "carbon_mud",
    "carbonaceous_mudstone": "carbon_mud",
    "carbonaceous mudstone": "carbon_mud",
    "silt": "silt",
    "siltstone": "silt",
    "mud": "mud",
    "mudstone": "mud",
    "overbank_mudstone": "mud",
    "overbank mudstone": "mud",
    "shale": "mud",
}


def normalize_facies_label(label: Any) -> str:
    """
    Normalizes a facies label string to its canonical technical identifier.

    Raises ValueError if label is unknown or null.
    """
    if pd.isna(label):
        raise ValueError("Encountered null or NaN facies label.")
    cleaned = str(label).strip().lower().replace("-", "_").replace(" ", "_")
    if cleaned in RAW_FACIES_ALIASES:
        return RAW_FACIES_ALIASES[cleaned]
    if cleaned in CANONICAL_FACIES_SCHEMA:
        return cleaned
    raise ValueError(f"Unknown facies label '{label}'. Cannot map to canonical SMALT classes.")


class LithologInspector:
    """
    Validates interval coverage, detects gaps and overlaps, and calculates
    rigorous descriptive statistics for continuous and discretized stratigraphy.
    """

    def __init__(self, raw_dir: Union[str, Path] = "data/raw_lithologs"):
        self.raw_dir = Path(raw_dir)

    def load_raw_litholog(self, litholog_id: str) -> pd.DataFrame:
        """Loads and normalizes raw litholog CSV file."""
        csv_path = self.raw_dir / f"{litholog_id}.csv"
        if not csv_path.exists():
            raise FileNotFoundError(f"Raw litholog file not found: {csv_path}")

        df = pd.read_csv(csv_path)
        cols_lower = {str(c).strip().lower(): c for c in df.columns}

        top_col = cols_lower.get("top") or cols_lower.get("from")
        bot_col = cols_lower.get("bottom") or cols_lower.get("to")
        fac_col = cols_lower.get("facies") or cols_lower.get("lithology")

        if not (top_col and bot_col and fac_col):
            raise ValueError(f"File {csv_path} lacks required (Top, Bottom, Facies) columns. Found: {list(df.columns)}")

        df_out = pd.DataFrame()
        df_out["Top"] = pd.to_numeric(df[top_col], errors="raise")
        df_out["Bottom"] = pd.to_numeric(df[bot_col], errors="raise")
        df_out["Raw_Facies"] = df[fac_col].astype(str)
        df_out["Facies"] = df_out["Raw_Facies"].apply(normalize_facies_label)
        df_out["Thickness"] = df_out["Bottom"] - df_out["Top"]

        if (df_out["Thickness"] <= 0).any():
            invalid = df_out[df_out["Thickness"] <= 0]
            raise ValueError(f"File {csv_path} contains non-positive thickness intervals:\n{invalid}")

        return df_out

    def inspect_litholog(self, litholog_id: str) -> Dict[str, Any]:
        """Performs rigorous interval-continuity and boundary inspection for a litholog."""
        df = self.load_raw_litholog(litholog_id)
        prov = PROVENANCE_METADATA.get(litholog_id, {
            "litholog_id": litholog_id,
            "provenance_category": "unknown",
            "provenance_description": "Unknown provenance",
            "independent_validation": "unknown",
            "benchmark_accuracy": None,
            "coordinates_available": False,
            "coordinates_status": "Unknown",
            "group": "unassigned",
        })

        min_depth = float(df["Top"].min())
        max_depth = float(df["Bottom"].max())
        represented_span = float(max_depth - min_depth)
        sum_thickness = float(df["Thickness"].sum())
        span_thickness_diff = float(represented_span - sum_thickness)

        gaps: List[Dict[str, Any]] = []
        overlaps: List[Dict[str, Any]] = []

        for i in range(len(df) - 1):
            curr_bottom = float(df.iloc[i]["Bottom"])
            next_top = float(df.iloc[i + 1]["Top"])
            delta = next_top - curr_bottom

            if delta > 1e-4:
                gaps.append({
                    "interval_index_before": i,
                    "depth_before": curr_bottom,
                    "depth_after": next_top,
                    "gap_size_m": round(delta, 3),
                    "facies_before": df.iloc[i]["Facies"],
                    "facies_after": df.iloc[i + 1]["Facies"],
                })
            elif delta < -1e-4:
                overlaps.append({
                    "interval_index_before": i,
                    "depth_before": curr_bottom,
                    "depth_after": next_top,
                    "overlap_size_m": round(-delta, 3),
                    "facies_before": df.iloc[i]["Facies"],
                    "facies_after": df.iloc[i + 1]["Facies"],
                })

        duplicates_count = int(df.duplicated(subset=["Top", "Bottom"]).sum())

        # Facies breakdown
        facies_counts: Dict[str, int] = {}
        facies_thickness: Dict[str, float] = {}
        facies_props: Dict[str, float] = {}

        for facies_key in CANONICAL_FACIES_SCHEMA.keys():
            subset = df[df["Facies"] == facies_key]
            count = len(subset)
            thick = float(subset["Thickness"].sum()) if count > 0 else 0.0
            prop = float(thick / sum_thickness) if sum_thickness > 0 else 0.0

            facies_counts[facies_key] = count
            facies_thickness[facies_key] = round(thick, 3)
            facies_props[facies_key] = round(prop, 4)

        # Net-to-gross ratios
        sand_thick = facies_thickness.get("sand", 0.0)
        silt_thick = facies_thickness.get("silt", 0.0)
        coal_thick = facies_thickness.get("coal", 0.0)

        ntg_pure = float(sand_thick / sum_thickness) if sum_thickness > 0 else 0.0
        ntg_coarse = float((sand_thick + silt_thick) / sum_thickness) if sum_thickness > 0 else 0.0

        # Coal bed statistics
        coal_subset = df[df["Facies"] == "coal"]
        coal_beds_count = len(coal_subset)
        coal_mean_thick = float(coal_subset["Thickness"].mean()) if coal_beds_count > 0 else 0.0
        coal_max_thick = float(coal_subset["Thickness"].max()) if coal_beds_count > 0 else 0.0

        # Sandstone bed statistics
        sand_subset = df[df["Facies"] == "sand"]
        sand_beds_count = len(sand_subset)
        sand_min_thick = float(sand_subset["Thickness"].min()) if sand_beds_count > 0 else 0.0
        sand_median_thick = float(sand_subset["Thickness"].median()) if sand_beds_count > 0 else 0.0
        sand_mean_thick = float(sand_subset["Thickness"].mean()) if sand_beds_count > 0 else 0.0
        sand_max_thick = float(sand_subset["Thickness"].max()) if sand_beds_count > 0 else 0.0
        sand_std_thick = float(sand_subset["Thickness"].std(ddof=1)) if sand_beds_count > 1 else 0.0

        return {
            "litholog_id": litholog_id,
            "provenance_category": prov["provenance_category"],
            "provenance_description": prov["provenance_description"],
            "independent_validation": prov["independent_validation"],
            "coordinates_available": prov["coordinates_available"],
            "coordinates_status": prov["coordinates_status"],
            "group": prov["group"],
            "num_intervals": len(df),
            "min_depth_m": min_depth,
            "max_depth_m": max_depth,
            "represented_span_m": represented_span,
            "sum_thickness_m": round(sum_thickness, 3),
            "span_thickness_diff_m": round(span_thickness_diff, 3),
            "has_gaps": len(gaps) > 0,
            "has_overlaps": len(overlaps) > 0,
            "gaps": gaps,
            "overlaps": overlaps,
            "duplicate_intervals": duplicates_count,
            "facies_counts": facies_counts,
            "facies_thickness_m": facies_thickness,
            "facies_proportions": facies_props,
            "ntg_pure": round(ntg_pure, 4),
            "ntg_coarse": round(ntg_coarse, 4),
            "coal_total_thickness_m": round(coal_thick, 3),
            "coal_beds_count": coal_beds_count,
            "coal_mean_bed_thickness_m": round(coal_mean_thick, 3),
            "coal_max_bed_thickness_m": round(coal_max_thick, 3),
            "sand_beds_count": sand_beds_count,
            "sand_min_bed_thickness_m": round(sand_min_thick, 3),
            "sand_median_bed_thickness_m": round(sand_median_thick, 3),
            "sand_mean_bed_thickness_m": round(sand_mean_thick, 3),
            "sand_max_bed_thickness_m": round(sand_max_thick, 3),
            "sand_std_bed_thickness_m": round(sand_std_thick, 3),
        }

    def discretize_litholog_1m(self, litholog_id: str) -> pd.DataFrame:
        """
        Discretizes a litholog into 1.0 m grid bins using midpoint interval evaluation.

        For each meter interval [d, d+1), evaluates the midpoint depth z = d + 0.5 m
        against the continuous stratigraphy. Tracks whether each cell is directly
        source-observed, assigned across an unrecorded gap, or resolved from overlapping intervals.
        """
        df = self.load_raw_litholog(litholog_id)
        min_depth = int(np.floor(df["Top"].min()))
        max_depth = int(np.ceil(df["Bottom"].max()))

        records = []
        for d in range(min_depth, max_depth):
            z_mid = d + 0.5
            # Find interval covering z_mid
            covering = df[(df["Top"] <= z_mid) & (df["Bottom"] > z_mid)]
            if len(covering) == 1:
                facies = covering.iloc[0]["Facies"]
                is_gap_filled = False
                is_overlap_resolved = False
                is_source_observed = True
            elif len(covering) > 1:
                # Overlap resolved by taking first valid covering interval
                facies = covering.iloc[0]["Facies"]
                is_gap_filled = False
                is_overlap_resolved = True
                is_source_observed = False
            else:
                # If z_mid falls in an unrecorded gap, nearest assignment
                covering_lo = df[df["Top"] <= z_mid]
                if len(covering_lo) > 0:
                    facies = covering_lo.iloc[-1]["Facies"]
                else:
                    facies = df.iloc[0]["Facies"]
                is_gap_filled = True
                is_overlap_resolved = False
                is_source_observed = False

            records.append({
                "litholog_id": litholog_id,
                "depth_m": float(d),
                "depth_mid_m": z_mid,
                "facies": facies,
                "facies_code": CANONICAL_FACIES_SCHEMA[facies]["code"],
                "is_source_observed": is_source_observed,
                "is_gap_filled": is_gap_filled,
                "is_overlap_resolved": is_overlap_resolved,
            })

        return pd.DataFrame(records)

    def build_embedded_bed_sequence(self, litholog_id: str) -> pd.DataFrame:
        """
        Constructs the sequence of consecutive distinct facies beds from the continuous
        stratigraphy, sorted in upward stratigraphic order (decreasing depth / bottom to top).

        Consecutive raw intervals of the same facies are merged into a single distinct bed
        to ensure the embedded sequence contains only genuine boundary-crossing transitions.

        Tracks whether transitions between beds cross unrecorded gaps or overlaps.
        """
        df = self.load_raw_litholog(litholog_id)
        # Sort descending by Bottom depth so iteration proceeds upward stratigraphically:
        # deepest bed (base of well) -> shallowest bed (top of well)
        df_sorted = df.sort_values(by="Bottom", ascending=False).reset_index(drop=True)

        beds: List[Dict[str, Any]] = []
        for _, r in df_sorted.iterrows():
            top = float(r["Top"])
            bot = float(r["Bottom"])
            fac = str(r["Facies"])
            code = int(CANONICAL_FACIES_SCHEMA[fac]["code"])

            if not beds:
                beds.append({
                    "bed_idx": 0,
                    "facies": fac,
                    "facies_code": code,
                    "top_m": top,
                    "bottom_m": bot,
                    "thickness_m": round(bot - top, 3),
                    "intervals_merged": 1,
                    "crosses_gap": False,
                    "crosses_overlap": False,
                })
            else:
                prev = beds[-1]
                prev_top = float(prev["top_m"])
                gap = prev_top - bot       # > 0 means unobserved interval between beds
                overlap = bot - prev_top   # > 0 means intervals overlap

                if fac == prev["facies"]:
                    # Consecutive intervals of identical facies: merge into single distinct bed
                    prev["top_m"] = min(prev["top_m"], top)
                    prev["thickness_m"] = round(prev["thickness_m"] + (bot - top), 3)
                    prev["intervals_merged"] += 1
                    if gap > 1e-4:
                        prev["crosses_gap"] = True
                else:
                    crosses_gap = bool(gap > 1e-4)
                    crosses_overlap = bool(overlap > 1e-4)
                    beds.append({
                        "bed_idx": len(beds),
                        "facies": fac,
                        "facies_code": code,
                        "top_m": top,
                        "bottom_m": bot,
                        "thickness_m": round(bot - top, 3),
                        "intervals_merged": 1,
                        "crosses_gap": crosses_gap,
                        "crosses_overlap": crosses_overlap,
                    })

        return pd.DataFrame(beds)

    def compare_continuous_vs_discretized(self, litholog_id: str) -> Dict[str, Any]:
        """Quantifies the difference between continuous interval and 1m discretized facies proportions."""
        insp = self.inspect_litholog(litholog_id)
        df_disc = self.discretize_litholog_1m(litholog_id)

        n_samples = len(df_disc)
        disc_props: Dict[str, float] = {}
        delta_props: Dict[str, float] = {}

        for facies_key in CANONICAL_FACIES_SCHEMA.keys():
            disc_count = (df_disc["facies"] == facies_key).sum()
            disc_prop = float(disc_count / n_samples) if n_samples > 0 else 0.0
            cont_prop = insp["facies_proportions"][facies_key]

            disc_props[facies_key] = round(disc_prop, 4)
            delta_props[facies_key] = round(disc_prop - cont_prop, 4)

        return {
            "litholog_id": litholog_id,
            "continuous_thickness_m": insp["sum_thickness_m"],
            "discretized_samples_1m": n_samples,
            "continuous_proportions": insp["facies_proportions"],
            "discretized_proportions": disc_props,
            "delta_proportions (disc - cont)": delta_props,
            "continuous_ntg_pure": insp["ntg_pure"],
            "discretized_ntg_pure": round(disc_props.get("sand", 0.0), 4),
            "delta_ntg_pure": round(disc_props.get("sand", 0.0) - insp["ntg_pure"], 4),
        }

    def compute_initial_state_distribution(
        self,
        litholog_ids: List[str],
        embedded: bool = False,
        use_discretized: bool = True,
        smoothing_alpha: float = 0.1,
    ) -> np.ndarray:
        """
        Computes the empirical basal (initial state) distribution across a set of training lithologs
        in upward stratigraphic order (the facies at the base / deepest point of each well),
        regularized with Laplace smoothing.
        """
        K = len(CANONICAL_FACIES_SCHEMA)
        counts = np.zeros(K, dtype=np.float64)
        for lid in litholog_ids:
            if embedded and not use_discretized:
                beds_df = self.build_embedded_bed_sequence(lid)
                if len(beds_df) > 0:
                    base_code = int(beds_df.iloc[0]["facies_code"])
                    counts[base_code] += 1.0
            else:
                disc_df = self.discretize_litholog_1m(lid)
                if len(disc_df) > 0:
                    sorted_disc = disc_df.sort_values(by="depth_m", ascending=False)
                    base_code = int(sorted_disc.iloc[0]["facies_code"])
                    counts[base_code] += 1.0

        total_wells = len(litholog_ids)
        p_base = (counts + smoothing_alpha) / (total_wells + K * smoothing_alpha)
        return p_base / p_base.sum()

    def compute_vertical_transitions(
        self,
        litholog_ids: List[str],
        embedded: bool = True,
        use_discretized: bool = False,
        smoothing_alpha: float = 0.1,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Estimates the vertical Markov transition probability matrix across a set of lithologs.

        Transitions are strictly tallied in upward stratigraphic order (decreasing depth).

        Args:
            litholog_ids: List of litholog IDs to fit.
            embedded: If True, suppresses self-transitions (P_ii = 0).
            use_discretized: If True, uses 1m discretized grids; if False, uses distinct continuous beds.
            smoothing_alpha: Additive Laplace smoothing prior (default 0.1).

        Returns:
            P: Row-stochastic K x K transition probability matrix.
            N: Raw K x K transition count matrix.
            pi: Stationary distribution vector (1 x K).
        """
        K = len(CANONICAL_FACIES_SCHEMA)
        N = np.zeros((K, K), dtype=np.float64)

        for lid in litholog_ids:
            if use_discretized:
                df = self.discretize_litholog_1m(lid)
                # Sort descending by depth to step upward stratigraphically
                df_sorted = df.sort_values(by="depth_m", ascending=False)
                codes = df_sorted["facies_code"].to_numpy()
                if len(codes) < 2:
                    continue
                for u, v in zip(codes[:-1], codes[1:]):
                    if embedded and u == v:
                        continue
                    N[u, v] += 1.0
            else:
                if embedded:
                    # Distinct bed sequence: consecutive identical facies already merged
                    beds_df = self.build_embedded_bed_sequence(lid)
                    codes = beds_df["facies_code"].to_numpy()
                    if len(codes) < 2:
                        continue
                    for u, v in zip(codes[:-1], codes[1:]):
                        N[u, v] += 1.0
                else:
                    df = self.load_raw_litholog(lid)
                    df_sorted = df.sort_values(by="Bottom", ascending=False)
                    codes = df_sorted["Facies"].map(lambda f: CANONICAL_FACIES_SCHEMA[f]["code"]).to_numpy()
                    if len(codes) < 2:
                        continue
                    for u, v in zip(codes[:-1], codes[1:]):
                        N[u, v] += 1.0

        # Construct transition probability matrix P
        P = np.zeros((K, K), dtype=np.float64)
        for i in range(K):
            row = N[i, :].copy()
            if embedded:
                row[i] = 0.0
                row_sum = np.sum(row)
                if row_sum == 0 and smoothing_alpha == 0:
                    P[i, :] = 1.0 / (K - 1)
                    P[i, i] = 0.0
                else:
                    P[i, :] = (row + smoothing_alpha) / (row_sum + (K - 1) * smoothing_alpha)
                    P[i, i] = 0.0
            else:
                row_sum = np.sum(row)
                if row_sum == 0 and smoothing_alpha == 0:
                    P[i, :] = 1.0 / K
                else:
                    P[i, :] = (row + smoothing_alpha) / (row_sum + K * smoothing_alpha)

        # Renormalize to ensure exact row stochasticity
        P = P / P.sum(axis=1, keepdims=True)

        # Compute stationary distribution pi
        pi = self._compute_stationary_vector(P)

        return P, N, pi

    @staticmethod
    def _compute_stationary_vector(P: np.ndarray) -> np.ndarray:
        """Solves invariant left eigenvector pi * P = pi, sum(pi) = 1.0."""
        K = P.shape[0]
        try:
            eigvals, eigvecs = np.linalg.eig(P.T)
            # Find eigenvalue closest to 1.0
            idx = np.argmin(np.abs(eigvals - 1.0))
            pi = np.real(eigvecs[:, idx])
            if np.all(pi < 0):
                pi = -pi
            if np.any(pi < 0):
                pi = np.clip(pi, 0.0, None)
            if pi.sum() > 0:
                pi = pi / pi.sum()
            else:
                pi = np.ones(K) / K
        except Exception:
            pi = np.ones(K) / K
        return pi
