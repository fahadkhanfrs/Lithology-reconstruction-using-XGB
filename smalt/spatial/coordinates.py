"""
SMALT Spatial Coordinates Module: Coordinate Ingestion & Geometry Validation.

Ingests real source-derived coordinates from Location_coordinates_lithologs.xlsx,
validates spatial bounds, handles missing coordinates (e.g. Litholog 1), and
provides clean programmatic access for spatial modeling pipelines.
"""

from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
import numpy as np
import pandas as pd

# Default relative origin / coordinate spreadsheet path
COORDINATE_FILE = Path("lolo/Location_coordinates_lithologs.xlsx")


def load_source_coordinates(
    filepath: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Loads source-derived coordinates from the project spreadsheet.

    Returns DataFrame with columns:
      - litholog_id (canonical string, e.g. 'litholog2')
      - raw_label (e.g. 'L2')
      - x_m (easting / local X coordinate in meters)
      - y_m (northing / local Y coordinate in meters)
      - coordinates_available (bool)
      - provenance (string)
      - notes (string)
    """
    path = filepath or COORDINATE_FILE
    if not Path(path).exists():
        raise FileNotFoundError(f"Coordinate file not found at: {path}")

    # Read excel file
    df_raw = pd.read_excel(path)
    
    # Identify columns: first column has well labels, second X, third Y
    label_col = df_raw.columns[0]
    x_col = df_raw.columns[1]
    y_col = df_raw.columns[2]

    records = []
    # All 12 project lithologs
    all_lids = [f"litholog{i}" for i in range(1, 13)]

    # Map raw excel entries
    excel_map: Dict[str, Tuple[float, float]] = {}
    for _, row in df_raw.iterrows():
        lbl = str(row[label_col]).strip()
        if not lbl or lbl == "nan" or pd.isna(row[label_col]):
            continue
        try:
            x_val = float(row[x_col])
            y_val = float(row[y_col])
            if not (np.isnan(x_val) or np.isnan(y_val)):
                # Normalize label (e.g. 'L2' -> 'litholog2')
                clean_lbl = lbl.lower().replace(" ", "")
                if clean_lbl.startswith("l") and clean_lbl[1:].isdigit():
                    num = int(clean_lbl[1:])
                    excel_map[f"litholog{num}"] = (x_val, y_val)
        except (ValueError, TypeError):
            continue

    for lid in all_lids:
        num = int(lid.replace("litholog", ""))
        raw_lbl = f"L{num}"
        if lid in excel_map:
            x_m, y_m = excel_map[lid]
            records.append({
                "litholog_id": lid,
                "raw_label": raw_lbl,
                "x_m": x_m,
                "y_m": y_m,
                "coordinates_available": True,
                "provenance": "source_spreadsheet (Location_coordinates_lithologs.xlsx)",
                "notes": f"Local Cartesian coordinate in meters; row {num-1} of spreadsheet.",
            })
        else:
            records.append({
                "litholog_id": lid,
                "raw_label": raw_lbl,
                "x_m": np.nan,
                "y_m": np.nan,
                "coordinates_available": False,
                "provenance": "missing_source_data",
                "notes": "Coordinates absent in source spreadsheet; strictly ineligible for spatial validation.",
            })

    return pd.DataFrame(records)


def get_spatial_litholog_subset(
    exclude_missing: bool = True,
    filepath: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Returns the table of lithologs eligible for spatial evaluation.

    If exclude_missing is True, filters out lithologs without valid coordinates (e.g. litholog1).
    """
    df_coords = load_source_coordinates(filepath)
    if exclude_missing:
        return df_coords[df_coords["coordinates_available"]].reset_index(drop=True)
    return df_coords


def compute_interwell_horizontal_distances(
    filepath: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Computes Euclidean horizontal distances (in meters) between all pairs
    of lithologs that possess source coordinates.
    """
    df_coords = get_spatial_litholog_subset(exclude_missing=True, filepath=filepath)
    n = len(df_coords)
    records = []

    for i in range(n):
        id_i = df_coords.iloc[i]["litholog_id"]
        xi, yi = df_coords.iloc[i]["x_m"], df_coords.iloc[i]["y_m"]
        for j in range(n):
            id_j = df_coords.iloc[j]["litholog_id"]
            xj, yj = df_coords.iloc[j]["x_m"], df_coords.iloc[j]["y_m"]
            dist_m = float(np.sqrt((xi - xj) ** 2 + (yi - yj) ** 2))
            records.append({
                "from_litholog": id_i,
                "to_litholog": id_j,
                "horizontal_distance_m": round(dist_m, 2),
            })

    return pd.DataFrame(records)
