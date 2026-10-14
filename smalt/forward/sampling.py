"""
SMALT Synthetic Borehole Sampling Pipeline (Sprint K).

Extracts sparse 1D synthetic borehole observations from complete 2D realizations.
Strictly separates sampled data from withheld subsurface ground truth:
- Downstream reconstruction baselines receive ONLY permitted borehole observations.
- Complete ground truth remains strictly withheld.
- Synthetic coordinates are quarantined with explicit 'synth_' prefixes to prevent
  accidental contamination of real-site spatial modeling tables.
"""

from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd

from smalt.descriptive.analyzer import CANONICAL_FACIES_SCHEMA
from smalt.forward.simulator import SyntheticRealization2D, CANONICAL_CODE_TO_NAME


class SyntheticBoreholeDataset:
    """
    Container for sparse borehole observations extracted from a synthetic realization.

    Guarantees:
    - Contains ONLY the 1D column data extracted at permitted borehole coordinates.
    - Withheld 2D ground truth, full grid arrays, and geobody bounding boxes are
      strictly excluded to prevent target or architectural leakage.
    - Coordinates are explicitly tagged with 'synthetic_' / 'synth_' to ensure they
      never contaminate real-site spatial databases.
    """

    def __init__(
        self,
        data: pd.DataFrame,
        borehole_ids: List[str],
        x_locations_m: List[float],
        z_coords_m: np.ndarray,
        spacing_m: Optional[float] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.data = data.copy()
        self.borehole_ids = list(borehole_ids)
        self.x_locations_m = list(x_locations_m)
        self.z_coords_m = z_coords_m.copy()
        self.spacing_m = spacing_m
        self.metadata = metadata or {}

    @property
    def num_boreholes(self) -> int:
        return len(self.borehole_ids)

    @property
    def num_observations(self) -> int:
        return len(self.data)

    def get_borehole(self, borehole_id: str) -> pd.DataFrame:
        """Returns observation DataFrame for a specific synthetic borehole."""
        if borehole_id not in self.borehole_ids:
            raise KeyError(f"Borehole '{borehole_id}' not found in dataset. Available: {self.borehole_ids}")
        return self.data[self.data["litholog_id"] == borehole_id].copy().reset_index(drop=True)

    def to_csv(self, filepath: str) -> None:
        """Exports sparse synthetic borehole observations to CSV."""
        self.data.to_csv(filepath, index=False)


class SyntheticBoreholeSampler:
    """
    Pipeline for extracting sparse 1D borehole observations from a SyntheticRealization2D.
    """

    def __init__(self, realization: SyntheticRealization2D) -> None:
        if not isinstance(realization, SyntheticRealization2D):
            raise TypeError("realization must be an instance of SyntheticRealization2D.")
        self.realization = realization

    def sample_regular(
        self,
        spacing_m: float = 400.0,
        offset_m: float = 0.0,
    ) -> SyntheticBoreholeDataset:
        """
        Extracts boreholes at regular lateral intervals across the synthetic domain.

        Args:
            spacing_m: Distance between adjacent boreholes in meters (must be > 0).
            offset_m: Initial horizontal offset in meters (0 <= offset < spacing_m).

        Returns:
            SyntheticBoreholeDataset containing sparse observations only.
        """
        if spacing_m <= 0.0:
            raise ValueError(f"spacing_m must be positive, got {spacing_m}")
        if offset_m < 0.0 or offset_m >= self.realization.config.grid.x_extent_m:
            raise ValueError(
                f"offset_m must be in [0, {self.realization.config.grid.x_extent_m}), got {offset_m}"
            )

        x_max = self.realization.config.grid.x_extent_m
        x_locs = []
        curr_x = offset_m
        while curr_x < x_max:
            x_locs.append(curr_x)
            curr_x += spacing_m

        if not x_locs:
            raise ValueError(f"Spacing {spacing_m}m with offset {offset_m}m yielded 0 boreholes.")

        return self.sample_locations(x_locs, spacing_m=spacing_m)

    def sample_locations(
        self,
        x_locations_m: List[float],
        spacing_m: Optional[float] = None,
    ) -> SyntheticBoreholeDataset:
        """
        Extracts 1D vertical borehole columns at explicitly specified lateral coordinates.

        Args:
            x_locations_m: List of horizontal coordinates in meters within domain [0, X_max].
            spacing_m: Optional nominal spacing value to attach as metadata.

        Returns:
            SyntheticBoreholeDataset containing sparse observations only.
        """
        if not x_locations_m:
            raise ValueError("x_locations_m cannot be empty.")

        x_max = self.realization.config.grid.x_extent_m
        z_max = self.realization.config.grid.z_extent_m
        dx = self.realization.config.grid.dx_m
        nx = self.realization.nx
        nz = self.realization.nz

        for x in x_locations_m:
            if not (0.0 <= x <= x_max):
                raise ValueError(f"Synthetic borehole location {x} m is outside domain [0, {x_max}].")

        records: List[Dict[str, Any]] = []
        borehole_ids: List[str] = []

        for idx, x_val in enumerate(x_locations_m):
            b_id = f"synth_well_{idx}"
            borehole_ids.append(b_id)

            # Map x_val to nearest grid column
            col_idx = int(np.clip(np.round(x_val / dx), 0, nx - 1))
            facies_column = self.realization.grid[:, col_idx]

            # In upward stratigraphic orientation: row 0 is base (z=0), row nz-1 is top (z=z_max)
            for r_idx in range(nz):
                z_height = float(self.realization.z_coords_m[r_idx])
                # Downward measured depth convention: depth = z_max - z_height
                depth_down = float(z_max - z_height)
                f_code = int(facies_column[r_idx])
                f_name = CANONICAL_CODE_TO_NAME[f_code]

                records.append({
                    "litholog_id": b_id,
                    "synthetic_x_m": round(x_val, 2),
                    "depth_m": round(depth_down, 2),
                    "strat_height_m": round(z_height, 2),
                    "facies_code": f_code,
                    "facies_name": f_name,
                    "grid_col_idx": col_idx,
                    "grid_row_idx": r_idx,
                })

        df_out = pd.DataFrame(records)

        metadata = {
            "source_seed": self.realization.seed,
            "domain_x_extent_m": x_max,
            "domain_z_extent_m": z_max,
            "num_boreholes": len(borehole_ids),
            "sampling_type": "regular" if spacing_m is not None else "custom_locations",
            "spacing_m": spacing_m,
            "coordinate_system": "project_local_synthetic_quarantined",
            "is_synthetic": True,
        }

        return SyntheticBoreholeDataset(
            data=df_out,
            borehole_ids=borehole_ids,
            x_locations_m=x_locations_m,
            z_coords_m=self.realization.z_coords_m,
            spacing_m=spacing_m,
            metadata=metadata,
        )
