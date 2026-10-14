"""
SMALT Fluvial Forward Simulator Engine (Sprint K).

Implements unconditioned 2D process-based forward simulation of fluvial stratigraphy:
- Reproducible complete 2D facies grids with known facies everywhere.
- Channel belt ribbon bodies with parabolic erosional scours and W/T scaling.
- Flanking crevasse splay lobes with ripple heterolithics and carbonaceous mud margins.
- Regional chronostratigraphic coal horizons.
- Overbank floodplain mudstone background.
- Deterministic erosional truncation (superposition).
- 100% seed reproducibility.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from smalt.descriptive.analyzer import CANONICAL_FACIES_SCHEMA
from smalt.forward.config import (
    ForwardSimulatorConfig,
    GridConfig,
    get_preset_config,
)

CANONICAL_CODE_TO_NAME: Dict[int, str] = {
    0: "sand",
    1: "p_sand",
    2: "ripples",
    3: "carbon_mud",
    4: "coal",
    5: "mud",
}


class SyntheticRealization2D:
    """
    Encapsulates a completed synthetic 2D geological realization.

    Attributes:
        grid: (nz, nx) 2D array of canonical integer facies codes (0 to 5).
        x_coords_m: (nx,) 1D array of horizontal coordinates in meters.
        z_coords_m: (nz,) 1D array of stratigraphic height coordinates in meters (0 at base).
        config: ForwardSimulatorConfig instance used for generation.
        seed: Random seed used for generation.
        objects_metadata: List of dictionaries describing placed geological objects.
        facies_proportions: Realized cell/volume proportions of all 6 facies.
        realized_net_to_gross: Realized channel, pure sand, and total coarse Net-to-Gross.
        vertical_diagnostics: Vertical transition counts and run-length statistics.
    """

    def __init__(
        self,
        grid: np.ndarray,
        x_coords_m: np.ndarray,
        z_coords_m: np.ndarray,
        config: ForwardSimulatorConfig,
        seed: int,
        objects_metadata: List[Dict[str, Any]],
        facies_proportions: Dict[str, float],
        realized_net_to_gross: Dict[str, float],
        vertical_diagnostics: Dict[str, Any],
    ) -> None:
        self.grid = grid
        self.x_coords_m = x_coords_m
        self.z_coords_m = z_coords_m
        self.config = config
        self.seed = seed
        self.objects_metadata = objects_metadata
        self.facies_proportions = facies_proportions
        self.realized_net_to_gross = realized_net_to_gross
        self.vertical_diagnostics = vertical_diagnostics

    @property
    def shape(self) -> Tuple[int, int]:
        return self.grid.shape

    @property
    def nz(self) -> int:
        return self.grid.shape[0]

    @property
    def nx(self) -> int:
        return self.grid.shape[1]

    def get_column(self, x_m: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extracts 1D vertical facies column at coordinate x_m.

        Returns:
            (z_coords_m, facies_codes): 1D arrays sorted bottom-to-top (increasing z).
        """
        x_idx = int(np.clip(np.round(x_m / self.config.grid.dx_m), 0, self.nx - 1))
        col_facies = self.grid[:, x_idx]
        return self.z_coords_m.copy(), col_facies.copy()

    def to_dict(self) -> Dict[str, Any]:
        """Returns JSON-serializable summary metadata."""
        return {
            "seed": self.seed,
            "preset_name": self.config.preset_name,
            "grid_dimensions": {
                "nx": self.nx,
                "nz": self.nz,
                "x_extent_m": self.config.grid.x_extent_m,
                "z_extent_m": self.config.grid.z_extent_m,
                "dx_m": self.config.grid.dx_m,
                "dz_m": self.config.grid.dz_m,
            },
            "num_objects": len(self.objects_metadata),
            "facies_proportions": self.facies_proportions,
            "realized_net_to_gross": self.realized_net_to_gross,
            "target_net_to_gross": self.config.target_net_to_gross,
            "ntg_discrepancy": round(
                self.realized_net_to_gross["ntg_channel"] - self.config.target_net_to_gross, 4
            ),
        }

    def plot_cross_section(
        self,
        output_path: Union[str, Path],
        borehole_x_locations: Optional[List[float]] = None,
        title: Optional[str] = None,
    ) -> None:
        """
        Renders a cross-section of the 2D realization with canonical facies legend.
        """
        fig, ax = plt.subplots(figsize=(14, 5.5))

        # Build custom discrete colormap using canonical colors
        colors = [CANONICAL_FACIES_SCHEMA[CANONICAL_CODE_TO_NAME[c]]["color_hex"] for c in range(6)]
        cmap = matplotlib.colors.ListedColormap(colors)
        bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

        extent = [0.0, self.config.grid.x_extent_m, 0.0, self.config.grid.z_extent_m]
        im = ax.imshow(
            self.grid,
            origin="lower",
            extent=extent,
            cmap=cmap,
            norm=norm,
            aspect="auto",
            interpolation="nearest",
        )

        # Plot synthetic borehole vertical tracks if provided
        if borehole_x_locations:
            for x_loc in borehole_x_locations:
                ax.axvline(x=x_loc, color="black", linestyle="--", linewidth=1.2, alpha=0.85)
                ax.scatter(
                    [x_loc],
                    [self.config.grid.z_extent_m],
                    marker="v",
                    color="black",
                    s=60,
                    zorder=5,
                )

        ax.set_xlabel("Lateral Distance X (m)", fontsize=11)
        ax.set_ylabel("Stratigraphic Height Z (m)", fontsize=11)

        display_title = title or (
            f"Synthetic 2D Fluvial Forward Realization (Seed={self.seed}, Preset='{self.config.preset_name}')\n"
            f"N/G Channel={self.realized_net_to_gross['ntg_channel']:.1%}, "
            f"Pure Sand={self.realized_net_to_gross['ntg_pure']:.1%}, "
            f"Coarse={self.realized_net_to_gross['ntg_coarse']:.1%}"
        )
        ax.set_title(display_title, fontsize=12, fontweight="bold")

        # Custom Legend
        legend_patches = [
            mpatches.Patch(
                color=CANONICAL_FACIES_SCHEMA[CANONICAL_CODE_TO_NAME[c]]["color_hex"],
                label=f"{CANONICAL_CODE_TO_NAME[c]} ({CANONICAL_FACIES_SCHEMA[CANONICAL_CODE_TO_NAME[c]]['canonical_name']})",
            )
            for c in range(6)
        ]
        if borehole_x_locations:
            legend_patches.append(
                plt.Line2D([0], [0], color="black", linestyle="--", linewidth=1.2, label="Synthetic Boreholes")
            )

        ax.legend(handles=legend_patches, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, fontsize=9)
        plt.tight_layout()

        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_p, dpi=200)
        plt.close(fig)


class FluvialForwardSimulator2D:
    """
    2D Fluvial Forward Simulator with process-like ribbon channels and splays.
    """

    def __init__(self, config: Optional[ForwardSimulatorConfig] = None) -> None:
        self.config = config or get_preset_config("upstream")
        self.config.validate()

    def generate(self, seed: Optional[int] = None) -> SyntheticRealization2D:
        """
        Executes a complete, reproducible 2D realization.

        Args:
            seed: Non-negative integer random seed (default 42 if None).

        Returns:
            SyntheticRealization2D instance.
        """
        actual_seed = 42 if seed is None else int(seed)
        rng = np.random.default_rng(actual_seed)

        nx = self.config.grid.nx
        nz = self.config.grid.nz
        dx = self.config.grid.dx_m
        dz = self.config.grid.dz_m
        x_extent = self.config.grid.x_extent_m
        z_extent = self.config.grid.z_extent_m

        # Coordinate arrays (cell center locations)
        x_coords_m = np.linspace(0.5 * dx, x_extent - 0.5 * dx, nx)
        z_coords_m = np.linspace(0.5 * dz, z_extent - 0.5 * dz, nz)

        # 1. Initialize background: Overbank Mudstone (Facies 5)
        grid = np.full((nz, nx), self.config.background.primary_facies_code, dtype=np.int32)
        objects_meta: List[Dict[str, Any]] = []

        # 2. Regional Chronostratigraphic Coal Horizons (Facies 4)
        if self.config.coal.enabled:
            coal_thick = self.config.coal.thickness_m
            coal_facies = self.config.coal.facies_code
            for i, rel_z in enumerate(self.config.coal.relative_elevations):
                z_center = rel_z * z_extent
                z_bot = max(0.0, z_center - 0.5 * coal_thick)
                z_top = min(z_extent, z_center + 0.5 * coal_thick)

                z_mask = (z_coords_m >= z_bot) & (z_coords_m <= z_top)
                grid[z_mask, :] = coal_facies

                objects_meta.append({
                    "object_type": "coal_horizon",
                    "object_id": f"coal_{i}",
                    "z_center_m": z_center,
                    "thickness_m": coal_thick,
                    "facies_code": coal_facies,
                })

        # 3. Stochastic Channel Belt and Splay Placement
        channel_cfg = self.config.channel
        splay_cfg = self.config.splay
        target_ntg = self.config.target_net_to_gross
        max_attempts = self.config.max_channel_attempts

        placed_channels = 0
        attempts = 0

        while attempts < max_attempts:
            attempts += 1

            # Sample channel thickness from empirical LogNormal distribution
            t_raw = float(rng.lognormal(mean=channel_cfg.thickness_mu, sigma=channel_cfg.thickness_sigma))
            thickness = float(np.clip(t_raw, channel_cfg.thickness_min_m, channel_cfg.thickness_max_m))

            # Sample width-to-thickness ratio from literature prior (Sahoo et al., 2016)
            wt_ratio = float(rng.normal(channel_cfg.width_to_thickness_ratio, channel_cfg.width_to_thickness_std))
            wt_ratio = float(np.clip(wt_ratio, 15.0, 70.0))
            width = thickness * wt_ratio

            # Sample spatial location
            xc = float(rng.uniform(0.0, x_extent))
            # Base elevation: placed vertically within column
            z_base = float(rng.uniform(0.0, max(0.1, z_extent - thickness)))

            # Pre-compute current sand proportion
            current_sand_cells = np.sum(grid == channel_cfg.core_facies_code)
            current_ntg = current_sand_cells / (nx * nz)

            # Stopping condition: if target reached after at least 3 channels, stop
            if placed_channels >= 3 and current_ntg >= (target_ntg - self.config.tolerance_net_to_gross):
                break

            # A. Place Flanking Crevasse Splays before channel scour (so channel erodes splay root)
            if splay_cfg.enabled:
                w_splay = float(rng.uniform(splay_cfg.width_min_m, splay_cfg.width_max_m))
                t_splay = float(rng.uniform(splay_cfg.thickness_min_m, splay_cfg.thickness_max_m))
                z_splay_base = min(z_extent - t_splay, z_base + 0.5 * thickness)
                z_splay_top = z_splay_base + t_splay

                # Left and right splay wings
                for direction, x_attach in [(-1, xc - 0.5 * width), (1, xc + 0.5 * width)]:
                    x_start = x_attach
                    x_end = x_attach + direction * w_splay
                    x_min_splay = min(x_start, x_end)
                    x_max_splay = max(x_start, x_end)

                    col_mask = (x_coords_m >= x_min_splay) & (x_coords_m <= x_max_splay)
                    row_mask = (z_coords_m >= z_splay_base) & (z_coords_m <= z_splay_top)

                    if np.any(col_mask) and np.any(row_mask):
                        col_indices = np.where(col_mask)[0]
                        for c_idx in col_indices:
                            x_val = x_coords_m[c_idx]
                            # Distance from channel margin
                            dist_from_chan = abs(x_val - x_attach)
                            frac = dist_from_chan / w_splay
                            # Proximal lobe: ripples (2); distal margin: carbon_mud (3)
                            splay_code = (
                                splay_cfg.margin_facies_code
                                if frac >= (1.0 - splay_cfg.margin_fraction)
                                else splay_cfg.lobe_facies_code
                            )
                            for r_idx in np.where(row_mask)[0]:
                                grid[r_idx, c_idx] = splay_code

            # B. Place Channel Belt (Parabolic Scour with basal planar sand)
            x_min_chan = xc - 0.5 * width
            x_max_chan = xc + 0.5 * width
            col_mask_chan = (x_coords_m >= x_min_chan) & (x_coords_m <= x_max_chan)

            if np.any(col_mask_chan):
                col_indices = np.where(col_mask_chan)[0]
                for c_idx in col_indices:
                    x_val = x_coords_m[c_idx]
                    # Parabolic top scour profile: z_scour(x) = z_base + T * (1 - 4 * (x - xc)^2 / W^2)
                    rel_x = (x_val - xc) / (0.5 * width)
                    local_thick = thickness * max(0.0, 1.0 - rel_x**2)
                    z_top_local = z_base + local_thick

                    # Basal planar sand thickness: f_basal * local_thick
                    t_basal = channel_cfg.basal_fraction * local_thick if channel_cfg.include_basal_planar_sand else 0.0
                    z_basal_cutoff = z_base + t_basal

                    # Apply erosional superposition to grid column
                    for r_idx, z_val in enumerate(z_coords_m):
                        if z_base <= z_val <= z_top_local:
                            if z_val <= z_basal_cutoff:
                                grid[r_idx, c_idx] = channel_cfg.basal_facies_code  # p_sand (1)
                            else:
                                grid[r_idx, c_idx] = channel_cfg.core_facies_code   # sand (0)

                objects_meta.append({
                    "object_type": "channel_belt",
                    "object_id": f"channel_{placed_channels}",
                    "xc_m": xc,
                    "z_base_m": z_base,
                    "thickness_m": thickness,
                    "width_m": width,
                    "wt_ratio": wt_ratio,
                    "x_bounds": [max(0.0, x_min_chan), min(x_extent, x_max_chan)],
                    "z_bounds": [z_base, min(z_extent, z_base + thickness)],
                })
                placed_channels += 1

        # 4. Compute Realized Metrics
        total_cells = nx * nz
        facies_props: Dict[str, float] = {}
        for code, name in CANONICAL_CODE_TO_NAME.items():
            count = int(np.sum(grid == code))
            facies_props[name] = round(count / total_cells, 4)

        sand_p = facies_props["sand"]
        p_sand_p = facies_props["p_sand"]
        ripples_p = facies_props["ripples"]

        realized_ntg = {
            "ntg_channel": round(sand_p, 4),
            "ntg_pure": round(sand_p + p_sand_p, 4),
            "ntg_coarse": round(sand_p + p_sand_p + ripples_p, 4),
        }

        # 5. Compute Vertical Column Succession Diagnostics
        vert_diag = self._compute_vertical_diagnostics(grid, dz)

        return SyntheticRealization2D(
            grid=grid,
            x_coords_m=x_coords_m,
            z_coords_m=z_coords_m,
            config=self.config,
            seed=actual_seed,
            objects_metadata=objects_meta,
            facies_proportions=facies_props,
            realized_net_to_gross=realized_ntg,
            vertical_diagnostics=vert_diag,
        )

    def _compute_vertical_diagnostics(self, grid: np.ndarray, dz: float) -> Dict[str, Any]:
        """
        Computes 1D vertical transition counts and run lengths across all columns (bottom-to-top).
        """
        nz, nx = grid.shape
        N_trans = np.zeros((6, 6), dtype=np.int64)
        run_lengths: Dict[int, List[float]] = {c: [] for c in range(6)}

        for col_idx in range(nx):
            col = grid[:, col_idx]
            # Transitions tallied upward: row 0 (base) -> row nz-1 (top)
            for from_s, to_s in zip(col[:-1], col[1:]):
                N_trans[from_s, to_s] += 1

            # Run lengths (distinct bed thicknesses in column)
            current_facies = col[0]
            current_count = 1
            for f in col[1:]:
                if f == current_facies:
                    current_count += 1
                else:
                    run_lengths[current_facies].append(current_count * dz)
                    current_facies = f
                    current_count = 1
            run_lengths[current_facies].append(current_count * dz)

        mean_run_lengths: Dict[str, float] = {}
        for c, name in CANONICAL_CODE_TO_NAME.items():
            lens = run_lengths[c]
            mean_run_lengths[name] = round(float(np.mean(lens)), 2) if lens else 0.0

        return {
            "transition_count_matrix": N_trans.tolist(),
            "mean_bed_thickness_m": mean_run_lengths,
        }
