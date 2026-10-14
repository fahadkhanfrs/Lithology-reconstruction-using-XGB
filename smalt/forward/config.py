"""
SMALT Forward Simulation Configuration Module (Sprint K).

Defines configuration dataclasses, parameter provenance tags, and presets
for reproducible 2D unconditioned fluvial forward simulation.

Parameter Provenance Categories:
- EMPIRICAL: Estimated from digitized lithologs (e.g. bed thicknesses, vertical Markov stats).
- LITERATURE_PRIOR: Derived from published literature (e.g. Sahoo et al., 2016 W/T ratios).
- MODELING_ASSUMPTION: Computational and geometric modeling choices (e.g. grid resolution, scour shape).
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional
import numpy as np

from smalt.descriptive.analyzer import CANONICAL_FACIES_SCHEMA


class ParameterProvenance(str, Enum):
    """Categorization of parameter evidence source."""
    EMPIRICAL = "EMPIRICAL"
    LITERATURE_PRIOR = "LITERATURE_PRIOR"
    MODELING_ASSUMPTION = "MODELING_ASSUMPTION"


@dataclass
class GridConfig:
    """
    2D Simulation Grid Domain Configuration.
    Stratigraphic height Z is oriented upward from base (Z=0) to top (Z=z_extent_m).
    """
    x_extent_m: float = 2000.0  # MODELING_ASSUMPTION: 2 km lateral transect
    z_extent_m: float = 50.0    # MODELING_ASSUMPTION: 50 m stratigraphic column
    dx_m: float = 10.0          # MODELING_ASSUMPTION: 10 m lateral resolution
    dz_m: float = 0.5           # MODELING_ASSUMPTION: 0.5 m vertical resolution

    @property
    def nx(self) -> int:
        return int(np.round(self.x_extent_m / self.dx_m))

    @property
    def nz(self) -> int:
        return int(np.round(self.z_extent_m / self.dz_m))

    def validate(self) -> None:
        if self.x_extent_m <= 0.0 or self.z_extent_m <= 0.0:
            raise ValueError(f"Grid extents must be positive. Got x={self.x_extent_m}, z={self.z_extent_m}")
        if self.dx_m <= 0.0 or self.dz_m <= 0.0:
            raise ValueError(f"Grid steps must be positive. Got dx={self.dx_m}, dz={self.dz_m}")
        if self.nx < 2 or self.nz < 2:
            raise ValueError(f"Grid dimensions must be at least 2x2. Got nx={self.nx}, nz={self.nz}")


@dataclass
class ChannelBeltConfig:
    """
    Channel Belt Geometry and Internal Facies Configuration.
    Channel thickness sampled from LogNormal distribution; width scaled via W/T ratio.
    """
    # Thickness distribution parameters (LogNormal: mu, sigma)
    # Default mu=1.7870, sigma=0.8823 is the MLE fit to 34 empirical Upstream channel beds
    thickness_mu: float = 1.7870          # EMPIRICAL: Upstream MLE log-thickness mean
    thickness_sigma: float = 0.8823       # EMPIRICAL: Upstream MLE log-thickness std
    thickness_min_m: float = 1.5          # MODELING_ASSUMPTION: Minimum channel thickness bound
    thickness_max_m: float = 25.0         # MODELING_ASSUMPTION: Maximum channel thickness bound

    # Width-to-thickness ratio (W/T)
    # Sahoo et al. (2016) reported mean W/T ~ 35 for fluvial channel belts in Ferron / Blackhawk
    width_to_thickness_ratio: float = 35.0  # LITERATURE_PRIOR: Sahoo et al. (2016)
    width_to_thickness_std: float = 5.0     # LITERATURE_PRIOR: Architectural variability

    # Internal facies assignment
    core_facies_code: int = 0             # EMPIRICAL: sand (Channel Sandstone, Facies 1)
    basal_facies_code: int = 1            # EMPIRICAL: p_sand (Planar Sandstone, Facies 2)
    basal_fraction: float = 0.20          # EMPIRICAL: Basal planar sand fraction (~20% of channel fill)
    include_basal_planar_sand: bool = True # Configurable modeling option

    def validate(self) -> None:
        if self.thickness_min_m <= 0.0 or self.thickness_max_m <= self.thickness_min_m:
            raise ValueError(f"Invalid thickness bounds: min={self.thickness_min_m}, max={self.thickness_max_m}")
        if self.thickness_sigma <= 0.0:
            raise ValueError(f"thickness_sigma must be positive. Got {self.thickness_sigma}")
        if self.width_to_thickness_ratio <= 0.0:
            raise ValueError(f"width_to_thickness_ratio must be positive. Got {self.width_to_thickness_ratio}")
        if not (0.0 <= self.basal_fraction <= 1.0):
            raise ValueError(f"basal_fraction must be in [0, 1]. Got {self.basal_fraction}")


@dataclass
class CrevasseSplayConfig:
    """
    Crevasse Splay Flanking Geometry Configuration.
    Extends laterally from channel margins into overbank floodplain.
    """
    enabled: bool = True
    # Sahoo et al. (2016) reported splay lateral widths of 10 to 130 m
    width_min_m: float = 20.0             # LITERATURE_PRIOR: Sahoo et al. (2016)
    width_max_m: float = 100.0            # LITERATURE_PRIOR: Sahoo et al. (2016)
    # Empirical bed thickness for ripples (mean 2.54 m) and carbon_mud (1.0 m)
    thickness_min_m: float = 0.5          # EMPIRICAL: Minimum splay thickness
    thickness_max_m: float = 2.5          # EMPIRICAL: Maximum splay thickness

    lobe_facies_code: int = 2             # EMPIRICAL: ripples (Rippled Heterolithics, Facies 3)
    margin_facies_code: int = 3           # EMPIRICAL: carbon_mud (Carbonaceous Mudstone, Facies 4)
    margin_fraction: float = 0.25         # EMPIRICAL: Distal muddy margin fraction

    def validate(self) -> None:
        if self.width_min_m <= 0.0 or self.width_max_m <= self.width_min_m:
            raise ValueError(f"Invalid splay width bounds: min={self.width_min_m}, max={self.width_max_m}")
        if self.thickness_min_m <= 0.0 or self.thickness_max_m <= self.thickness_min_m:
            raise ValueError(f"Invalid splay thickness bounds: min={self.thickness_min_m}, max={self.thickness_max_m}")
        if not (0.0 <= self.margin_fraction <= 1.0):
            raise ValueError(f"margin_fraction must be in [0, 1]. Got {self.margin_fraction}")


@dataclass
class CoalHorizonConfig:
    """
    Regional Coal Horizon Configuration.
    Represents extensive mire/swamp accumulation acting as chronostratigraphic markers.
    """
    enabled: bool = True
    num_horizons: int = 2                 # EMPIRICAL: Upstream coal frequency (~2 beds per 50m)
    thickness_m: float = 1.0              # EMPIRICAL: Upstream coal thickness (median and mode = 1.0 m)
    facies_code: int = 4                  # EMPIRICAL: coal (Coal, Facies 5)
    relative_elevations: List[float] = field(default_factory=lambda: [0.25, 0.75])

    def validate(self) -> None:
        if self.enabled:
            if self.num_horizons < 0:
                raise ValueError(f"num_horizons must be non-negative. Got {self.num_horizons}")
            if self.thickness_m <= 0.0:
                raise ValueError(f"Coal thickness must be positive. Got {self.thickness_m}")
            for z_rel in self.relative_elevations:
                if not (0.0 <= z_rel <= 1.0):
                    raise ValueError(f"Relative elevation must be in [0, 1]. Got {z_rel}")


@dataclass
class BackgroundConfig:
    """
    Background Floodplain Deposition Configuration.
    """
    primary_facies_code: int = 5          # EMPIRICAL: mud (Overbank Mudstone, Facies 6)

    def validate(self) -> None:
        if not (0 <= self.primary_facies_code <= 5):
            raise ValueError(f"Invalid background facies code: {self.primary_facies_code}")


@dataclass
class ForwardSimulatorConfig:
    """
    Master 2D Fluvial Forward Simulator Configuration.
    """
    preset_name: str = "upstream"
    grid: GridConfig = field(default_factory=GridConfig)
    channel: ChannelBeltConfig = field(default_factory=ChannelBeltConfig)
    splay: CrevasseSplayConfig = field(default_factory=CrevasseSplayConfig)
    coal: CoalHorizonConfig = field(default_factory=CoalHorizonConfig)
    background: BackgroundConfig = field(default_factory=BackgroundConfig)

    # Net-to-Gross target and stopping criteria
    target_net_to_gross: float = 0.45     # EMPIRICAL: Upstream reference (channel 42.9%, pure 50.5%)
    max_channel_attempts: int = 50        # MODELING_ASSUMPTION: Max placement attempts
    tolerance_net_to_gross: float = 0.08  # MODELING_ASSUMPTION: Stopping window

    def validate(self) -> None:
        self.grid.validate()
        self.channel.validate()
        self.splay.validate()
        self.coal.validate()
        self.background.validate()
        if not (0.0 < self.target_net_to_gross < 1.0):
            raise ValueError(f"target_net_to_gross must be in (0, 1). Got {self.target_net_to_gross}")
        if self.max_channel_attempts < 1:
            raise ValueError(f"max_channel_attempts must be >= 1. Got {self.max_channel_attempts}")


def get_preset_config(preset: str = "upstream") -> ForwardSimulatorConfig:
    """
    Returns a validated ForwardSimulatorConfig based on empirical group calibration.

    Presets:
    - 'upstream': Default. Calibrated to 8 Upstream logs (L2-L8, L10) and Sahoo et al. (2016) priors.
    - 'downstream_outcrop': Calibrated to 2 Downstream outcrop logs (L9, L11; higher ripples, wider splays).
    - 'pooled': Calibrated to all 12 lithologs (L1-L12).
    """
    preset_clean = preset.strip().lower()

    if preset_clean == "upstream":
        return ForwardSimulatorConfig(
            preset_name="upstream",
            grid=GridConfig(x_extent_m=2000.0, z_extent_m=50.0, dx_m=10.0, dz_m=0.5),
            channel=ChannelBeltConfig(
                thickness_mu=1.7870,
                thickness_sigma=0.8823,
                thickness_min_m=1.5,
                thickness_max_m=25.0,
                width_to_thickness_ratio=35.0,
                width_to_thickness_std=5.0,
                core_facies_code=0,
                basal_facies_code=1,
                basal_fraction=0.20,
                include_basal_planar_sand=True,
            ),
            splay=CrevasseSplayConfig(
                enabled=True,
                width_min_m=20.0,
                width_max_m=100.0,
                thickness_min_m=0.5,
                thickness_max_m=2.5,
                lobe_facies_code=2,
                margin_facies_code=3,
                margin_fraction=0.25,
            ),
            coal=CoalHorizonConfig(
                enabled=True,
                num_horizons=2,
                thickness_m=1.0,
                facies_code=4,
                relative_elevations=[0.25, 0.75],
            ),
            background=BackgroundConfig(primary_facies_code=5),
            target_net_to_gross=0.45,
            max_channel_attempts=50,
        )

    elif preset_clean == "downstream_outcrop":
        return ForwardSimulatorConfig(
            preset_name="downstream_outcrop",
            grid=GridConfig(x_extent_m=2000.0, z_extent_m=50.0, dx_m=10.0, dz_m=0.5),
            channel=ChannelBeltConfig(
                thickness_mu=2.20,      # Higher mean thickness in downstream outcrop
                thickness_sigma=0.60,
                thickness_min_m=2.0,
                thickness_max_m=20.0,
                width_to_thickness_ratio=35.0,
                width_to_thickness_std=5.0,
                core_facies_code=0,
                basal_facies_code=1,
                basal_fraction=0.18,
                include_basal_planar_sand=True,
            ),
            splay=CrevasseSplayConfig(
                enabled=True,
                width_min_m=40.0,       # Wider splays in downstream reach
                width_max_m=130.0,
                thickness_min_m=1.0,
                thickness_max_m=4.5,
                lobe_facies_code=2,
                margin_facies_code=3,
                margin_fraction=0.20,
            ),
            coal=CoalHorizonConfig(
                enabled=True,
                num_horizons=1,
                thickness_m=1.5,
                facies_code=4,
                relative_elevations=[0.50],
            ),
            background=BackgroundConfig(primary_facies_code=5),
            target_net_to_gross=0.42,
            max_channel_attempts=50,
        )

    elif preset_clean == "pooled":
        return ForwardSimulatorConfig(
            preset_name="pooled",
            grid=GridConfig(x_extent_m=2000.0, z_extent_m=50.0, dx_m=10.0, dz_m=0.5),
            channel=ChannelBeltConfig(
                thickness_mu=1.70,
                thickness_sigma=0.85,
                thickness_min_m=1.5,
                thickness_max_m=25.0,
                width_to_thickness_ratio=35.0,
                width_to_thickness_std=5.0,
                core_facies_code=0,
                basal_facies_code=1,
                basal_fraction=0.15,
                include_basal_planar_sand=True,
            ),
            splay=CrevasseSplayConfig(
                enabled=True,
                width_min_m=20.0,
                width_max_m=110.0,
                thickness_min_m=0.5,
                thickness_max_m=3.0,
                lobe_facies_code=2,
                margin_facies_code=3,
                margin_fraction=0.22,
            ),
            coal=CoalHorizonConfig(
                enabled=True,
                num_horizons=2,
                thickness_m=1.0,
                facies_code=4,
                relative_elevations=[0.30, 0.70],
            ),
            background=BackgroundConfig(primary_facies_code=5),
            target_net_to_gross=0.44,
            max_channel_attempts=50,
        )

    else:
        raise ValueError(f"Unknown preset '{preset}'. Expected one of ['upstream', 'downstream_outcrop', 'pooled'].")
