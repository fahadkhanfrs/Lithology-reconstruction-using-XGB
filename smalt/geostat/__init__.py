"""
SMALT Geostatistical Core Module: 1D Markov Chains, Transition Probabilities, and Spatial Geostatistics.
"""

from smalt.geostat.markov import StratigraphicMarkovChain
from smalt.geostat.spatial_markov import (
    SpatialMarkovTransitionAnalyzer,
    SpatialMarkovPredictor,
    DEFAULT_LATERAL_FACIES_LENGTHS_M,
)

__all__ = [
    "StratigraphicMarkovChain",
    "SpatialMarkovTransitionAnalyzer",
    "SpatialMarkovPredictor",
    "DEFAULT_LATERAL_FACIES_LENGTHS_M",
]

