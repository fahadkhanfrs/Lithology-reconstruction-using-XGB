"""
SMALT Visualization Module: Transition Matrices, Succession Graphs, and Geostatistical Diagnostics.
"""

from smalt.viz.markov_viz import (
    plot_transition_matrix,
    plot_facies_succession_network,
    plot_stationary_vs_empirical,
)

__all__ = [
    "plot_transition_matrix",
    "plot_facies_succession_network",
    "plot_stationary_vs_empirical",
]
