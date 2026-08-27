"""
SMALT Visualization Module: Diagnostic plots for 1D Markov Transition Analysis.

Publication-quality routines for:
1. Annotated transition probability heatmaps (plot_transition_matrix)
2. Directed facies succession networks (plot_facies_succession_network)
3. Stationary vs. empirical facies frequency diagnostics (plot_stationary_vs_empirical)
"""

from pathlib import Path
from typing import Optional, Union
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for headless execution
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

from smalt.geostat.markov import StratigraphicMarkovChain

# Standard geological facies color palette
FACIES_PALETTE = {
    0: "#2B2B2B",  # Coal (Dark Charcoal/Black)
    1: "#F1B738",  # Channel Sandstone (Warm Golden Amber)
    2: "#D97D27",  # Fine Sandstone / Splay (Terracotta/Orange)
    3: "#8DA365",  # Siltstone (Olive/Sage Green)
    4: "#4A6B82",  # Overbank Mudstone (Slate Blue/Gray)
}


def plot_transition_matrix(
    model: StratigraphicMarkovChain,
    title: str = "Stratigraphic Transition Probability Matrix",
    save_path: Optional[Union[str, Path]] = None,
) -> plt.Figure:
    """
    Renders an annotated heatmap of the transition probability matrix P.

    Args:
        model: Fitted StratigraphicMarkovChain instance.
        title: Figure title.
        save_path: Optional file path to save figure at 300 DPI.

    Returns:
        fig: Matplotlib Figure object.
    """
    if model.transition_matrix_ is None:
        raise ValueError("Model must be fitted before plotting transition matrix.")

    P = model.transition_matrix_
    K = model.num_classes
    labels = [model.facies_map.get(i, f"Facies {i}") for i in range(K)]

    fig, ax = plt.subplots(figsize=(8, 6.5), dpi=300)

    # Use a modern, high-contrast colormap
    im = ax.imshow(P, cmap="YlGnBu", vmin=0.0, vmax=np.max(P) if np.max(P) > 0 else 1.0)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Transition Probability $P_{ij}$", fontsize=11, fontweight="bold")
    cbar.ax.tick_params(labelsize=10)

    # Ticks and Labels
    ax.set_xticks(range(K))
    ax.set_yticks(range(K))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=10, fontweight="medium")
    ax.set_yticklabels(labels, fontsize=10, fontweight="medium")

    ax.set_xlabel("To Facies $S_t$ (Upper Bed)", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_ylabel("From Facies $S_{t-1}$ (Lower Bed)", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=15)

    # Cell Annotations
    thresh = np.max(P) / 2.0 if np.max(P) > 0 else 0.5
    for i in range(K):
        for j in range(K):
            val = P[i, j]
            count = model.count_matrix_[i, j] if model.count_matrix_ is not None else 0
            text_color = "white" if val > thresh else "black"
            ax.text(
                j,
                i,
                f"{val:.3f}\n(n={int(count)})",
                ha="center",
                va="center",
                color=text_color,
                fontsize=9,
                fontweight="bold" if val > 0.15 else "normal",
            )

    plt.tight_layout()

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return fig


def plot_facies_succession_network(
    model: StratigraphicMarkovChain,
    threshold: float = 0.10,
    save_path: Optional[Union[str, Path]] = None,
) -> plt.Figure:
    """
    Renders a directed transition network diagram showing dominant upward facies transitions.

    Args:
        model: Fitted StratigraphicMarkovChain instance.
        threshold: Minimum transition probability P_ij to render an edge.
        save_path: Optional file path to save figure at 300 DPI.

    Returns:
        fig: Matplotlib Figure object.
    """
    if model.transition_matrix_ is None:
        raise ValueError("Model must be fitted before plotting succession network.")

    P = model.transition_matrix_
    K = model.num_classes
    labels = [model.facies_map.get(i, f"State {i}") for i in range(K)]

    fig, ax = plt.subplots(figsize=(9, 7.5), dpi=300)
    ax.set_aspect("equal")
    ax.axis("off")

    # Place nodes in a regular circle
    angles = np.linspace(np.pi / 2, -3 * np.pi / 2, K, endpoint=False)
    radius = 2.5
    node_coords = {i: (radius * np.cos(angle), radius * np.sin(angle)) for i, angle in enumerate(angles)}

    # Draw directed edges
    for i in range(K):
        for j in range(K):
            if i == j:
                continue
            prob = P[i, j]
            if prob >= threshold:
                x1, y1 = node_coords[i]
                x2, y2 = node_coords[j]

                dx = x2 - x1
                dy = y2 - y1
                dist = np.sqrt(dx**2 + dy**2)
                if dist == 0:
                    continue

                # Shorten arrows to avoid overlapping node circles
                node_r = 0.45
                sx = x1 + (dx / dist) * node_r
                sy = y1 + (dy / dist) * node_r
                ex = x2 - (dx / dist) * node_r
                ey = y2 - (dy / dist) * node_r

                # Add curvature for bidirectional distinction
                rad = 0.18 if P[j, i] >= threshold else 0.08
                arrow = patches.FancyArrowPatch(
                    (sx, sy),
                    (ex, ey),
                    connectionstyle=f"arc3,rad={rad}",
                    color="#333333",
                    alpha=min(0.95, 0.4 + 0.6 * (prob / np.max(P))),
                    arrowstyle="-|>",
                    mutation_scale=14 + 10 * prob,
                    linewidth=1.2 + 4.0 * prob,
                    zorder=2,
                )
                ax.add_patch(arrow)

                # Label transition probability along arc
                mx = 0.5 * (sx + ex) - 0.25 * rad * (ey - sy)
                my = 0.5 * (sy + ey) + 0.25 * rad * (ex - sx)
                ax.text(
                    mx,
                    my,
                    f"{prob:.2f}",
                    fontsize=8.5,
                    fontweight="bold",
                    color="#111111",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85),
                    zorder=3,
                    ha="center",
                    va="center",
                )

    # Draw nodes
    for i in range(K):
        x, y = node_coords[i]
        node_color = FACIES_PALETTE.get(i, "#6C757D")

        # Determine text contrast
        brightness = sum(int(node_color.lstrip("#")[c : c + 2], 16) for c in (0, 2, 4)) / 3.0
        text_color = "white" if brightness < 140 else "black"

        circle = plt.Circle((x, y), 0.45, color=node_color, ec="#222222", lw=2.0, zorder=4)
        ax.add_patch(circle)

        # Multi-line label inside/around node
        label_text = labels[i].replace(" / ", "\n").replace(" Sandstone", "\nSand").replace(" Mudstone", "\nMud")
        ax.text(
            x,
            y,
            f"[{i}]\n{label_text}",
            ha="center",
            va="center",
            color=text_color,
            fontsize=8,
            fontweight="bold",
            zorder=5,
        )

    ax.set_xlim(-3.6, 3.6)
    ax.set_ylim(-3.6, 3.6)

    chain_type = "Embedded" if model.embedded else "Regular"
    ax.set_title(
        f"Facies Succession Network ({chain_type} Markov Chain, $P_{{ij}} \\geq {threshold:.2f}$)\n"
        f"Arrows indicate upward stratigraphic transitions (Lower Bed $\\to$ Upper Bed)",
        fontsize=12,
        fontweight="bold",
        pad=10,
    )

    plt.tight_layout()

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return fig


def plot_stationary_vs_empirical(
    model: StratigraphicMarkovChain,
    df: pd.DataFrame,
    facies_col: str = "facies_code",
    save_path: Optional[Union[str, Path]] = None,
) -> plt.Figure:
    """
    Bar chart comparing theoretical stationary distribution pi against empirical facies frequencies.

    Args:
        model: Fitted StratigraphicMarkovChain instance with computed stationary_dist_.
        df: Validated DataFrame containing digitized litholog observations.
        facies_col: Name of integer facies code column.
        save_path: Optional file path to save figure at 300 DPI.

    Returns:
        fig: Matplotlib Figure object.
    """
    if model.stationary_dist_ is None:
        raise ValueError("Model must have stationary_dist_ computed before plotting.")

    pi = model.stationary_dist_
    K = model.num_classes
    labels = [model.facies_map.get(i, f"State {i}") for i in range(K)]

    # Compute empirical frequency across all observations
    counts = df[facies_col].value_counts().to_dict()
    total_obs = len(df)
    empirical = np.array([counts.get(i, 0) / total_obs for i in range(K)], dtype=np.float64)

    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=300)

    x = np.arange(K)
    width = 0.35

    rects1 = ax.bar(
        x - width / 2,
        pi,
        width,
        label=r"Stationary Distribution ($\pi$)",
        color="#2B5B84",
        edgecolor="#1A3750",
        linewidth=1.2,
    )
    rects2 = ax.bar(
        x + width / 2,
        empirical,
        width,
        label=r"Empirical Observed Proportion ($p_{\mathrm{emp}}$)",
        color="#D97D27",
        edgecolor="#965214",
        linewidth=1.2,
    )

    # Annotate bar values
    def autolabel(rects):
        for rect in rects:
            h = rect.get_height()
            ax.annotate(
                f"{h:.3f}\n({h*100:.1f}%)",
                xy=(rect.get_x() + rect.get_width() / 2, h),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8.5,
                fontweight="bold",
            )

    autolabel(rects1)
    autolabel(rects2)

    ax.set_ylabel("Facies Proportion / Probability", fontsize=11, fontweight="bold")
    ax.set_title(
        r"Facies Proportion Diagnostic: Stationary Distribution ($\pi$) vs. Empirical Observed ($p_{\mathrm{emp}}$)",
        fontsize=12,
        fontweight="bold",
        pad=12,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=10, fontweight="medium")
    ax.set_ylim(0, max(np.max(pi), np.max(empirical)) * 1.25)
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    ax.legend(frameon=True, fontsize=10, loc="upper right")

    plt.tight_layout()

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return fig
