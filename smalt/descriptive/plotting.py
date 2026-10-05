"""
SMALT Stratigraphic Visualization Engine (Sprint C).

Generates clear, publishable figures suitable for undergraduate research presentation
and faculty review by Professor Hiranya Sahoo.

All figures strictly use color-blind accessible, standard sedimentological palettes,
clearly distinguish source-derived vs AI-reconstructed provenance, and annotate
geological caveats.
"""

from pathlib import Path
from typing import Dict, List, Optional, Union
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from smalt.descriptive.analyzer import (
    LithologInspector,
    CANONICAL_FACIES_SCHEMA,
    PROVENANCE_METADATA,
)


def plot_vertical_successions(
    inspector: LithologInspector,
    output_path: Union[str, Path],
    litholog_ids: Optional[List[str]] = None,
) -> None:
    """
    Plots side-by-side vertical stratigraphic facies columns for all lithologs.
    """
    if litholog_ids is None:
        litholog_ids = [f"litholog{i}" for i in range(1, 13)]

    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(litholog_ids),
        figsize=(22, 12),
        sharey=True,
    )

    if len(litholog_ids) == 1:
        axes = [axes]

    for ax, lid in zip(axes, litholog_ids):
        df = inspector.load_raw_litholog(lid)
        prov = PROVENANCE_METADATA.get(lid, {})

        # Plot each interval as a colored rectangle
        for _, row in df.iterrows():
            top = row["Top"]
            bot = row["Bottom"]
            facies = row["Facies"]
            color = CANONICAL_FACIES_SCHEMA[facies]["color_hex"]
            ax.bar(
                x=0.5,
                height=bot - top,
                bottom=top,
                width=0.8,
                color=color,
                edgecolor="black",
                linewidth=0.5,
            )

        ax.set_xlim(0, 1)
        ax.set_xticks([])
        # Depth downwards (cliff top 0 down to bottom)
        ax.invert_yaxis()

        # Category tag
        pcat = prov.get("provenance_category", "")
        if "source_derived" in pcat:
            tag = "Source"
            box_col = "#D4EFDF"
        elif "ai_reconstructed" in pcat:
            tag = "AI"
            box_col = "#FCF3CF"
        else:
            tag = "Core"
            box_col = "#E8DAEF"

        group_tag = prov.get("group", "").upper()
        if "UPSTREAM" in group_tag:
            gt = "UP"
        elif "DOWNSTREAM" in group_tag:
            gt = "DN"
        else:
            gt = "NA"

        ax.set_title(
            f"{lid.upper()}\n[{tag}|{gt}]\n{df['Bottom'].max():.0f}m",
            fontsize=10,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor=box_col, edgecolor="gray", alpha=0.9),
        )

    axes[0].set_ylabel("Depth Below Cliff Top (meters)", fontsize=12, fontweight="bold")
    axes[0].set_ylim(115, 0)

    # Legend
    legend_patches = [
        mpatches.Patch(
            facecolor=meta["color_hex"],
            edgecolor="black",
            label=f"{meta['canonical_name']} (Code {meta['code']})",
        )
        for meta in CANONICAL_FACIES_SCHEMA.values()
    ]
    fig.legend(
        handles=legend_patches,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.98),
        ncol=5,
        fontsize=11,
        frameon=True,
        shadow=True,
    )

    fig.suptitle(
        "SMALT Stratigraphic Transect: Vertical Facies Successions (Lithologs L1 - L12)\n"
        "[Note: Equal numerical depths do NOT represent a common chronostratigraphic datum; modern cliff tops are not aligned]",
        fontsize=13,
        fontweight="bold",
        y=1.03,
    )

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_facies_proportions(
    inspector: LithologInspector,
    output_path: Union[str, Path],
    litholog_ids: Optional[List[str]] = None,
) -> None:
    """
    Plots stacked horizontal bar chart of continuous facies proportions for each litholog.
    """
    if litholog_ids is None:
        litholog_ids = [f"litholog{i}" for i in range(1, 13)]

    data = []
    for lid in litholog_ids:
        insp = inspector.inspect_litholog(lid)
        row = {"litholog_id": lid, "group": insp["group"], "provenance": insp["provenance_category"]}
        row.update(insp["facies_proportions"])
        row["ntg_pure"] = insp["ntg_pure"]
        data.append(row)

    df_plot = pd.DataFrame(data)

    fig, ax = plt.subplots(figsize=(14, 8))

    facies_keys = list(CANONICAL_FACIES_SCHEMA.keys())
    bottom_vals = np.zeros(len(df_plot))

    for fkey in facies_keys:
        meta = CANONICAL_FACIES_SCHEMA[fkey]
        vals = df_plot[fkey].to_numpy() * 100.0
        ax.barh(
            y=df_plot["litholog_id"],
            width=vals,
            left=bottom_vals,
            color=meta["color_hex"],
            edgecolor="black",
            linewidth=0.5,
            label=f"{meta['canonical_name']}",
        )
        bottom_vals += vals

    # Annotate Net-to-Gross
    for i, (_, row) in enumerate(df_plot.iterrows()):
        ax.text(
            102,
            i,
            f"N/G: {row['ntg_pure']*100:.1f}%\n[{row['group'][:2].upper()}]",
            va="center",
            fontsize=9,
            fontweight="bold",
        )

    ax.set_xlim(0, 115)
    ax.set_xlabel("Continuous Stratigraphic Facies Proportion (%)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Litholog Identifier", fontsize=11, fontweight="bold")
    ax.set_title(
        "Facies Proportions and Net-to-Gross Across All Lithologs (Continuous Interval Thicknesses)",
        fontsize=13,
        fontweight="bold",
        pad=15,
    )
    ax.legend(loc="lower right", bbox_to_anchor=(0.95, 0.05), fontsize=10)
    ax.grid(axis="x", linestyle="--", alpha=0.5)

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_upstream_vs_downstream(
    inspector: LithologInspector,
    output_path: Union[str, Path],
) -> None:
    """
    Compares facies proportions and Net-to-Gross between Upstream and Downstream groups.
    """
    upstream_ids = [f"litholog{i}" for i in [2, 3, 4, 5, 6, 7, 8, 10]]
    downstream_all_ids = ["litholog9", "litholog11", "litholog12"]
    downstream_outcrop_ids = ["litholog9", "litholog11"]

    def aggregate_props(lids: List[str]) -> Tuple[Dict[str, float], float, float]:
        tot_thick = 0.0
        fac_thick = {k: 0.0 for k in CANONICAL_FACIES_SCHEMA.keys()}
        for lid in lids:
            insp = inspector.inspect_litholog(lid)
            tot_thick += insp["sum_thickness_m"]
            for k in CANONICAL_FACIES_SCHEMA.keys():
                fac_thick[k] += insp["facies_thickness_m"][k]

        props = {k: fac_thick[k] / tot_thick for k in CANONICAL_FACIES_SCHEMA.keys()}
        ntg_pure = fac_thick["sand"] / tot_thick
        ntg_coarse = (fac_thick["sand"] + fac_thick["silt"]) / tot_thick
        return props, ntg_pure, ntg_coarse

    up_props, up_ntg_p, up_ntg_c = aggregate_props(upstream_ids)
    dn_all_props, dn_all_ntg_p, dn_all_ntg_c = aggregate_props(downstream_all_ids)
    dn_out_props, dn_out_ntg_p, dn_out_ntg_c = aggregate_props(downstream_outcrop_ids)

    groups = [
        "Upstream (L2-L8, L10)\n[8 logs | 664.0m]",
        "Downstream Outcrop (L9, L11)\n[2 logs | 156.0m]",
        "Downstream Composite (L9, L11, L12)\n[3 logs | 267.0m]",
    ]

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 7), gridspec_kw={"width_ratios": [2.5, 1]})

    facies_keys = list(CANONICAL_FACIES_SCHEMA.keys())
    x = np.arange(len(groups))
    width = 0.15

    for idx, fkey in enumerate(facies_keys):
        meta = CANONICAL_FACIES_SCHEMA[fkey]
        vals = [
            up_props[fkey] * 100.0,
            dn_out_props[fkey] * 100.0,
            dn_all_props[fkey] * 100.0,
        ]
        rects = ax1.bar(
            x + (idx - 2) * width,
            vals,
            width,
            label=meta["canonical_name"],
            color=meta["color_hex"],
            edgecolor="black",
            linewidth=0.5,
        )
        for rect in rects:
            height = rect.get_height()
            if height >= 2.0:
                ax1.annotate(
                    f"{height:.1f}%",
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )

    ax1.set_xticks(x)
    ax1.set_xticklabels(groups, fontsize=10, fontweight="bold")
    ax1.set_ylabel("Facies Proportion (%)", fontsize=11, fontweight="bold")
    ax1.set_title("Facies Proportion Comparison: Upstream vs Downstream", fontsize=12, fontweight="bold")
    ax1.set_ylim(0, 60)
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(axis="y", linestyle="--", alpha=0.5)

    # Net-to-gross comparison
    ntg_p_vals = [up_ntg_p * 100, dn_out_ntg_p * 100, dn_all_ntg_p * 100]
    ntg_c_vals = [up_ntg_c * 100, dn_out_ntg_c * 100, dn_all_ntg_c * 100]

    ax2.bar(x - 0.15, ntg_p_vals, width=0.3, label="N/G Pure Sand", color="#F4D03F", edgecolor="black")
    ax2.bar(x + 0.15, ntg_c_vals, width=0.3, label="N/G Sand + Silt", color="#E59866", edgecolor="black")

    for i in range(len(groups)):
        ax2.text(i - 0.15, ntg_p_vals[i] + 1.5, f"{ntg_p_vals[i]:.1f}%", ha="center", fontsize=9, fontweight="bold")
        ax2.text(i + 0.15, ntg_c_vals[i] + 1.5, f"{ntg_c_vals[i]:.1f}%", ha="center", fontsize=9, fontweight="bold")

    ax2.set_xticks(x)
    ax2.set_xticklabels(["Upstream", "Downstream\n(Outcrop)", "Downstream\n(Composite)"], fontsize=10, fontweight="bold")
    ax2.set_ylabel("Net-to-Gross Ratio (%)", fontsize=11, fontweight="bold")
    ax2.set_title("Net-to-Gross Comparison", fontsize=12, fontweight="bold")
    ax2.set_ylim(0, 75)
    ax2.legend(loc="lower right", fontsize=9)
    ax2.grid(axis="y", linestyle="--", alpha=0.5)

    fig.suptitle(
        "Sedimentological Comparison: Proximal (Upstream) vs Distal (Downstream) Stratigraphy",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_sandstone_thickness_distributions(
    inspector: LithologInspector,
    output_path: Union[str, Path],
) -> None:
    """
    Plots sandstone interval bed thickness distributions across lithologs.
    """
    records = []
    for i in range(1, 13):
        lid = f"litholog{i}"
        df = inspector.load_raw_litholog(lid)
        sand_beds = df[df["Facies"] == "sand"]
        for _, r in sand_beds.iterrows():
            records.append({
                "litholog_id": lid,
                "thickness_m": float(r["Thickness"]),
                "group": PROVENANCE_METADATA[lid]["group"],
            })

    df_sand = pd.DataFrame(records)

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))

    # Boxplot by well
    wells = [f"litholog{i}" for i in range(1, 13)]
    box_data = [df_sand[df_sand["litholog_id"] == w]["thickness_m"].values for w in wells]

    bp = ax1.boxplot(box_data, tick_labels=[f"L{i}" for i in range(1, 13)], patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("#F4D03F")
        patch.set_alpha(0.8)

    ax1.set_xlabel("Litholog", fontsize=11, fontweight="bold")
    ax1.set_ylabel("Sandstone Bed Thickness (m)", fontsize=11, fontweight="bold")
    ax1.set_title("Sandstone Bed Thickness Distributions by Litholog", fontsize=12, fontweight="bold")
    ax1.grid(axis="y", linestyle="--", alpha=0.5)

    # Histogram / density overall
    ax2.hist(df_sand["thickness_m"], bins=15, color="#F39C12", edgecolor="black", alpha=0.7)
    mean_val = df_sand["thickness_m"].mean()
    median_val = df_sand["thickness_m"].median()
    ax2.axvline(mean_val, color="red", linestyle="--", linewidth=1.5, label=f"Mean: {mean_val:.2f} m")
    ax2.axvline(median_val, color="blue", linestyle="-.", linewidth=1.5, label=f"Median: {median_val:.2f} m")

    ax2.set_xlabel("Sandstone Bed Thickness (m)", fontsize=11, fontweight="bold")
    ax2.set_ylabel("Bed Count", fontsize=11, fontweight="bold")
    ax2.set_title(f"All Sandstone Intervals (N = {len(df_sand)} beds)", fontsize=12, fontweight="bold")
    ax2.legend(fontsize=10)
    ax2.grid(axis="y", linestyle="--", alpha=0.5)

    fig.suptitle(
        "Fluvial Sandstone Lithosome Thickness Architecture (Continuous Field Intervals)",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_transition_matrix_heatmaps(
    P_reg: np.ndarray,
    P_emb: np.ndarray,
    pi_reg: np.ndarray,
    pi_emb: np.ndarray,
    output_path: Union[str, Path],
) -> None:
    """
    Plots annotated heatmaps of Regular and Embedded transition probability matrices.
    """
    labels = [meta["canonical_name"].split(" ")[0] for meta in CANONICAL_FACIES_SCHEMA.values()]

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 7))

    # Regular matrix
    im1 = ax1.imshow(P_reg, cmap="YlGnBu", vmin=0.0, vmax=1.0)
    ax1.set_xticks(range(len(labels)))
    ax1.set_yticks(range(len(labels)))
    ax1.set_xticklabels(labels, rotation=45, ha="right", fontsize=10, fontweight="bold")
    ax1.set_yticklabels(labels, fontsize=10, fontweight="bold")
    ax1.set_xlabel("To Facies (Upward in Succession)", fontsize=11, fontweight="bold")
    ax1.set_ylabel("From Facies (Lower Bed)", fontsize=11, fontweight="bold")
    ax1.set_title("Regular Transition Matrix (1m Discretized, P_ii > 0)\n[Self-Transitions Retained]", fontsize=12, fontweight="bold")

    for i in range(len(labels)):
        for j in range(len(labels)):
            val = P_reg[i, j]
            text_col = "white" if val > 0.5 else "black"
            ax1.text(j, i, f"{val:.3f}", ha="center", va="center", color=text_col, fontsize=9, fontweight="bold")

    # Add stationary bar underneath or title
    stat_str_reg = " | ".join([f"{labels[k]}: {pi_reg[k]*100:.1f}%" for k in range(len(labels))])
    ax1.text(0.5, -0.25, f"Stationary Vector:\n{stat_str_reg}", ha="center", va="center", transform=ax1.transAxes, fontsize=9)

    # Embedded matrix
    im2 = ax2.imshow(P_emb, cmap="YlOrRd", vmin=0.0, vmax=1.0)
    ax2.set_xticks(range(len(labels)))
    ax2.set_yticks(range(len(labels)))
    ax2.set_xticklabels(labels, rotation=45, ha="right", fontsize=10, fontweight="bold")
    ax2.set_yticklabels(labels, fontsize=10, fontweight="bold")
    ax2.set_xlabel("To Facies (Upward in Succession)", fontsize=11, fontweight="bold")
    ax2.set_ylabel("From Facies (Lower Bed)", fontsize=11, fontweight="bold")
    ax2.set_title("Embedded Transition Matrix (Continuous Beds, P_ii = 0)\n[Boundary-Crossing Changes Only]", fontsize=12, fontweight="bold")

    for i in range(len(labels)):
        for j in range(len(labels)):
            val = P_emb[i, j]
            text_col = "white" if val > 0.5 else "black"
            ax2.text(j, i, f"{val:.3f}", ha="center", va="center", color=text_col, fontsize=9, fontweight="bold")

    stat_str_emb = " | ".join([f"{labels[k]}: {pi_emb[k]*100:.1f}%" for k in range(len(labels))])
    ax2.text(0.5, -0.25, f"Stationary Vector:\n{stat_str_emb}", ha="center", va="center", transform=ax2.transAxes, fontsize=9)

    fig.suptitle(
        "1D Vertical Stratigraphic Markov Transition Models (Tallied Upward: Deep -> Shallow)",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
