"""
SMALT Litholog Group Heterogeneity Analysis Module (Sprint K).

Computes rigorous statistical comparisons between established litholog groups:
- Upstream: L2-L8 and L10 (8 outcrop logs)
- Downstream outcrops: L9 and L11 (2 outcrop logs)
- Downstream composite: L9, L11, and L12 (retaining L12's distinct core provenance)
- L12 core alone: EM-137C (1 subsurface core log, 63 digitized intervals, 0-111 m)
- Pooled all: L1-L12 (12 logs)
- Pooled spatial: L2-L12 (11 logs with coordinates)

Calculates:
1. Facies thickness and interval proportions across canonical 6 states.
2. Continuous merged bed-thickness distributions and descriptive metrics.
3. Regular (1m discretized) and embedded (boundary-crossing) vertical transition matrices.
4. Quantitative discrepancy metrics (TV distance, Frobenius norm) and sample-size caveats.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from smalt.descriptive.analyzer import (
    LithologInspector,
    CANONICAL_FACIES_SCHEMA,
    PROVENANCE_METADATA,
)
from smalt.geostat.markov import StratigraphicMarkovChain, CANONICAL_TECHNICAL_MAP

# Established litholog groups
GROUPS_CONFIG: Dict[str, Dict[str, Any]] = {
    "upstream": {
        "name": "Upstream Group",
        "wells": ["litholog2", "litholog3", "litholog4", "litholog5", "litholog6", "litholog7", "litholog8", "litholog10"],
        "description": "8 outcrop sections from upstream depositional reach (L2-L8, L10).",
        "observation_type": "outcrop",
    },
    "downstream_outcrop": {
        "name": "Downstream Outcrop Group",
        "wells": ["litholog9", "litholog11"],
        "description": "2 measured outcrop sections from downstream reach (L9, L11).",
        "observation_type": "outcrop",
    },
    "downstream_composite": {
        "name": "Downstream Composite Group",
        "wells": ["litholog9", "litholog11", "litholog12"],
        "description": "Downstream outcrop sections (L9, L11) plus subsurface core EM-137C (L12).",
        "observation_type": "mixed_outcrop_core",
    },
    "l12_core": {
        "name": "L12 Core Alone",
        "wells": ["litholog12"],
        "description": "Single subsurface core log EM-137C (0-111 m, 63 digitized intervals).",
        "observation_type": "core",
    },
    "pooled_all": {
        "name": "Pooled All 12 Logs",
        "wells": [f"litholog{i}" for i in range(1, 13)],
        "description": "All 12 digitized lithologs in project collection (L1-L12).",
        "observation_type": "all_collection",
    },
    "pooled_spatial": {
        "name": "Pooled L2-L12",
        "wells": [f"litholog{i}" for i in range(2, 13)],
        "description": "11 lithologs with coordinates documented in spreadsheet (L2-L12).",
        "observation_type": "spatial_eligible",
    },
}

FACIES_KEYS: List[str] = ["sand", "p_sand", "ripples", "carbon_mud", "coal", "mud"]
FACIES_CODES: List[int] = [0, 1, 2, 3, 4, 5]


class LithologGroupAnalyzer:
    """
    Analyzes stratigraphic heterogeneity across defined litholog groups.
    """

    def __init__(
        self,
        raw_dir: str = "data/raw_lithologs",
        groups_config: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> None:
        self.inspector = LithologInspector(raw_dir=raw_dir)
        self.groups_config = groups_config or GROUPS_CONFIG.copy()
        self._cached_beds: Dict[str, pd.DataFrame] = {}
        self._cached_disc: Dict[str, pd.DataFrame] = {}

    def get_group_beds(self, group_key: str) -> pd.DataFrame:
        """Returns consolidated merged continuous beds for all wells in a group."""
        if group_key in self._cached_beds:
            return self._cached_beds[group_key].copy()

        wells = self.groups_config[group_key]["wells"]
        records = []
        for well in wells:
            well_beds = self.inspector.build_embedded_bed_sequence(well)
            well_beds["litholog_id"] = well
            well_beds["group"] = group_key
            records.append(well_beds)

        if not records:
            df = pd.DataFrame()
        else:
            df = pd.concat(records, ignore_index=True)

        self._cached_beds[group_key] = df
        return df.copy()

    def get_group_discretized(self, group_key: str) -> pd.DataFrame:
        """Returns consolidated 1m discretized profiles for all wells in a group."""
        if group_key in self._cached_disc:
            return self._cached_disc[group_key].copy()

        wells = self.groups_config[group_key]["wells"]
        records = []
        for well in wells:
            well_disc = self.inspector.discretize_litholog_1m(well)
            well_disc["group"] = group_key
            records.append(well_disc)

        if not records:
            df = pd.DataFrame()
        else:
            df = pd.concat(records, ignore_index=True)

        self._cached_disc[group_key] = df
        return df.copy()

    def compute_facies_proportions(self, group_key: str) -> Dict[str, Any]:
        """
        Computes sample counts, thickness sums, and proportions per facies for a group.
        """
        beds = self.get_group_beds(group_key)
        total_thickness = float(beds["thickness_m"].sum()) if not beds.empty else 0.0
        total_beds = len(beds)
        wells = self.groups_config[group_key]["wells"]

        prop_dict: Dict[str, float] = {}
        thick_dict: Dict[str, float] = {}
        bed_count_dict: Dict[str, int] = {}
        bed_count_prop: Dict[str, float] = {}

        for f in FACIES_KEYS:
            sub = beds[beds["facies"] == f]
            t = float(sub["thickness_m"].sum())
            c = len(sub)
            thick_dict[f] = round(t, 2)
            bed_count_dict[f] = c
            prop_dict[f] = round(t / total_thickness, 4) if total_thickness > 0 else 0.0
            bed_count_prop[f] = round(c / total_beds, 4) if total_beds > 0 else 0.0

        # Net-to-gross variations
        sand_t = thick_dict["sand"]
        p_sand_t = thick_dict["p_sand"]
        ripples_t = thick_dict["ripples"]

        ntg_channel = round(sand_t / total_thickness, 4) if total_thickness > 0 else 0.0
        ntg_pure = round((sand_t + p_sand_t) / total_thickness, 4) if total_thickness > 0 else 0.0
        ntg_coarse = round((sand_t + p_sand_t + ripples_t) / total_thickness, 4) if total_thickness > 0 else 0.0

        return {
            "group_key": group_key,
            "group_name": self.groups_config[group_key]["name"],
            "num_wells": len(wells),
            "total_beds": total_beds,
            "total_thickness_m": round(total_thickness, 2),
            "thickness_m": thick_dict,
            "thickness_proportions": prop_dict,
            "bed_counts": bed_count_dict,
            "bed_count_proportions": bed_count_prop,
            "ntg_channel": ntg_channel,
            "ntg_pure": ntg_pure,
            "ntg_coarse": ntg_coarse,
        }

    def compute_bed_thickness_distributions(self, group_key: str) -> Dict[str, Dict[str, Any]]:
        """
        Computes summary statistics of bed thickness for each facies in a group.
        """
        beds = self.get_group_beds(group_key)
        stats_by_facies: Dict[str, Dict[str, Any]] = {}

        for f in FACIES_KEYS:
            sub = beds[beds["facies"] == f]
            thick_vals = sub["thickness_m"].to_numpy()
            n = len(thick_vals)

            if n == 0:
                stats_by_facies[f] = {
                    "count": 0,
                    "mean_m": 0.0,
                    "std_m": 0.0,
                    "median_m": 0.0,
                    "min_m": 0.0,
                    "max_m": 0.0,
                    "p25_m": 0.0,
                    "p75_m": 0.0,
                    "iqr_m": 0.0,
                    "p90_m": 0.0,
                    "skewness": 0.0,
                }
            elif n == 1:
                val = float(thick_vals[0])
                stats_by_facies[f] = {
                    "count": 1,
                    "mean_m": val,
                    "std_m": 0.0,
                    "median_m": val,
                    "min_m": val,
                    "max_m": val,
                    "p25_m": val,
                    "p75_m": val,
                    "iqr_m": 0.0,
                    "p90_m": val,
                    "skewness": 0.0,
                }
            else:
                p25 = float(np.percentile(thick_vals, 25))
                p75 = float(np.percentile(thick_vals, 75))
                stats_by_facies[f] = {
                    "count": n,
                    "mean_m": round(float(np.mean(thick_vals)), 3),
                    "std_m": round(float(np.std(thick_vals, ddof=1)), 3),
                    "median_m": round(float(np.median(thick_vals)), 3),
                    "min_m": round(float(np.min(thick_vals)), 3),
                    "max_m": round(float(np.max(thick_vals)), 3),
                    "p25_m": round(p25, 3),
                    "p75_m": round(p75, 3),
                    "iqr_m": round(p75 - p25, 3),
                    "p90_m": round(float(np.percentile(thick_vals, 90)), 3),
                    "skewness": round(float(stats.skew(thick_vals)), 3) if n >= 3 else 0.0,
                }

        return stats_by_facies

    def compute_transition_matrices(
        self,
        group_key: str,
        smoothing_alpha: float = 0.0,
    ) -> Dict[str, Any]:
        """
        Computes both regular (1m discretized) and embedded (continuous boundary)
        Markov transition matrices for a group in upward stratigraphic order.
        """
        wells = self.groups_config[group_key]["wells"]

        # 1. Regular 1m discretized transition matrix
        disc_df = self.get_group_discretized(group_key)
        reg_mc = StratigraphicMarkovChain(
            num_classes=6,
            embedded=False,
            smoothing_alpha=smoothing_alpha,
        )
        reg_mc.fit(
            disc_df,
            depth_col="depth_m",
            facies_col="facies_code",
            well_col="litholog_id",
        )

        # 2. Embedded boundary-crossing transition matrix
        # Fit on bed sequences where self-transitions are strictly 0
        bed_df = self.get_group_beds(group_key)
        # Note: StratigraphicMarkovChain expects depth_col sorted descending.
        # bed_df has bottom_m and top_m. We can use bottom_m as depth proxy.
        emb_mc = StratigraphicMarkovChain(
            num_classes=6,
            embedded=True,
            smoothing_alpha=smoothing_alpha,
        )
        emb_mc.fit(
            bed_df,
            depth_col="bottom_m",
            facies_col="facies_code",
            well_col="litholog_id",
        )

        return {
            "group_key": group_key,
            "regular_count_matrix": reg_mc.count_matrix_.tolist(),
            "regular_transition_matrix": reg_mc.transition_matrix_.tolist(),
            "regular_stationary_dist": reg_mc.stationary_dist_.tolist() if reg_mc.stationary_dist_ is not None else [],
            "regular_total_transitions": reg_mc.n_transitions_,
            "embedded_count_matrix": emb_mc.count_matrix_.tolist(),
            "embedded_transition_matrix": emb_mc.transition_matrix_.tolist(),
            "embedded_stationary_dist": emb_mc.stationary_dist_.tolist() if emb_mc.stationary_dist_ is not None else [],
            "embedded_total_transitions": emb_mc.n_transitions_,
        }

    def compare_groups(
        self,
        group_a: str = "upstream",
        group_b: str = "downstream_outcrop",
    ) -> Dict[str, Any]:
        """
        Computes formal discrepancy metrics and statistical tests between two groups.
        """
        prop_a = self.compute_facies_proportions(group_a)
        prop_b = self.compute_facies_proportions(group_b)

        # 1. Total Variation (TV) distance between thickness proportions
        p_vec_a = np.array([prop_a["thickness_proportions"][f] for f in FACIES_KEYS])
        p_vec_b = np.array([prop_b["thickness_proportions"][f] for f in FACIES_KEYS])
        tv_dist = 0.5 * np.sum(np.abs(p_vec_a - p_vec_b))

        # 2. Chi-squared test on bed counts (contingency table)
        counts_a = np.array([prop_a["bed_counts"][f] for f in FACIES_KEYS])
        counts_b = np.array([prop_b["bed_counts"][f] for f in FACIES_KEYS])
        contingency = np.vstack([counts_a, counts_b])
        # Filter columns where both are 0
        valid_cols = (counts_a + counts_b) > 0
        if np.sum(valid_cols) >= 2:
            chi2_stat, p_val_chi2, dof, _ = stats.chi2_contingency(contingency[:, valid_cols])
        else:
            chi2_stat, p_val_chi2, dof = 0.0, 1.0, 0

        # 3. Transition matrix discrepancy
        tm_a = self.compute_transition_matrices(group_a, smoothing_alpha=0.1)
        tm_b = self.compute_transition_matrices(group_b, smoothing_alpha=0.1)
        P_reg_a = np.array(tm_a["regular_transition_matrix"])
        P_reg_b = np.array(tm_b["regular_transition_matrix"])
        frob_reg = float(np.linalg.norm(P_reg_a - P_reg_b, ord="fro"))
        mse_reg = float(np.mean((P_reg_a - P_reg_b) ** 2))

        P_emb_a = np.array(tm_a["embedded_transition_matrix"])
        P_emb_b = np.array(tm_b["embedded_transition_matrix"])
        frob_emb = float(np.linalg.norm(P_emb_a - P_emb_b, ord="fro"))
        mse_emb = float(np.mean((P_emb_a - P_emb_b) ** 2))

        # 4. Facies-specific bed thickness tests (Mann-Whitney U)
        beds_a = self.get_group_beds(group_a)
        beds_b = self.get_group_beds(group_b)
        bed_tests: Dict[str, Dict[str, Any]] = {}

        for f in FACIES_KEYS:
            vals_a = beds_a[beds_a["facies"] == f]["thickness_m"].to_numpy()
            vals_b = beds_b[beds_b["facies"] == f]["thickness_m"].to_numpy()
            n_a = len(vals_a)
            n_b = len(vals_b)

            if n_a >= 3 and n_b >= 3:
                u_stat, u_pval = stats.mannwhitneyu(vals_a, vals_b, alternative="two-sided")
                bed_tests[f] = {
                    "n_a": n_a,
                    "n_b": n_b,
                    "mean_a_m": round(float(np.mean(vals_a)), 2),
                    "mean_b_m": round(float(np.mean(vals_b)), 2),
                    "u_stat": float(u_stat),
                    "p_value": round(float(u_pval), 4),
                    "statistically_significant": bool(u_pval < 0.05),
                    "sample_size_adequate": bool(n_a >= 5 and n_b >= 5),
                }
            else:
                bed_tests[f] = {
                    "n_a": n_a,
                    "n_b": n_b,
                    "mean_a_m": round(float(np.mean(vals_a)), 2) if n_a > 0 else 0.0,
                    "mean_b_m": round(float(np.mean(vals_b)), 2) if n_b > 0 else 0.0,
                    "u_stat": None,
                    "p_value": None,
                    "statistically_significant": False,
                    "sample_size_adequate": False,
                    "note": f"Insufficient sample size for formal testing (n_a={n_a}, n_b={n_b}).",
                }

        return {
            "comparison": f"{group_a} vs {group_b}",
            "tv_distance_proportions": round(float(tv_dist), 4),
            "chi2_contingency": {
                "statistic": round(float(chi2_stat), 3),
                "p_value": round(float(p_val_chi2), 4),
                "dof": int(dof),
            },
            "transition_frobenius_reg": round(frob_reg, 4),
            "transition_mse_reg": round(mse_reg, 6),
            "transition_frobenius_emb": round(frob_emb, 4),
            "transition_mse_emb": round(mse_emb, 6),
            "bed_thickness_tests": bed_tests,
        }

    def generate_all_summaries(self, output_dir: str = "sprints/audit_sprint_k") -> Dict[str, Any]:
        """
        Executes full group analysis, exports CSV tables, and plots figures.
        """
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        # 1. Proportions table across all groups
        prop_rows = []
        for g_key in self.groups_config.keys():
            p_res = self.compute_facies_proportions(g_key)
            row = {
                "group_key": g_key,
                "group_name": p_res["group_name"],
                "num_wells": p_res["num_wells"],
                "total_beds": p_res["total_beds"],
                "total_thickness_m": p_res["total_thickness_m"],
                "ntg_channel": p_res["ntg_channel"],
                "ntg_pure": p_res["ntg_pure"],
                "ntg_coarse": p_res["ntg_coarse"],
            }
            for f in FACIES_KEYS:
                row[f"{f}_thick_m"] = p_res["thickness_m"][f]
                row[f"{f}_prop"] = p_res["thickness_proportions"][f]
                row[f"{f}_bed_count"] = p_res["bed_counts"][f]
                row[f"{f}_count_prop"] = p_res["bed_count_proportions"][f]
            prop_rows.append(row)

        df_props = pd.DataFrame(prop_rows)
        csv_props_path = out_path / "litholog_group_facies_proportions.csv"
        df_props.to_csv(csv_props_path, index=False)

        # 2. Bed thickness summary table
        thick_rows = []
        for g_key in self.groups_config.keys():
            t_res = self.compute_bed_thickness_distributions(g_key)
            for f in FACIES_KEYS:
                st = t_res[f]
                thick_rows.append({
                    "group_key": g_key,
                    "group_name": self.groups_config[g_key]["name"],
                    "facies": f,
                    "facies_code": CANONICAL_FACIES_SCHEMA[f]["code"],
                    "count": st["count"],
                    "mean_m": st["mean_m"],
                    "std_m": st["std_m"],
                    "median_m": st["median_m"],
                    "min_m": st["min_m"],
                    "max_m": st["max_m"],
                    "p25_m": st["p25_m"],
                    "p75_m": st["p75_m"],
                    "iqr_m": st["iqr_m"],
                    "p90_m": st["p90_m"],
                    "skewness": st["skewness"],
                })

        df_thick = pd.DataFrame(thick_rows)
        csv_thick_path = out_path / "litholog_group_bed_thickness_summary.csv"
        df_thick.to_csv(csv_thick_path, index=False)

        # 3. Transition matrices table
        trans_rows = []
        for g_key in self.groups_config.keys():
            t_mat = self.compute_transition_matrices(g_key, smoothing_alpha=0.1)
            reg_mat = np.array(t_mat["regular_transition_matrix"])
            emb_mat = np.array(t_mat["embedded_transition_matrix"])
            reg_counts = np.array(t_mat["regular_count_matrix"])
            emb_counts = np.array(t_mat["embedded_count_matrix"])

            for i, f_from in enumerate(FACIES_KEYS):
                for j, f_to in enumerate(FACIES_KEYS):
                    trans_rows.append({
                        "group_key": g_key,
                        "group_name": self.groups_config[g_key]["name"],
                        "from_facies": f_from,
                        "to_facies": f_to,
                        "from_code": i,
                        "to_code": j,
                        "reg_count": float(reg_counts[i, j]),
                        "reg_prob": round(float(reg_mat[i, j]), 4),
                        "emb_count": float(emb_counts[i, j]),
                        "emb_prob": round(float(emb_mat[i, j]), 4),
                    })

        df_trans = pd.DataFrame(trans_rows)
        csv_trans_path = out_path / "litholog_group_transition_matrices.csv"
        df_trans.to_csv(csv_trans_path, index=False)

        # 4. Group comparisons summary
        comp_pairs = [
            ("upstream", "downstream_outcrop"),
            ("upstream", "downstream_composite"),
            ("upstream", "l12_core"),
            ("downstream_outcrop", "l12_core"),
            ("upstream", "pooled_all"),
        ]
        comp_records = []
        for g1, g2 in comp_pairs:
            cmp_res = self.compare_groups(g1, g2)
            comp_records.append({
                "comparison": cmp_res["comparison"],
                "tv_distance_proportions": cmp_res["tv_distance_proportions"],
                "chi2_stat": cmp_res["chi2_contingency"]["statistic"],
                "chi2_p_value": cmp_res["chi2_contingency"]["p_value"],
                "frob_reg": cmp_res["transition_frobenius_reg"],
                "mse_reg": cmp_res["transition_mse_reg"],
                "frob_emb": cmp_res["transition_frobenius_emb"],
                "mse_emb": cmp_res["transition_mse_emb"],
                "sand_test_p_val": cmp_res["bed_thickness_tests"]["sand"].get("p_value"),
                "sand_n_a": cmp_res["bed_thickness_tests"]["sand"]["n_a"],
                "sand_n_b": cmp_res["bed_thickness_tests"]["sand"]["n_b"],
                "mud_test_p_val": cmp_res["bed_thickness_tests"]["mud"].get("p_value"),
                "mud_n_a": cmp_res["bed_thickness_tests"]["mud"]["n_a"],
                "mud_n_b": cmp_res["bed_thickness_tests"]["mud"]["n_b"],
            })

        df_comp = pd.DataFrame(comp_records)
        csv_comp_path = out_path / "litholog_group_comparison_metrics.csv"
        df_comp.to_csv(csv_comp_path, index=False)

        # 5. Raw consolidated beds for transparency
        raw_beds_list = []
        for g_key in ["upstream", "downstream_outcrop", "l12_core"]:
            b_df = self.get_group_beds(g_key)
            raw_beds_list.append(b_df)
        df_all_beds = pd.concat(raw_beds_list, ignore_index=True)
        csv_raw_beds_path = out_path / "litholog_group_raw_beds.csv"
        df_all_beds.to_csv(csv_raw_beds_path, index=False)

        # 6. Generate Figures
        self._plot_proportions(df_props, out_path / "fig_group_facies_proportions.png")
        self._plot_bed_thickness_distributions(out_path / "fig_group_bed_thickness_distributions.png")
        self._plot_transition_heatmaps(out_path / "fig_group_transition_heatmaps.png")

        return {
            "csv_props": str(csv_props_path),
            "csv_thick": str(csv_thick_path),
            "csv_trans": str(csv_trans_path),
            "csv_comp": str(csv_comp_path),
            "csv_raw_beds": str(csv_raw_beds_path),
            "fig_props": str(out_path / "fig_group_facies_proportions.png"),
            "fig_thick": str(out_path / "fig_group_bed_thickness_distributions.png"),
            "fig_trans": str(out_path / "fig_group_transition_heatmaps.png"),
            "proportions_df": df_props,
            "thickness_df": df_thick,
            "comparisons_df": df_comp,
        }

    def _plot_proportions(self, df_props: pd.DataFrame, output_file: Path) -> None:
        """Plots stacked facies proportion comparison across groups."""
        display_groups = ["upstream", "downstream_outcrop", "l12_core", "downstream_composite", "pooled_all"]
        sub_df = df_props[df_props["group_key"].isin(display_groups)].copy()
        sub_df = sub_df.set_index("group_key").reindex(display_groups).reset_index()

        fig, ax = plt.subplots(figsize=(10, 6))
        group_names = [self.groups_config[g]["name"] for g in display_groups]
        x = np.arange(len(display_groups))
        bottom = np.zeros(len(display_groups))

        for f in FACIES_KEYS:
            color = CANONICAL_FACIES_SCHEMA[f]["color_hex"]
            label = f"{f} ({CANONICAL_FACIES_SCHEMA[f]['canonical_name']})"
            props = sub_df[f"{f}_prop"].to_numpy()
            ax.bar(x, props, bottom=bottom, label=label, color=color, edgecolor="black", linewidth=0.5, width=0.55)
            bottom += props

        ax.set_ylabel("Facies Proportion (Thickness-Weighted)", fontsize=12)
        ax.set_title("Facies Proportion Heterogeneity Across Litholog Groups\n(SMALT Sprint K Stage 1 Audit)", fontsize=13, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(group_names, rotation=15, ha="right", fontsize=10)
        ax.set_ylim(0, 1.0)
        ax.grid(axis="y", linestyle="--", alpha=0.5)
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, fontsize=9)
        plt.tight_layout()
        plt.savefig(output_file, dpi=200)
        plt.close(fig)

    def _plot_bed_thickness_distributions(self, output_file: Path) -> None:
        """Plots boxplots of bed thicknesses for channel sand and mud across groups."""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        groups = ["upstream", "downstream_outcrop", "l12_core", "pooled_all"]
        group_labels = ["Upstream\n(L2-L8,L10)", "Downstream Outcrop\n(L9,L11)", "L12 Core Alone\n(EM-137C)", "Pooled All\n(L1-L12)"]

        # 1. Channel Sandstone
        sand_data = []
        for g in groups:
            b = self.get_group_beds(g)
            vals = b[b["facies"] == "sand"]["thickness_m"].to_numpy()
            sand_data.append(vals if len(vals) > 0 else np.array([0.0]))

        ax0 = axes[0]
        bp0 = ax0.boxplot(sand_data, tick_labels=group_labels, patch_artist=True, showmeans=True)
        for box in bp0["boxes"]:
            box.set_facecolor(CANONICAL_FACIES_SCHEMA["sand"]["color_hex"])
            box.set_alpha(0.8)
        ax0.set_title("Channel Sandstone Bed Thickness (Facies 0)", fontsize=11, fontweight="bold")
        ax0.set_ylabel("Bed Thickness (m)", fontsize=11)
        ax0.grid(axis="y", linestyle="--", alpha=0.5)

        # 2. Overbank Mudstone
        mud_data = []
        for g in groups:
            b = self.get_group_beds(g)
            vals = b[b["facies"] == "mud"]["thickness_m"].to_numpy()
            mud_data.append(vals if len(vals) > 0 else np.array([0.0]))

        ax1 = axes[1]
        bp1 = ax1.boxplot(mud_data, tick_labels=group_labels, patch_artist=True, showmeans=True)
        for box in bp1["boxes"]:
            box.set_facecolor(CANONICAL_FACIES_SCHEMA["mud"]["color_hex"])
            box.set_alpha(0.8)
        ax1.set_title("Overbank Mudstone Bed Thickness (Facies 5)", fontsize=11, fontweight="bold")
        ax1.set_ylabel("Bed Thickness (m)", fontsize=11)
        ax1.grid(axis="y", linestyle="--", alpha=0.5)

        plt.suptitle("Stratigraphic Bed Thickness Comparison: Upstream vs Downstream vs L12 Core", fontsize=13, fontweight="bold", y=1.02)
        plt.tight_layout()
        plt.savefig(output_file, dpi=200)
        plt.close(fig)

    def _plot_transition_heatmaps(self, output_file: Path) -> None:
        """Plots embedded transition probability heatmaps for Upstream vs Downstream Outcrop vs Pooled."""
        groups = ["upstream", "downstream_outcrop", "pooled_all"]
        titles = ["Upstream (L2-L8, L10)", "Downstream Outcrop (L9, L11)", "Pooled All (L1-L12)"]

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        for ax, g_key, title in zip(axes, groups, titles):
            t_res = self.compute_transition_matrices(g_key, smoothing_alpha=0.1)
            P_emb = np.array(t_res["embedded_transition_matrix"])
            im = ax.imshow(P_emb, cmap="YlGnBu", vmin=0, vmax=0.8)
            ax.set_title(f"{title}\nEmbedded Transitions ($P_{{ii}}=0$)", fontsize=11, fontweight="bold")
            ax.set_xticks(range(6))
            ax.set_yticks(range(6))
            ax.set_xticklabels(FACIES_KEYS, rotation=45, ha="right", fontsize=9)
            ax.set_yticklabels(FACIES_KEYS, fontsize=9)
            ax.set_xlabel("Upper Facies $S_{t+1}$", fontsize=10)
            ax.set_ylabel("Lower Facies $S_t$", fontsize=10)

            # Annotate cells
            for i in range(6):
                for j in range(6):
                    val = P_emb[i, j]
                    text_color = "white" if val > 0.4 else "black"
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=text_color, fontsize=8)

        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.7, label="Transition Probability")
        plt.savefig(output_file, dpi=200)
        plt.close(fig)


if __name__ == "__main__":
    analyzer = LithologGroupAnalyzer()
    res = analyzer.generate_all_summaries()
    print("Litholog Group Heterogeneity Analysis Completed Successfully.")
    print("Generated CSVs:")
    print(" -", res["csv_props"])
    print(" -", res["csv_thick"])
    print(" -", res["csv_trans"])
    print(" -", res["csv_comp"])
    print("Generated Plots:")
    print(" -", res["fig_props"])
    print(" -", res["fig_thick"])
    print(" -", res["fig_trans"])
