#!/usr/bin/env python3
"""
Module 7: Visualization for ReCoN Multicellular Network Analysis

Generates:
1. Heatmaps of cell type coordination strength
2. Top genes heatmaps
3. Condition comparison plots (disease vs normal)
4. Differential cascade analysis (volcano-style)
5. CCC network heatmaps
6. Fibrosis/seed marker gene bar plots
7. Summary multi-panel figures
8. Sankey diagrams for top ligand-TF-gene cascades (if seed_categories configured)

All paths, conditions, and parameters are driven by ReconConfig.
"""

import json
import hashlib
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

# Sankey-related imports
from recon.data import load_receptor_genes
from recon.explore import Celltype, Multicell, set_lambda
from recon.plot import plot_intracell_sankey, plot_ligand_sankey, plot_intercell_sankey
from recon.plot.sankey_paths import (
    get_cell_communication_layer,
    get_celltype_gene_layer,
    get_celltype_grn_receptor_bipartite,
    get_top_ligands,
    get_top_receptors,
    get_top_tfs,
    build_partial_networks,
    plot_3layer_sankey,
    plot_4layer_sankey,
)
import plotly.graph_objects as go

from config import ReconConfig, get_config

warnings.filterwarnings("ignore")

# Figure parameters
FIGSIZE_LARGE = (11, 8)
FIGSIZE_MEDIUM = (8, 6)
FIGSIZE_SMALL = (8, 6)
DPI = 150

# Known acronyms that should keep their casing (excluding SSC which is special-cased)
_KNOWN_ACRONYMS = {"IPF", "COPD", "ALS", "IBD", "CKD", "NASH", "HCC", "AML"}


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def format_condition_name(name: str) -> str:
    """
    Format condition name for display.

    Preserves known acronyms (IPF, COPD, etc.).
    Special case: 'ssc'/'SSC' always returns 'SSc' (standard medical abbreviation).
    Otherwise capitalizes first letter of each word.
    """
    upper = name.upper()
    # Special case: SSC -> SSc (standard medical abbreviation)
    if upper == "SSC":
        return "SSc"
    # Other known acronyms: keep uppercase
    if upper in _KNOWN_ACRONYMS:
        return upper
    # Unknown: title case
    return name.title()


def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def setup_directories(config: ReconConfig) -> Tuple[Path, Path, Path]:
    """Create and return (figures_dir, data_dir, sankey_dir)."""
    figures_dir = config.get_figures_dir()
    data_dir = figures_dir / "data"
    sankey_dir = figures_dir / "sankey"
    figures_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    sankey_dir.mkdir(parents=True, exist_ok=True)
    print(f"[{_timestamp()}] Output directory: {figures_dir}")
    print(f"[{_timestamp()}] Data export directory: {data_dir}")
    print(f"[{_timestamp()}] Sankey output directory: {sankey_dir}")
    return figures_dir, data_dir, sankey_dir


def save_data_table(
    df: pd.DataFrame, output_path: Path, data_dir: Path, description: str = ""
) -> None:
    """Save a data table as a companion to a figure."""
    data_path = data_dir / f"{output_path.stem}_data.csv"
    df.to_csv(data_path)
    print(f"  Data: {data_path} ({len(df)} rows)")


# =============================================================================
# DATA LOADING
# =============================================================================

def load_recon_results(config: ReconConfig) -> Dict:
    """
    Load ReCoN combined-effect CSVs for every condition.

    Returns dict with keys like '<condition>_combined' and
    'comparison_<disease>_vs_<normal>'.
    """
    print(f"\n[{_timestamp()}] Loading ReCoN results...")
    recon_dir = config.get_recon_dir()
    results: Dict[str, pd.DataFrame] = {}

    for cond in config.conditions:
        cond_lower = cond.lower()
        path = recon_dir / f"{cond_lower}_combined_effects.csv"
        if path.exists():
            results[f"{cond_lower}_combined"] = pd.read_csv(path, index_col=0)
            print(f"  Loaded {cond_lower} combined: {results[f'{cond_lower}_combined'].shape}")
        else:
            print(f"  Warning: {path} not found")

    # Comparisons: each disease vs normal
    normal_lower = config.normal_condition.lower()
    for disease in config.disease_conditions:
        d_lower = disease.lower()
        comp_path = recon_dir / f"comparison_{d_lower}_vs_{normal_lower}.csv"
        if comp_path.exists():
            key = f"comparison_{d_lower}"
            results[key] = pd.read_csv(comp_path)
            print(f"  Loaded comparison {d_lower} vs {normal_lower}: {results[key].shape}")

    return results


def load_ccc_data(config: ReconConfig) -> Dict:
    """Load CCC data based on config.ccc_source for every condition."""
    print(f"\n[{_timestamp()}] Loading CCC data (source: {config.ccc_source})...")

    ccc_dir = config.get_ccc_dir()
    ccc_data: Dict[str, pd.DataFrame] = {}

    if config.ccc_source == "merged":
        base_path = ccc_dir / "merged"
    elif config.ccc_source == "cellphonedb":
        base_path = ccc_dir / "cellphonedb"
    else:
        base_path = ccc_dir / "cellchat"

    for cond in config.conditions:
        cond_lower = cond.lower()
        ccc_path = base_path / f"{cond_lower}_ccc.csv"
        if ccc_path.exists():
            df = pd.read_csv(ccc_path)

            if config.ccc_source == "cellchat":
                df = df.rename(columns={
                    "source": "celltype_source",
                    "target": "celltype_target",
                    "ligand": "source",
                    "receptor": "target",
                })
                df["lr_means"] = df["prob"]

            ccc_data[cond_lower] = df
            print(f"  Loaded {cond_lower.upper()} CCC: {len(df)} interactions")

    # Differential files
    normal_lower = config.normal_condition.lower()
    for disease in config.disease_conditions:
        d_lower = disease.lower()
        diff_path = ccc_dir / f"differential_{d_lower}_vs_{normal_lower}.csv"
        if diff_path.exists():
            key = f"diff_{d_lower}"
            ccc_data[key] = pd.read_csv(diff_path)
            print(f"  Loaded differential {d_lower}: {len(ccc_data[key])} interactions")

    return ccc_data


# =============================================================================
# STATIC FIGURE FUNCTIONS
# =============================================================================

def plot_coordination_heatmap(
    combined_effects: pd.DataFrame,
    condition_name: str,
    output_path: Path,
    data_dir: Path,
) -> None:
    """Create heatmap of cell type coordination strength."""
    print(f"\n[{_timestamp()}] Creating coordination heatmap for {condition_name}...")

    if combined_effects.empty:
        print("  No data available for heatmap")
        return

    corr_matrix = combined_effects.corr()

    fig, ax = plt.subplots(figsize=FIGSIZE_MEDIUM)
    sns.heatmap(
        corr_matrix,
        annot=True,
        fmt=".2f",
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0.5,
        ax=ax,
    )
    ax.set_title(
        f"Cell Type Coordination - {format_condition_name(condition_name)}\n"
        f"(Correlation of Combined Effects)",
        fontsize=14,
    )
    ax.set_xlabel("Cell Type", fontsize=12)
    ax.set_ylabel("Cell Type", fontsize=12)
    ax.tick_params(axis="y", labelsize=13)
    ax.tick_params(axis="x", labelsize=13, rotation=90)

    plt.tight_layout()
    plt.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")

    save_data_table(corr_matrix, output_path, data_dir, "Cell type correlation matrix")


def plot_top_genes_heatmap(
    combined_effects: pd.DataFrame,
    condition_name: str,
    output_path: Path,
    data_dir: Path,
    top_n: int = 30,
) -> None:
    """Create heatmap of top affected genes across cell types."""
    print(f"\n[{_timestamp()}] Creating top genes heatmap for {condition_name}...")

    if combined_effects.empty:
        print("  No data available for heatmap")
        return

    mean_effects = combined_effects.mean(axis=1)
    top_genes = mean_effects.nlargest(top_n).index.tolist()
    data = combined_effects.loc[top_genes]

    data_norm = (data - data.min()) / (data.max() - data.min())

    fig, ax = plt.subplots(figsize=FIGSIZE_LARGE)
    sns.heatmap(data_norm, cmap="YlOrRd", linewidths=0.5, ax=ax)
    ax.set_title(
        f"Top {top_n} Affected Genes - {format_condition_name(condition_name)}",
        fontsize=14,
    )
    ax.set_xlabel("Cell Type", fontsize=12)
    ax.set_ylabel("Gene", fontsize=12)
    ax.tick_params(axis="y", labelsize=13)
    ax.tick_params(axis="x", labelsize=13, rotation=90)

    plt.tight_layout()
    plt.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")

    export_data = data.copy()
    export_data["mean_effect"] = mean_effects.loc[top_genes]
    save_data_table(export_data, output_path, data_dir, f"Top {top_n} genes by mean effect")


def plot_condition_comparison(
    disease_effects: pd.DataFrame,
    normal_effects: pd.DataFrame,
    output_path: Path,
    data_dir: Path,
    condition_name: str = "Fibrotic",
) -> None:
    """Create comparison plot of disease vs normal effects."""
    display_name = format_condition_name(condition_name)
    print(f"\n[{_timestamp()}] Creating condition comparison plot ({display_name} vs Normal)...")

    if disease_effects.empty or normal_effects.empty:
        print("  Missing data for comparison plot")
        return

    common_celltypes = list(set(disease_effects.columns) & set(normal_effects.columns))
    common_genes = list(set(disease_effects.index) & set(normal_effects.index))

    if not common_celltypes or not common_genes:
        print("  No common cell types or genes for comparison")
        return

    n_celltypes = len(common_celltypes)
    n_cols = 3
    n_rows = max(1, (n_celltypes + n_cols - 1) // n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for i, celltype in enumerate(common_celltypes):
        ax = axes[i]
        dis_scores = disease_effects.loc[common_genes, celltype]
        norm_scores = normal_effects.loc[common_genes, celltype]

        ax.scatter(norm_scores, dis_scores, alpha=0.3, s=10)

        max_val = max(norm_scores.max(), dis_scores.max())
        ax.plot([0, max_val], [0, max_val], "r--", alpha=0.5)

        combined_score = dis_scores + norm_scores
        top_genes = combined_score.nlargest(10).index
        texts = []
        for gene in top_genes:
            texts.append(
                ax.text(norm_scores[gene], dis_scores[gene], gene, fontsize=6, alpha=0.8)
            )
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", alpha=0.5))

        corr = dis_scores.corr(norm_scores)
        ax.set_title(f"{celltype}\nr = {corr:.3f}", fontsize=11)
        ax.set_xlabel("Normal Effect", fontsize=10)
        ax.set_ylabel(f"{display_name} Effect", fontsize=10)

    for i in range(n_celltypes, len(axes)):
        axes[i].set_visible(False)

    plt.suptitle(
        f"{display_name} vs Normal Combined Effects by Cell Type", fontsize=14, y=1.02
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")

    corr_data = []
    for celltype in common_celltypes:
        dis_scores = disease_effects.loc[common_genes, celltype]
        norm_scores = normal_effects.loc[common_genes, celltype]
        corr_data.append({
            "celltype": celltype,
            "correlation": dis_scores.corr(norm_scores),
            "disease_mean": dis_scores.mean(),
            "normal_mean": norm_scores.mean(),
            "n_genes": len(common_genes),
        })
    save_data_table(
        pd.DataFrame(corr_data), output_path, data_dir, "Condition comparison correlations"
    )


def plot_differential_effects(
    disease_effects: pd.DataFrame,
    normal_effects: pd.DataFrame,
    output_path: Path,
    data_dir: Path,
    condition_name: str = "Fibrotic",
    top_n: int = 20,
) -> None:
    """Create volcano-style plot of differential effects."""
    display_name = format_condition_name(condition_name)
    print(f"\n[{_timestamp()}] Creating differential effects plot ({display_name} vs Normal)...")

    if disease_effects.empty or normal_effects.empty:
        print("  Missing data for differential plot")
        return

    common_celltypes = list(set(disease_effects.columns) & set(normal_effects.columns))
    common_genes = list(set(disease_effects.index) & set(normal_effects.index))

    if not common_celltypes or not common_genes:
        print("  No common data for differential plot")
        return

    n_celltypes = len(common_celltypes)
    n_cols = 3
    n_rows = max(1, (n_celltypes + n_cols - 1) // n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for i, celltype in enumerate(common_celltypes):
        ax = axes[i]
        dis = disease_effects.loc[common_genes, celltype]
        norm = normal_effects.loc[common_genes, celltype]

        eps = 1e-10
        log2_fc = np.log2((dis + eps) / (norm + eps))
        mean_expr = (dis + norm) / 2

        ax.scatter(mean_expr, log2_fc, alpha=0.3, s=10, c="gray")

        sig_up = log2_fc > 1
        sig_down = log2_fc < -1
        ax.scatter(mean_expr[sig_up], log2_fc[sig_up], alpha=0.5, s=15, c="red", label="Up")
        ax.scatter(
            mean_expr[sig_down], log2_fc[sig_down], alpha=0.5, s=15, c="blue", label="Down"
        )

        ax.axhline(y=1, color="red", linestyle="--", alpha=0.3)
        ax.axhline(y=-1, color="blue", linestyle="--", alpha=0.3)
        ax.axhline(y=0, color="black", linestyle="-", alpha=0.3)

        ax.set_title(f"{celltype}", fontsize=11)
        ax.set_xlabel("Mean Effect", fontsize=10)
        ax.set_ylabel(f"log2({display_name}/Normal)", fontsize=10)

        top_up = log2_fc.nlargest(3)
        top_down = log2_fc.nsmallest(3)
        for gene in top_up.index:
            ax.annotate(gene, (mean_expr[gene], log2_fc[gene]), fontsize=6)
        for gene in top_down.index:
            ax.annotate(gene, (mean_expr[gene], log2_fc[gene]), fontsize=6)

    for i in range(n_celltypes, len(axes)):
        axes[i].set_visible(False)

    plt.suptitle(
        f"Differential Effects: {display_name} vs Normal\n(log2 fold change)",
        fontsize=14,
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")

    diff_data = []
    eps = 1e-10
    for celltype in common_celltypes:
        dis = disease_effects.loc[common_genes, celltype]
        norm = normal_effects.loc[common_genes, celltype]
        log2_fc = np.log2((dis + eps) / (norm + eps))
        for gene in common_genes:
            diff_data.append({
                "gene": gene,
                "celltype": celltype,
                "disease_effect": dis[gene],
                "normal_effect": norm[gene],
                "log2_fc": log2_fc[gene],
            })
    save_data_table(
        pd.DataFrame(diff_data),
        output_path,
        data_dir,
        "Differential effects by gene and cell type",
    )


def plot_ccc_network(
    ccc_df: pd.DataFrame,
    condition_name: str,
    output_path: Path,
    data_dir: Path,
    top_n: int = 50,
) -> None:
    """Create heatmap of cell-cell communication strength."""
    print(f"\n[{_timestamp()}] Creating CCC network plot for {condition_name}...")

    if ccc_df.empty:
        print("  No CCC data available")
        return

    ccc_agg = ccc_df.groupby(["celltype_source", "celltype_target"])["lr_means"].sum().reset_index()
    celltypes = sorted(list(set(ccc_agg["celltype_source"]) | set(ccc_agg["celltype_target"])))
    ccc_matrix = pd.DataFrame(0.0, index=celltypes, columns=celltypes)

    for _, row in ccc_agg.iterrows():
        ccc_matrix.loc[row["celltype_source"], row["celltype_target"]] = row["lr_means"]

    fig, ax = plt.subplots(figsize=FIGSIZE_MEDIUM)
    sns.heatmap(
        ccc_matrix,
        cmap="YlOrRd",
        square=True,
        linewidths=0.5,
        ax=ax,
        annot=True,
        fmt=".1f",
    )
    ax.set_title(
        f"Cell-Cell Communication Strength - {format_condition_name(condition_name)}",
        fontsize=14,
    )
    ax.set_xlabel("Receiver Cell Type", fontsize=12)
    ax.set_ylabel("Sender Cell Type", fontsize=12)
    ax.tick_params(axis="y", labelsize=13)
    ax.tick_params(axis="x", labelsize=13, rotation=90)

    plt.tight_layout()
    plt.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")

    save_data_table(ccc_matrix, output_path, data_dir, "Cell-cell communication matrix")


def plot_fibrosis_ligand_effects(
    combined_effects: pd.DataFrame,
    condition_name: str,
    output_path: Path,
    data_dir: Path,
) -> None:
    """Create bar plot of top fibrosis/seed ligand downstream effects."""
    print(f"\n[{_timestamp()}] Creating fibrosis ligand effects plot for {condition_name}...")

    if combined_effects.empty:
        print("  No data available")
        return

    fibrosis_genes = [
        "COL1A1", "COL1A2", "COL3A1", "FN1", "ELN", "VIM",
        "CTGF", "SERPINE1", "ACTA2", "TGFBI",
        "IL6", "IL1B", "TNF", "CCL2", "CXCL10",
        "PDGFRA", "PDGFRB", "FAP",
    ]

    available_genes = [g for g in fibrosis_genes if g in combined_effects.index]
    if not available_genes:
        print("  No fibrosis marker genes found in data")
        return

    data = combined_effects.loc[available_genes]

    fig, ax = plt.subplots(figsize=FIGSIZE_LARGE)
    x = np.arange(len(available_genes))
    n_celltypes = len(data.columns)
    width = 0.8 / max(n_celltypes, 1)

    for i, celltype in enumerate(data.columns):
        ax.bar(x + i * width, data[celltype], width, label=celltype)

    ax.set_xlabel("Gene", fontsize=12)
    ax.set_ylabel("Combined Effect Score", fontsize=12)
    ax.set_title(
        f"Fibrosis Marker Gene Effects - {format_condition_name(condition_name)}", fontsize=14
    )
    ax.set_xticks(x + width * (n_celltypes - 1) / 2)
    ax.set_xticklabels(available_genes, rotation=45, ha="right")
    ax.legend(title="Cell Type", bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.tight_layout()
    plt.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")

    save_data_table(data, output_path, data_dir, "Fibrosis marker gene effects")


def create_summary_figure(
    disease_effects: pd.DataFrame,
    normal_effects: pd.DataFrame,
    ccc_disease: pd.DataFrame,
    ccc_normal: pd.DataFrame,
    output_path: Path,
    data_dir: Path,
    condition_name: str = "Disease",
) -> None:
    """Create a summary figure with multiple panels."""
    display_name = format_condition_name(condition_name)
    print(f"\n[{_timestamp()}] Creating summary figure for {display_name}...")

    fig = plt.figure(figsize=(18, 12))

    # Panel A: Disease coordination
    if not disease_effects.empty:
        ax1 = fig.add_subplot(2, 3, 1)
        corr_dis = disease_effects.corr()
        sns.heatmap(corr_dis, cmap="RdBu_r", center=0, ax=ax1, cbar_kws={"shrink": 0.5})
        ax1.set_title(f"A) {display_name} Cell Coordination", fontsize=11)

    # Panel B: Normal coordination
    if not normal_effects.empty:
        ax2 = fig.add_subplot(2, 3, 2)
        corr_norm = normal_effects.corr()
        sns.heatmap(corr_norm, cmap="RdBu_r", center=0, ax=ax2, cbar_kws={"shrink": 0.5})
        ax2.set_title("B) Normal Cell Coordination", fontsize=11)

    # Panel C: Coordination difference
    if not disease_effects.empty and not normal_effects.empty:
        ax3 = fig.add_subplot(2, 3, 3)
        common_ct = list(set(disease_effects.columns) & set(normal_effects.columns))
        if common_ct:
            corr_diff = disease_effects[common_ct].corr() - normal_effects[common_ct].corr()
            sns.heatmap(corr_diff, cmap="RdBu_r", center=0, ax=ax3, cbar_kws={"shrink": 0.5})
        ax3.set_title(f"C) Coordination Change ({display_name} - Norm)", fontsize=11)

    # Panel D: Mean effects by cell type - Disease
    if not disease_effects.empty:
        ax4 = fig.add_subplot(2, 3, 4)
        mean_dis = disease_effects.mean()
        mean_dis.plot(kind="bar", ax=ax4, color="firebrick", alpha=0.7)
        ax4.set_title(f"D) Mean Effect by Cell Type ({display_name})", fontsize=11)
        ax4.set_ylabel("Mean Effect")
        ax4.tick_params(axis="x", rotation=90)

    # Panel E: Mean effects by cell type - Normal
    if not normal_effects.empty:
        ax5 = fig.add_subplot(2, 3, 5)
        mean_norm = normal_effects.mean()
        mean_norm.plot(kind="bar", ax=ax5, color="steelblue", alpha=0.7)
        ax5.set_title("E) Mean Effect by Cell Type (Normal)", fontsize=11)
        ax5.set_ylabel("Mean Effect")
        ax5.tick_params(axis="x", rotation=90)

    # Panel F: CCC interaction counts
    ax6 = fig.add_subplot(2, 3, 6)
    ccc_counts = pd.DataFrame(
        {
            display_name: [len(ccc_disease)] if not ccc_disease.empty else [0],
            "Normal": [len(ccc_normal)] if not ccc_normal.empty else [0],
        },
        index=["Interactions"],
    )
    ccc_counts.T.plot(kind="bar", ax=ax6, color=["firebrick", "steelblue"], alpha=0.7)
    ax6.set_title("F) CCC Interactions by Condition", fontsize=11)
    ax6.set_ylabel("Number of Interactions")
    ax6.tick_params(axis="x", rotation=0)
    ax6.legend().remove()

    plt.suptitle(
        f"{display_name} ReCoN Analysis Summary", fontsize=16, y=1.02
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")

    summary_data = []
    if not disease_effects.empty:
        for ct in disease_effects.columns:
            summary_data.append({
                "condition": condition_name,
                "celltype": ct,
                "mean_effect": disease_effects[ct].mean(),
                "std_effect": disease_effects[ct].std(),
            })
    if not normal_effects.empty:
        for ct in normal_effects.columns:
            summary_data.append({
                "condition": "normal",
                "celltype": ct,
                "mean_effect": normal_effects[ct].mean(),
                "std_effect": normal_effects[ct].std(),
            })
    if summary_data:
        save_data_table(
            pd.DataFrame(summary_data), output_path, data_dir, "Summary statistics"
        )


# =============================================================================
# CASCADE STATS LOADING FOR SANKEY HOVER
# =============================================================================

def load_cascade_stats_index(config: ReconConfig) -> Dict[str, Dict[str, dict]]:
    """
    Load edge-level differential statistics for Sankey hover.

    Returns dict keyed by disease condition (lowercase), each containing
    edge-key -> stats mappings.
    """
    cascade_dir = config.get_cascade_dir()
    indices: Dict[str, Dict[str, dict]] = {}

    for disease in config.disease_conditions:
        d_lower = disease.lower()
        edge_path = cascade_dir / d_lower / "edge_results.csv"
        cascade_path = cascade_dir / d_lower / "cascade_results.csv"

        if edge_path.exists():
            print(f"  Loading edge stats for {disease.upper()} from {edge_path.name}...")
            indices[d_lower] = _load_edge_results_index(edge_path)
            print(f"    Loaded {len(indices[d_lower]):,} edges with edge-level statistics")
        elif cascade_path.exists():
            print(f"  Warning: edge_results.csv not found, falling back to cascade_results.csv")
            print(f"  Loading cascade stats for {disease.upper()} from {cascade_path.name}...")
            indices[d_lower] = _build_edge_cascade_index_chunked(cascade_path)
            print(f"    Built index with {len(indices[d_lower]):,} unique edges (cascade-level)")
        else:
            print(f"  Warning: No differential results found for {disease}")
            indices[d_lower] = {}

    return indices


def _load_edge_results_index(csv_path: Path) -> Dict[str, dict]:
    """Load edge-level statistics from edge_results.csv."""
    df = pd.read_csv(csv_path)
    index = {}
    for _, row in df.iterrows():
        edge_key = row["edge_key"]
        index[edge_key] = {
            "diff": row["diff"],
            "pval": row["pval"],
            "padj": row["padj"],
            "weight_disease": row.get("weight_disease", float("nan")),
            "weight_normal": row.get("weight_normal", float("nan")),
            "edge_type": row.get("edge_type", "unknown"),
        }
    return index


def _build_edge_cascade_index_chunked(
    csv_path: Path, chunksize: int = 500000
) -> Dict[str, dict]:
    """Fallback: Build edge index from cascade_results.csv (cascade-level stats)."""
    index: Dict[str, dict] = {}
    cols = [
        "cascade_id", "ligand", "receptor", "tf", "gene",
        "cell_source", "cell_target", "diff", "pval", "padj",
    ]

    chunk_num = 0
    for chunk in pd.read_csv(csv_path, usecols=cols, chunksize=chunksize):
        chunk_num += 1
        if chunk_num % 5 == 0:
            print(f"    Processing chunk {chunk_num}...")

        for _, row in chunk.iterrows():
            ccc_key = f"{row['ligand']}::{row['cell_source']}|{row['receptor']}::{row['cell_target']}|ccc"
            rtf_key = f"{row['receptor']}::{row['cell_target']}|{row['tf']}::{row['cell_target']}|rtf"
            grn_key = f"{row['tf']}::{row['cell_target']}|{row['gene']}::{row['cell_target']}|grn"

            stats = {
                "diff": row["diff"],
                "pval": row["pval"],
                "padj": row["padj"],
                "cascade_id": row["cascade_id"],
            }

            for key in [ccc_key, rtf_key, grn_key]:
                if key not in index or row["pval"] < index[key]["pval"]:
                    index[key] = stats

    return index


# =============================================================================
# SANKEY DIAGRAM FUNCTIONS
# =============================================================================

def _get_all_seeds(config: ReconConfig) -> List[str]:
    """Get flattened list of all seed genes from config.seed_categories."""
    all_seeds: List[str] = []
    for genes in config.seed_categories.values():
        all_seeds.extend(genes)
    return list(set(all_seeds))


def _detect_cell_types(config: ReconConfig) -> List[str]:
    """Auto-detect cell types from GRN directory."""
    grn_dir = config.get_grn_dir()
    celltypes = set()
    for p in grn_dir.glob("*_rna_network.csv"):
        # filename: <celltype>_<condition>_rna_network.csv
        stem = p.stem  # e.g. fibroblast_ssc_rna_network
        # Remove known condition suffixes and _rna_network
        name_part = stem.replace("_rna_network", "")
        for cond in config.conditions:
            suffix = f"_{cond.lower()}"
            if name_part.endswith(suffix):
                ct_raw = name_part[: -len(suffix)]
                celltypes.add(ct_raw)
                break
    return sorted(celltypes)


def load_grn_for_sankey(condition: str, config: ReconConfig) -> pd.DataFrame:
    """Load and merge cell-type GRNs for Sankey visualization."""
    grn_dir = config.get_grn_dir()
    all_edges = []

    cell_types = _detect_cell_types(config)
    for ct_safe in cell_types:
        grn_path = grn_dir / f"{ct_safe}_{condition.lower()}_rna_network.csv"
        if grn_path.exists():
            grn = pd.read_csv(grn_path)
            grn = grn[grn["weight"] > config.min_grn_weight][["source", "target", "weight"]]
            grn["source"] = grn["source"] + "_TF"
            all_edges.append(grn)

    if not all_edges:
        return pd.DataFrame(columns=["source", "target", "weight"])

    combined = pd.concat(all_edges, ignore_index=True)
    return combined.groupby(["source", "target"], as_index=False).agg({"weight": "max"})


def load_ccc_for_sankey(condition: str, config: ReconConfig) -> pd.DataFrame:
    """Load CCC data based on config.ccc_source for Sankey visualization."""
    ccc_dir = config.get_ccc_dir()

    if config.ccc_source == "cellchat":
        path = ccc_dir / "cellchat" / f"{condition.lower()}_ccc.csv"
        ccc = pd.read_csv(path)
        ccc = ccc.rename(columns={
            "source": "celltype_source",
            "target": "celltype_target",
            "ligand": "source",
            "receptor": "target",
        })
        ccc["lr_means"] = ccc["prob"]
    elif config.ccc_source == "cellphonedb":
        path = ccc_dir / "cellphonedb" / f"{condition.lower()}_ccc.csv"
        ccc = pd.read_csv(path)
    else:  # merged
        path = ccc_dir / "merged" / f"{condition.lower()}_ccc.csv"
        ccc = pd.read_csv(path)

    if "lr_means" in ccc.columns:
        ccc["lr_means"] = ccc["lr_means"].fillna(0)

    return ccc


def build_multicell_with_seeds(
    celltypes: List[str],
    ccc: pd.DataFrame,
    grn: pd.DataFrame,
    receptor_grn: pd.DataFrame,
    seeds: List[str],
    restart_proba: float = 0.6,
) -> tuple:
    """Build Multicell object and run propagation for Sankey visualization."""
    celltype_objs = []
    for ct in celltypes:
        celltype_objs.append(
            Celltype(
                celltype_name=ct,
                grn_graph=grn,
                receptor_grn_bipartite=receptor_grn.copy(),
                receptor_graph=None,
                grn_graph_directed=False,
                grn_graph_weighted=True,
            )
        )

    celltype_dict = {ct.celltype_name: ct for ct in celltype_objs}

    starting_nodes = [f"{seed}-{ct}" for seed in seeds for ct in celltypes]

    multicell = Multicell(
        celltypes=celltype_dict,
        cell_communication_graph=ccc,
        cell_communication_graph_directed=False,
        cell_communication_graph_weighted=True,
        seeds=starting_nodes,
        verbose=False,
    )

    multicell.lamb = set_lambda(
        multicell=multicell,
        direction="downstream",
        strategy="intercell",
    )

    print("  Running MultiXrank propagation...")
    multixrank_obj = multicell.Multixrank(restart_proba=restart_proba, verbose=False)
    results_df = multixrank_obj.random_walk_rank()

    return multicell, results_df


def save_sankey_data_tables(
    multicell: Multicell,
    results: pd.DataFrame,
    cell_type: str,
    seeds: List[str],
    ligand_cells: List[str],
    condition: str,
    output_dir: Path,
    top_ligand_n: int = 50,
    top_receptor_n: int = 30,
    top_tf_n: int = 10,
) -> Dict[str, pd.DataFrame]:
    """Extract and save data tables underlying Sankey diagrams."""
    prefix = f"{condition}_{cell_type.lower().replace(' ', '_')}"
    tables: Dict[str, pd.DataFrame] = {}

    try:
        cc_df = get_cell_communication_layer(
            multicell, as_dataframe=True,
            ligand_cells=ligand_cells, receptor_cells=[cell_type],
        )
        cc_df.to_csv(output_dir / f"{prefix}_ccc_layer.csv", index=False)
        tables["ccc"] = cc_df
        print(f"    Saved: {prefix}_ccc_layer.csv ({len(cc_df)} rows)")
    except Exception as e:
        print(f"    Warning: Could not extract CCC layer: {e}")
        return tables

    try:
        top_ligands = get_top_ligands(results, cc_df, n=top_ligand_n, per_celltype=True)
        top_ligands.to_csv(output_dir / f"{prefix}_top_ligands.csv", index=False)
        tables["top_ligands"] = top_ligands
        print(f"    Saved: {prefix}_top_ligands.csv ({len(top_ligands)} rows)")
    except Exception as e:
        print(f"    Warning: Could not extract top ligands: {e}")

    try:
        top_receptors = get_top_receptors(results, cell_type=cell_type, n=top_receptor_n)
        top_receptors.to_csv(output_dir / f"{prefix}_top_receptors.csv", index=False)
        tables["top_receptors"] = top_receptors
        print(f"    Saved: {prefix}_top_receptors.csv ({len(top_receptors)} rows)")
    except Exception as e:
        print(f"    Warning: Could not extract top receptors: {e}")

    try:
        top_tfs = get_top_tfs(results, cell_type=cell_type, n=top_tf_n)
        top_tfs.to_csv(output_dir / f"{prefix}_top_tfs.csv", index=False)
        tables["top_tfs"] = top_tfs
        print(f"    Saved: {prefix}_top_tfs.csv ({len(top_tfs)} rows)")
    except Exception as e:
        print(f"    Warning: Could not extract top TFs: {e}")

    try:
        tf_gene_df = get_celltype_gene_layer(
            multicell, cell_type=cell_type, layer_name="gene", as_dataframe=True
        )
        tf_gene_df.to_csv(output_dir / f"{prefix}_tf_gene_layer.csv", index=False)
        tables["tf_gene"] = tf_gene_df
        print(f"    Saved: {prefix}_tf_gene_layer.csv ({len(tf_gene_df)} rows)")
    except Exception as e:
        print(f"    Warning: Could not extract TF-gene layer: {e}")

    try:
        receptor_tf_df = get_celltype_grn_receptor_bipartite(
            multicell, cell_type=cell_type, as_dataframe=True
        )
        receptor_tf_df.to_csv(output_dir / f"{prefix}_receptor_tf_layer.csv", index=False)
        tables["receptor_tf"] = receptor_tf_df
        print(f"    Saved: {prefix}_receptor_tf_layer.csv ({len(receptor_tf_df)} rows)")
    except Exception as e:
        print(f"    Warning: Could not extract receptor-TF layer: {e}")

    return tables


def fetch_ppi_partners(
    seed_genes: List[str],
    output_dir: Path,
    grn_tf_set: set,
    ppi_min_score: int = 400,
) -> pd.DataFrame:
    """
    Fetch PPI partners for seed genes from STRING API.

    Queries the STRING database for physical protein-protein interactions
    involving any of the seed genes. Results are cached to CSV.

    Args:
        seed_genes: List of gene symbols to query.
        output_dir: Directory for the cache file (ppi_partners.csv).
        grn_tf_set: Set of TF names from GRN (to classify partners).
        ppi_min_score: Minimum combined score for STRING interactions.

    Returns:
        DataFrame with columns [seed, partner, combined_score, is_tf].
    """
    cache_path = Path(output_dir) / "ppi_partners.csv"
    empty_cols = ["seed", "partner", "combined_score", "is_tf"]

    if cache_path.exists():
        print(f"  Loading cached PPI data from {cache_path}")
        return pd.read_csv(cache_path)

    if not seed_genes:
        print("  Warning: No seed genes provided for PPI query.")
        return pd.DataFrame(columns=empty_cols)

    print(f"  Fetching PPI partners for {len(seed_genes)} seed genes from STRING API...")
    try:
        import requests
    except ImportError:
        print("  Warning: 'requests' package not available. Skipping PPI fetch.")
        return pd.DataFrame(columns=empty_cols)

    try:
        url = "https://string-db.org/api/json/network"
        params = {
            "identifiers": "|".join(seed_genes),
            "species": 9606,
            "required_score": ppi_min_score,
            "network_type": "physical",
            "caller_identity": "recon_multinetwork",
        }
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        seed_set = set(seed_genes)
        rows = []
        for item in data:
            pA = item.get("preferredName_A", "")
            pB = item.get("preferredName_B", "")
            score = item.get("score", 0)

            if pA in seed_set and pB not in seed_set:
                rows.append({"seed": pA, "partner": pB, "combined_score": score})
            elif pB in seed_set and pA not in seed_set:
                rows.append({"seed": pB, "partner": pA, "combined_score": score})
            elif pA in seed_set and pB in seed_set:
                rows.append({"seed": pA, "partner": pB, "combined_score": score})
                rows.append({"seed": pB, "partner": pA, "combined_score": score})

        if not rows:
            print("  Warning: No PPI partners found. Returning empty DataFrame.")
            ppi_df = pd.DataFrame(columns=empty_cols)
        else:
            ppi_df = pd.DataFrame(rows)
            ppi_df = ppi_df.drop_duplicates(subset=["seed", "partner"])
            ppi_df["is_tf"] = ppi_df["partner"].isin(grn_tf_set)
            ppi_df = ppi_df.sort_values("combined_score", ascending=False)
            print(
                f"  Found {len(ppi_df)} PPI partner links "
                f"({ppi_df['is_tf'].sum()} partners are TFs)"
            )

        ppi_df.to_csv(cache_path, index=False)
        return ppi_df

    except Exception as e:
        print(f"  Warning: STRING API failed: {e}")
        print("  Returning empty PPI DataFrame.")
        return pd.DataFrame(columns=empty_cols)


def build_5layer_ppi_links(
    condition: str,
    cell_type: str,
    ppi_df: pd.DataFrame,
    config: ReconConfig,
) -> Optional[Dict[str, pd.DataFrame]]:
    """
    Build 5-layer PPI Sankey links:
    Ligand::SrcCell -> Seed::Cell -> PPI_Partner::Cell -> TF::Cell -> Gene::Cell

    Generalised from the CDH11 reference for any set of seed genes.

    Args:
        condition: Condition name (lowercase).
        cell_type: Focal cell type (display name matching CCC/GRN).
        ppi_df: DataFrame from fetch_ppi_partners (seed, partner, combined_score, is_tf).
        config: ReconConfig with seed_categories, min_sankey_grn_weight, etc.

    Returns:
        Dict with keys ligand_receptor, receptor_ppi, ppi_tf, tf_gene,
        or None if no links can be built.
    """
    if ppi_df is None or len(ppi_df) == 0:
        return None

    seeds = config.load_seeds()
    if not seeds:
        seeds = _get_all_seeds(config)
    seed_set = set(seeds)

    # Load CCC filtered to target cell
    ccc = load_ccc_for_sankey(condition, config)
    ccc_to_cell = ccc[ccc["celltype_target"] == cell_type]
    if len(ccc_to_cell) == 0:
        return None

    # Filter CCC to interactions where receptor is a seed gene
    ccc_seed = ccc_to_cell[ccc_to_cell["target"].isin(seed_set)]
    if len(ccc_seed) == 0:
        return None

    # Layer 1: Ligand::SrcCell -> Seed::Cell  (CCC)
    ligand_receptor = pd.DataFrame({
        "source": ccc_seed["source"] + "::" + ccc_seed["celltype_source"],
        "target": ccc_seed["target"] + "::" + cell_type,
        "weight": ccc_seed["lr_means"].values,
    })

    # Filter PPI to seeds that appear in this CCC
    active_seeds = set(ccc_seed["target"].unique())
    ppi_active = ppi_df[ppi_df["seed"].isin(active_seeds)].copy()
    if len(ppi_active) == 0:
        return None

    # Layer 2: Seed::Cell -> PPI_Partner::Cell
    receptor_ppi = pd.DataFrame({
        "source": ppi_active["seed"] + "::" + cell_type,
        "target": ppi_active["partner"] + "::" + cell_type,
        "weight": ppi_active["combined_score"].values / 1000.0,
    })

    # Load GRN for this cell type
    ct_safe = cell_type.lower().replace(" ", "_")
    grn_dir = config.get_grn_dir()
    grn_path = grn_dir / f"{ct_safe}_{condition}_rna_network.csv"
    if not grn_path.exists():
        return None

    grn = pd.read_csv(grn_path)
    grn = grn[grn["weight"] >= config.min_sankey_grn_weight]
    grn_tfs = set(grn["source"].unique())

    # Load NicheNet PKN for receptor -> TF links
    receptor_grn = load_receptor_genes("human_receptor_gene_from_NichenetPKN")

    # Layers 3 & 4: PPI partner -> TF -> Gene
    ppi_to_tf_rows = []
    tf_to_gene_rows = []

    for _, prow in ppi_active.iterrows():
        partner = prow["partner"]
        if partner in grn_tfs:
            # Partner IS a TF: connect partner -> partner_TF -> genes
            partner_genes = grn[grn["source"] == partner]
            partner_genes = partner_genes.sort_values(
                "weight", ascending=False
            ).head(config.top_grn_genes_sankey)
            for _, grow in partner_genes.iterrows():
                ppi_to_tf_rows.append({
                    "source": f"{partner}::{cell_type}",
                    "target": f"{partner}_TF::{cell_type}",
                    "weight": prow["combined_score"] / 1000.0,
                })
                tf_to_gene_rows.append({
                    "source": f"{partner}_TF::{cell_type}",
                    "target": f"{grow['target']}::{cell_type}",
                    "weight": grow["weight"],
                })
        else:
            # Partner is NOT a TF: check NicheNet PKN for partner -> TF links
            partner_rtf = receptor_grn[receptor_grn["source"] == partner]
            partner_rtf = partner_rtf[partner_rtf["target"].isin(grn_tfs)]
            partner_rtf = partner_rtf[
                partner_rtf["weight"] >= config.min_rtf_weight
            ]
            partner_rtf = partner_rtf.sort_values(
                "weight", ascending=False
            ).head(config.top_tfs_sankey)
            for _, rtf_row in partner_rtf.iterrows():
                tf = rtf_row["target"]
                ppi_to_tf_rows.append({
                    "source": f"{partner}::{cell_type}",
                    "target": f"{tf}_TF::{cell_type}",
                    "weight": rtf_row["weight"],
                })
                tf_genes = grn[grn["source"] == tf].sort_values(
                    "weight", ascending=False
                ).head(config.top_grn_genes_sankey)
                for _, grow in tf_genes.iterrows():
                    tf_to_gene_rows.append({
                        "source": f"{tf}_TF::{cell_type}",
                        "target": f"{grow['target']}::{cell_type}",
                        "weight": grow["weight"],
                    })

    if not ppi_to_tf_rows:
        return None

    ppi_tf = pd.DataFrame(ppi_to_tf_rows).drop_duplicates(
        subset=["source", "target"]
    )
    tf_gene = pd.DataFrame(tf_to_gene_rows).drop_duplicates(
        subset=["source", "target"]
    )

    # Filter: only keep PPI partners that connect downstream
    connected_ppi = set(ppi_tf["source"].unique())
    receptor_ppi = receptor_ppi[receptor_ppi["target"].isin(connected_ppi)]

    if len(receptor_ppi) == 0:
        return None

    # Filter: only keep seeds that connect to active PPI partners
    connected_seeds = set(receptor_ppi["source"].unique())
    ligand_receptor = ligand_receptor[
        ligand_receptor["target"].isin(connected_seeds)
    ]

    if len(ligand_receptor) == 0:
        return None

    return {
        "ligand_receptor": ligand_receptor,
        "receptor_ppi": receptor_ppi,
        "ppi_tf": ppi_tf,
        "tf_gene": tf_gene,
    }


def plot_5layer_ppi_sankey(
    links: Dict[str, pd.DataFrame],
    condition: str,
    cell_type: str,
    cascade_stats: Dict[str, Dict[str, dict]],
    config: ReconConfig,
    save_path: Path,
) -> None:
    """
    Render a 5-layer PPI Sankey diagram using Plotly.

    Layers: Ligands -> Receptors/Seeds -> PPI Partners -> TFs -> Genes

    Colors:
        - Ligand->Receptor (CCC): gray rgba(160,160,160,0.4)
        - Receptor->PPI (PPI): steelblue rgba(70,130,180,0.7)
        - PPI->TF: steelblue rgba(70,130,180,0.5)
        - TF->Gene (GRN): green rgba(100,200,100,0.6)

    Args:
        links: Dict with keys ligand_receptor, receptor_ppi, ppi_tf, tf_gene.
        condition: Condition name (lowercase).
        cell_type: Focal cell type.
        cascade_stats: Dict of {disease_lower: {edge_key: stats}}.
        config: ReconConfig instance.
        save_path: Output HTML path.
    """
    layer_order = ["ligand_receptor", "receptor_ppi", "ppi_tf", "tf_gene"]
    layer_colors = {
        "ligand_receptor": "rgba(160,160,160,0.4)",
        "receptor_ppi": "rgba(70,130,180,0.7)",
        "ppi_tf": "rgba(70,130,180,0.5)",
        "tf_gene": "rgba(100,200,100,0.6)",
    }
    layer_edge_types = {
        "ligand_receptor": "ccc",
        "receptor_ppi": "ppi",
        "ppi_tf": "ppi_tf",
        "tf_gene": "grn",
    }

    present_layers = [
        l for l in layer_order if l in links and len(links[l]) > 0
    ]
    if not present_layers:
        print(f"    Warning: No links for PPI Sankey {condition} {cell_type}")
        return

    def hex_to_rgba(hex_color, alpha=0.6):
        h = hex_color.lstrip("#")
        return f"rgba({int(h[0:2], 16)}, {int(h[2:4], 16)}, {int(h[4:6], 16)}, {alpha})"

    def string_to_color(string):
        h = hashlib.md5(string.encode()).hexdigest()
        return "#" + h[:6]

    all_link_dfs = []
    for layer_name in present_layers:
        df = links[layer_name].copy()
        color_spec = layer_colors.get(layer_name, "rgba(160,160,160,0.4)")
        # For CCC layer, color by source cell type
        if layer_name == "ligand_receptor":
            df["color"] = df["source"].apply(
                lambda x: hex_to_rgba(
                    string_to_color(_extract_celltype(x))
                )
            )
        else:
            df["color"] = color_spec
        df["layer_type"] = layer_edge_types[layer_name]

        total = df["weight"].sum()
        if total > 0:
            df["value"] = df["weight"] / total
        else:
            df["value"] = df["weight"]

        all_link_dfs.append(df)

    all_links = pd.concat(all_link_dfs, ignore_index=True)

    all_nodes = pd.unique(all_links[["source", "target"]].values.ravel())
    node_idx = {name: i for i, name in enumerate(all_nodes)}
    all_links["source_idx"] = all_links["source"].map(node_idx)
    all_links["target_idx"] = all_links["target"].map(node_idx)

    labels = [_format_label(n) for n in all_nodes]
    node_celltypes = [_extract_celltype(n) for n in all_nodes]

    disease_conditions = [d.lower() for d in config.disease_conditions]

    # Build customdata for hover
    link_customdata = []
    for _, row in all_links.iterrows():
        src_gene = _format_label(row["source"])
        src_ct = _extract_celltype(row["source"])
        tgt_gene = _format_label(row["target"])
        tgt_ct = _extract_celltype(row["target"])
        layer_type = row["layer_type"]

        src_norm = _normalize_for_lookup(row["source"])
        tgt_norm = _normalize_for_lookup(row["target"])
        edge_key = f"{src_norm}|{tgt_norm}|{layer_type}"

        entry = [src_gene, src_ct, tgt_gene, tgt_ct, layer_type]
        for d_lower in disease_conditions:
            d_stats = cascade_stats.get(d_lower, {}).get(edge_key, {})
            entry.extend([
                d_stats.get("diff", float("nan")),
                d_stats.get("pval", float("nan")),
                d_stats.get("padj", float("nan")),
            ])
        link_customdata.append(entry)

    # Build hover template
    hover_parts = [
        "<b>%{customdata[0]}</b> (%{customdata[1]}) -> "
        "<b>%{customdata[2]}</b> (%{customdata[3]})<br>",
        "Weight: %{value:.4f}<br>",
        "<br>",
        "<b>Edge Type:</b> %{customdata[4]}<br>",
    ]
    for i, d_lower in enumerate(disease_conditions):
        d_display = format_condition_name(d_lower)
        base_idx = 5 + i * 3
        hover_parts.append(
            f"<b>{d_display} vs Normal:</b> "
            f"diff=%{{customdata[{base_idx}]:.4f}}, "
            f"pval=%{{customdata[{base_idx + 1}]:.2e}}, "
            f"padj=%{{customdata[{base_idx + 2}]:.2e}}<br>"
        )
    hover_parts.append("<extra></extra>")
    link_hovertemplate = "".join(hover_parts)

    sankey_data = go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            customdata=node_celltypes,
            hovertemplate=(
                "<b>%{label}</b><br>Cell Type: %{customdata}<extra></extra>"
            ),
        ),
        link=dict(
            source=all_links["source_idx"],
            target=all_links["target_idx"],
            value=all_links["value"],
            color=all_links["color"],
            customdata=link_customdata,
            hovertemplate=link_hovertemplate,
        ),
        orientation="h",
    )

    fig = go.Figure(data=[sankey_data])

    cond_display = format_condition_name(condition)
    title = f"PPI 5-Layer Sankey | {cond_display} | {cell_type}"

    fig.update_layout(
        title_text=title,
        font_size=14,
        font_color="black",
        height=1200,
        width=2000,
        margin=dict(l=50, r=200, t=50, b=100),
    )

    # Layer labels at bottom
    layer_labels = ["Ligands", "Receptors", "PPI Partners", "TFs", "Genes"]
    x_positions = [i / (len(layer_labels) - 1) for i in range(len(layer_labels))]
    for x, name in zip(x_positions, layer_labels):
        fig.add_annotation(
            x=x,
            y=-0.08,
            text=f"<b>{name}</b>",
            showarrow=False,
            font=dict(size=16),
            xref="paper",
            yref="paper",
        )

    # Cell type legend
    all_cts = set()
    for df in all_link_dfs:
        for node_col in ["source", "target"]:
            all_cts.update(
                df[node_col].apply(_extract_celltype).unique()
            )
    all_cts.discard("N/A")

    for i, ct in enumerate(sorted(all_cts)):
        fig.add_annotation(
            x=1.02,
            y=1.0 - i * 0.05,
            text=f"<b>{ct}</b>",
            showarrow=False,
            font=dict(size=12),
            bgcolor=hex_to_rgba(string_to_color(ct)),
            bordercolor="black",
            borderwidth=0.5,
            align="left",
            xanchor="left",
        )

    fig.write_html(str(save_path))

    try:
        fig.show()
    except Exception:
        pass


def plot_6layer_sankey_with_hover(
    before_receptor_tf_df: pd.DataFrame,
    before_tf_ligand_df: pd.DataFrame,
    ligand_receptor_df: pd.DataFrame,
    receptor_tf_df: pd.DataFrame,
    gene_tf_df: pd.DataFrame,
    flow: str = "upstream",
    save_path=None,
    cascade_stats: Optional[Dict[str, Dict[str, dict]]] = None,
    disease_conditions: Optional[List[str]] = None,
):
    """
    Patched version of plot_6layer_sankey with hover tooltips showing cell type info
    and differential cascade statistics for all disease conditions.

    Args:
        cascade_stats: Dict of {disease_lower: {edge_key: stats_dict}}
        disease_conditions: List of disease condition names (lowercase)
    """
    if cascade_stats is None:
        cascade_stats = {}
    if disease_conditions is None:
        disease_conditions = sorted(cascade_stats.keys())

    def format_links(df, source_col, target_col):
        return df.loc[:, [source_col, target_col, "weight"]].rename(
            columns={source_col: "source", target_col: "target", "weight": "value"}
        )

    def hex_to_rgba(hex_color, alpha=0.6):
        h = hex_color.lstrip("#")
        return f"rgba({int(h[0:2], 16)}, {int(h[2:4], 16)}, {int(h[4:6], 16)}, {alpha})"

    def string_to_color(string):
        h = hashlib.md5(string.encode()).hexdigest()
        return "#" + h[:6]

    br_bt = format_links(before_receptor_tf_df, "receptor", "tf")
    bt_l = format_links(before_tf_ligand_df, "tf_clean", "gene")
    l_r = format_links(ligand_receptor_df, "ligand", "receptor_clean")
    r_t = format_links(receptor_tf_df, "receptor", "tf")
    t_g = format_links(gene_tf_df, "tf_clean", "gene")

    # Filter unconnected
    l_r = l_r[l_r["target"].isin(r_t["source"])]
    l_r = l_r[l_r["source"].isin(bt_l["target"])]
    r_t = r_t[r_t["source"].isin(l_r["target"])]
    t_g = t_g[t_g["source"].isin(r_t["target"])]
    bt_l = bt_l[bt_l["target"].isin(l_r["source"])]
    bt_l = bt_l[bt_l["source"].isin(br_bt["target"])]
    l_r = l_r[l_r["source"].isin(bt_l["target"])]
    br_bt = br_bt[br_bt["target"].isin(bt_l["source"])]

    if flow.lower() == "downstream":
        for df in [br_bt, bt_l, l_r, r_t, t_g]:
            df.loc[:, ["source", "target"]] = df.loc[:, ["target", "source"]]

    def assign_group_colors(df, column):
        unique_types = df[column].str.extract(r"::(.+)$")[0].fillna("Unknown")
        return unique_types.apply(lambda ct: hex_to_rgba(string_to_color(ct)))

    if len(br_bt) > 0:
        br_bt.loc[:, "color"] = assign_group_colors(br_bt, "source")
        br_bt.loc[:, "layer_type"] = "upstream_rtf"
    if len(bt_l) > 0:
        bt_l.loc[:, "color"] = assign_group_colors(bt_l, "source")
        bt_l.loc[:, "layer_type"] = "upstream_grn"
    if len(l_r) > 0:
        l_r.loc[:, "color"] = "rgba(160,160,160,0.4)"
        l_r.loc[:, "layer_type"] = "ccc"
    if len(r_t) > 0:
        r_t.loc[:, "color"] = "rgba(100,200,100,1)"
        r_t.loc[:, "layer_type"] = "rtf"
    if len(t_g) > 0:
        t_g.loc[:, "color"] = "rgba(100,200,100,1)"
        t_g.loc[:, "layer_type"] = "grn"

    for df in [br_bt, bt_l, l_r, r_t, t_g]:
        total = df["value"].sum()
        if total > 0:
            df["value"] /= total

    non_empty_dfs = [df for df in [br_bt, bt_l, l_r, r_t, t_g] if len(df) > 0]
    if not non_empty_dfs:
        print("  Warning: All networks are empty. Cannot generate Sankey plot.")
        return

    links = pd.concat(non_empty_dfs, ignore_index=True)
    all_nodes = pd.unique(links[["source", "target"]].values.ravel())
    node_idx = {name: i for i, name in enumerate(all_nodes)}
    links.loc[:, "source_idx"] = links.loc[:, "source"].map(node_idx)
    links.loc[:, "target_idx"] = links.loc[:, "target"].map(node_idx)

    def _format_label(x: str) -> str:
        parts = x.split("::", 1)
        return parts[0].split("_")[0] if len(parts) == 2 else x.split("_")[0]

    def _extract_celltype(x: str) -> str:
        parts = x.split("::", 1)
        return parts[1] if len(parts) == 2 else "N/A"

    labels = [_format_label(n) for n in all_nodes]
    node_celltypes = [_extract_celltype(n) for n in all_nodes]

    def _normalize_for_lookup(node: str) -> str:
        parts = node.split("::", 1)
        if len(parts) == 2:
            gene = parts[0]
            celltype = parts[1]
            for suffix in ["_TF", "_receptor", "_clean"]:
                if gene.endswith(suffix):
                    gene = gene[: -len(suffix)]
                    break
            return f"{gene}::{celltype}"
        return node

    # Build customdata: [src_gene, src_ct, tgt_gene, tgt_ct, edge_type,
    #                     d1_diff, d1_pval, d1_padj, d2_diff, d2_pval, d2_padj, ...]
    link_customdata = []
    for _, row in links.iterrows():
        src = row["source"]
        tgt = row["target"]
        src_gene = _format_label(src)
        src_ct = _extract_celltype(src)
        tgt_gene = _format_label(tgt)
        tgt_ct = _extract_celltype(tgt)

        layer_type = row.get("layer_type", "unknown")
        src_norm = _normalize_for_lookup(src)
        tgt_norm = _normalize_for_lookup(tgt)
        edge_key = f"{src_norm}|{tgt_norm}|{layer_type}"

        entry = [src_gene, src_ct, tgt_gene, tgt_ct, layer_type]
        for d_lower in disease_conditions:
            d_stats = cascade_stats.get(d_lower, {}).get(edge_key, {})
            edge_type = d_stats.get("edge_type", layer_type)
            entry.extend([
                d_stats.get("diff", float("nan")),
                d_stats.get("pval", float("nan")),
                d_stats.get("padj", float("nan")),
            ])
        link_customdata.append(entry)

    # Build hover template dynamically for all diseases
    hover_parts = [
        "<b>%{customdata[0]}</b> (%{customdata[1]}) -> <b>%{customdata[2]}</b> (%{customdata[3]})<br>",
        "Weight: %{value:.4f}<br>",
        "<br>",
        "<b>Edge Type:</b> %{customdata[4]}<br>",
    ]
    for i, d_lower in enumerate(disease_conditions):
        d_display = format_condition_name(d_lower)
        base_idx = 5 + i * 3
        hover_parts.append(
            f"<b>{d_display} vs Normal:</b> "
            f"diff=%{{customdata[{base_idx}]:.4f}}, "
            f"pval=%{{customdata[{base_idx + 1}]:.2e}}, "
            f"padj=%{{customdata[{base_idx + 2}]:.2e}}<br>"
        )
    hover_parts.append("<extra></extra>")
    link_hovertemplate = "".join(hover_parts)

    sankey_data = go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            customdata=node_celltypes,
            hovertemplate="<b>%{label}</b><br>Cell Type: %{customdata}<extra></extra>",
        ),
        link=dict(
            source=links.loc[:, "source_idx"],
            target=links.loc[:, "target_idx"],
            value=links.loc[:, "value"],
            color=links.loc[:, "color"],
            customdata=link_customdata,
            hovertemplate=link_hovertemplate,
        ),
        orientation="h",
    )

    fig = go.Figure(data=[sankey_data])

    color_map = (
        pd.concat([br_bt, bt_l, r_t, t_g]).loc[:, ["source", "color"]]
        .dropna()
        .drop_duplicates()
    )
    color_map["celltype"] = color_map.loc[:, "source"].str.extract(r"::(.+)$")[0]
    color_map = color_map.dropna(subset=["celltype"])
    color_map = dict(zip(color_map.loc[:, "celltype"], color_map.loc[:, "color"]))

    for i, (ct, color) in enumerate(sorted(color_map.items())):
        fig.add_annotation(
            x=1.02,
            y=1.0 - i * 0.05,
            text=f"<b>{ct}</b>",
            showarrow=False,
            font=dict(size=12),
            bgcolor=color,
            bordercolor="black",
            borderwidth=0.5,
            align="left",
            xanchor="left",
        )

    if flow.lower() == "upstream":
        title_text = "Upstream Receptor -> Upstream TF -> Ligand -> Downstream Receptor -> Downstream TF -> Gene"
        layer_names = [
            "Upstream Receptors", "Upstream TFs", "Ligands",
            "Receptors", "TFs", "Genes",
        ]
    else:
        title_text = "Gene -> First TF -> First Receptor -> Ligand -> Upstream TF -> Upstream Receptor"
        layer_names = [
            "Genes", "TFs", "Receptors",
            "Ligands", "Upstream TFs", "Upstream Receptors",
        ]

    fig.update_layout(
        title_text=title_text,
        font_size=14,
        font_color="black",
        height=1200,
        width=2000,
        margin=dict(l=50, r=200, t=50, b=100),
    )

    x_positions = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    for x, name in zip(x_positions, layer_names):
        fig.add_annotation(
            x=x,
            y=-0.08,
            text=f"<b>{name}</b>",
            showarrow=False,
            font=dict(size=16),
            xref="paper",
            yref="paper",
        )

    if save_path:
        fig.write_html(save_path)

    try:
        fig.show()
    except Exception:
        pass


def generate_sankey_for_condition(
    condition: str,
    output_dir: Path,
    config: ReconConfig,
    cascade_stats: Optional[Dict[str, Dict[str, dict]]] = None,
    ppi_df: Optional[pd.DataFrame] = None,
) -> None:
    """
    Generate Sankey diagrams and data tables for a condition.

    Args:
        condition: Condition name (e.g. 'ssc', 'ipf', 'normal')
        output_dir: Output directory for Sankey files
        config: ReconConfig instance
        cascade_stats: Pre-loaded dict of {disease_lower: {edge_key: stats}}
        ppi_df: Optional PPI partners DataFrame from fetch_ppi_partners()
    """
    if cascade_stats is None:
        cascade_stats = {}

    # Monkey-patch plot_6layer_sankey with our closure-wrapped version
    disease_conditions = [d.lower() for d in config.disease_conditions]

    def plot_6layer_sankey_with_stats(
        before_receptor_tf_df, before_tf_ligand_df, ligand_receptor_df,
        receptor_tf_df, gene_tf_df, flow="upstream", save_path=None,
    ):
        return plot_6layer_sankey_with_hover(
            before_receptor_tf_df, before_tf_ligand_df, ligand_receptor_df,
            receptor_tf_df, gene_tf_df, flow=flow, save_path=save_path,
            cascade_stats=cascade_stats,
            disease_conditions=disease_conditions,
        )

    import recon.plot.sankey_paths as sankey_module
    sankey_module.plot_6layer_sankey = plot_6layer_sankey_with_stats

    print(f"\n{'=' * 60}")
    print(f"GENERATING SANKEY DIAGRAMS FOR {condition.upper()}")
    print(f"{'=' * 60}")

    data_subdir = output_dir / "data"
    data_subdir.mkdir(exist_ok=True)

    seeds = _get_all_seeds(config)
    print(f"Seeds: {len(seeds)} genes from {len(config.seed_categories)} categories")
    for cat, genes in config.seed_categories.items():
        print(f"  {cat}: {len(genes)} genes")

    grn = load_grn_for_sankey(condition, config)
    ccc = load_ccc_for_sankey(condition, config)
    receptor_grn = load_receptor_genes("human_receptor_gene_from_NichenetPKN")

    print(f"GRN edges: {len(grn):,}")
    print(f"CCC interactions: {len(ccc):,}")

    # Detect cell types from GRN files
    cell_types_raw = _detect_cell_types(config)
    # Convert to display format (title case, matching CCC)
    # We use the raw names as-is since they come from filenames
    cell_types = cell_types_raw

    print(f"\nBuilding Multicell and running propagation...")
    multicell, results = build_multicell_with_seeds(
        cell_types, ccc, grn, receptor_grn, seeds,
        restart_proba=config.restart_proba,
    )
    print(f"  Results shape: {results.shape}")
    print(f"  Result columns: {list(results.columns)}")

    focal_celltypes = config.focal_celltypes
    ligand_source_cells = config.ligand_source_cells

    for cell_type in focal_celltypes:
        print(f"\n--- {cell_type} ---")
        ct_safe = cell_type.lower().replace(" ", "_")

        print("  Extracting data tables...")
        save_sankey_data_tables(
            multicell, results, cell_type, seeds,
            ligand_source_cells, condition, data_subdir,
        )

        # 3-layer intracellular Sankey
        try:
            save_path = output_dir / f"sankey_intracell_{condition}_{ct_safe}.html"
            plot_intracell_sankey(
                multicell_obj=multicell,
                results=results,
                cell_type=cell_type,
                seeds=seeds,
                top_receptor_n=30,
                top_tf_n=10,
                flow="upstream",
                save_path=str(save_path),
            )
            print(f"  Intracellular Sankey: {save_path.name}")
        except Exception as e:
            print(f"  Intracellular Sankey failed: {e}")

        # 4-layer ligand Sankey
        try:
            save_path = output_dir / f"sankey_ligand_{condition}_{ct_safe}.html"
            plot_ligand_sankey(
                multicell_obj=multicell,
                results=results,
                cell_type=cell_type,
                seeds=seeds,
                ligand_cells=ligand_source_cells,
                top_ligand_n=50,
                top_receptor_n=30,
                top_tf_n=10,
                per_celltype=True,
                flow="upstream",
                save_path=str(save_path),
            )
            print(f"  Ligand Sankey: {save_path.name}")
        except Exception as e:
            print(f"  Ligand Sankey failed: {e}")

        # 5-layer PPI Sankey
        if ppi_df is not None and len(ppi_df) > 0 and config.ppi_min_score > 0:
            try:
                ppi_links = build_5layer_ppi_links(condition, cell_type, ppi_df, config)
                if ppi_links:
                    save_path = output_dir / f"sankey_ppi_{condition}_{ct_safe}.html"
                    plot_5layer_ppi_sankey(
                        ppi_links, condition, cell_type,
                        cascade_stats, config, save_path,
                    )
                    print(f"  PPI Sankey: {save_path.name}")
            except Exception as e:
                print(f"  PPI Sankey failed: {e}")

        # 6-layer intercellular Sankey (full cascade)
        try:
            save_path = output_dir / f"sankey_intercell_{condition}_{ct_safe}.html"
            plot_intercell_sankey(
                multicell_obj=multicell,
                results=results,
                cell_type=cell_type,
                seeds=seeds,
                ligand_cells=ligand_source_cells,
                top_ligand_n=50,
                top_receptor_n=30,
                top_tf_n=10,
                before_top_n=5,
                per_celltype=True,
                flow="upstream",
                save_path=str(save_path),
            )
            print(f"  Intercellular Sankey: {save_path.name}")
        except Exception as e:
            print(f"  Intercellular Sankey failed: {e}")


# =============================================================================
# MAIN
# =============================================================================

def main(config: Optional[ReconConfig] = None) -> None:
    """Main execution function."""
    config = get_config(config)

    print("=" * 60)
    print("RECON VISUALIZATION")
    print("Module 7: Figures and Plots")
    print(f"CCC Source: {config.ccc_source}")
    print(f"Conditions: {config.conditions}")
    print(f"Disease conditions: {config.disease_conditions}")
    print(f"Normal condition: {config.normal_condition}")
    print("=" * 60)

    start_time = datetime.now()

    # Setup directories
    figures_dir, data_dir, sankey_dir = setup_directories(config)

    # Load data
    recon_results = load_recon_results(config)
    ccc_data = load_ccc_data(config)

    normal_lower = config.normal_condition.lower()
    normal_combined = recon_results.get(f"{normal_lower}_combined", pd.DataFrame())
    ccc_normal = ccc_data.get(normal_lower, pd.DataFrame())

    # Generate per-condition visualizations
    for cond in config.conditions:
        cond_lower = cond.lower()
        combined = recon_results.get(f"{cond_lower}_combined", pd.DataFrame())
        ccc_cond = ccc_data.get(cond_lower, pd.DataFrame())

        if not combined.empty:
            plot_coordination_heatmap(
                combined, cond,
                figures_dir / f"coordination_heatmap_{cond_lower}.png",
                data_dir,
            )
            plot_top_genes_heatmap(
                combined, cond,
                figures_dir / f"top_genes_heatmap_{cond_lower}.png",
                data_dir,
            )

        # Fibrosis markers only for disease conditions
        if cond_lower != normal_lower and not combined.empty:
            plot_fibrosis_ligand_effects(
                combined, cond,
                figures_dir / f"fibrosis_markers_{cond_lower}.png",
                data_dir,
            )

        # CCC network plot
        if not ccc_cond.empty:
            plot_ccc_network(
                ccc_cond, cond,
                figures_dir / f"ccc_network_{cond_lower}.png",
                data_dir,
            )

    # Generate disease vs normal comparisons
    for disease in config.disease_conditions:
        d_lower = disease.lower()
        disease_combined = recon_results.get(f"{d_lower}_combined", pd.DataFrame())

        if not disease_combined.empty and not normal_combined.empty:
            plot_condition_comparison(
                disease_combined, normal_combined,
                figures_dir / f"condition_comparison_{d_lower}_vs_{normal_lower}.png",
                data_dir,
                condition_name=disease,
            )
            plot_differential_effects(
                disease_combined, normal_combined,
                figures_dir / f"differential_effects_{d_lower}_vs_{normal_lower}.png",
                data_dir,
                condition_name=disease,
            )

        # Summary figure
        ccc_disease = ccc_data.get(d_lower, pd.DataFrame())
        create_summary_figure(
            disease_combined, normal_combined,
            ccc_disease, ccc_normal,
            figures_dir / f"analysis_summary_{d_lower}.png",
            data_dir,
            condition_name=disease,
        )

    # Sankey diagrams (only if seed_categories and focal_celltypes are configured)
    if not config.seed_categories:
        print("\n[ERROR] seed_categories is empty. Sankey diagrams require seed gene categories.")
        print("  Configure seed_categories in config JSON, e.g.:")
        print('    "seed_categories": {"fibrosis": ["TGFB1", "COL1A1"], "inflammatory": ["IL6", "TNF"]}')
        print("  Skipping Sankey generation.")
    if not config.focal_celltypes:
        print("\n[ERROR] focal_celltypes is empty. Sankey diagrams require focal cell types.")
        print("  Configure focal_celltypes in config JSON, e.g.:")
        print('    "focal_celltypes": ["Fibroblast", "Myeloid", "Endothelial"]')
        print("  Skipping Sankey generation.")
    if not config.ligand_source_cells:
        print("\n[WARNING] ligand_source_cells is empty. Using all cell types as ligand sources.")

    if config.seed_categories and config.focal_celltypes:
        print(f"\n{'=' * 60}")
        print("GENERATING SANKEY DIAGRAMS")
        print(f"{'=' * 60}")

        # Save setup.json with seeds and config
        all_seeds = _get_all_seeds(config)
        setup = {
            "timestamp": datetime.now().isoformat(),
            "ccc_source": config.ccc_source,
            "focal_celltypes": config.focal_celltypes,
            "ligand_source_cells": config.ligand_source_cells,
            "min_grn_weight": config.min_grn_weight,
            "seed_categories": config.seed_categories,
            "seed_count_by_category": {
                k: len(v) for k, v in config.seed_categories.items()
            },
            "total_unique_seeds": len(all_seeds),
        }
        setup_path = sankey_dir / "setup.json"
        with open(setup_path, "w") as f:
            json.dump(setup, f, indent=2)
        print(f"Saved: {setup_path}")

        # Save seeds summary
        seeds_summary = pd.DataFrame([
            {"category": cat, "genes": ", ".join(genes), "count": len(genes)}
            for cat, genes in config.seed_categories.items()
        ])
        seeds_summary.to_csv(sankey_dir / "sankey_seeds_summary.csv", index=False)
        print(f"Saved: sankey_seeds_summary.csv")

        # Fetch PPI data once for all conditions
        all_seeds = _get_all_seeds(config)
        all_grn_tfs = set()
        for cond in config.conditions:
            grn = load_grn_for_sankey(cond, config)
            all_grn_tfs.update(grn["source"].str.replace("_TF", "").unique())
        ppi_df = fetch_ppi_partners(all_seeds, sankey_dir, all_grn_tfs, config.ppi_min_score)

        # Load cascade stats once for all conditions
        print(f"\n[{_timestamp()}] Loading cascade differential statistics for Sankey hover tooltips...")
        cascade_indices = load_cascade_stats_index(config)
        for d_lower, idx in cascade_indices.items():
            print(f"  {d_lower} index: {len(idx):,} edges")

        # Generate for each condition
        for condition in config.conditions:
            generate_sankey_for_condition(
                condition.lower(), sankey_dir, config,
                cascade_stats=cascade_indices,
                ppi_df=ppi_df,
            )
    else:
        pass  # Error messages already printed above

    # Summary
    elapsed = datetime.now() - start_time
    print(f"\n{'=' * 60}")
    print("MODULE 7 COMPLETE")
    print(f"{'=' * 60}")
    print(f"[{_timestamp()}] Total time: {elapsed}")
    print(f"Output directory: {figures_dir}")
    if config.seed_categories and config.focal_celltypes:
        print(f"Sankey directory: {sankey_dir}")

    # List generated figures
    print("\nGenerated figures:")
    for fig_path in sorted(figures_dir.glob("*.png")):
        print(f"  {fig_path.name}")

    if config.seed_categories and config.focal_celltypes:
        print("\nGenerated Sankey diagrams:")
        sankey_count = len(list(sankey_dir.glob("*.html")))
        print(f"  {sankey_count} HTML files in {sankey_dir.name}/")


if __name__ == "__main__":
    main()
