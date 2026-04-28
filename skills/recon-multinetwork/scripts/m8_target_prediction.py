#!/usr/bin/env python3
"""
Module 8: Target Prediction for ReCoN Multicellular Network Analysis

Predicts cell-type-specific effects of seed gene perturbation across conditions
using ReCoN's multilayer RWR framework. Generalized from CDH11 treatment pipeline.

Outputs per condition:
- direct/indirect/combined effect tables (ALL genes × cell types)
- Gene rankings with target flag
- Top genes per cell type

Cross-condition outputs:
- Differential tables (ALL genes, log2FC)
- Interactive Plotly scatter plots (ALL genes, targets highlighted)
- Interactive Plotly Sankey plots (4-layer cascades)
- Interactive Plotly heatmaps (target gene log2FC)

Usage:
    # Via pipeline
    python run_pipeline.py --config config.json --start-from 8 --end-at 8

    # Standalone
    python m8_target_prediction.py --config config.json
    python m8_target_prediction.py --config config.json --conditions ssc ipf
    python m8_target_prediction.py --config config.json --sankey-only
    python m8_target_prediction.py --config config.json --skip-plots
"""

import argparse
import gc
import json
import warnings
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import re

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import gseapy as gp
from gseapy.plot import gseaplot
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from recon.data import load_receptor_genes
from recon.explore import multicell_targets, combine_effects

from config import ReconConfig, get_config
from m4_recon_analysis import _apply_multixrank_patch
from m7_visualization import format_condition_name

warnings.filterwarnings("ignore")


# =============================================================================
# TIMESTAMP UTILITY
# =============================================================================

def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


# =============================================================================
# CELL TYPE DISCOVERY
# =============================================================================

def _detect_cell_types(config: ReconConfig) -> List[str]:
    """
    Discover cell types from data_prep_metadata.json or GRN file glob.

    Canonical source is M1 metadata. Falls back to GRN filename parsing.
    """
    # Try metadata first
    meta_path = config.get_output_path("data_prep_metadata.json")
    if meta_path.exists():
        with open(meta_path) as f:
            meta = json.load(f)
        if "cell_types" in meta:
            return sorted(meta["cell_types"])

    # Fallback: discover from GRN file names
    grn_dir = config.get_grn_dir()
    celltypes = set()
    for p in grn_dir.glob("*_rna_network.csv"):
        name_part = p.stem.replace("_rna_network", "")
        for cond in config.conditions:
            suffix = f"_{cond.lower()}"
            if name_part.endswith(suffix):
                ct_raw = name_part[: -len(suffix)]
                celltypes.add(ct_raw)
                break
    return sorted(celltypes)


def _safe_celltype_name(ct: str) -> str:
    """Convert cell type to safe filename component."""
    return ct.lower().replace(" ", "_")


# =============================================================================
# NETWORK LOADING
# =============================================================================

def load_celltype_grns(
    condition: str, config: ReconConfig
) -> Dict[str, pd.DataFrame]:
    """Load cell-type-specific GRNs for a condition."""
    grn_dir = config.get_grn_dir()
    cell_types = _detect_cell_types(config)
    grns = {}

    print(f"\nLoading cell-type GRNs for {condition} (weight > {config.min_grn_weight})...")
    for ct in cell_types:
        ct_safe = _safe_celltype_name(ct)
        grn_path = grn_dir / f"{ct_safe}_{condition}_rna_network.csv"
        if grn_path.exists():
            grn = pd.read_csv(grn_path)
            original_count = len(grn)
            grn = grn[grn["weight"] > config.min_grn_weight]
            print(f"  {ct}: {original_count:,} -> {len(grn):,} edges")
            grns[ct] = grn
        else:
            print(f"  {ct}: No GRN at {grn_path}")
    return grns


def load_ccc(condition: str, config: ReconConfig) -> pd.DataFrame:
    """Load CCC data from configured source."""
    ccc_dir = config.get_ccc_dir()

    if config.ccc_source == "cellchat":
        path = ccc_dir / "cellchat" / f"{condition}_ccc.csv"
        ccc = pd.read_csv(path)
        ccc = ccc.rename(columns={
            "source": "celltype_source",
            "target": "celltype_target",
            "ligand": "source",
            "receptor": "target",
        })
        ccc["lr_means"] = ccc["prob"]
    elif config.ccc_source == "cellphonedb":
        path = ccc_dir / "cellphonedb" / f"{condition}_ccc.csv"
        ccc = pd.read_csv(path)
    else:
        path = ccc_dir / "merged" / f"{condition}_ccc.csv"
        ccc = pd.read_csv(path)

    if "lr_means" in ccc.columns:
        nan_count = ccc["lr_means"].isna().sum()
        if nan_count > 0:
            ccc["lr_means"] = ccc["lr_means"].fillna(0)
    print(f"  CCC ({config.ccc_source}): {len(ccc):,} interactions")
    return ccc


def load_receptor_gene_network() -> pd.DataFrame:
    """Load receptor-gene relationships from NicheNet PKN."""
    print("\nLoading receptor-gene network from NicheNet...")
    receptor_grn = load_receptor_genes("human_receptor_gene_from_NichenetPKN")
    print(f"  {len(receptor_grn):,} receptor-gene relationships")
    return receptor_grn


def merge_celltype_grns(grns: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Merge cell-type GRNs: union of edges, max weight for duplicates."""
    all_edges = []
    for grn in grns.values():
        if grn.empty:
            continue
        all_edges.append(grn[["source", "target", "weight"]].copy())

    if not all_edges:
        return pd.DataFrame(columns=["source", "target", "weight"])

    combined = pd.concat(all_edges, ignore_index=True)
    merged = combined.groupby(["source", "target"], as_index=False).agg({"weight": "max"})
    print(f"  Merged GRN: {len(merged):,} unique edges from {len(grns)} cell types")
    return merged


# =============================================================================
# PER-CONDITION PREDICTION
# =============================================================================

def run_condition_prediction(
    condition: str,
    receptor_grn: pd.DataFrame,
    config: ReconConfig,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict]:
    """Run target prediction for one condition using RWR."""
    seeds = config.load_seeds()
    seed_label = ", ".join(seeds[:5]) + ("..." if len(seeds) > 5 else "")

    print(f"\n{'=' * 60}")
    print(f"TARGET PREDICTION: {condition.upper()}")
    print(f"Seeds: {seed_label}")
    print(f"{'=' * 60}")

    start = datetime.now()

    grns = load_celltype_grns(condition, config)
    ccc = load_ccc(condition, config)
    merged_grn = merge_celltype_grns(grns)

    if merged_grn.empty:
        print("ERROR: No GRN edges!")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), {"error": "empty GRN"}

    # CRITICAL: copy receptor_grn to prevent in-place _receptor suffix mutation
    receptor_grn_copy = receptor_grn.copy()

    # CRITICAL: copy cell types list — ReCoN mutates it into Celltype objects
    all_cell_types = _detect_cell_types(config)
    celltypes_copy = list(all_cell_types)

    print(f"\nRunning multicell_targets(seeds={seeds}, celltypes={len(celltypes_copy)})...")
    direct, indirect = multicell_targets(
        seeds=seeds,
        celltypes=celltypes_copy,
        ccc=ccc,
        grn=merged_grn,
        receptor_grn=receptor_grn_copy,
        restart_proba=config.restart_proba,
        extend_seeds=config.extend_seeds,
        njobs=config.n_jobs,
        verbose=True,
    )

    print(f"  Direct shape: {direct.shape}")
    print(f"  Indirect shape: {indirect.shape}")

    # CRITICAL: copy before combine_effects (it normalizes in-place)
    direct_copy = direct.copy()
    indirect_copy = indirect.copy()

    combined = combine_effects(direct, indirect, alpha=config.alpha)
    print(f"  Combined shape: {combined.shape}")

    elapsed = datetime.now() - start
    stats = {
        "condition": condition,
        "n_grn_edges": int(len(merged_grn)),
        "n_ccc": int(len(ccc)),
        "direct_shape": list(direct_copy.shape),
        "indirect_shape": list(indirect_copy.shape),
        "combined_shape": list(combined.shape),
        "compute_time_seconds": elapsed.total_seconds(),
    }

    return direct_copy, indirect_copy, combined, stats


# =============================================================================
# PER-CONDITION TABLES
# =============================================================================

def make_gene_rankings(
    combined: pd.DataFrame, config: ReconConfig, out_dir: Path
) -> pd.DataFrame:
    """Create gene rankings for ALL genes across focus cell types."""
    focus = config.focus_cell_types if config.focus_cell_types else list(combined.columns)
    focus_cols = [ct for ct in focus if ct in combined.columns]
    df = combined[focus_cols].copy()

    for ct in focus_cols:
        df[f"{ct}_rank"] = df[ct].rank(ascending=False, method="min").astype(int)

    rank_cols = [f"{ct}_rank" for ct in focus_cols]
    df["mean_rank"] = df[rank_cols].mean(axis=1)
    df["is_target"] = df.index.isin(config.target_genes)
    df = df.sort_values("mean_rank")

    df.to_csv(out_dir / "gene_rankings.csv")
    print(f"  gene_rankings.csv: {len(df):,} genes x {len(focus_cols)} cell types")
    return df


def make_top_genes_per_celltype(
    combined: pd.DataFrame, config: ReconConfig, out_dir: Path
) -> pd.DataFrame:
    """Create ALL genes ranked per focus cell type."""
    focus = config.focus_cell_types if config.focus_cell_types else list(combined.columns)
    focus_cols = [ct for ct in focus if ct in combined.columns]
    rows = []
    for ct in focus_cols:
        scores = combined[ct].sort_values(ascending=False)
        for rank, (gene, score) in enumerate(scores.items(), 1):
            rows.append({
                "gene": gene,
                "celltype": ct,
                "rank": rank,
                "score": score,
                "is_target": gene in config.target_genes,
            })

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "top_genes_per_celltype.csv", index=False)
    print(f"  top_genes_per_celltype.csv: {len(df):,} rows")
    return df


# =============================================================================
# CROSS-CONDITION DIFFERENTIAL TABLES
# =============================================================================

def make_differential_tables(
    results: Dict[str, Dict],
    config: ReconConfig,
    diff_dir: Path,
) -> Dict[str, pd.DataFrame]:
    """Create differential tables (disease vs normal) for ALL genes."""
    normal_combined = results[config.normal_condition.lower()]["combined"]
    eps = 1e-10
    diff_tables = {}

    for disease in config.disease_conditions:
        d_lower = disease.lower()
        if d_lower not in results:
            continue
        disease_combined = results[d_lower]["combined"]

        common_genes = sorted(set(disease_combined.index) & set(normal_combined.index))
        common_cts = sorted(set(disease_combined.columns) & set(normal_combined.columns))

        rows = []
        for ct in common_cts:
            d_scores = disease_combined.loc[common_genes, ct]
            n_scores = normal_combined.loc[common_genes, ct]
            log2fc = np.log2((d_scores + eps) / (n_scores + eps))

            for gene in common_genes:
                rows.append({
                    "gene": gene,
                    "celltype": ct,
                    "disease_score": d_scores[gene],
                    "normal_score": n_scores[gene],
                    "log2FC": log2fc[gene],
                    "abs_log2FC": abs(log2fc[gene]),
                    "is_target": gene in config.target_genes,
                })

        df_long = pd.DataFrame(rows)
        fname = f"{d_lower}_vs_normal_all.csv"
        df_long.to_csv(diff_dir / fname, index=False)
        print(f"  {fname}: {len(df_long):,} rows")

        pivot = df_long.pivot_table(
            index="gene", columns="celltype", values="log2FC", aggfunc="first"
        )
        pivot["is_target"] = pivot.index.isin(config.target_genes)
        pivot_fname = f"{d_lower}_vs_normal_pivot.csv"
        pivot.to_csv(diff_dir / pivot_fname)
        print(f"  {pivot_fname}: {pivot.shape}")

        diff_tables[d_lower] = df_long

    return diff_tables


# =============================================================================
# SCATTER PLOTS
# =============================================================================

def make_scatter_pairwise(
    combined: pd.DataFrame, condition: str, config: ReconConfig, out_dir: Path
):
    """Create pairwise scatter plots for all focus cell type pairs."""
    focus = config.focus_cell_types if config.focus_cell_types else list(combined.columns)
    focus_cols = [ct for ct in focus if ct in combined.columns]
    pairs = list(combinations(focus_cols, 2))
    seeds = config.load_seeds()
    seed_label = ", ".join(seeds[:3]) + ("..." if len(seeds) > 3 else "")

    for ct1, ct2 in pairs:
        is_target = combined.index.isin(config.target_genes)
        bg = combined[~is_target]
        tg = combined[is_target]

        threshold1 = combined[ct1].quantile(0.998)
        threshold2 = combined[ct2].quantile(0.998)
        top_genes = combined[
            (combined[ct1] >= threshold1) | (combined[ct2] >= threshold2)
        ]

        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=bg[ct1], y=bg[ct2],
            mode="markers",
            marker=dict(size=3, color="rgba(180,180,180,0.4)"),
            text=bg.index,
            hovertemplate=(
                "<b>%{text}</b><br>"
                f"{ct1}: %{{x:.6f}}<br>"
                f"{ct2}: %{{y:.6f}}<br>"
                "<extra>gene</extra>"
            ),
            name="All genes",
        ))

        fig.add_trace(go.Scatter(
            x=tg[ct1], y=tg[ct2],
            mode="markers+text",
            marker=dict(size=8, color="red", line=dict(width=0.5, color="darkred")),
            text=tg.index,
            textposition="top right",
            textfont=dict(size=8, color="red"),
            hovertemplate=(
                "<b>%{text}</b> (target)<br>"
                f"{ct1}: %{{x:.6f}}<br>"
                f"{ct2}: %{{y:.6f}}<br>"
                "<extra>target</extra>"
            ),
            name="Target genes",
        ))

        top_non_target = top_genes[~top_genes.index.isin(config.target_genes)]
        if len(top_non_target) > 0:
            fig.add_trace(go.Scatter(
                x=top_non_target[ct1], y=top_non_target[ct2],
                mode="markers+text",
                marker=dict(size=5, color="blue", symbol="diamond"),
                text=top_non_target.index,
                textposition="bottom right",
                textfont=dict(size=7, color="blue"),
                hovertemplate=(
                    "<b>%{text}</b> (top quantile)<br>"
                    f"{ct1}: %{{x:.6f}}<br>"
                    f"{ct2}: %{{y:.6f}}<br>"
                    "<extra>top gene</extra>"
                ),
                name="Top 0.2% genes",
            ))

        ct1_safe = _safe_celltype_name(ct1)
        ct2_safe = _safe_celltype_name(ct2)
        cond_display = format_condition_name(condition)

        fig.update_layout(
            title=f"{cond_display}: {ct1} vs {ct2} — {seed_label} Treatment Effect",
            xaxis_title=f"{ct1} combined effect",
            yaxis_title=f"{ct2} combined effect",
            height=700, width=800,
            template="plotly_white",
        )

        fname = f"scatter_{ct1_safe}_vs_{ct2_safe}.html"
        fig.write_html(str(out_dir / fname))
        print(f"    {fname}")


def make_scatter_overview(
    combined: pd.DataFrame, condition: str, config: ReconConfig, out_dir: Path
):
    """Create overview grid of all pairwise comparisons."""
    focus = config.focus_cell_types if config.focus_cell_types else list(combined.columns)
    focus_cols = [ct for ct in focus if ct in combined.columns]
    pairs = list(combinations(focus_cols, 2))
    n_pairs = len(pairs)
    if n_pairs == 0:
        return

    ncols = min(4, n_pairs)
    nrows = (n_pairs + ncols - 1) // ncols

    fig = make_subplots(
        rows=nrows, cols=ncols,
        subplot_titles=[f"{a} vs {b}" for a, b in pairs],
        horizontal_spacing=0.06, vertical_spacing=0.08,
    )

    is_target = combined.index.isin(config.target_genes)

    for idx, (ct1, ct2) in enumerate(pairs):
        row = idx // ncols + 1
        col = idx % ncols + 1

        bg = combined[~is_target]
        tg = combined[is_target]

        fig.add_trace(go.Scatter(
            x=bg[ct1], y=bg[ct2],
            mode="markers",
            marker=dict(size=2, color="rgba(180,180,180,0.3)"),
            text=bg.index,
            hovertemplate="<b>%{text}</b><br>x=%{x:.4f}<br>y=%{y:.4f}<extra></extra>",
            showlegend=False,
        ), row=row, col=col)

        fig.add_trace(go.Scatter(
            x=tg[ct1], y=tg[ct2],
            mode="markers+text",
            marker=dict(size=5, color="red"),
            text=tg.index,
            textposition="top right",
            textfont=dict(size=6, color="red"),
            hovertemplate="<b>%{text}</b> (target)<br>x=%{x:.4f}<br>y=%{y:.4f}<extra></extra>",
            showlegend=(idx == 0),
            name="Target genes",
        ), row=row, col=col)

    seeds = config.load_seeds()
    seed_label = ", ".join(seeds[:3]) + ("..." if len(seeds) > 3 else "")
    cond_display = format_condition_name(condition)

    fig.update_layout(
        title=f"{cond_display}: {seed_label} Treatment Effect — All Cell Type Pairs",
        height=350 * nrows,
        width=350 * ncols,
        template="plotly_white",
        showlegend=True,
    )

    fig.write_html(str(out_dir / "scatter_all_celltypes.html"))
    print(f"    scatter_all_celltypes.html")


# =============================================================================
# SANKEY PLOTS
# =============================================================================

def build_sankey_links(
    condition: str,
    target_cell: str,
    config: ReconConfig,
) -> Dict[str, pd.DataFrame]:
    """
    Build 4-layer Sankey links: Ligand::SrcCell -> Seed::Cell -> TF::Cell -> Gene::Cell
    """
    seeds = config.load_seeds()

    # Layer 1: CCC ligands -> seeds
    ccc = load_ccc(condition, config)
    seed_ccc = ccc[
        (ccc["target"].isin(seeds)) &
        (ccc["celltype_target"] == target_cell)
    ].copy()

    if len(seed_ccc) == 0:
        return {}

    ligand_receptor_rows = []
    for _, row in seed_ccc.iterrows():
        ligand_receptor_rows.append({
            "source": f"{row['source']}::{row['celltype_source']}",
            "target": f"{row['target']}::{target_cell}",
            "weight": row["lr_means"],
        })
    ligand_receptor = pd.DataFrame(ligand_receptor_rows)

    # Layer 2: Seeds -> TFs (NicheNet PKN filtered to TFs in GRN)
    rtf_all = load_receptor_genes("human_receptor_gene_from_NichenetPKN")
    seed_targets = rtf_all[rtf_all["source"].isin(seeds)].copy()

    ct_safe = _safe_celltype_name(target_cell)
    grn_path = config.get_grn_dir() / f"{ct_safe}_{condition}_rna_network.csv"
    if not grn_path.exists():
        return {}

    grn = pd.read_csv(grn_path)
    grn_tfs = set(grn[grn["weight"] >= config.min_sankey_grn_weight]["source"].unique())

    seed_tfs = seed_targets[seed_targets["target"].isin(grn_tfs)].copy()
    seed_tfs = seed_tfs[seed_tfs["weight"] >= config.min_rtf_weight]

    # Post-filter: keep TFs that regulate at least one target gene
    target_genes_set = set(config.target_genes) if config.target_genes else None
    if target_genes_set:
        grn_filtered = grn[
            (grn["source"].isin(seed_tfs["target"])) &
            (grn["target"].isin(target_genes_set)) &
            (grn["weight"] >= config.min_sankey_grn_weight)
        ]
    else:
        grn_filtered = grn[
            (grn["source"].isin(seed_tfs["target"])) &
            (grn["weight"] >= config.min_sankey_grn_weight)
        ]

    connected_tfs = set(grn_filtered["source"].unique())
    seed_tfs = seed_tfs[seed_tfs["target"].isin(connected_tfs)]
    seed_tfs = seed_tfs.sort_values("weight", ascending=False).head(config.top_tfs_sankey)

    if len(seed_tfs) == 0:
        return {}

    receptor_tf_rows = []
    for _, row in seed_tfs.iterrows():
        receptor_tf_rows.append({
            "source": f"{row['source']}::{target_cell}",
            "target": f"{row['target']}::{target_cell}",
            "weight": row["weight"],
        })
    receptor_tf = pd.DataFrame(receptor_tf_rows)

    # Layer 3: TF -> target genes (from GRN)
    tf_list = seed_tfs["target"].tolist()
    if target_genes_set:
        tf_gene_edges = grn[
            (grn["source"].isin(tf_list)) &
            (grn["target"].isin(target_genes_set)) &
            (grn["weight"] >= config.min_sankey_grn_weight)
        ]
    else:
        tf_gene_edges = grn[
            (grn["source"].isin(tf_list)) &
            (grn["weight"] >= config.min_sankey_grn_weight)
        ].head(config.top_grn_genes_sankey * len(tf_list))

    tf_gene = pd.DataFrame({
        "source": tf_gene_edges["source"] + "::" + target_cell,
        "target": tf_gene_edges["target"] + "::" + target_cell,
        "weight": tf_gene_edges["weight"].values,
    })

    return {
        "ligand_receptor": ligand_receptor,
        "receptor_tf": receptor_tf,
        "tf_gene": tf_gene,
    }


def plot_treatment_sankey(
    links: Dict[str, pd.DataFrame],
    condition: str,
    target_cell: str,
    combined: pd.DataFrame,
    normal_combined: pd.DataFrame,
    config: ReconConfig,
    save_path: Path,
):
    """Render Sankey with treatment prediction stats in hover."""
    layer_order = ["ligand_receptor", "receptor_tf", "tf_gene"]
    present_layers = [l for l in layer_order if l in links and len(links[l]) > 0]
    if not present_layers:
        return

    layer_colors = {
        "ligand_receptor": "rgba(160,160,160,0.4)",
        "receptor_tf": "rgba(100,200,100,0.8)",
        "tf_gene": "rgba(100,200,100,0.6)",
    }

    all_link_dfs = []
    for layer_name in present_layers:
        df = links[layer_name].copy()
        df["color"] = layer_colors.get(layer_name, "rgba(200,200,200,0.4)")
        total = df["weight"].sum()
        df["value"] = df["weight"] / total if total > 0 else df["weight"]
        df["layer_type"] = layer_name
        all_link_dfs.append(df)

    all_links = pd.concat(all_link_dfs, ignore_index=True)

    all_nodes = pd.unique(all_links[["source", "target"]].values.ravel())
    node_idx = {name: i for i, name in enumerate(all_nodes)}
    all_links["source_idx"] = all_links["source"].map(node_idx)
    all_links["target_idx"] = all_links["target"].map(node_idx)

    def _fmt(x):
        return x.split("::", 1)[0] if "::" in x else x

    def _cell(x):
        return x.split("::", 1)[1] if "::" in x else "N/A"

    labels = [_fmt(n) for n in all_nodes]
    node_celltypes = [_cell(n) for n in all_nodes]

    node_customdata = []
    for n in all_nodes:
        gene = _fmt(n)
        ct = _cell(n)
        cond_score = float("nan")
        norm_score = float("nan")
        diff_score = float("nan")

        if ct in combined.columns and gene in combined.index:
            cond_score = combined.loc[gene, ct]
        if ct in normal_combined.columns and gene in normal_combined.index:
            norm_score = normal_combined.loc[gene, ct]
        if not np.isnan(cond_score) and not np.isnan(norm_score):
            diff_score = cond_score - norm_score

        node_customdata.append([ct, cond_score, norm_score, diff_score])

    link_customdata = []
    for _, row in all_links.iterrows():
        src_gene = _fmt(row["source"])
        src_ct = _cell(row["source"])
        tgt_gene = _fmt(row["target"])
        tgt_ct = _cell(row["target"])

        src_score = float("nan")
        if src_ct in combined.columns and src_gene in combined.index:
            src_score = combined.loc[src_gene, src_ct]

        tgt_score = float("nan")
        tgt_norm = float("nan")
        tgt_diff = float("nan")
        if tgt_ct in combined.columns and tgt_gene in combined.index:
            tgt_score = combined.loc[tgt_gene, tgt_ct]
        if tgt_ct in normal_combined.columns and tgt_gene in normal_combined.index:
            tgt_norm = normal_combined.loc[tgt_gene, tgt_ct]
        if not np.isnan(tgt_score) and not np.isnan(tgt_norm):
            tgt_diff = tgt_score - tgt_norm

        link_customdata.append([
            src_gene, src_ct, tgt_gene, tgt_ct,
            row["layer_type"], row["weight"],
            src_score, tgt_score, tgt_norm, tgt_diff,
        ])

    seeds = config.load_seeds()
    seed_label = ", ".join(seeds[:3]) + ("..." if len(seeds) > 3 else "")
    cond_display = format_condition_name(condition)

    sankey = go.Sankey(
        node=dict(
            pad=15, thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            customdata=node_customdata,
            hovertemplate=(
                "<b>%{label}</b><br>"
                "Cell: %{customdata[0]}<br>"
                f"{cond_display} score: %{{customdata[1]:.6f}}<br>"
                "Normal score: %{customdata[2]:.6f}<br>"
                "Diff (vs normal): %{customdata[3]:.6f}"
                "<extra></extra>"
            ),
        ),
        link=dict(
            source=all_links["source_idx"],
            target=all_links["target_idx"],
            value=all_links["value"],
            color=all_links["color"],
            customdata=link_customdata,
            hovertemplate=(
                "<b>%{customdata[0]}</b> (%{customdata[1]}) -> "
                "<b>%{customdata[2]}</b> (%{customdata[3]})<br>"
                "Layer: %{customdata[4]}<br>"
                "Edge weight: %{customdata[5]:.4f}<br>"
                "Normalized: %{value:.4f}<br><br>"
                f"<b>{cond_display} scores:</b><br>"
                "Source: %{customdata[6]:.6f}<br>"
                "Target: %{customdata[7]:.6f}<br>"
                "Target normal: %{customdata[8]:.6f}<br>"
                "Target diff: %{customdata[9]:.6f}"
                "<extra></extra>"
            ),
        ),
        orientation="h",
    )

    fig = go.Figure(data=[sankey])
    fig.update_layout(
        title_text=f"{seed_label} Treatment Sankey | {cond_display} | {target_cell}",
        font_size=14, font_color="black",
        height=1000, width=1800,
        margin=dict(l=50, r=150, t=50, b=80),
    )

    layer_labels = ["Ligands", "Receptors", "TFs", "Genes"]
    for i, name in enumerate(layer_labels):
        x = i / (len(layer_labels) - 1) if len(layer_labels) > 1 else 0.5
        fig.add_annotation(
            x=x, y=-0.06,
            text=f"<b>{name}</b>",
            showarrow=False, font=dict(size=15),
            xref="paper", yref="paper",
        )

    fig.write_html(str(save_path))
    print(f"    {save_path.name}")


# =============================================================================
# DIFFERENTIAL HEATMAPS
# =============================================================================

def make_differential_heatmap(
    results: Dict[str, Dict],
    disease: str,
    config: ReconConfig,
    diff_dir: Path,
):
    """Create interactive Plotly heatmap of target gene log2FC."""
    normal_combined = results[config.normal_condition.lower()]["combined"]
    disease_combined = results[disease.lower()]["combined"]
    eps = 1e-10

    focus = config.focus_cell_types if config.focus_cell_types else list(disease_combined.columns)
    focus_cols = [ct for ct in focus
                  if ct in disease_combined.columns and ct in normal_combined.columns]
    target_genes = [g for g in config.target_genes
                    if g in disease_combined.index and g in normal_combined.index]

    if not target_genes or not focus_cols:
        print(f"  Skipping heatmap for {disease}: insufficient data")
        return

    d_vals = disease_combined.loc[target_genes, focus_cols]
    n_vals = normal_combined.loc[target_genes, focus_cols]
    log2fc = np.log2((d_vals + eps) / (n_vals + eps))

    hover_text = []
    for gene in target_genes:
        row = []
        for ct in focus_cols:
            row.append(
                f"Gene: {gene}<br>"
                f"Cell: {ct}<br>"
                f"{format_condition_name(disease)} score: {d_vals.loc[gene, ct]:.6f}<br>"
                f"Normal score: {n_vals.loc[gene, ct]:.6f}<br>"
                f"log2FC: {log2fc.loc[gene, ct]:.4f}"
            )
        hover_text.append(row)

    max_abs = max(abs(log2fc.values.min()), abs(log2fc.values.max()))

    fig = go.Figure(data=go.Heatmap(
        z=log2fc.values,
        x=focus_cols,
        y=target_genes,
        colorscale="RdBu_r",
        zmid=0,
        zmin=-max_abs,
        zmax=max_abs,
        text=hover_text,
        hoverinfo="text",
        colorbar=dict(title="log2FC"),
    ))

    seeds = config.load_seeds()
    seed_label = ", ".join(seeds[:3]) + ("..." if len(seeds) > 3 else "")
    cond_display = format_condition_name(disease)

    fig.update_layout(
        title=f"{seed_label} Target Genes: {cond_display} vs Normal (log2FC)",
        xaxis_title="Cell Type",
        yaxis_title="Gene",
        height=500, width=700,
        template="plotly_white",
    )

    fname = f"heatmap_target_genes_{disease.lower()}_vs_normal.html"
    fig.write_html(str(diff_dir / fname))
    print(f"  {fname}")


# =============================================================================
# GSEA ENRICHMENT
# =============================================================================

_EXCLUDE_GENES = re.compile(
    r"^(AC|AL|AP|AF|AJ|BX)\d{5}"
    r"|^MIR\d"
    r"|^RNU\d"
)


def _detect_gsea_cell_types(df: pd.DataFrame) -> List[str]:
    """Auto-detect cell type score columns from gene_rankings.csv."""
    exclude_patterns = {"gene", "mean_rank"}
    cell_types = []
    for col in df.columns:
        if col in exclude_patterns:
            continue
        if col.endswith("_rank"):
            continue
        if col.startswith("is_"):
            continue
        if pd.api.types.is_numeric_dtype(df[col]):
            cell_types.append(col)
    return cell_types


def _detect_custom_gene_set(df: pd.DataFrame) -> Tuple[str, List[str]]:
    """Auto-detect custom gene set from is_*_target column."""
    target_cols = [c for c in df.columns if c.startswith("is_") and "target" in c.lower()]
    if not target_cols:
        return "", []
    col = target_cols[0]
    # Extract target name: is_cdh11_target -> CDH11, is_target -> target
    name_part = col.replace("is_", "").replace("_target", "").upper()
    if not name_part:
        name_part = "SEED"
    set_name = f"{name_part}_function_target"
    gene_col = "gene" if "gene" in df.columns else df.columns[0]
    target_genes = df.loc[df[col] == True, gene_col].tolist()
    return set_name, target_genes


def _fdr_to_stars(val: float) -> str:
    if pd.isna(val) or val >= 0.05:
        return ""
    if val < 0.0001:
        return "****"
    if val < 0.001:
        return "***"
    if val < 0.01:
        return "**"
    return "*"


def _make_gsea_heatmap(
    gsea_all: pd.DataFrame,
    fdr_threshold: float,
    title: str,
    save_path: Path,
) -> None:
    """Generate clustered NES heatmap with FDR asterisks."""
    nes_pivot = gsea_all.pivot_table(
        index="Term", columns="celltype", values="NES", aggfunc="first"
    )
    fdr_pivot = gsea_all.pivot_table(
        index="Term", columns="celltype", values="FDR q-val", aggfunc="first"
    )

    # Filter to pathways with at least one significant result
    sig_terms = fdr_pivot.index[fdr_pivot.min(axis=1) < fdr_threshold]
    if len(sig_terms) == 0:
        print(f"    No significant pathways for heatmap")
        return
    nes_pivot = nes_pivot.loc[sig_terms]
    fdr_pivot = fdr_pivot.loc[sig_terms]

    nes_filled = nes_pivot.fillna(0)
    annot_matrix = fdr_pivot.map(_fdr_to_stars)

    n_rows = len(nes_filled)
    n_cols = len(nes_filled.columns)
    fig_height = max(8, n_rows * 0.25 + 2)
    fig_width = max(5, n_cols * 1.0 + 3)

    # Dynamic colormap based on NES sign distribution
    has_pos = (nes_filled.values > 0).any()
    has_neg = (nes_filled.values < 0).any()

    if has_pos and has_neg:
        cmap = "bwr"
        max_abs = np.ceil(nes_filled.abs().max().max() * 10) / 10
        vmin, vmax, center = -max_abs, max_abs, 0
    elif has_neg:
        cmap = "YlGnBu"
        vmin, vmax, center = nes_filled.min().min(), nes_filled.max().max(), None
    else:
        cmap = "YlOrRd"
        vmin, vmax, center = nes_filled.min().min(), nes_filled.max().max(), None

    clustermap_kwargs = dict(
        data=nes_filled,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        figsize=(fig_width, fig_height),
        annot=annot_matrix.reindex_like(nes_filled),
        fmt="",
        annot_kws={"size": 8, "color": "black", "weight": "bold"},
        linewidths=0.5,
        linecolor="gray",
        dendrogram_ratio=(0.12, 0.08),
        cbar_kws={"label": "NES", "shrink": 0.6},
        xticklabels=True,
        yticklabels=True,
        method="ward",
        metric="euclidean",
    )
    if center is not None:
        clustermap_kwargs["center"] = center

    g = sns.clustermap(**clustermap_kwargs)
    g.cax.set_position([0.96, 0.3, 0.02, 0.4])
    g.ax_heatmap.set_xlabel("Cell Type", fontsize=12)
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.tick_params(axis="y", labelsize=13)
    g.ax_heatmap.tick_params(axis="x", labelsize=13, rotation=90)
    g.figure.suptitle(title, fontsize=14, fontweight="bold", y=1.02)

    g.savefig(str(save_path), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"    Saved heatmap: {save_path.name} ({n_rows} pathways, cmap={cmap})")


def run_gsea(
    condition: str, config: ReconConfig, out_dir: Path
) -> pd.DataFrame:
    """Run GSEA prerank on gene_rankings.csv for one condition.

    Returns concatenated results DataFrame with NES/FDR per pathway per cell type.
    """
    rankings_csv = out_dir / "gene_rankings.csv"
    if not rankings_csv.exists():
        print(f"    gene_rankings.csv not found in {out_dir}")
        return pd.DataFrame()

    df = pd.read_csv(rankings_csv)
    # Handle index-as-first-column format
    if "gene" not in df.columns and df.columns[0] == "Unnamed: 0":
        df = df.rename(columns={"Unnamed: 0": "gene"})
    elif "gene" not in df.columns:
        df = df.rename(columns={df.columns[0]: "gene"})

    # Filter non-protein-coding genes
    n_before = len(df)
    df = df[~df["gene"].str.match(_EXCLUDE_GENES)].reset_index(drop=True)
    print(f"    Filtered genes: {n_before} -> {len(df)}")

    # Auto-detect cell types and custom gene set
    cell_types = _detect_gsea_cell_types(df)
    set_name, target_genes = _detect_custom_gene_set(df)
    print(f"    Cell types: {cell_types}")
    if target_genes:
        print(f"    Custom gene set: {set_name} ({len(target_genes)} genes)")

    # Load gene set libraries
    combined_gene_sets = {}
    for lib_name in config.gsea_gene_sets:
        lib = gp.get_library(lib_name)
        combined_gene_sets.update(lib)
    if target_genes:
        combined_gene_sets[set_name] = target_genes
    print(f"    Total gene sets: {len(combined_gene_sets)}")

    # Output directory
    gsea_dir = out_dir / "target_gsea"
    gsea_dir.mkdir(exist_ok=True)

    all_results = []
    for ct in cell_types:
        if ct not in df.columns:
            continue
        ct_dir = gsea_dir / ct.replace(" ", "_")
        ct_dir.mkdir(exist_ok=True)

        ranking = (
            df[["gene", ct]]
            .drop_duplicates(subset="gene")
            .sort_values(ct, ascending=False)
            .set_index("gene")[ct]
        )

        pre_res = gp.prerank(
            rnk=ranking,
            gene_sets=combined_gene_sets,
            outdir=str(ct_dir),
            min_size=config.gsea_min_size,
            max_size=config.gsea_max_size,
            permutation_num=config.gsea_permutations,
            seed=42,
            no_plot=True,
            verbose=False,
        )

        res_df = pre_res.res2d.copy()
        res_df["celltype"] = ct
        all_results.append(res_df)
        n_sig = (res_df["FDR q-val"] < config.gsea_fdr_threshold).sum()
        print(f"    {ct}: {n_sig}/{len(res_df)} significant (FDR<{config.gsea_fdr_threshold})")

        # Enrichment plots for significant terms
        viz_metric = pd.Series(np.linspace(1, -1, len(ranking)), index=ranking.index)
        sig_terms = res_df.loc[
            res_df["FDR q-val"] < config.gsea_fdr_threshold, "Term"
        ].tolist()
        for term in sig_terms:
            r = pre_res.results[term]
            safe_term = term.replace("/", "-").replace(":", "_")
            gseaplot(
                term=term,
                hits=r["hits"],
                nes=r["nes"],
                pval=r["pval"],
                fdr=r["fdr"],
                RES=r["RES"],
                rank_metric=viz_metric,
                pheno_pos="High",
                pheno_neg="Low",
                ofname=str(ct_dir / f"{safe_term}.png"),
                figsize=(6, 5.5),
            )

    if not all_results:
        return pd.DataFrame()

    gsea_all = pd.concat(all_results, ignore_index=True)
    gsea_csv = gsea_dir / "gsea_all_celltypes.csv"
    gsea_all.to_csv(gsea_csv, index=False)
    print(f"    Saved: {gsea_csv.name} ({len(gsea_all)} rows)")

    # Per-condition heatmap
    cond_title = format_condition_name(condition)
    _make_gsea_heatmap(
        gsea_all,
        config.gsea_fdr_threshold,
        f"GSEA: {cond_title}\n(NES, {' + '.join(config.gsea_gene_sets)})",
        gsea_dir / "gsea_heatmap.png",
    )
    return gsea_all


def make_cross_condition_heatmap(
    all_gsea_results: Dict[str, pd.DataFrame],
    config: ReconConfig,
    prediction_dir: Path,
) -> None:
    """Generate cross-condition NES comparison heatmap."""
    frames = []
    for cond, df in all_gsea_results.items():
        if df.empty:
            continue
        df_copy = df.copy()
        cond_label = format_condition_name(cond)
        df_copy["celltype_cond"] = df_copy["celltype"] + f" ({cond_label})"
        frames.append(df_copy)

    if not frames:
        return

    combined = pd.concat(frames, ignore_index=True)

    # Replace celltype with celltype_cond for the heatmap
    cross_df = combined.drop(columns=["celltype"]).rename(
        columns={"celltype_cond": "celltype"}
    )

    _make_gsea_heatmap(
        cross_df,
        config.gsea_fdr_threshold,
        "GSEA: Cross-Condition Comparison\n(NES, all conditions)",
        prediction_dir / "gsea_cross_condition_heatmap.png",
    )


# =============================================================================
# UTILITIES
# =============================================================================

def save_sankey_data(
    links: Dict[str, pd.DataFrame],
    condition: str,
    target_cell: str,
    sankey_dir: Path,
) -> None:
    """Save each Sankey layer's DataFrame as CSV."""
    ct_safe = _safe_celltype_name(target_cell)
    data_dir = sankey_dir / "data"
    data_dir.mkdir(exist_ok=True)
    for layer_name, df in links.items():
        fname = f"sankey_{condition}_{ct_safe}_{layer_name}.csv"
        df.to_csv(data_dir / fname, index=False)
    print(f"    Saved {len(links)} layer CSVs to sankey/data/")


# =============================================================================
# MAIN
# =============================================================================

def main(config: Optional[ReconConfig] = None) -> None:
    """Main execution function for target prediction."""

    # Apply multixrank patch for sparse networks
    _apply_multixrank_patch()

    # Handle standalone CLI
    standalone = config is None
    if standalone:
        parser = argparse.ArgumentParser(
            description="M8: Target Prediction",
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        parser.add_argument("--config", required=True, help="Path to JSON config file")
        parser.add_argument("--conditions", nargs="+", type=str, help="Specific conditions to process")
        parser.add_argument("--sankey-only", action="store_true", help="Load existing CSVs, skip RWR")
        parser.add_argument("--skip-plots", action="store_true", help="Skip generating plots")
        parser.add_argument("--output-dir", type=str, help="Override prediction output directory")
        args = parser.parse_args()

        config = ReconConfig.from_json(args.config)
        conditions_override = args.conditions
        sankey_only = args.sankey_only
        skip_plots = args.skip_plots
        if args.output_dir:
            config.prediction_output_dir = args.output_dir
    else:
        conditions_override = None
        sankey_only = False
        skip_plots = False

    seeds = config.load_seeds()
    if not seeds:
        print("=" * 70)
        print("ERROR: No seed genes configured for target prediction.")
        print("  Set seeds in config JSON, e.g.:")
        print('    "seeds": ["CDH11"]')
        print("  Or provide a seeds file:")
        print('    "seeds_file": "path/to/seeds.txt"')
        print("=" * 70)
        raise ValueError("seeds is empty — cannot run target prediction without seed genes")

    seed_label = ", ".join(seeds[:5]) + ("..." if len(seeds) > 5 else "")

    print("=" * 70)
    print("TARGET PREDICTION")
    print(f"Seeds: {seed_label} ({config.seed_type})")
    print(f"CCC source: {config.ccc_source}")
    print(f"Conditions: {config.conditions}")
    print(f"Target genes: {len(config.target_genes)} configured")
    print(f"Focus cell types: {config.focus_cell_types or 'all'}")
    print("=" * 70)

    start_time = datetime.now()

    # Determine output directory
    out_dir = config.get_prediction_dir()
    diff_dir = out_dir / "differential"
    sankey_dir = out_dir / "sankey"
    out_dir.mkdir(parents=True, exist_ok=True)
    diff_dir.mkdir(exist_ok=True)
    sankey_dir.mkdir(exist_ok=True)

    print(f"Output: {out_dir}")

    # Determine which conditions to process
    conditions_to_run = conditions_override if conditions_override else [c.lower() for c in config.conditions]
    for cond in conditions_to_run:
        (out_dir / cond).mkdir(exist_ok=True)

    # Load receptor-gene network (shared)
    receptor_grn = load_receptor_gene_network()

    # Per-condition predictions
    results = {}
    all_stats = {}

    if sankey_only:
        print(f"\n[{_timestamp()}] --sankey-only: Loading existing CSVs...")
        for cond in conditions_to_run:
            cond_dir = out_dir / cond
            combined_path = cond_dir / "combined_effects.csv"
            if combined_path.exists():
                combined = pd.read_csv(combined_path, index_col=0)
                results[cond] = {"combined": combined}
                print(f"  Loaded {cond}: {combined.shape}")
            else:
                print(f"  Warning: {combined_path} not found")
    else:
        for condition in conditions_to_run:
            direct, indirect, combined, stats = run_condition_prediction(
                condition, receptor_grn, config
            )
            results[condition] = {
                "direct": direct,
                "indirect": indirect,
                "combined": combined,
            }
            all_stats[condition] = stats

            cond_dir = out_dir / condition
            if not direct.empty:
                direct.to_csv(cond_dir / "direct_effects.csv")
            if not indirect.empty:
                indirect.to_csv(cond_dir / "indirect_effects.csv")
            if not combined.empty:
                combined.to_csv(cond_dir / "combined_effects.csv")

            print(f"  Saved CSVs to {cond_dir}")
            gc.collect()

    # Auto-populate target_genes if empty (top 20 by mean combined score)
    if not config.target_genes:
        print(f"\n[{_timestamp()}] Auto-selecting top 20 target genes...")
        first_cond = next(
            (c for c in conditions_to_run if c in results and not results[c]["combined"].empty),
            None,
        )
        if first_cond:
            combined = results[first_cond]["combined"]
            mean_scores = combined.mean(axis=1)
            top_20 = mean_scores.nlargest(20).index.tolist()
            config.target_genes = top_20
            print(f"  Selected: {top_20}")

    # Per-condition rankings + top genes
    print(f"\n{'=' * 60}")
    print("PER-CONDITION TABLES")
    print(f"{'=' * 60}")

    for cond in conditions_to_run:
        if cond not in results or results[cond]["combined"].empty:
            continue
        combined = results[cond]["combined"]
        cond_dir = out_dir / cond
        print(f"\n  {cond.upper()}:")
        make_gene_rankings(combined, config, cond_dir)
        make_top_genes_per_celltype(combined, config, cond_dir)

    # GSEA enrichment per condition
    print(f"\n{'=' * 60}")
    print("GSEA ENRICHMENT")
    print(f"{'=' * 60}")

    gsea_results = {}
    for cond in conditions_to_run:
        cond_dir = out_dir / cond
        rankings_csv = cond_dir / "gene_rankings.csv"
        if not rankings_csv.exists():
            continue
        print(f"\n  {cond.upper()}:")
        gsea_df = run_gsea(cond, config, cond_dir)
        if not gsea_df.empty:
            gsea_results[cond] = gsea_df

    # Cross-condition GSEA heatmap
    if len(gsea_results) > 1:
        print(f"\n  Cross-condition heatmap:")
        make_cross_condition_heatmap(gsea_results, config, out_dir)

    # Cross-condition differential tables
    normal_lower = config.normal_condition.lower()
    if normal_lower in results and not results[normal_lower]["combined"].empty:
        print(f"\n{'=' * 60}")
        print("DIFFERENTIAL TABLES")
        print(f"{'=' * 60}")
        make_differential_tables(results, config, diff_dir)

    # Plots
    if not skip_plots:
        # Scatter plots
        print(f"\n{'=' * 60}")
        print("SCATTER PLOTS")
        print(f"{'=' * 60}")

        for cond in conditions_to_run:
            if cond not in results or results[cond]["combined"].empty:
                continue
            combined = results[cond]["combined"]
            cond_dir = out_dir / cond
            print(f"\n  {cond.upper()}:")
            make_scatter_pairwise(combined, cond, config, cond_dir)
            make_scatter_overview(combined, cond, config, cond_dir)

        # Sankey plots
        print(f"\n{'=' * 60}")
        print("SANKEY PLOTS")
        print(f"{'=' * 60}")

        normal_combined = results.get(normal_lower, {}).get("combined", pd.DataFrame())
        focus = config.focus_cell_types if config.focus_cell_types else _detect_cell_types(config)

        for cond in conditions_to_run:
            if cond not in results or results[cond]["combined"].empty:
                continue
            combined = results[cond]["combined"]

            for target_cell in focus:
                ct_safe = _safe_celltype_name(target_cell)
                print(f"\n  {cond.upper()} | {target_cell}:")

                links = build_sankey_links(cond, target_cell, config)
                if not links:
                    print(f"    No Sankey links found")
                    continue

                save_path = sankey_dir / f"sankey_{cond}_{ct_safe}.html"
                plot_treatment_sankey(
                    links, cond, target_cell,
                    combined, normal_combined, config, save_path,
                )
                save_sankey_data(links, cond, target_cell, sankey_dir)

        # Differential heatmaps
        print(f"\n{'=' * 60}")
        print("DIFFERENTIAL HEATMAPS")
        print(f"{'=' * 60}")

        if normal_lower in results and not results[normal_lower]["combined"].empty and config.target_genes:
            for disease in config.disease_conditions:
                d_lower = disease.lower()
                if d_lower in results and not results[d_lower]["combined"].empty:
                    make_differential_heatmap(results, disease, config, diff_dir)

    # Save setup.json
    setup = {
        "timestamp": datetime.now().isoformat(),
        "seeds": config.load_seeds(),
        "seed_type": config.seed_type,
        "conditions": conditions_to_run,
        "all_cell_types": _detect_cell_types(config),
        "focus_cell_types": config.focus_cell_types,
        "target_genes": config.target_genes,
        "restart_proba": config.restart_proba,
        "alpha": config.alpha,
        "min_grn_weight": config.min_grn_weight,
        "min_sankey_grn_weight": config.min_sankey_grn_weight,
        "ccc_source": config.ccc_source,
        "n_jobs": config.n_jobs,
        "per_condition_stats": all_stats,
    }

    with open(out_dir / "setup.json", "w") as f:
        json.dump(setup, f, indent=2, default=str)
    print(f"\nSaved setup.json")

    elapsed = datetime.now() - start_time
    print(f"\n{'=' * 70}")
    print(f"TARGET PREDICTION COMPLETE")
    print(f"Total time: {elapsed}")
    print(f"Output: {out_dir}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
