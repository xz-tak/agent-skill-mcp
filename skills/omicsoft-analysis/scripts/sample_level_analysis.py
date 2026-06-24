#!/usr/bin/env python3
"""
Sample-Level Expression Analysis (Python port)
===============================================
Main entry point for sample-level expression analysis across internal and
external bulk transcriptomic studies. Runs comparison boxplots, GSVA scoring
(via gseapy), and correlation modules.

Usage:
    conda run -n spatial python scripts/sample_level_analysis.py \
        --expr-uri "s3://..." \
        --deg-uri "s3://..." \
        --target-name <target> \
        --targets <genes> \
        --signatures <sigs> \
        --output-dir <output>/sample_level \
        --per-sample-studies <studies> \
        --modules comparison,gsva,corr \
        [--no-cache]
"""

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, ttest_ind
from statsmodels.regression.linear_model import OLS
from statsmodels.stats.multitest import multipletests
from statsmodels.tools import add_constant

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

SCRIPT_DIR = Path(__file__).resolve().parent
SKILL_DIR = SCRIPT_DIR.parent
CACHE_DIR = SKILL_DIR / "cache"
CONFIG_PATH = SKILL_DIR / "references" / "internal_study_configs.json"
EXTRACT_SCRIPT = SCRIPT_DIR / "soma_expr_extract.py"


# =============================================================================
# Argument Parsing
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="Sample-Level Expression Analysis")
    parser.add_argument("--expr-uri", required=True, help="SOMA/h5ad URI for expression data")
    parser.add_argument("--deg-uri", required=True, help="SOMA/h5ad URI for DEG statistics")
    parser.add_argument("--target-name", required=True, help="Target name for labeling")
    parser.add_argument("--targets", required=True, help="Comma-separated target gene symbols")
    parser.add_argument("--signatures", default=None,
                        help="Signatures: SigName:Gene1,Gene2;Sig2:GeneA,GeneB")
    parser.add_argument("--output-dir", required=True, help="Output directory path")
    parser.add_argument("--sources", default="internal,curated,omicsoft",
                        help="Comma-separated source categories")
    parser.add_argument("--per-sample-studies", default=None,
                        help="Comma-separated external study IDs")
    parser.add_argument("--modules", default="comparison,gsva,corr",
                        help="Comma-separated modules to run")
    parser.add_argument("--config-json", default=None, help="Path to internal_study_configs.json")
    parser.add_argument("--extract-script", default=None, help="Path to soma_expr_extract.py")
    parser.add_argument("--backend", default="soma", help="Data backend: soma, h5ad, or auto")
    parser.add_argument("--conda-env", default="spatial", help="Conda environment name")
    parser.add_argument("--no-cache", action="store_true", default=False,
                        help="Force re-extraction (skip cached data)")
    return parser.parse_args()


def parse_signatures(sig_str):
    if not sig_str:
        return {}
    sigs = {}
    for entry in sig_str.split(";"):
        entry = entry.strip()
        if ":" not in entry:
            continue
        parts = entry.split(":", 1)
        sig_name = parts[0].strip()
        genes = [g.strip() for g in parts[1].split(",") if g.strip()]
        if genes:
            sigs[sig_name] = genes
    return sigs


def sanitize_path(x):
    return re.sub(r"[^A-Za-z0-9_]", "_", x)


# =============================================================================
# Data Loading
# =============================================================================

def load_expr_data(uri, genes, studies, extract_script, conda_env, backend="soma"):
    study_str = ",".join(studies) if studies else None
    output_file = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
    try:
        cmd_parts = [
            "conda", "run", "-n", conda_env,
            "python3", str(extract_script),
            "--uri", uri, "--mode", "expr", "--genes", "ALL"
        ]
        if study_str:
            cmd_parts += ["--per-sample-studies", study_str]
        cmd = " ".join(cmd_parts) + f" > {output_file}"
        print(f"[Data] Loading EXPR via: {cmd[:200]}...", file=sys.stderr)
        status = os.system(cmd)
        if status != 0 or not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            print("ERROR: EXPR extraction failed or produced empty output", file=sys.stderr)
            return None
        df = pd.read_csv(output_file, sep="\t")
        print(f"[Data] EXPR loaded: {df.shape[0]} samples x {df.shape[1]} columns", file=sys.stderr)
        return df
    except Exception as e:
        print(f"ERROR loading EXPR data: {e}", file=sys.stderr)
        return None
    finally:
        if os.path.exists(output_file):
            os.unlink(output_file)


def load_deg_data(uri, genes, studies, extract_script, conda_env, backend="soma"):
    gene_str = ",".join(genes)
    study_str = ",".join(studies) if studies else None
    output_file = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
    try:
        cmd_parts = [
            "conda", "run", "-n", conda_env,
            "python3", str(extract_script),
            "--uri", uri, "--mode", "deg", "--genes", gene_str
        ]
        if study_str:
            cmd_parts += ["--per-sample-studies", study_str]
        cmd = " ".join(cmd_parts) + f" > {output_file}"
        print(f"[Data] Loading DEG via: {cmd[:200]}...", file=sys.stderr)
        status = os.system(cmd)
        if status != 0 or not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            print("WARNING: DEG extraction failed", file=sys.stderr)
            return None
        df = pd.read_csv(output_file, sep="\t")
        print(f"[Data] DEG loaded: {df.shape[0]} comparisons x {df.shape[1]} columns", file=sys.stderr)
        return df
    except Exception as e:
        print(f"WARNING: DEG loading failed: {e}", file=sys.stderr)
        return None
    finally:
        if os.path.exists(output_file):
            os.unlink(output_file)


# =============================================================================
# Group Derivation Functions
# =============================================================================

def derive_varsity_groups(df):
    df = df.copy()
    if "metadata_Visit.Type" not in df.columns:
        df["week_response"] = np.nan
        return df
    week_map = {"Visit 1": "Wk0", "Visit 9": "Wk14", "Visit 28": "Wk52"}
    df["varsity_week"] = df["metadata_Visit.Type"].map(week_map)
    df["varsity_treatment"] = df["metadata_Planned.Treatment.for.Period.01"].apply(
        lambda x: "Adalimumab" if pd.notna(x) and "Adalimumab" in str(x)
        else ("Vedolizumab" if pd.notna(x) and "Vedolizumab" in str(x) else np.nan))
    resp14 = df.get("metadata_Clinical.Response.at.Week.14", pd.Series(dtype=str))
    rem52 = df.get("metadata_Clinical.Remission.at.Week.52", pd.Series(dtype=str))
    df["varsity_response"] = np.where(
        (resp14 == "Yes") & (rem52 == "Yes"), "Yes.Yes",
        np.where((resp14 == "No") & (rem52 == "No"), "No.No", "_INVALID_"))
    df["week_response"] = df["varsity_week"].astype(str) + "_" + df["varsity_response"].astype(str)
    valid_levels = {"Wk0_No.No", "Wk0_Yes.Yes", "Wk14_No.No",
                    "Wk14_Yes.Yes", "Wk52_No.No", "Wk52_Yes.Yes"}
    df.loc[~df["week_response"].isin(valid_levels), "week_response"] = np.nan
    return df


def derive_yokohama_rna_fibrosis(df):
    df = df.copy()
    df["yokohama_rna_fibrosis"] = "F" + pd.to_numeric(df["metadata_stage"], errors="coerce").astype("Int64").astype(str)
    valid = {"F0", "F1", "F2", "F3", "F4"}
    df.loc[~df["yokohama_rna_fibrosis"].isin(valid), "yokohama_rna_fibrosis"] = np.nan
    rna_score_map = {"F0": 0, "F1": 1, "F2": 2, "F3": 3, "F4": 4}
    df["fibrosis_score"] = df["yokohama_rna_fibrosis"].map(rna_score_map)
    return df


def derive_yokohama_rna_nash(df):
    df = df.copy()
    nas = pd.to_numeric(df.get("metadata_nas_score"), errors="coerce")
    stage = pd.to_numeric(df.get("metadata_stage"), errors="coerce")
    df["yokohama_rna_nash"] = np.where(
        (nas >= 4) & (stage >= 2), "at_risk",
        np.where(nas.notna(), "control", "_INVALID_"))
    return df


def derive_yokohama_rna_diagnosis(df):
    df = df.copy()
    diag = df.get("metadata_diagnosis", pd.Series(dtype=str))
    df["yokohama_rna_diagnosis"] = np.where(diag.isin(["NAFL", "MASH"]), diag, "_INVALID_")
    return df


def derive_yokohama_prot_fibrosis(df):
    df = df.copy()
    fib = df.get("metadata_Fibrosis", pd.Series(dtype=str)).astype(str)
    valid = {"Healthy", "F0", "F1", "F2", "F3", "F4"}
    df["yokohama_prot_fibrosis"] = np.where(fib.isin(valid), fib, "_INVALID_")
    prot_score_map = {"Healthy": 0, "F0": 0.5, "F1": 1, "F2": 2, "F3": 3, "F4": 4}
    df["fibrosis_score"] = df["yokohama_prot_fibrosis"].map(prot_score_map)
    return df


def derive_yokohama_prot_nash(df):
    df = df.copy()
    nas = pd.to_numeric(df.get("metadata_NAS"), errors="coerce")
    fib = df.get("metadata_Fibrosis", pd.Series(dtype=str)).astype(str)
    at_risk = (nas >= 4) & fib.isin({"F2", "F3", "F4"})
    control = nas.notna() & ~at_risk
    df["yokohama_prot_nash"] = np.where(at_risk, "at_risk",
                                         np.where(control, "control", "_INVALID_"))
    return df


def derive_yokohama_prot_diagnosis(df):
    df = df.copy()
    diag = df.get("metadata_diagnosis", pd.Series(dtype=str))
    df["yokohama_prot_diagnosis"] = np.where(diag.isin(["NAFL", "MASH"]), diag, "_INVALID_")
    return df


def derive_sparc_tissue(df):
    df = df.copy()
    if "metadata_CHARACTERISTICS_BIO_MATERIAL" in df.columns:
        df["sparc_tissue"] = df["metadata_CHARACTERISTICS_BIO_MATERIAL"].apply(
            lambda x: "ileum" if pd.notna(x) and "ileum" in str(x).lower() else "nonileum")
    elif "meta_tissue" in df.columns:
        df["sparc_tissue"] = df["meta_tissue"]
    else:
        df["sparc_tissue"] = "unknown"
    return df


def subset_sparc_disease(df):
    if "metadata_DIAGNOSIS" not in df.columns:
        return {}
    diag = df["metadata_DIAGNOSIS"].astype(str)
    valid_mask = diag.notna() & ~diag.isin(["", "-1", "NA", "\\N", "nan"])
    valid_mask &= ~diag.str.contains("Unclassified", case=False, na=False)
    df_valid = df[valid_mask]
    label_map = {
        "Crohn's Disease": "CD", "crohn's disease (CD)": "CD",
        "Ulcerative Colitis": "UC", "ulcerative colitis (UC)": "UC"
    }
    subsets = {}
    for d in df_valid["metadata_DIAGNOSIS"].unique():
        if d in ["", "-1", "NA", "\\N", "IBD Unclassified"]:
            continue
        short = label_map.get(d, d)
        mask = df_valid["metadata_DIAGNOSIS"] == d
        if mask.sum() > 0:
            subsets[short] = df_valid[mask].copy()
    return subsets


# =============================================================================
# Significance Utilities
# =============================================================================

def _sorted_annotation_patterns(annotation_map):
    """Sort annotation_map patterns by length descending to avoid substring collisions.
    E.g., 'InfFibrvsH' must match before 'FibrvsH' since the latter is a substring."""
    if not annotation_map:
        return []
    return sorted(annotation_map.items(), key=lambda x: len(x[0]), reverse=True)


def get_asterisk(pval):
    if pd.isna(pval):
        return ""
    if pval < 0.0001:
        return "****"
    if pval < 0.001:
        return "***"
    if pval < 0.01:
        return "**"
    if pval < 0.05:
        return "*"
    return ""


def build_sig_annotations(deg_df, expr_long, target_genes, group_levels,
                           control_group, annotation_map, sig_threshold=0.05):
    annotations = []
    if deg_df is None or deg_df.empty or "padj" not in deg_df.columns:
        return pd.DataFrame(columns=["Gene", "group", "y_pos", "label"])

    gene_max = expr_long.groupby("Gene")["Expression"].max().to_dict()

    for _, row in deg_df.iterrows():
        gene = row.get("gene")
        padj_val = row.get("padj")
        if pd.isna(padj_val) or padj_val >= sig_threshold:
            continue
        if gene not in target_genes:
            continue
        asterisk = get_asterisk(padj_val)
        if not asterisk:
            continue

        matched_groups = None
        if annotation_map:
            comp_id = row.get("comparison_id", "")
            for pattern, target_grps in _sorted_annotation_patterns(annotation_map):
                if pattern in str(comp_id):
                    matched_groups = target_grps if isinstance(target_grps, list) else [target_grps]
                    break
        else:
            matched_groups = ["Case"]

        if not matched_groups:
            continue

        max_val = gene_max.get(gene, 0)
        for mg in matched_groups:
            if control_group and mg == control_group:
                continue
            if mg not in group_levels:
                continue
            annotations.append({
                "Gene": gene, "group": mg,
                "y_pos": max_val * 1.002, "label": asterisk
            })

    if not annotations:
        return pd.DataFrame(columns=["Gene", "group", "y_pos", "label"])
    annot_df = pd.DataFrame(annotations)
    annot_df = annot_df.drop_duplicates(subset=["Gene", "group"])
    return annot_df


def build_yokohama_sig_annotations(deg_df, expr_long, target_genes, group_levels,
                                    comparisons, sig_threshold=0.05):
    """Build colored midpoint annotations for Yokohama-style multi-comparison studies.
    ONE asterisk per comparison at the midpoint of its numerator stages."""
    annotations = []
    if deg_df is None or deg_df.empty or not comparisons:
        return pd.DataFrame(columns=["Gene", "x_pos", "y_pos", "label", "comp_label", "comp_color"])
    if "padj" not in deg_df.columns:
        return pd.DataFrame(columns=["Gene", "x_pos", "y_pos", "label", "comp_label", "comp_color"])

    gene_max = expr_long.groupby("Gene")["Expression"].max().to_dict()

    n_groups = len(group_levels)
    dodge_width = 0.75
    group_offsets = {g: dodge_width * (i - (n_groups - 1) / 2) / n_groups
                    for i, g in enumerate(group_levels)}

    annotation_count = {}
    n_total_comps = len(comparisons)
    stack_gap = min(0.025, 0.05 / max(n_total_comps, 1))

    for comp in comparisons:
        comp_id = comp.get("comparison_id", "")
        comp_label = comp.get("label", comp_id)
        comp_color = comp.get("color", "black")
        numerator_stages = comp.get("numerator_stages", [])

        comp_deg = deg_df[deg_df["comparison_id"].astype(str).str.contains(comp_id, na=False)]
        if comp_deg.empty:
            continue

        for gene in target_genes:
            gene_deg = comp_deg[comp_deg["gene"] == gene]
            if gene_deg.empty:
                continue
            padj_val = gene_deg["padj"].iloc[0]
            if pd.isna(padj_val) or padj_val >= sig_threshold:
                continue
            asterisk = get_asterisk(padj_val)
            if not asterisk:
                continue

            gene_idx = target_genes.index(gene)
            valid_stages = [s for s in numerator_stages if s in group_levels]
            if not valid_stages:
                continue

            stage_offsets = [group_offsets[s] for s in valid_stages]
            x_pos = gene_idx + np.mean(stage_offsets)

            gene_max_val = gene_max.get(gene, 0)
            if gene_max_val == 0:
                continue
            existing_count = annotation_count.get(gene, 0)
            y_pos = gene_max_val * 1.002 + existing_count * gene_max_val * stack_gap
            annotation_count[gene] = existing_count + 1

            annotations.append({
                "Gene": gene, "x_pos": x_pos, "y_pos": y_pos,
                "label": asterisk, "comp_label": comp_label, "comp_color": comp_color
            })

    if not annotations:
        return pd.DataFrame(columns=["Gene", "x_pos", "y_pos", "label", "comp_label", "comp_color"])
    annot_df = pd.DataFrame(annotations)
    n_cols_grid = 2
    for gene in annot_df["Gene"].unique():
        gene_mask = annot_df["Gene"] == gene
        n = gene_mask.sum()
        if n <= 1:
            continue
        actual_n_cols = min(n, n_cols_grid)
        x_offsets = []
        row_indices = []
        for idx in range(n):
            row = idx // actual_n_cols
            col = idx % actual_n_cols
            items_in_row = min(actual_n_cols, n - row * actual_n_cols)
            x_off = (col - (items_in_row - 1) / 2) * 0.15
            x_offsets.append(x_off)
            row_indices.append(row)
        gene_max_val = annot_df.loc[gene_mask, "y_pos"].iloc[0]
        annot_df.loc[gene_mask, "x_pos"] = annot_df.loc[gene_mask, "x_pos"].values + np.array(x_offsets)
        annot_df.loc[gene_mask, "y_pos"] = gene_max_val + np.array(row_indices) * gene_max_val * 0.03
    return annot_df


def build_yokohama_gsva_annotations(gsva_stats, sig_names, group_levels, comparisons, sig_threshold=0.05):
    """Build colored midpoint annotations for Yokohama-style GSVA plots.
    ONE asterisk per comparison×signature at the midpoint of its numerator stages."""
    annotations = []
    if gsva_stats is None or gsva_stats.empty or not comparisons:
        return pd.DataFrame(columns=["Signature", "x_pos", "label", "comp_label", "comp_color", "stack_idx"])

    n_groups = len(group_levels)
    dodge_width = 0.75
    group_offsets = {g: dodge_width * (i - (n_groups - 1) / 2) / n_groups
                    for i, g in enumerate(group_levels)}

    annotation_count = {}
    n_total_comps = len(comparisons)
    stack_gap = min(0.025, 0.05 / max(n_total_comps, 1))

    for comp in comparisons:
        comp_label = comp.get("label", comp.get("comparison_id", ""))
        comp_color = comp.get("color", "black")
        numerator_stages = comp.get("numerator_stages", [])
        valid_stages = [s for s in numerator_stages if s in group_levels]
        if not valid_stages:
            continue

        for sig in sig_names:
            best_padj = 1.0
            for stage in valid_stages:
                match = gsva_stats[(gsva_stats["signature"] == sig) & (gsva_stats["group"] == stage)]
                if not match.empty:
                    padj = match["adj.P.Val"].iloc[0]
                    if pd.notna(padj) and padj < best_padj:
                        best_padj = padj

            if best_padj >= sig_threshold:
                continue
            asterisk = get_asterisk(best_padj)
            if not asterisk:
                continue

            sig_idx = sig_names.index(sig)
            stage_offsets = [group_offsets[s] for s in valid_stages]
            x_pos = sig_idx + np.mean(stage_offsets)

            existing_count = annotation_count.get(sig, 0)
            annotation_count[sig] = existing_count + 1

            annotations.append({
                "Signature": sig, "x_pos": x_pos,
                "label": asterisk, "comp_label": comp_label,
                "comp_color": comp_color, "stack_idx": existing_count
            })

    if not annotations:
        return pd.DataFrame(columns=["Signature", "x_pos", "label", "comp_label", "comp_color", "stack_idx"])
    annot_df = pd.DataFrame(annotations)
    n_cols_grid = 2
    for sig in annot_df["Signature"].unique():
        sig_mask = annot_df["Signature"] == sig
        n = sig_mask.sum()
        if n <= 1:
            continue
        actual_n_cols = min(n, n_cols_grid)
        x_offsets = []
        row_indices = []
        for idx in range(n):
            row = idx // actual_n_cols
            col = idx % actual_n_cols
            items_in_row = min(actual_n_cols, n - row * actual_n_cols)
            x_off = (col - (items_in_row - 1) / 2) * 0.15
            x_offsets.append(x_off)
            row_indices.append(row)
        annot_df.loc[sig_mask, "x_pos"] = annot_df.loc[sig_mask, "x_pos"].values + np.array(x_offsets)
        annot_df.loc[sig_mask, "stack_idx"] = row_indices
    return annot_df


# =============================================================================
# Expression Boxplots (Module: comparison)
# =============================================================================

def plot_expression_boxplots(expr_long, group_col, group_levels, group_colors,
                              control_group, sig_annot, title, subtitle,
                              target_genes, facet_col=None, output_dir=".",
                              file_suffix="", deg_df=None, annotation_map=None,
                              show_comparison_legend=False, colored_sig_annot=None,
                              sig_threshold=0.05):
    os.makedirs(output_dir, exist_ok=True)
    import matplotlib.lines as mlines

    fig, ax = plt.subplots(figsize=(12, 8))
    expr_plot = expr_long[expr_long[group_col].isin(group_levels)].copy()
    expr_plot["Gene"] = pd.Categorical(expr_plot["Gene"], categories=target_genes, ordered=True)
    expr_plot[group_col] = pd.Categorical(expr_plot[group_col], categories=group_levels, ordered=True)

    color_list = [group_colors.get(g, "#999999") for g in group_levels]

    sns.boxplot(data=expr_plot, x="Gene", y="Expression", hue=group_col,
                hue_order=group_levels, palette=color_list, dodge=True,
                width=0.75, fliersize=0, ax=ax, linewidth=0.8)
    sns.stripplot(data=expr_plot, x="Gene", y="Expression", hue=group_col,
                  hue_order=group_levels, palette=color_list, dodge=True,
                  jitter=0.1, alpha=0.5, size=1.5, ax=ax, legend=False)

    comp_handles = []
    if colored_sig_annot is not None and not colored_sig_annot.empty:
        for _, row in colored_sig_annot.iterrows():
            ax.text(row["x_pos"], row["y_pos"], row["label"],
                    ha="center", va="bottom", fontsize=12, fontweight="bold",
                    color=row["comp_color"])
        seen = set()
        for _, row in colored_sig_annot.iterrows():
            if row["comp_label"] not in seen:
                comp_handles.append(mlines.Line2D([], [], color=row["comp_color"],
                                    marker="*", markersize=14, linestyle="None",
                                    label=f"* {row['comp_label']}"))
                seen.add(row["comp_label"])
    elif sig_annot is not None and not sig_annot.empty:
        n_groups = len(group_levels)
        dodge_width = 0.75
        offsets = {g: dodge_width * (i - (n_groups - 1) / 2) / n_groups
                   for i, g in enumerate(group_levels)}
        for _, row in sig_annot.iterrows():
            gene_idx = target_genes.index(row["Gene"]) if row["Gene"] in target_genes else None
            if gene_idx is None:
                continue
            offset = offsets.get(row["group"], 0)
            ax.text(gene_idx + offset, row["y_pos"], row["label"],
                    ha="center", va="bottom", fontsize=12, fontweight="bold",
                    color="black")

    annot_y_max = None
    if colored_sig_annot is not None and not colored_sig_annot.empty:
        annot_y_max = colored_sig_annot["y_pos"].max()
    elif sig_annot is not None and not sig_annot.empty:
        annot_y_max = sig_annot["y_pos"].max()
    if annot_y_max is not None:
        cur_bottom, cur_top = ax.get_ylim()
        ax.set_ylim(cur_bottom, max(cur_top, annot_y_max * 1.03))

    ax.set_title(title, fontsize=18, fontweight="bold")
    ax.set_xlabel("Gene", fontsize=16)
    ax.set_ylabel("Normalized Expression", fontsize=16)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    handles, labels = ax.get_legend_handles_labels()
    extra_artists = []
    legend1 = ax.legend(handles[:len(group_levels)], group_levels, title="Group",
                        fontsize=10, title_fontsize=11,
                        bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    ax.add_artist(legend1)
    extra_artists.append(legend1)

    if comp_handles:
        legend2 = ax.legend(handles=comp_handles, title="Comparison",
                            fontsize=9, title_fontsize=10,
                            bbox_to_anchor=(1.02, 0.5), loc="upper left", borderaxespad=0)
        ax.add_artist(legend2)
        extra_artists.append(legend2)

    fig.text(0.5, 0.01, subtitle, ha="center", fontsize=10, color="grey")
    plt.tight_layout(rect=[0, 0.03, 0.78, 1])

    base_name = f"target_expression{file_suffix}"
    fig.savefig(os.path.join(output_dir, f"{base_name}.png"), dpi=300,
                bbox_inches="tight", bbox_extra_artists=extra_artists)
    plt.close(fig)
    print(f"  Saved: {base_name}.png", file=sys.stderr)

    _save_expression_html(expr_long, group_col, group_levels, group_colors,
                          target_genes, title, subtitle, output_dir, base_name,
                          deg_df=deg_df, annotation_map=annotation_map,
                          control_group=control_group,
                          show_comparison_legend=show_comparison_legend,
                          colored_sig_annot=colored_sig_annot,
                          sig_threshold=sig_threshold)


def _save_expression_html(expr_long, group_col, group_levels, group_colors,
                           target_genes, title, subtitle, output_dir, base_name,
                           deg_df=None, annotation_map=None, control_group=None,
                           show_comparison_legend=False, colored_sig_annot=None,
                           sig_threshold=0.05):
    expr_plot = expr_long[expr_long[group_col].isin(group_levels)].copy()
    expr_plot["Gene"] = pd.Categorical(expr_plot["Gene"], categories=target_genes, ordered=True)
    expr_plot[group_col] = pd.Categorical(expr_plot[group_col], categories=group_levels, ordered=True)

    fig = px.box(expr_plot, x="Gene", y="Expression", color=group_col,
                 category_orders={"Gene": target_genes, group_col: group_levels},
                 color_discrete_map=group_colors,
                 title=f"<b>{title}</b><br><span style='font-size:11px;color:grey'>{subtitle}</span>",
                 labels={"Expression": "Normalized Expression"})
    fig.update_traces(boxpoints="all", jitter=0.1, pointpos=0, marker_opacity=0.5, marker_size=3)

    if deg_df is not None and not deg_df.empty:
        _enrich_plotly_hover(fig, expr_plot, deg_df, group_col, annotation_map,
                             control_group, show_comparison_legend,
                             colored_sig_annot=colored_sig_annot,
                             sig_threshold=sig_threshold)

    fig.update_layout(
        legend_title_text="Group", hovermode="closest",
        xaxis=dict(title_font_size=16, tickfont_size=14),
        yaxis=dict(title_font_size=16, tickfont_size=14),
        plot_bgcolor="white",
        paper_bgcolor="rgba(0,0,0,0)",
    )
    fig.write_html(os.path.join(output_dir, f"{base_name}.html"))
    print(f"  Saved: {base_name}.html", file=sys.stderr)


def _enrich_plotly_hover(fig, expr_plot, deg_df, group_col, annotation_map,
                          control_group, show_comparison_legend=False,
                          colored_sig_annot=None, sig_threshold=0.05):
    """Add DEG stats to hover and visual significance asterisks to Plotly boxplot."""
    if deg_df is None or deg_df.empty:
        return

    genes_in_plot = list(expr_plot["Gene"].cat.categories) if hasattr(expr_plot["Gene"], "cat") else sorted(expr_plot["Gene"].unique())
    group_levels_in_plot = list(expr_plot[group_col].cat.categories) if hasattr(expr_plot[group_col], "cat") else sorted(expr_plot[group_col].unique())

    deg_lookup = {}
    for _, row in deg_df.iterrows():
        gene = row.get("gene")
        if gene not in genes_in_plot:
            continue
        comp_id = str(row.get("comparison_id", ""))
        stats = {
            "log2fc": row.get("log2fc"),
            "pval": row.get("pval"),
            "padj": row.get("padj"),
            "comparison_contrast": row.get("comparison_contrast", comp_id),
        }
        if annotation_map:
            for pattern, target_grps in _sorted_annotation_patterns(annotation_map):
                if pattern in comp_id:
                    grps = target_grps if isinstance(target_grps, list) else [target_grps]
                    for g in grps:
                        if g in group_levels_in_plot:
                            key = (gene, g)
                            if key not in deg_lookup or (pd.notna(stats["padj"]) and
                                (pd.isna(deg_lookup[key]["padj"]) or stats["padj"] < deg_lookup[key]["padj"])):
                                deg_lookup[key] = stats
                    break
        else:
            for g in group_levels_in_plot:
                if g != control_group:
                    key = (gene, g)
                    if key not in deg_lookup or (pd.notna(stats["padj"]) and
                        (pd.isna(deg_lookup[key]["padj"]) or stats["padj"] < deg_lookup[key]["padj"])):
                        deg_lookup[key] = stats

    for trace in fig.data:
        group_name = trace.name
        if group_name == control_group:
            trace.hovertemplate = (
                "<b>%{x}</b><br>"
                f"Group: {control_group} (Reference)<br>"
                "Expression: %{y:.3f}<extra></extra>"
            )
        else:
            x_vals = trace.x if trace.x is not None else []
            customdata_list = []
            for gene_val in x_vals:
                key = (str(gene_val), group_name)
                if key in deg_lookup:
                    s = deg_lookup[key]
                    l2fc = f"{s['log2fc']:.3f}" if pd.notna(s.get("log2fc")) else "N/A"
                    pv = f"{s['pval']:.2e}" if pd.notna(s.get("pval")) else "N/A"
                    pa = f"{s['padj']:.2e}" if pd.notna(s.get("padj")) else "N/A"
                    customdata_list.append(f"log2fc: {l2fc} | pval: {pv} | padj: {pa}")
                else:
                    customdata_list.append("No DEG data")
            trace.customdata = np.array(customdata_list).reshape(-1, 1)
            trace.hovertemplate = (
                "<b>%{x}</b><br>"
                f"Group: {group_name}<br>"
                "Expression: %{y:.3f}<br>"
                "%{customdata[0]}<extra></extra>"
            )

    # Annotations: use colored_sig_annot (midpoint) if provided, else per-group
    if colored_sig_annot is not None and not colored_sig_annot.empty:
        for _, row in colored_sig_annot.iterrows():
            fig.add_annotation(
                x=row["x_pos"], y=row["y_pos"],
                text=f"<b>{row['label']}</b>",
                showarrow=False,
                font=dict(size=18, color=row["comp_color"]),
                xanchor="center", yanchor="bottom",
                xref="x", yref="y",
            )
        seen = set()
        for _, row in colored_sig_annot.iterrows():
            if row["comp_label"] not in seen:
                fig.add_trace(go.Scatter(
                    x=[None], y=[None], mode="markers",
                    marker=dict(symbol="star", color=row["comp_color"], size=12),
                    name=f"* {row['comp_label']}",
                    legendgroup="comparisons",
                    legendgrouptitle_text="Comparison",
                    showlegend=True,
                ))
                seen.add(row["comp_label"])
    else:
        gene_x_index = {g: i for i, g in enumerate(genes_in_plot)}
        n_groups = len(group_levels_in_plot)
        group_width = 0.7 / n_groups
        plotly_offsets = {g: (i - (n_groups - 1) / 2) * group_width
                         for i, g in enumerate(group_levels_in_plot)}

        gene_max = expr_plot.groupby("Gene")["Expression"].max().to_dict()
        for (gene, group), stats in deg_lookup.items():
            padj = stats.get("padj")
            if pd.isna(padj) or padj >= sig_threshold:
                continue
            asterisk = get_asterisk(padj)
            if not asterisk or gene not in genes_in_plot:
                continue
            max_val = gene_max.get(gene, 0)
            y_pos = max_val * 1.01
            fig.add_annotation(
                x=gene_x_index[gene] + plotly_offsets.get(group, 0),
                y=y_pos,
                text=f"<b>{asterisk}</b>",
                showarrow=False,
                font=dict(size=18, color="black"),
                xanchor="center", yanchor="bottom",
                xref="x", yref="y",
            )


# =============================================================================
# GSVA Scoring (Module: gsva)
# =============================================================================

def run_gsva_scoring(expr_mat, signatures):
    try:
        import gseapy
    except ImportError:
        print("  WARNING: gseapy not installed. Skipping GSVA.", file=sys.stderr)
        return None

    sigs_filtered = {}
    available_genes = set(expr_mat.index)
    for name, genes in signatures.items():
        found = [g for g in genes if g in available_genes]
        if len(found) >= 3:
            sigs_filtered[name] = found

    if not sigs_filtered:
        print("  WARNING: No signatures have >= 3 genes in data.", file=sys.stderr)
        return None

    for name, genes in sigs_filtered.items():
        print(f"  Signature {name}: {len(genes)} of {len(signatures[name])} genes in data; "
              f"{len(available_genes) - len(genes)} background genes", file=sys.stderr)

    try:
        result = gseapy.gsva(data=expr_mat, gene_sets=sigs_filtered,
                             method="gsva", kcdf="Gaussian",
                             min_size=3, threads=4, verbose=False)
        if hasattr(result, "res2d") and result.res2d is not None:
            df = result.res2d
            if "Term" in df.columns and "Name" in df.columns and "ES" in df.columns:
                return df.pivot(index="Term", columns="Name", values="ES")
            return df
        elif hasattr(result, "pivot_score"):
            return result.pivot_score
        else:
            return result
    except Exception as e:
        print(f"  WARNING: GSVA computation failed: {e}", file=sys.stderr)
        return None


def compare_gsva_scores(gsva_df, sample_groups, control_group, case_groups, group_levels):
    stats_rows = []
    groups_arr = sample_groups.values if hasattr(sample_groups, 'values') else np.array(sample_groups)

    if gsva_df.shape[1] != len(groups_arr):
        print(f"  WARNING: GSVA cols ({gsva_df.shape[1]}) != sample_groups ({len(groups_arr)}). "
              f"Aligning by position.", file=sys.stderr)
        min_len = min(gsva_df.shape[1], len(groups_arr))
        gsva_df = gsva_df.iloc[:, :min_len]
        groups_arr = groups_arr[:min_len]

    for sig_name in gsva_df.index:
        scores = gsva_df.loc[sig_name].values.astype(float)

        if np.all(np.isnan(scores)):
            print(f"  WARNING: GSVA scores all NaN for '{sig_name}'. Skipping.", file=sys.stderr)
            continue

        ctrl_mask = groups_arr == control_group
        ctrl_scores = scores[ctrl_mask]
        ctrl_scores = ctrl_scores[~np.isnan(ctrl_scores)]
        if len(ctrl_scores) < 2:
            continue

        for case_grp in case_groups:
            case_mask = groups_arr == case_grp
            case_scores = scores[case_mask]
            case_scores = case_scores[~np.isnan(case_scores)]
            if len(case_scores) < 2:
                continue
            if np.std(ctrl_scores) == 0 and np.std(case_scores) == 0:
                continue
            logfc = float(case_scores.mean() - ctrl_scores.mean())
            try:
                combined = np.concatenate([ctrl_scores, case_scores])
                binary = np.array([0] * len(ctrl_scores) + [1] * len(case_scores))
                model = OLS(combined, add_constant(binary)).fit()
                pval = float(model.pvalues[1])
            except Exception as e:
                print(f"  WARNING: OLS test failed for '{sig_name}' group '{case_grp}': {e}",
                      file=sys.stderr)
                pval = 1.0

            stats_rows.append({
                "signature": sig_name, "group": case_grp,
                "vs_control": control_group,
                "logFC": logfc, "P.Value": pval,
                "n_group": len(case_scores), "n_control": len(ctrl_scores)
            })

    if not stats_rows:
        print("  WARNING: No valid GSVA comparisons produced.", file=sys.stderr)
        return None
    stats_df = pd.DataFrame(stats_rows)
    _, padj, _, _ = multipletests(stats_df["P.Value"].fillna(1), method="fdr_bh")
    stats_df["adj.P.Val"] = padj
    return stats_df


def plot_gsva_boxplots(gsva_long, group_col, group_levels, group_colors,
                       control_group, gsva_stats, title, subtitle,
                       output_dir=".", file_suffix="", comparisons_list=None,
                       colored_gsva_annot=None):
    os.makedirs(output_dir, exist_ok=True)
    import matplotlib.lines as mlines
    sig_names = gsva_long["Signature"].unique().tolist()

    fig, ax = plt.subplots(figsize=(12, 8))
    color_list = [group_colors.get(g, "#999999") for g in group_levels]

    sns.boxplot(data=gsva_long, x="Signature", y="GSVA_Score", hue=group_col,
                hue_order=group_levels, palette=color_list, dodge=True,
                width=0.75, fliersize=0, ax=ax, linewidth=0.8)
    sns.stripplot(data=gsva_long, x="Signature", y="GSVA_Score", hue=group_col,
                  hue_order=group_levels, palette=color_list, dodge=True,
                  jitter=0.1, alpha=0.5, size=1.5, ax=ax, legend=False)

    sig_max = gsva_long.groupby("Signature")["GSVA_Score"].max().to_dict()

    annot_y_max = None
    if colored_gsva_annot is not None and not colored_gsva_annot.empty:
        n_colored = colored_gsva_annot["stack_idx"].max() + 1
        _gsva_gap = min(0.025, 0.05 / max(n_colored, 1))
        for _, row in colored_gsva_annot.iterrows():
            sig = row["Signature"]
            max_val = sig_max.get(sig, 0)
            y_pos = max_val * 1.002 + row["stack_idx"] * max_val * _gsva_gap
            ax.text(row["x_pos"], y_pos, row["label"],
                    ha="center", va="bottom", fontsize=12, fontweight="bold",
                    color=row["comp_color"])
            if annot_y_max is None or y_pos > annot_y_max:
                annot_y_max = y_pos
    elif gsva_stats is not None and not gsva_stats.empty:
        n_groups = len(group_levels)
        dodge_width = 0.75
        offsets = {g: dodge_width * (i - (n_groups - 1) / 2) / n_groups
                   for i, g in enumerate(group_levels)}
        for _, row in gsva_stats.iterrows():
            asterisk = get_asterisk(row.get("adj.P.Val", 1))
            if not asterisk:
                continue
            sig = row["signature"]
            grp = row["group"]
            if sig not in sig_names or grp not in group_levels:
                continue
            sig_idx = sig_names.index(sig)
            offset = offsets.get(grp, 0)
            y_pos = sig_max.get(sig, 0) * 1.002
            ax.text(sig_idx + offset, y_pos, asterisk,
                    ha="center", va="bottom", fontsize=12, fontweight="bold")
            if annot_y_max is None or y_pos > annot_y_max:
                annot_y_max = y_pos

    if annot_y_max is not None:
        cur_bottom, cur_top = ax.get_ylim()
        ax.set_ylim(cur_bottom, max(cur_top, annot_y_max * 1.03))

    ax.set_title(title, fontsize=18, fontweight="bold")
    ax.set_xlabel("Signature", fontsize=16)
    ax.set_ylabel("GSVA Score", fontsize=16)
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=16)
    handles, labels = ax.get_legend_handles_labels()
    extra_artists = []
    legend1 = ax.legend(handles[:len(group_levels)], group_levels, title="Group",
                        fontsize=10, title_fontsize=11,
                        bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    ax.add_artist(legend1)
    extra_artists.append(legend1)

    comp_handles = []
    if comparisons_list:
        for comp in comparisons_list:
            comp_handles.append(mlines.Line2D([], [], color=comp.get("color", "black"),
                                marker="*", markersize=14, linestyle="None",
                                label=f"* {comp.get('label', comp.get('comparison_id', ''))}"))
    if comp_handles:
        legend2 = ax.legend(handles=comp_handles, title="Comparison",
                            fontsize=9, title_fontsize=10,
                            bbox_to_anchor=(1.02, 0.5), loc="upper left", borderaxespad=0)
        ax.add_artist(legend2)
        extra_artists.append(legend2)

    fig.text(0.5, 0.01, subtitle, ha="center", fontsize=10, color="grey")
    plt.tight_layout(rect=[0, 0.03, 0.78, 1])

    base_name = f"signature_gsva{file_suffix}"
    fig.savefig(os.path.join(output_dir, f"{base_name}.png"), dpi=300,
                bbox_inches="tight", bbox_extra_artists=extra_artists)
    plt.close(fig)
    print(f"  Saved: {base_name}.png", file=sys.stderr)

    # HTML interactive plot with hover stats
    fig_html = px.box(gsva_long, x="Signature", y="GSVA_Score", color=group_col,
                      category_orders={"Signature": sig_names, group_col: group_levels},
                      color_discrete_map=group_colors,
                      title=f"<b>{title}</b><br><span style='font-size:11px;color:grey'>{subtitle}</span>")
    fig_html.update_traces(boxpoints="all", jitter=0.1, pointpos=0,
                           marker_opacity=0.5, marker_size=3)

    # Enrich hover with GSVA stats (logFC, P.Value, adj.P.Val)
    if gsva_stats is not None and not gsva_stats.empty:
        stats_lookup = {}
        for _, row in gsva_stats.iterrows():
            stats_lookup[(row["signature"], row["group"])] = {
                "logFC": row.get("logFC"),
                "P.Value": row.get("P.Value"),
                "adj.P.Val": row.get("adj.P.Val"),
            }
        for trace in fig_html.data:
            group_name = trace.name
            if group_name == control_group:
                trace.hovertemplate = (
                    "<b>%{x}</b><br>"
                    f"Group: {control_group} (Reference)<br>"
                    "GSVA Score: %{y:.3f}<extra></extra>"
                )
            else:
                x_vals = trace.x if trace.x is not None else []
                customdata_list = []
                for sig_val in x_vals:
                    key = (str(sig_val), group_name)
                    if key in stats_lookup:
                        s = stats_lookup[key]
                        lfc = f"{s['logFC']:.3f}" if pd.notna(s.get("logFC")) else "N/A"
                        pv = f"{s['P.Value']:.2e}" if pd.notna(s.get("P.Value")) else "N/A"
                        pa = f"{s['adj.P.Val']:.2e}" if pd.notna(s.get("adj.P.Val")) else "N/A"
                        customdata_list.append(f"logFC: {lfc} | P.Value: {pv} | adj.P.Val: {pa}")
                    else:
                        customdata_list.append("No stats")
                trace.customdata = np.array(customdata_list).reshape(-1, 1)
                trace.hovertemplate = (
                    "<b>%{x}</b><br>"
                    f"Group: {group_name}<br>"
                    "GSVA Score: %{y:.3f}<br>"
                    "%{customdata[0]}<extra></extra>"
                )

    # Add asterisk annotations to HTML plot
    if colored_gsva_annot is not None and not colored_gsva_annot.empty:
        for _, row in colored_gsva_annot.iterrows():
            sig = row["Signature"]
            max_val = sig_max.get(sig, 0)
            y_pos = max_val * 1.002 + row["stack_idx"] * max_val * 0.025
            fig_html.add_annotation(
                x=row["x_pos"], y=y_pos,
                text=f"<b>{row['label']}</b>",
                showarrow=False,
                font=dict(size=18, color=row["comp_color"]),
                xanchor="center", yanchor="bottom",
                xref="x", yref="y",
            )
    elif gsva_stats is not None and not gsva_stats.empty:
        sig_x_index = {s: i for i, s in enumerate(sig_names)}
        n_groups = len(group_levels)
        group_width = 0.7 / n_groups
        plotly_offsets = {g: (i - (n_groups - 1) / 2) * group_width
                         for i, g in enumerate(group_levels)}
        for _, row in gsva_stats.iterrows():
            asterisk = get_asterisk(row.get("adj.P.Val", 1))
            if not asterisk:
                continue
            sig = row["signature"]
            grp = row["group"]
            if sig not in sig_names or grp not in group_levels:
                continue
            max_val = sig_max.get(sig, 0)
            fig_html.add_annotation(
                x=sig_x_index[sig] + plotly_offsets.get(grp, 0),
                y=max_val * 1.01,
                text=f"<b>{asterisk}</b>",
                showarrow=False,
                font=dict(size=18, color="black"),
                xanchor="center", yanchor="bottom",
                xref="x", yref="y",
            )

    if comparisons_list:
        for comp in comparisons_list:
            fig_html.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers",
                marker=dict(symbol="star", color=comp.get("color", "black"), size=12),
                name=f"* {comp.get('label', comp.get('comparison_id', ''))}",
                legendgroup="comparisons",
                legendgrouptitle_text="Comparison",
                showlegend=True,
            ))

    fig_html.update_layout(
        legend_title_text="Group", hovermode="closest",
        xaxis=dict(title_font_size=16, tickfont_size=14),
        yaxis=dict(title_font_size=16, tickfont_size=14),
        plot_bgcolor="white",
        paper_bgcolor="rgba(0,0,0,0)",
    )
    fig_html.write_html(os.path.join(output_dir, f"{base_name}.html"))
    print(f"  Saved: {base_name}.html", file=sys.stderr)


# =============================================================================
# Correlation Modules (Module: corr)
# =============================================================================

def compute_pairwise_cor(mat, method="spearman"):
    """Compute pairwise correlation between rows of mat (genes x samples)."""
    mat = mat.apply(pd.to_numeric, errors="coerce")
    genes = list(mat.index)
    n = len(genes)
    if n < 2:
        return None, None

    r_mat = np.full((n, n), np.nan)
    p_mat = np.full((n, n), np.nan)
    np.fill_diagonal(r_mat, 1.0)
    np.fill_diagonal(p_mat, 0.0)

    for i in range(n):
        for j in range(i + 1, n):
            x = mat.iloc[i].values
            y = mat.iloc[j].values
            valid = ~np.isnan(x) & ~np.isnan(y)
            if valid.sum() < 5:
                continue
            xv, yv = x[valid], y[valid]
            if np.std(xv) == 0 or np.std(yv) == 0:
                continue
            if method == "spearman":
                r, p = spearmanr(xv, yv)
            else:
                r, p = pearsonr(xv, yv)
            r_mat[i, j] = r_mat[j, i] = r
            p_mat[i, j] = p_mat[j, i] = p

    r_df = pd.DataFrame(r_mat, index=genes, columns=genes)
    p_df = pd.DataFrame(p_mat, index=genes, columns=genes)
    return r_df, p_df


def format_cor_annotation(r_val, p_val):
    if pd.isna(r_val) or pd.isna(p_val):
        return ""
    asterisk = get_asterisk(p_val)
    return f"{r_val:.2f}{asterisk}"


def run_gene_gene_correlation(expr_df, target_genes, output_dir, study_name=""):
    print(f"\n--- Module: Gene-Gene Correlation --- {study_name}", file=sys.stderr)
    found_genes = [g for g in target_genes if g in expr_df.columns]
    if len(found_genes) < 2:
        print("  WARNING: < 2 target genes found. Skipping.", file=sys.stderr)
        return None

    gene_mat = expr_df[found_genes].T.astype(float)
    n_samples = gene_mat.shape[1]
    if n_samples < 5:
        print("  WARNING: < 5 samples. Skipping.", file=sys.stderr)
        return None

    print(f"  Genes: {len(found_genes)} | Samples: {n_samples}", file=sys.stderr)
    r_df, p_df = compute_pairwise_cor(gene_mat, method="spearman")
    if r_df is None:
        return None

    # BH correction on upper triangle p-values
    stats_rows = []
    p_values_for_correction = []
    for i in range(len(found_genes)):
        for j in range(i + 1, len(found_genes)):
            p_values_for_correction.append(p_df.iloc[i, j])
            stats_rows.append({
                "gene1": found_genes[i], "gene2": found_genes[j],
                "spearman_r": r_df.iloc[i, j], "pvalue": p_df.iloc[i, j],
                "n": n_samples
            })

    if not stats_rows:
        return None

    stats_df = pd.DataFrame(stats_rows)
    if len(p_values_for_correction) > 0:
        p_arr = np.array(p_values_for_correction)
        p_arr = np.where(np.isnan(p_arr), 1.0, p_arr)
        _, padj, _, _ = multipletests(p_arr, method="fdr_bh")
        stats_df["padj"] = padj
    else:
        stats_df["padj"] = np.nan

    # Annotation matrix for heatmap
    annot_mat = np.full((len(found_genes), len(found_genes)), "", dtype=object)
    for i in range(len(found_genes)):
        annot_mat[i, i] = "1.00"
        for j in range(i + 1, len(found_genes)):
            ann = format_cor_annotation(r_df.iloc[i, j], p_df.iloc[i, j])
            annot_mat[i, j] = ann
            annot_mat[j, i] = ann

    # Heatmap
    cor_dir = os.path.join(output_dir, "correlations", "gene_gene")
    os.makedirs(cor_dir, exist_ok=True)

    n_genes = len(found_genes)
    plot_size = max(8, n_genes * 1.8)
    fig, ax = plt.subplots(figsize=(plot_size, max(6, n_genes * 1.8)))
    sns.heatmap(r_df.values.astype(float), annot=annot_mat, fmt="",
                xticklabels=found_genes, yticklabels=found_genes,
                cmap="bwr", vmin=-1, vmax=1, center=0, ax=ax,
                linewidths=0.3, linecolor="#CCCCCC", annot_kws={"size": 8})
    ax.set_title(f"Gene-Gene Correlation - {study_name}", fontsize=14, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.tight_layout()
    png_path = os.path.join(cor_dir, f"{study_name}_gene_gene_heatmap.png")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {png_path}", file=sys.stderr)

    # Interactive HTML heatmap
    hover_text = []
    for i in range(len(found_genes)):
        row_hover = []
        for j in range(len(found_genes)):
            r_val = r_df.iloc[i, j]
            p_val = p_df.iloc[i, j]
            row_hover.append(
                f"Genes: {found_genes[i]} vs {found_genes[j]}<br>"
                f"rho = {r_val:.3f}<br>p = {p_val:.2e}" if not pd.isna(r_val) else "")
            row_hover.append("")
        hover_text.append(row_hover[:len(found_genes)])

    fig_html = go.Figure(data=go.Heatmap(
        z=r_df.values, x=found_genes, y=found_genes,
        colorscale="RdBu_r", zmin=-1, zmax=1,
        hovertext=hover_text, hoverinfo="text",
        text=annot_mat, texttemplate="%{text}"))
    fig_html.update_layout(
        title=f"Gene-Gene Correlation - {study_name}",
        xaxis_title="", yaxis_title="")
    html_path = os.path.join(cor_dir, f"{study_name}_gene_gene_heatmap.html")
    fig_html.write_html(html_path)
    print(f"  Saved: {html_path}", file=sys.stderr)

    # CSV
    csv_path = os.path.join(cor_dir, f"{study_name}_gene_gene_stats.csv")
    stats_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path} ({len(stats_df)} rows)", file=sys.stderr)
    return stats_df


def run_gsva_gsva_correlation(gsva_scores_df, output_dir, study_name=""):
    """GSVA-GSVA pairwise correlation if >= 2 signatures."""
    if gsva_scores_df is None or gsva_scores_df.shape[0] < 2:
        return None
    print(f"\n--- Module: GSVA-GSVA Correlation --- {study_name}", file=sys.stderr)
    n_samples = gsva_scores_df.shape[1]
    if n_samples < 5:
        print("  WARNING: < 5 samples. Skipping.", file=sys.stderr)
        return None

    sig_names = list(gsva_scores_df.index)
    r_df, p_df = compute_pairwise_cor(gsva_scores_df, method="spearman")
    if r_df is None:
        return None

    stats_rows = []
    for i in range(len(sig_names)):
        for j in range(i + 1, len(sig_names)):
            stats_rows.append({
                "sig1": sig_names[i], "sig2": sig_names[j],
                "spearman_r": r_df.iloc[i, j], "pvalue": p_df.iloc[i, j],
                "n": n_samples
            })
    if not stats_rows:
        return None
    stats_df = pd.DataFrame(stats_rows)

    cor_dir = os.path.join(output_dir, "correlations", "gsva_gsva")
    os.makedirs(cor_dir, exist_ok=True)

    # Heatmap
    n_sigs = len(sig_names)
    annot_mat = np.full((n_sigs, n_sigs), "", dtype=object)
    for i in range(n_sigs):
        annot_mat[i, i] = "1.00"
        for j in range(i + 1, n_sigs):
            ann = format_cor_annotation(r_df.iloc[i, j], p_df.iloc[i, j])
            annot_mat[i, j] = ann
            annot_mat[j, i] = ann

    fig, ax = plt.subplots(figsize=(max(8, n_sigs * 2), max(6, n_sigs * 2)))
    sns.heatmap(r_df.values.astype(float), annot=annot_mat, fmt="",
                xticklabels=sig_names, yticklabels=sig_names,
                cmap="bwr", vmin=-1, vmax=1, center=0, ax=ax,
                linewidths=0.3, linecolor="#CCCCCC")
    ax.set_title(f"GSVA Signature Correlation - {study_name}", fontsize=14, fontweight="bold")
    plt.tight_layout()
    png_path = os.path.join(cor_dir, f"{study_name}_gsva_gsva_heatmap.png")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {png_path}", file=sys.stderr)

    fig_html = go.Figure(data=go.Heatmap(
        z=r_df.values, x=sig_names, y=sig_names,
        colorscale="RdBu_r", zmin=-1, zmax=1,
        text=annot_mat, texttemplate="%{text}"))
    fig_html.update_layout(title=f"GSVA Signature Correlation - {study_name}")
    html_path = os.path.join(cor_dir, f"{study_name}_gsva_gsva_heatmap.html")
    fig_html.write_html(html_path)
    print(f"  Saved: {html_path}", file=sys.stderr)

    csv_path = os.path.join(cor_dir, f"{study_name}_gsva_gsva_stats.csv")
    stats_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}", file=sys.stderr)
    return stats_df


def run_target_vs_continuous(expr_df, target_genes, continuous_cols, output_dir, study_name=""):
    """Target gene expression vs continuous clinical variables."""
    if not continuous_cols:
        return None
    print(f"\n--- Module: Target vs Continuous --- {study_name}", file=sys.stderr)

    found_genes = [g for g in target_genes if g in expr_df.columns]
    if not found_genes:
        print("  WARNING: No target genes found. Skipping.", file=sys.stderr)
        return None

    valid_cols = [c for c in continuous_cols if c in expr_df.columns]
    if not valid_cols:
        print("  WARNING: No valid continuous columns. Skipping.", file=sys.stderr)
        return None

    cor_dir = os.path.join(output_dir, "correlations", "target_vs_clinical")
    os.makedirs(cor_dir, exist_ok=True)

    all_stats = []
    for var_name in valid_cols:
        var_vals = pd.to_numeric(expr_df[var_name], errors="coerce")
        if var_vals.dropna().std() == 0 or var_vals.notna().sum() < 3:
            continue

        scatter_data = []
        annot_data = []
        for gene in found_genes:
            expr_vals = pd.to_numeric(expr_df[gene], errors="coerce")
            valid = expr_vals.notna() & var_vals.notna()
            n_valid = valid.sum()
            if n_valid < 5:
                continue
            xv = var_vals[valid].values
            yv = expr_vals[valid].values
            if np.std(xv) == 0 or np.std(yv) == 0:
                continue
            rho, p_spear = spearmanr(xv, yv)
            r_pear, p_pear = pearsonr(xv, yv)
            all_stats.append({
                "gene": gene, "continuous_var": var_name,
                "spearman_rho": rho, "spearman_pval": p_spear,
                "pearson_r": r_pear, "pearson_pval": p_pear, "n": n_valid
            })
            for xi, yi in zip(xv, yv):
                scatter_data.append({"Gene": gene, "x_val": xi, "y_val": yi})
            annot_data.append({
                "Gene": gene,
                "label": f"Spearman: r={rho:.2f}, p={p_spear:.2e}\nPearson: r={r_pear:.2f}, p={p_pear:.2e}\nn={n_valid}"
            })

        if not scatter_data:
            continue
        scatter_df = pd.DataFrame(scatter_data)
        annot_df = pd.DataFrame(annot_data)

        n_panels = len(found_genes)
        n_cols_p = min(n_panels, 4)
        n_rows_p = int(np.ceil(n_panels / n_cols_p))

        fig, axes = plt.subplots(n_rows_p, n_cols_p,
                                 figsize=(max(10, n_cols_p * 4), max(6, n_rows_p * 3.5)),
                                 squeeze=False)
        from scipy.stats import t as t_dist
        for idx, gene in enumerate(found_genes):
            r_idx, c_idx = divmod(idx, n_cols_p)
            ax = axes[r_idx, c_idx]
            gdf = scatter_df[scatter_df["Gene"] == gene]
            if gdf.empty:
                ax.set_visible(False)
                continue
            ax.scatter(gdf["x_val"], gdf["y_val"], alpha=0.6, s=8, color="#2166AC")
            if len(gdf) >= 3:
                x = gdf["x_val"].values
                y = gdf["y_val"].values
                z = np.polyfit(x, y, 1)
                p_line = np.poly1d(z)
                x_range = np.linspace(x.min(), x.max(), 50)
                y_pred = p_line(x_range)
                ax.plot(x_range, y_pred, color="#B2182B", linewidth=1.2)
                n_pts = len(x)
                residuals = y - p_line(x)
                se = np.sqrt(np.sum(residuals**2) / (n_pts - 2))
                x_mean = x.mean()
                ss_x = np.sum((x - x_mean)**2)
                se_line = se * np.sqrt(1/n_pts + (x_range - x_mean)**2 / ss_x)
                t_val = t_dist.ppf(0.975, n_pts - 2)
                ax.fill_between(x_range, y_pred - t_val * se_line, y_pred + t_val * se_line,
                                alpha=0.15, color="#B2182B")
            ann = annot_df[annot_df["Gene"] == gene]
            if not ann.empty:
                ax.text(0.05, 0.95, ann.iloc[0]["label"], transform=ax.transAxes,
                        fontsize=6, va="top", ha="left", linespacing=1.4)
            ax.set_title(gene, fontsize=11, fontweight="bold")
            ax.set_xlabel(var_name, fontsize=9)
            ax.set_ylabel("Expression", fontsize=9)

        for idx in range(len(found_genes), n_rows_p * n_cols_p):
            r_idx, c_idx = divmod(idx, n_cols_p)
            axes[r_idx, c_idx].set_visible(False)

        fig.suptitle(f"Target Expression vs {var_name} - {study_name}", fontsize=14, fontweight="bold")
        plt.tight_layout()
        base = os.path.join(cor_dir, f"{study_name}_target_vs_{var_name}")
        fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {base}.png", file=sys.stderr)

        # HTML scatter with CI band
        from scipy.stats import t as t_dist_html
        scatter_df["hover_text"] = scatter_df.apply(
            lambda r: f"Gene: {r['Gene']}<br>{var_name}: {r['x_val']:.2f}<br>Expression: {r['y_val']:.3f}",
            axis=1)
        fig_html = px.scatter(scatter_df, x="x_val", y="y_val", facet_col="Gene",
                              facet_col_wrap=n_cols_p, hover_data=["hover_text"],
                              color_discrete_sequence=["#2166AC"],
                              trendline="ols",
                              title=f"Target Expression vs {var_name} - {study_name}",
                              labels={"x_val": var_name, "y_val": "Expression"})
        for i, gene in enumerate(found_genes):
            gdf = scatter_df[scatter_df["Gene"] == gene]
            if len(gdf) >= 3:
                x = gdf["x_val"].values
                y = gdf["y_val"].values
                z = np.polyfit(x, y, 1)
                p_line = np.poly1d(z)
                x_range = np.linspace(x.min(), x.max(), 50)
                y_pred = p_line(x_range)
                n_pts = len(x)
                residuals = y - p_line(x)
                se = np.sqrt(np.sum(residuals**2) / (n_pts - 2))
                x_mean = x.mean()
                ss_x = np.sum((x - x_mean)**2)
                se_line = se * np.sqrt(1/n_pts + (x_range - x_mean)**2 / ss_x)
                t_val = t_dist_html.ppf(0.975, n_pts - 2)
                ci_upper = y_pred + t_val * se_line
                ci_lower = y_pred - t_val * se_line
                axis_suffix = str(i + 1) if i > 0 else ""
                fig_html.add_trace(go.Scatter(
                    x=np.concatenate([x_range, x_range[::-1]]),
                    y=np.concatenate([ci_upper, ci_lower[::-1]]),
                    fill='toself', fillcolor='rgba(178,24,43,0.15)',
                    line=dict(color='rgba(0,0,0,0)'), showlegend=False,
                    hoverinfo='skip', xaxis=f"x{axis_suffix}", yaxis=f"y{axis_suffix}",
                ))
            ann = annot_df[annot_df["Gene"] == gene]
            if not ann.empty:
                fig_html.add_annotation(
                    text=ann.iloc[0]["label"].replace("\n", "<br>"),
                    x=0.05, y=0.95, xanchor="left", yanchor="top",
                    showarrow=False, font=dict(size=9),
                    xref=f"x{i+1} domain" if i > 0 else "x domain",
                    yref=f"y{i+1} domain" if i > 0 else "y domain",
                )
        fig_html.write_html(f"{base}.html")
        print(f"  Saved: {base}.html", file=sys.stderr)

    if all_stats:
        stats_combined = pd.DataFrame(all_stats)
        csv_path = os.path.join(cor_dir, f"{study_name}_target_vs_clinical_stats.csv")
        stats_combined.to_csv(csv_path, index=False)
        print(f"  Saved: {csv_path}", file=sys.stderr)
        return stats_combined
    return None


def run_gsva_vs_continuous(gsva_scores_df, expr_df, continuous_cols, output_dir, study_name=""):
    """GSVA scores vs continuous clinical variables."""
    if not continuous_cols or gsva_scores_df is None or gsva_scores_df.empty:
        return None
    print(f"\n--- Module: GSVA vs Continuous --- {study_name}", file=sys.stderr)

    valid_cols = [c for c in continuous_cols if c in expr_df.columns]
    if not valid_cols:
        return None

    cor_dir = os.path.join(output_dir, "correlations", "gsva_vs_clinical")
    os.makedirs(cor_dir, exist_ok=True)
    sig_names = list(gsva_scores_df.index)
    common_samples = [s for s in gsva_scores_df.columns if s in expr_df.index]
    if len(common_samples) < 5:
        return None

    all_stats = []
    for var_name in valid_cols:
        var_vals = pd.to_numeric(expr_df.loc[common_samples, var_name], errors="coerce")
        if var_vals.dropna().std() == 0 or var_vals.notna().sum() < 3:
            continue
        scatter_data = []
        for sig in sig_names:
            score_vals = gsva_scores_df.loc[sig, common_samples]
            valid = score_vals.notna() & var_vals.notna()
            n_valid = valid.sum()
            if n_valid < 5:
                continue
            xv = var_vals[valid].values
            yv = score_vals[valid].values.astype(float)
            if np.std(xv) == 0 or np.std(yv) == 0:
                continue
            rho, p_spear = spearmanr(xv, yv)
            r_pear, p_pear = pearsonr(xv, yv)
            all_stats.append({
                "signature": sig, "continuous_var": var_name,
                "spearman_rho": rho, "spearman_pval": p_spear,
                "pearson_r": r_pear, "pearson_pval": p_pear, "n": n_valid
            })
            for xi, yi in zip(xv, yv):
                scatter_data.append({"Signature": sig, "x_val": xi, "y_val": yi})

        if not scatter_data:
            continue
        scatter_df = pd.DataFrame(scatter_data)
        from scipy.stats import t as t_dist_gsva
        fig_html = px.scatter(scatter_df, x="x_val", y="y_val", facet_col="Signature",
                              facet_col_wrap=min(len(sig_names), 3),
                              color_discrete_sequence=["#2166AC"], trendline="ols",
                              title=f"GSVA Score vs {var_name} - {study_name}",
                              labels={"x_val": var_name, "y_val": "GSVA Score"})
        for i, sig in enumerate(sig_names):
            sdf = scatter_df[scatter_df["Signature"] == sig]
            if len(sdf) >= 3:
                x = sdf["x_val"].values
                y = sdf["y_val"].values
                z = np.polyfit(x, y, 1)
                p_line = np.poly1d(z)
                x_range = np.linspace(x.min(), x.max(), 50)
                y_pred = p_line(x_range)
                n_pts = len(x)
                residuals = y - p_line(x)
                se = np.sqrt(np.sum(residuals**2) / (n_pts - 2))
                x_mean = x.mean()
                ss_x = np.sum((x - x_mean)**2)
                se_line = se * np.sqrt(1/n_pts + (x_range - x_mean)**2 / ss_x)
                t_val = t_dist_gsva.ppf(0.975, n_pts - 2)
                ci_upper = y_pred + t_val * se_line
                ci_lower = y_pred - t_val * se_line
                axis_suffix = str(i + 1) if i > 0 else ""
                fig_html.add_trace(go.Scatter(
                    x=np.concatenate([x_range, x_range[::-1]]),
                    y=np.concatenate([ci_upper, ci_lower[::-1]]),
                    fill='toself', fillcolor='rgba(178,24,43,0.15)',
                    line=dict(color='rgba(0,0,0,0)'), showlegend=False,
                    hoverinfo='skip', xaxis=f"x{axis_suffix}", yaxis=f"y{axis_suffix}",
                ))
            sig_stat = next((s for s in all_stats if s["signature"] == sig and s["continuous_var"] == var_name), None)
            if sig_stat:
                ann_text = (f"Spearman: r={sig_stat['spearman_rho']:.2f}, p={sig_stat['spearman_pval']:.2e}<br>"
                            f"Pearson: r={sig_stat['pearson_r']:.2f}, p={sig_stat['pearson_pval']:.2e}<br>"
                            f"n={sig_stat['n']}")
                fig_html.add_annotation(
                    text=ann_text, x=0.05, y=0.95, xanchor="left", yanchor="top",
                    showarrow=False, font=dict(size=9),
                    xref=f"x{i+1} domain" if i > 0 else "x domain",
                    yref=f"y{i+1} domain" if i > 0 else "y domain",
                )
        base = os.path.join(cor_dir, f"{study_name}_gsva_vs_{var_name}")
        fig_html.write_html(f"{base}.html")
        print(f"  Saved: {base}.html", file=sys.stderr)

        # PNG scatter
        n_panels = len(sig_names)
        n_cols_p = min(n_panels, 3)
        n_rows_p = int(np.ceil(n_panels / n_cols_p))
        fig, axes = plt.subplots(n_rows_p, n_cols_p,
                                 figsize=(max(10, n_cols_p * 4.5), max(5, n_rows_p * 3.5)),
                                 squeeze=False)
        from scipy.stats import t as t_dist
        for idx, sig in enumerate(sig_names):
            r_idx, c_idx = divmod(idx, n_cols_p)
            ax = axes[r_idx, c_idx]
            sdf = scatter_df[scatter_df["Signature"] == sig]
            if sdf.empty:
                ax.set_visible(False)
                continue
            ax.scatter(sdf["x_val"], sdf["y_val"], alpha=0.6, s=8, color="#2166AC")
            if len(sdf) >= 3:
                x = sdf["x_val"].values
                y = sdf["y_val"].values
                z = np.polyfit(x, y, 1)
                p_line = np.poly1d(z)
                x_range = np.linspace(x.min(), x.max(), 50)
                y_pred = p_line(x_range)
                ax.plot(x_range, y_pred, color="#B2182B", linewidth=1.2)
                n_pts = len(x)
                residuals = y - p_line(x)
                se = np.sqrt(np.sum(residuals**2) / (n_pts - 2))
                x_mean = x.mean()
                ss_x = np.sum((x - x_mean)**2)
                se_line = se * np.sqrt(1/n_pts + (x_range - x_mean)**2 / ss_x)
                t_val = t_dist.ppf(0.975, n_pts - 2)
                ax.fill_between(x_range, y_pred - t_val * se_line, y_pred + t_val * se_line,
                                alpha=0.15, color="#B2182B")
            sig_stat = next((s for s in all_stats if s["signature"] == sig and s["continuous_var"] == var_name), None)
            if sig_stat:
                ann_text = (f"Spearman: r={sig_stat['spearman_rho']:.2f}, p={sig_stat['spearman_pval']:.2e}\n"
                            f"Pearson: r={sig_stat['pearson_r']:.2f}, p={sig_stat['pearson_pval']:.2e}\n"
                            f"n={sig_stat['n']}")
                ax.text(0.05, 0.95, ann_text, transform=ax.transAxes,
                        fontsize=6, va="top", ha="left", linespacing=1.4)
            ax.set_title(sig, fontsize=11, fontweight="bold")
            ax.set_xlabel(var_name, fontsize=9)
            ax.set_ylabel("GSVA Score", fontsize=9)
        for idx in range(len(sig_names), n_rows_p * n_cols_p):
            r_idx, c_idx = divmod(idx, n_cols_p)
            axes[r_idx, c_idx].set_visible(False)
        fig.suptitle(f"GSVA Score vs {var_name} - {study_name}", fontsize=14, fontweight="bold")
        plt.tight_layout()
        fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {base}.png", file=sys.stderr)

    if all_stats:
        stats_combined = pd.DataFrame(all_stats)
        csv_path = os.path.join(cor_dir, f"{study_name}_gsva_vs_clinical_stats.csv")
        stats_combined.to_csv(csv_path, index=False)
        print(f"  Saved: {csv_path}", file=sys.stderr)
        return stats_combined
    return None


# =============================================================================
# Internal Study Handler
# =============================================================================

def apply_group_derivation(df, grp_cfg, study_name):
    if grp_cfg.get("derive_fibrosis_from_stage"):
        return derive_yokohama_rna_fibrosis(df)
    elif grp_cfg.get("derive_nash_nas4_stage2"):
        return derive_yokohama_rna_nash(df)
    elif grp_cfg.get("derive_diagnosis_from_metadata"):
        return derive_yokohama_rna_diagnosis(df)
    elif grp_cfg.get("derive_fibrosis_from_metadata"):
        return derive_yokohama_prot_fibrosis(df)
    elif grp_cfg.get("derive_nash_nas4_fib2"):
        return derive_yokohama_prot_nash(df)
    elif grp_cfg.get("derive_diagnosis_from_prot_metadata"):
        return derive_yokohama_prot_diagnosis(df)
    return df


def run_internal_study(study_name, study_df, expr_mat, deg_df, config,
                       target_genes, signatures, gene_cols, output_dir):
    print(f"\n{'='*40}", file=sys.stderr)
    print(f"Processing Internal Study: {study_name}", file=sys.stderr)
    print(f"{'='*40}", file=sys.stderr)

    study_output = os.path.join(output_dir, sanitize_path(study_name))
    os.makedirs(study_output, exist_ok=True)
    print(f"Study: {study_name} | Samples: {len(study_df)}", file=sys.stderr)

    # Derive Varsity groups if needed
    if config.get("derive_groups") and study_name == "Varsity":
        study_df = derive_varsity_groups(study_df)
        varsity_mask = study_df["week_response"].isin(config["group_levels"]).values
        study_df = study_df.loc[varsity_mask]
        print(f"  Varsity: derived week_response, remaining samples: {len(study_df)}", file=sys.stderr)

    # Derive fibrosis_score for continuous correlation (Yokohama studies)
    if "fibrosis_score" in config.get("continuous_cols", []):
        if study_name == "Yokohama_RNA":
            rna_map = {"F0": 0, "F1": 1, "F2": 2, "F3": 3, "F4": 4}
            stage_col = "metadata_stage"
            if stage_col in study_df.columns:
                fib_label = "F" + pd.to_numeric(study_df[stage_col], errors="coerce").astype("Int64").astype(str)
                study_df["fibrosis_score"] = fib_label.map(rna_map)
        elif study_name == "Yokohama_Protein":
            prot_map = {"Healthy": 0, "F0": 0.5, "F1": 1, "F2": 2, "F3": 3, "F4": 4}
            fib_col = "metadata_Fibrosis"
            if fib_col in study_df.columns:
                study_df["fibrosis_score"] = study_df[fib_col].astype(str).map(prot_map)

    # Determine grouping configs
    if "groupings" in config:
        grouping_list = config["groupings"]
    else:
        grouping_list = {"default": {
            "group_col": config.get("group_col"),
            "group_levels": config.get("group_levels"),
            "control": config.get("control"),
            "group_colors": config.get("group_colors"),
            "annotation_map": config.get("annotation_map"),
            "facet_col": config.get("facet_col"),
            "disease_col": config.get("disease_col"),
        }}

    found_targets = [g for g in target_genes if g in study_df.columns]
    if not found_targets:
        found_targets = [g for g in target_genes if g in gene_cols and g in study_df.columns]
    if not found_targets:
        print(f"  WARNING: No target genes in data. Skipping {study_name}.", file=sys.stderr)
        return

    all_gsva_stats = []

    for grp_name, grp_cfg in grouping_list.items():
        group_col = grp_cfg.get("group_col")
        group_levels = grp_cfg.get("group_levels", [])
        control_group = grp_cfg.get("control")
        group_colors = grp_cfg.get("group_colors", {})
        annotation_map = grp_cfg.get("annotation_map", {})
        facet_col = grp_cfg.get("facet_col")
        disease_col = grp_cfg.get("disease_col")
        grp_suffix = "" if grp_name == "default" else f"_{grp_name}"

        view_df = apply_group_derivation(study_df, grp_cfg, study_name)

        if not group_col or group_col not in view_df.columns:
            print(f"    SKIP grouping {grp_name} - group_col not found: {group_col}", file=sys.stderr)
            continue

        grp_mask = view_df[group_col].isin(group_levels).values
        grp_expr = view_df.loc[grp_mask].copy().reset_index(drop=True)
        if grp_expr.empty:
            continue
        if len(grp_expr) < 4:
            print(f"  SKIP: {study_name} {grp_name} - too few samples", file=sys.stderr)
            continue

        print(f"  Grouping: {grp_name} | {' | '.join(group_levels)} | Samples: {len(grp_expr)}",
              file=sys.stderr)

        # Derive SPARC tissue facet
        if facet_col == "sparc_tissue":
            grp_expr = derive_sparc_tissue(grp_expr)

        # Disease subsets for SPARC
        disease_subsets = [{"label": None, "data": None}]
        if disease_col and disease_col in grp_expr.columns:
            subsets = subset_sparc_disease(grp_expr)
            if len(subsets) > 1:
                disease_subsets = [{"label": nm, "data": sdf} for nm, sdf in subsets.items()]

        for ds in disease_subsets:
            ds_expr = ds["data"] if ds["data"] is not None else grp_expr
            file_suffix = grp_suffix
            title_suffix = "" if grp_name == "default" else f" [{grp_name}]"

            if ds["label"]:
                file_suffix = f"{grp_suffix}_{ds['label']}"
                title_suffix = f"{title_suffix} ({ds['label']})"
                print(f"    Disease subset: {ds['label']} | Samples: {len(ds_expr)}", file=sys.stderr)
                if len(ds_expr) < 4:
                    continue

            available_targets = [g for g in target_genes if g in ds_expr.columns]
            if not available_targets:
                continue

            gene_cols_for_gsva = [c for c in gene_cols if c in ds_expr.columns
                                  and pd.api.types.is_numeric_dtype(ds_expr[c])]

            # Build expression long format
            expr_long = ds_expr[available_targets + [group_col]].melt(
                id_vars=[group_col], var_name="Gene", value_name="Expression")
            expr_long["Gene"] = pd.Categorical(expr_long["Gene"],
                                               categories=available_targets, ordered=True)

            # Filter DEG for disease subset
            ds_deg = deg_df
            if ds["label"] and ds_deg is not None and not ds_deg.empty:
                if "comparison_id" in ds_deg.columns:
                    mask = ds_deg["comparison_id"].str.contains(ds["label"], case=False, na=False)
                    if mask.any():
                        ds_deg = ds_deg[mask]

            # Significance annotations
            sig_threshold_val = config.get("sig_threshold", 0.05)
            show_comp_legend = grp_cfg.get("show_comparison_legend", False)
            comparisons_list = grp_cfg.get("comparisons", [])

            colored_sig_annot = None
            sig_annot = None
            if show_comp_legend and comparisons_list:
                colored_sig_annot = build_yokohama_sig_annotations(
                    ds_deg, expr_long, available_targets, group_levels,
                    comparisons_list, sig_threshold=sig_threshold_val)
            else:
                sig_annot = build_sig_annotations(
                    ds_deg, expr_long, available_targets, group_levels,
                    control_group, annotation_map, sig_threshold=sig_threshold_val)

            sig_col_name = config.get("sig_col", "padj") or "padj"
            sig_threshold_label = f"{sig_col_name}<{sig_threshold_val}"
            ctrl_label = control_group if control_group else "baseline"
            plot_title = f"Target Gene Expression - {study_name}{title_suffix}"
            plot_subtitle = f"* {sig_threshold_label} vs {ctrl_label}"

            # Expression boxplot
            plot_expression_boxplots(
                expr_long, group_col, group_levels, group_colors,
                control_group, sig_annot, plot_title, plot_subtitle,
                available_targets, facet_col=facet_col,
                output_dir=study_output, file_suffix=file_suffix,
                deg_df=ds_deg, annotation_map=annotation_map,
                show_comparison_legend=show_comp_legend,
                colored_sig_annot=colored_sig_annot,
                sig_threshold=sig_threshold_val)

            # GSVA
            gsva_scores_df = None
            gsva_stats = None
            _gsva_sample_index = None
            if signatures and len(gene_cols_for_gsva) >= 5:
                expr_mat_gsva = ds_expr[gene_cols_for_gsva].T.astype(float)
                _npad = len(str(expr_mat_gsva.shape[1] - 1))
                _gsva_sample_index = list(ds_expr.index)
                expr_mat_gsva.columns = [f"s{i:0{_npad}d}" for i in range(expr_mat_gsva.shape[1])]
                expr_mat_gsva = expr_mat_gsva.loc[expr_mat_gsva.abs().sum(axis=1) > 0]
                expr_mat_gsva = expr_mat_gsva[~expr_mat_gsva.index.duplicated()]
                print(f"  GSVA matrix: {expr_mat_gsva.shape[0]} genes x {expr_mat_gsva.shape[1]} samples",
                      file=sys.stderr)
                gsva_scores_df = run_gsva_scoring(expr_mat_gsva, signatures)

                if gsva_scores_df is not None and not gsva_scores_df.empty:
                    nan_count = gsva_scores_df.isna().sum().sum()
                    print(f"  GSVA scores: {gsva_scores_df.shape[0]} sigs x {gsva_scores_df.shape[1]} samples"
                          f", range [{gsva_scores_df.min().min():.4f}, {gsva_scores_df.max().max():.4f}]"
                          f", NaN={nan_count}", file=sys.stderr)

                    if control_group:
                        sample_groups = ds_expr[group_col].values
                        case_groups = [g for g in group_levels if g != control_group]
                        gsva_stats = compare_gsva_scores(
                            gsva_scores_df, sample_groups, control_group, case_groups, group_levels)
                        if gsva_stats is not None:
                            sig_count = (gsva_stats["adj.P.Val"] < 0.05).sum()
                            print(f"  GSVA stats: {sig_count}/{len(gsva_stats)} significant (adj.P.Val < 0.05)",
                                  file=sys.stderr)
                    elif facet_col and facet_col in ds_expr.columns:
                        gsva_stats_parts = []
                        for trt in ds_expr[facet_col].dropna().unique():
                            trt_mask = (ds_expr[facet_col] == trt).values
                            trt_groups = ds_expr.loc[trt_mask, group_col].values
                            trt_gsva = gsva_scores_df.iloc[:, trt_mask]
                            unique_trt_groups = set(trt_groups)
                            no_groups = [g for g in group_levels if "_No.No" in g and g in unique_trt_groups]
                            for no_grp in no_groups:
                                week = no_grp.split("_")[0]
                                yes_grp = f"{week}_Yes.Yes"
                                if yes_grp not in unique_trt_groups:
                                    continue
                                trt_part = compare_gsva_scores(
                                    trt_gsva, trt_groups, yes_grp, [no_grp], [yes_grp, no_grp])
                                if trt_part is not None:
                                    trt_part["treatment"] = trt
                                    gsva_stats_parts.append(trt_part)
                        if gsva_stats_parts:
                            gsva_stats = pd.concat(gsva_stats_parts, ignore_index=True)
                            sig_count = (gsva_stats["adj.P.Val"] < 0.05).sum()
                            print(f"  GSVA stats: {sig_count}/{len(gsva_stats)} significant (adj.P.Val < 0.05)",
                                  file=sys.stderr)

                    # Build GSVA long format
                    gsva_long = gsva_scores_df.reset_index().melt(
                        id_vars=["index"] if "index" in gsva_scores_df.reset_index().columns
                        else [gsva_scores_df.reset_index().columns[0]],
                        var_name="Sample", value_name="GSVA_Score")
                    if "index" in gsva_long.columns:
                        gsva_long = gsva_long.rename(columns={"index": "Signature"})
                    else:
                        gsva_long = gsva_long.rename(
                            columns={gsva_long.columns[0]: "Signature"})
                    # Map sample→group using same zero-padded keys as GSVA columns
                    sample_group_series = ds_expr[group_col].reset_index(drop=True)
                    _npad_map = len(str(len(sample_group_series) - 1))
                    sample_group_map = {f"s{i:0{_npad_map}d}": v for i, v in sample_group_series.items()}
                    gsva_long[group_col] = gsva_long["Sample"].map(sample_group_map)
                    gsva_long = gsva_long.dropna(subset=[group_col])

                    # Build colored GSVA annotations for Yokohama-style comparisons
                    colored_gsva_annot = None
                    if show_comp_legend and comparisons_list and gsva_stats is not None:
                        sig_names_list = gsva_long["Signature"].unique().tolist()
                        colored_gsva_annot = build_yokohama_gsva_annotations(
                            gsva_stats, sig_names_list, group_levels, comparisons_list)

                    gsva_title = f"GSVA Signature Scores - {study_name}{title_suffix}"
                    ctrl_label = control_group if control_group else "all groups"
                    gsva_subtitle = f"* adj.P.Val < 0.05 vs {ctrl_label} (OLS + BH)"
                    plot_gsva_boxplots(
                        gsva_long, group_col, group_levels, group_colors,
                        control_group, gsva_stats, gsva_title, gsva_subtitle,
                        output_dir=study_output, file_suffix=file_suffix,
                        comparisons_list=comparisons_list if show_comp_legend else None,
                        colored_gsva_annot=colored_gsva_annot)

                    if gsva_stats is not None:
                        gsva_stats["view"] = grp_name
                        all_gsva_stats.append(gsva_stats)

    # Save stats CSVs
    if deg_df is not None and not deg_df.empty:
        deg_df.to_csv(os.path.join(study_output, "target_comparison_stats.csv"), index=False)
        print(f"  Saved: target_comparison_stats.csv ({len(deg_df)} rows)", file=sys.stderr)
    if all_gsva_stats:
        combined = pd.concat(all_gsva_stats, ignore_index=True)
        combined.to_csv(os.path.join(study_output, "gsva_comparison_stats.csv"), index=False)
        print(f"  Saved: gsva_comparison_stats.csv", file=sys.stderr)

    # Correlation modules
    print(f"\n  --- Running Correlation Modules ---", file=sys.stderr)
    if len(study_df) >= 5 and len(found_targets) >= 2:
        run_gene_gene_correlation(study_df, found_targets, study_output, study_name)

    if gsva_scores_df is not None and gsva_scores_df.shape[0] >= 2:
        run_gsva_gsva_correlation(gsva_scores_df, study_output, study_name)

    continuous_cols = config.get("continuous_cols", [])
    if continuous_cols and found_targets:
        run_target_vs_continuous(study_df, found_targets, continuous_cols, study_output, study_name)
        if gsva_scores_df is not None and _gsva_sample_index is not None:
            _npad = len(str(len(_gsva_sample_index) - 1))
            col_map = {f"s{i:0{_npad}d}": idx for i, idx in enumerate(_gsva_sample_index)}
            gsva_remapped = gsva_scores_df.rename(columns=col_map)
            run_gsva_vs_continuous(gsva_remapped, study_df, continuous_cols, study_output, study_name)


# =============================================================================
# External Study Handler
# =============================================================================

def assign_comparison_roles(comparison_group_col, comparison_id):
    roles = []
    comp_escaped = re.escape(comparison_id)
    for cg in comparison_group_col:
        if pd.isna(cg) or str(cg) in ("", "\\N"):
            roles.append(np.nan)
        elif re.search(f"{comp_escaped}@case", str(cg)):
            roles.append("Case")
        elif re.search(f"{comp_escaped}@control", str(cg)):
            roles.append("Control")
        else:
            roles.append(np.nan)
    return roles


def run_external_study(study_id, study_df, expr_mat, deg_df,
                       target_genes, signatures, gene_cols, output_dir):
    print(f"\n{'='*40}", file=sys.stderr)
    print(f"Processing External Study: {study_id}", file=sys.stderr)
    print(f"{'='*40}", file=sys.stderr)

    study_output = os.path.join(output_dir, sanitize_path(study_id))
    os.makedirs(study_output, exist_ok=True)

    group_levels = ["Control", "Case"]
    control_group = "Control"
    group_colors = {"Control": "lightblue", "Case": "#F8766D"}
    ext_padj_threshold = 0.05

    if "meta_comparison_group" not in study_df.columns:
        print(f"  SKIP: {study_id} - no meta_comparison_group", file=sys.stderr)
        return
    if deg_df is None or deg_df.empty:
        print(f"  SKIP: {study_id} - no DEG data", file=sys.stderr)
        return

    found_targets = [g for g in target_genes if g in study_df.columns]
    if not found_targets:
        print(f"  SKIP: {study_id} - no target genes in data", file=sys.stderr)
        return

    gene_cols_for_gsva = [c for c in gene_cols if c in study_df.columns
                          and pd.api.types.is_numeric_dtype(study_df[c])]
    comparisons = deg_df["comparison_id"].unique() if "comparison_id" in deg_df.columns else []
    print(f"  Comparisons: {len(comparisons)}", file=sys.stderr)

    all_gsva_stats = []

    for comp_id in comparisons:
        comp_deg = deg_df[deg_df["comparison_id"] == comp_id]
        comp_contrast = comp_deg["comparison_contrast"].iloc[0] if "comparison_contrast" in comp_deg.columns else comp_id
        comp_label = re.sub(r"[^A-Za-z0-9_]", "_", str(comp_contrast))
        comp_label = re.sub(r"_+", "_", comp_label).strip("_")[:80]
        print(f"    Processing comparison: {comp_contrast}", file=sys.stderr)

        group_col = "query_group"
        comp_expr = study_df.copy()
        comp_expr[group_col] = assign_comparison_roles(
            comp_expr["meta_comparison_group"], comp_id)
        comp_expr = comp_expr[comp_expr[group_col].isin(group_levels)]

        if len(comp_expr) < 4:
            print(f"    SKIP comparison: {comp_contrast} - too few samples", file=sys.stderr)
            continue

        n_case = (comp_expr[group_col] == "Case").sum()
        n_ctrl = (comp_expr[group_col] == "Control").sum()
        print(f"    Samples: Case={n_case} Control={n_ctrl}", file=sys.stderr)

        # Expression long
        expr_long = comp_expr[found_targets + [group_col]].melt(
            id_vars=[group_col], var_name="Gene", value_name="Expression")
        expr_long["Gene"] = pd.Categorical(expr_long["Gene"],
                                           categories=found_targets, ordered=True)

        sig_annot = build_sig_annotations(
            comp_deg, expr_long, found_targets, group_levels,
            control_group, None, sig_threshold=ext_padj_threshold)

        plot_title = f"Target Gene Expression - {study_id}"
        plot_subtitle = f"{comp_contrast}\n* padj < {ext_padj_threshold} vs Control"

        plot_expression_boxplots(
            expr_long, group_col, group_levels, group_colors,
            control_group, sig_annot, plot_title, plot_subtitle,
            found_targets, output_dir=study_output,
            file_suffix=f"_{comp_label}", deg_df=comp_deg)

        # GSVA per comparison
        if signatures and len(gene_cols_for_gsva) >= 5:
            expr_mat_gsva = comp_expr[gene_cols_for_gsva].T.astype(float)
            _npad = len(str(expr_mat_gsva.shape[1] - 1))
            expr_mat_gsva.columns = [f"s{i:0{_npad}d}" for i in range(expr_mat_gsva.shape[1])]
            expr_mat_gsva = expr_mat_gsva.loc[expr_mat_gsva.abs().sum(axis=1) > 0]
            expr_mat_gsva = expr_mat_gsva[~expr_mat_gsva.index.duplicated()]
            gsva_scores_df = run_gsva_scoring(expr_mat_gsva, signatures)
            if gsva_scores_df is not None and not gsva_scores_df.empty:
                sample_groups = comp_expr[group_col]
                case_groups = ["Case"]
                gsva_stats = compare_gsva_scores(
                    gsva_scores_df, sample_groups, control_group, case_groups, group_levels)

                gsva_long = gsva_scores_df.reset_index().melt(
                    id_vars=[gsva_scores_df.reset_index().columns[0]],
                    var_name="Sample", value_name="GSVA_Score")
                gsva_long = gsva_long.rename(columns={gsva_long.columns[0]: "Signature"})
                _npad_map = len(str(len(comp_expr) - 1))
                sample_group_map = {f"s{i:0{_npad_map}d}": v for i, v in comp_expr[group_col].reset_index(drop=True).items()}
                gsva_long[group_col] = gsva_long["Sample"].map(sample_group_map)

                gsva_title = f"GSVA Signature Scores - {study_id}"
                gsva_subtitle = f"{comp_contrast}\n* adj.P.Val < 0.05 vs Control (OLS + BH)"
                plot_gsva_boxplots(
                    gsva_long, group_col, group_levels, group_colors,
                    control_group, gsva_stats, gsva_title, gsva_subtitle,
                    output_dir=study_output, file_suffix=f"_{comp_label}")
                if gsva_stats is not None:
                    all_gsva_stats.append(gsva_stats)

    # Stats at study level
    if all_gsva_stats:
        combined = pd.concat(all_gsva_stats, ignore_index=True)
        combined.to_csv(os.path.join(study_output, "gsva_comparison_stats.csv"), index=False)
        print(f"  Saved: gsva_comparison_stats.csv", file=sys.stderr)
    deg_df.to_csv(os.path.join(study_output, "target_comparison_stats.csv"), index=False)
    print(f"  Saved: target_comparison_stats.csv ({len(deg_df)} rows)", file=sys.stderr)

    # Correlation pooled per study
    print(f"\n  --- Running Correlation Modules ---", file=sys.stderr)
    if len(study_df) >= 5 and len(found_targets) >= 2:
        run_gene_gene_correlation(study_df, found_targets, study_output, study_id)


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()

    config_path = args.config_json or str(CONFIG_PATH)
    extract_script = args.extract_script or str(EXTRACT_SCRIPT)

    target_genes = [g.strip() for g in args.targets.split(",") if g.strip()]
    signatures = parse_signatures(args.signatures)
    modules = [m.strip() for m in args.modules.split(",")]
    external_studies = ([s.strip() for s in args.per_sample_studies.split(",") if s.strip()]
                        if args.per_sample_studies else [])

    print("=" * 58, file=sys.stderr)
    print(f"SAMPLE-LEVEL EXPRESSION ANALYSIS: {args.target_name}", file=sys.stderr)
    print("=" * 58, file=sys.stderr)
    print(f"  Targets: {', '.join(target_genes)}", file=sys.stderr)
    print(f"  Signatures: {len(signatures)}", file=sys.stderr)
    print(f"  Modules: {', '.join(modules)}", file=sys.stderr)
    print(f"  Output: {args.output_dir}", file=sys.stderr)

    os.makedirs(args.output_dir, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # Load config
    with open(config_path) as f:
        config = json.load(f)
    internal_studies = list(config.keys())
    print(f"  Internal studies: {', '.join(internal_studies)}", file=sys.stderr)

    all_genes = list(set(target_genes + [g for genes in signatures.values() for g in genes]))
    print(f"  Total unique genes: {len(all_genes)}", file=sys.stderr)

    # Cache keys
    no_cache = args.no_cache
    if no_cache:
        for f_path in CACHE_DIR.glob("*"):
            if f_path.is_file():
                f_path.unlink()
                print(f"[Cache] --no-cache: removed {f_path.name}", file=sys.stderr)

    expr_cache_key = hashlib.md5(
        str(sorted(internal_studies) + ["ALL_GENES"] + [args.backend]).encode()
    ).hexdigest()
    expr_cache_path = CACHE_DIR / f"expr_{expr_cache_key}.tsv"

    deg_cache_key = hashlib.md5(
        str(sorted(internal_studies) + sorted(target_genes) + ["deg"] + [args.backend]).encode()
    ).hexdigest()
    deg_cache_path = CACHE_DIR / f"deg_{deg_cache_key}.tsv"

    # Load internal EXPR
    internal_expr_df = None
    if not no_cache and expr_cache_path.exists():
        print(f"[Cache] Using skill cache: {expr_cache_path}", file=sys.stderr)
        internal_expr_df = pd.read_csv(expr_cache_path, sep="\t", index_col=0)
        if not any(c.startswith("metadata_") for c in internal_expr_df.columns):
            print("[Cache] Stale cache (no metadata_* columns). Re-extracting.", file=sys.stderr)
            internal_expr_df = None

    if internal_expr_df is None:
        internal_expr_df = load_expr_data(
            args.expr_uri, all_genes, internal_studies,
            extract_script, args.conda_env, args.backend)
        if internal_expr_df is not None and not internal_expr_df.empty:
            internal_expr_df.to_csv(expr_cache_path, sep="\t")
            print(f"[Cache] Saved EXPR: {expr_cache_path}", file=sys.stderr)

    # Load internal DEG
    internal_deg_df = None
    if not no_cache and deg_cache_path.exists():
        print(f"[Cache] Using skill DEG cache: {deg_cache_path}", file=sys.stderr)
        internal_deg_df = pd.read_csv(deg_cache_path, sep="\t")

    if internal_deg_df is None:
        internal_deg_df = load_deg_data(
            args.deg_uri, target_genes, internal_studies,
            extract_script, args.conda_env, args.backend)
        if internal_deg_df is not None and not internal_deg_df.empty:
            internal_deg_df.to_csv(deg_cache_path, sep="\t", index=False)
            print(f"[Cache] Saved DEG: {deg_cache_path}", file=sys.stderr)

    # Load external
    external_expr_df = None
    external_deg_df = None
    if external_studies:
        external_expr_df = load_expr_data(
            args.expr_uri, all_genes, external_studies,
            extract_script, args.conda_env, args.backend)
        external_deg_df = load_deg_data(
            args.deg_uri, target_genes, external_studies,
            extract_script, args.conda_env, args.backend)

    # Combine data
    expr_df = internal_expr_df
    if external_expr_df is not None and not external_expr_df.empty:
        if expr_df is not None:
            expr_df = pd.concat([expr_df, external_expr_df], ignore_index=True)
        else:
            expr_df = external_expr_df

    if expr_df is None or expr_df.empty:
        print("FATAL: Could not load expression data. Exiting.", file=sys.stderr)
        sys.exit(1)

    deg_df = internal_deg_df
    if external_deg_df is not None and not external_deg_df.empty:
        if deg_df is not None:
            common_cols = list(set(deg_df.columns) & set(external_deg_df.columns))
            deg_df = pd.concat([deg_df[common_cols], external_deg_df[common_cols]], ignore_index=True)
        else:
            deg_df = external_deg_df

    # Identify gene vs metadata columns
    meta_cols = [c for c in expr_df.columns if c.startswith("meta_") or c.startswith("metadata_")]
    gene_cols = [c for c in expr_df.columns if c not in meta_cols]

    if not gene_cols:
        print("FATAL: No gene expression columns found. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Determine project_id column
    pid_col = None
    for candidate in ["meta_project_id", "project_id", "study"]:
        if candidate in expr_df.columns:
            pid_col = candidate
            break
    if pid_col is None:
        print("FATAL: No project_id/study column found. Exiting.", file=sys.stderr)
        sys.exit(1)

    print(f"[Matrix] Expression: {len(gene_cols)} genes x {len(expr_df)} samples", file=sys.stderr)

    # Manifest
    manifest = {
        "status": "success",
        "timestamp": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
        "studies_processed": [],
        "studies_skipped": [],
        "modules_run": modules,
        "warnings": [],
        "package_versions": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "soma_uris": {"expr": args.expr_uri, "deg": args.deg_uri}
    }

    # Process Internal Studies
    print("\n\n=== INTERNAL STUDIES ===", file=sys.stderr)
    for study_name in internal_studies:
        study_mask = (expr_df[pid_col] == study_name).values
        if study_mask.sum() == 0:
            msg = f"{study_name}: 0 samples in data"
            print(f"  SKIP: {msg}", file=sys.stderr)
            manifest["studies_skipped"].append({"name": study_name, "reason": "0 samples"})
            manifest["warnings"].append(msg)
            continue

        study_df_sub = expr_df.loc[study_mask].copy().reset_index(drop=True)
        study_config = config[study_name]

        # Study DEG
        study_deg = None
        if deg_df is not None:
            deg_pid = "project_id" if "project_id" in deg_df.columns else "study"
            if deg_pid in deg_df.columns:
                study_deg = deg_df[deg_df[deg_pid] == study_name].copy()

        try:
            run_internal_study(
                study_name, study_df_sub, None, study_deg,
                study_config, target_genes, signatures,
                gene_cols, args.output_dir)
            manifest["studies_processed"].append(study_name)
        except Exception as e:
            msg = f"{study_name}: {str(e)}"
            print(f"  ERROR: {msg}", file=sys.stderr)
            manifest["studies_skipped"].append({"name": study_name, "reason": str(e)})
            manifest["warnings"].append(msg)

    # Process External Studies
    if external_studies:
        print("\n\n=== EXTERNAL STUDIES ===", file=sys.stderr)
        for study_id in external_studies:
            study_mask = (expr_df[pid_col] == study_id).values
            if study_mask.sum() == 0:
                msg = f"{study_id}: 0 samples in data"
                print(f"  SKIP: {msg}", file=sys.stderr)
                manifest["studies_skipped"].append({"name": study_id, "reason": "0 samples"})
                continue

            study_df_sub = expr_df.loc[study_mask].copy().reset_index(drop=True)
            study_deg = None
            if deg_df is not None:
                deg_pid = "project_id" if "project_id" in deg_df.columns else "study"
                if deg_pid in deg_df.columns:
                    study_deg = deg_df[deg_df[deg_pid] == study_id].copy()

            try:
                run_external_study(
                    study_id, study_df_sub, None, study_deg,
                    target_genes, signatures, gene_cols, args.output_dir)
                manifest["studies_processed"].append(study_id)
            except Exception as e:
                msg = f"{study_id}: {str(e)}"
                print(f"  ERROR: {msg}", file=sys.stderr)
                manifest["studies_skipped"].append({"name": study_id, "reason": str(e)})
                manifest["warnings"].append(msg)

    # Write manifest
    manifest_path = os.path.join(args.output_dir, "manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\n[Manifest] Saved: {manifest_path}", file=sys.stderr)

    print("\n" + "=" * 58, file=sys.stderr)
    print("SAMPLE-LEVEL ANALYSIS COMPLETE", file=sys.stderr)
    print("=" * 58, file=sys.stderr)
    print(f"  Studies processed: {len(manifest['studies_processed'])}", file=sys.stderr)
    print(f"  Studies skipped: {len(manifest['studies_skipped'])}", file=sys.stderr)
    print(f"  Warnings: {len(manifest['warnings'])}", file=sys.stderr)
    print(f"  Output: {args.output_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
