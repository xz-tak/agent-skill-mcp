#!/usr/bin/env python3
"""
Module 6: Differential Cascade Analysis

Performs statistical comparison of multicellular signaling cascades
(Ligand→Receptor→TF→Gene) between disease and normal conditions.

Cascade Structure:
    Ligand:Cell_A → Receptor:Cell_B → TF:Cell_B → Gene:Cell_B
         (CCC)         (Receptor-TF)       (GRN)

Analysis Approach:
- Analyzes ALL cascades in disease ∪ normal (not just shared)
- Uses simple difference (diff = score_disease - score_normal) as effect size
- Uses t-distribution with kurtosis-estimated df for p-values
- FDR correction via Benjamini-Hochberg method
- Edge-level statistics (CCC + GRN) with raw weights
- Cell-pair aggregation with Fisher's combined p-value

Analysis Levels:
1. Individual cascades - Full cascade statistics (diff, p-value, padj)
2. Cell-type pairs - Aggregated by sender→receiver (Fisher's combined p-value)
3. Individual edges - CCC and GRN edge statistics

Output Structure:
    {output_dir}/differential_cascades/{ccc_source}_ccc/
    ├── {disease_cond_1}/
    │   ├── cascade_results.csv
    │   ├── cellpair_results.csv
    │   └── edge_results.csv
    ├── {disease_cond_2}/
    │   └── (same files)
    └── summary_stats.json
"""

import gc
import json
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues, false_discovery_control, t as t_dist
from tqdm import tqdm

from config import ReconConfig, get_config
from recon.data import load_receptor_genes

warnings.filterwarnings("ignore")


def log(msg: str) -> None:
    """Print message with timestamp."""
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_ccc(condition: str, config: ReconConfig) -> pd.DataFrame:
    """
    Load CCC data for a condition (percentile-rank normalized).

    Returns DataFrame with columns: ligand, receptor, cell_source, cell_target, weight, condition
    """
    ccc_dir = config.get_ccc_dir()
    ccc_source = config.ccc_source

    if ccc_source == "cellchat":
        path = ccc_dir / "cellchat" / f"{condition}_ccc.csv"
        ccc = pd.read_csv(path)
        ccc = ccc.rename(columns={
            "source": "cell_source",
            "target": "cell_target",
        })
        ccc["weight"] = ccc["prob"]
    elif ccc_source == "cellphonedb":
        path = ccc_dir / "cellphonedb" / f"{condition}_ccc.csv"
        ccc = pd.read_csv(path)
        ccc = ccc.rename(columns={
            "source": "ligand",
            "target": "receptor",
            "celltype_source": "cell_source",
            "celltype_target": "cell_target",
        })
        ccc["weight"] = ccc["lr_means"]
    else:  # merged
        path = ccc_dir / "merged" / f"{condition}_ccc.csv"
        ccc = pd.read_csv(path)
        ccc = ccc.rename(columns={
            "source": "ligand",
            "target": "receptor",
            "celltype_source": "cell_source",
            "celltype_target": "cell_target",
        })
        ccc["weight"] = ccc["lr_means"]

    # Percentile-rank normalize weights to 0-1
    ccc["weight"] = ccc["weight"].rank(pct=True)

    ccc = ccc[["ligand", "receptor", "cell_source", "cell_target", "weight"]].copy()
    ccc["condition"] = condition
    return ccc


def _load_raw_ccc(condition: str, config: ReconConfig) -> pd.DataFrame:
    """
    Load RAW CCC data (not percentile-ranked) for edge-level statistics.

    Returns DataFrame with columns: ligand, receptor, cell_source, cell_target, weight
    """
    ccc_dir = config.get_ccc_dir()
    ccc_source = config.ccc_source

    if ccc_source == "cellchat":
        path = ccc_dir / "cellchat" / f"{condition}_ccc.csv"
        ccc = pd.read_csv(path)
        ccc = ccc.rename(columns={
            "source": "cell_source",
            "target": "cell_target",
        })
        ccc["weight"] = ccc["prob"]
    elif ccc_source == "cellphonedb":
        path = ccc_dir / "cellphonedb" / f"{condition}_ccc.csv"
        ccc = pd.read_csv(path)
        ccc = ccc.rename(columns={
            "source": "ligand",
            "target": "receptor",
            "celltype_source": "cell_source",
            "celltype_target": "cell_target",
        })
        ccc["weight"] = ccc["lr_means"]
    else:  # merged
        path = ccc_dir / "merged" / f"{condition}_ccc.csv"
        ccc = pd.read_csv(path)
        ccc = ccc.rename(columns={
            "source": "ligand",
            "target": "receptor",
            "celltype_source": "cell_source",
            "celltype_target": "cell_target",
        })
        ccc["weight"] = ccc["lr_means"]

    return ccc[["ligand", "receptor", "cell_source", "cell_target", "weight"]].copy()


def load_receptor_tf_network() -> pd.DataFrame:
    """
    Load receptor-TF relationships from NicheNet PKN.

    Returns DataFrame with columns: receptor, tf, weight
    """
    log("Loading receptor-TF network from NicheNet PKN...")
    receptor_grn = load_receptor_genes("human_receptor_gene_from_NichenetPKN")

    receptor_grn = receptor_grn.rename(columns={
        "source": "receptor",
        "target": "tf",
    })

    if "weight" not in receptor_grn.columns:
        receptor_grn["weight"] = 1.0

    log(f"  {len(receptor_grn):,} receptor-TF edges")
    log(f"  {receptor_grn['receptor'].nunique()} unique receptors")
    log(f"  {receptor_grn['tf'].nunique()} unique TFs")

    return receptor_grn


def load_grn(celltype: str, condition: str, config: ReconConfig) -> pd.DataFrame:
    """
    Load GRN for a specific cell type and condition.

    Returns DataFrame with columns: tf, gene, weight
    """
    grn_dir = config.get_grn_dir()
    min_weight = config.min_cascade_grn_weight
    ct_safe = celltype.lower().replace(" ", "_")
    grn_path = grn_dir / f"{ct_safe}_{condition}_rna_network.csv"

    if not grn_path.exists():
        return pd.DataFrame(columns=["tf", "gene", "weight"])

    grn = pd.read_csv(grn_path)
    grn = grn.rename(columns={"source": "tf", "target": "gene"})
    grn = grn[grn["weight"] >= min_weight].copy()

    return grn[["tf", "gene", "weight"]]


def load_all_grns(condition: str, config: ReconConfig) -> Dict[str, pd.DataFrame]:
    """Load GRNs for all available cell types for a condition."""
    grn_dir = config.get_grn_dir()

    # Auto-detect cell types from available GRN files
    pattern = f"*_{condition}_rna_network.csv"
    grn_files = sorted(grn_dir.glob(pattern))

    grns: Dict[str, pd.DataFrame] = {}
    for path in grn_files:
        ct_safe = path.stem.replace(f"_{condition}_rna_network", "")
        # Reconstruct display name from safe name
        ct = ct_safe
        grn = load_grn(ct, condition, config)
        if not grn.empty:
            grns[ct] = grn

    return grns


# ---------------------------------------------------------------------------
# Cascade enumeration
# ---------------------------------------------------------------------------

def enumerate_cascades(
    ccc: pd.DataFrame,
    receptor_tf: pd.DataFrame,
    grns: Dict[str, pd.DataFrame],
    condition: str,
    max_per_pair: int = 50000,
) -> pd.DataFrame:
    """
    Enumerate cascades using vectorized pandas merges.

    Cascade: Ligand:Cell_A -> Receptor:Cell_B -> TF:Cell_B -> Gene:Cell_B
    Limits cascades per cell-type pair to manage memory.
    """
    all_cascades = []

    # Pre-filter receptor-TF to receptors in CCC
    receptors_in_ccc = set(ccc["receptor"].unique())
    rtf_filtered = receptor_tf[receptor_tf["receptor"].isin(receptors_in_ccc)].copy()

    # Also filter to TFs that exist in ANY GRN
    all_tfs_in_grns = set()
    for grn in grns.values():
        all_tfs_in_grns.update(grn["tf"].unique())
    rtf_filtered = rtf_filtered[rtf_filtered["tf"].isin(all_tfs_in_grns)].copy()
    log(f"  RTF filtered to {len(rtf_filtered):,} edges (from {len(receptor_tf):,})")
    log(f"  (filtered to {len(receptors_in_ccc)} receptors, {len(all_tfs_in_grns)} TFs)")

    cell_targets = ccc["cell_target"].unique()

    for cell_target in tqdm(cell_targets, desc="Processing cell types"):
        if cell_target not in grns:
            continue

        grn = grns[cell_target]
        if grn.empty:
            continue

        ccc_ct = ccc[ccc["cell_target"] == cell_target].copy()

        # Filter RTF to TFs in this cell type's GRN
        tfs_in_grn = set(grn["tf"].unique())
        rtf_ct = rtf_filtered[rtf_filtered["tf"].isin(tfs_in_grn)]

        # Step 1: Merge CCC with Receptor-TF on receptor
        merged = ccc_ct.merge(
            rtf_ct, on="receptor", how="inner", suffixes=("_ccc", "_rtf")
        )
        if merged.empty:
            continue

        # Step 2: Merge with GRN on TF
        merged = merged.merge(grn, on="tf", how="inner", suffixes=("", "_grn"))
        if merged.empty:
            continue

        merged = merged.rename(columns={
            "weight_ccc": "ccc_weight",
            "weight_rtf": "rtf_weight",
            "weight": "grn_weight",
        })

        # Cascade score (normalized)
        merged["grn_weight_norm"] = np.minimum(merged["grn_weight"] / 10.0, 1.0)
        merged["cascade_score"] = (
            merged["ccc_weight"] * merged["rtf_weight"] * merged["grn_weight_norm"]
        )

        # Cascade ID
        merged["cascade_id"] = (
            merged["ligand"] + ":" + merged["cell_source"] + "→" +
            merged["receptor"] + ":" + merged["cell_target"] + "→" +
            merged["tf"] + ":" + merged["cell_target"] + "→" +
            merged["gene"] + ":" + merged["cell_target"]
        )

        # Limit cascades per cell-type pair
        for cell_source in merged["cell_source"].unique():
            pair_cascades = merged[merged["cell_source"] == cell_source]
            if len(pair_cascades) > max_per_pair:
                pair_cascades = pair_cascades.nlargest(max_per_pair, "cascade_score")
            all_cascades.append(pair_cascades)

    if not all_cascades:
        return pd.DataFrame()

    result = pd.concat(all_cascades, ignore_index=True)

    cols = [
        "cascade_id", "ligand", "receptor", "tf", "gene",
        "cell_source", "cell_target", "ccc_weight", "rtf_weight",
        "grn_weight", "cascade_score",
    ]
    result = result[cols].copy()
    result["condition"] = condition

    log(f"  Total cascades: {len(result):,}")
    return result


# ---------------------------------------------------------------------------
# Statistical analysis
# ---------------------------------------------------------------------------

def compute_effect_sizes(
    cascades_disease: pd.DataFrame,
    cascades_normal: pd.DataFrame,
    disease_name: str,
) -> pd.DataFrame:
    """
    Compute effect size (diff = score_disease - score_normal) for ALL cascades.

    Analyzes the union of disease and normal cascades. Cascades absent from
    one condition get score=0, making diff equal to the present score.
    """
    disease_scores = cascades_disease.set_index("cascade_id")["cascade_score"].to_dict()
    normal_scores = cascades_normal.set_index("cascade_id")["cascade_score"].to_dict()

    all_cascade_ids = set(disease_scores.keys()) | set(normal_scores.keys())

    log(f"  Total unique cascades (union): {len(all_cascade_ids):,}")
    log(f"  Disease cascades: {len(disease_scores):,}")
    log(f"  Normal cascades: {len(normal_scores):,}")

    # Build metadata lookup (disease first, then fill from normal)
    cascade_metadata: Dict[str, dict] = {}
    for source_df in (cascades_disease, cascades_normal):
        for _, row in source_df.iterrows():
            cid = row["cascade_id"]
            if cid not in cascade_metadata:
                cascade_metadata[cid] = {
                    "ligand": row["ligand"],
                    "receptor": row["receptor"],
                    "tf": row["tf"],
                    "gene": row["gene"],
                    "cell_source": row["cell_source"],
                    "cell_target": row["cell_target"],
                    "ccc_weight": row["ccc_weight"],
                    "rtf_weight": row["rtf_weight"],
                    "grn_weight": row["grn_weight"],
                }

    results = []
    for cascade_id in all_cascade_ids:
        score_d = disease_scores.get(cascade_id, 0)
        score_n = normal_scores.get(cascade_id, 0)
        meta = cascade_metadata.get(cascade_id)
        if meta is None:
            continue

        results.append({
            "cascade_id": cascade_id,
            **meta,
            "score_disease": score_d,
            "score_normal": score_n,
            "diff": score_d - score_n,
        })

    result_df = pd.DataFrame(results)

    n_disease_only = (result_df["score_normal"] == 0).sum()
    n_normal_only = (result_df["score_disease"] == 0).sum()
    n_shared = len(result_df) - n_disease_only - n_normal_only

    log(f"  Shared cascades (in both): {n_shared:,}")
    log(f"  Disease-specific (score_normal=0): {n_disease_only:,}")
    log(f"  Normal-specific (score_disease=0): {n_normal_only:,}")

    return result_df


def run_permutation_tests(
    effect_df: pd.DataFrame,
    disease_name: str,
    n_perm: int,
    n_jobs: int,
) -> pd.DataFrame:
    """
    Compute p-values using t-distribution with kurtosis-estimated df.

    Estimates degrees of freedom from excess kurtosis:
        df = 6 / kurtosis + 4  (for kurtosis > 0)
    """
    log("Computing t-distribution p-values...")

    all_diff = effect_df["diff"].values
    n_cascades = len(all_diff)

    mean_diff = np.mean(all_diff)
    std_diff = np.std(all_diff)
    z_scores = (all_diff - mean_diff) / std_diff if std_diff > 0 else np.zeros_like(all_diff)

    # Estimate df from excess kurtosis (Fisher's definition)
    kurtosis = pd.Series(all_diff).kurtosis()

    if kurtosis > 0:
        df_est = max(4.5, min(100, 6.0 / kurtosis + 4))
    else:
        df_est = 100

    log(f"  {n_cascades:,} cascades")
    log(f"  diff mean: {mean_diff:.6f}, std: {std_diff:.6f}")
    log(f"  kurtosis: {kurtosis:.4f}, estimated df: {df_est:.1f}")

    # Two-tailed p-value from t-distribution
    pvals = 2 * t_dist.sf(np.abs(z_scores), df=df_est)

    effect_df["pval"] = pvals
    effect_df["zscore"] = z_scores

    log(f"  Z-score range: [{z_scores.min():.2f}, {z_scores.max():.2f}]")
    log(f"  P-value range: [{pvals.min():.2e}, {pvals.max():.4f}]")
    log(f"  Cascades with p < 0.05: {(pvals < 0.05).sum():,}")
    log(f"  Cascades with p < 0.01: {(pvals < 0.01).sum():,}")
    log(f"  Cascades with p < 1e-5: {(pvals < 1e-5).sum():,}")

    return effect_df


def apply_fdr_correction(effect_df: pd.DataFrame, disease_name: str) -> pd.DataFrame:
    """Apply Benjamini-Hochberg FDR correction to p-values."""
    pvals = effect_df["pval"].values

    valid_mask = ~np.isnan(pvals)
    padj = np.full_like(pvals, np.nan)

    if valid_mask.sum() > 0:
        padj[valid_mask] = false_discovery_control(pvals[valid_mask], method="bh")

    effect_df["padj"] = padj
    return effect_df


# ---------------------------------------------------------------------------
# Edge-level statistics
# ---------------------------------------------------------------------------

def compute_edge_statistics(
    grn_disease: Dict[str, pd.DataFrame],
    grn_normal: Dict[str, pd.DataFrame],
    disease_name: str,
    config: ReconConfig,
) -> pd.DataFrame:
    """
    Compute edge-level diff/pval/padj for CCC and GRN edges.

    Compares RAW weights between disease and normal conditions
    for individual edges (not combined cascade scores).
    """
    log(f"Computing edge-level statistics for {disease_name.upper()} vs Normal...")
    edges = []

    # --- CCC EDGES ---
    log("  Loading RAW CCC data for edge statistics...")
    ccc_disease = _load_raw_ccc(disease_name, config)
    ccc_normal = _load_raw_ccc(config.normal_condition, config)
    log(f"    Disease CCC: {len(ccc_disease):,} edges")
    log(f"    Normal CCC: {len(ccc_normal):,} edges")

    log("  Processing CCC edges...")
    ccc_merged = ccc_disease.merge(
        ccc_normal,
        on=["ligand", "receptor", "cell_source", "cell_target"],
        how="outer",
        suffixes=("_disease", "_normal"),
    ).fillna(0)

    log(f"    {len(ccc_merged):,} CCC edges (union)")

    for _, row in ccc_merged.iterrows():
        edges.append({
            "edge_key": f"{row['ligand']}::{row['cell_source']}|{row['receptor']}::{row['cell_target']}|ccc",
            "edge_type": "ccc",
            "src": row["ligand"],
            "tgt": row["receptor"],
            "src_ct": row["cell_source"],
            "tgt_ct": row["cell_target"],
            "weight_disease": row["weight_disease"],
            "weight_normal": row["weight_normal"],
            "diff": row["weight_disease"] - row["weight_normal"],
        })

    # --- GRN EDGES ---
    log("  Processing GRN edges...")
    all_cell_types = set(grn_disease.keys()) | set(grn_normal.keys())

    for ct in all_cell_types:
        grn_d = grn_disease.get(ct, pd.DataFrame(columns=["tf", "gene", "weight"]))
        grn_n = grn_normal.get(ct, pd.DataFrame(columns=["tf", "gene", "weight"]))

        if grn_d.empty and grn_n.empty:
            continue

        grn_merged = grn_d.merge(
            grn_n, on=["tf", "gene"], how="outer", suffixes=("_disease", "_normal")
        ).fillna(0)

        for _, row in grn_merged.iterrows():
            edges.append({
                "edge_key": f"{row['tf']}::{ct}|{row['gene']}::{ct}|grn",
                "edge_type": "grn",
                "src": row["tf"],
                "tgt": row["gene"],
                "src_ct": ct,
                "tgt_ct": ct,
                "weight_disease": row["weight_disease"],
                "weight_normal": row["weight_normal"],
                "diff": row["weight_disease"] - row["weight_normal"],
            })

    log(f"    Total GRN edges: {sum(1 for e in edges if e['edge_type'] == 'grn'):,}")

    # --- STATISTICAL TESTING ---
    df = pd.DataFrame(edges)
    log(f"  Total edges: {len(df):,}")

    if len(df) == 0:
        return df

    all_diff = df["diff"].values
    mean_diff = np.mean(all_diff)
    std_diff = np.std(all_diff)

    z_scores = (all_diff - mean_diff) / std_diff if std_diff > 0 else np.zeros_like(all_diff)

    kurtosis = pd.Series(all_diff).kurtosis()
    df_est = max(4.5, min(100, 6.0 / kurtosis + 4)) if kurtosis > 0 else 100

    log(f"    diff mean: {mean_diff:.6f}, std: {std_diff:.6f}")
    log(f"    kurtosis: {kurtosis:.4f}, estimated df: {df_est:.1f}")

    pvals = 2 * t_dist.sf(np.abs(z_scores), df=df_est)

    df["pval"] = pvals
    df["zscore"] = z_scores

    # FDR correction
    valid_mask = ~np.isnan(pvals)
    padj = np.full_like(pvals, np.nan)
    if valid_mask.sum() > 0:
        padj[valid_mask] = false_discovery_control(pvals[valid_mask], method="bh")
    df["padj"] = padj

    log(f"    Edges with p < 0.05: {(pvals < 0.05).sum():,}")
    log(f"    Edges with padj < 0.05: {(df['padj'] < 0.05).sum():,}")

    return df


# ---------------------------------------------------------------------------
# Cell-pair aggregation
# ---------------------------------------------------------------------------

def aggregate_cellpair_level(
    effect_df: pd.DataFrame,
    disease_name: str,
    fdr_threshold: float,
) -> pd.DataFrame:
    """
    Aggregate cascades by cell-type pairs using Fisher's combined p-values.

    Columns:
    - mean_diff, median_diff, sum_diff
    - n_sig_positive, n_sig_negative, net_cascades
    - combined_pval, padj
    """
    log("Aggregating at cell-pair level...")

    cellpair_results = []
    grouped = effect_df.groupby(["cell_source", "cell_target"])

    for (cell_source, cell_target), group in tqdm(grouped, desc="Cell pairs"):
        n_cascades = len(group)

        pvals = group["pval"].dropna().values
        if len(pvals) >= 2:
            _, combined_pval = combine_pvalues(pvals, method="fisher")
        elif len(pvals) == 1:
            combined_pval = pvals[0]
        else:
            combined_pval = np.nan

        diffs = group["diff"].dropna()
        sig_positive = (group["padj"] < fdr_threshold) & (group["diff"] > 0)
        sig_negative = (group["padj"] < fdr_threshold) & (group["diff"] < 0)

        n_sig_pos = int(sig_positive.sum())
        n_sig_neg = int(sig_negative.sum())

        cellpair_results.append({
            "cell_source": cell_source,
            "cell_target": cell_target,
            "n_cascades": n_cascades,
            "mean_diff": diffs.mean() if len(diffs) > 0 else np.nan,
            "median_diff": diffs.median() if len(diffs) > 0 else np.nan,
            "sum_diff": diffs.sum() if len(diffs) > 0 else 0,
            "n_sig_positive": n_sig_pos,
            "n_sig_negative": n_sig_neg,
            "net_cascades": n_sig_pos - n_sig_neg,
            "combined_pval": combined_pval,
        })

    cellpair_df = pd.DataFrame(cellpair_results)

    # FDR on combined p-values
    if "combined_pval" in cellpair_df.columns and len(cellpair_df) > 0:
        pvals = cellpair_df["combined_pval"].values
        valid_mask = ~np.isnan(pvals)
        padj = np.full_like(pvals, np.nan)
        if valid_mask.sum() > 0:
            padj[valid_mask] = false_discovery_control(pvals[valid_mask], method="bh")
        cellpair_df["padj"] = padj

    return cellpair_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(config: Optional[ReconConfig] = None, disease: str = "both") -> Dict:
    """
    Run differential cascade analysis.

    Args:
        config: Pipeline configuration (falls back to CLI args if None)
        disease: Which disease condition(s) to analyze vs normal.
                 "both" processes all disease_conditions from config.
                 Otherwise must match one of config.disease_conditions.

    Returns:
        Summary dict with statistics for each disease comparison.
    """
    config = get_config(config)

    log("=" * 70)
    log("DIFFERENTIAL CASCADE ANALYSIS")
    log("=" * 70)
    log(f"CCC Source: {config.ccc_source}")
    log(f"Permutations: {config.n_permutations}")
    log(f"FDR threshold: {config.fdr_threshold}")
    log(f"Edge weight threshold: {config.edge_weight_threshold}")
    log(f"Min GRN weight: {config.min_cascade_grn_weight}")
    log(f"Max cascades/cell-pair: {config.max_cascades_per_cellpair}")
    log(f"Normal condition: {config.normal_condition}")
    log(f"Disease conditions: {config.disease_conditions}")
    log("=" * 70)

    start_time = datetime.now()

    base_output_dir = config.get_cascade_dir()
    base_output_dir.mkdir(parents=True, exist_ok=True)
    log(f"Output directory: {base_output_dir}")

    # Load receptor-TF network (shared)
    receptor_tf = load_receptor_tf_network()

    # Determine which conditions to load
    if disease == "both":
        diseases_to_process = list(config.disease_conditions)
    else:
        diseases_to_process = [disease]

    conditions_to_load = list(set(diseases_to_process + [config.normal_condition]))

    # Load data and enumerate cascades for each condition
    all_cascades: Dict[str, pd.DataFrame] = {}
    all_grns: Dict[str, Dict[str, pd.DataFrame]] = {}

    for condition in conditions_to_load:
        log(f"\n{'=' * 60}")
        log(f"PROCESSING: {condition.upper()}")
        log(f"{'=' * 60}")

        log(f"Loading CCC for {condition}...")
        ccc = load_ccc(condition, config)
        log(f"  {len(ccc):,} CCC edges")

        log(f"Loading GRNs for {condition}...")
        grns = load_all_grns(condition, config)
        total_grn_edges = sum(len(g) for g in grns.values())
        log(f"  {len(grns)} cell types, {total_grn_edges:,} total GRN edges")

        all_grns[condition] = grns

        log(f"Enumerating cascades for {condition}...")
        cascades = enumerate_cascades(
            ccc, receptor_tf, grns, condition,
            max_per_pair=config.max_cascades_per_cellpair,
        )
        log(f"  {len(cascades):,} cascades enumerated")

        all_cascades[condition] = cascades

    # Compute statistics for each disease vs normal
    all_results: Dict[str, pd.DataFrame] = {}
    summary: Dict = {
        "ccc_source": config.ccc_source,
        "n_permutations": config.n_permutations,
        "fdr_threshold": config.fdr_threshold,
        "edge_weight_threshold": config.edge_weight_threshold,
        "min_grn_weight": config.min_cascade_grn_weight,
        "n_cascades": {cond: len(df) for cond, df in all_cascades.items()},
        "compute_time_seconds": 0,
        "timestamp": datetime.now().isoformat(),
        "diseases": {},
    }

    normal_condition = config.normal_condition

    for disease_name in diseases_to_process:
        if disease_name not in all_cascades or normal_condition not in all_cascades:
            log(f"Skipping {disease_name}: missing data")
            continue

        log(f"\n{'=' * 60}")
        log(f"COMPUTING STATISTICS: {disease_name.upper()} vs {normal_condition.upper()}")
        log(f"{'=' * 60}")

        disease_dir = base_output_dir / disease_name
        disease_dir.mkdir(parents=True, exist_ok=True)
        log(f"Output directory: {disease_dir}")

        # Effect sizes
        log("Computing effect sizes...")
        effect_df = compute_effect_sizes(
            all_cascades[disease_name],
            all_cascades[normal_condition],
            disease_name,
        )

        if len(effect_df) == 0:
            log("WARNING: No cascades found!")
            continue

        log(f"  {len(effect_df):,} total cascades for statistical testing")

        # Permutation tests (t-distribution)
        effect_df = run_permutation_tests(
            effect_df, disease_name, config.n_permutations, config.n_jobs
        )

        # FDR correction
        log("Applying FDR correction...")
        effect_df = apply_fdr_correction(effect_df, disease_name)

        sig_positive = int(((effect_df["padj"] < config.fdr_threshold) & (effect_df["diff"] > 0)).sum())
        sig_negative = int(((effect_df["padj"] < config.fdr_threshold) & (effect_df["diff"] < 0)).sum())
        log(f"  Significant: {sig_positive} increased in disease, {sig_negative} decreased")

        # Diff statistics
        log(f"  diff statistics:")
        log(f"    Mean: {effect_df['diff'].mean():.6f}")
        log(f"    Median: {effect_df['diff'].median():.6f}")
        log(f"    Range: [{effect_df['diff'].min():.6f}, {effect_df['diff'].max():.6f}]")

        all_results[disease_name] = effect_df

        # Save cascade results
        cascade_path = disease_dir / "cascade_results.csv"
        effect_df.to_csv(cascade_path, index=False)
        log(f"\nSaved: {cascade_path}")
        log(f"  {len(effect_df):,} total cascades")

        # Cell-pair aggregation
        cellpair_df = aggregate_cellpair_level(
            effect_df, disease_name, config.fdr_threshold
        )
        cellpair_path = disease_dir / "cellpair_results.csv"
        cellpair_df.to_csv(cellpair_path, index=False)
        log(f"\nSaved: {cellpair_path}")
        log(f"  {len(cellpair_df):,} cell pairs")

        # Edge-level statistics
        log(f"\n{'=' * 50}")
        edge_df = compute_edge_statistics(
            grn_disease=all_grns[disease_name],
            grn_normal=all_grns[normal_condition],
            disease_name=disease_name,
            config=config,
        )
        edge_path = disease_dir / "edge_results.csv"
        edge_df.to_csv(edge_path, index=False)
        log(f"\nSaved: {edge_path}")
        log(f"  {len(edge_df):,} edges")
        n_ccc_edges = int((edge_df["edge_type"] == "ccc").sum())
        n_grn_edges = int((edge_df["edge_type"] == "grn").sum())
        log(f"  ({n_ccc_edges:,} CCC edges, {n_grn_edges:,} GRN edges)")

        # Summary for this disease
        n_disease_specific = int((effect_df["score_normal"] == 0).sum())
        n_normal_specific = int((effect_df["score_disease"] == 0).sum())
        n_shared_cascades = len(effect_df) - n_disease_specific - n_normal_specific
        n_p_lt_05 = int((effect_df["pval"] < 0.05).sum())
        n_p_lt_01 = int((effect_df["pval"] < 0.01).sum())

        summary["diseases"][disease_name] = {
            "n_total_cascades": len(effect_df),
            "n_shared_cascades": n_shared_cascades,
            "n_disease_specific": n_disease_specific,
            "n_normal_specific": n_normal_specific,
            "n_p_lt_05": n_p_lt_05,
            "n_p_lt_01": n_p_lt_01,
            "sig_positive": sig_positive,
            "sig_negative": sig_negative,
            "n_cellpairs": len(cellpair_df),
            "diff_mean": float(effect_df["diff"].mean()),
            "diff_median": float(effect_df["diff"].median()),
            "diff_range": [float(effect_df["diff"].min()), float(effect_df["diff"].max())],
        }

    # Save summary
    elapsed = datetime.now() - start_time
    summary["compute_time_seconds"] = elapsed.total_seconds()

    summary_path = base_output_dir / "summary_stats.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    log(f"\nSaved: {summary_path}")

    log("\n" + "=" * 70)
    log("DIFFERENTIAL CASCADE ANALYSIS COMPLETE")
    log("=" * 70)
    log(f"Total time: {elapsed}")
    log(f"Output directory: {base_output_dir}")

    log("\nCondition-specific outputs:")
    for d in diseases_to_process:
        disease_dir = base_output_dir / d
        if disease_dir.exists():
            log(f"\n  {d}/")
            for f_path in sorted(disease_dir.glob("*.csv")):
                log(f"    - {f_path.name}")
    log(f"\n  summary_stats.json")

    return summary


if __name__ == "__main__":
    main()
