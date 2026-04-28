#!/usr/bin/env python3
"""
Module 3: Cell-Cell Communication Analysis

Combines CellPhoneDB (via LIANA+), CellChat import, and merge logic into a
single generalized module.

Functions:
- run_cellphonedb: Compute CellPhoneDB via LIANA+ for one condition
- import_cellchat: Parse/standardize pre-computed CellChat R output
- merge_ccc_methods: Union with percentile-rank logic
- main: Orchestrate full CCC pipeline

Expects:
- CellPhoneDB: computed from scRNA-seq AnnData
- CellChat: pre-computed CSVs from R (one combined file or per-condition)
- Merged: union of both with OVERALL percentile rank as lr_means
"""

import gc
import json
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from config import ReconConfig, get_config

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Column conventions
# ---------------------------------------------------------------------------
# ReCoN expects: source (ligand gene), target (receptor gene),
#                celltype_source, celltype_target, lr_means
KEY_COLS = ["celltype_source", "celltype_target", "source", "target"]


# ---------------------------------------------------------------------------
# 1. CellPhoneDB via LIANA+
# ---------------------------------------------------------------------------

def run_cellphonedb(
    adata,
    condition_name: str,
    config: ReconConfig,
) -> Tuple[pd.DataFrame, Dict]:
    """
    Run CellPhoneDB via LIANA+ for a single condition.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object for this condition.
    condition_name : str
        Label for this condition (e.g. "ssc", "normal").
    config : ReconConfig
        Pipeline configuration.

    Returns
    -------
    ccc_results : pd.DataFrame
        CCC interactions with standardised columns.
    stats : dict
        Summary statistics for logging.
    """
    import anndata as ad
    import liana as li

    print(f"\n{'=' * 60}")
    print(f"[{datetime.now():%H:%M:%S}] RUNNING CELLPHONEDB FOR: {condition_name.upper()}")
    print(f"{'=' * 60}")

    start_time = datetime.now()
    celltype_col = config.celltype_col

    if celltype_col not in adata.obs.columns:
        raise ValueError(f"Cell type column '{celltype_col}' not found in obs")

    celltypes = adata.obs[celltype_col].unique().tolist()
    print(f"Cell types: {len(celltypes)}")
    for ct in celltypes:
        n_cells = int((adata.obs[celltype_col] == ct).sum())
        print(f"  {ct}: {n_cells:,} cells")

    print(f"\nRunning LIANA+ CellPhoneDB...")
    print(f"  Resource: {config.resource_name}")
    print(f"  Expression proportion threshold: {config.expr_prop}")

    li.method.cellphonedb(
        adata,
        resource_name=config.resource_name,
        expr_prop=config.expr_prop,
        groupby=celltype_col,
        key_added="cpdb_res",
        use_raw=False,
        verbose=True,
    )

    ccc_results = adata.uns["cpdb_res"].copy()
    print(f"CellPhoneDB complete: {len(ccc_results):,} interactions")

    # Standardise columns: LIANA uses ligand/receptor for genes,
    # source/target for cell types.  ReCoN expects the opposite.
    ccc_results = ccc_results.rename(columns={
        "ligand": "source",
        "receptor": "target",
        "source": "celltype_source",
        "target": "celltype_target",
    })

    # Filter by minimum lr_means
    initial_count = len(ccc_results)
    ccc_results = ccc_results[ccc_results["lr_means"] > config.min_lr_means]
    print(f"Filtered to {len(ccc_results):,} interactions (lr_means > {config.min_lr_means})")

    ccc_results["condition"] = condition_name

    elapsed = datetime.now() - start_time
    stats = {
        "condition": condition_name,
        "n_cells": int(adata.n_obs),
        "n_celltypes": len(celltypes),
        "celltypes": celltypes,
        "n_interactions_raw": int(initial_count),
        "n_interactions_filtered": int(len(ccc_results)),
        "n_unique_ligands": int(ccc_results["source"].nunique()) if len(ccc_results) > 0 else 0,
        "n_unique_receptors": int(ccc_results["target"].nunique()) if len(ccc_results) > 0 else 0,
        "compute_time_seconds": elapsed.total_seconds(),
    }

    print(f"\nCCC Statistics for {condition_name}:")
    for key, value in stats.items():
        if isinstance(value, list):
            continue
        elif isinstance(value, float):
            print(f"  {key}: {value:.2f}")
        else:
            print(f"  {key}: {value}")

    return ccc_results, stats


# ---------------------------------------------------------------------------
# 2. Import pre-computed CellChat R output
# ---------------------------------------------------------------------------

def import_cellchat(config: ReconConfig) -> Dict[str, pd.DataFrame]:
    """
    Parse and standardise pre-computed CellChat CSV(s).

    Handles two layouts:
    - Single combined file with a "Condition" column (split by condition)
    - Per-condition files already split ({condition}_ccc.csv)

    Column standardisation:
    - CellChat source/target (cell types) -> celltype_source / celltype_target
    - CellChat ligand_name/receptor_name or ligand/receptor -> source / target (genes)

    Parameters
    ----------
    config : ReconConfig

    Returns
    -------
    dict mapping condition name -> DataFrame with standardised columns.
    """
    print(f"\n{'=' * 60}")
    print(f"[{datetime.now():%H:%M:%S}] IMPORTING CELLCHAT DATA")
    print(f"{'=' * 60}")

    cellchat_dir = config.get_ccc_dir() / "cellchat"
    results: Dict[str, pd.DataFrame] = {}

    # Strategy 1: per-condition files already exist
    per_condition_found = False
    for condition in config.conditions:
        path = cellchat_dir / f"{condition}_ccc.csv"
        if path.exists():
            per_condition_found = True
            df = pd.read_csv(path)
            df = _standardise_cellchat_columns(df)
            df["condition"] = condition
            results[condition] = df
            print(f"  {condition}: {len(df):,} interactions (from {path.name})")

    if per_condition_found:
        return results

    # Strategy 2: single combined file (cellchat_data_path in config)
    src_path = None
    if config.cellchat_data_path:
        src_path = Path(config.cellchat_data_path)
    if src_path is None or not src_path.exists():
        print("  No CellChat data found. Skipping.")
        return results

    print(f"  Loading combined file: {src_path}")
    df = pd.read_csv(src_path)
    print(f"  Total rows: {len(df):,}")

    # Detect condition column (CellChat uses "Condition" by default)
    cond_col = None
    for candidate in ["Condition", "condition", "group"]:
        if candidate in df.columns:
            cond_col = candidate
            break

    if cond_col is None:
        raise ValueError(f"Cannot find condition column in CellChat file. Columns: {list(df.columns)}")

    print(f"  Conditions in file: {df[cond_col].unique().tolist()}")

    # Build a flexible condition mapping (case-insensitive matching)
    raw_conditions = df[cond_col].unique().tolist()
    cond_map = {}
    for raw in raw_conditions:
        lower = raw.strip().lower()
        for our_cond in config.conditions:
            if lower == our_cond.lower():
                cond_map[raw] = our_cond
                break

    # Split and save per-condition files
    cellchat_dir.mkdir(parents=True, exist_ok=True)
    for raw_cond, our_cond in cond_map.items():
        subset = df[df[cond_col] == raw_cond].copy()
        if len(subset) == 0:
            continue
        subset = _standardise_cellchat_columns(subset)
        subset["condition"] = our_cond
        out_path = cellchat_dir / f"{our_cond}_ccc.csv"
        subset.to_csv(out_path, index=False)
        results[our_cond] = subset
        print(f"  {our_cond}: {len(subset):,} rows -> {out_path.name}")

    return results


def _standardise_cellchat_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardise CellChat columns to ReCoN convention.

    CellChat uses:
    - source/target for cell types
    - ligand_name/receptor_name (or ligand/receptor) for genes

    ReCoN expects:
    - celltype_source/celltype_target for cell types
    - source/target for genes
    """
    df = df.copy()

    # Step 1: rename cell type columns first (before gene columns clash)
    rename_map = {}

    # Cell type columns: source/target -> celltype_source/celltype_target
    if "source" in df.columns and "ligand" in df.columns:
        rename_map["source"] = "celltype_source"
        rename_map["target"] = "celltype_target"
    elif "source" in df.columns and "ligand_name" in df.columns:
        rename_map["source"] = "celltype_source"
        rename_map["target"] = "celltype_target"

    if rename_map:
        df = df.rename(columns=rename_map)

    # Step 2: rename gene columns -> source/target
    gene_rename = {}
    if "ligand_name" in df.columns:
        gene_rename["ligand_name"] = "source"
        gene_rename["receptor_name"] = "target"
    elif "ligand" in df.columns:
        gene_rename["ligand"] = "source"
        gene_rename["receptor"] = "target"

    if gene_rename:
        df = df.rename(columns=gene_rename)

    return df


# ---------------------------------------------------------------------------
# 3. Merge CCC methods (percentile-rank union)
# ---------------------------------------------------------------------------

def merge_ccc_methods(
    cpdb_results: Dict[str, pd.DataFrame],
    cellchat_results: Dict[str, pd.DataFrame],
    config: ReconConfig,
) -> Dict[str, pd.DataFrame]:
    """
    Merge CellPhoneDB and CellChat results using percentile-rank logic.

    1. Compute OVERALL percentile ranks across ALL conditions per method
    2. Union both methods
    3. For duplicates (same celltype_source/celltype_target/source/target),
       keep the interaction with HIGHER overall percentile rank
    4. Use overall percentile rank as lr_means for downstream

    Parameters
    ----------
    cpdb_results : dict
        condition -> DataFrame from CellPhoneDB.
    cellchat_results : dict
        condition -> DataFrame from CellChat.
    config : ReconConfig

    Returns
    -------
    dict mapping condition -> merged DataFrame.
    """
    print(f"\n{'=' * 60}")
    print(f"[{datetime.now():%H:%M:%S}] MERGING CCC METHODS (PERCENTILE-RANK)")
    print(f"{'=' * 60}")

    # Concatenate all conditions per method
    cpdb_all = _concat_conditions(cpdb_results, "cellphonedb")
    cc_all = _concat_conditions(cellchat_results, "cellchat")

    if cpdb_all.empty and cc_all.empty:
        print("  No CCC data to merge.")
        return {}

    # Compute OVERALL percentile ranks across all conditions
    if not cpdb_all.empty:
        cpdb_all["lr_means_cellphonedb"] = cpdb_all["lr_means"]
        cpdb_all["pct_rank_overall"] = cpdb_all["lr_means_cellphonedb"].rank(pct=True)
        print(f"  CellPhoneDB: {len(cpdb_all):,} interactions, "
              f"lr_means range: {cpdb_all['lr_means_cellphonedb'].min():.3f} - "
              f"{cpdb_all['lr_means_cellphonedb'].max():.3f}")

    if not cc_all.empty:
        score_col = "prob" if "prob" in cc_all.columns else "lr_means"
        cc_all["prob_cellchat"] = cc_all[score_col]
        cc_all["pct_rank_overall"] = cc_all["prob_cellchat"].rank(pct=True)
        print(f"  CellChat: {len(cc_all):,} interactions, "
              f"prob range: {cc_all['prob_cellchat'].min():.6f} - "
              f"{cc_all['prob_cellchat'].max():.6f}")

    # Compute per-condition percentile ranks (for reference)
    for condition in config.conditions:
        if not cpdb_all.empty:
            mask = cpdb_all["condition"] == condition
            if mask.sum() > 0:
                cpdb_all.loc[mask, "pct_rank_percondition"] = (
                    cpdb_all.loc[mask, "lr_means_cellphonedb"].rank(pct=True)
                )
        if not cc_all.empty:
            mask = cc_all["condition"] == condition
            if mask.sum() > 0:
                cc_all.loc[mask, "pct_rank_percondition"] = (
                    cc_all.loc[mask, "prob_cellchat"].rank(pct=True)
                )

    # Merge per condition
    merged_dir = config.get_ccc_dir() / "merged"
    merged_dir.mkdir(parents=True, exist_ok=True)

    results: Dict[str, pd.DataFrame] = {}

    for condition in config.conditions:
        print(f"\n  Merging {condition.upper()}...")
        cpdb_cond = cpdb_all[cpdb_all["condition"] == condition].copy() if not cpdb_all.empty else pd.DataFrame()
        cc_cond = cc_all[cc_all["condition"] == condition].copy() if not cc_all.empty else pd.DataFrame()

        if cpdb_cond.empty and cc_cond.empty:
            print(f"    No data for {condition}")
            results[condition] = pd.DataFrame()
            continue

        # Identify overlapping interactions
        both_keys = set()
        if not cpdb_cond.empty and not cc_cond.empty:
            cpdb_keys = set(cpdb_cond[KEY_COLS].apply(tuple, axis=1))
            cc_keys = set(cc_cond[KEY_COLS].apply(tuple, axis=1))
            both_keys = cpdb_keys & cc_keys
            print(f"    CellPhoneDB-only: {len(cpdb_keys - cc_keys):,}")
            print(f"    CellChat-only: {len(cc_keys - cpdb_keys):,}")
            print(f"    Overlapping: {len(both_keys):,}")

        # Combine and deduplicate: keep HIGHER overall percentile rank
        combined = pd.concat([cpdb_cond, cc_cond], ignore_index=True)
        combined = combined.sort_values("pct_rank_overall", ascending=False)
        combined = combined.drop_duplicates(subset=KEY_COLS, keep="first")

        # Mark interactions from both sources
        combined["key"] = combined[KEY_COLS].apply(tuple, axis=1)
        combined.loc[combined["key"].isin(both_keys), "source_method"] = "both"
        combined = combined.drop(columns=["key"])

        # Use OVERALL pct_rank as lr_means for downstream
        combined["lr_means"] = combined["pct_rank_overall"]

        if "pval" in combined.columns:
            combined["pval"] = combined["pval"].fillna(1.0)

        # Select output columns
        output_cols = [
            "source", "target", "celltype_source", "celltype_target",
            "lr_means", "pval", "source_method", "condition",
            "lr_means_cellphonedb", "prob_cellchat",
            "pct_rank_percondition", "pct_rank_overall",
            "pathway_name", "annotation",
        ]
        output_cols = [c for c in output_cols if c in combined.columns]
        combined = combined[output_cols]

        print(f"    Total merged: {len(combined):,}")
        if "source_method" in combined.columns:
            print(f"    Method breakdown: {combined['source_method'].value_counts().to_dict()}")

        # Save
        out_path = merged_dir / f"{condition}_ccc.csv"
        combined.to_csv(out_path, index=False)
        print(f"    Saved: {out_path}")

        results[condition] = combined

    return results


def _concat_conditions(
    per_condition: Dict[str, pd.DataFrame],
    method_label: str,
) -> pd.DataFrame:
    """Concatenate per-condition DataFrames, adding source_method label."""
    frames = []
    for condition, df in per_condition.items():
        if df is not None and not df.empty:
            chunk = df.copy()
            chunk["source_method"] = method_label
            if "condition" not in chunk.columns:
                chunk["condition"] = condition
            frames.append(chunk)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# 4. Main orchestrator
# ---------------------------------------------------------------------------

def main(config: Optional[ReconConfig] = None) -> Dict[str, pd.DataFrame]:
    """
    Run the full CCC pipeline: compute CellPhoneDB, import CellChat, merge.

    Parameters
    ----------
    config : ReconConfig or None
        If None, builds from CLI args.

    Returns
    -------
    dict mapping condition -> final CCC DataFrame (source depends on config).
    """
    config = get_config(config)

    print("=" * 60)
    print(f"[{datetime.now():%H:%M:%S}] MODULE 3: CELL-CELL COMMUNICATION ANALYSIS")
    print(f"  Compute methods: {config.ccc_compute_methods}")
    print(f"  Downstream source: {config.ccc_source}")
    print("=" * 60)

    start_time = datetime.now()

    ccc_dir = config.get_ccc_dir()
    ccc_dir.mkdir(parents=True, exist_ok=True)

    all_stats: Dict = {}
    cpdb_results: Dict[str, pd.DataFrame] = {}
    cellchat_results: Dict[str, pd.DataFrame] = {}

    # --- CellPhoneDB ---
    if "cellphonedb" in config.ccc_compute_methods:
        import scanpy as sc
        import anndata as ad

        cpdb_dir = ccc_dir / "cellphonedb"
        cpdb_dir.mkdir(parents=True, exist_ok=True)

        for condition in config.conditions:
            print(f"\n{'#' * 60}")
            print(f"[{datetime.now():%H:%M:%S}] CellPhoneDB: {condition.upper()}")

            # Load condition-specific h5ad (convention: adata_{condition}.h5ad)
            adata_path = Path(config.output_dir) / f"adata_{condition}.h5ad"
            if not adata_path.exists():
                print(f"  Warning: {adata_path} not found, skipping {condition}")
                continue

            adata = sc.read_h5ad(adata_path)
            print(f"  Loaded {condition}: {adata.shape}")

            ccc, stats = run_cellphonedb(adata, condition, config)
            cpdb_results[condition] = ccc
            all_stats[f"cellphonedb_{condition}"] = stats

            # Save per-condition checkpoint
            out_path = cpdb_dir / f"{condition}_ccc.csv"
            ccc.to_csv(out_path, index=False)
            print(f"  Saved: {out_path}")

            del adata
            gc.collect()

    # --- CellChat ---
    if "cellchat" in config.ccc_compute_methods:
        cellchat_results = import_cellchat(config)

    # --- Merge ---
    merged_results: Dict[str, pd.DataFrame] = {}
    if cpdb_results or cellchat_results:
        merged_results = merge_ccc_methods(cpdb_results, cellchat_results, config)

    # Determine which results to return based on ccc_source
    if config.ccc_source == "merged":
        final = merged_results
    elif config.ccc_source == "cellphonedb":
        final = cpdb_results
    elif config.ccc_source == "cellchat":
        final = cellchat_results
    else:
        final = merged_results

    # Save stats
    stats_path = ccc_dir / "ccc_stats.json"
    with open(stats_path, "w") as f:
        json.dump(all_stats, f, indent=2, default=str)
    print(f"\n[{datetime.now():%H:%M:%S}] Saved stats: {stats_path}")

    elapsed = datetime.now() - start_time
    print(f"\n{'=' * 60}")
    print(f"MODULE 3 COMPLETE ({elapsed})")
    print(f"{'=' * 60}")
    for cond, df in final.items():
        n = len(df) if df is not None and not df.empty else 0
        print(f"  {cond}: {n:,} interactions")

    return final


if __name__ == "__main__":
    main()
