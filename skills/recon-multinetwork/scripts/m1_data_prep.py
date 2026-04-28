#!/usr/bin/env python3
"""
M1: Data Preparation for ReCoN Multicellular Network Analysis

Generalised from the SSc Lung Atlas pipeline. This module:
1. Loads an integrated scRNA h5ad file
2. Validates required columns (condition, cell type)
3. Subsets by condition (disease vs normal)
4. Saves per-condition AnnData objects + metadata

Usage:
    # With JSON config
    python m1_data_prep.py --config config.json

    # With CLI args
    python m1_data_prep.py --h5ad path.h5ad --condition-col condition \
        --celltype-col cluster_l2 --disease-conditions ssc ipf \
        --normal-condition normal --output-dir results/
"""

import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import anndata as ad
import scanpy as sc

from config import ReconConfig, get_config


def _log(msg: str) -> None:
    """Print timestamped message."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {msg}", flush=True)


def load_and_validate_data(config: ReconConfig) -> ad.AnnData:
    """
    Load the h5ad file and validate required columns.

    Args:
        config: Pipeline configuration.

    Returns:
        Loaded AnnData object.

    Raises:
        FileNotFoundError: If h5ad_path does not exist.
        ValueError: If required columns are missing.
    """
    h5ad_path = Path(config.h5ad_path)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"h5ad file not found: {h5ad_path}")

    _log(f"Loading data from: {h5ad_path}")
    adata = sc.read_h5ad(str(h5ad_path))

    _log(f"Total cells: {adata.n_obs:,}")
    _log(f"Total genes: {adata.n_vars:,}")

    # Validate cell type column
    if config.celltype_col not in adata.obs.columns:
        raise ValueError(
            f"Cell type column '{config.celltype_col}' not found in obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    celltypes = adata.obs[config.celltype_col].value_counts()
    _log(f"Cell types ({config.celltype_col}): {len(celltypes)} types")
    for ct, count in celltypes.items():
        _log(f"  {ct}: {count:,} ({count / adata.n_obs * 100:.1f}%)")

    # Validate condition column
    if config.condition_col not in adata.obs.columns:
        raise ValueError(
            f"Condition column '{config.condition_col}' not found in obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    conditions = adata.obs[config.condition_col].value_counts()
    _log(f"Conditions ({config.condition_col}):")
    for cond, count in conditions.items():
        _log(f"  {cond}: {count:,} ({count / adata.n_obs * 100:.1f}%)")

    # Validate requested conditions exist
    available_conditions = set(adata.obs[config.condition_col].unique())
    all_requested = set(config.disease_conditions) | {config.normal_condition}
    missing = all_requested - available_conditions
    if missing:
        raise ValueError(
            f"Requested conditions not found in data: {missing}. "
            f"Available: {available_conditions}"
        )

    # Report layers
    _log(f"Layers available: {list(adata.layers.keys())}")
    if "counts" not in adata.layers and "raw_counts" not in adata.layers:
        _log("WARNING: No raw counts layer found. GRN inference may need raw counts.")

    return adata


def subset_conditions(
    adata: ad.AnnData, config: ReconConfig
) -> Dict[str, ad.AnnData]:
    """
    Subset data by condition into separate AnnData objects.

    Args:
        adata: Full AnnData object.
        config: Pipeline configuration.

    Returns:
        Dict mapping condition name (lowercase) to subsetted AnnData.
    """
    _log("Subsetting by condition")

    condition_data: Dict[str, ad.AnnData] = {}

    # Process each condition (disease + normal)
    for condition in config.conditions:
        mask = adata.obs[config.condition_col] == condition
        n_cells = mask.sum()

        if n_cells == 0:
            _log(f"  WARNING: No cells found for condition '{condition}', skipping")
            continue

        subset = adata[mask].copy()
        key = condition.lower().replace(" ", "_")

        # Add comparison_group annotation for downstream analysis
        subset.obs['comparison_group'] = key

        condition_data[key] = subset
        _log(f"  {condition}: {n_cells:,} cells")

        # Cell type distribution
        ct_counts = subset.obs[config.celltype_col].value_counts()
        for ct, count in ct_counts.items():
            _log(f"    {ct}: {count:,}")

    return condition_data


def save_results(
    condition_data: Dict[str, ad.AnnData], config: ReconConfig
) -> None:
    """
    Save per-condition AnnData objects and metadata.

    Args:
        condition_data: Dict mapping condition key to AnnData.
        config: Pipeline configuration.
    """
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    _log(f"Saving results to: {output_dir}")

    for key, adata_cond in condition_data.items():
        out_path = output_dir / f"adata_{key}.h5ad"
        _log(f"  Saving {key}: {adata_cond.shape} -> {out_path}")
        adata_cond.write_h5ad(out_path)

    # Save metadata summary
    metadata = {
        "input_file": config.h5ad_path,
        "condition_col": config.condition_col,
        "celltype_col": config.celltype_col,
        "conditions": list(condition_data.keys()),
        "cells_per_condition": {
            k: int(v.n_obs) for k, v in condition_data.items()
        },
        "total_genes": int(next(iter(condition_data.values())).n_vars),
        "celltypes": sorted(
            list(
                next(iter(condition_data.values()))
                .obs[config.celltype_col]
                .unique()
            )
        ),
        "layers": list(next(iter(condition_data.values())).layers.keys()),
        "celltypes_per_condition": {
            k: v.obs[config.celltype_col].value_counts().to_dict()
            for k, v in condition_data.items()
        },
    }

    metadata_path = output_dir / "data_prep_metadata.json"
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)
    _log(f"  Saved metadata: {metadata_path}")


def main(config: Optional[ReconConfig] = None) -> Dict[str, ad.AnnData]:
    """
    Run the data preparation pipeline.

    Args:
        config: Optional ReconConfig. Falls back to CLI args if None.

    Returns:
        Dict mapping condition key to subsetted AnnData.
    """
    config = get_config(config)

    _log("=" * 60)
    _log("M1: DATA PREPARATION")
    _log("=" * 60)
    start_time = datetime.now()

    # Validate config
    if not config.h5ad_path:
        _log("ERROR: h5ad_path is required")
        sys.exit(1)
    if not config.disease_conditions:
        _log("ERROR: disease_conditions is required (at least one)")
        sys.exit(1)

    # Load and validate
    adata = load_and_validate_data(config)

    # Subset by condition
    condition_data = subset_conditions(adata, config)

    # Save
    save_results(condition_data, config)

    elapsed = datetime.now() - start_time
    _log("=" * 60)
    _log("M1 COMPLETE")
    _log("=" * 60)
    _log(f"Total time: {elapsed}")
    for key, adata_cond in condition_data.items():
        _log(f"  {key}: {adata_cond.n_obs:,} cells")
    _log(f"Output directory: {config.output_dir}")

    return condition_data


if __name__ == "__main__":
    main()
