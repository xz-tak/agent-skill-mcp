#!/usr/bin/env python3
"""
CELLxGENE-style cell type marker scoring for CyteType annotations.

Computes pairwise Cohen's d (pooled SD) between each CyteType-annotated cell type
and all others, then takes the bootstrapped 10th percentile as a specificity-weighted
marker score. Exports top markers per cell type as CSV + pkl.

This module is designed to be called by the CyteType skill (Phase 7) after annotation
is complete. It can also be run standalone via CLI.

Usage (standalone):
    python compute_markers.py \
        --input cytetype/integration_cytetype.h5ad \
        --output-dir cytetype/ \
        --top-n 200 \
        --bootstrap 1000

Usage (as module):
    from compute_markers import compute_cellxgene_markers
    markers = compute_cellxgene_markers(adata, output_dir="cytetype/")
"""

import argparse
import gc
import os
import pickle
import time

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_TOP_N = 200
DEFAULT_BOOTSTRAP_B = 1000
DEFAULT_PERCENTILE = 10
DEFAULT_SEED = 42


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _phase_timer(name):
    """Context manager that prints elapsed time for a phase."""
    class _Timer:
        def __enter__(self):
            self.t0 = time.time()
            print(f"\n{'='*60}")
            print(f"  {name}")
            print(f"{'='*60}")
            return self
        def __exit__(self, *_):
            elapsed = time.time() - self.t0
            m, s = divmod(elapsed, 60)
            print(f"  -> {name} done in {int(m)}m {s:.1f}s")
            self.elapsed = elapsed
    return _Timer()


def _compute_cluster_stats(X, obs_cluster, clusters):
    """Compute per-cluster mean, variance, and cell count from sparse X."""
    stats = {}
    for ct in tqdm(clusters, desc="Per-cluster stats"):
        mask = (obs_cluster == ct).values
        X_sub = X[mask]
        n = X_sub.shape[0]
        mean = np.asarray(X_sub.mean(axis=0)).ravel()
        var = np.asarray(X_sub.power(2).mean(axis=0)).ravel() - mean ** 2
        var = np.maximum(var, 0.0)
        stats[ct] = (mean, var, n)
    return stats


def _compute_markers_for_target(target, clusters, cluster_stats, rng,
                                n_genes, B, percentile):
    """Pairwise Cohen's d + bootstrapped P10 for one target cluster."""
    mean_T, var_T, _ = cluster_stats[target]
    others = [c for c in clusters if c != target]
    n_others = len(others)

    # Pairwise Cohen's d: (n_genes, n_others)
    effect_sizes = np.empty((n_genes, n_others), dtype=np.float32)
    for j, comp in enumerate(others):
        mean_C, var_C, _ = cluster_stats[comp]
        pooled_sd = np.sqrt((var_T + var_C) / 2.0)
        d = np.where(pooled_sd > 0, (mean_T - mean_C) / pooled_sd, 0.0)
        effect_sizes[:, j] = d

    # Bootstrap P10
    p10_accum = np.zeros(n_genes, dtype=np.float64)
    for _ in range(B):
        idx = rng.choice(n_others, size=n_others, replace=True)
        p10_accum += np.percentile(effect_sizes[:, idx], percentile, axis=1)
    marker_score = (p10_accum / B).astype(np.float32)

    # Specificity: fraction of pairwise d > 0
    specificity = (effect_sizes > 0).mean(axis=1).astype(np.float32)

    return marker_score, specificity


def _patch_h5ad_markers(h5ad_path, markers_df):
    """Patch uns['cytetype_marker'] directly via h5py without rewriting X.

    Using adata.write_h5ad() would promote float32→float64 and strip gzip
    compression, bloating a 23GB file to 139GB. This function opens the
    existing h5ad in append mode and only touches the uns group.
    """
    with h5py.File(h5ad_path, "a") as f:
        # Remove old marker data if present
        if "cytetype_marker" in f["uns"]:
            del f["uns"]["cytetype_marker"]

        mg = f["uns"].create_group("cytetype_marker")
        mg.attrs["encoding-type"] = "dataframe"
        mg.attrs["encoding-version"] = "0.2.0"
        mg.attrs["column-order"] = list(markers_df.columns)
        mg.attrs["_index"] = "_index"

        # Write index
        mg.create_dataset("_index", data=np.arange(len(markers_df)))
        mg["_index"].attrs["encoding-type"] = "array"
        mg["_index"].attrs["encoding-version"] = "0.2.0"

        # Write each column
        for col in markers_df.columns:
            vals = markers_df[col].values
            if vals.dtype == object or vals.dtype.kind == "U":
                str_list = [str(v) for v in vals]
                mg.create_dataset(col, data=str_list, dtype=h5py.string_dtype())
            else:
                mg.create_dataset(col, data=vals)
            mg[col].attrs["encoding-type"] = "array"
            mg[col].attrs["encoding-version"] = "0.2.0"

    print(f"  Patched uns['cytetype_marker']: {len(markers_df)} rows")


# ---------------------------------------------------------------------------
# Main API
# ---------------------------------------------------------------------------

def compute_cellxgene_markers(
    adata,
    output_dir=None,
    top_n=DEFAULT_TOP_N,
    bootstrap_b=DEFAULT_BOOTSTRAP_B,
    percentile=DEFAULT_PERCENTILE,
    seed=DEFAULT_SEED,
    update_h5ad=True,
    h5ad_path=None,
):
    """
    Compute CELLxGENE-style marker scores for CyteType-annotated cell types.

    Parameters
    ----------
    adata : AnnData
        Annotated AnnData with 'cytetype_cluster' and 'cytetype_leiden' in obs,
        and 'rank_genes_groups' with 'pts' in uns.
    output_dir : str or None
        Directory for CSV + pkl output. If None, skip file export.
    top_n : int
        Number of top marker genes per cell type (default: 200).
    bootstrap_b : int
        Number of bootstrap resamples (default: 1000).
    percentile : int
        Percentile of pairwise effect sizes (default: 10).
    seed : int
        Random seed for bootstrap reproducibility (default: 42).
    update_h5ad : bool
        Whether to update adata.uns['cytetype_marker'] and rewrite h5ad (default: True).
    h5ad_path : str or None
        Path to write updated h5ad. Required if update_h5ad=True.

    Returns
    -------
    dict
        Marker DataFrames keyed by cell type name.
    """
    t_total = time.time()

    X = adata.X
    gene_names = adata.var_names.tolist()
    n_genes = len(gene_names)
    obs_cluster = adata.obs["cytetype_cluster"]
    clusters = obs_cluster.cat.categories.tolist()
    print(f"  {len(clusters)} cell types, {n_genes} genes, top_n={top_n}, B={bootstrap_b}")

    # Reuse pts from rank_genes_groups
    pts_df = adata.uns["rank_genes_groups"]["pts"]

    # Build cluster -> leiden mapping
    cluster_to_leiden = (
        adata.obs[["cytetype_leiden", "cytetype_cluster"]]
        .drop_duplicates()
        .set_index("cytetype_cluster")["cytetype_leiden"]
        .to_dict()
    )

    # Phase 1: Per-cluster stats
    with _phase_timer("Precompute per-cluster mean/var"):
        cluster_stats = _compute_cluster_stats(X, obs_cluster, clusters)

    gc.collect()

    # Phase 2: Pairwise Cohen's d + bootstrap
    rng = np.random.default_rng(seed)
    markers = {}

    with _phase_timer(f"Pairwise Cohen's d + bootstrap P{percentile}"):
        for target in tqdm(clusters, desc="Marker scoring"):
            marker_score, specificity = _compute_markers_for_target(
                target, clusters, cluster_stats, rng,
                n_genes, bootstrap_b, percentile
            )

            leiden_id = str(cluster_to_leiden[target])
            mean_T = cluster_stats[target][0]

            if leiden_id in pts_df.columns:
                pct_cells = pts_df[leiden_id].values
            else:
                pct_cells = (mean_T > 0).astype(np.float32)

            df = pd.DataFrame({
                "cytetype_leiden": int(leiden_id),
                "cytetype_cluster": target,
                "Symbol": gene_names,
                "Effect Size": np.round(marker_score, 4),
                "Specificity": np.round(specificity, 4),
                "Mean Expression": np.round(mean_T, 4),
                "% of Cells": np.round(pct_cells * 100, 2),
            })
            df = (
                df.sort_values("Effect Size", ascending=False)
                .head(top_n)
                .reset_index(drop=True)
            )
            markers[target] = df

    # Phase 3: Export CSV + pkl
    all_markers = pd.concat(markers.values(), ignore_index=True)

    if output_dir is not None:
        with _phase_timer("Save CSV + pkl"):
            os.makedirs(output_dir, exist_ok=True)
            csv_path = os.path.join(output_dir, "cytetype_markers.csv")
            pkl_path = os.path.join(output_dir, "cytetype_markers.pkl")

            all_markers.to_csv(csv_path, index=False)
            print(f"  CSV: {csv_path} ({len(all_markers)} rows)")

            with open(pkl_path, "wb") as f:
                pickle.dump(markers, f)
            print(f"  PKL: {pkl_path}")

    # Phase 4: Update h5ad
    # IMPORTANT: Patch uns directly via h5py — NEVER use adata.write_h5ad().
    # anndata's read/write round-trip promotes float32→float64 and strips gzip
    # compression, bloating a 23GB file to 139GB. h5py patches only the uns group,
    # preserving X encoding and compression.
    #
    # Store as a single concatenated DataFrame, NOT dict-of-DataFrames.
    # Cell type names contain '/' (e.g., "ACTA2/TAGLN/...") which HDF5 interprets
    # as path separators, corrupting dict keys.
    # The pkl file retains the dict-of-DataFrames format (pickle handles '/' fine).
    if update_h5ad and h5ad_path is not None:
        with _phase_timer("Update h5ad (h5py patch)"):
            _patch_h5ad_markers(h5ad_path, all_markers)
            fsize = os.path.getsize(h5ad_path) / (1024 ** 3)
            print(f"  h5ad: {h5ad_path} ({fsize:.1f} GB)")
    elif update_h5ad:
        # Just update in-memory
        adata.uns["cytetype_marker"] = all_markers

    # Summary
    elapsed = time.time() - t_total
    m, s = divmod(elapsed, 60)
    print(f"\n{'='*60}")
    print(f"  DONE — {len(markers)} cell types, top {top_n} genes each")
    print(f"  Total time: {int(m)}m {s:.1f}s")
    print(f"{'='*60}")

    return markers


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="CELLxGENE-style marker scoring for CyteType annotations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python compute_markers.py --input cytetype/integration_cytetype.h5ad
    python compute_markers.py --input data.h5ad --top-n 500 --bootstrap 500
        """
    )
    parser.add_argument("--input", "-i", required=True, help="Path to annotated h5ad")
    parser.add_argument("--output-dir", "-o", default=None,
                        help="Output directory for CSV + pkl (default: same as input)")
    parser.add_argument("--top-n", type=int, default=DEFAULT_TOP_N,
                        help=f"Top genes per cell type (default: {DEFAULT_TOP_N})")
    parser.add_argument("--bootstrap", "-B", type=int, default=DEFAULT_BOOTSTRAP_B,
                        help=f"Bootstrap resamples (default: {DEFAULT_BOOTSTRAP_B})")
    parser.add_argument("--percentile", "-P", type=int, default=DEFAULT_PERCENTILE,
                        help=f"Effect size percentile (default: {DEFAULT_PERCENTILE})")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help=f"Random seed (default: {DEFAULT_SEED})")
    parser.add_argument("--no-update-h5ad", action="store_true",
                        help="Skip updating h5ad with markers")
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = args.output_dir or os.path.dirname(args.input) or "."

    with _phase_timer("Load h5ad"):
        adata = sc.read_h5ad(args.input)
        print(f"  Shape: {adata.shape}")
        print(f"  X dtype: {adata.X.dtype}, format: {type(adata.X).__name__}")

    compute_cellxgene_markers(
        adata,
        output_dir=output_dir,
        top_n=args.top_n,
        bootstrap_b=args.bootstrap,
        percentile=args.percentile,
        seed=args.seed,
        update_h5ad=not args.no_update_h5ad,
        h5ad_path=args.input if not args.no_update_h5ad else None,
    )


if __name__ == "__main__":
    main()
