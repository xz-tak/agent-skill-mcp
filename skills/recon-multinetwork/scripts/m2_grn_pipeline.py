#!/usr/bin/env python3
"""
M2: GRN Pipeline for ReCoN Multicellular Network Analysis

Combines scATAC preparation, CIRCE co-accessibility, and 5-layer GRN
construction into a single configurable module.

Sub-steps (--step):
    scatac  - Prepare scATAC object (filter, preprocess, convert peaks)
    circe   - Run CIRCE co-accessibility per cell type
    grn     - Build 5-layer GRN per cell type x condition
    all     - Run scatac -> circe -> grn sequentially

If config.scatac_path is None, scATAC/CIRCE steps are skipped and
RNA-only GRNs are built (layers 1 + 3 only).

Usage:
    # Full pipeline with config
    python m2_grn_pipeline.py --config config.json --step all

    # GRN only for one cell type
    python m2_grn_pipeline.py --config config.json --step grn --celltype Fibroblast

    # scATAC prep only
    python m2_grn_pipeline.py --config config.json --step scatac
"""

import argparse
import gc
import gzip
import json
import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from config import ReconConfig, get_config

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _log(msg: str) -> None:
    """Print timestamped message."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {msg}", flush=True)


def _celltype_key(name: str) -> str:
    """Normalise a cell type name to a filesystem-safe key."""
    return name.lower().replace(" ", "_")


# ---------------------------------------------------------------------------
# Sub-step 1: scATAC preparation
# ---------------------------------------------------------------------------

def _load_peak_files(peaks_dir: Path, peak_files: List[str]) -> set:
    """Load unique peaks from .bed.gz files."""
    all_peaks: list = []
    for peak_file in peak_files:
        peak_path = peaks_dir / peak_file
        if not peak_path.exists():
            _log(f"  WARNING: {peak_file} not found, skipping")
            continue
        try:
            peaks_df = pd.read_csv(
                peak_path, sep="\t", header=None, compression="gzip", usecols=[0, 1, 2]
            )
            peaks_df.columns = ["chr", "start", "end"]
            peak_names = peaks_df.apply(
                lambda x: f"{x['chr']}:{x['start']}-{x['end']}", axis=1
            ).tolist()
            all_peaks.extend(peak_names)
            _log(f"  {peak_file}: {len(peak_names):,} peaks")
        except Exception as e:
            _log(f"  ERROR loading {peak_file}: {e}")
    return set(all_peaks)


def _load_scatac_metadata(metadata_path: Path, scatac_celltype_mapping: Dict[str, List[str]]) -> pd.DataFrame:
    """Load cell metadata and filter to relevant cell types."""
    _log(f"Loading scATAC metadata from {metadata_path}")

    with gzip.open(str(metadata_path), "rt") as f:
        meta = pd.read_csv(f, sep="\t")
    _log(f"  Total cells in metadata: {len(meta):,}")

    # Collect all scATAC cell type names from the mapping values
    all_scatac_types: set = set()
    for types in scatac_celltype_mapping.values():
        all_scatac_types.update(types)

    # Filter to Adult cells matching mapped types
    if "Life stage" in meta.columns:
        adult_mask = meta["Life stage"] == "Adult"
    else:
        adult_mask = pd.Series(True, index=meta.index)

    ct_col = "cell type" if "cell type" in meta.columns else "cell_type"
    type_mask = meta[ct_col].isin(all_scatac_types)
    filtered = meta[adult_mask & type_mask].copy()

    _log(f"  Filtered to {len(filtered):,} cells matching mapping")
    return filtered


def _preprocess_scatac(scatac: ad.AnnData, config: ReconConfig) -> ad.AnnData:
    """Preprocess scATAC: binarize, filter, select variable features, normalize."""
    import episcanpy as epi

    _log("Preprocessing scATAC with episcanpy")

    max_val = np.max(scatac.X)
    if max_val > 1:
        _log("  Binarizing matrix...")
        epi.pp.binarize(scatac)

    # Remove empty features/barcodes
    epi.pp.filter_cells(scatac, min_features=1)
    epi.pp.filter_features(scatac, min_cells=1)

    # Add feature counts
    scatac.obs["nb_features"] = np.array((scatac.X > 0).sum(axis=1)).flatten()
    scatac.obs["log_nb_features"] = np.log10(scatac.obs["nb_features"] + 1)

    # Filter by min features and min cells
    _log(f"  Filtering: min_features={config.min_features_scatac}, min_cells={config.min_cells_scatac}")
    epi.pp.filter_cells(scatac, min_features=config.min_features_scatac)
    epi.pp.filter_features(scatac, min_cells=config.min_cells_scatac)
    _log(f"  After filtering: {scatac.n_obs:,} cells x {scatac.n_vars:,} peaks")

    # Calculate variability
    epi.pp.cal_var(scatac)

    # Save raw before normalization
    scatac.raw = scatac.copy()
    scatac.layers["binary"] = scatac.X.copy()

    # Select variable features
    if scatac.n_vars > config.nb_features_selected:
        _log(f"  Selecting top {config.nb_features_selected:,} variable features...")
        scatac = epi.pp.select_var_feature(
            scatac, nb_features=config.nb_features_selected, show=False, copy=True
        )

    # Normalize and log
    sc.pp.normalize_total(scatac)
    scatac.layers["normalized"] = scatac.X.copy()
    sc.pp.log1p(scatac)

    _log(f"  Final shape: {scatac.n_obs:,} cells x {scatac.n_vars:,} peaks")
    return scatac


def _convert_peak_names_for_circe(scatac: ad.AnnData) -> ad.AnnData:
    """Convert peak names from chr:start-end to chr_start_end for CIRCE."""
    new_names = [name.replace(":", "_").replace("-", "_") for name in scatac.var_names]
    scatac.var_names = new_names
    scatac.var_names_make_unique()
    return scatac


def prepare_scatac(config: ReconConfig) -> Optional[ad.AnnData]:
    """
    Prepare scATAC object: load, filter cells, preprocess, convert peaks.

    Args:
        config: Pipeline configuration with scatac_path, scatac_metadata_path,
                scatac_celltype_mapping, and preprocessing params.

    Returns:
        Preprocessed AnnData or None if scATAC is not configured.
    """
    if not config.scatac_path:
        _log("No scATAC path configured, skipping scATAC preparation")
        return None

    _log("=" * 60)
    _log("M2 / SCATAC PREP")
    _log("=" * 60)

    scatac_path = Path(config.scatac_path)
    if not scatac_path.exists():
        raise FileNotFoundError(f"scATAC file not found: {scatac_path}")

    grn_dir = config.get_grn_dir()

    # Check for cached result
    cached_path = grn_dir / "scatac_prep.h5ad"
    if cached_path.exists():
        _log(f"Loading cached scATAC from {cached_path}")
        return sc.read_h5ad(str(cached_path))

    # Load full scATAC
    _log(f"Loading scATAC from {scatac_path}")
    scatac_full = ad.read_h5ad(str(scatac_path))
    _log(f"  Full matrix shape: {scatac_full.shape}")

    # Filter cells using metadata if provided
    if config.scatac_metadata_path and config.scatac_celltype_mapping:
        meta_path = Path(config.scatac_metadata_path)
        if not meta_path.exists():
            raise FileNotFoundError(f"scATAC metadata not found: {meta_path}")

        meta = _load_scatac_metadata(meta_path, config.scatac_celltype_mapping)

        # Find common cells
        cell_id_col = "cellID" if "cellID" in meta.columns else meta.columns[0]
        common_cells = sorted(set(scatac_full.obs_names) & set(meta[cell_id_col]))
        _log(f"  Common cells: {len(common_cells):,}")

        if len(common_cells) == 0:
            _log("ERROR: No common cells between matrix and metadata")
            return None

        scatac = scatac_full[scatac_full.obs_names.isin(common_cells), :].copy()

        # Attach metadata
        meta_indexed = meta.set_index(cell_id_col).loc[scatac.obs_names]
        ct_col = "cell type" if "cell type" in meta_indexed.columns else "cell_type"
        scatac.obs["cell_type"] = meta_indexed[ct_col].values
    else:
        _log("  No metadata/mapping provided, using all cells")
        scatac = scatac_full.copy()

    del scatac_full
    gc.collect()

    # Ensure float64
    scatac.X = scatac.X.astype("float64")

    # Preprocess
    scatac = _preprocess_scatac(scatac, config)

    # Convert peak names for CIRCE
    scatac = _convert_peak_names_for_circe(scatac)

    # Add chromosome annotations
    _add_chromosome_annotations(scatac)

    # Save
    _log(f"Saving prepared scATAC to {cached_path}")
    scatac.write(cached_path)

    # Save peak list
    peaks_df = pd.DataFrame({"peak_id": scatac.var_names})
    peaks_df.to_csv(grn_dir / "peaks.csv", index=False)

    _log(f"scATAC prep complete: {scatac.shape}")
    return scatac


# ---------------------------------------------------------------------------
# Sub-step 2: CIRCE co-accessibility
# ---------------------------------------------------------------------------

def _add_chromosome_annotations(scatac: ad.AnnData) -> None:
    """Add chromosome/start/end columns to var from peak names."""
    if "chromosome" in scatac.var.columns:
        return
    chromosomes, starts, ends = [], [], []
    for name in scatac.var_names:
        if ":" in name:
            chrom, coords = name.split(":")
            start, end = coords.split("-")
        else:
            parts = name.split("_")
            chrom, start, end = parts[0], parts[1], parts[2]
        chromosomes.append(chrom)
        starts.append(int(start))
        ends.append(int(end))
    scatac.var["chromosome"] = chromosomes
    scatac.var["start"] = starts
    scatac.var["end"] = ends


def _get_cells_for_l2(
    scatac: ad.AnnData,
    l2_celltype: str,
    mapping: Dict[str, List[str]],
    min_cells: int,
) -> Optional[ad.AnnData]:
    """Get scATAC cells matching an L2 cell type via the mapping."""
    scatac_types = mapping.get(l2_celltype, [])
    if not scatac_types:
        return None

    # Build reverse map and assign L2 labels
    if "l2_celltype" not in scatac.obs.columns:
        reverse_map = {}
        for l2, types in mapping.items():
            for t in types:
                reverse_map[t] = l2
        scatac.obs["l2_celltype"] = scatac.obs.get("cell_type", pd.Series(dtype=str)).map(reverse_map)

    mask = scatac.obs["l2_celltype"] == l2_celltype
    n_cells = mask.sum()
    if n_cells < min_cells:
        _log(f"  {l2_celltype}: only {n_cells} cells (need {min_cells}), skipping")
        return None

    subset = scatac[mask].copy()
    _log(f"  {l2_celltype}: {n_cells} cells")
    return subset


def _run_circe_single(
    scatac_subset: ad.AnnData, celltype: str, window_size: int, n_cpus: int
) -> Optional[pd.DataFrame]:
    """Run CIRCE co-accessibility on a cell subset."""
    import circe as ci

    _log(f"  Running CIRCE (window={window_size / 1000:.0f}kb, cells={scatac_subset.n_obs})...")
    start_time = datetime.now()

    try:
        if hasattr(scatac_subset.X, "toarray"):
            scatac_subset.X = scatac_subset.X.toarray()

        ci.compute_atac_network(
            scatac_subset,
            window_size=window_size,
            unit_distance=1000,
            distance_constraint=None,
            njobs=n_cpus,
            verbose=1,
        )

        connections = ci.extract_atac_links(
            scatac_subset,
            key="atac_network",
            columns=("peak1", "peak2", "coaccess"),
        )

        elapsed = datetime.now() - start_time
        _log(f"  CIRCE complete in {elapsed}: {len(connections):,} connections")
        return connections

    except Exception as e:
        _log(f"  CIRCE error for {celltype}: {e}")
        return None


def run_circe_per_celltype(
    config: ReconConfig, scatac: Optional[ad.AnnData] = None
) -> Dict[str, pd.DataFrame]:
    """
    Run CIRCE co-accessibility for each cell type defined in the mapping.

    Args:
        config: Pipeline configuration.
        scatac: Pre-loaded scATAC AnnData. Loaded from cache if None.

    Returns:
        Dict mapping cell type to CIRCE connections DataFrame.
    """
    if not config.scatac_path:
        _log("No scATAC configured, skipping CIRCE")
        return {}

    _log("=" * 60)
    _log("M2 / CIRCE CO-ACCESSIBILITY")
    _log("=" * 60)

    grn_dir = config.get_grn_dir()

    # Load scATAC if not provided
    if scatac is None:
        cached_path = grn_dir / "scatac_prep.h5ad"
        if not cached_path.exists():
            _log("ERROR: scATAC not prepared. Run --step scatac first.")
            return {}
        scatac = sc.read_h5ad(str(cached_path))

    _add_chromosome_annotations(scatac)

    # Determine cell types from the mapping
    mapping = config.scatac_celltype_mapping
    if not mapping:
        _log("WARNING: scatac_celltype_mapping is empty, cannot run CIRCE")
        return {}

    results: Dict[str, pd.DataFrame] = {}
    summary: Dict[str, dict] = {}

    for l2_celltype in mapping:
        _log(f"\nProcessing: {l2_celltype}")
        ct_key = _celltype_key(l2_celltype)
        output_path = grn_dir / f"circe_{ct_key}.csv"

        # Check cache
        if output_path.exists():
            _log(f"  Loading cached CIRCE: {output_path}")
            df = pd.read_csv(output_path)
            results[l2_celltype] = df
            summary[l2_celltype] = {"status": "cached", "n_connections": len(df)}
            continue

        subset = _get_cells_for_l2(scatac, l2_celltype, mapping, config.min_cells_scatac)
        if subset is None:
            summary[l2_celltype] = {"status": "skipped_insufficient_cells"}
            continue

        connections = _run_circe_single(subset, l2_celltype, config.circe_window, config.n_cpus)
        del subset
        gc.collect()

        if connections is not None and len(connections) > 0:
            connections.to_csv(output_path, index=False)
            results[l2_celltype] = connections
            summary[l2_celltype] = {"status": "complete", "n_connections": len(connections)}
        else:
            summary[l2_celltype] = {"status": "failed"}

    # Save summary
    summary_path = grn_dir / "circe_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    _log(f"\nCIRCE summary saved: {summary_path}")

    return results


# ---------------------------------------------------------------------------
# Sub-step 3: 5-layer GRN construction
# ---------------------------------------------------------------------------

def _load_tfs() -> List[str]:
    """Load human TF list from hummuspy."""
    from hummuspy.loader import load_tfs as hummus_load_tfs

    tfs = hummus_load_tfs("human_tfs_r_hummus")
    _log(f"Loaded {len(tfs)} TFs")
    return tfs


def _load_scrna_condition(config: ReconConfig, condition: str) -> ad.AnnData:
    """Load scRNA data for a specific condition."""
    rna_path = Path(config.output_dir) / f"adata_{condition}.h5ad"
    if not rna_path.exists():
        raise FileNotFoundError(
            f"scRNA for condition '{condition}' not found: {rna_path}. "
            f"Run m1_data_prep.py first."
        )
    _log(f"Loading scRNA ({condition}) from {rna_path}")
    adata = sc.read_h5ad(str(rna_path))
    _log(f"  Shape: {adata.shape}")
    return adata


def _subset_scrna(
    adata: ad.AnnData, celltype: str, celltype_col: str, min_cells: int
) -> Optional[ad.AnnData]:
    """Subset scRNA to a specific cell type."""
    mask = adata.obs[celltype_col] == celltype
    n_cells = mask.sum()
    if n_cells < min_cells:
        _log(f"  Insufficient cells for {celltype}: {n_cells} < {min_cells}")
        return None
    _log(f"  Subset {celltype}: {n_cells} cells")
    return adata[mask].copy()


def _load_circe_network(grn_dir: Path, celltype: str) -> Optional[pd.DataFrame]:
    """Load CIRCE co-accessibility network for a cell type."""
    ct_key = _celltype_key(celltype)
    circe_path = grn_dir / f"circe_{ct_key}.csv"
    if not circe_path.exists():
        return None
    circe = pd.read_csv(circe_path)
    _log(f"  Loaded CIRCE for {celltype}: {len(circe):,} connections")
    return circe


def _compute_layer1_rna(
    rna_subset: ad.AnnData,
    tfs: List[str],
    grn_dir: Path,
    celltype: str,
    condition: str,
    n_cpus: int,
) -> pd.DataFrame:
    """Layer 1: TF -> Gene network via GRNBoost2."""
    from recon.infer_grn import compute_rna_network

    ct_key = _celltype_key(celltype)
    output_path = grn_dir / f"{ct_key}_{condition}_rna_network.csv"
    if output_path.exists():
        _log(f"  Layer 1: Loading cached RNA network")
        return pd.read_csv(output_path)

    _log(f"  Layer 1: Computing RNA network (GRNBoost2)...")
    available_tfs = [tf for tf in tfs if tf in rna_subset.var_names]
    _log(f"    TFs: {len(available_tfs)}/{len(tfs)} available")

    rna_network = compute_rna_network(
        df_exp_mtx=rna_subset, tf_names=available_tfs, method="GBM", n_cpu=n_cpus, seed=42
    )
    _log(f"    Edges: {len(rna_network):,}")
    rna_network.to_csv(output_path, index=False)
    return rna_network


def _compute_layer3_tf(
    rna_subset: ad.AnnData,
    tfs: List[str],
    grn_dir: Path,
    celltype: str,
    condition: str,
) -> pd.DataFrame:
    """Layer 3: TF <-> TF correlation network."""
    from recon.infer_grn import compute_tf_network

    ct_key = _celltype_key(celltype)
    output_path = grn_dir / f"{ct_key}_{condition}_tf_network.csv"
    if output_path.exists():
        _log(f"  Layer 3: Loading cached TF network")
        return pd.read_csv(output_path)

    _log(f"  Layer 3: Computing TF correlation network...")
    available_tfs = [tf for tf in tfs if tf in rna_subset.var_names]
    tf_network = compute_tf_network(rna=rna_subset, tfs_list=available_tfs)
    _log(f"    TF network edges: {len(tf_network):,}")
    tf_network.to_csv(output_path, index=False)
    return tf_network


def _compute_layer4_tf_atac(
    scatac: ad.AnnData,
    tfs: List[str],
    grn_dir: Path,
    celltype: str,
    config: ReconConfig,
) -> pd.DataFrame:
    """Layer 4: TF -> Peak links via motif scanning (shared across conditions)."""
    from recon.infer_grn import compute_tf_to_atac_links

    ct_key = _celltype_key(celltype)
    output_path = grn_dir / f"{ct_key}_tf_atac_links.csv"
    if output_path.exists():
        _log(f"  Layer 4: Loading cached TF-ATAC links")
        return pd.read_csv(output_path)

    _log(f"  Layer 4: Computing TF -> Peak links (motif scanning)...")
    tf_atac_links = compute_tf_to_atac_links(
        atac=scatac,
        ref_genome=config.ref_genome,
        tfs_list=tfs,
        fpr=config.motif_fpr,
        n_cpus=config.n_cpus,
        verbose=True,
    )
    _log(f"    TF-Peak links: {len(tf_atac_links):,}")
    tf_atac_links.to_csv(output_path, index=False)
    return tf_atac_links


def _compute_layer5_atac_rna(
    scatac: ad.AnnData,
    rna_subset: ad.AnnData,
    grn_dir: Path,
    celltype: str,
    config: ReconConfig,
) -> pd.DataFrame:
    """Layer 5: Peak -> Gene links via TSS proximity (shared across conditions)."""
    from recon.infer_grn import compute_atac_to_rna_links

    ct_key = _celltype_key(celltype)
    output_path = grn_dir / f"{ct_key}_atac_rna_links.csv"
    if output_path.exists():
        _log(f"  Layer 5: Loading cached ATAC-RNA links")
        return pd.read_csv(output_path)

    _log(f"  Layer 5: Computing Peak -> Gene links (TSS proximity)...")
    atac_rna_links = compute_atac_to_rna_links(
        atac=scatac, rna=rna_subset, ref_genome=config.ref_genome
    )
    _log(f"    Peak-Gene links: {len(atac_rna_links):,}")
    atac_rna_links.to_csv(output_path, index=False)
    return atac_rna_links


def _integrate_layers(
    rna_network: pd.DataFrame,
    atac_network: Optional[pd.DataFrame],
    tf_network: pd.DataFrame,
    tf_atac_links: Optional[pd.DataFrame],
    atac_rna_links: Optional[pd.DataFrame],
    grn_dir: Path,
    celltype: str,
    condition: str,
    n_cpus: int,
) -> pd.DataFrame:
    """Integrate all available layers into a final GRN."""
    from recon.infer_grn import generate_grn

    ct_key = _celltype_key(celltype)
    output_path = grn_dir / f"{ct_key}_{condition}_5layer_grn.csv"
    if output_path.exists():
        _log(f"  Integration: Loading cached GRN")
        return pd.read_csv(output_path)

    _log(f"  Integrating layers...")

    # Prepare ATAC network columns if present
    if atac_network is not None and "peak1" in atac_network.columns:
        atac_network = atac_network.rename(
            columns={"peak1": "source", "peak2": "target", "coaccess": "weight"}
        )

    # Use empty DataFrames for missing ATAC layers
    if atac_network is None:
        atac_network = pd.DataFrame(columns=["source", "target", "weight"])
    if tf_atac_links is None:
        tf_atac_links = pd.DataFrame()
    if atac_rna_links is None:
        atac_rna_links = pd.DataFrame()

    # RNA-only guard: if all ATAC layers are empty, skip generate_grn integration
    # and save the RNA network directly as the GRN (layers 1+3 only)
    atac_all_empty = (
        atac_network.empty
        and (tf_atac_links.empty if isinstance(tf_atac_links, pd.DataFrame) else True)
        and (atac_rna_links.empty if isinstance(atac_rna_links, pd.DataFrame) else True)
    )
    if atac_all_empty:
        _log("    RNA-only mode: skipping 5-layer integration (ATAC layers empty)")
        _log("    Saving RNA network (layer 1) + TF network (layer 3) as GRN")
        grn = pd.concat([rna_network, tf_network], ignore_index=True) if tf_network is not None else rna_network
        grn.to_csv(output_path, index=False)
        return grn

    grn = generate_grn(
        rna_network=rna_network,
        atac_network=atac_network,
        tf_network=tf_network,
        tf_to_atac_links=tf_atac_links,
        atac_to_rna_links=atac_rna_links,
        n_jobs=n_cpus,
    )

    _log(f"    Final GRN: {len(grn):,} edges")
    grn.to_csv(output_path, index=False)
    return grn


def _build_grn_single(
    celltype: str,
    condition: str,
    tfs: List[str],
    scatac: Optional[ad.AnnData],
    config: ReconConfig,
) -> Dict:
    """Build complete GRN for one cell type x condition."""
    ct_key = _celltype_key(celltype)
    grn_dir = config.get_grn_dir()
    has_atac = scatac is not None

    # Check cache
    final_path = grn_dir / f"{ct_key}_{condition}_5layer_grn.csv"
    if final_path.exists():
        _log(f"  GRN already exists: {ct_key} x {condition}")
        return {"celltype": celltype, "condition": condition, "status": "cached"}

    _log(f"\nBuilding GRN: {celltype} x {condition}")
    start_time = datetime.now()
    stats: Dict = {"celltype": celltype, "condition": condition}

    try:
        # Load and subset scRNA
        adata = _load_scrna_condition(config, condition)
        rna_subset = _subset_scrna(adata, celltype, config.celltype_col, config.min_cells_grn)
        if rna_subset is None:
            return {**stats, "status": "skipped_insufficient_cells"}

        stats["n_cells"] = rna_subset.n_obs

        # Layer 1: RNA network
        rna_network = _compute_layer1_rna(rna_subset, tfs, grn_dir, celltype, condition, config.n_cpus)
        stats["layer1_edges"] = len(rna_network)

        # Layer 3: TF network
        tf_network = _compute_layer3_tf(rna_subset, tfs, grn_dir, celltype, condition)
        stats["layer3_edges"] = len(tf_network)

        # ATAC-dependent layers (2, 4, 5)
        atac_network = None
        tf_atac_links = None
        atac_rna_links = None

        if has_atac:
            atac_network = _load_circe_network(grn_dir, celltype)
            if atac_network is not None:
                stats["layer2_edges"] = len(atac_network)

            tf_atac_links = _compute_layer4_tf_atac(scatac, tfs, grn_dir, celltype, config)
            stats["layer4_edges"] = len(tf_atac_links)

            atac_rna_links = _compute_layer5_atac_rna(scatac, rna_subset, grn_dir, celltype, config)
            stats["layer5_edges"] = len(atac_rna_links)

        # Integrate
        final_grn = _integrate_layers(
            rna_network, atac_network, tf_network, tf_atac_links, atac_rna_links,
            grn_dir, celltype, condition, config.n_cpus,
        )
        stats["final_edges"] = len(final_grn)

        elapsed = datetime.now() - start_time
        stats["status"] = "complete"
        stats["time_seconds"] = elapsed.total_seconds()
        _log(f"  Complete: {len(final_grn):,} edges in {elapsed}")

        del adata, rna_subset
        gc.collect()
        return stats

    except Exception as e:
        _log(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()
        return {**stats, "status": "error", "error": str(e)}


def build_grn(
    config: ReconConfig, celltype: Optional[str] = None
) -> Dict:
    """
    Build 5-layer GRNs for cell type x condition combinations.

    Args:
        config: Pipeline configuration.
        celltype: If provided, build only for this cell type. Otherwise all.

    Returns:
        Dict with build statistics per cell type x condition.
    """
    _log("=" * 60)
    _log("M2 / 5-LAYER GRN CONSTRUCTION")
    _log("=" * 60)

    grn_dir = config.get_grn_dir()

    # Load TFs
    tfs = _load_tfs()

    # Load scATAC if available
    scatac = None
    if config.scatac_path:
        cached_path = grn_dir / "scatac_prep.h5ad"
        if cached_path.exists():
            _log(f"Loading scATAC from {cached_path}")
            scatac = sc.read_h5ad(str(cached_path))
        else:
            _log("WARNING: scATAC configured but not prepared. Building RNA-only GRNs.")

    # Determine cell types
    if celltype:
        cell_types = [celltype]
    else:
        # Get cell types from M1 output metadata
        meta_path = Path(config.output_dir) / "data_prep_metadata.json"
        if meta_path.exists():
            with open(meta_path) as f:
                meta = json.load(f)
            cell_types = meta.get("celltypes", [])
        else:
            _log("WARNING: No data_prep_metadata.json found. Provide --celltype.")
            return {}

    # Determine conditions
    conditions = [c.lower().replace(" ", "_") for c in config.conditions]

    _log(f"Cell types: {cell_types}")
    _log(f"Conditions: {conditions}")

    all_stats = []
    for ct in cell_types:
        for condition in conditions:
            stats = _build_grn_single(ct, condition, tfs, scatac, config)
            all_stats.append(stats)
            gc.collect()

    # Save summary
    summary_path = grn_dir / "grn_build_summary.json"
    with open(summary_path, "w") as f:
        json.dump(all_stats, f, indent=2)
    _log(f"\nGRN build summary saved: {summary_path}")

    return {"stats": all_stats}


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main(config: Optional[ReconConfig] = None, step: str = "all") -> None:
    """
    Run the GRN pipeline.

    Args:
        config: Optional ReconConfig. Falls back to CLI args if None.
        step: Which sub-step to run: scatac, circe, grn, all.
    """
    # Parse CLI args if no config provided
    if config is None:
        parser = argparse.ArgumentParser(
            description="M2: GRN Pipeline",
            parents=[_build_m2_parser()],
            add_help=False,
        )
        args = parser.parse_args()
        step = args.step
        celltype = getattr(args, "celltype", None)
        config = get_config()
    else:
        celltype = None

    _log("=" * 60)
    _log("M2: GRN PIPELINE")
    _log(f"Step: {step}")
    _log("=" * 60)
    start_time = datetime.now()

    if step in ("scatac", "all"):
        prepare_scatac(config)

    if step in ("circe", "all"):
        run_circe_per_celltype(config)

    if step in ("grn", "all"):
        ct = celltype if step == "grn" else None
        build_grn(config, celltype=ct)

    elapsed = datetime.now() - start_time
    _log("=" * 60)
    _log(f"M2 COMPLETE (step={step}, elapsed={elapsed})")
    _log("=" * 60)


def _build_m2_parser() -> argparse.ArgumentParser:
    """Build M2-specific CLI argument parser."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--step",
        choices=["scatac", "circe", "grn", "all"],
        default="all",
        help="Which sub-step to run (default: all)",
    )
    parser.add_argument(
        "--celltype",
        type=str,
        default=None,
        help="Process single cell type (for --step grn only)",
    )
    parser.add_argument("--config", type=str, help="Path to JSON config file")
    return parser


if __name__ == "__main__":
    main()
