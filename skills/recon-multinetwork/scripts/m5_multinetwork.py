#!/usr/bin/env python3
"""
Module 5: Multinetwork Generation

Merges CCC + GRN + Receptor-Gene networks into a unified multinetwork
for ReCoN analysis. Creates expression-weighted bidirectional receptor-gene
edges with optional module annotations.

Output: Parquet files with standardized node format (GENE:CellType)
"""

import ast
import gc
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

from recon.data import load_receptor_genes

from config import ReconConfig, get_config

warnings.filterwarnings("ignore")


def log(msg: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def load_ccc(condition: str, config: ReconConfig) -> pd.DataFrame:
    """
    Load CCC data based on config.ccc_source, format nodes as gene:celltype.

    Score handling:
    - merged: Uses OVERALL percentile-rank (0-1) as lr_means
    - cellphonedb: Uses RAW lr_means (original score)
    - cellchat: Uses RAW prob as lr_means (original score)
    """
    ccc_dir = config.get_ccc_dir()
    ccc_source = config.ccc_source

    if ccc_source == "cellchat":
        path = ccc_dir / "cellchat" / f"{condition}_ccc.csv"
        log(f"  Loading CCC (CellChat) from {path}")
        df = pd.read_csv(path)
        log(f"    Raw edges: {len(df):,}")
        df = df.rename(columns={
            "source": "celltype_source",
            "target": "celltype_target",
            "ligand": "source",
            "receptor": "target",
        })
        # Handle alternative CellChat column names (ligand_name/receptor_name)
        if "ligand_name" in df.columns:
            df = df.rename(columns={"ligand_name": "source", "receptor_name": "target"})
        df["lr_means"] = df["prob"]
    elif ccc_source == "cellphonedb":
        path = ccc_dir / "cellphonedb" / f"{condition}_ccc.csv"
        log(f"  Loading CCC (CellPhoneDB) from {path}")
        df = pd.read_csv(path)
        log(f"    Raw edges: {len(df):,}")
    else:  # merged
        path = ccc_dir / "merged" / f"{condition}_ccc.csv"
        log(f"  Loading CCC (merged) from {path}")
        df = pd.read_csv(path)
        log(f"    Raw edges: {len(df):,}")

    # Format nodes as gene:celltype
    df["source_node"] = df["source"] + ":" + df["celltype_source"]
    df["target_node"] = df["target"] + ":" + df["celltype_target"]

    # Weight: Apply percentile-rank normalization to 0-1 range
    if ccc_source == "merged" and "pct_rank_overall" in df.columns:
        df["weight"] = df["pct_rank_overall"].fillna(0)
    elif ccc_source == "cellphonedb":
        df["weight"] = df["lr_means"].rank(pct=True).fillna(0)
    elif ccc_source == "cellchat":
        df["weight"] = df["prob"].rank(pct=True).fillna(0)
    else:
        df["weight"] = 1 - df["pval"].fillna(1.0)

    df["edge_type"] = "ccc"
    df["interaction"] = "ligand-receptor"

    result = df[["source_node", "target_node", "weight", "edge_type", "interaction"]].copy()
    log(f"    Final edges: {len(result):,}")
    return result


def load_grn(condition: str, config: ReconConfig) -> pd.DataFrame:
    """
    Load 5-layer GRN for all cell types, format nodes as gene:celltype.

    Auto-detects cell types from available GRN files.
    Filters edges by config.grn_score_threshold.
    """
    grn_dir = config.get_grn_dir()
    score_threshold = config.grn_score_threshold

    # Auto-detect cell types from available GRN files for this condition
    pattern = f"*_{condition}_5layer_grn.csv"
    grn_files = sorted(grn_dir.glob(pattern))
    log(f"  Loading GRN: found {len(grn_files)} files matching {pattern}")

    all_edges = []

    for path in grn_files:
        # Extract cell type from filename: {ct_safe}_{condition}_5layer_grn.csv
        ct_safe = path.stem.replace(f"_{condition}_5layer_grn", "")
        df = pd.read_csv(path)
        original_count = len(df)

        if score_threshold > 0:
            df = df[df["score"] > score_threshold]

        if df.empty:
            log(f"    {ct_safe}: {original_count:,} -> 0 edges (filtered)")
            continue

        # Reconstruct proper-case cell type name from filename
        ct = ct_safe.replace("_", " ").title()
        # Fix known single-letter/acronym cell types
        _ct_fixes = {"Nk": "NK", "Smooth Muscle": "Smooth muscle"}
        ct = _ct_fixes.get(ct, ct)

        # Clean gene names (remove _TF suffix from seed and target)
        df["source_gene"] = df["seed"].str.replace("_TF", "", regex=False)
        df["target_gene"] = (
            df["target"]
            .str.replace("_TF", "", regex=False)
            .str.replace("_target", "", regex=False)
        )

        # Format nodes as gene:celltype
        df["source_node"] = df["source_gene"] + ":" + ct
        df["target_node"] = df["target_gene"] + ":" + ct

        # Normalize weight to 0-1
        max_score = df["score"].max() if len(df) > 0 else 0.7
        df["weight"] = df["score"] / max(max_score, 0.7)
        df["edge_type"] = "grn"
        df["interaction"] = "tf-target"

        edges = df[["source_node", "target_node", "weight", "edge_type", "interaction"]].copy()
        all_edges.append(edges)
        log(f"    {ct_safe}: {original_count:,} -> {len(edges):,} edges")

        del df
        gc.collect()

    if not all_edges:
        return pd.DataFrame(columns=["source_node", "target_node", "weight", "edge_type", "interaction"])

    result = pd.concat(all_edges, ignore_index=True)
    log(f"    Total GRN edges: {len(result):,}")
    return result


def load_receptor_gene_base() -> pd.DataFrame:
    """Load raw NicheNet receptor-gene network (without cell type expansion)."""
    log("  Loading receptor-gene base network from NicheNet PKN")
    rg = load_receptor_genes("human_receptor_gene_from_NichenetPKN")
    log(f"    Raw edges: {len(rg):,}")

    rg["pkn_weight"] = rg["weight"] / rg["weight"].max()
    return rg[["source", "target", "pkn_weight"]]


def build_receptor_gene_edges(
    rg_base: pd.DataFrame,
    expr_dict: Dict[Tuple[str, str], float],
    cell_types: List[str],
) -> pd.DataFrame:
    """
    Build condition-specific receptor-gene edges with expression weighting.

    Creates bidirectional edges (receptor<->gene) with weights based on:
    weight = pkn_weight * sqrt(source_expression * target_expression)
    """
    log("  Building receptor-gene edges (expression-weighted, bidirectional)")
    all_edges = []

    for ct in cell_types:
        df = rg_base.copy()
        df["source_node"] = df["source"] + ":" + ct
        df["target_node"] = df["target"] + ":" + ct

        df["source_expr"] = df["source"].apply(lambda g: expr_dict.get((g, ct), 0))
        df["target_expr"] = df["target"].apply(lambda g: expr_dict.get((g, ct), 0))

        df["raw_weight"] = df["pkn_weight"] * np.sqrt(df["source_expr"] * df["target_expr"])

        df["edge_type"] = "receptor_gene"
        df["interaction"] = "intra-intra"

        # Forward edges (receptor -> gene)
        forward = df[["source_node", "target_node", "raw_weight", "edge_type", "interaction"]].copy()
        forward.columns = ["source_node", "target_node", "weight", "edge_type", "interaction"]

        # Reverse edges (gene -> receptor, bidirectional)
        reverse = forward.copy()
        reverse["source_node"], reverse["target_node"] = (
            forward["target_node"].values.copy(),
            forward["source_node"].values.copy(),
        )

        all_edges.extend([forward, reverse])
        log(f"    {ct}: {len(forward):,} forward + {len(reverse):,} reverse = {len(forward) * 2:,} edges")

    if not all_edges:
        return pd.DataFrame(columns=["source_node", "target_node", "weight", "edge_type", "interaction"])

    result = pd.concat(all_edges, ignore_index=True)

    if result["weight"].max() > 0:
        result["weight"] = result["weight"] / result["weight"].max()

    log(f"    Total receptor-gene edges (bidirectional): {len(result):,}")
    return result


def load_expression_data(condition: str, config: ReconConfig) -> Dict[Tuple[str, str], float]:
    """
    Compute mean expression per gene x cell type for a condition.

    Uses per-condition h5ad files: {output_dir}/adata_{condition}.h5ad
    """
    log(f"  Loading expression data for {condition}")

    expr_path = Path(config.output_dir) / f"adata_{condition}.h5ad"
    if not expr_path.exists():
        log(f"    Expression file not found: {expr_path}")
        return {}

    log(f"    Loading {expr_path.name}...")
    adata = sc.read_h5ad(expr_path)
    log(f"    Shape: {adata.shape}")

    celltype_col = config.celltype_col
    gene_names = adata.var_names.tolist()
    expr_dict: Dict[Tuple[str, str], float] = {}

    # Auto-detect cell types from data
    available_cts = adata.obs[celltype_col].unique().tolist()

    for ct in available_cts:
        ct_mask = adata.obs[celltype_col] == ct
        n_cells = ct_mask.sum()

        if n_cells == 0:
            continue

        ct_indices = np.where(ct_mask.values)[0]

        # Compute mean expression in chunks
        chunk_size = 10000
        mean_expr = np.zeros(adata.n_vars)

        for start in range(0, len(ct_indices), chunk_size):
            end = min(start + chunk_size, len(ct_indices))
            chunk_indices = ct_indices[start:end]

            if hasattr(adata.X, "toarray"):
                X_chunk = adata.X[chunk_indices, :].toarray()
            else:
                X_chunk = adata.X[chunk_indices, :]

            chunk_sum = X_chunk.sum(axis=0)
            mean_expr += chunk_sum.flatten() if hasattr(chunk_sum, "flatten") else chunk_sum
            del X_chunk
            gc.collect()

        mean_expr = mean_expr / n_cells

        for i, gene in enumerate(gene_names):
            expr_dict[(gene, ct)] = float(mean_expr[i])

        log(f"    {ct}: {n_cells:,} cells")

    del adata
    gc.collect()

    log(f"    Expression dict entries: {len(expr_dict):,}")
    return expr_dict


def load_module_weights(config: ReconConfig) -> Dict[str, str]:
    """
    Load module weights, return gene -> module_index mapping.

    Uses config.module_file; returns empty dict if not configured.
    """
    if config.module_file is None:
        log("  Module file not configured, skipping module annotations")
        return {}

    module_path = Path(config.module_file)
    if not module_path.exists():
        log(f"  Module file not found: {module_path}")
        return {}

    log(f"  Loading module weights from {module_path}")
    df = pd.read_csv(module_path, sep="\t")

    gene_to_module: Dict[str, str] = {}
    for idx, row in df.iterrows():
        genes = ast.literal_eval(row["sets"])
        for gene in genes:
            gene_to_module[gene] = str(idx)

    log(f"    Genes in modules: {len(gene_to_module):,}")
    log(f"    Modules: {len(df)}")
    return gene_to_module


def _parse_node(node: str) -> Tuple[str, Optional[str]]:
    """Parse gene:celltype format into (gene, celltype) tuple."""
    parts = node.rsplit(":", 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return node, None


def _detect_cell_types(config: ReconConfig, condition: str) -> List[str]:
    """Auto-detect cell types from available GRN files for a condition."""
    grn_dir = config.get_grn_dir()
    pattern = f"*_{condition}_5layer_grn.csv"
    grn_files = sorted(grn_dir.glob(pattern))
    cell_types = []
    for path in grn_files:
        ct_safe = path.stem.replace(f"_{condition}_5layer_grn", "")
        cell_types.append(ct_safe)
    return cell_types


def generate_multinetwork(
    condition: str,
    expr_dict: Dict[Tuple[str, str], float],
    gene_to_module: Dict[str, str],
    rg_base: pd.DataFrame,
    config: ReconConfig,
) -> pd.DataFrame:
    """
    Merge all network layers with expression and module annotations.

    Produces a unified multinetwork for a single condition.
    """
    log(f"\n{'=' * 60}")
    log(f"GENERATING MULTINETWORK FOR {condition.upper()}")
    log(f"{'=' * 60}")

    # Auto-detect cell types
    cell_types = _detect_cell_types(config, condition)
    if not cell_types:
        # Fallback: try to detect from expression data keys
        cell_types = sorted({ct for (_, ct) in expr_dict.keys()}) if expr_dict else []
    if not cell_types:
        log(f"  WARNING: No cell types detected for {condition}. Check GRN files or expression data.")
        log(f"  Returning empty multinetwork.")
        return pd.DataFrame()
    log(f"  Cell types: {cell_types}")

    # Load all network layers
    ccc = load_ccc(condition, config)
    grn = load_grn(condition, config)

    # Build expression-weighted bidirectional receptor-gene edges
    rg = build_receptor_gene_edges(rg_base, expr_dict, cell_types)

    # Merge networks
    log("\n  Merging networks...")
    multinetwork = pd.concat([ccc, grn, rg], ignore_index=True)
    multinetwork["condition"] = condition
    log(f"    Total edges: {len(multinetwork):,}")

    # Extract celltype from node names
    def get_celltype(node: str) -> str:
        parts = node.rsplit(":", 1)
        return parts[1] if len(parts) == 2 else ""

    multinetwork["source_celltype"] = multinetwork["source_node"].apply(get_celltype)
    multinetwork["target_celltype"] = multinetwork["target_node"].apply(get_celltype)

    # Add expression columns
    log("\n  Adding expression annotations...")

    def get_expression(node: str) -> float:
        gene, ct = _parse_node(node)
        if ct is None:
            return np.nan
        return expr_dict.get((gene, ct), np.nan)

    multinetwork["source_expression"] = multinetwork["source_node"].apply(get_expression)
    multinetwork["target_expression"] = multinetwork["target_node"].apply(get_expression)

    source_expr_coverage = multinetwork["source_expression"].notna().sum()
    target_expr_coverage = multinetwork["target_expression"].notna().sum()
    log(f"    Source expression coverage: {source_expr_coverage:,} / {len(multinetwork):,}")
    log(f"    Target expression coverage: {target_expr_coverage:,} / {len(multinetwork):,}")

    # Add module columns (if available)
    if gene_to_module:
        log("\n  Adding module annotations...")

        def get_module(node: str) -> str:
            gene, _ = _parse_node(node)
            return gene_to_module.get(gene, "na")

        multinetwork["source_module"] = multinetwork["source_node"].apply(get_module)
        multinetwork["target_module"] = multinetwork["target_node"].apply(get_module)

        source_mod_coverage = (multinetwork["source_module"] != "na").sum()
        target_mod_coverage = (multinetwork["target_module"] != "na").sum()
        log(f"    Source module coverage: {source_mod_coverage:,} / {len(multinetwork):,}")
        log(f"    Target module coverage: {target_mod_coverage:,} / {len(multinetwork):,}")

    # Summary by edge type
    log("\n  Edge counts by type:")
    for edge_type, count in multinetwork["edge_type"].value_counts().items():
        log(f"    {edge_type}: {count:,}")

    return multinetwork


def main(config: Optional[ReconConfig] = None) -> None:
    """Main execution function."""
    config = get_config(config)

    log("=" * 70)
    log("MULTINETWORK GENERATION")
    log(f"CCC Source: {config.ccc_source}")
    log(f"GRN score threshold: {config.grn_score_threshold}")
    log(f"Module file: {config.module_file}")
    log("=" * 70)

    output_dir = config.get_multinetwork_dir()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load module weights (shared across conditions)
    gene_to_module = load_module_weights(config)

    # Load receptor-gene base network (shared across conditions)
    rg_base = load_receptor_gene_base()

    # Process each condition
    for condition in config.conditions:
        log(f"\n{'#' * 70}")
        log(f"PROCESSING {condition.upper()}")
        log(f"{'#' * 70}")

        # Load expression data for this condition
        expr_dict = load_expression_data(condition, config)

        # Generate multinetwork
        multinetwork = generate_multinetwork(
            condition, expr_dict, gene_to_module, rg_base, config
        )

        # Save to parquet
        output_path = output_dir / f"multinetwork_{condition}.parquet"
        multinetwork.to_parquet(output_path, index=False)
        log(f"\n  Saved: {output_path}")
        log(f"  File size: {output_path.stat().st_size / (1024 * 1024):.1f} MB")

        del multinetwork, expr_dict
        gc.collect()

    # Final summary
    log(f"\n{'=' * 70}")
    log("MULTINETWORK GENERATION COMPLETE")
    log(f"{'=' * 70}")
    log(f"Output directory: {output_dir}")

    for condition in config.conditions:
        output_path = output_dir / f"multinetwork_{condition}.parquet"
        if output_path.exists():
            df = pd.read_parquet(output_path)
            log(f"\n{condition.upper()}:")
            log(f"  Total edges: {len(df):,}")
            log(f"  Edge types: {df['edge_type'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
