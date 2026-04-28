#!/usr/bin/env python3
"""
M4: ReCoN Multicellular Coordination Analysis

Generalised from the SSc Lung Atlas Phase 4 pipeline. This module runs ReCoN
multicellular coordination network analysis using cell-type-specific 5-layer GRNs.

Key features preserved from source:
- receptor_grn.copy() to prevent in-place mutation by ReCoN library
- direct.copy() and indirect.copy() before combine_effects() (prevents normalization)
- Fill NaN lr_means with 0 for CCC data
- CellChat gene alias mapping to standardize column names
- multixrank_patch for deeper network exploration

Usage:
    # With JSON config
    python m4_recon_analysis.py --config config.json

    # Specify single condition
    python m4_recon_analysis.py --config config.json --condition ssc

    # Use specific seed genes
    python m4_recon_analysis.py --config config.json --seeds TGFB1 IL6 TNF
"""

import argparse
import gc
import json
import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from config import ReconConfig, get_config

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Inline multixrank patch for sparse networks (norm=0 division by zero)
# ---------------------------------------------------------------------------

_PATCH_APPLIED = False


def _apply_multixrank_patch() -> None:
    """
    Monkey-patch multixrank.TransitionMatrix to handle sparse networks.

    When a node has no inter-layer connections (norm=0), the original code
    divides by zero at TransitionMatrix.py:177, causing NaN propagation.
    This patch uses epsilon=1e-10 when norm=0.

    Fixes sparse CCC networks (e.g., IPF with 14,513 vs SSc with 20,051
    interactions) where some nodes become isolated in the multilayer network.
    """
    global _PATCH_APPLIED
    if _PATCH_APPLIED:
        return

    try:
        import numpy
        import scipy.sparse
        from multixrank.TransitionMatrix import TransitionMatrix

        def _patched_normalization(self, alpha, beta, matrix, diago):
            EPSILON = 1e-10
            Tot_strength_nodes = matrix.sum(axis=0)
            adjacency = self.multiplex_all.supra_adj_matrix_list[beta].sum(axis=0)
            adjacency = list(adjacency.flat)
            diago_up = list(Tot_strength_nodes.flat)
            diago_down = list(Tot_strength_nodes.flat)
            self_L = [len(x.layer_tuple) for x in self.multiplex_all.multiplex_tuple]

            for k in range(len(diago_up)):
                if (self_L[beta] == 1) and (adjacency[k] == 0):
                    diago_up[k] = 0
                    list_value_diago = numpy.zeros(len(self_L))
                    for l in range(len(self_L)):
                        if diago[l][k] != 0:
                            list_value_diago[l] = 1
                    list_value_diago = list_value_diago * self.lamb[:, beta].T
                    norm = 0
                    for l in range(len(self_L)):
                        if l != beta:
                            norm += list_value_diago[l]
                    if diago_down[k] == 0:
                        diago_down[k] = 1
                    else:
                        if norm == 0:
                            diago_down[k] = EPSILON * (1 / diago_down[k])
                        else:
                            diago_down[k] = (self.lamb[alpha, beta] / norm) * (1 / diago_down[k])
                else:
                    diago_down[k] = 0
                    if diago_up[k] == 0:
                        diago_up[k] = 1
                    else:
                        diago_up[k] = self.lamb[alpha, beta] * (1 / diago_up[k])

            Normalization_matrix_up = scipy.sparse.diags(diago_up, format="coo")
            Normalization_matrix_down = scipy.sparse.diags(diago_down, format="coo")
            Transition_up = matrix.dot(Normalization_matrix_up)
            Transition_down = matrix.dot(Normalization_matrix_down)
            return Transition_up + Transition_down

        TransitionMatrix.get_normalization_bipartite_alpha_beta = _patched_normalization
        _PATCH_APPLIED = True
        _log("[PATCH] Applied multixrank TransitionMatrix patch for sparse networks")

    except ImportError:
        _log("  Warning: multixrank not installed, skipping sparse-network patch")


# ---------------------------------------------------------------------------
# CellChat gene alias -> HGNC mapping
# ---------------------------------------------------------------------------

CELLCHAT_ALIAS_TO_HGNC = {
    "CD45": "PTPRC",
    "NKG2D": "KLRK1",
    "VEGFR1": "FLT1",
    "VEGFR2": "KDR",
    "VEGFR1R2": "FLT1",
    "JAM1": "F11R",
    "PIRB": "LILRB3",
    "FASL": "FASLG",
    "ICOSL": "ICOSLG",
    "INTEGRIN": None,
    "HLA": None,
    "HOXB": None,
    "HOXD": None,
}


def _normalize_gene_name(gene: str) -> Optional[str]:
    """Convert CellChat alias to HGNC symbol. Returns None for generic names."""
    if gene in CELLCHAT_ALIAS_TO_HGNC:
        return CELLCHAT_ALIAS_TO_HGNC[gene]
    return gene


def _log(msg: str) -> None:
    """Print timestamped message."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {msg}", flush=True)


def load_seeds(config: ReconConfig) -> List[str]:
    """
    Load seed genes from config.

    Priority:
    1. CLI seeds (--seeds)
    2. config.seeds_file
    3. config.seeds list

    Args:
        config: Pipeline configuration.

    Returns:
        List of seed gene names (e.g., ["TGFB1", "IL6"])
    """
    seeds = config.load_seeds()

    if not seeds:
        _log("ERROR: No seeds provided. Use --seeds, --seeds-file, or config.seeds")
        sys.exit(1)

    _log(f"Loaded {len(seeds):,} seed genes")
    if len(seeds) <= 20:
        _log(f"  Seeds: {seeds}")
    else:
        _log(f"  Seeds (first 20): {seeds[:20]}...")

    return seeds


def load_celltype_grns(
    config: ReconConfig, condition: str, min_weight: float = 1.0
) -> Dict[str, pd.DataFrame]:
    """
    Load cell-type-specific GRNs for a condition.

    Uses RNA network files (GRNBoost2 direct edges) rather than 5-layer
    integrated networks which may be too large.

    Args:
        config: Pipeline configuration.
        condition: Condition name (lowercase, e.g., 'ssc', 'ipf', 'normal').
        min_weight: Minimum GRNBoost2 weight threshold.

    Returns:
        Dict mapping cell type (proper case) to GRN DataFrame.
    """
    # Load celltypes from M1 metadata
    meta_path = Path(config.output_dir) / "data_prep_metadata.json"
    if not meta_path.exists():
        _log("ERROR: data_prep_metadata.json not found. Run M1 first.")
        return {}

    with open(meta_path) as f:
        meta = json.load(f)
    celltypes = meta.get("celltypes", [])

    if not celltypes:
        _log("ERROR: No celltypes in metadata")
        return {}

    _log(f"Loading cell-type GRNs for {condition} (weight > {min_weight})...")
    grn_dir = config.get_grn_dir()
    grns = {}

    for celltype in celltypes:
        ct_key = celltype.lower().replace(" ", "_")
        grn_path = grn_dir / f"{ct_key}_{condition}_rna_network.csv"

        if grn_path.exists():
            grn = pd.read_csv(grn_path)
            original_count = len(grn)

            # Filter by weight
            grn = grn[grn["weight"] > min_weight]
            _log(f"  {celltype}: {original_count:,} -> {len(grn):,} edges (weight > {min_weight})")

            grns[celltype] = grn
        else:
            _log(f"  {celltype}: No GRN found at {grn_path}")

    return grns


def load_ccc(config: ReconConfig, condition: str) -> pd.DataFrame:
    """
    Load CCC data from specified source.

    Handles three sources with different score semantics:
    - merged: Uses OVERALL percentile-rank (0-1) as lr_means
    - cellphonedb: Uses RAW lr_means (original score)
    - cellchat: Uses RAW prob as lr_means (original score)

    IMPORTANT: Fills NaN lr_means with 0 (prevents downstream errors).

    Args:
        config: Pipeline configuration.
        condition: Condition name (lowercase).

    Returns:
        DataFrame with CCC interactions.
    """
    ccc_source = config.ccc_source
    ccc_dir = config.get_ccc_dir()

    _log(f"Loading CCC ({ccc_source}) for {condition}...")

    if ccc_source == "cellchat":
        # Per-condition CellChat file
        path = ccc_dir / "cellchat" / f"{condition}_ccc.csv"
        if not path.exists():
            raise FileNotFoundError(f"CellChat CCC not found: {path}")

        ccc = pd.read_csv(path)
        _log(f"  CCC (CellChat): {len(ccc):,} interactions")

        # Standardize CellChat column names
        ccc = ccc.rename(
            columns={
                "source": "celltype_source",
                "target": "celltype_target",
                "ligand": "source",
                "receptor": "target",
            }
        )
        # Handle alternative CellChat column names (ligand_name/receptor_name)
        if "ligand_name" in ccc.columns:
            ccc = ccc.rename(columns={"ligand_name": "source", "receptor_name": "target"})
        # Use RAW prob as lr_means (no percentile normalization)
        ccc["lr_means"] = ccc["prob"]

        # Apply CellChat gene alias normalization (CD45→PTPRC, VEGFR1→FLT1, etc.)
        alias_count = 0
        for col in ["source", "target"]:
            if col in ccc.columns:
                original = ccc[col].copy()
                ccc[col] = ccc[col].map(_normalize_gene_name)
                null_mask = ccc[col].isna()
                if null_mask.any():
                    ccc = ccc[~null_mask]
                alias_count += int((original != ccc[col]).sum())
        if alias_count > 0:
            _log(f"    Applied {alias_count} CellChat alias conversions")

    elif ccc_source == "cellphonedb":
        # Per-condition CellPhoneDB file
        path = ccc_dir / "cellphonedb" / f"{condition}_ccc.csv"
        if not path.exists():
            raise FileNotFoundError(f"CellPhoneDB CCC not found: {path}")

        ccc = pd.read_csv(path)
        _log(f"  CCC (CellPhoneDB): {len(ccc):,} interactions")
        # lr_means already exists with original values

    else:  # merged
        # Merged CCC (consensus from cellphonedb + cellchat)
        path = ccc_dir / "merged" / f"{condition}_ccc.csv"
        if not path.exists():
            raise FileNotFoundError(f"Merged CCC not found: {path}")

        ccc = pd.read_csv(path)
        _log(f"  CCC (merged): {len(ccc):,} interactions")
        if "source_method" in ccc.columns:
            method_counts = ccc["source_method"].value_counts().to_dict()
            _log(f"    Sources: {method_counts}")

    # PRESERVE: Fill NaN lr_means with 0
    if "lr_means" in ccc.columns:
        nan_count = ccc["lr_means"].isna().sum()
        if nan_count > 0:
            ccc["lr_means"] = ccc["lr_means"].fillna(0)
            _log(f"    Filled {nan_count:,} NaN lr_means with 0")

    return ccc


def load_receptor_gene_network() -> pd.DataFrame:
    """Load receptor-gene relationships from NicheNet PKN."""
    from recon.data import load_receptor_genes

    _log("Loading receptor-gene network from NicheNet...")
    receptor_grn = load_receptor_genes("human_receptor_gene_from_NichenetPKN")
    _log(f"  {len(receptor_grn):,} receptor-gene relationships")
    _log(f"  {receptor_grn['source'].nunique()} unique receptors")
    _log(f"  {receptor_grn['target'].nunique()} unique targets")
    return receptor_grn


def merge_celltype_grns(grns: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Merge cell-type-specific GRNs into combined network.

    Takes union of all edges, using max weight for duplicates.
    """
    if not grns:
        return pd.DataFrame(columns=["source", "target", "weight"])

    all_edges = []
    for celltype, grn in grns.items():
        if grn.empty:
            continue
        grn = grn.copy()

        if "source" in grn.columns and "target" in grn.columns and "weight" in grn.columns:
            grn = grn[["source", "target", "weight"]].copy()
            grn["celltype"] = celltype
            all_edges.append(grn)
        else:
            _log(f"  Warning: {celltype} GRN missing columns, skipping")
            continue

    if not all_edges:
        return pd.DataFrame(columns=["source", "target", "weight"])

    combined = pd.concat(all_edges, ignore_index=True)

    # Group by source-target, take max weight
    merged = combined.groupby(["source", "target"], as_index=False).agg({"weight": "max"})

    _log(f"  Merged GRN: {len(merged):,} unique edges from {len(grns)} cell types")
    return merged


def filter_seeds_to_network(
    seeds: List[str],
    grn: pd.DataFrame,
    receptor_grn: pd.DataFrame,
    ccc: pd.DataFrame,
) -> List[str]:
    """Filter seed genes to those present in the networks."""
    grn_genes = set(grn["source"]) | set(grn["target"])
    receptor_genes = set(receptor_grn["source"]) | set(receptor_grn["target"])
    ccc_genes = set(ccc["source"]) | set(ccc["target"])
    all_genes = grn_genes | receptor_genes | ccc_genes

    valid_seeds = []
    missing_genes = []

    for gene in seeds:
        if gene in all_genes:
            valid_seeds.append(gene)
        else:
            missing_genes.append(gene)

    if missing_genes:
        unique_missing = list(set(missing_genes))[:20]
        _log(f"  Warning: {len(set(missing_genes))} seeds not in network (first 20): {unique_missing}")

    return valid_seeds


def run_recon_analysis(
    condition_name: str,
    grns: Dict[str, pd.DataFrame],
    ccc: pd.DataFrame,
    receptor_grn: pd.DataFrame,
    seeds: List[str],
    config: ReconConfig,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    """
    Run ReCoN multicellular analysis for a condition.

    Args:
        condition_name: Condition identifier.
        grns: Cell-type-specific GRNs.
        ccc: Cell-cell communication network.
        receptor_grn: Receptor-gene relationships.
        seeds: Seed genes.
        config: Pipeline configuration.

    Returns:
        (direct_effect, indirect_effect, stats)
    """
    from recon.explore import combine_effects, multicell_targets

    _log("")
    _log("=" * 60)
    _log(f"RECON ANALYSIS FOR: {condition_name.upper()}")
    _log("=" * 60)

    start_time = datetime.now()

    # Merge cell-type GRNs
    merged_grn = merge_celltype_grns(grns)

    if merged_grn.empty:
        _log("ERROR: No GRN edges available!")
        return pd.DataFrame(), pd.DataFrame(), {"error": "No GRN edges"}

    celltypes = list(grns.keys())
    _log(f"Cell types: {celltypes}")

    # Filter seeds to network
    valid_seeds = filter_seeds_to_network(seeds, merged_grn, receptor_grn, ccc)
    _log(f"Seeds: {len(valid_seeds)}/{len(seeds)} valid")
    if len(valid_seeds) <= 20:
        _log(f"  Valid: {valid_seeds}")
    else:
        _log(f"  Valid (first 20): {valid_seeds[:20]}...")

    if len(valid_seeds) == 0:
        _log("ERROR: No valid seeds found in network!")
        return pd.DataFrame(), pd.DataFrame(), {"error": "No valid seeds"}

    # Run multicell_targets
    _log(f"\nRunning multicell_targets...")
    _log(f"  Cell types: {len(celltypes)}")
    _log(f"  GRN edges: {len(merged_grn):,}")
    _log(f"  CCC interactions: {len(ccc):,}")
    _log(f"  Receptor-gene edges: {len(receptor_grn):,}")
    _log(f"  Restart probability: {config.restart_proba}")
    _log(f"  Alpha (indirect weight): {config.alpha}")

    try:
        # PRESERVE: receptor_grn.copy() to prevent in-place mutation by ReCoN
        # The library adds "_receptor" suffix to gene names in-place, corrupting
        # the network for subsequent conditions without this copy
        receptor_grn_copy = receptor_grn.copy()

        direct_effect, indirect_effect = multicell_targets(
            seeds=valid_seeds,
            celltypes=celltypes,
            ccc=ccc,
            grn=merged_grn,
            receptor_grn=receptor_grn_copy,
            restart_proba=config.restart_proba,
            extend_seeds=True,  # Let ReCoN create gene-celltype combinations
            njobs=config.n_jobs,
            verbose=True,
        )

        _log(f"\nResults:")
        _log(f"  Direct effect shape: {direct_effect.shape}")
        _log(f"  Indirect effect shape: {indirect_effect.shape}")

        direct_nan = direct_effect.isna().sum().sum()
        if direct_nan > 0:
            _log(f"  Warning: {direct_nan} NaN values in direct effects")

    except Exception as e:
        _log(f"Error in multicell_targets: {e}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame(), pd.DataFrame(), {"error": str(e)}

    # Compute statistics
    elapsed = datetime.now() - start_time
    stats = {
        "condition": condition_name,
        "n_seeds": len(valid_seeds),
        "seeds": valid_seeds,
        "n_celltypes": len(celltypes),
        "celltypes": celltypes,
        "n_grn_edges": int(len(merged_grn)),
        "n_ccc_interactions": int(len(ccc)),
        "n_receptor_gene_edges": int(len(receptor_grn)),
        "restart_proba": config.restart_proba,
        "alpha": config.alpha,
        "direct_effect_shape": list(direct_effect.shape),
        "indirect_effect_shape": list(indirect_effect.shape)
        if hasattr(indirect_effect, "shape")
        else "MultiIndex",
        "compute_time_seconds": elapsed.total_seconds(),
    }

    return direct_effect, indirect_effect, stats


def combine_and_save_results(
    condition_name: str,
    direct_effect: pd.DataFrame,
    indirect_effect: pd.DataFrame,
    config: ReconConfig,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Combine effects and save results for one condition.

    PRESERVE: direct.copy() and indirect.copy() BEFORE combine_effects()
    because combine_effects modifies direct_effect in-place (normalizes columns).
    This prevents division-by-zero NaN corruption if column sums are ~0.

    Returns:
        (direct_copy, indirect_copy, combined)
    """
    from recon.explore import combine_effects

    recon_dir = config.get_recon_dir()

    # PRESERVE: Save copies BEFORE combine_effects modifies them in-place
    direct_copy = direct_effect.copy()
    indirect_copy = indirect_effect.copy()

    # Combine effects
    if not direct_effect.empty:
        combined = combine_effects(direct_effect, indirect_effect, alpha=config.alpha)
    else:
        combined = pd.DataFrame()

    # Save all three
    for result_type, df in [
        ("direct", direct_copy),
        ("indirect", indirect_copy),
        ("combined", combined),
    ]:
        if isinstance(df, pd.DataFrame) and not df.empty:
            path = recon_dir / f"{condition_name}_{result_type}_effects.csv"
            df.to_csv(path)
            _log(f"  Saved: {path}")

    return direct_copy, indirect_copy, combined


def compare_conditions(
    condition1: str,
    combined1: pd.DataFrame,
    condition2: str,
    combined2: pd.DataFrame,
) -> pd.DataFrame:
    """Compare combined effects between two conditions."""
    if combined1.empty or combined2.empty:
        _log(f"Cannot compare {condition1} vs {condition2}: missing data")
        return pd.DataFrame()

    # Common cell types and genes
    common_celltypes = sorted(set(combined1.columns) & set(combined2.columns))
    common_genes = sorted(set(combined1.index) & set(combined2.index))

    _log(f"  Common cell types: {len(common_celltypes)}")
    _log(f"  Common genes: {len(common_genes)}")

    if not common_celltypes or not common_genes:
        return pd.DataFrame()

    # Compute differential effects
    comparison_results = []
    eps = 1e-10

    for celltype in common_celltypes:
        scores1 = combined1.loc[common_genes, celltype]
        scores2 = combined2.loc[common_genes, celltype]

        log2_fc = np.log2((scores1 + eps) / (scores2 + eps))
        correlation = scores1.corr(scores2)

        comparison_results.append(
            {
                "celltype": celltype,
                f"mean_{condition1}": float(scores1.mean()),
                f"mean_{condition2}": float(scores2.mean()),
                "mean_log2_fc": float(log2_fc.mean()),
                "correlation": float(correlation),
                "n_upregulated": int((log2_fc > 1).sum()),
                "n_downregulated": int((log2_fc < -1).sum()),
            }
        )

        _log(f"  {celltype}: corr={correlation:.3f}, up={(log2_fc > 1).sum()}, down={(log2_fc < -1).sum()}")

    return pd.DataFrame(comparison_results)


def main(config: Optional[ReconConfig] = None, condition: Optional[str] = None) -> Dict:
    """
    Run ReCoN analysis pipeline.

    Args:
        config: Optional ReconConfig. Falls back to CLI args if None.
        condition: Optional single condition to process. If None, processes all.

    Returns:
        Dict with all results and statistics.
    """
    config = get_config(config)

    # Apply multixrank patch for sparse networks BEFORE any analysis
    _apply_multixrank_patch()

    # Validate config
    errors = config.validate()
    if errors:
        for e in errors:
            _log(f"Config error: {e}")
        sys.exit(1)

    _log("=" * 70)
    _log("M4: RECON MULTICELLULAR COORDINATION ANALYSIS")
    _log(f"CCC Source: {config.ccc_source}")
    _log("=" * 70)

    start_time = datetime.now()
    recon_dir = config.get_recon_dir()
    recon_dir.mkdir(parents=True, exist_ok=True)

    # Load seeds
    seeds = load_seeds(config)

    # Load receptor-gene network (shared across conditions)
    receptor_grn = load_receptor_gene_network()

    # Collect results
    all_stats = {}
    all_results = {}

    # Determine which conditions to process
    conditions_to_process = [condition] if condition else config.conditions

    # Run ReCoN for each condition
    for cond in conditions_to_process:
        cond_key = cond.lower().replace(" ", "_")

        _log(f"\n{'#' * 70}")
        _log(f"PROCESSING {cond.upper()}")
        _log(f"{'#' * 70}")

        # Load cell-type GRNs
        grns = load_celltype_grns(config, cond_key, min_weight=config.min_grn_weight)

        if not grns:
            _log(f"No GRNs available for {cond}, skipping")
            continue

        # Load CCC
        ccc = load_ccc(config, cond_key)

        # Run analysis
        direct, indirect, stats = run_recon_analysis(
            cond, grns, ccc, receptor_grn, seeds, config
        )
        all_stats[cond] = stats

        # Combine and save
        direct_copy, indirect_copy, combined = combine_and_save_results(
            cond_key, direct, indirect, config
        )

        all_results[cond] = {
            "direct": direct_copy,
            "indirect": indirect_copy,
            "combined": combined,
        }

        gc.collect()

    # Compare conditions
    if len(all_results) >= 2:
        _log(f"\n{'#' * 70}")
        _log("CONDITION COMPARISONS")
        _log(f"{'#' * 70}")

        # Compare each disease condition vs normal
        for disease_cond in config.disease_conditions:
            if disease_cond in all_results and config.normal_condition in all_results:
                if (not all_results[disease_cond]["combined"].empty
                    and not all_results[config.normal_condition]["combined"].empty):
                    _log(f"\n=== {disease_cond.title()} vs {config.normal_condition.title()} ===")
                    comparison = compare_conditions(
                        disease_cond,
                        all_results[disease_cond]["combined"],
                        config.normal_condition,
                        all_results[config.normal_condition]["combined"],
                    )

                    comparison_key = f"{disease_cond}_vs_{config.normal_condition}"
                    comparison_path = recon_dir / f"comparison_{comparison_key}.csv"
                    comparison.to_csv(comparison_path, index=False)
                    _log(f"  Saved: {comparison_path}")

    # Save statistics
    stats_path = recon_dir / "recon_stats.json"
    with open(stats_path, "w") as f:
        json.dump(all_stats, f, indent=2, default=str)
    _log(f"\nSaved statistics: {stats_path}")

    # Summary
    elapsed = datetime.now() - start_time
    _log(f"\n{'=' * 70}")
    _log("M4 COMPLETE")
    _log(f"{'=' * 70}")
    _log(f"Total time: {elapsed}")
    _log(f"Output directory: {recon_dir}")

    return {"results": all_results, "stats": all_stats}


def _build_m4_parser() -> argparse.ArgumentParser:
    """Build M4-specific CLI argument parser."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--condition",
        type=str,
        default=None,
        help="Process single condition (default: all)",
    )
    parser.add_argument("--config", type=str, help="Path to JSON config file")
    return parser


if __name__ == "__main__":
    # Parse CLI args
    parser = argparse.ArgumentParser(
        description="M4: ReCoN Multicellular Coordination Analysis",
        parents=[_build_m4_parser()],
        add_help=False,
    )
    args = parser.parse_args()

    # Run
    main(condition=args.condition)
