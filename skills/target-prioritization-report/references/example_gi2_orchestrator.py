"""
Master data loading orchestrator that coordinates all data source loaders.
"""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any

from .loaders import (
    load_cortellis_data,
    load_offx_safety_data,
    load_deg_data,
    load_biobridge_data,
    load_ultra_data,
    load_primekg_data,
    load_pathway_data,
    load_ppi_data,
    load_coexpression_data,
    load_competitive_intelligence,
    load_ci_dashboard
)

logger = logging.getLogger(__name__)


def load_all_data(
    genes: List[str],
    results_dir: str = "results",
    kgpred_dir: Optional[str] = None,
    is_combo: bool = False
) -> Dict[str, Any]:
    """
    Load all data sources for a target (individual gene or combination).

    This is the master orchestrator function that loads data from all sources
    and returns a unified dictionary with raw scores and metadata.

    Args:
        genes: List of gene names (single gene or combination)
        results_dir: Base results directory path
        kgpred_dir: Optional kgpred root path (e.g., `kgpred_ibd`). If not provided, loaders
            should resolve `kgpred_<indication>` internally (recommended) or gracefully degrade.
        is_combo: Whether this is a combination (True) or individual gene (False)

    Returns:
        Dictionary containing all loaded data with raw scores
    """
    results_path = Path(results_dir)
    kgpred_path = Path(kgpred_dir) if kgpred_dir else None
    gene_str = '-'.join(genes) if is_combo else genes[0]

    logger.info(f"Loading all data for: {gene_str} (combo={is_combo})")

    # Initialize result dictionary
    data = {
        'target': gene_str,
        'genes': genes,
        'is_combination': is_combo,
        'data_sources_loaded': [],
        'data_sources_failed': []
    }

    # Load CI dashboard once (shared across genes) to support CI family-level scoring.
    ci_dashboard = load_ci_dashboard(results_path)
    data['ci_dashboard'] = ci_dashboard

    # Load data for each gene (always needed for individuals and combos)
    for gene in genes:
        gene_prefix = f"{gene}_" if is_combo else ""

        # Cortellis
        try:
            cortellis = load_cortellis_data(gene, results_path)
            if cortellis['total_score_raw'] is not None:
                data[f'{gene_prefix}cortellis'] = cortellis
                data['data_sources_loaded'].append(f'cortellis_{gene}')
            else:
                data['data_sources_failed'].append(f'cortellis_{gene}')
        except Exception as e:
            logger.error(f"Failed to load Cortellis for {gene}: {e}")
            data['data_sources_failed'].append(f'cortellis_{gene}')

        # OFF-X Safety
        try:
            offx = load_offx_safety_data(gene, results_path)
            if offx['safety_breakdown_raw'] is not None:
                data[f'{gene_prefix}offx'] = offx
                data['data_sources_loaded'].append(f'offx_{gene}')
            else:
                data['data_sources_failed'].append(f'offx_{gene}')
        except Exception as e:
            logger.error(f"Failed to load OFF-X for {gene}: {e}")
            data['data_sources_failed'].append(f'offx_{gene}')

        # DEG
        try:
            deg = load_deg_data(gene, results_path)
            if deg['deg_score_raw'] is not None:
                data[f'{gene_prefix}deg'] = deg
                data['data_sources_loaded'].append(f'deg_{gene}')
            else:
                data['data_sources_failed'].append(f'deg_{gene}')
        except Exception as e:
            logger.error(f"Failed to load DEG for {gene}: {e}")
            data['data_sources_failed'].append(f'deg_{gene}')

        # BioBridge
        try:
            combo_name = gene_str if is_combo else None
            biobridge = load_biobridge_data(
                gene,
                results_path,
                kgpred_path=kgpred_path,
                is_combo=is_combo,
                combo_name=combo_name,
            )
            if biobridge['biobridge_percentile_raw'] is not None:
                data[f'{gene_prefix}biobridge'] = biobridge
                data['data_sources_loaded'].append(f'biobridge_{gene}')
            else:
                data['data_sources_failed'].append(f'biobridge_{gene}')
        except Exception as e:
            logger.error(f"Failed to load BioBridge for {gene}: {e}")
            data['data_sources_failed'].append(f'biobridge_{gene}')

        # ULTRA (optional)
        try:
            ultra = load_ultra_data(gene, results_path, kgpred_path=kgpred_path)
            if ultra['ultra_percentile_raw'] is not None:
                data[f'{gene_prefix}ultra'] = ultra
                data['data_sources_loaded'].append(f'ultra_{gene}')
            # Don't mark as failed if optional
        except Exception as e:
            logger.error(f"Failed to load ULTRA for {gene}: {e}")

        # PrimeKG
        try:
            primekg = load_primekg_data(gene, results_path, kgpred_path=kgpred_path)
            if primekg['primekg_connections_raw'] is not None:
                data[f'{gene_prefix}primekg'] = primekg
                data['data_sources_loaded'].append(f'primekg_{gene}')
            else:
                data['data_sources_failed'].append(f'primekg_{gene}')
        except Exception as e:
            logger.error(f"Failed to load PrimeKG for {gene}: {e}")
            data['data_sources_failed'].append(f'primekg_{gene}')

        # Competitive Intelligence
        try:
            ci = load_competitive_intelligence(gene, results_path, ci_dashboard=ci_dashboard)
            if ci['total_programs'] > 0:
                data[f'{gene_prefix}ci'] = ci
                data['data_sources_loaded'].append(f'ci_{gene}')
            else:
                data['data_sources_failed'].append(f'ci_{gene}')
        except Exception as e:
            logger.error(f"Failed to load CI for {gene}: {e}")
            data['data_sources_failed'].append(f'ci_{gene}')

    # Load combination-specific data (only if is_combo)
    if is_combo:
        # Pathway overlap
        try:
            pathway = load_pathway_data(genes, results_path)
            if pathway['shared_pathways'] is not None:
                data['pathway'] = pathway
                data['data_sources_loaded'].append('pathway')
            else:
                data['data_sources_failed'].append('pathway')
        except Exception as e:
            logger.error(f"Failed to load pathway data: {e}")
            data['data_sources_failed'].append('pathway')

        # PPI connections
        try:
            ppi = load_ppi_data(genes, results_path)
            if ppi['min_hops'] is not None:
                data['ppi'] = ppi
                data['data_sources_loaded'].append('ppi')
            else:
                data['data_sources_failed'].append('ppi')
        except Exception as e:
            logger.error(f"Failed to load PPI data: {e}")
            data['data_sources_failed'].append('ppi')

        # Coexpression
        try:
            coexpr = load_coexpression_data(genes, results_path)
            if coexpr['bulk_correlation'] is not None or coexpr['sc_correlation']:
                data['coexpression'] = coexpr
                data['data_sources_loaded'].append('coexpression')
            else:
                data['data_sources_failed'].append('coexpression')
        except Exception as e:
            logger.error(f"Failed to load coexpression data: {e}")
            data['data_sources_failed'].append('coexpression')

    # Log summary
    loaded_count = len(data['data_sources_loaded'])
    failed_count = len(data['data_sources_failed'])
    total_count = loaded_count + failed_count
    success_rate = (loaded_count / total_count * 100) if total_count > 0 else 0

    logger.info(f"Data loading complete for {gene_str}:")
    logger.info(f"  Loaded: {loaded_count}/{total_count} ({success_rate:.1f}%)")
    if failed_count > 0:
        logger.warning(f"  Failed: {data['data_sources_failed']}")

    return data


def load_multiple_targets(
    target_list: List[Dict[str, Any]],
    results_dir: str = "results"
) -> Dict[str, Dict[str, Any]]:
    """
    Load data for multiple targets (genes and combinations).

    Args:
        target_list: List of target dictionaries with 'genes' and optional 'is_combo'
        results_dir: Base results directory path

    Returns:
        Dictionary mapping target names to their data
    """
    all_data = {}

    for target_spec in target_list:
        genes = target_spec['genes']
        is_combo = target_spec.get('is_combo', len(genes) > 1)
        target_name = '-'.join(genes)

        logger.info(f"Loading target: {target_name}")

        try:
            data = load_all_data(genes, results_dir, is_combo)
            all_data[target_name] = data
        except Exception as e:
            logger.error(f"Failed to load target {target_name}: {e}")
            all_data[target_name] = {
                'target': target_name,
                'genes': genes,
                'is_combination': is_combo,
                'error': str(e)
            }

    return all_data
