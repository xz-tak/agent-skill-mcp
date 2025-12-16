#!/usr/bin/env python3
"""
Module Query and Perturbation Enrichment Script

This script provides functionality to:
1. Query genes in disease modules
2. Perform perturbation enrichment analysis (TF, Gene, LINCS-CRISPR, LINCS-CHEM)
3. Find potential regulators for target genes
4. Generate PPI network plots
"""

import os
import pickle
import argparse
import logging
from pathlib import Path
from datetime import datetime
import warnings
import pandas as pd

try:
    from ncbi_metadata import enrich_perturbation_tables
    NCBI_AVAILABLE = True
except ImportError:
    NCBI_AVAILABLE = False
    logging.warning("ncbi_metadata module not available. Install requests to enable metadata enrichment.")


def setup_logging(output_dir):
    """Setup logging configuration."""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    date = datetime.now().strftime('%Y%m%d')
    log_file = output_dir / f'module_query_{timestamp}.log'

    # Create logger
    logger = logging.getLogger('module_query')
    logger.setLevel(logging.INFO)

    # Remove existing handlers
    logger.handlers = []

    # File handler
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger, date


def query_genes_in_modules(module_path, genelist, output_dir=None):
    """
    Query which modules contain the specified genes.

    Parameters
    ----------
    module_path : str
        Path to the trained module pickle file
    genelist : list
        List of gene symbols to query
    output_dir : str, optional
        Directory to save output files. If None, uses same directory as module.

    Returns
    -------
    gene_module_df : DataFrame
        DataFrame showing which module each gene belongs to
    genes_not_found : list
        List of genes not found in the module
    """
    # Determine output directory
    if output_dir is None:
        output_dir = Path(module_path).parent
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    logger, date = setup_logging(output_dir)

    logger.info("=" * 60)
    logger.info("Module Gene Query")
    logger.info("=" * 60)
    logger.info(f"Module path: {module_path}")
    logger.info(f"Query genes: {len(genelist)}")
    logger.info(f"Output directory: {output_dir}")

    # Load module
    logger.info("Loading module...")
    with open(module_path, 'rb') as handle:
        module = pickle.load(handle)

    logger.info(f"Module has {module.n_modules} modules")

    # Create module dictionary
    module_dict = module.weights.to_df().apply(lambda x: x[x != 0].to_dict().keys()).T

    # Query genes
    genes_in_module = [k for k in set(genelist) if k in module_dict.index]
    genes_not_found = [k for k in set(genelist) if k not in module_dict.index]

    logger.info(f"Genes found in modules: {len(genes_in_module)}")
    logger.info(f"Genes NOT found: {len(genes_not_found)}")

    if genes_not_found:
        logger.warning(f"Genes not found in module: {', '.join(genes_not_found)}")

    # Get module assignments
    if genes_in_module:
        gene_module_df = module_dict.loc[genes_in_module].astype(int).sort_values(by=0)
        gene_module_df.columns = ['Module']

        # Save results
        output_file = output_dir / f'moduleweights_{date}.txt'
        gene_module_df.to_csv(output_file, sep='\t')
        logger.info(f"Gene-module mapping saved to {output_file}")

        # Print results
        logger.info("\nGene-Module Mapping:")
        logger.info(f"\n{gene_module_df}")
    else:
        gene_module_df = pd.DataFrame()
        logger.warning("No genes found in modules!")

    # Save genes not found
    if genes_not_found:
        not_found_file = output_dir / f'genes_not_found_{date}.txt'
        with open(not_found_file, 'w') as f:
            for gene in genes_not_found:
                f.write(f"{gene}\n")
        logger.info(f"Genes not found saved to {not_found_file}")

    logger.info("=" * 60)
    logger.info("Query complete!")
    logger.info("=" * 60)

    return gene_module_df, genes_not_found


def perturbation_enrichment(module_path, module_list, output_dir=None, species='human', enrich_metadata=False, email='user@example.com'):
    """
    Perform perturbation enrichment analysis for specified modules.

    Parameters
    ----------
    module_path : str
        Path to the trained module pickle file
    module_list : list
        List of module IDs to analyze
    output_dir : str, optional
        Directory to save output files. If None, uses same directory as module.
    species : str, optional
        Species for enrichment analysis. Default is 'human'.
    enrich_metadata : bool, optional
        Whether to enrich results with NCBI metadata. Default is False.
    email : str, optional
        Email for NCBI API (required if enrich_metadata=True). Default is 'user@example.com'.

    Returns
    -------
    results : dict
        Dictionary containing perturbation enrichment results
    """
    # Determine output directory
    if output_dir is None:
        output_dir = Path(module_path).parent
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    logger, date = setup_logging(output_dir)

    logger.info("=" * 60)
    logger.info("Perturbation Enrichment Analysis")
    logger.info("=" * 60)
    logger.info(f"Module path: {module_path}")
    logger.info(f"Modules to analyze: {module_list}")
    logger.info(f"Output directory: {output_dir}")

    # Load module
    logger.info("Loading module...")
    with open(module_path, 'rb') as handle:
        module = pickle.load(handle)

    # Create module dictionary
    module_dict = module.weights.to_df().apply(lambda x: x[x != 0].to_dict().keys()).T

    results = {}

    # TF regulation
    logger.info("\nRunning TF regulation analysis...")
    try:
        module._infer.tf_regulation(species=species)
        tf_regulation = module.weights.obs.TF_regulation.loc[module_list]
        logger.info(f"\n{tf_regulation}")

        tf_genes = module_dict.loc[
            [k for k in tf_regulation.keys() if k in module_dict.index]
        ].astype(int).sort_values(by=0)
        logger.info(f"\nTF genes in modules:\n{tf_genes}")

        # Save results
        tf_reg_file = output_dir / f'tf_regulation_{date}.txt'
        tf_regulation.to_frame().to_csv(tf_reg_file, sep='\t')
        logger.info(f"TF regulation saved to {tf_reg_file}")

        tf_genes_file = output_dir / f'tf_regulation_genes_{date}.txt'
        tf_genes.to_csv(tf_genes_file, sep='\t')
        logger.info(f"TF genes saved to {tf_genes_file}")

        results['tf_regulation'] = tf_regulation
        results['tf_genes'] = tf_genes
    except Exception as e:
        logger.warning(f"TF regulation analysis failed: {e}")

    # Kinase regulation
    logger.info("\nRunning kinase regulation analysis...")
    try:
        module._infer.kinase_regulation(species=species)
        kinase_regulation = module.weights.obs.kinase_regulation.loc[module_list]
        logger.info(f"\n{kinase_regulation}")

        kinase_genes = module_dict.loc[
            [k for k in kinase_regulation.keys() if k in module_dict.index]
        ].astype(int).sort_values(by=0)
        logger.info(f"\nKinase genes in modules:\n{kinase_genes}")

        # Save results
        kinase_reg_file = output_dir / f'kinase_regulation_{date}.txt'
        kinase_regulation.to_frame().to_csv(kinase_reg_file, sep='\t')
        logger.info(f"Kinase regulation saved to {kinase_reg_file}")

        kinase_genes_file = output_dir / f'kinase_regulation_genes_{date}.txt'
        kinase_genes.to_csv(kinase_genes_file, sep='\t')
        logger.info(f"Kinase genes saved to {kinase_genes_file}")

        results['kinase_regulation'] = kinase_regulation
        results['kinase_genes'] = kinase_genes
    except Exception as e:
        logger.warning(f"Kinase regulation analysis failed: {e}")

    # TF perturbation
    logger.info("\n" + "=" * 60)
    logger.info("Module non-restricted TF perturbation")
    logger.info("=" * 60)
    try:
        perturb_tf_df = module._infer.tf_pert(module_list, method='enrich_pert', tf_pert='tf', lookup=False)
        logger.info(f"\n{perturb_tf_df}")

        # Add Module_restricted column
        perturb_tf_df['Module_restricted'] = perturb_tf_df['pert_TF'].apply(
            lambda x: x.upper() in module_dict.index and
                     any(module_dict.loc[x.upper()].values[0] in module_list)
            if x.upper() in module_dict.index else False
        )

        perturb_reg_module = module_dict.loc[
            [k.upper() for k in perturb_tf_df.pert_TF if k.upper() in module_dict.index]
        ].astype(int).sort_values(by=0)
        perturb_reg_module = perturb_reg_module[~perturb_reg_module.index.duplicated(keep='first')]
        logger.info(f"\nPerturbed TF genes in modules:\n{perturb_reg_module}")

        # Enrich with NCBI metadata if requested
        if enrich_metadata and NCBI_AVAILABLE:
            logger.info("Enriching TF perturbations with NCBI metadata...")
            perturb_tf_df = enrich_perturbation_tables(perturb_tf_df, pert_type='tf', email=email)

        # Save results
        output_file = output_dir / f'perturb_tf_{date}.txt'
        perturb_tf_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved to {output_file}")

        # Save module mapping
        module_file = output_dir / f'perturb_tf_modules_{date}.txt'
        perturb_reg_module.to_csv(module_file, sep='\t')
        logger.info(f"Module mapping saved to {module_file}")

        results['perturb_tf'] = perturb_tf_df
        results['perturb_tf_modules'] = perturb_reg_module
    except Exception as e:
        logger.warning(f"TF perturbation analysis failed: {e}")

    # Gene perturbation
    logger.info("\n" + "=" * 60)
    logger.info("Module non-restricted Gene perturbation")
    logger.info("=" * 60)
    try:
        perturb_gene_df = module._infer.tf_pert(module_list, method='enrich_pert', tf_pert='gene', lookup=False)
        logger.info(f"\n{perturb_gene_df}")

        # Add Module_restricted column
        perturb_gene_df['Module_restricted'] = perturb_gene_df['pert_GENE'].apply(
            lambda x: x.upper() in module_dict.index and
                     any(module_dict.loc[x.upper()].values[0] in module_list)
            if x.upper() in module_dict.index else False
        )

        perturb_reg_module = module_dict.loc[
            [k.upper() for k in perturb_gene_df.pert_GENE if k.upper() in module_dict.index]
        ].astype(int).sort_values(by=0)
        perturb_reg_module = perturb_reg_module[~perturb_reg_module.index.duplicated(keep='first')]
        logger.info(f"\nPerturbed genes in modules:\n{perturb_reg_module}")

        # Enrich with NCBI metadata if requested
        if enrich_metadata and NCBI_AVAILABLE:
            logger.info("Enriching Gene perturbations with NCBI metadata...")
            perturb_gene_df = enrich_perturbation_tables(perturb_gene_df, pert_type='gene', email=email)

        # Save results
        output_file = output_dir / f'perturb_gene_{date}.txt'
        perturb_gene_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved to {output_file}")

        # Save module mapping
        module_file = output_dir / f'perturb_gene_modules_{date}.txt'
        perturb_reg_module.to_csv(module_file, sep='\t')
        logger.info(f"Module mapping saved to {module_file}")

        results['perturb_gene'] = perturb_gene_df
        results['perturb_gene_modules'] = perturb_reg_module
    except Exception as e:
        logger.warning(f"Gene perturbation analysis failed: {e}")

    # LINCS-CRISPR perturbation
    logger.info("\n" + "=" * 60)
    logger.info("Module non-restricted LINCS-CRISPR perturbation")
    logger.info("=" * 60)
    try:
        perturb_crispr_df = module._infer.tf_pert(module_list, method='enrich_pert', tf_pert='LINCS_CRISPR', lookup=False)
        logger.info(f"\n{perturb_crispr_df}")

        # Add Module_restricted column
        perturb_crispr_df['Module_restricted'] = perturb_crispr_df['pert_GENE'].apply(
            lambda x: x.upper() in module_dict.index and
                     any(module_dict.loc[x.upper()].values[0] in module_list)
            if x.upper() in module_dict.index else False
        )

        perturb_reg_module = module_dict.loc[
            [k.upper() for k in perturb_crispr_df.pert_GENE if k.upper() in module_dict.index]
        ].astype(int).sort_values(by=0)
        perturb_reg_module = perturb_reg_module[~perturb_reg_module.index.duplicated(keep='first')]
        logger.info(f"\nPerturbed CRISPR genes in modules:\n{perturb_reg_module}")

        # Enrich with NCBI metadata if requested
        if enrich_metadata and NCBI_AVAILABLE:
            logger.info("Enriching CRISPR perturbations with NCBI metadata...")
            perturb_crispr_df = enrich_perturbation_tables(perturb_crispr_df, pert_type='crispr', email=email)

        # Save results
        output_file = output_dir / f'perturb_crispr_{date}.txt'
        perturb_crispr_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved to {output_file}")

        # Save module mapping
        module_file = output_dir / f'perturb_crispr_modules_{date}.txt'
        perturb_reg_module.to_csv(module_file, sep='\t')
        logger.info(f"Module mapping saved to {module_file}")

        results['perturb_crispr'] = perturb_crispr_df
        results['perturb_crispr_modules'] = perturb_reg_module
    except Exception as e:
        logger.warning(f"LINCS-CRISPR perturbation analysis failed: {e}")

    # LINCS-CHEM perturbation
    logger.info("\n" + "=" * 60)
    logger.info("Module non-restricted LINCS-CHEM perturbation")
    logger.info("=" * 60)
    try:
        perturb_chem_df = module._infer.tf_pert(module_list, method='enrich_pert', tf_pert='LINCS_CHEM', lookup=False)
        logger.info(f"\n{perturb_chem_df}")

        # Add Module_restricted column
        perturb_chem_df['Module_restricted'] = perturb_chem_df['pert_CHEM'].apply(
            lambda x: x.upper() in module_dict.index and
                     any(module_dict.loc[x.upper()].values[0] in module_list)
            if x.upper() in module_dict.index else False
        )

        perturb_reg_module = module_dict.loc[
            [k.upper() for k in perturb_chem_df.pert_CHEM if k.upper() in module_dict.index]
        ].astype(int).sort_values(by=0)
        perturb_reg_module = perturb_reg_module[~perturb_reg_module.index.duplicated(keep='first')]
        logger.info(f"\nPerturbed chemicals in modules:\n{perturb_reg_module}")

        # Save results
        output_file = output_dir / f'perturb_chem_{date}.txt'
        perturb_chem_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved to {output_file}")

        # Save module mapping
        module_file = output_dir / f'perturb_chem_modules_{date}.txt'
        perturb_reg_module.to_csv(module_file, sep='\t')
        logger.info(f"Module mapping saved to {module_file}")

        results['perturb_chem'] = perturb_chem_df
        results['perturb_chem_modules'] = perturb_reg_module
    except Exception as e:
        logger.warning(f"LINCS-CHEM perturbation analysis failed: {e}")

    # PPI network plot
    logger.info("\nGenerating PPI network plot...")
    try:
        module.plot.ppi_network(module_list, save=str(output_dir / f'ppi_network_{date}'))
        logger.info(f"PPI network plot saved")
    except Exception as e:
        logger.warning(f"PPI network plot failed: {e}")

    logger.info("=" * 60)
    logger.info("Perturbation enrichment complete!")
    logger.info("=" * 60)

    return results


def find_regulators_for_targets(module_path, target_genes, module_list=None, output_dir=None, species='human', enrich_metadata=False, email='user@example.com'):
    """
    Find potential regulators for target genes using perturbation data.

    Parameters
    ----------
    module_path : str
        Path to the trained module pickle file
    target_genes : list
        List of target gene symbols
    module_list : list, optional
        List of module IDs to restrict analysis. If None, uses all modules.
    output_dir : str, optional
        Directory to save output files. If None, uses same directory as module.
    species : str, optional
        Species for enrichment analysis. Default is 'human'.
    enrich_metadata : bool, optional
        Whether to enrich results with NCBI metadata. Default is False.
    email : str, optional
        Email for NCBI API (required if enrich_metadata=True). Default is 'user@example.com'.

    Returns
    -------
    regulators : dict
        Dictionary containing regulator results for each perturbation type
    """
    # Determine output directory
    if output_dir is None:
        output_dir = Path(module_path).parent
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    logger, date = setup_logging(output_dir)

    logger.info("=" * 60)
    logger.info("Find Regulators for Target Genes")
    logger.info("=" * 60)
    logger.info(f"Module path: {module_path}")
    logger.info(f"Target genes: {target_genes}")
    logger.info(f"Output directory: {output_dir}")

    # Load module
    logger.info("Loading module...")
    with open(module_path, 'rb') as handle:
        module = pickle.load(handle)

    # Create module dictionary
    module_dict = module.weights.to_df().apply(lambda x: x[x != 0].to_dict().keys()).T

    # Use all modules if not specified
    if module_list is None:
        module_list = module.weights.obs.index.tolist()
        logger.info(f"Using all {len(module_list)} modules")
    else:
        logger.info(f"Using specified modules: {module_list}")

    regulators = {}

    # TF perturbation
    logger.info("\n" + "=" * 60)
    logger.info(f"Module non-restricted TF perturbation with targets: {target_genes}")
    logger.info("=" * 60)
    try:
        perturb_tf_df = module._infer.tf_pert(module_list, method='enrich_pert', tf_pert='tf', lookup=False)

        msk = perturb_tf_df['overlap_genes'].apply(lambda x: not set(x).isdisjoint(target_genes))
        target2reg = perturb_tf_df.loc[msk, :]

        if len(target2reg) > 0:
            # Add Module_restricted column
            target2reg['Module_restricted'] = target2reg['pert_TF'].apply(
                lambda x: x.upper() in module_dict.index and
                         any(module_dict.loc[x.upper()].values[0] in module_list)
                if x.upper() in module_dict.index else False
            )

            logger.info(f"\nTF regulators found:\n{target2reg}")

            perturb_reg_module = module_dict.loc[
                [k.upper() for k in target2reg.pert_TF if k.upper() in module_dict.index]
            ].astype(int).sort_values(by=0)
            perturb_reg_module = perturb_reg_module[~perturb_reg_module.index.duplicated(keep='first')]
            logger.info(f"\nRegulator TFs in modules:\n{perturb_reg_module}")

            # Enrich with NCBI metadata if requested
            if enrich_metadata and NCBI_AVAILABLE:
                logger.info("Enriching TF regulators with NCBI metadata...")
                target2reg = enrich_perturbation_tables(target2reg, pert_type='tf', email=email)

            # Save results
            output_file = output_dir / f'regulators_tf_{date}.txt'
            target2reg.to_csv(output_file, sep='\t', index=False)
            logger.info(f"Saved to {output_file}")

            # Save module mapping
            module_file = output_dir / f'regulators_tf_modules_{date}.txt'
            perturb_reg_module.to_csv(module_file, sep='\t')
            logger.info(f"Module mapping saved to {module_file}")

            regulators['tf'] = target2reg
            regulators['tf_modules'] = perturb_reg_module
        else:
            logger.info("No TF regulators found")
    except Exception as e:
        logger.warning(f"TF regulator analysis failed: {e}")

    # Gene perturbation
    logger.info("\n" + "=" * 60)
    logger.info(f"Module non-restricted GENE perturbation with targets: {target_genes}")
    logger.info("=" * 60)
    try:
        perturb_gene_df = module._infer.tf_pert(module_list, method='enrich_pert', tf_pert='gene', lookup=False)

        msk = perturb_gene_df['overlap_genes'].apply(lambda x: not set(x).isdisjoint(target_genes))
        target2reg = perturb_gene_df.loc[msk, :]

        if len(target2reg) > 0:
            # Add Module_restricted column
            target2reg['Module_restricted'] = target2reg['pert_GENE'].apply(
                lambda x: x.upper() in module_dict.index and
                         any(module_dict.loc[x.upper()].values[0] in module_list)
                if x.upper() in module_dict.index else False
            )

            logger.info(f"\nGene regulators found:\n{target2reg}")

            perturb_reg_module = module_dict.loc[
                [k.upper() for k in target2reg.pert_GENE if k.upper() in module_dict.index]
            ].astype(int).sort_values(by=0)
            perturb_reg_module = perturb_reg_module[~perturb_reg_module.index.duplicated(keep='first')]
            logger.info(f"\nRegulator genes in modules:\n{perturb_reg_module}")

            # Enrich with NCBI metadata if requested
            if enrich_metadata and NCBI_AVAILABLE:
                logger.info("Enriching Gene regulators with NCBI metadata...")
                target2reg = enrich_perturbation_tables(target2reg, pert_type='gene', email=email)

            # Save results
            output_file = output_dir / f'regulators_gene_{date}.txt'
            target2reg.to_csv(output_file, sep='\t', index=False)
            logger.info(f"Saved to {output_file}")

            # Save module mapping
            module_file = output_dir / f'regulators_gene_modules_{date}.txt'
            perturb_reg_module.to_csv(module_file, sep='\t')
            logger.info(f"Module mapping saved to {module_file}")

            regulators['gene'] = target2reg
            regulators['gene_modules'] = perturb_reg_module
        else:
            logger.info("No gene regulators found")
    except Exception as e:
        logger.warning(f"Gene regulator analysis failed: {e}")

    # LINCS-CRISPR perturbation
    logger.info("\n" + "=" * 60)
    logger.info(f"Module non-restricted LINCS-CRISPR perturbation with targets: {target_genes}")
    logger.info("=" * 60)
    try:
        perturb_crispr_df = module._infer.tf_pert(module_list, method='enrich_pert', tf_pert='LINCS_CRISPR', lookup=False)

        msk = perturb_crispr_df['overlap_genes'].apply(lambda x: not set(x).isdisjoint(target_genes))
        target2reg = perturb_crispr_df.loc[msk, :]

        if len(target2reg) > 0:
            # Add Module_restricted column
            target2reg['Module_restricted'] = target2reg['pert_GENE'].apply(
                lambda x: x.upper() in module_dict.index and
                         any(module_dict.loc[x.upper()].values[0] in module_list)
                if x.upper() in module_dict.index else False
            )

            logger.info(f"\nLINCS-CRISPR regulators found:\n{target2reg}")

            perturb_reg_module = module_dict.loc[
                [k.upper() for k in target2reg.pert_GENE if k.upper() in module_dict.index]
            ].astype(int).sort_values(by=0)
            perturb_reg_module = perturb_reg_module[~perturb_reg_module.index.duplicated(keep='first')]
            logger.info(f"\nRegulator CRISPR genes in modules:\n{perturb_reg_module}")

            # Enrich with NCBI metadata if requested
            if enrich_metadata and NCBI_AVAILABLE:
                logger.info("Enriching CRISPR regulators with NCBI metadata...")
                target2reg = enrich_perturbation_tables(target2reg, pert_type='crispr', email=email)

            # Save results
            output_file = output_dir / f'regulators_crispr_{date}.txt'
            target2reg.to_csv(output_file, sep='\t', index=False)
            logger.info(f"Saved to {output_file}")

            # Save module mapping
            module_file = output_dir / f'regulators_crispr_modules_{date}.txt'
            perturb_reg_module.to_csv(module_file, sep='\t')
            logger.info(f"Module mapping saved to {module_file}")

            regulators['crispr'] = target2reg
            regulators['crispr_modules'] = perturb_reg_module
        else:
            logger.info("No LINCS-CRISPR regulators found")
    except Exception as e:
        logger.warning(f"LINCS-CRISPR regulator analysis failed: {e}")

    # LINCS-CHEM perturbation
    logger.info("\n" + "=" * 60)
    logger.info(f"Module non-restricted LINCS-CHEM perturbation with targets: {target_genes}")
    logger.info("=" * 60)
    try:
        perturb_chem_df = module._infer.tf_pert(module_list, method='enrich_pert', tf_pert='LINCS_CHEM', lookup=False)

        msk = perturb_chem_df['overlap_genes'].apply(lambda x: not set(x).isdisjoint(target_genes))
        target2reg = perturb_chem_df.loc[msk, :]

        if len(target2reg) > 0:
            # Add Module_restricted column
            target2reg['Module_restricted'] = target2reg['pert_CHEM'].apply(
                lambda x: x.upper() in module_dict.index and
                         any(module_dict.loc[x.upper()].values[0] in module_list)
                if x.upper() in module_dict.index else False
            )

            logger.info(f"\nLINCS-CHEM regulators found:\n{target2reg}")

            perturb_reg_module = module_dict.loc[
                [k.upper() for k in target2reg.pert_CHEM if k.upper() in module_dict.index]
            ].astype(int).sort_values(by=0)
            perturb_reg_module = perturb_reg_module[~perturb_reg_module.index.duplicated(keep='first')]
            logger.info(f"\nRegulator chemicals in modules:\n{perturb_reg_module}")

            # Save results
            output_file = output_dir / f'regulators_chem_{date}.txt'
            target2reg.to_csv(output_file, sep='\t', index=False)
            logger.info(f"Saved to {output_file}")

            # Save module mapping
            module_file = output_dir / f'regulators_chem_modules_{date}.txt'
            perturb_reg_module.to_csv(module_file, sep='\t')
            logger.info(f"Module mapping saved to {module_file}")

            regulators['chem'] = target2reg
            regulators['chem_modules'] = perturb_reg_module
        else:
            logger.info("No LINCS-CHEM regulators found")
    except Exception as e:
        logger.warning(f"LINCS-CHEM regulator analysis failed: {e}")

    logger.info("=" * 60)
    logger.info("Regulator search complete!")
    logger.info("=" * 60)

    return regulators


def main():
    parser = argparse.ArgumentParser(
        description='Query genes in disease modules and perform perturbation enrichment',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Query genes in modules
  python module_query.py --module module.pickle --query-genes genes.txt --mode query

  # Perturbation enrichment for modules
  python module_query.py --module module.pickle --modules 0 1 2 --mode perturb

  # Find regulators for target genes
  python module_query.py --module module.pickle --target-genes targets.txt --mode regulators

  # Find regulators in specific modules
  python module_query.py --module module.pickle --target-genes targets.txt --modules 0 1 2 --mode regulators
        """
    )

    parser.add_argument('--module', required=True, help='Path to trained module pickle file')
    parser.add_argument('--mode', required=True, choices=['query', 'perturb', 'regulators'],
                        help='Analysis mode: query genes, perturbation enrichment, or find regulators')
    parser.add_argument('--query-genes', help='Path to file with genes to query (one per line)')
    parser.add_argument('--target-genes', help='Path to file with target genes (one per line)')
    parser.add_argument('--modules', nargs='+', help='Module IDs to analyze')
    parser.add_argument('--output', help='Output directory (default: same as module)')
    parser.add_argument('--species', default='human', help='Species for enrichment (default: human)')
    parser.add_argument('--enrich-metadata', action='store_true',
                        help='Enrich results with NCBI metadata (species, tissue, cell type, disease)')
    parser.add_argument('--email', default='user@example.com',
                        help='Email for NCBI API (required if using --enrich-metadata)')

    args = parser.parse_args()

    if args.mode == 'query':
        if not args.query_genes:
            parser.error("--query-genes required for query mode")

        # Load gene list
        with open(args.query_genes, 'r') as f:
            genelist = [line.strip() for line in f if line.strip()]

        query_genes_in_modules(
            module_path=args.module,
            genelist=genelist,
            output_dir=args.output,
        )

    elif args.mode == 'perturb':
        if not args.modules:
            parser.error("--modules required for perturb mode")

        perturbation_enrichment(
            module_path=args.module,
            module_list=[str(m) for m in args.modules],
            output_dir=args.output,
            species=args.species,
            enrich_metadata=args.enrich_metadata,
            email=args.email,
        )

    elif args.mode == 'regulators':
        if not args.target_genes:
            parser.error("--target-genes required for regulators mode")

        # Load target genes
        with open(args.target_genes, 'r') as f:
            target_genes = [line.strip() for line in f if line.strip()]

        module_list = [str(m) for m in args.modules] if args.modules else None

        find_regulators_for_targets(
            module_path=args.module,
            target_genes=target_genes,
            module_list=module_list,
            output_dir=args.output,
            species=args.species,
            enrich_metadata=args.enrich_metadata,
            email=args.email,
        )


if __name__ == '__main__':
    main()
