#!/usr/bin/env python3
"""
Disease Module Analysis Script

This script provides functionality to:
1. Train study-specific disease modules from single-cell data
2. Transform data using pre-trained modules
3. Generate dotplots and enrichment analysis
4. Save results with timestamps and logging
"""

import os
import pickle
import argparse
import logging
from pathlib import Path
from datetime import datetime
import warnings

import scanpy as sc
import xzsc_module.disease_modules as ci


def setup_logging(output_dir, mode='train'):
    """Setup logging configuration."""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = output_dir / f'disease_module_{mode}_{timestamp}.log'

    # Create logger
    logger = logging.getLogger('disease_module')
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

    return logger, timestamp


def train_and_transform(
    scanpy_h5ad_path,
    genelist_module=None,
    output_dir=None,
    use_hvg=True,
    resolution=5,
    n_neighbors=15,
    dotplot_groupby=None,
):
    """
    Train a new disease module from single-cell data and transform it.

    Parameters
    ----------
    scanpy_h5ad_path : str
        Path to the input h5ad file
    genelist_module : list, optional
        List of genes to force include in module training. Default is None.
    output_dir : str, optional
        Directory to save output files. If None, uses same directory as input h5ad.
    use_hvg : bool or int, optional
        Whether to use highly variable genes. If int, specifies number of HVGs. Default is True.
    resolution : float, optional
        Resolution parameter for Leiden clustering. Default is 5.
    n_neighbors : int, optional
        Number of neighbors for correlation calculation. Default is 15.
    dotplot_groupby : list, optional
        List of obs columns to group by for dotplots. Default is ['cluster', 'condition'].

    Returns
    -------
    module : CorrModules
        Trained module object
    adata_mod : AnnData
        Transformed anndata with module scores
    """
    # Determine output directory
    if output_dir is None:
        output_dir = Path(scanpy_h5ad_path).parent
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    logger, timestamp = setup_logging(output_dir, mode='train')

    logger.info("=" * 60)
    logger.info("Disease Module Training")
    logger.info("=" * 60)
    logger.info(f"Input h5ad: {scanpy_h5ad_path}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Resolution: {resolution}")
    logger.info(f"N neighbors: {n_neighbors}")
    logger.info(f"Use HVG: {use_hvg}")

    # Read data
    logger.info(f"Reading data from {scanpy_h5ad_path}...")
    adata = sc.read_h5ad(scanpy_h5ad_path)
    logger.info(f"Data shape: {adata.shape}")
    logger.info(f"Observations: {adata.n_obs}, Variables: {adata.n_vars}")

    # Fix log1p base if needed
    if 'log1p' in adata.uns:
        adata.uns['log1p']["base"] = None

    # Prepare module genes
    if genelist_module is not None:
        logger.info(f"Using provided gene list with {len(genelist_module)} genes")
        # Filter to genes present in data
        if hasattr(adata, 'raw') and adata.raw is not None:
            genelist_module = [k for k in genelist_module if k in adata.raw.var_names]
        else:
            genelist_module = [k for k in genelist_module if k in adata.var_names]
        logger.info(f"Found {len(genelist_module)} genes in dataset")

        # Combine with HVGs if requested
        if use_hvg:
            if 'highly_variable' not in adata.var.columns:
                n_hvg = 2000 if use_hvg is True else use_hvg
                sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
            hvg_genes = list(adata.var[adata.var.highly_variable].index)
            module_gene = list(set(hvg_genes) | set(genelist_module))
            logger.info(f"Combined HVGs ({len(hvg_genes)}) with gene list")
        else:
            module_gene = list(set(genelist_module))
    else:
        logger.info("No gene list provided, using HVGs only")
        if 'highly_variable' not in adata.var.columns:
            n_hvg = 2000 if use_hvg is True else use_hvg
            sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
        module_gene = list(adata.var[adata.var.highly_variable].index)

    logger.info(f"Training modules with {len(module_gene)} genes...")

    # Initialize and train module
    module = ci.CorrModules(adata[:, module_gene])
    module.preprocess(use_hvg=False, include=module_gene, landmarks=False)
    module.train(resolution=resolution, n_neighbors=n_neighbors)

    logger.info(f"Training complete. Number of modules: {module.n_modules}")

    # Save trained module
    module_path = output_dir / f'module_{timestamp}.pickle'
    logger.info(f"Saving trained module to {module_path}...")
    with open(module_path, 'wb') as handle:
        pickle.dump(module, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Module composition
    logger.info("Module composition statistics:")
    logger.info(f"\n{module.weights.obs[['n_genes']].describe()}")

    # Enrichment analysis
    logger.info("Running enrichment analysis...")
    try:
        module.enrich()
        weights_path = output_dir / f'moduleweights_{timestamp}.txt'
        module.weights.obs.to_csv(weights_path, sep='\t')
        logger.info(f"Module weights saved to {weights_path}")
        logger.info("\nTop 20 modules by enrichment:")
        logger.info(f"\n{module.weights.obs.head(20)}")
    except Exception as e:
        logger.warning(f"Enrichment analysis failed: {e}")

    # Transform data
    logger.info("Transforming data with trained modules...")
    adata_mod = module.transform(adata)

    # Save transformed data
    mod_path = output_dir / f'module_mod_{timestamp}.pickle'
    logger.info(f"Saving transformed data to {mod_path}...")
    with open(mod_path, 'wb') as handle:
        pickle.dump(adata_mod, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Generate dotplots
    if dotplot_groupby is None:
        # Auto-detect groupby columns
        dotplot_groupby = []
        for col in ['cluster', 'condition', 'comb_condition', 'disease']:
            if col in adata_mod.obs.columns:
                dotplot_groupby.append(col)
        # Also check for any column containing 'condition'
        for col in adata_mod.obs.columns:
            if 'condition' in col.lower() and col not in dotplot_groupby:
                dotplot_groupby.append(col)

    logger.info(f"Generating dotplots for: {dotplot_groupby}")

    # Create figures directory
    figures_dir = output_dir / 'figures'
    figures_dir.mkdir(parents=True, exist_ok=True)

    import matplotlib.pyplot as plt

    for groupby in dotplot_groupby:
        if groupby in adata_mod.obs.columns:
            try:
                logger.info(f"  - Plotting by {groupby}...")

                # Create dotplot with standard scale='var' (default)
                fig_path = figures_dir / f'dotplot_{groupby}_{timestamp}.png'
                sc.pl.dotplot(
                    adata_mod,
                    groupby=groupby,
                    var_names=adata_mod.var_names,
                    standard_scale='var',
                    show=False,
                )
                plt.savefig(fig_path, dpi=300, bbox_inches='tight')
                plt.close()
                logger.info(f"    Saved to {fig_path}")

                # Create dotplot without standard scale
                fig_path_noscale = figures_dir / f'dotplot_{groupby}_noscale_{timestamp}.png'
                sc.pl.dotplot(
                    adata_mod,
                    groupby=groupby,
                    var_names=adata_mod.var_names,
                    show=False,
                )
                plt.savefig(fig_path_noscale, dpi=300, bbox_inches='tight')
                plt.close()
                logger.info(f"    Saved unscaled version to {fig_path_noscale}")

            except Exception as e:
                logger.warning(f"Could not create dotplot for {groupby}: {e}")
        else:
            logger.warning(f"Column '{groupby}' not found in obs, skipping")

    logger.info("=" * 60)
    logger.info("Training complete!")
    logger.info("=" * 60)

    return module, adata_mod


def transform_with_pretrained(
    scanpy_h5ad_path,
    module_path,
    output_dir=None,
    dotplot_groupby=None,
):
    """
    Transform data using a pre-trained disease module.

    Parameters
    ----------
    scanpy_h5ad_path : str
        Path to the input h5ad file
    module_path : str
        Path to the pre-trained module pickle file
    output_dir : str, optional
        Directory to save output files. If None, uses same directory as input h5ad.
    dotplot_groupby : list, optional
        List of obs columns to group by for dotplots. Default is ['cluster', 'condition'].

    Returns
    -------
    adata_mod : AnnData
        Transformed anndata with module scores
    """
    # Determine output directory
    if output_dir is None:
        output_dir = Path(scanpy_h5ad_path).parent
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    logger, timestamp = setup_logging(output_dir, mode='transform')

    logger.info("=" * 60)
    logger.info("Disease Module Transformation")
    logger.info("=" * 60)
    logger.info(f"Input h5ad: {scanpy_h5ad_path}")
    logger.info(f"Pre-trained module: {module_path}")
    logger.info(f"Output directory: {output_dir}")

    # Read data
    logger.info(f"Reading data from {scanpy_h5ad_path}...")
    adata = sc.read_h5ad(scanpy_h5ad_path)
    logger.info(f"Data shape: {adata.shape}")
    logger.info(f"Observations: {adata.n_obs}, Variables: {adata.n_vars}")

    # Fix log1p base if needed
    if 'log1p' in adata.uns:
        adata.uns['log1p']["base"] = None

    # Load pre-trained module
    logger.info(f"Loading pre-trained module from {module_path}...")
    with open(module_path, 'rb') as handle:
        module = pickle.load(handle)

    logger.info(f"Pre-trained module has {module.n_modules} modules")

    # Enrichment analysis on pre-trained module
    logger.info("Running enrichment analysis on pre-trained module...")
    try:
        module.enrich()
        weights_path = output_dir / f'moduleweights_{timestamp}.txt'
        module.weights.obs.to_csv(weights_path, sep='\t')
        logger.info(f"Module weights saved to {weights_path}")
        logger.info("\nTop 20 modules by enrichment:")
        logger.info(f"\n{module.weights.obs.head(20)}")
    except Exception as e:
        logger.warning(f"Enrichment analysis failed: {e}")

    # Transform data
    logger.info("Transforming data with pre-trained modules...")
    adata_mod = module.transform(adata)

    # Save transformed data
    mod_path = output_dir / f'module_mod_{timestamp}.pickle'
    logger.info(f"Saving transformed data to {mod_path}...")
    with open(mod_path, 'wb') as handle:
        pickle.dump(adata_mod, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Generate dotplots
    if dotplot_groupby is None:
        # Auto-detect groupby columns
        dotplot_groupby = []
        for col in ['cluster', 'condition', 'comb_condition', 'disease']:
            if col in adata_mod.obs.columns:
                dotplot_groupby.append(col)
        # Also check for any column containing 'condition'
        for col in adata_mod.obs.columns:
            if 'condition' in col.lower() and col not in dotplot_groupby:
                dotplot_groupby.append(col)

    logger.info(f"Generating dotplots for: {dotplot_groupby}")

    # Create figures directory
    figures_dir = output_dir / 'figures'
    figures_dir.mkdir(parents=True, exist_ok=True)

    import matplotlib.pyplot as plt

    for groupby in dotplot_groupby:
        if groupby in adata_mod.obs.columns:
            try:
                logger.info(f"  - Plotting by {groupby}...")

                # Create dotplot with standard scale='var' (default)
                fig_path = figures_dir / f'dotplot_{groupby}_{timestamp}.png'
                sc.pl.dotplot(
                    adata_mod,
                    groupby=groupby,
                    var_names=adata_mod.var_names,
                    standard_scale='var',
                    show=False,
                )
                plt.savefig(fig_path, dpi=300, bbox_inches='tight')
                plt.close()
                logger.info(f"    Saved to {fig_path}")

                # Create dotplot without standard scale
                fig_path_noscale = figures_dir / f'dotplot_{groupby}_noscale_{timestamp}.png'
                sc.pl.dotplot(
                    adata_mod,
                    groupby=groupby,
                    var_names=adata_mod.var_names,
                    show=False,
                )
                plt.savefig(fig_path_noscale, dpi=300, bbox_inches='tight')
                plt.close()
                logger.info(f"    Saved unscaled version to {fig_path_noscale}")

            except Exception as e:
                logger.warning(f"Could not create dotplot for {groupby}: {e}")
        else:
            logger.warning(f"Column '{groupby}' not found in obs, skipping")

    logger.info("=" * 60)
    logger.info("Transformation complete!")
    logger.info("=" * 60)

    return adata_mod


def main():
    parser = argparse.ArgumentParser(
        description='Train or apply disease modules to single-cell data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Train new modules
  python disease_module.py --h5ad data.h5ad --train

  # Train with specific genes
  python disease_module.py --h5ad data.h5ad --train --genelist genes.txt

  # Use pre-trained module
  python disease_module.py --h5ad data.h5ad --pretrained module_20231215_143022.pickle

  # Specify output directory and groupby columns
  python disease_module.py --h5ad data.h5ad --train --output results/ --groupby cluster condition
        """
    )

    parser.add_argument('--h5ad', required=True, help='Path to input h5ad file')
    parser.add_argument('--train', action='store_true', help='Train new modules')
    parser.add_argument('--pretrained', help='Path to pre-trained module pickle file')
    parser.add_argument('--genelist', help='Path to file with gene list (one gene per line)')
    parser.add_argument('--output', help='Output directory (default: same as h5ad)')
    parser.add_argument('--groupby', nargs='+', help='Columns to group by for dotplots')
    parser.add_argument('--resolution', type=float, default=5, help='Leiden resolution (default: 5)')
    parser.add_argument('--n-neighbors', type=int, default=15, help='Number of neighbors (default: 15)')
    parser.add_argument('--no-hvg', action='store_true', help='Do not use HVGs')

    args = parser.parse_args()

    # Validate inputs
    if not args.train and not args.pretrained:
        parser.error("Must specify either --train or --pretrained")

    if args.train and args.pretrained:
        parser.error("Cannot specify both --train and --pretrained")

    # Load gene list if provided
    genelist = None
    if args.genelist:
        print(f"Loading gene list from {args.genelist}...")
        with open(args.genelist, 'r') as f:
            genelist = [line.strip() for line in f if line.strip()]
        print(f"Loaded {len(genelist)} genes")

    # Run analysis
    if args.train:
        train_and_transform(
            scanpy_h5ad_path=args.h5ad,
            genelist_module=genelist,
            output_dir=args.output,
            use_hvg=not args.no_hvg,
            resolution=args.resolution,
            n_neighbors=args.n_neighbors,
            dotplot_groupby=args.groupby,
        )
    else:
        transform_with_pretrained(
            scanpy_h5ad_path=args.h5ad,
            module_path=args.pretrained,
            output_dir=args.output,
            dotplot_groupby=args.groupby,
        )


if __name__ == '__main__':
    main()
