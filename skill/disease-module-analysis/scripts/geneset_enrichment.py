#!/usr/bin/env python3
"""
Gene Set Enrichment Analysis Script

Perform enrichment analysis on gene sets/lists in a dataframe.
Supports multiple enrichment methods and automatic PPI enrichment.
"""

import os
import argparse
import logging
import pickle
from pathlib import Path
from datetime import datetime
import warnings
import pandas as pd
import xzsc_module as xz


def setup_logging(output_dir):
    """Setup logging configuration."""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    date = datetime.now().strftime('%Y%m%d')
    log_file = output_dir / f'geneset_enrichment_{timestamp}.log'

    # Create logger
    logger = logging.getLogger('geneset_enrichment')
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


def geneset_annotation(
    df,
    set_column,
    species='human',
    method='gprofiler',
    pval=0.05,
    output_dir=None,
):
    """
    Perform enrichment analysis on gene sets.

    Follows the module_annotation pattern from xzsc_module.

    Parameters
    ----------
    df : DataFrame
        Input dataframe containing gene sets
    set_column : str
        Column name containing gene sets (lists)
    species : str, optional
        Species for enrichment analysis. Default is 'human'.
    method : str, optional
        Enrichment method to use. Default is 'gprofiler'.
        Options: 'gprofiler', 'enrichr', 'x2k'
    pval : float, optional
        P-value threshold for enrichment. Default is 0.05.
    output_dir : str, optional
        Directory to save output files. If None, uses current directory.

    Returns
    -------
    df_annotated : DataFrame
        DataFrame with added enrichment annotation columns
    """
    # Determine output directory
    if output_dir is None:
        output_dir = Path.cwd()
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    logger, date = setup_logging(output_dir)

    logger.info("=" * 60)
    logger.info("Gene Set Enrichment Analysis")
    logger.info("=" * 60)
    logger.info(f"Input dataframe shape: {df.shape}")
    logger.info(f"Gene set column: {set_column}")
    logger.info(f"Species: {species}")
    logger.info(f"Enrichment method: {method}")
    logger.info(f"P-value threshold: {pval}")

    # Validate input
    if set_column not in df.columns:
        raise ValueError(f"Column '{set_column}' not found in dataframe")

    # Check if column contains lists
    if not isinstance(df[set_column].iloc[0], list):
        logger.warning(f"Column '{set_column}' does not contain lists. Attempting to convert...")
        # Try to convert if it's a string representation of a list
        try:
            df[set_column] = df[set_column].apply(eval)
        except:
            raise ValueError(f"Column '{set_column}' must contain lists of genes")

    # Convert to dictionary for enrichment
    genesets = df[set_column].to_dict()
    logger.info(f"Number of gene sets: {len(genesets)}")

    # Get enrichment method
    try:
        met_func = getattr(xz.enrich, method)
    except AttributeError:
        logger.error(f"Enrichment method '{method}' not found")
        raise ValueError(f"Unknown enrichment method: {method}")

    # Perform enrichment for each gene set
    logger.info(f"\nRunning {method} enrichment analysis...")
    res = {}
    for k, v in genesets.items():
        if not isinstance(v, list):
            logger.warning(f"Gene set {k} is not a list, skipping")
            continue
        if len(v) == 0:
            logger.warning(f"Gene set {k} is empty, skipping")
            continue

        logger.info(f"  Processing gene set {k} ({len(v)} genes)...")
        try:
            res[k] = met_func(v, species=species, pval=pval, as_dict=True)
        except Exception as e:
            logger.warning(f"  Enrichment failed for gene set {k}: {e}")
            res[k] = {}

    # Extract enrichment terms
    logger.info("\nExtracting enrichment terms...")
    df_anno = xz.disease_modules._utils._modulexterm(res, prefix=method)

    # PPI enrichment
    logger.info("\nRunning PPI enrichment analysis...")
    try:
        ppi = xz.enrich.PPI(genesets, species=species)
        ppi_results = ppi.ppi_enrichment()
        df_anno['PPI'] = ppi_results['p_value'].astype(float).tolist()
        logger.info("PPI enrichment complete")
    except Exception as e:
        logger.warning(f"PPI enrichment failed: {e}")
        df_anno['PPI'] = None

    # Merge with original dataframe
    logger.info("\nMerging enrichment results with input dataframe...")
    df_annotated = df.merge(df_anno, left_index=True, right_index=True, how='left')

    # Save results as moduleweights_{date}.txt
    output_file = output_dir / f'moduleweights_{date}.txt'
    df_annotated.to_csv(output_file, sep='\t')
    logger.info(f"Results saved to {output_file}")

    # Print summary
    logger.info("\n" + "=" * 60)
    logger.info("Enrichment Summary")
    logger.info("=" * 60)
    enrichment_cols = [col for col in df_annotated.columns if col not in df.columns]
    logger.info(f"Added {len(enrichment_cols)} enrichment columns:")
    for col in enrichment_cols:
        non_null = df_annotated[col].notna().sum()
        logger.info(f"  - {col}: {non_null}/{len(df_annotated)} gene sets enriched")

    logger.info("\n" + "=" * 60)
    logger.info("Analysis complete!")
    logger.info("=" * 60)

    return df_annotated


def load_geneset_data(input_file, set_column):
    """
    Load gene set data from various formats.

    Parameters
    ----------
    input_file : str
        Path to input file (CSV, TSV, or pickle)
    set_column : str
        Column name containing gene sets

    Returns
    -------
    df : DataFrame
        Loaded dataframe
    """
    input_path = Path(input_file)

    if input_path.suffix == '.pickle':
        with open(input_path, 'rb') as f:
            df = pickle.load(f)
            # Check if it's a module object
            if hasattr(df, 'weights'):
                df = df.weights.obs
    elif input_path.suffix == '.csv':
        df = pd.read_csv(input_path, index_col=0)
    elif input_path.suffix in ['.txt', '.tsv']:
        df = pd.read_csv(input_path, sep='\t', index_col=0)
    else:
        # Try to auto-detect
        try:
            df = pd.read_csv(input_path, index_col=0)
        except:
            try:
                df = pd.read_csv(input_path, sep='\t', index_col=0)
            except:
                raise ValueError(f"Cannot read file: {input_file}")

    # Check if set_column exists
    if set_column not in df.columns:
        raise ValueError(f"Column '{set_column}' not found in {input_file}")

    return df


def main():
    parser = argparse.ArgumentParser(
        description='Perform enrichment analysis on gene sets',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with module pickle file
  python geneset_enrichment.py --input module.pickle --set-column sets

  # With custom parameters
  python geneset_enrichment.py \
    --input genesets.csv \
    --set-column gene_list \
    --species human \
    --method gprofiler \
    --pval 0.01 \
    --output results/

  # Using enrichr method
  python geneset_enrichment.py \
    --input module.pickle \
    --set-column sets \
    --method enrichr \
    --species mouse

Available enrichment methods:
  - gprofiler (default): Comprehensive pathway and ontology enrichment
  - enrichr: Enrichr database enrichment
  - x2k: Expression2Kinases upstream regulator analysis

Available species:
  - human (default)
  - mouse
  - rat
        """
    )

    parser.add_argument('--input', required=True,
                        help='Input file (CSV, TSV, or pickle with gene sets)')
    parser.add_argument('--set-column', required=True,
                        help='Column name containing gene sets (lists)')
    parser.add_argument('--species', default='human',
                        help='Species for enrichment analysis (default: human)')
    parser.add_argument('--method', default='gprofiler',
                        choices=['gprofiler', 'enrichr', 'x2k'],
                        help='Enrichment method (default: gprofiler)')
    parser.add_argument('--pval', type=float, default=0.05,
                        help='P-value threshold (default: 0.05)')
    parser.add_argument('--output', help='Output directory (default: current directory)')

    args = parser.parse_args()

    print("=" * 70)
    print("Gene Set Enrichment Analysis")
    print("=" * 70)
    print()

    # Load data
    print(f"Loading data from {args.input}...")
    df = load_geneset_data(args.input, args.set_column)
    print(f"Loaded dataframe with {len(df)} gene sets")
    print()

    # Run enrichment
    df_annotated = geneset_annotation(
        df=df,
        set_column=args.set_column,
        species=args.species,
        method=args.method,
        pval=args.pval,
        output_dir=args.output,
    )

    # Display results
    print("\nEnriched dataframe preview:")
    print(df_annotated.head())
    print()
    print(f"Output shape: {df_annotated.shape}")
    print()


if __name__ == '__main__':
    main()
