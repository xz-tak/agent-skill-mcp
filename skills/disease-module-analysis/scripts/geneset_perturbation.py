#!/usr/bin/env python3
"""
Gene Set Perturbation Enrichment Script

Perform perturbation enrichment analysis on gene sets/lists in a dataframe.
Does NOT require a trained module - works directly on gene sets.
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
    log_file = output_dir / f'geneset_perturbation_{timestamp}.log'

    # Create logger
    logger = logging.getLogger('geneset_perturbation')
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


def perturbation_enrichment_df(
    df,
    set_column,
    pert_types=['tf', 'gene'],
    species='human',
    enrich_metadata=False,
    email='user@example.com',
    output_dir=None,
):
    """
    Perform perturbation enrichment analysis on gene sets in dataframe.

    Does NOT require a module - works directly on gene sets.

    Parameters
    ----------
    df : DataFrame
        Input dataframe containing gene sets
    set_column : str
        Column name containing gene sets (lists)
    pert_types : list, optional
        List of perturbation types to test. Options: 'tf', 'gene', 'LINCS_CRISPR', 'LINCS_CHEM'
        Default is ['tf', 'gene'].
    species : str, optional
        Species for enrichment analysis. Default is 'human'.
    enrich_metadata : bool, optional
        Whether to enrich results with NCBI metadata. Default is False.
    email : str, optional
        Email for NCBI API (required if enrich_metadata=True). Default is 'user@example.com'.
    output_dir : str, optional
        Directory to save output files. If None, uses current directory.

    Returns
    -------
    results : dict
        Dictionary containing perturbation enrichment results for each type
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
    logger.info("Gene Set Perturbation Enrichment")
    logger.info("=" * 60)
    logger.info(f"Input dataframe shape: {df.shape}")
    logger.info(f"Gene set column: {set_column}")
    logger.info(f"Perturbation types: {pert_types}")
    logger.info(f"Species: {species}")
    logger.info(f"Enrich metadata: {enrich_metadata}")

    # Validate input
    if set_column not in df.columns:
        raise ValueError(f"Column '{set_column}' not found in dataframe")

    # Check if column contains lists
    if not isinstance(df[set_column].iloc[0], list):
        logger.warning(f"Column '{set_column}' does not contain lists. Attempting to convert...")
        try:
            df[set_column] = df[set_column].apply(eval)
        except:
            raise ValueError(f"Column '{set_column}' must contain lists of genes")

    # Convert to dictionary for enrichment
    genesets = df[set_column].to_dict()
    logger.info(f"Number of gene sets: {len(genesets)}")

    results = {}

    # Process each perturbation type
    for pert_type in pert_types:
        logger.info(f"\n{'=' * 60}")
        logger.info(f"Processing {pert_type} perturbation enrichment")
        logger.info(f"{'=' * 60}")

        try:
            # Get perturbation enrichment function
            if pert_type == 'tf':
                pert_func = xz.enrich.pert_TF
                pert_col = 'pert_TF'
            elif pert_type == 'gene':
                pert_func = xz.enrich.pert_GENE
                pert_col = 'pert_GENE'
            elif pert_type == 'LINCS_CRISPR':
                pert_func = xz.enrich.pert_LINCS_CRISPR
                pert_col = 'pert_GENE'
            elif pert_type == 'LINCS_CHEM':
                pert_func = xz.enrich.pert_LINCS_CHEM
                pert_col = 'pert_CHEM'
            else:
                logger.warning(f"Unknown perturbation type: {pert_type}, skipping")
                continue

            # Perform enrichment for each gene set
            all_results = []
            for idx, geneset in genesets.items():
                if not isinstance(geneset, list) or len(geneset) == 0:
                    logger.warning(f"Gene set {idx} is invalid, skipping")
                    continue

                logger.info(f"  Processing gene set {idx} ({len(geneset)} genes)...")
                try:
                    pert_result = pert_func(geneset, species=species, lookup=False, as_dataframe=True)

                    if pert_result is not None and len(pert_result) > 0:
                        # Add gene set identifier
                        pert_result['geneset_id'] = idx
                        all_results.append(pert_result)
                        logger.info(f"    Found {len(pert_result)} significant perturbations")
                except Exception as e:
                    logger.warning(f"  Perturbation enrichment failed for gene set {idx}: {e}")

            # Combine results
            if all_results:
                combined_df = pd.concat(all_results, ignore_index=True)
                logger.info(f"\nTotal {pert_type} perturbations found: {len(combined_df)}")

                # Enrich with NCBI metadata if requested
                if enrich_metadata and NCBI_AVAILABLE:
                    logger.info(f"Enriching {pert_type} perturbations with NCBI metadata...")
                    combined_df = enrich_perturbation_tables(combined_df, pert_type=pert_type, email=email)

                # Save results
                output_file = output_dir / f'perturbation_{pert_type}_{date}.txt'
                combined_df.to_csv(output_file, sep='\t', index=False)
                logger.info(f"Saved to {output_file}")

                results[pert_type] = combined_df
            else:
                logger.warning(f"No {pert_type} perturbations found")
                results[pert_type] = None

        except Exception as e:
            logger.error(f"Failed to process {pert_type} perturbation: {e}")
            results[pert_type] = None

    logger.info("\n" + "=" * 60)
    logger.info("Perturbation enrichment complete!")
    logger.info("=" * 60)

    return results


def find_regulators_df(
    df,
    set_column,
    target_genes,
    pert_types=['tf', 'gene'],
    species='human',
    enrich_metadata=False,
    email='user@example.com',
    output_dir=None,
):
    """
    Find regulators for target genes using perturbation enrichment on gene sets.

    Does NOT require a module - works directly on gene sets.

    Parameters
    ----------
    df : DataFrame
        Input dataframe containing gene sets
    set_column : str
        Column name containing gene sets (lists)
    target_genes : list
        List of target gene symbols to find regulators for
    pert_types : list, optional
        List of perturbation types to test. Options: 'tf', 'gene', 'LINCS_CRISPR', 'LINCS_CHEM'
        Default is ['tf', 'gene'].
    species : str, optional
        Species for enrichment analysis. Default is 'human'.
    enrich_metadata : bool, optional
        Whether to enrich results with NCBI metadata. Default is False.
    email : str, optional
        Email for NCBI API (required if enrich_metadata=True). Default is 'user@example.com'.
    output_dir : str, optional
        Directory to save output files. If None, uses current directory.

    Returns
    -------
    regulators : dict
        Dictionary containing regulator results for each perturbation type
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
    logger.info("Find Regulators for Target Genes")
    logger.info("=" * 60)
    logger.info(f"Input dataframe shape: {df.shape}")
    logger.info(f"Gene set column: {set_column}")
    logger.info(f"Target genes: {target_genes}")
    logger.info(f"Perturbation types: {pert_types}")
    logger.info(f"Species: {species}")

    # Validate input
    if set_column not in df.columns:
        raise ValueError(f"Column '{set_column}' not found in dataframe")

    # Check if column contains lists
    if not isinstance(df[set_column].iloc[0], list):
        logger.warning(f"Column '{set_column}' does not contain lists. Attempting to convert...")
        try:
            df[set_column] = df[set_column].apply(eval)
        except:
            raise ValueError(f"Column '{set_column}' must contain lists of genes")

    # Convert to dictionary for enrichment
    genesets = df[set_column].to_dict()
    logger.info(f"Number of gene sets: {len(genesets)}")

    target_genes_set = set(target_genes)
    regulators = {}

    # Process each perturbation type
    for pert_type in pert_types:
        logger.info(f"\n{'=' * 60}")
        logger.info(f"Finding {pert_type} regulators for target genes")
        logger.info(f"{'=' * 60}")

        try:
            # Get perturbation enrichment function
            if pert_type == 'tf':
                pert_func = xz.enrich.pert_TF
                pert_col = 'pert_TF'
            elif pert_type == 'gene':
                pert_func = xz.enrich.pert_GENE
                pert_col = 'pert_GENE'
            elif pert_type == 'LINCS_CRISPR':
                pert_func = xz.enrich.pert_LINCS_CRISPR
                pert_col = 'pert_GENE'
            elif pert_type == 'LINCS_CHEM':
                pert_func = xz.enrich.pert_LINCS_CHEM
                pert_col = 'pert_CHEM'
            else:
                logger.warning(f"Unknown perturbation type: {pert_type}, skipping")
                continue

            # Perform enrichment for each gene set
            all_regulators = []
            for idx, geneset in genesets.items():
                if not isinstance(geneset, list) or len(geneset) == 0:
                    logger.warning(f"Gene set {idx} is invalid, skipping")
                    continue

                logger.info(f"  Processing gene set {idx} ({len(geneset)} genes)...")
                try:
                    pert_result = pert_func(geneset, species=species, lookup=False, as_dataframe=True)

                    if pert_result is not None and len(pert_result) > 0:
                        # Filter for target genes in overlap_genes
                        if 'overlap_genes' in pert_result.columns:
                            mask = pert_result['overlap_genes'].apply(
                                lambda x: not set(x).isdisjoint(target_genes_set)
                            )
                            target_regs = pert_result[mask].copy()

                            if len(target_regs) > 0:
                                # Add gene set identifier
                                target_regs['geneset_id'] = idx
                                all_regulators.append(target_regs)
                                logger.info(f"    Found {len(target_regs)} regulators affecting target genes")
                except Exception as e:
                    logger.warning(f"  Regulator search failed for gene set {idx}: {e}")

            # Combine results
            if all_regulators:
                combined_df = pd.concat(all_regulators, ignore_index=True)
                logger.info(f"\nTotal {pert_type} regulators found: {len(combined_df)}")

                # Enrich with NCBI metadata if requested
                if enrich_metadata and NCBI_AVAILABLE:
                    logger.info(f"Enriching {pert_type} regulators with NCBI metadata...")
                    combined_df = enrich_perturbation_tables(combined_df, pert_type=pert_type, email=email)

                # Save results
                output_file = output_dir / f'regulators_{pert_type}_{date}.txt'
                combined_df.to_csv(output_file, sep='\t', index=False)
                logger.info(f"Saved to {output_file}")

                regulators[pert_type] = combined_df
            else:
                logger.warning(f"No {pert_type} regulators found for target genes")
                regulators[pert_type] = None

        except Exception as e:
            logger.error(f"Failed to process {pert_type} regulators: {e}")
            regulators[pert_type] = None

    logger.info("\n" + "=" * 60)
    logger.info("Regulator search complete!")
    logger.info("=" * 60)

    return regulators


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
            obj = pickle.load(f)
            # Check if it's a module object
            if hasattr(obj, 'weights'):
                df = obj.weights.obs
            else:
                df = obj
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
        description='Perform perturbation enrichment on gene sets (no module required)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Perturbation enrichment on gene sets
  python geneset_perturbation.py \\
    --input genesets.csv \\
    --set-column gene_list \\
    --mode perturb \\
    --pert-types tf gene

  # With NCBI metadata enrichment
  python geneset_perturbation.py \\
    --input genesets.csv \\
    --set-column gene_list \\
    --mode perturb \\
    --pert-types tf gene LINCS_CRISPR \\
    --enrich-metadata \\
    --email your.email@example.com

  # Find regulators for target genes
  python geneset_perturbation.py \\
    --input genesets.csv \\
    --set-column gene_list \\
    --mode regulators \\
    --target-genes targets.txt \\
    --pert-types tf gene

Available perturbation types:
  - tf: Transcription factor perturbations
  - gene: Gene perturbations
  - LINCS_CRISPR: LINCS CRISPR knockout perturbations
  - LINCS_CHEM: LINCS chemical perturbations
        """
    )

    parser.add_argument('--input', required=True,
                        help='Input file (CSV, TSV, or pickle with gene sets)')
    parser.add_argument('--set-column', required=True,
                        help='Column name containing gene sets (lists)')
    parser.add_argument('--mode', required=True, choices=['perturb', 'regulators'],
                        help='Analysis mode: perturb or regulators')
    parser.add_argument('--pert-types', nargs='+', default=['tf', 'gene'],
                        choices=['tf', 'gene', 'LINCS_CRISPR', 'LINCS_CHEM'],
                        help='Perturbation types to test (default: tf gene)')
    parser.add_argument('--target-genes', help='Path to file with target genes (one per line) for regulators mode')
    parser.add_argument('--species', default='human',
                        help='Species for enrichment (default: human)')
    parser.add_argument('--enrich-metadata', action='store_true',
                        help='Enrich results with NCBI metadata')
    parser.add_argument('--email', default='user@example.com',
                        help='Email for NCBI API (required if using --enrich-metadata)')
    parser.add_argument('--output', help='Output directory (default: current directory)')

    args = parser.parse_args()

    print("=" * 70)
    print("Gene Set Perturbation Enrichment (No Module Required)")
    print("=" * 70)
    print()

    # Load data
    print(f"Loading data from {args.input}...")
    df = load_geneset_data(args.input, args.set_column)
    print(f"Loaded dataframe with {len(df)} gene sets")
    print()

    if args.mode == 'perturb':
        # Perturbation enrichment
        results = perturbation_enrichment_df(
            df=df,
            set_column=args.set_column,
            pert_types=args.pert_types,
            species=args.species,
            enrich_metadata=args.enrich_metadata,
            email=args.email,
            output_dir=args.output,
        )

        # Display summary
        print("\nResults summary:")
        for pert_type, result_df in results.items():
            if result_df is not None:
                print(f"  {pert_type}: {len(result_df)} perturbations")
            else:
                print(f"  {pert_type}: No results")

    elif args.mode == 'regulators':
        if not args.target_genes:
            parser.error("--target-genes required for regulators mode")

        # Load target genes
        with open(args.target_genes, 'r') as f:
            target_genes = [line.strip() for line in f if line.strip()]

        # Find regulators
        regulators = find_regulators_df(
            df=df,
            set_column=args.set_column,
            target_genes=target_genes,
            pert_types=args.pert_types,
            species=args.species,
            enrich_metadata=args.enrich_metadata,
            email=args.email,
            output_dir=args.output,
        )

        # Display summary
        print("\nRegulators summary:")
        for pert_type, result_df in regulators.items():
            if result_df is not None:
                print(f"  {pert_type}: {len(result_df)} regulators")
            else:
                print(f"  {pert_type}: No regulators found")

    print()


if __name__ == '__main__':
    main()
