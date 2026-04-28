#!/usr/bin/env python3
"""
DEG Multi-Dataset Analysis - Table Generator

This script generates TSV tables with log2FC values and significance annotations
from multiple data sources. Tables can be used for downstream heatmap visualization.

Usage:
    python generate_tables.py --config config.json
    python generate_tables.py --mode gene_list --genes "GENE1,GENE2" --config config.json

The script supports:
- Discovery mode: Find top N genes meeting reversal criteria
- Gene List mode: Analyze specific genes regardless of significance
- Multiple data sources with different significance criteria
- Flexible score calculation based on directional logic
"""

import argparse
import json
import pandas as pd
import numpy as np
from pathlib import Path
import re
import sys


# =============================================================================
# CONFIGURATION SECTION
# =============================================================================
# These values are typically loaded from a config file, but can be customized
# directly here for specific analyses.

# Default output directory name (created in working directory)
DEFAULT_OUTPUT_DIR = "deg_multidata_output"

DEFAULT_CONFIG = {
    "mode": "gene_list",  # "discovery" or "gene_list"
    "genes": [],  # List of gene symbols (for gene_list mode)
    "top_n": 50,  # Number of top genes to return (for discovery mode)
    "data_sources": [],  # List of data source configurations
    "column_groups": {},  # Grouping of columns for organization
    "score_logic": {
        "UP": {},  # Direction contributing to score for UP table
        "DOWN": {}  # Direction contributing to score for DOWN table
    },
    "output": {
        "prefix": "analysis",
        "directory": DEFAULT_OUTPUT_DIR
    }
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def format_with_stars(log2fc, padj, log2fc_cutoff=None, missing=False):
    """
    Format log2FC value with significance stars.

    Args:
        log2fc: log2 fold change value
        padj: adjusted p-value
        log2fc_cutoff: optional |log2FC| threshold for significance
        missing: if True, gene was not found in data source

    Returns:
        Formatted string with log2FC and significance stars

    Star notation:
        **** = padj < 0.0001
        ***  = padj < 0.001
        **   = padj < 0.01
        *    = padj < 0.05
        .    = padj < 0.05 but |log2FC| below cutoff
        (none) = not significant
    """
    if missing:
        return '0'

    if pd.isna(log2fc) or pd.isna(padj):
        return '0'

    log2fc_str = f'{log2fc:.2f}'

    # Not significant
    if padj >= 0.05:
        return log2fc_str

    # Significant but below log2FC cutoff
    if log2fc_cutoff is not None and abs(log2fc) <= log2fc_cutoff:
        return f'{log2fc_str}.'

    # Significant - add stars based on p-value
    if padj < 0.0001:
        stars = '****'
    elif padj < 0.001:
        stars = '***'
    elif padj < 0.01:
        stars = '**'
    else:
        stars = '*'

    return f'{log2fc_str}{stars}'


def is_significant(formatted_value, direction=None):
    """
    Check if a formatted value is significant (has stars, not just dot).

    Args:
        formatted_value: String formatted by format_with_stars()
        direction: Optional 'up' or 'down' to check direction

    Returns:
        True if significant (and in specified direction if given)
    """
    if pd.isna(formatted_value) or formatted_value == '' or formatted_value == '0':
        return False

    # Must have stars (not just dot)
    if not re.search(r'\*+$', str(formatted_value)):
        return False

    # Extract numeric value
    numeric_str = re.sub(r'[*.]', '', str(formatted_value))
    try:
        value = float(numeric_str)
    except ValueError:
        return False

    # Check direction if specified
    if direction == 'up':
        return value > 0
    elif direction == 'down':
        return value < 0
    else:
        return True


def load_data_source(config):
    """
    Load and process a single data source.

    Args:
        config: Dictionary with data source configuration:
            - file: Path to data file
            - sheet: Sheet name (for Excel files)
            - gene_col: Column containing gene symbols
            - log2fc_col: Column containing log2FC values
            - padj_col: Column containing adjusted p-values
            - contrast_filter: Optional dict with column/value for filtering

    Returns:
        Dictionary mapping gene symbols to {log2fc, padj}
    """
    file_path = Path(config['file'])

    # Load based on file type
    if file_path.suffix in ['.xlsx', '.xls']:
        df = pd.read_excel(file_path, sheet_name=config.get('sheet'))
    elif file_path.suffix == '.tsv' or file_path.suffix == '.txt':
        df = pd.read_csv(file_path, sep='\t')
    else:
        df = pd.read_csv(file_path)

    # Apply contrast filter if specified
    if 'contrast_filter' in config and config['contrast_filter']:
        filter_config = config['contrast_filter']
        filter_col = filter_config['column']
        filter_val = filter_config['value']
        df = df[df[filter_col].str.contains(filter_val, case=False, na=False)]

    # Extract gene data
    gene_col = config['gene_col']
    log2fc_col = config['log2fc_col']
    padj_col = config['padj_col']

    gene_data = {}
    for _, row in df.iterrows():
        gene = row[gene_col]
        if pd.notna(gene):
            gene_data[gene] = {
                'log2fc': row[log2fc_col],
                'padj': row[padj_col]
            }

    return gene_data


def calculate_score(row, table_type, column_groups, score_logic, all_columns):
    """
    Calculate score based on configured logic.

    Args:
        row: DataFrame row with formatted values
        table_type: 'UP' or 'DOWN'
        column_groups: Dict mapping group names to column lists
        score_logic: Dict mapping groups to expected direction
        all_columns: List of all data columns

    Returns:
        Integer score (count of direction-matching significant values)
    """
    score = 0
    logic = score_logic.get(table_type, {})

    for group_name, columns in column_groups.items():
        expected_direction = logic.get(group_name, 'up')
        for col in columns:
            if col in row and col in all_columns:
                if is_significant(row[col], expected_direction):
                    score += 1

    return score


def count_by_group(row, table_type, column_groups, score_logic):
    """
    Count significant values by group in expected direction.

    Returns:
        Dict mapping group names to counts
    """
    counts = {}
    logic = score_logic.get(table_type, {})

    for group_name, columns in column_groups.items():
        expected_direction = logic.get(group_name, 'up')
        count = sum(1 for col in columns if col in row and is_significant(row[col], expected_direction))
        counts[group_name] = count

    return counts


def create_table(genes, data_sources, column_groups, score_logic, table_type, all_columns):
    """
    Create biomarker table for specified genes.

    Args:
        genes: List of gene symbols to include
        data_sources: Dict mapping column labels to {data, log2fc_cutoff}
        column_groups: Dict mapping group names to column lists
        score_logic: Dict with UP/DOWN score configuration
        table_type: 'UP' or 'DOWN'
        all_columns: Ordered list of all data columns

    Returns:
        DataFrame with Gene, Score, data columns, and count columns
    """
    rows = []

    for gene in genes:
        row = {'Gene': gene}

        # Add data columns
        for col_label in all_columns:
            source_data = data_sources[col_label]
            gene_data = source_data['data']
            log2fc_cutoff = source_data['log2fc_cutoff']

            if gene in gene_data:
                data = gene_data[gene]
                formatted = format_with_stars(data['log2fc'], data['padj'], log2fc_cutoff)
            else:
                formatted = '0'
            row[col_label] = formatted

        # Calculate score
        row['Score'] = calculate_score(row, table_type, column_groups, score_logic, all_columns)

        # Count by group
        counts = count_by_group(row, table_type, column_groups, score_logic)
        for group_name, count in counts.items():
            direction = score_logic[table_type].get(group_name, 'up')
            row[f'n_{group_name}_{direction}'] = count

        rows.append(row)

    # Build column order
    count_cols = [f'n_{g}_{score_logic[table_type].get(g, "up")}' for g in column_groups.keys()]
    columns = ['Gene', 'Score'] + all_columns + count_cols

    df = pd.DataFrame(rows, columns=columns)

    # Sort by Score (descending), then Gene name
    df = df.sort_values(['Score', 'Gene'], ascending=[False, True])

    return df


def discover_top_genes(data_sources, column_groups, score_logic, table_type, all_columns, top_n=50):
    """
    Find top N genes meeting criteria (Discovery mode).

    Returns:
        List of gene symbols ranked by score
    """
    # Collect all genes from all sources
    all_genes = set()
    for source_data in data_sources.values():
        all_genes.update(source_data['data'].keys())

    # Calculate scores for all genes
    gene_scores = []
    for gene in all_genes:
        row = {'Gene': gene}
        for col_label in all_columns:
            source_data = data_sources[col_label]
            gene_data = source_data['data']
            log2fc_cutoff = source_data['log2fc_cutoff']

            if gene in gene_data:
                data = gene_data[gene]
                formatted = format_with_stars(data['log2fc'], data['padj'], log2fc_cutoff)
            else:
                formatted = '0'
            row[col_label] = formatted

        score = calculate_score(row, table_type, column_groups, score_logic, all_columns)
        gene_scores.append((gene, score))

    # Sort by score and return top N
    gene_scores.sort(key=lambda x: (-x[1], x[0]))
    return [g[0] for g in gene_scores[:top_n]]


# =============================================================================
# MAIN FUNCTIONS
# =============================================================================

def load_config(config_path):
    """Load configuration from JSON file."""
    with open(config_path, 'r') as f:
        config = json.load(f)

    # Merge with defaults
    full_config = DEFAULT_CONFIG.copy()
    full_config.update(config)
    return full_config


def run_analysis(config):
    """
    Run the full analysis based on configuration.

    Args:
        config: Configuration dictionary

    Returns:
        Tuple of (up_df, down_df) DataFrames
    """
    print('DEG Multi-Dataset Analysis - Table Generator')
    print('=' * 50)

    # Load all data sources
    print('\nLoading data sources...')
    data_sources = {}
    all_columns = []

    for source_config in config['data_sources']:
        label = source_config['label']
        print(f'  Loading {label}...')

        gene_data = load_data_source(source_config)
        print(f'    {len(gene_data)} genes loaded')

        data_sources[label] = {
            'data': gene_data,
            'log2fc_cutoff': source_config.get('log2fc_cutoff'),
        }
        all_columns.append(label)

    # Get column groups
    column_groups = config['column_groups']
    score_logic = config['score_logic']

    # Determine genes to analyze
    if config['mode'] == 'gene_list':
        genes = config['genes']
        print(f'\nGene list mode: {len(genes)} genes')
    else:
        print(f'\nDiscovery mode: Finding top {config["top_n"]} genes...')
        genes = discover_top_genes(
            data_sources, column_groups, score_logic, 'UP', all_columns, config['top_n']
        )
        print(f'  Found {len(genes)} genes')

    # Create output directory
    output_dir = Path(config['output']['directory'])
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = config['output']['prefix']

    # Create UP table
    print('\nCreating UP table...')
    up_df = create_table(genes, data_sources, column_groups, score_logic, 'UP', all_columns)
    up_file = output_dir / f'{prefix}_UP_by_stim_DOWN_by_treatment.tsv'
    up_df.to_csv(up_file, sep='\t', index=False)
    print(f'  Saved {len(up_df)} genes to {up_file.name}')
    print(f'  Score range: {up_df["Score"].min()} - {up_df["Score"].max()}')

    # Create DOWN table
    print('\nCreating DOWN table...')
    down_df = create_table(genes, data_sources, column_groups, score_logic, 'DOWN', all_columns)
    down_file = output_dir / f'{prefix}_DOWN_by_stim_UP_by_treatment.tsv'
    down_df.to_csv(down_file, sep='\t', index=False)
    print(f'  Saved {len(down_df)} genes to {down_file.name}')
    print(f'  Score range: {down_df["Score"].min()} - {down_df["Score"].max()}')

    print('\nDone!')
    return up_df, down_df


def main():
    parser = argparse.ArgumentParser(
        description='Generate DEG tables from multiple data sources'
    )
    parser.add_argument('--config', required=True, help='Path to JSON config file')
    parser.add_argument('--mode', choices=['discovery', 'gene_list'],
                        help='Analysis mode (overrides config)')
    parser.add_argument('--genes', help='Comma-separated gene list (overrides config)')
    parser.add_argument('--top_n', type=int, help='Number of top genes (discovery mode)')
    parser.add_argument('--output_dir', help='Output directory (overrides config)')
    parser.add_argument('--prefix', help='Output file prefix (overrides config)')

    args = parser.parse_args()

    # Load configuration
    config = load_config(args.config)

    # Apply command-line overrides
    if args.mode:
        config['mode'] = args.mode
    if args.genes:
        config['genes'] = [g.strip() for g in args.genes.split(',')]
    if args.top_n:
        config['top_n'] = args.top_n
    if args.output_dir:
        config['output']['directory'] = args.output_dir
    if args.prefix:
        config['output']['prefix'] = args.prefix

    # Run analysis
    run_analysis(config)


if __name__ == '__main__':
    main()
