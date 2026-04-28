#!/usr/bin/env python3
"""
DEG Multi-Dataset Analysis - Pathway Table Generator

This script generates TSV tables with NES values and significance annotations
from GSEA results across multiple data sources. Tables can be used for
downstream heatmap visualization.

Usage:
    python generate_pathway_tables.py --config config.json
    python generate_pathway_tables.py --gsea_file gsea_all.txt --config config.json

The script supports:
- Reading combined GSEA results from run_gsea.R
- Scoring pathways by reversal pattern
- Formatting NES with significance stars
- Generating UP/DOWN tables with column groups
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

# Default output directory name (created in working directory)
DEFAULT_OUTPUT_DIR = "deg_multidata_output"

DEFAULT_CONFIG = {
    "mode": "pathway_analysis",
    "top_n": 100,
    "column_groups": {},
    "score_logic": {
        "UP": {},
        "DOWN": {}
    },
    "output": {
        "prefix": "analysis",
        "directory": DEFAULT_OUTPUT_DIR
    }
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def format_nes_with_stars(nes, padj, missing=False):
    """
    Format NES value with significance stars.

    Args:
        nes: Normalized Enrichment Score
        padj: adjusted p-value
        missing: if True, pathway was not found for this data source

    Returns:
        Formatted string with NES and significance stars

    Star notation:
        **** = padj < 0.0001
        ***  = padj < 0.001
        **   = padj < 0.01
        *    = padj < 0.05
        (none) = not significant
    """
    if missing:
        return '0'

    if pd.isna(nes) or pd.isna(padj):
        return '0'

    nes_str = f'{nes:.2f}'

    # Not significant
    if padj >= 0.05:
        return nes_str

    # Significant - add stars based on p-value
    if padj < 0.0001:
        stars = '****'
    elif padj < 0.001:
        stars = '***'
    elif padj < 0.01:
        stars = '**'
    else:
        stars = '*'

    return f'{nes_str}{stars}'


def is_significant(formatted_value, direction=None):
    """
    Check if a formatted value is significant (has stars).

    Args:
        formatted_value: String formatted by format_nes_with_stars()
        direction: Optional 'up' (positive NES) or 'down' (negative NES)

    Returns:
        True if significant (and in specified direction if given)
    """
    if pd.isna(formatted_value) or formatted_value == '' or formatted_value == '0':
        return False

    # Must have stars
    if not re.search(r'\*+$', str(formatted_value)):
        return False

    # Extract numeric value
    numeric_str = re.sub(r'\*+$', '', str(formatted_value))
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


def load_gsea_results(gsea_file):
    """
    Load combined GSEA results from run_gsea.R output.

    Args:
        gsea_file: Path to GSEA results file (TSV format)

    Returns:
        Dictionary mapping source labels to {pathway: {NES, p.adjust}}
    """
    df = pd.read_csv(gsea_file, sep='\t')

    # Group by source
    source_data = {}
    for source in df['source'].unique():
        source_df = df[df['source'] == source]
        pathway_data = {}
        for _, row in source_df.iterrows():
            pathway_id = row['ID']
            pathway_data[pathway_id] = {
                'NES': row['NES'],
                'padj': row['p.adjust']
            }
        source_data[source] = pathway_data

    return source_data


def calculate_score(row, table_type, column_groups, score_logic, all_columns):
    """
    Calculate score based on configured logic.

    Args:
        row: Dictionary with formatted values
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


def discover_top_pathways(source_data, column_groups, score_logic, table_type, all_columns, top_n=100):
    """
    Find top N pathways meeting criteria.

    Returns:
        List of pathway IDs ranked by score
    """
    # Collect all pathways from all sources
    all_pathways = set()
    for source_pathways in source_data.values():
        all_pathways.update(source_pathways.keys())

    # Calculate scores for all pathways
    pathway_scores = []
    for pathway in all_pathways:
        row = {'Pathway': pathway}
        for col_label in all_columns:
            if col_label in source_data:
                pathway_data = source_data[col_label]
                if pathway in pathway_data:
                    data = pathway_data[pathway]
                    formatted = format_nes_with_stars(data['NES'], data['padj'])
                else:
                    formatted = '0'
            else:
                formatted = '0'
            row[col_label] = formatted

        score = calculate_score(row, table_type, column_groups, score_logic, all_columns)
        pathway_scores.append((pathway, score))

    # Sort by score and return top N
    pathway_scores.sort(key=lambda x: (-x[1], x[0]))
    return [p[0] for p in pathway_scores[:top_n]]


def create_pathway_table(pathways, source_data, column_groups, score_logic, table_type, all_columns):
    """
    Create biomarker table for specified pathways.

    Args:
        pathways: List of pathway IDs to include
        source_data: Dict mapping source labels to {pathway: {NES, padj}}
        column_groups: Dict mapping group names to column lists
        score_logic: Dict with UP/DOWN score configuration
        table_type: 'UP' or 'DOWN'
        all_columns: Ordered list of all data columns

    Returns:
        DataFrame with Pathway, Score, data columns, and count columns
    """
    rows = []

    for pathway in pathways:
        row = {'Pathway': pathway}

        # Add data columns
        for col_label in all_columns:
            if col_label in source_data:
                pathway_data = source_data[col_label]
                if pathway in pathway_data:
                    data = pathway_data[pathway]
                    formatted = format_nes_with_stars(data['NES'], data['padj'])
                else:
                    formatted = '0'
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
    columns = ['Pathway', 'Score'] + all_columns + count_cols

    df = pd.DataFrame(rows, columns=columns)

    # Sort by Score (descending), then Pathway name
    df = df.sort_values(['Score', 'Pathway'], ascending=[False, True])

    return df


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


def run_analysis(config, gsea_file=None):
    """
    Run the full pathway analysis based on configuration.

    Args:
        config: Configuration dictionary
        gsea_file: Optional path to GSEA results file

    Returns:
        Tuple of (up_df, down_df) DataFrames
    """
    print('DEG Multi-Dataset Analysis - Pathway Table Generator')
    print('=' * 55)

    # Determine GSEA file path
    if gsea_file is None:
        output_dir = Path(config['output']['directory'])
        prefix = config['output']['prefix']
        gsea_file = output_dir / 'gsea' / f'{prefix}_gsea_all.txt'

    print(f'\nLoading GSEA results: {gsea_file}')
    source_data = load_gsea_results(gsea_file)
    print(f'  Loaded {len(source_data)} data sources')

    # Get column order from column_groups
    column_groups = config['column_groups']
    all_columns = []
    for group_cols in column_groups.values():
        all_columns.extend(group_cols)

    # Verify all columns are in source_data
    available_sources = set(source_data.keys())
    missing_sources = [c for c in all_columns if c not in available_sources]
    if missing_sources:
        print(f'  Warning: Missing sources in GSEA data: {missing_sources}')
        all_columns = [c for c in all_columns if c in available_sources]

    score_logic = config['score_logic']
    top_n = config.get('top_n', 100)

    # Discover top pathways for UP table
    print(f'\nDiscovering top {top_n} pathways for UP table...')
    up_pathways = discover_top_pathways(
        source_data, column_groups, score_logic, 'UP', all_columns, top_n
    )
    print(f'  Found {len(up_pathways)} pathways')

    # Discover top pathways for DOWN table
    print(f'\nDiscovering top {top_n} pathways for DOWN table...')
    down_pathways = discover_top_pathways(
        source_data, column_groups, score_logic, 'DOWN', all_columns, top_n
    )
    print(f'  Found {len(down_pathways)} pathways')

    # Create output directory
    output_dir = Path(config['output']['directory'])
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = config['output']['prefix']

    # Create UP table
    print('\nCreating UP table...')
    up_df = create_pathway_table(up_pathways, source_data, column_groups, score_logic, 'UP', all_columns)
    up_file = output_dir / f'{prefix}_UP_by_stim_DOWN_by_treatment.tsv'
    up_df.to_csv(up_file, sep='\t', index=False)
    print(f'  Saved {len(up_df)} pathways to {up_file.name}')
    print(f'  Score range: {up_df["Score"].min()} - {up_df["Score"].max()}')

    # Create DOWN table
    print('\nCreating DOWN table...')
    down_df = create_pathway_table(down_pathways, source_data, column_groups, score_logic, 'DOWN', all_columns)
    down_file = output_dir / f'{prefix}_DOWN_by_stim_UP_by_treatment.tsv'
    down_df.to_csv(down_file, sep='\t', index=False)
    print(f'  Saved {len(down_df)} pathways to {down_file.name}')
    print(f'  Score range: {down_df["Score"].min()} - {down_df["Score"].max()}')

    print('\nDone!')
    return up_df, down_df


def main():
    parser = argparse.ArgumentParser(
        description='Generate pathway tables from GSEA results'
    )
    parser.add_argument('--config', required=True, help='Path to JSON config file')
    parser.add_argument('--gsea_file', help='Path to GSEA results file (overrides config)')
    parser.add_argument('--top_n', type=int, help='Number of top pathways (overrides config)')
    parser.add_argument('--output_dir', help='Output directory (overrides config)')
    parser.add_argument('--prefix', help='Output file prefix (overrides config)')

    args = parser.parse_args()

    # Load configuration
    config = load_config(args.config)

    # Apply command-line overrides
    if args.top_n:
        config['top_n'] = args.top_n
    if args.output_dir:
        config['output']['directory'] = args.output_dir
    if args.prefix:
        config['output']['prefix'] = args.prefix

    # Run analysis
    run_analysis(config, args.gsea_file)


if __name__ == '__main__':
    main()
