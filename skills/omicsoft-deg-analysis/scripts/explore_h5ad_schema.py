#!/usr/bin/env python3
"""
Explore H5AD Schema
Generates a comprehensive report of available filter values in an Omicsoft h5ad file

Usage:
    python explore_h5ad_schema.py --file <path_to_h5ad> [--output schema_report.json] [--format json]
"""

import argparse
import sys
import os
import json

import scanpy as sc
import pandas as pd
import numpy as np


def load_h5ad_safe(file_path):
    """Load h5ad file with automatic cleanup of expression data"""
    print(f"Loading h5ad file: {file_path}")

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    adata = sc.read_h5ad(file_path)
    print(f"✓ Loaded AnnData: {adata.shape[0]:,} observations × {adata.shape[1]:,} variables")

    # Clean up expression data from uns if present
    base_name = os.path.basename(file_path).replace('.h5ad', '')

    if base_name.endswith('_deg'):
        potential_expr_key = base_name.replace('_deg', '_expr')

        if potential_expr_key in adata.uns:
            print(f"ℹ Removing expression data from uns['{potential_expr_key}']...")
            del adata.uns[potential_expr_key]
            print(f"✓ Memory optimized")

    return adata


def get_column_stats(series, column_name):
    """Get statistics for a column"""
    stats = {
        'column': column_name,
        'total_values': len(series),
        'unique_values': series.nunique(),
        'null_count': series.isna().sum(),
        'dtype': str(series.dtype)
    }
    return stats


def format_value_counts(series, max_display=50):
    """Format value counts with counts and percentages"""
    value_counts = series.value_counts()
    total = len(series)

    lines = []
    for i, (value, count) in enumerate(value_counts.items()):
        if i >= max_display:
            remaining = len(value_counts) - max_display
            lines.append(f"    ... and {remaining} more values")
            break

        pct = (count / total) * 100
        lines.append(f"    - {value}: {count:,} ({pct:.1f}%)")

    return '\n'.join(lines)


def search_values(series, search_term, column_name):
    """Search for values matching a search term"""
    search_term_lower = search_term.lower()
    matching_values = [v for v in series.unique() if pd.notna(v) and search_term_lower in str(v).lower()]

    if matching_values:
        result = f"\n  Matches for '{search_term}' in {column_name}:\n"
        for val in sorted(matching_values):
            count = (series == val).sum()
            result += f"    - {val}: {count:,} observations\n"
        return result
    else:
        return f"\n  No matches for '{search_term}' in {column_name}\n"


def get_all_value_counts(series):
    """Get complete value counts for a column"""
    value_counts = series.value_counts()
    total = len(series)

    result = []
    for value, count in value_counts.items():
        pct = (count / total) * 100
        result.append({
            'value': str(value) if pd.notna(value) else 'NULL',
            'count': int(count),
            'percentage': round(pct, 2)
        })

    return result


def explore_schema_json(file_path):
    """Explore schema and return JSON structure"""

    # Load data
    adata = load_h5ad_safe(file_path)

    # Main schema structure
    schema = {
        'file_info': {
            'file_path': file_path,
            'n_observations': int(adata.n_obs),
            'n_variables': int(adata.n_vars),
            'observation_description': 'Each observation represents a differential expression comparison (e.g., disease vs. control in a study)',
            'variable_description': 'Each variable represents a gene'
        },
        'obs_columns': {},
        'var_info': {
            'n_genes': int(adata.n_vars),
            'gene_list_sample': list(adata.var_names[:100]),
            'all_genes_available': True
        },
        'layers': {},
        'filtering_guide': {
            'exact_match_filters': ['comparison_category', 'study', 'case_treatment', 'control_treatment', 'case_treatment_status', 'control_treatment_status'],
            'partial_match_filters': ['tissue', 'disease', 'comparison'],
            'note': 'Exact match filters require exact string matching. Partial match filters are case-insensitive substring matching.'
        }
    }

    # Define column categories
    column_categories = {
        'key_filtering': [
            'comparison_category', 'tissue', 'disease', 'disease_category',
            'case_treatment', 'control_treatment',
            'case_treatment_status', 'control_treatment_status',
            'comparison', 'study', 'database'
        ],
        'demographic': [
            'case_age_category', 'case_gender', 'case_ethnicity',
            'control_age_category', 'control_gender', 'control_ethnicity'
        ],
        'case_disease': [
            'case_disease_state', 'case_disease_subtype',
            'case_disease_group', 'case_disease_location'
        ],
        'control_disease': [
            'control_disease_state', 'control_disease_subtype',
            'control_disease_group', 'control_disease_location'
        ],
        'treatment': [
            'case_dosage', 'case_treatment_group', 'case_treat_time',
            'control_dosage', 'control_treatment_group', 'control_treat_time'
        ],
        'response': [
            'case_response', 'control_response'
        ],
        'sample': [
            'case_tissue', 'case_sample_material',
            'control_tissue', 'control_sample_material',
            'sample', 'project_id'
        ],
        'comparison_details': [
            'comparison_id', 'comparison_contrast'
        ]
    }

    # Process all obs columns
    for col in adata.obs.columns:
        col_data = {
            'total_values': int(len(adata.obs[col])),
            'unique_values': int(adata.obs[col].nunique()),
            'null_count': int(adata.obs[col].isna().sum()),
            'dtype': str(adata.obs[col].dtype),
            'all_values': get_all_value_counts(adata.obs[col])
        }

        # Add category
        col_data['category'] = 'other'
        for cat_name, cols in column_categories.items():
            if col in cols:
                col_data['category'] = cat_name
                break

        schema['obs_columns'][col] = col_data

    # Process layers
    for layer_name in adata.layers.keys():
        layer = adata.layers[layer_name]

        if hasattr(layer, 'toarray'):
            layer_dense = layer.toarray()
        else:
            layer_dense = layer

        non_zero = int(np.count_nonzero(layer_dense))
        total_elements = int(layer_dense.size)
        sparsity = round((1 - non_zero / total_elements) * 100, 2)

        schema['layers'][layer_name] = {
            'shape': list(layer.shape),
            'type': type(layer).__name__,
            'non_zero_elements': non_zero,
            'total_elements': total_elements,
            'sparsity_percent': sparsity,
            'min': float(layer_dense.min()),
            'max': float(layer_dense.max()),
            'mean': float(layer_dense.mean())
        }

    # Add categorized column lists
    schema['column_categories'] = {}
    for cat_name, cols in column_categories.items():
        existing_cols = [col for col in cols if col in adata.obs.columns]
        schema['column_categories'][cat_name] = existing_cols

    return schema


def explore_schema(file_path, output_file=None, search_terms=None, output_format='text'):
    """Main exploration function"""

    if output_format == 'json':
        # Generate JSON schema
        schema = explore_schema_json(file_path)

        if output_file:
            with open(output_file, 'w') as f:
                json.dump(schema, f, indent=2)
            print(f"\n✓ JSON schema report saved to: {output_file}")
        else:
            print(json.dumps(schema, indent=2))

        return schema

    # Load data for text format
    adata = load_h5ad_safe(file_path)

    # Prepare output
    lines = []
    lines.append("=" * 100)
    lines.append(f"H5AD FILE SCHEMA REPORT")
    lines.append("=" * 100)
    lines.append(f"\nFile: {file_path}")
    lines.append(f"Observations (comparisons): {adata.n_obs:,}")
    lines.append(f"Variables (genes): {adata.n_vars:,}")

    # Observation metadata columns
    lines.append("\n" + "=" * 100)
    lines.append("OBSERVATION METADATA (obs)")
    lines.append("=" * 100)
    lines.append(f"\nTotal columns: {len(adata.obs.columns)}")
    lines.append(f"Available columns: {', '.join(adata.obs.columns)}")

    # Key filtering columns
    key_columns = [
        'comparison_category',
        'tissue',
        'disease',
        'disease_category',
        'case_treatment',
        'control_treatment',
        'case_treatment_status',
        'control_treatment_status',
        'comparison',
        'study',
        'database'
    ]

    # Filter to columns that exist
    existing_key_columns = [col for col in key_columns if col in adata.obs.columns]

    lines.append("\n" + "-" * 100)
    lines.append("KEY FILTERING COLUMNS")
    lines.append("-" * 100)

    for col in existing_key_columns:
        lines.append(f"\n[{col}]")
        stats = get_column_stats(adata.obs[col], col)
        lines.append(f"  Total: {stats['total_values']:,} | Unique: {stats['unique_values']:,} | Null: {stats['null_count']:,} | Type: {stats['dtype']}")
        lines.append(f"\n  All values:")
        lines.append(format_value_counts(adata.obs[col]))

    # Additional metadata columns
    other_columns = [col for col in adata.obs.columns if col not in existing_key_columns]

    if other_columns:
        lines.append("\n" + "-" * 100)
        lines.append("OTHER METADATA COLUMNS")
        lines.append("-" * 100)

        for col in sorted(other_columns):
            stats = get_column_stats(adata.obs[col], col)
            lines.append(f"\n[{col}]")
            lines.append(f"  Total: {stats['total_values']:,} | Unique: {stats['unique_values']:,} | Null: {stats['null_count']:,} | Type: {stats['dtype']}")

            # Only show value counts for categorical columns with reasonable number of unique values
            if stats['unique_values'] <= 100:
                lines.append(f"\n  Values:")
                lines.append(format_value_counts(adata.obs[col], max_display=30))

    # Variable metadata
    lines.append("\n" + "=" * 100)
    lines.append("VARIABLE METADATA (var)")
    lines.append("=" * 100)
    lines.append(f"\nTotal genes: {adata.n_vars:,}")
    lines.append(f"Variable columns: {', '.join(adata.var.columns) if len(adata.var.columns) > 0 else 'None'}")

    if len(adata.var.columns) > 0:
        lines.append("\nFirst 10 genes:")
        for i, gene in enumerate(adata.var_names[:10]):
            lines.append(f"  {i+1}. {gene}")

    # Layers
    lines.append("\n" + "=" * 100)
    lines.append("LAYERS (expression matrices)")
    lines.append("=" * 100)
    lines.append(f"\nAvailable layers: {', '.join(adata.layers.keys())}")

    for layer_name in adata.layers.keys():
        layer = adata.layers[layer_name]
        lines.append(f"\n[{layer_name}]")
        lines.append(f"  Shape: {layer.shape}")
        lines.append(f"  Type: {type(layer).__name__}")

        # Get some statistics
        if hasattr(layer, 'toarray'):
            layer_dense = layer.toarray()
        else:
            layer_dense = layer

        non_zero = np.count_nonzero(layer_dense)
        total_elements = layer_dense.size
        sparsity = (1 - non_zero / total_elements) * 100

        lines.append(f"  Non-zero elements: {non_zero:,} / {total_elements:,} ({100-sparsity:.1f}%)")
        lines.append(f"  Min: {layer_dense.min():.4f} | Max: {layer_dense.max():.4f} | Mean: {layer_dense.mean():.4f}")

    # Search functionality
    if search_terms:
        lines.append("\n" + "=" * 100)
        lines.append("SEARCH RESULTS")
        lines.append("=" * 100)

        for search_term in search_terms:
            lines.append(f"\nSearching for: '{search_term}'")
            lines.append("-" * 100)

            for col in existing_key_columns:
                result = search_values(adata.obs[col], search_term, col)
                lines.append(result)

    # Quick reference guide
    lines.append("\n" + "=" * 100)
    lines.append("QUICK REFERENCE: HOW TO USE THESE VALUES")
    lines.append("=" * 100)
    lines.append("""
When running deg_analysis.py, use the exact values shown above for filtering:

1. --comparison-category: Exact match, comma-separated
   Example: --comparison-category "Disease vs. Normal,Treatment vs. Control"

2. --tissue: Partial match (case-insensitive), comma-separated
   Example: --tissues "skin,blood,lung"

3. --diseases: Partial match (case-insensitive), comma-separated
   Example: --diseases "scleroderma,sclerosis,crohn,colitis"

4. --case-treatment: Exact match, comma-separated
   Example: --case-treatment "none,NA"

5. --comparison: Partial match (case-insensitive), comma-separated
   Example: --comparison "response vs no response"

6. --studies: Exact match, comma-separated
   Example: --studies "GSE130955,GSE181549"

Note: All filters are optional. If not specified, all data will be included.
""")

    # Output results
    report = '\n'.join(lines)

    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)
        print(f"\n✓ Schema report saved to: {output_file}")
    else:
        print("\n" + report)

    return report


def main():
    parser = argparse.ArgumentParser(
        description='Explore h5ad schema and generate filter reference',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--file',
        required=True,
        help='Path to h5ad file'
    )

    parser.add_argument(
        '--output',
        default=None,
        help='Output file for schema report (default: print to stdout)'
    )

    parser.add_argument(
        '--format',
        choices=['text', 'json'],
        default='text',
        help='Output format: text or json (default: text)'
    )

    parser.add_argument(
        '--search',
        nargs='+',
        default=None,
        help='Search terms to find matching values (e.g., --search scleroderma ibd skin) [text format only]'
    )

    args = parser.parse_args()

    try:
        explore_schema(args.file, args.output, args.search, args.format)
    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
