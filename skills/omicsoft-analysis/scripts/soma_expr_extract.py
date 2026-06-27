#!/usr/bin/env python3
"""
SOMA Expression/DEG Extraction Script

Efficient Python SOMA extraction that outputs TSV to stdout for R pipe consumption.
Uses tiledbsoma AxisQuery for server-side filtering.

Usage:
    python soma_expr_extract.py --uri <S3_URI> --mode expr --genes GENE1,GENE2 --per-sample-studies STUDY1,STUDY2
    python soma_expr_extract.py --uri <S3_URI> --mode deg --genes GENE1,GENE2 --per-sample-studies STUDY1,STUDY2

Output: TSV to stdout (piped to R via read.delim(pipe(...)))

Dependencies:
    - tiledbsoma
    - pandas
    - numpy
"""

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

try:
    import tiledbsoma
    SOMA_AVAILABLE = True
except ImportError:
    SOMA_AVAILABLE = False


def get_soma_context(region=None):
    """Build SOMA context for S3 access."""
    if not SOMA_AVAILABLE:
        raise ImportError("tiledbsoma is not installed")

    region = region or os.environ.get("AWS_DEFAULT_REGION", "us-east-1")
    config = {"vfs.s3.region": region}

    if hasattr(tiledbsoma, 'SOMAContext'):
        return tiledbsoma.SOMAContext(tiledb_config=config)
    else:
        context = tiledbsoma.SOMATileDBContext()
        return context.replace(tiledb_config=config)


def build_obs_filter(studies=None, sources=None, diseases=None, tissues=None):
    """Build SOMA obs value_filter string from filter params."""
    clauses = []

    if studies:
        escaped = [f"'{s.replace(chr(39), chr(39)+chr(39))}'" for s in studies]
        clauses.append(f"project_id in [{', '.join(escaped)}]")

    if sources:
        escaped = [f"'{s}'" for s in sources]
        clauses.append(f"source in [{', '.join(escaped)}]")

    if diseases:
        escaped = [f"'{d.replace(chr(39), chr(39)+chr(39))}'" for d in diseases]
        clauses.append(f"disease in [{', '.join(escaped)}]")

    if tissues:
        escaped = [f"'{t.replace(chr(39), chr(39)+chr(39))}'" for t in tissues]
        clauses.append(f"tissue in [{', '.join(escaped)}]")

    return " and ".join(clauses) if clauses else None


def extract_expr(uri, genes, studies=None, sources=None, diseases=None,
                 tissues=None, context=None):
    """
    Extract sample-level expression data from EXPR SOMA.

    Returns DataFrame with metadata columns + gene expression columns.
    """
    print(f"[soma_expr_extract] Opening EXPR SOMA: {uri}", file=sys.stderr)

    with tiledbsoma.open(uri, mode="r", context=context) as exp:
        ms = exp.ms["RNA"]

        # Get var (gene) info to filter by gene symbols
        var_df = ms.var.read().concat().to_pandas()
        if 'gene_symbol' in var_df.columns:
            gene_col = 'gene_symbol'
        elif 'feature_name' in var_df.columns:
            gene_col = 'feature_name'
        elif 'gene_id' in var_df.columns:
            gene_col = 'gene_id'
        elif 'gene_name' in var_df.columns:
            gene_col = 'gene_name'
        else:
            gene_col = var_df.columns[0]

        # Filter to requested genes (or extract all if genes == ['ALL'])
        if genes == ['ALL']:
            gene_mask = pd.Series([True] * len(var_df), index=var_df.index)
            matched_genes = var_df[gene_col].tolist()
            print(f"[soma_expr_extract] Extracting ALL {len(matched_genes)} genes",
                  file=sys.stderr)
        else:
            gene_mask = var_df[gene_col].isin(genes)
            matched_genes = var_df[gene_mask][gene_col].tolist()
            missing_genes = set(genes) - set(matched_genes)

            if missing_genes:
                print(f"[soma_expr_extract] WARNING: genes not found in var: "
                      f"{sorted(missing_genes)}", file=sys.stderr)

            if not matched_genes:
                print("[soma_expr_extract] ERROR: No genes found in SOMA var",
                      file=sys.stderr)
                sys.exit(1)

        # Get soma_joinid for matched genes
        var_ids = var_df[gene_mask]["soma_joinid"].tolist()

        n_total = len(genes) if genes != ['ALL'] else 'ALL'
        print(f"[soma_expr_extract] Matched {len(matched_genes)}/{n_total} genes",
              file=sys.stderr)

        # Build obs filter
        obs_filter = build_obs_filter(studies, sources, diseases, tissues)
        obs_query = None
        if obs_filter:
            obs_query = tiledbsoma.AxisQuery(value_filter=obs_filter)
            print(f"[soma_expr_extract] obs filter: {obs_filter[:200]}",
                  file=sys.stderr)

        var_query = tiledbsoma.AxisQuery(coords=(var_ids,))

        # Execute query
        query = exp.axis_query(
            measurement_name="RNA",
            obs_query=obs_query,
            var_query=var_query
        )

        # Get obs metadata
        obs_df = query.obs().concat().to_pandas()
        print(f"[soma_expr_extract] Retrieved {len(obs_df)} samples",
              file=sys.stderr)

        if len(obs_df) == 0:
            print("[soma_expr_extract] ERROR: No samples after filtering",
                  file=sys.stderr)
            sys.exit(1)

        # Get expression matrix - use first available layer
        x_layers = list(ms.X.keys())
        # Prefer normalized data
        layer_priority = ['data', 'normalized', 'tpm', 'fpkm', 'raw_counts']
        x_layer = None
        for lp in layer_priority:
            if lp in x_layers:
                x_layer = lp
                break
        if x_layer is None:
            x_layer = x_layers[0]

        print(f"[soma_expr_extract] Using X layer: '{x_layer}'", file=sys.stderr)

        # Convert to AnnData for easy matrix access
        adata = query.to_anndata(X_name=x_layer)
        query.close()

        # Build output DataFrame: obs metadata + gene expression columns
        result_df = obs_df.copy()
        result_df.index = range(len(result_df))

        # Add gene expression values as columns (vectorized, no per-column insertion)
        from scipy import sparse
        X = adata.X
        if sparse.issparse(X):
            X = X.toarray()

        query_var_df = var_df[gene_mask].reset_index(drop=True)
        gene_names = query_var_df[gene_col].values[:X.shape[1]]
        expr_df = pd.DataFrame(X[:, :len(gene_names)], columns=gene_names, index=result_df.index)
        result_df = pd.concat([result_df, expr_df], axis=1)

        # Drop soma_joinid if present (internal ID)
        if 'soma_joinid' in result_df.columns:
            result_df = result_df.drop(columns=['soma_joinid'])

        print(f"[soma_expr_extract] Output: {result_df.shape[0]} rows x "
              f"{result_df.shape[1]} cols", file=sys.stderr)

        return result_df


def extract_deg(uri, genes, studies=None, sources=None, diseases=None,
                tissues=None, context=None):
    """
    Extract DEG statistics from DEG SOMA.

    Returns DataFrame with comparison metadata + per-gene log2fc/padj columns.
    """
    print(f"[soma_expr_extract] Opening DEG SOMA: {uri}", file=sys.stderr)

    with tiledbsoma.open(uri, mode="r", context=context) as exp:
        ms = exp.ms["RNA"]

        # Get var info
        var_df = ms.var.read().concat().to_pandas()
        if 'gene_symbol' in var_df.columns:
            gene_col = 'gene_symbol'
        elif 'feature_name' in var_df.columns:
            gene_col = 'feature_name'
        elif 'gene_id' in var_df.columns:
            gene_col = 'gene_id'
        elif 'gene_name' in var_df.columns:
            gene_col = 'gene_name'
        else:
            gene_col = var_df.columns[0]

        gene_mask = var_df[gene_col].isin(genes)
        matched_genes = var_df[gene_mask][gene_col].tolist()
        missing_genes = set(genes) - set(matched_genes)

        if missing_genes:
            print(f"[soma_expr_extract] WARNING: genes not in DEG var: "
                  f"{sorted(missing_genes)}", file=sys.stderr)

        if not matched_genes:
            print("[soma_expr_extract] ERROR: No genes found in DEG SOMA var",
                  file=sys.stderr)
            sys.exit(1)

        var_ids = var_df[gene_mask]["soma_joinid"].tolist()

        print(f"[soma_expr_extract] DEG: matched {len(matched_genes)}/{len(genes)} genes",
              file=sys.stderr)

        # Build obs filter
        obs_filter = build_obs_filter(studies, sources, diseases, tissues)
        obs_query = None
        if obs_filter:
            obs_query = tiledbsoma.AxisQuery(value_filter=obs_filter)

        var_query = tiledbsoma.AxisQuery(coords=(var_ids,))

        query = exp.axis_query(
            measurement_name="RNA",
            obs_query=obs_query,
            var_query=var_query
        )

        obs_df = query.obs().concat().to_pandas()
        print(f"[soma_expr_extract] DEG: {len(obs_df)} comparisons",
              file=sys.stderr)

        if len(obs_df) == 0:
            print("[soma_expr_extract] WARNING: No DEG comparisons after filtering",
                  file=sys.stderr)
            query.close()
            return pd.DataFrame()

        # Get available layers
        x_layers = list(ms.X.keys())
        deg_layers = ['log2fc', 'pval', 'padj', 'sig_score']
        available_layers = [l for l in deg_layers if l in x_layers]

        if not available_layers:
            print(f"[soma_expr_extract] WARNING: No DEG layers found. "
                  f"Available: {x_layers}", file=sys.stderr)
            query.close()
            return pd.DataFrame()

        # Build long-format result: one row per (gene × comparison)
        from scipy import sparse
        query_var_df = var_df[gene_mask].reset_index(drop=True)

        X_arrays = {}
        for layer_name in available_layers:
            adata = query.to_anndata(X_name=layer_name)
            X = adata.X
            if sparse.issparse(X):
                X = X.toarray()
            X_arrays[layer_name] = X

        rows = []
        for i in range(obs_df.shape[0]):
            obs_row = obs_df.iloc[i]
            for j, gene_name in enumerate(query_var_df[gene_col].values):
                row_data = {
                    "comparison_id": obs_row.get("comparison_id", ""),
                    "comparison_contrast": obs_row.get("comparison_contrast", ""),
                    "project_id": obs_row.get("project_id", ""),
                    "tissue": obs_row.get("tissue", ""),
                    "disease": obs_row.get("disease", obs_row.get("disease_state", "")),
                    "gene": gene_name,
                }
                for layer_name in available_layers:
                    if j < X_arrays[layer_name].shape[1]:
                        row_data[layer_name] = X_arrays[layer_name][i, j]
                rows.append(row_data)

        query.close()

        result_df = pd.DataFrame(rows)

        print(f"[soma_expr_extract] DEG output (long): {result_df.shape[0]} rows x "
              f"{result_df.shape[1]} cols", file=sys.stderr)

        return result_df


def main():
    parser = argparse.ArgumentParser(
        description='Extract expression/DEG data from SOMA to TSV stdout')
    parser.add_argument('--uri', required=True,
                        help='SOMA experiment URI (S3 or local)')
    parser.add_argument('--mode', required=True, choices=['expr', 'deg'],
                        help='Extraction mode: expr (expression) or deg (statistics)')
    parser.add_argument('--genes', required=True,
                        help='Comma-separated gene symbols')
    parser.add_argument('--per-sample-studies', default=None,
                        help='Comma-separated project_id values for per-sample analysis')
    parser.add_argument('--sources', default=None,
                        help='Comma-separated source values (internal,curated,omicsoft)')
    parser.add_argument('--diseases', default=None,
                        help='Comma-separated disease exact values')
    parser.add_argument('--tissues', default=None,
                        help='Comma-separated tissue exact values')
    parser.add_argument('--region', default=None,
                        help='AWS region (default: AWS_DEFAULT_REGION or us-east-1)')

    args = parser.parse_args()

    if not SOMA_AVAILABLE:
        print("ERROR: tiledbsoma is not installed", file=sys.stderr)
        sys.exit(1)

    # Parse list arguments
    genes = [g.strip() for g in args.genes.split(',') if g.strip()]
    studies = [s.strip() for s in args.per_sample_studies.split(',') if s.strip()] if args.per_sample_studies else None
    sources = [s.strip() for s in args.sources.split(',') if s.strip()] if args.sources else None
    diseases = [d.strip() for d in args.diseases.split(',') if d.strip()] if args.diseases else None
    tissues = [t.strip() for t in args.tissues.split(',') if t.strip()] if args.tissues else None

    # Get SOMA context
    context = get_soma_context(args.region)

    # Extract data
    if args.mode == 'expr':
        result_df = extract_expr(
            args.uri, genes, studies=studies, sources=sources,
            diseases=diseases, tissues=tissues, context=context
        )
    else:
        result_df = extract_deg(
            args.uri, genes, studies=studies, sources=sources,
            diseases=diseases, tissues=tissues, context=context
        )

    # Output TSV to stdout
    if result_df is not None and len(result_df) > 0:
        result_df.to_csv(sys.stdout, sep='\t', index=False)
    else:
        print("ERROR: No data extracted", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
