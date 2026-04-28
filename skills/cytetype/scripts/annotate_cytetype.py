#!/usr/bin/env python3
"""
CyteType Cell Type Annotation Script

Performs automated cell type annotation on single-cell RNA-seq data using CyteType.
Loads LLM configuration from environment variables via dotenv.

Usage:
    python annotate_cytetype.py \
        --input data.h5ad \
        --group-key cluster \
        --study-context "Human PBMC from healthy donor" \
        --output-dir ./results
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import scanpy as sc
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="CyteType cell type annotation for single-cell RNA-seq data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python annotate_cytetype.py --input pbmc.h5ad --study-context "Human PBMC"

    # Custom parameters
    python annotate_cytetype.py \\
        --input data.h5ad \\
        --group-key leiden \\
        --study-context "Mouse liver from NASH model" \\
        --n-top-genes 100 \\
        --output-dir ./cytetype_results
        """
    )

    # Required arguments
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to input h5ad file"
    )
    parser.add_argument(
        "--study-context", "-s",
        required=True,
        help="Biological context (organism, tissue, disease, etc.)"
    )

    # CyteType initialization parameters
    parser.add_argument(
        "--group-key", "-g",
        default="cluster",
        help="Column in adata.obs containing cluster labels (default: cluster)"
    )
    parser.add_argument(
        "--rank-key",
        default="rank_genes_groups",
        help="Key in adata.uns with DE results (default: rank_genes_groups)"
    )
    parser.add_argument(
        "--gene-symbols-column",
        default="gene_symbols",
        help="Column in adata.var with gene symbols (default: gene_symbols)"
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=50,
        help="Number of top marker genes per cluster (default: 50)"
    )
    parser.add_argument(
        "--aggregate-metadata",
        action="store_true",
        default=True,
        help="Aggregate metadata from AnnData (default: True)"
    )
    parser.add_argument(
        "--no-aggregate-metadata",
        action="store_false",
        dest="aggregate_metadata",
        help="Disable metadata aggregation"
    )
    parser.add_argument(
        "--min-percentage",
        type=int,
        default=10,
        help="Min percentage for cluster context (default: 10)"
    )
    parser.add_argument(
        "--coordinates-key",
        default="X_umap",
        help="Coordinates key for visualization (default: X_umap)"
    )
    parser.add_argument(
        "--max-cells-per-group",
        type=int,
        default=1000,
        help="Max cells per group for visualization (default: 1000)"
    )

    # CyteType run parameters
    parser.add_argument(
        "--n-parallel",
        type=int,
        default=2,
        help="Number of parallel cluster requests (default: 2, max: 50)"
    )
    parser.add_argument(
        "--results-prefix",
        default="cytetype",
        help="Prefix for result columns (default: cytetype)"
    )
    parser.add_argument(
        "--poll-interval",
        type=int,
        default=10,
        help="Polling interval in seconds (default: 10)"
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=7200,
        help="Timeout in seconds (default: 7200 = 2 hours)"
    )

    # LLM configuration
    parser.add_argument(
        "--provider",
        default="openai",
        choices=["anthropic", "bedrock", "google", "groq", "mistral", "openai", "openrouter"],
        help="LLM provider (default: openai)"
    )
    parser.add_argument(
        "--model",
        default="gpt-5.2",
        help="LLM model name (default: gpt-5.2)"
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="API key (default: from OPENAI_API_KEY env var)"
    )
    parser.add_argument(
        "--base-url",
        default=None,
        help="API base URL (default: from OPENAI_API_BASE env var)"
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=0.0,
        help="LLM temperature (default: 0.0)"
    )
    parser.add_argument(
        "--max-tokens",
        type=int,
        default=4096,
        help="LLM max tokens (default: 4096)"
    )

    # Output options
    parser.add_argument(
        "--output-dir", "-o",
        default=".",
        help="Output directory (default: current directory)"
    )
    parser.add_argument(
        "--save-query",
        action="store_true",
        default=True,
        help="Save query JSON (default: True)"
    )
    parser.add_argument(
        "--no-save-query",
        action="store_false",
        dest="save_query",
        help="Don't save query JSON"
    )

    return parser.parse_args()


def get_llm_config(args):
    """Build LLM configuration from arguments and environment variables."""
    # Get API key from args or environment
    api_key = args.api_key or os.getenv("OPENAI_API_KEY")
    if not api_key:
        raise ValueError(
            "API key not found. Set OPENAI_API_KEY in .env file or pass --api-key"
        )

    # Get base URL from args or environment
    base_url = args.base_url or os.getenv("OPENAI_API_BASE")

    config = {
        "provider": args.provider,
        "name": args.model,
        "apiKey": api_key,
        "modelSettings": {
            "temperature": args.temperature,
            "max_tokens": args.max_tokens,
        }
    }

    if base_url:
        config["baseUrl"] = base_url

    return [config]


def validate_adata(adata, args):
    """Validate AnnData object has required components."""
    errors = []
    warnings = []

    # Check group_key exists
    if args.group_key not in adata.obs.columns:
        errors.append(f"group_key '{args.group_key}' not found in adata.obs")
        # Suggest alternatives
        categorical_cols = [c for c in adata.obs.columns if adata.obs[c].dtype.name == 'category']
        if categorical_cols:
            warnings.append(f"Available categorical columns: {categorical_cols[:10]}")

    # Check rank_key exists
    if args.rank_key not in adata.uns:
        errors.append(f"rank_key '{args.rank_key}' not found in adata.uns")
        warnings.append(
            f"Run: sc.tl.rank_genes_groups(adata, groupby='{args.group_key}') first"
        )

    # Check gene symbols column
    if args.gene_symbols_column not in adata.var.columns:
        # Try using var_names as fallback
        if adata.var_names.str.match(r'^[A-Z][A-Z0-9]+$').any():
            warnings.append(
                f"gene_symbols_column '{args.gene_symbols_column}' not found, "
                "but var_names appear to be gene symbols"
            )
        else:
            errors.append(
                f"gene_symbols_column '{args.gene_symbols_column}' not found in adata.var"
            )

    # Check coordinates
    if args.coordinates_key not in adata.obsm:
        warnings.append(
            f"coordinates_key '{args.coordinates_key}' not found in adata.obsm. "
            "Visualization may be limited."
        )

    return errors, warnings


def print_summary(adata, args, llm_config):
    """Print configuration summary for user confirmation."""
    print("\n" + "=" * 60)
    print("CyteType Configuration Summary")
    print("=" * 60)

    print(f"\nInput File: {args.input}")
    print(f"Output Directory: {args.output_dir}")

    print(f"\nStudy Context: {args.study_context}")

    print("\nData Summary:")
    print(f"  - Cells: {adata.n_obs:,}")
    print(f"  - Genes: {adata.n_vars:,}")
    if args.group_key in adata.obs.columns:
        n_clusters = adata.obs[args.group_key].nunique()
        print(f"  - Clusters ({args.group_key}): {n_clusters}")

    print("\nInitialization Parameters:")
    print(f"  - group_key: {args.group_key}")
    print(f"  - rank_key: {args.rank_key}")
    print(f"  - gene_symbols_column: {args.gene_symbols_column}")
    print(f"  - n_top_genes: {args.n_top_genes}")
    print(f"  - aggregate_metadata: {args.aggregate_metadata}")
    print(f"  - coordinates_key: {args.coordinates_key}")

    print("\nRun Parameters:")
    print(f"  - n_parallel_clusters: {args.n_parallel}")
    print(f"  - results_prefix: {args.results_prefix}")
    print(f"  - timeout_seconds: {args.timeout}")

    print("\nLLM Configuration:")
    print(f"  - provider: {llm_config[0]['provider']}")
    print(f"  - model: {llm_config[0]['name']}")
    if llm_config[0].get('baseUrl'):
        print(f"  - baseUrl: {llm_config[0]['baseUrl']}")
    print(f"  - temperature: {llm_config[0]['modelSettings']['temperature']}")

    print("\n" + "=" * 60)


def main():
    """Main function to run CyteType annotation."""
    args = parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate timestamp for output files
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    print(f"\nLoading h5ad file: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Validate AnnData
    errors, warnings = validate_adata(adata, args)

    if warnings:
        print("\nWarnings:")
        for w in warnings:
            print(f"  - {w}")

    if errors:
        print("\nErrors:")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)

    # Build LLM config
    try:
        llm_config = get_llm_config(args)
    except ValueError as e:
        print(f"\nError: {e}")
        sys.exit(1)

    # Print summary
    print_summary(adata, args, llm_config)

    # Import CyteType (after validation to fail fast if data issues)
    try:
        from cytetype import CyteType
    except ImportError:
        print("\nError: CyteType not installed. Run: pip install cytetype")
        sys.exit(1)

    print("\nInitializing CyteType...")
    annotator = CyteType(
        adata,
        group_key=args.group_key,
        rank_key=args.rank_key,
        gene_symbols_column=args.gene_symbols_column,
        n_top_genes=args.n_top_genes,
        aggregate_metadata=args.aggregate_metadata,
        min_percentage=args.min_percentage,
        coordinates_key=args.coordinates_key,
        max_cells_per_group=args.max_cells_per_group,
    )

    print("\nRunning CyteType annotation...")
    print(f"This may take several minutes depending on the number of clusters.")
    print(f"Timeout set to {args.timeout} seconds ({args.timeout/60:.1f} minutes)\n")

    # Prepare metadata
    metadata = {
        "input_file": str(args.input),
        "timestamp": timestamp,
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
    }

    # Query filename
    query_filename = output_dir / f"{args.results_prefix}_query_{timestamp}.json"

    # Run annotation
    adata = annotator.run(
        study_context=args.study_context,
        llm_configs=llm_config,
        metadata=metadata,
        n_parallel_clusters=args.n_parallel,
        results_prefix=args.results_prefix,
        poll_interval_seconds=args.poll_interval,
        timeout_seconds=args.timeout,
        save_query=args.save_query,
        query_filename=str(query_filename) if args.save_query else "query.json",
        show_progress=True,
    )

    # Save annotated h5ad
    output_h5ad = output_dir / f"{args.results_prefix}_annotated_{timestamp}.h5ad"
    print(f"\nSaving annotated h5ad to: {output_h5ad}")
    adata.write_h5ad(output_h5ad)

    # Save results JSON
    if f"{args.results_prefix}_results" in adata.uns:
        results_json = output_dir / f"{args.results_prefix}_results_{timestamp}.json"
        with open(results_json, "w") as f:
            json.dump(adata.uns[f"{args.results_prefix}_results"], f, indent=2)
        print(f"Saved results JSON to: {results_json}")

    # Print annotation summary
    annotation_col = f"{args.results_prefix}_annotation_{args.group_key}"
    if annotation_col in adata.obs.columns:
        print("\nAnnotation Summary:")
        print("-" * 40)
        for ct, count in adata.obs[annotation_col].value_counts().items():
            pct = count / adata.n_obs * 100
            print(f"  {ct}: {count:,} cells ({pct:.1f}%)")

    print("\n" + "=" * 60)
    print("CyteType annotation complete!")
    print("=" * 60)
    print(f"\nOutput files:")
    print(f"  - Annotated h5ad: {output_h5ad}")
    if args.save_query:
        print(f"  - Query JSON: {query_filename}")
    print(f"  - Results JSON: {results_json}")
    print(f"\nNext: Run visualize_cytetype.py to generate UMAP and dotplot")

    return adata


if __name__ == "__main__":
    main()
