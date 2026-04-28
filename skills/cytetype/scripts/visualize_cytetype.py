#!/usr/bin/env python3
"""
CyteType Visualization Script

Generates UMAP and dotplot visualizations for CyteType-annotated data.

Dotplot specifications:
- X-axis: CyteType annotation
- Y-axis: Original cluster annotation (group_key)
- Dot color: Confidence score (bwr colormap, centered at 0)
- Dot size: Fraction of cells shared between cluster and annotation

Usage:
    python visualize_cytetype.py \
        --input annotated.h5ad \
        --group-key cluster \
        --output-dir ./results
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scanpy as sc


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate CyteType visualizations (UMAP and dotplot)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python visualize_cytetype.py --input annotated.h5ad

    # Custom options
    python visualize_cytetype.py \\
        --input annotated.h5ad \\
        --group-key leiden \\
        --results-prefix cytetype \\
        --output-dir ./figures
        """
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to CyteType-annotated h5ad file"
    )
    parser.add_argument(
        "--group-key", "-g",
        default="cluster",
        help="Original cluster column in adata.obs (default: cluster)"
    )
    parser.add_argument(
        "--results-prefix",
        default="cytetype",
        help="CyteType results prefix (default: cytetype)"
    )
    parser.add_argument(
        "--output-dir", "-o",
        default=".",
        help="Output directory (default: current directory)"
    )
    parser.add_argument(
        "--figsize",
        default="12,10",
        help="Figure size as width,height (default: 12,10)"
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for output figures (default: 300)"
    )
    parser.add_argument(
        "--umap-only",
        action="store_true",
        help="Only generate UMAP plot"
    )
    parser.add_argument(
        "--dotplot-only",
        action="store_true",
        help="Only generate dotplot"
    )

    return parser.parse_args()


def get_confidence_scores(adata, results_prefix, group_key):
    """Extract confidence scores from CyteType results."""
    results_key = f"{results_prefix}_results"

    if results_key not in adata.uns:
        print(f"Warning: {results_key} not found in adata.uns")
        return None

    # Parse results JSON
    results = adata.uns[results_key]
    if isinstance(results, dict) and "result" in results:
        result_data = json.loads(results["result"])
    else:
        result_data = results

    # Extract confidence per cluster
    confidence_map = {}
    if "annotations" in result_data:
        for ann in result_data["annotations"]:
            cluster_id = ann.get("clusterId")
            confidence = ann.get("confidence", 0.5)
            if cluster_id is not None:
                confidence_map[str(cluster_id)] = confidence

    return confidence_map


def create_annotation_comparison_df(adata, group_key, results_prefix, confidence_map):
    """Create DataFrame comparing original clusters vs CyteType annotations."""
    annotation_col = f"{results_prefix}_annotation_{group_key}"

    if annotation_col not in adata.obs.columns:
        raise ValueError(f"Annotation column '{annotation_col}' not found in adata.obs")

    # Create cross-tabulation
    crosstab = pd.crosstab(
        adata.obs[group_key],
        adata.obs[annotation_col],
        normalize='index'  # Normalize by row (cluster)
    )

    # Calculate fraction shared for each cluster-annotation pair
    data_rows = []
    for cluster in crosstab.index:
        for annotation in crosstab.columns:
            fraction = crosstab.loc[cluster, annotation]
            if fraction > 0:  # Only include non-zero entries
                # Get confidence for this annotation
                # Map cluster to CyteType cluster ID (1-indexed)
                cluster_id = str(list(adata.obs[group_key].cat.categories).index(cluster) + 1)
                confidence = confidence_map.get(cluster_id, 0.5) if confidence_map else 0.5

                data_rows.append({
                    'cluster': str(cluster),
                    'annotation': annotation,
                    'fraction': fraction,
                    'confidence': confidence,
                    # Confidence centered at 0.5, map to -0.5 to 0.5 for bwr colormap
                    'confidence_centered': confidence - 0.5,
                })

    return pd.DataFrame(data_rows)


def plot_umap(adata, group_key, results_prefix, output_path, figsize, dpi):
    """Generate UMAP plot colored by CyteType annotation."""
    annotation_col = f"{results_prefix}_annotation_{group_key}"

    if annotation_col not in adata.obs.columns:
        print(f"Error: {annotation_col} not found in adata.obs")
        return False

    # Check for UMAP coordinates
    if "X_umap" not in adata.obsm:
        print("Warning: X_umap not found. Attempting to compute UMAP...")
        try:
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
        except Exception as e:
            print(f"Error computing UMAP: {e}")
            return False

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Original clusters
    sc.pl.umap(
        adata,
        color=group_key,
        ax=axes[0],
        show=False,
        title=f"Original Clusters ({group_key})",
        legend_loc="right margin" if adata.obs[group_key].nunique() <= 20 else "on data",
    )

    # CyteType annotations
    sc.pl.umap(
        adata,
        color=annotation_col,
        ax=axes[1],
        show=False,
        title="CyteType Annotations",
        legend_loc="right margin" if adata.obs[annotation_col].nunique() <= 20 else "on data",
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()

    print(f"Saved UMAP plot to: {output_path}")
    return True


def plot_dotplot(df, output_path, figsize, dpi):
    """
    Generate custom dotplot showing annotation comparison.

    X-axis: CyteType annotation
    Y-axis: Original cluster
    Dot color: Confidence (bwr colormap, centered at 0)
    Dot size: Fraction shared
    """
    if df.empty:
        print("Error: No data for dotplot")
        return False

    # Get unique values and sort
    clusters = sorted(df['cluster'].unique(), key=lambda x: int(x) if x.isdigit() else x)
    annotations = sorted(df['annotation'].unique())

    # Create figure
    fig_width, fig_height = figsize
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Create position mappings
    annotation_pos = {ann: i for i, ann in enumerate(annotations)}
    cluster_pos = {clust: i for i, clust in enumerate(clusters)}

    # Normalize sizes for plotting
    # Size should be proportional to fraction, with reasonable min/max
    min_size = 20
    max_size = 500
    size_scale = (max_size - min_size)

    # Create colormap: bwr centered at 0
    # Confidence ranges 0-1, centered at 0.5
    # Map to bwr: 0->blue, 0.5->white, 1->red
    norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=1)
    cmap = plt.cm.bwr

    # Plot dots
    for _, row in df.iterrows():
        x = annotation_pos[row['annotation']]
        y = cluster_pos[row['cluster']]
        size = min_size + row['fraction'] * size_scale
        color = cmap(norm(row['confidence']))

        ax.scatter(x, y, s=size, c=[color], edgecolors='black', linewidths=0.5)

    # Set axis labels and ticks
    ax.set_xticks(range(len(annotations)))
    ax.set_xticklabels(annotations, rotation=45, ha='right', fontsize=10)
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels(clusters, fontsize=10)

    ax.set_xlabel("CyteType Annotation", fontsize=12)
    ax.set_ylabel(f"Original Cluster", fontsize=12)
    ax.set_title("CyteType Annotation vs Original Clusters", fontsize=14)

    # Add grid
    ax.set_axisbelow(True)
    ax.grid(True, linestyle='--', alpha=0.3)

    # Adjust limits
    ax.set_xlim(-0.5, len(annotations) - 0.5)
    ax.set_ylim(-0.5, len(clusters) - 0.5)

    # Create legend for size
    size_legend_values = [0.25, 0.5, 0.75, 1.0]
    size_legend_handles = []
    for val in size_legend_values:
        size = min_size + val * size_scale
        handle = ax.scatter([], [], s=size, c='gray', edgecolors='black', linewidths=0.5)
        size_legend_handles.append(handle)

    legend1 = ax.legend(
        size_legend_handles,
        [f"{v:.0%}" for v in size_legend_values],
        title="Fraction Shared",
        loc='upper left',
        bbox_to_anchor=(1.02, 1),
        frameon=True,
    )
    ax.add_artist(legend1)

    # Create colorbar for confidence
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, aspect=20, pad=0.15)
    cbar.set_label("Confidence", fontsize=10)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar.set_ticklabels(['0.0', '0.25', '0.5', '0.75', '1.0'])

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()

    print(f"Saved dotplot to: {output_path}")
    return True


def main():
    """Main function to generate visualizations."""
    args = parse_args()

    # Parse figsize
    try:
        figsize = tuple(map(float, args.figsize.split(',')))
    except ValueError:
        print(f"Error: Invalid figsize format '{args.figsize}'. Use 'width,height'")
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate timestamp
    timestamp = datetime.now().strftime("%Y%m%d")

    print(f"\nLoading annotated h5ad: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Check for CyteType annotations
    annotation_col = f"{args.results_prefix}_annotation_{args.group_key}"
    if annotation_col not in adata.obs.columns:
        print(f"\nError: CyteType annotation column '{annotation_col}' not found.")
        print(f"Available columns: {list(adata.obs.columns)}")
        print("\nMake sure the h5ad file was annotated with CyteType using:")
        print(f"  --group-key {args.group_key}")
        print(f"  --results-prefix {args.results_prefix}")
        sys.exit(1)

    print(f"\nFound CyteType annotations in: {annotation_col}")
    print(f"Annotation counts:")
    for ct, count in adata.obs[annotation_col].value_counts().items():
        print(f"  - {ct}: {count:,}")

    # Get confidence scores
    confidence_map = get_confidence_scores(adata, args.results_prefix, args.group_key)
    if confidence_map:
        print(f"\nExtracted confidence scores for {len(confidence_map)} clusters")
    else:
        print("\nWarning: Could not extract confidence scores, using default 0.5")

    # Generate UMAP
    if not args.dotplot_only:
        umap_path = output_dir / f"{args.results_prefix}_umap_{timestamp}.png"
        plot_umap(adata, args.group_key, args.results_prefix, umap_path, figsize, args.dpi)

    # Generate dotplot
    if not args.umap_only:
        try:
            comparison_df = create_annotation_comparison_df(
                adata, args.group_key, args.results_prefix, confidence_map
            )
            dotplot_path = output_dir / f"{args.results_prefix}_dotplot_{timestamp}.png"
            plot_dotplot(comparison_df, dotplot_path, figsize, args.dpi)

            # Also save the comparison data
            csv_path = output_dir / f"{args.results_prefix}_comparison_{timestamp}.csv"
            comparison_df.to_csv(csv_path, index=False)
            print(f"Saved comparison data to: {csv_path}")

        except Exception as e:
            print(f"Error generating dotplot: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 60)
    print("Visualization complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
