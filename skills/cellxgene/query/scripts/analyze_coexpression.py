#!/usr/bin/env python3
"""
Gene Coexpression Analysis Script

Analyze pairwise gene coexpression from AnnData objects (local or from CELLxGENE queries).
Computes correlation matrices and generates heatmap visualizations.
"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Dict, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
from scipy.stats import pearsonr, spearmanr
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
import scanpy as sc

# Optional plotly import for interactive visualizations
try:
    import plotly.graph_objects as go
    import plotly.figure_factory as ff
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("Warning: plotly not installed. Interactive HTML heatmaps will not be generated.")


def load_gene_list(gene_list_arg: str) -> List[str]:
    """
    Load gene list from file or comma-separated string.

    Args:
        gene_list_arg: File path or comma-separated gene names

    Returns:
        List of gene symbols
    """
    # Check if it's a file
    if Path(gene_list_arg).is_file():
        with open(gene_list_arg, 'r') as f:
            genes = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        print(f"Loaded {len(genes)} genes from file: {gene_list_arg}")
        return genes
    else:
        # Treat as comma-separated list
        genes = [g.strip() for g in gene_list_arg.split(',') if g.strip()]
        print(f"Loaded {len(genes)} genes from argument")
        return genes


def filter_adata(
    adata,
    cell_type: Optional[str] = None,
    tissue: Optional[str] = None,
    disease: Optional[str] = None,
    metadata_filter: Optional[str] = None
):
    """
    Filter AnnData by metadata criteria.

    Args:
        adata: AnnData object
        cell_type: Filter by cell type (partial match)
        tissue: Filter by tissue (partial match)
        disease: Filter by disease (partial match)
        metadata_filter: Custom filter expression (e.g., "sex == 'female'")

    Returns:
        Filtered AnnData object
    """
    original_cells = adata.n_obs

    if cell_type:
        if 'cell_type' in adata.obs.columns:
            mask = adata.obs['cell_type'].astype(str).str.contains(cell_type, case=False, na=False)
            adata = adata[mask].copy()
            print(f"Filtered by cell_type containing '{cell_type}': {adata.n_obs} cells")
        else:
            print(f"Warning: 'cell_type' column not found in metadata")

    if tissue:
        tissue_col = None
        for col in ['tissue_general', 'tissue']:
            if col in adata.obs.columns:
                tissue_col = col
                break

        if tissue_col:
            mask = adata.obs[tissue_col].astype(str).str.contains(tissue, case=False, na=False)
            adata = adata[mask].copy()
            print(f"Filtered by {tissue_col} containing '{tissue}': {adata.n_obs} cells")
        else:
            print(f"Warning: tissue columns not found in metadata")

    if disease:
        if 'disease' in adata.obs.columns:
            mask = adata.obs['disease'].astype(str).str.contains(disease, case=False, na=False)
            adata = adata[mask].copy()
            print(f"Filtered by disease containing '{disease}': {adata.n_obs} cells")
        else:
            print(f"Warning: 'disease' column not found in metadata")

    if metadata_filter:
        try:
            adata = adata[adata.obs.eval(metadata_filter)].copy()
            print(f"Filtered by custom expression '{metadata_filter}': {adata.n_obs} cells")
        except Exception as e:
            print(f"Error applying metadata filter: {e}")
            sys.exit(1)

    print(f"Final cell count: {adata.n_obs} / {original_cells} ({100*adata.n_obs/original_cells:.1f}%)")

    return adata


def compute_correlation_matrix(
    adata,
    genes: List[str],
    method: str = 'pearson',
    use_raw: bool = False,
    min_cells: int = 10
) -> tuple:
    """
    Compute pairwise correlation matrix and p-values for gene list.

    Args:
        adata: AnnData object
        genes: List of gene symbols
        method: 'pearson' or 'spearman'
        use_raw: Use .raw.X instead of .X
        min_cells: Minimum cells expressing each gene

    Returns:
        Tuple of (correlation_df, pvalue_df)
    """
    print(f"\nComputing {method.capitalize()} correlation matrix...")

    # Get gene names from var
    if use_raw and adata.raw is not None:
        var_names = adata.raw.var_names
        X = adata.raw.X
    else:
        var_names = adata.var_names
        X = adata.X

    # Find genes in dataset
    genes_found = []
    genes_missing = []

    for gene in genes:
        if gene in var_names:
            genes_found.append(gene)
        else:
            genes_missing.append(gene)

    if genes_missing:
        print(f"Warning: {len(genes_missing)} genes not found in dataset: {genes_missing[:10]}" +
              (f" ... and {len(genes_missing)-10} more" if len(genes_missing) > 10 else ""))

    if len(genes_found) < 2:
        print(f"Error: Need at least 2 genes in dataset. Only found: {genes_found}")
        sys.exit(1)

    print(f"Found {len(genes_found)} genes in dataset")

    # Extract expression matrix for these genes
    gene_indices = [var_names.get_loc(g) for g in genes_found]

    if sparse.issparse(X):
        expr_matrix = X[:, gene_indices].toarray()
    else:
        expr_matrix = X[:, gene_indices]

    # Check expression levels
    n_cells_expressing = (expr_matrix > 0).sum(axis=0)
    low_expression_genes = [genes_found[i] for i, count in enumerate(n_cells_expressing) if count < min_cells]

    if low_expression_genes:
        print(f"Warning: {len(low_expression_genes)} genes expressed in fewer than {min_cells} cells: {low_expression_genes}")

    # Compute correlation matrix and p-values
    n_genes = len(genes_found)
    corr_matrix = np.zeros((n_genes, n_genes))
    pval_matrix = np.ones((n_genes, n_genes))  # Initialize with 1s

    for i in range(n_genes):
        for j in range(i, n_genes):
            if i == j:
                corr_matrix[i, j] = 1.0
                pval_matrix[i, j] = 0.0  # Diagonal p-value = 0
            else:
                x = expr_matrix[:, i]
                y = expr_matrix[:, j]

                # Compute correlation and p-value
                if method == 'pearson':
                    if np.std(x) > 0 and np.std(y) > 0:
                        corr, pval = pearsonr(x, y)
                    else:
                        corr = 0.0
                        pval = 1.0
                elif method == 'spearman':
                    if len(np.unique(x)) > 1 and len(np.unique(y)) > 1:
                        corr, pval = spearmanr(x, y)
                    else:
                        corr = 0.0
                        pval = 1.0
                else:
                    raise ValueError(f"Unknown method: {method}")

                corr_matrix[i, j] = corr
                corr_matrix[j, i] = corr
                pval_matrix[i, j] = pval
                pval_matrix[j, i] = pval

    # Create DataFrames
    corr_df = pd.DataFrame(corr_matrix, index=genes_found, columns=genes_found)
    pval_df = pd.DataFrame(pval_matrix, index=genes_found, columns=genes_found)

    # Print summary statistics
    print(f"\nCorrelation Matrix Summary:")
    print(f"  Dimensions: {corr_df.shape[0]} × {corr_df.shape[1]}")

    # Get upper triangle (excluding diagonal)
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)
    upper_triangle_corr = corr_matrix[mask]
    upper_triangle_pval = pval_matrix[mask]

    print(f"  Mean correlation: {upper_triangle_corr.mean():.3f}")
    print(f"  Median correlation: {np.median(upper_triangle_corr):.3f}")
    print(f"  Std correlation: {upper_triangle_corr.std():.3f}")
    print(f"  Min correlation: {upper_triangle_corr.min():.3f}")
    print(f"  Max correlation: {upper_triangle_corr.max():.3f}")

    # P-value statistics
    print(f"\nP-value Statistics:")
    print(f"  Mean p-value: {upper_triangle_pval.mean():.4f}")
    print(f"  Median p-value: {np.median(upper_triangle_pval):.4f}")
    n_sig_001 = (upper_triangle_pval < 0.001).sum()
    n_sig_01 = (upper_triangle_pval < 0.01).sum()
    n_sig_05 = (upper_triangle_pval < 0.05).sum()
    n_total = len(upper_triangle_pval)
    print(f"  Significant at p < 0.001: {n_sig_001} / {n_total} ({100*n_sig_001/n_total:.1f}%)")
    print(f"  Significant at p < 0.01:  {n_sig_01} / {n_total} ({100*n_sig_01/n_total:.1f}%)")
    print(f"  Significant at p < 0.05:  {n_sig_05} / {n_total} ({100*n_sig_05/n_total:.1f}%)")

    # Find strongest correlations with p-values
    print(f"\nTop 10 strongest correlations (with p-values):")
    corr_pairs = []
    for i in range(n_genes):
        for j in range(i+1, n_genes):
            corr_pairs.append((genes_found[i], genes_found[j], corr_matrix[i, j], pval_matrix[i, j]))

    corr_pairs.sort(key=lambda x: abs(x[2]), reverse=True)
    for gene1, gene2, corr, pval in corr_pairs[:10]:
        sig_marker = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
        print(f"  {gene1} <-> {gene2}: {corr:.3f} (p={pval:.2e}) {sig_marker}")

    return corr_df, pval_df


def plot_interactive_heatmap(
    corr_df: pd.DataFrame,
    pval_df: pd.DataFrame,
    output_path: Path,
    method: str = 'pearson',
    title: Optional[str] = None,
    cluster: bool = True
):
    """
    Create interactive HTML heatmap with hierarchical clustering and significance markers.

    Args:
        corr_df: Correlation matrix DataFrame
        pval_df: P-value matrix DataFrame
        output_path: Output HTML file path
        method: Correlation method (for title)
        title: Custom plot title
        cluster: Perform hierarchical clustering
    """
    if not PLOTLY_AVAILABLE:
        print("Plotly not available. Skipping interactive heatmap generation.")
        return

    print(f"\nCreating interactive HTML heatmap...")

    genes = corr_df.index.tolist()
    n_genes = len(genes)

    # Create text annotations with correlation values and p-value significance
    text_annot = []
    for i in range(n_genes):
        row_annot = []
        for j in range(n_genes):
            r = corr_df.iloc[i, j]
            p = pval_df.iloc[i, j]

            # Add significance stars
            sig = ""
            if i != j:  # Skip diagonal
                if p < 0.001:
                    sig = "***"
                elif p < 0.01:
                    sig = "**"
                elif p < 0.05:
                    sig = "*"

            # Format annotation
            if i == j:
                row_annot.append(f"{r:.2f}")
            else:
                row_annot.append(f"{r:.3f}{sig}")
        text_annot.append(row_annot)

    # Perform hierarchical clustering if requested and more than 2 genes
    if cluster and n_genes > 2:
        # Convert correlation to distance (1 - abs(correlation))
        dist_matrix = 1 - np.abs(corr_df.values)
        np.fill_diagonal(dist_matrix, 0)

        # Convert to condensed distance matrix
        dist_condensed = squareform(dist_matrix, checks=False)

        # Perform hierarchical clustering
        linkage = hierarchy.linkage(dist_condensed, method='average')
        dendro = hierarchy.dendrogram(linkage, no_plot=True)
        gene_order = dendro['leaves']

        # Reorder matrices
        corr_ordered = corr_df.iloc[gene_order, gene_order]
        text_ordered = [[text_annot[i][j] for j in gene_order] for i in gene_order]
        genes_ordered = [genes[i] for i in gene_order]

        print(f"Applied hierarchical clustering (average linkage)")
    else:
        corr_ordered = corr_df
        text_ordered = text_annot
        genes_ordered = genes

    # Create interactive heatmap using plotly
    fig = go.Figure(data=go.Heatmap(
        z=corr_ordered.values,
        x=genes_ordered,
        y=genes_ordered,
        text=text_ordered,
        texttemplate='%{text}',
        textfont={"size": 10 if n_genes <= 20 else 8},
        colorscale='RdBu_r',  # Red-Blue reversed (red=positive, blue=negative)
        zmid=0,  # Center colorscale at 0
        zmin=-1,
        zmax=1,
        colorbar=dict(
            title=f"{method.capitalize()}<br>Correlation",
            tickmode="linear",
            tick0=-1,
            dtick=0.5
        ),
        hovertemplate='<b>%{y} vs %{x}</b><br>Correlation: %{z:.3f}<extra></extra>'
    ))

    # Update layout
    if title is None:
        title = f"Gene Coexpression ({method.capitalize()} correlation)"

    fig.update_layout(
        title={
            'text': title,
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 16, 'family': 'Arial Black'}
        },
        xaxis={'title': '', 'side': 'bottom', 'tickangle': -45},
        yaxis={'title': '', 'autorange': 'reversed'},
        width=max(600, n_genes * 40),
        height=max(600, n_genes * 40),
        plot_bgcolor='white',
        paper_bgcolor='white'
    )

    # Add annotations for significance levels
    fig.add_annotation(
        text="* p<0.05, ** p<0.01, *** p<0.001",
        xref="paper", yref="paper",
        x=0.5, y=-0.10,
        showarrow=False,
        font=dict(size=10, color="gray"),
        xanchor='center'
    )

    # Save as HTML
    html_path = Path(str(output_path).replace('.png', '.html'))
    fig.write_html(html_path)
    print(f"Interactive heatmap saved to: {html_path}")


def plot_correlation_heatmap(
    corr_df: pd.DataFrame,
    pval_df: pd.DataFrame,
    output_path: Path,
    method: str = 'pearson',
    title: Optional[str] = None,
    figsize: tuple = (10, 8),
    cmap: str = 'bwr',
    vmin: float = -1.0,
    vmax: float = 1.0,
    cluster: bool = True,
    annot: bool = True,
    fmt: str = '.2f',
    generate_html: bool = True
):
    """
    Create correlation heatmap with hierarchical clustering and significance markers.
    Generates both static PNG and interactive HTML versions.

    Args:
        corr_df: Correlation matrix DataFrame
        pval_df: P-value matrix DataFrame
        output_path: Output file path (PNG)
        method: Correlation method (for title)
        title: Custom plot title
        figsize: Figure size (width, height)
        cmap: Colormap (default: bwr - blue-white-red)
        vmin: Minimum value for colormap
        vmax: Maximum value for colormap
        cluster: Perform hierarchical clustering
        annot: Annotate cells with correlation values and significance
        fmt: Format string for annotations
        generate_html: Also generate interactive HTML heatmap (default: True)
    """
    print(f"\nCreating correlation heatmap...")

    # Create significance markers based on p-values
    def get_sig_marker(pval):
        """Get significance marker based on p-value."""
        if pval < 0.0001:
            return '****'
        elif pval < 0.001:
            return '***'
        elif pval < 0.01:
            return '**'
        elif pval < 0.05:
            return '*'
        else:
            return ''

    # Create custom annotation matrix with correlation values and significance stars
    annot_matrix = []
    for i in range(len(corr_df)):
        row = []
        for j in range(len(corr_df)):
            corr_val = corr_df.iloc[i, j]
            pval = pval_df.iloc[i, j]
            sig = get_sig_marker(pval)
            # Format: correlation value with significance stars
            if i == j:
                text = f"{corr_val:.2f}"  # Diagonal without stars
            else:
                text = f"{corr_val:.2f}{sig}"
            row.append(text)
        annot_matrix.append(row)

    annot_matrix = np.array(annot_matrix)

    # Set up plot
    fig, ax = plt.subplots(figsize=figsize)

    # Create heatmap
    if cluster and corr_df.shape[0] > 2:
        # Perform hierarchical clustering
        linkage = hierarchy.linkage(corr_df.values, method='average')
        dendro = hierarchy.dendrogram(linkage, no_plot=True)
        order = dendro['leaves']

        corr_df_ordered = corr_df.iloc[order, order]
        annot_matrix_ordered = annot_matrix[order, :][:, order]

        print(f"Applied hierarchical clustering (average linkage)")
    else:
        corr_df_ordered = corr_df
        annot_matrix_ordered = annot_matrix

    # Determine annotation based on matrix size
    if annot and corr_df.shape[0] > 20:
        annot = False
        annot_matrix_ordered = False
        print(f"Disabled annotations (>20 genes)")
    elif annot:
        annot = annot_matrix_ordered

    # Create heatmap
    sns.heatmap(
        corr_df_ordered,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        center=0,
        square=True,
        linewidths=0.5,
        cbar_kws={"shrink": 0.8, "label": f"{method.capitalize()} correlation"},
        annot=annot,
        fmt='',  # Use empty format since we're providing custom strings
        ax=ax,
        annot_kws={'fontsize': 8}
    )

    # Set title
    if title is None:
        title = f"Gene Coexpression ({method.capitalize()} correlation)"

    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

    # Adjust layout
    plt.tight_layout()

    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to: {output_path}")

    plt.close()

    # Generate interactive HTML version if requested
    if generate_html:
        plot_interactive_heatmap(
            corr_df=corr_df,
            pval_df=pval_df,
            output_path=output_path,
            method=method,
            title=title,
            cluster=cluster
        )


def generate_summary_report(
    adata,
    genes_requested: List[str],
    corr_df: pd.DataFrame,
    pval_df: pd.DataFrame,
    method: str,
    filters: Dict[str, Any],
    output_dir: Path,
    prefix: str
):
    """
    Generate summary report of coexpression analysis.

    Args:
        adata: AnnData object used
        genes_requested: Original gene list requested
        corr_df: Correlation matrix DataFrame
        pval_df: P-value matrix DataFrame
        method: Correlation method used
        filters: Dictionary of filters applied
        output_dir: Output directory
        prefix: Output file prefix
    """
    report_path = output_dir / f"{prefix}_summary.txt"

    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("GENE COEXPRESSION ANALYSIS SUMMARY\n")
        f.write("=" * 80 + "\n\n")

        f.write("Analysis Parameters:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Correlation method: {method}\n")
        f.write(f"Genes requested: {len(genes_requested)}\n")
        f.write(f"Genes analyzed: {len(corr_df)}\n")
        f.write(f"Genes missing: {len(genes_requested) - len(corr_df)}\n")
        f.write(f"Total cells: {adata.n_obs:,}\n")
        f.write(f"Total genes in dataset: {adata.n_vars:,}\n")

        f.write("\nFilters Applied:\n")
        f.write("-" * 40 + "\n")
        if filters['cell_type']:
            f.write(f"Cell type: {filters['cell_type']}\n")
        if filters['tissue']:
            f.write(f"Tissue: {filters['tissue']}\n")
        if filters['disease']:
            f.write(f"Disease: {filters['disease']}\n")
        if filters['metadata_filter']:
            f.write(f"Custom filter: {filters['metadata_filter']}\n")
        if not any(filters.values()):
            f.write("None\n")

        f.write("\nCorrelation Statistics:\n")
        f.write("-" * 40 + "\n")

        mask = np.triu(np.ones_like(corr_df.values, dtype=bool), k=1)
        upper_triangle_corr = corr_df.values[mask]
        upper_triangle_pval = pval_df.values[mask]

        f.write(f"Mean: {upper_triangle_corr.mean():.4f}\n")
        f.write(f"Median: {np.median(upper_triangle_corr):.4f}\n")
        f.write(f"Std: {upper_triangle_corr.std():.4f}\n")
        f.write(f"Min: {upper_triangle_corr.min():.4f}\n")
        f.write(f"Max: {upper_triangle_corr.max():.4f}\n")

        f.write("\nP-value Statistics:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Mean p-value: {upper_triangle_pval.mean():.4f}\n")
        f.write(f"Median p-value: {np.median(upper_triangle_pval):.4f}\n")
        n_sig_0001 = (upper_triangle_pval < 0.0001).sum()
        n_sig_001 = (upper_triangle_pval < 0.001).sum()
        n_sig_01 = (upper_triangle_pval < 0.01).sum()
        n_sig_05 = (upper_triangle_pval < 0.05).sum()
        n_total = len(upper_triangle_pval)
        f.write(f"Significant at p < 0.0001: {n_sig_0001} / {n_total} ({100*n_sig_0001/n_total:.1f}%)\n")
        f.write(f"Significant at p < 0.001:  {n_sig_001} / {n_total} ({100*n_sig_001/n_total:.1f}%)\n")
        f.write(f"Significant at p < 0.01:   {n_sig_01} / {n_total} ({100*n_sig_01/n_total:.1f}%)\n")
        f.write(f"Significant at p < 0.05:   {n_sig_05} / {n_total} ({100*n_sig_05/n_total:.1f}%)\n")

        f.write("\nGenes Analyzed:\n")
        f.write("-" * 40 + "\n")
        for gene in corr_df.index:
            f.write(f"  {gene}\n")

        if len(genes_requested) > len(corr_df):
            f.write("\nGenes Not Found:\n")
            f.write("-" * 40 + "\n")
            missing = set(genes_requested) - set(corr_df.index)
            for gene in sorted(missing):
                f.write(f"  {gene}\n")

        f.write("\nTop 20 Strongest Correlations (with p-values):\n")
        f.write("-" * 40 + "\n")
        f.write("Significance: **** p<0.0001, *** p<0.001, ** p<0.01, * p<0.05\n\n")

        corr_pairs = []
        n = len(corr_df)
        for i in range(n):
            for j in range(i+1, n):
                corr_pairs.append((corr_df.index[i], corr_df.index[j],
                                 corr_df.values[i, j], pval_df.values[i, j]))

        corr_pairs.sort(key=lambda x: abs(x[2]), reverse=True)
        for rank, (gene1, gene2, corr, pval) in enumerate(corr_pairs[:20], 1):
            sig_marker = "****" if pval < 0.0001 else "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
            f.write(f"{rank:3d}. {gene1:15s} <-> {gene2:15s}: {corr:7.4f} (p={pval:.2e}) {sig_marker}\n")

        f.write("\n" + "=" * 80 + "\n")

    print(f"Summary report saved to: {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze gene coexpression from AnnData objects',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis with gene list file
  python analyze_coexpression.py \\
    --input data.h5ad \\
    --genes my_genes.txt \\
    --output coexpression

  # With comma-separated gene list
  python analyze_coexpression.py \\
    --input data.h5ad \\
    --genes "APOE,APOC1,APOC2,APOC3" \\
    --output apoe_coexp

  # Filter by cell type
  python analyze_coexpression.py \\
    --input lung_data.h5ad \\
    --genes immune_genes.txt \\
    --cell-type "T cell" \\
    --output tcell_coexp

  # Spearman correlation with clustering
  python analyze_coexpression.py \\
    --input data.h5ad \\
    --genes genes.txt \\
    --method spearman \\
    --cluster \\
    --output spearman_coexp

  # Custom metadata filter
  python analyze_coexpression.py \\
    --input data.h5ad \\
    --genes genes.txt \\
    --metadata-filter "sex == 'female' and disease == 'normal'" \\
    --output female_normal_coexp
        """
    )

    # Input/output arguments
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Input AnnData file (.h5ad)')
    parser.add_argument('--genes', '-g', type=str, required=True,
                        help='Gene list file (one per line) or comma-separated gene names')
    parser.add_argument('--output', '-o', type=str, default='coexpression',
                        help='Output prefix for files')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Output directory')

    # Analysis parameters
    parser.add_argument('--method', type=str, default='pearson',
                        choices=['pearson', 'spearman'],
                        help='Correlation method (default: pearson)')
    parser.add_argument('--use-raw', action='store_true',
                        help='Use .raw.X instead of .X')
    parser.add_argument('--min-cells', type=int, default=10,
                        help='Minimum cells expressing each gene (default: 10)')

    # Filtering parameters
    parser.add_argument('--cell-type', type=str,
                        help='Filter by cell type (partial match)')
    parser.add_argument('--tissue', type=str,
                        help='Filter by tissue (partial match)')
    parser.add_argument('--disease', type=str,
                        help='Filter by disease (partial match)')
    parser.add_argument('--metadata-filter', type=str,
                        help='Custom metadata filter expression (e.g., "sex == \'female\'")')

    # Visualization parameters
    parser.add_argument('--figsize', type=float, nargs=2, default=[10, 8],
                        help='Figure size (width height) (default: 10 8)')
    parser.add_argument('--cmap', type=str, default='RdBu_r',
                        help='Colormap for heatmap (default: RdBu_r)')
    parser.add_argument('--cluster', action='store_true',
                        help='Perform hierarchical clustering')
    parser.add_argument('--no-annot', action='store_true',
                        help='Disable correlation value annotations')
    parser.add_argument('--title', type=str,
                        help='Custom plot title')

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load AnnData
    print("=" * 80)
    print("GENE COEXPRESSION ANALYSIS")
    print("=" * 80)
    print(f"\nLoading AnnData from: {args.input}")

    try:
        adata = sc.read_h5ad(args.input)
        print(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    except Exception as e:
        print(f"Error loading AnnData: {e}")
        sys.exit(1)

    # Load gene list
    genes = load_gene_list(args.genes)

    if len(genes) < 2:
        print(f"Error: Need at least 2 genes for coexpression analysis. Found: {len(genes)}")
        sys.exit(1)

    # Filter AnnData
    if args.cell_type or args.tissue or args.disease or args.metadata_filter:
        print("\nApplying filters...")
        adata = filter_adata(
            adata,
            cell_type=args.cell_type,
            tissue=args.tissue,
            disease=args.disease,
            metadata_filter=args.metadata_filter
        )

        if adata.n_obs == 0:
            print("Error: No cells remaining after filtering")
            sys.exit(1)

    # Compute correlation matrix and p-values
    corr_df, pval_df = compute_correlation_matrix(
        adata,
        genes,
        method=args.method,
        use_raw=args.use_raw,
        min_cells=args.min_cells
    )

    # Save correlation matrix
    corr_path = output_dir / f"{args.output}_correlation_matrix.csv"
    corr_df.to_csv(corr_path)
    print(f"\nCorrelation matrix saved to: {corr_path}")

    # Save p-value matrix
    pval_path = output_dir / f"{args.output}_pvalue_matrix.csv"
    pval_df.to_csv(pval_path)
    print(f"P-value matrix saved to: {pval_path}")

    # Create heatmap with date in filename
    from datetime import datetime
    date_str = datetime.now().strftime("%Y%m%d")
    heatmap_path = output_dir / f"{args.output}_heatmap_{date_str}.png"
    plot_correlation_heatmap(
        corr_df,
        pval_df,
        heatmap_path,
        method=args.method,
        title=args.title,
        figsize=tuple(args.figsize),
        cmap=args.cmap,
        cluster=args.cluster,
        annot=not args.no_annot
    )

    # Generate summary report
    filters = {
        'cell_type': args.cell_type,
        'tissue': args.tissue,
        'disease': args.disease,
        'metadata_filter': args.metadata_filter
    }

    generate_summary_report(
        adata,
        genes,
        corr_df,
        pval_df,
        args.method,
        filters,
        output_dir,
        args.output
    )

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nOutput files:")
    print(f"  1. Correlation matrix: {corr_path}")
    print(f"  2. P-value matrix: {pval_path}")
    print(f"  3. Heatmap: {heatmap_path}")
    print(f"  4. Summary report: {output_dir / f'{args.output}_summary.txt'}")
    print()


if __name__ == "__main__":
    main()
