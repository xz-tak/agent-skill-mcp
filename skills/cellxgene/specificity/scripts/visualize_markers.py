#!/usr/bin/env python3
"""
Marker Gene Visualization Tool

Create dotplots and other visualizations for marker genes extracted from CellGuide.
Dotplot features:
- Color intensity indicates effect size (marker_score or specificity)
- Dot size indicates percentage of cells expressing (pc)

Usage:
    python visualize_markers.py --input markers.csv --output dotplot.png
    python visualize_markers.py --input markers.csv --top-n 20 --group-by organism
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import sys
from pathlib import Path

# Optional Plotly support for interactive visualizations
try:
    import plotly.graph_objects as go
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("Note: Plotly not available. Install with 'pip install plotly' for interactive HTML visualizations.")


def load_marker_data(input_file):
    """Load marker gene data from CSV file"""
    try:
        df = pd.read_csv(input_file)
        return df
    except Exception as e:
        print(f"Error loading file {input_file}: {e}")
        sys.exit(1)


def create_dotplot_interactive(df, output_file, top_n=None, group_by=None,
                               color_by='marker_score', size_by='pc',
                               title=None, filter_organism=None,
                               filter_tissue=None):
    """Create interactive HTML dotplot using Plotly"""
    if not PLOTLY_AVAILABLE:
        return

    # Apply filters
    plot_df = df.copy()

    if filter_organism:
        if isinstance(filter_organism, str):
            filter_organism = [filter_organism]
        plot_df = plot_df[plot_df['organism_ontology_term_label'].isin(filter_organism)]

    if filter_tissue:
        if isinstance(filter_tissue, str):
            filter_tissue = [filter_tissue]
        plot_df = plot_df[plot_df['tissue_ontology_term_label'].isin(filter_tissue)]

    if len(plot_df) == 0:
        return

    # Get top N genes
    if top_n:
        plot_df = plot_df.nlargest(top_n, 'marker_score')

    # Ensure required columns exist
    required_cols = ['symbol', color_by, size_by]
    missing_cols = [col for col in required_cols if col not in plot_df.columns]
    if missing_cols:
        return

    # Sort by marker_score for better visualization
    plot_df = plot_df.sort_values('marker_score', ascending=True)

    # Create hover text
    hover_text = []
    for _, row in plot_df.iterrows():
        text = f"Gene: {row['symbol']}<br>"
        text += f"Marker Score: {row.get('marker_score', 'N/A'):.3f}<br>"
        text += f"Specificity: {row.get('specificity', 'N/A'):.3f}<br>"
        text += f"% Cells: {row.get('pc', 'N/A'):.1f}<br>"
        if group_by and group_by in row:
            text += f"{group_by}: {row[group_by]}"
        hover_text.append(text)

    # Create figure
    fig = go.Figure()

    if group_by and group_by in plot_df.columns:
        # Group by specified column
        groups = plot_df[group_by].unique()
        for i, grp in enumerate(groups):
            grp_df = plot_df[plot_df[group_by] == grp]
            grp_hover = [hover_text[j] for j, idx in enumerate(plot_df.index) if idx in grp_df.index]

            fig.add_trace(go.Scatter(
                x=grp_df[color_by],
                y=grp_df['symbol'],
                mode='markers',
                name=str(grp),
                marker=dict(
                    size=grp_df[size_by] * 2,  # Scale size for visibility
                    color=grp_df[color_by],
                    colorscale='Reds',
                    showscale=(i == 0),
                    colorbar=dict(title=color_by) if i == 0 else None,
                    line=dict(width=0.5, color='white')
                ),
                hovertext=grp_hover,
                hoverinfo='text'
            ))
    else:
        fig.add_trace(go.Scatter(
            x=plot_df[color_by],
            y=plot_df['symbol'],
            mode='markers',
            marker=dict(
                size=plot_df[size_by] * 2,  # Scale size for visibility
                color=plot_df[color_by],
                colorscale='Reds',
                showscale=True,
                colorbar=dict(title=color_by),
                line=dict(width=0.5, color='white')
            ),
            hovertext=hover_text,
            hoverinfo='text'
        ))

    # Update layout
    fig.update_layout(
        title=title or "Marker Gene Expression",
        xaxis_title=color_by,
        yaxis_title="Gene Symbol",
        height=max(400, len(plot_df) * 20),
        hovermode='closest',
        template='plotly_white'
    )

    # Save HTML
    html_output = Path(output_file).with_suffix('.html')
    fig.write_html(str(html_output))
    print(f"✓ Interactive dotplot saved to {html_output}")


def create_heatmap_interactive(df, output_file, top_n=None, title=None):
    """Create interactive HTML heatmap using Plotly"""
    if not PLOTLY_AVAILABLE:
        return

    # Get top N genes
    if top_n:
        df = df.nlargest(top_n, 'marker_score')

    # Check required columns
    if 'symbol' not in df.columns or 'marker_score' not in df.columns:
        return

    # Pivot data for heatmap (if we have grouping columns)
    # For simple case, create a 1D heatmap
    genes = df['symbol'].values
    scores = df['marker_score'].values.reshape(1, -1)

    # Create hover text
    hover_text = []
    for gene, score in zip(genes, df['marker_score']):
        hover_text.append(f"Gene: {gene}<br>Marker Score: {score:.3f}")
    hover_array = np.array(hover_text).reshape(1, -1)

    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=scores,
        x=genes,
        y=['Marker Score'],
        colorscale='Reds',
        text=hover_array,
        hoverinfo='text',
        colorbar=dict(title="Marker Score")
    ))

    # Update layout
    fig.update_layout(
        title=title or "Marker Gene Heatmap",
        xaxis_title="Gene Symbol",
        height=300,
        template='plotly_white'
    )

    # Save HTML
    html_output = Path(output_file).with_suffix('.html')
    fig.write_html(str(html_output))
    print(f"✓ Interactive heatmap saved to {html_output}")


def create_barplot_interactive(df, output_file, top_n=None, title=None):
    """Create interactive HTML barplot using Plotly"""
    if not PLOTLY_AVAILABLE:
        return

    # Get top N genes
    if top_n:
        df = df.nlargest(top_n, 'marker_score')

    # Check required columns
    if 'symbol' not in df.columns or 'marker_score' not in df.columns:
        return

    # Sort by marker_score
    df = df.sort_values('marker_score', ascending=True)

    # Create hover text
    hover_text = []
    for _, row in df.iterrows():
        text = f"Gene: {row['symbol']}<br>"
        text += f"Marker Score: {row['marker_score']:.3f}<br>"
        text += f"Specificity: {row.get('specificity', 'N/A'):.3f}"
        hover_text.append(text)

    # Create barplot
    fig = go.Figure(data=go.Bar(
        x=df['marker_score'],
        y=df['symbol'],
        orientation='h',
        marker=dict(
            color=df.get('specificity', df['marker_score']),
            colorscale='Reds',
            showscale=True,
            colorbar=dict(title="Specificity")
        ),
        hovertext=hover_text,
        hoverinfo='text'
    ))

    # Update layout
    fig.update_layout(
        title=title or "Marker Gene Scores",
        xaxis_title="Marker Score",
        yaxis_title="Gene Symbol",
        height=max(400, len(df) * 20),
        template='plotly_white'
    )

    # Save HTML
    html_output = Path(output_file).with_suffix('.html')
    fig.write_html(str(html_output))
    print(f"✓ Interactive barplot saved to {html_output}")


def create_dotplot(df, output_file, top_n=None, group_by=None,
                   color_by='marker_score', size_by='pc',
                   figsize=None, title=None, filter_organism=None,
                   filter_tissue=None, generate_static=False):
    """
    Create dotplot visualization of marker genes.

    By default, generates interactive HTML visualization. Set generate_static=True to also create PNG.

    Parameters:
    -----------
    df : pd.DataFrame
        Marker gene data
    output_file : str
        Output file path
    top_n : int, optional
        Show only top N genes by marker_score
    group_by : str, optional
        Group genes by this column (e.g., 'organism_ontology_term_label', 'tissue_ontology_term_label')
    color_by : str
        Column to use for color intensity (default: 'marker_score')
    size_by : str
        Column to use for dot size (default: 'pc')
    figsize : tuple, optional
        Figure size (width, height)
    title : str, optional
        Plot title
    filter_organism : str or list, optional
        Filter to specific organism(s)
    filter_tissue : str or list, optional
        Filter to specific tissue(s)
    generate_static : bool, optional
        If True, also generate static PNG image (default: False)
    """
    # Apply filters
    plot_df = df.copy()

    if filter_organism:
        if isinstance(filter_organism, str):
            filter_organism = [filter_organism]
        plot_df = plot_df[plot_df['organism_ontology_term_label'].isin(filter_organism)]

    if filter_tissue:
        if isinstance(filter_tissue, str):
            filter_tissue = [filter_tissue]
        plot_df = plot_df[plot_df['tissue_ontology_term_label'].isin(filter_tissue)]

    if len(plot_df) == 0:
        print("Error: No data remaining after filtering")
        sys.exit(1)

    # Get top N genes
    if top_n:
        plot_df = plot_df.nlargest(top_n, 'marker_score')

    # Ensure required columns exist
    required_cols = ['symbol', color_by, size_by]
    missing_cols = [col for col in required_cols if col not in plot_df.columns]
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}")
        sys.exit(1)

    # Calculate color and size ranges from filtered data
    color_min = plot_df[color_by].min()
    color_max = plot_df[color_by].max()
    size_min = plot_df[size_by].min()
    size_max = plot_df[size_by].max()

    # Normalize sizes for better visual distinction
    # Map the size range to a visual range (e.g., 50 to 800 points)
    size_range = size_max - size_min
    if size_range > 0:
        # Normalize to 0-1, then scale to visual range
        def normalize_size(val):
            normalized = (val - size_min) / size_range
            return 50 + normalized * 750  # Range: 50 to 800
    else:
        # All values are the same, use medium size
        def normalize_size(val):
            return 400

    # Generate interactive HTML version (default)
    create_dotplot_interactive(
        df=df,
        output_file=output_file,
        top_n=top_n,
        group_by=group_by,
        color_by=color_by,
        size_by=size_by,
        title=title,
        filter_organism=filter_organism,
        filter_tissue=filter_tissue
    )

    # Only generate static PNG if requested
    if not generate_static:
        return

    # Group data if requested
    if group_by and group_by in plot_df.columns:
        # Create a multi-level plot
        groups = plot_df[group_by].unique()
        groups = [g for g in groups if pd.notna(g)]

        if not figsize:
            figsize = (12, max(8, len(plot_df) * 0.3))

        fig, ax = plt.subplots(figsize=figsize)

        # Prepare data for grouped dotplot
        y_pos = 0
        y_labels = []
        y_positions = []
        group_boundaries = []

        for group in sorted(groups):
            group_data = plot_df[plot_df[group_by] == group].nlargest(
                min(len(plot_df[plot_df[group_by] == group]), top_n if top_n else float('inf')),
                'marker_score'
            )

            for _, row in group_data.iterrows():
                y_labels.append(f"{row['symbol']}")
                y_positions.append(y_pos)

                # Plot the dot
                color_val = row[color_by]
                size_val = normalize_size(row[size_by])

                scatter = ax.scatter(0, y_pos, s=size_val, c=[color_val],
                                   cmap='Reds', vmin=color_min, vmax=color_max,
                                   alpha=0.8, edgecolors='black', linewidth=0.5)

                y_pos += 1

            group_boundaries.append(y_pos)
            y_pos += 0.5  # Add space between groups

        # Add group separators and labels
        for i, (group, boundary) in enumerate(zip(sorted(groups), group_boundaries)):
            start = group_boundaries[i-1] + 0.5 if i > 0 else 0
            ax.axhline(y=boundary - 0.5, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            mid_point = (start + boundary - 0.5) / 2
            ax.text(0.5, mid_point, str(group), fontsize=9, color='blue',
                   rotation=0, ha='left', va='center', weight='bold')

        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=8)
        ax.set_xticks([0])
        ax.set_xticklabels(['Markers'])
        ax.set_xlim(-0.3, 0.8)

    else:
        # Simple dotplot without grouping
        if not figsize:
            figsize = (10, max(6, len(plot_df) * 0.25))

        fig, ax = plt.subplots(figsize=figsize)

        # Sort by marker score
        plot_df = plot_df.sort_values('marker_score', ascending=True)

        # Create scatter plot
        y_pos = range(len(plot_df))
        sizes = [normalize_size(val) for val in plot_df[size_by]]
        colors = plot_df[color_by]

        scatter = ax.scatter([0] * len(plot_df), y_pos, s=sizes, c=colors,
                           cmap='Reds', vmin=color_min, vmax=color_max,
                           alpha=0.8, edgecolors='black', linewidth=0.5)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot_df['symbol'], fontsize=8)
        ax.set_xticks([0])
        ax.set_xticklabels(['Markers'])
        ax.set_xlim(-0.3, 0.3)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label(color_by.replace('_', ' ').title(), rotation=270, labelpad=20)

    # Add size legend based on actual data range
    # Create 4 evenly spaced values from min to max
    if size_range > 0:
        size_values_data = [
            size_min,
            size_min + size_range * 0.33,
            size_min + size_range * 0.67,
            size_max
        ]
    else:
        size_values_data = [size_min]

    # Create legend with normalized sizes
    size_legend = [plt.scatter([], [], s=normalize_size(s), c='gray', alpha=0.6,
                              edgecolors='black', linewidth=0.5)
                  for s in size_values_data]

    # Format labels based on the size_by column
    if size_by == 'pc':
        labels = [f'{s:.2f}' for s in size_values_data]
    else:
        labels = [f'{s:.1f}' for s in size_values_data]

    legend1 = ax.legend(size_legend, labels, title=size_by.replace('_', ' ').title(),
                       loc='upper right', frameon=True, fontsize=8, title_fontsize=9)
    ax.add_artist(legend1)

    # Set title
    if title:
        ax.set_title(title, fontsize=14, weight='bold', pad=20)
    else:
        ax.set_title('Marker Gene Expression', fontsize=14, weight='bold', pad=20)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Dotplot saved to {output_file}")
    plt.close()


def create_heatmap(df, output_file, top_n=20, cluster=True, figsize=None, title=None, generate_static=False):
    """
    Create heatmap visualization of marker genes across contexts.

    By default, generates interactive HTML visualization. Set generate_static=True to also create PNG.

    Parameters:
    -----------
    df : pd.DataFrame
        Marker gene data
    output_file : str
        Output file path
    top_n : int
        Show top N genes by marker_score
    cluster : bool
        Whether to cluster genes
    figsize : tuple, optional
        Figure size
    title : str, optional
        Plot title
    generate_static : bool, optional
        If True, also generate static PNG image (default: False)
    """
    # Generate interactive HTML version (default)
    create_heatmap_interactive(
        df=df,
        output_file=output_file,
        top_n=top_n,
        title=title
    )

    # Only generate static PNG if requested
    if not generate_static:
        return

    # Get top genes
    top_genes = df.nlargest(top_n, 'marker_score')

    # Create pivot table for heatmap
    if 'tissue_ontology_term_label' in df.columns:
        # Create matrix: genes x tissues
        pivot_data = top_genes.pivot_table(
            values='marker_score',
            index='symbol',
            columns='tissue_ontology_term_label',
            aggfunc='max'
        )
    elif 'organism_ontology_term_label' in df.columns:
        # Create matrix: genes x organisms
        pivot_data = top_genes.pivot_table(
            values='marker_score',
            index='symbol',
            columns='organism_ontology_term_label',
            aggfunc='max'
        )
    else:
        print("Error: Cannot create heatmap without grouping column")
        return

    # Fill NaN with 0
    pivot_data = pivot_data.fillna(0)

    if not figsize:
        figsize = (max(8, len(pivot_data.columns) * 0.8), max(6, len(pivot_data) * 0.4))

    fig, ax = plt.subplots(figsize=figsize)

    # Create heatmap
    sns.heatmap(pivot_data, cmap='Reds', annot=False, fmt='.2f',
               cbar_kws={'label': 'Marker Score'}, ax=ax,
               linewidths=0.5, linecolor='lightgray')

    if title:
        ax.set_title(title, fontsize=14, weight='bold', pad=20)
    else:
        ax.set_title('Marker Gene Expression Heatmap', fontsize=14, weight='bold', pad=20)

    ax.set_xlabel('')
    ax.set_ylabel('Gene Symbol', fontsize=10)
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(fontsize=8)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to {output_file}")
    plt.close()


def create_barplot(df, output_file, top_n=20, figsize=None, title=None, generate_static=False):
    """
    Create barplot of top marker genes by score.

    By default, generates interactive HTML visualization. Set generate_static=True to also create PNG.

    Parameters:
    -----------
    df : pd.DataFrame
        Marker gene data
    output_file : str
        Output file path
    top_n : int
        Show top N genes
    figsize : tuple, optional
        Figure size
    title : str, optional
        Plot title
    generate_static : bool, optional
        If True, also generate static PNG image (default: False)
    """
    # Generate interactive HTML version (default)
    create_barplot_interactive(
        df=df,
        output_file=output_file,
        top_n=top_n,
        title=title
    )

    # Only generate static PNG if requested
    if not generate_static:
        return

    # Get top genes
    top_genes = df.nlargest(top_n, 'marker_score').sort_values('marker_score', ascending=True)

    if not figsize:
        figsize = (10, max(6, len(top_genes) * 0.35))

    fig, ax = plt.subplots(figsize=figsize)

    # Create barplot
    y_pos = range(len(top_genes))
    bars = ax.barh(y_pos, top_genes['marker_score'], color='steelblue', alpha=0.8)

    # Color bars by specificity if available
    if 'specificity' in top_genes.columns:
        colors = plt.cm.Reds(top_genes['specificity'])
        for bar, color in zip(bars, colors):
            bar.set_color(color)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_genes['symbol'], fontsize=8)
    ax.set_xlabel('Marker Score', fontsize=10)
    ax.set_ylabel('Gene Symbol', fontsize=10)

    if title:
        ax.set_title(title, fontsize=14, weight='bold', pad=20)
    else:
        ax.set_title(f'Top {top_n} Marker Genes', fontsize=14, weight='bold', pad=20)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Barplot saved to {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize marker genes from CellGuide data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic dotplot (default: marker_score >= 0.5)
  python visualize_markers.py --input markers.csv --output dotplot.png

  # Top N genes with grouping (disables marker_score cutoff)
  python visualize_markers.py --input markers.csv --top-n 20 --group-by organism_ontology_term_label --output dotplot.png

  # Filter and visualize with custom cutoff
  python visualize_markers.py --input markers.csv --filter-organism "Homo sapiens" --marker-score-cutoff 1.0 --output human_dotplot.png

  # Create heatmap with top N genes
  python visualize_markers.py --input markers.csv --plot-type heatmap --top-n 25 --output heatmap.png

  # Create barplot with custom cutoff
  python visualize_markers.py --input markers.csv --plot-type barplot --marker-score-cutoff 0.8 --output barplot.png
        """
    )

    # Input/output
    parser.add_argument('--input', type=str, required=True,
                       help='Input CSV file with marker data')
    parser.add_argument('--output', type=str, required=True,
                       help='Output file path for plot')

    # Plot options
    parser.add_argument('--plot-type', type=str, default='dotplot',
                       choices=['dotplot', 'heatmap', 'barplot'],
                       help='Type of plot to generate (default: dotplot)')
    parser.add_argument('--marker-score-cutoff', type=float,
                       help='Filter genes with marker_score above this cutoff (default: 0.5 if --top-n not specified)')
    parser.add_argument('--top-n', type=int,
                       help='Number of top genes to show (disables marker_score cutoff filtering)')
    parser.add_argument('--group-by', type=str,
                       choices=['organism_ontology_term_label', 'tissue_ontology_term_label'],
                       help='Group genes by this attribute')

    # Filtering
    parser.add_argument('--filter-organism', type=str, nargs='+',
                       help='Filter to specific organism(s)')
    parser.add_argument('--filter-tissue', type=str, nargs='+',
                       help='Filter to specific tissue(s)')

    # Appearance
    parser.add_argument('--title', type=str,
                       help='Plot title')
    parser.add_argument('--figsize', type=float, nargs=2,
                       help='Figure size as width height (e.g., 12 8)')
    parser.add_argument('--color-by', type=str, default='marker_score',
                       help='Column for color intensity (default: marker_score)')
    parser.add_argument('--size-by', type=str, default='pc',
                       help='Column for dot size (default: pc)')

    args = parser.parse_args()

    # Load data
    print(f"Loading marker data from {args.input}...")
    df = load_marker_data(args.input)
    print(f"Loaded {len(df)} markers for {df['symbol'].nunique()} unique genes")

    # Apply top-N or marker_score cutoff filtering (mutually exclusive)
    if args.top_n:
        # User explicitly requested top-N: use it and disable cutoff
        df = df.nlargest(args.top_n, 'marker_score')
        print(f"Selected top {len(df)} markers by marker_score")
    else:
        # Default behavior: filter by marker_score cutoff
        cutoff = args.marker_score_cutoff if args.marker_score_cutoff is not None else 0.5
        original_count = len(df)
        df = df[df['marker_score'] >= cutoff]
        print(f"Applied marker_score cutoff >= {cutoff}: {len(df)} markers (removed {original_count - len(df)})")

    # Set figure size
    figsize = tuple(args.figsize) if args.figsize else None

    # Create visualization
    print(f"Creating {args.plot_type}...")

    if args.plot_type == 'dotplot':
        create_dotplot(
            df, args.output,
            top_n=None,  # Filtering already applied above
            group_by=args.group_by,
            color_by=args.color_by,
            size_by=args.size_by,
            figsize=figsize,
            title=args.title,
            filter_organism=args.filter_organism,
            filter_tissue=args.filter_tissue
        )
    elif args.plot_type == 'heatmap':
        create_heatmap(df, args.output, top_n=None, figsize=figsize, title=args.title)
    elif args.plot_type == 'barplot':
        create_barplot(df, args.output, top_n=None, figsize=figsize, title=args.title)

    print("✓ Done!")


if __name__ == '__main__':
    main()
