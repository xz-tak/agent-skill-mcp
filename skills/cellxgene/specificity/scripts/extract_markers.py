#!/usr/bin/env python3
"""
CellGuide Marker Gene Extraction Tool

Extract computational and canonical marker genes from CellGuide database
with optional filtering by organism and tissue.

Usage examples:
    python extract_markers.py --cell-type "B cell" --output bcell_markers
    python extract_markers.py --cell-type "T cell" --organism "Homo sapiens" --output tcell_human
    python extract_markers.py --ontology-id CL_0000236 --tissue spleen --format excel --output markers
"""

import argparse
import requests
import pandas as pd
import sys
import os

# CellGuide API configuration
CELL_GUIDE_BASE_URI = "https://cellguide.cellxgene.cziscience.com"
LATEST_SNAPSHOT = requests.get(f"{CELL_GUIDE_BASE_URI}/latest_snapshot_identifier").text


def _get_cellguide_file(relpth: str, snapshot=LATEST_SNAPSHOT):
    """Internal function to fetch CellGuide files"""
    req = requests.get(f"{CELL_GUIDE_BASE_URI}/{snapshot}/{relpth}")
    if req.text == "":
        raise ValueError(f"No record found for {snapshot}/{relpth}")
    return req


def get_computational_marker_genes(ontology_id, *, snapshot=LATEST_SNAPSHOT):
    """
    Fetch computational marker genes for a cell type.

    Parameters:
    -----------
    ontology_id : str
        Cell Ontology ID in format 'CL_XXXXXXX'
    snapshot : str, optional
        CellGuide snapshot version

    Returns:
    --------
    pd.DataFrame
        DataFrame with marker genes and their statistics
    """
    resp = _get_cellguide_file(f"computational_marker_genes/{ontology_id}.json", snapshot=snapshot)
    res_df = pd.DataFrame.from_records(resp.json())
    expanded_cols = pd.json_normalize(res_df['groupby_dims'])
    res_df = res_df.join(expanded_cols).drop(columns=['groupby_dims'])
    return res_df


def get_canonical_marker_genes(ontology_id, *, snapshot=LATEST_SNAPSHOT):
    """
    Fetch canonical (literature-curated) marker genes for a cell type.

    Parameters:
    -----------
    ontology_id : str
        Cell Ontology ID in format 'CL_XXXXXXX'
    snapshot : str, optional
        CellGuide snapshot version

    Returns:
    --------
    pd.DataFrame
        DataFrame with canonical marker genes
    """
    resp = _get_cellguide_file(f"canonical_marker_genes/{ontology_id}.json", snapshot=snapshot)
    res_df = pd.DataFrame.from_records(resp.json())
    return res_df


def get_cell_ontology_id(cell_name):
    """
    Convert cell type name to Cell Ontology ID using OLS API.

    Parameters:
    -----------
    cell_name : str
        Cell type name (e.g., 'B cell', 'T cell', 'neuron')

    Returns:
    --------
    str
        Cell Ontology ID in format 'CL:XXXXXXX'
    """
    base_url = "https://www.ebi.ac.uk/ols/api/search"
    params = {
        "q": cell_name,
        "ontology": "cl",
        "exact": "true",
        "rows": 1
    }

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()

        if data["response"]["numFound"] > 0:
            result = data["response"]["docs"][0]
            ontology_id = result.get("obo_id")
            if ontology_id:
                return ontology_id
        return None
    except requests.exceptions.RequestException as e:
        print(f"Error making API request: {e}")
        raise
    except (KeyError, IndexError) as e:
        print(f"Error parsing API response: {e}")
        return None


def filter_markers(df, organism=None, tissue=None):
    """
    Filter computational marker genes by organism and/or tissue.

    Parameters:
    -----------
    df : pd.DataFrame
        Computational marker genes DataFrame
    organism : str or list, optional
        Organism name(s) to filter (e.g., 'Homo sapiens', 'Mus musculus')
    tissue : str or list, optional
        Tissue name(s) to filter

    Returns:
    --------
    pd.DataFrame
        Filtered DataFrame
    """
    filtered_df = df.copy()

    if organism is not None:
        if isinstance(organism, str):
            organism = [organism]
        filtered_df = filtered_df[filtered_df['organism_ontology_term_label'].isin(organism)]

    if tissue is not None:
        if isinstance(tissue, str):
            tissue = [tissue]
        filtered_df = filtered_df[filtered_df['tissue_ontology_term_label'].isin(tissue)]

    return filtered_df


def export_markers(df, output_path, file_format='csv', index=False):
    """
    Export marker genes to file.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame to export
    output_path : str
        Output file path (without extension)
    file_format : str
        Output format: 'csv', 'tsv', or 'excel'
    index : bool
        Whether to include row index
    """
    if file_format == 'csv':
        df.to_csv(f"{output_path}.csv", index=index)
        print(f"Exported to {output_path}.csv")
    elif file_format == 'tsv':
        df.to_csv(f"{output_path}.tsv", sep='\t', index=index)
        print(f"Exported to {output_path}.tsv")
    elif file_format == 'excel':
        df.to_excel(f"{output_path}.xlsx", index=index)
        print(f"Exported to {output_path}.xlsx")
    else:
        raise ValueError(f"Unsupported format: {file_format}")


def get_marker_summary(df):
    """
    Generate summary statistics for marker genes.

    Parameters:
    -----------
    df : pd.DataFrame
        Computational marker genes DataFrame

    Returns:
    --------
    dict
        Summary statistics
    """
    summary = {
        'total_markers': len(df),
        'unique_genes': df['symbol'].nunique(),
        'mean_marker_score': df['marker_score'].mean(),
        'mean_specificity': df['specificity'].mean(),
        'mean_expression': df['me'].mean(),
        'organisms': df['organism_ontology_term_label'].value_counts().to_dict(),
    }

    tissue_counts = df['tissue_ontology_term_label'].value_counts(dropna=False).to_dict()
    summary['tissues'] = tissue_counts

    return summary


def main():
    parser = argparse.ArgumentParser(
        description='Extract marker genes from CellGuide database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract B cell markers (default: marker_score >= 0.5)
  python extract_markers.py --cell-type "B cell" --output bcell_markers

  # Filter by human organism (default: marker_score >= 0.5)
  python extract_markers.py --cell-type "T cell" --organism "Homo sapiens" --output tcell_human

  # Use custom marker_score cutoff
  python extract_markers.py --cell-type "neuron" --marker-score-cutoff 1.0 --output neuron_high_score

  # Filter by tissue
  python extract_markers.py --cell-type "B cell" --tissue "spleen" --output bcell_spleen

  # Combined filtering with summary
  python extract_markers.py --cell-type "macrophage" --organism "Homo sapiens" --tissue "lung" --summary --output macro_human_lung

  # Get top N markers (disables marker_score cutoff)
  python extract_markers.py --cell-type "dendritic cell" --top-n 100 --output dc_top100

  # Use Cell Ontology ID directly
  python extract_markers.py --ontology-id CL_0000236 --output CL0000236_markers

  # Export both computational and canonical markers
  python extract_markers.py --cell-type "monocyte" --canonical --output monocyte_all
        """
    )

    # Cell type specification
    cell_group = parser.add_mutually_exclusive_group(required=True)
    cell_group.add_argument('--cell-type', type=str,
                           help='Cell type name (e.g., "B cell", "neuron")')
    cell_group.add_argument('--ontology-id', type=str,
                           help='Cell Ontology ID (e.g., CL_0000236 or CL:0000236)')

    # Filtering options
    parser.add_argument('--organism', type=str, nargs='+',
                       help='Filter by organism(s) (e.g., "Homo sapiens" "Mus musculus")')
    parser.add_argument('--tissue', type=str, nargs='+',
                       help='Filter by tissue(s)')

    # Output options
    parser.add_argument('--output', type=str, required=True,
                       help='Output file base name (without extension)')
    parser.add_argument('--format', type=str, default='csv',
                       choices=['csv', 'tsv', 'excel'],
                       help='Output file format (default: csv)')
    parser.add_argument('--marker-score-cutoff', type=float,
                       help='Filter genes with marker_score above this cutoff (default: 0.5 if --top-n not specified)')
    parser.add_argument('--top-n', type=int,
                       help='Export only top N markers by marker_score (disables marker_score cutoff filtering)')

    # Marker type options
    parser.add_argument('--canonical', action='store_true',
                       help='Also export canonical markers')
    parser.add_argument('--summary', action='store_true',
                       help='Print summary statistics')

    args = parser.parse_args()

    try:
        # Get ontology ID
        if args.cell_type:
            print(f"Looking up Cell Ontology ID for '{args.cell_type}'...")
            cell_id = get_cell_ontology_id(args.cell_type)
            if not cell_id:
                print(f"Error: Could not find Cell Ontology ID for '{args.cell_type}'")
                sys.exit(1)
            ontology_id = cell_id.replace(':', '_')
            print(f"Found: {cell_id} ({ontology_id})")
        else:
            ontology_id = args.ontology_id.replace(':', '_')
            print(f"Using Cell Ontology ID: {ontology_id}")

        # Fetch computational markers
        print("\nFetching computational marker genes...")
        comp_markers = get_computational_marker_genes(ontology_id)
        print(f"Retrieved {len(comp_markers)} computational markers")

        # Apply filters
        if args.organism or args.tissue:
            print("Applying filters...")
            comp_markers = filter_markers(
                comp_markers,
                organism=args.organism,
                tissue=args.tissue
            )
            print(f"After filtering: {len(comp_markers)} markers")

        # Apply top-N or marker_score cutoff filtering (mutually exclusive)
        if args.top_n:
            # User explicitly requested top-N: use it and disable cutoff
            comp_markers = comp_markers.nlargest(args.top_n, 'marker_score')
            print(f"Selected top {len(comp_markers)} markers by marker_score")
        else:
            # Default behavior: filter by marker_score cutoff
            cutoff = args.marker_score_cutoff if args.marker_score_cutoff is not None else 0.5
            original_count = len(comp_markers)
            comp_markers = comp_markers[comp_markers['marker_score'] >= cutoff]
            print(f"Applied marker_score cutoff >= {cutoff}: {len(comp_markers)} markers (removed {original_count - len(comp_markers)})")

        # Export computational markers
        output_name = f"{args.output}_computational"
        export_markers(comp_markers, output_name, file_format=args.format)

        # Summary
        if args.summary:
            summary = get_marker_summary(comp_markers)
            print("\n=== Summary Statistics ===")
            print(f"Total markers: {summary['total_markers']}")
            print(f"Unique genes: {summary['unique_genes']}")
            print(f"Mean marker score: {summary['mean_marker_score']:.3f}")
            print(f"Mean specificity: {summary['mean_specificity']:.3f}")
            print(f"Mean expression: {summary['mean_expression']:.3f}")
            print("\nOrganism distribution:")
            for org, count in sorted(summary['organisms'].items()):
                print(f"  {org}: {count}")

            # Show top tissues if any
            tissues = {k: v for k, v in summary['tissues'].items() if pd.notna(k)}
            if tissues:
                print("\nTop 5 tissues:")
                for tissue, count in sorted(tissues.items(), key=lambda x: x[1], reverse=True)[:5]:
                    print(f"  {tissue}: {count}")

        # Fetch canonical markers if requested
        if args.canonical:
            print("\nFetching canonical marker genes...")
            canon_markers = get_canonical_marker_genes(ontology_id)
            print(f"Retrieved {len(canon_markers)} canonical markers")

            # Export canonical markers
            output_name = f"{args.output}_canonical"
            export_markers(canon_markers, output_name, file_format=args.format)

        print("\n✓ Done!")

    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
