#!/usr/bin/env python3
"""
Multi-gene pathway analysis with UpSet plot visualization.

Query multiple genes across pathway databases (KEGG, Reactome, MSigDB),
visualize shared pathways using UpSet plots, and export results to a table.

Usage:
    python multi_gene_analysis.py TP53 BRCA1 EGFR
    python multi_gene_analysis.py TP53 BRCA1 --output results --msigdb-collections H C2
    python multi_gene_analysis.py gene1 gene2 gene3 --min-intersection 2

Requires: pandas, matplotlib, upsetplot, requests, gseapy
"""

import json
import logging
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Set
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt

try:
    from upsetplot import UpSet, from_contents
except ImportError:
    print("Error: upsetplot not installed. Run: pip install upsetplot")
    sys.exit(1)

# Import pathway query modules
try:
    from kegg_api import get_pathways_for_gene
    from reactome_api import get_reactome_terms_for_gene
    from msigdb_api import get_gene_sets_across_collections
except ImportError as e:
    print(f"Error importing pathway modules: {e}")
    sys.exit(1)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def query_gene_across_databases(
    gene: str,
    kegg_organism: str = "hsa",
    reactome_species: int = 9606,
    msigdb_collections: List[str] = None,
    msigdb_version: str = "2025.1.Hs"
) -> Dict[str, Set[str]]:
    """
    Query a single gene across all databases and return pathway sets.

    Returns:
        Dictionary with keys 'kegg', 'reactome', 'msigdb' containing sets of pathway names
    """
    results = {
        "kegg": set(),
        "reactome": set(),
        "msigdb": set()
    }

    # Query KEGG
    logger.info(f"Querying KEGG for {gene}")
    try:
        kegg_result = get_pathways_for_gene(gene, kegg_organism)
        for gene_data in kegg_result.get("genes_found", []):
            for pathway in gene_data.get("pathways", []):
                pathway_name = pathway.get("pathway_name", "")
                pathway_id = pathway.get("pathway_id", "")
                # Prefer pathway name; use ID only if name is empty
                if pathway_name:
                    results["kegg"].add(pathway_name)
                elif pathway_id:
                    # Use ID as last resort (shouldn't happen with fix)
                    logger.warning(f"Using pathway ID instead of name: {pathway_id}")
                    results["kegg"].add(pathway_id)
        logger.info(f"Found {len(results['kegg'])} KEGG pathways for {gene}")
    except Exception as e:
        logger.error(f"Failed to query KEGG for {gene}: {e}")

    # Query Reactome
    logger.info(f"Querying Reactome for {gene}")
    try:
        reactome_result = get_reactome_terms_for_gene(gene, species=reactome_species)
        for term in reactome_result.get("terms", []):
            term_name = term.get("name", "")
            if term_name:
                results["reactome"].add(term_name)
        logger.info(f"Found {len(results['reactome'])} Reactome pathways for {gene}")
    except Exception as e:
        logger.error(f"Failed to query Reactome for {gene}: {e}")

    # Query MSigDB
    logger.info(f"Querying MSigDB for {gene}")
    try:
        if msigdb_collections is None:
            msigdb_collections = ["H", "C2"]  # Default to Hallmark and Curated

        msigdb_result = get_gene_sets_across_collections(
            gene,
            collections=msigdb_collections,
            dbver=msigdb_version
        )

        for coll_name, coll_data in msigdb_result.get("collections", {}).items():
            if "error" not in coll_data:
                for gene_set in coll_data.get("gene_sets", []):
                    gene_set_name = gene_set.get("gene_set_name", "")
                    if gene_set_name:
                        results["msigdb"].add(gene_set_name)
        logger.info(f"Found {len(results['msigdb'])} MSigDB gene sets for {gene}")
    except Exception as e:
        logger.error(f"Failed to query MSigDB for {gene}: {e}")

    return results


def query_multiple_genes(
    genes: List[str],
    kegg_organism: str = "hsa",
    reactome_species: int = 9606,
    msigdb_collections: List[str] = None,
    msigdb_version: str = "2025.1.Hs"
) -> Dict[str, Dict[str, Set[str]]]:
    """
    Query multiple genes across all databases.

    Returns:
        Nested dictionary: {gene: {database: set of pathways}}
    """
    all_results = {}

    for gene in genes:
        logger.info(f"Processing gene: {gene}")
        all_results[gene] = query_gene_across_databases(
            gene,
            kegg_organism=kegg_organism,
            reactome_species=reactome_species,
            msigdb_collections=msigdb_collections,
            msigdb_version=msigdb_version
        )

    return all_results


def create_upset_plot(
    gene_results: Dict[str, Dict[str, Set[str]]],
    output_prefix: str = "pathway_upset",
    min_subset_size: int = None,
    max_subsets: int = 20
):
    """
    Create UpSet plots for pathway overlap across genes and databases.

    Creates separate plots for each database and one combined plot.
    """
    genes = list(gene_results.keys())

    # Create plots for each database
    for database in ["kegg", "reactome", "msigdb"]:
        logger.info(f"Creating UpSet plot for {database}")

        # Collect pathways for each gene in this database
        gene_pathways = {}
        for gene in genes:
            pathways = gene_results[gene].get(database, set())
            if pathways:
                gene_pathways[gene] = pathways

        if not gene_pathways:
            logger.warning(f"No pathways found for {database}, skipping plot")
            continue

        # Create UpSet plot
        try:
            upset_data = from_contents(gene_pathways)

            # Calculate min_subset_size to limit to max_subsets bars
            if min_subset_size is None:
                subset_counts = upset_data.value_counts()
                if len(subset_counts) > max_subsets:
                    # Sort by count and take the value at position max_subsets
                    min_subset_size_calc = sorted(subset_counts.values(), reverse=True)[max_subsets - 1]
                else:
                    min_subset_size_calc = 1
            else:
                min_subset_size_calc = min_subset_size

            fig = plt.figure(figsize=(14, 10))
            upset = UpSet(
                upset_data,
                subset_size='count',
                intersection_plot_elements=3,
                show_counts='%d',
                min_subset_size=min_subset_size_calc,
                sort_by='cardinality',
                sort_categories_by='cardinality'
            )
            upset.plot(fig=fig)

            fig.suptitle(
                f'Pathway Overlap Across Genes - {database.upper()} (intersections with ≥{min_subset_size_calc} pathways)',
                fontsize=10,
                y=0.98,
                fontweight='bold'
            )

            # Reserve 20% of the figure height free at the top for the title
            fig.tight_layout(rect=[0, 0, 1, 0.8])

            output_file = f"{output_prefix}_{database}.png"
            fig.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.3)
            logger.info(f"Saved {database} UpSet plot to {output_file}")
            plt.close(fig)

        except Exception as e:
            logger.error(f"Failed to create UpSet plot for {database}: {e}")

    # Create combined plot across all databases
    logger.info("Creating combined UpSet plot across all databases")
    try:
        combined_pathways = {}
        for gene in genes:
            all_gene_pathways = set()
            for database in ["kegg", "reactome", "msigdb"]:
                all_gene_pathways.update(gene_results[gene].get(database, set()))
            if all_gene_pathways:
                combined_pathways[gene] = all_gene_pathways

        if combined_pathways:
            upset_data = from_contents(combined_pathways)

            # Calculate min_subset_size for combined plot
            if min_subset_size is None:
                subset_counts = upset_data.value_counts()
                if len(subset_counts) > max_subsets:
                    min_subset_size_calc = sorted(subset_counts.values(), reverse=True)[max_subsets - 1]
                else:
                    min_subset_size_calc = 1
            else:
                min_subset_size_calc = min_subset_size

            fig = plt.figure(figsize=(14, 10))
            upset = UpSet(
                upset_data,
                subset_size='count',
                intersection_plot_elements=3,
                show_counts='%d',
                min_subset_size=min_subset_size_calc,
                sort_by='cardinality',
                sort_categories_by='cardinality'
            )
            upset.plot(fig=fig)

            fig.suptitle(
                f'Pathway Overlap Across Genes - All Databases Combined (intersections with ≥{min_subset_size_calc} pathways)',
                fontsize=10,
                y=0.98,
                fontweight='bold'
            )

            # Reserve 20% of the figure height free at the top for the title
            fig.tight_layout(rect=[0, 0, 1, 0.8])

            output_file = f"{output_prefix}_combined.png"
            fig.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.3)
            logger.info(f"Saved combined UpSet plot to {output_file}")
            plt.close(fig)

    except Exception as e:
        logger.error(f"Failed to create combined UpSet plot: {e}")


def export_to_table(
    gene_results: Dict[str, Dict[str, Set[str]]],
    output_prefix: str = "pathway_table"
):
    """
    Export results to CSV and Excel tables.

    Creates:
    1. Wide format table with separate columns for each database
    2. Long format table with all pathways
    """
    genes = list(gene_results.keys())

    # Wide format table
    logger.info("Creating wide format table")
    wide_data = []
    for gene in genes:
        kegg_pathways = list(gene_results[gene].get("kegg", set()))
        reactome_pathways = list(gene_results[gene].get("reactome", set()))
        msigdb_pathways = list(gene_results[gene].get("msigdb", set()))

        wide_data.append({
            "Gene": gene,
            "KEGG_Pathways": "; ".join(sorted(kegg_pathways)),
            "KEGG_Count": len(kegg_pathways),
            "Reactome_Pathways": "; ".join(sorted(reactome_pathways)),
            "Reactome_Count": len(reactome_pathways),
            "MSigDB_GeneSets": "; ".join(sorted(msigdb_pathways)),
            "MSigDB_Count": len(msigdb_pathways),
            "Total_Count": len(kegg_pathways) + len(reactome_pathways) + len(msigdb_pathways)
        })

    wide_df = pd.DataFrame(wide_data)
    wide_csv = f"{output_prefix}_wide.csv"
    wide_df.to_csv(wide_csv, index=False)
    logger.info(f"Saved wide format table to {wide_csv}")

    # Excel format with separate sheets
    try:
        excel_file = f"{output_prefix}_wide.xlsx"
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            wide_df.to_excel(writer, sheet_name='Summary', index=False)

            # Create separate sheets for each database
            for database in ["kegg", "reactome", "msigdb"]:
                db_data = []
                for gene in genes:
                    pathways = sorted(gene_results[gene].get(database, set()))
                    for pathway in pathways:
                        db_data.append({"Gene": gene, "Pathway": pathway})

                if db_data:
                    db_df = pd.DataFrame(db_data)
                    db_df.to_excel(writer, sheet_name=database.upper(), index=False)

        logger.info(f"Saved Excel table to {excel_file}")
    except Exception as e:
        logger.warning(f"Failed to create Excel file: {e}")

    # Long format table
    logger.info("Creating long format table")
    long_data = []
    for gene in genes:
        for database in ["kegg", "reactome", "msigdb"]:
            pathways = gene_results[gene].get(database, set())
            for pathway in sorted(pathways):
                long_data.append({
                    "Gene": gene,
                    "Database": database.upper(),
                    "Pathway": pathway
                })

    long_df = pd.DataFrame(long_data)
    long_csv = f"{output_prefix}_long.csv"
    long_df.to_csv(long_csv, index=False)
    logger.info(f"Saved long format table to {long_csv}")

    # Create summary statistics
    summary_data = []
    for gene in genes:
        for database in ["kegg", "reactome", "msigdb"]:
            count = len(gene_results[gene].get(database, set()))
            summary_data.append({
                "Gene": gene,
                "Database": database.upper(),
                "Pathway_Count": count
            })

    summary_df = pd.DataFrame(summary_data)
    summary_csv = f"{output_prefix}_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    logger.info(f"Saved summary statistics to {summary_csv}")

    return wide_df, long_df, summary_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Multi-gene pathway analysis with UpSet plot visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s TP53 BRCA1 EGFR
  %(prog)s TP53 BRCA1 --output results
  %(prog)s TP53 BRCA1 EGFR --msigdb-collections H C2 C6
  %(prog)s gene1 gene2 gene3 --min-intersection 3 --max-bars 15

Output Files:
  {prefix}_kegg.png          - UpSet plot for KEGG pathways
  {prefix}_reactome.png      - UpSet plot for Reactome pathways
  {prefix}_msigdb.png        - UpSet plot for MSigDB gene sets
  {prefix}_combined.png      - UpSet plot for all databases combined
  {prefix}_wide.csv          - Wide format table (one row per gene)
  {prefix}_wide.xlsx         - Excel file with separate sheets
  {prefix}_long.csv          - Long format table (one row per gene-pathway)
  {prefix}_summary.csv       - Summary statistics
        """
    )

    parser.add_argument(
        "genes",
        nargs="+",
        help="Gene symbols to query (e.g., TP53 BRCA1 EGFR)"
    )
    parser.add_argument(
        "-o", "--output",
        default="pathway_analysis",
        help="Output file prefix (default: pathway_analysis)"
    )
    parser.add_argument(
        "--kegg-organism",
        default="hsa",
        help="KEGG organism code (default: hsa)"
    )
    parser.add_argument(
        "--reactome-species",
        type=int,
        default=9606,
        help="Reactome species taxonomy ID (default: 9606)"
    )
    parser.add_argument(
        "--msigdb-collections",
        nargs="+",
        default=["H", "C2"],
        help="MSigDB collections to query (default: H C2)"
    )
    parser.add_argument(
        "--msigdb-version",
        default="2025.1.Hs",
        help="MSigDB version (default: 2025.1.Hs)"
    )
    parser.add_argument(
        "--min-intersection",
        type=int,
        default=None,
        help="Minimum intersection size to show in UpSet plot"
    )
    parser.add_argument(
        "--max-bars",
        type=int,
        default=20,
        help="Maximum number of bars in UpSet plot (default: 20)"
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip UpSet plot generation, only export tables"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Validate input
    if len(args.genes) < 2:
        logger.error("Please provide at least 2 genes for comparison")
        sys.exit(1)

    logger.info(f"Analyzing {len(args.genes)} genes: {', '.join(args.genes)}")

    # Query all genes
    try:
        gene_results = query_multiple_genes(
            genes=args.genes,
            kegg_organism=args.kegg_organism,
            reactome_species=args.reactome_species,
            msigdb_collections=args.msigdb_collections,
            msigdb_version=args.msigdb_version
        )

        # Export to tables
        logger.info("Exporting results to tables")
        wide_df, long_df, summary_df = export_to_table(gene_results, args.output)

        # Print summary
        print("\n" + "="*80)
        print("ANALYSIS SUMMARY")
        print("="*80)
        print(f"\nGenes analyzed: {', '.join(args.genes)}")
        print(f"\nPathway counts by gene:")
        print(summary_df.pivot(index='Gene', columns='Database', values='Pathway_Count'))

        # Create UpSet plots
        if not args.no_plot:
            logger.info("Creating UpSet plots")
            create_upset_plot(
                gene_results,
                output_prefix=args.output,
                min_subset_size=args.min_intersection,
                max_subsets=args.max_bars
            )
            print(f"\nUpSet plots saved with prefix: {args.output}")

        print(f"\nResults exported to:")
        print(f"  - {args.output}_wide.csv")
        print(f"  - {args.output}_wide.xlsx")
        print(f"  - {args.output}_long.csv")
        print(f"  - {args.output}_summary.csv")

        if not args.no_plot:
            print(f"  - {args.output}_kegg.png")
            print(f"  - {args.output}_reactome.png")
            print(f"  - {args.output}_msigdb.png")
            print(f"  - {args.output}_combined.png")

        print("\n" + "="*80)

    except Exception as e:
        logger.error(f"Analysis failed: {e}", exc_info=args.verbose)
        sys.exit(1)
