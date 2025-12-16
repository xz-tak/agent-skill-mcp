#!/usr/bin/env python3
"""
Unified interface to query multiple pathway databases for a given gene.

This script provides a convenient way to query KEGG, MSigDB, and Reactome
databases simultaneously to find all pathways/gene sets containing a target gene.

Usage:
    python pathway_query.py TP53
    python pathway_query.py TP53 --databases kegg msigdb
    python pathway_query.py TP53 --output json
    python pathway_query.py TP53 --all --export results.json

Requires: requests, gseapy
"""

import json
import logging
import sys
from typing import List, Dict, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

# Import individual API modules
try:
    from kegg_api import get_pathways_for_gene as query_kegg
except ImportError:
    query_kegg = None

try:
    from msigdb_api import get_gene_sets_for_gene as query_msigdb
except ImportError:
    query_msigdb = None

try:
    from reactome_api import get_reactome_terms_for_gene as query_reactome
except ImportError:
    query_reactome = None

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def query_all_databases(
    gene_symbol: str,
    databases: Optional[List[str]] = None,
    kegg_organism: str = "hsa",
    msigdb_collection: str = "H",
    msigdb_version: str = "2025.1.Hs",
    reactome_species: int = 9606,
    reactome_resource: str = "UniProt",
    parallel: bool = True
) -> Dict[str, any]:
    """
    Query multiple pathway databases for a given gene.

    Args:
        gene_symbol: Gene symbol to search for
        databases: List of databases to query (kegg, msigdb, reactome). If None, queries all.
        kegg_organism: KEGG organism code (default: 'hsa' for human)
        msigdb_collection: MSigDB collection (default: 'H' for Hallmark)
        msigdb_version: MSigDB version (default: '2025.1.Hs')
        reactome_species: Reactome species taxonomy ID (default: 9606 for human)
        reactome_resource: Reactome resource type (default: 'UniProt')
        parallel: Run queries in parallel (default: True)

    Returns:
        Dictionary with results from each database
    """
    if databases is None:
        databases = ["kegg", "msigdb", "reactome"]

    databases = [db.lower() for db in databases]
    results = {
        "gene_symbol": gene_symbol,
        "databases": {}
    }

    def query_database(db_name: str) -> tuple:
        """Query a single database and return results."""
        try:
            if db_name == "kegg" and query_kegg:
                logger.info(f"Querying KEGG for {gene_symbol}")
                result = query_kegg(gene_symbol, kegg_organism)
                return db_name, result, None
            elif db_name == "msigdb" and query_msigdb:
                logger.info(f"Querying MSigDB for {gene_symbol}")
                result = query_msigdb(gene_symbol, msigdb_collection, msigdb_version)
                return db_name, result, None
            elif db_name == "reactome" and query_reactome:
                logger.info(f"Querying Reactome for {gene_symbol}")
                result = query_reactome(gene_symbol, reactome_resource, reactome_species)
                return db_name, result, None
            else:
                error_msg = f"{db_name} API not available or not supported"
                logger.warning(error_msg)
                return db_name, None, error_msg
        except Exception as e:
            error_msg = f"Error querying {db_name}: {e}"
            logger.error(error_msg)
            return db_name, None, error_msg

    # Query databases in parallel or sequentially
    if parallel:
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = {executor.submit(query_database, db): db for db in databases}
            for future in as_completed(futures):
                db_name, result, error = future.result()
                if error:
                    results["databases"][db_name] = {"error": error}
                else:
                    results["databases"][db_name] = result
    else:
        for db in databases:
            db_name, result, error = query_database(db)
            if error:
                results["databases"][db_name] = {"error": error}
            else:
                results["databases"][db_name] = result

    return results


def summarize_results(results: Dict[str, any]) -> Dict[str, any]:
    """Generate a summary of results across all databases."""
    summary = {
        "gene_symbol": results["gene_symbol"],
        "total_databases_queried": len(results["databases"]),
        "databases_with_results": 0,
        "total_pathways": 0,
        "breakdown": {}
    }

    for db_name, db_result in results["databases"].items():
        if "error" in db_result:
            summary["breakdown"][db_name] = {"status": "error", "count": 0}
            continue

        # Extract counts based on database structure
        count = 0
        if db_name == "kegg":
            count = db_result.get("total_pathways", 0)
        elif db_name == "msigdb":
            count = db_result.get("total_gene_sets", 0)
        elif db_name == "reactome":
            count = db_result.get("total_terms", 0)

        summary["breakdown"][db_name] = {"status": "success", "count": count}
        summary["total_pathways"] += count
        if count > 0:
            summary["databases_with_results"] += 1

    return summary


def export_to_tables(results: Dict[str, any], output_prefix: str) -> None:
    """Export results to CSV and Excel tables."""
    gene_symbol = results["gene_symbol"]

    # Create unified long format table (gene-database-pathway)
    long_data = []

    for db_name, db_result in results["databases"].items():
        if "error" in db_result:
            continue

        # Extract pathways/terms from each database
        if db_name == "kegg":
            for gene_data in db_result.get("genes_found", []):
                for pathway in gene_data.get("pathways", []):
                    long_data.append({
                        "Gene": gene_symbol,
                        "Database": "KEGG",
                        "Database_ID": pathway.get("pathway_id", ""),
                        "Pathway_Name": pathway.get("pathway_name", ""),
                        "Gene_ID": gene_data.get("kegg_gene_id", "")
                    })

        elif db_name == "msigdb":
            for gene_set in db_result.get("gene_sets", []):
                long_data.append({
                    "Gene": gene_symbol,
                    "Database": "MSIGDB",
                    "Database_ID": gene_set.get("gene_set_name", ""),
                    "Pathway_Name": gene_set.get("gene_set_name", ""),
                    "Gene_ID": gene_symbol
                })

        elif db_name == "reactome":
            for term in db_result.get("terms", []):
                long_data.append({
                    "Gene": gene_symbol,
                    "Database": "REACTOME",
                    "Database_ID": term.get("stId", ""),
                    "Pathway_Name": term.get("name", ""),
                    "Gene_ID": gene_symbol
                })

    if long_data:
        long_df = pd.DataFrame(long_data)

        # Save as CSV
        csv_file = f"{output_prefix}_pathways.csv"
        long_df.to_csv(csv_file, index=False)
        logger.info(f"Saved pathways to {csv_file}")

        # Save as Excel with separate sheets per database
        try:
            excel_file = f"{output_prefix}_pathways.xlsx"
            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                # All pathways sheet
                long_df.to_excel(writer, sheet_name='All_Pathways', index=False)

                # Database-specific sheets
                for db in ["KEGG", "MSIGDB", "REACTOME"]:
                    db_df = long_df[long_df["Database"] == db]
                    if not db_df.empty:
                        db_df.to_excel(writer, sheet_name=db, index=False)

                # Summary sheet
                summary = summarize_results(results)
                summary_data = []
                for db_name, info in summary["breakdown"].items():
                    summary_data.append({
                        "Database": db_name.upper(),
                        "Status": info["status"],
                        "Pathway_Count": info["count"]
                    })
                summary_df = pd.DataFrame(summary_data)
                summary_df.to_excel(writer, sheet_name='Summary', index=False)

            logger.info(f"Saved pathways to {excel_file}")
        except Exception as e:
            logger.warning(f"Failed to create Excel file: {e}")
    else:
        logger.warning("No pathways found across any database, skipping table export")


def format_output(results: Dict[str, any], output_format: str = "text", summary_only: bool = False) -> str:
    """Format output in the specified format."""
    if output_format == "json":
        return json.dumps(results, indent=2)

    # Text format
    output = []
    output.append(f"=" * 80)
    output.append(f"Pathway Database Query Results for: {results['gene_symbol']}")
    output.append(f"=" * 80)

    # Generate summary
    summary = summarize_results(results)
    output.append(f"\nSummary:")
    output.append(f"  Total databases queried: {summary['total_databases_queried']}")
    output.append(f"  Databases with results: {summary['databases_with_results']}")
    output.append(f"  Total pathways/gene sets found: {summary['total_pathways']}")

    for db_name, info in summary["breakdown"].items():
        status_icon = "✓" if info["status"] == "success" else "✗"
        output.append(f"    {status_icon} {db_name.upper()}: {info['count']} results")

    if summary_only:
        return "\n".join(output)

    # Detailed results for each database
    for db_name, db_result in results["databases"].items():
        output.append(f"\n{'-' * 80}")
        output.append(f"{db_name.upper()} Results")
        output.append(f"{'-' * 80}")

        if "error" in db_result:
            output.append(f"Error: {db_result['error']}")
            continue

        # Format based on database
        if db_name == "kegg":
            output.append(f"Organism: {db_result['query']['organism']}")
            output.append(f"Total pathways: {db_result['total_pathways']}\n")
            for gene_data in db_result["genes_found"]:
                output.append(f"Gene: {gene_data['kegg_gene_id']}")
                for pathway in gene_data["pathways"][:10]:  # Show first 10
                    output.append(f"  {pathway['pathway_id']}: {pathway['pathway_name']}")
                if len(gene_data["pathways"]) > 10:
                    output.append(f"  ... and {len(gene_data['pathways']) - 10} more")

        elif db_name == "msigdb":
            output.append(f"Collection: {db_result['query']['collection']}")
            output.append(f"Version: {db_result['query']['dbver']}")
            output.append(f"Total gene sets: {db_result['total_gene_sets']}\n")
            for gs in db_result["gene_sets"][:10]:  # Show first 10
                output.append(f"  {gs['gene_set_name']} (n={gs['gene_count']})")
            if len(db_result["gene_sets"]) > 10:
                output.append(f"  ... and {len(db_result['gene_sets']) - 10} more")

        elif db_name == "reactome":
            output.append(f"Species: {db_result['query']['species']}")
            output.append(f"Resource: {db_result['query']['resource']}")
            output.append(f"Total terms: {db_result['total_terms']}\n")
            for term in db_result["terms"][:10]:  # Show first 10
                output.append(f"  {term['stId']}: {term['name']}")
            if len(db_result["terms"]) > 10:
                output.append(f"  ... and {len(db_result['terms']) - 10} more")

    output.append(f"\n{'=' * 80}")
    return "\n".join(output)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Query multiple pathway databases for a given gene",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s TP53
  %(prog)s TP53 --databases kegg msigdb
  %(prog)s TP53 --all --output json
  %(prog)s TP53 --summary-only
  %(prog)s TP53 --export results.json

Databases:
  kegg     - KEGG pathway database
  msigdb   - Molecular Signatures Database (MSigDB)
  reactome - Reactome pathway database
        """
    )

    parser.add_argument("gene", help="Gene symbol (e.g., TP53)")
    parser.add_argument(
        "-d", "--databases",
        nargs="+",
        choices=["kegg", "msigdb", "reactome"],
        help="Databases to query (default: all)"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Query all databases (default behavior)"
    )
    parser.add_argument(
        "-o", "--output",
        choices=["text", "json", "csv"],
        default="csv",
        help="Output format (default: csv)"
    )
    parser.add_argument(
        "--output-prefix", "-p",
        default=None,
        help="Output file prefix for CSV/Excel (default: pathway_query_{gene})"
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Show only summary, not detailed results"
    )
    parser.add_argument(
        "--export",
        metavar="FILE",
        help="Export full results to JSON file"
    )
    parser.add_argument(
        "--kegg-organism",
        default="hsa",
        help="KEGG organism code (default: hsa)"
    )
    parser.add_argument(
        "--msigdb-collection",
        default="H",
        help="MSigDB collection (default: H)"
    )
    parser.add_argument(
        "--msigdb-version",
        default="2025.1.Hs",
        help="MSigDB version (default: 2025.1.Hs)"
    )
    parser.add_argument(
        "--reactome-species",
        type=int,
        default=9606,
        help="Reactome species taxonomy ID (default: 9606)"
    )
    parser.add_argument(
        "--no-parallel",
        action="store_true",
        help="Run queries sequentially instead of in parallel"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Determine which databases to query
    databases = args.databases if args.databases else None

    try:
        results = query_all_databases(
            gene_symbol=args.gene,
            databases=databases,
            kegg_organism=args.kegg_organism,
            msigdb_collection=args.msigdb_collection,
            msigdb_version=args.msigdb_version,
            reactome_species=args.reactome_species,
            parallel=not args.no_parallel
        )

        # Export to CSV/Excel by default
        if args.output in ["csv", "text"]:
            output_prefix = args.output_prefix or f"pathway_query_{args.gene}"
            export_to_tables(results, output_prefix)

            # Print summary to console
            summary = summarize_results(results)
            print(f"\n{'='*80}")
            print(f"Pathway Database Query: {args.gene}")
            print(f"{'='*80}")
            print(f"Total databases queried: {summary['total_databases_queried']}")
            print(f"Databases with results: {summary['databases_with_results']}")
            print(f"Total pathways found: {summary['total_pathways']}")
            print(f"\nBreakdown by database:")
            for db_name, info in summary["breakdown"].items():
                status_icon = "✓" if info["status"] == "success" else "✗"
                print(f"  {status_icon} {db_name.upper()}: {info['count']} results")
            print(f"\nResults exported to:")
            print(f"  - {output_prefix}_pathways.csv")
            print(f"  - {output_prefix}_pathways.xlsx")
            print(f"{'='*80}\n")

        # Export to JSON file if requested
        if args.export:
            with open(args.export, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results exported to {args.export}")

        # JSON output format
        if args.output == "json":
            output = format_output(results, args.output, args.summary_only)
            print(output)

    except Exception as e:
        logger.error(f"Error: {e}", exc_info=args.verbose)
        sys.exit(1)
