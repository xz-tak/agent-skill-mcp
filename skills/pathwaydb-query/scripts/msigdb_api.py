#!/usr/bin/env python3
"""
Query MSigDB for gene sets/signatures containing a given gene.

For a given gene and collection (H or C1–C8), retrieve all gene sets
that contain that gene as a member.

Usage:
    python msigdb_api.py TP53
    python msigdb_api.py TP53 --collection H --output json
    python msigdb_api.py TP53 -c C2 -v 2025.1.Hs

Requires: gseapy >= 1.1.8
"""

import json
import logging
import sys
from functools import lru_cache
from typing import List, Dict, Optional
import pandas as pd
from gseapy import Msigdb

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Mapping from MSigDB collection code to Msigdb "category" name
COLLECTION_TO_CATEGORY = {
    "H": "h.all",   # Hallmark
    "C1": "c1.all",  # Positional gene sets
    "C2": "c2.all",  # Curated gene sets
    "C3": "c3.all",  # Regulatory target gene sets
    "C4": "c4.all",  # Computational gene sets
    "C5": "c5.all",  # Ontology gene sets
    "C6": "c6.all",  # Oncogenic signature gene sets
    "C7": "c7.all",  # Immunologic signature gene sets
    "C8": "c8.all",  # Cell type signature gene sets
}


@lru_cache(maxsize=32)
def _load_gmt_cached(category: str, dbver: str) -> Dict[str, List[str]]:
    """
    Load and cache GMT files from MSigDB to avoid repeated downloads.

    Args:
        category: MSigDB category (e.g., 'h.all', 'c2.all')
        dbver: MSigDB version (e.g., '2025.1.Hs')

    Returns:
        Dictionary mapping gene set names to gene lists
    """
    logger.info(f"Loading MSigDB GMT: category={category}, version={dbver}")
    try:
        msig = Msigdb()
        gmt = msig.get_gmt(category=category, dbver=dbver)
        logger.info(f"Loaded {len(gmt)} gene sets from {category}")
        return gmt
    except Exception as e:
        logger.error(f"Failed to load GMT for {category} ({dbver}): {e}")
        raise


def get_gene_sets_for_gene(
    gene_symbol: str,
    collection: str = "H",
    dbver: str = "2025.1.Hs"
) -> Dict[str, any]:
    """
    Query MSigDB for gene sets containing a given gene.

    Args:
        gene_symbol: Gene symbol to search for (case-insensitive)
        collection: MSigDB collection code (H, C1-C8)
        dbver: MSigDB database version (e.g., '2025.1.Hs')

    Returns:
        Dictionary with query info and matching gene sets

    Raises:
        ValueError: If collection is not supported

    Example:
        >>> result = get_gene_sets_for_gene("TP53", "H")
        >>> print(result["total_gene_sets"])
        15
    """
    gene_symbol = gene_symbol.strip().upper()
    coll = collection.upper()

    if coll not in COLLECTION_TO_CATEGORY:
        raise ValueError(
            f"Unsupported collection '{collection}'. "
            f"Use one of: {', '.join(COLLECTION_TO_CATEGORY.keys())}"
        )

    category = COLLECTION_TO_CATEGORY[coll]
    logger.info(f"Querying MSigDB for gene={gene_symbol}, collection={coll}")

    # Load GMT with caching
    try:
        gmt = _load_gmt_cached(category, dbver)
    except Exception as e:
        logger.error(f"Failed to query MSigDB: {e}")
        return {
            "query": {
                "gene_symbol": gene_symbol,
                "collection": coll,
                "dbver": dbver
            },
            "gene_sets": [],
            "total_gene_sets": 0
        }

    # Search for gene sets containing the query gene
    hits = []
    for term, genes in gmt.items():
        genes_upper = {g.upper() for g in genes}
        if gene_symbol in genes_upper:
            hits.append({
                "gene_set_name": term,
                "gene_count": len(genes)
            })

    logger.info(f"Found {len(hits)} gene sets containing {gene_symbol}")

    return {
        "query": {
            "gene_symbol": gene_symbol,
            "collection": coll,
            "dbver": dbver
        },
        "gene_sets": sorted(hits, key=lambda x: x["gene_set_name"]),
        "total_gene_sets": len(hits)
    }


def get_gene_sets_across_collections(
    gene_symbol: str,
    collections: List[str] = None,
    dbver: str = "2025.1.Hs"
) -> Dict[str, any]:
    """
    Query MSigDB for gene sets containing a given gene across multiple collections.

    Args:
        gene_symbol: Gene symbol to search for (case-insensitive)
        collections: List of MSigDB collection codes (H, C1-C8). If None, queries all.
        dbver: MSigDB database version (e.g., '2025.1.Hs')

    Returns:
        Dictionary with results from all collections

    Example:
        >>> result = get_gene_sets_across_collections("TP53", ["H", "C2"])
        >>> print(result["total_gene_sets_all_collections"])
        150
    """
    if collections is None:
        collections = list(COLLECTION_TO_CATEGORY.keys())

    gene_symbol = gene_symbol.strip().upper()
    results = {
        "gene_symbol": gene_symbol,
        "dbver": dbver,
        "collections": {},
        "total_gene_sets_all_collections": 0,
        "total_collections_queried": len(collections)
    }

    for coll in collections:
        logger.info(f"Querying collection {coll} for {gene_symbol}")
        try:
            result = get_gene_sets_for_gene(gene_symbol, coll, dbver)
            results["collections"][coll] = result
            results["total_gene_sets_all_collections"] += result["total_gene_sets"]
        except Exception as e:
            logger.error(f"Failed to query collection {coll}: {e}")
            results["collections"][coll] = {
                "error": str(e),
                "total_gene_sets": 0
            }

    return results


def export_to_tables(data: Dict[str, any], output_prefix: str) -> None:
    """Export results to CSV and Excel tables."""
    # Check if single or multi-collection
    if "collections" in data:
        # Multi-collection format
        long_data = []
        for coll_name, coll_data in data["collections"].items():
            if "error" not in coll_data:
                for gs in coll_data.get("gene_sets", []):
                    long_data.append({
                        "Gene_Symbol": data["gene_symbol"],
                        "Collection": coll_name,
                        "Gene_Set_Name": gs["gene_set_name"],
                        "Gene_Count": gs["gene_count"]
                    })

        if long_data:
            long_df = pd.DataFrame(long_data)

            # Save CSV
            csv_file = f"{output_prefix}_genesets.csv"
            long_df.to_csv(csv_file, index=False)
            logger.info(f"Saved gene sets to {csv_file}")

            # Save Excel
            try:
                excel_file = f"{output_prefix}_genesets.xlsx"
                with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                    long_df.to_excel(writer, sheet_name='Gene_Sets', index=False)

                    # Add summary sheet
                    summary = long_df.groupby('Collection').size().reset_index(name='Gene_Set_Count')
                    summary.to_excel(writer, sheet_name='Summary', index=False)

                logger.info(f"Saved gene sets to {excel_file}")
            except Exception as e:
                logger.warning(f"Failed to create Excel file: {e}")
        else:
            logger.warning("No gene sets found, skipping table export")
    else:
        # Single collection format
        long_data = []
        query = data["query"]
        for gs in data.get("gene_sets", []):
            long_data.append({
                "Gene_Symbol": query["gene_symbol"],
                "Collection": query["collection"],
                "Gene_Set_Name": gs["gene_set_name"],
                "Gene_Count": gs["gene_count"]
            })

        if long_data:
            long_df = pd.DataFrame(long_data)

            # Save CSV
            csv_file = f"{output_prefix}_genesets.csv"
            long_df.to_csv(csv_file, index=False)
            logger.info(f"Saved gene sets to {csv_file}")

            # Save Excel
            try:
                excel_file = f"{output_prefix}_genesets.xlsx"
                long_df.to_excel(excel_file, sheet_name='Gene_Sets', index=False, engine='openpyxl')
                logger.info(f"Saved gene sets to {excel_file}")
            except Exception as e:
                logger.warning(f"Failed to create Excel file: {e}")
        else:
            logger.warning("No gene sets found, skipping table export")


def format_output(data: Dict[str, any], output_format: str = "text") -> str:
    """Format output in the specified format."""
    if output_format == "json":
        return json.dumps(data, indent=2)

    # Text format
    output = []

    # Check if this is multi-collection or single-collection result
    if "collections" in data:
        # Multi-collection format
        output.append(f"=" * 80)
        output.append(f"MSigDB Query Results for: {data['gene_symbol']}")
        output.append(f"=" * 80)
        output.append(f"Version: {data['dbver']}")
        output.append(f"Total collections queried: {data['total_collections_queried']}")
        output.append(f"Total gene sets found: {data['total_gene_sets_all_collections']}\n")

        for coll_name, coll_data in data["collections"].items():
            output.append(f"{'-' * 80}")
            output.append(f"Collection {coll_name}: {COLLECTION_TO_CATEGORY.get(coll_name, coll_name)}")
            output.append(f"{'-' * 80}")

            if "error" in coll_data:
                output.append(f"Error: {coll_data['error']}\n")
                continue

            output.append(f"Found {coll_data['total_gene_sets']} gene set(s)")

            if coll_data["gene_sets"]:
                # Show first 20 gene sets per collection
                for gs in coll_data["gene_sets"][:20]:
                    output.append(f"  {gs['gene_set_name']} (n={gs['gene_count']} genes)")
                if len(coll_data["gene_sets"]) > 20:
                    output.append(f"  ... and {len(coll_data['gene_sets']) - 20} more")
            else:
                output.append("  (no gene sets found)")
            output.append("")

        output.append(f"{'=' * 80}")
    else:
        # Single collection format (backward compatible)
        query = data["query"]
        output.append(f"Query: {query['gene_symbol']} in collection {query['collection']} (MSigDB {query['dbver']})")
        output.append(f"Found {data['total_gene_sets']} gene set(s)\n")

        if data["gene_sets"]:
            for gs in data["gene_sets"]:
                output.append(f"{gs['gene_set_name']} (n={gs['gene_count']} genes)")
        else:
            output.append("(no gene sets found)")

    return "\n".join(output)


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Query MSigDB for gene sets containing a given gene",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s TP53
  %(prog)s TP53 --collection C2
  %(prog)s TP53 --collections H C2 C6
  %(prog)s TP53 --all
  %(prog)s TP53 -c H -v 2025.1.Hs --output json
  %(prog)s TP53 --verbose

Collections:
  H  - Hallmark gene sets
  C1 - Positional gene sets
  C2 - Curated gene sets (pathway databases)
  C3 - Regulatory target gene sets
  C4 - Computational gene sets
  C5 - Ontology gene sets (GO terms)
  C6 - Oncogenic signature gene sets
  C7 - Immunologic signature gene sets
  C8 - Cell type signature gene sets
        """
    )
    parser.add_argument("gene", help="Gene symbol (e.g., TP53)")
    parser.add_argument(
        "-c", "--collection",
        default=None,
        help="Single MSigDB collection: H, C1-C8 (default: H if no --collections/--all)"
    )
    parser.add_argument(
        "--collections",
        nargs="+",
        choices=list(COLLECTION_TO_CATEGORY.keys()),
        help="Multiple MSigDB collections to query (e.g., H C2 C6)"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Query all MSigDB collections (H, C1-C8)"
    )
    parser.add_argument(
        "-v", "--version",
        default="2025.1.Hs",
        help="MSigDB version (default: 2025.1.Hs)"
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
        help="Output file prefix for CSV/Excel (default: msigdb_{gene})"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    try:
        # Determine which collections to query
        if args.all:
            # Query all collections
            result = get_gene_sets_across_collections(
                gene_symbol=args.gene,
                collections=None,
                dbver=args.version
            )
        elif args.collections:
            # Query multiple specified collections
            result = get_gene_sets_across_collections(
                gene_symbol=args.gene,
                collections=args.collections,
                dbver=args.version
            )
        else:
            # Query single collection (default or specified)
            collection = args.collection if args.collection else "H"
            result = get_gene_sets_for_gene(
                gene_symbol=args.gene,
                collection=collection,
                dbver=args.version,
            )

        # Export to CSV/Excel by default
        if args.output in ["csv", "text"]:
            output_prefix = args.output_prefix or f"msigdb_{args.gene}"
            export_to_tables(result, output_prefix)

            # Print summary
            print(f"\n{'='*80}")
            print(f"MSigDB Query: {args.gene}")
            print(f"{'='*80}")
            if "collections" in result:
                print(f"Collections queried: {result['total_collections_queried']}")
                print(f"Total gene sets: {result['total_gene_sets_all_collections']}")
            else:
                print(f"Collection: {result['query']['collection']}")
                print(f"Total gene sets: {result['total_gene_sets']}")
            print(f"\nResults exported to:")
            print(f"  - {output_prefix}_genesets.csv")
            print(f"  - {output_prefix}_genesets.xlsx")
            print(f"{'='*80}\n")

        # JSON output
        if args.output == "json":
            output = format_output(result, args.output)
            print(output)

    except ValueError as e:
        logger.error(f"Invalid input: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error querying MSigDB: {e}", exc_info=args.verbose)
        sys.exit(1)
