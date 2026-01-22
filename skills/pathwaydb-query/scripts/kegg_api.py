#!/usr/bin/env python3
"""
Query KEGG database for pathways containing a given gene.

Given a gene symbol (e.g. TP53) and KEGG organism code (e.g. hsa),
retrieve all KEGG pathway IDs and names that contain this gene.

Usage:
    python kegg_api.py TP53 hsa
    python kegg_api.py --gene TP53 --organism hsa --output json
"""

import sys
import json
import logging
import time
from functools import lru_cache
from typing import List, Dict, Tuple, Optional
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

KEGG_REST = "https://rest.kegg.jp"
REQUEST_TIMEOUT = 30
MAX_RETRIES = 3
BACKOFF_FACTOR = 0.5


def get_session_with_retry() -> requests.Session:
    """Create a requests session with retry logic."""
    session = requests.Session()
    retry_strategy = Retry(
        total=MAX_RETRIES,
        backoff_factor=BACKOFF_FACTOR,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


def find_kegg_gene_ids(symbol: str, org: str, session: Optional[requests.Session] = None) -> List[Tuple[str, str]]:
    """
    Use KEGG 'find' to get KEGG gene IDs for a gene symbol in an organism.

    Args:
        symbol: Gene symbol to search for
        org: KEGG organism code (e.g., 'hsa' for human)
        session: Optional requests session with retry logic

    Returns:
        List of tuples (kegg_gene_id, description)

    Example:
        >>> find_kegg_gene_ids("TP53", "hsa")
        [("hsa:7157", "TP53, tumor protein p53")]
    """
    if session is None:
        session = get_session_with_retry()

    url = f"{KEGG_REST}/find/{org}/{symbol}"
    logger.debug(f"Querying KEGG find: {url}")

    try:
        resp = session.get(url, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to query KEGG find for {symbol} in {org}: {e}")
        return []

    hits = []
    for line in resp.text.strip().splitlines():
        if not line.strip():
            continue
        try:
            entry_id, desc = line.split("\t", 1)
            hits.append((entry_id, desc))
        except ValueError:
            logger.warning(f"Failed to parse line: {line}")
            continue

    logger.info(f"Found {len(hits)} KEGG gene IDs for {symbol}")
    return hits


@lru_cache(maxsize=128)
def load_pathway_names(org: str, session: Optional[requests.Session] = None) -> Dict[str, str]:
    """
    Load all organism-specific pathway IDs and names as a dict.
    Results are cached to avoid repeated API calls.

    Args:
        org: KEGG organism code (e.g., 'hsa' for human)
        session: Optional requests session with retry logic

    Returns:
        Dictionary mapping pathway IDs to pathway names

    Example:
        >>> load_pathway_names("hsa")
        {"path:hsa04115": "p53 signaling pathway - Homo sapiens (human)", ...}
    """
    if session is None:
        session = get_session_with_retry()

    url = f"{KEGG_REST}/list/pathway/{org}"
    logger.debug(f"Loading pathway names: {url}")

    try:
        resp = session.get(url, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to load pathway names for {org}: {e}")
        return {}

    mapping = {}
    for line in resp.text.strip().splitlines():
        if not line.strip():
            continue
        try:
            pid, name = line.split("\t", 1)
            mapping[pid] = name
        except ValueError:
            logger.warning(f"Failed to parse line: {line}")
            continue

    logger.info(f"Loaded {len(mapping)} pathway names for {org}")
    return mapping


def pathways_for_gene(kegg_gene_id: str, session: Optional[requests.Session] = None) -> List[str]:
    """
    Get all KEGG pathways linked to a given KEGG gene ID.

    Args:
        kegg_gene_id: KEGG gene identifier (e.g., 'hsa:7157')
        session: Optional requests session with retry logic

    Returns:
        List of KEGG pathway IDs

    Example:
        >>> pathways_for_gene("hsa:7157")
        ["path:hsa04115", "path:hsa04151", ...]
    """
    if session is None:
        session = get_session_with_retry()

    url = f"{KEGG_REST}/link/pathway/{kegg_gene_id}"
    logger.debug(f"Querying pathways for gene: {url}")

    try:
        resp = session.get(url, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to get pathways for {kegg_gene_id}: {e}")
        return []

    path_ids = []
    for line in resp.text.strip().splitlines():
        if not line.strip():
            continue
        try:
            gene_id, path_id = line.split("\t", 1)
            path_ids.append(path_id)
        except ValueError:
            logger.warning(f"Failed to parse line: {line}")
            continue

    logger.info(f"Found {len(path_ids)} pathways for {kegg_gene_id}")
    return path_ids


def get_pathways_for_gene(symbol: str, org: str = "hsa") -> Dict[str, any]:
    """
    Main function to get all pathways containing a given gene.

    Args:
        symbol: Gene symbol to search for
        org: KEGG organism code (default: 'hsa' for human)

    Returns:
        Dictionary with query info and pathway results
    """
    session = get_session_with_retry()

    # Step 1: find KEGG gene IDs for the symbol
    gene_hits = find_kegg_gene_ids(symbol, org, session=session)
    if not gene_hits:
        logger.warning(f"No KEGG genes found for {symbol!r} in organism {org!r}")
        return {
            "query": {"gene_symbol": symbol, "organism": org},
            "genes_found": [],
            "pathways": []
        }

    # Step 2: load all pathway names once
    pathway_names = load_pathway_names(org, session=session)

    # Step 3: for each gene hit, collect its pathways
    results = []
    for kegg_gene_id, desc in gene_hits:
        path_ids = pathways_for_gene(kegg_gene_id, session=session)

        pathways = []
        for path_id in sorted(set(path_ids)):
            # Normalize pathway ID by removing 'path:' prefix for lookup
            lookup_id = path_id.replace("path:", "") if path_id.startswith("path:") else path_id
            pathway_name = pathway_names.get(lookup_id, pathway_names.get(path_id, ""))

            pathways.append({
                "pathway_id": path_id,
                "pathway_name": pathway_name
            })

        results.append({
            "kegg_gene_id": kegg_gene_id,
            "description": desc,
            "pathways": pathways
        })

    return {
        "query": {"gene_symbol": symbol, "organism": org},
        "genes_found": results,
        "total_genes": len(results),
        "total_pathways": sum(len(r["pathways"]) for r in results)
    }


def export_to_tables(data: Dict[str, any], output_prefix: str) -> None:
    """Export results to CSV and Excel tables."""
    # Create long format table (gene-pathway pairs)
    long_data = []
    for gene_data in data["genes_found"]:
        gene_id = gene_data["kegg_gene_id"]
        gene_desc = gene_data["description"]
        for pathway in gene_data["pathways"]:
            long_data.append({
                "KEGG_Gene_ID": gene_id,
                "Gene_Description": gene_desc,
                "Pathway_ID": pathway["pathway_id"],
                "Pathway_Name": pathway["pathway_name"]
            })

    if long_data:
        long_df = pd.DataFrame(long_data)

        # Save as CSV
        csv_file = f"{output_prefix}_pathways.csv"
        long_df.to_csv(csv_file, index=False)
        logger.info(f"Saved pathways to {csv_file}")

        # Save as Excel
        try:
            excel_file = f"{output_prefix}_pathways.xlsx"
            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                long_df.to_excel(writer, sheet_name='Pathways', index=False)

                # Add summary sheet
                summary_data = []
                for gene_data in data["genes_found"]:
                    summary_data.append({
                        "KEGG_Gene_ID": gene_data["kegg_gene_id"],
                        "Gene_Description": gene_data["description"],
                        "Pathway_Count": len(gene_data["pathways"])
                    })
                summary_df = pd.DataFrame(summary_data)
                summary_df.to_excel(writer, sheet_name='Summary', index=False)

            logger.info(f"Saved pathways to {excel_file}")
        except Exception as e:
            logger.warning(f"Failed to create Excel file: {e}")
    else:
        logger.warning("No pathways found, skipping table export")


def format_output(data: Dict[str, any], output_format: str = "text") -> str:
    """Format output in the specified format."""
    if output_format == "json":
        return json.dumps(data, indent=2)

    # Text format (for display only)
    output = []
    query = data["query"]
    output.append(f"Query: {query['gene_symbol']} in organism {query['organism']}")
    output.append(f"Found {data['total_genes']} gene(s) with {data['total_pathways']} pathway(s)\n")

    for gene_data in data["genes_found"]:
        output.append(f"Gene: {gene_data['kegg_gene_id']}\t{gene_data['description']}")
        if not gene_data["pathways"]:
            output.append("  (no pathways found)")
        else:
            for pathway in gene_data["pathways"]:
                output.append(f"  {pathway['pathway_id']}\t{pathway['pathway_name']}")
        output.append("")

    return "\n".join(output)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Query KEGG database for pathways containing a given gene",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s TP53
  %(prog)s TP53 hsa
  %(prog)s --gene TP53 --organism hsa --output json
  %(prog)s TP53 --verbose
        """
    )
    parser.add_argument("gene", nargs="?", help="Gene symbol (e.g., TP53)")
    parser.add_argument("organism", nargs="?", default="hsa",
                        help="KEGG organism code (default: hsa)")
    parser.add_argument("--gene", "-g", dest="gene_flag",
                        help="Gene symbol (alternative to positional arg)")
    parser.add_argument("--organism", "-o", dest="organism_flag",
                        help="KEGG organism code (alternative to positional arg)")
    parser.add_argument("--output", "-f", choices=["text", "json", "csv"],
                        default="csv", help="Output format (default: csv)")
    parser.add_argument("--output-prefix", "-p", default=None,
                        help="Output file prefix for CSV/Excel (default: kegg_{gene})")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable verbose logging")

    args = parser.parse_args()

    # Handle both positional and flag-based arguments
    gene_symbol = args.gene_flag or args.gene
    organism = args.organism_flag or args.organism

    if not gene_symbol:
        parser.print_help()
        sys.exit(1)

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    try:
        result = get_pathways_for_gene(gene_symbol, organism)

        # Export to CSV/Excel by default
        if args.output in ["csv", "text"]:
            output_prefix = args.output_prefix or f"kegg_{gene_symbol}_{organism}"
            export_to_tables(result, output_prefix)

            # Print summary to console
            print(f"\n{'='*80}")
            print(f"KEGG Query: {gene_symbol} in organism {organism}")
            print(f"{'='*80}")
            print(f"Total genes found: {result['total_genes']}")
            print(f"Total pathways: {result['total_pathways']}")
            print(f"\nResults exported to:")
            print(f"  - {output_prefix}_pathways.csv")
            print(f"  - {output_prefix}_pathways.xlsx")
            print(f"{'='*80}\n")

        # JSON output
        if args.output == "json":
            output = format_output(result, args.output)
            print(output)

    except Exception as e:
        logger.error(f"Error querying KEGG: {e}", exc_info=args.verbose)
        sys.exit(1)
