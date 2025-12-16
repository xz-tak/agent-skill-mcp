#!/usr/bin/env python3
"""
Query Reactome database for pathways/reactions containing a given gene.

Given a gene identifier and resource type, retrieve all Reactome pathways
or reactions that contain this gene.

Usage:
    python reactome_api.py TP53
    python reactome_api.py TP53 --resource UniProt --species 9606
    python reactome_api.py TP53 --map-to reactions --output json

Requires: requests
"""

import json
import logging
import sys
from typing import List, Dict, Union, Optional
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

BASE_URL = "https://reactome.org/ContentService"
REQUEST_TIMEOUT = 30
MAX_RETRIES = 3
BACKOFF_FACTOR = 0.5

# Common species mapping
SPECIES_MAP = {
    "human": 9606,
    "mouse": 10090,
    "rat": 10116,
    "zebrafish": 7955,
    "fly": 7227,
    "worm": 6239,
}


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


def get_reactome_terms_for_gene(
    gene_id: str,
    resource: str = "UniProt",
    species: Union[str, int] = 9606,
    map_to: str = "pathways",
    session: Optional[requests.Session] = None
) -> Dict[str, any]:
    """
    Query Reactome for all pathways/reactions that contain a given gene.

    Args:
        gene_id: Gene identifier (gene symbol or accession)
        resource: Identifier type ('UniProt', 'NCBI', 'ENSEMBL', etc.)
        species: Species name or NCBI taxonomy ID (default: 9606 for human)
        map_to: Query type - 'pathways' or 'reactions'
        session: Optional requests session with retry logic

    Returns:
        Dictionary with query info and pathway/reaction results

    Raises:
        ValueError: If map_to is not 'pathways' or 'reactions'

    Example:
        >>> result = get_reactome_terms_for_gene("TP53", species=9606)
        >>> print(result["total_terms"])
        25
    """
    if map_to not in {"pathways", "reactions"}:
        raise ValueError("map_to must be 'pathways' or 'reactions'")

    # Convert species name to taxonomy ID if needed
    if isinstance(species, str):
        species = SPECIES_MAP.get(species.lower(), species)

    if session is None:
        session = get_session_with_retry()

    url = f"{BASE_URL}/data/mapping/{resource}/{gene_id}/{map_to}"
    params = {"species": species}

    logger.info(f"Querying Reactome: gene={gene_id}, resource={resource}, species={species}, map_to={map_to}")

    try:
        r = session.get(url, params=params, timeout=REQUEST_TIMEOUT)
        r.raise_for_status()
        data = r.json()
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to query Reactome for {gene_id}: {e}")
        return {
            "query": {
                "gene_id": gene_id,
                "resource": resource,
                "species": species,
                "map_to": map_to
            },
            "terms": [],
            "total_terms": 0
        }

    # Parse and structure the results
    terms = []
    for item in data:
        terms.append({
            "stId": item.get("stId"),
            "dbId": item.get("dbId"),
            "name": item.get("displayName"),
            "species": item.get("speciesName"),
            "type": item.get("type", ""),
        })

    logger.info(f"Found {len(terms)} {map_to} for {gene_id}")

    return {
        "query": {
            "gene_id": gene_id,
            "resource": resource,
            "species": species,
            "map_to": map_to
        },
        "terms": terms,
        "total_terms": len(terms)
    }


def export_to_tables(data: Dict[str, any], output_prefix: str) -> None:
    """Export results to CSV and Excel tables."""
    query = data["query"]

    # Create long format table
    long_data = []
    for term in data["terms"]:
        long_data.append({
            "Gene_ID": query["gene_id"],
            "Resource": query["resource"],
            "Species": query["species"],
            "Map_To": query["map_to"],
            "Reactome_StId": term["stId"],
            "Reactome_DbId": term["dbId"],
            "Term_Name": term["name"],
            "Term_Species": term["species"],
            "Term_Type": term["type"]
        })

    if long_data:
        long_df = pd.DataFrame(long_data)

        # Save as CSV
        csv_file = f"{output_prefix}_{query['map_to']}.csv"
        long_df.to_csv(csv_file, index=False)
        logger.info(f"Saved {query['map_to']} to {csv_file}")

        # Save as Excel
        try:
            excel_file = f"{output_prefix}_{query['map_to']}.xlsx"
            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                long_df.to_excel(writer, sheet_name='Terms', index=False)

                # Add summary sheet
                summary_data = [{
                    "Gene_ID": query["gene_id"],
                    "Resource": query["resource"],
                    "Species": query["species"],
                    "Map_To": query["map_to"],
                    "Term_Count": len(data["terms"])
                }]
                summary_df = pd.DataFrame(summary_data)
                summary_df.to_excel(writer, sheet_name='Summary', index=False)

            logger.info(f"Saved {query['map_to']} to {excel_file}")
        except Exception as e:
            logger.warning(f"Failed to create Excel file: {e}")
    else:
        logger.warning(f"No {query['map_to']} found, skipping table export")


def format_output(data: Dict[str, any], output_format: str = "text") -> str:
    """Format output in the specified format."""
    if output_format == "json":
        return json.dumps(data, indent=2)

    # Text format
    output = []
    query = data["query"]
    output.append(f"Query: {query['gene_id']} ({query['resource']}) in species {query['species']}")
    output.append(f"Found {data['total_terms']} {query['map_to']}\n")

    if data["terms"]:
        for term in data["terms"]:
            output.append(f"{term['stId']}\t{term['name']} ({term['species']})")
    else:
        output.append(f"(no {query['map_to']} found)")

    return "\n".join(output)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Query Reactome database for pathways/reactions containing a given gene",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s TP53
  %(prog)s TP53 --resource UniProt --species 9606
  %(prog)s TP53 --map-to reactions
  %(prog)s PTEN --species mouse --output json

Species (common names or taxonomy IDs):
  human (9606), mouse (10090), rat (10116), zebrafish (7955)

Resources:
  UniProt, NCBI, ENSEMBL, etc.
        """
    )
    parser.add_argument("gene", help="Gene identifier (symbol or accession)")
    parser.add_argument(
        "-r", "--resource",
        default="UniProt",
        help="Identifier resource type (default: UniProt)"
    )
    parser.add_argument(
        "-s", "--species",
        default="9606",
        help="Species name or taxonomy ID (default: 9606 for human)"
    )
    parser.add_argument(
        "-m", "--map-to",
        choices=["pathways", "reactions"],
        default="pathways",
        help="Query type: pathways or reactions (default: pathways)"
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
        help="Output file prefix for CSV/Excel (default: reactome_{gene})"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Convert species to int if it's a number
    species = args.species
    if species.isdigit():
        species = int(species)

    try:
        result = get_reactome_terms_for_gene(
            gene_id=args.gene,
            resource=args.resource,
            species=species,
            map_to=args.map_to
        )

        # Export to CSV/Excel by default
        if args.output in ["csv", "text"]:
            output_prefix = args.output_prefix or f"reactome_{args.gene}"
            export_to_tables(result, output_prefix)

            # Print summary to console
            print(f"\n{'='*80}")
            print(f"Reactome Query: {args.gene}")
            print(f"{'='*80}")
            print(f"Resource: {result['query']['resource']}")
            print(f"Species: {result['query']['species']}")
            print(f"Map to: {result['query']['map_to']}")
            print(f"Total {result['query']['map_to']}: {result['total_terms']}")
            print(f"\nResults exported to:")
            print(f"  - {output_prefix}_{result['query']['map_to']}.csv")
            print(f"  - {output_prefix}_{result['query']['map_to']}.xlsx")
            print(f"{'='*80}\n")

        # JSON output
        if args.output == "json":
            output = format_output(result, args.output)
            print(output)

    except ValueError as e:
        logger.error(f"Invalid input: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error querying Reactome: {e}", exc_info=args.verbose)
        sys.exit(1)
