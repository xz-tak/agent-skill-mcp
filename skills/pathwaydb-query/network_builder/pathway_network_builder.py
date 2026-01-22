#!/usr/bin/env python3
"""
Build a pathway similarity network using Jaccard indices.

This script:
1. Collects all pathways from KEGG, MSigDB, and Reactome with their associated genes
2. Computes pairwise Jaccard similarity (intersection/union) between all pathways
3. Outputs an edge list where nodes are pathways and edges are weighted by Jaccard index

Usage:
    python pathway_network_builder.py --output pathway_network.csv
    python pathway_network_builder.py --kegg-org hsa --msigdb-collection H --min-jaccard 0.1
"""

import sys
import logging
import argparse
import pandas as pd
import requests
from typing import Dict, List, Set, Tuple
from itertools import combinations
from functools import lru_cache
import time

# Import from existing modules
sys.path.append('/home/sagemaker-user/.claude/skills/pathwaydb-query/scripts')
from kegg_api import get_session_with_retry, load_pathway_names, KEGG_REST, REQUEST_TIMEOUT
from msigdb_api import _load_gmt_cached, COLLECTION_TO_CATEGORY
from reactome_api import BASE_URL as REACTOME_BASE_URL

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def get_genes_for_kegg_pathway(pathway_id: str, organism: str, session: requests.Session = None) -> Set[str]:
    """
    Get all genes for a KEGG pathway.

    Args:
        pathway_id: KEGG pathway ID (e.g., 'hsa04115' or 'path:hsa04115')
        organism: KEGG organism code (e.g., 'hsa')
        session: Optional requests session

    Returns:
        Set of gene IDs (e.g., 'hsa:7157')
    """
    if session is None:
        session = get_session_with_retry()

    # Remove 'path:' prefix if present
    pid = pathway_id.replace('path:', '')

    # Use link/{organism}/{pathway} endpoint
    url = f"{KEGG_REST}/link/{organism}/{pid}"
    logger.debug(f"Fetching genes for pathway: {url}")

    try:
        resp = session.get(url, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
    except requests.exceptions.RequestException as e:
        logger.warning(f"Failed to get genes for pathway {pathway_id}: {e}")
        return set()

    genes = set()
    for line in resp.text.strip().splitlines():
        if not line.strip():
            continue
        try:
            _, gene_id = line.split("\t", 1)
            # gene_id format: 'hsa:7157'
            genes.add(gene_id)
        except ValueError:
            logger.warning(f"Failed to parse line: {line}")
            continue

    return genes


def collect_kegg_pathways(organism: str = "hsa", max_pathways: int = None) -> Dict[str, Set[str]]:
    """
    Collect all KEGG pathways for an organism with their gene sets.

    Args:
        organism: KEGG organism code (default: 'hsa' for human)
        max_pathways: Maximum number of pathways to collect (None = all)

    Returns:
        Dictionary mapping pathway names to gene sets
        Key format: "pathway_name (KEGG:pathway_id)"
    """
    logger.info(f"Collecting KEGG pathways for organism: {organism}")
    session = get_session_with_retry()

    # Get all pathway names
    pathway_names = load_pathway_names(organism, session=session)
    logger.info(f"Found {len(pathway_names)} KEGG pathways")

    if max_pathways:
        logger.info(f"Limiting to first {max_pathways} pathways for testing")
        pathway_names = dict(list(pathway_names.items())[:max_pathways])

    pathways = {}
    for i, (pathway_id, pathway_name) in enumerate(pathway_names.items(), 1):
        if i % 10 == 0:
            logger.info(f"Processing KEGG pathway {i}/{len(pathway_names)}")

        # Get genes for this pathway
        genes = get_genes_for_kegg_pathway(pathway_id, organism=organism, session=session)

        # Format: "pathway_name (KEGG:pathway_id)"
        key = f"{pathway_name} (KEGG:{pathway_id})"
        pathways[key] = genes

        # Rate limiting
        time.sleep(0.1)

    logger.info(f"Collected {len(pathways)} KEGG pathways with genes")
    return pathways


def collect_msigdb_pathways(collections: List[str] = None, dbver: str = "2025.1.Hs") -> Dict[str, Set[str]]:
    """
    Collect all MSigDB gene sets from one or more collections.

    Args:
        collections: List of MSigDB collection codes (H, C1-C8). If None, collects all.
        dbver: MSigDB database version

    Returns:
        Dictionary mapping gene set names to gene sets
        Key format: "gene_set_name (MSigDB:collection)"
    """
    # Default to all collections if not specified
    if collections is None:
        collections = list(COLLECTION_TO_CATEGORY.keys())

    # Ensure collections is a list
    if isinstance(collections, str):
        collections = [collections]

    all_pathways = {}

    for collection in collections:
        logger.info(f"Collecting MSigDB gene sets from collection: {collection}")

        collection_upper = collection.upper()
        if collection_upper not in COLLECTION_TO_CATEGORY:
            logger.warning(f"Invalid MSigDB collection: {collection}, skipping")
            continue

        category = COLLECTION_TO_CATEGORY[collection_upper]

        # Load GMT file
        try:
            gmt = _load_gmt_cached(category, dbver)
            logger.info(f"Loaded {len(gmt)} gene sets from MSigDB {collection}")

            for gene_set_name, gene_list in gmt.items():
                # Format: "gene_set_name (MSigDB:collection)"
                key = f"{gene_set_name} (MSigDB:{collection_upper})"
                all_pathways[key] = set(gene_list)
        except Exception as e:
            logger.error(f"Failed to load collection {collection}: {e}")
            continue

    logger.info(f"Collected {len(all_pathways)} MSigDB gene sets from {len(collections)} collection(s)")
    return all_pathways


def collect_reactome_pathways(species: int = 9606) -> Dict[str, Set[str]]:
    """
    Collect all Reactome pathways with their genes.

    Args:
        species: NCBI taxonomy ID (default: 9606 for human)

    Returns:
        Dictionary mapping pathway names to gene sets
        Key format: "pathway_name (Reactome:stId)"
    """
    logger.info(f"Collecting Reactome pathways for species: {species}")
    session = get_session_with_retry()

    # Get all top-level pathways for species
    url = f"{REACTOME_BASE_URL}/data/pathways/top/{species}"
    logger.debug(f"Fetching Reactome top pathways: {url}")

    try:
        resp = session.get(url, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        top_pathways = resp.json()
    except Exception as e:
        logger.error(f"Failed to get Reactome top pathways: {e}")
        return {}

    logger.info(f"Found {len(top_pathways)} top-level Reactome pathways")

    # For each pathway, get its genes
    pathways = {}
    for i, pathway in enumerate(top_pathways, 1):
        st_id = pathway.get('stId')
        name = pathway.get('displayName', st_id)

        if i % 10 == 0:
            logger.info(f"Processing Reactome pathway {i}/{len(top_pathways)}")

        # Get entities (genes) for this pathway
        try:
            entity_url = f"{REACTOME_BASE_URL}/data/pathway/{st_id}/participatingPhysicalEntities"
            resp = session.get(entity_url, timeout=REQUEST_TIMEOUT)
            resp.raise_for_status()
            entities = resp.json()

            # Extract gene names
            genes = set()
            for entity in entities:
                # Try to get gene names from various fields
                gene_name = entity.get('geneName')
                if gene_name:
                    if isinstance(gene_name, list):
                        genes.update(gene_name)
                    else:
                        genes.add(gene_name)

            # Format: "pathway_name (Reactome:stId)"
            key = f"{name} (Reactome:{st_id})"
            pathways[key] = genes

            time.sleep(0.1)  # Rate limiting

        except Exception as e:
            logger.warning(f"Failed to get genes for Reactome pathway {st_id}: {e}")
            continue

    logger.info(f"Collected {len(pathways)} Reactome pathways with genes")
    return pathways


def jaccard_similarity(set1: Set[str], set2: Set[str]) -> float:
    """
    Calculate Jaccard similarity between two sets.

    Jaccard = |intersection| / |union|

    Args:
        set1: First gene set
        set2: Second gene set

    Returns:
        Jaccard index (0 to 1)
    """
    if not set1 or not set2:
        return 0.0

    intersection = len(set1 & set2)
    union = len(set1 | set2)

    if union == 0:
        return 0.0

    return intersection / union


def build_pathway_network(
    pathways: Dict[str, Set[str]],
    min_jaccard: float = 0.0
) -> List[Tuple[str, str, float]]:
    """
    Build a pathway similarity network based on Jaccard indices.

    Args:
        pathways: Dictionary mapping pathway names to gene sets
        min_jaccard: Minimum Jaccard index to include an edge (default: 0.0)

    Returns:
        List of edges as (pathway1, pathway2, jaccard_weight)
    """
    logger.info(f"Building pathway network from {len(pathways)} pathways")
    logger.info(f"Total possible edges: {len(pathways) * (len(pathways) - 1) // 2}")

    edges = []
    pathway_names = list(pathways.keys())
    total_pairs = len(pathway_names) * (len(pathway_names) - 1) // 2

    for i, (pw1, pw2) in enumerate(combinations(pathway_names, 2), 1):
        if i % 100000 == 0:
            logger.info(f"Processed {i}/{total_pairs} pathway pairs ({i/total_pairs*100:.1f}%)")

        jaccard = jaccard_similarity(pathways[pw1], pathways[pw2])

        if jaccard >= min_jaccard:
            edges.append((pw1, pw2, jaccard))

    logger.info(f"Generated {len(edges)} edges with Jaccard >= {min_jaccard}")
    return edges


def export_edge_list(edges: List[Tuple[str, str, float]], output_file: str) -> None:
    """
    Export edge list to parquet file (or CSV if specified).

    Args:
        edges: List of (pathway1, pathway2, weight) tuples
        output_file: Output filename (.parquet or .csv)
    """
    df = pd.DataFrame(edges, columns=['Pathway1', 'Pathway2', 'Jaccard_Index'])
    df = df.sort_values('Jaccard_Index', ascending=False)

    # Determine format from file extension
    if output_file.endswith('.csv'):
        df.to_csv(output_file, index=False)
        logger.info(f"Saved edge list to {output_file} (CSV format)")
    else:
        # Default to parquet
        if not output_file.endswith('.parquet'):
            output_file = output_file + '.parquet'
        df.to_parquet(output_file, engine='pyarrow', compression='snappy', index=False)
        logger.info(f"Saved edge list to {output_file} (Parquet format)")

    # Print summary statistics
    logger.info(f"\nNetwork Summary:")
    logger.info(f"  Total edges: {len(edges)}")
    logger.info(f"  Mean Jaccard: {df['Jaccard_Index'].mean():.4f}")
    logger.info(f"  Median Jaccard: {df['Jaccard_Index'].median():.4f}")
    logger.info(f"  Max Jaccard: {df['Jaccard_Index'].max():.4f}")
    logger.info(f"  Min Jaccard: {df['Jaccard_Index'].min():.4f}")


def main():
    parser = argparse.ArgumentParser(
        description="Build pathway similarity network using Jaccard indices",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --output pathway_network.csv
  %(prog)s --kegg-org hsa --msigdb-collections H --min-jaccard 0.1
  %(prog)s --databases kegg msigdb --msigdb-collections H C2 --output network.csv
  %(prog)s --databases reactome --species 9606 --min-jaccard 0.05
  %(prog)s --databases msigdb --msigdb-collections H --output hallmark_only.csv
        """
    )

    parser.add_argument(
        '--databases',
        nargs='+',
        choices=['kegg', 'msigdb', 'reactome'],
        default=['kegg', 'msigdb', 'reactome'],
        help='Databases to include (default: all)'
    )
    parser.add_argument(
        '--kegg-org',
        default='hsa',
        help='KEGG organism code (default: hsa)'
    )
    parser.add_argument(
        '--msigdb-collections',
        nargs='+',
        choices=list(COLLECTION_TO_CATEGORY.keys()),
        default=None,
        help='MSigDB collections to include (default: all collections H C1-C8). Example: --msigdb-collections H C2'
    )
    parser.add_argument(
        '--msigdb-version',
        default='2025.1.Hs',
        help='MSigDB version (default: 2025.1.Hs)'
    )
    parser.add_argument(
        '--species',
        type=int,
        default=9606,
        help='Species taxonomy ID for Reactome (default: 9606 for human)'
    )
    parser.add_argument(
        '--min-jaccard',
        type=float,
        default=0.0,
        help='Minimum Jaccard index to include edge (default: 0.0)'
    )
    parser.add_argument(
        '--max-pathways-per-db',
        type=int,
        default=None,
        help='Maximum pathways to collect per database (for testing, default: all)'
    )
    parser.add_argument(
        '--output',
        default='pathway_network_edges.parquet',
        help='Output file (default: pathway_network_edges.parquet). Use .csv extension for CSV format, .parquet for Parquet (default).'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Collect pathways from selected databases
    all_pathways = {}

    try:
        if 'kegg' in args.databases:
            kegg_pathways = collect_kegg_pathways(
                organism=args.kegg_org,
                max_pathways=args.max_pathways_per_db
            )
            all_pathways.update(kegg_pathways)

        if 'msigdb' in args.databases:
            msigdb_pathways = collect_msigdb_pathways(
                collections=args.msigdb_collections,
                dbver=args.msigdb_version
            )
            all_pathways.update(msigdb_pathways)

        if 'reactome' in args.databases:
            reactome_pathways = collect_reactome_pathways(species=args.species)
            all_pathways.update(reactome_pathways)

        if not all_pathways:
            logger.error("No pathways collected from any database")
            sys.exit(1)

        logger.info(f"\nTotal pathways collected: {len(all_pathways)}")

        # Build network
        edges = build_pathway_network(
            pathways=all_pathways,
            min_jaccard=args.min_jaccard
        )

        # Export edge list
        export_edge_list(edges, args.output)

        print(f"\n{'='*80}")
        print(f"Pathway Network Built Successfully!")
        print(f"{'='*80}")
        print(f"Total pathways: {len(all_pathways)}")
        print(f"Total edges: {len(edges)}")
        print(f"Minimum Jaccard threshold: {args.min_jaccard}")
        print(f"\nOutput saved to: {args.output}")
        print(f"{'='*80}\n")

    except Exception as e:
        logger.error(f"Error building pathway network: {e}", exc_info=args.verbose)
        sys.exit(1)


if __name__ == "__main__":
    main()
