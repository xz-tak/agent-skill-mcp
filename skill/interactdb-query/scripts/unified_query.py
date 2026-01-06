#!/usr/bin/env python3
"""
Unified query interface for all three interaction databases.

This module provides wrapper functions that query STRING, IntAct, and BioGRID
by default, and always export results even when empty.

Features
--------
- Query all databases with a single function call
- Automatic result export (CSV for single-gene, JSON for paths)
- Handle empty results gracefully
- Support both single-gene and shortest-path queries

Usage
-----
    # Single gene query
    results = query_single_gene_all_databases("TP53", export_dir="./results")

    # Shortest path query
    results = query_shortest_paths_all_databases(["GREM1", "MRC2", "IL11"], export_dir="./results")
"""

import os
import json
import pandas as pd
from typing import Union, List, Dict, Any, Optional
from pathlib import Path

# Import from each database module
import sys
# Add the scripts directory to the Python path for imports
scripts_dir = Path(__file__).parent
if str(scripts_dir) not in sys.path:
    sys.path.insert(0, str(scripts_dir))

from string_api import StringClient
from intact_api import get_neighbors_multihop, find_shortest_paths_intact
try:
    from biogrid_api import BioGRIDClient
    BIOGRID_AVAILABLE = True
except ImportError:
    BIOGRID_AVAILABLE = False


def query_single_gene_all_databases(
    gene: str,
    species: Union[int, str] = 9606,
    top_n: int = 100,
    export_results: bool = True,
    output_dir: Optional[str] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Query a single gene across all three databases (STRING, IntAct, BioGRID).

    Parameters
    ----------
    gene : str
        Gene symbol to query (e.g., "TP53")
    species : int or str
        Species identifier: 9606 or "human" for Homo sapiens
    top_n : int
        Maximum number of neighbors to return per database
    export_results : bool
        Always export results to CSV, even if empty (default True)
    output_dir : Optional[str]
        Directory to save result files (default: current working directory)
    **kwargs : dict
        Additional database-specific parameters:
        - min_combined_score, min_experimental_score, etc. for STRING
        - min_miscore, organism_filter for IntAct
        - min_score for BioGRID

    Returns
    -------
    Dict[str, Any]
        {
            'string': DataFrame or error message,
            'intact': DataFrame or error message,
            'biogrid': List[NeighborRecord] or error message,
            'exports': {
                'string': file path or None,
                'intact': file path or None,
                'biogrid': file path or None
            }
        }
    """
    results = {}
    exports = {}

    # Use current working directory if no output_dir specified
    if output_dir is None:
        output_dir = os.getcwd()

    # Convert species string to appropriate format
    species_int = 9606 if species in (9606, "human", "homo sapiens") else int(species)
    species_str = "human" if species_int == 9606 else str(species_int)

    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Query STRING (using multi-hop BFS expansion)
    print(f"Querying STRING for {gene}...")
    try:
        client = StringClient()
        string_neighbors = client.get_neighbors_multihop(
            identifier=gene,
            species=species_int,
            top_n=top_n,
            max_hops=kwargs.get('max_hops', 5),
            min_combined_score=kwargs.get('min_combined_score', kwargs.get('min_score', 0)),
            min_experimental_score=kwargs.get('min_experimental_score', 0),
            min_database_score=kwargs.get('min_database_score', 0),
            min_textmining_score=kwargs.get('min_textmining_score', 0),
            min_coexpression_score=kwargs.get('min_coexpression_score', 0),
            network_type=kwargs.get('network_type', 'functional')
        )
        # Convert to DataFrame
        import pandas as pd
        from dataclasses import asdict
        string_df = pd.DataFrame([asdict(n) for n in string_neighbors])
        results['string'] = string_df

        if export_results:
            exports['string'] = export_neighbor_results(
                string_df, gene, "string", output_dir
            )
            print(f"✓ STRING: {len(string_df)} neighbors → {exports['string']}")
        else:
            exports['string'] = None
            print(f"✓ STRING: {len(string_df)} neighbors")
    except Exception as e:
        results['string'] = f"Error: {e}"
        exports['string'] = None
        print(f"✗ STRING: {e}")

    # Query IntAct (using multi-hop BFS expansion)
    print(f"Querying IntAct for {gene}...")
    try:
        intact_df = get_neighbors_multihop(
            gene=gene,
            species=species_str,
            top_n=top_n,
            max_hops=kwargs.get('max_hops', 5),
            min_miscore=kwargs.get('min_miscore', 0.0),
            organism_filter=kwargs.get('organism_filter', None),
            miql_max_results=kwargs.get('miql_max_results', 20000)
        )
        results['intact'] = intact_df

        if export_results:
            exports['intact'] = export_neighbor_results(
                intact_df, gene, "intact", output_dir
            )
            print(f"✓ IntAct: {len(intact_df)} neighbors → {exports['intact']}")
        else:
            exports['intact'] = None
            print(f"✓ IntAct: {len(intact_df)} neighbors")
    except Exception as e:
        results['intact'] = f"Error: {e}"
        exports['intact'] = None
        print(f"✗ IntAct: {e}")

    # Query BioGRID (if API key available)
    print(f"Querying BioGRID for {gene}...")
    if not BIOGRID_AVAILABLE:
        results['biogrid'] = "BioGRID module not available"
        exports['biogrid'] = None
        print("✗ BioGRID: Module not available")
    else:
        api_key = os.environ.get('BIOGRID_API_KEY')
        if not api_key:
            results['biogrid'] = "BioGRID API key not found (set BIOGRID_API_KEY environment variable)"
            exports['biogrid'] = None
            print("✗ BioGRID: No API key")
        else:
            try:
                client = BioGRIDClient(api_key)  # Default timeout: 5 minutes
                biogrid_neighbors = client.get_neighbors(
                    seed_gene=gene,
                    tax_id=str(species_int),
                    max_hops=kwargs.get('max_hops', 5),  # Multi-hop enabled (may take several minutes)
                    max_neighbors=top_n,
                    min_score=kwargs.get('min_score', None)
                )
                results['biogrid'] = biogrid_neighbors

                if export_results:
                    # Convert to DataFrame for export
                    from dataclasses import asdict
                    biogrid_df = pd.DataFrame([asdict(n) for n in biogrid_neighbors])
                    exports['biogrid'] = export_neighbor_results(
                        biogrid_df, gene, "biogrid", output_dir
                    )
                    print(f"✓ BioGRID: {len(biogrid_neighbors)} neighbors → {exports['biogrid']}")
                else:
                    exports['biogrid'] = None
                    print(f"✓ BioGRID: {len(biogrid_neighbors)} neighbors")
            except Exception as e:
                results['biogrid'] = f"Error: {e}"
                exports['biogrid'] = None
                print(f"✗ BioGRID: {e}")

    results['exports'] = exports
    return results


def query_shortest_paths_all_databases(
    gene_list: List[str],
    species: Union[int, str] = 9606,
    max_distance: int = 50,
    export_results: bool = True,
    output_dir: Optional[str] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Find shortest paths between genes across all three databases.

    Parameters
    ----------
    gene_list : List[str]
        List of gene symbols to find paths between
    species : int or str
        Species identifier: 9606 or "human" for Homo sapiens
    max_distance : int
        Maximum path length (default 50)
    export_results : bool
        Always export results to JSON, even if empty (default True)
    output_dir : Optional[str]
        Directory to save result files (default: current working directory)
    **kwargs : dict
        Additional database-specific parameters:
        - max_network_expansion: BFS expansion depth (default 5)
        - min_combined_score, min_experimental_score for STRING
        - min_miscore, organism_filter for IntAct
        - min_score for BioGRID

    Returns
    -------
    Dict[str, Any]
        {
            'string': Dict of paths or error message,
            'intact': Dict of paths or error message,
            'biogrid': Dict of paths or error message,
            'exports': {
                'string': file path or None,
                'intact': file path or None,
                'biogrid': file path or None
            }
        }
    """
    if len(gene_list) < 2:
        raise ValueError("gene_list must contain at least 2 genes")

    results = {}
    exports = {}

    # Use current working directory if no output_dir specified
    if output_dir is None:
        output_dir = os.getcwd()

    # Convert species string to appropriate format
    species_int = 9606 if species in (9606, "human", "homo sapiens") else int(species)
    species_str = "human" if species_int == 9606 else str(species_int)

    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Query STRING
    print(f"Finding paths in STRING for {gene_list}...")
    try:
        client = StringClient()
        string_paths = client.find_shortest_paths(
            gene_list=gene_list,
            species=species_int,
            max_distance=max_distance,
            max_network_expansion=kwargs.get('max_network_expansion', 5),
            min_combined_score=kwargs.get('min_combined_score', 400),
            min_experimental_score=kwargs.get('min_experimental_score', 0),
            min_database_score=kwargs.get('min_database_score', 0),
            min_textmining_score=kwargs.get('min_textmining_score', 0),
            min_coexpression_score=kwargs.get('min_coexpression_score', 0),
            network_type=kwargs.get('network_type', 'functional')
        )
        results['string'] = string_paths

        if export_results:
            exports['string'] = export_path_results(
                string_paths, gene_list, "string", output_dir
            )
            print(f"✓ STRING: {len(string_paths)} paths → {exports['string']}")
        else:
            exports['string'] = None
            print(f"✓ STRING: {len(string_paths)} paths")
    except Exception as e:
        results['string'] = f"Error: {e}"
        exports['string'] = None
        print(f"✗ STRING: {e}")

    # Query IntAct
    print(f"Finding paths in IntAct for {gene_list}...")
    try:
        intact_paths = find_shortest_paths_intact(
            gene_list=gene_list,
            species=species_str,
            max_distance=max_distance,
            min_miscore=kwargs.get('min_miscore', 0.4),
            organism_filter=kwargs.get('organism_filter', None),
            miql_max_results=kwargs.get('miql_max_results', 50000)
        )
        results['intact'] = intact_paths

        if export_results:
            exports['intact'] = export_path_results(
                intact_paths, gene_list, "intact", output_dir
            )
            print(f"✓ IntAct: {len(intact_paths)} paths → {exports['intact']}")
        else:
            exports['intact'] = None
            print(f"✓ IntAct: {len(intact_paths)} paths")
    except Exception as e:
        results['intact'] = f"Error: {e}"
        exports['intact'] = None
        print(f"✗ IntAct: {e}")

    # Query BioGRID (if API key available)
    print(f"Finding paths in BioGRID for {gene_list}...")
    if not BIOGRID_AVAILABLE:
        results['biogrid'] = "BioGRID module not available"
        exports['biogrid'] = None
        print("✗ BioGRID: Module not available")
    else:
        api_key = os.environ.get('BIOGRID_API_KEY')
        if not api_key:
            results['biogrid'] = "BioGRID API key not found (set BIOGRID_API_KEY environment variable)"
            exports['biogrid'] = None
            print("✗ BioGRID: No API key")
        else:
            try:
                client = BioGRIDClient(api_key)
                biogrid_paths = client.find_shortest_paths(
                    gene_list=gene_list,
                    tax_id=str(species_int),
                    max_distance=max_distance,
                    min_score=kwargs.get('min_score', 0)
                )
                results['biogrid'] = biogrid_paths

                if export_results:
                    exports['biogrid'] = export_path_results(
                        biogrid_paths, gene_list, "biogrid", output_dir
                    )
                    print(f"✓ BioGRID: {len(biogrid_paths)} paths → {exports['biogrid']}")
                else:
                    exports['biogrid'] = None
                    print(f"✓ BioGRID: {len(biogrid_paths)} paths")
            except Exception as e:
                results['biogrid'] = f"Error: {e}"
                exports['biogrid'] = None
                print(f"✗ BioGRID: {e}")

    results['exports'] = exports
    return results


def export_neighbor_results(
    df: pd.DataFrame,
    gene: str,
    database: str,
    output_dir: str
) -> str:
    """
    Export neighbor results to CSV, even if DataFrame is empty.

    Parameters
    ----------
    df : pd.DataFrame
        Neighbor results (may be empty)
    gene : str
        Query gene name
    database : str
        Database name (string, intact, biogrid)
    output_dir : str
        Output directory path

    Returns
    -------
    str
        Path to exported CSV file
    """
    output_file = os.path.join(output_dir, f"{gene}_{database}_neighbors.csv")

    # Export even if empty
    df.to_csv(output_file, index=False)

    return output_file


def export_path_results(
    paths: Dict,
    genes: List[str],
    database: str,
    output_dir: str
) -> str:
    """
    Export shortest path results to JSON, even if no paths found.

    Parameters
    ----------
    paths : Dict
        Path results (may be empty)
    genes : List[str]
        Query gene list
    database : str
        Database name (string, intact, biogrid)
    output_dir : str
        Output directory path

    Returns
    -------
    str
        Path to exported JSON file
    """
    genes_str = "-".join(genes)
    output_file = os.path.join(output_dir, f"{genes_str}_{database}_paths.json")

    # Create result structure
    if not paths or len(paths) == 0:
        result = {
            "query_genes": genes,
            "database": database,
            "paths": {},
            "message": "No paths found between query genes"
        }
    else:
        # Convert tuple keys to strings for JSON serialization
        serializable_paths = {
            f"{a}-{b}": {
                'path': info['path'],
                'hops': info['hops'],
                'distance': float(info['distance']),
                'scores': [float(s) if s is not None else None for s in info['scores']],
                'algorithm': info.get('algorithm', 'Dijkstra'),
                'weight_formula': info.get('weight_formula', 'N/A')
            }
            for (a, b), info in paths.items()
        }

        result = {
            "query_genes": genes,
            "database": database,
            "paths": serializable_paths,
            "num_paths": len(serializable_paths)
        }

    # Export to JSON
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2)

    return output_file


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Query all interaction databases with unified interface."
    )

    subparsers = parser.add_subparsers(dest='command', help='Query mode')

    # Single gene query
    single_parser = subparsers.add_parser('single', help='Query single gene')
    single_parser.add_argument('gene', help='Gene symbol (e.g., TP53)')
    single_parser.add_argument('--top-n', type=int, default=100, help='Max neighbors per database')
    single_parser.add_argument('--output-dir', default='.', help='Output directory')
    single_parser.add_argument('--species', default=9606, help='Species (9606 or "human")')
    single_parser.add_argument('--min-score', type=int, default=0, help='Minimum score filter')

    # Shortest paths query
    path_parser = subparsers.add_parser('paths', help='Find shortest paths')
    path_parser.add_argument('genes', nargs='+', help='Gene symbols (e.g., TP53 MDM2 ATM)')
    path_parser.add_argument('--max-distance', type=int, default=50, help='Max path length')
    path_parser.add_argument('--output-dir', default='.', help='Output directory')
    path_parser.add_argument('--species', default=9606, help='Species (9606 or "human")')
    path_parser.add_argument('--min-score', type=int, default=400, help='Minimum score filter')

    args = parser.parse_args()

    if args.command == 'single':
        print(f"\n{'='*80}")
        print(f"Single Gene Query: {args.gene}")
        print(f"{'='*80}\n")

        results = query_single_gene_all_databases(
            gene=args.gene,
            species=args.species,
            top_n=args.top_n,
            output_dir=args.output_dir,
            min_score=args.min_score
        )

        print(f"\n{'='*80}")
        print("SUMMARY")
        print(f"{'='*80}")
        for db, export_path in results['exports'].items():
            if export_path:
                print(f"{db.upper()}: {export_path}")

    elif args.command == 'paths':
        print(f"\n{'='*80}")
        print(f"Shortest Paths Query: {' ↔ '.join(args.genes)}")
        print(f"{'='*80}\n")

        results = query_shortest_paths_all_databases(
            gene_list=args.genes,
            species=args.species,
            max_distance=args.max_distance,
            output_dir=args.output_dir,
            min_combined_score=args.min_score
        )

        print(f"\n{'='*80}")
        print("SUMMARY")
        print(f"{'='*80}")
        for db, export_path in results['exports'].items():
            if export_path:
                print(f"{db.upper()}: {export_path}")

    else:
        parser.print_help()
