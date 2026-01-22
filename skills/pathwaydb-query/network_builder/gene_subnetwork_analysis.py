#!/usr/bin/env python3
"""
Extract gene-specific pathway subnetwork and compute centrality metrics.

Given a gene and its associated pathways, this script:
1. Extracts a subnetwork from the full pathway network
2. Computes comprehensive centrality metrics
3. Identifies hubs, authorities, and leaf nodes
4. Exports results to CSV/Excel with detailed statistics

Usage:
    python gene_subnetwork_analysis.py --gene TP53 \
        --network data/all_pathway_network_12052025.parquet \
        --pathways kegg_TP53_pathways.csv \
        --output tp53_subnetwork
"""

import sys
import logging
import argparse
import pandas as pd
import networkx as nx
from pathlib import Path
from typing import Dict, Set, List, Tuple
import warnings
warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_pathway_network(network_file: str, min_jaccard: float = 0.0) -> nx.Graph:
    """
    Load pathway network from Parquet or CSV file.

    Args:
        network_file: Path to network file (.parquet or .csv)
        min_jaccard: Minimum Jaccard index threshold

    Returns:
        NetworkX graph
    """
    logger.info(f"Loading pathway network from {network_file}")

    if network_file.endswith('.parquet'):
        df = pd.read_parquet(network_file)
    else:
        df = pd.read_csv(network_file)

    # Filter by Jaccard threshold
    if min_jaccard > 0:
        original_edges = len(df)
        df = df[df['Jaccard_Index'] >= min_jaccard]
        logger.info(f"Filtered edges: {original_edges:,} → {len(df):,} (Jaccard >= {min_jaccard})")

    # Add distance attribute as inverse of similarity
    # distance = 1/Jaccard (higher Jaccard = lower distance)
    df['Distance'] = 1.0 / df['Jaccard_Index']

    logger.info(f"Creating graph from {len(df):,} edges")
    G = nx.from_pandas_edgelist(
        df,
        source='Pathway1',
        target='Pathway2',
        edge_attr=['Jaccard_Index', 'Distance']
    )

    logger.info(f"Graph created: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")
    return G


def extract_pathways_from_query_results(pathway_files: List[str]) -> Set[str]:
    """
    Extract pathway names from pathway query result files.

    Args:
        pathway_files: List of CSV/Excel files from pathway queries

    Returns:
        Set of pathway names
    """
    pathways = set()

    for file_path in pathway_files:
        logger.info(f"Reading pathways from {file_path}")

        if file_path.endswith('.xlsx'):
            # Try reading different sheets
            try:
                df = pd.read_excel(file_path, sheet_name='Pathways')
            except:
                df = pd.read_excel(file_path)
        else:
            df = pd.read_csv(file_path)

        # Look for pathway name columns
        pathway_cols = [col for col in df.columns if 'pathway' in col.lower() or 'name' in col.lower()]

        if pathway_cols:
            for col in pathway_cols:
                pathways.update(df[col].dropna().unique())
        else:
            # Try first column
            pathways.update(df.iloc[:, 0].dropna().unique())

    logger.info(f"Extracted {len(pathways)} unique pathways")
    return pathways


def extract_subnetwork(G: nx.Graph, seed_pathways: Set[str],
                       include_neighbors: bool = True,
                       max_distance: int = 1) -> nx.Graph:
    """
    Extract subnetwork containing seed pathways and optionally their neighbors.

    Args:
        G: Full pathway network
        seed_pathways: Set of pathway names to include
        include_neighbors: Include direct neighbors of seed pathways
        max_distance: Maximum distance for neighbor inclusion

    Returns:
        Subnetwork as NetworkX graph
    """
    # Find seed pathways that exist in the network
    available_seeds = seed_pathways & set(G.nodes())
    missing = seed_pathways - available_seeds

    if missing:
        logger.warning(f"{len(missing)} pathways not found in network")
        if len(missing) <= 10:
            for p in list(missing)[:10]:
                logger.warning(f"  Missing: {p}")

    logger.info(f"Found {len(available_seeds)} seed pathways in network")

    if not available_seeds:
        raise ValueError("No seed pathways found in network")

    # Collect nodes to include
    nodes_to_include = set(available_seeds)

    if include_neighbors:
        for seed in available_seeds:
            # Get neighbors within max_distance
            if max_distance == 1:
                neighbors = set(G.neighbors(seed))
            else:
                neighbors = set()
                for node in nx.single_source_shortest_path_length(G, seed, cutoff=max_distance):
                    neighbors.add(node)
            nodes_to_include.update(neighbors)

        logger.info(f"Including {len(nodes_to_include) - len(available_seeds)} neighbors")

    # Extract subnetwork
    subG = G.subgraph(nodes_to_include).copy()
    logger.info(f"Subnetwork: {subG.number_of_nodes()} nodes, {subG.number_of_edges()} edges")

    return subG


def compute_centrality_metrics(G: nx.Graph, seed_pathways: Set[str] = None) -> pd.DataFrame:
    """
    Compute comprehensive centrality metrics for all nodes.

    Args:
        G: NetworkX graph
        seed_pathways: Optional set of seed pathways to mark

    Returns:
        DataFrame with centrality metrics
    """
    logger.info("Computing centrality metrics...")

    results = []

    # Basic metrics
    logger.info("  - Degree centrality")
    degree_cent = nx.degree_centrality(G)

    logger.info("  - Closeness centrality")
    try:
        # Use Distance for closeness (distance-based metric)
        closeness_cent = nx.closeness_centrality(G, distance='Distance', wf_improved=True)
    except:
        closeness_cent = {node: 0 for node in G.nodes()}

    logger.info("  - Betweenness centrality")
    # Use Distance for betweenness (distance-based, shortest path metric)
    betweenness_cent = nx.betweenness_centrality(G, weight='Distance')

    logger.info("  - Eigenvector centrality")
    try:
        eigen_cent = nx.eigenvector_centrality(G, weight='Jaccard_Index', max_iter=1000)
    except:
        logger.warning("  Eigenvector centrality failed, using zeros")
        eigen_cent = {node: 0 for node in G.nodes()}

    logger.info("  - PageRank")
    pagerank = nx.pagerank(G, weight='Jaccard_Index')

    # Hub and authority scores
    logger.info("  - HITS (Hubs and Authorities)")
    try:
        hits = nx.hits(G, max_iter=1000)
        hubs = hits[0]
        authorities = hits[1]
    except:
        logger.warning("  HITS failed, using zeros")
        hubs = {node: 0 for node in G.nodes()}
        authorities = {node: 0 for node in G.nodes()}

    # Clustering coefficient
    logger.info("  - Clustering coefficient")
    clustering = nx.clustering(G, weight='Jaccard_Index')

    # Compile results
    for node in G.nodes():
        degree = G.degree(node)

        # Classify node type
        if degree == 1:
            node_type = "Leaf"
        elif degree >= degree_cent[node] * len(G) * 0.5:  # High degree relative to average
            node_type = "Hub"
        elif betweenness_cent[node] > 0.1:
            node_type = "Bridge"
        else:
            node_type = "Regular"

        # Extract database from pathway name
        # Format: "pathway_name (DATABASE:id)"
        database = "Unknown"
        if '(' in node and ':' in node:
            try:
                database = node.split('(')[1].split(':')[0]
            except:
                pass

        # Check if seed pathway
        is_seed = node in seed_pathways if seed_pathways else False

        results.append({
            'Pathway': node,
            'Database': database,
            'Is_Seed': is_seed,
            'Node_Type': node_type,
            'Degree': degree,
            'Degree_Centrality': degree_cent[node],
            'Closeness_Centrality': closeness_cent[node],
            'Betweenness_Centrality': betweenness_cent[node],
            'Eigenvector_Centrality': eigen_cent[node],
            'PageRank': pagerank[node],
            'Hub_Score': hubs[node],
            'Authority_Score': authorities[node],
            'Clustering_Coefficient': clustering[node]
        })

    df = pd.DataFrame(results)
    df = df.sort_values('Degree_Centrality', ascending=False)

    logger.info(f"Computed metrics for {len(df)} pathways")
    return df


def compute_network_statistics(G: nx.Graph, centrality_df: pd.DataFrame) -> Dict:
    """
    Compute overall network statistics.

    Args:
        G: NetworkX graph
        centrality_df: DataFrame with centrality metrics

    Returns:
        Dictionary of statistics
    """
    stats = {
        'Nodes': G.number_of_nodes(),
        'Edges': G.number_of_edges(),
        'Density': nx.density(G),
        'Average_Degree': sum(dict(G.degree()).values()) / G.number_of_nodes(),
        'Average_Clustering': nx.average_clustering(G, weight='Jaccard_Index'),
        'N_Components': nx.number_connected_components(G),
    }

    # Largest component
    largest_cc = max(nx.connected_components(G), key=len)
    stats['Largest_Component_Size'] = len(largest_cc)
    stats['Largest_Component_Fraction'] = len(largest_cc) / G.number_of_nodes()

    # Centrality distributions
    stats['Mean_Degree_Centrality'] = centrality_df['Degree_Centrality'].mean()
    stats['Mean_Closeness_Centrality'] = centrality_df['Closeness_Centrality'].mean()
    stats['Mean_Betweenness_Centrality'] = centrality_df['Betweenness_Centrality'].mean()
    stats['Mean_Eigenvector_Centrality'] = centrality_df['Eigenvector_Centrality'].mean()

    # Node types
    stats['N_Hubs'] = (centrality_df['Node_Type'] == 'Hub').sum()
    stats['N_Bridges'] = (centrality_df['Node_Type'] == 'Bridge').sum()
    stats['N_Leaves'] = (centrality_df['Node_Type'] == 'Leaf').sum()
    stats['N_Regular'] = (centrality_df['Node_Type'] == 'Regular').sum()

    return stats


def export_results(centrality_df: pd.DataFrame, network_stats: Dict,
                   output_prefix: str, G: nx.Graph = None):
    """
    Export results to CSV, Excel, and optionally network file.

    Args:
        centrality_df: Centrality metrics DataFrame
        network_stats: Network statistics dictionary
        output_prefix: Output file prefix
        G: Optional NetworkX graph to export
    """
    # Export centrality table
    csv_file = f"{output_prefix}_centrality.csv"
    centrality_df.to_csv(csv_file, index=False)
    logger.info(f"Saved centrality table to {csv_file}")

    # Export to Excel with multiple sheets
    excel_file = f"{output_prefix}_analysis.xlsx"
    with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
        # Centrality metrics
        centrality_df.to_excel(writer, sheet_name='Centrality', index=False)

        # Network statistics
        stats_df = pd.DataFrame([network_stats]).T
        stats_df.columns = ['Value']
        stats_df.to_excel(writer, sheet_name='Network_Stats')

        # Top hubs
        top_hubs = centrality_df.nlargest(20, 'Degree_Centrality')[
            ['Pathway', 'Database', 'Node_Type', 'Degree', 'Degree_Centrality']
        ]
        top_hubs.to_excel(writer, sheet_name='Top_Hubs', index=False)

        # Top bridges
        top_bridges = centrality_df.nlargest(20, 'Betweenness_Centrality')[
            ['Pathway', 'Database', 'Node_Type', 'Degree', 'Betweenness_Centrality']
        ]
        top_bridges.to_excel(writer, sheet_name='Top_Bridges', index=False)

        # Seed pathways (if any)
        seed_pathways = centrality_df[centrality_df['Is_Seed']]
        if len(seed_pathways) > 0:
            seed_pathways.to_excel(writer, sheet_name='Seed_Pathways', index=False)

    logger.info(f"Saved comprehensive analysis to {excel_file}")

    # Export network for visualization
    if G is not None:
        edgelist_file = f"{output_prefix}_edgelist.csv"
        nx.to_pandas_edgelist(G).to_csv(edgelist_file, index=False)
        logger.info(f"Saved network edge list to {edgelist_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract gene-specific pathway subnetwork and compute centrality metrics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # With pathway query results
  %(prog)s --gene TP53 \\
    --network data/all_pathway_network_12052025.parquet \\
    --pathways kegg_TP53_pathways.csv reactome_TP53_pathways.csv \\
    --output tp53_subnetwork

  # With pre-extracted pathway list
  %(prog)s --gene BRCA1 \\
    --network data/all_pathway_network_12052025.parquet \\
    --pathway-list "p53 signaling,DNA repair,BRCA1 pathway" \\
    --output brca1_subnetwork

  # Include 2-hop neighbors
  %(prog)s --gene EGFR \\
    --network data/all_pathway_network_12052025.parquet \\
    --pathways egfr_pathways.csv \\
    --include-neighbors \\
    --max-distance 2 \\
    --output egfr_subnetwork
        """
    )

    parser.add_argument(
        '--gene',
        required=True,
        help='Gene symbol for labeling output'
    )
    parser.add_argument(
        '--network',
        required=True,
        help='Path to full pathway network (.parquet or .csv)'
    )
    parser.add_argument(
        '--pathways',
        nargs='+',
        help='Pathway query result files (CSV/Excel)'
    )
    parser.add_argument(
        '--pathway-list',
        help='Comma-separated list of pathway names'
    )
    parser.add_argument(
        '--min-jaccard',
        type=float,
        default=0.0,
        help='Minimum Jaccard index for network edges (default: 0.0)'
    )
    parser.add_argument(
        '--include-neighbors',
        action='store_true',
        help='Include neighboring pathways in subnetwork'
    )
    parser.add_argument(
        '--max-distance',
        type=int,
        default=1,
        help='Maximum distance for neighbor inclusion (default: 1)'
    )
    parser.add_argument(
        '--output',
        default='gene_subnetwork',
        help='Output file prefix (default: gene_subnetwork)'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    try:
        # Load full network
        G = load_pathway_network(args.network, min_jaccard=args.min_jaccard)

        # Extract seed pathways
        if args.pathways:
            seed_pathways = extract_pathways_from_query_results(args.pathways)
        elif args.pathway_list:
            seed_pathways = set([p.strip() for p in args.pathway_list.split(',')])
        else:
            raise ValueError("Must provide either --pathways or --pathway-list")

        logger.info(f"Gene: {args.gene}")
        logger.info(f"Seed pathways: {len(seed_pathways)}")

        # Extract subnetwork
        subG = extract_subnetwork(
            G,
            seed_pathways,
            include_neighbors=args.include_neighbors,
            max_distance=args.max_distance
        )

        # Compute centrality metrics
        centrality_df = compute_centrality_metrics(subG, seed_pathways)

        # Compute network statistics
        network_stats = compute_network_statistics(subG, centrality_df)
        network_stats['Gene'] = args.gene
        network_stats['Seed_Pathways'] = len(seed_pathways & set(subG.nodes()))

        # Export results
        export_results(centrality_df, network_stats, args.output, subG)

        # Print summary
        print(f"\n{'='*80}")
        print(f"Gene-Specific Pathway Subnetwork Analysis: {args.gene}")
        print(f"{'='*80}")
        print(f"Seed pathways: {network_stats['Seed_Pathways']}")
        print(f"Subnetwork nodes: {network_stats['Nodes']}")
        print(f"Subnetwork edges: {network_stats['Edges']}")
        print(f"Network density: {network_stats['Density']:.4f}")
        print(f"Average degree: {network_stats['Average_Degree']:.2f}")
        print(f"\nNode classification:")
        print(f"  Hubs: {network_stats['N_Hubs']}")
        print(f"  Bridges: {network_stats['N_Bridges']}")
        print(f"  Leaves: {network_stats['N_Leaves']}")
        print(f"  Regular: {network_stats['N_Regular']}")
        print(f"\nTop 5 most central pathways (by degree):")
        for idx, row in centrality_df.head(5).iterrows():
            print(f"  {row['Pathway'][:70]}")
            print(f"    Degree: {row['Degree']}, Centrality: {row['Degree_Centrality']:.4f}")
        print(f"\nResults saved to:")
        print(f"  - {args.output}_centrality.csv")
        print(f"  - {args.output}_analysis.xlsx")
        print(f"  - {args.output}_edgelist.csv")
        print(f"{'='*80}\n")

    except Exception as e:
        logger.error(f"Error: {e}", exc_info=args.verbose)
        sys.exit(1)


if __name__ == "__main__":
    main()
