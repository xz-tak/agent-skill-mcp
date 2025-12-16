#!/usr/bin/env python3
"""
Example script showing how to analyze the pathway similarity network.

This demonstrates basic network analysis tasks on the edge list output.
"""

import pandas as pd
import sys

def analyze_network(edge_list_file: str):
    """Analyze pathway network from edge list file (CSV or Parquet)."""

    # Load edge list
    print(f"Loading network from {edge_list_file}...")

    # Auto-detect format from extension
    if edge_list_file.endswith('.parquet'):
        df = pd.read_parquet(edge_list_file)
        print(f"  Format: Parquet")
    elif edge_list_file.endswith('.csv'):
        df = pd.read_csv(edge_list_file)
        print(f"  Format: CSV")
    else:
        # Try parquet first, then CSV
        try:
            df = pd.read_parquet(edge_list_file)
            print(f"  Format: Parquet (auto-detected)")
        except:
            df = pd.read_csv(edge_list_file)
            print(f"  Format: CSV (auto-detected)")

    print(f"\n{'='*80}")
    print("NETWORK STATISTICS")
    print(f"{'='*80}")

    # Basic stats
    print(f"\nTotal edges: {len(df)}")

    # Get unique pathways
    pathways = set(df['Pathway1'].tolist() + df['Pathway2'].tolist())
    print(f"Total nodes (pathways): {len(pathways)}")

    # Jaccard statistics
    print(f"\nJaccard Index Statistics:")
    print(f"  Mean:   {df['Jaccard_Index'].mean():.4f}")
    print(f"  Median: {df['Jaccard_Index'].median():.4f}")
    print(f"  Std:    {df['Jaccard_Index'].std():.4f}")
    print(f"  Min:    {df['Jaccard_Index'].min():.4f}")
    print(f"  Max:    {df['Jaccard_Index'].max():.4f}")

    # Extract database from pathway names
    # Format: "pathway_name (DATABASE:id)"
    df['DB1'] = df['Pathway1'].str.extract(r'\(([^:]+):')[0]
    df['DB2'] = df['Pathway2'].str.extract(r'\(([^:]+):')[0]

    # Count edges by database
    print(f"\n{'='*80}")
    print("EDGES BY DATABASE")
    print(f"{'='*80}")

    for db1 in sorted(df['DB1'].unique()):
        for db2 in sorted(df['DB2'].unique()):
            if db1 <= db2:  # Only count once for undirected edges
                mask = ((df['DB1'] == db1) & (df['DB2'] == db2)) | \
                       ((df['DB1'] == db2) & (df['DB2'] == db1))
                count = mask.sum()
                if count > 0:
                    avg_jaccard = df[mask]['Jaccard_Index'].mean()
                    print(f"{db1:10s} - {db2:10s}: {count:6d} edges (avg Jaccard: {avg_jaccard:.4f})")

    # Top 10 most similar pathway pairs
    print(f"\n{'='*80}")
    print("TOP 10 MOST SIMILAR PATHWAY PAIRS")
    print(f"{'='*80}\n")

    top_pairs = df.nlargest(10, 'Jaccard_Index')
    for idx, row in top_pairs.iterrows():
        pw1 = row['Pathway1']
        pw2 = row['Pathway2']
        jaccard = row['Jaccard_Index']

        # Shorten pathway names for display
        pw1_short = pw1[:60] + '...' if len(pw1) > 60 else pw1
        pw2_short = pw2[:60] + '...' if len(pw2) > 60 else pw2

        print(f"Jaccard: {jaccard:.4f}")
        print(f"  → {pw1_short}")
        print(f"  → {pw2_short}")
        print()

    # Degree distribution
    print(f"\n{'='*80}")
    print("DEGREE DISTRIBUTION (Top 10 most connected pathways)")
    print(f"{'='*80}\n")

    # Count connections per pathway
    pathway_degree = {}
    for _, row in df.iterrows():
        pathway_degree[row['Pathway1']] = pathway_degree.get(row['Pathway1'], 0) + 1
        pathway_degree[row['Pathway2']] = pathway_degree.get(row['Pathway2'], 0) + 1

    # Sort by degree
    top_pathways = sorted(pathway_degree.items(), key=lambda x: x[1], reverse=True)[:10]

    for pathway, degree in top_pathways:
        pathway_short = pathway[:70] + '...' if len(pathway) > 70 else pathway
        print(f"{degree:4d} connections: {pathway_short}")

    print(f"\n{'='*80}\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analyze_network_example.py <edge_list_file>")
        print("\nSupports both CSV and Parquet formats:")
        print("  python analyze_network_example.py pathway_network.csv")
        print("  python analyze_network_example.py pathway_network.parquet")
        sys.exit(1)

    edge_list_file = sys.argv[1]
    analyze_network(edge_list_file)
