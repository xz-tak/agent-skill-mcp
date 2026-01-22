#!/usr/bin/env python3
"""
Summarize and interpret gene-specific pathway subnetwork centrality analysis.

This script reads centrality analysis results and provides:
1. Executive summary of network structure
2. Interpretation of top pathways
3. Biological insights based on centrality metrics
4. Recommendations for further analysis

Usage:
    python summarize_centrality_analysis.py <centrality_csv> [--top N]
    python summarize_centrality_analysis.py examples/tp53_centrality.csv --top 10
"""

import sys
import pandas as pd
import argparse
from pathlib import Path


def load_centrality_data(csv_file: str) -> pd.DataFrame:
    """Load centrality analysis results."""
    return pd.read_csv(csv_file)


def summarize_network_structure(df: pd.DataFrame) -> dict:
    """Summarize overall network structure."""
    summary = {
        'total_pathways': len(df),
        'seed_pathways': (df['Is_Seed'] == True).sum(),
        'neighbor_pathways': (df['Is_Seed'] == False).sum(),
        'avg_degree': df['Degree'].mean(),
        'max_degree': df['Degree'].max(),
        'hub_count': (df['Node_Type'] == 'Hub').sum(),
        'bridge_count': (df['Node_Type'] == 'Bridge').sum(),
        'leaf_count': (df['Node_Type'] == 'Leaf').sum(),
        'regular_count': (df['Node_Type'] == 'Regular').sum()
    }

    # Database distribution
    summary['databases'] = df['Database'].value_counts().to_dict()

    return summary


def interpret_centrality_metrics(df: pd.DataFrame, top_n: int = 10) -> dict:
    """Interpret centrality metrics and identify key pathways."""
    interpretation = {}

    # Top pathways by different metrics
    interpretation['top_by_degree'] = df.nlargest(top_n, 'Degree_Centrality')[
        ['Pathway', 'Degree', 'Degree_Centrality', 'Database', 'Is_Seed']
    ].to_dict('records')

    interpretation['top_by_betweenness'] = df.nlargest(top_n, 'Betweenness_Centrality')[
        ['Pathway', 'Betweenness_Centrality', 'Node_Type', 'Database', 'Is_Seed']
    ].to_dict('records')

    interpretation['top_by_pagerank'] = df.nlargest(top_n, 'PageRank')[
        ['Pathway', 'PageRank', 'Degree', 'Database', 'Is_Seed']
    ].to_dict('records')

    interpretation['top_by_closeness'] = df.nlargest(top_n, 'Closeness_Centrality')[
        ['Pathway', 'Closeness_Centrality', 'Database', 'Is_Seed']
    ].to_dict('records')

    # Seed pathway metrics
    seed_df = df[df['Is_Seed'] == True].copy()
    if len(seed_df) > 0:
        interpretation['seed_pathways'] = seed_df.sort_values('Degree_Centrality', ascending=False)[
            ['Pathway', 'Degree', 'Degree_Centrality', 'Betweenness_Centrality', 'PageRank']
        ].to_dict('records')
    else:
        interpretation['seed_pathways'] = []

    return interpretation


def generate_biological_insights(summary: dict, interpretation: dict) -> list:
    """Generate biological insights based on network metrics."""
    insights = []

    # Network density insights
    if summary['total_pathways'] > 1:
        if summary['hub_count'] / summary['total_pathways'] > 0.7:
            insights.append({
                'type': 'Network Structure',
                'finding': 'Highly interconnected pathway module',
                'interpretation': 'Most pathways are hubs (>70%), suggesting the gene participates in a tightly integrated biological process or pathway module with extensive cross-talk.',
                'biological_significance': 'High degree of functional redundancy and coordination among pathways.'
            })
        elif summary['leaf_count'] / summary['total_pathways'] > 0.5:
            insights.append({
                'type': 'Network Structure',
                'finding': 'Sparse, loosely connected network',
                'interpretation': 'Many leaf nodes (>50%) indicate pathways with minimal connections, suggesting the gene is involved in diverse, independent biological processes.',
                'biological_significance': 'Gene has pleiotropic effects across multiple distinct pathways.'
            })

    # Bridge pathway insights
    if summary['bridge_count'] > 0:
        top_bridge = interpretation['top_by_betweenness'][0]
        insights.append({
            'type': 'Pathway Integration',
            'finding': f"Identified {summary['bridge_count']} bridge pathway(s)",
            'interpretation': f"Top bridge: '{top_bridge['Pathway'][:60]}...' connects different pathway modules.",
            'biological_significance': 'Bridge pathways are critical for information flow between distinct biological processes. Disruption may affect multiple downstream pathways.'
        })

    # Seed pathway centrality
    if interpretation['seed_pathways']:
        top_seed = interpretation['seed_pathways'][0]
        if top_seed['Degree_Centrality'] > 0.5:
            insights.append({
                'type': 'Gene Pathway Importance',
                'finding': f"Highly central seed pathway (centrality={top_seed['Degree_Centrality']:.3f})",
                'interpretation': f"'{top_seed['Pathway'][:60]}...' is a hub with {top_seed['Degree']} connections.",
                'biological_significance': 'This pathway is highly similar to many other pathways, suggesting it represents a core biological function with broad downstream effects.'
            })

    # Database diversity
    if len(summary['databases']) > 1:
        dominant_db = max(summary['databases'], key=summary['databases'].get)
        pct = summary['databases'][dominant_db] / summary['total_pathways'] * 100
        insights.append({
            'type': 'Database Coverage',
            'finding': f"{dominant_db} pathways dominate ({pct:.1f}%)",
            'interpretation': f"Network primarily contains {dominant_db} annotations.",
            'biological_significance': f"Results reflect {dominant_db}'s curation focus and granularity. Consider querying other databases for complementary perspectives."
        })

    return insights


def generate_recommendations(summary: dict, df: pd.DataFrame) -> list:
    """Generate recommendations for further analysis."""
    recommendations = []

    # Check if neighbors were included
    if summary['neighbor_pathways'] > 0:
        recommendations.append({
            'action': 'Compare with seed-only analysis',
            'rationale': 'You included neighbors. Compare with seed-only results to distinguish direct gene involvement from broader pathway context.'
        })
    else:
        if summary['total_pathways'] < 10:
            recommendations.append({
                'action': 'Consider adding --include-neighbors',
                'rationale': f'Only {summary["total_pathways"]} pathways in network. Including neighbors provides broader functional context.'
            })

    # Check for bridges
    if summary['bridge_count'] > 0:
        recommendations.append({
            'action': 'Investigate bridge pathways in detail',
            'rationale': 'Bridge pathways connect distinct modules. Understanding their function explains how the gene integrates multiple biological processes.'
        })

    # Visualization recommendation
    if summary['total_pathways'] >= 10:
        recommendations.append({
            'action': 'Visualize network in Cytoscape',
            'rationale': 'Import edgelist CSV to visualize pathway relationships. Color nodes by Database, size by Degree_Centrality.'
        })

    # Functional enrichment
    top_pathways = df.nlargest(20, 'Degree_Centrality')
    unique_dbs = top_pathways['Database'].nunique()
    if unique_dbs > 1:
        recommendations.append({
            'action': 'Perform cross-database functional analysis',
            'rationale': f'Top pathways span {unique_dbs} databases. Look for common themes across KEGG, MSigDB, and Reactome annotations.'
        })

    return recommendations


def print_summary(summary: dict, interpretation: dict, insights: list, recommendations: list, top_n: int):
    """Print formatted summary report."""
    print("=" * 80)
    print("PATHWAY CENTRALITY ANALYSIS SUMMARY")
    print("=" * 80)
    print()

    # Network structure
    print("## Network Structure")
    print("-" * 80)
    print(f"Total pathways:        {summary['total_pathways']:,}")
    print(f"  Seed pathways:       {summary['seed_pathways']} (directly associated with gene)")
    print(f"  Neighbor pathways:   {summary['neighbor_pathways']} (functionally related)")
    print()
    print(f"Node classification:")
    print(f"  Hubs:                {summary['hub_count']} ({summary['hub_count']/summary['total_pathways']*100:.1f}%)")
    print(f"  Bridges:             {summary['bridge_count']}")
    print(f"  Leaves:              {summary['leaf_count']}")
    print(f"  Regular:             {summary['regular_count']}")
    print()
    print(f"Connectivity:")
    print(f"  Average degree:      {summary['avg_degree']:.1f}")
    print(f"  Maximum degree:      {summary['max_degree']}")
    print()
    print(f"Database distribution:")
    for db, count in sorted(summary['databases'].items(), key=lambda x: x[1], reverse=True):
        pct = count / summary['total_pathways'] * 100
        print(f"  {db:15s}      {count:4d} ({pct:5.1f}%)")
    print()

    # Top pathways
    print("## Top Pathways by Centrality Metrics")
    print("-" * 80)

    print(f"\n### Top {top_n} by Degree Centrality (Most Connected)")
    print("  (High degree = many similar pathways, suggests core biological function)\n")
    for i, pw in enumerate(interpretation['top_by_degree'][:top_n], 1):
        seed_marker = " [SEED]" if pw['Is_Seed'] else ""
        pathway_short = pw['Pathway'][:65] + "..." if len(pw['Pathway']) > 65 else pw['Pathway']
        print(f"{i:2d}. {pathway_short}{seed_marker}")
        print(f"    Degree: {pw['Degree']:3d} | Centrality: {pw['Degree_Centrality']:.4f} | DB: {pw['Database']}")

    if interpretation['top_by_betweenness'][0]['Betweenness_Centrality'] > 0:
        print(f"\n### Top {min(5, top_n)} by Betweenness Centrality (Bridges)")
        print("  (High betweenness = lies on many shortest paths, connects different modules)\n")
        for i, pw in enumerate(interpretation['top_by_betweenness'][:5], 1):
            if pw['Betweenness_Centrality'] > 0:
                seed_marker = " [SEED]" if pw['Is_Seed'] else ""
                pathway_short = pw['Pathway'][:65] + "..." if len(pw['Pathway']) > 65 else pw['Pathway']
                print(f"{i:2d}. {pathway_short}{seed_marker}")
                print(f"    Betweenness: {pw['Betweenness_Centrality']:.4f} | Type: {pw['Node_Type']} | DB: {pw['Database']}")

    print(f"\n### Top {min(5, top_n)} by PageRank (Overall Importance)")
    print("  (High PageRank = connected to other important pathways, high influence)\n")
    for i, pw in enumerate(interpretation['top_by_pagerank'][:5], 1):
        seed_marker = " [SEED]" if pw['Is_Seed'] else ""
        pathway_short = pw['Pathway'][:65] + "..." if len(pw['Pathway']) > 65 else pw['Pathway']
        print(f"{i:2d}. {pathway_short}{seed_marker}")
        print(f"    PageRank: {pw['PageRank']:.6f} | Degree: {pw['Degree']:3d} | DB: {pw['Database']}")

    # Seed pathway analysis
    if interpretation['seed_pathways']:
        print(f"\n### Seed Pathways (Gene-Associated)")
        print("  (Pathways directly associated with the query gene)\n")
        for i, pw in enumerate(interpretation['seed_pathways'], 1):
            pathway_short = pw['Pathway'][:65] + "..." if len(pw['Pathway']) > 65 else pw['Pathway']
            print(f"{i:2d}. {pathway_short}")
            print(f"    Degree: {pw['Degree']:3d} | Centrality: {pw['Degree_Centrality']:.4f} | "
                  f"Betweenness: {pw['Betweenness_Centrality']:.4f} | PageRank: {pw['PageRank']:.6f}")

    print()

    # Biological insights
    if insights:
        print("## Biological Insights")
        print("-" * 80)
        print()
        for i, insight in enumerate(insights, 1):
            print(f"### {i}. {insight['type']}: {insight['finding']}")
            print(f"**Interpretation:** {insight['interpretation']}")
            print(f"**Biological Significance:** {insight['biological_significance']}")
            print()

    # Recommendations
    if recommendations:
        print("## Recommendations for Further Analysis")
        print("-" * 80)
        print()
        for i, rec in enumerate(recommendations, 1):
            print(f"{i}. **{rec['action']}**")
            print(f"   Rationale: {rec['rationale']}")
            print()

    print("=" * 80)
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Summarize and interpret pathway centrality analysis results",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'centrality_csv',
        help='Path to centrality CSV file (output from gene_subnetwork_analysis.py)'
    )
    parser.add_argument(
        '--top',
        type=int,
        default=10,
        help='Number of top pathways to show (default: 10)'
    )

    args = parser.parse_args()

    # Validate file exists
    if not Path(args.centrality_csv).exists():
        print(f"Error: File not found: {args.centrality_csv}", file=sys.stderr)
        sys.exit(1)

    # Load data
    df = load_centrality_data(args.centrality_csv)

    # Analyze
    summary = summarize_network_structure(df)
    interpretation = interpret_centrality_metrics(df, args.top)
    insights = generate_biological_insights(summary, interpretation)
    recommendations = generate_recommendations(summary, df)

    # Print report
    print_summary(summary, interpretation, insights, recommendations, args.top)


if __name__ == "__main__":
    main()
