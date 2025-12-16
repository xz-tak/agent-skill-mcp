#!/usr/bin/env python3
"""
Test script for shortest path functions across all three interaction databases.

Tests:
1. STRING API - find_shortest_paths()
2. IntAct API - find_shortest_paths_intact()
3. BioGRID API - find_shortest_paths() (requires API key)
"""

import os
from string_api import StringClient
from intact_api import find_shortest_paths_intact
from biogrid_api import BioGRIDClient


def format_path_result(gene_a, gene_b, info):
    """Format a single path result for display."""
    path_str = " → ".join(info['path'])
    scores_str = ", ".join([f"{s:.3f}" if isinstance(s, float) else str(s) for s in info['scores']])

    output = []
    output.append(f"\n{gene_a} ↔ {gene_b}:")
    output.append(f"  Path: {path_str}")
    output.append(f"  Hops: {info['hops']}")
    output.append(f"  Distance: {info['distance']:.3f}")
    output.append(f"  Edge scores: [{scores_str}]")

    # Show algorithm info if available
    if 'algorithm' in info:
        output.append(f"  Algorithm: {info['algorithm']}")
        output.append(f"  Weight formula: {info['weight_formula']}")

    return "\n".join(output)


def test_string():
    """Test STRING API shortest paths."""
    print("="*70)
    print("TEST 1: STRING API - Shortest Paths")
    print("="*70)

    genes = ["TP53", "MDM2", "ATM"]
    print(f"\nQuery genes: {', '.join(genes)}")
    print(f"Parameters: species=9606 (human), max_distance=3, min_score=400")

    client = StringClient()

    try:
        paths = client.find_shortest_paths(
            gene_list=genes,
            species=9606,
            max_distance=3,
            min_combined_score=400,
            network_type="functional"
        )

        if paths:
            print(f"\n✓ Found {len(paths)} pairwise paths:")
            for (gene_a, gene_b), info in sorted(paths.items()):
                print(format_path_result(gene_a, gene_b, info))
        else:
            print("\n✗ No paths found")

    except Exception as e:
        print(f"\n✗ Error: {e}")

    print()


def test_intact():
    """Test IntAct API shortest paths."""
    print("="*70)
    print("TEST 2: IntAct API - Shortest Paths")
    print("="*70)

    genes = ["TP53", "MDM2", "EP300"]
    print(f"\nQuery genes: {', '.join(genes)}")
    print(f"Parameters: species='human', max_distance=3, min_miscore=0.4")

    try:
        paths = find_shortest_paths_intact(
            gene_list=genes,
            species="human",
            max_distance=3,
            min_miscore=0.4,
            miql_max_results=50000
        )

        if paths:
            print(f"\n✓ Found {len(paths)} pairwise paths:")
            for (gene_a, gene_b), info in sorted(paths.items()):
                print(format_path_result(gene_a, gene_b, info))
        else:
            print("\n✗ No paths found")

    except Exception as e:
        print(f"\n✗ Error: {e}")

    print()


def test_biogrid():
    """Test BioGRID API shortest paths."""
    print("="*70)
    print("TEST 3: BioGRID API - Shortest Paths")
    print("="*70)

    api_key = os.environ.get("BIOGRID_API_KEY")

    if not api_key:
        print("\n⚠ Skipping BioGRID test - No API key found")
        print("  Set BIOGRID_API_KEY environment variable to test")
        print()
        return

    genes = ["TP53", "MDM2", "ATM"]
    print(f"\nQuery genes: {', '.join(genes)}")
    print(f"Parameters: tax_id='9606', max_distance=3, min_score=0.5")

    client = BioGRIDClient(api_key)

    try:
        paths = client.find_shortest_paths(
            gene_list=genes,
            tax_id="9606",
            max_distance=3,
            min_score=0.5,
            experimental_system_types=None,
            throughput_tag="any"
        )

        if paths:
            print(f"\n✓ Found {len(paths)} pairwise paths:")
            for (gene_a, gene_b), info in sorted(paths.items()):
                print(format_path_result(gene_a, gene_b, info))
        else:
            print("\n✗ No paths found")

    except Exception as e:
        print(f"\n✗ Error: {e}")

    print()


def test_larger_set():
    """Test with a larger set of genes."""
    print("="*70)
    print("TEST 4: STRING API - Larger Gene Set")
    print("="*70)

    genes = ["TP53", "MDM2", "ATM", "CHEK2", "BRCA1"]
    print(f"\nQuery genes: {', '.join(genes)}")
    print(f"Finding shortest paths between all {len(genes) * (len(genes)-1) // 2} pairs...")

    client = StringClient()

    try:
        paths = client.find_shortest_paths(
            gene_list=genes,
            species=9606,
            max_distance=3,
            min_combined_score=500,
            network_type="functional"
        )

        if paths:
            print(f"\n✓ Found paths for {len(paths)} gene pairs:")

            # Group by hop count
            by_hops = {}
            for (gene_a, gene_b), info in paths.items():
                hops = info['hops']
                if hops not in by_hops:
                    by_hops[hops] = []
                by_hops[hops].append((gene_a, gene_b, info))

            for hops in sorted(by_hops.keys()):
                print(f"\n  {hops}-hop connections: {len(by_hops[hops])}")
                for gene_a, gene_b, info in by_hops[hops][:3]:  # Show first 3
                    path_str = " → ".join(info['path'])
                    print(f"    {gene_a} ↔ {gene_b}: {path_str}")

                if len(by_hops[hops]) > 3:
                    print(f"    ... and {len(by_hops[hops]) - 3} more")
        else:
            print("\n✗ No paths found")

    except Exception as e:
        print(f"\n✗ Error: {e}")

    print()


if __name__ == "__main__":
    print("\n" + "="*70)
    print("SHORTEST PATH TESTING - ALL DATABASES")
    print("="*70)
    print("\nThis script tests the shortest path finding functionality")
    print("for STRING, IntAct, and BioGRID interaction databases.")
    print()

    # Run tests
    test_string()
    test_intact()
    test_biogrid()
    test_larger_set()

    print("="*70)
    print("TESTING COMPLETE")
    print("="*70)
    print()
