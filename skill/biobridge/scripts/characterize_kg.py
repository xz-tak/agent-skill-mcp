#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Knowledge Graph Characterization Tool

This script analyzes the BioBridge knowledge graph to provide statistics on:
- Entity types and their counts
- Relation types and their frequencies
- Entity type combinations per relation
- Overall graph statistics

Usage:
    python characterize_kg.py [OPTIONS]

Examples:
    python characterize_kg.py
    python characterize_kg.py --kg-path /path/to/kg.csv
    python characterize_kg.py --export kg_stats.json
    python characterize_kg.py --show-combinations
"""

import argparse
import json
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Tuple

import pandas as pd


# Default paths
DEFAULT_KG_PATH = "/home/sagemaker-user/biobridge/bbridge/data/PrimeKG/kg.csv"


def load_kg(kg_path: str) -> pd.DataFrame:
    """Load knowledge graph."""
    print(f"Loading knowledge graph from: {kg_path}", file=sys.stderr)
    if not os.path.exists(kg_path):
        print(f"Error: KG file not found: {kg_path}", file=sys.stderr)
        sys.exit(1)

    kg = pd.read_csv(kg_path, low_memory=False)
    print(f"Loaded {len(kg):,} edges", file=sys.stderr)
    return kg


def characterize_entities(kg: pd.DataFrame) -> Dict:
    """Characterize entity types from KG edges."""
    entity_stats = {}

    # Collect all entity types from both head and tail
    head_types = kg["x_type"].tolist() if "x_type" in kg.columns else []
    tail_types = kg["y_type"].tolist() if "y_type" in kg.columns else []
    all_types = head_types + tail_types

    type_counts = Counter(all_types)

    entity_stats["unique_entity_types"] = len(type_counts)
    entity_stats["entity_type_occurrences"] = dict(type_counts)

    return entity_stats


def characterize_relations(kg: pd.DataFrame) -> Dict:
    """Characterize relations in the knowledge graph."""
    relation_stats = {}

    relation_stats["total_edges"] = len(kg)

    # Display relation types
    if "display_relation" in kg.columns:
        rel_counts = kg["display_relation"].value_counts()
        relation_stats["unique_relations"] = len(rel_counts)
        relation_stats["relation_distribution"] = rel_counts.to_dict()

    # Raw relation column
    if "relation" in kg.columns:
        raw_rel_counts = kg["relation"].value_counts()
        relation_stats["raw_relation_types"] = len(raw_rel_counts)
        relation_stats["raw_relation_distribution"] = raw_rel_counts.to_dict()

    return relation_stats


def characterize_entity_relation_combinations(kg: pd.DataFrame) -> Dict:
    """Characterize valid entity type combinations per relation."""
    combinations = defaultdict(set)

    if "display_relation" in kg.columns and "x_type" in kg.columns and "y_type" in kg.columns:
        for _, row in kg.iterrows():
            rel = row["display_relation"]
            x_type = row["x_type"]
            y_type = row["y_type"]
            combinations[rel].add((x_type, y_type))

    # Convert sets to lists for JSON serialization
    combinations_dict = {
        rel: [{"head_type": x, "tail_type": y} for x, y in sorted(combos)]
        for rel, combos in combinations.items()
    }

    return combinations_dict


def print_summary(entity_stats: Dict, relation_stats: Dict, combinations: Dict = None):
    """Print human-readable summary."""
    print("\n" + "="*80)
    print("BIOBRIDGE KNOWLEDGE GRAPH CHARACTERIZATION")
    print("="*80)
    print()

    # Entity statistics
    print("ENTITY TYPE STATISTICS:")
    print(f"  Unique Entity Types: {entity_stats.get('unique_entity_types', 0)}")
    print()

    print("  Entity Type Occurrences (in edges):")
    type_dist = entity_stats.get('entity_type_occurrences', {})
    total_occurrences = sum(type_dist.values())
    for etype, count in sorted(type_dist.items(), key=lambda x: x[1], reverse=True):
        pct = 100 * count / total_occurrences if total_occurrences > 0 else 0
        print(f"    {etype:<30} {count:>10,} ({pct:>5.1f}%)")
    print()

    # Relation statistics
    print("="*80)
    print("RELATION STATISTICS:")
    print(f"  Total Edges: {relation_stats.get('total_edges', 0):,}")
    print(f"  Unique Relation Types (display): {relation_stats.get('unique_relations', 0)}")
    if relation_stats.get('raw_relation_types'):
        print(f"  Unique Relation Types (raw): {relation_stats.get('raw_relation_types', 0)}")
    print()

    print("  Relation Distribution (Top 30):")
    rel_dist = relation_stats.get('relation_distribution', {})
    for idx, (rel, count) in enumerate(sorted(rel_dist.items(), key=lambda x: x[1], reverse=True)[:30], 1):
        pct = 100 * count / relation_stats.get('total_edges', 1)
        print(f"    {idx:>2}. {rel:<40} {count:>9,} ({pct:>5.1f}%)")
    print()

    # Entity-relation combinations
    if combinations:
        print("="*80)
        print("VALID ENTITY-RELATION COMBINATIONS:")
        print()

        for rel in sorted(combinations.keys()):
            combos = combinations[rel]
            print(f"  {rel}: ({len(combos)} combinations)")
            for combo in combos[:15]:  # Show first 15 combinations
                print(f"    {combo['head_type']:<30} → {combo['tail_type']}")
            if len(combos) > 15:
                print(f"    ... and {len(combos) - 15} more")
            print()

    print("="*80)


def main():
    parser = argparse.ArgumentParser(
        description="Characterize BioBridge knowledge graph",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic characterization
  python characterize_kg.py

  # Custom KG path
  python characterize_kg.py --kg-path /path/to/kg.csv

  # Export to JSON
  python characterize_kg.py --export kg_stats.json

  # Show entity-relation combinations
  python characterize_kg.py --show-combinations
        """
    )

    parser.add_argument(
        "--kg-path",
        default=DEFAULT_KG_PATH,
        help=f"Path to knowledge graph CSV (default: {DEFAULT_KG_PATH})"
    )
    parser.add_argument(
        "--export",
        help="Export statistics to JSON file"
    )
    parser.add_argument(
        "--show-combinations",
        action="store_true",
        help="Show entity-relation combinations (verbose)"
    )
    parser.add_argument(
        "--format",
        choices=["text", "json"],
        default="text",
        help="Output format (default: text)"
    )

    args = parser.parse_args()

    # Load data
    kg = load_kg(args.kg_path)
    print("Analyzing...\n", file=sys.stderr)

    # Characterize
    entity_stats = characterize_entities(kg)
    relation_stats = characterize_relations(kg)
    combinations = None
    if args.show_combinations or args.format == "json" or args.export:
        print("Computing entity-relation combinations...", file=sys.stderr)
        combinations = characterize_entity_relation_combinations(kg)

    # Output
    if args.format == "json" or args.export:
        output_data = {
            "kg_path": args.kg_path,
            "entity_statistics": entity_stats,
            "relation_statistics": relation_stats,
        }
        if combinations:
            output_data["entity_relation_combinations"] = combinations

        if args.export:
            with open(args.export, 'w') as f:
                json.dump(output_data, f, indent=2)
            print(f"\nStatistics exported to: {args.export}", file=sys.stderr)

        if args.format == "json":
            print(json.dumps(output_data, indent=2))
    else:
        # Text format
        print_summary(
            entity_stats,
            relation_stats,
            combinations if args.show_combinations else None
        )


if __name__ == "__main__":
    main()
