#!/usr/bin/env python3
"""
Build complete PrimeKG graph schema from all edges.

This script processes ALL edges (train, valid, test) to build a complete
mapping of (head_type, relation) -> set of valid tail_types.

The resulting schema file is used for query validation and expected tail type
computation in the UltraQuery inference server.
"""
import sys
import os
import pickle
import json
from collections import defaultdict
from datetime import datetime

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

import polars as pl
from easydict import EasyDict
from ultra import util, datasets

# Configuration
DEFAULT_CONFIG = {
    "dataset": {"class": "PrimeKG1", "root": os.path.join(os.path.dirname(__file__), "data")},
}


def build_complete_schema(dataset_root: str, dataset_name: str):
    """
    Build complete graph schema from all edges in the dataset.

    Args:
        dataset_root: Path to data directory
        dataset_name: Name of dataset (e.g., "primekg1")

    Returns:
        Dictionary with:
        - schema: Dict mapping (h_type, relation) -> set of t_types
        - relation_to_types: Dict mapping relation -> list of (h_type, t_type) pairs
        - entity_type_counts: Dict mapping entity_type -> count
        - relation_counts: Dict mapping relation -> count
        - statistics: Overall statistics
    """
    print("=" * 80)
    print("Building Complete PrimeKG Graph Schema")
    print("=" * 80)

    # Load dataset
    print(f"\n[1/5] Loading dataset: {dataset_name}")
    cfg = EasyDict(DEFAULT_CONFIG)
    dataset = util.build_dataset(cfg)
    print(f"  ✓ Dataset loaded: {len(dataset[0].edge_index[0])} edges in training graph")

    # Load dictionaries
    print("\n[2/5] Loading entity/relation dictionaries...")
    with open(os.path.join(dataset_root, "id2ent_dict.pkl"), "rb") as f:
        id2ent_dict = pickle.load(f)
    with open(os.path.join(dataset_root, "id2rel_dict.pkl"), "rb") as f:
        id2rel_dict = pickle.load(f)
    print(f"  ✓ Loaded {len(id2ent_dict)} entities, {len(id2rel_dict)} relations")

    # Load nodes for type information
    print("\n[3/5] Loading node metadata...")
    dataset_dir = os.path.join(dataset_root, dataset_name.lower())
    nodes_path = os.path.join(dataset_dir, "raw", "nodes.txt")
    nodes_df = pl.read_csv(
        nodes_path,
        separator="\t",
        has_header=True,
        schema={
            "source_id": pl.Utf8,
            "name": pl.Utf8,
            "type": pl.Utf8,
            "source": pl.Utf8,
            "source_label": pl.Utf8,
        },
    )
    print(f"  ✓ Loaded {len(nodes_df)} nodes")

    # Create entity_id -> type mapping for fast lookup
    entity_to_type = dict(zip(nodes_df["source_label"], nodes_df["type"]))

    # Build schema from ALL edges
    print("\n[4/5] Processing all edges to build complete schema...")
    edges = dataset._data.target_edge_index.t()  # Shape: (num_edges, 2)
    edge_types = dataset._data.target_edge_type  # Shape: (num_edges,)
    num_edges = len(edges)

    print(f"  Processing {num_edges:,} edges...")

    # Schema: (h_type, relation) -> set of t_types
    schema = defaultdict(set)
    # Inverse map: relation -> list of (h_type, t_type) pairs
    relation_to_types = defaultdict(set)
    # Statistics
    entity_type_counts = defaultdict(int)
    relation_counts = defaultdict(int)

    # Track progress
    edges_processed = 0
    edges_with_types = 0
    report_interval = num_edges // 20  # Report every 5%

    for i in range(num_edges):
        h_idx, t_idx = edges[i].tolist()
        rel_idx = edge_types[i].item()

        # Get entity IDs and relation
        h_id = id2ent_dict[h_idx]
        t_id = id2ent_dict[t_idx]
        rel_label = id2rel_dict[rel_idx]

        # Get types
        h_type = entity_to_type.get(h_id)
        t_type = entity_to_type.get(t_id)

        if h_type and t_type:
            # Add to schema
            schema[(h_type, rel_label)].add(t_type)
            relation_to_types[rel_label].add((h_type, t_type))

            # Update counts
            entity_type_counts[h_type] += 1
            entity_type_counts[t_type] += 1
            relation_counts[rel_label] += 1

            edges_with_types += 1

        edges_processed += 1

        # Progress report
        if edges_processed % report_interval == 0:
            progress = (edges_processed / num_edges) * 100
            print(f"  Progress: {progress:.1f}% ({edges_processed:,}/{num_edges:,} edges)")

    print(f"  ✓ Processed {edges_processed:,} edges ({edges_with_types:,} with type information)")
    print(f"  ✓ Built schema with {len(schema)} (head_type, relation) pairs")
    print(f"  ✓ Found {len(relation_counts)} unique relations")
    print(f"  ✓ Found {len(entity_type_counts)} unique entity types")

    # Convert sets to sorted lists for JSON serialization
    schema_serializable = {
        f"{h_type}|{rel}": sorted(list(t_types))
        for (h_type, rel), t_types in schema.items()
    }

    relation_to_types_serializable = {
        rel: sorted([f"{h_type}|{t_type}" for h_type, t_type in pairs])
        for rel, pairs in relation_to_types.items()
    }

    # Compile statistics
    statistics = {
        "total_edges": num_edges,
        "edges_with_types": edges_with_types,
        "num_schema_pairs": len(schema),
        "num_relations": len(relation_counts),
        "num_entity_types": len(entity_type_counts),
        "generated_at": datetime.now().isoformat(),
    }

    result = {
        "schema": schema_serializable,
        "relation_to_types": relation_to_types_serializable,
        "entity_type_counts": dict(entity_type_counts),
        "relation_counts": dict(relation_counts),
        "statistics": statistics,
    }

    return result


def save_schema(schema_data: dict, output_path: str):
    """Save schema to JSON file."""
    print(f"\n[5/5] Saving schema to {output_path}...")

    with open(output_path, "w") as f:
        json.dump(schema_data, f, indent=2)

    # Calculate file size
    file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f"  ✓ Schema saved ({file_size_mb:.2f} MB)")


def print_schema_summary(schema_data: dict):
    """Print summary of the schema."""
    print("\n" + "=" * 80)
    print("Schema Summary")
    print("=" * 80)

    stats = schema_data["statistics"]
    print(f"\nStatistics:")
    print(f"  Total edges processed: {stats['total_edges']:,}")
    print(f"  Edges with type info: {stats['edges_with_types']:,}")
    print(f"  Schema pairs (head_type, relation): {stats['num_schema_pairs']:,}")
    print(f"  Unique relations: {stats['num_relations']}")
    print(f"  Unique entity types: {stats['num_entity_types']}")

    print(f"\nEntity Types ({len(schema_data['entity_type_counts'])}):")
    sorted_types = sorted(
        schema_data['entity_type_counts'].items(),
        key=lambda x: x[1],
        reverse=True
    )
    for entity_type, count in sorted_types:
        print(f"  {entity_type:30s} {count:>10,}")

    print(f"\nRelations ({len(schema_data['relation_counts'])}):")
    sorted_relations = sorted(
        schema_data['relation_counts'].items(),
        key=lambda x: x[1],
        reverse=True
    )
    for relation, count in sorted_relations[:20]:  # Top 20
        print(f"  {relation:40s} {count:>10,}")

    if len(sorted_relations) > 20:
        print(f"  ... and {len(sorted_relations) - 20} more relations")

    print("\n" + "=" * 80)


def main():
    """Main execution."""
    # Paths - script is now in src/, data is at ../data/
    script_dir = os.path.dirname(os.path.abspath(__file__))
    mcp_root = os.path.dirname(script_dir)  # Go up one level from src/
    dataset_root = os.path.join(mcp_root, "data")
    dataset_name = "primekg1"
    # Save schema to data directory
    data_dir = os.path.join(mcp_root, "data")
    os.makedirs(data_dir, exist_ok=True)
    output_path = os.path.join(data_dir, "primekg_schema.json")

    # Build schema
    schema_data = build_complete_schema(dataset_root, dataset_name)

    # Save schema
    save_schema(schema_data, output_path)

    # Print summary
    print_schema_summary(schema_data)

    print(f"\n✓ Complete! Schema saved to: {output_path}")


if __name__ == "__main__":
    main()
