#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Entity Exploration Tool for BioBridge Knowledge Graph

This script provides a quick way to explore entity names in the KG.
For actual predictions, use predict_link.py which uses LLM-based entity matching.

Usage:
    python match_entity.py [ENTITY_NAME] [OPTIONS]

Examples:
    python match_entity.py "IL11" --type "gene/protein"
    python match_entity.py "Crohn disease" --type disease
    python match_entity.py "gremlin" --limit 5
"""

import argparse
import json
import os
import re
import sys
from typing import List, Optional

import pandas as pd

# BioBridge data paths
BIOBRIDGE_DATA_DIR = "/home/sagemaker-user/biobridge/bbridge/data"
NODES_CSV_PATH = os.path.join(BIOBRIDGE_DATA_DIR, "PrimeKG", "nodes.csv")

ENTITY_TYPES = [
    "biological_process", "molecular_function", "cellular_component",
    "gene/protein", "disease", "drug", "pathway", "effect/phenotype",
    "anatomy", "exposure", "biologics_drug",
]


def _norm(s: str) -> str:
    """Normalize string to alphanumeric lowercase."""
    return re.sub(r"[^a-z0-9]+", "", str(s or "").lower())


def load_nodes() -> pd.DataFrame:
    """Load nodes.csv."""
    if not os.path.exists(NODES_CSV_PATH):
        print(f"Error: {NODES_CSV_PATH} not found", file=sys.stderr)
        sys.exit(1)
    return pd.read_csv(NODES_CSV_PATH, low_memory=False)


def quick_search(nodes: pd.DataFrame, query: str, entity_type: Optional[str], limit: int) -> pd.DataFrame:
    """Quick substring search."""
    df = nodes if not entity_type else nodes[nodes["node_type"] == entity_type]
    query_lower = query.lower()

    # Try exact match first
    exact = df[df["node_name"].astype(str).str.lower() == query_lower]
    if not exact.empty:
        return exact.head(limit)

    # Then contains
    contains = df[df["node_name"].astype(str).str.lower().str.contains(re.escape(query_lower), na=False)]
    return contains.head(limit)


def main():
    parser = argparse.ArgumentParser(description="Quick entity exploration (for full matching, use predict_link.py)")
    parser.add_argument("entity", help="Entity name")
    parser.add_argument("--type", choices=ENTITY_TYPES, help="Filter by type")
    parser.add_argument("--limit", type=int, default=10, help="Max results")
    parser.add_argument("--format", choices=["text", "json"], default="text")
    args = parser.parse_args()

    nodes = load_nodes()
    matches = quick_search(nodes, args.entity, args.type, args.limit)

    if matches.empty:
        print(f"No matches for '{args.entity}'", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(matches)} match(es):\n", file=sys.stderr)

    if args.format == "json":
        results = matches[["node_name", "node_type", "node_index"]].to_dict("records")
        print(json.dumps(results, indent=2))
    else:
        for _, row in matches.iterrows():
            print(f"- {row['node_name']} (type: {row['node_type']}, node_index: {row['node_index']})")


if __name__ == "__main__":
    main()
