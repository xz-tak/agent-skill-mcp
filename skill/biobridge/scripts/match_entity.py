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
from difflib import SequenceMatcher
from typing import List, Optional, Tuple

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


def similarity_score(s1: str, s2: str) -> float:
    """Calculate similarity score between two strings (0-1)."""
    return SequenceMatcher(None, s1.lower(), s2.lower()).ratio()


def load_nodes() -> pd.DataFrame:
    """Load nodes.csv."""
    if not os.path.exists(NODES_CSV_PATH):
        print(f"Error: {NODES_CSV_PATH} not found", file=sys.stderr)
        sys.exit(1)
    return pd.read_csv(NODES_CSV_PATH, low_memory=False)


def fuzzy_search(nodes: pd.DataFrame, query: str, entity_type: Optional[str], limit: int) -> Tuple[pd.DataFrame, List[float]]:
    """
    Enhanced search with fuzzy matching.
    Returns matches and their similarity scores.
    """
    df = nodes if not entity_type else nodes[nodes["node_type"] == entity_type]
    query_lower = query.lower()

    # Strategy 1: Try exact match first
    exact = df[df["node_name"].astype(str).str.lower() == query_lower]
    if not exact.empty:
        scores = [1.0] * len(exact.head(limit))
        return exact.head(limit), scores

    # Strategy 2: Try contains match
    contains = df[df["node_name"].astype(str).str.lower().str.contains(re.escape(query_lower), na=False)]
    if not contains.empty and len(contains) <= limit:
        # Calculate similarity scores for ranking
        contains_list = contains.copy()
        contains_list['similarity'] = contains_list['node_name'].apply(
            lambda x: similarity_score(query, str(x))
        )
        contains_sorted = contains_list.sort_values('similarity', ascending=False)
        scores = contains_sorted['similarity'].head(limit).tolist()
        return contains_sorted.head(limit), scores

    # Strategy 3: Fuzzy matching on all entities (more expensive)
    if df.empty:
        return pd.DataFrame(), []

    # Calculate similarity scores for all entities
    df_copy = df.copy()
    df_copy['similarity'] = df_copy['node_name'].apply(
        lambda x: similarity_score(query, str(x))
    )

    # Sort by similarity and take top matches
    df_sorted = df_copy.sort_values('similarity', ascending=False)
    top_matches = df_sorted.head(limit * 3)  # Get more candidates

    # Filter to only keep reasonable matches (similarity > 0.4)
    good_matches = top_matches[top_matches['similarity'] > 0.4]

    if good_matches.empty:
        # If no good fuzzy matches, return top matches anyway
        scores = top_matches['similarity'].head(limit).tolist()
        return top_matches.head(limit), scores

    scores = good_matches['similarity'].head(limit).tolist()
    return good_matches.head(limit), scores


def quick_search(nodes: pd.DataFrame, query: str, entity_type: Optional[str], limit: int) -> pd.DataFrame:
    """Quick substring search (legacy function)."""
    matches, _ = fuzzy_search(nodes, query, entity_type, limit)
    return matches


def format_match_text(matches: pd.DataFrame, scores: List[float], query: str, show_details: bool = False) -> str:
    """Format matches as numbered list with details."""
    lines = []

    if matches.empty:
        return f"No matches found for '{query}'"

    # Check if exact match
    is_exact = any(score == 1.0 for score in scores)

    if is_exact:
        lines.append(f"Found exact match(es) for '{query}':\n")
    else:
        lines.append(f"No exact match for '{query}'. Here are the closest matches:\n")

    for idx, (_, row) in enumerate(matches.iterrows(), 1):
        score = scores[idx - 1] if idx - 1 < len(scores) else 0.0
        match_quality = ""
        if score == 1.0:
            match_quality = " [EXACT]"
        elif score > 0.8:
            match_quality = f" [similarity: {score:.0%}]"
        elif score > 0.6:
            match_quality = f" [similarity: {score:.0%}]"
        elif score > 0.4:
            match_quality = f" [weak match: {score:.0%}]"

        lines.append(f"{idx}. {row['node_name']} (type: {row['node_type']}, node_index: {row['node_index']}){match_quality}")

        if show_details:
            # Show node_id if available
            if pd.notna(row.get('node_id')):
                lines.append(f"   ID: {row['node_id']}")

            # Show synonyms if available
            if 'node_synonyms' in row and pd.notna(row['node_synonyms']):
                synonyms = str(row['node_synonyms'])
                if synonyms and synonyms.lower() not in ['nan', 'none', '']:
                    # Truncate long synonym lists
                    if len(synonyms) > 100:
                        synonyms = synonyms[:100] + "..."
                    lines.append(f"   Aliases: {synonyms}")

        lines.append("")  # Blank line between matches

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Entity exploration tool with fuzzy matching",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Find exact or similar matches
  python match_entity.py "IL11" --type "gene/protein"

  # Show detailed information including aliases
  python match_entity.py "Crohn disease" --type disease --show-details

  # Get more matches
  python match_entity.py "gremlin" --limit 5

  # JSON output for programmatic use
  python match_entity.py "TP53" --format json
        """
    )
    parser.add_argument("entity", help="Entity name to search for")
    parser.add_argument("--type", choices=ENTITY_TYPES, help="Filter by entity type")
    parser.add_argument("--limit", type=int, default=5, help="Maximum number of matches to return (default: 5)")
    parser.add_argument("--show-details", action="store_true", help="Show node IDs and aliases")
    parser.add_argument("--format", choices=["text", "json"], default="text", help="Output format")
    args = parser.parse_args()

    nodes = load_nodes()
    matches, scores = fuzzy_search(nodes, args.entity, args.type, args.limit)

    if matches.empty:
        print(f"No matches found for '{args.entity}'", file=sys.stderr)
        if args.type:
            print(f"Try searching without --type filter or with a different entity type.", file=sys.stderr)
        sys.exit(1)

    if args.format == "json":
        # JSON output for programmatic use
        results = []
        for idx, (_, row) in enumerate(matches.iterrows()):
            result = {
                "rank": idx + 1,
                "node_name": row["node_name"],
                "node_type": row["node_type"],
                "node_index": int(row["node_index"]),
                "similarity_score": float(scores[idx]) if idx < len(scores) else 0.0
            }
            if pd.notna(row.get('node_id')):
                result["node_id"] = row["node_id"]
            if args.show_details and 'node_synonyms' in row and pd.notna(row['node_synonyms']):
                result["aliases"] = str(row['node_synonyms'])
            results.append(result)

        output = {
            "query": args.entity,
            "entity_type": args.type,
            "total_matches": len(results),
            "matches": results
        }
        print(json.dumps(output, indent=2))
    else:
        # Human-readable text output
        output_text = format_match_text(matches, scores, args.entity, args.show_details)
        print(output_text)

        # Add helpful message
        if not any(score == 1.0 for score in scores):
            print("\nWhich entity did you mean? Please:")
            print("- Reply with the number (e.g., '1' for the first match)")
            print("- Provide a more precise entity name")
            print("- Try a different search with --type filter")


if __name__ == "__main__":
    main()
