#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
BioBridge Link Prediction Tool

This script uses the BioBridge MCP server's predict_associations function,
which includes LLM-based entity matching and neural link prediction.

Usage:
    python predict_link.py --head [ENTITY] --head-type [TYPE] --tail-type [TYPE] [OPTIONS]

Examples:
    # Predict diseases for a gene
    python predict_link.py --head GREM1 --head-type "gene/protein" --tail-type disease

    # Validate specific pair
    python predict_link.py --head IL11 --head-type "gene/protein" \
        --tail "Crohn disease" --tail-type disease

    # With export
    python predict_link.py --head TP53 --head-type "gene/protein" \
        --tail-type disease --export results.csv --format csv
"""

import argparse
import json
import os
import sys
from datetime import datetime
from typing import Any, Dict, List, Optional

# Add BioBridge source to path
BIOBRIDGE_SRC_DIR = "/home/sagemaker-user/biobridge/bbridge"
if BIOBRIDGE_SRC_DIR not in sys.path:
    sys.path.append(BIOBRIDGE_SRC_DIR)

# Import MCP server module
sys.path.insert(0, "/home/sagemaker-user/biobridge/bbridge/experiments/interpretation/agent_interpretation/tools")
import biobridge_mcp_server
from biobridge_mcp_server import predict_associations, _get_app

import pandas as pd


# Entity types
ENTITY_TYPES = [
    "biological_process", "molecular_function", "cellular_component",
    "gene/protein", "disease", "drug", "pathway", "effect/phenotype",
    "anatomy", "exposure", "biologics_drug",
]

# Common relations
RELATIONS = [
    "associated with", "treats", "interacts with", "side effect",
    "regulates", "participates in", "targets",
]


def set_custom_paths(
    kg_path: Optional[str] = None,
    nodes_path: Optional[str] = None,
    model_ckpt_path: Optional[str] = None,
    embedding_dir: Optional[str] = None
) -> None:
    """Override default paths in biobridge_mcp_server module."""
    if kg_path:
        biobridge_mcp_server.KG_CSV_PATH = kg_path
        print(f"Using custom KG path: {kg_path}", file=sys.stderr)

    if nodes_path:
        biobridge_mcp_server.NODES_CSV_PATH = nodes_path
        print(f"Using custom nodes path: {nodes_path}", file=sys.stderr)

    if model_ckpt_path:
        biobridge_mcp_server.CKPT_PATH = model_ckpt_path
        # Also update model config path (assume it's in same directory)
        ckpt_dir = os.path.dirname(model_ckpt_path)
        biobridge_mcp_server.MODEL_CONFIG_PATH = os.path.join(ckpt_dir, "model_config.json")
        print(f"Using custom model checkpoint: {model_ckpt_path}", file=sys.stderr)

    if embedding_dir:
        biobridge_mcp_server.EMB_SUBDIR = embedding_dir
        print(f"Using custom embedding directory: {embedding_dir}", file=sys.stderr)


def build_context(
    head: str,
    head_type: str,
    tail: Optional[str],
    tail_type: str,
    relation: Optional[str],
    custom_context: Optional[str]
) -> str:
    """Build context string for LLM-based entity matching."""
    if custom_context:
        return custom_context

    parts = []

    if tail:
        parts.append(f"Is {head} ({head_type})")
        if relation:
            parts.append(relation)
        parts.append(f"{tail} ({tail_type})?")
    else:
        if relation:
            parts.append(f"Find {tail_type} {relation} {head} ({head_type}).")
        else:
            parts.append(f"Find {tail_type} associated with {head} ({head_type}).")

    return " ".join(parts)


def print_resolved_mapping(result: Dict[str, Any]) -> None:
    """Print the resolved entity mapping for user confirmation."""
    resolved = result.get("resolved", {})

    print("\n" + "="*60)
    print("ENTITY MAPPING (from LLM-based matching)")
    print("="*60)

    print(f"\n**Head Entity:**")
    print(f"  Name: {resolved.get('head_name')}")
    print(f"  Type: {resolved.get('head_type')}")

    print(f"\n**Tail Entity:**")
    if resolved.get('tail_name'):
        print(f"  Name: {resolved.get('tail_name')}")
    else:
        print(f"  Type: {resolved.get('tail_type')} (exploring all entities)")

    print(f"\n**Relation:**")
    print(f"  {resolved.get('relation_family')}")

    chosen_rels = resolved.get('chosen_relations', [])
    if chosen_rels:
        print(f"  Relation IDs: {[r.get('id') for r in chosen_rels]}")

    print(f"\n**LLM Used:** {resolved.get('llm_used', False)}")

    print("\n" + "="*60 + "\n")


def deduplicate_results(results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Deduplicate results by entity, keeping highest score."""
    seen = {}
    for r in results:
        name = r.get("node_name") or r.get("mondo_name", "Unknown")
        node_index = r.get("node_index")
        key = (name, node_index)

        if key not in seen or r.get("cos_sim", 0) > seen[key].get("cos_sim", 0):
            seen[key] = r

    # Sort by score descending
    deduped = sorted(seen.values(), key=lambda x: x.get("cos_sim", 0), reverse=True)
    return deduped


def print_results(result: Dict[str, Any], topk: int = 20) -> None:
    """Print prediction results in readable format."""
    results = result.get("results", [])

    if not results:
        print("No predictions found.")
        return

    # Deduplicate results
    deduped_results = deduplicate_results(results)
    result["results"] = deduped_results  # Update for export

    print(f"\nTop {min(len(deduped_results), topk)} Predictions:")
    print("-" * 80)
    print(f"{'#':<4} {'Entity':<30} {'Cosine Sim':<12} {'Percentile':<12} {'ID':<10}")
    print("-" * 80)

    for idx, r in enumerate(deduped_results[:topk], 1):
        name = r.get("node_name") or r.get("mondo_name", "Unknown")
        cos_sim = r.get("cos_sim", 0.0)
        pct_rank = r.get("pct_rank", 0.0)
        entity_id = r.get("mondo_id") or r.get("node_id", "N/A")

        print(f"{idx:<4} {name[:28]:<30} {cos_sim:<12.4f} {pct_rank:<12.3f} {entity_id:<10}")

    print("-" * 80)
    print(f"\nTotal unique predictions: {len(deduped_results)}")


def export_results(
    result: Dict[str, Any],
    filepath: str,
    format: str = "csv"
) -> None:
    """Export results to file."""
    results = result.get("results", [])
    resolved = result.get("resolved", {})

    if not results:
        print("No results to export.", file=sys.stderr)
        return

    if format == "json":
        # Export as JSON
        export_data = {
            "query": resolved,
            "results": results,
            "metadata": {
                "timestamp": datetime.now().isoformat(),
                "total_results": len(results)
            }
        }
        with open(filepath, 'w') as f:
            json.dump(export_data, f, indent=2)
        print(f"Exported {len(results)} results to {filepath} (JSON)")

    elif format == "csv":
        # Export as CSV
        rows = []
        for r in results:
            rows.append({
                "Head_Entity": resolved.get("head_name"),
                "Head_Type": resolved.get("head_type"),
                "Tail_Entity": r.get("node_name") or r.get("mondo_name"),
                "Tail_Type": resolved.get("tail_type"),
                "Relation": resolved.get("relation_family"),
                "Cosine_Similarity": r.get("cos_sim"),
                "Percentile_Rank": r.get("pct_rank"),
                "Node_Index": r.get("node_index"),
                "Entity_ID": r.get("mondo_id") or r.get("node_id"),
            })

        df = pd.DataFrame(rows)
        df.to_csv(filepath, index=False)
        print(f"Exported {len(results)} results to {filepath} (CSV)")

    elif format == "excel":
        # Export as Excel with multiple sheets
        with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
            # Results sheet
            rows = []
            for r in results:
                rows.append({
                    "Head_Entity": resolved.get("head_name"),
                    "Head_Type": resolved.get("head_type"),
                    "Tail_Entity": r.get("node_name") or r.get("mondo_name"),
                    "Tail_Type": resolved.get("tail_type"),
                    "Relation": resolved.get("relation_family"),
                    "Cosine_Similarity": r.get("cos_sim"),
                    "Percentile_Rank": r.get("pct_rank"),
                    "Node_Index": r.get("node_index"),
                    "Entity_ID": r.get("mondo_id") or r.get("node_id"),
                })
            df = pd.DataFrame(rows)
            df.to_excel(writer, sheet_name="Results", index=False)

            # Summary sheet
            summary_data = {
                "Query Parameter": ["Head Entity", "Head Type", "Tail Type", "Relation", "Total Results"],
                "Value": [
                    resolved.get("head_name"),
                    resolved.get("head_type"),
                    resolved.get("tail_type"),
                    resolved.get("relation_family"),
                    len(results)
                ]
            }
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_excel(writer, sheet_name="Summary", index=False)

        print(f"Exported {len(results)} results to {filepath} (Excel)")


def main():
    parser = argparse.ArgumentParser(
        description="Predict biomedical entity associations using BioBridge",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Predict diseases for a gene
  python predict_link.py --head GREM1 --head-type "gene/protein" --tail-type disease

  # Validate specific pair
  python predict_link.py --head IL11 --head-type "gene/protein" \\
      --tail "Crohn disease" --tail-type disease --topk 1

  # Export results
  python predict_link.py --head TP53 --head-type "gene/protein" \\
      --tail-type disease --export tp53_diseases.csv --format csv
        """
    )

    # Required arguments
    parser.add_argument("--head", required=True, help="Head entity name (e.g., 'GREM1', 'TP53')")
    parser.add_argument("--head-type", required=True, choices=ENTITY_TYPES, help="Head entity type")
    parser.add_argument("--tail-type", required=True, choices=ENTITY_TYPES, help="Tail entity type to predict")

    # Optional arguments
    parser.add_argument("--tail", help="Specific tail entity name (for pair validation)")
    parser.add_argument("--relation", choices=RELATIONS, help="Relation type hint")
    parser.add_argument("--context", help="Custom context string for LLM matching")
    parser.add_argument("--topk", type=int, default=25, help="Number of top predictions (default: 25)")

    # Output options
    parser.add_argument("--export", help="Export results to file")
    parser.add_argument("--format", choices=["csv", "json", "excel"], default="csv", help="Export format")
    parser.add_argument("--show-all", action="store_true", help="Show all results (not just top 20)")

    # Advanced options
    parser.add_argument("--slidewindow", action="store_true", help="Use slidewindow embeddings")
    parser.add_argument("--no-slidewindow", action="store_true", help="Disable slidewindow embeddings")
    parser.add_argument("--debug", action="store_true", help="Show debug information")

    # Custom paths
    parser.add_argument("--kg-path", help="Custom path to knowledge graph CSV file")
    parser.add_argument("--nodes-path", help="Custom path to nodes CSV file")
    parser.add_argument("--model-ckpt", help="Custom path to model checkpoint (.bin file)")
    parser.add_argument("--embedding-dir", help="Custom path to embeddings directory")

    args = parser.parse_args()

    # Set custom paths if provided
    if any([args.kg_path, args.nodes_path, args.model_ckpt, args.embedding_dir]):
        print("Setting custom paths...", file=sys.stderr)
        set_custom_paths(
            kg_path=args.kg_path,
            nodes_path=args.nodes_path,
            model_ckpt_path=args.model_ckpt,
            embedding_dir=args.embedding_dir
        )

    # Build context for LLM-based entity matching
    context = build_context(
        args.head,
        args.head_type,
        args.tail,
        args.tail_type,
        args.relation,
        args.context
    )

    print(f"Query Context: {context}\n")
    print("Running prediction with LLM-based entity matching...")

    # Determine slidewindow setting
    slidewindow = None
    if args.slidewindow:
        slidewindow = True
    elif args.no_slidewindow:
        slidewindow = False

    # Call MCP predict_associations
    try:
        # Get app state (this initializes the MCP server state)
        app = _get_app()

        # Create a dummy context object (required by MCP tool signature)
        class DummyContext:
            pass

        ctx = DummyContext()

        # Call predict_associations
        result = predict_associations(
            ctx,
            context=context,
            topk=args.topk,
            override_head_name=args.head,
            override_head_type=args.head_type,
            override_tail_name=args.tail,
            override_tail_type=args.tail_type,
            relation_hint=args.relation,
            slidewindow=slidewindow,
            include_relation_catalog=False,
            include_debug=args.debug
        )

        # Check for errors
        if "error" in result:
            print(f"Error: {result['error']}", file=sys.stderr)
            sys.exit(1)

        # Print resolved mapping for user confirmation
        print_resolved_mapping(result)

        # Print results
        display_topk = None if args.show_all else 20
        print_results(result, topk=display_topk)

        # Export if requested
        if args.export:
            export_results(result, args.export, args.format)

        # Print debug info if requested
        if args.debug and "debug" in result:
            print("\nDebug Information:")
            print(json.dumps(result["debug"], indent=2))

    except Exception as e:
        print(f"Error during prediction: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
