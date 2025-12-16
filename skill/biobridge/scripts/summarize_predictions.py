#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
BioBridge Prediction Summary Generator

This script generates human-readable summaries from BioBridge prediction results.

Usage:
    python summarize_predictions.py [RESULTS_FILE] [OPTIONS]

Examples:
    python summarize_predictions.py results.csv
    python summarize_predictions.py results.json --threshold 0.1 --top 15
    python summarize_predictions.py results.csv --output summary.txt
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import pandas as pd


def load_results(filepath: str) -> Dict:
    """Load results from CSV or JSON file."""
    path = Path(filepath)

    if not path.exists():
        print(f"Error: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)

    if path.suffix == ".json":
        with open(filepath, 'r') as f:
            data = json.load(f)

        # Extract query and results
        if "query" in data and "results" in data:
            return data
        else:
            print("Error: JSON file missing 'query' or 'results' keys", file=sys.stderr)
            sys.exit(1)

    elif path.suffix == ".csv":
        df = pd.read_csv(filepath)

        # Reconstruct query and results from CSV
        if df.empty:
            print("Error: CSV file is empty", file=sys.stderr)
            sys.exit(1)

        # Query info from first row
        query = {
            "head_name": df.iloc[0]["Head_Entity"],
            "head_type": df.iloc[0]["Head_Type"],
            "tail_type": df.iloc[0]["Tail_Type"],
            "relation_family": df.iloc[0]["Relation"]
        }

        # Results
        results = []
        for _, row in df.iterrows():
            results.append({
                "node_name": row.get("Tail_Entity"),
                "cos_sim": row.get("Cosine_Similarity"),
                "pct_rank": row.get("Percentile_Rank"),
                "node_index": row.get("Node_Index"),
                "entity_id": row.get("Entity_ID")
            })

        return {"query": query, "results": results}

    else:
        print(f"Error: Unsupported file format: {path.suffix}", file=sys.stderr)
        print("Supported formats: .json, .csv", file=sys.stderr)
        sys.exit(1)


def interpret_score(cos_sim: float, pct_rank: float) -> str:
    """Interpret association scores."""
    if cos_sim >= 0.2 or pct_rank >= 0.9:
        return "Strong"
    elif cos_sim >= 0.1 or pct_rank >= 0.75:
        return "Moderate"
    elif cos_sim >= 0.05 or pct_rank >= 0.5:
        return "Weak"
    else:
        return "Very Weak"


def filter_by_threshold(results: List[Dict], threshold: float) -> List[Dict]:
    """Filter results by minimum cosine similarity."""
    return [r for r in results if r.get("cos_sim", 0) >= threshold]


def generate_summary(data: Dict, threshold: float, top_n: int) -> str:
    """Generate human-readable summary."""
    query = data["query"]
    results = data["results"]

    # Filter by threshold
    filtered_results = filter_by_threshold(results, threshold)

    # Sort by score
    sorted_results = sorted(filtered_results, key=lambda x: x.get("cos_sim", 0), reverse=True)

    # Build summary
    lines = []
    lines.append("=" * 80)
    lines.append("BIOBRIDGE LINK PREDICTION SUMMARY")
    lines.append("=" * 80)
    lines.append("")

    # Query information
    lines.append("QUERY:")
    lines.append(f"  Head Entity: {query.get('head_name')} ({query.get('head_type')})")
    if query.get('tail_name'):
        lines.append(f"  Tail Entity: {query.get('tail_name')} ({query.get('tail_type')})")
    else:
        lines.append(f"  Tail Type: {query.get('tail_type')} (exploring all)")
    lines.append(f"  Relation: {query.get('relation_family')}")
    lines.append("")

    # Results summary
    lines.append("RESULTS SUMMARY:")
    lines.append(f"  Total Predictions: {len(results)}")
    lines.append(f"  Above Threshold ({threshold}): {len(sorted_results)}")
    lines.append("")

    if sorted_results:
        # Score distribution
        scores = [r.get("cos_sim", 0) for r in sorted_results]
        lines.append("  Score Distribution:")
        lines.append(f"    Maximum: {max(scores):.4f}")
        lines.append(f"    Minimum: {min(scores):.4f}")
        lines.append(f"    Mean:    {sum(scores)/len(scores):.4f}")
        lines.append("")

        # Strength breakdown
        strong = len([r for r in sorted_results if interpret_score(r.get("cos_sim", 0), r.get("pct_rank", 0)) == "Strong"])
        moderate = len([r for r in sorted_results if interpret_score(r.get("cos_sim", 0), r.get("pct_rank", 0)) == "Moderate"])
        weak = len([r for r in sorted_results if interpret_score(r.get("cos_sim", 0), r.get("pct_rank", 0)) in ["Weak", "Very Weak"]])

        lines.append("  Evidence Strength:")
        lines.append(f"    Strong:   {strong} ({100*strong/len(sorted_results):.1f}%)")
        lines.append(f"    Moderate: {moderate} ({100*moderate/len(sorted_results):.1f}%)")
        lines.append(f"    Weak:     {weak} ({100*weak/len(sorted_results):.1f}%)")
        lines.append("")

    # Top predictions
    lines.append("=" * 80)
    lines.append(f"TOP {min(top_n, len(sorted_results))} PREDICTIONS:")
    lines.append("=" * 80)
    lines.append("")

    if not sorted_results:
        lines.append("No predictions above the specified threshold.")
    else:
        for idx, result in enumerate(sorted_results[:top_n], 1):
            name = result.get("node_name", "Unknown")
            cos_sim = result.get("cos_sim", 0.0)
            pct_rank = result.get("pct_rank", 0.0)
            strength = interpret_score(cos_sim, pct_rank)

            lines.append(f"{idx}. {name}")
            lines.append(f"   Cosine Similarity: {cos_sim:.4f}")
            lines.append(f"   Percentile Rank:   {pct_rank:.3f} ({pct_rank*100:.1f}%ile)")
            lines.append(f"   Evidence Strength: {strength}")
            lines.append("")

    # Interpretation
    lines.append("=" * 80)
    lines.append("INTERPRETATION:")
    lines.append("=" * 80)
    lines.append("")

    if not sorted_results:
        lines.append(f"No significant associations found above threshold {threshold}.")
        lines.append("This could indicate:")
        lines.append("  - The entities have minimal connection in the knowledge graph")
        lines.append("  - The relation type may not be appropriate")
        lines.append("  - Consider lowering the threshold or exploring different entity types")
    elif sorted_results and sorted_results[0].get("cos_sim", 0) > 0.2:
        lines.append("Strong associations detected:")
        lines.append(f"  The top prediction ({sorted_results[0].get('node_name')}) shows robust")
        lines.append(f"  evidence with a cosine similarity of {sorted_results[0].get('cos_sim', 0):.4f}.")
        lines.append("")
        lines.append(f"  This ranks in the {sorted_results[0].get('pct_rank', 0)*100:.1f}th percentile,")
        lines.append("  indicating it's among the strongest associations for this entity")
        lines.append("  in the knowledge graph.")
    else:
        lines.append("Moderate to weak associations detected:")
        lines.append("  The predictions show some evidence of connection, but scores are")
        lines.append("  relatively modest. These associations may be:")
        lines.append("    - Less well-established in the literature")
        lines.append("    - Indirect or mediated by other entities")
        lines.append("    - Worthy of further investigation")
        lines.append("")
        lines.append("  Consider comparing these scores to other entity pairs as a baseline.")

    lines.append("")
    lines.append("=" * 80)
    lines.append("Note: These predictions are based on neural retrieval over a knowledge")
    lines.append("graph and should be validated with literature or experimental evidence.")
    lines.append("=" * 80)

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Generate summary from BioBridge prediction results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python summarize_predictions.py results.csv
  python summarize_predictions.py results.json --threshold 0.1 --top 15
  python summarize_predictions.py results.csv --output summary.txt
        """
    )

    parser.add_argument("results_file", help="Input file (CSV or JSON from predict_link.py)")
    parser.add_argument("--threshold", type=float, default=0.05,
                        help="Minimum cosine similarity threshold (default: 0.05)")
    parser.add_argument("--top", type=int, default=10,
                        help="Number of top predictions to show (default: 10)")
    parser.add_argument("--output", help="Save summary to file (default: print to stdout)")

    args = parser.parse_args()

    # Load results
    print("Loading results...", file=sys.stderr)
    data = load_results(args.results_file)

    # Generate summary
    print("Generating summary...", file=sys.stderr)
    summary = generate_summary(data, args.threshold, args.top)

    # Output
    if args.output:
        with open(args.output, 'w') as f:
            f.write(summary)
        print(f"\nSummary saved to: {args.output}", file=sys.stderr)
    else:
        print()  # Blank line before summary
        print(summary)


if __name__ == "__main__":
    main()
