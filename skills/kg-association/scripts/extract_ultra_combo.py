#!/usr/bin/env python3
"""
ULTRA Combo/Intersection Score Extraction
Extracts scores using the same pipeline as ULTRA individual inference:
  1. Filter to entity_type == "disease"
  2. Compute probability = sigmoid(score) = exp(score) / (1 + exp(score))
  3. Compute pct_rank = rank(probability, ascending=True, ordinal) / total

UltraQuery Combo Parquet Column Schema (may vary):
    entity_id       - Entity ID (e.g., "MONDO:5011")
    entity_name     - Entity name (e.g., "Crohn disease")
    entity_type     - Entity type (e.g., "disease")
    score           - Model prediction score (raw logit, float64)
    rank            - Global rank
    percentile_rank - Global percentile
    schema_match    - Whether entity matches expected schema

IMPORTANT: Native rank/percentile_rank are across ALL entity types.
           We MUST recalculate within target entity type (e.g., disease only).
           Method matches ultra_inference_server.py lines 926-936 exactly.

Usage:
    python extract_ultra_combo.py --data-dir ./kgpred_IBD_2025-01-05/data/ultra/combo \\
                                  --diseases "crohn,ulcerative colitis,inflammatory bowel disease" \\
                                  --output ./kgpred_IBD_2025-01-05/data/ultra/combo_scores.json

Example:
    # From working directory with parquet files
    python extract_ultra_combo.py --data-dir data/ultra/combo \\
                                  --diseases "crohn,ulcerative colitis,IBD" \\
                                  --output data/ultra/combo_scores.json
"""

import numpy as np
import pandas as pd
import json
import argparse
from pathlib import Path
from typing import List, Dict, Optional


def infer_tail_type(df: pd.DataFrame, target_names: List[str],
                    name_col: str, type_col: str) -> Optional[str]:
    """
    Infer the tail_type by finding which entity_type contains the target names.

    Searches each entity_type partition for matches against target_names.
    Returns the type with the most matches. Returns None if no matches found
    in any type (caller should ask user or include all types).
    """
    available_types = df[type_col].unique().tolist()

    if len(available_types) == 1:
        return available_types[0]

    type_match_counts = {}
    for etype in available_types:
        subset = df[df[type_col] == etype]
        matches = 0
        for name in target_names:
            if subset[name_col].str.lower().str.contains(name.lower()).any():
                matches += 1
        type_match_counts[etype] = matches

    best_type = max(type_match_counts, key=type_match_counts.get)
    if type_match_counts[best_type] > 0:
        return best_type

    # Cannot infer — return None to signal caller
    return None


def extract_combo_scores(
    parquet_path: str,
    target_names: List[str],
    tail_type: Optional[str] = None
) -> Dict[str, Optional[dict]]:
    """
    Extract entity-specific scores from ULTRA combo parquet file.

    Matches the ULTRA individual inference pipeline (ultra_inference_server.py:926-936):
      1. Filter to tail_type — inferred from target_names if not provided
      2. Compute probability = sigmoid(score) = exp(score) / (1 + exp(score))
      3. Compute pct_rank = rank(probability, ascending=True, ordinal) / total
         → highest probability gets pct_rank ≈ 1.0

    Note: Combo files may use different column names than individual files:
        - 'entity_name' vs 't_pred_name'
        - 'entity_id' vs 't_pred_label'
        - 'entity_type' vs 't_pred_type'
        - 'score' vs 't_pred_score'

    Args:
        parquet_path: Path to the filtered parquet file
        target_names: List of entity name patterns to search for (case-insensitive)
        tail_type: Entity type to filter to. If None, inferred from target_names.

    Returns:
        Dictionary mapping each target name to its score data or None
    """
    df = pd.read_parquet(parquet_path)

    # Detect column names (combo files may differ from individual files)
    if 'entity_name' in df.columns:
        name_col = 'entity_name'
        id_col = 'entity_id'
        type_col = 'entity_type'
        score_col = 'score'
    else:
        name_col = 't_pred_name'
        id_col = 't_pred_label'
        type_col = 't_pred_type'
        score_col = 't_pred_score'

    # Step 1: Determine tail_type — infer from target_names if not explicit
    if tail_type is None:
        tail_type = infer_tail_type(df, target_names, name_col, type_col)

    # If inference failed (None), include all entity types
    if tail_type is None:
        tail_type = "__all__"
        df_filtered = df.copy()
    else:
        df_filtered = df[df[type_col] == tail_type].copy()

    if len(df_filtered) == 0:
        return {name: None for name in target_names}

    # Step 2: Compute probability = sigmoid(score)
    scores = df_filtered[score_col].values
    df_filtered['probability'] = np.where(
        scores >= 0,
        1.0 / (1.0 + np.exp(-scores)),
        np.exp(scores) / (1.0 + np.exp(scores))
    )

    # Step 3: pct_rank = rank(probability, ascending=True, ordinal) / total
    total_entities = len(df_filtered)
    df_filtered['pct_rank'] = df_filtered['probability'].rank(
        ascending=True, method='first'
    ) / total_entities

    # entity_rank: 1 = highest probability
    df_filtered['entity_rank'] = df_filtered['probability'].rank(
        ascending=False, method='first'
    ).astype(int)

    results = {}
    for name in target_names:
        match = df_filtered[df_filtered[name_col].str.lower().str.contains(name.lower())]
        if len(match) > 0:
            row = match.iloc[0]
            results[name] = {
                'entity_id': row[id_col],
                'entity_name': row[name_col],
                'score': float(row[score_col]),
                'probability': float(row['probability']),
                'entity_rank': int(row['entity_rank']),
                'pct_rank': float(row['pct_rank']),
                'total_entities': total_entities,
                'tail_type': tail_type
            }
        else:
            results[name] = None

    return results


def process_all_combos(
    data_dir: str,
    target_names: List[str],
    output_path: Optional[str] = None,
    tail_type: Optional[str] = None
) -> Dict[str, dict]:
    """
    Process all combo parquet files in a directory.

    Args:
        data_dir: Directory containing combo parquet files
        target_names: List of entity name patterns to search for
        output_path: Optional path to save JSON results
        tail_type: Entity type to filter to. If None, inferred per file from target_names.

    Returns:
        Dictionary mapping combo name to entity scores
    """
    data_path = Path(data_dir)

    # Look for filtered files first, then all predictions
    parquet_files = list(data_path.glob("*_filtered.parquet"))
    if not parquet_files:
        parquet_files = list(data_path.glob("*_predictions.parquet"))
    if not parquet_files:
        parquet_files = list(data_path.glob("*.parquet"))

    if not parquet_files:
        print(f"No parquet files found in {data_dir}")
        return {}

    results = {}
    for pq_file in sorted(parquet_files):
        combo_name = pq_file.stem
        for suffix in ["_predictions_filtered", "_filtered", "_predictions", "_all"]:
            combo_name = combo_name.replace(suffix, "")

        print(f"Processing {combo_name}...")
        results[combo_name] = extract_combo_scores(str(pq_file), target_names, tail_type)

    # Print summary
    print("\n" + "=" * 80)
    print("ULTRA Combo Intersection Scores (sigmoid→pct_rank pipeline)")
    print("=" * 80)
    for combo, combo_results in results.items():
        print(f"\n{combo}:")
        for target, data in combo_results.items():
            if data:
                print(f"  {target}: pct_rank={data['pct_rank']:.4f}, "
                      f"prob={data['probability']:.6f}, "
                      f"rank={data['entity_rank']}/{data['total_entities']} "
                      f"[{data['tail_type']}]")
            else:
                print(f"  {target}: NOT FOUND")

    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to {output_path}")

    return results


def calculate_synergy(
    combo_pct_rank: float,
    individual_ranks: List[Optional[float]],
    total_entities: int,
    threshold: float = 0.02
) -> Dict[str, any]:
    """
    Compute synergy using geometric mean of raw ranks (log-space for stability).

    Args:
        combo_pct_rank: The combo's pct_rank score
        individual_ranks: List of individual component raw ranks (1-indexed, 1=best)
        total_entities: Total number of entities in the ranking
        threshold: Delta threshold for classification (default: 0.02)

    Returns:
        Dictionary with combo_pct_rank, geo_rank_mean, geo_pct_rank, delta, classification
    """
    ranks = [r for r in individual_ranks if r is not None]
    if len(ranks) != len(individual_ranks):
        return {
            'combo_pct_rank': combo_pct_rank,
            'geo_rank_mean': None,
            'geo_pct_rank': None,
            'delta': None,
            'classification': 'INCOMPLETE'
        }

    geo_rank_mean = np.exp(np.mean(np.log(np.array(ranks, dtype=np.float64))))
    geo_pct_rank = 1.0 - (geo_rank_mean - 1) / total_entities
    delta = combo_pct_rank - geo_pct_rank

    if delta > threshold:
        classification = "SYNERGY"
    elif delta < -threshold:
        classification = "DILUTION"
    else:
        classification = "NEAR-ADDITIVE"

    return {
        'combo_pct_rank': round(combo_pct_rank, 6),
        'geo_rank_mean': round(geo_rank_mean, 2),
        'geo_pct_rank': round(geo_pct_rank, 6),
        'delta': round(delta, 6),
        'classification': classification
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract ULTRA combo scores from parquet files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Disease targets (tail_type inferred as "disease")
    python extract_ultra_combo.py --data-dir data/ultra/combo \\
                                  --targets "crohn,ulcerative colitis,IBD"

    # Gene targets (tail_type inferred as "gene/protein")
    python extract_ultra_combo.py --data-dir data/ultra/combo \\
                                  --targets "GREM1,IL23A,TYK2"

    # Explicit tail_type override
    python extract_ultra_combo.py --data-dir data/ultra/combo \\
                                  --targets "crohn,ulcerative colitis" \\
                                  --tail-type disease

    # With output file
    python extract_ultra_combo.py --data-dir data/ultra/combo \\
                                  --targets "crohn,ulcerative colitis,IBD" \\
                                  --output results/combo_scores.json
        """
    )
    parser.add_argument("--data-dir", required=True, help="Directory with parquet files")
    parser.add_argument("--targets", "--diseases", required=True,
                        help="Comma-separated entity name patterns to search for")
    parser.add_argument("--output", help="Output JSON file path")
    parser.add_argument("--tail-type", default=None,
                        help="Entity type to filter to. If omitted, inferred from targets.")

    args = parser.parse_args()
    targets = [t.strip() for t in args.targets.split(",")]

    process_all_combos(args.data_dir, targets, args.output, args.tail_type)
