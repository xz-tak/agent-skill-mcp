#!/usr/bin/env python3
"""
ULTRA Combo/Intersection Score Extraction
Extracts and recalculates pct_rank within disease type from filtered predictions.

UltraQuery Combo Parquet Column Schema (may vary):
    entity_id       - Entity ID (e.g., "MONDO:5011")
    entity_name     - Entity name (e.g., "Crohn disease")
    entity_type     - Entity type (e.g., "disease")
    score           - Model prediction score (float64)
    rank            - Global rank
    percentile_rank - Global percentile
    schema_match    - Whether entity matches expected schema

IMPORTANT: Native rank/percentile_rank are across ALL entity types.
           We MUST recalculate within target entity type (e.g., disease only).

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

import pandas as pd
import json
import argparse
from pathlib import Path
from typing import List, Dict, Optional


def extract_combo_scores(
    parquet_path: str,
    disease_list: List[str],
    entity_type: str = "disease"
) -> Dict[str, Optional[dict]]:
    """
    Extract disease-specific scores from ULTRA combo parquet file.

    Note: Combo files may use different column names than individual files:
        - 'entity_name' vs 't_pred_name'
        - 'entity_id' vs 't_pred_label'
        - 'entity_type' vs 't_pred_type'
        - 'score' vs 't_pred_score'

    Args:
        parquet_path: Path to the filtered parquet file
        disease_list: List of disease name patterns to search for
        entity_type: Entity type to filter to (default: "disease")

    Returns:
        Dictionary with disease scores
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

    # Filter to target entity type only
    df_filtered = df[df[type_col] == entity_type].copy()

    if len(df_filtered) == 0:
        return {disease: None for disease in disease_list}

    # Recalculate rank within entity type
    df_filtered = df_filtered.sort_values(score_col, ascending=False).reset_index(drop=True)
    df_filtered['entity_rank'] = range(1, len(df_filtered) + 1)
    total_entities = len(df_filtered)
    df_filtered['pct_rank'] = 1.0 - (df_filtered['entity_rank'] / total_entities)

    results = {}
    for disease in disease_list:
        match = df_filtered[df_filtered[name_col].str.lower().str.contains(disease.lower())]
        if len(match) > 0:
            row = match.iloc[0]
            results[disease] = {
                'entity_id': row[id_col],
                'entity_name': row[name_col],
                'score': float(row[score_col]),
                'entity_rank': int(row['entity_rank']),
                'pct_rank': float(row['pct_rank']),
                'total_entities': total_entities
            }
        else:
            results[disease] = None

    return results


def process_all_combos(
    data_dir: str,
    diseases: List[str],
    output_path: Optional[str] = None,
    entity_type: str = "disease"
) -> Dict[str, dict]:
    """
    Process all combo parquet files in a directory.

    Args:
        data_dir: Directory containing combo parquet files
        diseases: List of disease patterns to search for
        output_path: Optional path to save JSON results
        entity_type: Entity type to filter to

    Returns:
        Dictionary mapping combo name to disease scores
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
        # Extract combo name from filename
        combo_name = pq_file.stem
        for suffix in ["_predictions_filtered", "_filtered", "_predictions", "_all"]:
            combo_name = combo_name.replace(suffix, "")

        print(f"Processing {combo_name}...")
        results[combo_name] = extract_combo_scores(str(pq_file), diseases, entity_type)

    # Print summary
    print("\n" + "=" * 80)
    print(f"ULTRA Combo Intersection Scores (within {entity_type} type only)")
    print("=" * 80)
    for combo, combo_results in results.items():
        print(f"\n{combo}:")
        for disease, data in combo_results.items():
            if data:
                print(f"  {disease}: pct_rank={data['pct_rank']:.4f}, "
                      f"rank={data['entity_rank']}/{data['total_entities']}")
            else:
                print(f"  {disease}: NOT FOUND")

    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to {output_path}")

    return results


def calculate_synergy(
    combo_score: float,
    individual_scores: List[float]
) -> Dict[str, any]:
    """
    Calculate synergy/dilution for a combo vs individual components.

    Args:
        combo_score: The combo's pct_rank score
        individual_scores: List of individual component pct_rank scores

    Returns:
        Dictionary with combo_score, individual_mean, delta, classification
    """
    individual_mean = sum(individual_scores) / len(individual_scores)
    delta = combo_score - individual_mean

    if delta > 0.02:
        classification = "SYNERGY"
    elif delta < -0.02:
        classification = "DILUTION"
    else:
        classification = "NEAR-ADDITIVE"

    return {
        'combo_score': combo_score,
        'individual_mean': individual_mean,
        'delta': delta,
        'classification': classification
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract ULTRA combo scores from parquet files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python extract_ultra_combo.py --data-dir data/ultra/combo \\
                                  --diseases "crohn,ulcerative colitis,IBD"

    # With output file
    python extract_ultra_combo.py --data-dir data/ultra/combo \\
                                  --diseases "crohn,ulcerative colitis,IBD" \\
                                  --output results/combo_scores.json
        """
    )
    parser.add_argument("--data-dir", required=True, help="Directory with parquet files")
    parser.add_argument("--diseases", required=True, help="Comma-separated disease patterns")
    parser.add_argument("--output", help="Output JSON file path")
    parser.add_argument("--entity-type", default="disease", help="Entity type to filter (default: disease)")

    args = parser.parse_args()
    diseases = [d.strip() for d in args.diseases.split(",")]

    process_all_combos(args.data_dir, diseases, args.output, args.entity_type)
