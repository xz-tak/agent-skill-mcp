#!/usr/bin/env python3
"""
ULTRA Individual Gene-Disease Score Extraction
Extracts and recalculates pct_rank within disease type only.

ULTRA Parquet Column Schema:
    h_label         - Head entity ID (e.g., "NCBI:7297")
    h_name          - Head entity name (e.g., "TYK2")
    h_type          - Head entity type (e.g., "gene/protein")
    r_label         - Relation label (e.g., "associated_with")
    t_pred_label    - Tail entity ID (e.g., "MONDO:5011")
    t_pred_name     - Tail entity name (e.g., "Crohn disease")
    t_pred_score    - Model prediction score (float64)
    t_pred_type     - Tail entity type (e.g., "disease")
    edge_in_primekg - Whether edge exists in PrimeKG (bool)
    rank            - Global rank across ALL entity types (int64)
    percentile_rank - Global percentile across ALL types (float64)

IMPORTANT: Native rank/percentile_rank are across ALL entity types.
           We MUST recalculate within target entity type (e.g., disease only).

Usage:
    python extract_ultra_individual.py --data-dir ./kgpred_IBD_2025-01-05/data/ultra/individual \\
                                       --diseases "crohn,ulcerative colitis,inflammatory bowel disease" \\
                                       --output ./kgpred_IBD_2025-01-05/data/ultra/individual_scores.json

Example:
    # From working directory with parquet files
    python extract_ultra_individual.py --data-dir data/ultra/individual \\
                                       --diseases "crohn,ulcerative colitis,IBD" \\
                                       --output data/ultra/scores.json
"""

import numpy as np
import pandas as pd
import json
import argparse
from pathlib import Path
from typing import List, Dict, Optional


def extract_disease_scores(
    parquet_path: str,
    disease_list: List[str],
    entity_type: str = "disease"
) -> Dict[str, Optional[dict]]:
    """
    Extract disease-specific scores from ULTRA parquet file.

    Args:
        parquet_path: Path to the parquet file
        disease_list: List of disease name patterns to search for (case-insensitive)
        entity_type: Entity type to filter to (default: "disease")

    Returns:
        Dictionary mapping disease pattern to score data or None if not found
    """
    df = pd.read_parquet(parquet_path)

    # Filter to target entity type only
    df_filtered = df[df['t_pred_type'] == entity_type].copy()

    if len(df_filtered) == 0:
        return {disease: None for disease in disease_list}

    # Handle both old schema (t_pred_score) and new schema (t_pred_logit + t_pred_probability)
    if 't_pred_probability' in df_filtered.columns:
        df_filtered['probability'] = df_filtered['t_pred_probability']
        score_col = 't_pred_logit' if 't_pred_logit' in df_filtered.columns else 't_pred_probability'
    else:
        score_col = 't_pred_score'
        scores = df_filtered[score_col].values
        df_filtered['probability'] = np.where(
            scores >= 0,
            1.0 / (1.0 + np.exp(-scores)),
            np.exp(scores) / (1.0 + np.exp(scores))
        )

    total_entities = len(df_filtered)
    df_filtered['pct_rank'] = df_filtered['probability'].rank(
        ascending=True, method='first'
    ) / total_entities

    # disease_rank: 1 = highest probability
    df_filtered['disease_rank'] = df_filtered['probability'].rank(
        ascending=False, method='first'
    ).astype(int)

    results = {}
    for disease in disease_list:
        # Case-insensitive partial match
        match = df_filtered[df_filtered['t_pred_name'].str.lower().str.contains(disease.lower())]
        if len(match) > 0:
            row = match.iloc[0]
            results[disease] = {
                'entity_id': row['t_pred_label'],
                'entity_name': row['t_pred_name'],
                'score': float(row[score_col]),
                'probability': float(row['probability']),
                'disease_rank': int(row['disease_rank']),
                'pct_rank': float(row['pct_rank']),
                'total_entities': total_entities
            }
        else:
            results[disease] = None

    return results


def process_all_genes(
    data_dir: str,
    diseases: List[str],
    output_path: Optional[str] = None,
    entity_type: str = "disease"
) -> Dict[str, dict]:
    """
    Process all gene parquet files in a directory.

    Args:
        data_dir: Directory containing <gene>_predictions.parquet files
        diseases: List of disease patterns to search for
        output_path: Optional path to save JSON results
        entity_type: Entity type to filter to

    Returns:
        Dictionary mapping gene name to disease scores
    """
    data_path = Path(data_dir)

    # Support both flat ({GENE}_predictions.parquet) and nested ({ID}/relation/predictions.parquet)
    parquet_files = list(data_path.glob("*_predictions.parquet"))
    nested_files = list(data_path.glob("*/*/predictions.parquet"))

    if not parquet_files and not nested_files:
        print(f"No parquet files found in {data_dir}")
        return {}

    results = {}
    for pq_file in sorted(parquet_files):
        gene_name = pq_file.stem.replace("_predictions", "")
        print(f"Processing {gene_name}...")
        results[gene_name] = extract_disease_scores(str(pq_file), diseases, entity_type)

    for pq_file in sorted(nested_files):
        # Read to get gene name from h_name column
        df_peek = pd.read_parquet(pq_file, columns=['h_name'])
        gene_name = df_peek['h_name'].iloc[0] if len(df_peek) > 0 else pq_file.parent.parent.name
        if gene_name not in results:
            print(f"Processing {gene_name} (nested)...")
            results[gene_name] = extract_disease_scores(str(pq_file), diseases, entity_type)

    # Print summary
    print("\n" + "=" * 80)
    print(f"ULTRA Individual Gene-Disease Scores (within {entity_type} type only)")
    print("=" * 80)
    for gene, gene_results in results.items():
        print(f"\n{gene}:")
        for disease, data in gene_results.items():
            if data:
                print(f"  {disease}: pct_rank={data['pct_rank']:.4f}, "
                      f"rank={data['disease_rank']}/{data['total_entities']}, "
                      f"entity={data['entity_name']}")
            else:
                print(f"  {disease}: NOT FOUND")

    # Save to JSON if output path provided
    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to {output_path}")

    return results


def generate_score_matrix(results: Dict[str, dict], diseases: List[str]) -> pd.DataFrame:
    """
    Generate a gene x disease matrix of pct_rank scores.

    Args:
        results: Output from process_all_genes()
        diseases: List of disease patterns

    Returns:
        DataFrame with genes as rows and diseases as columns
    """
    matrix = []
    for gene, gene_results in results.items():
        row = {'Gene': gene}
        for disease in diseases:
            if gene_results.get(disease):
                row[disease] = gene_results[disease]['pct_rank']
            else:
                row[disease] = None
        matrix.append(row)

    return pd.DataFrame(matrix).set_index('Gene')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract ULTRA disease scores from parquet files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python extract_ultra_individual.py --data-dir data/ultra/individual \\
                                       --diseases "crohn,ulcerative colitis,IBD"

    # With output file
    python extract_ultra_individual.py --data-dir data/ultra/individual \\
                                       --diseases "crohn,ulcerative colitis,IBD" \\
                                       --output results/scores.json

    # Different entity type
    python extract_ultra_individual.py --data-dir data/ultra/individual \\
                                       --diseases "schizophrenia,depression" \\
                                       --entity-type disease
        """
    )
    parser.add_argument("--data-dir", required=True, help="Directory with parquet files")
    parser.add_argument("--diseases", required=True, help="Comma-separated disease patterns")
    parser.add_argument("--output", help="Output JSON file path")
    parser.add_argument("--entity-type", default="disease", help="Entity type to filter (default: disease)")

    args = parser.parse_args()
    diseases = [d.strip() for d in args.diseases.split(",")]

    process_all_genes(args.data_dir, diseases, args.output, args.entity_type)
