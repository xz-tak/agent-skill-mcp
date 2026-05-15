#!/usr/bin/env python3
"""
Parquet Utility Functions for KG Association Analysis
Helper functions for score extraction and analysis.

Usage:
    # Import as module
    from parquet_utils import get_ultra_disease_score, inspect_parquet, generate_score_matrix

    # Or run inspection from command line
    python parquet_utils.py inspect path/to/file.parquet
    python parquet_utils.py search path/to/file.parquet "crohn"
"""

import pandas as pd
import json
import sys
from pathlib import Path
from typing import List, Dict, Optional, Union


# =============================================================================
# ULTRA Parquet Column Reference
# =============================================================================

ULTRA_INDIVIDUAL_COLUMNS = {
    'h_label': 'Head entity ID (e.g., "NCBI:7297")',
    'h_name': 'Head entity name (e.g., "TYK2")',
    'h_type': 'Head entity type (e.g., "gene/protein")',
    'r_label': 'Relation label (e.g., "associated_with")',
    't_pred_label': 'Tail entity ID (e.g., "MONDO:5011")',
    't_pred_name': 'Tail entity name (e.g., "Crohn disease")',
    't_pred_score': 'Model prediction score (float64)',
    't_pred_type': 'Tail entity type (e.g., "disease")',
    'edge_in_primekg': 'Whether edge exists in PrimeKG (bool)',
    'rank': 'Global rank across ALL entity types (int64)',
    'percentile_rank': 'Global percentile across ALL types (float64)'
}

ULTRA_COMBO_COLUMNS = {
    'entity_id': 'Entity ID (e.g., "MONDO:5011")',
    'entity_name': 'Entity name (e.g., "Crohn disease")',
    'entity_type': 'Entity type (e.g., "disease")',
    'score': 'Model prediction score (float64)',
    'rank': 'Global rank',
    'percentile_rank': 'Global percentile',
    'schema_match': 'Whether entity matches expected schema (bool)'
}


# =============================================================================
# Quick Score Lookup Functions
# =============================================================================

def get_ultra_disease_score(
    parquet_path: str,
    disease_pattern: str,
    entity_type: str = "disease"
) -> Optional[dict]:
    """
    Quick lookup of a single disease score from ULTRA parquet.

    Args:
        parquet_path: Path to parquet file
        disease_pattern: Disease name pattern (case-insensitive partial match)
        entity_type: Filter to this entity type

    Returns:
        dict with pct_rank, rank, entity_name, etc. or None if not found

    Example:
        >>> score = get_ultra_disease_score("TYK2_predictions.parquet", "crohn")
        >>> print(f"TYK2-Crohn: pct_rank={score['pct_rank']:.4f}")
    """
    df = pd.read_parquet(parquet_path)

    # Detect column names
    if 't_pred_type' in df.columns:
        type_col, name_col, id_col, score_col = 't_pred_type', 't_pred_name', 't_pred_label', 't_pred_score'
    else:
        type_col, name_col, id_col, score_col = 'entity_type', 'entity_name', 'entity_id', 'score'

    df_filtered = df[df[type_col] == entity_type].copy()

    if len(df_filtered) == 0:
        return None

    # Recalculate within type
    df_filtered = df_filtered.sort_values(score_col, ascending=False).reset_index(drop=True)
    df_filtered['rank'] = range(1, len(df_filtered) + 1)
    total = len(df_filtered)
    df_filtered['pct_rank'] = 1.0 - (df_filtered['rank'] / total)

    # Find match
    match = df_filtered[df_filtered[name_col].str.lower().str.contains(disease_pattern.lower())]

    if len(match) == 0:
        return None

    row = match.iloc[0]
    return {
        'entity_id': row[id_col],
        'entity_name': row[name_col],
        'raw_score': float(row[score_col]),
        'rank': int(row['rank']),
        'pct_rank': float(row['pct_rank']),
        'total': total
    }


def get_all_disease_scores(
    parquet_path: str,
    disease_patterns: List[str],
    entity_type: str = "disease"
) -> Dict[str, Optional[dict]]:
    """
    Get scores for multiple diseases from a single parquet file.

    Args:
        parquet_path: Path to parquet file
        disease_patterns: List of disease patterns to search
        entity_type: Filter to this entity type

    Returns:
        Dictionary mapping disease pattern to score data
    """
    return {
        disease: get_ultra_disease_score(parquet_path, disease, entity_type)
        for disease in disease_patterns
    }


# =============================================================================
# Inspection Functions
# =============================================================================

def inspect_parquet(parquet_path: str, n_rows: int = 5) -> dict:
    """
    Inspect a parquet file's structure and content.

    Args:
        parquet_path: Path to parquet file
        n_rows: Number of rows to preview

    Returns:
        Dictionary with columns, dtypes, shape, and sample data
    """
    df = pd.read_parquet(parquet_path)

    info = {
        'path': parquet_path,
        'shape': df.shape,
        'columns': df.columns.tolist(),
        'dtypes': {col: str(dtype) for col, dtype in df.dtypes.items()},
        'sample': df.head(n_rows).to_dict('records')
    }

    # Detect file type
    if 't_pred_type' in df.columns:
        info['file_type'] = 'ULTRA Individual (predict_tail_entities)'
        info['entity_type_counts'] = df['t_pred_type'].value_counts().to_dict()
    elif 'entity_type' in df.columns:
        info['file_type'] = 'ULTRA Combo (answer_complex_query)'
        info['entity_type_counts'] = df['entity_type'].value_counts().to_dict()
    else:
        info['file_type'] = 'Unknown'

    return info


def search_parquet(
    parquet_path: str,
    pattern: str,
    n_results: int = 10
) -> pd.DataFrame:
    """
    Search for entities matching a pattern in parquet file.

    Args:
        parquet_path: Path to parquet file
        pattern: Search pattern (case-insensitive)
        n_results: Max results to return

    Returns:
        DataFrame with matching rows
    """
    df = pd.read_parquet(parquet_path)

    # Detect name column
    name_col = 't_pred_name' if 't_pred_name' in df.columns else 'entity_name'

    matches = df[df[name_col].str.lower().str.contains(pattern.lower())]
    return matches.head(n_results)


# =============================================================================
# Matrix Generation Functions
# =============================================================================

def generate_score_matrix(
    data_dir: str,
    genes: List[str],
    diseases: List[str],
    entity_type: str = "disease"
) -> pd.DataFrame:
    """
    Generate a gene x disease matrix of pct_rank scores.

    Args:
        data_dir: Directory containing <gene>_predictions.parquet files
        genes: List of gene names
        diseases: List of disease patterns
        entity_type: Entity type to filter to

    Returns:
        DataFrame with genes as rows and diseases as columns
    """
    matrix = []
    for gene in genes:
        pq_path = Path(data_dir) / f"{gene}_predictions.parquet"
        if not pq_path.exists():
            print(f"Warning: {pq_path} not found")
            row = {'Gene': gene}
            row.update({d: None for d in diseases})
        else:
            scores = get_all_disease_scores(str(pq_path), diseases, entity_type)
            row = {'Gene': gene}
            for disease in diseases:
                if scores[disease]:
                    row[disease] = scores[disease]['pct_rank']
                else:
                    row[disease] = None
        matrix.append(row)

    return pd.DataFrame(matrix).set_index('Gene')


def generate_combo_matrix(
    data_dir: str,
    combos: List[str],
    diseases: List[str],
    entity_type: str = "disease"
) -> pd.DataFrame:
    """
    Generate a combo x disease matrix of pct_rank scores.

    Args:
        data_dir: Directory containing combo parquet files
        combos: List of combo names (matching filenames without _predictions.parquet)
        diseases: List of disease patterns
        entity_type: Entity type to filter to

    Returns:
        DataFrame with combos as rows and diseases as columns
    """
    matrix = []
    data_path = Path(data_dir)

    for combo in combos:
        # Try different filename patterns
        possible_files = [
            data_path / f"{combo}_predictions_filtered.parquet",
            data_path / f"{combo}_filtered.parquet",
            data_path / f"{combo}_predictions.parquet",
            data_path / f"{combo}.parquet"
        ]

        pq_path = None
        for f in possible_files:
            if f.exists():
                pq_path = f
                break

        if pq_path is None:
            print(f"Warning: No file found for {combo}")
            row = {'Combo': combo}
            row.update({d: None for d in diseases})
        else:
            scores = get_all_disease_scores(str(pq_path), diseases, entity_type)
            row = {'Combo': combo}
            for disease in diseases:
                if scores[disease]:
                    row[disease] = scores[disease]['pct_rank']
                else:
                    row[disease] = None
        matrix.append(row)

    return pd.DataFrame(matrix).set_index('Combo')


# =============================================================================
# Synergy/Dilution Calculation
# =============================================================================

def calculate_synergy(
    combo_score: float,
    individual_scores: Optional[List[float]] = None,
    individual_ranks: Optional[List[Optional[float]]] = None,
    total_entities: Optional[int] = None,
    threshold: float = 0.02
) -> Dict[str, Union[float, str, None]]:
    """
    Compute synergy/dilution.
    Preferred: pass individual_ranks + total_entities for geometric mean.
    Fallback: pass individual_scores for legacy arithmetic mean.
    """
    import numpy as np

    if individual_ranks is not None and total_entities is not None:
        ranks = [r for r in individual_ranks if r is not None]
        if len(ranks) != len(individual_ranks):
            return {
                'combo_score': combo_score, 'baseline': None,
                'geo_rank_mean': None, 'delta': None, 'classification': 'INCOMPLETE'
            }
        geo_rank_mean = np.exp(np.mean(np.log(np.array(ranks, dtype=np.float64))))
        baseline = 1.0 - (geo_rank_mean - 1) / total_entities
    else:
        baseline = sum(individual_scores) / len(individual_scores)
        geo_rank_mean = None

    delta = combo_score - baseline

    if delta > threshold:
        classification = "SYNERGY"
    elif delta < -threshold:
        classification = "DILUTION"
    else:
        classification = "NEAR-ADDITIVE"

    return {
        'combo_score': round(combo_score, 6),
        'baseline': round(baseline, 6),
        'geo_rank_mean': round(geo_rank_mean, 2) if geo_rank_mean else None,
        'delta': round(delta, 6),
        'classification': classification
    }


# =============================================================================
# Command Line Interface
# =============================================================================

def main():
    """Command line interface for parquet utilities."""
    if len(sys.argv) < 3:
        print("Usage:")
        print("  python parquet_utils.py inspect <parquet_file>")
        print("  python parquet_utils.py search <parquet_file> <pattern>")
        print("  python parquet_utils.py score <parquet_file> <disease_pattern>")
        sys.exit(1)

    command = sys.argv[1]
    parquet_path = sys.argv[2]

    if command == "inspect":
        info = inspect_parquet(parquet_path)
        print(f"\nFile: {info['path']}")
        print(f"Shape: {info['shape']}")
        print(f"Type: {info['file_type']}")
        print(f"\nColumns: {info['columns']}")
        if 'entity_type_counts' in info:
            print(f"\nEntity type counts:")
            for k, v in info['entity_type_counts'].items():
                print(f"  {k}: {v}")
        print(f"\nFirst {len(info['sample'])} rows:")
        for row in info['sample']:
            print(f"  {row}")

    elif command == "search":
        if len(sys.argv) < 4:
            print("Usage: python parquet_utils.py search <parquet_file> <pattern>")
            sys.exit(1)
        pattern = sys.argv[3]
        results = search_parquet(parquet_path, pattern)
        print(f"\nResults for '{pattern}':")
        print(results.to_string())

    elif command == "score":
        if len(sys.argv) < 4:
            print("Usage: python parquet_utils.py score <parquet_file> <disease_pattern>")
            sys.exit(1)
        disease = sys.argv[3]
        score = get_ultra_disease_score(parquet_path, disease)
        if score:
            print(f"\nScore for '{disease}':")
            print(f"  Entity: {score['entity_name']} ({score['entity_id']})")
            print(f"  pct_rank: {score['pct_rank']:.4f}")
            print(f"  Rank: {score['rank']}/{score['total']}")
            print(f"  Raw score: {score['raw_score']:.4f}")
        else:
            print(f"\nNo match found for '{disease}'")

    else:
        print(f"Unknown command: {command}")
        sys.exit(1)


if __name__ == "__main__":
    main()
