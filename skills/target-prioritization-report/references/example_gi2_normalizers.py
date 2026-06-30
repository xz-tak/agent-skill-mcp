"""
Normalizer functions for target prioritization scoring.

These functions transform raw values into the 0-100 scoring scale.
"""

from typing import Dict, List, Optional, Any

# Clinical validation max (based on Cortellis data)
MAX_CLINICAL_SCORE = 134.5

# Safety risk category base scores
SAFETY_CATEGORY_SCORES = {
    'very_high': 0,
    'high': 10,
    'medium': 20,
    'low': 60,
    'very_low': 80,
    'not_assoc': 100,
}


def normalize_clinical(raw: float, max_score: float = MAX_CLINICAL_SCORE) -> float:
    """
    Normalize clinical validation score to 0-100.

    Args:
        raw: Raw Cortellis clinical score
        max_score: Maximum clinical score for normalization

    Returns:
        Normalized score (0-100)
    """
    if raw is None:
        return 50.0
    normalized = (float(raw) / max_score) * 100.0
    return max(0.0, min(100.0, normalized))


def normalize_deg(
    raw: float,
    max_score: Optional[float] = None,
    approach: str = "antagonist",
    all_scores: Optional[List[float]] = None,
) -> float:
    """
    Normalize DEG score to 0-100 with approach-based transformation.

    NEW BEHAVIOR (dynamic max, no clamping, approach support):
    1. Apply approach-based transformation:
       - antagonist (default): use raw score as-is
       - agonist: negate raw score (-raw)
    2. Normalize using dynamic max (computed from data, not hardcoded)
    3. No clamping - scores reflect actual distribution

    Args:
        raw: Raw DEG score (can be positive or negative)
        max_score: Dynamic max score for normalization. If None and all_scores
                   is provided, computes from all_scores. If both None, returns 50.0.
        approach: Therapeutic approach ("antagonist" or "agonist"). Default: "antagonist"
        all_scores: Optional list of all adjusted scores for dynamic max computation

    Returns:
        Normalized score (not clamped - can exceed 0-100 range)
    """
    if raw is None:
        return 50.0

    # Step 1: Apply approach-based transformation
    if approach.lower() == "agonist":
        adjusted_raw = -float(raw)
    else:
        adjusted_raw = float(raw)

    # Step 2: Determine max_score dynamically
    if max_score is None:
        if all_scores is not None and len(all_scores) > 0:
            max_score = max(abs(s) for s in all_scores)
        else:
            # No max available - return neutral
            return 50.0

    if max_score == 0:
        return 50.0

    # Step 3: Normalize WITHOUT clamping
    normalized = (adjusted_raw / max_score) * 100.0

    return normalized


def compute_deg_max_score(
    deg_scores: Dict[str, float],
    approaches: Optional[Dict[str, str]] = None,
) -> float:
    """
    Compute the dynamic max DEG score from a dataset.

    Args:
        deg_scores: Dict mapping gene names to raw DEG scores
        approaches: Optional dict mapping gene names to approach ("antagonist" or "agonist")
                    Default approach is "antagonist" if not specified.

    Returns:
        Maximum absolute value of adjusted scores
    """
    if not deg_scores:
        return 1.0  # Avoid division by zero

    approaches = approaches or {}

    adjusted_scores = []
    for gene, raw in deg_scores.items():
        if raw is None:
            continue
        approach = approaches.get(gene, "antagonist")
        if approach.lower() == "agonist":
            adjusted_scores.append(-float(raw))
        else:
            adjusted_scores.append(float(raw))

    if not adjusted_scores:
        return 1.0

    return max(abs(s) for s in adjusted_scores)


def normalize_biobridge(raw: float) -> float:
    """
    Normalize BioBridge percentile to 0-100.

    BioBridge percentiles are already in 0-1 range (pct_rank).
    Convert to 0-100 scale.

    Args:
        raw: BioBridge pct_rank (0-1)

    Returns:
        Normalized score (0-100)
    """
    if raw is None:
        return 50.0
    # pct_rank is 0-1, convert to 0-100
    normalized = float(raw) * 100.0
    return max(0.0, min(100.0, normalized))


def normalize_ultra(raw: float) -> float:
    """
    Normalize ULTRA percentile to 0-100.

    ULTRA percentiles should use pct_rank (0-1 range, recalculated
    within disease entities only).

    Args:
        raw: ULTRA pct_rank (0-1)

    Returns:
        Normalized score (0-100)
    """
    if raw is None:
        return 50.0
    # pct_rank is 0-1, convert to 0-100
    normalized = float(raw) * 100.0
    return max(0.0, min(100.0, normalized))


def normalize_primekg(
    raw: Optional[float] = None,
    max_connections: Optional[int] = None,
    path_length: Optional[int] = None,
    connectivity_score: Optional[float] = None,
) -> float:
    """
    Normalize PrimeKG connectivity to 0-100.

    Preferred inputs (in order):
    1. connectivity_score: Path-based score in [0, 1] (preferred)
    2. path_length: Integer path length >= 1, mapped via 0.9^(len-1)
    3. raw (connections count): Legacy proxy, mapped to capped 0-100

    Args:
        raw: Legacy connection count proxy
        max_connections: Max connections for legacy normalization
        path_length: Shortest path length (>= 1)
        connectivity_score: Path-based connectivity score in [0, 1]

    Returns:
        Normalized connectivity score (0-100)
    """
    # Preferred: connectivity_score (0-1)
    if connectivity_score is not None:
        normalized = float(connectivity_score) * 100.0
        return max(0.0, min(100.0, normalized))

    # Alternative: path_length
    if path_length is not None and int(path_length) >= 1:
        score = 0.9 ** (int(path_length) - 1)
        normalized = score * 100.0
        return max(0.0, min(100.0, normalized))

    # Legacy: connection count
    if raw is not None:
        if max_connections is not None and max_connections > 0:
            normalized = (float(raw) / float(max_connections)) * 100.0
        else:
            # Default legacy normalization
            normalized = min((float(raw) / 50.0) * 100.0, 100.0)
        return max(0.0, min(100.0, normalized))

    return 50.0  # Neutral default


def calculate_risk_weighted_safety(breakdown: Dict[str, int]) -> float:
    """
    Calculate risk-weighted safety score from OFF-X breakdown.

    Formula: Safety = sum(count[class] * score[class]) / sum(count[class])

    Categories and base scores:
    - very_high: 0
    - high: 10
    - medium: 20
    - low: 60
    - very_low: 80
    - not_assoc: 100

    NA rows are excluded (treated as unknown).

    Args:
        breakdown: Dict mapping category names to counts

    Returns:
        Safety score (0-100, higher = safer)
    """
    if not breakdown:
        return 50.0

    total_count = 0
    weighted_sum = 0.0

    for category, count in breakdown.items():
        if count is None or count == 0:
            continue

        # Normalize category name
        cat_key = category.lower().replace(' ', '_').replace('-', '_')

        # Skip NA/unknown
        if cat_key in ('na', 'unknown', 'not_available'):
            continue

        score = SAFETY_CATEGORY_SCORES.get(cat_key)
        if score is None:
            continue

        weighted_sum += count * score
        total_count += count

    if total_count == 0:
        return 50.0

    safety = weighted_sum / total_count
    return max(0.0, min(100.0, safety))
