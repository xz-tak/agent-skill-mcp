"""
Subscore calculation functions for the 5 dimensions:
1. Clinical Validation
2. Disease Association
3. Safety
4. Opportunity
5. Novelty
"""
import logging
from typing import Dict, List, Optional, Any, Tuple

from .normalizers import (
    normalize_clinical,
    normalize_deg,
    normalize_biobridge,
    normalize_ultra,
    normalize_primekg,
    calculate_risk_weighted_safety
)

logger = logging.getLogger(__name__)

# Weights for Disease Association subscore
DISEASE_ASSOC_WEIGHTS = {
    'deg': 0.40,
    'biobridge': 0.25,
    'ultra': 0.25,
    'primekg': 0.10
}

# Weights for Opportunity subscore
OPPORTUNITY_WEIGHTS = {
    'disease_assoc': 0.40,
    'clinical_novelty': 0.30,
    'ci_score': 0.30
}

# Weights for Novelty subscore
NOVELTY_WEIGHTS = {
    'clinical_novelty': 0.70,
    'literature_novelty': 0.30
}


def calculate_clinical_validation(data: Dict[str, Any], is_combo: bool = False) -> Dict[str, Any]:
    """
    Calculate Clinical Validation subscore.

    Args:
        data: Data dictionary from loader
        is_combo: Whether this is a combination

    Returns:
        Dictionary with raw score, normalized score, and components
    """
    result = {
        'raw_score': None,
        'normalized_score': 50.0,
        'components': {}
    }

    if is_combo:
        # Average individual gene scores
        genes = data.get('genes', [])
        gene_scores = []

        for gene in genes:
            gene_key = f"{gene}_cortellis"
            if gene_key in data:
                raw = data[gene_key].get('total_score_raw')
                if raw is not None:
                    gene_scores.append(raw)
                    result['components'][gene] = {
                        'raw': raw,
                        'normalized': normalize_clinical(raw)
                    }

        if gene_scores:
            result['raw_score'] = sum(gene_scores) / len(gene_scores)
            result['normalized_score'] = normalize_clinical(result['raw_score'])
    else:
        # Individual gene
        if 'cortellis' in data:
            raw = data['cortellis'].get('total_score_raw')
            if raw is not None:
                result['raw_score'] = raw
                result['normalized_score'] = normalize_clinical(raw)

    return result


def calculate_disease_association(
    data: Dict[str, Any],
    is_combo: bool = False,
    max_primekg: Optional[int] = None
) -> Dict[str, Any]:
    """
    Calculate Disease Association subscore.

    Components:
    - DEG: 40%
    - BioBridge: 25%
    - ULTRA: 25%
    - PrimeKG: 10%

    Args:
        data: Data dictionary from loader
        is_combo: Whether this is a combination
        max_primekg: Maximum PrimeKG connectivity proxy for normalization (legacy)

    Returns:
        Dictionary with subscore and component breakdown
    """
    result = {
        'score': 50.0,
        'components': {}
    }

    components_normalized = {}
    weight_sum = sum(DISEASE_ASSOC_WEIGHTS.values()) or 1.0

    if is_combo:
        genes = data.get('genes', [])

        # DEG: Average individual genes
        deg_scores = []
        for gene in genes:
            gene_key = f"{gene}_deg"
            if gene_key in data:
                raw = data[gene_key].get('deg_score_raw')
                if raw is not None:
                    norm = normalize_deg(raw)
                    deg_scores.append(norm)
                    result['components'].setdefault('deg', {})[gene] = {'raw': raw, 'normalized': norm}

        if deg_scores:
            components_normalized['deg'] = sum(deg_scores) / len(deg_scores)
        else:
            components_normalized['deg'] = 50.0

        # BioBridge: Use combo percentile directly
        if 'biobridge' in data:
            raw = data['biobridge'].get('biobridge_percentile_raw')
            if raw is not None:
                norm = normalize_biobridge(raw)
                components_normalized['biobridge'] = norm
                result['components']['biobridge'] = {'raw': raw, 'normalized': norm}
        elif f"{genes[0]}_biobridge" in data:
            # Fallback to first gene if combo not available
            raw = data[f"{genes[0]}_biobridge"].get('biobridge_percentile_raw')
            if raw is not None:
                norm = normalize_biobridge(raw)
                components_normalized['biobridge'] = norm
                result['components']['biobridge'] = {'raw': raw, 'normalized': norm}

        if 'biobridge' not in components_normalized:
            components_normalized['biobridge'] = 50.0

        # ULTRA: Average individual genes
        ultra_scores = []
        for gene in genes:
            gene_key = f"{gene}_ultra"
            if gene_key in data:
                raw = data[gene_key].get('ultra_percentile_raw')
                if raw is not None:
                    norm = normalize_ultra(raw)
                    ultra_scores.append(norm)
                    result['components'].setdefault('ultra', {})[gene] = {'raw': raw, 'normalized': norm}

        if ultra_scores:
            components_normalized['ultra'] = sum(ultra_scores) / len(ultra_scores)
        else:
            components_normalized['ultra'] = 50.0

        # PrimeKG: Average individual genes
        primekg_scores = []
        for gene in genes:
            gene_key = f"{gene}_primekg"
            if gene_key in data:
                raw = data[gene_key].get('primekg_connections_raw')
                if raw is not None:
                    norm = normalize_primekg(raw, max_primekg)
                    primekg_scores.append(norm)
                    result['components'].setdefault('primekg', {})[gene] = {'raw': raw, 'normalized': norm}

        if primekg_scores:
            components_normalized['primekg'] = sum(primekg_scores) / len(primekg_scores)
        else:
            components_normalized['primekg'] = 50.0

    else:
        # Individual gene
        # DEG
        if 'deg' in data:
            raw = data['deg'].get('deg_score_raw')
            if raw is not None:
                norm = normalize_deg(raw)
                components_normalized['deg'] = norm
                result['components']['deg'] = {'raw': raw, 'normalized': norm}

        if 'deg' not in components_normalized:
            components_normalized['deg'] = 50.0

        # BioBridge
        if 'biobridge' in data:
            raw = data['biobridge'].get('biobridge_percentile_raw')
            if raw is not None:
                norm = normalize_biobridge(raw)
                components_normalized['biobridge'] = norm
                result['components']['biobridge'] = {'raw': raw, 'normalized': norm}

        if 'biobridge' not in components_normalized:
            components_normalized['biobridge'] = 50.0

        # ULTRA
        if 'ultra' in data:
            raw = data['ultra'].get('ultra_percentile_raw')
            if raw is not None:
                norm = normalize_ultra(raw)
                components_normalized['ultra'] = norm
                result['components']['ultra'] = {'raw': raw, 'normalized': norm}

        if 'ultra' not in components_normalized:
            components_normalized['ultra'] = 50.0

        # PrimeKG
        if 'primekg' in data:
            raw = data['primekg'].get('primekg_connections_raw')
            if raw is not None:
                norm = normalize_primekg(raw, max_primekg)
                components_normalized['primekg'] = norm
                result['components']['primekg'] = {'raw': raw, 'normalized': norm}

        if 'primekg' not in components_normalized:
            components_normalized['primekg'] = 50.0

    # Calculate weighted score
    final_score = 0.0
    for component, weight in DISEASE_ASSOC_WEIGHTS.items():
        # Normalize weights to sum to 1.0 (keeps disease score on 0-100 scale).
        final_score += components_normalized.get(component, 50.0) * (weight / weight_sum)

    result['score'] = final_score
    result['components_normalized'] = components_normalized

    return result


def calculate_safety(data: Dict[str, Any], is_combo: bool = False) -> Dict[str, Any]:
    """
    Calculate Safety subscore using risk-weighted calculation.

    Args:
        data: Data dictionary from loader
        is_combo: Whether this is a combination

    Returns:
        Dictionary with safety score and breakdown
    """
    result = {
        'score': 50.0,
        'breakdown': None
    }

    if is_combo:
        # Try combo OFF-X data first
        if 'offx' in data:
            breakdown = data['offx'].get('safety_breakdown_raw')
            if breakdown:
                result['score'] = calculate_risk_weighted_safety(breakdown)
                result['breakdown'] = breakdown
                return result

        # Fallback: Average individual gene scores
        genes = data.get('genes', [])
        gene_scores = []

        for gene in genes:
            gene_key = f"{gene}_offx"
            if gene_key in data:
                breakdown = data[gene_key].get('safety_breakdown_raw')
                if breakdown:
                    score = calculate_risk_weighted_safety(breakdown)
                    gene_scores.append(score)

        if gene_scores:
            result['score'] = sum(gene_scores) / len(gene_scores)

    else:
        # Individual gene
        if 'offx' in data:
            breakdown = data['offx'].get('safety_breakdown_raw')
            if breakdown:
                result['score'] = calculate_risk_weighted_safety(breakdown)
                result['breakdown'] = breakdown

    return result


def calculate_opportunity(
    data: Dict[str, Any],
    disease_assoc_score: float,
    clinical_validation_score: float,
    is_combo: bool = False,
    pathway_synergy: float = 0.0
) -> Dict[str, Any]:
    """
    Calculate Opportunity subscore.

    Components:
    - Disease Association: 50%
    - Clinical Novelty: 30%
    - Competitive Intelligence: 20%

    Args:
        data: Data dictionary from loader
        disease_assoc_score: Pre-calculated disease association score
        clinical_validation_score: Pre-calculated clinical validation score
        is_combo: Whether this is a combination
        pathway_synergy: Pathway synergy score (for combos)

    Returns:
        Dictionary with opportunity score and components
    """
    result = {
        'score': 50.0,
        'components': {}
    }

    # Clinical Novelty (inverse of clinical validation)
    clinical_novelty = 100.0 - clinical_validation_score
    result['components']['clinical_novelty'] = clinical_novelty

    # Competitive Intelligence Score
    ci_details = calculate_ci_details(data, is_combo=is_combo)
    ci_score = float(ci_details.get('ci_blend', 50.0))
    result['components']['ci_score'] = ci_score
    result['components']['ci_gene'] = float(ci_details.get('ci_gene', 50.0))
    result['components']['ci_family'] = float(ci_details.get('ci_family', 50.0))
    if ci_details.get('family'):
        result['components']['ci_family_label'] = ci_details.get('family')

    if is_combo:
        # Calculate individual opportunities
        genes = data.get('genes', [])
        # For combos: mean individual opportunities (60%) + novel mechanism bonus (40%)

        # Mean individual opportunities (simplified - would need individual scores)
        mean_individual_opps = (disease_assoc_score * 0.4 + clinical_novelty * 0.3 + ci_score * 0.3)

        # Novel mechanism bonus (inverse of synergy)
        combined_synergy = pathway_synergy  # Simplified
        novel_mechanism_bonus = 100.0 - (combined_synergy * 100.0)

        result['score'] = (mean_individual_opps * 0.6) + (novel_mechanism_bonus * 0.4)
        result['components']['mean_individual_opps'] = mean_individual_opps
        result['components']['novel_mechanism_bonus'] = novel_mechanism_bonus

    else:
        # Individual gene: weighted combination
        result['score'] = (
            disease_assoc_score * OPPORTUNITY_WEIGHTS['disease_assoc'] +
            clinical_novelty * OPPORTUNITY_WEIGHTS['clinical_novelty'] +
            ci_score * OPPORTUNITY_WEIGHTS['ci_score']
        )

    return result


def _ci_phase_weight(phase_label: str) -> float:
    phase = (phase_label or "").strip()
    if "Marketed" in phase:
        return 1.0
    if "Phase III" in phase:
        return 0.7
    if "Phase II" in phase:
        return 0.4
    if "Phase I" in phase:
        return 0.2
    if "Preclinical" in phase:
        return 0.1
    return 0.0


def _calculate_ci_gene_specific(programs_by_phase: Dict[str, Any], total_programs: float) -> Optional[float]:
    if not programs_by_phase or not total_programs or total_programs <= 0:
        return None
    weighted = 0.0
    weighted += float(programs_by_phase.get('Marketed', 0) or 0) * 1.0
    weighted += float(programs_by_phase.get('Phase III', 0) or 0) * 0.7
    weighted += float(programs_by_phase.get('Phase II', 0) or 0) * 0.4
    weighted += float(programs_by_phase.get('Phase I', 0) or 0) * 0.2
    weighted += float(programs_by_phase.get('Preclinical', 0) or 0) * 0.1
    ci = 100.0 - ((weighted / float(total_programs)) * 100.0)
    return max(0.0, min(ci, 100.0))


def _calculate_family_weighted_totals(ci_dashboard: Dict[str, Any]) -> Tuple[Dict[str, float], float]:
    """
    Return (family_weighted_totals, all_weighted_total) using CI dashboard entries.
    Only entries with non-empty `ibdTags` are counted (IBD-relevant programs).
    """
    entries = ci_dashboard.get('entries', []) if isinstance(ci_dashboard, dict) else []
    if not isinstance(entries, list):
        return {}, 0.0

    fam_totals: Dict[str, float] = {}
    all_total = 0.0
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        ibd_tags = entry.get('ibdTags', [])
        if not ibd_tags:
            continue
        fam = entry.get('targetFamilyPrimary') or 'Other'
        if not isinstance(fam, str) or not fam.strip():
            fam = 'Other'
        weight = _ci_phase_weight(str(entry.get('ibdPhase', '') or ''))
        if weight <= 0:
            continue
        fam_totals[fam] = fam_totals.get(fam, 0.0) + weight
        all_total += weight

    return fam_totals, all_total


def _calculate_ci_family_whitespace(
    family: Optional[str],
    fam_totals: Dict[str, float],
    all_total: float
) -> Optional[float]:
    if not family or not fam_totals or all_total <= 0:
        return None
    fam_weighted = fam_totals.get(family, 0.0)
    frac = fam_weighted / all_total
    ci = 100.0 - (frac * 100.0)
    return max(0.0, min(ci, 100.0))


def calculate_ci_details(data: Dict[str, Any], is_combo: bool = False) -> Dict[str, Any]:
    """
    Calculate blended CI details:
      - CI_gene_specific: per-target (phase-weighted) whitespace proxy
      - CI_family_whitespace: family crowding whitespace proxy
      - CI_blended: 0.5*gene + 0.5*family

    For combinations: compute per-gene blended CI and return the mean.
    """
    result: Dict[str, Any] = {
        'ci_gene': 50.0,
        'ci_family': 50.0,
        'ci_blend': 50.0,
        'family': None,
    }

    ci_dashboard = data.get('ci_dashboard') if isinstance(data, dict) else None
    fam_totals, all_total = _calculate_family_weighted_totals(ci_dashboard) if isinstance(ci_dashboard, dict) else ({}, 0.0)

    if is_combo:
        genes = data.get('genes', [])
        per_gene = []
        for gene in genes:
            gene_key = f"{gene}_ci"
            if gene_key not in data:
                continue
            gci = data[gene_key]
            if not isinstance(gci, dict):
                continue
            gene_specific = _calculate_ci_gene_specific(
                gci.get('programs_by_phase', {}) or {},
                float(gci.get('total_programs', 0) or 0),
            )
            family = gci.get('target_family_primary')
            family_ci = _calculate_ci_family_whitespace(family, fam_totals, all_total)
            if gene_specific is None:
                gene_specific = 50.0
            if family_ci is None:
                family_ci = 50.0
            per_gene.append((gene_specific, family_ci, family))

        if not per_gene:
            return result

        result['ci_gene'] = sum(x[0] for x in per_gene) / len(per_gene)
        result['ci_family'] = sum(x[1] for x in per_gene) / len(per_gene)
        # Not meaningful for combos, but keep the most common family label for display/debug.
        fam_counts: Dict[str, int] = {}
        for _, _, fam in per_gene:
            if isinstance(fam, str) and fam.strip():
                fam_counts[fam] = fam_counts.get(fam, 0) + 1
        if fam_counts:
            result['family'] = sorted(fam_counts.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)[0][0]

    else:
        gci = data.get('ci', {}) if isinstance(data.get('ci', {}), dict) else {}
        gene_specific = _calculate_ci_gene_specific(
            gci.get('programs_by_phase', {}) or {},
            float(gci.get('total_programs', 0) or 0),
        )
        family = gci.get('target_family_primary')
        family_ci = _calculate_ci_family_whitespace(family, fam_totals, all_total)

        result['family'] = family
        result['ci_gene'] = gene_specific if gene_specific is not None else 50.0
        result['ci_family'] = family_ci if family_ci is not None else 50.0

    result['ci_blend'] = 0.5 * float(result['ci_gene']) + 0.5 * float(result['ci_family'])
    return result


def calculate_ci_score(data: Dict[str, Any], is_combo: bool = False, pathway_synergy: float = 0.0) -> float:
    """
    Calculate Competitive Intelligence score.

    For individuals: % market left = 100 - (marketed / total * 100)
    For combos: average of member-gene % market left (no scaling / synergy adjustment)

    Args:
        data: Data dictionary
        is_combo: Whether this is a combination
        pathway_synergy: Pathway synergy (for combo adjustment)

    Returns:
        CI score (0-100, higher = more opportunity)
    """
    details = calculate_ci_details(data, is_combo=is_combo)
    return float(details.get('ci_blend', 50.0))


def calculate_novelty(
    clinical_validation_score: float,
    primekg_connectivity_score: Optional[float] = None,
    primekg_path_length: Optional[int] = None,
    primekg_connections: Optional[int] = None,
    is_combo: bool = False,
) -> Dict[str, Any]:
    """Calculate Novelty subscore.

    Components:
    - Clinical Novelty: 70%
    - Literature Novelty: 30%

    Literature novelty is defined as the inverse of a PrimeKG *connectivity* signal:
    higher connectivity => lower novelty.

    Supported raw inputs (prefer in this order):
    - primekg_connectivity_score: expected in [0, 1]
    - primekg_path_length: integer path length (>= 1), mapped to score via 0.9^(len-1)
    - primekg_connections: legacy proxy (connection count), mapped to a capped 0-100 percentile

    Args:
        clinical_validation_score: Pre-calculated clinical validation score
        primekg_connectivity_score: PrimeKG connectivity score in [0, 1] (preferred)
        primekg_path_length: PrimeKG shortest path length (optional)
        primekg_connections: Legacy PrimeKG connection count proxy (optional)
        is_combo: Whether this is a combination

    Returns:
        Dictionary with novelty score and components
    """

    result: Dict[str, Any] = {
        'score': 50.0,
        'components': {},
    }

    # Clinical Novelty (inverse of clinical validation)
    clinical_novelty = 100.0 - clinical_validation_score
    result['components']['clinical_novelty'] = clinical_novelty

    # PrimeKG connectivity (0-100), then invert to novelty.
    if primekg_connectivity_score is not None:
        connectivity_100 = max(0.0, min(float(primekg_connectivity_score) * 100.0, 100.0))
    elif primekg_path_length is not None and int(primekg_path_length) >= 1:
        score = 0.9 ** (int(primekg_path_length) - 1)
        connectivity_100 = max(0.0, min(float(score) * 100.0, 100.0))
    elif primekg_connections is not None:
        # Legacy proxy: more connections => more connectivity.
        connectivity_100 = min((float(primekg_connections) / 50.0) * 100.0, 100.0)
    else:
        connectivity_100 = 50.0

    result['components']['primekg_connectivity'] = float(connectivity_100)
    literature_novelty = 100.0 - float(connectivity_100)
    result['components']['literature_novelty'] = literature_novelty

    # Calculate weighted score
    result['score'] = (
        clinical_novelty * NOVELTY_WEIGHTS['clinical_novelty'] +
        literature_novelty * NOVELTY_WEIGHTS['literature_novelty']
    )

    if is_combo:
        # Apply 10% combination premium
        result['score'] *= 1.10
        result['combo_premium'] = True

    return result
