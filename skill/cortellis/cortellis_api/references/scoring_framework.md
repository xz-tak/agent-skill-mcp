# Flexible Scoring Framework for Target Analysis

This document describes the scoring methodology for evaluating clinical validation and competitive intelligence of drug targets. The framework is designed to be flexible and context-dependent.

## Overview

The scoring system quantifies the level of clinical validation for a drug target based on:
1. **Development Phase** - How far drugs have progressed
2. **Indication Relevance** - Whether drugs target the disease of interest
3. **Context-Specific Weights** - Adjustable based on therapeutic area and analysis goals

## Default Scoring Matrix

### Standard Clinical Validation Scores

| Phase | Disease-Specific | Non-Specific | Rationale |
|-------|------------------|--------------|-----------|
| **FDA Approved / On Market** | 7 | 4 | Highest validation; disease-specific indicates proven efficacy |
| **Phase 3 Clinical** | 3 | 2 | Late-stage validation; high probability of approval |
| **Phase 2 Clinical** | 2 | 1 | Mid-stage validation; proof-of-concept established |
| **Phase 1 Clinical** | 1 | 0.5 | Early validation; safety and PK/PD characterized |
| **Preclinical** | 0.1 | 0.1 | Minimal validation; target engagement shown |
| **Discontinued / Suspended** | 0 | 0 | No active development; may indicate safety/efficacy issues |
| **No Development Reported** | 0 | 0 | No validation data |

**Disease-Specific Indications** are scored higher because they represent direct evidence of efficacy in the target disease.

## Customizing Scoring Weights

### Context-Specific Adjustments

The scoring weights should be adjusted based on:

#### 1. Therapeutic Area Risk Tolerance

**High-Risk Tolerance (Oncology, Rare Diseases):**
```python
SCORING_WEIGHTS = {
    'approved': {'disease_specific': 10, 'non_specific': 5},
    'phase_3': {'disease_specific': 5, 'non_specific': 3},
    'phase_2': {'disease_specific': 3, 'non_specific': 2},
    'phase_1': {'disease_specific': 2, 'non_specific': 1},
    'preclinical': {'disease_specific': 0.5, 'non_specific': 0.5}
}
```
*Rationale:* Higher weights for early stages reflect value of any clinical data in high-risk areas.

**Conservative Tolerance (Common Diseases, Crowded Markets):**
```python
SCORING_WEIGHTS = {
    'approved': {'disease_specific': 10, 'non_specific': 2},
    'phase_3': {'disease_specific': 4, 'non_specific': 1},
    'phase_2': {'disease_specific': 2, 'non_specific': 0.5},
    'phase_1': {'disease_specific': 0.5, 'non_specific': 0.1},
    'preclinical': {'disease_specific': 0.05, 'non_specific': 0.05}
}
```
*Rationale:* Emphasizes approved drugs; de-emphasizes early-stage programs in competitive markets.

#### 2. Analysis Purpose

**Competitive Intelligence (identify crowded targets):**
```python
# Count all programs equally to show competitive intensity
SCORING_WEIGHTS = {
    'approved': {'disease_specific': 1, 'non_specific': 1},
    'phase_3': {'disease_specific': 1, 'non_specific': 1},
    'phase_2': {'disease_specific': 1, 'non_specific': 1},
    'phase_1': {'disease_specific': 1, 'non_specific': 1},
    'preclinical': {'disease_specific': 1, 'non_specific': 1}
}
```

**Clinical Validation (biological confidence):**
```python
# Weight by clinical maturity only
SCORING_WEIGHTS = {
    'approved': {'disease_specific': 10, 'non_specific': 10},
    'phase_3': {'disease_specific': 5, 'non_specific': 5},
    'phase_2': {'disease_specific': 2, 'non_specific': 2},
    'phase_1': {'disease_specific': 0.5, 'non_specific': 0.5},
    'preclinical': {'disease_specific': 0.1, 'non_specific': 0.1}
}
```

**Repurposing Opportunities:**
```python
# Higher weight for approved drugs in adjacent indications
SCORING_WEIGHTS = {
    'approved': {'disease_specific': 10, 'non_specific': 8},  # High non-specific weight
    'phase_3': {'disease_specific': 5, 'non_specific': 4},
    'phase_2': {'disease_specific': 2, 'non_specific': 1.5},
    'phase_1': {'disease_specific': 1, 'non_specific': 0.5},
    'preclinical': {'disease_specific': 0.1, 'non_specific': 0.05}
}
```

#### 3. Disease Similarity

For related diseases (e.g., Crohn's and Ulcerative Colitis in IBD):

```python
# Define indication groups
IBD_INDICATIONS = [
    'inflammatory bowel disease', 'ibd',
    "crohn's disease", 'crohn',
    'ulcerative colitis', 'uc', 'colitis'
]

AUTOIMMUNE_INDICATIONS = [
    'rheumatoid arthritis', 'psoriasis', 'lupus',
    'multiple sclerosis', 'autoimmune disease'
]

# Score with partial credit for related diseases
def calculate_indication_score(indication):
    if matches_primary_disease(indication):
        return 1.0  # Full credit
    elif matches_related_disease(indication):
        return 0.5  # Partial credit
    else:
        return 0.0  # No credit
```

## Implementation Example

### Configurable Scoring Function

```python
def calculate_drug_score(phase, indication, scoring_config):
    """
    Calculate score for a drug with flexible weights.

    Args:
        phase: Development phase string
        indication: Indication text
        scoring_config: Dictionary with scoring weights

    Returns:
        float: Calculated score
    """
    # Parse phase to standard format
    phase_normalized = normalize_phase(phase)

    # Check indication relevance
    is_disease_specific = check_indication_match(
        indication,
        scoring_config['target_indications']
    )

    # Get weight from config
    phase_key = phase_normalized.lower().replace(' ', '_')
    indication_key = 'disease_specific' if is_disease_specific else 'non_specific'

    score = scoring_config['weights'].get(phase_key, {}).get(indication_key, 0)

    # Apply modifiers
    if scoring_config.get('apply_company_weight'):
        company_multiplier = get_company_weight(drug_company)
        score *= company_multiplier

    if scoring_config.get('apply_recency_weight'):
        recency_multiplier = get_recency_weight(last_update_date)
        score *= recency_multiplier

    return score
```

### Configuration Examples

**Example 1: IBD Analysis (Default)**
```python
IBD_CONFIG = {
    'target_indications': [
        'inflammatory bowel disease', 'ibd',
        "crohn's disease", 'ulcerative colitis'
    ],
    'weights': {
        'approved': {'disease_specific': 7, 'non_specific': 4},
        'phase_3_clinical': {'disease_specific': 3, 'non_specific': 2},
        'phase_2_clinical': {'disease_specific': 2, 'non_specific': 1},
        'phase_1_clinical': {'disease_specific': 1, 'non_specific': 0.5},
        'preclinical': {'disease_specific': 0.1, 'non_specific': 0.1},
        'discontinued': {'disease_specific': 0, 'non_specific': 0}
    },
    'apply_company_weight': False,
    'apply_recency_weight': False
}
```

**Example 2: Oncology Analysis (High Risk Tolerance)**
```python
ONCOLOGY_CONFIG = {
    'target_indications': [
        'cancer', 'tumor', 'carcinoma', 'leukemia',
        'lymphoma', 'melanoma', 'sarcoma'
    ],
    'weights': {
        'approved': {'disease_specific': 10, 'non_specific': 5},
        'phase_3_clinical': {'disease_specific': 5, 'non_specific': 3},
        'phase_2_clinical': {'disease_specific': 3, 'non_specific': 2},
        'phase_1_clinical': {'disease_specific': 2, 'non_specific': 1},
        'preclinical': {'disease_specific': 0.5, 'non_specific': 0.5},
        'discontinued': {'disease_specific': -1, 'non_specific': -1}  # Penalty for failures
    },
    'apply_company_weight': True,  # Weight by company quality
    'apply_recency_weight': True   # Prefer recent programs
}
```

**Example 3: Competitive Intensity (Equal Weights)**
```python
COMPETITION_CONFIG = {
    'target_indications': ['all'],  # Count all
    'weights': {
        'approved': {'disease_specific': 1, 'non_specific': 1},
        'phase_3_clinical': {'disease_specific': 1, 'non_specific': 1},
        'phase_2_clinical': {'disease_specific': 1, 'non_specific': 1},
        'phase_1_clinical': {'disease_specific': 1, 'non_specific': 1},
        'preclinical': {'disease_specific': 1, 'non_specific': 1},
        'discontinued': {'disease_specific': 0, 'non_specific': 0}
    },
    'apply_company_weight': False,
    'apply_recency_weight': False
}
```

## Advanced Scoring Modifiers

### 1. Company Quality Weighting

```python
def get_company_weight(company_name):
    """Weight by company development success rate."""
    COMPANY_TIERS = {
        'tier_1': 1.5,  # Top pharma with high success rates
        'tier_2': 1.2,  # Mid-tier pharma
        'tier_3': 1.0,  # Biotech / smaller companies
    }

    tier = classify_company(company_name)
    return COMPANY_TIERS.get(tier, 1.0)
```

### 2. Recency Weighting

```python
from datetime import datetime, timedelta

def get_recency_weight(last_update_date, decay_years=3):
    """Apply exponential decay for old programs."""
    if not last_update_date:
        return 1.0

    years_old = (datetime.now() - last_update_date).days / 365
    if years_old <= 1:
        return 1.0
    elif years_old > decay_years:
        return 0.5
    else:
        # Linear decay
        return 1.0 - (0.5 * (years_old - 1) / (decay_years - 1))
```

### 3. Mechanism Diversity Bonus

```python
def calculate_mechanism_bonus(drugs_list):
    """Bonus for diversity of mechanisms."""
    unique_mechanisms = set(drug['mechanism'] for drug in drugs_list)

    # Bonus points for each unique mechanism beyond the first
    diversity_bonus = (len(unique_mechanisms) - 1) * 0.5
    return diversity_bonus
```

## Interpreting Scores

### Score Ranges and Interpretation

| Score Range | Interpretation | Strategic Guidance |
|-------------|----------------|-------------------|
| **> 100** | Extremely crowded, extensively validated | High competition; need strong differentiation |
| **50-100** | Well-validated, high competition | Proven biology; focus on novelty |
| **20-50** | Moderate validation, active development | Balanced risk-reward; good opportunity |
| **5-20** | Emerging target, limited validation | Early stage; first-mover advantage possible |
| **1-5** | Minimal activity, exploratory | High risk; potential breakthrough |
| **0** | No development activity | Blue ocean or undruggable target |

## Best Practices

1. **Define Context First:** Choose weights based on analysis purpose before scoring
2. **Document Assumptions:** Clearly state which scoring config was used
3. **Sensitivity Analysis:** Test multiple scoring schemes to ensure robust conclusions
4. **Normalize by Target Class:** Compare scores within similar target types (kinases vs GPCRs)
5. **Combine with Qualitative Assessment:** Scores guide but don't replace expert judgment

## Reference Implementation

See complete working example:
- **Analysis Script:** `/home/sagemaker-user/claude_code/mockup/cortellis/analyze_ibd_targets.py`
- **Report Generator:** `/home/sagemaker-user/claude_code/mockup/cortellis/generate_report.py`

These scripts demonstrate:
- Flexible scoring implementation
- Phase parsing
- Indication matching
- Score calculation and aggregation
- Result interpretation and ranking
