# Score Interpretation Guide

## BioBridge Scores

### Cosine Similarity (cos_sim)

Measures embedding space proximity between head entity and tail entity.

| Score Range | Classification | Interpretation |
|-------------|----------------|----------------|
| > 0.7 | **Very Strong** | High-confidence association, likely validated |
| 0.5 - 0.7 | **Strong** | Substantial evidence, well-supported |
| 0.3 - 0.5 | **Moderate** | Moderate evidence, worth investigating |
| 0.1 - 0.3 | **Weak** | Low evidence, suggestive but not conclusive |
| < 0.1 | **Very Weak** | Minimal evidence, may be noise |

### Percentile Rank (pct_rank)

Position relative to all possible associations for this entity type combination.

| Percentile | Interpretation |
|------------|----------------|
| > 99% | Top 1% of associations — exceptionally strong |
| 95-99% | Top 5% — very strong association |
| 90-95% | Top 10% — strong association |
| 75-90% | Top 25% — above average |
| 50-75% | Average range |
| < 50% | Below average — weak signal |

### Combined Interpretation

| cos_sim | pct_rank | Overall Assessment |
|---------|----------|-------------------|
| > 0.7 | > 95% | **HIGH CONFIDENCE** — validated target |
| 0.5-0.7 | > 90% | **GOOD CONFIDENCE** — strong candidate |
| 0.3-0.5 | > 75% | **MODERATE** — worth following up |
| 0.1-0.3 | > 50% | **LOW** — requires additional evidence |
| < 0.1 | < 50% | **MINIMAL** — likely not meaningful |

---

## ULTRA/UltraQuery Scores

### Model Score

Raw prediction score from the ULTRA foundation model.

| Score Range | Interpretation |
|-------------|----------------|
| > 0.8 | Very high confidence prediction |
| 0.6 - 0.8 | High confidence |
| 0.4 - 0.6 | Moderate confidence |
| 0.2 - 0.4 | Low confidence |
| < 0.2 | Very low confidence |

### Percentile Rank

Similar interpretation to BioBridge — position among all predictions.

### Filtered vs Unfiltered

- **Filtered**: Only schema-valid predictions (recommended)
- **Unfiltered**: All predictions including invalid type combinations

Always use filtered results for interpretation.

---

## PrimeKG Scores (Shortest Path)

### Path Length Score

Exponential decay based on graph distance: `score = 0.9^(path_length - 1)`

| Path Length | Score | Classification | Interpretation |
|-------------|-------|----------------|----------------|
| 1 (direct) | 1.00 | **Direct** | Direct edge in KG - very strong |
| 2 | 0.90 | **Strong** | 1 intermediate node |
| 3 | 0.81 | **Moderate-Strong** | 2 intermediate nodes |
| 4 | 0.729 | **Moderate** | 3 intermediate nodes |
| 5 | 0.656 | **Weak** | 4 intermediate nodes |
| 6 | 0.590 | **Weak** | 5 intermediate nodes |
| 7+ | <0.53 | **Very Weak** | Distant in graph |
| No path | 0.00 | **None** | Disconnected entities |

### Path Quality Assessment

| Score Range | Assessment | Action |
|-------------|------------|--------|
| 1.0 | Direct association | High confidence, check edge type |
| 0.81-0.99 | Strong proximity | Examine intermediate nodes |
| 0.59-0.80 | Moderate proximity | May be indirect relationship |
| 0.01-0.58 | Weak proximity | Distant, investigate path quality |
| 0.0 | No connection | Entities not connected in PrimeKG |

### Combo Score (Average)

PrimeKG combo scores are arithmetic averages of component scores:

```
combo_score = mean(component_scores)
```

**Important**: Unlike BioBridge/ULTRA, PrimeKG combo scores:
- Always equal the arithmetic mean (by definition)
- Cannot detect synergy/dilution
- Delta is always 0

Focus interpretation on:
- Which component contributes highest score
- Range of component scores (consistency)
- Path structure differences between components

---

## Comparing Scores Across Methods

### Direct Comparison

BioBridge, ULTRA, and PrimeKG scores are NOT directly comparable due to:
- Different underlying models
- Different scoring methodologies (embeddings vs queries vs paths)
- Different score calibrations and ranges

### Relative Comparison

Compare **relative rankings** not absolute scores:
- If all 3 methods rank Entity A > Entity B, high confidence
- If methods disagree on ranking, investigate further

### Three-Method Concordance Assessment

| Scenario | Confidence |
|----------|------------|
| All 3 methods agree (same top entities) | **HIGH** |
| 2 of 3 methods agree | **MODERATE** |
| All 3 methods disagree significantly | **LOW** — requires investigation |

### Cross-Method Agreement Matrix

| BioBridge | ULTRA | PrimeKG | Interpretation |
|-----------|-------|---------|----------------|
| Strong | Strong | Short path | **HIGH CONFIDENCE** - validated across all methods |
| Strong | Strong | Long path | Semantic similarity, indirect mechanism |
| Strong | Weak | Short path | BioBridge may capture broader context |
| Weak | Strong | Short path | UltraQuery intersection is specific |
| Weak | Weak | Short path | Direct edge but weak annotation/literature |
| Strong | Strong | No path | Not in PrimeKG but validated by other methods |
| Weak | Weak | Long/No path | **LOW CONFIDENCE** - weak signal across all |

---

## Biological Interpretation Guidelines

### Strong Association Signals

When cos_sim > 0.5 or pct_rank > 90%:

1. **Check for therapeutic validation**
   - Is this a known drug target?
   - Are there clinical trials?
   - FDA-approved therapies?

2. **Check for genetic validation**
   - GWAS hits?
   - Mendelian genetics?
   - Known disease genes?

3. **Check mechanism plausibility**
   - Pathway involvement?
   - Expression in relevant tissue?
   - Known biological function?

### Weak Association Signals

When cos_sim < 0.3 or pct_rank < 75%:

1. **Consider indirect associations**
   - May be downstream effect
   - May require pathway context
   - May be tissue-specific

2. **Check for confounders**
   - General inflammation marker?
   - Housekeeping gene?
   - Cell proliferation signal?

3. **Evaluate in context**
   - Compare to known positive controls
   - Compare to known negatives
   - Check consistency across diseases

---

## Synergy/Dilution Interpretation

### Synergy (Δ > +0.02)

**Biological Meaning**: Combined signature captures more disease-relevant biology than individual components.

**Common Causes**:
- Sequential pathway (A→B in disease)
- Complementary mechanisms
- Shared literature co-citation
- Tissue-specific co-expression

**Example**: TNFRSF25 + PCOLCE
- TNFRSF25: Acute inflammation
- PCOLCE: Chronic fibrosis
- Combined: Captures inflammation→fibrosis progression

### Near-Additive (Δ ≈ 0)

**Biological Meaning**: Signature represents average of components, no emergent property.

**Common Causes**:
- Related but independent pathways
- No direct biological coupling
- Both individually relevant

### Dilution (Δ < -0.02)

**Biological Meaning**: Combined signature loses signal compared to individual components.

**Common Causes**:
- Unrelated pathways (distant embeddings)
- Tissue specificity mismatch
- One strong + one weak component
- No literature co-citation

**Example**: CDKN2D + PCOLCE
- CDKN2D: Cell cycle/cancer
- PCOLCE: ECM/fibrosis
- Combined: Falls into "null space" between cancer and fibrosis

---

## Clinical Relevance Assessment

### Biomarker Utility

| Score Level | Biomarker Recommendation |
|-------------|-------------------------|
| Very Strong (>0.7) | Primary biomarker candidate |
| Strong (0.5-0.7) | Good biomarker, combine with others |
| Moderate (0.3-0.5) | Secondary biomarker, needs validation |
| Weak (<0.3) | Not recommended as standalone |

### Therapeutic Target Assessment

| Evidence | Recommendation |
|----------|----------------|
| Strong association + known mechanism | High priority target |
| Strong association + unknown mechanism | Investigate mechanism first |
| Moderate association + known mechanism | Secondary target |
| Weak association | Low priority unless strong rationale |

### Combination Therapy Prediction

| Combo Effect | Recommendation |
|--------------|----------------|
| Synergy | Strong candidate for combination |
| Near-additive | Possible combination, no synergy benefit |
| Dilution | Avoid combination, use individually |

---

## Reporting Best Practices

### Always Report

1. **Raw scores**: cos_sim and pct_rank
2. **Context**: Entity types and relation
3. **Comparison**: To known positives/negatives
4. **Caveats**: Model limitations, validation needed

### Avoid

1. Overinterpreting small differences (<0.05)
2. Claiming causation from association
3. Ignoring biological plausibility
4. Presenting without validation context

### Recommended Phrasing

**Strong signal**: "BioBridge predicts a strong association (cos_sim=0.72, 98th percentile), consistent with [validation evidence]"

**Weak signal**: "BioBridge shows a weak association (cos_sim=0.18, 65th percentile), suggesting possible indirect relationship requiring further investigation"

**Synergy**: "The [A+B] signature shows synergy (Δ=+0.04), potentially reflecting [biological mechanism]"

**Dilution**: "The [A+B] signature shows dilution (Δ=-0.12), suggesting the genes have unrelated biological contexts"

**PrimeKG Direct Path**: "PrimeKG shows a direct edge (score=1.0, 1 hop) between [Gene] and [Disease], confirming the association exists in the knowledge graph"

**PrimeKG Indirect Path**: "PrimeKG shows an indirect path (score=0.81, 3 hops via [Intermediate1]→[Intermediate2]), suggesting a mechanistic relationship through [pathway/process]"

**PrimeKG No Path**: "No path found in PrimeKG (score=0.0), indicating these entities are not directly connected in the knowledge graph, though BioBridge/ULTRA may still show semantic associations"

---

## Cross-Method Validation

### High-Confidence Findings

A finding has **HIGH CONFIDENCE** when:
1. BioBridge shows pct_rank > 0.95
2. ULTRA shows pct_rank > 0.95 (within entity type)
3. PrimeKG shows score > 0.81 (path length ≤ 3)

### Method-Specific Insights

| Method | Unique Insight |
|--------|----------------|
| **BioBridge** | Captures semantic/literature similarity even without direct KG edge |
| **ULTRA** | Identifies entities satisfying multiple conditions (intersection) |
| **PrimeKG** | Reveals mechanistic path with intermediate nodes |

### When Methods Disagree

1. **BioBridge strong, PrimeKG weak**: Association exists in literature but not curated in KG
2. **PrimeKG strong, BioBridge weak**: Direct KG edge but under-annotated in literature
3. **ULTRA strong, others weak**: Intersection-specific finding, one component dominates

### Recommended Interpretation Workflow

1. Start with BioBridge/ULTRA for initial ranking
2. Validate top candidates with PrimeKG path analysis
3. Examine PrimeKG intermediate nodes for mechanism insights
4. Flag discordant results for manual review
5. Prioritize findings where all 3 methods agree
