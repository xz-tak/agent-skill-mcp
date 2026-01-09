# BioBridge Association Analysis

## Overview

BioBridge is a neural knowledge graph link prediction system that uses multimodal embeddings to predict associations between biomedical entities. It leverages pre-trained embeddings from the BioBridge model to compute semantic similarity between entities in the knowledge graph.

**Key Characteristics:**
- Uses embedding-based similarity (cosine similarity)
- Supports mean embeddings for gene signatures/combos
- Returns both raw similarity scores and percentile ranks
- Can detect synergy/dilution effects in combinations

## When to Use BioBridge

| Use Case | Recommendation |
|----------|----------------|
| Single gene-disease association | ✓ Good choice |
| Gene signature/combo analysis | ✓ Excellent - uses mean embeddings |
| Semantic similarity queries | ✓ Best choice |
| Synergy/dilution detection | ✓ Supported via pct_rank delta |
| Complex multi-hop reasoning | Use ULTRA instead |
| Network topology analysis | Use PrimeKG instead |

## Scoring Methodology

### Primary Score: Percentile Rank (pct_rank)

BioBridge returns predictions ranked by cosine similarity. The `pct_rank` score indicates where a prediction ranks among all possible tail entities:

| pct_rank | Interpretation | Assessment |
|----------|---------------|------------|
| ≥ 0.99 | Top 1% of all predictions | **Very Strong** |
| 0.95-0.99 | Top 5% | **Strong** |
| 0.90-0.95 | Top 10% | Moderate-Strong |
| 0.80-0.90 | Top 20% | Moderate |
| 0.50-0.80 | Top 50% | Weak |
| < 0.50 | Below median | Very Weak |

### Secondary Score: Cosine Similarity (cos_sim)

The raw cosine similarity between entity embeddings is also reported for reference:
- Range: -1 to 1 (typically 0.3-0.9 for related entities)
- Higher values indicate stronger semantic similarity
- Useful for comparing relative strength within the same query

---

## API Requirements

### CRITICAL: Explicit Entity Overrides

When calling `mcp__biobridge__predict_associations`, you **MUST ALWAYS** provide explicit entity overrides. Never rely on context inference alone.

**Why explicit overrides are required:**
1. **Determinism**: Ensures reproducible results across runs
2. **Accuracy**: Prevents LLM inference errors on entity names (e.g., "TYK2" vs "Tyrosine Kinase 2")
3. **Matching**: Guarantees exact entity matching in the knowledge graph
4. **Scoring**: Required for accurate pct_rank extraction

---

## Individual Entity Analysis

### Required Parameters

| Parameter | Required | Description | Example |
|-----------|----------|-------------|---------|
| `context` | Yes | Natural language description | "Find association between TYK2 and Crohn disease" |
| `override_head_name` | **MANDATORY** | Exact head entity name | "TYK2" |
| `override_head_type` | **MANDATORY** | Head entity type | "gene/protein" |
| `override_tail_name` | **MANDATORY** | Exact tail entity name | "Crohn disease" |
| `override_tail_type` | **MANDATORY** | Tail entity type | "disease" |
| `relation_hint` | Recommended | Relation family | "associated with" |
| `topk` | Recommended | Number of results | 1 (for pair validation) |

### Individual Query Code Example

```python
# Query each gene-disease pair individually
for gene in genes:
    for disease in diseases:
        result = mcp__biobridge__predict_associations(
            context=f"Find association between {gene} and {disease}",
            override_head_name=gene,           # REQUIRED
            override_head_type="gene/protein", # REQUIRED
            override_tail_name=disease,        # REQUIRED
            override_tail_type="disease",      # REQUIRED
            relation_hint="associated with",
            topk=1
        )

        # Extract scores from result
        pct_rank = result['pct_rank']   # PRIMARY score (0-1)
        cos_sim = result['cos_sim']     # Reference score

        # Store results
        results[gene][disease] = {
            'pct_rank': pct_rank,
            'cos_sim': cos_sim,
            'assessment': get_assessment(pct_rank)
        }

def get_assessment(pct_rank):
    if pct_rank >= 0.99: return "Very Strong"
    elif pct_rank >= 0.95: return "Strong"
    elif pct_rank >= 0.90: return "Moderate-Strong"
    elif pct_rank >= 0.80: return "Moderate"
    elif pct_rank >= 0.50: return "Weak"
    else: return "Very Weak"
```

### Individual Analysis Report Format

**Part 1: Individual Analysis Table**

```markdown
## BioBridge Individual Gene-Disease Analysis

### Crohn Disease

| Gene | Score (pct_rank) | Cos_Sim | Assessment |
|------|------------------|---------|------------|
| ITGA4 | 0.9959 | 0.776 | Very Strong |
| TYK2 | 0.9940 | 0.752 | Very Strong |
| JAK1 | 0.9953 | 0.761 | Very Strong |
| IL17A | 0.9821 | 0.698 | Strong |

**Key Findings:**
- ITGA4 shows the strongest association (top 0.4%)
- All JAK-STAT pathway genes show very strong associations
- IL17A shows strong but slightly lower association
```

---

## Combo/Signature Analysis

### Mean Embedding Approach

BioBridge computes combo associations by taking the **mean of component embeddings**:

```
combo_embedding = mean(embed(gene1), embed(gene2), ...)
similarity = cosine_similarity(combo_embedding, disease_embedding)
```

This approach can detect **synergy** (combo stronger than expected) or **dilution** (combo weaker than expected).

### Required Parameters for Combos

| Parameter | Required | Description | Example |
|-----------|----------|-------------|---------|
| `context` | Yes | Natural language description | "Find diseases for TYK2+JAK1 signature" |
| `override_head_names` | **MANDATORY** | List of head entity names | ["TYK2", "JAK1"] |
| `override_head_type` | **MANDATORY** | Head entity type | "gene/protein" |
| `override_tail_name` | **MANDATORY** | Exact tail entity name | "Crohn disease" |
| `override_tail_type` | **MANDATORY** | Tail entity type | "disease" |
| `relation_hint` | Recommended | Relation family | "associated with" |
| `topk` | Recommended | Number of results | 1 |

**Note:** Use `override_head_names` (plural) for combos, not `override_head_name`.

### Combo Query Code Example

```python
# Query each combo-disease pair
for combo in combos:  # e.g., [["TYK2", "JAK1"], ["ITGA4", "ITGB7"]]
    combo_name = "+".join(combo)

    for disease in diseases:
        result = mcp__biobridge__predict_associations(
            context=f"Find association between {combo_name} signature and {disease}",
            override_head_names=combo,         # REQUIRED: list of genes
            override_head_type="gene/protein", # REQUIRED
            override_tail_name=disease,        # REQUIRED
            override_tail_type="disease",      # REQUIRED
            relation_hint="associated with",
            topk=1
        )

        # Extract combo scores
        combo_pct_rank = result['pct_rank']
        combo_cos_sim = result['cos_sim']

        # Calculate synergy/dilution
        individual_scores = [
            individual_results[gene][disease]['pct_rank']
            for gene in combo
        ]
        individual_mean = sum(individual_scores) / len(individual_scores)
        delta = combo_pct_rank - individual_mean

        # Classify effect
        if delta > 0.02:
            classification = "SYNERGY"
        elif delta < -0.02:
            classification = "DILUTION"
        else:
            classification = "NEAR-ADDITIVE"

        results[combo_name][disease] = {
            'combo_score': combo_pct_rank,
            'cos_sim': combo_cos_sim,
            'individual_mean': individual_mean,
            'delta': delta,
            'classification': classification
        }
```

### Synergy/Dilution Interpretation

| Delta (Δ) | Classification | Interpretation |
|-----------|----------------|----------------|
| > +0.02 | **SYNERGY** | Combo is stronger than the average of its components. The genes may act on complementary pathways that reinforce the disease association. |
| -0.02 to +0.02 | **NEAR-ADDITIVE** | Combo performs as expected based on component averages. No emergent synergy or interference detected. |
| < -0.02 | **DILUTION** | Combo is weaker than expected. One gene may be "diluting" the signal of a stronger component, or the genes may have conflicting associations. |

### Combo Analysis Report Format

**Part 2: Combo Analysis Table**

```markdown
## BioBridge Combo Analysis

### Crohn Disease

| Combo | Combo Score | Cos_Sim | Individual Mean | Delta | Classification |
|-------|-------------|---------|-----------------|-------|----------------|
| TYK2+JAK1 | 0.9951 | 0.565 | 0.9947 | +0.0004 | NEAR-ADDITIVE |
| ITGA4+ITGB7 | 0.9972 | 0.612 | 0.9945 | +0.0027 | SYNERGY |
| IL17A+IL17F | 0.9756 | 0.521 | 0.9821 | -0.0065 | DILUTION |

**Synergy Analysis:**
- **ITGA4+ITGB7**: Shows synergy (+0.27%). These integrins form a functional heterodimer,
  and their combined signal is stronger than individual components.
- **IL17A+IL17F**: Shows dilution (-0.65%). While both are IL-17 family members, IL17F
  has weaker disease associations that dilute the IL17A signal.
```

---

## Comparative Analysis (Part 3)

Part 3 provides deep comparative analysis for each combo:

### Report Structure

```markdown
## BioBridge Part 3: Comparative Analysis

### TYK2+JAK1 Signature Analysis

#### Performance Comparison
| Entity | Score | vs Combo |
|--------|-------|----------|
| TYK2+JAK1 (combo) | 0.9951 | baseline |
| TYK2 (individual) | 0.9940 | -0.0011 |
| JAK1 (individual) | 0.9953 | +0.0002 |

#### Biological Interpretation
- TYK2 and JAK1 are both JAK family kinases involved in cytokine signaling
- Their combined embedding captures shared JAK-STAT pathway biology
- Near-additive effect suggests overlapping rather than complementary mechanisms

#### Clinical Implications
- Dual JAK inhibition (TYK2+JAK1) may not provide synergistic efficacy
- Consider alternative combinations with complementary mechanisms
- Monitor for overlapping safety signals from shared pathway inhibition
```

---

## Output Files

BioBridge analysis generates the following files in `<output_dir>/biobridge/`:

| File | Generated When | Contents |
|------|----------------|----------|
| `Part1_Individual_Analysis.md` | Always | Individual gene-disease scores and assessments |
| `Part2_Combo_Analysis.md` | If combos provided | Combo scores with synergy/dilution classification |
| `Part3_Comparative_Analysis.md` | If combos provided | Deep comparative analysis with biological interpretation |

---

## Troubleshooting

### Common Issues

**Issue: Low or zero pct_rank scores**
- Check entity names match exactly (case-sensitive)
- Verify entity type is correct ("gene/protein" not "gene")
- Ensure entity exists in BioBridge knowledge graph

**Issue: Unexpected synergy/dilution**
- Verify individual scores were calculated correctly
- Check that all combo components have valid scores
- Consider biological plausibility of the result

**Issue: MCP call fails**
- Ensure all MANDATORY override parameters are provided
- Check that `override_head_names` (plural) is used for combos
- Verify entity names don't contain special characters

### Validation Checklist

- [ ] All override parameters provided for every call
- [ ] Using `override_head_names` (plural) for combos
- [ ] Entity types match expected values
- [ ] pct_rank used as primary score (not cos_sim)
- [ ] Synergy/dilution calculated from pct_rank delta
