# Combo Methods: Mean Embeddings vs Intersection vs Shortest Path

## Overview

This skill supports three methods for analyzing entity combinations (signatures):

1. **Mean Embeddings** (BioBridge): Average embedding space representation
2. **Intersection Queries** (UltraQuery): Logical AND operation in query space
3. **Shortest Path** (PrimeKG): Graph distance with exponential decay scoring

## Method 1: Mean Embeddings (BioBridge)

### Concept

Combines multiple entity embeddings by computing their arithmetic mean:

```
signature_embedding = (embedding_A + embedding_B + ... + embedding_N) / N
```

Then predicts associations using the mean embedding as the query.

### Implementation

```python
mcp__biobridge__predict_associations(
    override_head_names=["GENE_A", "GENE_B", "GENE_C"],
    override_head_type="gene/protein",
    override_tail_name="Crohn disease",
    override_tail_type="disease",
    topk=1
)
```

### Characteristics

| Aspect | Description |
|--------|-------------|
| **Operation** | Geometric averaging in embedding space |
| **Effect** | Creates "centroid" representing combined entity |
| **Synergy Detection** | Can detect if centroid lands in high-association region |
| **Dilution Risk** | Weak entities can pull centroid away from strong region |

### When Synergy Occurs

Synergy (combo > mean of individuals) occurs when:
1. Entities share pathway biology (embeddings are nearby)
2. Combined centroid lands in disease-relevant region
3. Literature co-cites both entities with disease

### When Dilution Occurs

Dilution (combo < mean of individuals) occurs when:
1. Entities have unrelated biology (distant embeddings)
2. Centroid lands in low-association "null space"
3. Weak entity pulls centroid away from strong associations

### Example

```
TNFRSF25 embedding: Near "inflammation" region
PCOLCE embedding: Near "fibrosis" region
Mean embedding: Between inflammation and fibrosis → "fibrostenotic IBD" region
Result: SYNERGY (captures inflammation→fibrosis pathway)

CDKN2D embedding: Near "cancer/cell cycle" region
PCOLCE embedding: Near "ECM/fibrosis" region
Mean embedding: Between cancer and ECM → low-association null space
Result: SEVERE DILUTION
```

---

## Method 2: Intersection Queries (UltraQuery)

### Concept

Finds entities that satisfy ALL specified conditions (logical AND):

```
Query: Which diseases are associated with BOTH Gene_A AND Gene_B?
Result: Diseases that have strong associations with ALL genes
```

### Implementation

**2-Gene Intersection (2i)**:
```python
mcp__ultraquery-inference__answer_complex_query(
    query_structure=[
        ["GENE_A", ["associated with"]],
        ["GENE_B", ["associated with"]]
    ],
    top_k=25
)
```

**3-Gene Intersection (3i)**:
```python
mcp__ultraquery-inference__answer_complex_query(
    query_structure=[
        ["GENE_A", ["associated with"]],
        ["GENE_B", ["associated with"]],
        ["GENE_C", ["associated with"]]
    ],
    top_k=25
)
```

### Characteristics

| Aspect | Description |
|--------|-------------|
| **Operation** | Logical intersection (AND) in query space |
| **Effect** | Returns entities satisfying ALL conditions |
| **Synergy Detection** | Finds shared associations naturally |
| **Dilution Risk** | More entities = stricter filter = fewer results |

### When to Use Intersection

Best for:
- Finding diseases affected by MULTIPLE genes/pathways
- Identifying shared therapeutic targets
- Bispecific/trispecific antibody target validation
- Pathway convergence analysis

### Comparison Table

| Aspect | Mean Embeddings | Intersection |
|--------|-----------------|--------------|
| **Question Answered** | "What is the combined signature associated with?" | "What is associated with ALL entities?" |
| **Mathematical** | Centroid calculation | Set intersection |
| **More Genes** | Averages toward center | Stricter filtering |
| **Synergy** | Detectable via Δ calculation | Implicit in query structure |
| **Best For** | Gene signatures, biomarkers | Shared targets, combinations |

---

## Calculating Synergy/Dilution

### From Part 1 Individual Ranks

For each entity in combo, extract individual RAW RANKS from Part 1 (1-indexed, 1=best):
```python
individual_ranks = [rank(GENE_A, DISEASE), rank(GENE_B, DISEASE)]
# Log-space geometric mean (numerically stable for large ranks × many genes):
geo_rank_mean = exp(mean(log(individual_ranks)))
# Convert to pct_rank space:
geo_pct_rank = 1.0 - (geo_rank_mean - 1) / total_entities
```

For N-ary combos, this generalizes naturally: `(rank_1 × ... × rank_N)^(1/N)`.

If any gene is missing from results, flag as `classification = INCOMPLETE` (don't guess or penalize).

### From Part 2 Combo Scores

Get the combo/signature score:
```
combo_pct_rank = pct_rank from signature query result
```

### Delta Calculation

```
delta = combo_pct_rank - geo_pct_rank
```

### Probability Mean (per-MCP only)

```python
# Mean in raw score space first, THEN transform:
probability_mean = sigmoid(mean(logit_A, logit_B, ...))  # ULTRA/UltraQuery
probability_mean = sigmoid(mean(cos_sim_A, cos_sim_B, ...))  # BioBridge
# Never combine probabilities across different MCPs (different raw score scales)
```

### Classification

| Delta | Classification | Interpretation |
|-------|----------------|----------------|
| > +0.02 | **SYNERGY** | Combo enhances association beyond components |
| -0.02 to +0.02 | **NEAR-ADDITIVE** | Combo ≈ geometric mean of components |
| < -0.02 | **DILUTION** | Combo weakens association vs components |
| N/A | **INCOMPLETE** | One or more genes missing from results |

Note: PrimeKG combo scoring cannot detect synergy (always equals mean by construction).

### Severity Scale

| Delta | Severity |
|-------|----------|
| > +0.05 | Strong synergy |
| +0.02 to +0.05 | Mild synergy |
| -0.02 to +0.02 | Near-additive |
| -0.05 to -0.02 | Mild dilution |
| -0.10 to -0.05 | Moderate dilution |
| < -0.10 | Severe dilution |

---

## Method Selection Guide

### Use Mean Embeddings (BioBridge) When:
- Creating gene expression signatures
- Building multi-gene biomarkers
- Want to capture "average" biology of gene set
- Interested in pathway centroid effects

### Use Intersection (UltraQuery) When:
- Finding shared therapeutic targets
- Validating bispecific/trispecific combinations
- Identifying diseases affected by multiple pathways
- Need logical AND relationships

### Use BOTH When:
- Comparing methods for robustness
- Comprehensive combo analysis
- Research/publication quality results

---

## Method 3: Shortest Path (PrimeKG)

### Concept

Measures structural proximity in the knowledge graph using BFS (Breadth-First Search):

```
score = 0.9^(path_length - 1)
```

Where:
- 1-hop (direct edge) = 1.0 (strongest)
- 2-hop = 0.9
- 3-hop = 0.81
- No path = 0.0 (disconnected)

### Implementation

**Individual Queries**:
```python
from query_primekg import PrimeKGData

primekg = PrimeKGData()

# Single query
result = primekg.get_shortest_path(
    source_entity="TYK2",          # Gene name or ID
    target_entity="Crohn disease", # Disease name or ID
    source_type="gene/protein",
    target_type="disease"
)

# Batch query (parallel processing)
batch_result = primekg.get_shortest_paths_batch(
    source_entities=["TYK2", "JAK1", "ITGA4"],
    target_entities=["Crohn disease", "ulcerative colitis"],
    source_type="gene/protein",
    target_type="disease",
    parallel=True,
    max_workers=4
)
```

**Combo Queries**:
```python
# Combo score = average of component scores from individual queries
combo_scores = []
for gene in ["TYK2", "JAK1"]:
    score = individual_scores[(gene, disease)]
    combo_scores.append(score)

combo_score = sum(combo_scores) / len(combo_scores)
```

### Characteristics

| Aspect | Description |
|--------|-------------|
| **Operation** | BFS traversal on unweighted graph |
| **Effect** | Measures structural/topological proximity |
| **Combo Method** | Average of component path scores |
| **Synergy Detection** | **Not applicable** - combo always equals arithmetic mean |

### Score Interpretation

| Path Length | Score | Interpretation |
|-------------|-------|----------------|
| 1 (direct) | 1.00 | Direct edge - very strong association |
| 2 | 0.90 | 1 intermediate - strong |
| 3 | 0.81 | 2 intermediates - moderate-strong |
| 4 | 0.729 | 3 intermediates - moderate |
| 5 | 0.656 | 4 intermediates - weak |
| 6+ | <0.59 | 5+ intermediates - very weak |
| No path | 0.00 | Disconnected - no association |

### When to Use Shortest Path

Best for:
- Validating direct vs indirect associations
- Understanding mechanistic paths between entities
- Complementing embedding-based predictions with graph structure
- Identifying intermediate nodes in pathways
- Network topology analysis

### Limitations

1. **No Synergy/Dilution**: Since combo = average, delta is always 0 by definition
2. **Path Quality**: Shortest path may not be most biologically relevant path
3. **Edge Weighting**: Currently unweighted - all edges treated equally
4. **Graph Coverage**: Limited to entities and relations in PrimeKG

---

## Three-Method Comparison Table

| Aspect | Mean Embeddings (BioBridge) | Intersection (UltraQuery) | Shortest Path (PrimeKG) |
|--------|---------------------------|--------------------------|------------------------|
| **What it measures** | Semantic similarity in embedding space | Logical AND satisfaction | Graph distance (hops) |
| **Combo calculation** | Mean of entity embeddings | Intersection query | Average of path scores |
| **Synergy detection** | Yes (Δ from mean) | Yes (Δ from mean) | **No** (combo = mean) |
| **Output type** | cos_sim, pct_rank | score, pct_rank | path_length, score |
| **Best for** | Gene signatures, biomarkers | Shared targets, combinations | Path validation, topology |
| **Captures** | Literature/annotation similarity | Multi-entity requirements | Structural connectivity |
| **Computational** | Fast (embedding lookup) | Moderate (query inference) | Variable (BFS traversal) |

---

## Cross-Method Agreement

### When All 3 Methods Agree

**HIGH CONFIDENCE** finding when:
- BioBridge shows strong cos_sim/pct_rank
- UltraQuery ranks entity highly in intersection results
- PrimeKG shows short path (1-3 hops)

### When Methods Disagree

| BioBridge | UltraQuery | PrimeKG | Interpretation |
|-----------|------------|---------|----------------|
| Strong | Strong | Weak (long path) | Semantic relationship, indirect mechanism |
| Weak | Weak | Strong (short path) | Direct edge exists but low annotation |
| Strong | Weak | Strong | One entity dominates, other has weaker evidence |
| Weak | Strong | Weak | May be artifact of intersection query |

### Recommended Workflow

For highest confidence results:

1. **Run all 3 methods** on same entity pairs
2. **Compare rankings** across methods
3. **Flag concordant findings** as high priority
4. **Investigate discordant findings** for biological insight
5. **Use PrimeKG paths** to understand mechanism of concordant findings

---

## Method Selection Summary

| Analysis Goal | Primary Method | Supporting Method |
|--------------|----------------|-------------------|
| Gene signature biomarker | BioBridge | UltraQuery |
| Multi-target drug | UltraQuery | BioBridge |
| Pathway mechanism | PrimeKG | BioBridge |
| Comprehensive analysis | **All 3** | - |
| Quick validation | BioBridge | - |
| Strict AND requirement | UltraQuery | - |
| Direct edge check | PrimeKG | - |
