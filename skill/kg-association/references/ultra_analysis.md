# ULTRA Association Analysis

# ULTRA Association Analysis

## Overview

ULTRA (Unified, Learnable, and TRAnsferable) is a foundation model for knowledge graph reasoning that can generalize to any knowledge graph without retraining. It uses several GNN-based architecture layers to learn relational patterns and predict missing links.

**Key Characteristics:**
- Foundation model approach (zero-shot generalization)
- Supports complex logical queries (intersection, union, negation)
- Two MCP tools: `ultra-inference` (individual) and `ultraquery-inference` (combos)
- Returns predictions ranked by model confidence scores


## When to Use ULTRA

| Use Case | Recommendation |
|----------|----------------|
| Single gene-disease association | ✓ Good choice (use ultra-inference) |
| Gene combo analysis | ✓ Excellent - uses intersection queries (ultraquery-inference) |
| Logical AND queries | ✓ Best choice |
| Complex multi-hop reasoning | ✓ Best choice |
| Multi-way intersections (3+ genes) | ✓ Supported |
| Semantic similarity | Use BioBridge instead |
| Network path analysis | Use PrimeKG instead |

---

## Scoring Methodology

### Primary Score
Use Percentile Rank (percentile_rank) in ULTRA and Filtered percentile rank (filtered_percentile_rank) in ULTRA Query

| pct_rank | Interpretation | Assessment |
|----------|---------------|------------|
| ≥ 0.99 | Top 1% among target type | **Very Strong** |
| 0.95-0.99 | Top 5% | **Strong** |
| 0.90-0.95 | Top 10% | Moderate-Strong |
| 0.80-0.90 | Top 20% | Moderate |
| 0.50-0.80 | Top 50% | Weak |
| < 0.50 | Below median | Very Weak |

---

## Parquet Column Schemas

### ultra-inference Parquet Schema (Individual Queries)

| Column | Description | Example |
|--------|-------------|---------|
| `h_label` | Head entity ID | "NCBI:7297" |
| `h_name` | Head entity name | "TYK2" |
| `t_pred_label` | Tail entity ID | "MONDO:5011" |
| `t_pred_name` | Tail entity name | "Crohn disease" |
| `t_pred_score` | Model prediction score | 0.847 |
| `t_pred_type` | Tail entity type | "disease" |
| `rank` | Rank within schema-matched (USE) | 1523 |
| `percentile_rank` | Percentile within schema-matched (USE) | 0.988 |

### ultraquery-inference Parquet Schema (Combo Queries)

| Column | Description | Example |
|--------|-------------|---------|
| `entity_id` | Entity ID | "MONDO:5011" |
| `entity_name` | Entity name | "Crohn disease" |
| `entity_type` | Entity type | "disease" |
| `score` | Model prediction score | 0.923 |
| `rank` | **Global** rank (DO NOT USE) | 847 |
| `percentile_rank` | **Global** percentile (DO NOT USE) | 0.993 |
| `filtered_rank` | Rank within schema-matched (USE) | 12 |
| `filtered_percentile_rank` (USE) | Percentile within schema-matched | 0.9996 |
---

## Individual Entity Analysis (ultra-inference)

### MCP Tool

Use `mcp__ultra-inference__predict_tail_entities` for individual gene queries.

### Required Parameters

| Parameter | Required | Description | Example |
|-----------|----------|-------------|---------|
| `head_entity` | Yes | Gene name or ID | "TYK2" or "NCBI:7297" |
| `relation` | Yes | Relation type | "associated with" |
| `top_k` | Optional | Limit results (None = all) | None |

### Individual Query Code Example

```python
import pandas as pd
import shutil
from pathlib import Path

working_dir = Path("./kgpred_IBD_2025-01-08")
data_dir = working_dir / "ultra" / "individual"
data_dir.mkdir(parents=True, exist_ok=True)

results = {}

for gene in genes:
    # Step 1: Run MCP query
    result = mcp__ultra-inference__predict_tail_entities(
        head_entity=gene,
        relation="associated with",
        top_k=None  # Get ALL predictions for proper pct_rank calculation
    )

    # Step 2: IMMEDIATELY move parquet to working directory
    src_file = result['output_file']
    dest_file = data_dir / f"{gene}_predictions.parquet"
    shutil.move(src_file, dest_file)

    # Step 3: Extract scores for target diseases
    results[gene] = {}
    for disease in target_diseases:
        match = df_disease[df_disease['t_pred_name'].str.lower().str.contains(disease.lower())]
        if not match.empty:
            results[gene][disease] = {
                'pct_rank': match['pct_rank'].iloc[0],
                'rank': match['rank'].iloc[0],
                'raw_score': match['t_pred_score'].iloc[0]
            }
```

### Individual Analysis Report Format

```markdown
## ULTRA Individual Gene-Disease Analysis

### Crohn Disease

| Gene | Score (pct_rank) | Rank (disease-only) | Raw Score | Assessment |
|------|------------------|---------------------|-----------|------------|
| TYK2 | 0.9985 | 26 | 0.847 | Very Strong |
| JAK1 | 0.9992 | 14 | 0.869 | Very Strong |
| ITGA4 | 0.9997 | 5 | 0.912 | Very Strong |

```

---

## Combo Analysis (ultraquery-inference)

### MCP Tool

Use `mcp__ultraquery-inference__answer_complex_query` for combo/intersection queries.

### Intersection Query Structure

ULTRA uses the BetaE query format for logical operations:

**2-way intersection (2i):**
```python
query = [
    ["TYK2", ["associated with"]],
    ["JAK1", ["associated with"]]
]
# Finds entities associated with BOTH TYK2 AND JAK1
```

**3-way intersection (3i):**
```python
query = [
    ["ITGA4", ["associated with"]],
    ["ITGB7", ["associated with"]],
    ["TYK2", ["associated with"]]
]
# Finds entities associated with ALL THREE genes
```

**N-way intersection (ni):**
```python
# Supports any number of branches
query = [
    ["Gene1", ["associated with"]],
    ["Gene2", ["associated with"]],
    ["Gene3", ["associated with"]],
    ["Gene4", ["associated with"]]
]
```

### Combo Query Code Example

```python
import pandas as pd
import shutil
from pathlib import Path

working_dir = Path("./kgpred_IBD_2025-01-08")
data_dir = working_dir / "ultra" / "combo"
data_dir.mkdir(parents=True, exist_ok=True)

combo_results = {}

for combo in combos:  # e.g., [["TYK2", "JAK1"], ["ITGA4", "ITGB7", "TYK2"]]
    combo_name = "+".join(combo)

    # Step 1: Build intersection query
    query = [[gene, ["associated with"]] for gene in combo]

    # Step 2: Run MCP query
    result = mcp__ultraquery-inference__answer_complex_query(
        query_structure=query,
        top_k=25  # Only affects JSON response; parquet has ALL predictions
    )

    # Step 3: IMMEDIATELY move parquet files to working directory
    combo_filename = "_".join(combo)
    shutil.move(result['output_file_filtered'],
                data_dir / f"{combo_filename}_filtered.parquet")
    shutil.move(result['output_file_all'],
                data_dir / f"{combo_filename}_all.parquet")

    # Step 4: Read and process the filtered predictions
    df = pd.read_parquet(data_dir / f"{combo_filename}_filtered.parquet")

    # Step 5: Extract scores and calculate synergy/dilution
    combo_results[combo_name] = {}
    for disease in target_diseases:
        match = df_disease[df_disease['entity_name'].str.lower().str.contains(disease.lower())]
        if not match.empty:
            combo_score = match['pct_rank'].iloc[0]

            # Get individual scores from Part 1 results
            individual_scores = [
                individual_results[gene].get(disease, {}).get('pct_rank', 0)
                for gene in combo
            ]
            individual_mean = sum(individual_scores) / len(individual_scores) if individual_scores else 0
            delta = combo_score - individual_mean

            # Classify
            if delta > 0.02:
                classification = "SYNERGY"
            elif delta < -0.02:
                classification = "DILUTION"
            else:
                classification = "NEAR-ADDITIVE"

            combo_results[combo_name][disease] = {
                'combo_score': combo_score,
                'rank': match['rank'].iloc[0],
                'individual_mean': individual_mean,
                'delta': delta,
                'classification': classification
            }
```

### Important: top_k Parameter Behavior

The `top_k` parameter in ultraquery-inference:
- **Only affects the JSON API response** (limits preview)
- **Does NOT affect parquet output files** (always contain ALL predictions)
- For ultra-inference always use default null to get complete data

### Combo Analysis Report Format

```markdown
## ULTRA Combo Analysis (Intersection Queries)

### Crohn Disease

| Combo | Query Type | Combo Score | Rank | Individual Mean | Delta | Classification |
|-------|------------|-------------|------|-----------------|-------|----------------|
| TYK2+JAK1 | 2i | 0.9990 | 17 | 0.9989 | +0.0001 | NEAR-ADDITIVE |
| ITGA4+ITGB7 | 2i | 0.9995 | 8 | 0.9982 | +0.0013 | NEAR-ADDITIVE |
| ITGA4+ITGB7+TYK2 | 3i | 0.9988 | 21 | 0.9991 | -0.0003 | NEAR-ADDITIVE |

**Interpretation:**
- 2i query finds diseases at the INTERSECTION of both gene associations
- 3i query is more restrictive (must be associated with ALL three genes)
- Near-additive results suggest overlapping pathway biology
```

---

## Move-Then-Process Workflow

**CRITICAL: Follow this exact pattern for every ULTRA query.**

```
FOR EACH gene/combo:
    1. Run MCP query
    2. IMMEDIATELY move parquet file(s) to working directory
    3. Read the MOVED file (not original path)
    4. Filter to target entity type
    5. Extract scores
    6. Proceed to next query
```

### Why Immediate Move Matters

1. **Volatility**: MCP output files may be overwritten by next query
2. **Accessibility**: Working directory is always accessible to user
3. **Traceability**: Each query's output is preserved with descriptive name
4. **Reproducibility**: All data files in one predictable location

---

## Synergy/Dilution Calculation

Same methodology as BioBridge:

```python
individual_mean = mean(pct_rank scores of combo components from Part 1)
combo_score = pct_rank from intersection query (recalculated)
delta = combo_score - individual_mean

if delta > 0.02: "SYNERGY"
elif delta < -0.02: "DILUTION"
else: "NEAR-ADDITIVE"
```

### Interpretation for ULTRA

| Delta (Δ) | Classification | Interpretation |
|-----------|----------------|----------------|
| > +0.02 | **SYNERGY** | Intersection query finds stronger disease association than individual queries. The genes may converge on the same disease pathways. |
| -0.02 to +0.02 | **NEAR-ADDITIVE** | Intersection performs as expected. Genes have overlapping but not synergistic associations. |
| < -0.02 | **DILUTION** | Intersection is more restrictive and reduces the score. One gene may be filtering out relevant diseases. |

---

## Output Files

ULTRA analysis generates the following files in `<output_dir>/ultra/`:

| File | Generated When | Contents |
|------|----------------|----------|
| `Part1_Individual_Analysis.md` | Always | Individual gene-disease scores (recalculated) |
| `Part2_Combo_Analysis.md` | If combos provided | Intersection query scores with synergy/dilution |
| `Part3_Comparative_Analysis.md` | If combos provided | Deep comparative analysis |

Data files in `<output_dir>/ultra/`:

| File Pattern | Contents |
|--------------|----------|
| `individual/<gene>_predictions.parquet` | Full predictions for each gene |
| `combo/<genes>_filtered.parquet` | Filtered predictions for intersection |
| `combo/<genes>_all.parquet` | All predictions for intersection |

---

## Available Scripts

### scripts/extract_ultra_individual.py

Extract scores from individual gene prediction parquets:

```bash
python scripts/extract_ultra_individual.py \
    --data-dir ./kgpred_IBD/ultra/individual \
    --diseases "crohn,ulcerative colitis,IBD" \
    --output ./individual_scores.json
```

### scripts/extract_ultra_combo.py

Extract scores from combo/intersection parquets:

```bash
python scripts/extract_ultra_combo.py \
    --data-dir ./kgpred_IBD/ultra/combo \
    --diseases "crohn,ulcerative colitis,IBD" \
    --output ./combo_scores.json
```

### scripts/parquet_utils.py

Utility functions for inspection and quick lookups:

```bash
# Inspect structure
python scripts/parquet_utils.py inspect ultra/TYK2_predictions.parquet

# Search for disease
python scripts/parquet_utils.py search ultra/TYK2_predictions.parquet "crohn"

# Get single score
python scripts/parquet_utils.py score ultra/TYK2_predictions.parquet "crohn"
```

---

## Troubleshooting

### Common Issues

**Issue: pct_rank values seem wrong (too high or too low)**
- Ensure ultra-inference is using `percentile_rank`
- Ensure ultraquery-inference uses `filtered_percentile_rank`

**Issue: Disease not found in predictions**
- Check disease name spelling (case-insensitive search recommended)
- Disease may use different naming (e.g., "Crohn disease" vs "Crohn's disease")
- Try partial matching with `.str.contains()`
- Query using the primekg identifiers (e.g. NCBI:5247)

**Issue: Parquet file not found**
- Ensure you moved the file IMMEDIATELY after query
- Check the correct output_file key (`output_file` vs `output_file_filtered`)
- Verify destination path exists

**Issue: Intersection query returns empty or very few results**
- More genes = more restrictive query
- Some gene combinations may have no overlapping associations
- Try 2-way intersections before 3+ way

### Validation Checklist

- [ ] parquet files moved IMMEDIATELY after each query
- [ ] Using correct column names for each MCP tool
- [ ] Using correct name or identifier for each MCP tool
- [ ] Synergy/dilution calculated from percentile_rank/filtered_percentile_rank for ULTRA/ULTRA Query
- [ ] top_k=None for individual queries (need all predictions)
- [ ] Query structure correct for intersection queries