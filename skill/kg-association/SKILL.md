---
name: kg-association
description: Analyze biomedical entity associations using multiple knowledge graph methods (BioBridge, ULTRA/ULTRAQuery, PrimeKG). Use this skill when users request multi-entity association analysis, gene-disease predictions, combo/signature analysis, or cross-method comparison. Generates structured markdown reports with biological interpretations for individual entities, combos, and comparative analysis across all three methods.
---

# Knowledge Graph Association Analysis

## Overview

This skill performs comprehensive biomedical entity association analysis using three complementary knowledge graph methods:

| Method | Approach | Best For | Detailed Docs |
|--------|----------|----------|---------------|
| **BioBridge** | Neural embeddings with cosine similarity | Semantic similarity, mean embeddings for combos | [`references/biobridge_analysis.md`](references/biobridge_analysis.md) |
| **ULTRA** | Foundation model for simple and complex logical queries on knowledge graphs | simple and complex reasoning queries | [`references/ultra_analysis.md`](references/ultra_analysis.md) |
| **PrimeKG** | Graph traversal (BFS shortest path) | Network topology, mechanistic paths | [`references/primekg_analysis.md`](references/primekg_analysis.md) |

**Output**: Structured markdown reports in `kgpred_<context>_<date>/` directory

---

## When to Use Each Method

### Decision Matrix

| Use Case | BioBridge | ULTRA | PrimeKG |
|----------|-----------|-------|---------|
| Single gene-disease association | ✓ Good | ✓ Good | ✓ Good |
| Gene combo/signature analysis | ✓ **Best** (mean embeddings) | ✓ Good (intersection) | Average only |
| Synergy/dilution detection | ✓ Yes | ✓ Yes | No |
| Complex multi-hop reasoning | Limited | ✓ **Best** | Limited |
| Network topology analysis | No | No | ✓ **Best** |
| Mechanistic path discovery | No | No | ✓ **Best** |
| Interactive visualization | No | No | ✓ Yes |

### Quick Selection Guide

```
"Find diseases associated with gene X"
  → Run all three methods for cross-validation

"Analyze gene combo for synergy"
  → BioBridge (mean embeddings) + ULTRA (intersection)

"How does gene X connect to disease Y?"
  → PrimeKG (shows intermediate nodes)

"Find diseases associated with BOTH gene X AND gene Y"
  → ULTRA (native intersection query)
```

---

## Required Inputs

### 1. Target Entities (entity1_list)

Entities to predict associations WITH (typically diseases):

```yaml
entity1_list:
  - Crohn's disease
  - Ulcerative colitis
  - Inflammatory bowel disease
entity1_type: disease  # Default
```

### 2. Source Entities (entity2_list)

Source entities in one of two formats:

**Format A: Flat list (individual entities only)**
```yaml
entity2_list:
  - TYK2
  - JAK1
  - GREM1
# → Only Part 1 (Individual Analysis) is generated
```

**Format B: Nested list (combos)**
```yaml
entity2_list:
  - [TYK2, JAK1]
  - [TNFRSF25, GREM1]
  - [CDKN2D, ITGA4, ITGB7]
# → Full analysis: Part 1, Part 2, Part 3
# → Individuals derived from union of combo elements
```

### Input Detection Logic

```
IF entity2_list contains nested lists:
    mode = "COMBO"
    → Generate Part 1, Part 2, Part 3 for each method
ELSE:
    mode = "INDIVIDUAL_ONLY"
    → Generate Part 1 only
```

---

## Output Structure

### Directory Layout

```
kgpred_<context_disease>/
├── KG_Association_Report.md           # Main report with pct_rank summary
├── biobridge/
│   ├── Part1_Individual_Analysis.md
│   ├── Part2_Combo_Analysis.md        # If combos provided
│   └── Part3_Comparative_Analysis.md  # If combos provided
├── ultra/
│   ├── Part1_Individual_Analysis.md
│   ├── Part2_Combo_Analysis.md        # If combos provided
│   ├── Part3_Comparative_Analysis.md  # If combos provided
│   ├── individual/*.parquet           # Individual query results
│   └── combo/*.parquet                # Combo query results
├── primekg/
│   ├── Part1_Individual_Analysis.md
│   ├── Part2_Combo_Analysis.md        # If combos provided
│   ├── Part3_Comparative_Analysis.md  # If combos provided
│   ├── {genes}_{disease}_network.html # Full interactive network (REQUIRED)
│   ├── {genes}_{disease}_subgraph.html # Per-combo subgraph (REQUIRED for combos)
│   └── primekg_shortest_paths.json    # Shortest path data
└── Cross_Method_Comparison.md
```

**REQUIRED Visualizations:**
- Always generate `primekg/{genes}_{disease}_network.html` with full gene-disease network
- For each combo, generate `primekg/{genes}_{disease}_subgraph.html`
- Use `pixi run python scripts/primekg_visualization.py` to create visualizations
- Example: `TYK2_JAK1_IBD_network.html`, `TYK2_JAK1_IBD_subgraph.html`

### Report Contents

| Report | Contents |
|--------|----------|
| **Part 1** | Individual entity pct_rank scores and assessments |
| **Part 2** | Combo pct_rank scores with synergy/dilution classification |
| **Part 3** | Deep comparative analysis with biological interpretation |
| **Cross_Method_Comparison** | Unified pct_rank comparison across all three methods |

### Required Report Structure

**Executive Summary MUST include pct_rank table:**
```markdown
| Method | Target1 → Disease (pct_rank) | Target2 → Disease (pct_rank) | Combo pct_rank |
|--------|------------------------------|------------------------------|----------------|
| BioBridge | 0.XXXX | 0.XXXX | 0.XXXX |
| ULTRA | 0.XXXX | 0.XXXX | 0.XXXX |
| PrimeKG | 0.XXXX | 0.XXXX | 0.XXXX |
```

**All tables MUST report pct_rank consistently across methods.**

---

## Scoring Quick Reference

### Score Interpretation (All Methods)

| Score Range | Assessment |
|-------------|------------|
| ≥ 0.99 | **Very Strong** |
| 0.95-0.99 | **Strong** |
| 0.90-0.95 | Moderate-Strong |
| 0.80-0.90 | Moderate |
| 0.50-0.80 | Weak |
| < 0.50 | Very Weak |

### Synergy/Dilution Thresholds

| Delta (Δ) | Classification |
|-----------|----------------|
| > +0.02 | **SYNERGY** |
| -0.02 to +0.02 | NEAR-ADDITIVE |
| < -0.02 | **DILUTION** |

### Cross-Method Agreement

| Agreement | Confidence |
|-----------|------------|
| All 3 methods agree | **HIGH** |
| 2 of 3 agree | **MODERATE** |
| All disagree | **LOW** - investigate |

---

## Workflow Summary

### Phase 0: Setup Output Directory

```bash
mkdir -p ./kgpred_<context_disease>/{biobridge,ultra,primekg}
```

### Phase 1: Individual Analysis (All Methods)

Query each gene-disease pair individually for BioBridge, ULTRA, and PrimeKG.

**See detailed workflows:**
- BioBridge: [`references/biobridge_analysis.md`](references/biobridge_analysis.md#individual-entity-analysis)
- ULTRA: [`references/ultra_analysis.md`](references/ultra_analysis.md#individual-entity-analysis-ultra-inference)
- PrimeKG: [`references/primekg_analysis.md`](references/primekg_analysis.md#running-primekg-analysis)

### Phase 2: Combo Analysis (If Combos Provided)

Query each combo-disease pair for synergy/dilution assessment.

**See detailed workflows:**
- BioBridge: [`references/biobridge_analysis.md`](references/biobridge_analysis.md#combosignature-analysis)
- ULTRA: [`references/ultra_analysis.md`](references/ultra_analysis.md#combointersection-analysis-ultraquery-inference)
- PrimeKG: [`references/primekg_analysis.md`](references/primekg_analysis.md#report-formats)

### Phase 3: Cross-Method Comparison

Generate unified comparison report highlighting agreements and discordances.

---

## Key Principles

1. **Run ALL THREE methods** for comprehensive cross-validated results
2. **Save outputs to working directory** - never to skill or MCP folders
3. **Move files IMMEDIATELY** after each MCP query
4. **BioBridge: ALWAYS use explicit overrides** - `override_head_name`, `override_tail_name`, etc.
5. **ULTRA: use percentile rank** - single queries using ultra-inference, use `percentile_rank`, complex queries (like intersections) use `filtered_percentile_rank`
6. **PrimeKG: No synergy detection** - combo = average of components; pct_rank = score (0.9^(path_length-1))
7. **Report pct_rank consistently** - ALL summary tables must include pct_rank for individual AND combo analyses
8. **Generate PrimeKG visualizations** - create network HTML files for all analyses using `scripts/primekg_visualization.py`
9. **Identify agreement** - high confidence when all methods concur
10. **Flag discordance** - investigate when methods disagree
11. **PrimeKG MCP: Use `limit=null`** - always get complete connection data (see [`references/primekg_analysis.md#mcp-query-best-practices`](references/primekg_analysis.md#mcp-query-best-practices))
12. **Prioritize shared connections** - shared disease genes and shared intermediary nodes get visualization priority boost

---

## Available Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `scripts/extract_ultra_individual.py` | Extract scores from individual parquets | See [`references/ultra_analysis.md`](references/ultra_analysis.md#helper-scripts) |
| `scripts/extract_ultra_combo.py` | Extract scores from combo parquets | See [`references/ultra_analysis.md`](references/ultra_analysis.md#helper-scripts) |
| `scripts/parquet_utils.py` | Parquet inspection and utilities | See [`references/ultra_analysis.md`](references/ultra_analysis.md#helper-scripts) |
| `scripts/primekg_visualization.py` | Network visualization | See [`references/primekg_analysis.md`](references/primekg_analysis.md#network-visualization) |

---

## Reference Documentation

| Document | Contents |
|----------|----------|
| [`references/biobridge_analysis.md`](references/biobridge_analysis.md) | BioBridge API, override requirements, mean embeddings, synergy detection |
| [`references/ultra_analysis.md`](references/ultra_analysis.md) | ULTRA tools, parquet schemas, intersection queries |
| [`references/primekg_analysis.md`](references/primekg_analysis.md) | Shortest path scoring, visualization, path interpretation |
| [`references/mcp_usage.md`](references/mcp_usage.md) | MCP tool documentation |
| [`references/report_templates.md`](references/report_templates.md) | Markdown report templates |
| [`references/combo_methods.md`](references/combo_methods.md) | Mean embeddings vs intersection comparison |
| [`references/interpretation_guide.md`](references/interpretation_guide.md) | Score interpretation guidelines |

---

## Quick Start Example

```python
# 1. Define inputs
genes = ["TYK2", "JAK1"]
diseases = ["Crohn disease", "ulcerative colitis", "inflammatory bowel disease"]
combos = [["TYK2", "JAK1"]]

# 2. Run BioBridge (see references/biobridge_analysis.md)
# Individual queries
for gene in genes:
    for disease in diseases:
        result = mcp__biobridge__predict_associations(
            context=f"Association: {gene} - {disease}",
            override_head_name=gene,
            override_head_type="gene/protein",
            override_tail_name=disease,
            override_tail_type="disease",
            topk=1
        )
        # Record: pct_rank = result['results'][0]['pct_rank']

# Combo query (mean embedding)
for combo in combos:
    result = mcp__biobridge__predict_associations(
        context=f"Combo association: {'+'.join(combo)} - {disease}",
        override_head_names=combo,  # LIST of genes
        override_head_type="gene/protein",
        override_tail_name=disease,
        override_tail_type="disease",
        topk=1
    )
    # Record: combo_pct_rank = result['results'][0]['pct_rank']

# 3. Run ULTRA (see references/ultra_analysis.md)
# Individual queries
for gene in genes:
    result = mcp__ultra-inference__predict_tail_entities(
        head_entity=gene,
        relation="associated with"
    )
    # Move parquet, read percentile_rank for target diseases

# Combo query (2i intersection)
for combo in combos:
    result = mcp__ultraquery-inference__answer_complex_query(
        query_structure=[[combo[0], ["associated with"]], [combo[1], ["associated with"]]],
        top_k=50
    )
    # Record: filtered_percentile_rank for target diseases

# 4. Run PrimeKG (see references/primekg_analysis.md)
# Use MCP tools for shortest paths
# Save results to primekg/primekg_shortest_paths.json

# 5. Generate PrimeKG visualizations (REQUIRED)
# Full network:
pixi run python scripts/primekg_visualization.py \
    --json-results ./primekg/primekg_shortest_paths.json \
    --output ./primekg/TYK2_JAK1_IBD_network.html \
    --title "TYK2+JAK1 IBD Network"

# Combo subgraph:
pixi run python scripts/primekg_visualization.py \
    --json-results ./primekg/primekg_shortest_paths.json \
    --combo "TYK2+JAK1" \
    --output ./primekg/TYK2_JAK1_IBD_subgraph.html

# 6. Generate KG_Association_Report.md with pct_rank summary table
# 7. Generate Cross_Method_Comparison.md
```
