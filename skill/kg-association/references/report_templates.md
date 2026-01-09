# Report Templates

## Executive Summary Template (REQUIRED for all reports)

The main report (`KG_Association_Report.md`) MUST begin with this pct_rank summary:

```markdown
# KG Association Analysis Report: {Combo} for {Disease}

**Date:** {YYYY-MM-DD}
**Target Combination:** {Target1} + {Target2}
**Disease Indication:** {Disease}

---

## Executive Summary

### Key Findings Summary

| Method | {Target1} → {Disease} (pct_rank) | {Target2} → {Disease} (pct_rank) | Combo pct_rank |
|--------|----------------------------------|----------------------------------|----------------|
| **BioBridge** | {0.XXXX} | {0.XXXX} | {0.XXXX} |
| **ULTRA** | {0.XXXX} | {0.XXXX} | {0.XXXX} |
| **PrimeKG** | {0.XXXX} ({path_info}) | {0.XXXX} ({path_info}) | {0.XXXX} (avg) |

**Verdict:** {1-2 sentence summary of overall support}
```

### pct_rank Source by Method

| Method | Individual pct_rank Source | Combo pct_rank Source |
|--------|---------------------------|----------------------|
| BioBridge | `result['results'][0]['pct_rank']` | `result['results'][0]['pct_rank']` (mean embedding) |
| ULTRA | `percentile_rank` column in parquet | `filtered_percentile_rank` from 2i query |
| PrimeKG | `score = 0.9^(path_length-1)` | Average of component scores |

---

## Part 1: Individual Entity Analysis Template

```markdown
# Part 1: Individual {entity2_type}-{entity1_type} Association Analysis
## {MCP_name} Predictions

**Date**: {YYYY-MM-DD}
**Context**: {context description}
**Method**: {method description - e.g., "Neural link prediction with BioBridge embeddings"}

---

## Executive Summary

{2-4 bullet points summarizing key findings:}
- Strongest association: {entity2} with {entity1} (score: X.XXX)
- Pattern observed: {e.g., "UC > CD > IBD across all genes"}
- Notable findings: {any standout results}

---

## Summary Table: All Association Scores

| {entity2_type} | {entity1[0]} | {entity1[1]} | {entity1[2]} |
|----------------|--------------|--------------|--------------|
| | cos_sim / pct_rank | cos_sim / pct_rank | cos_sim / pct_rank |
| **{entity2[0]}** | {score} / {rank}% | {score} / {rank}% | {score} / {rank}% |
| **{entity2[1]}** | {score} / {rank}% | {score} / {rank}% | {score} / {rank}% |
{... repeat for all entity2}

---

## Rankings by Association Strength

### {entity1[0]} (e.g., "Ulcerative Colitis")
1. **{top_entity2}**: {score} ({rank}%) — {brief interpretation}
2. **{second_entity2}**: {score} ({rank}%)
{... top 5-8}

### {entity1[1]} (e.g., "Crohn's Disease")
{same format}

---

## Individual {entity2_type} Analysis

### 1. {entity2[0]} {optional: emoji indicator}

**Association Scores:**
- **{entity1[0]}**: {score} ({rank}%)
- **{entity1[1]}**: {score} ({rank}%)
- **{entity1[2]}**: {score} ({rank}%)

**Interpretation:**
{2-3 sentences explaining the scores in biological context}

**Biological Context:**
- {bullet point about mechanism}
- {bullet point about pathway}
- {bullet point about therapeutic relevance}

**Clinical Relevance:** {1 sentence}

---

{... repeat for all entity2}

---

## Cross-{entity2_type} Pattern Analysis

### Pattern 1: {pattern name}
{Description of pattern observed across entities}

| {entity2_type} | {metric} | Interpretation |
|----------------|----------|----------------|
| {entity} | {value} | {interpretation} |

### Pattern 2: {pattern name}
{...}

---

## Key Biological Insights

### 1. {Insight title}
{Explanation}

### 2. {Insight title}
{Explanation}

---

## Clinical Recommendations

### For {use case 1}:
- {recommendation}

### For {use case 2}:
- {recommendation}

---

*Report generated: {date}*
```

---

## Part 2: Combo/Signature Analysis Template

```markdown
# Part 2: {entity2_type} Signature (Multi-{entity2_type}) Association Analysis
## {MCP_name} {method} Predictions

**Date**: {YYYY-MM-DD}
**Method**: {e.g., "Mean embeddings from multiple genes per signature"}
**Signatures Tested**: {N} combinations ({size} gene signatures)
**Comparison**: Signature vs. individual {entity2_type}s to identify synergy or dilution

---

## Executive Summary

{Key findings with synergy/dilution classification}

### Key Results:
- **[{combo1}]**: {classification} ({entity1[0]}: {score})
- **[{combo2}]**: {classification} ({entity1[0]}: {score})
{... for each combo}

---

## Summary Table: Signature Association Scores

| Signature | Genes | {entity1[0]} | {entity1[1]} | {entity1[2]} |
|-----------|-------|--------------|--------------|--------------|
| | | cos_sim / pct_rank | cos_sim / pct_rank | cos_sim / pct_rank |
| **[{combo1}]** | {N} | {score} / {rank}% | {score} / {rank}% | {score} / {rank}% |
{... for each combo}

---

## Signature {N}: [{combo}] - {descriptive name}

### Signature Scores
- **{entity1[0]}**: {score} ({rank}%)
- **{entity1[1]}**: {score} ({rank}%)
- **{entity1[2]}**: {score} ({rank}%)

### Comparison to Individual {entity2_type}s

| Metric | {component1} | {component2} | Mean | Signature | Δ from Mean |
|--------|--------------|--------------|------|-----------|-------------|
| **{entity1[0]}** | {score} | {score} | {mean} | **{sig_score}** | **{delta}** |
| **{entity1[1]}** | {score} | {score} | {mean} | **{sig_score}** | **{delta}** |
| **{entity1[2]}** | {score} | {score} | {mean} | **{sig_score}** | **{delta}** |

### Effect: **{SYNERGY/NEAR-ADDITIVE/DILUTION}**

**Interpretation:**
{Explanation of why this effect occurs}

**Biological Mechanism:**
{Pathway/mechanism explanation}

**Clinical Relevance:**
{Clinical implications}

---

{... repeat for each signature}

---

## Cross-Signature Performance Analysis

### Rankings by Association Strength

#### {entity1[0]} (Best Overall)
1. **[{combo}]**: **{score}** {emoji} {classification}
2. **[{combo}]**: {score}
{...}

---

### Synergy vs. Dilution Summary

| Signature | {entity1[0]} Effect | {entity1[1]} Effect | {entity1[2]} Effect | Overall |
|-----------|---------------------|---------------------|---------------------|---------|
| **[{combo}]** | {delta} {emoji} | {delta} | {delta} | **{classification}** |
{...}

---

## Key Findings

### 1. {Finding title}
{Explanation}

### 2. {Finding title}
{Explanation}

---

## Design Principles for Signatures

### DO:
1. {recommendation}
2. {recommendation}

### DON'T:
1. {anti-pattern}
2. {anti-pattern}

---

## Clinical Recommendations

### RECOMMENDED SIGNATURES:
**1. [{combo}] — {use case}**
- **Use for**: {application}
- **Score**: {entity1[0]} {score} ({classification})
- **Advantage**: {why it works}

### NOT RECOMMENDED:
**{combo}**
- {reason to avoid}

---

*Report generated: {date}*
```

---

## Part 3: Comparative Analysis Template

```markdown
# Part 3: Deep Comparative Analysis
## Mechanistic Insights on Synergy, Dilution, and Design Principles

**Date**: {YYYY-MM-DD}
**Focus**: Biological explanations for signature performance vs. individual {entity2_type}s

---

## Executive Summary

This analysis reveals **{N} critical principles** for signature design:

1. **{Principle 1}**: {brief description}
2. **{Principle 2}**: {brief description}
3. **{Principle 3}**: {brief description}

---

## Performance Matrix

### Comparison to Best Individual {entity2_type}

| Signature | Best Individual | Ind. Score | Sig. Score | Δ | Effect |
|-----------|----------------|------------|------------|---|--------|
| **[{combo}]** | {entity} | {score} | **{score}** | **{delta}** | {emoji} {class} |
{...}

### Comparison to Arithmetic Mean

| Signature | Mean | Signature | Δ from Mean | Interpretation |
|-----------|------|-----------|-------------|----------------|
| **[{combo}]** | {mean} | **{score}** | **{delta}** | **{class}** |
{...}

---

## Deep Dive {N}: [{combo}] — {title}

### {Synergy/Dilution} Evidence

| Indication | {comp1} | {comp2} | Mean | Signature | Δ | Severity |
|-----------|---------|---------|------|-----------|---|----------|
| **{entity1[0]}** | {score} | {score} | {mean} | **{sig}** | **{delta}** | {class} |
{...}

### Why {Synergy/Dilution} Occurs

#### 1. {Mechanism title}
{Detailed explanation with pathway diagram if helpful}

#### 2. {Mechanism title}
{Explanation}

### Clinical Translation

#### For Biomarker Development:
{Recommendations}

#### For Therapeutic Strategy:
{Recommendations}

---

{... repeat for key signatures}

---

## Cross-Signature Mechanistic Insights

### Principle 1: {title}
**Successful Example:** [{combo}]
{Explanation}

**Failed Example:** [{combo}]
{Explanation}

**Rule:** {Actionable rule}

---

### Principle 2: {title}
{Same structure}

---

## Therapeutic Implications

### 1. {Strategy title}
{Detailed recommendations}

### 2. {Strategy title}
{Detailed recommendations}

---

## Design Principles for Future Signatures

### Successful Signature Checklist
{checklist items}

### Warning Signs for Failed Signatures
{anti-patterns}

### Recommended Signatures to Test
{Future directions}

---

## Final Clinical Recommendations

### FOR BIOMARKER DEVELOPMENT:
**Use:**
1. [{combo}] — {use case}
2. [{combo}] — {use case}

**Avoid:**
- [{combo}] — {reason}

### FOR THERAPEUTIC DEVELOPMENT:
{Recommendations}

---

*Analysis complete. Date: {date}*
```

---

## PrimeKG Part 1: Individual Analysis Template

```markdown
# Part 1: Individual {entity2_type}-{entity1_type} Shortest Path Analysis
## PrimeKG Graph Predictions

**Date**: {YYYY-MM-DD}
**Context**: {context description}
**Method**: Shortest path via BFS on PrimeKG knowledge graph
**Scoring**: Exponential decay: score = 0.9^(path_length - 1)

---

## Executive Summary

{2-4 bullet points summarizing key findings:}
- Strongest association: {entity2} with {entity1} (score: X.XX, path: N hops)
- Direct associations found: N of M pairs have 1-hop path
- Notable findings: {any standout results}

---

## Summary Table: All Shortest Path Scores

| {entity2_type} | {entity1[0]} | {entity1[1]} | {entity1[2]} |
|----------------|--------------|--------------|--------------|
| | path / score | path / score | path / score |
| **{entity2[0]}** | {path} / {score} | {path} / {score} | {path} / {score} |
| **{entity2[1]}** | {path} / {score} | {path} / {score} | {path} / {score} |
{... repeat for all entity2}

---

## Path Length Distribution

| Path Length | Count | Percentage | Assessment |
|-------------|-------|------------|------------|
| 1 (direct) | N | X% | Direct association |
| 2 | N | X% | Strong |
| 3 | N | X% | Moderate |
| 4+ | N | X% | Weak |
| No path | N | X% | Disconnected |

---

## Individual {entity2_type} Analysis

### 1. {entity2[0]}

**Path Scores:**
- **{entity1[0]}**: {score} (path: {path_length} hops)
  - Path: {entity2[0]} → {intermediate} → {entity1[0]}
- **{entity1[1]}**: {score} (path: {path_length} hops)
- **{entity1[2]}**: {score} (path: {path_length} hops)

**Interpretation:**
{2-3 sentences explaining the scores in biological context}

**Clinical Relevance:** {1 sentence}

---

{... repeat for all entity2}

---

*Report generated: {date}*
```

---

## PrimeKG Part 2: Combo Analysis Template

```markdown
# Part 2: {entity2_type} Combo Shortest Path Analysis
## PrimeKG Average Score Predictions

**Date**: {YYYY-MM-DD}
**Method**: Average of individual shortest path scores
**Combos Tested**: {N} combinations

---

## Executive Summary

{Key findings}

---

## Summary Table: Combo Scores

| Combo | {entity1[0]} | {entity1[1]} | {entity1[2]} |
|-------|--------------|--------------|--------------|
| | avg score | avg score | avg score |
| **[{combo1}]** | {score} | {score} | {score} |
{... for each combo}

---

## Combo Analysis: [{combo}]

### Combo Score
- **{entity1[0]}**: {combo_score}

### Component Breakdown

| Component | Path Length | Score | Contribution |
|-----------|-------------|-------|--------------|
| {comp1} | {path} | {score} | {above/below avg} |
| {comp2} | {path} | {score} | {above/below avg} |

**Max Component**: {component} ({score})
**Min Component**: {component} ({score})
**Range**: {max - min}

**Interpretation:**
{Explanation of component contributions}

---

{... repeat for each combo}

---

*Report generated: {date}*
```

---

## PrimeKG Part 3: Comparative Analysis Template

```markdown
# Part 3: PrimeKG Comparative Analysis
## Path Structure Insights

**Date**: {YYYY-MM-DD}

---

## Combo vs Component Comparison

| Combo | Disease | Combo Score | Max Component | Min Component | Range |
|-------|---------|-------------|---------------|---------------|-------|
| [{combo}] | {disease} | {score} | {max_comp}: {score} | {min_comp}: {score} | {range} |

---

## Path Structure Analysis

### Direct Associations (1-hop)
{List of entity pairs with direct edges}

### Strong Associations (2-3 hops)
{List with interpretation}

### Weak/No Associations
{List with explanation}

---

## Network Topology Insights

{Patterns observed in path structures}

---

*Report generated: {date}*
```

---

## Cross-Method Comparison Template

```markdown
# Cross-Method Comparison Report
## BioBridge vs ULTRA vs PrimeKG Analysis

**Date**: {YYYY-MM-DD}
**Methods Compared**: BioBridge (embeddings), ULTRA (intersection), PrimeKG (shortest path)
**Analysis Type**: Comprehensive cross-method pct_rank comparison

---

## Executive Summary

{Overview of how all 3 methods compare}

---

## SECTION A: pct_rank Summary Tables (REQUIRED)

### A1. Individual Target pct_rank Comparison

| Disease | {Target1} BioBridge | {Target1} ULTRA | {Target1} PrimeKG | {Target2} BioBridge | {Target2} ULTRA | {Target2} PrimeKG |
|---------|---------------------|-----------------|-------------------|---------------------|-----------------|-------------------|
| {disease1} | **{pct_rank}** | **{pct_rank}** | **{pct_rank}** | **{pct_rank}** | **{pct_rank}** | **{pct_rank}** |
| {disease2} | {pct_rank} | {pct_rank} | {pct_rank} | {pct_rank} | {pct_rank} | {pct_rank} |

### A2. Combo pct_rank Comparison

| Disease | BioBridge Combo | ULTRA Combo (2i) | PrimeKG Combo (avg) |
|---------|-----------------|------------------|---------------------|
| {disease1} | **{pct_rank}** | **{pct_rank}** | **{pct_rank}** |
| {disease2} | {pct_rank} | {pct_rank} | {pct_rank} |

---

## SECTION B: Cross-Method Agreement Analysis

### B1. Individual Analysis Cross-Comparison

| {entity2} | {entity1} | BioBridge pct_rank | ULTRA pct_rank | PrimeKG pct_rank | Path | Agreement |
|-----------|-----------|-------------------|----------------|------------------|------|-----------|
| {gene} | {disease} | {pct_rank} | {pct_rank} | {pct_rank} | {hops} | {HIGH/MOD/LOW} |
{...}

**Concordant Findings** (all 3 agree):
{list}

**Discordant Findings**:
{list with explanations}

### B2. Combo Analysis Cross-Comparison

| Combo | Disease | BioBridge | ULTRA | PrimeKG | BB Class | ULTRA Class | Agreement |
|-------|---------|-----------|-------|---------|----------|-------------|-----------|
| [{combo}] | {disease} | {score} | {score} | {score} | {class} | {class} | {level} |
{...}

Note: PrimeKG combo = average (always NEAR-ADDITIVE)

### B3. Design Principles Cross-Comparison

**Unified Principles** (supported by all methods):
{list}

**Method-Specific Insights**:
- BioBridge: {insights}
- ULTRA: {insights}
- PrimeKG: {insights}

---

## SECTION C: Final Cross-Method Analysis

### C1. Confidence Assessment Matrix

| Finding | BioBridge | ULTRA | PrimeKG | Confidence |
|---------|-----------|-------|---------|------------|
| {finding} | {evidence} | {evidence} | {evidence} | {HIGH/MOD/LOW} |

### C2. Method Agreement Statistics

- All 3 agree: {N} findings ({X}%)
- 2 of 3 agree: {N} findings ({X}%)
- All disagree: {N} findings ({X}%)

### C3. Final Recommendations

**Target Prioritization** (by cross-method agreement):
1. {target} - HIGH confidence (all 3 methods)
2. {target} - MODERATE confidence (2/3 methods)

**Biomarker Combinations** (validated by multiple methods):
{list}

### C4. Methodological Recommendations

| Analysis Goal | Recommended Method | Reason |
|---------------|-------------------|--------|
| Semantic similarity | BioBridge | Embedding proximity |
| Logical AND queries | ULTRA | Intersection |
| Path analysis | PrimeKG | Graph structure |
| Highest confidence | All 3 | Cross-validation |

---

*Comparison generated: {date}*
```
