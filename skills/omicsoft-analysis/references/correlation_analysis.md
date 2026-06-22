# Correlation Analysis Module Documentation

## Overview

Three correlation sub-modules that quantify relationships between target genes,
GSVA signatures, and continuous clinical metadata.

## Statistical Methods

### Spearman Rank Correlation (Primary)

- Non-parametric: no assumption of linearity or normal distribution
- Robust to outliers
- Reports: rho (effect size), p-value (significance)
- Used for: gene-gene heatmaps, annotation text

### Pearson Linear Correlation (Secondary)

- Reports: r (linear correlation coefficient), p-value
- Used for: target vs continuous scatter plots (both reported in CSV)
- Linear regression line + 95% CI shown on plots

## Minimum Sample Requirements

| Analysis | Minimum n | Rationale |
|----------|-----------|-----------|
| Gene-gene heatmap | n >= 5 per group | Statistical power for rank test |
| Target vs continuous | n >= 5 overall | Meaningful regression |
| GSVA vs continuous | n >= 5 overall | Meaningful regression |

## Zero-Variance Handling

Before computing correlation:
1. Check `sd(x) > 0` and `sd(y) > 0`
2. If either variable has zero variance (all identical values): return NA, skip pair
3. Log reason to analysis_log.txt

This prevents `cor.test()` from producing warnings or nonsensical results.

## Multiple Testing

- **No FDR correction** on correlation p-values
- Rationale: these are descriptive/exploratory, not inferential
- Asterisk notation indicates raw p-value significance level
- Users should interpret patterns rather than individual p-values

## Asterisk Notation

| Symbol | p-value range |
|--------|---------------|
| `*` | < 0.05 |
| `**` | < 0.01 |
| `***` | < 0.001 |
| `****` | < 0.0001 |

## Plot Aesthetics

### Gene-Gene Heatmap

- Color scale: blue-white-red (bwr), centered at 0
- Range: fixed -1 to +1
- Cell annotations: "0.72***" format (r value + asterisk)
- Faceted by condition group (and optionally by treatment)
- Grey borders between cells for readability

### Scatter Plots (Module 2 & 3)

- Points: blue (#2166AC), alpha=0.6, size=2
- Regression line: red (#B2182B), linewidth=0.8
- 95% CI band: light red (#FDDBC7)
- Annotation text: top-left of each panel (rho, p, n)
- Faceted by gene or signature (free_y scales)

## Edge Cases

| Case | Behavior |
|------|----------|
| All values identical in one group | Skip that group, log |
| < 2 target genes available | Skip gene-gene module entirely |
| Continuous column >50% non-numeric | Skip with warning |
| Single sample in a group | Skip correlation for that group |
| Gene not in expression matrix | Skip silently, report in summary |
| NaN/Inf in expression | Treated as NA, excluded from computation |

## Output Files

### Per-Study Correlation Outputs

```
correlations/
├── gene_gene/
│   ├── {study}_gene_gene_correlation.png
│   ├── {study}_gene_gene_correlation.html
│   └── {study}_gene_gene_correlation.csv
├── target_vs_clinical/
│   ├── {study}_target_vs_{var}.png
│   ├── {study}_target_vs_{var}.html
│   └── {study}_target_vs_clinical_stats.csv
└── gsva_vs_clinical/
    ├── {study}_gsva_vs_{var}.png
    ├── {study}_gsva_vs_{var}.html
    └── {study}_gsva_vs_clinical_stats.csv
```

### CSV Column Descriptions

**gene_gene_correlation.csv:**
- gene1, gene2: gene pair
- condition: group label
- treatment: facet label (NA if no faceting)
- spearman_r: Spearman correlation coefficient
- pvalue: raw p-value
- n: number of samples used

**target_vs_clinical_stats.csv:**
- gene: target gene
- continuous_var: metadata column name
- spearman_rho, spearman_pval: Spearman results
- pearson_r, pearson_pval: Pearson results
- n: number of valid pairs

**gsva_vs_clinical_stats.csv:**
- signature: GSVA signature name
- continuous_var: metadata column name
- spearman_rho, spearman_pval, pearson_r, pearson_pval, n
