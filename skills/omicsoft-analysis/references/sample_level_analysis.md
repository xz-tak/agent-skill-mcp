# Sample-Level Expression Analysis

## Overview

Three analysis modules that operate on sample-level expression data from SOMA experiments.
Produces per-study outputs with interactive plots and statistical summaries.

## Data Flow

1. `soma_expr_extract.py` queries SOMA (server-side filtered by genes + studies)
2. TSV piped to R via stdout
3. R partitions by study, applies configs, runs modules
4. Outputs: PNG + HTML (plotly) + CSV per study

## Module: comparison (Expression Boxplots)

### Internal Studies

- Grouped boxplot of target gene expression (one boxplot per gene)
- Groups defined by study-specific config (metadata column + levels)
- Significance asterisks from DEG SOMA stats (study-specific thresholds)
- Multi-view support (Yokohama: fibrosis view, nash view, diagnosis view)
- Faceted plots (Varsity: by treatment; SPARC: by tissue)

### External Studies (Omicsoft/Curated)

- One plot PER comparison (Control vs Case from comparison_group)
- `assign_comparison_roles()` parses `@case`/`@control` patterns
- External threshold: padj < 0.05

### Functions

- `build_expr_plot()` — ggplot2 boxplot with jitter, custom colors, faceting
- `build_standard_sig_annotations()` — black asterisks for Engitix/SPARC
- `build_varsity_sig_annotations()` — per-facet dual-threshold annotations
- `build_yokohama_sig_annotations()` — colored midpoint asterisks per comparison
- `save_plot_pair()` — PNG (300dpi) + interactive HTML (plotly with hover)
- `enrich_hover_text()` — adds log2FC, pvalue, FDR, comparison info to plotly hover

## Module: gsva (GSVA Signature Scoring)

### Method

- GSVA package: `gsvaParam()` with `kcdf="Gaussian"` (for normalized expression)
- Limma: design matrix + contrasts + eBayes
- Significance: adj.P.Val < 0.05 (Benjamini-Hochberg FDR)
- All studies use same GSVA threshold (unified)

### Internal Studies

- Standard: all groups vs control (Engitix, SPARC, Yokohama)
- Varsity: per-treatment limma (`run_varsity_gsva`)
- Study-specific signatures (from config `study_signatures`) added to user signatures

### External Studies

- Case vs Control per comparison
- `run_gsva_analysis()` with limma

### Functions

- `run_gsva_analysis()` — standard GSVA + limma for single-control designs
- `build_gsva_plot()` — ggplot2 boxplot with signature-level annotations

## Module: corr (Correlation Analysis)

### Module 1: Gene-Gene Correlation

- Pairwise Spearman correlation of target genes within each group
- Faceted heatmap (bwr colormap, -1 to 1 scale)
- r value + asterisk annotation per cell
- Minimum: 5 samples per group

### Module 2: Target vs Continuous Metadata

- Scatter plot of each target gene vs each `continuous_cols` variable
- Linear regression + 95% CI (geom_smooth)
- Spearman rho annotation
- Skip if: continuous_cols empty, column all NA, column zero variance, <3 non-NA values

### Module 3: GSVA vs Continuous Metadata

- Same as Module 2 but using GSVA scores instead of gene expression
- Requires GSVA module to have run first (uses its output)

### Functions

- `compute_pairwise_cor()` — cor.test() with zero-variance check
- `format_cor_annotation()` — "0.72***" style
- `run_gene_gene_correlation()` — per-group heatmaps
- `run_target_vs_continuous()` — scatter + lm fit
- `run_gsva_vs_continuous()` — scatter + lm fit for GSVA scores
- `save_cor_plot()` — PNG + HTML

## Graceful Failure Rules

| Condition | Action |
|-----------|--------|
| Gene not in EXPR var | Warn, skip gene, continue |
| Group n<3 samples | Skip limma for that group, log |
| Continuous col all NA | Skip that column, log |
| Continuous col zero variance | Skip that column, log |
| GSVA signature <3 genes matched | Warn, skip signature |
| Study 0 samples | Skip entire study, log |
| DEG comparison not found | Plot without asterisks, note in log |
| metadata_* column missing for derivation | Skip that grouping view, log |

## Caveats

- External studies: no cross-study pooling; per-comparison only
- Correlation only meaningful with sufficient samples and variance
- Different platforms use different normalization — absolute values not comparable across studies
- Varsity derived columns depend on exact metadata_* column names from SOMA
- Yokohama NAS thresholds: RNA uses NAS>=4 + stage>=2; Protein uses NAS>=4 + Fibrosis in F2-F4

## Internal Study Configurations

Configuration is externalized to `references/internal_study_configs.json`. Each study has:
- `group_col`: metadata column defining sample groups
- `group_levels`: ordered list of group labels
- `control`: which group is the baseline (null for Varsity)
- `group_colors`: named map of group -> color
- `annotation_map`: comparison_id -> target group for DEG asterisks
- `facet_col`: optional column for plot faceting
- `continuous_cols`: metadata columns for correlation analysis
- `sig_col`: which DEG stat column to use for significance (padj, pval)
- `sig_threshold`: cutoff for significance annotation

Multi-view studies (Yokohama_RNA, Yokohama_Protein) use a `groupings` key containing
nested views (fibrosis, nash, diagnosis), each with their own group definitions and
derivation logic.

## Output Structure

```
sample_level/
├── manifest.json
├── Engitix_FFPE/
│   ├── target_expression.png / .html
│   ├── signature_gsva.png / .html
│   ├── target_comparison_stats.csv
│   ├── gsva_comparison_stats.csv
│   ├── analysis_log.txt
│   └── correlations/
│       ├── gene_gene/
│       ├── target_vs_clinical/
│       └── gsva_vs_clinical/
├── SPARC/
│   ├── CD/
│   └── UC/
├── Varsity/
├── Yokohama_RNA/
│   ├── fibrosis/
│   ├── nash/
│   └── diagnosis/
├── Yokohama_Protein/
│   ├── fibrosis/
│   └── nash/
└── [external_study]/
    └── [comparison_name]/
```
