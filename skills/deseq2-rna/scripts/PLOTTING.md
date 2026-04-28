# Standalone Visualization Functions Module

Four standalone plotting functions for ad-hoc visualization of RNA-seq results. These are **not** part of the main DESeq2 workflow - they can be invoked separately when users request specific visualizations.

## When to Use

Use these functions when users request:
- Heatmaps comparing signatures across conditions
- Gene expression visualizations for specific gene sets
- GSVA/ssGSEA score comparisons with statistical testing
- Custom visualizations of DEG or enrichment results

## Available Functions

| Function | Purpose | Primary Inputs |
|----------|---------|----------------|
| `plot_signature_heatmap()` | NES heatmap (signatures × comparisons) | GSEA results file OR rerun GSEA |
| `plot_deg_heatmap()` | log2FC heatmap (genes × comparisons) | DEG summary stats file |
| `plot_expression_heatmap()` | Z-scored expression (genes × samples) | Expression matrix, metadata |
| `plot_gsva_boxplot()` | GSVA scores by treatment group | GSVA matrix, metadata OR rerun GSVA |

## Script Location

```
scripts/plotting_functions.R
```

## Quick Usage

```r
# Load functions
source("scripts/plotting_functions.R")

# Example: Signature NES heatmap
plot_signature_heatmap(
  gsea_data = "deg/RNAseq_gsea_all.txt",
  signatures = c("sig1", "sig2", "sig3"),
  output_file = "figures/custom/signature_heatmap.png",
  source_filter = "CUSTOM"
)
```

## Rerun GSEA/GSVA Option

Both `plot_signature_heatmap()` and `plot_gsva_boxplot()` support rerunning GSEA/GSVA analysis instead of using pre-computed results. This is useful when:
- You want to analyze new custom gene sets not in the original analysis
- You want to recalculate scores with different parameters
- You're working with data from a different source

**IMPORTANT:** GSEA pipeline saves ALL results (including non-significant) to output files. Dotplots filter for significance internally.

### Rerun GSEA Parameters

```r
plot_signature_heatmap(
  gsea_data = NULL,           # NULL when rerunning
  rerun_gsea = TRUE,          # Enable rerun mode
  deg_data = "path/to/summstats.txt",  # DEG results
  gene_sets = list(           # Named list of gene sets
    sig1 = c("GENE1", "GENE2"),
    sig2 = readLines("sig.txt")
  ),
  pval_cutoff = 1             # Include all results (default)
)
```

### Rerun GSVA Parameters

```r
plot_gsva_boxplot(
  gsva_data = NULL,           # NULL when rerunning
  rerun_gsva = TRUE,          # Enable rerun mode
  counts_data = "path/to/counts.rds",  # Expression/counts matrix
  gene_sets = list(           # Named list of gene sets
    sig1 = c("GENE1", "GENE2"),
    sig2 = readLines("sig.txt")
  ),
  kcdf = "auto"               # Auto-detect: "Poisson" for integers, "Gaussian" for normalized
)
```

**kcdf auto-detection:** The function automatically detects whether input is integer counts (uses "Poisson") or normalized values (uses "Gaussian") based on whether all values equal their floor.

## Required Inputs (User-Specified)

For each function, Claude should collect from user:

### plot_signature_heatmap
- **Use existing data OR rerun analysis** (ask first!)
  - If existing: GSEA results file path
  - If rerun: DEG summary stats file path + gene_sets (named list)
- Signature names to include
- Source filter (CUSTOM, GO:BP, H, C2, etc.)
- Output path

### plot_deg_heatmap
- DEG summary stats file path
- Signature file(s) OR gene list
- Output path

### plot_expression_heatmap
- Expression data (RDS or matrix file)
- Metadata file
- Signature file(s) OR gene list
- Group column name
- Sample ID column name
- Output path

### plot_gsva_boxplot
- **Use existing data OR rerun analysis** (ask first!)
  - If existing: GSVA data (RDS or matrix file)
  - If rerun: Counts data file path + gene_sets (named list)
- Metadata file
- Signature names
- Group column name
- Statistical method (adaptive/kruskal/none)
- Max pairs to display
- Output path

## AskUserQuestion Workflow for Rerun Options

When Claude is asked to generate GSVA boxplots or NES heatmaps, Claude MUST use `AskUserQuestion` tool:

**Question 1: Data Source**
```
header: "Data source"
question: "Use existing pre-computed results or rerun analysis?"
options:
  - label: "Use existing results (Recommended)"
    description: "Use pre-computed GSVA/GSEA scores from analysis RDS file"
  - label: "Rerun analysis"
    description: "Recompute GSVA/GSEA scores with new gene sets or matrix"
```

**If "Use existing results":**
- Read gsva_scores/gsea results from the RDS/text file
- No additional questions needed

**If "Rerun analysis" - Question 2 (GSVA only): Input Matrix**
```
header: "Input matrix"
question: "Which counts/expression matrix file should be used for GSVA?"
options: [list available .rds/.txt files in data directory]
```

**Note:** For GSVA, kcdf is auto-detected - no manual question needed.

## Output Files

All output files (static PNG, interactive HTML, and statistics TSV) are saved to subfolders under `figures/`:
- **Static PNG**: `figures/expr_comp_heatmap_boxplot/`
- **Interactive HTML**: `figures/interactive/expr_comp_heatmap_boxplot/`
- **Statistics TSV**: `figures/expr_comp_heatmap_boxplot/` (same as PNG)

| Function | PNG Output | HTML Output |
|----------|------------|-------------|
| plot_signature_heatmap | `*_NES_heatmap.png` | `*_NES_heatmap_interactive.html` |
| plot_deg_heatmap | `*_[signature]_deg_heatmap.png` | `*_[signature]_deg_heatmap_interactive.html` |
| plot_expression_heatmap | `*_[signature]_expression_heatmap.png` | `*_[signature]_expression_heatmap_interactive.html` |
| plot_gsva_boxplot | `*_[signature]_boxplot.png` | `*_[signature]_boxplot_interactive.html` |

**Additional Output for `plot_gsva_boxplot()`:**
- Pairwise statistics TSV: `figures/expr_comp_heatmap_boxplot/*_pairwise_stats.tsv` (when `export_stats = TRUE`)

## Interactive Plots

All plotting functions generate interactive HTML versions by default (`interactive = TRUE`). Interactive plots include hover tooltips with comprehensive statistical information:

| Function | Hover Information | Cell Annotation |
|----------|-------------------|-----------------|
| `plot_signature_heatmap()` | Signature name, Comparison, NES, p-value, FDR | `*` significance stars |
| `plot_deg_heatmap()` | Gene symbol, Comparison, log2FC, p-value (if available), FDR | `*` significance stars |
| `plot_expression_heatmap()` | Gene symbol, Sample, Treatment group, Z-score/expression value | None (colors show z-score) |
| `plot_gsva_boxplot()` | Treatment group, GSVA score, individual data points | N/A (boxplot) |

**Significance annotations (displayed in heatmap cells):**
- `*` p < 0.05
- `**` p < 0.01
- `***` p < 0.001
- `****` p < 0.0001

To disable interactive output, set `interactive = FALSE` in the function call.

**Note:** Interactive HTML plots require a web browser. They do NOT work in PDF, PowerPoint, or static documents.

## Dependencies

Additional R packages required:
- ComplexHeatmap
- ggsignif
- rstatix
- car
- heatmaply (for interactive heatmaps)
- plotly (for interactive plots)
- htmlwidgets (for saving HTML output)
- GSVA (for GSVA analysis)
- limma (for differential expression)
- ggrepel (for volcano plot labels)

---

## 5. run_gsva_limma_de() - GSVA Signature DE Analysis

### Purpose

Compute GSVA scores for custom gene signatures, run limma-based differential
expression across treatment groups, generate volcano plots and summary heatmap.

### When to Use

- Sample-level signature activity analysis (vs rank-based GSEA)
- Comparing signature scores across multiple treatment conditions
- Generating publication-ready volcano + heatmap visualizations
- When you need statistical testing between groups for signature scores

### Required Inputs

| Input | Description |
|-------|-------------|
| `counts_file` | Raw counts (genes × samples), CSV/TSV with Gene.name column or gene symbols as rownames |
| `metadata_file` | Sample metadata (TSV) with sample IDs as rownames and treatment column |
| `gmt_file` | GMT file with gene signatures (or named list of gene vectors) |
| `comparisons` | List of comparisons, each with `name`, `numerator`, `denominator` |

### Quick Usage

```r
source("~/.claude/skills/deseq2-rna/scripts/plotting_functions.R")

comparisons <- list(
  list(name = "TGFb_vs_Control", numerator = "TGFb", denominator = "Control"),
  list(name = "Drug_vs_TGFb", numerator = "TGFb+Drug", denominator = "TGFb")
)

run_gsva_limma_de(
  counts_file = "raw_counts.txt",
  metadata_file = "processed_metadata.txt",
  gmt_file = "custom_signatures.gmt",
  comparisons = comparisons,
  group_col = "Treatment",
  output_dir = "deseq2_output"
)
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `counts_file` | (required*) | Path to raw counts file (*not required in replot mode) |
| `metadata_file` | (required*) | Path to metadata file (*not required in replot mode) |
| `gmt_file` | (required*) | Path to GMT file or named list of gene vectors (*not required in replot mode) |
| `comparisons` | (required*) | List of comparison definitions (*not required in replot mode) |
| `group_col` | `"Treatment_std"` | Column in metadata for treatment groups |
| `output_dir` | `"."` | Base output directory |
| `padj_cutoff` | `0.05` | Adjusted p-value threshold for reporting |
| `pval_cutoff` | `0.05` | Nominal p-value threshold for volcano/heatmap significance stars |
| `cluster_rows` | `TRUE` | Cluster rows in heatmap (hierarchical clustering) |
| `cluster_columns` | `TRUE` | Cluster columns in heatmap |
| `column_order` | `NULL` | Custom column order (vector of comparison names) |
| `comparisons_filter` | `NULL` | Filter to include only these comparisons in heatmap |
| `row_order` | `NULL` | Custom row order (vector of signature names) |
| `column_labels` | `NULL` | Named vector mapping comparison names to display labels |
| `limma_results_file` | `NULL` | Path to existing all_limma_combined.txt for replot mode |
| `interactive` | `TRUE` | Generate interactive HTML plots |
| `width` | `14` | Heatmap width in inches |
| `height` | `10` | Heatmap height in inches |
| `volcano_width` | `10` | Volcano plot width in inches |
| `volcano_height` | `8` | Volcano plot height in inches |

### Replot Mode

Regenerate heatmap from existing limma results without re-running GSVA + limma:

```r
run_gsva_limma_de(
  limma_results_file = "output/deg/gsva/all_limma_combined.txt",
  output_dir = "output",
  column_order = c("TGFb_vs_Control", "Drug_vs_TGFb"),  # Optional: custom order
  column_labels = c("TGFb_vs_Control" = "Stim", "Drug_vs_TGFb" = "Drug"),  # Optional: custom labels
  cluster_rows = FALSE,  # Override default clustering
  cluster_columns = FALSE
)
```

### Output Files

```
output_dir/
├── deg/gsva/
│   ├── gsva_scores.txt              # GSVA matrix (signatures × samples)
│   ├── {comparison}_limma.txt       # Per-comparison limma results
│   └── all_limma_combined.txt       # Combined results (all comparisons)
└── figures/gsva/
    ├── volcano_{comparison}.png     # Per-comparison volcano plots
    ├── signatures_heatmap.png       # Summary heatmap (logFC + asterisks)
    └── interactive/
        ├── volcano_{comparison}_interactive.html
        └── signatures_heatmap_interactive.html
```

### Comparison Definition Format

Each comparison in the list must have three elements:

```r
list(
  name = "TGFb_vs_Control",     # Output file name prefix
  numerator = "TGFb",           # Treatment group (numerator in contrast)
  denominator = "Control"       # Reference group (denominator in contrast)
)
```

The `numerator` and `denominator` values must exactly match values in the `group_col` column of your metadata.

### Statistical Notes

- **GSVA**: Uses Poisson kernel (`kcdf = "Poisson"`) for raw count data
- **Limma**: Cell-means model (`~0 + Treatment`) with explicit contrasts
- **Significance**: Uses **nominal p-value** (not adjusted) for heatmap asterisks and volcano y-axis. This is appropriate for small signature panels (<100 signatures) where FDR correction is overly conservative.
- **Asterisk thresholds**: `*` p<0.05, `**` p<0.01, `***` p<0.001, `****` p<0.0001

### Return Value

Returns a list (invisibly) with:
- `gsva_matrix`: GSVA score matrix (signatures × samples)
- `limma_results`: List of limma results per comparison
- `combined_results`: Combined data frame of all results
- `output_files`: Paths to all output files

---

## Reference Documentation

See `references/plotting_guide.md` for detailed documentation including:
- Full parameter tables for each function
- Statistical method details
- Data compatibility
- Troubleshooting guide
