# RNA-seq Visualization Functions Reference

## Overview

Four standalone plotting functions for **ad-hoc visualization of DESeq2 analysis pipeline outputs**.

These functions are designed for on-demand visualization requests after the main DESeq2 differential expression analysis is complete. They are **not** part of the 10-step DESeq2 workflow - instead, they serve as supplementary tools for:

- Creating custom heatmaps from GSEA/DEG results
- Visualizing gene expression patterns for specific signatures
- Generating publication-ready GSVA/ssGSEA boxplots with statistics
- Exploring results interactively based on user questions

While optimized for DESeq2 workflow outputs (`*_analysis_data.rds`, `deg/*.txt`), these functions also accept compatible data from external sources.

## Function 1: plot_signature_heatmap()

### Purpose
Creates clustered heatmap with GSEA/GSVA signatures as rows and comparisons as columns.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| gsea_data | char/df | required | GSEA results file or data frame |
| signatures | char[] | NULL | Signature names (NULL = all) |
| output_file | char | "figures/signature_NES_heatmap.png" | Output path |
| source_filter | char | "CUSTOM" | Filter by source column |
| title | char | NULL | Plot title (auto if NULL) |
| width | num | 12 | Width in inches |
| height | num | 8 | Height in inches |
| cluster_rows | bool | TRUE | Cluster signatures |
| cluster_cols | bool | TRUE | Cluster comparisons |

### Input Format
Tab-delimited with columns: ID, comparison, NES, p.adjust, source

### Output
- PNG heatmap with significance annotations (*, **, ***, ****) in cells
- Interactive HTML file saved to `figures/interactive/expr_comp_heatmap_boxplot/`
- Blue-white-red color scale centered at 0

### Interactive Features
- Hover tooltips: Signature name, Comparison, NES, p-value, FDR
- Cell annotations: `*` significance stars displayed in heatmap cells
- Set `interactive = FALSE` to disable HTML generation

### Example
```r
plot_signature_heatmap(
  gsea_data = "deg/RNAseq_gsea_all.txt",
  signatures = c("iaf_sig_Doug75", "tgfb_sig_Jinjin"),
  output_file = "figures/custom/NES_heatmap.png",
  source_filter = "CUSTOM"
)
```

---

## Function 2: plot_deg_heatmap()

### Purpose
Creates clustered heatmap with genes as rows and comparisons as columns, showing log2FC.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| deg_data | char/df | required | DEG stats file or data frame |
| signature_files | char[] | NULL | Paths to gene list files |
| genes | char[] | NULL | Gene symbols (if no files) |
| output_file | char | "figures/deg_gene_heatmap.png" | Output path |
| export_individual | bool | TRUE | Separate plot per signature |
| title | char | NULL | Plot title (auto if NULL) |
| width | num | 12 | Width in inches |
| height | num | 10 | Height in inches |

### Input Format
- DEG file: Tab-delimited with columns: symbol, comparison, log2FoldChange, padj
- Signature files: Plain text, one gene per line

### Output
- PNG heatmap(s) with significance annotations (*, **, ***, ****) in cells
- Interactive HTML file saved to `figures/interactive/expr_comp_heatmap_boxplot/`
- Individual plots per signature if export_individual=TRUE

### Interactive Features
- Hover tooltips: Gene symbol, Comparison, log2FC, p-value, FDR
- Cell annotations: `*` significance stars displayed in heatmap cells
- Set `interactive = FALSE` to disable HTML generation

### Example
```r
plot_deg_heatmap(
  deg_data = "deg/RNAseq_summstats_all.txt",
  signature_files = c("signatures/iaf_sig.txt", "signatures/tgfb_sig.txt"),
  output_file = "figures/custom/deg_heatmap.png"
)
```

---

## Function 3: plot_expression_heatmap()

### Purpose
Creates clustered heatmap with genes as rows and samples as columns, grouped by treatment.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| expr_data | char/matrix/RDS | required | Expression matrix |
| metadata | char/df | required | Sample metadata |
| signature_files | char[] | NULL | Paths to gene list files |
| genes | char[] | NULL | Gene symbols (if no files) |
| group_col | char | "Treatment" | Grouping column |
| sample_col | char | "Sample.ID" | Sample ID column |
| output_file | char | "figures/expression_heatmap.png" | Output path |
| export_individual | bool | TRUE | Separate plot per signature |
| zscore | bool | TRUE | Z-score normalize rows |
| width | num | 14 | Width in inches |
| height | num | 10 | Height in inches |

### Input Format
- Expression: genes as rows, samples as columns (or RDS with vsd/normalized_counts)
- Metadata: Tab-delimited with sample IDs and group column
- Signature files: Plain text, one gene per line

### Output
- PNG heatmap(s) with treatment group color bar
- Interactive HTML file saved to `figures/interactive/expr_comp_heatmap_boxplot/`
- Z-scored expression values

### Interactive Features
- Hover tooltips: Gene symbol, Sample, Treatment group, Z-score/expression value
- No cell annotations (colors represent expression values)
- Set `interactive = FALSE` to disable HTML generation

### Example
```r
plot_expression_heatmap(
  expr_data = "RNAseq_analysis_data.rds",
  metadata = "processed_metadata.txt",
  signature_files = c("signatures/iaf_sig.txt"),
  group_col = "Treatment",
  sample_col = "Sample ID",
  output_file = "figures/custom/expression_heatmap.png"
)
```

---

## Function 4: plot_gsva_boxplot()

### Purpose
Creates boxplots of GSVA/ssGSEA scores by treatment with statistical comparisons.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| gsva_data | char/matrix/RDS | required | GSVA score matrix |
| metadata | char/df | required | Sample metadata |
| signatures | char[] | required | Signature names to plot |
| group_col | char | "Treatment" | Grouping column |
| sample_col | char | "Sample.ID" | Sample ID column |
| output_file | char | "figures/gsva_boxplot.png" | Output path |
| stat_compare | char | "adaptive" | Statistical method |
| max_pairs | num | 10 | Max pairs to display on plot |
| export_stats | bool | TRUE | Export statistics TSV |
| width | num | 12 | Width in inches |
| height | num | 8 | Height in inches |

### Statistical Methods
- **adaptive**: Levene's test -> Tukey HSD (equal variance) or Games-Howell (unequal)
- **kruskal**: Kruskal-Wallis global test
- **wilcox_ref**: Wilcoxon vs reference group
- **wilcox_all**: All pairwise Wilcoxon
- **none**: No statistics

### Input Format
- GSVA: signatures as rows, samples as columns (or RDS with gsva_scores)
- Metadata: Tab-delimited with sample IDs and group column

### Output
- PNG boxplot(s) saved to `figures/expr_comp_heatmap_boxplot/`
- Interactive HTML file saved to `figures/interactive/expr_comp_heatmap_boxplot/`
- TSV file with all pairwise statistics saved to `figures/expr_comp_heatmap_boxplot/*_pairwise_stats.tsv`
- Subtitle shows "Showing X of Y significant pairs" if truncated

### Interactive Features
- Hover tooltips: Treatment group, GSVA score, individual data points
- Set `interactive = FALSE` to disable HTML generation

### Example
```r
plot_gsva_boxplot(
  gsva_data = "RNAseq_analysis_data.rds",
  metadata = "processed_metadata.txt",
  signatures = c("iaf_sig_Doug75", "tgfb_sig_Jinjin"),
  group_col = "Treatment",
  sample_col = "Sample ID",
  output_file = "figures/custom/gsva_boxplot.png",
  stat_compare = "adaptive",
  max_pairs = 10,
  export_stats = TRUE
)
```

---

## Common Features

### Significance Annotations
- `*` p < 0.05
- `**` p < 0.01
- `***` p < 0.001
- `****` p < 0.0001

### Color Scheme
- Blue-White-Red (BWR) centered at 0
- Blue = negative (downregulated/negative NES)
- Red = positive (upregulated/positive NES)

### Dependencies
```r
# Core plotting
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(circlize)
library(rstatix)
library(car)

# Interactive output
library(heatmaply)
library(plotly)
library(htmlwidgets)
```

---

## Example Workflow

```r
# Set working directory to DESeq2 output
setwd("/path/to/deseq2_output")

# Load functions
source("path/to/plotting_functions.R")

# 1. Signature NES heatmap
plot_signature_heatmap(
  gsea_data = "deg/RNAseq_gsea_all.txt",
  signatures = c("iaf_sig_Doug75", "tgfb_sig_Jinjin"),
  output_file = "figures/custom/NES_heatmap.png",
  source_filter = "CUSTOM"
)

# 2. DEG gene heatmap
plot_deg_heatmap(
  deg_data = "deg/RNAseq_summstats_all.txt",
  signature_files = c("signatures/iaf_sig.txt", "signatures/tgfb_sig.txt"),
  output_file = "figures/custom/deg_heatmap.png"
)

# 3. Expression heatmap
plot_expression_heatmap(
  expr_data = "RNAseq_analysis_data.rds",
  metadata = "processed_metadata.txt",
  signature_files = c("signatures/iaf_sig.txt"),
  group_col = "Treatment",
  sample_col = "Sample ID",
  output_file = "figures/custom/expression_heatmap.png"
)

# 4. GSVA boxplot with statistics
plot_gsva_boxplot(
  gsva_data = "RNAseq_analysis_data.rds",
  metadata = "processed_metadata.txt",
  signatures = c("iaf_sig_Doug75", "tgfb_sig_Jinjin"),
  group_col = "Treatment",
  sample_col = "Sample ID",
  output_file = "figures/custom/gsva_boxplot.png",
  stat_compare = "adaptive",
  max_pairs = 10,
  export_stats = TRUE
)
```

---

## Data Compatibility

### From DESeq2 Workflow
These functions work seamlessly with outputs from the main DESeq2 analysis:

| Function | Compatible Files |
|----------|-----------------|
| plot_signature_heatmap | `deg/{prefix}_gsea_all.txt` |
| plot_deg_heatmap | `deg/{prefix}_summstats_all.txt` |
| plot_expression_heatmap | `{prefix}_analysis_data.rds` (vsd slot) |
| plot_gsva_boxplot | `{prefix}_analysis_data.rds` (gsva_scores slot) |

### From External Sources
Functions also accept:
- Tab-delimited text files with appropriate columns
- R data frames
- Expression matrices (genes x samples)
- Any GSEA results with ID, NES, p.adjust, comparison columns

---

## Troubleshooting

### Common Issues

**"Gene not found in expression matrix"**
- Check that gene symbols match between signature file and expression data
- Verify gene symbol column name in expression data

**"Sample ID not found in metadata"**
- Ensure sample_col parameter matches actual column name in metadata
- Check for whitespace in column names

**"No significant pairs to display"**
- All comparisons may be non-significant
- Try `stat_compare = "none"` to show plot without statistics
- Check `figures/expr_comp_heatmap_boxplot/*_pairwise_stats.tsv` for all p-values

**Heatmap too crowded**
- Reduce number of signatures or genes
- Increase width/height parameters
- Use export_individual=TRUE for separate plots
