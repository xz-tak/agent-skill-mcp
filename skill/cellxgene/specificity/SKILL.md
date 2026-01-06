---
name: cellxgene-specificity
description: Analyze gene expression specificity across tissues and cell types using the CELLxGENE Census CellGuide API. Use this for extracting computational and canonical marker genes for any cell type, with optional filtering by organism and tissue, and visualization through dotplots, heatmaps, and barplots. This subskill applies when users request cell type markers, specificity analysis, or marker visualizations from single-cell data.
---

# Gene Expression Specificity - CellGuide

Extract and visualize cell type marker genes from the CellGuide database with optional filtering by organism and tissue.

## Overview

This subskill provides tools to:
1. Extract computational marker genes (derived from CELLxGENE Census data)
2. Extract canonical marker genes (literature-curated)
3. Filter markers by organism (e.g., human, mouse), tissue, and marker_score quality
4. Visualize markers with dotplots, heatmaps, and barplots
5. Generate summary statistics and export results

**Default Behavior**: By default, markers are filtered to include only genes with marker_score >= 0.5 to ensure high-quality results. Users can customize this threshold or request top N markers instead.

## When to Use This Subskill

Use this subskill when users request:
- "Find markers for [cell type]"
- "What are the top markers for [cell type] in [organism]?"
- "Show me [cell type] markers specific to [tissue]"
- "Visualize marker genes for [cell type]"
- "Get canonical markers for [cell type]"
- "Compare [cell type] markers across tissues/organisms"

## Core Concepts

### Marker Types

**Computational Markers**: Identified through statistical analysis of CELLxGENE Census data
- Includes metrics: marker_score, specificity, mean expression (me), percent cells (pc)
- Available for specific organism-tissue combinations
- Quantitative and data-driven

**Canonical Markers**: Curated from scientific literature
- Includes tissue context and publication references
- Qualitative and expert-validated
- Generally smaller sets of well-established markers

### Key Metrics

- **marker_score**: Overall score indicating marker quality (higher = better marker)
- **specificity**: How specific the gene is to the cell type (0-1 scale)
- **me**: Mean expression level across cells
- **pc**: Percentage of cells expressing the gene (0-1 scale)

## Workflow

### Step 1: Extract Markers

Use `scripts/extract_markers.py` to query CellGuide and extract marker genes.

**Basic Extraction** (default: marker_score >= 0.5):
```bash
python scripts/extract_markers.py \
  --cell-type "B cell" \
  --output bcell_markers \
  --summary
```

**With Organism Filter** (default: marker_score >= 0.5):
```bash
python scripts/extract_markers.py \
  --cell-type "T cell" \
  --organism "Homo sapiens" \
  --output tcell_human \
  --summary
```

**With Custom Marker Score Cutoff**:
```bash
python scripts/extract_markers.py \
  --cell-type "macrophage" \
  --organism "Homo sapiens" \
  --marker-score-cutoff 1.0 \
  --output macro_human_high_quality
```

**With Combined Filters**:
```bash
python scripts/extract_markers.py \
  --cell-type "macrophage" \
  --organism "Homo sapiens" \
  --tissue "lung" \
  --output macro_human_lung
```

**Using Cell Ontology ID**:
```bash
python scripts/extract_markers.py \
  --ontology-id CL_0000236 \
  --output markers
```

**Top N Markers** (disables marker_score cutoff):
```bash
python scripts/extract_markers.py \
  --cell-type "dendritic cell" \
  --organism "Homo sapiens" \
  --top-n 50 \
  --output dc_top50
```

**Include Canonical Markers**:
```bash
python scripts/extract_markers.py \
  --cell-type "B cell" \
  --canonical \
  --output bcell_all_markers
```

### Step 2: Visualize Markers

Use `scripts/visualize_markers.py` to create interactive visualizations.

**Default Output**: Interactive HTML files with Plotly visualizations
- **Hover details**: Gene name, marker scores, specificity, and expression metrics
- **Interactive features**: Zoom, pan, and export capabilities
- **Color scheme**: Uses 'Reds' colormap for consistency
- **Saved to**: Working directory (or user-specified path)
- **Requires**: `plotly` package (install with: `pip install plotly`)

**Optional Static PNG**: Set `generate_static=True` in function call to also create high-resolution PNG (300 DPI)

**Dotplot Visualization**:

The dotplot is the primary visualization for marker genes:
- **Color intensity**: Indicates effect size (marker_score or specificity)
- **Dot size**: Indicates percentage of cells expressing (pc)
- **Default**: Shows genes with marker_score >= 0.5

```bash
# Default: generates interactive HTML with markers having marker_score >= 0.5
# Output: bcell_dotplot.html
python scripts/visualize_markers.py \
  --input bcell_markers_computational.csv \
  --output bcell_dotplot.html \
  --title "B Cell Marker Genes"
```

**Dotplot with Top N Genes**:
```bash
# Output: dotplot_top.html
python scripts/visualize_markers.py \
  --input markers_computational.csv \
  --output dotplot_top.html \
  --top-n 30 \
  --title "Top 30 Marker Genes"
```

**Dotplot with Grouping**:
```bash
# Output: dotplot_grouped.html
python scripts/visualize_markers.py \
  --input markers_computational.csv \
  --output dotplot_grouped.html \
  --top-n 40 \
  --group-by organism_ontology_term_label
```

**Filtered Dotplot with Custom Cutoff**:
```bash
# Output: human_dotplot.html
python scripts/visualize_markers.py \
  --input markers_computational.csv \
  --output human_dotplot.html \
  --filter-organism "Homo sapiens" \
  --marker-score-cutoff 1.0
```

**Heatmap Visualization**:
```bash
# Output: markers_heatmap.html
python scripts/visualize_markers.py \
  --input markers_computational.csv \
  --plot-type heatmap \
  --top-n 30 \
  --output markers_heatmap.html
```

**Barplot Visualization**:
```bash
# Output: markers_barplot.html
python scripts/visualize_markers.py \
  --input markers_computational.csv \
  --plot-type barplot \
  --marker-score-cutoff 0.8 \
  --output markers_barplot.html
```

**Custom Color and Size Metrics**:
```bash
# Output: dotplot.html
python scripts/visualize_markers.py \
  --input markers_computational.csv \
  --output dotplot.html \
  --color-by specificity \
  --size-by me \
  --top-n 30
```

### Step 3: Interpret Results

After extraction and visualization:

1. **Review marker_score**: Higher scores indicate more robust markers
2. **Check specificity**: Values close to 1.0 indicate cell type-specific expression
3. **Evaluate pc**: Higher percentage means more consistent expression across cells
4. **Consider organism/tissue context**: Markers may vary across contexts

## Common Usage Patterns

### Pattern 1: Quick Marker Lookup

User asks: "What are the top markers for B cells?"

```bash
# Extract high-quality markers (default: marker_score >= 0.5)
python scripts/extract_markers.py \
  --cell-type "B cell" \
  --output bcell \
  --summary

# Visualize (will apply default marker_score >= 0.5 filter)
python scripts/visualize_markers.py \
  --input bcell_computational.csv \
  --output bcell_dotplot.png
```

Alternative with top N:
```bash
# Extract top 30 markers
python scripts/extract_markers.py \
  --cell-type "B cell" \
  --top-n 30 \
  --output bcell_top30 \
  --summary

# Visualize top 20
python scripts/visualize_markers.py \
  --input bcell_top30_computational.csv \
  --output bcell_dotplot.png \
  --top-n 20
```

### Pattern 2: Species-Specific Analysis

User asks: "Show me human-specific T cell markers"

```bash
# Extract human T cell markers (default: marker_score >= 0.5)
python scripts/extract_markers.py \
  --cell-type "T cell" \
  --organism "Homo sapiens" \
  --output tcell_human \
  --summary

# Visualize with default cutoff
python scripts/visualize_markers.py \
  --input tcell_human_computational.csv \
  --output tcell_human_dotplot.png \
  --title "Human T Cell Markers"
```

### Pattern 3: Tissue-Specific Markers

User asks: "Find macrophage markers in human lung tissue"

```bash
# Extract with default marker_score >= 0.5
python scripts/extract_markers.py \
  --cell-type "macrophage" \
  --organism "Homo sapiens" \
  --tissue "lung" \
  --output macro_lung

# Visualize with default cutoff
python scripts/visualize_markers.py \
  --input macro_lung_computational.csv \
  --output macro_lung_dotplot.png \
  --title "Macrophage Markers in Human Lung"
```

### Pattern 4: Comparative Analysis

User asks: "Compare B cell markers across human and mouse"

```bash
# Extract markers (no organism filter)
python scripts/extract_markers.py \
  --cell-type "B cell" \
  --output bcell_multi \
  --summary

# Visualize with grouping
python scripts/visualize_markers.py \
  --input bcell_multi_computational.csv \
  --output bcell_comparison_dotplot.png \
  --top-n 40 \
  --group-by organism_ontology_term_label \
  --title "B Cell Markers: Human vs Mouse"
```

### Pattern 5: Literature-Validated Markers

User asks: "Get canonical markers for dendritic cells"

```bash
python scripts/extract_markers.py \
  --cell-type "dendritic cell" \
  --canonical \
  --output dc_canonical

# Canonical markers don't have computational metrics
# Review the CSV directly or create a simple visualization
```

### Pattern 6: Multiple Cell Types

User asks: "Get markers for several immune cell types"

```bash
# Extract markers for each cell type
for cell_type in "B cell" "T cell" "monocyte" "natural killer cell"; do
  python scripts/extract_markers.py \
    --cell-type "$cell_type" \
    --organism "Homo sapiens" \
    --top-n 30 \
    --output "${cell_type// /_}_human"
done

# Visualize each
for file in *_human_computational.csv; do
  base="${file%_computational.csv}"
  python scripts/visualize_markers.py \
    --input "$file" \
    --output "${base}_dotplot.png" \
    --top-n 20
done
```

## Script Reference

### extract_markers.py

**Purpose**: Query CellGuide API and extract marker genes

**Key Arguments**:
- `--cell-type`: Cell type name (e.g., "B cell", "neuron")
- `--ontology-id`: Cell Ontology ID (e.g., CL_0000236)
- `--organism`: Filter by organism(s)
- `--tissue`: Filter by tissue(s)
- `--marker-score-cutoff`: Filter genes with marker_score above this cutoff (default: 0.5 if --top-n not specified)
- `--top-n`: Limit to top N markers by marker_score (disables marker_score cutoff filtering)
- `--canonical`: Also extract canonical markers
- `--summary`: Print summary statistics
- `--output`: Output file base name (required)
- `--format`: Output format (csv, tsv, excel)

**Default Behavior**:
- Filters markers to include only genes with marker_score >= 0.5
- This ensures high-quality, well-validated markers
- Use `--marker-score-cutoff` to customize the threshold
- Use `--top-n` to get top N markers instead (disables cutoff filter)

**Outputs**:
- `{output}_computational.csv`: Computational markers
- `{output}_canonical.csv`: Canonical markers (if --canonical flag used)

### visualize_markers.py

**Purpose**: Create visualizations from marker data

**Key Arguments**:
- `--input`: Input CSV file with marker data (required)
- `--output`: Output file path for plot (required)
- `--plot-type`: dotplot, heatmap, or barplot (default: dotplot)
- `--marker-score-cutoff`: Filter genes with marker_score above this cutoff (default: 0.5 if --top-n not specified)
- `--top-n`: Number of top genes to show (disables marker_score cutoff filtering)
- `--group-by`: Group by organism or tissue
- `--filter-organism`: Filter to specific organism(s)
- `--filter-tissue`: Filter to specific tissue(s)
- `--color-by`: Column for color intensity (default: marker_score)
- `--size-by`: Column for dot size (default: pc)
- `--title`: Custom plot title
- `--figsize`: Figure dimensions (width height)

**Default Behavior**:
- Filters to genes with marker_score >= 0.5 before plotting
- Ensures visualizations show only high-quality markers
- Use `--top-n` to show specific number of top genes (disables cutoff)
- Use `--marker-score-cutoff` to customize the threshold

**Outputs**:
- High-resolution PNG image (300 DPI)

## Data Structure

### Computational Markers CSV Columns

- `me`: Mean expression across expressing cells
- `pc`: Percentage of cells expressing (0-1)
- `marker_score`: Overall marker quality score
- `specificity`: Cell type specificity (0-1)
- `gene_ontology_term_id`: Ensembl gene ID
- `symbol`: Gene symbol (e.g., CD79A, MS4A1)
- `name`: Full gene name
- `organism_ontology_term_label`: Organism (e.g., "Homo sapiens")
- `tissue_ontology_term_label`: Tissue context (may be NaN)

### Canonical Markers CSV Columns

- `tissue`: Tissue context
- `symbol`: Gene symbol
- `name`: Full gene name
- `publication`: Publication reference
- `publication_titles`: Publication titles

## Tips and Best Practices

1. **Default filtering is enabled**: By default, results are filtered to marker_score >= 0.5 for quality
2. **Adjust cutoff for stricter results**: Use `--marker-score-cutoff 1.0` or higher for highly specific markers
3. **Use --summary flag**: Get overview statistics to understand the data before visualization
4. **Specify organism**: Many queries benefit from organism filtering (human or mouse)
5. **Top N vs cutoff**:
   - Use `--marker-score-cutoff` (default) for quality-based filtering
   - Use `--top-n` when you need a specific number of markers
   - These options are mutually exclusive
6. **Combine with canonical**: Cross-reference computational with canonical markers for validation
7. **Group visualizations**: Use --group-by to compare across organisms or tissues
8. **Export format**: Use Excel format for sharing with non-computational collaborators
9. **Cell Ontology IDs**: Use when exact cell type matching is critical

## Troubleshooting

**Issue**: "Could not find Cell Ontology ID for [cell type]"
**Solution**: Try alternative names (e.g., "B lymphocyte" instead of "B cell") or find the Cell Ontology ID manually and use --ontology-id

**Issue**: "No data remaining after filtering"
**Solution**: Filters may be too restrictive. Remove tissue filter or check available tissues in summary output

**Issue**: Empty visualization
**Solution**: Ensure CSV has data and check column names match expectations (marker_score, pc, symbol)

**Issue**: Dotplot too crowded
**Solution**: Reduce --top-n parameter or increase --figsize

## Example Session

```bash
# User request: "Find and visualize high-quality markers for human B cells"

# Step 1: Extract markers (default: marker_score >= 0.5)
python scripts/extract_markers.py \
  --cell-type "B cell" \
  --organism "Homo sapiens" \
  --output bcell_human \
  --summary

# Output:
# Looking up Cell Ontology ID for 'B cell'...
# Found: CL:0000236 (CL_0000236)
# Fetching computational marker genes...
# Retrieved 3067 computational markers
# Applying filters...
# After filtering: 2016 markers
# Applied marker_score cutoff >= 0.5: 145 markers (removed 1871)
# Exported to bcell_human_computational.csv
#
# === Summary Statistics ===
# Total markers: 145
# Unique genes: 102
# Mean marker score: 1.523
# Mean specificity: 0.892
# Mean expression: 2.841

# Step 2: Create dotplot visualization (default: marker_score >= 0.5)
python scripts/visualize_markers.py \
  --input bcell_human_computational.csv \
  --output bcell_human_dotplot.png \
  --title "Human B Cell Markers"

# Output:
# Loading marker data from bcell_human_computational.csv...
# Loaded 145 markers for 102 unique genes
# Applied marker_score cutoff >= 0.5: 145 markers (removed 0)
# Creating dotplot...
# Dotplot saved to bcell_human_dotplot.png
# ✓ Done!

# Result: User receives CSV file with high-quality markers and publication-quality dotplot

# Alternative: Request specific top N markers
python scripts/extract_markers.py \
  --cell-type "B cell" \
  --organism "Homo sapiens" \
  --top-n 25 \
  --output bcell_human_top25 \
  --summary

# This will extract exactly 25 markers, bypassing the marker_score cutoff
```

## Census-based Specificity Analysis

For analyzing cell type specificity of **specific genes** rather than finding markers for a cell type, use `compute_specificity.py`. This tool:

1. Queries raw expression data from CELLxGENE Census
2. Filters by tissue/disease/species but includes **all cell types**
3. **Optionally target specific cell type(s)** via name, substring, regex, or Cell Ontology ID
4. Runs scanpy's `rank_genes_groups` with wilcoxon test
5. Computes custom specificity scores for target cell types

### When to Use compute_specificity.py

Use this tool when you want to:
- Find which cell types express specific genes (reverse lookup)
- Analyze cell type specificity for a gene signature
- **Calculate specificity for specific cell type(s)** (e.g., fibroblasts, plasma cells)
- Get detailed statistical metrics (fold change, p-values, % expressing)
- Filter by tissue/disease context while comparing all cell types

### compute_specificity.py Usage

**Target-specific Analysis (Recommended):**
```bash
# Analyze specificity for fibroblasts in intestinal tissue
python scripts/compute_specificity.py \
  --genes "IL11,GREM1,TYK2,JAK1" \
  --tissue intestine \
  --target-cell-type "fibroblast" \
  --output fibroblast_specificity
```

**Substring Matching:**
```bash
# Match all cell types containing "plasma"
python scripts/compute_specificity.py \
  --genes "IL11,GREM1" \
  --tissue intestine \
  --target-cell-type "plasma" \
  --match-mode substring \
  --output plasma_specificity
```

**Regex Pattern Matching:**
```bash
# Match fibroblasts, myofibroblasts, and smooth muscle cells
python scripts/compute_specificity.py \
  --genes "IL11,GREM1" \
  --tissue intestine \
  --target-cell-type "fibroblast|myofibroblast|smooth muscle" \
  --match-mode regex \
  --output stromal_specificity
```

**Cell Ontology ID Matching:**
```bash
# Use Cell Ontology ID (CL:0000057 = fibroblast)
python scripts/compute_specificity.py \
  --genes "APOE,APOC1" \
  --tissue liver \
  --target-cell-type "CL:0000057" \
  --match-mode id \
  --output hepatocyte_specificity
```

**All Cell Types (no target filter):**
```bash
python scripts/compute_specificity.py \
  --genes "IL11,GREM1" \
  --tissue intestine \
  --output all_cells_specificity
```

**With Disease Filter:**
```bash
python scripts/compute_specificity.py \
  --genes "IL11,GREM1,TYK2,JAK1" \
  --tissue intestine \
  --disease "Crohn disease,ulcerative colitis,normal" \
  --target-cell-type "fibroblast" \
  --output ibd_fibroblast_specificity
```

**Using Gene File:**
```bash
# genes.txt contains one gene symbol per line
python scripts/compute_specificity.py \
  --genes genes.txt \
  --tissue liver \
  --target-cell-type "hepatocyte" \
  --output liver_gene_specificity
```

### compute_specificity.py Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--genes` | Gene list (comma-separated) or file path | Required |
| `--output` | Output prefix for result files | Required |
| `--tissue` | Tissue filter(s), comma-separated | None |
| `--disease` | Disease filter(s), comma-separated | None |
| `--species` | Species (human/mouse) | human |
| `--development-stage` | Development stage filter | adult |
| `--target-cell-type` | Target cell type pattern for specificity calculation | None (all) |
| `--match-mode` | Cell type matching: exact/substring/regex/id | substring |
| `--group-by` | Column for cell grouping | cell_type |
| `--pval-cutoff` | P-value cutoff for markers | 0.05 |
| `--log2fc-min` | Minimum log2 fold change | 1.0 |
| `--output-dir` | Output directory | . |
| `--no-interactive` | Skip cell type matching confirmation | False |

### Cell Type Matching Modes

| Mode | Description | Example |
|------|-------------|---------|
| `exact` | Case-insensitive exact match | "B cell" matches only "B cell" |
| `substring` | Match if pattern is substring | "plasma" matches "plasma cell", "plasmablast" |
| `regex` | Regular expression pattern | "fibroblast\|myo" matches fibroblast, myofibroblast |
| `id` | Match by Cell Ontology ID | "CL:0000057" matches fibroblast |

**Interactive Confirmation:** By default, matched cell types are displayed and you can:
- `yes` - proceed with all matches
- `no` - cancel
- `select` - choose specific matches by number

Use `--no-interactive` to skip confirmation (auto-proceed with all matches).

### Specificity Score Methodology

The scoring formula combines multiple metrics:

```
lfc_logp = logfoldchanges × -log(pvals_adj + 1e-323) ×
           (pct_nz_group - pct_nz_reference) × (1 - pct_nz_reference)
```

This composite score rewards genes that:
- Have high fold change (strong differential expression)
- Have low p-values (statistically significant)
- Are expressed in a higher fraction of target cells vs reference
- Are NOT widely expressed in reference cells (specificity)

Then, specificity ratio is calculated:
```
mean_specificity_ratio = cell_type_lfc_logp / total_lfc_logp_across_all_cell_types
```

Finally, normalized as percentile rank:
```
score = percentile_rank(mean_specificity_ratio)
```

### Output Columns

| Column | Description |
|--------|-------------|
| `names` | Gene symbol |
| `group` | Cell type (from groupby column) |
| `logfoldchanges` | Log2 fold change vs other cell types |
| `pvals` | Raw p-values |
| `pvals_adj` | Adjusted p-values (Benjamini-Hochberg) |
| `pct_nz_group` | % cells expressing in this cell type |
| `pct_nz_reference` | % cells expressing in other cell types |
| `lfc_logp` | Raw composite score |
| `total_lfc_logp` | Sum of lfc_logp across all cell types |
| `mean_specificity_ratio` | Cell type specificity ratio |
| `score` | Normalized percentile rank (0-1) |

### Output Files

- `{output}_specificity.csv` - Full results with all scores
- `{output}_summary.json` - Query parameters and summary statistics

### Example: Finding Cell Types for IBD Targets

```bash
# Query: Which cell types express IL11 and GREM1 in IBD intestinal tissue?
python scripts/compute_specificity.py \
  --genes "IL11,GREM1" \
  --tissue intestine \
  --disease "Crohn disease,ulcerative colitis" \
  --output ibd_targets

# Results show specificity scores for each gene × cell type combination
# High scores indicate cell types where the gene is specifically expressed
```

### Comparison: extract_markers.py vs compute_specificity.py

| Feature | extract_markers.py | compute_specificity.py |
|---------|-------------------|----------------------|
| Data source | CellGuide API (pre-computed) | Census API (raw data) |
| Input | Cell type name | Gene list |
| Output | Markers for that cell type | Cell types for those genes |
| Use case | "Find markers for B cells" | "Which cells express IL11?" |
| Speed | Fast (pre-computed) | Slower (computes on-the-fly) |
| Customization | Limited filters | Full Census filters |

## Advanced Usage

### Custom Analysis Pipeline

For complex analyses, combine scripts with data processing:

```bash
# Extract markers for multiple tissues
for tissue in "spleen" "bone marrow" "lymph node"; do
  python scripts/extract_markers.py \
    --cell-type "B cell" \
    --organism "Homo sapiens" \
    --tissue "$tissue" \
    --output "bcell_${tissue// /_}"
done

# Combine and visualize
# (User can merge CSVs and create comparative visualizations)
```

### Integration with Analysis Pipelines

Marker data can be integrated with downstream analyses:
- Gene set enrichment analysis (GSEA)
- Pathway analysis
- Regulatory network inference
- Cross-species comparisons
- Multi-omics integration

## Environment Requirements

The scripts require Python 3.7+ with the following packages:
- pandas
- requests
- matplotlib
- seaborn
- numpy

Install via:
```bash
conda install pandas requests matplotlib seaborn numpy
# or
pip install pandas requests matplotlib seaborn numpy
```

## Data Source

All marker data comes from the CellGuide database, which aggregates and analyzes single-cell RNA-seq data from the CELLxGENE Census. The database is continuously updated with new data and improved algorithms.

CellGuide API: https://cellguide.cellxgene.cziscience.com
