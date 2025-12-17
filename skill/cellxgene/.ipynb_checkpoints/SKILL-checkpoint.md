---
name: gene-expression-specificity-cellxgene
description: Analyze gene expression specificity across tissues and cell types using the CELLxGENE Census CellGuide API. Use this skill for extracting computational and canonical marker genes for any cell type, with optional filtering by organism and tissue, and visualization through dotplots, heatmaps, and barplots. This skill applies when users request cell type markers, specificity analysis, or marker visualizations from single-cell data.
---

# Gene Expression Specificity - CellGuide

Extract and visualize cell type marker genes from the CellGuide database with optional filtering by organism and tissue.

## Overview

This skill provides tools to:
1. Extract computational marker genes (derived from CELLxGENE Census data)
2. Extract canonical marker genes (literature-curated)
3. Filter markers by organism (e.g., human, mouse) and tissue
4. Visualize markers with dotplots, heatmaps, and barplots
5. Generate summary statistics and export results

## When to Use This Skill

Use this skill when users request:
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

**Basic Extraction**:
```bash
python scripts/extract_markers.py \
  --cell-type "B cell" \
  --output bcell_markers \
  --summary
```

**With Organism Filter**:
```bash
python scripts/extract_markers.py \
  --cell-type "T cell" \
  --organism "Homo sapiens" \
  --output tcell_human \
  --summary
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

**Top N Markers**:
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

Use `scripts/visualize_markers.py` to create publication-quality visualizations.

**Dotplot Visualization**:

The dotplot is the primary visualization for marker genes:
- **Color intensity**: Indicates effect size (marker_score or specificity)
- **Dot size**: Indicates percentage of cells expressing (pc)

```bash
python scripts/visualize_markers.py \
  --input bcell_markers_computational.csv \
  --output bcell_dotplot.png \
  --top-n 30 \
  --title "B Cell Marker Genes"
```

**Dotplot with Grouping**:
```bash
python scripts/visualize_markers.py \
  --input markers_computational.csv \
  --output dotplot_grouped.png \
  --top-n 40 \
  --group-by organism_ontology_term_label
```

**Filtered Dotplot**:
```bash
python scripts/visualize_markers.py \
  --input markers_computational.csv \
  --output human_dotplot.png \
  --filter-organism "Homo sapiens" \
  --top-n 25
```

**Heatmap Visualization**:
```bash
python scripts/visualize_markers.py \
  --input markers_computational.csv \
  --plot-type heatmap \
  --top-n 30 \
  --output markers_heatmap.png
```

**Barplot Visualization**:
```bash
python scripts/visualize_markers.py \
  --input markers_computational.csv \
  --plot-type barplot \
  --top-n 20 \
  --output markers_barplot.png
```

**Custom Color and Size Metrics**:
```bash
python scripts/visualize_markers.py \
  --input markers_computational.csv \
  --output dotplot.png \
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
# Extract and visualize in one workflow
python scripts/extract_markers.py \
  --cell-type "B cell" \
  --top-n 30 \
  --output bcell \
  --summary

python scripts/visualize_markers.py \
  --input bcell_computational.csv \
  --output bcell_dotplot.png \
  --top-n 20
```

### Pattern 2: Species-Specific Analysis

User asks: "Show me human-specific T cell markers"

```bash
python scripts/extract_markers.py \
  --cell-type "T cell" \
  --organism "Homo sapiens" \
  --top-n 50 \
  --output tcell_human \
  --summary

python scripts/visualize_markers.py \
  --input tcell_human_computational.csv \
  --output tcell_human_dotplot.png \
  --top-n 30 \
  --title "Human T Cell Markers"
```

### Pattern 3: Tissue-Specific Markers

User asks: "Find macrophage markers in human lung tissue"

```bash
python scripts/extract_markers.py \
  --cell-type "macrophage" \
  --organism "Homo sapiens" \
  --tissue "lung" \
  --output macro_lung

python scripts/visualize_markers.py \
  --input macro_lung_computational.csv \
  --output macro_lung_dotplot.png \
  --top-n 25 \
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
- `--top-n`: Limit to top N markers by marker_score
- `--canonical`: Also extract canonical markers
- `--summary`: Print summary statistics
- `--output`: Output file base name (required)
- `--format`: Output format (csv, tsv, excel)

**Outputs**:
- `{output}_computational.csv`: Computational markers
- `{output}_canonical.csv`: Canonical markers (if --canonical flag used)

### visualize_markers.py

**Purpose**: Create visualizations from marker data

**Key Arguments**:
- `--input`: Input CSV file with marker data (required)
- `--output`: Output file path for plot (required)
- `--plot-type`: dotplot, heatmap, or barplot (default: dotplot)
- `--top-n`: Number of top genes to show (default: 20)
- `--group-by`: Group by organism or tissue
- `--filter-organism`: Filter to specific organism(s)
- `--filter-tissue`: Filter to specific tissue(s)
- `--color-by`: Column for color intensity (default: marker_score)
- `--size-by`: Column for dot size (default: pc)
- `--title`: Custom plot title
- `--figsize`: Figure dimensions (width height)

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

1. **Start broad, then filter**: Extract all markers first, then apply filters as needed
2. **Use --summary flag**: Get overview statistics before diving into details
3. **Specify organism**: Many queries benefit from organism filtering (human or mouse)
4. **Top N selection**: Use --top-n to focus on highest-scoring markers (typically 20-50)
5. **Combine with canonical**: Cross-reference computational with canonical markers for validation
6. **Group visualizations**: Use --group-by to compare across organisms or tissues
7. **Export format**: Use Excel format for sharing with non-computational collaborators
8. **Cell Ontology IDs**: Use when exact cell type matching is critical

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
# User request: "Find and visualize the top 25 markers for human B cells"

# Step 1: Extract markers
python scripts/extract_markers.py \
  --cell-type "B cell" \
  --organism "Homo sapiens" \
  --top-n 25 \
  --output bcell_human \
  --summary

# Output:
# Looking up Cell Ontology ID for 'B cell'...
# Found: CL:0000236 (CL_0000236)
# Fetching computational marker genes...
# Retrieved 3067 computational markers
# Applying filters...
# After filtering: 2016 markers
# Selected top 25 markers by marker_score
# Exported to bcell_human_computational.csv
#
# === Summary Statistics ===
# Total markers: 25
# Unique genes: 25
# Mean marker score: 2.145
# Mean specificity: 0.967
# Mean expression: 3.421

# Step 2: Create dotplot visualization
python scripts/visualize_markers.py \
  --input bcell_human_computational.csv \
  --output bcell_human_dotplot.png \
  --top-n 25 \
  --title "Top 25 Human B Cell Markers"

# Output:
# Loading marker data from bcell_human_computational.csv...
# Loaded 25 markers for 25 unique genes
# Creating dotplot...
# Dotplot saved to bcell_human_dotplot.png
# ✓ Done!

# Result: User receives CSV file and publication-quality dotplot visualization
```

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
