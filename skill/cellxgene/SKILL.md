---
name: cellxgene
description: Comprehensive toolkit for working with CELLxGENE Census single-cell RNA-seq data. Use this skill for (1) querying and downloading expression data via CLI tools or programmatic API with flexible filtering by species, tissue, cell type, disease, and other metadata, including out-of-core processing and PyTorch integration for ML workflows, or (2) analyzing gene expression specificity and extracting cell type marker genes. This skill applies when users request Census data retrieval, marker gene analysis, specificity visualization, cohort building, or integration with analysis pipelines.
---

# CELLxGENE Census Toolkit

Comprehensive toolkit for querying and analyzing single-cell RNA-seq data from the CELLxGENE Census database.

## Overview

This skill provides two main capabilities:

### 1. **Query Subskill** (`query/`)
Download and filter single-cell RNA-seq data from CELLxGENE Census using CLI tools or programmatic API:
- **CLI Tools**: Interactive query scripts with confirmation prompts
- **Programmatic API**: Direct Python access for pipelines and ML workflows
- Query with flexible filters (species, tissue, cell type, disease, sex, development stage)
- Download expression matrices and comprehensive metadata
- Build custom cohorts for analysis
- Export in AnnData format (.h5ad)
- **28 cell metadata fields** including quality metrics
- **7 gene metadata fields** with Ensembl IDs and biotypes
- Analyze gene coexpression with correlation matrices and heatmaps
- Out-of-core processing for large-scale queries
- PyTorch integration for machine learning

### 2. **Specificity Subskill** (`specificity/`)
Analyze gene expression specificity and extract cell type markers:
- Extract computational marker genes (data-driven from Census)
- Extract canonical marker genes (literature-curated)
- Filter by organism, tissue, and marker quality
- Visualize with dotplots, heatmaps, and barplots
- Compare markers across organisms and tissues

## When to Use This Skill

Use this skill when users request:

**Data Retrieval** (→ use `query/` subskill):
- "Query single-cell data from CELLxGENE for [conditions]"
- "Download [cell type] from [tissue] tissue"
- "Get expression data for [genes] in [disease] samples"
- "Build a cohort of [cell type] matching [filters]"
- "What data is available for [tissue/cell type/disease]?"
- "Export cells with specific metadata criteria"
- "Integrate Census data into [analysis pipeline/ML workflow]"
- "Process large-scale Census queries out-of-core"
- "Train a model on Census data with PyTorch"

**Marker Analysis** (→ use `specificity/` subskill):
- "Find markers for [cell type]"
- "What are the top marker genes for [cell type]?"
- "Show me [cell type] markers specific to [tissue]"
- "Visualize marker genes for [cell type]"
- "Get canonical markers for [cell type]"
- "Compare [cell type] markers across organisms/tissues"

## Decision Guide

**Ask yourself**: What is the user trying to do?

### → Download/Query Data
User wants to:
- Get expression matrices
- Build a dataset/cohort
- Download cells matching criteria
- Export metadata
- Filter Census by multiple conditions

**Use**: `query/` subskill
**Tools**: `query_cellxgene.py`, `inspect_metadata_fields.py`

### → Find Marker Genes
User wants to:
- Identify cell type markers
- Visualize specificity
- Compare markers across contexts
- Get literature-validated markers
- Analyze marker quality

**Use**: `specificity/` subskill
**Tools**: `extract_markers.py`, `visualize_markers.py`

### → Both
Some requests require both:
1. Use `query/` to download relevant data
2. Use `specificity/` to analyze markers in that data

Example: "Download intestinal epithelial cells and find their marker genes"
1. `query_cellxgene.py --tissue intestine --cell-type epithelial`
2. `extract_markers.py --cell-type "intestinal epithelial cell"`

## Quick Start

### Query Workflow

```bash
# Download lung T cells
python query/scripts/query_cellxgene.py \
  --tissue lung \
  --cell-type "T cell" \
  --output lung_tcells

# Output: lung_tcells_data.h5ad, lung_tcells_metadata.csv,
#         lung_tcells_summary.json, lung_tcells_log.json
```

### Specificity Workflow

```bash
# Extract B cell markers
python specificity/scripts/extract_markers.py \
  --cell-type "B cell" \
  --organism "Homo sapiens" \
  --output bcell_markers \
  --summary

# Visualize markers
python specificity/scripts/visualize_markers.py \
  --input bcell_markers_computational.csv \
  --output bcell_dotplot.png \
  --title "B Cell Markers"
```

## Subskill Details

### Query Subskill (`query/`)

**Purpose**: Download expression data from CELLxGENE Census and analyze gene coexpression

**Two Approaches**:
1. **CLI Tools**: User-friendly scripts for interactive queries
2. **Programmatic API**: Python API for pipelines, ML workflows, and large-scale processing

**Key Features**:
- Flexible filtering by species, tissue, cell type, disease, sex, development stage
- Gene-specific queries to reduce data size
- Interactive download with summary preview (CLI)
- Out-of-core processing for large queries (API)
- PyTorch integration for machine learning (API)
- Comprehensive metadata export (28 cell fields, 7 gene fields)
- Quality metrics for cell filtering
- Gene coexpression analysis with correlation matrices
- Hierarchical clustering and heatmap visualization

**CLI Scripts**:
- `query_cellxgene.py` - Main query and download tool
- `inspect_metadata_fields.py` - Explore available metadata fields
- `example_query.py` - Example usage demonstrations
- `analyze_coexpression.py` - Gene coexpression analysis

**Programmatic API**:
- `cellxgene_census.open_soma()` - Context manager for Census access
- `cellxgene_census.get_anndata()` - Load data into AnnData objects
- `cellxgene_census.get_obs()`/`get_var()` - Metadata-only queries
- `axis_query()` - Out-of-core processing for large queries
- `experiment_dataloader()` - PyTorch integration

**Documentation**: See `query/SKILL.md` for detailed usage of both approaches

**Common Filters**:
- `--tissue`: lung, intestine, liver, brain, blood, etc.
- `--cell-type`: T cell, B cell, macrophage, epithelial cell, etc.
- `--disease`: normal, COVID-19, breast cancer, etc.
- `--species`: human (default), mouse
- `--development-stage`: adult (default), embryonic, fetal

**Example Outputs**:
```
query_data.h5ad          # AnnData with expression + metadata
query_metadata.csv       # 28 cell metadata fields
query_summary.json       # Aggregated statistics
query_log.json          # Query parameters
```

### Specificity Subskill (`specificity/`)

**Purpose**: Extract and visualize cell type marker genes

**Key Features**:
- Computational markers from Census data analysis
- Canonical markers from literature curation
- Quality filtering (marker_score >= 0.5 default)
- Multiple visualization types (dotplot, heatmap, barplot)
- Organism and tissue filtering

**Scripts**:
- `extract_markers.py` - Extract marker genes from CellGuide API
- `visualize_markers.py` - Create publication-quality visualizations

**Documentation**: See `specificity/SKILL.md` for detailed usage

**Key Metrics**:
- `marker_score`: Overall marker quality
- `specificity`: Cell type specificity (0-1)
- `me`: Mean expression level
- `pc`: Percentage of cells expressing

**Example Outputs**:
```
markers_computational.csv    # Computational markers with metrics
markers_canonical.csv        # Literature-curated markers
markers_dotplot.png          # Visualization
```

## Common Patterns

### Pattern 1: Tissue Survey

User: "What cell types are in lung tissue and what are their markers?"

```bash
# Step 1: Query to see what's available
python query/scripts/query_cellxgene.py \
  --tissue lung \
  --output lung_survey

# Check summary.json for cell type distribution

# Step 2: Extract markers for major cell types
for cell_type in "T cell" "B cell" "macrophage" "epithelial cell"; do
  python specificity/scripts/extract_markers.py \
    --cell-type "$cell_type" \
    --organism "Homo sapiens" \
    --tissue "lung" \
    --output "lung_${cell_type// /_}_markers"
done
```

### Pattern 2: Disease Cohort Analysis

User: "Build a COVID-19 lung dataset and identify macrophage markers"

```bash
# Step 1: Query disease cohort
python query/scripts/query_cellxgene.py \
  --tissue lung \
  --disease "COVID-19" \
  --cell-type macrophage \
  --output covid_macrophages

# Step 2: Extract markers for validation
python specificity/scripts/extract_markers.py \
  --cell-type "macrophage" \
  --organism "Homo sapiens" \
  --tissue "lung" \
  --output macrophage_markers \
  --summary
```

### Pattern 3: Gene-Specific Query

User: "Get expression of APOE family genes in liver cells"

```bash
# Query with gene set
python query/scripts/query_cellxgene.py \
  --tissue liver \
  --genes "APOE,APOC1,APOC2,APOC3,APOA1" \
  --output liver_apoe

# Analyze in Python
import scanpy as sc
adata = sc.read_h5ad("liver_apoe_data.h5ad")
sc.pl.violin(adata, keys=['APOE', 'APOC1'], groupby='cell_type')
```

### Pattern 4: Comparative Marker Analysis

User: "Compare B cell markers between human and mouse"

```bash
# Extract markers without organism filter to get both
python specificity/scripts/extract_markers.py \
  --cell-type "B cell" \
  --output bcell_multi \
  --summary

# Visualize with grouping
python specificity/scripts/visualize_markers.py \
  --input bcell_multi_computational.csv \
  --output bcell_comparison.png \
  --top-n 40 \
  --group-by organism_ontology_term_label \
  --title "B Cell Markers: Human vs Mouse"
```

### Pattern 5: Quality-Filtered Cohort

User: "Get high-quality epithelial cells from normal intestine"

```bash
# Step 1: Query data
python query/scripts/query_cellxgene.py \
  --tissue intestine \
  --cell-type epithelial \
  --disease normal \
  --output intestine_epithelial_normal

# Step 2: Filter by quality in Python
import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("intestine_epithelial_normal_data.h5ad")
metadata = pd.read_csv("intestine_epithelial_normal_metadata.csv")

# Filter high quality
high_quality = metadata[
    (metadata['n_measured_vars'] > 10000) &
    (metadata['raw_sum'] > 5000)
]

# Subset AnnData
adata = adata[high_quality.index]
sc.write("intestine_epithelial_high_quality.h5ad", adata)
```

## Reference Documentation

Detailed documentation available in `references/` directory:

**Query Subskill - CLI Tools**:
- `CELLXGENE_QUERY_QUICKSTART.md` - Quick reference with examples
- `QUERY_README.md` - Comprehensive usage guide
- `METADATA_REFERENCE.md` - Metadata field descriptions and usage
- `METADATA_FIELDS_ACTUAL.md` - Complete field listing (28 cell + 7 gene)
- `COEXPRESSION_ANALYSIS.md` - Gene coexpression analysis guide

**Query Subskill - Programmatic API**:
- `census_schema.md` - Census data structure, metadata fields, and filter syntax
- `common_patterns.md` - API usage patterns, best practices, and code examples

**Specificity Subskill**:
- See `specificity/SKILL.md` for complete documentation

## Environment Requirements

Python 3.7+ with:
- **For Query**: cellxgene-census, anndata, pandas, numpy, tiledbsoma
- **For Specificity**: pandas, requests, matplotlib, seaborn, numpy

Install:
```bash
# For query functionality
pip install cellxgene-census anndata pandas numpy

# For specificity functionality
pip install pandas requests matplotlib seaborn numpy

# Or install all
pip install cellxgene-census anndata pandas numpy requests matplotlib seaborn
```

## Data Source

All data comes from **CELLxGENE Census**, a comprehensive single-cell RNA-seq database:
- Aggregates data from multiple studies
- Standardized cell type annotations (Cell Ontology)
- Harmonized metadata across experiments
- Quality-controlled expression matrices
- Regular updates with new data

**CellGuide API**: Provides marker gene annotations and specificity metrics
**Census API**: Provides raw expression data and metadata

Census Documentation: https://chanzuckerberg.github.io/cellxgene-census/
CellGuide: https://cellguide.cellxgene.cziscience.com

## Tips and Best Practices

1. **Start with inspection**: Use `inspect_metadata_fields.py` to see available values before querying
2. **Use specific filters**: Tissue + cell type filters keep queries manageable
3. **Check summaries first**: Review summary before downloading large datasets
4. **Gene-specific queries**: Use `--genes` to reduce data size for targeted analyses
5. **Default filtering**: Specificity tools default to marker_score >= 0.5 for quality
6. **Validate markers**: Cross-reference computational and canonical markers
7. **Quality metrics**: Filter cells by `n_measured_vars` and `raw_sum` post-download
8. **Save query logs**: Track analyses with log.json files
9. **Organism defaults**: Query defaults to human, specify `--species mouse` for mouse
10. **Development stage**: Query defaults to adult, specify `--development-stage` for other stages

## Troubleshooting

**Issue**: Don't know which subskill to use
**Solution**: Use `query/` for data download, `specificity/` for marker analysis

**Issue**: Query returns no results
**Solution**: Filters too restrictive. Remove some filters or use broader terms

**Issue**: Can't find cell type for markers
**Solution**: Try alternative names or find Cell Ontology ID manually

**Issue**: Download too large
**Solution**: Add more filters or use `--genes` to limit to specific genes

**Issue**: Marker visualization empty
**Solution**: Check that CSV has data and marker_score values

**Issue**: Out of memory
**Solution**: Query smaller subsets by adding more restrictive filters

## Integration

Both subskills integrate with standard analysis tools:

**Scanpy** (Python):
```python
import scanpy as sc
adata = sc.read_h5ad("query_data.h5ad")
# Standard Scanpy workflow
```

**Seurat** (R):
```r
library(Seurat)
library(SeuratDisk)
Convert("query_data.h5ad", dest = "h5seurat")
```

**Downstream analyses**:
- Differential expression
- Cell type annotation
- Trajectory inference
- Gene set enrichment
- Regulatory networks

## Example Full Workflow

```bash
# User: "Analyze human intestinal epithelial cells and their markers"

# 1. Inspect available data
python query/scripts/inspect_metadata_fields.py --sample-size 1000

# 2. Query intestinal epithelial cells
python query/scripts/query_cellxgene.py \
  --tissue intestine \
  --cell-type epithelial \
  --organism human \
  --disease normal \
  --output intestine_epithelial

# Reviews summary, confirms download
# Output: intestine_epithelial_data.h5ad + metadata + summary + log

# 3. Extract marker genes
python specificity/scripts/extract_markers.py \
  --cell-type "intestinal epithelial cell" \
  --organism "Homo sapiens" \
  --output intestine_epi_markers \
  --summary

# 4. Visualize markers
python specificity/scripts/visualize_markers.py \
  --input intestine_epi_markers_computational.csv \
  --output intestine_epi_dotplot.png \
  --top-n 30 \
  --title "Intestinal Epithelial Cell Markers"

# 5. Analyze in Python
import scanpy as sc
adata = sc.read_h5ad("intestine_epithelial_data.h5ad")
markers = pd.read_csv("intestine_epi_markers_computational.csv")

# Check marker expression in dataset
marker_genes = markers['symbol'].head(20).tolist()
sc.pl.dotplot(adata, var_names=marker_genes, groupby='cell_type')

# Result: User has dataset, marker genes, visualizations, and can proceed with analysis
```

## Support

For detailed usage of each subskill:
- **Query**: Read `query/SKILL.md` and `references/CELLXGENE_QUERY_QUICKSTART.md`
- **Specificity**: Read `specificity/SKILL.md`

For CELLxGENE Census issues:
- Census Documentation: https://chanzuckerberg.github.io/cellxgene-census/
- CellGuide API: https://cellguide.cellxgene.cziscience.com
