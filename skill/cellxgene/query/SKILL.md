---
name: cellxgene-query
description: Query and download single-cell RNA-seq data from CELLxGENE Census with flexible filtering by species, tissue, cell type, disease, sex, development stage, and gene sets. Returns AnnData objects with expression matrices and comprehensive cell metadata. Use this for data retrieval, cohort building, and exploratory analysis of Census data.
---

# CELLxGENE Census Query

Query and download single-cell RNA-seq data from the CELLxGENE Census database with comprehensive filtering options.

## Overview

This subskill provides tools to:
1. Query CELLxGENE Census database with flexible filters
2. Download expression matrices and cell metadata
3. Generate summary statistics before download
4. Export data in AnnData format (.h5ad)
5. Extract and inspect metadata fields
6. Analyze gene coexpression patterns with correlation matrices and heatmaps

**Key Features**:
- Filter by species, tissue, cell type, disease, sex, development stage
- Case-insensitive partial matching for flexible queries
- Interactive download confirmation with summary preview
- Comprehensive metadata export (28 cell fields, 7 gene fields)
- Quality metrics for cell filtering
- Pairwise gene correlation analysis (Pearson/Spearman)
- Hierarchical clustering and heatmap visualization

## When to Use This Skill

Use this skill when users request:
- "Query single-cell data from CELLxGENE for [conditions]"
- "Download [cell type] from [tissue] tissue"
- "Get expression data for [genes] in [condition]"
- "Build a cohort of [cell type] from [disease] samples"
- "What data is available for [tissue/cell type/disease]?"
- "Download cells matching [multiple filters]"
- "Analyze coexpression of [gene list]"
- "Generate correlation matrix for [genes] in [cell type]"
- "Show me how [genes] correlate in [condition]"

## Core Concepts

### CELLxGENE Census

The CELLxGENE Census is a comprehensive collection of single-cell RNA-seq data:
- Aggregates data from multiple studies and datasets
- Standardized cell type annotations (Cell Ontology)
- Harmonized metadata across experiments
- Quality-controlled expression matrices
- Regular updates with new data

### Filtering Logic

**String Input** (single value):
- Partial case-insensitive match
- Example: `--tissue lung` matches "lung", "left lung", "lung parenchyma"

**List Input** (comma-separated):
- Any match (OR logic)
- Example: `--tissue "lung,intestine"` matches cells from either tissue

**Defaults**:
- Species: `human` (if not specified)
- Development stage: `adult` (if not specified)
- All other filters: None (no filtering)

### Output Files

Each query produces 4 files:
1. **`{prefix}_data.h5ad`** - AnnData object with expression matrix and metadata
2. **`{prefix}_metadata.csv`** - Cell metadata (28 fields) in CSV format
3. **`{prefix}_summary.json`** - Aggregated statistics and counts
4. **`{prefix}_log.json`** - Query parameters and timestamps

## Workflow

### Step 1: Query Data

Use `scripts/query_cellxgene.py` to query and download data from Census.

**Basic Query** (human adult cells):
```bash
python scripts/query_cellxgene.py \
  --tissue lung \
  --output lung_adult
```

**With Cell Type Filter**:
```bash
python scripts/query_cellxgene.py \
  --tissue intestine \
  --cell-type "epithelial cell" \
  --output intestine_epithelial
```

**With Multiple Filters**:
```bash
python scripts/query_cellxgene.py \
  --tissue lung \
  --cell-type "T cell" \
  --disease "COVID-19" \
  --sex male \
  --output covid_lung_tcells
```

**With Gene Set** (reduces data size):
```bash
python scripts/query_cellxgene.py \
  --tissue liver \
  --genes "APOE,APOC1,APOC2,APOC3" \
  --output liver_apoe_genes
```

**Multi-Tissue Query**:
```bash
python scripts/query_cellxgene.py \
  --tissue "lung,intestine,liver" \
  --cell-type macrophage \
  --output macrophage_multi_tissue
```

**Multiple Species**:
```bash
python scripts/query_cellxgene.py \
  --species "human,mouse" \
  --tissue brain \
  --cell-type neuron \
  --output brain_neurons
```

**Development Stage**:
```bash
python scripts/query_cellxgene.py \
  --tissue brain \
  --development-stage embryonic \
  --output embryonic_brain
```

**Non-Interactive Mode** (skip confirmation):
```bash
python scripts/query_cellxgene.py \
  --tissue lung \
  --output lung_data \
  --no-interactive
```

### Step 2: Inspect Metadata Fields

Use `scripts/inspect_metadata_fields.py` to explore available metadata fields.

**Inspect Human Data**:
```bash
python scripts/inspect_metadata_fields.py --sample-size 1000
```

**Inspect Mouse Data**:
```bash
python scripts/inspect_metadata_fields.py \
  --organism mus_musculus \
  --sample-size 1000
```

This displays:
- All observation (cell) metadata fields
- All variable (gene) metadata fields
- Example values for each field
- Value distributions for categorical fields

### Step 3: Analyze Results

After download, work with the data:

**Load AnnData**:
```python
import scanpy as sc
adata = sc.read_h5ad("lung_tcells_data.h5ad")
print(adata)
# AnnData object with n_obs × n_vars = 50000 × 36000
```

**Load Metadata CSV**:
```python
import pandas as pd
metadata = pd.read_csv("lung_tcells_metadata.csv")
print(metadata.columns)
# cell_type, tissue_general, disease, sex, donor_id, ...
```

**Load Summary JSON**:
```python
import json
with open("lung_tcells_summary.json") as f:
    summary = json.load(f)
print(f"Total cells: {summary['n_cells']:,}")
print(f"Cell types: {summary['cell_types']}")
```

### Step 4: Analyze Gene Coexpression

Use `scripts/analyze_coexpression.py` to assess pairwise gene coexpression patterns.

**Basic Coexpression Analysis**:
```bash
python scripts/analyze_coexpression.py \
  --input lung_tcells_data.h5ad \
  --genes "CD3D,CD3E,CD4,CD8A,CD8B" \
  --output tcell_markers_coexp
```

**With Gene List File**:
```bash
# genes.txt contains one gene per line
python scripts/analyze_coexpression.py \
  --input data.h5ad \
  --genes my_genes.txt \
  --output gene_coexp
```

**Filter by Cell Type**:
```bash
python scripts/analyze_coexpression.py \
  --input lung_data.h5ad \
  --genes immune_genes.txt \
  --cell-type "T cell" \
  --output tcell_coexp
```

**Spearman Correlation with Clustering**:
```bash
python scripts/analyze_coexpression.py \
  --input data.h5ad \
  --genes genes.txt \
  --method spearman \
  --cluster \
  --output spearman_coexp
```

**Custom Metadata Filtering**:
```bash
python scripts/analyze_coexpression.py \
  --input data.h5ad \
  --genes genes.txt \
  --tissue lung \
  --disease "COVID-19" \
  --metadata-filter "sex == 'male'" \
  --output covid_male_coexp
```

**Output Files**:
- `{prefix}_correlation_matrix.csv` - Pairwise correlation matrix
- `{prefix}_pvalue_matrix.csv` - Pairwise p-value matrix
- `{prefix}_heatmap_YYYYMMDD.png` - Correlation heatmap with significance markers (date-stamped)
- `{prefix}_summary.txt` - Summary report with correlation and p-value statistics

**Heatmap Features**:
- **BWR colormap**: Blue (negative correlation) - White (no correlation) - Red (positive correlation)
- **Significance markers**: Annotations show correlation values with significance stars
  - `****` p < 0.0001
  - `***` p < 0.001
  - `**` p < 0.01
  - `*` p < 0.05
- **Hierarchical clustering**: Optional dendrogram-based gene grouping

## Common Usage Patterns

### Pattern 1: Quick Tissue Query

User asks: "Get lung single-cell data"

```bash
python scripts/query_cellxgene.py \
  --tissue lung \
  --output lung_adult
```

### Pattern 2: Specific Cell Type

User asks: "Download intestinal epithelial cells"

```bash
python scripts/query_cellxgene.py \
  --tissue intestine \
  --cell-type epithelial \
  --output intestine_epithelial
```

### Pattern 3: Disease Cohort

User asks: "Get COVID-19 lung samples"

```bash
python scripts/query_cellxgene.py \
  --tissue lung \
  --disease "COVID-19" \
  --output covid_lung
```

### Pattern 4: Gene-Specific Query

User asks: "Get liver expression data for APOE and related genes"

```bash
python scripts/query_cellxgene.py \
  --tissue liver \
  --genes "APOE,APOC1,APOC2,APOC3,APOA1" \
  --output liver_apoe_family
```

### Pattern 5: Comparative Analysis Setup

User asks: "Get T cells from multiple tissues for comparison"

```bash
python scripts/query_cellxgene.py \
  --tissue "lung,intestine,spleen,blood" \
  --cell-type "T cell" \
  --output tcells_multi_tissue
```

### Pattern 6: Quality-Filtered Cohort

User asks: "Get high-quality epithelial cells from normal lung tissue"

```bash
# Step 1: Query data
python scripts/query_cellxgene.py \
  --tissue lung \
  --cell-type epithelial \
  --disease normal \
  --output lung_epithelial_normal

# Step 2: Filter by quality in Python
import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("lung_epithelial_normal_data.h5ad")
metadata = pd.read_csv("lung_epithelial_normal_metadata.csv")

# Filter high quality cells
high_quality = metadata[
    (metadata['n_measured_vars'] > 10000) &
    (metadata['raw_sum'] > 5000)
]

print(f"High quality cells: {len(high_quality)} / {len(metadata)}")
```

### Pattern 7: Gene Coexpression Analysis

User asks: "Analyze coexpression of APOE family genes in liver"

```bash
# Step 1: Query liver data
python scripts/query_cellxgene.py \
  --tissue liver \
  --disease normal \
  --output liver_normal

# Step 2: Analyze coexpression of APOE family
python scripts/analyze_coexpression.py \
  --input liver_normal_data.h5ad \
  --genes "APOE,APOC1,APOC2,APOC3,APOA1" \
  --method pearson \
  --cluster \
  --title "APOE Family Coexpression in Liver" \
  --output apoe_coexp

# Outputs:
# - apoe_coexp_correlation_matrix.csv (pairwise correlations)
# - apoe_coexp_pvalue_matrix.csv (pairwise p-values)
# - apoe_coexp_heatmap_20251215.png (clustered heatmap with significance stars)
# - apoe_coexp_summary.txt (statistics, p-values, and top correlations)
```

Alternative with cell type filtering:

```bash
# Analyze only in hepatocytes
python scripts/analyze_coexpression.py \
  --input liver_normal_data.h5ad \
  --genes "APOE,APOC1,APOC2,APOC3,APOA1" \
  --cell-type hepatocyte \
  --method pearson \
  --cluster \
  --output apoe_hepatocyte_coexp
```

## Script Reference

### query_cellxgene.py

**Purpose**: Query CELLxGENE Census and download data

**Key Arguments**:
- `--species`: Species name(s), comma-separated (default: human)
- `--tissue`: Tissue name(s), comma-separated
- `--cell-type`: Cell type name(s), comma-separated
- `--disease`: Disease name(s), comma-separated
- `--sex`: Sex, comma-separated
- `--development-stage`: Development stage(s), comma-separated (default: adult)
- `--genes`: Gene symbols, comma-separated (reduces data size)
- `--output`: Output file prefix (required)
- `--output-dir`: Output directory (default: current directory)
- `--no-interactive`: Skip confirmation prompt

**Workflow**:
1. Applies all filters to Census metadata
2. Shows summary statistics (cell count, tissues, cell types, etc.)
3. Asks for download confirmation (unless --no-interactive)
4. Downloads expression matrix and metadata
5. Saves 4 files: .h5ad, metadata.csv, summary.json, log.json

**Example**:
```bash
python scripts/query_cellxgene.py \
  --tissue lung \
  --cell-type "epithelial cell" \
  --disease normal \
  --output lung_epithelial
```

### inspect_metadata_fields.py

**Purpose**: Inspect available metadata fields in Census

**Key Arguments**:
- `--organism`: Organism to inspect (default: homo_sapiens, also: mus_musculus)
- `--sample-size`: Number of cells to sample (default: 100)

**Output**:
- Lists all observation (cell) metadata fields
- Lists all variable (gene) metadata fields
- Shows example values and data types
- Displays value distributions for key categorical fields

**Example**:
```bash
python scripts/inspect_metadata_fields.py --sample-size 1000
```

### example_query.py

**Purpose**: Demonstrate API usage with 5 examples

**Examples**:
1. Simple tissue query
2. Cell type across multiple tissues
3. Gene set query
4. Disease query
5. Multi-species query

**Usage**:
```bash
# Run specific example
python scripts/example_query.py --example 1

# Run all examples (warning: large downloads!)
python scripts/example_query.py
```

### analyze_coexpression.py

**Purpose**: Analyze pairwise gene coexpression from AnnData objects

**Key Arguments**:
- `--input`: Input AnnData file (.h5ad) (required)
- `--genes`: Gene list file (one per line) or comma-separated names (required)
- `--output`: Output prefix for files (default: coexpression)
- `--output-dir`: Output directory (default: current directory)
- `--method`: Correlation method - pearson or spearman (default: pearson)
- `--use-raw`: Use .raw.X instead of .X
- `--min-cells`: Minimum cells expressing each gene (default: 10)

**Filtering Arguments**:
- `--cell-type`: Filter by cell type (partial match)
- `--tissue`: Filter by tissue (partial match)
- `--disease`: Filter by disease (partial match)
- `--metadata-filter`: Custom filter expression (e.g., "sex == 'female'")

**Visualization Arguments**:
- `--figsize`: Figure size (width height) (default: 10 8)
- `--cmap`: Colormap for heatmap (default: RdBu_r)
- `--cluster`: Perform hierarchical clustering
- `--no-annot`: Disable correlation value annotations
- `--title`: Custom plot title

**Workflow**:
1. Loads AnnData from .h5ad file
2. Applies optional metadata filters
3. Extracts expression for specified genes
4. Computes pairwise correlations (Pearson or Spearman)
5. Generates correlation matrix CSV
6. Creates heatmap with optional hierarchical clustering
7. Produces summary report with statistics

**Outputs**:
- `{prefix}_correlation_matrix.csv` - Pairwise correlation matrix
- `{prefix}_pvalue_matrix.csv` - Pairwise p-value matrix
- `{prefix}_heatmap_YYYYMMDD.png` - Correlation heatmap with significance stars (300 DPI, date-stamped)
- `{prefix}_summary.txt` - Summary report with correlation and p-value statistics

**Example**:
```bash
python scripts/analyze_coexpression.py \
  --input lung_tcells_data.h5ad \
  --genes "CD3D,CD3E,CD4,CD8A,CD8B,GZMB,PRF1" \
  --method pearson \
  --cluster \
  --output tcell_markers_coexp
```

**Advanced Example with Filtering**:
```bash
python scripts/analyze_coexpression.py \
  --input lung_data.h5ad \
  --genes immune_genes.txt \
  --cell-type "T cell" \
  --disease "COVID-19" \
  --method spearman \
  --cluster \
  --title "T Cell Coexpression in COVID-19" \
  --output covid_tcell_coexp
```

## Data Structure

### Cell Metadata (28 fields)

**Cell Identity**:
- `cell_type` - Cell type annotation
- `cell_type_ontology_term_id` - Cell Ontology ID
- `soma_joinid` - Unique cell ID
- `donor_id` - Patient/donor identifier
- `dataset_id` - Source dataset ID

**Tissue/Anatomy**:
- `tissue` - Detailed tissue name
- `tissue_general` - General tissue category
- `tissue_ontology_term_id` - UBERON ID
- `suspension_type` - "cell" or "nucleus"

**Donor Characteristics**:
- `sex` - female, male, unknown
- `development_stage` - e.g., "29-year-old stage", "adult"
- `self_reported_ethnicity` - Ethnicity
- `disease` - Disease annotation

**Experimental**:
- `assay` - Sequencing method (e.g., "10x 3' v3")
- `assay_ontology_term_id` - EFO ID

**Quality Metrics**:
- `n_measured_vars` - Number of genes measured
- `nnz` - Non-zero gene counts
- `raw_sum` - Total UMI counts
- `raw_mean_nnz` - Mean expression
- `raw_variance_nnz` - Expression variance
- `is_primary_data` - True/False

### Gene Metadata (7 fields)

- `feature_id` - Ensembl gene ID
- `feature_name` - Gene symbol
- `feature_type` - Gene biotype (protein_coding, lncRNA, etc.)
- `feature_length` - Gene length in bp
- `n_measured_obs` - Number of cells with this gene measured
- `nnz` - Non-zero observations

### Summary Statistics Structure

```json
{
  "timestamp": "2025-12-05T10:30:45",
  "filters_applied": "(tissue_general == \"lung\")",
  "n_cells": 125432,
  "n_genes": 36601,
  "donors": 45,
  "organisms": {"Homo sapiens": 125432},
  "tissues": {"lung": 98234, "lung parenchyma": 15678},
  "cell_types": {"T cell": 58234, "macrophage": 35678},
  "diseases": {"normal": 98234, "COVID-19": 15678},
  "sexes": {"male": 65432, "female": 60000},
  "development_stages": {"adult": 120000},
  "assays": {"10x 3' v3": 85432}
}
```

## Tips and Best Practices

1. **Start restrictive**: Use tissue + cell type filters to keep data manageable
2. **Review summary first**: Check cell counts before downloading
3. **Use gene sets**: Specify genes with `--genes` to reduce data size significantly
4. **Check quality metrics**: Filter cells by `n_measured_vars`, `raw_sum` after download
5. **Save query logs**: The log.json file helps reproduce analyses
6. **Species default**: Defaults to human; specify `--species mouse` for mouse data
7. **Development stage**: Defaults to adult; use `--development-stage embryonic` for developmental data
8. **Partial matching**: Leverage partial matching for flexible queries (e.g., "epithelial" matches "intestinal epithelial cell")
9. **Non-interactive mode**: Use `--no-interactive` for automated pipelines
10. **Metadata exploration**: Use `inspect_metadata_fields.py` to see what's available before querying

## Troubleshooting

**Issue**: Query returns no results
**Solution**: Filters may be too restrictive. Remove some filters or check filter spelling. Try broader terms (e.g., "epithelial" instead of "intestinal epithelial cell")

**Issue**: Download is very slow
**Solution**: Add more restrictive filters (tissue, cell type, development stage) or use `--genes` to limit to specific genes

**Issue**: Out of memory errors
**Solution**: Query smaller subsets by adding more filters, or query specific genes only

**Issue**: "Organism not found"
**Solution**: Use `homo_sapiens` or `human` for human data, `mus_musculus` or `mouse` for mouse data

**Issue**: Empty summary statistics
**Solution**: Check that filters are correctly spelled and match available values. Use `inspect_metadata_fields.py` to see available values

## Reference Documentation

Detailed reference documentation is available in the `references/` directory:

- **CELLXGENE_QUERY_QUICKSTART.md** - Quick reference with copy-paste examples
- **QUERY_README.md** - Comprehensive usage guide with all options
- **METADATA_REFERENCE.md** - Detailed metadata field descriptions and usage
- **METADATA_FIELDS_ACTUAL.md** - Complete list of all 28 cell and 7 gene metadata fields

To access these files:
```bash
# Read in terminal
cat references/CELLXGENE_QUERY_QUICKSTART.md

# Or open in editor
```

## Example Session

```bash
# User request: "Download human lung T cells for analysis"

# Query data
python scripts/query_cellxgene.py \
  --tissue lung \
  --cell-type "T cell" \
  --organism human \
  --output lung_tcells

# Output shows:
# ================================================================================
# QUERY RESULTS SUMMARY
# ================================================================================
# Total Cells: 125,432
# Total Genes: 36,601
# Unique Donors: 45
#
# Cell Types:
#   T cell: 98,234
#   CD4-positive, alpha-beta T cell: 15,678
#   CD8-positive, alpha-beta T cell: 11,520
# ...
#
# Do you want to download and save this data? (yes/no): yes
#
# Saving AnnData object to: lung_tcells_data.h5ad
# Saving metadata to: lung_tcells_metadata.csv
# Saving summary to: lung_tcells_summary.json
# Saving query log to: lung_tcells_log.json
# All results saved successfully!

# Analyze in Python
import scanpy as sc
adata = sc.read_h5ad("lung_tcells_data.h5ad")
sc.pp.filter_cells(adata, min_genes=1000)
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['cell_type', 'disease', 'sex'])
```

## Integration with Analysis Workflows

Query results integrate seamlessly with standard single-cell analysis tools:

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
seurat_obj <- LoadH5Seurat("query_data.h5seurat")
```

**Downstream analyses**:
- Differential expression
- Cell type annotation
- Trajectory inference
- Gene set enrichment
- Regulatory network analysis

## Environment Requirements

Python 3.7+ with:
- cellxgene-census
- anndata
- pandas
- numpy
- tiledbsoma

Install via:
```bash
conda activate your_env
pip install cellxgene-census anndata pandas numpy
```

## Data Source

All data comes from CELLxGENE Census, which aggregates and standardizes single-cell RNA-seq data from multiple sources with:
- Quality control and validation
- Standardized cell type ontologies
- Harmonized metadata
- Regular updates

Census Documentation: https://chanzuckerberg.github.io/cellxgene-census/
