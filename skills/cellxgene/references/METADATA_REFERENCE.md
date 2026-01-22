# CELLxGENE Query Metadata & Summary Reference

## Overview

When you run a CELLxGENE query, 4 files are generated:
1. **`{prefix}_data.h5ad`** - AnnData object with expression matrix and all metadata
2. **`{prefix}_metadata.csv`** - Cell metadata in CSV format (easier to browse)
3. **`{prefix}_summary.json`** - Aggregated statistics about the query results
4. **`{prefix}_log.json`** - Query parameters and timestamps

This document explains what's in the metadata CSV and summary JSON files.

---

## Cell Metadata CSV (`{prefix}_metadata.csv`)

The metadata CSV contains **one row per cell** with detailed annotations. The exact fields depend on the CELLxGENE Census version, but typically include:

### Core Identification Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `soma_joinid` | int64 | Unique cell identifier | 123456789 |
| `obs_id` | string | Original observation ID from dataset | AAACCTGAGCGCTCCA-1 |
| `dataset_id` | string | Unique dataset identifier | abc123-def456 |
| `donor_id` | string | Unique donor/patient identifier | D001 |
| `cell_type` | string | Cell type annotation | CD4-positive, alpha-beta T cell |
| `cell_type_ontology_term_id` | string | Cell Ontology ID | CL:0000624 |

### Biological Context Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `organism` | string | Species name | Homo sapiens |
| `organism_ontology_term_id` | string | NCBITaxon ID | NCBITaxon:9606 |
| `tissue` | string | Detailed tissue name | lung parenchyma |
| `tissue_general` | string | General tissue category | lung |
| `tissue_ontology_term_id` | string | UBERON tissue ID | UBERON:0002048 |
| `assay` | string | Experimental technique | 10x 3' v3 |
| `assay_ontology_term_id` | string | EFO assay ID | EFO:0009922 |

### Donor/Sample Metadata

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `sex` | string | Biological sex | female, male, unknown |
| `sex_ontology_term_id` | string | PATO sex ID | PATO:0000383 |
| `ethnicity` | string | Self-reported ethnicity | European |
| `ethnicity_ontology_term_id` | string | HANCESTRO ethnicity ID | HANCESTRO:0005 |
| `development_stage` | string | Developmental stage | adult |
| `development_stage_ontology_term_id` | string | HsapDv stage ID | HsapDv:0000087 |

### Disease & Health Status

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `disease` | string | Disease annotation | normal, COVID-19 |
| `disease_ontology_term_id` | string | MONDO/PATO disease ID | PATO:0000461 |
| `self_reported_ethnicity` | string | Self-reported ethnicity | White |

### Data Provenance

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `is_primary_data` | bool | Whether data is primary or derived | True, False |
| `suspension_type` | string | Sample preparation type | cell, nucleus |
| `observation_joinid` | int64 | Join ID for linking | 123456789 |

### Dataset Information

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `dataset_title` | string | Dataset publication title | Single-cell atlas of... |
| `dataset_h5ad_path` | string | Path to original dataset | datasets/abc123.h5ad |
| `dataset_total_cell_count` | int64 | Total cells in source dataset | 50000 |

### Additional Fields (when available)

| Field | Type | Description |
|-------|------|-------------|
| `cell_culture` | string | Whether cells were cultured |
| `tissue_type` | string | Tissue type classification |
| `bmi` | float | Body mass index (if available) |
| `age` | string | Age or age range |
| `anatomical_information` | string | Additional anatomical details |

### Example Row

```csv
soma_joinid,cell_type,tissue_general,organism,sex,disease,development_stage,assay,donor_id
12345,CD4-positive alpha-beta T cell,lung,Homo sapiens,male,normal,adult,10x 3' v3,D001
```

---

## Summary Statistics JSON (`{prefix}_summary.json`)

The summary JSON provides **aggregated statistics** across all cells in your query. This gives you a quick overview without loading the full dataset.

### Structure

```json
{
  "timestamp": "2025-12-05T10:30:45.123456",
  "filters_applied": "(tissue_general == \"lung\") and (cell_type == \"T cell\")",
  "n_cells": 125432,
  "n_genes": 36601,
  "donors": 45,
  "organisms": {
    "Homo sapiens": 125432
  },
  "tissues": {
    "lung": 98234,
    "lung parenchyma": 15678,
    "left lung": 11520
  },
  "cell_types": {
    "T cell": 58234,
    "CD4-positive, alpha-beta T cell": 35678,
    "CD8-positive, alpha-beta T cell": 20120,
    "regulatory T cell": 8400,
    "memory T cell": 3000
  },
  "diseases": {
    "normal": 98234,
    "COVID-19": 15678,
    "ARDS": 11520
  },
  "sexes": {
    "male": 65432,
    "female": 60000
  },
  "development_stages": {
    "adult": 120000,
    "child": 5432
  },
  "assays": {
    "10x 3' v3": 85432,
    "10x 5' v2": 25000,
    "Smart-seq2": 15000
  }
}
```

### Fields Explained

| Field | Type | Description |
|-------|------|-------------|
| `timestamp` | string | When the query was run (ISO format) |
| `filters_applied` | string | SOMA filter string used in query |
| `n_cells` | int | Total number of cells retrieved |
| `n_genes` | int | Total number of genes in expression matrix |
| `donors` | int | Number of unique donors/patients |
| `organisms` | dict | Cell count per species |
| `tissues` | dict | Cell count per tissue type |
| `cell_types` | dict | Cell count per cell type annotation |
| `diseases` | dict | Cell count per disease status |
| `sexes` | dict | Cell count per sex |
| `development_stages` | dict | Cell count per developmental stage |
| `assays` | dict | Cell count per experimental assay |

---

## Query Log JSON (`{prefix}_log.json`)

Tracks all queries performed, useful for reproducing analyses.

```json
[
  {
    "timestamp": "2025-12-05T10:30:45.123456",
    "filters": {
      "species": ["homo sapiens"],
      "tissue": ["lung"],
      "cell_type": ["t cell"],
      "disease": null,
      "sex": null,
      "development_stage": ["adult"],
      "drug_treatment": null,
      "gene_set": null
    }
  }
]
```

---

## Inspecting Available Fields

To see all available metadata fields in your environment:

```bash
conda activate claude_test
cd /home/sagemaker-user/claude_code/gene-expression-specificity-cellxgene

# Inspect fields (connects to Census and shows all columns)
conda run -n claude_test python scripts/inspect_metadata_fields.py

# Inspect for different organism
conda run -n claude_test python scripts/inspect_metadata_fields.py --organism "Mus musculus"

# Sample more cells for better overview
conda run -n claude_test python scripts/inspect_metadata_fields.py --sample-size 1000
```

This will display:
- All observation (cell) metadata columns
- All variable (gene) metadata columns
- Example values for each field
- Value distributions for key categorical fields

---

## Loading and Using Metadata

### Load Metadata CSV
```python
import pandas as pd

# Read metadata
metadata = pd.read_csv("lung_tcells_metadata.csv")

# Explore
print(f"Total cells: {len(metadata)}")
print(f"Columns: {list(metadata.columns)}")

# Group by cell type
cell_type_counts = metadata['cell_type'].value_counts()
print(cell_type_counts)

# Filter for specific donor
donor_cells = metadata[metadata['donor_id'] == 'D001']

# Check disease distribution
disease_dist = metadata.groupby(['disease', 'sex']).size()
print(disease_dist)
```

### Load Summary JSON
```python
import json

# Read summary
with open("lung_tcells_summary.json", 'r') as f:
    summary = json.load(f)

# Access statistics
print(f"Total cells: {summary['n_cells']:,}")
print(f"Total genes: {summary['n_genes']:,}")
print(f"Unique donors: {summary['donors']}")

# Top cell types
for cell_type, count in sorted(summary['cell_types'].items(),
                                key=lambda x: x[1],
                                reverse=True)[:5]:
    print(f"{cell_type}: {count:,}")
```

### Load Full AnnData
```python
import scanpy as sc

# Read AnnData (includes expression + all metadata)
adata = sc.read_h5ad("lung_tcells_data.h5ad")

# Metadata is in .obs
print(adata.obs.columns)
print(adata.obs['cell_type'].value_counts())

# Gene metadata is in .var
print(adata.var.columns)

# Expression matrix
print(adata.X)  # Sparse matrix (cells × genes)
```

---

## Common Metadata Queries

### 1. Count Cells per Donor
```python
metadata.groupby('donor_id').size().sort_values(ascending=False)
```

### 2. Disease × Sex Cross-tabulation
```python
pd.crosstab(metadata['disease'], metadata['sex'])
```

### 3. Find Rare Cell Types
```python
cell_type_counts = metadata['cell_type'].value_counts()
rare_types = cell_type_counts[cell_type_counts < 100]
print(rare_types)
```

### 4. Filter by Multiple Conditions
```python
# Adult females with normal tissue
subset = metadata[
    (metadata['development_stage'] == 'adult') &
    (metadata['sex'] == 'female') &
    (metadata['disease'] == 'normal')
]
```

### 5. Unique Datasets
```python
print(f"Data from {metadata['dataset_id'].nunique()} datasets")
metadata['dataset_title'].unique()
```

---

## Field Completeness

Not all fields are available for all cells. Some fields may be:
- **Missing** (NaN) - not collected for that sample
- **"unknown"** - explicitly marked as unknown
- **"na"** - not applicable

Always check for missing data:
```python
# Check completeness
metadata.info()

# Count missing values
metadata.isnull().sum()

# Filter out unknowns
clean_metadata = metadata[metadata['sex'] != 'unknown']
```

---

## Tips

1. **Start with Summary**: Check the summary JSON first to understand your data without loading large files

2. **Use Metadata for Filtering**: The metadata CSV is much smaller than the full AnnData and easier to explore

3. **Check Ontology IDs**: Ontology term IDs provide standardized identifiers that link to external databases (Cell Ontology, UBERON, etc.)

4. **Validate Filters**: The `filters_applied` field in summary shows exactly what filters were used

5. **Track Provenance**: Use `dataset_id` and `donor_id` to track cell origins

6. **Mind the Size**: Large queries can produce huge metadata files. Use filters to keep queries manageable

---

## Additional Resources

- **CELLxGENE Census Documentation**: https://chanzuckerberg.github.io/cellxgene-census/
- **Cell Ontology**: https://www.ebi.ac.uk/ols/ontologies/cl
- **UBERON Anatomy**: https://www.ebi.ac.uk/ols/ontologies/uberon
- **Experimental Factor Ontology**: https://www.ebi.ac.uk/ols/ontologies/efo

---

## Quick Reference

| File | Content | Use Case |
|------|---------|----------|
| `*_data.h5ad` | Expression matrix + metadata | Full analysis in Python |
| `*_metadata.csv` | Cell annotations only | Quick filtering, exploration |
| `*_summary.json` | Aggregated statistics | Quick overview, QC |
| `*_log.json` | Query parameters | Reproducibility, tracking |
