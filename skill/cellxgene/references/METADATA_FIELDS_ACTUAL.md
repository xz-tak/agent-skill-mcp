# Actual CELLxGENE Census Metadata Fields

Based on inspection of the CELLxGENE Census (stable release 2025-01-30).

---

## Cell Metadata (28 fields total)

These fields are included in the `{prefix}_metadata.csv` file for each cell.

### 1. Cell Identity & Annotations (6 fields)

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `cell_type` | category | Cell type annotation | endothelial cell, T cell, macrophage |
| `cell_type_ontology_term_id` | category | Cell Ontology ID | CL:0000115, CL:0001064 |
| `soma_joinid` | int64 | Unique cell identifier | 0, 1, 2, ... |
| `observation_joinid` | object | Join ID for linking | qsW0>t$X%X |
| `donor_id` | category | Unique donor/patient ID | HTAPP-330 |
| `dataset_id` | category | Dataset identifier | d7476ae2-e320-4703-8304-da5c42 |

### 2. Tissue & Anatomy (6 fields)

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `tissue` | category | Detailed tissue name | liver, lung parenchyma, left lung |
| `tissue_general` | category | General tissue category | liver, lung, intestine |
| `tissue_ontology_term_id` | category | UBERON tissue ID | UBERON:0002107 |
| `tissue_general_ontology_term_id` | category | UBERON general tissue ID | UBERON:0002107 |
| `tissue_type` | category | Tissue type classification | tissue |
| `suspension_type` | category | Sample preparation | cell, nucleus |

### 3. Donor Characteristics (8 fields)

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `sex` | category | Biological sex | female, male, unknown |
| `sex_ontology_term_id` | category | PATO sex ID | PATO:0000383 |
| `development_stage` | category | Developmental stage | 29-year-old stage, adult |
| `development_stage_ontology_term_id` | category | HsapDv stage ID | HsapDv:0000123 |
| `self_reported_ethnicity` | category | Self-reported ethnicity | European, Asian |
| `self_reported_ethnicity_ontology_term_id` | category | HANCESTRO ethnicity ID | HANCESTRO:0005 |
| `disease` | category | Disease annotation | breast cancer, COVID-19, normal |
| `disease_ontology_term_id` | category | MONDO/PATO disease ID | MONDO:0007254 |

### 4. Experimental Details (2 fields)

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `assay` | category | Experimental technique | 10x 3' v3, 10x 5' v2, Smart-seq2 |
| `assay_ontology_term_id` | category | EFO assay ID | EFO:0009922 |

### 5. Data Quality Metrics (6 fields)

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `is_primary_data` | bool | Primary vs derived data | True, False |
| `n_measured_vars` | int64 | Number of genes measured | 12641 |
| `nnz` | int64 | Number of non-zero gene counts | 7157, 5388 |
| `raw_sum` | float64 | Total UMI counts | 19641.0, 17251.0 |
| `raw_mean_nnz` | float64 | Mean expression of detected genes | 2.744306 |
| `raw_variance_nnz` | float64 | Variance of detected genes | 696.13164 |

---

## Gene Metadata (7 fields total)

These fields describe each gene in the expression matrix (stored in AnnData `.var`).

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `feature_id` | object | Ensembl gene ID | ENSG00000237491, ENSG00000188976 |
| `feature_name` | object | Gene symbol | LINC01409, NOC2L, PERM1 |
| `feature_type` | object | Gene biotype | lncRNA, protein_coding, pseudogene |
| `feature_length` | int64 | Gene length (bp) | 1059, 1244, 2765 |
| `soma_joinid` | int64 | Unique gene identifier | 0, 1, 2, ... |
| `n_measured_obs` | int64 | Number of cells with this gene | 92252850 |
| `nnz` | int64 | Number of non-zero observations | 7958785 |

---

## Summary Statistics JSON Structure

The `{prefix}_summary.json` file contains aggregated statistics:

```json
{
  "timestamp": "2025-12-05T10:30:45.123456",
  "filters_applied": "(tissue_general == \"lung\")",
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
    "macrophage": 35678,
    "endothelial cell": 20120
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

---

## Common Field Values

### Available Tissues (examples)
- `lung`, `liver`, `intestine`, `brain`, `blood`
- `adipose tissue`, `breast`, `kidney`, `heart`, `spleen`
- `bone marrow`, `lymph node`, `thymus`, `pancreas`

### Available Cell Types (examples)
- **Immune**: T cell, B cell, macrophage, monocyte, NK cell, dendritic cell
- **Epithelial**: epithelial cell, goblet cell, enterocyte, pneumocyte
- **Stromal**: fibroblast, endothelial cell, pericyte, smooth muscle cell
- **Other**: hepatocyte, neuron, adipocyte, erythrocyte

### Available Assays (examples)
- **10x Genomics**: 10x 3' v1/v2/v3, 10x 5' v1/v2, 10x gene expression flex
- **Other**: Smart-seq2, BD Rhapsody, CITE-seq, sci-RNA-seq

### Available Diseases (examples)
- **Normal**: normal (healthy tissue)
- **Cancer**: breast cancer, lung adenocarcinoma, glioblastoma
- **Infectious**: COVID-19, HIV, tuberculosis
- **Autoimmune**: Crohn disease, ulcerative colitis, rheumatoid arthritis
- **Other**: Alzheimer disease, diabetes, COPD

### Development Stages (examples)
- **Age-based**: 1-year-old stage, 29-year-old stage, 65-year-old stage
- **Life stage**: embryonic, fetal, postnatal, adult
- **Detailed**: 10th week post-fertilization stage

---

## Example Usage

### Load and Explore Metadata
```python
import pandas as pd

# Load metadata CSV
metadata = pd.read_csv("lung_tcells_metadata.csv")

# Show all available columns
print(metadata.columns.tolist())

# Basic statistics
print(f"Total cells: {len(metadata)}")
print(f"Unique donors: {metadata['donor_id'].nunique()}")
print(f"Unique datasets: {metadata['dataset_id'].nunique()}")

# Cell type distribution
print(metadata['cell_type'].value_counts())

# Filter by quality
high_quality = metadata[metadata['n_measured_vars'] > 10000]
print(f"High quality cells (>10k genes): {len(high_quality)}")
```

### Analyze by Donor
```python
# Cells per donor
donor_counts = metadata.groupby('donor_id').size()
print(donor_counts.describe())

# Disease status per donor
donor_disease = metadata.groupby('donor_id')['disease'].first()
print(donor_disease.value_counts())
```

### Quality Control
```python
# Check UMI counts
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 4))

plt.subplot(1, 2, 1)
plt.hist(metadata['raw_sum'], bins=50)
plt.xlabel('Total UMI counts')
plt.ylabel('Number of cells')
plt.title('UMI Count Distribution')

plt.subplot(1, 2, 2)
plt.hist(metadata['n_measured_vars'], bins=50)
plt.xlabel('Genes detected')
plt.ylabel('Number of cells')
plt.title('Gene Detection Distribution')

plt.tight_layout()
plt.savefig('qc_metrics.png')
```

### Filter by Multiple Criteria
```python
# Adult females with normal lung tissue, high quality
filtered = metadata[
    (metadata['tissue_general'] == 'lung') &
    (metadata['sex'] == 'female') &
    (metadata['disease'] == 'normal') &
    (metadata['development_stage'].str.contains('adult')) &
    (metadata['n_measured_vars'] > 10000) &
    (metadata['raw_sum'] > 5000)
]

print(f"Filtered cells: {len(filtered)}")
print(f"Cell types: {filtered['cell_type'].unique()}")
```

---

## Inspect Fields in Your Environment

Run this command to see current field values:

```bash
conda activate claude_test
cd /home/sagemaker-user/claude_code/gene-expression-specificity-cellxgene

# Inspect human data
conda run -n claude_test python scripts/inspect_metadata_fields.py --sample-size 1000

# Inspect mouse data
conda run -n claude_test python scripts/inspect_metadata_fields.py --organism mus_musculus --sample-size 1000
```

This will show:
- All available field names and types
- Example values for each field
- Value distributions for key categorical fields
- Current Census version being used

---

## Key Takeaways

1. **28 cell metadata fields** covering identity, tissue, donor characteristics, experimental details, and quality metrics

2. **7 gene metadata fields** with Ensembl IDs, gene symbols, biotypes, and statistics

3. **Quality metrics** (`n_measured_vars`, `nnz`, `raw_sum`) help filter high-quality cells

4. **Ontology IDs** link to standardized databases (Cell Ontology, UBERON, MONDO, etc.)

5. **Development stages** are specific (e.g., "29-year-old stage") not general (e.g., "adult")

6. **Summary JSON** provides quick overview without loading full dataset

7. **All fields** are available in both the CSV metadata and the AnnData object

---

## Additional Resources

- **Inspect script**: `scripts/inspect_metadata_fields.py`
- **Full documentation**: `METADATA_REFERENCE.md`
- **Query guide**: `scripts/QUERY_README.md`
- **Quick start**: `CELLXGENE_QUERY_QUICKSTART.md`
