# Disease Module Analysis API Reference

Comprehensive API documentation for disease module analysis workflows.

## Table of Contents

1. [Module-Based Analysis](#module-based-analysis)
2. [Dataframe-Based Analysis](#dataframe-based-analysis)
3. [NCBI Metadata Enrichment](#ncbi-metadata-enrichment)
4. [Common Parameters](#common-parameters)
5. [Output Files](#output-files)

---

## Module-Based Analysis

Functions that require a trained module object.

### disease_module.py

#### `train_and_transform()`

Train disease modules and transform data in one step.

**Parameters:**
- `scanpy_h5ad_path` (str, required): Path to single-cell AnnData H5AD file
- `genelist_module` (list, optional): List of gene symbols to force include in training. Default: None
- `output_dir` (str, optional): Output directory. Default: same directory as H5AD file
- `use_hvg` (bool, optional): Use highly variable genes. Default: True
- `resolution` (int, optional): Leiden clustering resolution. Default: 5
- `n_neighbors` (int, optional): Number of neighbors for KNN. Default: 15
- `dotplot_groupby` (list, optional): List of obs columns for dotplot grouping. Default: ['cluster', 'condition'] or columns containing 'condition'
- `standard_scale` (str, optional): Standard scaling for dotplots ('var', 'obs', or None). Default: 'var'

**Returns:**
- `module`: Trained CorrModules object
- `adata_mod`: Transformed AnnData with module scores

**Output Files:**
- `module_YYYYMMDD_HHMMSS.pickle`: Trained module
- `transformed_YYYYMMDD_HHMMSS.h5ad`: Transformed data
- `dotplot_SCALED_{groupby}_YYYYMMDD_HHMMSS.png`: Scaled dotplots
- `dotplot_UNSCALED_{groupby}_YYYYMMDD_HHMMSS.png`: Unscaled dotplots
- `module_annotation_YYYYMMDD_HHMMSS.txt`: Enrichment results
- `disease_module_YYYYMMDD_HHMMSS.log`: Log file

**Command Line:**
```bash
python disease_module.py \
    --scanpy-h5ad data.h5ad \
    --genelist-module genes.txt \
    --resolution 5 \
    --n-neighbors 15 \
    --output results/
```

#### `transform_with_pretrained()`

Transform data using a pretrained module without training.

**Parameters:**
- `scanpy_h5ad_path` (str, required): Path to single-cell AnnData H5AD file
- `module_path` (str, required): Path to trained module pickle file
- `output_dir` (str, optional): Output directory. Default: same directory as H5AD file
- `dotplot_groupby` (list, optional): List of obs columns for dotplot grouping
- `standard_scale` (str, optional): Standard scaling for dotplots. Default: 'var'

**Returns:**
- `module`: Loaded CorrModules object
- `adata_mod`: Transformed AnnData with module scores

**Output Files:**
- `transformed_YYYYMMDD_HHMMSS.h5ad`: Transformed data
- `dotplot_SCALED_{groupby}_YYYYMMDD_HHMMSS.png`: Scaled dotplots
- `dotplot_UNSCALED_{groupby}_YYYYMMDD_HHMMSS.png`: Unscaled dotplots
- `module_annotation_YYYYMMDD_HHMMSS.txt`: Enrichment results
- `disease_module_YYYYMMDD_HHMMSS.log`: Log file

**Command Line:**
```bash
python disease_module.py \
    --scanpy-h5ad data.h5ad \
    --pretrained-module module.pickle \
    --output results/
```

### module_query.py

Functions for querying and analyzing trained modules.

#### `query_genes_in_modules()`

Query which modules contain specified genes.

**Parameters:**
- `module_path` (str, required): Path to trained module pickle file
- `genelist` (list, required): List of gene symbols to query
- `output_dir` (str, optional): Output directory. Default: same as module

**Returns:**
- `gene_module_df`: DataFrame with gene-module mappings
- `genes_not_found`: List of genes not found in modules

**Output Files:**
- `moduleweights_YYYYMMDD.txt`: Gene-module mapping (tab-separated)
- `genes_not_found_YYYYMMDD.txt`: Genes not found in modules
- `module_query_YYYYMMDD_HHMMSS.log`: Log file

**Command Line:**
```bash
python module_query.py \
    --module module.pickle \
    --query-genes genes.txt \
    --mode query \
    --output results/
```

#### `perturbation_enrichment()`

Perform perturbation enrichment for specified modules.

**Parameters:**
- `module_path` (str, required): Path to trained module pickle file
- `module_list` (list, required): List of module IDs to analyze
- `output_dir` (str, optional): Output directory. Default: same as module
- `species` (str, optional): Species for enrichment ('human', 'mouse', 'rat'). Default: 'human'
- `enrich_metadata` (bool, optional): Add NCBI metadata columns. Default: False
- `email` (str, optional): Email for NCBI API (required if enrich_metadata=True)

**Returns:**
- `results` (dict): Dictionary containing perturbation results for different types

**Output Files (per perturbation type):**
- `perturb_{type}_YYYYMMDD.txt`: Perturbation results (tab-separated)
- `perturb_{type}_modules_YYYYMMDD.txt`: Module mapping (tab-separated)
- `ppi_network_YYYYMMDD.pdf`: PPI network plot

**Command Line:**
```bash
python module_query.py \
    --module module.pickle \
    --modules 0 1 2 \
    --mode perturb \
    --species human \
    --enrich-metadata \
    --email user@example.com \
    --output results/
```

#### `find_regulators_for_targets()`

Find potential regulators for target genes.

**Parameters:**
- `module_path` (str, required): Path to trained module pickle file
- `target_genes` (list, required): List of target gene symbols
- `module_list` (list, optional): List of module IDs to restrict analysis
- `output_dir` (str, optional): Output directory
- `species` (str, optional): Species for enrichment
- `enrich_metadata` (bool, optional): Add NCBI metadata columns
- `email` (str, optional): Email for NCBI API

**Returns:**
- `regulators` (dict): Dictionary containing regulator results

**Output Files:**
- `regulators_{type}_YYYYMMDD.txt`: Regulator results (tab-separated)
- `regulators_{type}_modules_YYYYMMDD.txt`: Module mapping

**Command Line:**
```bash
python module_query.py \
    --module module.pickle \
    --target-genes targets.txt \
    --mode regulators \
    --enrich-metadata \
    --email user@example.com
```

---

## Dataframe-Based Analysis

Functions that work directly on dataframes without requiring a trained module.

### geneset_enrichment.py

#### `geneset_annotation()`

Perform enrichment analysis on gene sets in a dataframe.

**Parameters:**
- `df` (DataFrame, required): Input dataframe containing gene sets
- `set_column` (str, required): Column name containing gene sets (lists)
- `species` (str, optional): Species for enrichment. Default: 'human'
- `method` (str, optional): Enrichment method ('gprofiler', 'enrichr', 'x2k'). Default: 'gprofiler'
- `pval` (float, optional): P-value threshold. Default: 0.05
- `output_dir` (str, optional): Output directory

**Returns:**
- `df_annotated`: DataFrame with added enrichment columns

**Output Files:**
- `moduleweights_YYYYMMDD.txt`: Annotated results (tab-separated)

**Programmatic Usage:**
```python
from geneset_enrichment import geneset_annotation
import pandas as pd

df = pd.DataFrame({
    'set_name': ['Fibrosis', 'Inflammation'],
    'gene_list': [
        ['COL1A1', 'ACTA2', 'FN1', 'TGFB1'],
        ['IL6', 'TNF', 'IL1B', 'CCL2'],
    ],
})

df_enriched = geneset_annotation(
    df=df,
    set_column='gene_list',
    species='human',
    method='gprofiler',
    pval=0.05,
)
```

### geneset_perturbation.py

#### `perturbation_enrichment_df()`

Perform perturbation enrichment on gene sets without a module.

**Parameters:**
- `df` (DataFrame, required): Input dataframe containing gene sets
- `set_column` (str, required): Column name containing gene sets
- `pert_types` (list, optional): Perturbation types. Default: ['tf', 'gene']
- `species` (str, optional): Species. Default: 'human'
- `enrich_metadata` (bool, optional): Add NCBI metadata. Default: False
- `email` (str, optional): Email for NCBI API

**Returns:**
- `results` (dict): Dictionary with perturbation DataFrames

**Output Files:**
- `perturbation_{type}_YYYYMMDD.txt`: Results (tab-separated)

**Programmatic Usage:**
```python
from geneset_perturbation import perturbation_enrichment_df

results = perturbation_enrichment_df(
    df=df,
    set_column='genes',
    pert_types=['tf', 'gene'],
    species='human',
    enrich_metadata=True,
    email='user@example.com',
)
```

#### `find_regulators_df()`

Find regulators for target genes using gene sets.

**Parameters:**
- `df` (DataFrame, required): Input dataframe containing gene sets
- `set_column` (str, required): Column name containing gene sets
- `target_genes` (list, required): Target gene symbols
- `pert_types` (list, optional): Perturbation types
- `species` (str, optional): Species
- `enrich_metadata` (bool, optional): Add NCBI metadata
- `email` (str, optional): Email for NCBI API

**Returns:**
- `regulators` (dict): Dictionary with regulator DataFrames

**Output Files:**
- `regulators_{type}_YYYYMMDD.txt`: Results (tab-separated)

---

## NCBI Metadata Enrichment

### ncbi_metadata.py

#### `NCBIMetadataExtractor`

Extract metadata from NCBI GEO datasets.

**Constructor:**
```python
extractor = NCBIMetadataExtractor(email='user@example.com', delay=0.4)
```

**Methods:**
- `extract_geo_id(pert_name)`: Extract GEO ID from perturbation name
- `query_geo_metadata(geo_id)`: Query NCBI for dataset metadata
- `extract_metadata_fields(metadata)`: Extract structured fields

**Extracted Fields:**
- Species: Human, Mouse, Rat
- Tissues: 20+ supported (Liver, Heart, Brain, etc.)
- Cell Types: 20+ supported (Fibroblast, Macrophage, etc.)
- Diseases: 30+ supported (Cancer, Fibrosis, etc.)

#### `enrich_perturbation_tables()`

Add metadata columns to perturbation DataFrames.

**Parameters:**
- `df` (DataFrame, required): Perturbation results
- `pert_type` (str, required): Perturbation type ('tf', 'gene', 'crispr', 'chem')
- `email` (str, required): Email for NCBI API

**Returns:**
- `df_enriched`: DataFrame with metadata columns

**Added Columns:**
- `geo_id`: GEO series ID
- `species`: Species
- `tissue`: Tissue/organ
- `cell_type`: Cell type
- `disease`: Disease/condition

---

## Common Parameters

### Species
- `'human'` (default): Homo sapiens
- `'mouse'`: Mus musculus
- `'rat'`: Rattus norvegicus

### Perturbation Types
- `'tf'`: Transcription factor perturbations
- `'gene'`: Gene perturbations
- `'LINCS_CRISPR'`: LINCS CRISPR knockouts
- `'LINCS_CHEM'`: LINCS chemical perturbations

### Enrichment Methods
- `'gprofiler'` (default): Pathway and ontology enrichment
- `'enrichr'`: Enrichr database enrichment
- `'x2k'`: Expression2Kinases upstream analysis

---

## Output Files

### File Naming Convention
- Format: `{prefix}_YYYYMMDD.txt` for main results
- Format: `{prefix}_YYYYMMDD_HHMMSS.{ext}` for logs and pickles
- Example: `perturb_tf_20231215.txt`

### File Formats
- **Tab-separated (.txt)**: All result tables
- **H5AD (.h5ad)**: Transformed single-cell data
- **Pickle (.pickle)**: Trained module objects
- **PNG (.png)**: Dotplot figures
- **PDF (.pdf)**: PPI network plots

---

## Performance Considerations

### NCBI API Rate Limits
- Limit: 3 requests per second
- Default delay: 0.4 seconds (2.5 requests/sec)

### Processing Time Estimates
- Module training: ~5-10 minutes
- Perturbation enrichment (with metadata):
  - 50 perturbations: ~20 seconds
  - 100 perturbations: ~40 seconds
  - 500 perturbations: ~3 minutes

---

## Dependencies

Required packages:
- Python 3.7+
- xzsc_module
- pandas
- scanpy (for disease_module.py)
- requests (for ncbi_metadata.py)
