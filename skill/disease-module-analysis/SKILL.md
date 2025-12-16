---
name: disease-module-analysis
description: This skill provides specialized workflows for disease module analysis using single-cell RNA-seq data and gene sets. Use this skill when users request disease module training/transformation, gene set enrichment analysis (pathway/ontology), perturbation enrichment (TF/Gene/LINCS), regulator identification, or NCBI metadata extraction for biological context. Supports both module-based and dataframe-based approaches for analyzing gene expression patterns, identifying regulatory relationships, and enriching results with biological annotations.
---

# Disease Module Analysis

## Overview

This skill enables comprehensive disease module analysis workflows for single-cell RNA-seq data and gene sets. Perform module training and transformation, gene set enrichment, perturbation enrichment, regulator discovery, and NCBI metadata enrichment. Supports both module-based analysis (with trained CorrModules objects) and dataframe-based analysis (working directly on gene sets without requiring modules).

## When to Use This Skill

Invoke this skill when users request:

- **Module training**: "Train disease modules from single-cell data"
- **Module transformation**: "Transform data with pretrained modules"
- **Gene queries**: "Which modules contain COL1A1 and ACTA2?"
- **Gene set enrichment**: "Perform pathway enrichment on these gene sets"
- **Perturbation enrichment**: "Find TF/Gene perturbations for these modules/gene sets"
- **Regulator discovery**: "Find regulators for TGFB1, COL1A1, ACTA2"
- **Metadata enrichment**: "Add tissue/cell type information from NCBI"

## Decision Tree: Module-Based vs Dataframe-Based

Choose the appropriate workflow based on the user's data and goals:

```
Do you have a trained module or single-cell data?
├─ YES: Use Module-Based Analysis
│   ├─ Need to train? → disease_module.py (train_and_transform)
│   ├─ Have trained module? → disease_module.py (transform_with_pretrained)
│   ├─ Query genes? → module_query.py (query mode)
│   ├─ Find perturbations? → module_query.py (perturb mode)
│   └─ Find regulators? → module_query.py (regulators mode)
│
└─ NO: Use Dataframe-Based Analysis
    ├─ Pathway/ontology enrichment? → geneset_enrichment.py
    ├─ Perturbation enrichment? → geneset_perturbation.py (perturb mode)
    └─ Find regulators? → geneset_perturbation.py (regulators mode)
```

**Key Difference:**
- **Module-based**: Requires trained module, provides module context, tracks module-restricted effects
- **Dataframe-based**: Works on any gene set table, faster, no module training needed

## Module-Based Analysis

Use when working with single-cell data or trained modules.

### 1. Train Disease Modules

Train modules from single-cell AnnData and transform data in one step.

**Command:**
```bash
python scripts/disease_module.py \
    --scanpy-h5ad data.h5ad \
    --genelist-module genes.txt \
    --resolution 5 \
    --n-neighbors 15 \
    --output results/
```

**Programmatic:**
```python
from scripts.disease_module import train_and_transform

module, adata_mod = train_and_transform(
    scanpy_h5ad_path='data.h5ad',
    genelist_module=['COL1A1', 'ACTA2', 'FN1'],  # Optional force-include genes
    output_dir='results/',
    use_hvg=True,
    resolution=5,
    n_neighbors=15,
)
```

**Parameters:**
- `scanpy_h5ad_path`: Path to H5AD file
- `genelist_module`: Optional list of genes to force include (useful for domain-specific genes)
- `resolution`: Leiden clustering resolution (higher = more modules)
- `n_neighbors`: KNN neighbors (affects module connectivity)
- `dotplot_groupby`: Obs columns for dotplot (default: ['cluster', 'condition'])
- `standard_scale`: Scaling method ('var', 'obs', or None; default: 'var')

**Outputs:**
- `module_YYYYMMDD_HHMMSS.pickle`: Trained module object
- `transformed_YYYYMMDD_HHMMSS.h5ad`: Transformed data with module scores
- `dotplot_SCALED_*.png` and `dotplot_UNSCALED_*.png`: Visualization
- `module_annotation_YYYYMMDD_HHMMSS.txt`: Enrichment results
- Log file with all parameters

**When to use genelist_module:**
- Include domain-specific marker genes
- Ensure key genes are not filtered out
- Force important pathways into analysis

### 2. Transform with Pretrained Module

Use existing trained module to transform new data without retraining.

**Command:**
```bash
python scripts/disease_module.py \
    --scanpy-h5ad new_data.h5ad \
    --pretrained-module module.pickle \
    --output results/
```

**Programmatic:**
```python
from scripts.disease_module import transform_with_pretrained

module, adata_mod = transform_with_pretrained(
    scanpy_h5ad_path='new_data.h5ad',
    module_path='module.pickle',
    output_dir='results/',
)
```

**Use cases:**
- Apply modules to validation datasets
- Compare module activity across conditions
- Batch processing multiple datasets with same modules

### 3. Query Genes in Modules

Find which modules contain specific genes.

**Command:**
```bash
python scripts/module_query.py \
    --module module.pickle \
    --query-genes genes.txt \
    --mode query \
    --output results/
```

**Programmatic:**
```python
from scripts.module_query import query_genes_in_modules

gene_module_df, genes_not_found = query_genes_in_modules(
    module_path='module.pickle',
    genelist=['COL1A1', 'ACTA2', 'FN1', 'TGFB1'],
    output_dir='results/',
)
```

**Outputs:**
- `moduleweights_YYYYMMDD.txt`: Gene-to-module mapping
- `genes_not_found_YYYYMMDD.txt`: Genes not in modules

**Important:** Always check genes_not_found to identify:
- Genes filtered during preprocessing
- Misspelled gene symbols
- Genes not expressed in dataset

### 4. Perturbation Enrichment (Module-Based)

Find TF/Gene/LINCS perturbations for specific modules.

**Command:**
```bash
python scripts/module_query.py \
    --module module.pickle \
    --modules 0 1 2 \
    --mode perturb \
    --species human \
    --enrich-metadata \
    --email user@example.com \
    --output results/
```

**Programmatic:**
```python
from scripts.module_query import perturbation_enrichment

results = perturbation_enrichment(
    module_path='module.pickle',
    module_list=['0', '1', '2'],
    output_dir='results/',
    species='human',
    enrich_metadata=True,
    email='user@example.com',
)

# Access results
tf_df = results['perturb_tf']
gene_df = results['perturb_gene']
crispr_df = results['perturb_crispr']
chem_df = results['perturb_chem']
```

**Outputs (per perturbation type):**
- `perturb_{type}_YYYYMMDD.txt`: Main results with p-values, overlap genes, Module_restricted flag
- `perturb_{type}_modules_YYYYMMDD.txt`: Module membership of perturbed genes
- `ppi_network_YYYYMMDD.pdf`: PPI network visualization
- Optional metadata columns: geo_id, species, tissue, cell_type, disease

**Key Column: Module_restricted**
- `True`: Perturbed gene is IN the same module it regulates (in-module effect)
- `False`: Perturbed gene is OUTSIDE the module it regulates (trans effect)

**Interpretation:**
- Module_restricted=True → likely direct/local regulation
- Module_restricted=False → likely upstream/trans regulation

**When to use enrich_metadata:**
- Need biological context (tissue, cell type, disease)
- Filter perturbations by experimental system
- Prioritize disease-relevant perturbations

### 5. Find Regulators (Module-Based)

Identify regulators for specific target genes using module context.

**Command:**
```bash
python scripts/module_query.py \
    --module module.pickle \
    --target-genes targets.txt \
    --modules 0 1 2 \
    --mode regulators \
    --enrich-metadata \
    --email user@example.com \
    --output results/
```

**Programmatic:**
```python
from scripts.module_query import find_regulators_for_targets

regulators = find_regulators_for_targets(
    module_path='module.pickle',
    target_genes=['TGFB1', 'COL1A1', 'ACTA2'],
    module_list=['0', '1', '2'],  # Optional: restrict to specific modules
    species='human',
    enrich_metadata=True,
    email='user@example.com',
)

# Access regulators
tf_regulators = regulators['tf']
gene_regulators = regulators['gene']
```

**Outputs:**
- `regulators_{type}_YYYYMMDD.txt`: Regulators affecting target genes
- `regulators_{type}_modules_YYYYMMDD.txt`: Module membership

**Use cases:**
- Drug target discovery: Find upstream regulators
- Mechanism investigation: Identify key control points
- Pathway analysis: Map regulatory networks

## Dataframe-Based Analysis

Use when working with gene sets without trained modules (faster, more flexible).

### 1. Gene Set Enrichment

Perform pathway/ontology enrichment on gene sets using gprofiler, enrichr, or x2k.

**Command:**
```bash
python scripts/geneset_enrichment.py \
    --input genesets.csv \
    --set-column gene_list \
    --method gprofiler \
    --species human \
    --pval 0.05 \
    --output results/
```

**Programmatic:**
```python
from scripts.geneset_enrichment import geneset_annotation
import pandas as pd

df = pd.DataFrame({
    'set_name': ['Fibrosis', 'Inflammation'],
    'gene_list': [
        ['COL1A1', 'ACTA2', 'FN1', 'TGFB1', 'CTGF'],
        ['IL6', 'TNF', 'IL1B', 'CCL2', 'CXCL8'],
    ],
})

df_enriched = geneset_annotation(
    df=df,
    set_column='gene_list',
    species='human',
    method='gprofiler',
    pval=0.05,
    output_dir='results/',
)
```

**Input Requirements:**
- DataFrame with a column containing lists of gene symbols
- Gene symbols as Python lists (not strings)
- Minimum 3-5 genes per set recommended

**Methods:**
- `gprofiler` (default): Most comprehensive, includes GO, KEGG, Reactome, WikiPathways
- `enrichr`: Enrichr database collection
- `x2k`: Expression2Kinases upstream regulator analysis

**Outputs:**
- `moduleweights_YYYYMMDD.txt`: Original data + enrichment columns
- Columns: `{method}_{term}` for each significant pathway/term
- Column: `PPI` for PPI enrichment p-value

**Enrichment Column Interpretation:**
- Each column = one significant pathway/term
- Value = enrichment score or p-value
- NaN = not significantly enriched for that term

### 2. Perturbation Enrichment (Dataframe-Based)

Find TF/Gene/LINCS perturbations for gene sets without requiring a module.

**Command:**
```bash
python scripts/geneset_perturbation.py \
    --input genesets.csv \
    --set-column gene_list \
    --mode perturb \
    --pert-types tf gene LINCS_CRISPR \
    --species human \
    --enrich-metadata \
    --email user@example.com \
    --output results/
```

**Programmatic:**
```python
from scripts.geneset_perturbation import perturbation_enrichment_df

df = pd.DataFrame({
    'pathway': ['TGF-beta', 'ECM'],
    'genes': [
        ['TGFB1', 'SMAD2', 'SMAD3', 'SMAD4'],
        ['COL1A1', 'FN1', 'SPARC', 'MMP2'],
    ],
})

results = perturbation_enrichment_df(
    df=df,
    set_column='genes',
    pert_types=['tf', 'gene', 'LINCS_CRISPR', 'LINCS_CHEM'],
    species='human',
    enrich_metadata=True,
    email='user@example.com',
    output_dir='results/',
)

# Access results by type
tf_perturbations = results['tf']
gene_perturbations = results['gene']
```

**Perturbation Types:**
- `tf`: Transcription factor perturbations
- `gene`: Gene knockdown/overexpression
- `LINCS_CRISPR`: LINCS CRISPR knockout screens
- `LINCS_CHEM`: LINCS chemical perturbations

**Outputs:**
- `perturbation_{type}_YYYYMMDD.txt`: Results with p-values and overlap genes
- Column `geneset_id`: Links back to input DataFrame index
- Optional metadata: geo_id, species, tissue, cell_type, disease

**Key Difference from Module-Based:**
- No Module_restricted column (no module context)
- Results include `geneset_id` to track which input gene set produced each result
- Faster: No need to load/query module structure

### 3. Find Regulators (Dataframe-Based)

Identify regulators for target genes using gene sets (no module required).

**Command:**
```bash
python scripts/geneset_perturbation.py \
    --input genesets.csv \
    --set-column gene_list \
    --mode regulators \
    --target-genes targets.txt \
    --pert-types tf gene \
    --enrich-metadata \
    --email user@example.com \
    --output results/
```

**Programmatic:**
```python
from scripts.geneset_perturbation import find_regulators_df

regulators = find_regulators_df(
    df=df,
    set_column='genes',
    target_genes=['TGFB1', 'COL1A1', 'ACTA2'],
    pert_types=['tf', 'gene'],
    species='human',
    enrich_metadata=True,
    email='user@example.com',
)
```

**Outputs:**
- `regulators_{type}_YYYYMMDD.txt`: Regulators affecting target genes

**Use cases:**
- Quick regulator discovery without module training
- Exploratory analysis of regulatory relationships
- Validate regulators across multiple gene sets

## NCBI Metadata Enrichment

Add biological context (species, tissue, cell type, disease) by querying NCBI GEO databases.

### Automatic Enrichment

Use `--enrich-metadata` flag in perturbation enrichment workflows:

```bash
# Module-based
python scripts/module_query.py ... --enrich-metadata --email user@example.com

# Dataframe-based
python scripts/geneset_perturbation.py ... --enrich-metadata --email user@example.com
```

### Manual Enrichment

For custom workflows:

```python
from scripts.ncbi_metadata import enrich_perturbation_tables

# Assume perturb_df has perturbation results
perturb_enriched = enrich_perturbation_tables(
    df=perturb_df,
    pert_type='gene',  # or 'tf', 'crispr', 'chem'
    email='user@example.com'
)

# Now has columns: geo_id, species, tissue, cell_type, disease
```

### Metadata Fields

**Extracted from NCBI GEO:**
- `geo_id`: GEO series ID (e.g., GSE12345)
- `species`: Human, Mouse, Rat
- `tissue`: Liver, Heart, Brain, Lung, Kidney (20+ tissues supported)
- `cell_type`: Fibroblast, Macrophage, Hepatocyte (20+ types supported)
- `disease`: Cancer, Fibrosis, Inflammation, Diabetes (30+ diseases supported)

**Multiple values:** Semicolon-separated (e.g., "Liver; Kidney")
**Missing data:** None or NaN

### NCBI API Guidelines

**Required:**
- Valid email address (NCBI requirement)
- Maximum 3 requests/second (automatically enforced with delay=0.4s)

**Best Practices:**
- Run during off-peak hours for large analyses
- Don't run multiple instances simultaneously
- Spot-check metadata accuracy
- Cache is session-specific (not persisted across runs)

**Processing Time:**
- 50 perturbations: ~20 seconds
- 100 perturbations: ~40 seconds
- 500 perturbations: ~3 minutes

### Use Cases for Metadata

**1. Filter by experimental system:**
```python
# Filter for liver studies
liver_df = results[results['tissue'].str.contains('Liver', na=False)]

# Filter for human fibroblasts
human_fib = results[
    (results['species'] == 'Human') &
    (results['cell_type'].str.contains('Fibroblast', na=False))
]
```

**2. Prioritize disease-relevant perturbations:**
```python
# Fibrosis-related perturbations
fibrosis = results[results['disease'].str.contains('Fibrosis', na=False)]
```

**3. Summarize by tissue:**
```python
# Count perturbations by tissue
tissue_counts = results['tissue'].value_counts()
```

## Common Parameters

### Species
- `human` (default): Homo sapiens
- `mouse`: Mus musculus
- `rat`: Rattus norvegicus

### Perturbation Types
- `tf`: Transcription factor perturbations
- `gene`: Gene knockdown/overexpression
- `LINCS_CRISPR`: LINCS CRISPR knockout screens
- `LINCS_CHEM`: LINCS chemical perturbations

### Enrichment Methods (geneset_enrichment.py only)
- `gprofiler`: GO, KEGG, Reactome, WikiPathways (most comprehensive)
- `enrichr`: Enrichr database collection
- `x2k`: Expression2Kinases upstream regulator analysis

### Standard Scale (disease_module.py only)
- `var` (default): Scale by variable (gene)
- `obs`: Scale by observation (cell/sample)
- `None`: No scaling

## Output Files

All scripts generate timestamped output files with consistent naming:

**Main results:** `{prefix}_YYYYMMDD.txt` (tab-separated)
- Example: `perturb_tf_20231215.txt`
- Example: `moduleweights_20231215.txt`

**Module objects:** `module_YYYYMMDD_HHMMSS.pickle`
**Transformed data:** `transformed_YYYYMMDD_HHMMSS.h5ad`
**Plots:** `dotplot_*.png`, `ppi_network_*.pdf`
**Logs:** `{script}_YYYYMMDD_HHMMSS.log`

**Log files contain:**
- All input parameters
- Progress updates
- Warnings and errors
- Summary statistics

## Best Practices

### Module Training
1. **Use HVG selection** for large datasets (>20k cells)
2. **Adjust resolution** based on expected biology (higher = more specific modules)
3. **Include domain genes** with genelist_module for key pathways
4. **Check gene filtering** by querying expected genes

### Gene Set Enrichment
1. **Minimum 3-5 genes** per set for meaningful enrichment
2. **Use official gene symbols** (uppercase: COL1A1 not col1a1)
3. **Start with gprofiler** (most comprehensive)
4. **Check PPI enrichment** for network context

### Perturbation Analysis
1. **Test multiple perturbation types** (TF, Gene, LINCS)
2. **Use enrich_metadata** for biological context
3. **Filter by significance** (p_value < 0.05 or 0.01)
4. **Check Module_restricted** (module-based only) to distinguish local vs trans effects

### NCBI Metadata
1. **Always provide valid email** (NCBI requirement)
2. **Don't modify delay** (respect NCBI limits)
3. **Spot-check results** for accuracy
4. **Process large datasets overnight**

### File Management
1. **Use dedicated output directories** to organize results
2. **Keep log files** for reproducibility
3. **Date-based naming** allows version tracking
4. **Tab-separated .txt** for easy import into Excel/R

## Troubleshooting

### Common Errors

**"Column not found"**
- Check `set_column` matches actual column name
- Use `df.columns` to list available columns

**"Column must contain lists"**
- Ensure gene sets are Python lists, not strings
- Use `df['col'].apply(eval)` if needed to convert strings to lists

**"No genes found in module"**
- Genes may be filtered during preprocessing
- Check gene symbol case (should be uppercase)
- Verify genes are expressed in dataset

**"No perturbations found"**
- Gene sets may be too small (<3 genes)
- Try different perturbation types
- Check species matches data

**NCBI API errors**
- Rate limit exceeded: Wait and retry
- Invalid email: Provide valid email address
- Network issues: Check internet connection

### Debugging Steps

1. **Check log files** for detailed error messages
2. **Verify input format** (gene lists as Python lists)
3. **Test with small example** before full analysis
4. **Check dependencies** are installed (pandas, scanpy, requests)

## Resources

### scripts/

All analysis scripts with both command-line and programmatic interfaces:

- `disease_module.py`: Module training and transformation
- `module_query.py`: Module querying and perturbation enrichment (module-based)
- `geneset_enrichment.py`: Gene set enrichment analysis (dataframe-based)
- `geneset_perturbation.py`: Perturbation enrichment (dataframe-based)
- `ncbi_metadata.py`: NCBI GEO metadata extraction

**Note:** Scripts can be executed without loading into context, but may need to be read for environment-specific adjustments.

### references/

Detailed API documentation with comprehensive function signatures, parameters, return values, and examples:

- `api_reference.md`: Complete API documentation for all functions

**When to load:** Reference when users need detailed parameter information, want to understand function signatures, or need programmatic usage examples.

**Search patterns for large files:**
```bash
# Find specific function documentation
grep -n "def train_and_transform" references/api_reference.md
grep -n "perturbation_enrichment" references/api_reference.md
```
