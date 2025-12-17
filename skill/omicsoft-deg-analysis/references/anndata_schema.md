# Omicsoft AnnData Object Schema

This document describes the structure of AnnData objects used in Omicsoft DEG analysis.

## DEG AnnData Object (`*_deg.h5ad`)

The primary object for differential expression analysis containing pre-computed statistics.

### Structure
```python
AnnData object with n_obs × n_vars = {n_comparisons} × {n_genes}
    obs: comparison metadata (study, tissue, disease, etc.)
    var: gene annotations (gene_name)
    uns: may contain embedded expression data (automatically removed)
    layers: 'log2fc', 'padj', 'sig_score'
    X: copy of sig_score layer (sparse CSR matrix)
```

### Key Components

**obs (Observations - Comparisons)**
Each row = one differential expression comparison (disease vs control in a study)
- Index: `sample_id` or `comparison_index`
- Core columns: `database`, `comparison_id`, `study`, `tissue`, `disease`, `disease_category`, `comparison`
- Case/Control metadata: `case_tissue`, `case_disease_state`, `control_disease_state`, etc.
- Optional columns vary by dataset (treatment, demographics, etc.)

**var (Variables - Genes)**
Each row = one gene
- Index: `gene_id`
- Column: `gene_name` (UPPERCASE for human, Title case for mouse)

**layers (Expression Data)**
All stored as sparse CSR matrices:
- `log2fc`: Log2 fold change (default=0)
- `padj`: Adjusted p-values (default=1)
- `sig_score`: Significance score = `log2fc * -log(padj)` where `padj < 0.05` (default=0)

**uns (Unstructured Annotations)**
- May contain `{query_name}_expr` key with embedded expression AnnData
- Example: `uns['ibd_mash_fibro_derm_rheum_10222025_expr']`
- **Automatically removed by `load_h5ad()`** to save memory

## Expression AnnData Object (`*_expr.h5ad`)

Companion object with raw expression values (optional, may be embedded in DEG object's `uns`).

### Structure
```python
AnnData object with n_obs × n_vars = {n_samples} × {n_genes}
    obs: sample metadata (disease_state, tissue, treatment, etc.)
    var: gene annotations (gene_name)
    layers: 'tpm', 'fpkm', 'raw_counts'
    X: TPM + microarray values (sparse CSR matrix)
```

### Key Components

**obs (Observations - Samples)**
Each row = one biological sample
- Core columns: `database`, `project_id`, `sample_id`, `tissue`, `disease_state`, `experiment_type`
- Demographics: `gender`, `age_summary`, `ethnicity`
- Treatment: `treatment`, `response`, `sampling_time`

**layers (Expression Values)**
- `tpm`: Transcripts per million (RNA-seq)
- `fpkm`: Fragments per kilobase million (RNA-seq)
- `raw_counts`: Raw read counts (RNA-seq)
- For microarray: all layers contain normalized values

## Automatic Cleanup

When loading `*_deg.h5ad` files, the `load_h5ad()` function:
1. Extracts the base filename (e.g., `ibd_mash_fibro_derm_rheum_10222025_deg`)
2. Constructs the potential expr key: `ibd_mash_fibro_derm_rheum_10222025_expr`
3. Checks if this key exists in `adata.uns`
4. If found, deletes it to reduce memory usage
5. Prints confirmation message

**Why remove expression data?**
- Expression data is 10-100x larger than DEG statistics
- Not needed for DEG analysis (only pre-computed log2fc/padj)
- Can always load `*_expr.h5ad` separately if needed
- Prevents memory issues

## Gene Symbol Conventions

**Critical**: Gene symbols must match organism:
- **Human**: ALL UPPERCASE (e.g., `CGAS`, `TMEM173`, `TBK1`)
- **Mouse**: Title case (e.g., `Cgas`, `Tmem173`, `Tbk1`)

See `signature_format.md` for detailed gene symbol guidelines.

## Exploring H5AD Schema

To discover available metadata values before running analysis, use the schema exploration tools:

```bash
# Generate comprehensive schema report
python scripts/explore_h5ad_schema.py --file data.h5ad --output schema.json --format json

# Create interactive HTML viewer
python scripts/generate_schema_viewer.py --json schema.json --output viewer.html
```

The schema exploration tools help identify:
- Available diseases, tissues, treatments for filtering
- Demographics (gender, ethnicity, age) when present
- Study identifiers and comparison types
- All unique values for every metadata column

See `schema_exploration.md` for complete documentation on using these tools.
