---
name: anndata-seurat-conversion
description: Convert AnnData h5ad files to Seurat RDS objects. This skill should be used when users request to convert single-cell data from Python (scanpy/AnnData h5ad format) to R (Seurat RDS format), or when working with h5ad files that need to be analyzed in R/Seurat workflows.
---

# AnnData to Seurat Conversion

Convert single-cell RNA-seq data from AnnData h5ad format (Python/scanpy) to Seurat RDS format (R) with full preservation of:
- Expression data (counts and normalized)
- Cell metadata (obs)
- Dimensional reductions (PCA, UMAP, t-SNE)
- Layer data

## When to Use

This skill applies when:
- User provides an `.h5ad` file path and wants to convert to Seurat
- User mentions converting from AnnData/scanpy to Seurat/R
- User needs to use R-based tools (CellChat, Seurat, etc.) with h5ad data

## Workflow

### Step 1: Validate the Input File

Before conversion, validate the h5ad file structure:

```bash
python @scripts/validate_h5ad.py <input.h5ad>
```

This returns JSON with:
- `valid`: Whether file is a valid AnnData object
- `error`: Error message if invalid
- `structure`: Object structure (cells, genes, layers, obsm keys, etc.)

**If validation fails:** Inform the user that the file is not a valid AnnData h5ad object and request a proper single-cell h5ad file.

### Step 2: Run Conversion

Execute the R conversion script:

```bash
Rscript @scripts/convert_h5ad_to_seurat.R <input.h5ad> [output.rds] [output.log]
```

**Arguments:**
- `input.h5ad` (required): Path to input h5ad file
- `output.rds` (optional): Path to output RDS file. Defaults to same name with `.rds` extension
- `output.log` (optional): Path to log file. Defaults to same name with `.log` extension

**Examples:**
```bash
# Basic conversion (output: data.rds, data.log)
Rscript @scripts/convert_h5ad_to_seurat.R data.h5ad

# Specify output location
Rscript @scripts/convert_h5ad_to_seurat.R data.h5ad /output/converted.rds

# Specify all paths
Rscript @scripts/convert_h5ad_to_seurat.R data.h5ad /output/converted.rds /output/conversion.log
```

### Step 3: Review Log and Summarize

After conversion, read the log file and provide a summary to the user:

**Summary should include:**
1. **Status**: SUCCESS or FAILED
2. **Output files**: RDS file path and size, log file path
3. **Conversion details**:
   - Number of cells and genes
   - Counts source (raw.X, layer, or X fallback)
   - Layers preserved
   - Reductions added (PCA, UMAP, etc.)
   - Metadata columns count
4. **Warnings**: Any issues encountered during conversion

## What Gets Converted

| AnnData Component | Seurat Component | Notes |
|-------------------|------------------|-------|
| `X` | `RNA$data` | Log-normalized expression |
| `raw.X` or `layers['counts']` | `RNA$counts` | Integer counts (validated) |
| `obs` | `meta.data` | Cell metadata |
| `obsm['X_pca']` | `pca` reduction | PCA embeddings |
| `obsm['X_umap']` | `umap` reduction | UMAP embeddings |
| `obsm['X_tsne']` | `tsne` reduction | t-SNE embeddings |
| `layers` | `RNA` layers | Additional layers |

## Counts Source Selection

The script automatically selects the best counts source:
1. **raw.X** - Preferred if contains integer/non-negative values
2. **layers['raw_counts']** or **layers['counts']** - If raw.X fails validation
3. **X (fallback)** - Last resort, warns user if X appears normalized

## Skipped Components

The following are intentionally skipped:
- **Non-dimensional obsm**: Pathway scores, activity estimates (have named columns instead of dimension numbers)
- **Duplicate layers**: Layer used as counts is not added again as separate layer

## Error Handling

| Error | Cause | Solution |
|-------|-------|----------|
| "File not found" | Invalid path | Verify file path exists |
| "Not a valid h5ad" | Corrupted or wrong format | Provide valid AnnData h5ad file |
| "Seurat v5+ required" | Old Seurat version | Upgrade Seurat package |
| "C stack overflow" | Memory issue (avoided) | Script handles this automatically |

## Requirements

**R packages:**
- Seurat (v5+)
- anndata
- Matrix
- reticulate

**Python packages:**
- anndata (for validation script)
