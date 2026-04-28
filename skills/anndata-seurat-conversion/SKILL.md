---
name: anndata-seurat-conversion
description: Bidirectional conversion between AnnData h5ad (Python/scanpy) and Seurat RDS (R) formats for single-cell RNA-seq data. Use when users request to convert between h5ad and RDS in either direction, preserving expression data, metadata, reductions, and layers.
---

# AnnData / Seurat Bidirectional Conversion

Convert single-cell RNA-seq data between AnnData h5ad (Python/scanpy) and Seurat RDS (R) formats with full preservation of:
- Expression data (X: normalized, layers: counts/raw_counts)
- Cell metadata (obs / meta.data)
- Gene metadata (var / feature metadata)
- Dimensional reductions (obsm: PCA, UMAP, t-SNE)
- Gene loadings (varm: PCA loadings)
- Additional layers

## When to Use

This skill applies when:
- User provides an `.h5ad` file and wants to convert to Seurat `.rds`
- User provides a Seurat `.rds`/`.Rds`/`.RDS` file and wants to convert to `.h5ad`
- User mentions converting between AnnData/scanpy and Seurat/R in either direction
- User needs to use R-based tools (CellChat, Seurat, etc.) with h5ad data
- User needs to use Python-based tools (scanpy, etc.) with Seurat data

## Direction Detection

Detect the conversion direction from the input file extension:
- `.h5ad` input → **h5ad to Seurat** (Section A)
- `.rds` / `.Rds` / `.RDS` input → **Seurat to h5ad** (Section B)

---

# Section A: H5AD to Seurat Conversion

## Workflow

### Step A1: Validate the Input File

```bash
python @scripts/validate_h5ad.py <input.h5ad>
```

Returns JSON with `valid`, `error`, and `structure` (cells, genes, layers, obsm keys).

**If validation fails:** Inform the user and request a valid h5ad file.

### Step A2: Run Conversion

```bash
conda run -n r_env Rscript @scripts/convert_h5ad_to_seurat.R <input.h5ad> [output.rds] [output.log]
```

**Arguments:**
- `input.h5ad` (required): Path to input h5ad file
- `output.rds` (optional): Defaults to same name with `.rds` extension
- `output.log` (optional): Defaults to same name with `.log` extension

### Step A3: Review Log and Summarize

Read the log file and summarize: status, cells/genes, counts source, layers, reductions, metadata columns, warnings.

## What Gets Converted (h5ad → Seurat)

| AnnData Component | Seurat Component | Notes |
|-------------------|------------------|-------|
| `X` | `RNA$data` | Log-normalized expression |
| `raw.X` or `layers['counts']` | `RNA$counts` | Integer counts (validated) |
| `obs` | `meta.data` | Cell metadata |
| `obsm['X_pca']` | `pca` reduction | PCA embeddings |
| `obsm['X_umap']` | `umap` reduction | UMAP embeddings |
| `obsm['X_tsne']` | `tsne` reduction | t-SNE embeddings |
| `layers` | `RNA` layers | Additional layers |

## Counts Source Selection (h5ad → Seurat)

1. **raw.X** - Preferred if integer/non-negative
2. **layers['raw_counts']** or **layers['counts']** - Fallback
3. **X** - Last resort (warns if normalized)

---

# Section B: Seurat to H5AD Conversion

## Workflow

### Step B1: Validate the Input File

```bash
conda run -n r_env Rscript @scripts/validate_rds.R <input.rds>
```

Returns JSON with `valid`, `class`, `n_cells`, assay structure (layers, features, var columns), reductions (dims, loadings), and metadata columns.

**If validation fails:** Inform the user and request a valid Seurat RDS file.

### Step B2: Run Conversion

**IMPORTANT:** The `RETICULATE_PYTHON` environment variable MUST be set to a Python with `anndata` installed before running. Ask the user which conda environment to use if unknown.

```bash
RETICULATE_PYTHON=/path/to/conda/envs/<env>/bin/python3 \
  conda run -n r_env Rscript @scripts/convert_seurat_to_h5ad.R <input.rds> [output.h5ad] [output.log] [assay]
```

**Arguments:**
- `input.rds` (required): Path to input Seurat RDS file
- `output.h5ad` (optional): Defaults to same basename with `.h5ad` in the **working directory**
- `output.log` (optional): Defaults to same basename with `.log`
- `assay` (optional): Assay to convert (default: `RNA`)

**Examples:**
```bash
# Basic conversion (output: ./data.h5ad in working dir)
RETICULATE_PYTHON=/home/user/.conda/envs/spatial/bin/python3 \
  conda run -n r_env Rscript @scripts/convert_seurat_to_h5ad.R /path/to/data.Rds

# Specify output and assay
RETICULATE_PYTHON=/home/user/.conda/envs/spatial/bin/python3 \
  conda run -n r_env Rscript @scripts/convert_seurat_to_h5ad.R data.Rds output.h5ad output.log SCT
```

### Step B3: Review Log and Summarize

Read the log file and summarize: status, cells/genes, X source, layers, obs/var columns, obsm keys, varm keys, warnings.

## What Gets Converted (Seurat → h5ad)

| Seurat Component | AnnData Component | Notes |
|------------------|-------------------|-------|
| `RNA$data` (log-normalized) | `X` | Cells x genes sparse matrix |
| `RNA$counts` (raw integers) | `layers["counts"]` | Preserved as-is |
| Additional assay layers | `layers[name]` | Excluding data/counts/scale.data |
| `meta.data` | `obs` | Factors → character, logicals → character |
| Feature metadata (`[[]]`) | `var` | From source assay + other assays |
| `Embeddings(pca)` | `obsm["X_pca"]` | All reductions mapped to X_ prefix |
| `Embeddings(umap)` | `obsm["X_umap"]` | |
| `Embeddings(tsne)` | `obsm["X_tsne"]` | |
| `Loadings(pca)` | `varm["PCs"]` | Zero-padded to full gene set |
| `Loadings(other)` | `varm["name_loadings"]` | Zero-padded to full gene set |

## Key Details (Seurat → h5ad)

- **Transposition**: Seurat (genes x cells) → AnnData (cells x genes) handled automatically
- **Factor columns**: Converted to character strings for h5ad compatibility
- **Logical columns**: Converted to character to avoid NA edge cases
- **varm zero-padding**: Only PCA-contributing genes have non-zero loadings; remaining genes padded with zeros for shape compatibility
- **var metadata**: Collected from source assay first, then fills in columns from other assays (matched by gene name, NA for missing)
- **scale.data layer**: Skipped by default (dense, large, and easily recomputed)

---

## Large Dataset Support (BPCells On-Disk)

For h5ad files with **>1M cells**, the standard reticulate-based conversion crashes (segfault during Python→R sparse matrix transfer when nnz exceeds R's 2.1B int32 limit). The script auto-detects this and switches to **BPCells on-disk mode**.

### How It Works

1. **Auto-detection**: After loading h5ad dimensions, if `n_cells > 1,000,000` and BPCells is installed, the BPCells path activates automatically
2. **BPCells reads h5ad directly** via HDF5 — no reticulate, no Python→R matrix transfer
3. **On-disk matrix**: Counts stored in a `bpcells_counts/` directory, not in RAM
4. **Metadata via hdf5r**: Cell metadata read directly from HDF5 categorical encoding, with Python fallback for edge cases

### Output (BPCells path)

| File | Description |
|------|-------------|
| `output.rds` | Seurat v5 object (~100-200MB, references BPCells dir) |
| `bpcells_counts/` | On-disk matrix directory (size ≈ compressed counts) |

**IMPORTANT**: The RDS and `bpcells_counts/` directory must stay together. Moving the RDS without the directory breaks the matrix reference.

### Downstream Compatibility

Some R/Bioconductor tools don't support BPCells IterableMatrix. Use helper functions from `bpcells_helpers.R`:

| Function | Replaces | Use When |
|----------|----------|----------|
| `aggregate_bpcells_per_cluster()` | `scran::aggregateAcrossCells()` | Pseudobulk aggregation |
| `bpcells_find_all_markers()` | `presto::wilcoxauc()` / `scran::findMarkers()` | Marker gene detection |
| `is_bpcells_matrix()` | — | Check if matrix is BPCells-backed |
| `estimate_h5ad_size()` | — | Pre-check h5ad size before loading |
| `read_h5ad_metadata()` | — | Read obs without loading expression data |
| `validate_bpcells_rds()` | — | Verify BPCells dir exists alongside RDS |

```r
source("path/to/bpcells_helpers.R")

# Per-cluster aggregation (materializes each cluster subset, fits in RAM)
sce_agg <- aggregate_bpcells_per_cluster(sce, "celltype", agg_cols)

# Marker detection (uses Seurat v5 FindAllMarkers, supports BPCells natively)
markers <- bpcells_find_all_markers(seurat_obj, "celltype")
```

### Thresholds

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Cell threshold | 1,000,000 | Reticulate handles <1M cells (~1B nnz) safely |
| nnz threshold | 1,500,000,000 | 70% of R's 2.1B int32 limit, safety margin |
| Per-cluster nnz limit | ~500,000,000 | Largest clusters (~500K cells) fit per-cluster |

---

## Error Handling (Both Directions)

| Error | Cause | Solution |
|-------|-------|----------|
| "File not found" | Invalid path | Verify file path exists |
| "Not a valid h5ad" | Corrupted or wrong format | Provide valid h5ad file |
| "Object is class X, not Seurat" | RDS does not contain Seurat object | Provide valid Seurat RDS |
| "Seurat v5+ required" | Old Seurat version | Upgrade Seurat package |
| "RETICULATE_PYTHON not set" | Missing Python path | Set env var to Python with anndata |
| "No 'data' or 'counts' layer" | Empty assay | Check assay name argument |

## Requirements

**R packages (in r_env conda environment):**
- Seurat (v5+)
- anndata (R)
- Matrix
- reticulate
- BPCells (for large datasets >1M cells)
- hdf5r (for large dataset metadata reading)

**Python packages (in RETICULATE_PYTHON environment):**
- anndata

**Python packages (for h5ad validation / metadata fallback):**
- anndata
