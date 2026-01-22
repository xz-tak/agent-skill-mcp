# Parameter Guide

## Pipeline Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input` | Path to Seurat RDS file | `/path/to/file.rds` |
| `--condition` | Column for condition/disease status | `sample_disease` |
| `--sample` | Column for sample identifier | `orig.ident` |
| `--celltype` | Column for cell type annotation | `scArches_Cell_Annotation` |
| `--comparisons` | Comparison pairs (semicolon-separated) | `"Healthy,Disease1;Healthy,Disease2"` |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output` | `./pseudobulk_results` | Output directory |
| `--batch` | `NULL` | Batch column for correction |
| `--covariate` | `NULL` | Additional covariate column |
| `--custom_signatures` | `NULL` | Path to custom gene signatures |
| `--min_count` | `3` | filterByExpr min.count |
| `--min_total_count` | `5` | filterByExpr min.total.count |
| `--min_prop` | `0.2` | filterByExpr min.prop |

---

## filterByExpr Parameters

The `pseudoBulkDGE` function uses edgeR's `filterByExpr` to remove lowly expressed genes before testing. Understanding these parameters is crucial for your analysis.

### `min.count` (default: 10, our override: 3)

**What it does:** A gene must have CPM (counts per million) ≥ `min.count` in at least `n` samples.

**The formula:**
```
CPM threshold = min.count / median_library_size * 1,000,000
```

**Effect of lowering:**
- Lower values (e.g., 3) → More genes tested → More multiple testing correction
- Higher values (e.g., 10) → Fewer genes tested → More stringent filtering

**When to lower:**
- When you have low sequencing depth
- When studying rare transcripts
- When cell types have low cell counts

### `min.total.count` (default: 15, our override: 5)

**What it does:** A gene must have ≥ `min.total.count` total reads across ALL samples.

**Effect:**
- This is a hard minimum - genes below this are always filtered
- Acts as a safety net for genes with sporadic expression

**When to adjust:**
- Lower for low-depth samples
- Higher if you want very strict filtering

### `min.prop` (default: 0.7, our override: 0.2)

**What it does:** Defines what fraction of the smallest group determines `n` (the number of samples where a gene must be expressed).

**The formula:**
```
n = min.prop × (size of smallest condition group)
```

**Example:**
- If smallest group has 10 samples and min.prop = 0.7: n = 7 samples
- If smallest group has 10 samples and min.prop = 0.2: n = 2 samples

**Effect of lowering:**
- Lower values (e.g., 0.2) → Gene only needs expression in 20% of smallest group
- Higher values (e.g., 0.7) → Gene needs expression in 70% of smallest group

**When to lower:**
- When cell types are highly variable
- When some samples may be missing certain genes
- When you expect heterogeneous expression patterns

### Recommended Settings by Scenario

| Scenario | min.count | min.total.count | min.prop |
|----------|-----------|-----------------|----------|
| Default (stringent) | 10 | 15 | 0.7 |
| Our recommendation | 3 | 5 | 0.2 |
| Very permissive | 1 | 3 | 0.1 |
| High-depth data | 15 | 25 | 0.7 |

---

## Custom Signatures

### Supported Formats

#### 1. GMT File (.gmt)
Standard Gene Matrix Transposed format used by MSigDB.

```
Signature_Name<TAB>Description<TAB>Gene1<TAB>Gene2<TAB>Gene3...
```

Example:
```
FIBROSIS_SIGNATURE	Fibrosis-related genes	COL1A1	COL3A1	FN1	ACTA2
INFLAMMATION_SIGNATURE	Inflammatory genes	IL6	TNF	IL1B	CCL2
```

#### 2. Tab-separated File
```
Signature_Name<TAB>Gene1,Gene2,Gene3
```

Example:
```
FIBROSIS_SIGNATURE	COL1A1,COL3A1,FN1,ACTA2
INFLAMMATION_SIGNATURE	IL6,TNF,IL1B,CCL2
```

#### 3. Simple Gene List (one gene per row)
```
Gene1
Gene2
Gene3
```

This creates a single signature called "custom_signature".

#### 4. Comma-separated Gene String
```
COL1A1,COL3A1,FN1,ACTA2,IL6,TNF
```

Pass directly as the `--custom_signatures` argument value.

### Important Notes

- **Custom signatures use GENE SYMBOLS** (same as MSigDB)
- They are analyzed with the same ranked gene list as MSigDB pathways
- Leading edge genes in results are also gene symbols
- Results saved to `fgsea_custom_signatures.tsv`

---

## Report Generation Parameters

These parameters are for the `generate_report.py` script:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fdr_cutoff` | 0.05 | FDR threshold for significance |
| `--logfc_cutoff` | 0.5 | Absolute log2FC threshold |
| `--nes_cutoff` | 0 | Absolute NES threshold (0 = any) |
| `--top_k` | 30 | Number of top results per comparison |

### Interpreting Results

**FDR (False Discovery Rate):**
- < 0.05: Conventionally significant
- < 0.01: Highly significant
- < 0.001: Very highly significant

**log2FC (Log2 Fold Change):**
- > 0: Upregulated in group2 vs group1
- < 0: Downregulated in group2 vs group1
- |log2FC| > 1: 2-fold change
- |log2FC| > 2: 4-fold change

**NES (Normalized Enrichment Score):**
- > 0: Pathway enriched in upregulated genes
- < 0: Pathway enriched in downregulated genes
- |NES| > 1.5: Typically considered meaningful

---

## Column Naming Conventions

The pipeline uses two versions of column names:

### Original Names (.obj versions)
Used for **data access** via `colData(sce)[[column]]`:
```r
condition.obj <- 'sample_disease'
celltype.obj <- 'scArches_Cell_Annotation'
```

### Safe Names (make.names versions)
Used for **formulas and coefficients**:
```r
condition <- make.names(condition.obj)  # 'sample_disease'
celltype <- make.names(celltype.obj)    # 'scArches_Cell_Annotation'
```

### Why This Matters

R formulas and edgeR coefficient names require syntactically valid R names:
- Spaces become dots: `"T cells"` → `"T.cells"`
- Special characters removed: `"CD4+ T"` → `"CD4..T"`

Always use:
- `.obj` versions when accessing data: `sce[[condition.obj]]`
- `make.names()` versions in formulas: `paste0(condition, make.names(group2))`
