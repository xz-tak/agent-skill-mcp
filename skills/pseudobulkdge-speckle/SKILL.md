---
name: pseudobulkdge-speckle
description: "Run single-cell RNA-seq pseudobulk differential expression analysis using R/Bioconductor (scran, edgeR, fgsea, speckle). Performs cell type-specific DEG analysis, pathway enrichment (Reactome, MSigDB), and cell composition differential analysis. Use when user mentions: pseudobulk analysis, DEG analysis, differential expression, cell type DEG, scRNA-seq comparison, cell composition analysis, speckle analysis, or cell proportion changes. Accepts Seurat RDS or h5ad input. Also invocable via /pseudobulkdge-speckle command."
version: "1.0.0"
---

# Pseudobulk DEG Analysis with Speckle

## Overview

This skill performs comprehensive single-cell RNA-seq pseudobulk differential expression analysis using R/Bioconductor packages. It handles:

1. **Pseudobulk DEG analysis** - Cell type-specific differential expression using scran/edgeR
2. **Pathway enrichment** - fgsea with Reactome, MSigDB (C2, C8, Hallmark), and custom signatures
3. **Cell composition analysis** - Differential cell proportions using speckle

## Performance Estimates

Runtime scales with number of cells, cell types, and comparisons. Based on benchmark with **172K cells, 49 cell types, 1 comparison**:

| Step | Description | Time | Output Size |
|------|-------------|------|-------------|
| 1-2 | Load + Aggregate | ~2 min | 1-2KB QC |
| 3 | DEG Analysis | ~1 min | 80-200MB |
| 4a | Reactome fgsea (groups) | ~15 min | 10-15MB |
| 4b | MSigDB (C8+C2+H) groups | ~20 min | 40-50MB |
| 4c | Custom signatures (groups) | ~1 min | <1MB |
| 5a | Cluster-level Reactome | ~15-80 min* | 20-30MB |
| 5b | Cluster-level MSigDB C8 | ~20 min | 15-20MB |
| 5c | Cluster-level Custom | ~3 min | <1MB |
| 6 | Cluster markers | ~4 min | 150-200MB |
| 7 | Speckle | ~1 min | <1MB |
| **Total** | **Full pipeline** | **~80-150 min** | |

*Step 5 timing varies significantly:
- With presto::wilcoxauc (fast): ~15 min
- With Seurat::FindAllMarkers fallback: ~80 min (used when presto fails on Seurat 5.0+)

**Scaling factors:**
- Each additional comparison: +35 min for Steps 3-4
- More cell types: linear increase in fgsea time
- More cells: moderate increase in DEG time, significant increase in Step 5

**Important:** Do NOT interrupt the pipeline during fgsea steps. Progress is logged but outputs are written at the end of each step.

## Workflow (9 Steps)

### Step 1: Input Handling

Accept input file and validate:

```
Supported formats:
- Seurat RDS (.rds, .RDS) - PREFERRED, loads directly
- AnnData h5ad (.h5ad) - Converts using anndata-seurat-conversion skill

If user provides h5ad:
1. Check if anndata-seurat-conversion skill is available
2. Convert h5ad to Seurat RDS
3. Continue with converted file
```

**Validation checks:**
- File exists and is readable
- Contains raw counts in appropriate slot
- Has metadata columns available

### Step 2: Metadata Inspection

Run the metadata inspection script to show available columns:

```bash
conda activate r_env && Rscript ${SKILL_DIR}/scripts/inspect_metadata.R <input_file>
```

The output will show grouped metadata summary:
```
=== Metadata Column Summary ===

Column: sample_disease
  - Unique values: 3
  - Examples: Healthy, IPF, PSS
  - Interpretation: Condition/disease status (likely comparison column)

Column: orig.ident
  - Unique values: 20
  - Examples: Sample1, Sample2, Sample3...
  - Interpretation: Sample identifier

Column: scArches_Cell_Annotation
  - Unique values: 15
  - Examples: T cells, B cells, Fibroblasts...
  - Interpretation: Cell type annotation
```

### Step 3: Parameter Alignment

Use AskUserQuestion to align on parameters ONE AT A TIME:

**Question 1: Comparison Column**
```
Which column contains the condition/disease status for comparison?
Options based on metadata inspection (show columns with 2-10 unique values)
```

**Question 2: Sample Column**
```
Which column identifies individual samples?
Options based on metadata inspection (show columns with many unique values)
```

**Question 3: Cell Type Column**
```
Which column contains cell type annotations?
Options based on metadata inspection (show columns with cell type-like values)
```

**Question 4: Batch Column (optional)**
```
Do you want to include batch correction?
Options: [No batch correction, <detected batch columns>]
```

**Question 5: Covariate Columns (optional)**
```
Do you want to include additional covariates?
Options: [No covariates, <potential covariate columns>]
```

### Step 4: Comparison Design

Use AskUserQuestion for comparison pairs:

**Question 1: Comparison pairs**
```
What conditions do you want to compare?
Based on the condition column, available values are: [list values]

Please specify comparison pairs as: group1 vs group2
Example: "Healthy vs Disease" or "Control vs Treatment1, Control vs Treatment2"
```

**Question 2: Custom Signatures (optional)**
```
Do you want to include custom gene signatures for enrichment analysis?
Options:
- No custom signatures
- GMT file (.gmt format)
- Gene list (comma-separated)
- Text file (one gene per row)
```

If custom signatures selected, ask for file path or gene list.

### Step 5: Interview & Confirmation

Summarize all parameters and confirm:

```
=== Analysis Configuration ===

Input file: /path/to/file.rds
Output directory: ./pseudobulk_results/

Columns:
  - Condition: sample_disease
  - Sample: orig.ident
  - Cell type: scArches_Cell_Annotation
  - Batch: None
  - Covariates: None

Comparisons:
  1. Healthy vs IPF
  2. Healthy vs PSS
  3. IPF vs PSS

Custom signatures: None

filterByExpr parameters:
  - min.count: 3
  - min.total.count: 5
  - min.prop: 0.2

Proceed with analysis? [Yes/No]
```

### Step 6: Execution (Long-running)

**CRITICAL: Do NOT interrupt the pipeline unless there is an error or user explicitly requests.**

Run the main R pipeline:

```bash
conda activate r_env && Rscript ${SKILL_DIR}/scripts/run_pseudobulk_pipeline.R \
  --input <input_file> \
  --output <output_dir> \
  --condition <condition_col> \
  --sample <sample_col> \
  --celltype <celltype_col> \
  --batch <batch_col> \
  --covariate <covariate_col> \
  --comparisons "<group1>,<group2>;<group1>,<group2>" \
  --custom_signatures <path_to_signatures> \
  --volcano_fdr <fdr_cutoff> \
  --volcano_log2fc <log2fc_cutoff>
```

**Volcano Plot Cutoffs (ask user before running):**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--volcano_fdr` | 0.05 | FDR/adjusted p-value cutoff for significance |
| `--volcano_log2fc` | 1 | Absolute log2 fold change cutoff |

Before running, confirm cutoffs with user:
```
Volcano plot significance cutoffs:
- FDR < 0.05 (default) - adjust?
- |log2FC| > 1 (default) - adjust?
```

The pipeline performs:

1. **Pseudobulk aggregation** - Aggregate cells by sample, condition, cell type
2. **pseudoBulkDGE** - Differential expression per comparison/cluster
3. **fgsea enrichment**:
   - Reactome pathways
   - MSigDB C2 (curated gene sets)
   - MSigDB C8 (cell type signatures)
   - MSigDB Hallmark
   - **Custom signatures** (if provided)
4. **Cluster markers** - findMarkers for cell type characterization
5. **Speckle cell composition** - propeller.ttest for differential proportions

**Custom Signature Integration:**
When custom signatures are provided, the pipeline runs additional fgsea:

```r
# Load custom signatures
custom_signatures <- load_custom_signatures(custom_sig_path)

# Run fgsea with custom gene sets
custom.fgsea <- fgsea(
  pathways = custom_signatures,
  stats = ranked_genes,
  minSize = 1,
  maxSize = Inf
)

# Save to separate output file
write.table(custom.fgsea,
  file = file.path(output_dir, "fgsea_custom_signatures.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
```

### Step 7: Sanity Check

Verify all outputs generated:

```
Expected outputs:
- [ ] pseudobulk_*_DE.tsv
- [ ] volcano_plots/*.png (static volcano plots)
- [ ] volcano_plots/*.html (interactive volcano plots with hover)
- [ ] fgsea_groups_clusters/*.tsv (pathway enrichment per comparison)
- [ ] fgsea_clusters/*.tsv + *.pdf (cluster-level enrichment)
- [ ] combined_markers_clusters.tsv
- [ ] speckle_diffprop/*.tsv (cell composition results)
- [ ] fgsea_custom_signatures.tsv (if custom signatures provided)
```

**Interactive Volcano Plots (HTML):**
Self-contained HTML files with interactive hover. Hover shows:
- **Gene name** (bold title)
- Log2 Fold Change
- P-value
- FDR/Adjusted P-value

Color coding based on user-defined cutoffs:
- Red: Up-regulated (FDR < cutoff, log2FC > cutoff)
- Blue: Down-regulated (FDR < cutoff, log2FC < -cutoff)
- Orange: Significant but small fold change
- Gray: Not significant

If any outputs are missing, report the issue and suggest troubleshooting steps.

### Step 8: Report Generation

**8a. Align with user on filter cutoffs:**

Use AskUserQuestion for cutoffs:

```
What filter cutoffs would you like for the report?

FDR/adjusted p-value cutoff: [Default: 0.05]
|log2FC| cutoff: [Default: 0.5]
NES cutoff for pathways: [Default: 0, meaning any non-zero]
Top k results per comparison: [Default: 30]
```

**8b. Filter and extract top results:**

Run the Python report generation script:

```bash
python ${SKILL_DIR}/scripts/generate_report.py \
  --output_dir <output_dir> \
  --fdr_cutoff <fdr> \
  --logfc_cutoff <logfc> \
  --nes_cutoff <nes> \
  --top_k <k>
```

This generates filtered result files for each comparison:
- `{comparison}_top_deg.tsv`
- `{comparison}_top_pathways.tsv`
- `{comparison}_speckle.tsv`

**8c. Claude Code Ultrathink Interpretation (REQUIRED):**

For EACH comparison, Claude Code MUST:

1. **Read the filtered results files:**
   ```
   {comparison}_top_deg.tsv
   {comparison}_top_pathways.tsv
   {comparison}_speckle.tsv
   ```

2. **Use ultrathink to generate interpretation:**
   - Analyze top upregulated genes: function, known disease associations
   - Analyze top downregulated genes: biological significance
   - Interpret pathway enrichment: mechanistic insights
   - Interpret cell composition changes: cellular basis of disease
   - Cross-reference with literature using WebFetch when needed

3. **Append the ultrathink interpretation to the report markdown**

4. **NEVER skip ultrathink interpretation - it is REQUIRED for each comparison**

### Step 9: Report Structure

The final report follows this structure:

```markdown
# Pseudobulk DEG Analysis Report

## Executive Summary
- Key findings across all comparisons (2-3 sentences per comparison)
- Most significant cell types affected
- Top pathways/biological processes

## Input Data
- Input file path and type
- Total cells, genes, cell types
- Sample distribution per condition

## Methods & Parameters
- All configuration parameters
- Filter thresholds (min.count, min.total.count, min.prop)
- Comparison design matrix
- Databases used (Reactome, MSigDB collections)

## Results

### Comparison 1: {Group2} vs {Group1}

#### Differential Expression Summary
| Gene | log2FC | FDR | Direction |
|------|--------|-----|-----------|
(Top 15 up, Top 15 down)

**Volcano Plot:** ![Volcano]({output_dir}/volcano_{comparison}.png)

#### Pathway Enrichment Analysis

##### Reactome Pathways
| Pathway | NES | FDR | Leading Edge |
(Top 10 up, Top 10 down)

##### MSigDB Hallmark
| Gene Set | NES | FDR |
(Top 10)

##### Custom Signatures (if provided)
| Signature | NES | FDR | Leading Edge |
(All significant)

#### Cell Composition Changes (Speckle)
| Cell Type | PropRatio | FDR | Interpretation |
(Significant only, FDR < 0.05)

#### Biological Interpretation
[ULTRATHINK SECTION - Deep analysis of:]
- Gene function and pathway connections
- Disease/condition relevance
- Therapeutic implications
- Cell type-specific biology

[WEBFETCH SECTION - Literature support:]
- Query PubMed/literature for top genes
- Known disease associations
- Drug target potential

### Comparison 2: ...
(repeat structure)

## Discussion

### Summary of Key Findings
- Cross-comparison patterns
- Consistent vs unique findings

### Biological Significance
[ULTRATHINK - Professional biologist interpretation]

### Clinical/Translational Relevance
- Biomarker potential
- Drug target candidates
- Patient stratification implications

### Limitations
- Sample size considerations
- Technical limitations

### Future Directions
- Validation experiments
- Follow-up analyses
```

## R Environment Detection

The skill auto-detects the R environment:

```r
detect_r_env <- function() {
  # Check if r_env conda environment exists
  conda_envs <- system("conda env list", intern = TRUE)
  if (any(grepl("r_env", conda_envs))) {
    return("conda activate r_env && Rscript")
  }
  # Check if required packages available in current R
  required <- c("Seurat", "scran", "edgeR", "fgsea", "speckle")
  missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  if (length(missing) == 0) {
    return("Rscript")
  }
  stop("Missing packages: ", paste(missing, collapse = ", "))
}
```

## Custom Signature Support

Supports multiple input formats:

```r
load_custom_signatures <- function(input) {
  if (file.exists(input) && grepl("\\.gmt$", input)) {
    # GMT file format
    return(fgsea::gmtPathways(input))
  } else if (file.exists(input)) {
    # Text file (one gene per row)
    genes <- readLines(input)
    return(list("custom_signature" = genes))
  } else {
    # Direct gene list (comma or newline separated)
    genes <- unlist(strsplit(input, "[,\n]+"))
    return(list("custom_signature" = trimws(genes)))
  }
}
```

## Resources

### scripts/
- `run_pseudobulk_pipeline.R` - Main R pipeline script
- `inspect_metadata.R` - Metadata inspection helper
- `generate_report.py` - Report generation with filtering

### references/
- `parameter_guide.md` - Detailed parameter documentation
- `output_schema.md` - Output file format specifications

## Interview Pattern

Follow these interview guidelines:
- Ask questions ONE AT A TIME
- Prefer multiple choice when possible
- Validate after each section
- Focus on: purpose, constraints, success criteria
- Use AskUserQuestion tool for all user interactions

## Long-running Process Handling

- Use `message()` for progress updates in R
- Checkpoint after each major step
- Do NOT interrupt unless error or user request
- If interrupted, report progress and suggest resume options

## Auto-trigger Keywords

This skill triggers when user mentions:
- "pseudobulk" / "pseudo-bulk"
- "DEG analysis" / "differential expression"
- "cell type specific DEG"
- "scRNA-seq comparison"
- "cell composition" / "cell proportion"
- "speckle analysis"
- "condition comparison scRNA"

## Required R Packages

```r
# Core analysis
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(scuttle)
library(edgeR)

# Pathway analysis
library(fgsea)
library(msigdbr)
library(reactome.db)
library(org.Hs.eg.db)

# Cell composition
library(speckle)
library(limma)

# Utilities
library(dplyr)
library(tibble)
library(glue)
library(presto)
```
