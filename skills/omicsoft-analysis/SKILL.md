---
name: omicsoft-analysis
description: Analyze pre-computed differential gene expression statistics from Omicsoft h5ad files or TileDB-SOMA experiments on S3 to investigate target gene expression patterns and signature enrichment, and single-cell RNA-seq target expression analysis with AUCell scoring across predefined scRNA-seq datasets. This skill should be used when users request to analyze bulk DEG data (log2fc, padj) from h5ad files or SOMA experiments, filter for specific diseases or studies, examine target genes or pathways with custom signatures, perform GSEA enrichment analysis with enhanced summaries including leading edge annotation, and generate interactive visualizations. Supports auto-detection of SOMA vs h5ad format for S3 URIs.
---

# Omicsoft DEG Analysis

Analyze pre-computed differential gene expression statistics from Omicsoft bulk transcriptomics data stored in h5ad format to investigate target gene expression patterns and signature enrichment.

## Overview

Analyze pre-computed DEG statistics (log2fc and padj) to filter h5ad files for specific diseases and/or studies, investigate custom gene signatures and target pathway expression patterns, generate interactive visualizations (scatter plots and hierarchically-clustered heatmaps), perform GSEA signature enrichment analysis with enhanced summaries (including significance counts, NES direction separation, and leading edge annotation), and export detailed results with configurable significance thresholds. **Step 4** (optional) extends analysis to sample-level expression data with per-study boxplots, GSVA scoring, and correlation modules.

**Automatic Data Cleanup**: When loading `*_deg.h5ad` files, embedded expression data in `uns` (e.g., `uns['xxx_expr']`) is automatically removed to reduce memory usage and focus on DEG analysis. See `references/anndata_schema.md` for details on the AnnData object structure.

## Workflow

### Step 0: Explore H5AD Schema (Optional but Recommended)

**When to use**: When users are unfamiliar with the h5ad file contents, need to discover available filter values, or want to explore metadata before analysis.

Generate an interactive schema browser to help users discover available filter values:

```bash
# Generate comprehensive JSON schema
conda run -n <env_name> python scripts/explore_h5ad_schema.py \
  --file <path_to_h5ad_file> \
  --output schema_report.json \
  --format json

# Create standalone HTML viewer (works without web server)
conda run -n <env_name> python scripts/generate_schema_viewer.py \
  --json schema_report.json \
  --output schema_viewer.html
```

The schema viewer provides:
- **All unique values** for every metadata column (no truncation)
- **Demographics**: gender, ethnicity, age categories
- **Filtering columns**: diseases (178 values), tissues (131 values), treatments, studies
- **Interactive search** and category filtering
- **Count and percentage** for each value
- **Visual distributions** with progress bars

Users can open `schema_viewer.html` directly in a browser to explore all available filter values before running the analysis.

**Example Files**: The skill includes real example outputs for reference:
- `references/schema_report_example.json`: Complete schema JSON from IBD/MASH/Fibro/Derm/Rheum dataset (7,604 observations × 60,071 genes) showing all 40 metadata columns with full value lists
- `assets/schema_viewer_example.html`: Standalone interactive HTML viewer with embedded data - open directly in any browser to explore the example dataset
- `references/schema_quickstart.md`: Quick start guide with troubleshooting and usage tips

For comprehensive documentation on schema exploration tools and usage, refer to `references/schema_exploration.md`.

### Step 1: Gather User Requirements

Ask the user for the following information:

1. **H5ad file location**: Path to the Omicsoft h5ad file containing pre-computed DEG statistics
   - File must contain log2fc and padj layers
   - Can be local path or S3 URI (s3://bucket/key)
   - S3 URIs are loaded directly if accessible (no download needed)

2. **Target name**: Descriptive name for the analysis (e.g., "CGAS", "Fibrosis", "JAK_STAT")

3. **Gene signatures** (optional, paired with targets):
   - Format: `"SignatureName:Gene1,Gene2,Gene3;AnotherSignature:GeneA,GeneB"`
   - Example: `"CGAS_STING:CGAS,TMEM173,TBK1,IRF3;TGFB:TGFB1,TGFB2,TGFB3"`
   - For SC analysis: signatures are positionally paired with target sets — their genes are MERGED with the paired target genes for AUCell scoring
   - **IMPORTANT**: Gene symbols must match organism:
     - Human: ALL UPPERCASE (e.g., `CGAS`, `TMEM173`)
     - Mouse: First letter uppercase, rest lowercase (e.g., `Cgas`, `Tmem173`)
   - Refer to `references/signature_format.md` for detailed formatting guide

3b. **Addon signatures** (optional, NOT paired with targets):
   - Additional gene signatures scored independently via AUCell — genes are NOT merged with target-paired signatures
   - Format: same as signatures `"Name:Gene1,Gene2;Name2:GeneA,GeneB"`
   - Use for pathway/module signatures that should be analyzed alongside but are conceptually separate from targets
   - Example: `--addon-signatures "Fibrosis_module:COL1A1,COL3A1,FN1;Inflammation:TNF,IL6,IL1B"`
   - For SC script: passed directly via `--addon-signatures` (not positionally paired)

   **When gathering requirements, use AskUserQuestion to ask:**
   - "Do you have target-signature pairs? (targets paired positionally with signatures for AUCell scoring)"
   - "Do you have additional standalone signatures? (scored independently, not merged with targets)"

4. **Optional filtering parameters** (all filters are optional and default to None):
   - `--diseases`: Comma-separated disease keywords (e.g., "scleroderma,sclerosis", "crohn,colitis,IBD"). Partial matches supported (case-insensitive)
   - `--exclude-diseases`: Comma-separated disease keywords to exclude (e.g., "ALS,amyotrophic lateral sclerosis"). Applied after --diseases filter. Useful for including broad disease terms while excluding specific subtypes
   - `--studies`: Comma-separated list of study names (e.g., "GSE12345,GSE67890")
   - `--tissues`: Comma-separated tissue keywords (e.g., "skin,blood,lung"). Partial matches supported (case-insensitive)
   - `--comparison-category`: Comma-separated comparison categories (e.g., "Disease vs. Normal")
   - `--case-treatment`: Comma-separated treatment keywords (e.g., "none,NA")
   - `--comparison`: Comma-separated comparison keywords for fuzzy/partial matching (e.g., "response vs no response")
   - If no filters are specified, all data in the h5ad file will be included

5. **Optional target genes**:
   - `--targets`: Comma-separated list of target genes for tracking in leading edge and plotting (e.g., "GENE1,GENE2,GENE3")
   - These targets will be plotted like signatures but NOT included in GSEA enrichment
   - The enhanced summary will track only these targets in the leading edge

6. **Optional thresholds** (use defaults if not specified):
   - `--lfc-threshold`: Absolute log2 fold change threshold (default: 0.0)
   - `--padj-threshold`: Adjusted p-value threshold (default: 0.05)

7. **Optional flags**:
   - `--run-gsea`: Enable GSEA enrichment analysis (**ON by default**, includes MSigDB_Hallmark_2020)
   - `--no-gsea`: Disable GSEA analysis
   - `--pathwaydb-gsea-only`: Run ONLY PathwayDB GSEA (indication-centric). Skips target/signature analysis, default MSigDB GSEA, and report generation. Uses combined_pathways.csv from `/pathwaydb-query` skill. Does NOT require `--signatures`. **DISABLED by default — never include this flag unless the user explicitly requests PathwayDB-only/indication-centric GSEA mode.**
   - `--pathwaydb-file`: Optional explicit path to pathwaydb CSV file (auto-detected if not specified)
   - `--indication`: Indication name for output directory naming (e.g., "IBD", "Fibrosis"). Used in output path: `{target}_{indication}_omicsoft_YYYYMMDD/`
   - `--n-workers`: Number of parallel workers for GSEA (default: min(4, cpu_count)). Increase for faster processing on multi-core systems.
   - `--output-dir`: Output directory name (default: `{target_name}_{indication}_omicsoft_YYYYMMDD/` in current working directory)

> **IMPORTANT**: The `--pathwaydb-gsea-only` flag must NEVER be added by default. Only include it when the user explicitly requests "PathwayDB only", "indication-centric GSEA", or specifically asks to skip target/signature analysis. The default workflow is target mode with `--signatures`.

8. **Single-cell RNA-seq analysis** (default: YES):
   - Ask user: "Include single-cell RNA-seq target analysis? Runs AUCell scoring across 8 scRNA-seq datasets (~2-4 hours). Skip with 'no-scrna'."
   - Uses the SAME --targets and --signatures from above
   - **Important**: For sc analysis, targets and signatures are PAIRED positionally (1st target set with 1st signature, etc.)
   - When collecting targets/signatures in Step 1, ask: "How many target-signature sets?" then collect each pair iteratively
   - Format for sc script: `--targets "set1;set2;set3" --signatures "sig1:genes;sig2:genes;sig3:genes"`
   - Addon signatures (if any): `--addon-signatures "Name1:genes;Name2:genes"` (passed directly, not paired with targets)
   - Both `--signatures` and `--addon-signatures` are optional — at minimum only `--targets` is required
   - For bulk DEG (Step 3), flatten all target sets into one `--targets` comma-separated list
   - Available studies: jhu, rieder, umcg_fibroblast, otar, tacolny, pf_lung, ssc_lung, hs_skin

### Step 2: Verify Environment

Before running analysis, ensure the conda environment has required packages:

```bash
conda run -n <env_name> pip install scanpy s3fs anndata gseapy plotly upsetplot kaleido boto3
```

Replace `<env_name>` with the user's environment name.

**Note**: `s3fs` and `boto3` are required for loading files from S3 URIs. If you're only using local files, these packages are optional.

### Step 2.5: Validate Filters Step-by-Step (Highly Recommended)

**When to use**: Before running the full analysis, validate all filters step-by-step to ensure they produce meaningful results and avoid running a full analysis that returns 0 observations.

**Why this matters**: Filter validation tests each filter incrementally (disease → tissue → studies → comparison_category, etc.) and shows:
- Which values match your query terms
- How many observations remain after each filter
- Whether any filter results in 0 observations
- Suggestions for fixing problematic filters

This prevents wasting time on full analyses that fail due to incorrect filter values or overly restrictive combinations.

#### Running Filter Validation

```bash
conda run -n <env_name> python scripts/validate_filters.py \
  --file <path_to_h5ad_file> \
  --target-name <target_name> \
  --signatures <signature_string> \
  [--diseases <disease_keywords>] \
  [--exclude-diseases <disease_keywords_to_exclude>] \
  [--tissues <tissue_keywords>] \
  [--studies <study_names>] \
  [--comparison-category <categories>] \
  [--case-treatment <treatments>] \
  [--comparison <comparison_keywords>] \
  [--targets <target_genes>] \
  [--lfc-threshold <value>] \
  [--padj-threshold <value>] \
  [--no-gsea]  # GSEA runs by default; use this to disable
```

**Example validation:**

```bash
conda run -n claude_test python scripts/validate_filters.py \
  --file /path/to/ibd_deg.h5ad \
  --target-name IBD_Targets \
  --diseases "crohn's disease,ulcerative colitis,inflammatory bowel disease" \
  --tissues "intestine,colon,rectum,ileum,sigmoid" \
  --comparison-category "Disease vs. Normal,Responder vs. Non-Responder" \
  --signatures "IBD_Targets:JAK1,TYK2,ITGA4,TNFRSF25,ITGB1,PCOLCE,GREM1,CDKN2D" \
  --targets "JAK1,TYK2,ITGA4,TNFRSF25,ITGB1,PCOLCE,GREM1,CDKN2D" \
  --lfc-threshold 0.0 \
  --padj-threshold 0.1
  # GSEA runs by default; add --no-gsea to disable
```

#### Understanding Validation Output

The validation script shows step-by-step progress:

```
[STEP 0] Loading h5ad file...
  ✓ Loaded: (7,604, 60,071) (obs x vars)
  Total observations: 7,604

[STEP 1] Disease Filter
  Query: crohn's disease,ulcerative colitis,inflammatory bowel disease
  Filter type: Substring match (case-insensitive)
  Search terms: ["crohn's disease", "ulcerative colitis", "inflammatory bowel disease"]

  ✓ Found 4 matching disease(s):
    - crohn's disease (CD): 698 obs
    - inflammatory bowel disease (IBD): 40 obs
    - ulcerative colitis (UC): 692 obs
    - primary sclerosing cholangitis (PSC);ulcerative colitis (UC): 3 obs

  → Observations after disease filter: 1,433

[STEP 2] Tissue Filter
  Query: intestine,colon,rectum,ileum,sigmoid
  Filter type: Substring match (case-insensitive)

  ✓ Found 14 matching tissue(s):
    - colon: 59 obs
    - colonic mucosa: 209 obs
    - ileum: 131 obs
    - sigmoid colon: 80 obs
    ... (10 more tissues)

  → Observations after tissue filter: 647

[STEP 4] Comparison Category Filter
  Query: Disease vs. Normal,Responder vs. Non-Responder
  Filter type: Exact match

  Available comparison categories in current dataset: 8
    - Disease vs. Normal: 35 obs
    - Responder vs. Non-Responder: 35 obs
    - Other Comparisons: 129 obs
    ... (5 more categories)

  ✓ Found 2 matching category(ies):
    - Disease vs. Normal: 35 obs
    - Responder vs. Non-Responder: 35 obs

  → Observations after comparison category filter: 70

FILTER VALIDATION SUMMARY
✓ All filters validated successfully!
  Initial observations: 7,604
  Final observations: 70
  Reduction: 99.1%

✓ Ready to proceed with full analysis
```

#### Key Features of Validation

**1. Filter Type Identification**
- **Substring match (case-insensitive)**: Matches partial strings (diseases, tissues, comparison)
- **Exact match**: Requires exact value (studies, comparison_category, case_treatment)
- **Fuzzy match**: Flexible substring matching (comparison field)

**2. Intermediate Results**
- Shows observations remaining after each filter
- Helps identify which filters are too restrictive
- Displays matched values before applying filter

**3. Zero Observation Detection**
If any filter results in 0 observations, validation stops and provides:
- Which filter caused the issue
- Available values in the filtered dataset
- Suggestions for alternative filter terms
- Reference to schema viewer for valid values

**Example failure scenario:**

```
[STEP 2] Tissue Filter
  Query: brain

  ✗ No matching tissues found in current filtered dataset!

  Available tissue keywords (showing first 20):
    - colon
    - colonic mucosa
    - ileum
    - sigmoid colon
    ...

VALIDATION FAILED
✗ Failed at filter: tissue

Suggestions:
  - No tissues match 'brain' in the disease-filtered dataset.
  - Please check available tissues in schema viewer.
  - Try broader tissue terms or check tissue-disease combinations.
```

#### When Validation Succeeds

The script outputs the complete command to run the full analysis:

```bash
conda run -n <env_name> python scripts/deg_analysis.py \
  --file <path> \
  --target-name <name> \
  --signatures <sigs> \
  --diseases "<terms>" \
  --tissues "<terms>" \
  --comparison-category "<categories>" \
  --targets "<genes>" \
  --lfc-threshold 0.0 \
  --padj-threshold 0.1 \
  --run-gsea
```

Copy this command and proceed to Step 3 to run the full analysis.

#### Troubleshooting Filter Issues

If validation fails:

1. **Check available values**: Run schema exploration (Step 0) to see all valid values
2. **Use broader terms**: Try more general keywords (e.g., "colitis" instead of "ulcerative colitis (UC)")
3. **Remove restrictive filters**: Some filter combinations may not exist (e.g., specific study + specific tissue)
4. **Verify exact matches**: Comparison categories require exact matches - check schema for precise category names
5. **Check tissue-disease combinations**: Not all diseases have data for all tissues

### Step 3: Execute Analysis Script

Run the analysis script `scripts/deg_analysis.py` with the collected parameters:

```bash
conda run -n <env_name> python scripts/deg_analysis.py \
  --file <path_to_h5ad_file> \
  --target-name <target_name> \
  --signatures <signature_string> \
  [--diseases <disease_keywords>] \
  [--exclude-diseases <disease_keywords_to_exclude>] \
  [--studies <study_names>] \
  [--tissues <tissue_keywords>] \
  [--comparison-category <categories>] \
  [--case-treatment <treatments>] \
  [--comparison <comparison_keywords>] \
  [--targets <target_genes>] \
  [--lfc-threshold <value>] \
  [--padj-threshold <value>] \
  [--output-dir <output_directory>] \
  [--run-gsea]
```

**Example command:**

```bash
conda run -n claude_test python scripts/deg_analysis.py \
  --file data/omicsoft_deg.h5ad \
  --target-name CGAS_STING \
  --diseases "scleroderma,sclerosis" \
  --signatures "CGAS_STING:CGAS,MB21D1,TMEM173,STING1,TBK1,IKBKE,IRF3,IRF7,IFNB1;TGFB:TGFB1,TGFB2,TGFB3,TGFBR1,TGFBR2" \
  --lfc-threshold 0.5 \
  --padj-threshold 0.05 \
  --run-gsea
```

### Step 4: Sample-Level Expression Analysis

**Runs by default** after Step 3 (DEG analysis). Provides per-study expression boxplots, GSVA scoring, and correlations.

**What it does**: Queries the EXPR SOMA for sample-level normalized expression,
queries the DEG SOMA for significance annotations, then runs three modules per study:
- **comparison**: Grouped boxplots of target gene expression with DEG significance asterisks
- **gsva**: GSVA signature scoring with limma-based group comparison (adj.P.Val < 0.05)
- **corr**: Gene-gene correlation heatmaps, target vs continuous metadata scatter plots

**Study behavior**:
- Internal studies (Engitix_FFPE, SPARC, Varsity, Yokohama_RNA, Yokohama_Protein): ALWAYS run
- Curated studies: ALWAYS run (when available)
- Omicsoft studies: Ask user which studies to include

**No disease/tissue/comparison filters are applied** — each study runs with its full sample set.
Filters from Step 3 (DEG analysis) do NOT carry over.

#### Running Sample-Level Analysis

After Step 3 completes, ask the user which omicsoft studies to include (if any).
Then construct the Python command:
```bash
conda run -n spatial python ~/.claude/skills/omicsoft-analysis/scripts/sample_level_analysis.py \
  --expr-uri "s3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/ibd_mash_fibro_derm_rheum_internal_05302026/ibd_mash_fibro_derm_rheum_internal_05302026_expr/" \
  --deg-uri "s3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/ibd_mash_fibro_derm_rheum_internal_05302026/ibd_mash_fibro_derm_rheum_internal_05302026_deg/" \
  --target-name <target> \
  --targets <target_genes_from_step1> \
  --signatures <signatures_from_step1> \
  --output-dir <deg_output_dir>/sample_level \
  --per-sample-studies <any_additional_omicsoft_studies> \
  --modules comparison,gsva,corr
```

> **Note**: Python version uses gseapy for GSVA (replaces R GSVA package), statsmodels OLS+BH for group comparison (replaces limma), and scipy for correlations. Config JSON auto-resolves from the skill directory.

**Legacy R command** (deprecated, kept for reference):
```bash
conda run -n claude_test Rscript ~/.claude/skills/omicsoft-analysis/scripts/sample_level_analysis.R \
  --expr-uri "s3://..." --deg-uri "s3://..." \
  --target-name <target> --targets <genes> --signatures <sigs> \
  --output-dir <output>/sample_level --modules comparison,gsva,corr --conda-env claude_test
```

#### Output Structure

Nested under the DEG output directory:
```
{target}_omicsoft_YYYYMMDD/sample_level/
├── manifest.json
├── Engitix_FFPE/
│   ├── target_expression.png / .html
│   ├── signature_gsva.png / .html
│   ├── target_comparison_stats.csv
│   ├── gsva_comparison_stats.csv
│   ├── analysis_log.txt
│   └── correlations/{gene_gene,target_vs_clinical,gsva_vs_clinical}/
├── SPARC/{CD,UC}/...
├── Varsity/...  (faceted by treatment)
├── Yokohama_RNA/{fibrosis,nash,diagnosis}/...
├── Yokohama_Protein/{fibrosis,nash}/...
└── [optional_omicsoft_study]/{comparison_name}/...
```

#### Prerequisites

conda spatial env must have: numpy, pandas, scipy, statsmodels, gseapy, matplotlib, seaborn, plotly, boto3, tiledbsoma (optional)

Legacy R prerequisites (deprecated): conda r_env must have: R, ggplot2, dplyr, tidyr, tibble, plotly, htmlwidgets, GSVA, limma, RColorBrewer, jsonlite, argparse

#### References

- `references/sample_level_analysis.md` — detailed module documentation and caveats
- `references/internal_study_configs.json` — internal study configuration
- `references/correlation_analysis.md` — correlation module specifics

### Step 4.5: Single-Cell RNA-seq Target Analysis

**Runs by default** after Step 4 (sample-level analysis), unless the user declined in Step 1 (said "no single cell", "no-scrna", "skip sc", etc.).

**What it does**: Processes 8 predefined scRNA-seq datasets from S3, performing:
- AUCell signature enrichment scoring
- Embedding plots (UMAP or study-specific embeddings) colored by condition, cell type, gene expression, and signature scores
- Dotplots per condition and cell type
- Stacked violin plots for AUCell scores
- Correlation heatmaps (target gene and signature correlations, global and per-cell-type)
- Population composition analysis (if xzsc installed)
- Per-cell-type sub-analyses (dotplots, correlations)

**Command**:
```bash
conda run -n spatial python ~/.claude/skills/omicsoft-analysis/scripts/sc_target_analysis.py \
  --targets <targets_from_step1> \
  --signatures <signatures_from_step1> \
  --output-dir <deg_output_dir>/sc_target_result \
  [--studies jhu,rieder,umcg_fibroblast,otar,tacolny,pf_lung,ssc_lung,hs_skin] \
  [--max-workers 1] \
  [--temp-dir /mnt/sagemaker-nvme/tmp]
```

**Parameter translation** from Step 1:
- `--targets`: Target genes, semicolon-separated per set for paired matching: `"IL11,IL11RA;TNFRSF25,TNFSF15"`
- `--signatures`: Signatures paired positionally with target sets: `"IL11_sig:IL13RA2,CEMIP;TL1A_sig:CHI3L1,MMP7"`
- `--output-dir`: Nested inside the bulk DEG output directory as `sc_target_result/`

**Target-signature pairing** (critical for AUCell scoring):
- Each target set (separated by `;`) is unioned into its positionally-matched signature for AUCell scoring
- Example: `--targets "IL11,IL11RA;TNFRSF25,TNFSF15" --signatures "IL11_sig:IL13RA2,CEMIP;TL1A_sig:CHI3L1,MMP7"`
  - IL11_sig scores: IL13RA2, CEMIP + IL11, IL11RA (1st target set → 1st signature)
  - TL1A_sig scores: CHI3L1, MMP7 + TNFRSF25, TNFSF15 (2nd target set → 2nd signature)
- If only 1 target set (no semicolons): all targets unioned into ALL signatures (backward compatible with `deg_analysis.py` format)

**Environment**: `conda run -n spatial` (requires scanpy, decoupler, plotly, boto3; xzsc optional)

**Runtime**: ~15-30 minutes per study, ~2-4 hours total for all 8 studies

**Available studies**:
| Study | Disease Area | Key Features |
|-------|-------------|--------------|
| jhu | IBD | Patient-level scRNA-seq |
| rieder | IBD | Uses MDE scANVI embedding (not UMAP) |
| umcg_fibroblast | IBD | Stromal/fibroblast focus |
| otar | MASH/NASH | Liver snuc-seq |
| tacolny | MASH/NASH | Liver snuc-seq |
| pf_lung | Pulmonary fibrosis | Multi-study lung atlas |
| ssc_lung | Systemic sclerosis | Lung atlas integration |
| hs_skin | Hidradenitis suppurativa | Skin atlas |

**Output structure**:
```
{deg_output_dir}/sc_target_result/
├── jhu/
│   ├── jhu_report.pdf
│   ├── umap_interactive.html
│   ├── heatmap_interactive.html
│   └── all_plots_interactive.html
├── rieder/
│   └── ... (same structure, uses MDE embedding)
├── umcg_fibroblast/
├── otar/
├── tacolny/
├── pf_lung/
├── ssc_lung/
└── hs_skin/
```

**Error handling**:
- Individual study failures print `[FAIL] study_name: error` but don't halt pipeline
- Missing embeddings: falls back to X_umap, then skips embedding plots with info message
- Single-item correlations: skipped silently with info message if <2 targets or <2 signatures
- xzsc not installed: population composition skipped with info message

**Notes**:
- Memory-intensive: each study loads a full h5ad from S3 into memory
- Studies are processed sequentially by default (`--max-workers 1`) to avoid OOM
- The `rieder` study uses `X_mde_scANVI_240923` embedding instead of standard UMAP
- All other studies use `X_umap`

### Step 5: Interpret Results

After all steps complete, check outputs:
- Check `{output_dir}/sc_target_result/` for per-study single-cell results (if Step 4.5 ran)
- Report any `[FAIL]` lines from sc analysis stdout
- Note which studies completed successfully vs failed

The analysis generates multiple output files in the output directory (default: `{target_name}_omicsoft/` in current working directory):

1. **study_summary.csv** / **study_summary.json**: List of studies included after filtering
   - Columns: study, tissue, disease, disease_category
   - JSON version provided for programmatic access

2. **target_summary.csv** / **target_summary.json**: Target gene expression patterns
   - Shows which genes from `--targets` are up/down-regulated per tissue, disease, and comparison category
   - Grouped by: Gene, tissue, disease_category, comparison_category, sig (up/dn)
   - Columns: Gene, tissue, disease_category, comparison_category, sig (up/dn), count, study, total_count, threshold
   - JSON version provided for programmatic access
   - Note: Renamed from `signature_summary.csv` to clarify that it contains target genes, not custom signatures

3. **target_signature_YYYYMMDD.html**: Interactive scatter plot of gene expression
   - Shows log2fc for each gene across studies
   - Dropdown menu to switch between signatures (includes target genes if provided)
   - Bubble size indicates -log(padj)
   - Green dashed lines mark threshold boundaries
   - Hover data includes: log2fc, padj, tissue, disease_category, comparison, case_treatment_status, control_treatment_status

4. **comparison_signature_YYYYMMDD.html**: Interactive clustered heatmap of gene expression
   - Shows log2fc heatmap for each signature with hierarchical clustering
   - Both rows (comparisons) and columns (genes) are clustered using average linkage with Euclidean distance
   - Clustering helps reveal expression patterns by grouping similar samples and genes together
   - Y-axis labels display: {study}_{disease_category}_{comparison} for full context
   - Asterisks (*) mark significant genes (padj < threshold)
   - Dropdown menu to switch between signatures (includes target genes if provided)

5. **detailed_gene_table.csv** / **detailed_gene_table.json**: Comprehensive gene-level data export
   - Long-format table with all metadata for signature and target genes
   - Each row represents one gene in one comparison
   - Contains all columns: log2fc, padj, tissue, disease, comparison, study, source, weight, etc.
   - Weight column: internal=2.0, curated=1.5, omicsoft=1.0 (applied when source column exists)
   - Useful for custom downstream analysis and filtering
   - JSON version provided for programmatic access

6. **internal_vs_external_summary.csv** / **internal_vs_external_summary.json**: Source-stratified summary
   - Per-gene statistics stratified by source (internal, curated, omicsoft)
   - Columns: Gene, source, n_comparisons, n_significant, mean_log2fc, median_log2fc, weight
   - Includes WEIGHTED_AGGREGATE rows with weighted mean log2fc across all sources
   - Useful for comparing internal study findings against external validation

7. **gsea_all_results.csv** / **gsea_all_results.json**: Complete GSEA enrichment results (if --run-gsea used)
   - All signature enrichment results across studies
   - Includes both custom signatures and MSigDB_Hallmark_2020
   - Contains comparison_category for tracking enrichment context
   - JSON version provided for programmatic access

8. **gsea_enhanced_summary.csv** / **gsea_enhanced_summary.json**: Enhanced GSEA summary with detailed statistics (if --run-gsea used)
   - **Groupby format**: Each Term has separate rows for each Comparison_Category
   - Allows comparison of enrichment patterns between "Disease vs. Normal" and "Responder vs. Non-Responder"
   - Number of significant comparisons (FDR < threshold) per term-category combination
   - Studies with positive NES (upregulated) and their NES values
   - Studies with negative NES (downregulated) and their NES values
   - Input target genes found in leading edge (tracks only --targets genes if provided, otherwise all signature genes)
   - Columns include: Term, Comparison_Category, N_Significant_Comparisons, Total_Comparisons, Percent_Significant, N_Positive_NES, N_Negative_NES, Positive_NES_Studies, Negative_NES_Studies, Input_Targets_In_Leading_Edge, N_Input_Targets_In_LE
   - Note: Target genes (--targets) are plotted but NOT included in GSEA enrichment analysis
   - JSON version provided for programmatic access

9. **gsea_score_pathway_YYYYMMDD.html**: GSEA visualization (if --run-gsea used)
   - Interactive scatter plot of normalized enrichment scores (NES)
   - Bubble size indicates significance

10. **gsea_{target}_{indication}_pathwaydb.csv** / **gsea_{target}_{indication}_pathwaydb.json**: Pathway database GSEA results (if --pathwaydb-gsea-only used)
   - GSEA results using pathways from pathwaydb-query skill output
   - Includes Database, Database_ID columns for traceability
   - Contains tissue, disease metadata columns
   - Columns: Term, Database, Database_ID, Pathway_Name, NES, NOM p-val, FDR q-val, FWER p-val, Tag %, Gene %, Lead_genes, study, comparison_category, tissue, disease
   - JSON version provided for programmatic access

11. **gsea_{target}_{indication}_pathwaydb_fdr0.05.csv** / **gsea_{target}_{indication}_pathwaydb_fdr0.05.json**: Filtered pathway GSEA results (if --pathwaydb-gsea-only used)
   - Only significant pathways with FDR < 0.05
   - Same columns as the full pathwaydb GSEA results
   - Use for downstream biomedical report generation
   - JSON version provided for programmatic access

### Step 6: Summarize Findings

After analysis completes, provide the user with:

1. **Dataset summary**:
   - Number of observations (comparisons) before and after filtering
   - Diseases found and filtered
   - Number of unique studies

2. **Key findings**:
   - Which signature genes show significant differential expression
   - Tissue, disease, and comparison category expression patterns
   - Top enriched signatures (if GSEA was run)
   - Enrichment statistics by comparison category: number of significant comparisons, NES direction, and input targets in leading edge
   - Patterns across "Disease vs. Normal" and "Responder vs. Non-Responder" comparisons

3. **Output file locations**:
   - List all generated files with brief descriptions
   - Note that HTML files are interactive and can be opened in a browser

## Important Notes

### Pre-computed DEG Statistics

This analysis works with h5ad files that contain **pre-computed** differential expression statistics:
- **log2fc layer**: Log2 fold changes from disease vs. control comparisons
- **padj layer**: Adjusted p-values (FDR corrected)
- Each observation represents a comparison (e.g., disease vs. normal in a specific study)

The analysis does NOT compute differential expression - it analyzes existing DEG results to identify target gene expression patterns and signature enrichment.

### AnnData Object Structure and Automatic Cleanup

**DEG Files**: `*_deg.h5ad` files may contain embedded expression data in `uns` using keys like `{filename}_expr`. For example:
- File: `ibd_mash_fibro_derm_rheum_10222025_deg.h5ad`
- Embedded data: `uns['ibd_mash_fibro_derm_rheum_10222025_expr']`

**Automatic Memory Optimization**: The `load_h5ad()` function automatically detects and removes this embedded expression data because:
1. Expression data is typically 10-100x larger than DEG statistics
2. DEG analysis only needs pre-computed log2fc and padj values
3. Removing it prevents memory issues and speeds up processing
4. Original expression data can be accessed from separate `*_expr.h5ad` files if needed

**Data Structure**: See `references/anndata_schema.md` for complete documentation of:
- DEG AnnData schema (obs, var, layers, uns structure)
- Expression AnnData schema (sample-level data)
- Metadata column descriptions
- Sparse matrix storage format
- Memory optimization strategies

### Gene Symbol Formatting

**Critical**: Gene symbols must match the organism in the h5ad file:

- **Human data**: Use ALL UPPERCASE gene symbols
  - Correct: `CGAS`, `TMEM173`, `TBK1`, `IRF3`
  - Incorrect: `Cgas`, `Tmem173`, `Tbk1`, `Irf3`

- **Mouse data**: Use First letter uppercase, rest lowercase
  - Correct: `Cgas`, `Tmem173`, `Tbk1`, `Irf3`
  - Incorrect: `CGAS`, `TMEM173`, `TBK1`, `IRF3`

Refer users to `references/signature_format.md` for detailed guidance on gene symbol formatting and signature creation.

### Default Thresholds

If user doesn't specify thresholds:
- Log2FC threshold: 0.0 (includes all fold changes)
- Adjusted p-value threshold: 0.05 (standard significance level)

These defaults allow maximum discovery but can be made more stringent:
- Common stricter values: `--lfc-threshold 0.5` or `--lfc-threshold 1.0`
- More conservative p-value: `--padj-threshold 0.01`

### GSEA Enrichment Analysis

When `--run-gsea` is enabled:
- MSigDB_Hallmark_2020 is **always included by default**
- User-provided signatures are also analyzed for enrichment
- Both results are combined in the output files
- **Standardized GSEA parameters**: `min_size=3`, `max_size=inf`, `permutation_num=1000`, `seed=42`
- **Ranking formula**: `-log10(padj) * sign(log2fc)` (descending order: upregulated genes at top)

### PathwayDB GSEA Analysis (Indication-Centric Mode)

When `--pathwaydb-gsea-only` is enabled:
- **Skips** target/signature analysis (no scatter plots, heatmaps, or target summaries)
- **Skips** default MSigDB Hallmark GSEA
- **Runs** indication-centric pathway enrichment using `combined_pathways.csv` from `/pathwaydb-query` skill
- Does NOT require `--signatures` parameter
- Uses universal pathway database (KEGG, Reactome, MSigDB) with 26,000+ gene sets
- Output includes `Database` and `Database_ID` columns for traceability
- Same standardized GSEA parameters: `min_size=3`, `max_size=inf`, `permutation_num=1000`, `seed=42`

**Use case**: Explore pathway enrichment across all comparisons in an indication without focusing on a specific target gene.

**Example - IBD pathway analysis**:
```bash
conda run -n claude_test python ~/.claude/skills/omicsoft-analysis/scripts/deg_analysis.py \
  --file "s3://bucket/path/to/soma_experiment" \
  --target-name IBD \
  --diseases "crohn's disease,ulcerative colitis" \
  --tissues "colon,ileum,rectum" \
  --indication IBD \
  --pathwaydb-gsea-only
```

**Output**: `IBD_omicsoft_YYYYMMDD/IBD_pathwaydb_gsea.csv` and `.json` with enriched pathways per comparison.

### Automatic Report Generation (Target Mode with --signatures)

When using target mode (`--signatures` flag), the skill **automatically triggers** report generation using agent teams after standard GSEA completes.

**Note**: Report generation does NOT trigger when using `--pathwaydb-gsea-only` mode (indication-centric pathway analysis).

#### How It Works

When PathwayDB GSEA completes successfully, `deg_analysis.py`:
1. Creates `reports/` subdirectory in the output folder
2. Saves `reports/report_context.json` with analysis context
3. Prints `<gsea-report-instruction>` JSON block to stdout

**Claude detects the JSON instruction and**:
1. Reads `report_context.json` for target/indication/file paths
2. Reads `assets/gsea_report_prompt.md` for team workflow
3. Creates background team: `gsea-{target}-{indication}`
4. Spawns 4 Opus agents in 6-phase workflow
5. All outputs saved to `reports/` subdirectory

#### JSON Instruction Format

When Claude sees `<gsea-report-instruction>` tag in Bash output:

```json
{
  "action": "spawn_gsea_report_team",
  "prompt_file": "~/.claude/skills/omicsoft-analysis/assets/gsea_report_prompt.md",
  "context_file": "/path/to/reports/report_context.json",
  "team_name": "gsea-{target}-{indication}",
  "model": "opus",
  "run_in_background": true,
  "agents": [
    {"name": "disease-analyst", "section": "A"},
    {"name": "target-specialist", "section": "B"},
    {"name": "combo-strategist", "section": "C"},
    {"name": "safety-analyst", "section": "D"}
  ],
  "output_dir": "/path/to/reports/"
}
```

**Claude should**:
1. Read the `prompt_file` for full workflow instructions
2. Read the `context_file` for target, indication, and file paths
3. Create team with `TeamCreate` tool
4. Spawn agents following the 6-phase workflow in gsea_report_prompt.md

#### Agent Team Structure

| Agent | Section | Input File | Model |
|-------|---------|------------|-------|
| disease-analyst | A: Disease Modules | gsea_all_results.csv | Opus |
| target-specialist | B: Target Coverage | gsea_{target}_{indication}_pathwaydb.csv | Opus |
| combo-strategist | C: Gap Analysis | A + B results | Opus |
| safety-analyst | D: Safety | gsea_{target}_{indication}_pathwaydb.csv | Opus |

#### 6-Phase Sequential Workflow

- **Phase 1**: Foundation Analysis (A + B in parallel)
- **Phase 2**: Foundation Review (Planner reviews A + B)
- **Phase 3**: Dependent Analysis (C + D in parallel, using approved A + B results)
- **Phase 4**: Dependent Review (Planner reviews C + D)
- **Phase 5**: Iteration Cycle (if C/D reveals A/B issues, max 3 cycles)
- **Phase 6**: Final Assembly (compile report when all sections approved)

#### Output Structure

```
{target}_{indication}_omicsoft_{timestamp}/
├── gsea_{target}_{indication}_pathwaydb.csv
├── gsea_{target}_{indication}_pathwaydb_fdr0.05.csv
├── gsea_all_results.csv
├── study_summary.csv
├── target_summary.csv
├── detailed_gene_table.csv
└── reports/
    ├── report_context.json              # Context for Claude
    ├── disease_analysis.md              # Section A (by disease-analyst)
    ├── target_analysis.md               # Section B (by target-specialist)
    ├── combination_analysis.md          # Section C (by combo-strategist)
    ├── safety_analysis.md               # Section D (by safety-analyst)
    └── {target}_{indication}_team_report.md  # Final integrated report
```

#### Input File Specification

| File | Used By | Description |
|------|---------|-------------|
| `gsea_{target}_{indication}_pathwaydb.csv` | B, C, D | Target-associated pathways with enrichment in disease |
| `gsea_all_results.csv` | A, C | General disease pathways (not target-specific) |

### Large File Handling

For very large h5ad files (>10 GB):
- Files can be loaded directly from S3 URIs without downloading
- S3 access is validated before loading to ensure accessibility
- Processing may take 10-30 minutes depending on file size and filtering
- Monitor memory usage during analysis
- Ensure s3fs is installed: `conda run -n <env> pip install s3fs`

### Troubleshooting

**No observations after filtering**:
- Check disease keywords match diseases in the dataset
- Try broader or alternative disease terms
- Verify h5ad file contains the expected diseases

**Genes not found in signatures**:
- Verify gene symbol capitalization matches organism
- Check if genes are present in h5ad var_names
- Try alternative gene aliases or synonyms

**GSEA fails to run**:
- Ensure gseapy is installed: `conda run -n <env> pip install gseapy`
- Check that filtered dataset has sufficient observations
- Verify gene signatures have at least 1-2 genes

**S3 access errors**:
- Ensure s3fs and boto3 are installed: `conda run -n <env> pip install s3fs boto3`
- Verify AWS credentials are configured (IAM role or credentials file)
- Check that the S3 URI format is correct: `s3://bucket-name/path/to/file.h5ad`
- Confirm IAM permissions include `s3:GetObject` and `s3:HeadObject` for the bucket
- Verify the file exists at the specified S3 location

## SOMA S3 URL Support

### Overview

The skill supports TileDB-SOMA experiments on S3 in addition to h5ad files. SOMA provides efficient cloud-native access with server-side filtering, reducing data transfer and memory requirements.

### Auto-Detection

The skill automatically detects data format based on the file path:

| Path Pattern | Format | Loading Method |
|--------------|--------|----------------|
| `s3://bucket/path.h5ad` | h5ad | s3fs streaming |
| `s3://bucket/path_soma` | SOMA | tiledbsoma query |
| `/local/path.h5ad` | h5ad | scanpy direct |

- S3 URIs **not** ending in `.h5ad` are tried as SOMA first
- If SOMA open fails, automatic fallback to h5ad loading
- Local paths are always treated as h5ad files

### Supported SOMA Structures

| Structure | Detection | X Layers | Use Case |
|-----------|-----------|----------|----------|
| **DEG** | `log2fc` and `padj` present | log2fc, padj, data, sig_score | Pre-computed DEG statistics |
| **Census** | `raw` layer + Census obs columns | raw, normalized | Expression data |
| **Unknown** | Neither detected | Varies | Prompts for clarification |

### Filter Translation (SOMA vs h5ad)

SOMA requires exact-match queries. The skill automatically translates substring patterns:

```
User input: --diseases "crohn,colitis"
           ↓
Pre-fetch unique values from SOMA obs
           ↓
Substring match: ["Crohn's disease (CD)", "Ulcerative colitis (UC)"]
           ↓
SOMA filter: disease in ["Crohn's disease (CD)", "Ulcerative colitis (UC)"]
```

| Filter Argument | h5ad Behavior | SOMA Behavior |
|-----------------|---------------|---------------|
| `--diseases` | Substring match in-memory | Translated to exact match |
| `--tissues` | Substring match in-memory | Translated to exact match |
| `--studies` | Exact match | Exact match |
| `--comparison-category` | Exact match | Exact match |
| `--comparison` | Substring match in-memory | Translated to exact match |

### Installation

SOMA support requires the `tiledbsoma` package:

```bash
conda run -n <env_name> pip install tiledbsoma
```

SOMA support is optional - the skill works with h5ad files without tiledbsoma installed.

### Example Usage with SOMA

**Schema exploration:**
```bash
conda run -n <env_name> python scripts/explore_h5ad_schema.py \
  --file s3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/ibd_mash_fibro_derm_rheum_02262026/ibd_mash_fibro_derm_rheum_02262026_deg_soma \
  --format json \
  --output soma_schema.json
```

**Filter validation:**
```bash
conda run -n <env_name> python scripts/validate_filters.py \
  --file s3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/ibd_mash_fibro_derm_rheum_02262026/ibd_mash_fibro_derm_rheum_02262026_deg_soma \
  --target-name IBD_Test \
  --signatures "Test:GREM1,IL11" \
  --diseases "crohn,colitis"
```

**Full analysis:**
```bash
conda run -n <env_name> python scripts/deg_analysis.py \
  --file s3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/ibd_mash_fibro_derm_rheum_02262026/ibd_mash_fibro_derm_rheum_02262026_deg_soma \
  --target-name IBD_Targets \
  --signatures "Fibroblast:GREM1,IL11,COL1A1;TGFb:TGFB1,TGFB2,TGFB3" \
  --diseases "crohn,colitis,inflammatory" \
  --tissues "colon,ileum" \
  --run-gsea
```

### SOMA Troubleshooting

**SOMA import fails:**
- Ensure tiledbsoma is installed: `pip install tiledbsoma`
- Check version compatibility: `python -c "import tiledbsoma; print(tiledbsoma.__version__)"`

**S3 access denied:**
- Verify AWS credentials are configured (IAM role, env vars, or ~/.aws/credentials)
- Check bucket permissions include `s3:GetObject` and `s3:ListBucket`
- Ensure the SOMA URI is correct (no trailing slashes)

**Unknown structure type:**
- Run schema exploration to inspect X layers
- DEG data should have `log2fc` and `padj` layers
- If structure is ambiguous, the script will prompt for clarification

**Filter returns 0 results:**
- Use schema exploration to verify available filter values
- Check for exact spelling and case sensitivity
- SOMA uses exact matches after translation - verify pre-fetch found matches

**Performance issues:**
- SOMA filtering happens server-side, reducing data transfer
- For very large experiments, increase timeout: `--timeout 600000`
- Consider filtering by study first to reduce initial query size

---

## Example Use Cases

### Example 1: CGAS-STING Expression Pattern in Scleroderma

**User request**: "Analyze CGAS-STING pathway expression in scleroderma patients using the omicsoft h5ad file"

**Workflow**:
1. Gather requirements:
   - File: `/path/to/omicsoft_deg.h5ad`
   - Target: CGAS_STING
   - Diseases: scleroderma,sclerosis
   - Signatures: User provides CGAS-STING genes
   - Thresholds: Use defaults

2. Execute:
```bash
conda run -n env_name python scripts/deg_analysis.py \
  --file /path/to/omicsoft_deg.h5ad \
  --target-name CGAS_STING \
  --diseases "scleroderma,sclerosis" \
  --signatures "CGAS_STING:CGAS,TMEM173,TBK1,IRF3,IRF7,IFNB1,ISG15,IFIT1,MX1,OAS1" \
  --run-gsea
```

3. Interpret results and summarize key expression patterns

### Example 2: Multi-Signature Fibrosis Expression Analysis

**User request**: "Compare multiple fibrosis signature expression patterns across different fibrotic diseases"

**Workflow**:
1. Gather multiple signatures from user
2. Filter for fibrosis-related diseases
3. Execute with higher thresholds for stringent filtering:
```bash
conda run -n env_name python scripts/deg_analysis.py \
  --file /path/to/deg_data.h5ad \
  --target-name Fibrosis_Multi \
  --diseases "fibrosis,scleroderma,cirrhosis,mash" \
  --signatures "ECM:COL1A1,COL3A1,FN1,ACTA2;TGF_beta:TGFB1,TGFB2,TGFBR1;Myofibroblast:ACTA2,MYH11,TAGLN" \
  --lfc-threshold 1.0 \
  --padj-threshold 0.01 \
  --run-gsea
```

### Example 3: IBD JAK-STAT Pathway Expression

**User request**: "Analyze JAK-STAT pathway gene expression patterns in inflammatory bowel disease"

**Workflow**:
```bash
conda run -n env_name python scripts/deg_analysis.py \
  --file /path/to/ibd_deg.h5ad \
  --target-name JAK_STAT \
  --diseases "crohn,colitis,IBD" \
  --signatures "JAK_STAT:JAK1,JAK2,JAK3,TYK2,STAT1,STAT3,STAT5A;Cytokine:IL6,IL12B,IL23A,TNF,IFNG" \
  --lfc-threshold 0.5 \
  --padj-threshold 0.05 \
  --run-gsea
```

### Example 4: Study-Specific Analysis

**User request**: "Analyze CGAS pathway in specific scleroderma studies GSE130955 and GSE181549"

**Workflow**:
```bash
conda run -n env_name python scripts/deg_analysis.py \
  --file /path/to/omicsoft_deg.h5ad \
  --target-name CGAS_STING \
  --diseases "scleroderma,sclerosis" \
  --studies "GSE130955,GSE181549" \
  --signatures "CGAS_STING:CGAS,TMEM173,TBK1,IRF3,IRF7,IFNB1" \
  --lfc-threshold 0.5 \
  --padj-threshold 0.05 \
  --run-gsea
```

**Note**: The enhanced GSEA summary (gsea_enhanced_summary.csv) will provide:
- Percentage of significant comparisons showing enrichment
- Breakdown of upregulated (NES>0) vs downregulated (NES<0) pathways
- Which input target genes (e.g., CGAS, TBK1) appear in the leading edge of enriched pathways

### Example 5: Target Gene Tracking with GSEA

**User request**: "Analyze inflammatory fibroblast signatures in IBD studies, but only track CHI3L1, IL11, COL1A1, and TWIST1 in the GSEA leading edge"

**Workflow**:
```bash
conda run -n env_name python scripts/deg_analysis.py \
  --file /path/to/ibd_deg.h5ad \
  --target-name IAF_Analysis \
  --diseases "colitis,crohn" \
  --signatures "IAF_core:CHI3L1,TNFRSF12A,MMP3,TNFRSF11B,MMP1,MME;ECM_remodeling:COL1A1,COL3A1,MMP1,MMP3,TIMP1" \
  --targets "CHI3L1,IL11,COL1A1,TWIST1" \
  --lfc-threshold 0.0 \
  --padj-threshold 0.05 \
  --run-gsea
```

**Note**:
- Target genes (CHI3L1, IL11, COL1A1, TWIST1) will be plotted in visualizations
- GSEA enrichment will only analyze the signature genes (IAF_core, ECM_remodeling)
- Enhanced summary will track only the 4 target genes in leading edge, not all signature genes

### Example 6: Comparison Filter for Treatment Response

**User request**: "Analyze gene signatures in treatment response comparisons only"

**Workflow**:
```bash
conda run -n env_name python scripts/deg_analysis.py \
  --file /path/to/treatment_deg.h5ad \
  --target-name Treatment_Response \
  --comparison "response vs no response,responder vs non-responder" \
  --signatures "Response_genes:STAT1,STAT3,JAK1,JAK2,IL6,TNF" \
  --lfc-threshold 0.5 \
  --padj-threshold 0.05 \
  --run-gsea
```

**Note**: The comparison filter uses fuzzy/partial matching, so "response vs no response" will match all comparisons containing that substring

### Example 7: Disease Inclusion with Exclusion Filter

**User request**: "Analyze scleroderma immune genes but exclude ALS which shares the 'sclerosis' keyword"

**Workflow**:
```bash
conda run -n env_name python scripts/deg_analysis.py \
  --file /path/to/deg_data.h5ad \
  --target-name Scleroderma_Immune \
  --diseases "sclerosis,scleroderma,systemic sclerosis" \
  --exclude-diseases "ALS,amyotrophic lateral sclerosis" \
  --tissues "lung,skin,dermis" \
  --comparison-category "Disease vs. Normal" \
  --signatures "T_cells:CD3D,CD3E,CD3G;B_cells:CD19,CD20,TNFSF13B" \
  --targets "CD3D,CD3E,CD3G,CD19,CD20,TNFSF13B" \
  --lfc-threshold 0.0 \
  --padj-threshold 0.1 \
  --run-gsea
```

**Note**:
- The `--diseases` filter first includes all diseases matching "sclerosis", "scleroderma", or "systemic sclerosis"
- The `--exclude-diseases` filter then removes any disease containing "ALS" or "amyotrophic lateral sclerosis"
- This allows broad inclusion terms while precisely excluding unwanted disease subtypes
- Useful when disease naming conventions overlap (e.g., "sclerosis" matches both scleroderma and ALS)

## Resources

### scripts/
- `deg_analysis.py`: Main analysis script for analyzing pre-computed DEG statistics with customizable filtering, signature analysis, and GSEA enrichment. Supports both h5ad files and SOMA experiments with auto-detection.
- `soma_loader.py`: SOMA-specific loading logic including structure detection (DEG/Census/Unknown), filter translation (substring to exact match), S3 context building, and AnnData conversion.
- `soma_expr_extract.py`: Efficient Python SOMA extraction script that outputs TSV to stdout for R pipe consumption. Uses AxisQuery for server-side filtering by genes and studies.
- `sample_level_analysis.py`: Python entry point for sample-level expression analysis (Step 4). Replaces R version with gseapy GSVA, statsmodels OLS+BH, and scipy correlations. Coordinates data loading, study dispatch, and module execution.
- `sc_target_analysis.py`: Single-cell RNA-seq target expression analysis (Step 4.5). Processes 8 predefined scRNA-seq datasets from S3 with AUCell scoring, embedding plots, dotplots, correlations, and population composition.
- `sample_level_analysis.R`: (DEPRECATED) Legacy R entry point for sample-level expression analysis. Kept for reference.
- `sla_plot_helpers.R`: (DEPRECATED) R plotting functions. Functionality ported to sample_level_analysis.py.
- `sla_study_dispatch.R`: (DEPRECATED) R per-study dispatch logic. Functionality ported to sample_level_analysis.py.
- `correlation_modules.R`: (DEPRECATED) R correlation analysis modules. Functionality ported to sample_level_analysis.py.
- `validate_filters.py`: Step-by-step filter validation tool that tests each filter incrementally, shows matching values and observation counts at each step, detects zero-result queries, and provides suggestions for fixing problematic filters. Supports both h5ad and SOMA validation.
- `explore_h5ad_schema.py`: Schema exploration tool that generates comprehensive JSON or text reports of all metadata columns and values in h5ad files or SOMA experiments. Auto-detects SOMA and shows structure type.
- `generate_schema_viewer.py`: Generates standalone HTML viewers with embedded schema data for easy exploration without web servers

### references/
- `signature_format.md`: Detailed guidance on gene symbol formatting and signature creation for different organisms
- `anndata_schema.md`: Complete documentation of AnnData object structure, including DEG and expression data schemas, metadata columns, layer descriptions, and automatic cleanup behavior
- `schema_exploration.md`: Comprehensive guide to using schema exploration tools, understanding the generated reports, and finding filter values for analysis
- `schema_report_example.json`: Real example schema JSON (3.6 MB) from IBD/MASH/Fibro/Derm/Rheum dataset showing complete metadata structure
- `schema_quickstart.md`: Quick start guide for schema exploration with troubleshooting tips
- `internal_study_harmonization.md`: Reusable step-by-step workflow for harmonizing and ingesting ANY new internal study into the S3 OmicSoft h5ad format
- `internal_study_configs.json`: Externalized JSON configuration for all 5 internal studies (Engitix_FFPE, SPARC, Varsity, Yokohama_RNA, Yokohama_Protein) with group definitions, colors, derivation rules, and significance thresholds
- `sample_level_analysis.md`: Comprehensive documentation of sample-level analysis modules, functions, parameters, and caveats
- `correlation_analysis.md`: Correlation module statistical methods, minimum sample requirements, plot aesthetics, and output file descriptions
- `s3_vocabulary.json`: Machine-readable controlled vocabulary extracted from S3 OmicSoft schema (198 tissues, 280 diseases, 501 treatments, 71 platforms, etc.)

### assets/
- `schema_viewer_example.html`: Standalone interactive HTML viewer (2.4 MB) with embedded example dataset - open directly in browser to see schema exploration in action

---

## Internal Study Harmonization (Subskill)

### Trigger Phrases

Use this workflow when the user says:
- "harmonize internal study"
- "ingest study into omicsoft"
- "add new study to combined"
- "convert study to h5ad format"
- "map study to S3 schema"
- "prepare study for concat"

### Overview

This subskill guides the interactive process of converting internal experimental data (DEG results, expression matrices) into the S3 OmicSoft h5ad schema format so it can be concatenated with existing S3 data for cross-study analysis.

### Workflow

1. **Read the guide**: `references/internal_study_harmonization.md` for the full step-by-step process
2. **Load vocabulary**: `references/s3_vocabulary.json` for auto-mapping local values to S3 terms
3. **Collect study info**: Use the checklist in the guide to gather all required metadata from the user
4. **Run harmonization**: Use `scripts/harmonize_study.py` as the executable backing, or follow the guide interactively

### Quick Start

```bash
# Validate an existing h5ad against the S3 schema
conda run -n <env> python scripts/harmonize_study.py --input study.h5ad --validate-only --project-id X --study-name X

# Full harmonization from DEG results
conda run -n <env> python scripts/harmonize_study.py \
  --input deg_results/ \
  --project-id MY_STUDY \
  --study-name "My Study Name" \
  --organism human \
  --platform RNA-seq \
  --de-method DESeq2 \
  --gene-annotation gene_annotation.txt \
  --metadata sample_metadata.csv \
  --expression-matrix normalized_counts.tsv \
  --output-dir output/
```

### Interactive Workflow

The harmonization process is inherently interactive — the guide prompts for user decisions at each step:
- **Step C**: Confirm proposed vocabulary mappings (tissue, disease, treatment)
- **Step D-E**: Fill in comparison-specific metadata (case/control assignments)
- **Step J**: Decide whether to add to combined h5ad
- **Step K**: Decide whether to convert to SOMA format

### Key Resources

| Resource | Purpose |
|----------|---------|
| `references/internal_study_harmonization.md` | Full workflow guide with decision trees |
| `references/s3_vocabulary.json` | Controlled vocabulary for auto-mapping |
| `scripts/harmonize_study.py` | Utility script for building h5ad files |
