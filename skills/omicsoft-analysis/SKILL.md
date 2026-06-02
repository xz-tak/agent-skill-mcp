---
name: omicsoft-analysis
description: Analyze pre-computed differential gene expression statistics from Omicsoft h5ad files to investigate target gene expression patterns and signature enrichment. This skill should be used when users request to analyze bulk DEG data (log2fc, padj) from h5ad files, filter for specific diseases or studies, examine target genes or pathways with custom signatures, perform GSEA enrichment analysis with enhanced summaries including leading edge annotation, and generate interactive visualizations.
---

# Omicsoft DEG Analysis

Analyze pre-computed differential gene expression statistics from Omicsoft bulk transcriptomics data stored in h5ad format to investigate target gene expression patterns and signature enrichment.

## Overview

Analyze pre-computed DEG statistics (log2fc and padj) to filter h5ad files for specific diseases and/or studies, investigate custom gene signatures and target pathway expression patterns, generate interactive visualizations (scatter plots and hierarchically-clustered heatmaps), perform GSEA signature enrichment analysis with enhanced summaries (including significance counts, NES direction separation, and leading edge annotation), and export detailed results with configurable significance thresholds.

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

3. **Gene signatures**: User-defined gene sets in the format:
   - Format: `"SignatureName:Gene1,Gene2,Gene3;AnotherSignature:GeneA,GeneB"`
   - Example: `"CGAS_STING:CGAS,TMEM173,TBK1,IRF3;TGFB:TGFB1,TGFB2,TGFB3"`
   - **IMPORTANT**: Gene symbols must match organism:
     - Human: ALL UPPERCASE (e.g., `CGAS`, `TMEM173`)
     - Mouse: First letter uppercase, rest lowercase (e.g., `Cgas`, `Tmem173`)
   - Refer to `references/signature_format.md` for detailed formatting guide

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
   - `--run-gsea`: Enable GSEA enrichment analysis (always includes MSigDB_Hallmark_2020 by default)
   - `--output-dir`: Output directory name (default: deg_results)

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
  [--run-gsea]
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
  --padj-threshold 0.1 \
  --run-gsea
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
  --output-dir cgas_results \
  --run-gsea
```

### Step 4: Interpret Results

The analysis generates multiple output files in the specified output directory:

1. **study_summary.csv**: List of studies included after filtering
   - Columns: study, tissue, disease, disease_category

2. **signature_summary.csv**: Target gene expression patterns for signatures
   - Shows which genes are up/down-regulated per tissue, disease, and comparison category
   - Grouped by: Gene, tissue, disease_category, comparison_category, sig (up/dn)
   - Columns: Gene, tissue, disease_category, comparison_category, sig (up/dn), count, study, total_count, threshold

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

5. **detailed_gene_table.csv**: Comprehensive gene-level data export
   - Long-format table with all metadata for signature and target genes
   - Each row represents one gene in one comparison
   - Contains all columns: log2fc, padj, tissue, disease, comparison, study, etc.
   - Useful for custom downstream analysis and filtering

6. **gsea_all_results.csv**: Complete GSEA enrichment results (if --run-gsea used)
   - All signature enrichment results across studies
   - Includes both custom signatures and MSigDB_Hallmark_2020
   - Contains comparison_category for tracking enrichment context

7. **gsea_enhanced_summary.csv**: Enhanced GSEA summary with detailed statistics (if --run-gsea used)
   - **Groupby format**: Each Term has separate rows for each Comparison_Category
   - Allows comparison of enrichment patterns between "Disease vs. Normal" and "Responder vs. Non-Responder"
   - Number of significant comparisons (FDR < threshold) per term-category combination
   - Studies with positive NES (upregulated) and their NES values
   - Studies with negative NES (downregulated) and their NES values
   - Input target genes found in leading edge (tracks only --targets genes if provided, otherwise all signature genes)
   - Columns include: Term, Comparison_Category, N_Significant_Comparisons, Total_Comparisons, Percent_Significant, N_Positive_NES, N_Negative_NES, Positive_NES_Studies, Negative_NES_Studies, Input_Targets_In_Leading_Edge, N_Input_Targets_In_LE
   - Note: Target genes (--targets) are plotted but NOT included in GSEA enrichment analysis

8. **gsea_score_pathway_YYYYMMDD.html**: GSEA visualization (if --run-gsea used)
   - Interactive scatter plot of normalized enrichment scores (NES)
   - Bubble size indicates significance

### Step 5: Summarize Findings

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
- GSEA parameters: `min_size=1`, `max_size=8000`, `permutation_num=1000`
- Uses pre-ranked GSEA based on sig_score (log2fc * -log(padj))

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
- `deg_analysis.py`: Main analysis script for analyzing pre-computed DEG statistics with customizable filtering, signature analysis, and GSEA enrichment
- `validate_filters.py`: Step-by-step filter validation tool that tests each filter incrementally, shows matching values and observation counts at each step, detects zero-result queries, and provides suggestions for fixing problematic filters before running full analysis
- `explore_h5ad_schema.py`: Schema exploration tool that generates comprehensive JSON or text reports of all metadata columns and values in h5ad files
- `generate_schema_viewer.py`: Generates standalone HTML viewers with embedded schema data for easy exploration without web servers

### references/
- `signature_format.md`: Detailed guidance on gene symbol formatting and signature creation for different organisms
- `anndata_schema.md`: Complete documentation of AnnData object structure, including DEG and expression data schemas, metadata columns, layer descriptions, and automatic cleanup behavior
- `schema_exploration.md`: Comprehensive guide to using schema exploration tools, understanding the generated reports, and finding filter values for analysis
- `schema_report_example.json`: Real example schema JSON (3.6 MB) from IBD/MASH/Fibro/Derm/Rheum dataset showing complete metadata structure
- `schema_quickstart.md`: Quick start guide for schema exploration with troubleshooting tips

### assets/
- `schema_viewer_example.html`: Standalone interactive HTML viewer (2.4 MB) with embedded example dataset - open directly in browser to see schema exploration in action
