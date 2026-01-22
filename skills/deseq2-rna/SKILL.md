---
name: deseq2-rna
description: Interactive bulk RNA-seq differential expression analysis using DESeq2 in R. This skill should be used when users request RNA-seq analysis, differential expression analysis, DEG analysis, or bulk transcriptomics analysis with raw count data and metadata files.
---

# DESeq2 RNA-seq Analysis Skill

Interactive workflow for bulk RNA-seq differential expression analysis using DESeq2.

## When to Use This Skill

This skill applies when users:
- Request bulk RNA-seq differential expression analysis
- Provide raw count data and sample metadata files
- Ask for DEG (differentially expressed genes) analysis
- Need pathway enrichment analysis (GSEA, g:Profiler)
- Want volcano plots, PCA, or DE summary visualizations

## Input Requirements

### Counts File
- **Format**: CSV or TSV with genes as rows, samples as columns
- **Values**: Non-negative integers (raw counts, not normalized)
- **Gene IDs**: Row names as gene symbols or ENSEMBL IDs (with optional Gene.name column)

### Metadata File
- **Format**: TSV with samples as rows
- **Required**: Sample names matching counts column names
- **Columns**: At least one categorical column for comparison (e.g., Treatment, Condition)

## Interactive Workflow (10 Steps)

### Step 1: Check R Environment & Input Files

Run `scripts/check_r_env.R` to validate:
- Conda r_env exists with required packages
- Counts file has valid format (non-negative integers)
- Metadata file matches counts columns
- Species-specific annotation packages available

**Required R Packages:**
- Core: DESeq2, clusterProfiler, msigdbr, gprofiler2, enrichplot, PCAtools, ggplot2, pheatmap, glue, dplyr, GSVA, BiocParallel, limma, tidyr, tibble, ggrepel
- Species: org.Hs.eg.db (human), org.Mm.eg.db (mouse), org.Rn.eg.db (rat)

### Step 2: Summarize Metadata

Present metadata summary table to user:

```
| Column | Type | Unique Values | Sample Values | Rationale |
|--------|------|---------------|---------------|-----------|
```

For each column, show:
- Data type (categorical/numeric)
- Number of unique values
- Example values (first 3-5)
- Suggested role (comparison factor, covariate, identifier)

### Step 3: Align Comparison & Covariate Columns

Use AskUserQuestion to collect:
- **Primary comparison column**: Main factor for DESeq2 contrasts
- **Covariate columns** (optional): Batch effects to include in design

### Step 4: Align Design Formula & Test Collinearity

Construct design formula: `~Covariate1 + Covariate2 + ComparisonColumn`

Test design matrix for collinearity:
```r
model_matrix <- model.matrix(design_formula, data = coldata)
if (qr(model_matrix)$rank < ncol(model_matrix)) {
  # Report error and ask for new column selection
}
```

If collinear, report which columns are problematic and ask user to revise.

### Step 5: Align Comparison Pairs & Reference

Use AskUserQuestion to collect:
- **Reference level**: Baseline group (e.g., "Control", "Non-treatment")
- **Contrast method**: Choose based on comparison complexity

#### Option A: Simple Contrast (Default for pairwise comparisons)
Standard DESeq2 contrast vector `c(factor, numerator, denominator)`:
```json
"comparisons": "all_pairwise"
// OR
"comparisons": [["TreatmentA", "Control"], ["TreatmentB", "Control"]]
```
**Use for:** Direct A vs B comparisons within the same factor level.

#### Option B: List Contrast (REQUIRED for interaction/complex comparisons)
DESeq2 list contrast with `listValues` for proper statistical inference:
```json
"advanced_comparisons": [
  {
    "name": "DrugEffect_in_Stimulated",
    "type": "list",
    "contrast": ["Treatment_DrugA_vs_Control", "Treatment_DMSO_vs_Control"],
    "listValues": [1, -1]
  },
  {
    "name": "Interaction_Drug_x_Stimulus",
    "type": "list",
    "contrast": [["Treatment_DrugA_vs_Control", "Stimulus_Stim_vs_None"], "Treatment_DMSO_vs_Control"],
    "listValues": [1, -2]
  }
]
```

**When to use List Contrast (ENFORCED):**
- Interaction terms / difference-in-differences
- Comparing drug effects across different stimulus conditions
- Any contrast requiring weighted coefficient combinations
- Reversal analyses (e.g., does drug reverse stimulus effect?)

**listValues patterns:**
| Pattern | Formula | Use Case |
|---------|---------|----------|
| `[1, -1]` | coef1 - coef2 | Compare two model coefficients |
| `[1, -2]` | (coef1+coef2) - 2×coef3 | Interaction/reversal effects |

**Coefficient names:** Available via `resultsNames(dds)` after model fitting.
Format: `{Factor}_{Level}_vs_{Reference}` (e.g., `Treatment_DrugA_vs_Control`)

### Step 6: Align Species & Signature Options

Use AskUserQuestion to collect:
- **Species**: hsapiens, mmusculus, or rnorvegicus
- **MSigDB collections**: Comma-separated (H, C2, C5, C7, etc.)
- **Custom signatures** (optional): File path or "none"

Signature file formats supported (auto-detected):
- GMT: `SET_NAME<tab>DESC<tab>GENE1<tab>GENE2...`
- Text: `>SET_NAME` followed by one gene per line

### Step 7: Align FDR & Log2FC Cutoffs

Use AskUserQuestion to collect:
- **padj threshold**: Default 0.05
- **|log2FC| threshold**: Default 1.0

### Step 8: Plan Summary & Confirmation

Present complete analysis plan:
```
=== DESeq2 Analysis Plan ===
Counts file: {counts_path}
Metadata file: {metadata_path}
Samples: {n_samples}
Genes: {n_genes}
Species: {species}

Design: ~{design_formula}
Reference: {reference_level}
Comparisons: {n_comparisons} pairs

Signatures: MSigDB {collections}
Custom signatures: {n_custom} sets
Cutoffs: padj < {padj}, |log2FC| > {lfc}

Output directory: {output_dir}
```

Confirm with user before proceeding. After confirmation (or any changes):
- Save the plan to `{output_dir}/PLAN.md` for documentation
- Include all parameters, comparison list, and timestamp

### Step 9: Run Analysis

Execute `scripts/deseq2_analysis.R` with parameters via command-line or JSON config.

**The script will automatically:**
1. Run DESeq2 differential expression analysis
2. Generate all plots and enrichment analyses
3. **Trigger draft report generation** upon successful completion

**Outputs generated:**
- `figures/`: PCA plots, volcano plots, GSEA dotplots, g:Profiler plots, DE summary barplot
- `deg/`: Per-comparison DEG tables, enrichment results, combined statistics
- `{prefix}_analysis_data.rds`: All R objects for downstream analysis
- `{prefix}_Draft_Report.md`: Draft report with data summaries and interpretation placeholders

### Step 10: Generate Biological Interpretations

---
#### **⚠️ ULTRATHINK REQUIRED ⚠️**

**This step REQUIRES extended thinking (ultrathink) for deep biological analysis.**

Claude MUST use extended thinking to provide expert-level interpretations.
---

#### Claude's Role: Professional Biologist

**For this step, Claude assumes the role of a professional biologist with broad expertise in:**
- **Molecular Biology**: Gene regulation, signaling pathways, protein function, transcriptional networks
- **Disease Biology**: Broad expertise across therapeutic areas (oncology, immunology, fibrosis, metabolic, neurological, cardiovascular, infectious disease, rare diseases, etc.)
- **Computational Biology**: Transcriptomics, pathway analysis, systems biology, biomarker discovery

**Context-Specific Focus:**
While Claude has broad biological expertise, interpretations MUST be focused on:
- The specific target/pathway being studied (infer from comparison names, gene lists, metadata)
- The relevant disease biology and therapeutic context
- The cell type/model system used in the experiment
- The mechanism of action relevant to the treatments tested

**Interpretation Standards:**
- Provide mechanistic explanations, not just descriptions
- Connect findings to established biological knowledge specific to the study context
- Identify clinically or therapeutically relevant insights for the indication under study
- Generate testable hypotheses grounded in the data and relevant to the disease/target

---

#### **🔬 ULTRATHINK: Deep Biological Analysis Protocol**

**Before writing ANY interpretation, Claude MUST:**

1. **ULTRATHINK - Context Identification Phase:**
   - Identify the target/pathway being studied from file names, comparison names, and top genes
   - Determine the disease/indication context from metadata and experimental design
   - Recognize the cell type/model system and its relevance
   - Understand what biological question each comparison addresses

2. **ULTRATHINK - Extended Analysis Phase:**
   - Deeply analyze the DEG patterns, pathway enrichments, and top genes
   - Consider multiple biological mechanisms that could explain the observations
   - Cross-reference with known biology of the SPECIFIC target/pathway under study
   - Identify unexpected findings that warrant special attention
   - Formulate mechanistic hypotheses based on the data AND disease context

3. **WebFetch - Adaptive Knowledge Retrieval:**
   - Use WebFetch tool to retrieve relevant biological context from the MOST APPROPRIATE sources based on task context
   - Automatically select relevant databases/resources which may include but are not limited to:
     - **Literature**: PubMed, Google Scholar, preprint servers for recent findings
     - **Gene/Protein**: UniProt, GeneCards, NCBI Gene, Ensembl
     - **Pathways**: KEGG, Reactome, WikiPathways, STRING
     - **Disease**: OMIM, DisGeNET, ClinVar, GWAS Catalog
     - **Drug/Target**: DrugBank, ChEMBL, OpenTargets
     - **Cell/Tissue**: Human Protein Atlas, GTEx, CellMarker
     - **Any other relevant scientific databases** based on the specific study context
   - Focus retrieval on the SPECIFIC target, disease, and mechanisms relevant to this study
   - Prioritize recent publications and authoritative sources
   - Integrate retrieved knowledge into context-specific interpretations

---

#### Workflow Steps

**1. Read key data files:**
   - `deg/{prefix}_DE_lfc*_count.txt` - DEG counts summary
   - `deg/{prefix}_summstats_all.txt` - Combined DEG statistics (for top genes)
   - `deg/{prefix}_gsea_all.txt` - Combined GSEA results (GO:BP + MSigDB + Custom) with `source` column
   - `deg/{prefix}_ora_all.txt` - Combined ORA results (g:Profiler + MSigDB + Custom) with `source` column
   - `figures/*.png` - Volcano plots, PCA, GSEA/ORA dotplots by source (view for visual patterns)

**2. For each comparison, generate ULTRATHINK interpretations:**

   **🧬 ULTRATHINK Required for Each Section:**

   - **Biological Context**: What biological question does this comparison address?
   - **DEG Pattern Analysis**:
     - What does the up/down gene ratio suggest about cellular state?
     - Are there signs of activation, suppression, stress, or differentiation?
   - **Top Upregulated Genes** (with ULTRATHINK):
     - Research each gene's function using WebFetch if needed
     - Explain biological roles and relevance to the experimental context
     - Identify gene families or functional clusters
   - **Top Downregulated Genes** (with ULTRATHINK):
     - Same deep analysis as upregulated genes
     - Consider what suppression of these pathways means
   - **Pathway Enrichment Insights** (with ULTRATHINK):
     - Explain biological significance of each enriched pathway
     - Connect pathways to the experimental treatment/condition
     - Identify pathway crosstalk and regulatory relationships
   - **Mechanistic Summary**:
     - Synthesize all findings into a coherent biological narrative
     - Propose molecular mechanisms explaining the observations

**3. For the Discussion section, generate ULTRATHINK analysis:**

   **🔬 ULTRATHINK Required for Discussion:**

   - **Summary of Major Findings**:
     - Rank comparisons by biological significance
     - Identify the most important discoveries
   - **Treatment Effects Analysis**:
     - Characterize each treatment's transcriptional signature
     - Compare and contrast treatment effects
   - **Pathway Convergence and Divergence**:
     - Which pathways are shared across treatments?
     - Which are treatment-specific?
     - What does this reveal about mechanism of action?
   - **Biological Implications**:
     - What do findings mean for the biological system under study?
     - Relevance to the disease/indication and therapeutic implications
   - **Testable Hypotheses**:
     - Generate specific, experimentally testable predictions
     - Suggest validation experiments
   - **Limitations**:
     - Acknowledge study constraints
     - Suggest improvements for future studies

**4. Update the draft report:**
   - Replace all `<!-- CLAUDE_INTERPRETATION_START -->` placeholders
   - Ensure ULTRATHINK depth is evident in all interpretations
   - Save as `{prefix}_Analysis_Report.md`

---

#### Quality Checklist for Interpretations

Before finalizing, verify each interpretation includes:
- [ ] Mechanistic explanation (not just gene/pathway lists)
- [ ] Connection to known biology (cite specific knowledge)
- [ ] Clinical or therapeutic relevance where applicable
- [ ] Specific, testable hypotheses
- [ ] Integration of multiple data types (DEGs + pathways + plots)

---

**Output:** Final report with comprehensive ULTRATHINK biological interpretations at:
`{output_dir}/{prefix}_Analysis_Report.md`

## Script Usage

### Command-line Arguments
```bash
conda run -n r_env Rscript scripts/deseq2_analysis.R \
  --counts /path/to/counts.txt \
  --metadata /path/to/metadata.txt \
  --design "~Batch + Treatment" \
  --comparison_col Treatment \
  --reference Control \
  --comparisons '[["TreatA","Control"],["TreatB","Control"]]' \
  --msigdb "H,C2" \
  --species hsapiens \
  --padj 0.05 \
  --lfc 1 \
  --output_dir /path/to/output \
  --conda_env r_env
```

### JSON Config File
```bash
conda run -n r_env Rscript scripts/deseq2_analysis.R --config /path/to/config.json
```

Config format:
```json
{
  "counts": "/path/to/counts.txt",
  "metadata": "/path/to/metadata.txt",
  "design": "~Batch + Treatment",
  "comparison_col": "Treatment",
  "reference": "Control",
  "comparisons": [["TreatA", "Control"], ["TreatB", "Control"]],
  "advanced_comparisons": [
    {
      "name": "Interaction_Drug_x_Stimulus",
      "contrast": [["Treatment_DrugA_vs_Control", "Stimulus_Stim_vs_None"], "Treatment_DMSO_vs_Control"],
      "listValues": [1, -2]
    },
    {
      "name": "DrugA_vs_DrugB_in_Stim",
      "contrast": ["Treatment_DrugA_vs_Control", "Treatment_DrugB_vs_Control"],
      "listValues": [1, -1]
    }
  ],
  "msigdb": ["H", "C2"],
  "custom_sigs": null,
  "species": "hsapiens",
  "padj": 0.05,
  "lfc": 1,
  "output_dir": "/path/to/output",
  "conda_env": "r_env"
}
```

**Note:** `advanced_comparisons` uses DESeq2 list contrasts with proper statistical inference.
Coefficient names are available via `resultsNames(dds)` after model fitting.

## Reference Documentation

- `references/workflow_guide.md`: Detailed workflow documentation
- `references/signature_format.md`: Custom signature file format guide

## Output Files

### Figures Directory
| File Pattern | Description |
|--------------|-------------|
| `{prefix}_PCA_basic_*.png` | PCA colored by comparison column |
| `{prefix}_PCA_screeplot.png` | Variance explained per PC |
| `{prefix}_PCA_eigencorplot.png` | PC correlations with metadata + CUSTOM signatures only (excludes MSigDB) |
| `{prefix}_PCA_biplot_{column}.png` | PCA biplots with gene loadings, colored by each categorical metadata column |
| `GSEA_{prefix}_PC1_loadings.png` | Pathways driving PC1 variance (GO:BP enrichment of PC1 loadings) |
| `GSEA_{prefix}_PC2_loadings.png` | Pathways driving PC2 variance (GO:BP enrichment of PC2 loadings) |
| `GSEA_{prefix}_PC3_loadings.png` | Pathways driving PC3 variance (GO:BP enrichment of PC3 loadings) |
| `volcano_{comparison}.png` | Volcano plot per comparison |
| `GSEA_{comparison}_{source}.png` | GSEA enrichment by source (GO_BP, H, C2, CUSTOM, etc.) |
| `ORA_{comparison}_{source}.png` | ORA enrichment by source (GO_BP, KEGG, REAC, H, C2, CUSTOM, etc.) |
| `DE_genes_summary_barplot.png` | DEG counts across comparisons |

### Interactive Figures Directory (figures/interactive/)

Interactive HTML plots are automatically generated alongside static PNG plots. These provide hover tooltips with detailed statistical information.

| File Pattern | Description | Hover Information | Cell Annotation |
|--------------|-------------|-------------------|-----------------|
| `eigencorplot_interactive.html` | PC-metadata correlations | Correlation (r), p-value, FDR | N/A |
| `volcano_{comparison}_interactive.html` | Volcano plots per comparison | Gene symbol, log2FC, p-value, padj, baseMean | N/A |
| `nes_heatmap_{source}_interactive.html` | NES heatmaps by source (CUSTOM, GO:BP, Hallmark) | Signature name, NES, p-value, FDR | `*` significance stars |

**Note:** Interactive HTML plots require a web browser to view. They do NOT work when embedded in PDF, PowerPoint, or static markdown documents. Use static PNG versions for reports and presentations.

### DEG Directory
| File Pattern | Description |
|--------------|-------------|
| `{prefix}_{comparison}.txt` | Full DEG statistics per comparison |
| `{prefix}_{comparison}_gsea.txt` | Combined GSEA results (GO:BP + MSigDB + Custom) |
| `{prefix}_{comparison}_ora.txt` | Combined ORA results (g:Profiler + enricher) |
| `{prefix}_summstats_all.txt` | Combined DEG statistics across all comparisons |
| `{prefix}_gsea_all.txt` | Combined GSEA results with `source` column |
| `{prefix}_ora_all.txt` | Combined ORA results with harmonized columns |
| `{prefix}_PC_gsea.txt` | PC1/2/3 GSEA results with `PC` column (pathways driving variance) |
| `{prefix}_DE_lfc{lfc}padj{padj}_count.txt` | DEG counts summary |

### GSEA Output Columns
| Column | Description |
|--------|-------------|
| ID | Gene set name |
| Description | Gene set description |
| NES | Normalized enrichment score |
| pvalue | p-value |
| p.adjust | Adjusted p-value |
| **source** | Origin: "GO:BP", "H", "C2", "CUSTOM", etc. |
| comparison | Comparison name |

### ORA Output Columns (Harmonized)
| Column | Description |
|--------|-------------|
| term_id | Term/pathway identifier |
| term_name | Pathway/term name |
| p_value | p-value |
| p_adjust | Adjusted p-value (FDR) |
| **source** | Origin: "GO:BP", "KEGG", "REAC", "H", "C2", "CUSTOM" |
| DE | Direction: "up" or "dn" |
| genes | Gene symbols (semicolon-separated) |
| gene_count | Number of DEGs in term |
| term_size | Total genes in term |
| comparison | Comparison name |

### Analysis Objects
| File | Description |
|------|-------------|
| `{prefix}_analysis_data.rds` | R objects with progressive checkpoints (stage: deseq2_gsva_complete → comparisons_complete → complete) |
| `{prefix}_analysis.RData` | Full R session |
| `{prefix}_Analysis_Report.md` | Generated markdown report |
| `{prefix}_reproducible.R` | Reproducible R script with all parameters filled in |
| `PLAN.md` | Analysis plan saved after user confirmation (Step 8) |

## Reproducibility

After analysis completes, a reproducible R script is automatically saved:

- **File**: `{output_dir}/{prefix}_reproducible.R`
- **Usage**: `conda run -n {conda_env} Rscript {prefix}_reproducible.R`

This script contains all parameters used in the analysis and can be:
1. Re-executed to reproduce results exactly
2. Modified and re-run with different parameters
3. Shared with collaborators for verification

### Conda Environment

The `--conda_env` parameter (default: `r_env`) specifies which conda environment to use for running the analysis. The reproducible script will preserve this setting.

```bash
# Using default r_env
conda run -n r_env Rscript scripts/deseq2_analysis.R --counts ...

# Using custom environment
conda run -n my_rna_env Rscript scripts/deseq2_analysis.R --conda_env my_rna_env --counts ...
```

### Reproducing an Analysis

```bash
# Navigate to the output directory
cd /path/to/output

# View the command (default - prints but doesn't execute)
Rscript RNAseq_reproducible.R

# To execute, edit the script and uncomment the system(cmd) line
```

---

## Standalone Visualization Functions

The skill includes 4 standalone plotting functions for ad-hoc visualization of RNA-seq results. These are **not** part of the main DESeq2 workflow - they can be invoked separately when users request specific visualizations.

### When to Use

Use these functions when users request:
- Heatmaps comparing signatures across conditions
- Gene expression visualizations for specific gene sets
- GSVA/ssGSEA score comparisons with statistical testing
- Custom visualizations of DEG or enrichment results

### Available Functions

| Function | Purpose | Primary Inputs |
|----------|---------|----------------|
| `plot_signature_heatmap()` | NES heatmap (signatures × comparisons) | GSEA results file |
| `plot_deg_heatmap()` | log2FC heatmap (genes × comparisons) | DEG summary stats file |
| `plot_expression_heatmap()` | Z-scored expression (genes × samples) | Expression matrix, metadata |
| `plot_gsva_boxplot()` | GSVA scores by treatment group | GSVA matrix, metadata |

### Script Location

```
scripts/plotting_functions.R
```

### Usage

Source the functions and call with user-specified parameters:

```r
# Load functions
source("scripts/plotting_functions.R")

# Example: Signature NES heatmap
plot_signature_heatmap(
  gsea_data = "deg/RNAseq_gsea_all.txt",
  signatures = c("sig1", "sig2", "sig3"),
  output_file = "figures/custom/signature_heatmap.png",
  source_filter = "CUSTOM"
)
```

### Required Inputs (User-Specified)

For each function, Claude should collect from user:

**plot_signature_heatmap:**
- GSEA results file path
- Signature names to include
- Source filter (CUSTOM, GO:BP, H, C2, etc.)
- Output path

**plot_deg_heatmap:**
- DEG summary stats file path
- Signature file(s) OR gene list
- Output path

**plot_expression_heatmap:**
- Expression data (RDS or matrix file)
- Metadata file
- Signature file(s) OR gene list
- Group column name
- Sample ID column name
- Output path

**plot_gsva_boxplot:**
- GSVA data (RDS or matrix file)
- Metadata file
- Signature names
- Group column name
- Statistical method (adaptive/kruskal/none)
- Max pairs to display
- Output path

### Output Files

All output files (static PNG, interactive HTML, and statistics TSV) are saved to subfolders under `figures/`:
- **Static PNG**: `figures/expr_comp_heatmap_boxplot/`
- **Interactive HTML**: `figures/interactive/expr_comp_heatmap_boxplot/`
- **Statistics TSV**: `figures/expr_comp_heatmap_boxplot/` (same as PNG)

| Function | PNG Output (in `expr_comp_heatmap_boxplot/`) | HTML Output (in `interactive/expr_comp_heatmap_boxplot/`) |
|----------|----------------------------------------------|-----------------------------------------------------------|
| plot_signature_heatmap | `*_NES_heatmap.png` | `*_NES_heatmap_interactive.html` |
| plot_deg_heatmap | `*_[signature]_deg_heatmap.png` | `*_[signature]_deg_heatmap_interactive.html` |
| plot_expression_heatmap | `*_[signature]_expression_heatmap.png` | `*_[signature]_expression_heatmap_interactive.html` |
| plot_gsva_boxplot | `*_[signature]_boxplot.png` | `*_[signature]_boxplot_interactive.html` |

**Additional Output for `plot_gsva_boxplot()`:**
- Pairwise statistics TSV: `figures/expr_comp_heatmap_boxplot/*_pairwise_stats.tsv` (when `export_stats = TRUE`)

### Interactive Plots

All plotting functions generate interactive HTML versions by default (`interactive = TRUE`). Interactive plots include hover tooltips with comprehensive statistical information:

| Function | Hover Information | Cell Annotation |
|----------|-------------------|-----------------|
| `plot_signature_heatmap()` | Signature name, Comparison, NES, p-value, FDR | `*` significance stars |
| `plot_deg_heatmap()` | Gene symbol, Comparison, log2FC, p-value (if available), FDR | `*` significance stars |
| `plot_expression_heatmap()` | Gene symbol, Sample, Treatment group, Z-score/expression value | None (colors show z-score) |
| `plot_gsva_boxplot()` | Treatment group, GSVA score, individual data points | N/A (boxplot) |

**Significance annotations (displayed in heatmap cells):**
- `*` p < 0.05
- `**` p < 0.01
- `***` p < 0.001
- `****` p < 0.0001

To disable interactive output, set `interactive = FALSE` in the function call.

**Note:** Interactive HTML plots require a web browser. They do NOT work in PDF, PowerPoint, or static documents.

### Dependencies

Additional R packages required:
- ComplexHeatmap
- ggsignif
- rstatix
- car
- heatmaply (for interactive heatmaps)
- plotly (for interactive plots)
- htmlwidgets (for saving HTML output)

See `references/plotting_guide.md` for detailed documentation.
