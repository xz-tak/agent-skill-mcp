# DESeq2 RNA-seq Workflow Guide

Detailed documentation for the interactive bulk RNA-seq differential expression analysis workflow.

## Table of Contents

1. [Input File Requirements](#input-file-requirements)
2. [Interactive Workflow Steps](#interactive-workflow-steps)
3. [Design Formula Guide](#design-formula-guide)
4. [Collinearity Troubleshooting](#collinearity-troubleshooting)
5. [Comparison Configuration](#comparison-configuration)
6. [Advanced Comparisons (List Contrasts)](#advanced-comparisons-list-contrasts)
7. [Output Files Reference](#output-files-reference)

---

## Input File Requirements

### Counts File

The counts file contains raw gene expression counts (not normalized values).

#### Format Requirements

| Requirement | Description |
|-------------|-------------|
| **Values** | Non-negative integers (raw counts) |
| **Rows** | Genes (one per row) |
| **Columns** | Samples (one per column) |
| **Row names** | Gene identifiers (ENSEMBL IDs or gene symbols) |
| **Column names** | Sample identifiers |
| **Delimiter** | Tab-separated (.txt/.tsv) or comma-separated (.csv) |

#### Example Counts File (TSV)

```
Gene_ID	Sample1	Sample2	Sample3	Sample4
ENSG00000000003	1250	1340	1180	1420
ENSG00000000005	0	2	1	0
ENSG00000000419	850	920	780	910
ENSG00000000457	320	280	350	310
```

#### Example Counts File (CSV)

```
Gene_ID,Sample1,Sample2,Sample3,Sample4
ENSG00000000003,1250,1340,1180,1420
ENSG00000000005,0,2,1,0
ENSG00000000419,850,920,780,910
ENSG00000000457,320,280,350,310
```

#### Optional Gene.name Column

If your counts file includes a `Gene.name` column for gene symbols, it will be automatically extracted and used for labeling:

```
Gene_ID	Gene.name	Sample1	Sample2	Sample3
ENSG00000000003	TSPAN6	1250	1340	1180
ENSG00000000005	TNMD	0	2	1
ENSG00000000419	DPM1	850	920	780
```

#### Common Issues

| Issue | Solution |
|-------|----------|
| Negative values | Counts must be non-negative; check for data corruption |
| Non-integer values | Values will be rounded; consider if data is already normalized |
| Missing values | Impute with 0 or remove genes/samples with missing data |
| Duplicate gene names | Keep ENSEMBL IDs as row names for uniqueness |

---

### Metadata File

The metadata file contains sample information and experimental design factors.

#### Format Requirements

| Requirement | Description |
|-------------|-------------|
| **Rows** | Samples (one per row) |
| **Columns** | Sample attributes |
| **Sample column** | One column must contain sample names matching counts columns |
| **Delimiter** | Tab-separated |

#### Required Columns

At minimum, metadata must include:
1. **Sample identifier** - Matches counts column names
2. **Comparison column** - Factor for DESeq2 contrasts (e.g., Treatment, Condition)

#### Example Metadata File

```
Sample_ID	Treatment	Batch	Replicate
Sample1	Control	Batch1	1
Sample2	Control	Batch1	2
Sample3	Treatment_A	Batch1	1
Sample4	Treatment_A	Batch1	2
Sample5	Treatment_B	Batch2	1
Sample6	Treatment_B	Batch2	2
```

#### Column Types

| Type | Description | Use in Design |
|------|-------------|---------------|
| **Categorical** | Discrete groups (Treatment, Batch) | Comparison factor or covariate |
| **Numeric** | Continuous values (Age, RIN) | Covariate only |
| **Identifier** | Unique sample names | Not used in design |

---

## Interactive Workflow Steps

### Step 1: Environment Validation

The `check_r_env.R` script validates:

1. **R Version** - Reports current R version
2. **Core Packages** - Checks all required packages are installed
3. **Species Packages** - Validates annotation package for selected species
4. **Counts File** - Validates format and data integrity
5. **Metadata File** - Validates structure and sample matching

#### Running Validation

```bash
conda run -n r_env Rscript scripts/check_r_env.R \
  --counts /path/to/counts.txt \
  --metadata /path/to/metadata.txt \
  --species hsapiens
```

### Step 2: Metadata Summary

Present a summary table to the user:

| Column | Type | Unique Values | Sample Values | Suggested Role |
|--------|------|---------------|---------------|----------------|
| Sample_ID | identifier | 6 | Sample1, Sample2, ... | Identifier |
| Treatment | categorical | 3 | Control, Treatment_A, Treatment_B | Comparison factor |
| Batch | categorical | 2 | Batch1, Batch2 | Covariate |
| Replicate | numeric | 2 | 1, 2 | Not recommended |

### Step 3: Column Selection

Collect from user:
- **Primary comparison column**: Main factor for DESeq2 contrasts
- **Covariate columns** (optional): Additional factors for batch correction

### Step 4: Design Formula Validation

Test for collinearity before proceeding. See [Collinearity Troubleshooting](#collinearity-troubleshooting).

### Step 5: Comparison Configuration

Collect from user:
- **Reference level**: Baseline group for comparisons
- **Comparison pairs**: Specific pairs or all pairwise

### Step 6: Species and Signatures

Collect from user:
- **Species**: Human, mouse, or rat
- **MSigDB collections**: H, C2, C5, C7, etc.
- **Custom signatures**: Optional file path

### Step 7: Significance Thresholds

Collect from user:
- **padj threshold**: Default 0.05
- **|log2FC| threshold**: Default 1.0

### Step 8: Plan Confirmation

Present complete analysis plan and confirm before execution.

### Step 9: Execute Analysis

Run the parameterized `deseq2_analysis.R` script.

### Step 10: Generate Report

Run `generate_report.py` to create markdown report.

---

## Design Formula Guide

### Basic Design Formula

The simplest design includes only the comparison factor:

```r
design = ~Treatment
```

This tests for differences between treatment groups.

### Design with Covariates

Include covariates to account for batch effects or confounders:

```r
design = ~Batch + Treatment
```

**Order matters:** Covariates should come before the comparison factor.

### Multiple Covariates

```r
design = ~Batch + Sex + Age + Treatment
```

### Interaction Terms

To test if treatment effects differ by another factor:

```r
design = ~Batch + Sex + Treatment + Sex:Treatment
```

### Nested Designs

For samples nested within subjects:

```r
design = ~Subject + Treatment
```

### Examples by Study Type

| Study Type | Design Formula |
|------------|----------------|
| Simple two-group | `~Condition` |
| With batch correction | `~Batch + Condition` |
| Paired samples | `~Patient + Treatment` |
| Time course | `~Subject + Time` |
| Two factors | `~FactorA + FactorB` |
| Interaction | `~FactorA + FactorB + FactorA:FactorB` |

---

## Collinearity Troubleshooting

### What is Collinearity?

Collinearity occurs when two or more columns in the design matrix are linearly dependent. DESeq2 cannot estimate coefficients when the design matrix is not full rank.

### Testing for Collinearity

```r
# Create model matrix
coldata <- read.table("metadata.txt", header=TRUE, row.names=1)
coldata$Treatment <- factor(coldata$Treatment)
coldata$Batch <- factor(coldata$Batch)

# Test design
model_matrix <- model.matrix(~Batch + Treatment, data = coldata)
qr_result <- qr(model_matrix)

if (qr_result$rank < ncol(model_matrix)) {
  cat("WARNING: Design matrix is rank deficient (collinear)\n")
  cat("Rank:", qr_result$rank, "vs Expected:", ncol(model_matrix), "\n")
}
```

### Common Causes

| Cause | Example | Solution |
|-------|---------|----------|
| **Perfect confounding** | All Batch1 samples are Control | Remove Batch from design or balance experiment |
| **Nested factors** | Subject nested in Treatment | Use `~Treatment` (Subject absorbed) |
| **Empty factor levels** | Treatment has level with no samples | Drop empty levels: `droplevels()` |
| **Redundant columns** | Both Sex and Gender included | Remove one column |

### Identifying Problematic Columns

```r
# Check for aliased coefficients
alias_check <- alias(lm(~ Batch + Treatment, data = coldata))
if (!is.null(alias_check$Complete)) {
  print("Aliased (collinear) terms:")
  print(alias_check$Complete)
}
```

### Solutions

1. **Remove confounded covariate**
   ```r
   # If Batch is confounded with Treatment
   design = ~Treatment  # Remove Batch
   ```

2. **Combine factor levels**
   ```r
   # If Batch1 = Control and Batch2 = Treatment
   # Data is fundamentally confounded; cannot separate effects
   ```

3. **Drop empty levels**
   ```r
   coldata$Treatment <- droplevels(coldata$Treatment)
   ```

4. **Simplify design**
   ```r
   # Instead of ~Batch + Sex + Treatment
   design = ~Batch + Treatment  # If Sex is confounded
   ```

---

## Comparison Configuration

### Reference Level

The reference level is the baseline group for all comparisons. Choose a biologically meaningful control:

| Study Type | Good Reference |
|------------|----------------|
| Treatment study | "Control" or "Vehicle" |
| Time course | "Day0" or "Baseline" |
| Disease study | "Healthy" or "Normal" |
| Dose response | "Untreated" or "0mg" |

### Comparison Pairs

#### All Pairwise

Compare every group to every other group:

```json
"all_pairwise"
```

For 3 groups (A, B, C), this generates:
- B vs A
- C vs A
- C vs B

#### Specific Pairs

Define exact comparisons:

```json
[["Treatment_A", "Control"], ["Treatment_B", "Control"]]
```

Each pair is `[Numerator, Denominator]`:
- Positive log2FC = higher in Numerator
- Negative log2FC = higher in Denominator

### Comparison Naming

Comparisons are named automatically:
```
{Numerator}_vs_{Denominator}
```

Special characters are replaced with underscores:
- `TGFβ+ALK5i vs TGFβ` → `TGF__ALK5i_vs_TGF_`

---

## Advanced Comparisons (List Contrasts)

For complex comparisons involving interaction terms, reversal effects, or coefficient comparisons, use list contrasts instead of simple contrasts.

### When to Use List Contrasts

| Scenario | Use Simple Contrast | Use List Contrast |
|----------|---------------------|-------------------|
| A vs B (direct comparison) | ✅ Yes | ❌ No |
| Treatment effect in one condition | ❌ No | ✅ Yes |
| Interaction effects | ❌ No | ✅ Yes (Required) |
| Reversal/rescue experiments | ❌ No | ✅ Yes |
| Comparing coefficients | ❌ No | ✅ Yes |

### Statistical Validity

**Simple contrasts** use `contrast = c(factor, numerator, denominator)`:
- Appropriate for direct A vs B comparisons
- DESeq2 computes proper Wald statistics

**List contrasts** use `contrast = list(...), listValues = c(...)`:
- Proper Wald test for coefficient comparisons
- Statistically valid for complex hypotheses
- **Required** for interaction/reversal comparisons

### JSON Configuration Format

```json
{
  "comparisons": [["TreatA", "Control"]],
  "advanced_comparisons": [
    {
      "name": "Descriptive_Comparison_Name",
      "type": "list",
      "contrast": ["coef1", "coef2"],
      "listValues": [1, -1]
    }
  ]
}
```

### Finding Coefficient Names

DESeq2 coefficient names follow the pattern `{factor}_{level}_vs_{reference}`:

```r
resultsNames(dds)
# Example output:
# [1] "Intercept"
# [2] "Treatment_DrugA_vs_DMSO"
# [3] "Treatment_DrugB_vs_DMSO"
# [4] "Stimulation_Stim_vs_Unstim"
# [5] "Treatment_DrugA_vs_DMSO.Stimulation_Stim_vs_Unstim"  # interaction
```

The script automatically prints available coefficients before processing advanced comparisons.

### listValues Patterns

| Pattern | Meaning | Use Case |
|---------|---------|----------|
| `[1, -1]` | coef1 - coef2 | Compare two coefficients |
| `[1, -2]` | (coef1 + coef2)/2 - coef3 | Average vs single (scaled) |
| `[-1, 1]` | coef2 - coef1 | Reverse direction |

### Example: Drug Effect in Stimulated Condition

**Biological question:** What is the effect of DrugA specifically in stimulated cells?

**Design:** `~Treatment + Stimulation + Treatment:Stimulation`

**Coefficients available:**
- `Treatment_DrugA_vs_DMSO` (drug effect in unstimulated)
- `Treatment_DrugA_vs_DMSO.Stimulation_Stim_vs_Unstim` (interaction)

**List contrast configuration:**
```json
{
  "name": "DrugA_effect_in_Stimulated",
  "type": "list",
  "contrast": [
    ["Treatment_DrugA_vs_DMSO", "Treatment_DrugA_vs_DMSO.Stimulation_Stim_vs_Unstim"],
    "Treatment_DMSO_vs_DMSO"
  ],
  "listValues": [1, -1]
}
```

This computes: (DrugA_Unstim + Interaction) - DMSO_baseline = DrugA effect in stimulated cells

### Example: Reversal/Rescue Experiment

**Biological question:** Does DrugA reverse the TGFβ-induced gene expression changes?

**Comparisons needed:**
1. TGFβ_vs_Control (what TGFβ does)
2. TGFβ+DrugA_vs_TGFβ (what adding DrugA does)
3. Reversal effect (does DrugA reverse TGFβ changes?)

**List contrast for reversal:**
```json
{
  "name": "DrugA_reverses_TGFb",
  "type": "list",
  "contrast": [
    "Treatment_TGFb_vs_Control",
    "Treatment_TGFb_DrugA_vs_Control"
  ],
  "listValues": [1, -1]
}
```

Interpretation:
- Positive NES/log2FC: TGFβ effect is reversed by DrugA
- Negative NES/log2FC: DrugA amplifies TGFβ effect

### Example: Multiple Coefficient Comparison

**Biological question:** Average effect of two treatments vs control

```json
{
  "name": "Average_DrugAB_vs_Control",
  "type": "list",
  "contrast": [
    ["Treatment_DrugA_vs_Control", "Treatment_DrugB_vs_Control"]
  ],
  "listValues": [0.5]
}
```

### Complete Advanced Comparisons Example

```json
{
  "counts": "counts.txt",
  "metadata": "metadata.txt",
  "design": "~Batch + Treatment + Stimulation + Treatment:Stimulation",
  "comparison_col": "Treatment",
  "reference": "DMSO",
  "comparisons": [
    ["DrugA", "DMSO"],
    ["DrugB", "DMSO"]
  ],
  "advanced_comparisons": [
    {
      "name": "DrugA_in_Stimulated",
      "type": "list",
      "contrast": [
        ["Treatment_DrugA_vs_DMSO", "Treatment_DrugA_vs_DMSO.Stimulation_Stim_vs_Unstim"]
      ],
      "listValues": [1]
    },
    {
      "name": "DrugA_reversal_of_Stim",
      "type": "list",
      "contrast": [
        "Stimulation_Stim_vs_Unstim",
        ["Stimulation_Stim_vs_Unstim", "Treatment_DrugA_vs_DMSO.Stimulation_Stim_vs_Unstim"]
      ],
      "listValues": [1, -1]
    }
  ],
  "species": "hsapiens",
  "msigdb": ["H", "C2"],
  "padj": 0.05,
  "lfc": 1
}
```

### Troubleshooting List Contrasts

| Error | Cause | Solution |
|-------|-------|----------|
| "Coefficient not found" | Typo in coefficient name | Check `resultsNames(dds)` output |
| "listValues length mismatch" | Wrong number of values | Match length to contrast elements |
| Unexpected results | Wrong coefficient combination | Verify mathematical formula |

### When NOT to Use List Contrasts

- Simple A vs B comparisons → use regular `comparisons`
- When reference level handles the comparison → use regular `comparisons`
- When you just need all pairwise → use `"all_pairwise"`

---

## Output Files Reference

### Directory Structure

```
output_dir/
├── figures/
│   ├── {prefix}_PCA_basic_{comparison_col}.png
│   ├── {prefix}_PCA_screeplot.png
│   ├── {prefix}_PCA_eigencorplot.png
│   ├── volcano_{comparison}.png
│   ├── GSEA_dotplot_{comparison}.png
│   └── DE_genes_summary_barplot.png
├── deg/
│   ├── {prefix}_{comparison}.txt
│   ├── {prefix}_{comparison}_GOBP.txt
│   ├── {prefix}_{comparison}_gp.txt
│   ├── {prefix}_summstats_all.txt
│   ├── {prefix}_gseaGObp_all.txt
│   ├── {prefix}_defisher_all.txt
│   └── {prefix}_DE_lfc{lfc}padj{padj}_count.txt
├── {prefix}_analysis_data.rds
├── {prefix}_analysis.RData
└── {prefix}_Analysis_Report.md
```

### DEG Results Files

| Column | Description |
|--------|-------------|
| gene | Gene identifier (row name) |
| symbol | Gene symbol (if available) |
| baseMean | Mean normalized count |
| log2FoldChange | Effect size |
| lfcSE | Standard error of LFC |
| stat | Wald statistic |
| pvalue | Raw p-value |
| padj | Adjusted p-value (BH) |
| comparison | Comparison name |

### GSEA Results Files

| Column | Description |
|--------|-------------|
| ID | GO term ID |
| Description | GO term name |
| setSize | Genes in term |
| enrichmentScore | Raw ES |
| NES | Normalized ES |
| pvalue | Nominal p-value |
| p.adjust | Adjusted p-value |
| core_enrichment | Leading edge genes |

### g:Profiler Results Files

| Column | Description |
|--------|-------------|
| query | Input gene list |
| source | Database (GO:BP, KEGG, REAC) |
| term_id | Term identifier |
| term_name | Term description |
| p_value | Raw p-value |
| term_size | Total genes in term |
| intersection_size | Query genes in term |
| DE | Direction (up/dn) |

---

## Best Practices

### Sample Size

| Scenario | Minimum Samples | Recommended |
|----------|-----------------|-------------|
| Simple comparison | 3 per group | 5+ per group |
| With covariates | 5 per group | 8+ per group |
| Complex design | 8 per group | 10+ per group |

### Gene Filtering

- Default: Remove genes with < 10 total counts
- Alternative: Remove genes with < 5 counts in < 2 samples

### Multiple Testing

- Always use adjusted p-values (padj) for significance
- padj < 0.05 is standard threshold
- Consider stricter threshold (0.01) for large studies

### Effect Size

- |log2FC| > 1 corresponds to 2-fold change
- |log2FC| > 0.585 corresponds to 1.5-fold change
- Consider biological context when choosing threshold
