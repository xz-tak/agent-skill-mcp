# Internal Study Harmonization Guide

## Overview

This guide provides a step-by-step workflow for harmonizing and ingesting ANY new internal study into the S3 OmicSoft h5ad format. It produces per-study h5ad files that are concat-ready with the existing S3 OmicSoft dataset.

**When to use**: When you have internal experimental data (RNA-seq DEGs, proteomics DEPs, expression matrices) that needs to be integrated with the OmicSoft S3 h5ad reference for cross-study analysis.

**Accepted input formats**: TSV/CSV, Excel (.xlsx), RDS (Seurat/limma objects), h5ad, raw text matrices

**Output**: Per-study `*_deg.h5ad` and/or `*_expr.h5ad` files matching the S3 schema exactly.

---

## Information to Collect Per Study (Checklist)

Before starting harmonization, gather the following from the user:

- [ ] **Study name / project_id**: Short identifier (e.g., "SPARC", "Yokohama_RNA")
- [ ] **Organism**: human or mouse (determines gene symbol capitalization)
- [ ] **Platform**: RNA-seq, SomaScan, Olink, microarray (Affymetrix, Illumina)
- [ ] **DE method**: DESeq2, limma, edgeR, Wilcoxon (determines available statistics)
- [ ] **Tissue type(s)**: What tissue(s) were sampled
- [ ] **Disease/indication**: What disease or condition is being studied
- [ ] **Treatment arms** (if any): Drug names, dosages, vehicle controls
- [ ] **Timepoints** (if any): Baseline, week X, post-treatment
- [ ] **Response categories** (if any): Responder/non-responder definitions
- [ ] **Control definition**: What constitutes the control group (healthy, baseline, vehicle)
- [ ] **Significance cutoffs**: padj threshold, pval threshold, log2fc threshold
- [ ] **Sample metadata file**: Location of sample-level metadata (demographics, clinical)
- [ ] **DEG results file(s)**: Location of differential expression output files
- [ ] **Expression matrix** (if building EXPR h5ad): Normalized counts and/or raw counts

---

## Step-by-Step Workflow

### Step A: Identify Input Format and Gene Identifier Type

**Goal**: Determine file format and gene ID system to plan the mapping strategy.

| Input Format | Typical Gene ID | Detection Method |
|-------------|----------------|------------------|
| DESeq2 output (TSV/CSV) | Ensembl (ENSG) or HGNC symbol | Check first column |
| limma output (TSV/RDS) | Ensembl, HGNC, or probe IDs | Check rownames/ID column |
| edgeR output (TSV/CSV) | Ensembl or HGNC symbol | Check first column |
| SomaScan (CSV/Excel) | SeqId (e.g., 10000-28) | Numeric dash pattern |
| Olink (CSV/Excel) | UniProt or HGNC symbol | Check OlinkID column |
| Microarray (TSV) | Probe IDs (Affy/Illumina) | Platform-specific pattern |

**Decision tree for gene ID type**:
```
Does the ID start with "ENSG" or "ENSMUSG"?
  → YES: Ensembl gene ID (may have version suffix like .13)
  → NO: Continue...

Does the ID match [0-9]+-[0-9]+ pattern?
  → YES: SomaScan SeqId
  → NO: Continue...

Does the ID start with a number or contain "_at" suffix?
  → YES: Microarray probe ID
  → NO: Continue...

Is the ID all uppercase letters (human) or Title Case (mouse)?
  → YES: Likely HGNC/MGI gene symbol (ready to use)
  → NO: Check platform-specific annotation file
```

### Step B: Gene Symbol Harmonization

**Goal**: Map all gene identifiers to the S3 gene symbol vocabulary (60,084 DEG genes / 60,604 EXPR genes).

**Reference gene set**: The S3 h5ad var index defines the authoritative symbol set.

| Source ID Type | Mapping Strategy |
|---------------|-----------------|
| Versioned Ensembl (ENSG*.13) | Strip version → lookup in gene_annotation.txt → HGNC symbol |
| Unversioned Ensembl (ENSG*) | Direct lookup in gene_annotation.txt → gene_name |
| SeqId (SomaScan) | Lookup in SomaScan annotation → EntrezGeneSymbol |
| Olink ID | Lookup in Olink annotation → gene symbol |
| Microarray probe | Map via platform GPL annotation → gene symbol |
| HGNC symbol | Verify against S3 var_names; fix case if needed |

**Duplicate handling** (when multiple IDs map to one gene):
1. Rank by statistical significance (lowest average FDR across all comparisons)
2. Best-ranked gets the base symbol (e.g., "AKT3")
3. Subsequent get numbered suffixes: "AKT3_2", "AKT3_3"
4. Store original ID in var column (e.g., `var['SeqId']`, `var['probe_id']`)

**Unmapped genes**: Keep raw identifier as var index. Add to `gene_map.json` "unmapped" list.

**Output**: `gene_map.json` with structure:
```json
{
  "mapped": {"ENSG00000187634": "SAMD11", ...},
  "unmapped": ["ENSG00000239945", ...],
  "duplicates": {"AKT3": ["10000-28", "10001-7"], ...}
}
```

### Step C: Metadata Collection

**Goal**: Identify all available metadata and map to S3 vocabulary.

**Ask the user for**:
1. Sample-level metadata file (CSV/TSV/Excel with columns like patient_id, tissue, disease, sex, age)
2. Comparison definitions (which samples are case, which are control, for each comparison)

**Auto-mapping from `s3_vocabulary.json`**:
- For each metadata value, search the vocabulary for the closest match
- Present proposed mapping to user for confirmation
- Flag values with no obvious match for interactive resolution

**Tissue mapping decision tree**:
```
Input contains "colon" or "large intestine"? → "colon"
Input contains "ileum" or "small intestine"? → "ileum"
Input contains "liver" or "hepat"? → "liver"
Input contains "blood" or "serum" or "plasma"? → "peripheral blood"
Input contains "skin" or "derm"? → "skin"
Input contains "lung" or "pulmon"? → "lung"
Otherwise → Search s3_vocabulary.json tissue list for substring match
No match → Ask user to select from vocabulary or propose new term
```

**Disease mapping decision tree**:
```
Input contains "crohn" or "CD"? → "crohn's disease (CD)"
Input contains "ulcerative colitis" or "UC"? → "ulcerative colitis (UC)"
Input contains "IBD"? → "inflammatory bowel disease (IBD)"
Input contains "NASH" or "MASH" or "steatohepatitis"? → "non-alcoholic steatohepatitis (NASH)"
Input contains "normal" or "healthy" or "control"? → "normal control"
Otherwise → Search s3_vocabulary.json disease_state list
No match → Ask user to confirm closest match or propose term
```

### Step D: DEG Obs Schema Population (42 columns)

**Goal**: Populate all 42 obs columns for the DEG h5ad. Each row = one comparison.

The 42 columns (in order):
```
database, comparison_id, comparison_contrast, comparison_category,
case_tissue, case_sample_material, case_disease_state, case_disease_subtype,
case_disease_group, case_disease_location, case_response,
control_tissue, control_sample_material, control_disease_state, control_disease_subtype,
control_disease_group, control_disease_location, control_response,
case_dosage, case_treatment, case_treatment_group, case_treatment_status, case_treat_time,
control_dosage, control_treatment, control_treatment_group, control_treatment_status, control_treat_time,
case_age_category, case_gender, case_ethnicity,
control_age_category, control_gender, control_ethnicity,
case_sample_ids, control_sample_ids,
project_id, study, tissue, sample, disease, comparison
```

**Decision tree per column**:

| Column | Rule |
|--------|------|
| `database` | Always `"internal"` |
| `comparison_id` | Unique comparison name (e.g., filename stem or contrast name) |
| `comparison_contrast` | Free-text: `"case vs control"` or derived from group names |
| `comparison_category` | See comparison_category decision tree below |
| `case_tissue` / `control_tissue` | Map via tissue decision tree |
| `case_sample_material` / `control_sample_material` | Match to S3 vocabulary (biopsy, FFPE, fresh frozen, serum, etc.) |
| `case_disease_state` / `control_disease_state` | Map via disease decision tree |
| `case_disease_subtype` / `control_disease_subtype` | Specific subtype or `"NA"` |
| `case_disease_group` / `control_disease_group` | Grouping label or `"NA"` |
| `case_disease_location` / `control_disease_location` | Anatomical location or `"NA"` |
| `case_response` / `control_response` | `"response"`, `"no response"`, or `"NA"` |
| `case_dosage` / `control_dosage` | Drug dose or `"NA"` |
| `case_treatment` / `control_treatment` | Drug name (lowercase) or `"NA"` |
| `case_treatment_group` / `control_treatment_group` | Treatment arm label or `"NA"` |
| `case_treatment_status` / `control_treatment_status` | Pre/post-treatment or `"NA"` |
| `case_treat_time` / `control_treat_time` | Timepoint label or `"NA"` |
| `case_age_category` / `control_age_category` | `"adult"`, `"child"`, `"fetus"`, or `"NA"` |
| `case_gender` / `control_gender` | `"male"`, `"female"`, or `"NA"` |
| `case_ethnicity` / `control_ethnicity` | From S3 vocabulary or `"NA"` |
| `case_sample_ids` / `control_sample_ids` | Semicolon-separated EXPR sample IDs |
| `project_id` | Short project identifier |
| `study` | Descriptive study name |
| `tissue` | Summary tissue (same as case_tissue typically) |
| `sample` | Derived sample description |
| `disease` | Summary disease (same as case_disease_state typically) |
| `comparison` | Full human-readable comparison label |

### Step E: EXPR Obs Schema Population (30 columns)

**Goal**: Populate all 30 obs columns for the EXPR h5ad. Each row = one sample.

The 30 columns (in order):
```
database, tissue, project_id, disease_state, disease_stage,
ethnicity, gender, age_summary, treatment, sampling_time,
response, subject_id, cell_type, symptom, infection,
transfection, sample_integration_id, sample_pathology, sample_source,
sample_type, collection, data_source, organism, description,
title, platform_name, experiment_type, project_description,
project_title, comparison_group
```

**Key rules**:
- `database` = `"internal"`
- `organism` = `"Homo sapiens"` (human) or `"Mus musculus"` (mouse)
- `platform_name`: Use S3 vocabulary (e.g., `"RNA-seq"` for bulk sequencing, actual platform for arrays)
- `experiment_type`: `"Expression profiling by high throughput sequencing"` or `"Expression profiling by array"`
- `comparison_group`: Format is `{project_id}.{platform}.{method}.{comparison_id}@{case|control}` (semicolon-separated for multiple)
- `age_summary`: Format as `"X years"` or `"adult"`/`"child"` category
- `gender`: `"male"`, `"female"`, or `"NA"`
- All unavailable fields → `"NA"` or `""` (empty string)

### Step F: Comparison Structure

**Goal**: Define case/control assignments and comparison_category.

**Comparison category decision tree**:
```
Is it disease samples vs healthy/normal controls?
  → "Disease vs. Normal"

Is it drug-treated vs untreated/placebo?
  → "Treatment vs. Control"

Is it responders vs non-responders?
  → "Responder vs. Non-Responder"
  (NOTE: If internal study is NR vs R, use "Non-Responder vs. Responder"
   to preserve log2fc sign direction)

Is it one disease vs another disease?
  → "Disease1 vs. Disease2"

Is it one treatment vs another treatment?
  → "Treatment1 vs. Treatment2"

Is it one cell type vs another?
  → "CellType1 vs. CellType2"

Is it tissue comparisons?
  → "Tissue1 vs. Tissue2"

None of the above?
  → "Other Comparisons"
```

**Direction convention**: case = numerator, control = denominator in log2fc.
- Positive log2fc = upregulated in case vs control
- Negative log2fc = downregulated in case vs control

### Step G: Treatment/Response/Timepoint Mapping

**Goal**: Map treatment, response, and timepoint metadata to S3 vocabulary.

**Treatment mapping**:
- Use lowercase for drug names: `"adalimumab"`, `"vedolizumab"`, `"infliximab"`
- Match against `s3_vocabulary.json` treatment list
- For combination treatments, use semicolons: `"adalimumab;methotrexate"`
- No treatment = `"NA"`

**Response mapping**:
- `"response"` = responder (clinical response achieved)
- `"no response"` = non-responder
- `"partial response"` = partial responder
- Map common abbreviations: R → "response", NR → "no response"

**Timepoint mapping** (for `case_treat_time` / `sampling_time`):
- Use study-specific labels: `"Baseline"`, `"Week 14"`, `"Week 52"`
- Or generic: `"pre-treatment"`, `"post-treatment"`
- Match S3 vocabulary sampling_time format where possible

### Step H: sig_score Computation + Zero P-value Floor

**Goal**: Compute the sig_score layer for DEG h5ad.

**Formula**:
```python
sig_score = log2fc * -log10(padj)   # when significance criterion met
sig_score = 0                        # when NOT significant
```

**Pre-processing (CRITICAL)**:
1. Compute per-study global floor for zero p-values:
   ```python
   f32_tiny = np.finfo(np.float32).tiny  # ~1.175e-38
   pval_floor = max(pval_mat[pval_mat > 0].min() / 10.0, f32_tiny)
   padj_floor = max(padj_mat[padj_mat > 0].min() / 10.0, f32_tiny)
   ```
2. Replace exact-zero pval/padj with floor values (NOT NaN, NOT near-zero)
3. Apply significance gating (configurable per study):
   - Default: `padj < 0.05`
   - Store threshold in obs: `cutoff_padj`, `cutoff_pval`, `cutoff_log2fc`

**Significance gating**:
```python
is_significant = padj_mat < cutoff_padj  # or pval_mat < cutoff_pval
sig_score_mat = np.where(is_significant, log2fc_mat * -np.log10(padj_mat), 0.0)
```

### Step I: Validate H5AD Schema (Concat-Ready Checks)

**Goal**: Verify the per-study h5ad can cleanly concatenate with S3 data.

**Validation checklist**:
- [ ] `obs.index.name == "sample_id"`
- [ ] `var.index.name == "gene_id"`
- [ ] All 42 DEG obs columns present (or 30 EXPR obs columns)
- [ ] Column order matches S3 schema exactly
- [ ] Correct layers present:
  - DEG: `X` (log2fc), `pval`, `padj`, `sig_score`
  - EXPR: `X` (normalized), `raw_counts`
- [ ] No numeric var_names (all gene symbols)
- [ ] Unique obs indices within study
- [ ] `source` column = `"internal"`
- [ ] `.uns['schema_version']` = `"1.0"`
- [ ] `.uns['schema_type']` = `"deg"` or `"expr"`
- [ ] `.uns['available_layers']` lists all layers present

**NA fill conventions**:
| Layer | NA fill value | Rationale |
|-------|---------------|-----------|
| log2fc | 0.0 | No measured change |
| pval | 1.0 | Non-significant |
| padj | 1.0 | Non-significant |
| sig_score | 0.0 | Not significant |
| X (normalized expr) | 0.0 | No expression |
| raw_counts | 0.0 | No counts |

### Step J: Build Combined H5AD (Optional)

**Goal**: Merge new study with existing S3 + internal data.

**Ask user**: "Do you want to add this study to the combined h5ad?"

If yes:
```python
import anndata as ad

# Load existing combined
combined_deg = ad.read_h5ad("combined_deg.h5ad")

# Add source tag to new study
new_deg.obs['source'] = 'internal'

# Concatenate
updated = ad.concat([combined_deg, new_deg], join='outer')

# Verify no index collisions
assert updated.obs.index.is_unique

# Save
updated.write_h5ad("combined_deg_updated.h5ad")
```

### Step K: SOMA Conversion (Optional)

**Goal**: Convert h5ad to TileDB-SOMA format for cloud-native access.

**Ask user**: "Do you want to convert the updated combined h5ad to SOMA format for S3?"

If yes:
```python
import tiledbsoma

tiledbsoma.io.from_h5ad(
    "combined_deg_updated.h5ad",
    uri="s3://bucket/path/experiment_soma",
    measurement_name="RNA"
)
```

---

## Data Quality Rules

### NA Fill Conventions

| H5AD | Layer/Matrix | NA fill value | Rationale |
|------|-------------|---------------|-----------|
| DEG | log2fc (X) | 0.0 | No measured change |
| DEG | pval | 1.0 | Non-significant |
| DEG | padj | 1.0 | Non-significant |
| DEG | sig_score | 0.0 | Not significant |
| EXPR | X (normalized) | 0.0 | No expression detected |
| EXPR | raw_counts | 0.0 | No counts |

### Zero P-value Replacement Protocol

Statistical tools sometimes output exact-zero p-values when the true value underflows float64.

**Strategy**: Per-study global floor
```python
f32_tiny = np.finfo(np.float32).tiny  # ~1.175e-38
pval_floor = max(pval_mat[pval_mat > 0].min() / 10.0, f32_tiny)
padj_floor = max(padj_mat[padj_mat > 0].min() / 10.0, f32_tiny)
```

- Computed once per study (global across all comparisons and genes)
- Replaces exact zeros ONLY (not NaN, not near-zero)
- Floor is 10x smaller than the smallest non-zero value in the study
- Minimum floor: float32 tiny (~1.175e-38) to survive float32 storage

### Metadata Preservation Policy

ALL original metadata columns are preserved with `metadata_` prefix:
- `gender` = "female" (harmonized) + `metadata_SEX` = "F" (original)
- `tissue` = "ileum" (harmonized) + `metadata_CHARACTERISTICS_BIO_MATERIAL` = "ileum biopsy" (original)

This enables audit trails, recovery from incorrect mappings, and access to study-specific metadata.

### raw_counts Layer (EXPR)

The `raw_counts` layer must ALWAYS be populated:
- If raw count file exists and dimensions align → use true raw counts
- If raw counts unavailable → use normalized expression as fallback
- Track provenance: `adata.uns["raw_counts_source"]` = `"raw"` or `"normalized_expression"`

---

## Validation Checklist

Run after building each h5ad:

```python
def validate_concat_ready(adata, schema_type='deg'):
    """Validate h5ad is ready for concat with S3 OmicSoft data."""
    errors = []

    # Index names
    if adata.obs.index.name != 'sample_id':
        errors.append(f"obs.index.name = '{adata.obs.index.name}', expected 'sample_id'")
    if adata.var.index.name != 'gene_id':
        errors.append(f"var.index.name = '{adata.var.index.name}', expected 'gene_id'")

    # Column presence
    if schema_type == 'deg':
        required = [
            'database', 'comparison_id', 'comparison_contrast', 'comparison_category',
            'case_tissue', 'case_sample_material', 'case_disease_state', 'case_disease_subtype',
            'case_disease_group', 'case_disease_location', 'case_response',
            'control_tissue', 'control_sample_material', 'control_disease_state',
            'control_disease_subtype', 'control_disease_group', 'control_disease_location',
            'control_response', 'case_dosage', 'case_treatment', 'case_treatment_group',
            'case_treatment_status', 'case_treat_time', 'control_dosage', 'control_treatment',
            'control_treatment_group', 'control_treatment_status', 'control_treat_time',
            'case_age_category', 'case_gender', 'case_ethnicity',
            'control_age_category', 'control_gender', 'control_ethnicity',
            'case_sample_ids', 'control_sample_ids',
            'project_id', 'study', 'tissue', 'sample', 'disease', 'comparison'
        ]
        required_layers = ['pval', 'padj', 'sig_score']
    else:
        required = [
            'database', 'tissue', 'project_id', 'disease_state', 'disease_stage',
            'ethnicity', 'gender', 'age_summary', 'treatment', 'sampling_time',
            'response', 'subject_id', 'cell_type', 'symptom', 'infection',
            'transfection', 'sample_integration_id', 'sample_pathology', 'sample_source',
            'sample_type', 'collection', 'data_source', 'organism', 'description',
            'title', 'platform_name', 'experiment_type', 'project_description',
            'project_title', 'comparison_group'
        ]
        required_layers = ['raw_counts']

    missing_cols = [c for c in required if c not in adata.obs.columns]
    if missing_cols:
        errors.append(f"Missing obs columns: {missing_cols}")

    missing_layers = [l for l in required_layers if l not in adata.layers]
    if missing_layers:
        errors.append(f"Missing layers: {missing_layers}")

    # Unique indices
    if not adata.obs.index.is_unique:
        dups = adata.obs.index[adata.obs.index.duplicated()].tolist()
        errors.append(f"Duplicate obs indices: {dups[:5]}")

    # No numeric var_names
    numeric_vars = [v for v in adata.var_names[:100] if v.isdigit()]
    if numeric_vars:
        errors.append(f"Numeric var_names found: {numeric_vars[:5]}")

    # .uns metadata
    if 'schema_version' not in adata.uns:
        errors.append("Missing uns['schema_version']")
    if 'schema_type' not in adata.uns:
        errors.append("Missing uns['schema_type']")

    if errors:
        print("VALIDATION FAILED:")
        for e in errors:
            print(f"  ✗ {e}")
        return False
    else:
        print("✓ Validation passed - concat-ready!")
        return True
```

---

## Category Vocabulary Reference

All controlled vocabulary values are in `s3_vocabulary.json` (same directory as this file).

Key fields and their cardinality:
- `tissue`: 198 unique values
- `disease_state`: 280 unique values
- `comparison_category`: 8 values (fixed set)
- `treatment`: 501 unique values
- `response`: 14 values
- `gender`: 5 values
- `ethnicity`: 70 values
- `sample_material`: 58 values
- `platform_name`: 71 values
- `sampling_time`: 174 values

**Usage**: When mapping a local value, search the relevant vocabulary list for the best match. Present uncertain mappings to the user for confirmation.

---

## Additional Columns (Beyond S3 Schema)

Internal studies may add extra columns that will be NA for S3 rows after concat:

| Column | Type | Location | Value | Purpose |
|--------|------|----------|-------|---------|
| `source` | obs | DEG + EXPR | `"internal"` | Filter internal vs S3 after concat |
| `cutoff_pval` | obs | DEG | float | Per-comparison significance threshold |
| `cutoff_padj` | obs | DEG | float | Per-comparison significance threshold |
| `cutoff_log2fc` | obs | DEG | float | Per-comparison fold-change threshold |
| `metadata_*` | obs | EXPR | varies | All original study metadata |

---

## Example Reference

For a concrete worked example with 5 internal studies (Engitix FFPE, SPARC, Varsity, Yokohama RNA, Yokohama Protein), see:
`data_prepare/HARMONIZATION_EXAMPLE_IBD_NASH.md` in the project directory.
