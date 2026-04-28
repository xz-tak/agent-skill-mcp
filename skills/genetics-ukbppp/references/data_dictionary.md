# UKB-PPP Disease Signatures — Data Dictionary

## Source

UK Biobank Pharma Proteomics Project (UKB-PPP) plasma proteomics data.
Association analysis between plasma protein levels and diseases coded in ICD-10.

- **S3 Bucket:** `s3://tec-dev-usvga-11158-ukb-sumstats-share/UKBPPP_disease_signatures/`
- **AWS Profile:** `cmp-dev`

## Files

| File | Description | Size |
|------|-------------|------|
| `final.signatures.csv.gz` | Main results (~2M rows, 2924 genes × 685 ICD10 codes) | ~162 MB |
| `protein_info.tsv` | Protein annotations (2958 entries) | ~228 KB |
| `ICD10_info.tsv` | ICD-10 code definitions (19190 entries) | ~976 KB |

## Statistical Models

| Model | Prefix | Description | Diagnosis Timeframe | Data Coverage |
|-------|--------|-------------|---------------------|---------------|
| Logistic regression | `logit_` | logit Y ~ Covariates + Protein X | All patients regardless of diagnosis time | ~100% of rows |
| Automatic Relevance Determination | `ARD_` | Protein X ~ Covariates + Disease1 + ... + Disease_n | Within ±5 years of sample collection | ~21% of rows |
| Cox Proportional Hazard | `CoxPH_` | Survival analysis | After sample collection only | ~54% of rows |

**Note on ARD:** Coefficients for some diseases will be zero by nature of the ARD method (automatic sparsity).

## Column Descriptions — Main Data (final.signatures.csv.gz)

### Common Columns

| Column | Description |
|--------|-------------|
| `UKBPPP_ProteinID` | Unique protein identifier (format: GENE:UniProt:OlinkID:version:Panel) |
| `ICD10` | 3-character ICD-10 code (e.g., K50, D50) |
| `Assay` | Gene symbol (HGNC) for the measured protein |

### Per-Model Columns

Each model has 6 associated columns. Replace `{prefix}` with `logit_`, `ARD_`, or `CoxPH_`:

| Column | Description |
|--------|-------------|
| `{prefix}beta` | Effect size (log-odds ratio for logit, hazard ratio log for CoxPH, regression coefficient for ARD) |
| `{prefix}se` | Standard error of the estimated effect size |
| `{prefix}pval` | Raw p-value from the statistical test |
| `{prefix}Ncontrol` | Number of control subjects in the analysis |
| `{prefix}Ncase` | Number of case subjects (diagnosed with the disease) |
| `{prefix}FDR` | False discovery rate (Benjamini-Hochberg correction) |

**Column name inconsistency:** The source data uses `ARD_fdr` (lowercase) while logit and CoxPH use uppercase `FDR`. The query script normalizes this to `ARD_FDR` in exported output.

## Protein Info Columns (protein_info.tsv)

| Column | Description |
|--------|-------------|
| `UKBPPP_ProteinID` | Matches the main data protein identifier |
| `olink_target_fullname` | Full protein name from Olink assay |
| `HGNC.symbol` | Official gene symbol (matches `Assay` in main data) |
| `ensembl_id` | Ensembl gene identifier |

## ICD-10 Info Columns (ICD10_info.tsv)

| Column | Description |
|--------|-------------|
| `coding` | ICD-10 code (includes 3-char and sub-codes) |
| `meaning` | Disease/condition name |

**Note:** Only 3-character ICD-10 codes (e.g., K50, E11) exist in the main signatures data. Sub-codes (K500, K501) are available in ICD10_info.tsv for reference/search but map to the same parent 3-char code in results.

## Interpretation Guide

- **Beta > 0:** Higher protein levels associated with higher disease risk (logit/CoxPH) or disease presence (ARD)
- **Beta < 0:** Higher protein levels associated with lower disease risk
- **FDR < 0.05:** Conventionally considered statistically significant after multiple testing correction
- **Empty/NaN values:** The model was not applicable for that protein-disease pair (common for ARD and CoxPH)
- **ARD sign convention:** ARD models the protein as the response variable, so the sign interpretation may differ from logit/CoxPH
