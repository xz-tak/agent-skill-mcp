---
name: genetics-ukbppp
description: Query UKB-PPP plasma proteomics disease association data. Use this skill when users mention UKB-PPP, UKBPPP, plasma proteomics disease associations, protein-disease signatures, ICD10 proteomics results, logistic regression/ARD/CoxPH proteomics, or want gene-disease associations from UK Biobank Pharma Proteomics Project. Also triggers for "is [gene] associated with [disease] in proteomics" or "proteomics disease signature for [gene]".
---

# UKB-PPP Disease Signatures Query

Query plasma proteomics protein-disease associations from the UK Biobank Pharma Proteomics Project.

**Data:** ~2M rows, 2,924 genes, 685 ICD-10 codes, 3 statistical models (logit, ARD, CoxPH).
**Source:** `s3://tec-dev-usvga-11158-ukb-sumstats-share/UKBPPP_disease_signatures/`
**AWS Profile:** `cmp-dev` (SSO)

## Environment

All Python commands run in conda: `conda run -n claude_test python3 ...`

Script location: `${CLAUDE_SKILL_DIR}/scripts/ukbppp_query.py`

For reference documentation, read: `${CLAUDE_SKILL_DIR}/references/data_dictionary.md`

---

## Workflow

Follow these steps **in order**. Do not skip any step.

### Step 1: SSO Auth Check

**Always run this first.** Check if the AWS SSO session is valid:

```bash
aws s3 ls s3://tec-dev-usvga-11158-ukb-sumstats-share/UKBPPP_disease_signatures/ --profile cmp-dev 2>&1
```

**If it fails** (error contains "SSO session" or exit code non-zero):

```bash
aws sso login --profile cmp-dev 2>&1
```

Capture the output. It will contain a URL and a user code. Present both to the user clearly:

> **SSO authentication required.** Please open this URL and enter the code:
> - URL: `<url from output>`
> - Code: `<code from output>`

Wait for the user to confirm they have authenticated before proceeding.

### Step 2: Download Data

Download reference files and (if not cached) the main signatures file to `~/tmp/`:

```bash
conda run -n claude_test python3 ${CLAUDE_SKILL_DIR}/scripts/ukbppp_query.py --action download
```

The script:
- Caches `final.signatures.csv.gz` at `~/tmp/final.signatures.csv.gz` (reuses if exists, ~162MB)
- Always re-downloads `protein_info.tsv` and `ICD10_info.tsv` (small files)

### Step 3: Parse User Input

Extract from the user's request:
- **Gene list:** Comma-separated gene symbols. Uppercase them (e.g., `il11` → `IL11`).
- **Disease input:** Can be an ICD-10 code (e.g., `K50`) or a disease name/keyword (e.g., `Crohn`, `diabetes`).

### Step 4: Resolve Disease

Run the disease resolver:

```bash
conda run -n claude_test python3 ${CLAUDE_SKILL_DIR}/scripts/ukbppp_query.py --action resolve-disease --disease "<user_input>"
```

The script returns JSON with matching 3-character ICD-10 codes.

**You MUST use AskUserQuestion** to present matches and let the user pick:

```
"The following ICD-10 codes match your query '<input>':
Options:
  - K50 (K50 Crohn's disease [regional enteritis])
  - K51 (K51 Ulcerative colitis)
  - ... (additional matches)
```

If the user provided a valid 3-char ICD-10 code directly and it matches, you may skip the question and confirm: "Using ICD-10 code K50 (Crohn's disease)."

**If no matches found:** Inform the user and ask for an alternative disease name or ICD-10 code.

### Step 5: Confirm Statistical Model

**You MUST use AskUserQuestion** to confirm model selection:

```
"Which statistical model(s) should be exported?"
Options:
  - All 3 models (Recommended) — logit + ARD + CoxPH columns side by side
  - logit only — Logistic regression (all patients, regardless of diagnosis time)
  - ARD only — Automatic Relevance Determination (±5 years of sample collection)
  - CoxPH only — Cox Proportional Hazard (patients diagnosed after sample collection)
```

**Note:** ARD data is available for ~21% of rows and CoxPH for ~54%. Many gene-disease pairs will have empty ARD/CoxPH values. Mention this when presenting options.

### Step 5b: Combo Genes (Optional)

If the user mentions combos, pairs, or combinations of genes, use **AskUserQuestion** to confirm which combos to evaluate. Format: `GENE1+GENE2` separated by commas.

Example: User says "evaluate IL11+OSM and OSM+GREM1 combos" → `--combos "IL11+OSM,OSM+GREM1"`

If the user's intent is unclear (e.g., "compare these genes"), ask which specific combo groupings they want.

### Step 6: Run Query & Export

Execute the query with all confirmed parameters:

```bash
conda run -n claude_test python3 ${CLAUDE_SKILL_DIR}/scripts/ukbppp_query.py \
  --action query \
  --genes "GENE1,GENE2,GENE3" \
  --icd10 "K50,K51" \
  --models all \
  --combos "GENE1+GENE2,GENE2+GENE3" \
  --output-dir .
```

Replace `--models all` with `logit`, `ARD`, or `CoxPH` if the user selected a single model. Omit `--combos` if no combos requested.

The script produces:
- **`{GENE}_{ICD10_1}_{ICD10_2}_ukbppp.tsv`** per gene — filtered results with enriched columns (disease_name, protein name, Ensembl ID)
- **`ukbppp_query.log`** — query config, column descriptions, row counts
- **`ukbppp_summary.md`** — scored summary with individual gene rankings, per-indication detail tables, and combo rankings (if combos provided)

**After export:** Read `ukbppp_summary.md` and present the key findings to the user.

---

## Output Format

### Per-Gene TSV (`{GENE}_{ICD10_1}_{ICD10_2}_ukbppp.tsv`)

Default "all models" columns:

| Column | Source |
|--------|--------|
| Assay | Gene symbol |
| UKBPPP_ProteinID | Protein identifier |
| ICD10 | 3-char ICD-10 code |
| disease_name | From ICD10_info.tsv |
| logit_beta, logit_se, logit_pval, logit_Ncontrol, logit_Ncase, logit_FDR | Logistic regression |
| ARD_beta, ARD_se, ARD_pval, ARD_Ncontrol, ARD_Ncase, ARD_FDR | ARD model |
| CoxPH_beta, CoxPH_se, CoxPH_pval, CoxPH_Ncontrol, CoxPH_Ncase, CoxPH_FDR | Cox PH model |
| olink_target_fullname | From protein_info.tsv |
| ensembl_id | From protein_info.tsv |

### Summary (`ukbppp_summary.md`)

Sections:
1. **Scoring Method** — explanation of the direction-based scoring rules
2. **Individual Gene Rankings** — table ranked by total score (descending), with per-indication scores
3. **Gene Details** — per gene: total score, per-indication model tables with concordance, sorted by rank
4. **Combo Rankings** (if combos provided) — table ranked by total combo score
5. **Combo Details** — per combo: component gene scores per indication, combo averages

**Scoring:** Only FDR<0.05 models count. Net positive direction: +1→30, +2→60, +3→100. Net negative: -1→-30, -2→-60, -3→-100. Mixed cancel out. Gene total = mean across indications. Missing genes = 0. Combo = mean of component genes' scores.

### Log (`ukbppp_query.log`)

Query parameters, column descriptions, model notes, results summary.

---

## Statistical Models

| Model | Formula | Timeframe | Coverage |
|-------|---------|-----------|----------|
| **logit** | Y ~ Covariates + Protein X | All patients | ~100% |
| **ARD** | Protein X ~ Covariates + Diseases | ±5 years of sample | ~21% |
| **CoxPH** | Survival analysis | After sample collection | ~54% |

---

## Edge Cases

| Scenario | Handling |
|----------|----------|
| Gene not in data | Reported in summary as "not found" with score 0; no TSV file exported |
| ICD-10 not in data | Inform user, ask for alternative |
| Disease name no matches | Inform user, suggest alternative input |
| Disease name many matches | AskUserQuestion to let user pick |
| Empty ARD/CoxPH values | Shown as "—" / "N/A" in summary; NaN in TSV |
| All models empty for a pair | Noted in summary as "No data available" |
| S3 download failure | Script reports error; check SSO auth |
