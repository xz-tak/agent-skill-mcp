---
name: genetics-gsp
description: Primary entry point for genetics analysis. Query aggregated human genetics data from GSP (Genetics Support Profiler), including GWAS, Mendelian genetics (OMIM), and gene burden, for target safety concerns, disease associations, and evidence supporting drug targets risk assessment.
category: workflow
version: 0.4.0
---

# Genetics GSP Analysis Skill

Comprehensive tools for analyzing GSP (Genetics Support Profiler) data for drug target validation.

**Workflows:**
1. **Biological Risk Assessment** - Evaluate genetic evidence for target-disease associations (primary)
2. **Safety Risk Assessment** - Identify potential safety concerns when targeting a gene
3. **Combo / Multi-Gene Scoring** - Score gene combinations across phenotypes with ranked reports

**Environment:** All Python commands must run in `conda run -n claude_test` (openpyxl dependency).

---

# Dataset Registry

On every invocation, read the dataset registry at `references/datasets.md` (YAML frontmatter). This maps dataset names to direct S3 file paths and available indications.

## Available Datasets

| Dataset | Indications | EFO IDs |
|---------|-------------|---------|
| **IBD** | Inflammatory bowel disease, Ulcerative colitis, Crohn's disease | EFO_0003767, EFO_0000729, EFO_0000384 |
| **SSc** | Idiopathic pulmonary fibrosis, Systemic sclerosis, Interstitial lung disease | EFO_0000768, EFO_0000717, EFO_0004244 |
| **HS** | Hidradenitis suppurativa | EFO_1000710 |
| **AtD** | Atopic dermatitis, Recalcitrant atopic dermatitis | EFO_0000274, EFO_1000651 |

## Resolution Workflow

When the user mentions a disease or indication, follow this 3-step resolution:

### Step 1: Programmatic match
Call `resolve_dataset(user_input)` — matches dataset names, indication names, aliases, and EFO IDs.

```python
from gsp_tools import resolve_dataset
ds = resolve_dataset("Crohn's")  # Returns the IBD dataset entry
```

### Step 2: Claude interpretation (if Step 1 returns None)
If `resolve_dataset()` returns None, use your biomedical knowledge to interpret the user's intent:

- **Synonyms/related terms**: "bowel disease" -> IBD, "fibrosis" -> SSc/IPF, "skin condition with abscesses" -> HS, "eczema" -> AtD, "colitis" -> UC/IBD
- **Abbreviations you recognize**: "PSO" -> Psoriasis (not in registry), "RA" -> Rheumatoid arthritis (not in registry)
- **Parent/child disease relationships**: "autoimmune GI" -> likely IBD

If you can confidently map the user's text to a registry dataset, proceed. If you are unsure (>1 possible dataset, or the mapping is a stretch), go to Step 3.

### Step 3: Validate with user (AskUserQuestion)
Use **AskUserQuestion** when:
- **Ambiguous**: user says "inflammation" (could be IBD, SSc, HS, or AtD) -> present dataset options
- **No match and unsure**: user mentions a disease not in registry and you can't confidently map it
- **Multiple interpretations**: "fibrosis" could be IPF or ILD within the SSc dataset -> ask which indication(s)

Example AskUserQuestion options:
```
"Which dataset matches your query '{user_text}'?"
Options:
  - IBD (Inflammatory bowel disease, UC, Crohn's)
  - SSc (IPF, Systemic sclerosis, ILD)  
  - HS (Hidradenitis suppurativa)
  - AtD (Atopic dermatitis, Recalcitrant atopic dermatitis)
  - None of these (specify manually)
```

**Key principle**: Never silently skip or guess wrong. If there's any doubt, ask.

## Mandatory Phenotype & Indication Alignment (BEFORE Data Download)

> **CRITICAL — ALL WORKFLOWS:** After resolving the dataset but **before** downloading data or running any analysis, you **MUST** use AskUserQuestion to confirm phenotype/indication coverage with the user. Never skip this step, even when the user's intent seems obvious.

Once the dataset is resolved, present the available indications and ask the user to confirm:

**Q1: Indication/Phenotype coverage** (always ask)
```
"The {DATASET} dataset contains the following indications:

{For each indication in the dataset registry:}
  - {name} ({efo_id}) {aliases if any}

Which indications should be included in this analysis?"
Options:
  - All indications (Recommended)
  - {List each individual indication as a selectable option}
  - Custom subset (specify)
```

**Q2: Aggregation method** (ask when multiple phenotypes are selected)
```
"How should scores be aggregated across the selected phenotypes?"
Options:
  - Average across all selected (Recommended)
  - Best (maximum) across selected
```

**Single-phenotype datasets** (e.g., HS with only Hidradenitis suppurativa): Still confirm coverage — do not silently assume.
```
"The {DATASET} dataset contains one indication:
  - {name} ({efo_id})

Proceed with this indication?"
Options:
  - Yes, proceed
  - No, I need a different dataset
```

**Only after the user confirms** should you proceed to data download (S3) and downstream analysis.

## Adding/Updating Datasets

When user provides a new S3 path:
1. **Launch a subagent** to read the Study Info sheet from S3 (pipe via `aws s3 cp ... - | gunzip`) and extract unique EFO IDs + trait names
2. Compare discovered indications with existing registry (if updating)
3. If the dataset name doesn't cover all indications -> **AskUserQuestion** to suggest a name change
4. Update `references/datasets.md` with new/updated entry

---

# Quick Start

> **Reminder:** Before downloading data or running any analysis, you **MUST** first align with the user on which phenotypes/indications to include via AskUserQuestion. See "Mandatory Phenotype & Indication Alignment" section above.

```python
# Discover indications in a dataset
from gsp_tools import list_indications
indications = list_indications("tmp/IBD_20260407")

# Generate a full report for a gene across all indications
from gsp_tools import generate_biological_report
report = generate_biological_report("tmp/IBD_20260407", "OSMR", dataset_name="IBD")
# -> gsp_IBD/OSMR_IBD_gsp_risk.md
```

**CLI:**
```bash
conda run -n claude_test workflow-gsp indications tmp/IBD_20260407
conda run -n claude_test workflow-gsp report tmp/IBD_20260407 OSMR --dataset-name IBD
conda run -n claude_test workflow-gsp safety IBD OSMR Medium
conda run -n claude_test workflow-gsp biological IBD OSMR EFO_0003767
```

---

# Data Loading

## S3 Download

The registry stores **direct S3 paths** to GSP.pkl.gz and GSP.xlsx - no S3 listing needed:

```python
from gsp_tools import download_gsp_from_s3
data_dir = download_gsp_from_s3(
    s3_path="s3://bucket/gsp/IBD_20260407/",
    data_dir="tmp/IBD_20260407",
    aws_profile="cmp-dev"
)
```

**Default data directory:** When `data_dir` is omitted, files download to `~/tmp` (user home directory). This keeps large data files in a stable, shared location across projects and avoids polluting the skill installation directory. The `local_dir` field in the dataset registry provides the recommended path for each dataset (e.g., `~/tmp/AtD_20260410`).

## Loading Priority (with caching)

1. **In-memory cache** -> instant (same process)
2. **GSP_v2.pkl.gz** -> 2-5s (self-healed compatible pickle)
3. **GSP.pkl.gz** -> may fail (pandas version mismatch)
4. **GSP.xlsx** -> 30-60s (fallback; auto-writes GSP_v2.pkl.gz for next time)

After first xlsx load, a compatible pickle is written automatically. All subsequent loads are fast.

## GSP Data Sources

| Source | Sheet | Description |
|--------|-------|-------------|
| **GWAS** | GWAS Summary | Common variant associations from GWAS Catalog, UK Biobank |
| **OMIM** | OMIM | Mendelian disease phenotypes with direction annotations |
| **Gene Burden** | OT Gene Burden | Rare variant burden associations from exome sequencing |

---

# Workflow 1: Biological Risk Assessment (Primary)

## End-to-End: Single Gene Report

```
┌─────────────────────────────────────────────────────────────────┐
│ BIOLOGICAL RISK ASSESSMENT                                      │
├─────────────────────────────────────────────────────────────────┤
│  Step 0: RESOLVE DATASET                                        │
│  ├── Read references/datasets.md registry                       │
│  ├── Match user's disease text to dataset + EFO IDs             │
│  └── If ambiguous -> AskUserQuestion                            │
│                                                                 │
│  Step 1: ALIGN ON PHENOTYPES & INDICATIONS (MANDATORY)          │
│  ├── Present available indications from registry to user        │
│  ├── AskUserQuestion: which phenotypes to include?              │
│  ├── For single-phenotype datasets, confirm coverage            │
│  └── Only proceed after user confirms                           │
│                                                                 │
│  Step 2: DOWNLOAD DATA (if not local)                           │
│  ├── Use s3_pickle/s3_xlsx paths from registry (no S3 search)   │
│  └── download_gsp_from_s3(s3_path, data_dir, aws_profile)      │
│                                                                 │
│  Step 3: GENERATE REPORT                                        │
│  ├── generate_biological_report(data_dir, gene, diseases,       │
│  │     dataset_name=name)                                       │
│  ├── Uses only user-confirmed indications (not auto-discover)   │
│  ├── Runs expand_gene_disease_pair per confirmed disease        │
│  ├── Applies biological_risk_score (7-rule algorithm)           │
│  └── Writes {workdir}/gsp_{dataset}/{GENE}_{dataset}_gsp_risk.md│
└─────────────────────────────────────────────────────────────────┘
```

## Batch Processing: Multiple Genes (Parallel Agents)

When the user provides a **list of genes**, spawn parallel agents:

```
User: "Assess OSMR, IL23R, NOD2 against IBD dataset"

Orchestrator:
  1. Resolve dataset -> IBD, EFOs = [EFO_0003767, EFO_0000729, EFO_0000384]
  2. Download data once (if not local)
  3. For EACH gene, spawn a parallel Agent:
       conda run -n claude_test python3 -c "
       import sys; sys.path.insert(0, '...')
       from gsp_tools import generate_biological_report
       generate_biological_report('tmp/IBD_20260407', '{GENE}',
           dataset_name='IBD')
       "
  4. Each agent writes: gsp_IBD/{GENE}_IBD_gsp_risk.md
  5. Collect and summarize all reports
```

No shared state needed - each agent writes its own file. Cache is per-process (first call may be slow; subsequent fast via GSP_v2.pkl.gz).

## Biological Risk Scoring (Automated)

The 7-rule scoring algorithm is now codified in `biological_risk_score()`:

```python
from gsp_tools import expand_gene_disease_pair, biological_risk_score

result = expand_gene_disease_pair("IBD", "OSMR", "EFO_0003767", data_dir="tmp/IBD_20260407")
score = biological_risk_score(result['results'])

print(score['risk_level'])        # "High risk" / "Medium risk" / "Low risk"
print(score['risk_label'])        # Full description
print(score['cautionary_notes'])  # Any caveats
```

### Scoring Rules Reference

1. **High-quality threshold**: max_sample_size * 0.1 (binary traits use nCases; quantitative use nSamples)
2. All HQ studies "Not Sig./Low/Very low" -> **High risk**: Low-to-no genetic evidence
3. Multiple HQ "Medium" (no Strong) -> **Medium risk**: preliminary associations
4. Multiple HQ with >= 1 "Strong/Very strong" -> **Low risk**
5. Single HQ "Medium" -> **High risk** + cautionary note (needs expert review)
6. Single HQ "Strong/Very strong" -> **Low risk** + cautionary note (needs confirmation)
7. Check `study_ldPopulationStructure` for ancestry discrepancies

---

# Workflow 2: Safety Risk Assessment

Generate a 1-5 sentence narrative summarizing potential safety concerns when targeting a gene.

```python
from gsp_tools import phenotypes_for_gene

result = phenotypes_for_gene("IBD", "OSMR", cutoff="Medium", data_dir="tmp/IBD_20260407")
# Returns: gwas_results, omim_results, burden_results
```

## Interpreting Direction of Effect

| Direction | Meaning | Implication for Inhibitors |
|-----------|---------|---------------------------|
| **LOF** (Loss-of-function) | Reduced gene function causes phenotype | Inhibiting the target may cause this phenotype |
| **GOF** (Gain-of-function) | Increased gene function causes phenotype | Inhibiting the target may protect against this phenotype |

---

# Workflow 3: Combo / Multi-Gene Scoring

Score gene combinations across multiple phenotypes with ranked report generation.

## End-to-End Flow

```
┌──────────────────────────────────────────────────────────────────────┐
│ COMBO / MULTI-GENE SCORING                                           │
├──────────────────────────────────────────────────────────────────────┤
│  Step 0: RESOLVE DATASET (same as Workflow 1)                        │
│  ├── Read references/datasets.md registry                            │
│  ├── Match user's disease text to dataset + EFO IDs                  │
│  └── If ambiguous -> AskUserQuestion                                 │
│                                                                      │
│  Step 1: ALIGN ON PHENOTYPES & INDICATIONS (MANDATORY)               │
│  ├── AskUserQuestion: present available indications from registry    │
│  ├── Q1: Which phenotypes to include? (multi-select from dataset)    │
│  ├── Q2: How to aggregate? (average vs best-across)                  │
│  └── Only proceed after user confirms                                │
│                                                                      │
│  Step 2: AUTH + DOWNLOAD DATA                                        │
│  ├── aws sso login --profile cmp-dev (AskUserQuestion if auth needed)│
│  └── download_gsp_from_s3(s3_path, data_dir, aws_profile)           │
│                                                                      │
│  Step 3: PARSE INPUT (genes + combos)                                │
│  ├── Inline: extract from user prompt                                │
│  └── File: load_combos_from_file(path) -- supports TSV + JSON       │
│                                                                      │
│  Step 4: SCORE GENES                                                 │
│  ├── discover_dataset_genes(data_dir) -- pre-check existence         │
│  ├── For each gene: score_gene_across_phenotypes(...)                │
│  └── Uses biological_risk_score() per gene-phenotype pair            │
│                                                                      │
│  Step 5: SCORE COMBOS                                                │
│  ├── compute_combo_scores(gene_results, combos, phenotypes)          │
│  └── Combo score = mean of member gene final_scores                  │
│                                                                      │
│  Step 6: EXPORT                                                      │
│  ├── export_individual_tsv(...) -> gsp_{dataset}_individual.tsv      │
│  ├── export_combo_tsv(...)     -> gsp_{dataset}_combo.tsv            │
│  └── generate_combo_report_template(...) -> two .md templates        │
│                                                                      │
│  Step 7: ENHANCE REPORTS (MANDATORY -- see "Step 7" section below)   │
│  ├── 7a: Fill <!-- Claude: ... --> placeholders in both .md files    │
│  │   ├── Per-gene biomedical interpretation paragraphs               │
│  │   ├── Absent-gene explanation table                               │
│  │   ├── Combo strategic analysis by tier                            │
│  │   └── Manual adjustment recommendations                          │
│  └── 7b: Generate standalone gsp_{dataset}_combo_summary_report.md   │
│      └── Consolidated interpretive document for stakeholders         │
└──────────────────────────────────────────────────────────────────────┘
```

## Step 1: Phenotype & Indication Alignment (AskUserQuestion — BEFORE Download)

> **MANDATORY:** This step must complete **before** downloading data or running any analysis. See the universal "Mandatory Phenotype & Indication Alignment" section in the Dataset Registry for full details and templates.

When a dataset has multiple phenotypes, **always** ask the user before proceeding:

**Q1: Phenotype selection** (multi-select)
```
"Which phenotypes should be included in scoring?"
Options (example for IBD):
  - All three (IBD + UC + CD) (Recommended)
  - IBD + Crohn's only
  - IBD + UC only
  - IBD umbrella only
```

**Q2: Aggregation method**
```
"How should scores be aggregated across phenotypes?"
Options:
  - Average across all selected (Recommended)
  - Best (maximum) across selected
```

For single-phenotype datasets, still confirm coverage (do not silently assume).

## Step 3: Input Parsing

**Inline:** Claude extracts gene lists and combo definitions from the user's prompt.

**File input:** Use `load_combos_from_file()`:
```python
from gsp_tools import load_combos_from_file

# TSV format (columns: combo_id, gene)
genes, combos = load_combos_from_file("combos.tsv")

# JSON format (list of lists)
genes, combos = load_combos_from_file("combos.json")
```

## Steps 4-6: Scoring + Export

```python
import sys
sys.path.insert(0, '/home/sagemaker-user/.claude/skills/genetics-gsp')
from gsp_tools import (
    discover_dataset_genes,
    score_gene_across_phenotypes,
    compute_combo_scores,
    export_individual_tsv,
    export_combo_tsv,
    generate_combo_report_template,
)

DATA_DIR = 'tmp/IBD_20260407'
DATASET = 'IBD'

# Phenotypes as (efo_id, short_name, full_name) tuples
PHENOTYPES = [
    ('EFO_0003767', 'IBD', 'Inflammatory bowel disease'),
    ('EFO_0000729', 'UC', 'Ulcerative colitis'),
    ('EFO_0000384', 'CD', "Crohn's disease"),
]

GENES = ['TYK2', 'JAK1', 'OSMR', 'TNFSF15']
COMBOS = [['TYK2', 'JAK1'], ['OSMR', 'TNFSF15']]

# Step 4: Score genes
dataset_genes = discover_dataset_genes(DATA_DIR)
gene_results = {}
for gene in GENES:
    gene_results[gene.upper()] = score_gene_across_phenotypes(
        gene, PHENOTYPES, DATA_DIR, DATASET,
        dataset_genes, aggregation='average'
    )

# Step 5: Score combos
combo_results = compute_combo_scores(gene_results, COMBOS, PHENOTYPES)

# Step 6: Export
export_individual_tsv(gene_results, PHENOTYPES, 'gsp_IBD_individual_scores.tsv')
export_combo_tsv(combo_results, gene_results, PHENOTYPES, 'gsp_IBD_combo_scores.tsv')
reports = generate_combo_report_template(
    gene_results, combo_results, PHENOTYPES,
    dataset_name=DATASET, output_dir='.',
    aggregation='average', include_disclaimer=True
)
# reports = {'individual': 'gsp_IBD_individual_scores.md',
#            'combo': 'gsp_IBD_combo_scores.md'}
```

## Score Mapping

| GSP `risk_level` | Evidence Meaning | User Label | Score |
|---|---|---|---|
| Low risk | Strong genetic evidence | High | 100 |
| Medium risk | Moderate genetic evidence | Medium | 50 |
| High risk | Weak/no genetic evidence | Low | 25 |
| Gene absent from dataset | No data | None | 0 |
| Error/ambiguous | Cannot determine | Unclear | 0 |

**Rationale:** GSP "Low risk" = strong evidence = high confidence in target. The inversion maps program-risk to evidence-strength.

## Output Files

```
{workdir}/
├── gsp_{dataset}_individual_scores.md   # Gene-level report (Claude enhances in-place)
├── gsp_{dataset}_individual_scores.tsv  # Gene-level detail (machine-readable)
├── gsp_{dataset}_combo_scores.md        # Combo-level report (Claude enhances in-place)
├── gsp_{dataset}_combo_scores.tsv       # Combo-level detail (machine-readable)
└── gsp_{dataset}_combo_summary_report.md  # Standalone summary (Claude generates)
```

## Step 7: Enhance Reports (MANDATORY)

After `generate_combo_report_template()` produces the template MD files, Claude **MUST** perform two actions:

### 7a. Fill in template placeholders (in-place)

Both `gsp_{dataset}_individual_scores.md` and `gsp_{dataset}_combo_scores.md` contain `<!-- Claude: ... -->` placeholder comments. Claude must replace every placeholder with substantive content:

> **CRITICAL: Computed Scores Only**
>
> All final rankings and scores in reports MUST use the **computed/reported scores** from the scoring pipeline. Never replace computed scores with suggested or adjusted values in ranking tables. Interpretation paragraphs may note discrepancies and flag genes for expert review, but the score and rank columns must reflect the actual pipeline output.

**In `gsp_{dataset}_individual_scores.md`:**
- Replace `<!-- Claude: Add biomedical interpretation for each gene below -->` with per-gene interpretation paragraphs. For each gene present in the dataset, add a **Biomedical context** paragraph above its evidence table covering:
  - Gene function and pathway role in the disease context
  - Therapeutic programs targeting this gene (approved drugs, clinical trials)
  - Interpretation of why the GSP score may over/underestimate true evidence
  - Flags for expert review where discrepancies exist (do NOT propose adjusted scores)
- Add a section for genes absent from the dataset with a table explaining why each gene is missing (locus mapping, expression-only target, etc.)
- Replace `<!-- Claude: Add comprehensive interpretation narrative here -->` with:
  - Summary of key findings (dataset coverage, score distribution patterns)
  - Genes flagged for expert review with rationale (using computed scores only)
  - Ancestry and study design caveats

**In `gsp_{dataset}_combo_scores.md`:**
- Replace `<!-- Claude: Add combo-level strategic interpretation here -->` with:
  - Tier-by-tier strategic analysis of combinations (group by score tier)
  - For top-ranked combos: biological rationale for dual/triple targeting
  - For zero-score combos: explain whether zero reflects true lack of evidence or dataset coverage gaps
  - A "Notes for Expert Review" section flagging combos where known discrepancies exist (do NOT propose adjusted combo scores)

### 7b. Generate standalone summary report

After enhancing the template files, Claude must also generate a standalone `gsp_{dataset}_combo_summary_report.md` that combines and synthesizes both individual and combo results into a single interpretive document. This file should contain:

1. **Scoring methodology** -- the score mapping table and aggregation method
2. **Individual gene scores** -- summary table (ranked by computed score) + per-gene biomedical interpretation (same content as 7a but in consolidated form)
3. **Combination scores** -- ranked table (by computed score) with tier classification and strategic analysis
4. **Key observations and caveats** -- dataset coverage gaps, systematic scoring patterns, GWAS locus-vs-gene mapping limitations
5. **Notes for expert review** -- genes/combos flagged with rationale for discrepancies; interpretation may discuss why scores may under/overestimate, but all ranking tables and score columns must use computed values only

> **CRITICAL:** The standalone report must use **computed/reported scores** for all final rankings and score values. Interpretation sections may note where scores may not fully reflect published evidence or clinical validation, but must never present alternative/adjusted scores in ranking tables or as "recommended" replacements.

The standalone report provides a self-contained document for stakeholders who do not need the detailed per-phenotype supplementary tables in the template files.

**TSV exports** include per-phenotype columns:
- Score, Label, RiskLevel, RiskLabel, TotalStudies, HQStudies, MaxSampleSize
- EvidenceDistribution, CautionaryNotes, Ancestry, HQStudyIDs

## Disclaimer (Configurable)

By default, reports include a prominent disclaimer:

> **RAW SCORES -- MANUAL REVIEW REQUIRED**

Set `include_disclaimer=False` in `generate_combo_report_template()` to suppress (e.g., after manual review has been completed).

---

# Tool Reference

## Discovery & Registry

| Function | Purpose |
|----------|---------|
| `load_dataset_registry()` | Load all datasets from `references/datasets.md` |
| `resolve_dataset(user_input)` | Match user text to a dataset (case-insensitive) |
| `list_indications(data_dir)` | List EFO IDs, trait names, study/gene counts for a dataset |
| `list_available_targets(data_dir)` | List available GSP files/directories |

## Data Access

| Function | Purpose |
|----------|---------|
| `download_gsp_from_s3(s3_path, data_dir, aws_profile)` | Download GSP.pkl.gz + GSP.xlsx from S3 |
| `clear_cache()` | Clear in-memory report cache |

## Assessment

| Function | Purpose |
|----------|---------|
| `expand_gene_disease_pair(gsp_prefix, gene, disease, data_dir)` | Get study-level evidence |
| `phenotypes_for_gene(gsp_prefix, gene, cutoff, data_dir)` | Get all phenotypes for safety |
| `biological_risk_score(results)` | Apply 7-rule scoring to study results |
| `generate_biological_report(data_dir, gene, diseases, output_dir, dataset_name)` | Full MD report |

## Combo / Multi-Gene Scoring

| Function | Purpose |
|----------|---------|
| `discover_dataset_genes(data_dir)` | Set of gene symbols in GWAS Summary |
| `load_combos_from_file(file_path)` | Parse gene combos from TSV or JSON file |
| `score_gene_across_phenotypes(gene, phenotypes, data_dir, dataset_name, dataset_genes, aggregation)` | Score one gene across phenotypes |
| `compute_combo_scores(gene_results, combos, phenotypes)` | Mean of member gene scores per combo |
| `export_individual_tsv(gene_results, phenotypes, output_path)` | TSV with per-phenotype gene detail |
| `export_combo_tsv(combo_results, gene_results, phenotypes, output_path)` | TSV with per-phenotype combo detail |
| `generate_combo_report_template(gene_results, combo_results, phenotypes, dataset_name, output_dir, aggregation, include_disclaimer)` | Structured MD reports for Claude to enhance |

---

# Evidence Strength Levels

| Level | Score | Meaning |
|-------|-------|---------|
| **Very strong** | 5.0 | Very high confidence genetic evidence |
| **Strong** | 4.0 | High confidence genetic evidence |
| **Medium** | 3.0 | Moderate genetic evidence |
| **Low** | 2.0 | Weak genetic evidence |
| **Very low** | 1.0 | Very weak signal |
| **None** | 0.0 | No significant association |

# Abbreviations for population/ancestry/ethnicity

| Abbreviation | Meaning |
|--------------|---------|
| nfe | Non-Finnish European |
| eas | East Asian |
| afr | African |
| csa | Central and Southern Asian |
| amr | Native American |

---

# Setup & Environment

**Required:** Run all Python commands in the `claude_test` conda environment:
```bash
conda run -n claude_test python3 -c "from gsp_tools import ..."
```

Dependencies: `pandas`, `openpyxl`, `requests` (all pre-installed in claude_test).

# File Structure

```
genetics-gsp/
├── SKILL.md                    # This file
├── pyproject.toml              # Project config
├── gsp_tools.py                # Tool functions
├── gsp_schema.json             # Data schema
├── references/
│   └── datasets.md             # Dataset registry (YAML frontmatter)
└── data/
    ├── IBD_20260407/output/    # Downloaded GSP data
    └── ...
```

# Output Structure

```
{workdir}/
└── gsp_{dataset}/
    ├── GENE1_{dataset}_gsp_risk.md
    ├── GENE2_{dataset}_gsp_risk.md
    └── ...
```
