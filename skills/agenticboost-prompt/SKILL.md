---
name: agenticboost-prompt
version: "3.0.0"
description: Generate target evaluation documents from template. Fetches Cortellis API data, pathway analysis (KEGG, Reactome, MSigDB), PPI data (STRING, IntAct, BioGRID), and auto-populates placeholders via codex gpt xhigh. Outputs markdown prompt. Generates summary.txt with validated metadata (TA, modality, disease, MoA) from reference catalogs. Requires user to provide target, indication, AND related_diseases.
allowed-tools: Read, Write, Edit, Bash, Task, AskUserQuestion, Skill
---

# Target Evaluation Prompt Generator (Enhanced)

Generate comprehensive target evaluation documents by combining:
- **Cortellis API** data for drug/disease annotations
- **Pathway databases** (KEGG, Reactome, MSigDB) via `/pathwaydb-query`
- **PPI databases** (STRING, IntAct, BioGRID) via `/interactdb-query` (combo targets only)
- **GPT-powered placeholder population** via codex

## Invocation

```
/agenticboost-prompt
```

## Templates & Reference Files

| File | Description |
|------|-------------|
| `templates/prompt_template_v2.md` | Comprehensive 9-section output template with Reasoning Pack |
| `templates/info_template.txt` | 12 core placeholders schema |
| `templates/gpt_prompts.json` | GPT prompt definitions for auto-populated fields |
| `reference/catalog.json` | Canonical catalog for therapeutic areas, modalities, MoA, and diseases |

## Required vs Auto-Populated Fields

| Field | Source | Required? | GPT Prompt Context |
|-------|--------|-----------|-------------------|
| `target` | User | **Yes** | N/A |
| `indication` | User | **Yes** | N/A |
| `related_diseases` | User | **Yes** | N/A |
| `therapeutic_area` | summary.txt (catalog) | **Yes** (validated) | Default: GI2; user confirms from catalog |
| `modality` | summary.txt (catalog) | **Yes** (validated) | Default: SM+BIO for combo, SM for single; user confirms |
| `disease_association` | summary.txt (catalog) | **Yes** (validated) | Inferred from indication, matched to catalog |
| `moa` | summary.txt (optional) | Optional | Auto-infer then confirm; multi-valued for combos |
| `modality_primary` | summary.txt override OR GPT | Partial override | Uses summary.txt value if set, else GPT |
| `target_full_name` | GPT | Auto | "{targets}" |
| `function_summary` | GPT | Auto | "{targets} in {indication}" |
| `canonical_pathway` | GPT | Auto | "{targets} in {indication}" |
| `expression_cells` | GPT | Auto | "{targets} in {indication}" |
| `effector_cells` | GPT | Auto | "{targets} in {indication}" |
| `moa_rationale` | GPT | Auto | "{targets} in {indication} + Cortellis context" |
| `combo_opportunities` | GPT | Auto | "{targets} in {indication}" |
| `provided_sources` | Auto | Auto | Directory listing |

**Partial Override Rule**: Fields populated in `summary.txt` (therapeutic_area, modality, disease_association, moa) take precedence over GPT auto-population for corresponding template placeholders. If a summary.txt field is blank, GPT fills it.

**GPT Prompt Pattern**: All GPT prompts include both `{targets}` and `{indication}` to ensure indication-specific synthesis:
```
"For {targets} **specifically in the context of {indication}**, [task description].
Focus on {indication} pathophysiology, disease mechanisms, and relevant cell types."
```

## Overview

This skill automates the target evaluation document generation process:
1. **Interactive Input**: Prompts for target(s), indication, and related diseases
2. **Output Directory**: Creates sanitized output folder
3. **Summary Metadata (summary.txt)**: Validates therapeutic area, modality, disease association, and MoA against `reference/catalog.json` via two-pass user confirmation
4. **Parallel Data Fetching**: Queries Cortellis API, pathway databases, and PPI databases
5. **Pathway Analysis**: Identifies shared/unique pathways across targets
6. **Interaction Analysis** (combo only): Finds shortest paths and bridge proteins
7. **Combo Effect Classification**: Determines synergistic/complementary/additive effects
8. **Placeholder Population**: Auto-fills template using codex gpt xhigh (summary.txt values override where set)
9. **Document Generation**: Produces final populated **Markdown** file (not DOCX)

## Orchestration Flow

### Step 1: Gather Inputs (All Required)

Use `AskUserQuestion` to collect ALL THREE required inputs before proceeding:

```markdown
**Question 1**: "Enter target gene(s) (comma-separated for combo, e.g., 'TYK2,JAK1'):"
**Question 2**: "Enter primary indication (e.g., 'Crohn's disease'):"
**Question 3**: "Enter related diseases (comma-separated, e.g., 'ulcerative colitis, MASH, PSC'):"
```

**Input Validation**:
- If any required input is missing, use `AskUserQuestion` to prompt before proceeding
- Do NOT start data queries until all three inputs are collected
- `related_diseases` is user-defined context, NOT auto-populated from Cortellis

Parse the response:
- Single target: `["GREM1"]`
- Combo targets: `["TYK2", "JAK1"]` (split by comma)
- Related diseases: passed directly to template (user knows best which diseases share relevant pathophysiology)

### Step 2: Create Output Directory

Create folder in current working directory:
```bash
# Single target
mkdir -p ./{TARGET}_{indication}/

# Combo targets
mkdir -p ./{TARGET1}_{TARGET2}_{indication}/
```

Sanitize indication name: replace spaces with `_`, remove special characters.

### Step 2.5: Generate summary.txt (Metadata Validation)

This step validates and generates `summary.txt` in the output directory using a **two-pass confirmation flow** against `reference/catalog.json`.

**Reference catalog location**: `~/.claude/skills/agenticboost-prompt/reference/catalog.json`

#### Phase A: Auto-Infer Defaults

1. **Load catalog**: Read `reference/catalog.json` for therapeutic_areas, modalities, diseases, moa
2. **Therapeutic Area**: Default to `GI2 - Gastrointestinal and Inflammation`
3. **Disease Association**: Filter catalog diseases by the selected TA, then fuzzy-match user's `indication` to find the closest `input_name`. If no match, show filtered list for user to pick
4. **Modality**: Default `SM - Small Molecule` for single targets; default both `SM - Small Molecule` AND `BIO - Biologics` for combo targets
5. **MoA**: Optional — attempt auto-inference from target/indication context. For combos, can be multi-valued (e.g., `BsAb - Bispecific Antibody, Inh - Inhibitor`)

#### Phase B: Two-Pass User Validation

**Pass 1 — Field Table Review**: Present all inferred values as a structured table via `AskUserQuestion`:

```
| Field                | Inferred Value                              | Source    |
|----------------------|---------------------------------------------|-----------|
| Target               | CDKN2D, PCOLCE                              | User      |
| Therapeutic Area     | GI2 - Gastrointestinal and Inflammation     | Default   |
| Disease Association  | IBD - Inflammatory Bowel Disease             | Catalog   |
| Modality             | SM - Small Molecule, BIO - Biologics         | Default   |
| MoA                  | BsAb - Bispecific Antibody                   | Inferred  |
```

Ask: "Review the inferred metadata. Which fields need changes? (Select all that apply, or confirm all are correct)"

**Pass 2 — Re-ask Flagged Fields Only**: For each field the user flags:
- Show available catalog options filtered by context (e.g., diseases filtered by TA)
- Allow user to pick from catalog OR type a custom value
- If custom value is **out-of-schema**: warn with closest catalog suggestion, then require explicit confirmation before accepting

#### Out-of-Schema Handling

When a user provides a value not in `catalog.json`:
1. Warn: `"⚠ '{value}' is not in the current catalog. Closest match: '{closest_match}'."`
2. Ask via `AskUserQuestion`: accept the suggestion, keep custom value, or enter a different value
3. If user confirms custom value, accept it (no `[custom]` flag needed — user has explicitly confirmed)

#### Writing summary.txt

After validation, write the file to the output directory:

**Single target format**:
```
## CARD9

- Therapeutic area: GI2 - Gastrointestinal and Inflammation (DDU)

- Disease association: IBD - Inflammatory Bowel Disease

- Modality: SM - Small Molecule

- MoA:
```

**Combo target format** (join target names with `_`):
```
## CDKN2D_PCOLCE

- Therapeutic area: GI2 - Gastrointestinal and Inflammation (DDU)

- Disease association: IBD - Inflammatory Bowel Disease

- Modality: BIO - Biologics

- MoA: BsAb - Bispecific Antibody
```

**Multi-valued MoA** (comma-separated for combos):
```
- MoA: BsAb - Bispecific Antibody, Inh - Inhibitor
```

**Do NOT proceed to Step 3 until summary.txt is written and user has confirmed all fields.**

### Step 3: Run Data Queries

#### Phase 1: Cortellis API Queries (parallel)

For EACH target, launch a separate query:

```python
python ~/.claude/skills/cortellis/cortellis_api/scripts/cortellis_gene_query.py \
    {TARGET} --all --excel --output-dir {output_dir}
```

This produces:
- `{TARGET}_cortellis_data.json` - Complete API response
- `{TARGET}_cortellis_data.xlsx` - Excel workbook

#### Phase 2: Pathway Database Queries

Query KEGG, Reactome, MSigDB for all targets:

```python
# Uses pathwaydb-query skill
from multi_gene_analysis_async import query_multiple_genes_async

pathway_results = asyncio.run(query_multiple_genes_async(
    genes=targets,
    kegg_organism="hsa",
    msigdb_collections=["H", "C2"],
    reactome_species=9606,
    network_compatible_names=True
))
```

This produces:
- `pathways_long.csv` - Gene-pathway pairs
- `pathways_wide.csv` - One row per gene
- `pathways_summary.csv` - Pathway counts
- `pathways_pathway_centric.csv` - One row per pathway

#### Phase 3: Interaction Database Queries (Combo Only)

For multi-target evaluations, find shortest paths between targets:

```python
# Uses interactdb-query skill
from unified_query import query_shortest_paths_all_databases

interaction_results = query_shortest_paths_all_databases(
    gene_list=targets,
    species=9606,
    max_distance=50,
    export_results=True,
    output_dir=str(output_dir)
)
```

This produces:
- `{TARGET1}-{TARGET2}_string_paths.json`
- `{TARGET1}-{TARGET2}_intact_paths.json`
- `{TARGET1}-{TARGET2}_biogrid_paths.json`

### Step 4: Analyze Combo Effects (Combo Only)

For combo targets, analyze the relationship:

| Effect | Criteria |
|--------|----------|
| **Synergistic** | >5 shared pathways + direct interaction |
| **Complementary** | 2-5 shared pathways OR bridge proteins present |
| **Additive** | <2 shared pathways, no direct interaction |

Generate `combo_analysis.json` with enhanced structure:

```json
{
    "targets": ["TYK2", "JAK1"],
    "indication": "Crohn's disease",
    "effect_classification": {
        "data_effect_type": "synergistic",
        "final_effect_type": "synergistic",
        "evidence": {
            "shared_pathway_count": 8,
            "has_direct_interaction": true,
            "bridge_protein_count": 3
        }
    },
    "pathway_analysis": {
        "shared_pathways": ["JAK-STAT signaling pathway (KEGG:hsa04630)", "..."],
        "unique_pathways": {
            "TYK2": ["Type I interferon signaling (Reactome:R-HSA-909733)", "..."],
            "JAK1": ["IL-6 signaling (KEGG:hsa04066)", "..."]
        },
        "pathway_counts": {"TYK2": 45, "JAK1": 62}
    },
    "interaction_analysis": {
        "direct_interaction": true,
        "shortest_path": {
            "path": ["TYK2", "STAT1", "JAK1"],
            "hops": 2,
            "bridge_proteins": ["STAT1"]
        },
        "interaction_scores": {"STRING": 850, "IntAct": 0.72}
    },
    "gpt_interpretation": {
        "rationale": "TYK2 and JAK1 exhibit synergistic potential...",
        "mechanistic_insight": "Targeting both upstream and downstream nodes..."
    }
}
```

### Step 5: Populate Placeholders

Auto-populate via codex gpt xhigh:

```bash
codex exec --skip-git-repo-check --sandbox read-only <<EOF
For target {TARGET} in {INDICATION}, provide a concise value for: {PLACEHOLDER_NAME}
Context: {CORTELLIS_SUMMARY}
EOF
```

Save to info.txt in TSV format:
```
${target}	{value}
${indication}	{value}
...
```

### Step 6: Generate Markdown Document

Replace placeholders in template and output `.md` file:

```python
# Load template (v2 by default)
template = Path("~/.claude/skills/agenticboost-prompt/templates/prompt_template_v2.md").read_text()

# Replace placeholders
content = re.sub(r'\$\{([^}]+)\}', lambda m: replacements.get(m.group(1), f"[{m.group(1)}]"), template)

# Save
output_path.write_text(content)
```

## Placeholder Reference

### Core Placeholders (19 total in info.txt)

| Placeholder | Source | Description |
|-------------|--------|-------------|
| `${target}` | **User input** | Target gene symbol(s) |
| `${target_full_name}` | GPT | Full descriptive name |
| `${indication}` | **User input** | Primary indication |
| `${modality_primary}` | GPT | Drug modality approach |
| `${related_diseases}` | **User input** | Related diseases (user-defined context) |
| `${function_summary}` | Cortellis + GPT | Mechanism summary |
| `${canonical_pathway}` | Cortellis + GPT | Signaling pathway |
| `${expression_cells}` | Cortellis + GPT | Cell types expressing target |
| `${effector_cells}` | Cortellis + GPT | Downstream effector cells |
| `${moa_rationale}` | GPT (xhigh) | MOA rationale synthesis |
| `${combo_opportunities}` | User/GPT | Combo therapy options |
| `${provided_sources}` | Auto | Data sources reference |
| `${shared_pathways}` | pathwaydb-query | Top shared pathways (comma-separated) |
| `${unique_pathways_T1}` | pathwaydb-query | Target 1 unique pathways |
| `${unique_pathways_T2}` | pathwaydb-query | Target 2 unique pathways |
| `${interaction_path}` | interactdb-query | Shortest path between targets |
| `${bridge_proteins}` | interactdb-query | Intermediate proteins |
| `${combo_effect_type}` | Analysis | synergistic/complementary/additive |
| `${pathway_count_shared}` | Analysis | Number of shared pathways |

### Extended Placeholders (populated via GPT)

The template contains ~70 total placeholders including:
- Decision framework: `${claim_1}`, `${evidence_1a}`, `${counter_1}`
- Risk assessment: `${failure_mode_1}`, `${mitigation_1}`
- Scoring: `${w1}` through `${w5}`, `${s1}` through `${s5}`
- Analysis: `${option_A}`, `${logic_1}`, `${decision_statement}`

All remaining placeholders auto-filled via codex gpt xhigh.

## CLI Usage

```bash
python ~/.claude/skills/agenticboost-prompt/scripts/target_evaluation_gen.py \
    --targets "TYK2,JAK1" \
    --indication "Crohn's disease" \
    --related-diseases "ulcerative colitis, MASH, PSC" \
    --output-dir ./output
```

### CLI Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--targets`, `-t` | **Yes** | Comma-separated target list |
| `--indication`, `-i` | **Yes** | Primary indication |
| `--related-diseases`, `-r` | **Yes** | Related diseases (user-defined context for cross-indication synthesis) |
| `--output-dir`, `-o` | No | Output directory (default: ./{targets}_{indication}/) |
| `--template` | No | Custom markdown template path (default: prompt_template_v2.md) |

## Output Structure

```
{TARGET1}_{TARGET2}_{indication}/
├── summary.txt                             # Validated metadata (TA, modality, disease, MoA)
├── {TARGET1}_cortellis_data.json           # Cortellis API response
├── {TARGET1}_cortellis_data.xlsx           # Excel workbook
├── {TARGET2}_cortellis_data.json           # (if combo)
├── {TARGET2}_cortellis_data.xlsx
├── pathways_long.csv                       # All pathways (long format)
├── pathways_wide.csv                       # Pathways per gene
├── pathways_summary.csv                    # Pathway stats
├── pathways_pathway_centric.csv            # Pathways with gene counts
├── {TARGET1}-{TARGET2}_string_paths.json   # STRING PPI paths (combo only)
├── {TARGET1}-{TARGET2}_intact_paths.json   # IntAct PPI paths
├── {TARGET1}-{TARGET2}_biogrid_paths.json  # BioGRID PPI paths
├── combo_analysis.json                     # Enhanced combo analysis (combo only)
├── {TARGET1}_{TARGET2}_info.txt            # All placeholder values
└── {TARGET1}_{TARGET2}_({indication})_Target_Evaluation_Prompt.md  # Output
```

## Integration Points

| Skill | Script | Function |
|-------|--------|----------|
| `/pathwaydb-query` | `multi_gene_analysis_async.py` | `query_multiple_genes_async()` |
| `/interactdb-query` | `unified_query.py` | `query_shortest_paths_all_databases()` |
| `/cortellis` | `cortellis_gene_query.py` | CLI with `--all --excel` |
| `/skill-codex` | codex CLI | `codex exec` |

## Dependencies

- Existing `/cortellis` skill for API access
- Existing `/pathwaydb-query` skill for pathway data
- Existing `/interactdb-query` skill for PPI data
- `/skill-codex` for GPT integration
- pandas for data manipulation

## Examples

### Single Target

```
/agenticboost-prompt
> Enter targets: GLP2R
> Enter indication: IBD
> Enter related diseases: Crohn's disease, ulcerative colitis, short bowel syndrome

Output: GLP2R_IBD/GLP2R_(IBD)_Target_Evaluation_Prompt.md
```

### Combo Targets

```
/agenticboost-prompt
> Enter targets: TYK2,JAK1
> Enter indication: Crohn's disease
> Enter related diseases: ulcerative colitis, MASH, PSC, ankylosing spondylitis

Output: TYK2_JAK1_Crohns_disease/TYK2_JAK1_(Crohns_disease)_Target_Evaluation_Prompt.md
```

## Changelog

### v3.0.0 (2026-04-29)
- **New**: Added Step 2.5 — `summary.txt` generation with two-pass metadata validation
- **New**: Created `reference/catalog.json` with canonical schemas for therapeutic areas, modalities, MoA, and diseases
- **New**: Therapeutic area, modality, disease association matched against catalog; out-of-schema values supported with warning + confirmation
- **New**: MoA field is optional, auto-inferred then confirmed; supports multi-valued entries for combos (comma-separated)
- **New**: Partial override rule — summary.txt values override GPT auto-population for modality, disease, TA fields
- **New**: Disease catalog includes TA auto-assignment (GI2, ONC, NS) with disease filtering by selected TA
- **New**: MoA catalog with 21 active entries across 6 groups (Receptor Pharmacology, Multi-specific Antibody, Antibody Conjugate, TPD, Nucleic Acid, Standalone)
- **Breaking**: Workflow now gates data queries behind summary.txt confirmation (Step 2.5 must complete before Step 3)
- **Removed**: Deprecated MoA entry "Bispecific - Bi" (replaced by "BsAb - Bispecific Antibody")

### v2.1.0 (2026-03-16)
- **Breaking**: `--related-diseases` is now a **required** CLI argument (not auto-populated from Cortellis)
- **New**: Added Templates Reference section documenting v2 templates
- **New**: Added Required vs Auto-Populated Fields table
- **Updated**: Default template is now `prompt_template_v2.md` (comprehensive 9-section format)
- **Updated**: Step 1 now collects ALL THREE required inputs (target, indication, related_diseases)
- **Improved**: GPT prompts include indication-specific context for better synthesis

### v2.0.0 (2026-02-22)
- **Breaking**: Output changed from DOCX to Markdown
- **Breaking**: Removed `--skip-cortellis` and `--disable-gpt` debug flags
- **New**: Integrated `/pathwaydb-query` for comprehensive pathway analysis
- **New**: Integrated `/interactdb-query` for PPI data (combo targets)
- **New**: Enhanced `combo_analysis.json` with real data + GPT interpretation
- **New**: Added 7 new placeholders for pathway/interaction data
- **Improved**: Effect classification now based on quantitative pathway/interaction evidence
