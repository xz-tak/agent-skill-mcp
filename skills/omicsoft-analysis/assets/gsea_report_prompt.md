# GSEA Team-Based Coverage Analysis Prompt

## Purpose

This prompt orchestrates a **team of specialized agents** using a **sequential 6-phase workflow** to generate a **coverage-focused** pathway analysis report. Unlike individual pathway deep dives, this prompt focuses on:

1. **Module coverage analysis** - What % of disease signatures are covered by target?
2. **Target-disease overlap** - Which target-associated pathways are significantly enriched in disease?
3. **Gap identification** - Target pathways NOT significantly enriched in disease = combo opportunities
4. **Safety prediction** - Target-associated pathways overlapping with safety-critical modules

---

## Context Loading (FIRST STEP)

**Before starting any analysis**, read the context file to obtain paths and parameters:

1. Read `reports/report_context.json` in the output directory
2. Extract:
   - `target`: Target gene name (e.g., "GREM1", "IL11")
   - `indication`: Disease indication (e.g., "IBD", "Fibrosis")
   - `output_dir`: Base output directory containing GSEA results
   - `reports_dir`: Where to save all report outputs
   - `input_files`: Paths to GSEA result files
     - `pathwaydb_gsea`: Target-associated pathways (e.g., `gsea_GREM1_IBD_pathwaydb.csv`)
     - `pathwaydb_gsea_fdr`: Filtered pathways FDR < 0.05
     - `gsea_all_results`: General disease pathways

**Example context file:**

```json
{
  "target": "GREM1",
  "indication": "IBD",
  "output_dir": "/path/to/GREM1_IBD_omicsoft_20260227",
  "reports_dir": "/path/to/GREM1_IBD_omicsoft_20260227/reports",
  "input_files": {
    "pathwaydb_gsea": "gsea_GREM1_IBD_pathwaydb.csv",
    "pathwaydb_gsea_fdr": "gsea_GREM1_IBD_pathwaydb_fdr0.05.csv",
    "gsea_all_results": "gsea_all_results.csv"
  }
}
```

All subsequent analysis uses these values. File paths in `input_files` are relative to `output_dir`.

---

## Input Files (CRITICAL)

### Required Files

| File | Used By | Description |
|------|---------|-------------|
| `{indication}_pathwaydb_gsea.csv` | A, B, C, D | **Universal pathway GSEA**: All pathways from combined_pathways.csv enriched across disease comparisons |
| `target_summary.csv` | E | **Target gene summary**: DEG status for target gene across all comparisons |
| `{target}_pathwaydb.csv` | B | **Target-associated pathways**: Pathways containing {target} (from /pathwaydb-query output) |

### File Locations

```
{indication}_omicsoft_{timestamp}/
├── {indication}_pathwaydb_gsea.csv               # Universal pathway GSEA (Sections A, B, C, D)
├── target_summary.csv                            # Target gene DEG summary (Section E)
└── reports/                                      # Output directory
    ├── report_context.json                       # Context for Claude (target, indication, paths)
    ├── disease_analysis.md                       # Section A output
    ├── target_analysis.md                        # Section B output
    ├── combination_analysis.md                   # Section C output
    ├── safety_analysis.md                        # Section D output
    ├── target_ranking_report.md                  # Section E output
    └── {target}_{indication}_team_report.md      # Final integrated report
```

### Target-Associated Pathways File

**Auto-detect**: Look for `{target}_pathwaydb_*/{target}_pathwaydb.csv` in workdir.
If not found: Use AskUserQuestion to prompt user for path.

### Key Distinction

- **`{indication}_pathwaydb_gsea.csv`**: Contains enrichment for ALL pathways (KEGG, Reactome, MSigDB) from the universal pathway database. Used to characterize overall disease biology (Section A) and calculate target-pathway overlap (Section B).

- **`{target}_pathwaydb.csv`**: From /pathwaydb-query skill output. Contains pathways that include {target} gene. Used to identify which pathways in the universal GSEA are target-associated.

- **`target_summary.csv`**: DEG statistics for target gene across all comparisons. Used for target scoring (Section E).

---

## Key Definitions (CRITICAL)

### Target-Associated Pathway

A pathway is **"target-associated"** if it contains the target gene in its gene set (i.e., it appears in `gsea_{target}_pathwaydb.csv`). This is about pathway membership, NOT about UP/DN direction.

### Enrichment Direction (UP/DN)

Enrichment direction refers to how a pathway behaves in a **disease comparison**:

- **UP** (NES > 0): Pathway is enriched/activated in the comparison (e.g., higher in Disease vs. Normal)
- **DN** (NES < 0): Pathway is suppressed in the comparison (e.g., lower in Disease vs. Normal)

**IMPORTANT**: UP/DN describes the pathway's behavior in disease, not the target's effect on the pathway.

### Significance Threshold

A pathway is only considered **"significantly enriched"** when **FDR/padj < 0.05**

### Coverage Calculation

```
Coverage % = (Target pathways significantly enriched in disease with FDR < 0.05) / (ALL target pathways) × 100
```

**Example**: "93% coverage" = 93% of target-associated pathways are significantly enriched (FDR < 0.05) in Disease vs. Normal

**ALWAYS separate UP vs DN for interpretation**:
- **UP coverage** = Target pathways significantly UP enriched in disease (NES > 0, FDR < 0.05)
- **DN coverage** = Target pathways significantly DN enriched in disease (NES < 0, FDR < 0.05)
- Report BOTH separately for interpretation

### Gap Definition

**Gap** = Target-associated pathways that are **NOT significantly enriched** in disease (FDR ≥ 0.05)

These represent modules where the target may not have strong disease relevance, or where combination therapy could add value.

### Safety Risk Definition

**Safety Risk** = Target pathway **significantly enriched** (UP or DN, FDR < 0.05) in disease **AND** the pathway/module is **associated with known safety concerns** (Hematologic, Hepatic, Cardiac, etc.)

---

### Output

All reports saved to: `{folder}/reports/`

| Output File | Description |
|-------------|-------------|
| `disease_analysis.md` | Disease module characterization (Section A) |
| `target_analysis.md` | Target-disease coverage analysis (Section B) |
| `target_overlap_summary.csv` | Target-pathway overlap metrics by group (Section B) |
| `target_overlap_detailed.csv` | Pathway-level overlap flags (Section B) |
| `combination_analysis.md` | Gap analysis and combo opportunities (Section C) |
| `safety_analysis.md` | Safety module overlap assessment (Section D) |
| `target_ranking_report.md` | Target gene scoring summary (Section E - parallel) |
| `target_scores.csv` | Raw target scores for downstream analysis (Section E) |
| `{target}_{indication}_team_report.md` | Final integrated report |

---

## Comparison Category Configuration

### Included Categories (3 Interpretable Types)

**FILTER**: Only process the following 3 comparison categories:

| Category | Include | Interpretation Context |
|----------|---------|------------------------|
| **Disease vs. Normal** | YES | Core disease signature - pathways UP/DN in disease state |
| **Treatment vs. Control** | YES | Treatment response - pathways affected by drug treatment |
| **Responder vs. Non-Responder** | YES | Response prediction - pathways that predict treatment response |

### Excluded Categories

| Category | Exclude | Reason |
|----------|---------|--------|
| Other Comparisons | EXCLUDE | Too vague, no clear biological context |
| CellType1 vs. CellType2 | EXCLUDE | Cell type comparisons, not disease/treatment relevant |
| Tissue1 vs. Tissue2 | EXCLUDE | Tissue comparisons, not disease/treatment relevant |

### Critical Requirement

**ALL analyses MUST**:
1. **Filter** to only the 3 included comparison categories
2. **Stratify** by `comparison_category`, `tissue`, `disease`
3. **Apply context-aware interpretation** based on category type

Never aggregate across categories without explicit breakdown.

---

## Context-Aware Interpretation by Comparison Category

### Disease vs. Normal Interpretation

**What UP/DN means**:
- UP (NES > 0): Pathway is MORE active in disease than normal
- DN (NES < 0): Pathway is LESS active in disease than normal

**Target Effect Interpretation**:
- Target pathways UP in disease → Target is involved in disease pathology
  - If therapeutic goal is to INHIBIT target → May reduce disease activity
  - If therapeutic goal is to ACTIVATE target → May worsen disease
- Target pathways DN in disease → Target is involved in protective mechanisms
  - If therapeutic goal is to ACTIVATE target → May restore protective function
  - If therapeutic goal is to INHIBIT target → May worsen suppression

**Combo Opportunity Interpretation**:
- Gap in UP modules: Disease pathways activated but target doesn't cover → Combo needed to suppress
- Gap in DN modules: Disease pathways suppressed but target doesn't restore → Combo needed to activate

**Safety Risk Interpretation**:
- Target pathways UP in disease + Safety module: May worsen safety concern if target activated in disease
- Target pathways DN in disease + Safety module: May worsen safety concern if target suppressed in disease

### Responder vs. Non-Responder Interpretation

**What UP/DN means**:
- UP (NES > 0): Pathway is MORE active in responders than non-responders
- DN (NES < 0): Pathway is MORE active in NON-responders than responders

**Target Effect Interpretation**:
- Target pathways UP in responders → Target activity predicts response
  - Patients with high target pathway activity → More likely to respond
  - Use as patient selection biomarker
- Target pathways DN (UP in non-responders) → Target activity predicts non-response
  - Pathway is a resistance mechanism
  - Consider combo to overcome resistance

**Combo Opportunity Interpretation**:
- Target pathways DN (UP in non-responders) → Resistance mechanisms
  - Combo target should address these pathways to improve response rate
- Target pathways not enriched in non-responders → Not a resistance mechanism

**Safety Risk Interpretation**:
- Target + Safety module both UP in responders: Safety concern correlates with response
  - Responders may have MORE toxicity (therapeutic index consideration)
- Target + Safety module both DN (UP in non-responders): Non-responders may have safety signals

### Treatment vs. Control Interpretation

**What UP/DN means**:
- UP (NES > 0): Pathway is ACTIVATED by treatment
- DN (NES < 0): Pathway is SUPPRESSED by treatment

**Target Effect Interpretation**:
- Target pathways UP after treatment → Treatment activates target pathways
  - Mechanism of action (MOA) evidence
  - On-target pharmacodynamic effect
- Target pathways DN after treatment → Treatment suppresses target pathways
  - MOA evidence for inhibitory mechanism
  - On-target effect

**Combo Opportunity Interpretation**:
- Disease pathways NOT affected by treatment → Treatment doesn't address these
  - Combo needed to cover disease pathways unaffected by treatment
- Disease UP pathways still UP after treatment → Incomplete treatment response
  - Combo may enhance effect

**Safety Risk Interpretation**:
- Safety module UP after treatment → Treatment activates safety-associated pathways
  - On-target toxicity signal
  - Monitor for this adverse effect
- Safety module DN after treatment → Treatment suppresses safety-associated pathways
  - May cause deficiency-related toxicity

---

## Team Architecture

### Agent Roles

| Agent | Role | Focus Area | Input File | Model |
|-------|------|------------|------------|-------|
| **Planner/Reviewer** | Orchestrates team, reviews outputs | Quality control, integration | All files | Opus |
| **Disease Module Analyst (A)** | Module-level disease characterization | Coverage summaries, NOT pathway details | {indication}_pathwaydb_gsea.csv | Opus |
| **Target Coverage Specialist (B)** | Target-disease overlap analysis | Coverage %, overlap matrices | {indication}_pathwaydb_gsea.csv + {target}_pathwaydb.csv | Opus |
| **Combination Strategist (C)** | Gap analysis for combinations | Modules NOT covered by target | A + B results | Opus |
| **Safety Analyst (D)** | Safety module overlap assessment | Target ↔ safety-critical module overlap | {indication}_pathwaydb_gsea.csv | Opus |
| **Target Scoring Specialist (E)** | Target gene scoring from DEG | Rank targets by disease relevance | target_summary.csv | Opus |

### Team Configuration

```json
{
  "team_name": "gsea-{target}-{indication}",
  "description": "Coverage-focused pathway analysis for {target} in {indication}",
  "members": [
    {
      "name": "planner",
      "subagent_type": "general-purpose",
      "model": "opus",
      "role": "Orchestrator and quality reviewer"
    },
    {
      "name": "disease-analyst",
      "subagent_type": "general-purpose",
      "model": "opus",
      "role": "Disease module coverage analysis (Agent A)"
    },
    {
      "name": "target-specialist",
      "subagent_type": "general-purpose",
      "model": "opus",
      "role": "Target-disease overlap and coverage analysis (Agent B)"
    },
    {
      "name": "combo-strategist",
      "subagent_type": "general-purpose",
      "model": "opus",
      "role": "Gap analysis and combination identification (Agent C)"
    },
    {
      "name": "safety-analyst",
      "subagent_type": "general-purpose",
      "model": "opus",
      "role": "Safety module overlap assessment (Agent D)"
    },
    {
      "name": "target-scorer",
      "subagent_type": "general-purpose",
      "model": "opus",
      "role": "Target gene scoring from DEG (Agent E - parallel track)"
    }
  ]
}
```

---

## 6-Phase Sequential Workflow with Parallel Target Scoring

### Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  TRACK 1 (Sequential): A → B → C → D                                        │
├─────────────────────────────────────────────────────────────────────────────┤
│                    PHASE 1: FOUNDATION                                      │
│         ┌─────────────┐       ┌─────────────┐                               │
│         │  Agent A    │       │  Agent B    │  (Parallel)                   │
│         │  Disease    │       │  Coverage   │                               │
│         └──────┬──────┘       └──────┬──────┘                               │
│                │                     │                                      │
│                ▼                     ▼                                      │
│                    PHASE 2: FOUNDATION REVIEW                               │
│                      Planner Reviews A + B                                  │
│                 ┌──────────────────────────┐                                │
│                 │  A APPROVED? B APPROVED? │                                │
│                 └────────────┬─────────────┘                                │
│              ┌───────────────┴───────────────┐                              │
│              ▼                               ▼                              │
│         [NO: Revise]                   [YES: Proceed]                       │
│                                               │                             │
│                                               ▼                             │
│                    PHASE 3: DEPENDENT ANALYSIS                              │
│         ┌─────────────┐       ┌─────────────┐                               │
│         │  Agent C    │       │  Agent D    │  (Parallel)                   │
│         │  Combo/Gap  │       │  Safety     │                               │
│         │ (uses A+B)  │       │ (uses A+B)  │                               │
│         └──────┬──────┘       └──────┬──────┘                               │
│                │                     │                                      │
│                ▼                     ▼                                      │
│                    PHASE 4: DEPENDENT REVIEW                                │
│                      Planner Reviews C + D                                  │
│    ┌─────────────────────────────────────────────────────┐                  │
│    │  C/D reveal A/B issues?  │  C/D APPROVED?           │                  │
│    └───────────┬──────────────┴───────────┬──────────────┘                  │
│                ▼                          ▼                                 │
│         [YES: Phase 5]            [YES: Phase 6]                            │
│                │                                                            │
│                ▼                                                            │
│                    PHASE 5: ITERATION CYCLE                                 │
│    1. Update A/B → Review A/B → Approve A/B                                 │
│    2. Update C/D with new A/B → Review C/D                                  │
│    3. Full team review (A+B+C+D consistency)                                │
│    4. Repeat until all APPROVED (max 3 cycles)                              │
│                              │                                              │
│                              ▼                                              │
│                    [All APPROVED → Phase 6]                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│  TRACK 2 (Parallel): E                                                      │
├─────────────────────────────────────────────────────────────────────────────┤
│                    RUNS INDEPENDENTLY                                       │
│                    ┌─────────────┐                                          │
│                    │  Agent E    │                                          │
│                    │  Target     │                                          │
│                    │  Scoring    │                                          │
│                    └──────┬──────┘                                          │
│                           │                                                 │
│                    Uses: target_summary.csv                                 │
│                    Output: target_ranking_report.md                         │
├─────────────────────────────────────────────────────────────────────────────┤
│  PHASE 6: FINAL ASSEMBLY (Track 1 + Track 2)                                │
├─────────────────────────────────────────────────────────────────────────────┤
│    Compile: Executive Summary + A + B + C + D + E + Appendix                │
│    Save: {target}_{indication}_team_report.md                               │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

### Phase 1: Foundation Analysis (A + B in Parallel)

**Planner spawns 2 agents in parallel:**
- **Disease Analyst (Agent A)** → Section A: Disease Module Characterization
- **Target Specialist (Agent B)** → Section B: Target Coverage Analysis

These are foundational analyses that C and D depend on.

**CRITICAL**: Wait for BOTH A and B to complete before proceeding to Phase 2.

---

#### 1.1 Disease Module Analyst Task (Agent A)

```markdown
## Task: Disease Module Coverage Analysis

**Input File**: {indication}_pathwaydb_gsea.csv (universal pathway GSEA from combined_pathways.csv)

**Instructions**:

### Step 0: Filter to Interpretable Comparison Categories (CRITICAL)

**BEFORE any analysis**, filter data to ONLY these 3 categories:
1. Disease vs. Normal
2. Treatment vs. Control
3. Responder vs. Non-Responder

**EXCLUDE**: Other Comparisons, CellType1 vs. CellType2, Tissue1 vs. Tissue2

```python
# Filter code
included_categories = [
    'Disease vs. Normal',
    'Treatment vs. Control',
    'Responder vs. Non-Responder'
]
df_filtered = df[df['comparison_category'].isin(included_categories)]

# Only include FDR < 0.05 pathways
df_significant = df_filtered[df_filtered['FDR q-val'] < 0.05]
```

### Step 1: Module Classification

Classify ALL significant pathways (FDR < 0.05) into biological modules:

| Module | Example Pathways | Typical Sources |
|--------|------------------|-----------------|
| Immune/Inflammatory | IL-6 signaling, TNF pathway, Complement | HALLMARK, KEGG, Reactome |
| Cell Cycle | G1/S transition, DNA replication | KEGG, Reactome |
| ECM/Fibrosis | Collagen formation, TGF-beta | Reactome, KEGG |
| Metabolic | Oxidative phosphorylation, Glycolysis | KEGG, Reactome |
| Apoptosis/Cell Death | p53 pathway, Caspase cascade | HALLMARK, Reactome |
| Signaling | MAPK, PI3K-AKT, Wnt | KEGG |
| DNA Damage/Repair | Mismatch repair, Base excision | Reactome |

### Step 2: Generate A1 - Global Module Enrichment Summary

**CRITICAL**: Full breakdown by Module × Comparison × Tissue × Disease with N_UP and N_DN counts (FDR < 0.05 only)

**Output Table A1: Global Module Enrichment by Comparison × Tissue × Disease**

| Module | comparison_category | Tissue | Disease | N_UP (FDR<0.05) | N_DN (FDR<0.05) | Avg_NES_UP | Avg_NES_DN |
|--------|---------------------|--------|---------|-----------------|-----------------|------------|------------|
| Immune/Inflammatory | Disease vs. Normal | Blood | UC | 156 | 23 | +2.4 | -1.8 |
| Immune/Inflammatory | Disease vs. Normal | Colon | UC | 89 | 45 | +1.9 | -2.1 |
| ECM/Fibrosis | Responder vs. Non-Responder | Colon | CD | 12 | 67 | +1.5 | -2.8 |
| ... | ... | ... | ... | ... | ... | ... | ... |

**Interpretation Guide for A1**:
- **Disease vs. Normal**: N_UP = pathways activated in disease; N_DN = pathways suppressed in disease
  - High N_UP in Immune module → inflammation is UP in disease
  - High N_DN in DNA Repair → DNA repair capacity DOWN in disease
- **Responder vs. Non-Responder**: N_UP = pathways higher in responders; N_DN = pathways higher in non-responders
  - High N_DN in ECM → non-responders have MORE fibrosis (ECM UP in non-responders)
- **Treatment vs. Control**: N_UP = pathways activated by treatment; N_DN = pathways suppressed by treatment

### Step 3: Generate A2 - Global UP-Enriched Clusters/Networks

**Output Table A2: Pathways Significantly UP-Enriched Across All Conditions**

| Module | Total UP Pathways | Key Biological Themes | Dominant Comparison | Dominant Tissue |
|--------|-------------------|----------------------|---------------------|-----------------|
| Immune/Inflammatory | 1,234 | Cytokines, Complement, T-cell activation | Disease vs. Normal | Blood, Colon |
| Cell Cycle | 456 | Proliferation, G1/S, DNA replication | Disease vs. Normal | Colon |
| Apoptosis | 234 | Cell death, Caspases, p53 signaling | Treatment vs. Control | Colon |

**Interpretation for UP Clusters**:
- **Target Effect**: If target pathways overlap with UP clusters → target modulation affects activated disease pathways
- **Combo Opportunity**: UP modules NOT covered by target → potential combo targets
- **Safety Risk**: UP pathways in safety modules (Hematologic, Hepatic) → potential toxicity from target activation

### Step 4: Generate A3 - Global DN-Enriched Clusters/Networks

**Output Table A3: Pathways Significantly DN-Enriched Across All Conditions**

| Module | Total DN Pathways | Key Biological Themes | Dominant Comparison | Dominant Tissue |
|--------|-------------------|----------------------|---------------------|-----------------|
| DNA Repair | 567 | BRCA1, ATM, Mismatch repair | Disease vs. Normal | Blood |
| Metabolic | 345 | OXPHOS, Glycolysis, Lipid metabolism | Disease vs. Normal | Colon |
| Treg/Immune regulation | 123 | FOXP3, Regulatory T-cells | Disease vs. Normal | Blood |

**Interpretation for DN Clusters**:
- **Target Effect**: If target pathways overlap with DN clusters → target modulation affects suppressed disease pathways
- **Combo Opportunity**: DN modules NOT covered by target → restore suppressed functions
- **Safety Risk**: DN pathways in safety modules (Immune) → further suppression may worsen immune deficiency

### Step 5: Context-Aware Interpretation by Category

**CRITICAL**: Interpret findings based on comparison category:

| comparison_category | Module UP Interpretation | Module DN Interpretation |
|---------------------|--------------------------|--------------------------|
| Disease vs. Normal | Active in disease pathology | Suppressed in disease |
| Treatment vs. Control | Activated by treatment | Suppressed by treatment |
| Responder vs. Non-Responder | Marker of response | Marker of non-response (UP in non-responders) |

**DO NOT**:
- List individual pathways by name
- Provide deep dives on specific pathways
- Include top 15/10/5 pathway rankings
```

---

#### 1.2 Target Coverage Specialist Task (Agent B)

```markdown
## Task: Target-Disease Coverage Analysis (Overlap-Based)

**Input Files**:
- `{indication}_pathwaydb_gsea.csv` (universal pathway GSEA results)
- `{target}_pathwaydb.csv` (target-associated pathways from /pathwaydb-query)

**Instructions**:

### Step 0: Locate Target Pathways File (AUTO-DETECT)

1. Search for: `{target}_pathwaydb_*/{target}_pathwaydb.csv` in workdir
2. If found: Load and proceed
3. If NOT found: Use AskUserQuestion to ask user for path:
   - "Target pathways file not found. Please provide path to {target}_pathwaydb.csv"

```python
import glob
target_file_pattern = f"{target}_pathwaydb_*/{target}_pathwaydb.csv"
matches = glob.glob(target_file_pattern)
if matches:
    target_pathways_file = matches[0]
else:
    # Prompt user for path
    pass
```

### Step 1: Filter to Interpretable Comparison Categories (CRITICAL)

**BEFORE any analysis**, filter data to ONLY these 3 categories:
1. Disease vs. Normal
2. Treatment vs. Control
3. Responder vs. Non-Responder

**EXCLUDE**: Other Comparisons, CellType1 vs. CellType2, Tissue1 vs. Tissue2

```python
# Filter code
included_categories = [
    'Disease vs. Normal',
    'Treatment vs. Control',
    'Responder vs. Non-Responder'
]
df_gsea = df_gsea[df_gsea['comparison_category'].isin(included_categories)]
```

### Step 2: Load and Match Target Pathways

Load target pathways and match with universal GSEA results:

```python
# Load target pathways
df_target_pathways = pd.read_csv(target_pathways_file)
target_pathway_names = df_target_pathways['Pathway_Name'].unique()

# Mark target-associated pathways in universal GSEA
df_gsea['is_target_pathway'] = df_gsea['Term'].isin(target_pathway_names)
```

### Step 3: Calculate Overlap Metrics

**% UP Overlap** = (Target pathways significantly UP in disease) / (All target pathways) × 100
**% DN Overlap** = (Target pathways significantly DN in disease) / (All target pathways) × 100

```python
# Calculate overlap by group
def calc_overlap(group):
    target_only = group[group['is_target_pathway']]
    total_target = len(target_only)
    sig_up = len(target_only[(target_only['FDR q-val'] < 0.05) & (target_only['NES'] > 0)])
    sig_dn = len(target_only[(target_only['FDR q-val'] < 0.05) & (target_only['NES'] < 0)])
    return {
        'Total_Target_Pathways': total_target,
        'Sig_UP': sig_up,
        'Sig_DN': sig_dn,
        'UP_Overlap_%': (sig_up / total_target * 100) if total_target > 0 else 0,
        'DN_Overlap_%': (sig_dn / total_target * 100) if total_target > 0 else 0
    }

overlap_summary = df_gsea.groupby(['comparison_category', 'tissue', 'disease']).apply(calc_overlap)
```

### Step 4: Understand Coverage Definition

**CRITICAL**: Coverage is about target pathways being significantly enriched in disease, NOT about target UP/DN direction.

**Coverage Formula**:
```
Coverage % = (Target pathways significantly enriched in disease with FDR < 0.05) / (ALL target pathways) × 100
```

**Separate by enrichment direction**:
- **UP Coverage** = (Target pathways with FDR < 0.05 AND NES > 0 in disease) / (ALL target pathways)
- **DN Coverage** = (Target pathways with FDR < 0.05 AND NES < 0 in disease) / (ALL target pathways)

### Step 2: Calculate Coverage Metrics by Module

**Output Table B1: Target Coverage by Module**

| Module | Total Target Pathways | Sig UP (FDR<0.05, NES>0) | Sig DN (FDR<0.05, NES<0) | UP Coverage % | DN Coverage % | Total Coverage % |
|--------|----------------------|--------------------------|--------------------------|---------------|---------------|------------------|
| Immune/Inflammatory | 200 | 156 | 23 | 78% | 12% | 90% |
| Cell Cycle | 100 | 45 | 12 | 45% | 12% | 57% |
| ECM/Fibrosis | 80 | 8 | 4 | 10% | 5% | 15% |
| Metabolic | 60 | 2 | 1 | 3% | 2% | 5% |
| **TOTAL** | **440** | **211** | **40** | **48%** | **9%** | **57%** |

**Interpretation**:
- High UP coverage (>70%) in a module → Target pathways are significantly activated in disease for this module
- High DN coverage (>70%) in a module → Target pathways are significantly suppressed in disease for this module
- Low coverage (<30%) → Target pathways not strongly enriched in disease → potential GAP for combo

### Step 3: Calculate Coverage by Comparison Category

**Output Table B2: Coverage by Comparison Category**

| comparison_category | Total Target Pathways | Sig UP | Sig DN | UP Coverage % | DN Coverage % | Total Coverage % | Clinical Implication |
|---------------------|----------------------|--------|--------|---------------|---------------|------------------|----------------------|
| Disease vs. Normal | 440 | 211 | 40 | 48% | 9% | 57% | Target addresses core disease |
| Responder vs. Non-Responder | 440 | 45 | 89 | 10% | 20% | 30% | Target pathways relate to response |
| Treatment vs. Control | 440 | 123 | 67 | 28% | 15% | 43% | Treatment modulates target pathways |

### Step 4: Calculate Coverage by Tissue

**Output Table B3: Tissue-Specific Coverage**

| Tissue | Total Target Pathways | Sig UP | Sig DN | UP Coverage % | DN Coverage % | Total Coverage % |
|--------|----------------------|--------|--------|---------------|---------------|------------------|
| Blood | 440 | 189 | 34 | 43% | 8% | 51% |
| Colon | 440 | 156 | 45 | 35% | 10% | 45% |
| Ileum | 440 | 78 | 23 | 18% | 5% | 23% |

### Step 5: Summarize Coverage Implications

**Coverage Interpretation Matrix**:

| Coverage Pattern | Meaning | Strategic Implication |
|-----------------|---------|----------------------|
| High UP coverage (>70%) | Target pathways mostly activated in disease | Target inhibition may reduce disease activity |
| High DN coverage (>70%) | Target pathways mostly suppressed in disease | Target activation may restore function |
| Low coverage (<30%) | Target pathways not strongly enriched | Target may not be primary disease driver |
| Mixed (UP ~ DN) | Target pathways show bidirectional changes | Need module-specific analysis |

### Step 6: Context-Aware Coverage Interpretation

**CRITICAL**: Interpret coverage based on comparison category:

| comparison_category | High UP Coverage Meaning | High DN Coverage Meaning |
|---------------------|-------------------------|-------------------------|
| **Disease vs. Normal** | Target pathways activated in disease → Target inhibition = benefit | Target pathways suppressed in disease → Target activation = benefit |
| **Treatment vs. Control** | Treatment activates target pathways | Treatment suppresses target pathways |
| **Responder vs. Non-Responder** | Target pathway = response marker | Target pathway = resistance marker |

**DO NOT**:
- Use "target UP/DN" language (target doesn't have UP/DN direction)
- List individual pathways by name
- Provide detailed pathway mechanisms
```

---

### Phase 2: Foundation Review (Planner Reviews A + B)

**Planner reviews Section A and Section B using review checklist.**

**Review Checklist for A:**
- [ ] A1 has full breakdown: Module × Comparison × Tissue × Disease
- [ ] A1 includes N_UP and N_DN columns (FDR < 0.05 only)
- [ ] A2 covers global UP-enriched clusters with themes
- [ ] A3 covers global DN-enriched clusters with themes

**Review Checklist for B:**
- [ ] Coverage formula: (Target pathways FDR < 0.05) / (ALL target pathways)
- [ ] UP coverage and DN coverage reported separately
- [ ] Breakdown by module and comparison category

**If A or B FAILS review:**
1. Send critique to agent with specific issues
2. Agent revises and resubmits
3. Planner reviews again
4. Repeat until APPROVED

**CRITICAL**: Only proceed to Phase 3 when BOTH A and B are APPROVED.

#### Critique Template (Team Lead Feedback)

```markdown
## Section [A/B] Review

**Status**: NEEDS REVISION / APPROVED

### Issues Found:
1. [Specific issue 1]
2. [Specific issue 2]

### Missing Elements:
- [ ] [Missing element 1]
- [ ] [Missing element 2]

### Required Fixes:
1. [What to fix and how]
2. [What to add]

### Next Step:
Revise and resubmit. Do NOT proceed until issues resolved.
```

---

### Phase 3: Dependent Analysis (C + D in Parallel)

**Prerequisites:** Section A and Section B must be APPROVED.

**Planner spawns 2 agents in parallel, providing approved A/B results:**

**Agent Spawn Prompt MUST Include:**
1. Summary of approved A findings (modules, N_UP, N_DN by comparison/tissue/disease)
2. Summary of approved B findings (coverage % by module, UP vs DN coverage)
3. Any planner feedback or context from Phase 2 review

---

#### 3.1 Combination Strategist Task (Agent C)

```markdown
## Task: Gap Analysis for Combination Opportunities

**Input Files**:
- Approved Section A results (disease modules from gsea_all_results.csv)
- Approved Section B results (target coverage from gsea_{target}_{indication}_pathwaydb.csv)
- Planner feedback/context from Phase 2

**Instructions**:

### Step 0: Filter to Interpretable Comparison Categories (CRITICAL)

**BEFORE any analysis**, ensure data is filtered to ONLY these 3 categories:
1. Disease vs. Normal
2. Treatment vs. Control
3. Responder vs. Non-Responder

**EXCLUDE**: Other Comparisons, CellType1 vs. CellType2, Tissue1 vs. Tissue2

### Step 1: Understand Gap Definition

**CRITICAL**: Gap = Target pathways NOT significantly enriched in disease (FDR ≥ 0.05)

**What is "Covered"?**
A target-associated pathway is "covered" in disease when it is **significantly enriched (FDR < 0.05)** in disease.

**Gap Categories**:
1. **Within-module gaps**: Target pathways in a module NOT significantly enriched (FDR ≥ 0.05) in disease
2. **Cross-module gaps**: Disease modules where target has very low coverage (<30%)
3. **Disease-specific gaps**: Disease pathways NOT overlapped with any target pathway

### Step 2: Identify Gap Modules

**Output Table C1: Gap Analysis - Target Pathways NOT Significantly Enriched in Disease**

| Gap Module | Total Target Pathways | Covered (FDR<0.05) | Gap (FDR≥0.05) | Gap % | Dominant Disease Direction | Combo Opportunity |
|------------|----------------------|--------------------| ---------------|-------|---------------------------|-------------------|
| ECM/Fibrosis | 80 | 12 | 68 | **85%** | UP in disease | Anti-TGFbeta, Pirfenidone |
| Metabolic | 60 | 3 | 57 | **95%** | DN in disease | AMPK activator, Metformin |
| DNA Damage | 50 | 8 | 42 | **84%** | DN in disease | PARP modulator |
| ... | ... | ... | ... | ... | ... | ... |

### Step 3: Non-Responder Specific Gaps

**Output Table C2: Non-Responder Modules NOT Covered by Target**

| Module | Enrichment in Non-Responders | Target Coverage % | Gap % | Combo to Improve Response |
|--------|------------------------------|-------------------|-------|---------------------------|
| ECM/Fibrosis | UP (NES > 0, more fibrosis) | 12% | 88% | Anti-fibrotic combination |
| Immune checkpoint | UP (NES > 0, more checkpoint) | 8% | 92% | Checkpoint inhibitor |
| Metabolic | DN (NES < 0, less metabolism) | 5% | 95% | Metabolic modulator |

**Interpretation**: Non-responder-specific gaps indicate **resistance mechanisms** - combos should address these to improve response rate.

### Step 4: Gap → Combo Translation

**Combo Opportunity Logic**:

| Gap Type | Combo Strategy | Example |
|----------|---------------|---------|
| Residual immune gaps | Add immune-targeting combo | JAK inhibitor + anti-IL12/23 |
| ECM/Fibrosis gaps | Add anti-fibrotic combo | Target X + anti-TGFbeta |
| Metabolic gaps | Add metabolic modulator | Target X + PPARgamma agonist |
| DNA repair gaps | Add DNA repair enhancer | Target X + ATM/ATR modulator |
| Non-responder resistance | Address resistance mechanism | Target X + checkpoint inhibitor |

### Step 5: Prioritize Combination Targets

**Output Table C3: Prioritized Combination Targets**

| Priority | Combo Target | Gap Module(s) Addressed | Coverage Increase | Clinical Stage | Key Drugs |
|----------|--------------|-------------------------|-------------------|----------------|-----------|
| 1 | JAK inhibitor | Immune (residual) | +15% → 72% total | Approved | Upadacitinib, Tofacitinib |
| 2 | Anti-TGFbeta | ECM/Fibrosis | +18% → 75% total | Phase 2 | Pirfenidone, Fresolimumab |
| 3 | IL-23 inhibitor | Th17/IL-17 | +8% → 65% total | Approved | Risankizumab, Guselkumab |

### Step 6: Context-Aware Gap Interpretation

**CRITICAL**: Interpret gaps based on comparison category:

| comparison_category | Gap Interpretation | Combo Opportunity Type |
|---------------------|-------------------|------------------------|
| **Disease vs. Normal** | Target pathways NOT enriched in disease | Fill core disease gaps |
| **Treatment vs. Control** | Target pathways NOT affected by treatment | Add complementary mechanism |
| **Responder vs. Non-Responder** | Target pathways enriched in non-responders = resistance | Overcome resistance |

**DO NOT**:
- List individual pathway names in gap analysis
- Provide detailed pathway mechanisms
- Deep dive on specific pathways

**VERIFY**: Coverage numbers MUST be consistent with Section B findings.
```

---

#### 3.2 Safety Analyst Task (Agent D)

```markdown
## Task: Safety Module Overlap Assessment

**Input Files**:
- gsea_{target}_{indication}_pathwaydb.csv (target-associated pathways)
- Approved Section A results (disease modules)
- Approved Section B results (target coverage)
- Planner feedback/context from Phase 2

**Instructions**:

### Step 0: Filter to Interpretable Comparison Categories (CRITICAL)

**BEFORE any analysis**, filter data to ONLY these 3 categories:
1. Disease vs. Normal
2. Treatment vs. Control
3. Responder vs. Non-Responder

**EXCLUDE**: Other Comparisons, CellType1 vs. CellType2, Tissue1 vs. Tissue2

```python
# Filter code
included_categories = [
    'Disease vs. Normal',
    'Treatment vs. Control',
    'Responder vs. Non-Responder'
]
df_filtered = df[df['comparison_category'].isin(included_categories)]
```

### Step 1: Understand Safety Risk Definition

**CRITICAL**: Safety Risk requires BOTH conditions:
1. Target pathway **significantly enriched** (UP or DN, FDR < 0.05) in disease/comparison
2. AND the pathway/module is **associated with known safety concerns** (Hematologic, Hepatic, Cardiac, etc.)

**Safety Risk = Enriched Target Pathway ∩ Safety Module**

### Step 2: Define Safety-Critical Modules

| Safety Module | Pathways Included | Associated Clinical Risk |
|---------------|-------------------|--------------------------|
| Hematologic | Cell cycle, Apoptosis in blood, Hematopoiesis | Cytopenia, Anemia, Thrombosis |
| Hepatic | Drug metabolism, CYP450, Detoxification | Hepatotoxicity, LFT elevation |
| Immune Suppression | Inflammatory response, Cytokines, Complement | Infection, Opportunistic pathogens |
| Cardiac | Ion channels, Cardiac muscle, Conduction | Arrhythmia, QT prolongation |
| Renal | Transport, Filtration, Organic anion | Nephrotoxicity, Cr elevation |
| Neurological | Neurotransmitter, Synaptic, Neuro-signaling | CNS effects, Neuropathy |
| GI Toxicity | Epithelial turnover, Mucus production | Diarrhea, Nausea |

### Step 3: Calculate Target-Safety Module Overlap

**Output Table D1: Target Pathways Enriched in Disease × Safety Modules**

| Safety Module | Target Pathways in Module | Sig UP (FDR<0.05) | Sig DN (FDR<0.05) | Enrichment Direction | Risk Level | Interpretation |
|---------------|--------------------------|-------------------|-------------------|---------------------|------------|----------------|
| Hematologic | 45 | 23 | 5 | Mostly UP | MODERATE | Target pathways activated in disease may affect blood counts |
| Immune Suppression | 89 | 67 | 12 | Mostly UP | HIGH | Target pathways activated in disease = therapeutic but infection risk |
| Hepatic | 12 | 2 | 1 | Minimal | LOW | Few target pathways enriched in hepatic module |
| Cardiac | 8 | 0 | 0 | NONE | LOW | No significant enrichment in cardiac module |
| Renal | 15 | 3 | 8 | Mostly DN | LOW-MODERATE | Target pathways suppressed in disease |

### Step 4: Safety Risk Assessment by Enrichment Direction

**Safety Risk Scoring by Direction**:

| Enrichment Direction | Safety Module | Interpretation | Risk Example |
|----------------------|---------------|----------------|--------------|
| **UP** (NES > 0) | Hematologic | Target pathways activated in disease + safety concern | May exacerbate blood count abnormalities |
| **UP** (NES > 0) | Apoptosis | Target pathways activated in disease + cell death | May cause tissue damage |
| **DN** (NES < 0) | Immune | Target pathways suppressed in disease + safety concern | May worsen immune suppression |
| **DN** (NES < 0) | Hepatic | Target pathways suppressed in disease + safety concern | May worsen liver function |

### Step 5: Risk Summary by Comparison Category

**Output Table D2: Risk Summary**

| Risk Category | Safety Module | Enrichment Direction | comparison_category | Severity | Frequency Estimate | Monitoring |
|---------------|---------------|---------------------|---------------------|----------|-------------------|------------|
| Immunosuppression | Immune | UP in disease | Disease vs. Normal | Therapeutic but risk | Common (10-15%) | Infection screen, TB/HBV |
| Cytopenia | Hematologic | UP in disease | Disease vs. Normal | Moderate | Uncommon (2-5%) | CBC q2wk × 12 wk |
| Response-associated risk | Immune | UP in responders | Responder vs. Non-Responder | Correlates with response | Variable | Monitor in responders |

### Step 6: Monitoring Protocol

**Output Table D3: Recommended Monitoring Protocol**

| Timepoint | Tests | Purpose | Action Threshold |
|-----------|-------|---------|------------------|
| Baseline | CBC, CMP, LFTs, TB, HBV/HCV | Risk assessment | Screen out high-risk |
| Week 2 | CBC | Early hematologic toxicity | ANC <1000: hold |
| Week 4 | CBC, LFTs | Confirm tolerability | ALT >3x ULN: hold |
| Week 8 | CBC, LFTs, Cr | Ongoing safety | Per standard criteria |
| Month 3+ | CBC, LFTs q3mo | Maintenance monitoring | Per standard criteria |

### Step 7: Context-Aware Safety Interpretation

**CRITICAL**: Interpret safety risks based on comparison category:

| comparison_category | Safety Interpretation |
|---------------------|----------------------|
| **Disease vs. Normal** | Baseline safety signal - target pathways enriched in disease + safety module |
| **Treatment vs. Control** | On-target toxicity - safety modules affected by treatment |
| **Responder vs. Non-Responder** | Response-associated safety - different profiles for responders vs. non-responders |

**VERIFY**: Safety assessment must be consistent with A/B findings (disease modules, enrichment patterns).

**DO NOT**:
- Deep dive on specific pathway toxicity mechanisms
- List individual pathways by name
- Provide exhaustive pathway descriptions
```

---

### Parallel Track: Target Scoring Specialist Task (Agent E)

**RUNS IN PARALLEL** with Agents A-D (independent track from Track 1)

```markdown
## Task: Target Gene Scoring and Summary

**Input File**: target_summary.csv

**Instructions**:

### Step 1: Load Target Summary Data

```python
df_target = pd.read_csv('target_summary.csv')
# Expected columns: Gene, study, comparison, comparison_category, log2fc, padj
```

### Step 2: Score Each Target Gene

Apply scoring based on comparison category and DEG direction:

| comparison_category | DEG Direction | Score | Interpretation |
|---------------------|---------------|-------|----------------|
| Disease vs. Normal | Positive (log2fc > 0, padj < 0.05) | +1 | Target upregulated in disease |
| Disease vs. Normal | Negative (log2fc < 0, padj < 0.05) | +1 | Target downregulated in disease |
| Responder vs. Non-Responder | Negative (log2fc < 0 = higher in non-responders, padj < 0.05) | +1 | Target predicts non-response |
| Responder vs. Non-Responder | Positive (log2fc > 0 = higher in responders, padj < 0.05) | +1 | Target predicts response |
| Treatment vs. Control | Any significant (padj < 0.05) | +0.5 | Target affected by treatment |

```python
def score_target(row):
    score = 0
    if row['padj'] < 0.05:
        if row['comparison_category'] == 'Disease vs. Normal':
            score += 1  # Positive for upregulated in disease
        elif row['comparison_category'] == 'Responder vs. Non-Responder':
            if row['log2fc'] < 0:  # Higher in non-responders
                score += 1  # Resistance marker
        elif row['comparison_category'] == 'Treatment vs. Control':
            score += 0.5
    return score

df_target['score'] = df_target.apply(score_target, axis=1)
```

### Step 3: Aggregate Scores by Target

```python
target_scores = df_target.groupby('Gene').agg({
    'score': 'sum',
    'comparison_category': lambda x: list(x.unique()),
    'study': 'count'
}).rename(columns={'study': 'n_comparisons'})

target_scores = target_scores.sort_values('score', ascending=False)
```

### Step 4: Generate Target Ranking Output

**Output Table E1: Target Ranking Summary**

| Rank | Target | Total Score | N Comparisons | Key Categories | Interpretation |
|------|--------|-------------|---------------|----------------|----------------|
| 1 | GREM1 | 12.5 | 45 | Disease vs. Normal, Responder | Strong disease + resistance signal |
| 2 | IL11 | 8.0 | 32 | Disease vs. Normal | Strong disease signal |
| 3 | TGFb1 | 6.5 | 28 | Disease vs. Normal, Treatment | Disease + treatment response |

### Step 5: Generate Markdown Summary

Save to: `reports/target_ranking_report.md`

```markdown
# Target Gene Ranking Report

## Summary

| Metric | Value |
|--------|-------|
| Total targets analyzed | N |
| Top scoring target | {target_name} (score: X) |
| Average score | Y |

## Top Ranked Targets

{Table E1}

## Interpretation

- **Disease vs. Normal contributions**: Targets with high scores in this category are strongly associated with disease pathology
- **Responder vs. Non-Responder contributions**: High scores indicate potential as biomarkers or resistance mechanisms
- **Treatment response contributions**: Scores indicate pharmacodynamic relevance

## Per-Target Analysis

### {Target 1}
- Total score: X
- Disease vs. Normal: Y comparisons significant
- Clinical implication: {interpretation}

### {Target 2}
...
```

**Output Files**:
- `reports/target_ranking_report.md` - Full target scoring report
- `reports/target_scores.csv` - Raw scores for downstream analysis
```

---

### Phase 4: Dependent Review (Planner Reviews C + D)

**Planner reviews Section C and Section D.**

**Review Checklist for C:**
- [ ] Gap = FDR ≥ 0.05 (not significantly enriched)
- [ ] Gaps categorized by module
- [ ] Non-responder gaps identified as resistance mechanisms
- [ ] Combo recommendations align with A/B findings
- [ ] Coverage numbers consistent with Section B

**Review Checklist for D:**
- [ ] Safety risk = enriched (FDR < 0.05) + safety module
- [ ] Enrichment direction (UP/DN) reported
- [ ] Safety assessment consistent with A/B findings
- [ ] Monitoring protocol provided

**Decision Point:**

```
┌─ C/D review reveals A/B issues → Go to Phase 5 (Iteration)
│
└─ C/D APPROVED and consistent with A/B → Go to Phase 6 (Final Assembly)
```

---

### Phase 5: Iteration Cycle (If Needed)

**Trigger:** C/D review reveals issues with A/B foundation

**Iteration Protocol:**

**Step 1: Update A/B First**
- Send critique to A and/or B agents with specific issues found during C/D review
- A/B agents revise their sections
- Planner reviews revised A/B
- Repeat until A/B APPROVED

**Step 2: Update C/D with Revised A/B**
- Send updated A/B results to C and D agents
- Include specific feedback from iteration review
- C/D agents update their sections based on revised A/B
- Planner reviews updated C/D

**Step 3: Full Team Review**
- Planner reviews ALL sections (A, B, C, D) together
- Check cross-section consistency:
  - Coverage numbers in B match gap calculations in C
  - Safety modules in D align with disease modules in A
  - Combo targets in C address gaps identified from B
- Verify no contradictions between sections

**Step 4: Repeat if Needed**
- If any section fails review, return to Step 1 or Step 2
- Maximum 3 iteration cycles before escalating to user

**Exit Condition:** All 4 sections APPROVED and mutually consistent

---

### Phase 6: Final Assembly (All Sections Approved)

**Prerequisites:** A, B, C, D all APPROVED and cross-section consistent

**Planner compiles final report:**

1. **Executive Summary**
   - Overall coverage metrics from Section B
   - Key gaps from Section C
   - Safety concerns from Section D

2. **Assemble Sections**
   - Section A: Disease Module Characterization
   - Section B: Target-Disease Coverage
   - Section C: Gap Analysis (Combination Opportunities)
   - Section D: Safety Analysis

3. **Add Cross-References**
   - Link coverage gaps (C) to module findings (A)
   - Link safety risks (D) to enrichment patterns (B)

4. **Appendix**
   - Key definitions
   - Data sources
   - Module classification criteria

5. **Save Final Report**
   - Output: `{target}_{indication}_team_report.md`
   - Location: `{folder}/reports/`

---

## Final Report Structure

```markdown
# {Target} Coverage Analysis for {Indication}

**Analysis Date**: {DATE}
**Input Folder**: {target}_{indication}_omicsoft_{timestamp}/
**Output Location**: {folder}/reports/

---

## Executive Summary

### Overall Coverage Metrics

| Metric | Value | Interpretation |
|--------|-------|----------------|
| Total target pathway coverage | X% | [Good/Moderate/Poor] target alignment |
| UP coverage (Disease vs. Normal) | Y% | Target pathways activated in disease |
| DN coverage (Disease vs. Normal) | Z% | Target pathways suppressed in disease |
| Responder vs Non-Responder coverage | W% | [Predicts/Does not predict] response |
| Key gaps | [List modules] | Combination opportunities |
| Key safety concerns | [List modules] | Monitoring required |

### Key Findings

1. **Disease Coverage**: {TARGET} pathways show X% total coverage in disease (Y% UP, Z% DN)
2. **Best Coverage**: [Module] at W% - target pathways significantly enriched
3. **Largest Gap**: [Module] at V% - target pathways NOT significantly enriched = combo opportunity
4. **Safety Signal**: [Risk] predicted from [safety module] × [enrichment direction] overlap

---

## Section A: Disease Module Characterization

### A1) Global Module Enrichment by Comparison × Tissue × Disease

| Module | comparison_category | Tissue | Disease | N_UP (FDR<0.05) | N_DN (FDR<0.05) | Avg_NES_UP | Avg_NES_DN |
|--------|---------------------|--------|---------|-----------------|-----------------|------------|------------|
| [Data from Disease Analyst] |

**Interpretation**: N_UP = pathways significantly activated; N_DN = pathways significantly suppressed

### A2) Global UP-Enriched Clusters/Networks

| Module | Total UP Pathways | Key Biological Themes | Dominant Comparison | Dominant Tissue |
|--------|-------------------|----------------------|---------------------|-----------------|
| [Data from Disease Analyst] |

### A3) Global DN-Enriched Clusters/Networks

| Module | Total DN Pathways | Key Biological Themes | Dominant Comparison | Dominant Tissue |
|--------|-------------------|----------------------|---------------------|-----------------|
| [Data from Disease Analyst] |

---

## Section B: Target-Disease Coverage

### B1) Module Coverage Matrix

| Module | Total Target Pathways | Sig UP (FDR<0.05) | Sig DN (FDR<0.05) | UP Coverage % | DN Coverage % | Total Coverage % |
|--------|----------------------|-------------------|-------------------|---------------|---------------|------------------|
| [Data from Target Specialist] |

**Coverage Formula**: Coverage % = (Target pathways with FDR < 0.05 in disease) / (ALL target pathways) × 100

### B2) Coverage by Comparison Category

| comparison_category | Total Target Pathways | Sig UP | Sig DN | UP Coverage % | DN Coverage % | Total Coverage % | Clinical Implication |
|---------------------|----------------------|--------|--------|---------------|---------------|------------------|----------------------|
| [Data from Target Specialist] |

### B3) Tissue-Specific Coverage

| Tissue | Total Target Pathways | Sig UP | Sig DN | UP Coverage % | DN Coverage % | Total Coverage % |
|--------|----------------------|--------|--------|---------------|---------------|------------------|
| [Data from Target Specialist] |

---

## Section C: Gap Analysis (Combination Opportunities)

### C1) Target Pathways NOT Significantly Enriched in Disease (Gaps)

| Gap Module | Total Target Pathways | Covered (FDR<0.05) | Gap (FDR≥0.05) | Gap % | Dominant Disease Direction | Combo Opportunity |
|------------|----------------------|--------------------| ---------------|-------|---------------------------|-------------------|
| [Data from Combo Strategist] |

**Gap Definition**: Gap = Target pathways with FDR ≥ 0.05 (not significantly enriched in disease)

### C2) Non-Responder Specific Gaps (Resistance Mechanisms)

| Module | Enrichment in Non-Responders | Target Coverage % | Gap % | Combo to Improve Response |
|--------|------------------------------|-------------------|-------|---------------------------|
| [Data from Combo Strategist] |

### C3) Prioritized Combination Targets

| Priority | Combo Target | Gap Module(s) | Coverage Increase | Clinical Stage |
|----------|--------------|---------------|-------------------|----------------|
| [Data from Combo Strategist] |

---

## Section D: Safety Analysis

### D1) Target-Safety Module Overlap

| Safety Module | Target Pathways in Module | Sig UP | Sig DN | Enrichment Direction | Risk Level | Interpretation |
|---------------|--------------------------|--------|--------|---------------------|------------|----------------|
| [Data from Safety Analyst] |

**Safety Risk Definition**: Target pathway significantly enriched (FDR < 0.05) in disease AND overlaps with safety-critical module

### D2) Risk Summary

| Risk Category | Safety Module | Enrichment Direction | Severity | Frequency | Monitoring |
|---------------|---------------|---------------------|----------|-----------|------------|
| [Data from Safety Analyst] |

### D3) Monitoring Protocol

| Timepoint | Tests | Purpose | Action Threshold |
|-----------|-------|---------|------------------|
| [Data from Safety Analyst] |

---

## Appendix: Data Sources

### Input Files Used
- `{indication}_pathwaydb_gsea.csv` (universal pathway GSEA from combined_pathways.csv)
- `{target}_pathwaydb.csv` (target-associated pathways from /pathwaydb-query skill)
- `target_summary.csv` (target gene DEG summary across comparisons)

### Key Definitions

| Term | Definition |
|------|------------|
| Target-associated pathway | Pathway containing target gene in its gene set |
| UP enriched | NES > 0 AND FDR < 0.05 |
| DN enriched | NES < 0 AND FDR < 0.05 |
| Coverage % | (Target pathways with FDR < 0.05) / (ALL target pathways) × 100 |
| Gap | Target pathways with FDR ≥ 0.05 (not significantly enriched) |
| Safety risk | Enriched target pathway (FDR < 0.05) ∩ Safety module |

### Module Definitions
[Standardized module classification used in this analysis]
```

---

## Quality Control Criteria

### Report QC Checklist

- [ ] **FILTERED to 3 categories ONLY**: Disease vs. Normal, Treatment vs. Control, Responder vs. Non-Responder
- [ ] **EXCLUDED categories**: Other Comparisons, CellType1 vs. CellType2, Tissue1 vs. Tissue2
- [ ] **Key Definitions section** includes: target-associated pathway, UP/DN, coverage formula, gap, safety risk
- [ ] **A1** shows full breakdown: Module × Comparison × Tissue × Disease with N_UP, N_DN (FDR < 0.05)
- [ ] **A2** covers global UP-enriched clusters/networks
- [ ] **A3** covers global DN-enriched clusters/networks
- [ ] **Coverage formula** = (Target pathways FDR < 0.05 in disease) / (ALL target pathways)
- [ ] **UP and DN coverage** always reported separately
- [ ] **Gap** = target pathways NOT significantly enriched (FDR ≥ 0.05)
- [ ] **Safety risk** = target pathway enriched in disease (FDR < 0.05) + safety module
- [ ] **NO "target UP/DN" language** - uses "enriched UP/DN in disease"
- [ ] **Context-aware interpretation** applied based on comparison category
- [ ] NO individual pathway deep dives (no top 15/10/5 pathway details)
- [ ] File paths: input from target_indication_omicsoft_timestamp/, output to reports/

### Cross-Section Consistency Checks

- [ ] Coverage numbers in B match gap calculations in C
- [ ] Safety modules in D align with disease modules in A
- [ ] Combo targets in C address gaps identified from B
- [ ] No contradictions between sections

### Validation Commands

```bash
# Verify output files exist
ls {folder}/reports/*_analysis.md {folder}/reports/*_team_report.md

# Check for forbidden patterns (individual pathway lists)
grep -c "Top 15\|Top 10\|Top 5\|Pathway 1:\|Pathway 2:" {folder}/reports/*.md
# Should return 0

# Check for required patterns (coverage tables)
grep -c "Coverage %\|Gap %\|N_UP\|N_DN" {folder}/reports/*.md
# Should return > 10

# Check for correct language (no "target UP/DN")
grep -c "target UP\|target DN" {folder}/reports/*.md
# Should return 0 (use "enriched UP/DN" instead)
```

---

## Usage

### Automatic Team Triggering

This team workflow is automatically triggered when `--pathwaydb-gsea` flag is used in `/omicsoft-analysis`.

### Workflow Summary

1. **GSEA Analysis Completes** → Output files generated:
   - `gsea_{target}_{indication}_pathwaydb.csv`
   - `gsea_all_results.csv`

2. **Report Generation Auto-Starts**:
   - Load prompt from: `~/.claude/skills/omicsoft-analysis/assets/gsea_report_prompt.md`
   - Create agent team: `gsea-{target}-{indication}`
   - Spawn agents following 6-phase workflow

3. **Agent Team Structure**:
   | Agent | Section | Input File | Model |
   |-------|---------|------------|-------|
   | disease-analyst | A: Disease Modules | gsea_all_results.csv | Opus |
   | target-specialist | B: Target Coverage | gsea_{target}_{indication}_pathwaydb.csv | Opus |
   | combo-strategist | C: Gap Analysis | A + B results | Opus |
   | safety-analyst | D: Safety | gsea_{target}_{indication}_pathwaydb.csv | Opus |

4. **Output**: `{target}_{indication}_team_report.md` in `reports/` folder

### Manual Execution

```bash
# Create team and run analysis
claude --team gsea-{target}-{indication} \
  --context "gsea_{target}_{indication}_pathwaydb.csv" \
  --prompt "$(cat ~/.claude/skills/omicsoft-analysis/assets/gsea_report_prompt.md)"
```

### Team Commands

```
# Monitor team progress
/team status gsea-{target}-{indication}

# View agent outputs
/team output disease-analyst
/team output target-specialist

# Compile final report
/team compile gsea-{target}-{indication} --output {folder}/reports/{target}_{indication}_team_report.md
```
