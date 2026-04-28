---
name: kg-association
description: Analyze biomedical entity associations using multiple knowledge graph methods (BioBridge, ULTRA/UltraQuery, PrimeKG). Use this skill when users request multi-entity association analysis, gene-disease predictions, combo/signature analysis, or cross-method comparison. Generates structured markdown reports with biological interpretations for individual entities, combos, and comparative analysis across all three methods.
---

# Knowledge Graph Association Analysis

## Overview

This skill performs comprehensive biomedical entity association analysis using knowledge graph MCP tools:
- **BioBridge**: Neural KG link prediction with mean embeddings for combos
- **ULTRA**: Foundation model for predictions (ultra-inference for individuals, ultraquery-inference for combos)
- **PrimeKG**: Shortest path analysis using graph traversal (via primekg skill script)

**Output**: Structured markdown reports in `kgpred_<context>_<date>/` directory

## Output Location Policy

**CRITICAL: All results, reports, and data files MUST be saved in the user's working directory.**

### Default Behavior
- All output goes to: `<working_dir>/kgpred_<context>_<date>/`
- Working directory = where the user invoked the skill (current directory)
- **NEVER** save outputs to skill folders, MCP folders, or other system directories

### User Override
- If user specifies a custom output directory, use that instead
- Example: User says "save to /path/to/my/results" → use `/path/to/my/results/kgpred_<context>_<date>/`

### What Gets Saved to Working Directory
1. **Markdown Reports**: All Part1, Part2, Part3, Cross_MCP_Comparison.md files
2. **Data Tables**: Any CSV/TSV exports of results
3. **Parquet Files**: Copy from MCP output folders to working directory
4. **Summary Files**: Any additional analysis outputs

### MCP Output File Handling

**CRITICAL: MOVE all MCP-generated data files to the working directory IMMEDIATELY after each query.**

MCP tools write outputs to their internal directories. These files MUST be moved to the user's working directory:
- **IMMEDIATELY** after each query completes (not batched at end)
- For user accessibility (MCP folders may not be accessible)
- For result reproducibility and archiving
- For downstream analysis by users

#### Direct Output with output_dir Parameter (Recommended)

**ULTRA MCP tools now support `output_dir` parameter for direct file output:**

1. **Run MCP query with output_dir** → Files saved directly to your working directory
2. **Read and process the file** → Extract scores (no move needed!)
3. **Proceed to next query**

```python
from pathlib import Path

# Define destination in working directory
working_dir = Path("./kgpred_<context>_<date>")
data_dir = working_dir / "data" / "ultra"
data_dir.mkdir(parents=True, exist_ok=True)

# MCP writes directly to your output_dir - no move needed!
result = mcp__ultra-inference__predict_tail_entities(
    head_entity="TYK2",
    relation="associated with",
    output_dir=str(data_dir)  # Files saved directly here
)

# Read directly from the output path
df = pd.read_parquet(result['output_file'])
```

#### Legacy: Move-Then-Process Workflow (Backward Compatible)

If `output_dir` is not specified, files are saved to MCP default directory and must be moved:

```python
import shutil
from pathlib import Path

# MCP returns output_file path (in MCP's internal directory)
mcp_output_file = result['output_file']  # or result['output_file_filtered']

# Define destination in working directory
working_dir = Path("./kgpred_<context>_<date>")
dest_subdir = working_dir / "data" / "<mcp_name>"  # e.g., data/ultra or data/biobridge
dest_subdir.mkdir(parents=True, exist_ok=True)

# Move file to working directory
dest_file = dest_subdir / Path(mcp_output_file).name
shutil.move(mcp_output_file, dest_file)
```

#### Directory Structure with Data Files

```
kgpred_<context>_<date>/
├── biobridge/
│   ├── Part1_Individual_Analysis.md
│   ├── Part2_Combo_Analysis.md
│   └── Part3_Comparative_Analysis.md
├── ultra/
│   ├── Part1_Individual_Analysis.md
│   ├── Part2_Combo_Analysis.md
│   └── Part3_Comparative_Analysis.md
├── data/
│   ├── ultra/
│   │   ├── TYK2_JAK1_predictions.parquet
│   │   ├── CDKN2D_ITGA4_ITGB7_predictions.parquet
│   │   └── ...
│   └── biobridge/
│       └── (if any parquet outputs)
└── Cross_MCP_Comparison.md
```

#### File Naming Convention for Moved Files

When moving, use descriptive names:
- Individual queries: `<gene>_predictions.parquet`
- Combo queries: `<gene1>_<gene2>_predictions.parquet` or `<gene1>_<gene2>_<gene3>_predictions.parquet`
- Filtered outputs: `<combo>_predictions_filtered.parquet`

#### Bash Alternative

```bash
# After ULTRA query, move output files
mv /path/to/mcp/output/predictions.parquet ./kgpred_IBD_2025-01-03/data/ultra/TYK2_JAK1_predictions.parquet
```

### MCP Folder Cleanup
- MCP tools write raw outputs to their own directories
- Files are MOVED (not copied) to working directory for user access
- MCP output folders will be emptied after move
- User accesses all results from working directory only

---

## Scoring Methodology

### Primary Score: Percentile Rank (pct_rank)

**All scores in reports use percentile rank (0-1 scale, where 1 = top ranked)**

This ensures consistent interpretation across both MCPs:
- **pct_rank = 0.99** means top 1% of all predictions
- **pct_rank = 0.95** means top 5% of all predictions
- Higher pct_rank = stronger association

### BioBridge Scoring
- **Primary score**: `pct_rank` (0-1, higher = stronger)
- **Also report**: `cos_sim` (cosine similarity) for reference
- pct_rank indicates where this prediction ranks among all possible tail entities

### ULTRA Scoring (ultra-inference AND ultraquery-inference)

**Primary score**: `pct_rank` recalculated within entity1_type only (e.g., disease only)
**Also report**: `rank` recalculated within entity1_type only

#### Why Recalculation is Required

**CRITICAL**: Both `ultra-inference` and `ultraquery-inference` return native percentile ranks that are NOT suitable for direct use:

| MCP Tool | Native Field | What It Measures | Problem |
|----------|--------------|------------------|---------|
| ultra-inference | `percentile_rank` | Rank across ALL entity types (~129K) | Includes genes, drugs, phenotypes, etc. |
| ultraquery-inference | `percentile_rank` | Rank across ALL entity types (~129K) | Same issue |
| ultraquery-inference | `filtered_percentile_rank` | Rank within schema-matched types | Combines ALL schema-matched entity types |

**Example**: For a 2i query (TYK2 ∩ JAK1 → diseases):
- `filtered_predictions` = 32,391 (17,080 diseases + 15,311 phenotypes)
- A phenotype could be rank 1, pushing diseases down
- `filtered_percentile_rank` does NOT give disease-only ranking

**Solution**: Always filter to target entity type and recalculate.

#### Calculation Method (pandas)

```python
import pandas as pd

# Read ULTRA output parquet file (works for both ultra-inference and ultraquery-inference)
df = pd.read_parquet(output_file)

# STEP 1: Filter to entity1_type only (e.g., disease only, NOT disease+phenotype)
# For ultra-inference: column is 't_pred_type'
# For ultraquery-inference: column is 'entity_type'
df_filtered = df[df['entity_type'] == entity1_type].copy()  # e.g., 'disease'

# STEP 2: Recalculate rank within entity type only (rank 1 = highest score = best)
df_filtered['rank'] = df_filtered['score'].rank(ascending=False, method='first').astype(int)

# STEP 3: Recalculate pct_rank within entity type only
# Formula: 1.0 - (rank / total) where total = count of entity1_type only
# Result: rank 1 → pct_rank ≈ 1.0 (best), rank N → pct_rank ≈ 0.0 (worst)
df_filtered['pct_rank'] = 1.0 - (df_filtered['rank'] / len(df_filtered))
```

#### Example: Disease-Only Recalculation

```python
# Native filtered_percentile_rank (disease + phenotype combined):
# periapical periodontitis: filtered_rank=2, filtered_pct_rank=0.9999

# After disease-only recalculation (17,080 diseases):
# periapical periodontitis: disease_rank=1, disease_pct_rank=0.9999
# (Now correctly ranked #1 among diseases only)
```

#### Result Interpretation
- `pct_rank ≈ 1.0`: TOP ranked among entity1_type = STRONG association
- `pct_rank ≈ 0.0`: BOTTOM ranked among entity1_type = WEAK association
- `rank`: position among entity1_type only (1 = best, N = worst)
- Higher pct_rank = better rank = stronger association

## When to Use This Skill

Invoke this skill when users request:

- Multi-entity association analysis (genes → diseases, drugs → phenotypes)
- Gene signature/combo analysis with synergy/dilution assessment
- Comparative predictions across multiple MCPs
- Structured reports with biological interpretations

## Required Inputs

### 1. entity1_list (Target Entities)
List of target entities to predict associations with (e.g., diseases):
```
- Crohn's disease
- Ulcerative colitis
- Inflammatory bowel disease
```

### 2. entity2_list (Source Entities)
Source entities in one of two formats:

**Format A: Flat list (individual entities only)**
```
- TYK2
- JAK1
- GREM1
```
Or: `[TYK2, JAK1, GREM1]` or `[a+b, c+d]` (single items)

→ **Only Part 1 (Individual Analysis) is generated**

**Format B: Nested list (combos)**
```
- [TYK2, JAK1]
- [TNFRSF25, GREM1]
- [TNFRSF25, PCOLCE]
- [CDKN2D, ITGA4, ITGB7]
- [CDKN2D, PCOLCE]
```

→ **Full analysis: Part 1, Part 2, Part 3**
→ Individual entities derived from union of all combo elements

### Optional Parameters
```yaml
entity1_type: disease          # Default: disease
entity2_type: "gene/protein"   # Default: gene/protein
relation_hint: "associated with"
context: "derived from entity1 or user-provided"
```

## Input Detection Logic

```
IF entity2_list contains nested lists (e.g., [[a,b], [c,d]]):
    mode = "COMBO"
    entity2_individuals = union of all combo elements
    entity2_combos = the nested lists
    → Generate Part 1, Part 2, Part 3
ELSE:
    mode = "INDIVIDUAL_ONLY"
    entity2_individuals = entity2_list
    entity2_combos = None
    → Generate Part 1 only
```

## Workflow

### Phase 0: Setup Output Directory

**Location**: User's current working directory (or user-specified directory)

Create output directory using naming convention:
```
<working_dir>/kgpred_<context>_<YYYY-MM-DD>/
```

Example: If working dir is `/home/user/projects/`, output goes to:
`/home/user/projects/kgpred_IBD_2025-01-03/`

```bash
# Use current working directory - DO NOT use skill or MCP folders
mkdir -p ./kgpred_<context>_$(date +%Y-%m-%d)/{biobridge,ultra}
```

**Important**:
- Use relative path `./` or absolute path to working directory
- Never hardcode paths to skill folders or MCP service directories

---

## Part 1: Individual Entity Analysis (ALWAYS GENERATED)

**For each MCP**, query each individual entity2 → entity1 pair:

### BioBridge Individual Queries

**CRITICAL: BioBridge Override Requirements**

When calling `mcp__biobridge__predict_associations`, you MUST ALWAYS provide explicit entity overrides:

| Parameter | Required | Description |
|-----------|----------|-------------|
| `override_head_name` | **MANDATORY** | Exact head entity name (e.g., "TYK2", "JAK1") |
| `override_head_type` | **MANDATORY** | Head entity type (e.g., "gene/protein") |
| `override_tail_name` | **MANDATORY** | Exact tail entity name (e.g., "hidradenitis suppurativa") |
| `override_tail_type` | **MANDATORY** | Tail entity type (e.g., "disease") |
| `relation_hint` | Recommended | Relation family (e.g., "associated with") |
| `topk` | Recommended | Set to 1 for pair validation queries |

**Why explicit overrides are required:**
- Ensures deterministic, reproducible results
- Prevents LLM inference errors on entity names
- Guarantees exact entity matching in the knowledge graph
- Required for accurate pct_rank extraction

```python
for entity2 in entity2_individuals:
    for entity1 in entity1_list:
        # MANDATORY: Always use explicit override parameters
        result = mcp__biobridge__predict_associations(
            context=f"Find association between {entity2} and {entity1}",
            override_head_name=entity2,      # REQUIRED: exact entity name
            override_head_type=entity2_type, # REQUIRED: entity type
            override_tail_name=entity1,      # REQUIRED: exact entity name
            override_tail_type=entity1_type, # REQUIRED: entity type
            relation_hint=relation_hint,
            topk=1
        )
        # Extract BOTH: cos_sim and pct_rank
        # Use pct_rank as primary SCORE
        score = result['pct_rank']
        cos_sim = result['cos_sim']  # Report for reference
```

### ULTRA Individual Queries (ultra-inference)

**Column names for ultra-inference parquet**: `t_pred_type`, `t_pred_name`, `t_pred_score`

```python
import pandas as pd
import shutil
from pathlib import Path

working_dir = Path("./kgpred_<context>_<date>")
data_dir = working_dir / "data" / "ultra"
data_dir.mkdir(parents=True, exist_ok=True)

for entity2 in entity2_individuals:
    result = mcp__ultra-inference__predict_tail_entities(
        head_entity=entity2,
        relation="associated with",
        top_k=None,  # Get all predictions for proper pct_rank calculation
        output_dir=str(data_dir)  # NEW: Direct output to working directory
    )

    # File is now saved directly to working directory - no move needed!
    output_file = result['output_file']

    # Read output parquet and filter to entity1_type ONLY (e.g., disease only)
    df = pd.read_parquet(output_file)
    df_filtered = df[df['t_pred_type'] == entity1_type].copy()  # e.g., 'disease'

    # CRITICAL: Recalculate rank and pct_rank within entity1_type only
    # Native percentile_rank includes ALL entity types - DO NOT USE directly
    df_filtered['rank'] = df_filtered['t_pred_score'].rank(ascending=False, method='first').astype(int)
    df_filtered['pct_rank'] = 1.0 - (df_filtered['rank'] / len(df_filtered))

    # For each entity1, extract pct_rank and rank
    for entity1 in entity1_list:
        match = df_filtered[df_filtered['t_pred_name'].str.contains(entity1, case=False)]
        if not match.empty:
            pct_rank = match['pct_rank'].iloc[0]  # Score (within entity1_type ONLY)
            rank = match['rank'].iloc[0]  # Rank among entity1_type only
```

**Report Format** (Part 1 tables):
```markdown
| Gene | Score (pct_rank) | Rank | Cos_Sim | Assessment |
|------|------------------|------|---------|------------|
| ITGA4 | 0.9959 | - | 0.776 | Very Strong |
```

**Generate**: `Part1_Individual_Analysis.md` for each MCP

---

## Part 2: Combo/Signature Analysis (ONLY IF COMBOS PROVIDED)

**Skip this phase if mode = "INDIVIDUAL_ONLY"**

### BioBridge Combo Queries (Mean Embeddings)

**CRITICAL: BioBridge Combo Override Requirements**

For combo/signature queries, use `override_head_names` (plural) with a list of entities:

| Parameter | Required | Description |
|-----------|----------|-------------|
| `override_head_names` | **MANDATORY** | List of head entity names (e.g., ["TYK2", "JAK1"]) |
| `override_head_type` | **MANDATORY** | Head entity type (e.g., "gene/protein") |
| `override_tail_name` | **MANDATORY** | Exact tail entity name (e.g., "hidradenitis suppurativa") |
| `override_tail_type` | **MANDATORY** | Tail entity type (e.g., "disease") |

```python
for combo in entity2_combos:
    for entity1 in entity1_list:
        # MANDATORY: Always use explicit override parameters for combos
        result = mcp__biobridge__predict_associations(
            context=f"Find diseases associated with {'+'.join(combo)} signature",
            override_head_names=combo,       # REQUIRED: list of entity names
            override_head_type=entity2_type, # REQUIRED: entity type
            override_tail_name=entity1,      # REQUIRED: exact entity name
            override_tail_type=entity1_type, # REQUIRED: entity type
            relation_hint=relation_hint,
            topk=1
        )
        # Extract BOTH: cos_sim and pct_rank
        # Use pct_rank as primary SCORE
        score = result['pct_rank']
        cos_sim = result['cos_sim']  # Report for reference
```

### ULTRA Combo Queries (ultraquery-inference with Intersection)

**Column names for ultraquery-inference parquet**: `entity_type`, `entity_name`, `score`

**Note on top_k**: The `top_k` parameter only limits the JSON API response, NOT the parquet output files.
Parquet files always contain ALL predictions regardless of `top_k` value.

```python
import pandas as pd
import shutil
from pathlib import Path

working_dir = Path("./kgpred_<context>_<date>")
data_dir = working_dir / "data" / "ultra"
data_dir.mkdir(parents=True, exist_ok=True)

for combo in entity2_combos:
    # Build intersection query based on combo size
    if len(combo) == 2:
        query = [
            [combo[0], ["associated with"]],
            [combo[1], ["associated with"]]
        ]
    elif len(combo) == 3:
        query = [
            [combo[0], ["associated with"]],
            [combo[1], ["associated with"]],
            [combo[2], ["associated with"]]
        ]

    result = mcp__ultraquery-inference__answer_complex_query(
        query_structure=query,
        top_k=5,  # Only affects API response; parquet files always contain ALL predictions
        output_dir=str(data_dir)  # NEW: Direct output to working directory
    )

    # Files are now saved directly to working directory - no move needed!
    output_file_filtered = result['output_file_filtered']
    output_file_all = result['output_file_all']

    # Read output parquet and filter to entity1_type ONLY (e.g., disease only)
    # CRITICAL: filtered_percentile_rank combines ALL schema-matched types - DO NOT USE directly
    df = pd.read_parquet(output_file_filtered)
    df_filtered = df[df['entity_type'] == entity1_type].copy()  # e.g., 'disease'

    # CRITICAL: Recalculate rank and pct_rank within entity1_type only
    # Native filtered_percentile_rank includes ALL schema-matched entity types
    df_filtered['rank'] = df_filtered['score'].rank(ascending=False, method='first').astype(int)
    df_filtered['pct_rank'] = 1.0 - (df_filtered['rank'] / len(df_filtered))

    # For each entity1, extract pct_rank and rank
    for entity1 in entity1_list:
        match = df_filtered[df_filtered['entity_name'].str.contains(entity1, case=False)]
        if not match.empty:
            pct_rank = match['pct_rank'].iloc[0]  # Score (within entity1_type ONLY)
            rank = match['rank'].iloc[0]  # Rank among entity1_type only
```

**Calculate Synergy/Dilution** (using pct_rank scores):
```
individual_mean = mean(pct_rank scores of combo components from Part 1)
combo_score = pct_rank from signature query
delta = combo_score - individual_mean

if delta > 0.02: "SYNERGY"
elif delta < -0.02: "DILUTION"
else: "NEAR-ADDITIVE"
```

**Report Format** (Part 2 tables):
```markdown
| Gene Combo | Score (pct_rank) | Cos_Sim | Individual Mean | Delta | Classification |
|------------|------------------|---------|-----------------|-------|----------------|
| TYK2+JAK1 | 0.9951 | 0.565 | 0.9947 | +0.0004 | NEAR-ADDITIVE |
```

**Generate**: `Part2_Combo_Analysis.md` for each MCP

---

## Part 3: Comparative Analysis (ONLY IF COMBOS PROVIDED)

**Skip this phase if mode = "INDIVIDUAL_ONLY"**

For each combo, generate deep comparative analysis:

1. **Signature vs Individual Performance**: Table comparing combo scores (pct_rank) to each component
2. **Synergy/Dilution Classification**: Categorize each combo based on pct_rank delta
3. **Biological Mechanisms**: Explain WHY synergy/dilution occurs
4. **Clinical Implications**: Biomarker utility, therapeutic strategies
5. **Design Principles**: What makes successful signatures

**Table Format**:
```markdown
| Combo | Combo Score | Component 1 | Component 2 | Component 3 | Mean | Δ | Status |
|-------|-------------|-------------|-------------|-------------|------|---|--------|
| TYK2+JAK1 | 0.9951 | TYK2: 0.9940 | JAK1: 0.9953 | - | 0.9947 | +0.0004 | ≈ |
```

**Generate**: `Part3_Comparative_Analysis.md` for each MCP

---

## PrimeKG Shortest Path Analysis (ALWAYS RUN)

**Run alongside BioBridge and ULTRA for comprehensive cross-method comparison.**

### PrimeKG Scoring Methodology

- **Method**: Shortest path using BFS (unweighted graph)
- **Score formula**: `score = 0.9^(path_length - 1)`
  - 1 hop (direct): 1.0
  - 2 hops: 0.9
  - 3 hops: 0.81
  - 4 hops: 0.729
  - No path: 0.0
- **Combo method**: Average of component scores (no synergy/dilution detection)
- **Path extraction**: Full intermediate nodes are extracted for each path

### PrimeKG Script Location

**Script**: `/home/sagemaker-user/.claude/skills/primekg/scripts/kg_association_shortest_path.py`

This script is specifically designed to integrate seamlessly with the kg-association workflow.

### Running PrimeKG Analysis

#### Option 1: Command Line (Recommended)

```bash
cd /home/sagemaker-user/.claude/skills/primekg && pixi run python scripts/kg_association_shortest_path.py \
    --genes "TYK2,JAK1,IL1RAP,IL17A,IL17F,IL11" \
    --diseases "hidradenitis suppurativa" \
    --combos "TYK2+JAK1,IL1RAP+IL17A+IL17F,IL1RAP+IL11" \
    --output-dir "./kgpred_<context>_<date>" \
    --data-dir "/home/sagemaker-user/github/ai-sci-mcp-services/ultra-inference/data"
```

**Arguments**:
- `--genes`: Comma-separated list of individual genes
- `--diseases`: Comma-separated list of target diseases
- `--combos`: Plus-separated gene combos, comma-separated (e.g., "A+B,C+D+E")
- `--output-dir`: Working directory for kg-association output
- `--data-dir`: Path to PrimeKG data directory containing `primekg.csv`

#### Option 2: Python API

```python
import sys
sys.path.append('/home/sagemaker-user/.claude/skills/primekg/scripts')
from kg_association_shortest_path import PrimeKGAssociationAnalysis

# Initialize with data directory
analysis = PrimeKGAssociationAnalysis(
    data_dir="/home/sagemaker-user/github/ai-sci-mcp-services/ultra-inference/data"
)

# Run full analysis
results = analysis.run_analysis(
    genes=["TYK2", "JAK1", "IL1RAP", "IL17A", "IL17F", "IL11"],
    combos=[["TYK2", "JAK1"], ["IL1RAP", "IL17A", "IL17F"], ["IL1RAP", "IL11"]],
    diseases=["hidradenitis suppurativa"],
    output_dir="./kgpred_<context>_<date>"
)
```

### Script Output

The script generates:

1. **JSON Data File**: `data/primekg/primekg_shortest_paths.json`
   - Complete results with all scores, paths, and metadata
   - Includes individual and combo analyses

2. **Markdown Reports**:
   - `primekg/Part1_Individual_Analysis.md`
   - `primekg/Part2_Combo_Analysis.md`
   - `primekg/Part3_Comparative_Analysis.md`

### JSON Output Structure

```json
{
  "genes": ["TYK2", "JAK1", ...],
  "diseases": ["hidradenitis suppurativa"],
  "combos": ["TYK2+JAK1", ...],
  "individual": {
    "TYK2": {
      "hidradenitis suppurativa": {
        "score": 0.81,
        "path_length": 3,
        "path_found": true,
        "path": ["TYK2", "brain", "NCSTN", "hidradenitis suppurativa"],
        "source_id": "NCBI:7297",
        "target_id": "MONDO:6559"
      }
    }
  },
  "combo": {
    "TYK2+JAK1": {
      "hidradenitis suppurativa": {
        "combo_score": 0.81,
        "individual_mean": 0.81,
        "delta": 0.0,
        "classification": "NEAR-ADDITIVE",
        "component_scores": {
          "TYK2": {"score": 0.81, "path_length": 3, "path": [...]},
          "JAK1": {"score": 0.81, "path_length": 3, "path": [...]}
        }
      }
    }
  },
  "summary": {
    "total_pairs": 6,
    "paths_found": 6,
    "avg_score": 0.81,
    "avg_path_length": 3.0
  }
}
```

### Report Format (Part 1 Individual Table)

```markdown
| Gene | Disease | Path Length | Score | Path | Assessment |
|------|---------|-------------|-------|------|------------|
| TYK2 | HS | 3 | 0.81 | TYK2 -> brain -> NCSTN -> HS | Moderate-Strong |
| IL17A | HS | 3 | 0.81 | IL17A -> extracellular region -> NLRP3 -> HS | Moderate-Strong |
```

### Report Format (Part 2 Combo Table)

```markdown
| Combo | Disease | Combo Score | Component Scores | Max | Min | Classification |
|-------|---------|-------------|------------------|-----|-----|----------------|
| TYK2+JAK1 | HS | 0.81 | TYK2: 0.81, JAK1: 0.81 | 0.81 | 0.81 | NEAR-ADDITIVE |
```

### Part 3 - PrimeKG Comparative Analysis

- Compare combo average to individual components
- Identify which component contributes most to combo score
- Highlight path structure insights (direct vs indirect associations)
- Analyze intermediate nodes (e.g., NCSTN, NLRP3) for mechanistic insights
- Note: No synergy/dilution possible with averaging method (delta always 0)

**Generate**: `Part1_Individual_Analysis.md`, `Part2_Combo_Analysis.md`, `Part3_Comparative_Analysis.md` in `primekg/` folder

---

## Cross-Method Comparison Report (ALWAYS GENERATED)

**This is the FINAL COMPREHENSIVE REPORT that summarizes and compares ALL THREE methods: BioBridge, ULTRA, and PrimeKG.**

### Structure of Cross_Method_Comparison.md

---

### SECTION A: Method Summaries

#### A1. BioBridge Summary
**Part 1 - Individual Analysis Summary**
- Table of all gene-disease scores (pct_rank, cos_sim)
- Top performers ranking
- Key findings

**Part 2 - Combo Analysis Summary**
- Table of all combo scores with synergy/dilution classification
- Best performing combos

**Part 3 - Comparative Analysis Summary**
- Design principles identified
- Mechanistic insights

---

#### A2. ULTRA Summary
**Part 1 - Individual Analysis Summary**
- Table of all gene-disease scores (pct_rank, rank)
- Top performers ranking
- Key findings

**Part 2 - Combo Analysis Summary**
- Table of all combo intersection scores
- Classification summary

**Part 3 - Comparative Analysis Summary**
- Design principles from intersection queries

---

#### A3. PrimeKG Summary
**Part 1 - Individual Analysis Summary**
- Table of all gene-disease shortest path scores
- Path length distribution
- Key findings (direct vs indirect associations)

**Part 2 - Combo Analysis Summary**
- Table of combo average scores
- Component contribution analysis

**Part 3 - Comparative Analysis Summary**
- Path structure insights
- Network topology patterns

---

### SECTION B: Cross-Method Comparisons

#### B1. Individual Analysis Cross-Comparison (Part 1 from all methods)

| Gene | Disease | BioBridge | ULTRA | PrimeKG | Path | Agreement |
|------|---------|-----------|-------|---------|------|-----------|
| ITGA4 | CD | 0.996 | 0.985 | 0.90 | 2 | HIGH |
| TYK2 | UC | 0.994 | 0.992 | 1.00 | 1 | HIGH |

**Concordant Findings** (all 3 methods agree - HIGH confidence):
- List findings where all methods show consistent strong/weak

**Discordant Findings** (methods disagree - investigate):
- List with possible explanations

---

#### B2. Combo Analysis Cross-Comparison (Part 2 from all methods)

| Combo | Disease | BioBridge | ULTRA | PrimeKG | BB Class | ULTRA Class | Agreement |
|-------|---------|-----------|-------|---------|----------|-------------|-----------|
| TYK2+JAK1 | CD | 0.995 | 0.990 | 0.855 | NEAR-ADD | NEAR-ADD | HIGH |

**Synergy Agreement** (BioBridge and ULTRA both show synergy)
**Dilution Agreement** (both show dilution)
**Discordant Classifications**

Note: PrimeKG combo = average of component scores (always NEAR-ADDITIVE)

---

#### B3. Comparative Analysis Cross-Comparison (Part 3 from all methods)

**Unified Design Principles** (supported by all methods):
1. Principle X - BioBridge, ULTRA, PrimeKG all support
2. ...

**Method-Specific Insights**:
- BioBridge-only: Embedding proximity patterns
- ULTRA-only: Intersection query insights
- PrimeKG-only: Path structure/topology insights

---

### SECTION C: Final Cross-Method Analysis

#### C1. Confidence Assessment Matrix

| Finding | BioBridge | ULTRA | PrimeKG | Confidence |
|---------|-----------|-------|---------|------------|
| ITGA4 top for CD | Strong | Strong | Direct | **HIGH** |
| TYK2+JAK1 synergy | Near-Add | Near-Add | 0.855 | **MODERATE** |

**Confidence Levels**:
- **HIGH**: All 3 methods agree
- **MODERATE**: 2 of 3 methods agree
- **LOW**: All methods disagree - requires investigation

#### C2. Method Agreement Statistics
- Proportion of findings where all 3 agree
- Proportion of 2/3 agreement
- Proportion of disagreement

#### C3. Final Recommendations
Based on cross-method consensus:
- Target prioritization (ranked by agreement)
- Biomarker combinations (validated by multiple methods)
- Design principles (unified from all methods)

#### C4. Methodological Recommendations
- When to use BioBridge: Semantic similarity, embedding-based predictions
- When to use ULTRA: Logical AND queries, complex multi-hop reasoning
- When to use PrimeKG: Direct/indirect path analysis, network topology
- Combined workflow: Use all 3 for highest confidence results

---

## Output Structure

### Mode: COMBO (nested entity2 input)
```
kgpred_<context>_<date>/
├── biobridge/
│   ├── Part1_Individual_Analysis.md
│   ├── Part2_Combo_Analysis.md
│   └── Part3_Comparative_Analysis.md
├── ultra/
│   ├── Part1_Individual_Analysis.md
│   ├── Part2_Combo_Analysis.md
│   └── Part3_Comparative_Analysis.md
├── primekg/
│   ├── Part1_Individual_Analysis.md
│   ├── Part2_Combo_Analysis.md
│   └── Part3_Comparative_Analysis.md
├── data/
│   ├── biobridge/
│   ├── ultra/
│   └── primekg/
└── Cross_Method_Comparison.md   # COMPREHENSIVE: includes ALL method summaries + comparisons
```

### Mode: INDIVIDUAL_ONLY (flat entity2 input)
```
kgpred_<context>_<date>/
├── biobridge/
│   └── Part1_Individual_Analysis.md
├── ultra/
│   └── Part1_Individual_Analysis.md
├── primekg/
│   └── Part1_Individual_Analysis.md
├── data/
│   ├── biobridge/
│   ├── ultra/
│   └── primekg/
└── Cross_Method_Comparison.md   # COMPREHENSIVE: includes ALL method summaries + comparisons
```

---

## Score Interpretation Guide

### BioBridge/ULTRA: Percentile Rank (pct_rank) Thresholds

| pct_rank | Interpretation | Assessment |
|----------|---------------|------------|
| ≥ 0.99 | Top 1% | **Very Strong** |
| 0.95-0.99 | Top 5% | **Strong** |
| 0.90-0.95 | Top 10% | Moderate-Strong |
| 0.80-0.90 | Top 20% | Moderate |
| 0.50-0.80 | Top 50% | Weak |
| < 0.50 | Below median | Very Weak |

### PrimeKG: Shortest Path Score Thresholds

| Path Length | Score (0.9^(n-1)) | Assessment |
|-------------|-------------------|------------|
| 1 (direct) | 1.00 | **Direct** - very strong |
| 2 | 0.90 | **Strong** - 1 intermediate |
| 3 | 0.81 | **Moderate-Strong** - 2 intermediates |
| 4 | 0.729 | **Moderate** - 3 intermediates |
| 5-6 | 0.59-0.66 | **Weak** - 4-5 intermediates |
| 7+ | < 0.53 | Very Weak - distant |
| No path | 0.00 | **None** - disconnected |

### Synergy/Dilution Thresholds (pct_rank delta)

| Delta (Δ) | Classification | Meaning |
|-----------|----------------|---------|
| > +0.02 | **SYNERGY** | Combo stronger than expected |
| -0.02 to +0.02 | NEAR-ADDITIVE | Combo = average of components |
| < -0.02 | **DILUTION** | Combo weaker than expected |

Note: PrimeKG combos use averaging, so delta = 0 by definition (always NEAR-ADDITIVE)

### Cross-Method Agreement Thresholds

| Agreement | Definition | Confidence |
|-----------|------------|------------|
| All 3 methods agree | Strong signal in same direction | **HIGH** |
| 2 of 3 methods agree | Majority consensus | **MODERATE** |
| All 3 methods disagree | No consensus | **LOW** - investigate |

---

## Reference Files

- `references/mcp_usage.md` - MCP tool documentation
- `references/report_templates.md` - Markdown report templates
- `references/combo_methods.md` - Mean embeddings vs intersection
- `references/interpretation_guide.md` - Score interpretation

## Key Principles

1. **Save outputs to working directory** - All reports/results go to user's current directory, NOT skill or MCP folders
2. **MOVE files IMMEDIATELY after each query** - After EACH MCP query, IMMEDIATELY move parquet/data files to `<working_dir>/data/<mcp_name>/` BEFORE processing or proceeding
3. **BioBridge: ALWAYS use explicit entity overrides** - MANDATORY: `override_head_name`/`override_head_names`, `override_head_type`, `override_tail_name`, `override_tail_type` must be provided for every BioBridge MCP call. Never rely on context inference alone.
4. **Use pct_rank as primary score** - For BioBridge and ULTRA; PrimeKG uses exponential decay score
5. **Report both metrics** - pct_rank (score) AND cos_sim/rank (reference)
6. **Filter by entity type** - ULTRA pct_rank uses only entity1_type predictions
7. **Detect input format** - Check if entity2 is nested (combos) or flat (individuals only)
8. **Run ALL THREE methods** - BioBridge, ULTRA, and PrimeKG for comprehensive comparison
9. **Calculate deltas** - Compare combo scores to individual means (if combos provided)
10. **Identify consistency** - Highlight where all 3 methods agree (high confidence)
11. **Flag discordance** - Note where methods disagree (needs investigation)
12. **Comprehensive final report** - Cross_Method_Comparison includes ALL method summaries
13. **Biological interpretation** - Provide mechanistic explanations
14. **Clinical relevance** - Include therapeutic implications
15. **PrimeKG batch processing** - Use parallel shortest path queries for efficiency

---

## Recommended Workflow: Direct Output with output_dir

**Use `output_dir` parameter to save files directly to your working directory.**

### Pattern for Each Query (Recommended)

```
FOR EACH entity/combo:
    1. Run MCP query with output_dir parameter
    2. Read the output file directly (no move needed!)
    3. Process/extract scores
    4. Proceed to next entity/combo
```

### Example: ULTRA Individual Gene Queries

```python
data_dir = f"./kgpred_{context}_{date}/data/ultra/individual"
os.makedirs(data_dir, exist_ok=True)

for gene in genes:
    # Step 1: Run query with output_dir - files saved directly!
    result = mcp__ultra-inference__predict_tail_entities(
        head_entity=gene,
        relation="associated with",
        output_dir=data_dir  # NEW: Direct output
    )

    # Step 2: Read directly - no move needed!
    df = pd.read_parquet(result['output_file'])
    scores = extract_disease_scores(df, disease_list)

    # Step 3: Record
    results[gene] = scores

    # Step 4: Loop continues to next gene
```

### Example: ULTRA Combo Queries

```python
data_dir = f"./kgpred_{context}_{date}/data/ultra/combo"
os.makedirs(data_dir, exist_ok=True)

for combo in combos:
    # Step 1: Run query with output_dir - files saved directly!
    result = mcp__ultraquery-inference__answer_complex_query(
        query_structure=build_intersection_query(combo),
        top_k=5,  # Only affects API response; parquet files contain ALL predictions
        output_dir=data_dir  # NEW: Direct output
    )

    # Step 2: Read directly - no move needed!
    df = pd.read_parquet(result['output_file_filtered'])

    # CRITICAL: Filter to entity1_type only and recalculate pct_rank
    df_disease = df[df['entity_type'] == 'disease'].copy()
    df_disease['rank'] = df_disease['score'].rank(ascending=False, method='first').astype(int)
    df_disease['pct_rank'] = 1.0 - (df_disease['rank'] / len(df_disease))

    scores = extract_disease_scores(df_disease, disease_list)

    # Step 3: Record
    combo_name = "_".join(combo)
    results[combo_name] = scores

    # Step 4: Loop continues to next combo
```

### Benefits of output_dir Parameter

1. **No file moves needed** - Files saved directly to working directory
2. **Simpler workflow** - One less step per query
3. **No volatility risk** - Files are in user-controlled location from start
4. **Backward compatible** - Omit output_dir to use legacy behavior

---

## Parquet Extraction Scripts

**Ready-to-use scripts for seamless score extraction from MCP-generated parquet files.**

Scripts are located in: `scripts/` folder

### ULTRA Parquet Column Schema

**ultra-inference** and **ultraquery-inference** use DIFFERENT column names:

#### ultra-inference Parquet Schema (Individual Queries)

| Column | Description |
|--------|-------------|
| `h_label` | Head entity ID (e.g., "NCBI:7297") |
| `h_name` | Head entity name (e.g., "TYK2") |
| `t_pred_label` | Tail entity ID (e.g., "MONDO:5011") |
| `t_pred_name` | Tail entity name (e.g., "Crohn disease") |
| `t_pred_score` | Model prediction score |
| `t_pred_type` | Tail entity type (e.g., "disease") |
| `rank` | **Global** rank across ALL entity types |
| `percentile_rank` | **Global** percentile across ALL types |

#### ultraquery-inference Parquet Schema (Combo/Intersection Queries)

| Column | Description |
|--------|-------------|
| `entity_id` | Entity ID (e.g., "MONDO:5011") |
| `entity_name` | Entity name (e.g., "Crohn disease") |
| `entity_type` | Entity type (e.g., "disease") |
| `score` | Model prediction score |
| `rank` | **Global** rank across ALL entity types |
| `percentile_rank` | **Global** percentile across ALL types |
| `filtered_rank` | Rank within schema-matched types |
| `filtered_percentile_rank` | Percentile within schema-matched types |

**CRITICAL**: For BOTH tools:
- Native `rank` and `percentile_rank` are across ALL entity types
- Native `filtered_percentile_rank` (ultraquery) combines ALL schema-matched entity types
- You MUST filter to your specific target entity type and recalculate pct_rank

### Available Scripts

#### 1. `scripts/extract_ultra_individual.py`
Extract scores from individual gene prediction parquets.

```bash
python scripts/extract_ultra_individual.py \
    --data-dir ./kgpred_IBD_2025-01-05/data/ultra/individual \
    --diseases "crohn,ulcerative colitis,inflammatory bowel disease" \
    --output ./kgpred_IBD_2025-01-05/data/ultra/individual_scores.json
```

#### 2. `scripts/extract_ultra_combo.py`
Extract scores from combo/intersection parquets.

```bash
python scripts/extract_ultra_combo.py \
    --data-dir ./kgpred_IBD_2025-01-05/data/ultra/combo \
    --diseases "crohn,ulcerative colitis,inflammatory bowel disease" \
    --output ./kgpred_IBD_2025-01-05/data/ultra/combo_scores.json
```

#### 3. `scripts/parquet_utils.py`
Utility functions for inspection and quick lookups.

```bash
# Inspect parquet structure
python scripts/parquet_utils.py inspect data/ultra/individual/TYK2_predictions.parquet

# Search for entities
python scripts/parquet_utils.py search data/ultra/individual/TYK2_predictions.parquet "crohn"

# Get single score
python scripts/parquet_utils.py score data/ultra/individual/TYK2_predictions.parquet "crohn"
```

### Inline Python Usage

```python
# Import utility functions
import sys
sys.path.append('/home/sagemaker-user/.claude/skills/kg-association/scripts')
from parquet_utils import get_ultra_disease_score, generate_score_matrix, calculate_synergy

# Quick score lookup
score = get_ultra_disease_score("TYK2_predictions.parquet", "crohn")
print(f"TYK2-Crohn: pct_rank={score['pct_rank']:.4f}")

# Generate gene x disease matrix
matrix = generate_score_matrix(
    data_dir="data/ultra/individual",
    genes=["TYK2", "JAK1", "ITGA4"],
    diseases=["crohn", "ulcerative colitis", "IBD"]
)
print(matrix.to_markdown())

# Calculate synergy/dilution
result = calculate_synergy(combo_score=0.9951, individual_scores=[0.9940, 0.9953])
print(f"Classification: {result['classification']}, Delta: {result['delta']:.4f}")
```

### Bash One-Liners

```bash
# View parquet schema and first rows
python3 -c "import pandas as pd; df=pd.read_parquet('$FILE'); print('Columns:', df.columns.tolist()); print(df.head(3))"

# Count entities by type
python3 -c "import pandas as pd; df=pd.read_parquet('$FILE'); print(df['t_pred_type'].value_counts())"

# Search for disease
python3 -c "import pandas as pd; df=pd.read_parquet('$FILE'); print(df[df['t_pred_name'].str.lower().str.contains('crohn')][['t_pred_name','t_pred_score','rank']].head(10))"
```
