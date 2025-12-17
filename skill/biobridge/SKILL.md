---
name: biobridge
description: Predict biomedical entity associations using BioBridge multimodal knowledge graph. Use this skill when users request to predict links between genes, diseases, drugs, phenotypes, or other biomedical entities, query association scores, match entities to the knowledge graph, or need to export and summarize prediction results. The skill supports flexible entity/relation inputs, entity matching with user confirmation, result export, and findings summarization.
---

# BioBridge Link Prediction

## Overview

This skill enables comprehensive prediction of biomedical entity associations using the BioBridge multimodal knowledge graph and neural link prediction models. The skill queries pre-trained embeddings to predict associations between genes, diseases, drugs, phenotypes, pathways, and other biomedical entities.

**Key Capabilities:**
- Predict associations between biomedical entities (gene-disease, drug-phenotype, etc.)
- **Gene signature analysis**: Query multiple entities simultaneously with mean embedding
- Match user-provided entity names to knowledge graph nodes with user confirmation
- Query with flexible parameters (head/tail entities, relation types, embeddings)
- Calculate association scores using neural retrieval (cosine similarity, percentile rank)
- Export results in multiple formats (CSV, JSON, Excel)
- Generate summaries interpreting prediction findings
- Support custom embedding models and relation types
- **S3 and local storage**: Read data from S3 URIs or local paths

## Workflow: Environment Check Before Tasks

**IMPORTANT:** Before executing any prediction, characterization, or analysis task, Claude should ensure the environment is ready:

### Quick Environment Check

```bash
# Check if environment is ready (silent, fast)
python scripts/ensure_env.py --quiet
```

Exit code 0 = ready, non-zero = needs setup

### Auto-Setup if Needed

If the environment check fails, automatically trigger setup:

```bash
bash scripts/setup_env.sh
```

This is a one-time process (2-5 minutes) that:
- Installs `uv` for fast package installation
- Creates `biobridge` conda environment with Python 3.11
- Installs all dependencies from `pyproject.toml`
- Verifies installation

### Claude's Task Execution Pattern

**For every user request:**
1. Run: `python scripts/ensure_env.py --quiet`
2. Check exit code
3. If 0 → proceed with task (e.g., `conda run -n biobridge python scripts/predict_link.py ...`)
4. If non-zero → run setup, then proceed with task

**Example:**
```bash
# Step 1: Check environment
python scripts/ensure_env.py --quiet

# Step 2: If ready, run task
conda run -n biobridge python scripts/predict_link.py \
  --head IL11 --head-type "gene/protein" --tail-type disease
```

## When to Use This Skill

Invoke this skill when users request:

**Link Prediction Queries:**
- "Predict diseases associated with [GENE]"
- "What is the association score between [GENE] and [DISEASE]?"
- "Find drugs that target [PROTEIN]"
- "Which phenotypes are linked to [DRUG]?"
- "Predict top associations for [ENTITY]"

**Gene Signature Analysis (NEW):**
- "Find diseases associated with the IL11 and GREM1 gene signature"
- "What pathways are enriched for genes [GENE1, GENE2, GENE3]?"
- "Predict associations for multiple genes: TP53, BRCA1, EGFR"

**Entity Matching:**
- "Match GREM1 to knowledge graph nodes"
- "Find the KG entity for Crohn's disease"
- "What node types exist for [ENTITY]?"

**Comparative Analysis:**
- "Compare association scores between IL11 and GREM1 for Crohn's disease"
- "Rank genes by their association with [DISEASE]"
- "Which has stronger evidence: [ENTITY1] or [ENTITY2] for [TARGET]?"

**Export and Reporting:**
- "Export prediction results to Excel"
- "Summarize the findings for [ENTITY] associations"
- "Generate a report of top predicted links"

## Critical Workflow: Entity Mapping Confirmation

**IMPORTANT:** Before executing any prediction, ALWAYS show the mapped entities and relations to the user and get their confirmation.

### Mandatory Confirmation Steps

**Step 1: Map entities to KG nodes**
```bash
conda run -n biobridge python scripts/match_entity.py [ENTITY_NAME] --type [TYPE]
```

**Step 2: Present mappings to user**
Show the user:
- Matched head entity: name, type, node_index, aliases
- Matched tail entity (if specified): name, type, node_index, aliases
- Selected relation: relation family, compatible with these types
- Any warnings: ambiguous matches, multiple candidates

**Step 3: Get user confirmation**
Ask the user explicitly:
"I've mapped your query as follows:
- Head: [ENTITY_NAME] (type: [TYPE], node_index: [INDEX])
- Tail: [ENTITY_NAME] (type: [TYPE], node_index: [INDEX])
- Relation: [RELATION]

Does this look correct? Reply 'yes' to proceed, or tell me how to adjust the mapping."

**Step 4a: If user confirms → proceed with prediction**
```bash
conda run -n biobridge python scripts/predict_link.py [confirmed parameters]
```

**Step 4b: If user disagrees → remap based on feedback**
- Ask clarifying questions about which entity/relation to use
- Run match_entity.py again with adjusted parameters
- Return to Step 2 with new mappings

### Example Confirmation Workflow

**User request:** "Predict association between IL11 and Crohn's disease"

**Claude actions:**
1. Map entities:
```bash
conda run -n biobridge python scripts/match_entity.py "IL11" --type "gene/protein"
# Returns: IL11 (node_index: 12345, type: gene/protein)

conda run -n biobridge python scripts/match_entity.py "Crohn disease" --type disease
# Returns: Crohn disease (mondo_id: 5011, node_index: 37784, type: disease)
```

2. Present to user:
"I've found these entities in the knowledge graph:
- **Head:** IL11 (gene/protein, node_index: 12345)
- **Tail:** Crohn disease (disease, MONDO:5011, node_index: 37784)
- **Relation:** associated with (gene/protein → disease)

Shall I proceed with this mapping to predict the association score?"

3. Wait for user confirmation before running predict_link.py

## MCP Tool Usage (Primary Method)

When the BioBridge MCP server is available, use the `mcp__biobridge_github__predict_associations` tool as the primary method for predictions. This tool provides:

- Automatic entity matching with LLM assistance
- Fast-path optimization when entities are explicit
- Gene signature analysis with multiple head entities
- Comprehensive validation and error handling

### Single Entity Query

```python
# Example: Predict diseases for GREM1
mcp__biobridge_github__predict_associations(
    context="Find diseases associated with GREM1 (gremlin 1). Focus on inflammatory bowel disease, Crohn disease variants.",
    override_head_name="GREM1",
    override_head_type="gene/protein",
    relation_hint="associated with",
    topk=25
)
```

### Gene Signature Query (NEW Feature)

For gene signature analysis, use `override_head_names` to pass multiple entities:

```python
# Example: Predict diseases for IL11 + GREM1 signature
mcp__biobridge_github__predict_associations(
    context="Find diseases associated with IL11 and GREM1 gene signature. Focus on IBD and inflammatory conditions.",
    override_head_names=["IL11", "GREM1"],
    override_head_type="gene/protein",
    tail_type="disease",
    relation_hint="associated with",
    topk=50
)
```

**How gene signatures work:**
- Computes embeddings for each entity separately
- Calculates mean embedding across all entities
- Uses the mean embedding for prediction
- Useful for pathway analysis, disease mechanism studies, and multi-gene biomarkers

**When to use gene signatures:**
- User provides multiple genes/proteins (e.g., "IL11, GREM1, and TP53")
- Analyzing gene sets or pathways
- Comparing multi-gene biomarkers
- Studying disease mechanisms involving multiple factors

### Pair Validation Query

```python
# Example: Validate specific association
mcp__biobridge_github__predict_associations(
    context="Is GREM1 associated with Crohn disease?",
    override_head_name="GREM1",
    override_head_type="gene/protein",
    override_tail_name="Crohn disease",
    override_tail_type="disease",
    relation_hint="associated with",
    topk=1
)
```

### MCP Tool Parameters

**Required:**
- `context` (str): Natural language description of the query

**Optional:**
- `override_head_name` (str): Single head entity name
- `override_head_names` (list[str]): **Multiple head entities for signature analysis**
- `override_head_type` (str): Head entity type (**required** when using `override_head_names`)
- `override_tail_name` (str): Specific tail entity name
- `override_tail_type` (str): Tail entity type
- `relation_hint` (str): Relation family (e.g., "associated with", "treats")
- `topk` (int): Number of results (default: 25, max: 100)
- `slidewindow` (bool): Use slidewindow embeddings for proteins
- `include_relation_catalog` (bool): Include relation IDs (default: True)
- `include_debug` (bool): Include debug info (default: False)

**Note:** `override_head_name` and `override_head_names` are mutually exclusive.

## Knowledge Graph Overview

The BioBridge knowledge graph (PrimeKG) contains:
- **8.1M edges** (relationships)
- **10 entity types**: gene/protein, drug, disease, anatomy, effect/phenotype, biological_process, molecular_function, cellular_component, pathway, exposure
- **18 relation types**: associated with, treats, interacts with, ppi, side effect, contraindication, etc.

To characterize your knowledge graph:

```bash
conda run -n biobridge python scripts/characterize_kg.py
conda run -n biobridge python scripts/characterize_kg.py --show-combinations
conda run -n biobridge python scripts/characterize_kg.py --export kg_stats.json
```

## Quick Start

### Automatic Environment Setup

The skill automatically checks and sets up the `biobridge` conda environment when needed. On first use, it will:

1. Check if conda is available
2. Check if the `biobridge` environment exists
3. If missing, install `pixi` (fast cross-platform package manager)
4. Create the environment with Python 3.11
5. Install all required packages from `pixi.toml`

**To manually set up or verify the environment:**

```bash
# Automatic setup (recommended)
bash scripts/setup_env.sh

# Or check environment status
python scripts/ensure_env.py

# Check without setting up
python scripts/ensure_env.py --check-only
```

**Required dependencies** (automatically installed):
- pytorch >= 2.9.1
- numpy >= 2.3.5
- pandas >= 2.3.3
- scipy >= 1.16.3
- transformers >= 4.57.3
- fastmcp >= 2.13.3
- scikit-learn >= 1.8.0
- boto3 >= 1.42.10 (for S3 support)
- python-dotenv (for environment configuration)

All dependencies are specified in `pixi.toml` and installed automatically via the pixi package manager.

### Manual Environment Activation

If the environment already exists, activate it manually:

```bash
conda activate biobridge
```

### Basic Prediction Patterns

**Pattern 1: Map entities first (RECOMMENDED)**
```bash
# Step 1: Match head entity
conda run -n biobridge python scripts/match_entity.py GREM1 --type "gene/protein"

# Step 2: Match tail type (or specific tail)
conda run -n biobridge python scripts/match_entity.py "Crohn disease" --type disease

# Step 3: Show user and get confirmation

# Step 4: Run prediction with confirmed parameters
conda run -n biobridge python scripts/predict_link.py \
  --head GREM1 \
  --head-type "gene/protein" \
  --tail "Crohn disease" \
  --tail-type disease \
  --relation "associated with"
```

**Pattern 2: Explore matches interactively**
```bash
# Show all matches for ambiguous entity
conda run -n biobridge python scripts/match_entity.py "gremlin"
# Present options to user, let them choose

# Proceed with selected match
conda run -n biobridge python scripts/predict_link.py --head "GREM1" ...
```

## Entity Types and Relations

### Supported Entity Types

The knowledge graph supports the following canonical entity types:
- **gene/protein** - Genes and proteins
- **disease** - Diseases and disorders
- **drug** - Chemical compounds and drugs
- **biologics_drug** - Biologic therapeutic agents
- **effect/phenotype** - Effects and phenotypes
- **pathway** - Biological pathways
- **biological_process** - GO biological processes
- **molecular_function** - GO molecular functions
- **cellular_component** - GO cellular components
- **anatomy** - Anatomical terms
- **exposure** - Environmental exposures

### Common Relations

- **associated with** - General association (gene-disease, protein-phenotype)
- **treats** - Therapeutic relationship (drug-disease)
- **interacts with** - Protein-protein interactions
- **side effect** - Drug adverse effects
- **regulates** - Regulatory relationships
- **participates in** - Pathway/process participation

### Valid Relation Combinations

Not all entity type combinations support all relations. Common valid combinations:

- gene/protein → disease: **associated with**
- gene/protein → pathway: **participates in**, **regulates**
- gene/protein → gene/protein: **interacts with**
- drug → disease: **treats**, **associated with**
- drug → effect/phenotype: **side effect**, **associated with**

The match_entity.py and predict_link.py scripts will validate combinations and suggest alternatives if invalid.

## Prediction Queries

### Predict Top Associations

Use `scripts/predict_link.py` to predict top associations for a head entity.

**Basic usage:**
```bash
conda run -n biobridge python scripts/predict_link.py \
  --head [ENTITY_NAME] \
  --head-type [TYPE] \
  --tail-type [TYPE] \
  [OPTIONS]
```

**Required parameters:**
- `--head [NAME]` - Head entity name (e.g., "GREM1", "TP53")
- `--head-type [TYPE]` - Head entity type (e.g., "gene/protein", "drug")
- `--tail-type [TYPE]` - Target entity type to predict

**Optional parameters:**
- `--tail [NAME]` - Specific tail entity to validate (for pair queries)
- `--relation [REL]` - Relation type hint (e.g., "associated with")
- `--context [TEXT]` - Additional context to guide entity matching
- `--topk [N]` - Number of top predictions to return (default: 25)
- `--export [FILE]` - Export results to file
- `--format [FMT]` - Export format: csv, json, or excel (default: csv)
- `--slidewindow` - Use slidewindow embeddings for proteins
- `--no-slidewindow` - Disable slidewindow embeddings
- `--confirm` - Skip interactive confirmation (use after manual entity matching)

**Example: Full workflow with confirmation**
```bash
# Step 1: Match entity
conda run -n biobridge python scripts/match_entity.py GREM1 --type "gene/protein"
# Output: GREM1 (gene/protein, node_index: 45678)

# Step 2: Present to user
# "Found: GREM1 (gene/protein). Proceed? [yes/no]"

# Step 3: Run prediction (after user confirms)
conda run -n biobridge python scripts/predict_link.py \
  --head GREM1 \
  --head-type "gene/protein" \
  --tail-type disease \
  --relation "associated with" \
  --topk 50 \
  --export grem1_diseases.csv \
  --confirm
```

**Example: Validate specific pair**
```bash
# Step 1: Match both entities
conda run -n biobridge python scripts/match_entity.py IL11 --type "gene/protein"
conda run -n biobridge python scripts/match_entity.py "Crohn disease" --type disease

# Step 2: Confirm with user

# Step 3: Validate pair
conda run -n biobridge python scripts/predict_link.py \
  --head IL11 \
  --head-type "gene/protein" \
  --tail "Crohn disease" \
  --tail-type disease \
  --context "inflammatory bowel disease" \
  --topk 1 \
  --confirm
```

**Output interpretation:**
- **cos_sim** - Cosine similarity score (0-1, higher is stronger)
- **pct_rank** - Percentile rank (0-1, higher means top percentile)
- **node_index** - Internal KG node identifier
- **node_name/mondo_name** - Entity name
- **node_id/mondo_id** - External database identifier

### Entity Matching and Validation

Use `scripts/match_entity.py` to find knowledge graph nodes matching a query.

**Usage:**
```bash
conda run -n biobridge python scripts/match_entity.py [ENTITY_NAME] [OPTIONS]
```

**Parameters:**
- `[ENTITY_NAME]` - Entity to search for
- `--type [TYPE]` - Filter by entity type
- `--limit [N]` - Maximum matches to return (default: 10)
- `--show-details` - Show synonyms, node_index, and IDs
- `--format json` - JSON output format

**Example:**
```bash
# Find all matches for "Crohn"
conda run -n biobridge python scripts/match_entity.py "Crohn"

# Find disease matches only
conda run -n biobridge python scripts/match_entity.py "Crohn disease" --type disease --show-details

# Export matches for programmatic use
conda run -n biobridge python scripts/match_entity.py "IL11" --type "gene/protein" --format json > il11_matches.json
```

**When to use:**
- **ALWAYS** before running predict_link.py (mandatory workflow step)
- User provides ambiguous entity names
- Need to verify entity exists in knowledge graph
- Want to see alternative names/synonyms
- Exploring entity types available for a term

### Presenting Entity Mappings to Users

**Good presentation format:**
```
I've mapped your entities to the knowledge graph:

**Head Entity:**
- Name: IL11
- Type: gene/protein
- Node Index: 12345
- Aliases: Interleukin 11, IL-11

**Tail Entity:**
- Name: Crohn disease
- Type: disease
- Node Index: 37784
- IDs: MONDO:5011
- Aliases: Crohn's disease, Regional enteritis

**Relation:**
- associated with (gene/protein → disease)

This mapping will query: "IL11 → associated with → Crohn disease"

Does this look correct? Reply 'yes' to proceed or tell me how to adjust.
```

## Export and Summarization

### Export Prediction Results

Export predictions to various formats for downstream analysis.

**Supported formats:**
- **CSV** - Tabular format, one row per prediction
- **JSON** - Structured format with full metadata
- **Excel** - Multi-sheet workbook with summary statistics

**CSV output columns:**
- Head_Entity, Head_Type, Head_Node_Index
- Tail_Entity, Tail_Type, Tail_Node_Index
- Relation
- Cosine_Similarity, Percentile_Rank
- Entity_ID (MONDO_ID, etc.)

**JSON output structure:**
```json
{
  "query": {
    "head": {
      "name": "GREM1",
      "type": "gene/protein",
      "node_index": 45678
    },
    "tail_type": "disease",
    "relation": "associated with"
  },
  "results": [
    {
      "tail_name": "Crohn disease",
      "tail_node_index": 37784,
      "mondo_id": 5011,
      "cos_sim": 0.039,
      "pct_rank": 0.586
    }
  ],
  "metadata": {
    "timestamp": "2025-12-11T19:30:00Z",
    "model": "BioBridge-6layer",
    "embedding": "esm2b_unimo_pubmedbert"
  }
}
```

### Generate Summary Reports

Use `scripts/summarize_predictions.py` to generate human-readable summaries.

**Usage:**
```bash
conda run -n biobridge python scripts/summarize_predictions.py [RESULTS_FILE] [OPTIONS]
```

**Parameters:**
- `[RESULTS_FILE]` - Input file (CSV or JSON from predict_link.py)
- `--threshold [SCORE]` - Minimum cosine similarity threshold (default: 0.1)
- `--top [N]` - Show top N predictions (default: 10)
- `--output [FILE]` - Save summary to file

**Example:**
```bash
# Generate summary from predictions
conda run -n biobridge python scripts/summarize_predictions.py grem1_diseases.csv \
  --threshold 0.05 \
  --top 20 \
  --output grem1_summary.txt
```

**Summary includes:**
- Query details (entities, types, relation)
- Entity mapping confirmation status
- Total predictions found
- Distribution of association scores
- Top predictions with interpretation
- Confidence assessment (based on score thresholds)

## Understanding Prediction Results

### Interpreting Scores

**Cosine Similarity (cos_sim):**
- Range: -1 to 1 (typically 0 to 1 for positive associations)
- **Strong (>0.2):** High confidence association
- **Moderate (0.1-0.2):** Moderate evidence, worth investigating
- **Weak (0.05-0.1):** Low evidence, suggestive but not conclusive
- **Very Weak (<0.05):** Minimal evidence or no clear association

**Percentile Rank (pct_rank):**
- Range: 0 to 1 (0 to 100th percentile)
- **Top 10% (>0.9):** Among strongest associations for this entity
- **Top 25% (>0.75):** Well-established connection
- **Middle (0.4-0.6):** Average association strength
- **Bottom (<0.3):** Weaker than most associations

**Example interpretation:**
```
IL11 → Crohn disease
cos_sim: 0.254 (25.4%)
pct_rank: 0.845 (84.5 percentile)

Interpretation: Strong association. IL11 shows substantial evidence
of connection to Crohn disease, ranking in the top 15% of associations.
This is a well-established link supported by the knowledge graph.
```

### Comparing Predictions

When comparing multiple entities:

**GREM1 → Crohn disease:**
- cos_sim: 0.039, pct_rank: 0.586
- **Interpretation:** Weak to moderate signal, middle of distribution

**IL11 → Crohn disease:**
- cos_sim: 0.254, pct_rank: 0.845
- **Interpretation:** Strong signal, top 15% of associations

**Comparison:** IL11 has ~6.5x stronger association score with Crohn disease
compared to GREM1, suggesting IL11 has more substantial evidence linking
it to this condition in the knowledge graph.

## Custom Paths and Model Configuration

### Using Custom Data and Models

The skill supports custom paths for all data files and model checkpoints, enabling:
- Testing different knowledge graph versions
- Using alternative embeddings
- Switching between model checkpoints
- Working with custom entity/relation databases

**Command-line options:**
```bash
--kg-path [PATH]          # Custom knowledge graph CSV file
--nodes-path [PATH]       # Custom nodes CSV file
--model-ckpt [PATH]       # Custom model checkpoint (.bin file)
--embedding-dir [PATH]    # Custom embeddings directory
```

### Example: Using Custom Knowledge Graph

```bash
conda run -n biobridge python scripts/predict_link.py \
  --head IL11 \
  --head-type "gene/protein" \
  --tail-type disease \
  --kg-path /path/to/custom_kg.csv \
  --nodes-path /path/to/custom_nodes.csv
```

### Example: Using Custom Model Checkpoint

```bash
conda run -n biobridge python scripts/predict_link.py \
  --head TP53 \
  --head-type "gene/protein" \
  --tail-type disease \
  --model-ckpt /path/to/checkpoints/model_epoch100.bin \
  --embedding-dir /path/to/embeddings/slidewindow
```

### Default Paths

If custom paths are not specified, the skill uses these defaults:

**Knowledge Graph:**
- KG: `/home/sagemaker-user/biobridge/bbridge/data/PrimeKG/kg.csv`
- Nodes: `/home/sagemaker-user/biobridge/bbridge/data/PrimeKG/nodes.csv`

**Model:**
- Checkpoint: `/home/sagemaker-user/biobridge/bbridge/checkpoints/biobridge_6layer_100epochs_slidewindow.bin`
- Config: `/home/sagemaker-user/biobridge/bbridge/checkpoints/model_config.json`

**Embeddings:**
- Directory: `/home/sagemaker-user/biobridge/bbridge/data/slidewindow_esm2b_unimo_pubmedbert`

### S3 Storage Support (NEW)

The MCP server now supports reading data from Amazon S3 in addition to local paths. This enables:
- Cloud-based data storage and sharing
- Scalable deployment across multiple instances
- Version control of knowledge graphs and embeddings

**S3 Configuration:**

Set the `BIOBRIDGE_SRC_DIR` environment variable to an S3 URI:

```bash
export BIOBRIDGE_SRC_DIR="s3://your-bucket/biobridge/"
```

**S3 Path Format:**
- All paths support S3 URIs: `s3://bucket-name/path/to/file`
- The server automatically detects S3 URIs and uses boto3 for access
- Requires AWS credentials configured via environment or IAM role

**Example S3 Structure:**
```
s3://your-bucket/biobridge/
├── data/
│   ├── PrimeKG/
│   │   ├── kg.csv
│   │   └── nodes.csv
│   ├── embeddings/
│   │   └── esm2b_unimo_pubmedbert/
│   │       ├── protein.pkl
│   │       ├── disease.pkl
│   │       └── ...
│   └── Processed/
│       ├── protein.csv
│       ├── disease.csv
│       └── ...
└── ckpt/
    └── model_6layer_100epoch/
        ├── model.bin
        └── model_config.json
```

**Benefits:**
- Centralized data management
- No need to replicate large files across instances
- Easy version control and rollback
- Automatic failover and redundancy

### Characterizing Custom Knowledge Graphs

When using a custom KG, characterize it first to understand entity and relation types:

```bash
# Characterize custom KG
conda run -n biobridge python scripts/characterize_kg.py \
  --kg-path /path/to/custom_kg.csv

# Export statistics
conda run -n biobridge python scripts/characterize_kg.py \
  --kg-path /path/to/custom_kg.csv \
  --export custom_kg_stats.json \
  --show-combinations
```

**Key statistics provided:**
- Total edges (relationships)
- Entity type distribution
- Relation type frequencies
- Valid entity-relation combinations

Use this information to:
1. Verify entity types match your queries
2. Confirm relation types are available
3. Understand coverage and completeness
4. Validate data quality

## Troubleshooting

### Entity Not Found

**Symptom:** "Entity not found in knowledge graph"

**Solutions:**
1. Try alternative names/synonyms (e.g., "CDKN1A" vs "P21")
2. Check entity type is correct (gene/protein vs drug)
3. Run match_entity.py to explore available matches
4. Provide additional context with `--context` flag
5. Check spelling and capitalization

### Multiple Entity Matches

**Symptom:** match_entity.py returns multiple candidates

**Action:**
1. Present ALL matches to the user with details
2. Ask user to select the correct match
3. Use the selected match's node_index or exact name in predict_link.py

**Example user interaction:**
```
Claude: "I found multiple matches for 'IL':
1. IL11 (gene/protein, node_index: 12345)
2. IL6 (gene/protein, node_index: 23456)
3. IL1B (gene/protein, node_index: 34567)

Which entity did you mean?"

User: "IL11"

Claude: "Got it, using IL11 (node_index: 12345). Proceeding with prediction..."
```

### Incompatible Entity-Relation Combination

**Symptom:** Prediction returns empty or error about invalid combination

**Solutions:**
1. Check if the entity types support the specified relation
2. Try removing --relation flag (let system auto-select)
3. Consult the "Valid Relation Combinations" section
4. Inform user the combination isn't supported and suggest alternatives

## Key Principles

1. **ALWAYS match entities first and get user confirmation** before prediction
2. **Present mappings clearly** with entity names, types, node indices, and aliases
3. **Wait for user agreement** before proceeding with predict_link.py
4. **Remap on disagreement** - don't assume the first match is correct
5. **Provide entity types explicitly** for better matching accuracy
6. **Use context when available** to disambiguate entity names
7. **Interpret scores in context** - compare to baseline/other entities
8. **Export results for analysis** - scores alone don't tell the full story
9. **Generate summaries for presentations** - raw predictions need interpretation
10. **Validate important findings** - high scores warrant literature verification

## Resources

This skill includes:

### scripts/

Executable Python tools for BioBridge predictions:

- **characterize_kg.py** - Analyze knowledge graph structure (entity types, relations, combinations)
- **match_entity.py** - Quick entity name exploration tool
- **predict_link.py** - Main prediction tool using MCP predict_associations with LLM-based entity matching
- **summarize_predictions.py** - Generate human-readable summary reports from prediction results

Execute these scripts via conda environment. Scripts use the existing BioBridge infrastructure without requiring additional package installations.

---

## Technical Details

### Model Architecture

**Model:** BioBridge multimodal KG with neural link prediction (6-layer transformer, 100 epochs)

**Embeddings:** Domain-specific pre-trained models for optimal biomedical performance:
- **ESM2** (proteins): Large-scale protein language model from Meta AI
  - Trained on 65M protein sequences
  - Captures structural and functional protein properties
  - Superior for protein-protein interactions and gene associations
- **UniMol** (drugs): Unified molecular representation model
  - Trained on molecular structures and properties
  - Optimized for drug-drug and drug-target interactions
- **PubMedBERT** (text/entities): BERT pre-trained on PubMed abstracts
  - Domain-adapted for biomedical text
  - Better understanding of disease, phenotype, and pathway names

**Why domain-specific over general large models?**
1. **Accuracy**: Domain-specific models capture biomedical nuances better than general-purpose LLMs
2. **Efficiency**: Smaller, focused models run faster with lower computational cost
3. **Proven performance**: These embeddings have been validated on biomedical benchmarks
4. **Multimodal integration**: Different entity types need different embedding strategies (sequences, structures, text)

**Knowledge Graph:** PrimeKG with slidewindow augmentation
- 8.1M edges across 10 entity types
- Slidewindow embeddings for proteins capture local sequence context

**Prediction Method:** Neural retrieval via entity projection and transformation
- Entity embeddings → latent space projection
- Relation-specific transformations
- Cosine similarity ranking

**Scoring:** Cosine similarity in embedding space, percentile ranking

**Entity Matching:** LLM-based (OpenAI) for intelligent entity name resolution
**Validation:** Automatic type compatibility checking for entity-relation combinations
