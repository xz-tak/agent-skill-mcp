---
name: drugnomeai
description: Run, configure, and analyze DrugnomeAI gene druggability predictions. Use this skill when users mention drugnomeai, drugnome, gene druggability, druggability prediction, PU learning for druggability, assess_gene, cross-run assessment, druggability score, consensus tier, composite score, druggable targets, run drugnomeai, target druggability, novel drug targets, modality prediction, druggability rank, druggability probability, or request druggability info/recommendations for specific gene targets. Also invoke for inspecting the master druggability table, querying top-N druggable candidates, filtering by modality or consensus tier, or any request involving predicting which genes are druggable. Even if the user just asks "is X druggable?" or "druggability of X", use this skill.
---

# DrugnomeAI Skill

## Overview

DrugnomeAI is an ensemble PU (Positive-Unlabelled) learning framework that predicts gene druggability across 20,080 human genes using 435 features from 34 data sources. It trains 8 classifiers in parallel across 11 complementary run configurations spanning 3 assessment dimensions:

| Dimension | Runs | Purpose |
|-----------|------|---------|
| PHAROS knowledge level | `pharos_tclin`, `pharos_tchem`, `pharos_tclin_tchem` | Drug-target maturity |
| Tractability tier | `tier_1`, `tier_1_2`, `tier_1_2_3A`, `tier_1_2_3A_3B`, `tclin_tier1_intersect` | Structural druggability |
| Therapeutic modality | `modality_small_mol`, `modality_antibody`, `modality_protac` | Best treatment approach |

**8 classifiers:** Extra Trees, Random Forest, SVC, Gradient Boosting, XGBoost, DNN (PyTorch), Stacking, Gaussian Naive Bayes

The cross-run assessment merges all 11 runs into a **master druggability table** with composite scores, consensus tiers, modality predictions, and novelty flags.

**Biological annotation reference:** See `references/annotations.md` in this skill directory for full definitions of PHAROS TDL categories, druggability tiers, therapeutic modalities, PU learning terms, and per-run interpretation patterns.

**Codebase:** `~/claude_code/DrugnomeAI`
**Conda env:** `drugnome`
**Pre-built master table:** `~/.claude/skills/drugnomeai/references/druggability_master_table.csv`
**Master table column interpretation:** `~/.claude/skills/drugnomeai/references/drugnomeai_mastertable_interpretation.md` — read this reference when interpreting master table columns, understanding what each run/metric means, or explaining results to the user
**Scoring scripts:** `~/.claude/skills/drugnomeai/scripts/drugnome_combo_score.py`, `write_report.py`

---

## Environment Check

Before any operation, verify the environment:

```bash
# 1. Check conda env
conda env list | grep drugnome

# 2. Check codebase
ls ~/claude_code/DrugnomeAI/drugnome_ai/modules/main/__main__.py

# 3. Check master table exists (for query operations)
ls ~/.claude/skills/drugnomeai/references/druggability_master_table.csv
```

If the conda env or codebase is missing, use AskUserQuestion:
> "DrugnomeAI environment not found. Install from https://github.com/xz-tak/DrugnomeAI-xz ?"

If approved:
```bash
conda create -n drugnome python=3.10 -y
conda activate drugnome
git clone https://github.com/xz-tak/DrugnomeAI-xz.git ~/claude_code/DrugnomeAI
cd ~/claude_code/DrugnomeAI && pip install -e .
```

If the master table is missing at the skill references path, use AskUserQuestion:
> "Master table not found at ~/.claude/skills/drugnomeai/references/druggability_master_table.csv. Copy from DrugnomeAI output or regenerate via Section F."

---

## A: Configure & Run Training

### 11 Run Configurations

| # | Run Name | CLI Flags | Positive Genes |
|---|----------|-----------|----------------|
| 1 | `pharos_tclin` | `-p tclin` | 613 |
| 2 | `pharos_tchem` | `-p tchem` | 1,598 |
| 3 | `pharos_tclin_tchem` | `-p tclin tchem` | 2,211 |
| 4 | `tier_1` | `-t 1` | 1,427 |
| 5 | `tier_1_2` | `-t 1 2` | 2,109 |
| 6 | `tier_1_2_3A` | `-t 1 2 3A` | 2,979 |
| 7 | `tier_1_2_3A_3B` | `-t 1 2 3A 3B` | 4,479 |
| 8 | `tclin_tier1_intersect` | `-k misc/gene_lists/tclin_tier1_intersect.txt` | 525 |
| 9 | `modality_small_mol` | `-k misc/gene_lists/small_moelcules_genes.txt` | 885 |
| 10 | `modality_antibody` | `-k misc/gene_lists/antibody_genes.txt` | 251 |
| 11 | `modality_protac` | `-k misc/gene_lists/protac.txt` | 266 |

### Launch a Single Run

Always use tmux for long-running training jobs. The DNN classifier uses GPU automatically via `find_free_gpu()` — it picks the highest-index free GPU with <500 MiB used and 0% utilization. Override with `DRUGNOME_DNN_DEVICE` if needed.

```bash
# Single run example (pharos_tclin, all 8 classifiers, 10 iterations, 8 threads)
tmux new-session -d -s drugnome_tclin "cd ~/claude_code/DrugnomeAI && \
  conda run -n drugnome \
  drugnomeai -o ./output03182026/pharos_tclin -p tclin \
  -s et rf gb svc xgb dnn stack nb -n 8 -i 10 -r all \
  2>&1 | tee ./output03182026/pharos_tclin.log"
```

### Launch All 11 Runs (batch script)

```bash
tmux new-session -d -s drugnome_all "cd ~/claude_code/DrugnomeAI && \
  conda run -n drugnome bash run_all_iterations.sh \
  2>&1 | tee ./output03182026/run_all.log"
```

### GPU Configuration

```bash
# Check GPU availability
nvidia-smi

# Override GPU selection (before launching drugnomeai)
export DRUGNOME_DNN_DEVICE=cuda:0   # Specific GPU
export DRUGNOME_DNN_DEVICE=cuda:2   # Different GPU
export DRUGNOME_DNN_DEVICE=cpu      # Force CPU
```

The auto-detection logic (`dnn_pytorch.py:find_free_gpu()`) queries `nvidia-smi` for GPUs with memory <500 MiB and 0% utilization, then selects the highest-index free GPU. If no GPU is available, it falls back to CPU.

---

## B: Monitor Running Jobs

### Check tmux sessions
```bash
tmux ls
```

### Tail logs
```bash
# Follow a specific run's log
tmux send-keys -t drugnome_tclin "" ""   # Check if alive
tail -f ~/claude_code/DrugnomeAI/output03182026/run_all_iterations.log
```

### Detect current stage
Grep the log for stage markers:
```bash
# Which run is currently active?
grep -E "^(Starting|COMPLETED|FAILED):" ~/claude_code/DrugnomeAI/output03182026/run_all_iterations.log | tail -5

# Which classifier iteration?
grep "Iteration:" ~/claude_code/DrugnomeAI/output03182026/pharos_tclin.log | tail -3
```

### Parse AUC results
```bash
conda run -n drugnome python3 -c "
import glob, pandas as pd
files = glob.glob('output03182026/pharos_tclin/supervised-learning/PU_*.evaluation_metrics.tsv')
for f in sorted(files):
    clf = f.split('PU_')[1].split('.')[0]
    df = pd.read_csv(f, sep='\t', index_col=0)
    print(f'{clf:30s} median AUC = {df.AUC.median():.4f}')
"
```

---

## C: Analyze Single-Run Results

After a run completes, its predictions are in `<output_dir>/<run_name>/Gene-Predictions/`.

### Read merged predictions
```bash
conda run -n drugnome python3 -c "
import pandas as pd
df = pd.read_csv('output03182026/pharos_tclin/Gene-Predictions/AllClassifiers.Merged.drugnome_ai_predictions.csv')
print(f'Genes: {len(df)}, Columns: {list(df.columns)}')
print(df.sort_values('drugnome_ai_proba', ascending=False).head(20).to_string(index=False))
"
```

### Compare classifiers within a run
```bash
conda run -n drugnome python3 -c "
import glob, pandas as pd
run = 'output03182026/pharos_tclin'
files = glob.glob(f'{run}/supervised-learning/PU_*.evaluation_metrics.tsv')
rows = []
for f in sorted(files):
    clf = f.split('PU_')[1].split('.')[0]
    df = pd.read_csv(f, sep='\t', index_col=0)
    rows.append({'Classifier': clf, 'Median_AUC': df.AUC.median(), 'Std_AUC': df.AUC.std()})
print(pd.DataFrame(rows).sort_values('Median_AUC', ascending=False).to_string(index=False))
"
```

### Top-N novel genes from a single run
```bash
conda run -n drugnome python3 -c "
import pandas as pd
df = pd.read_csv('output03182026/pharos_tclin/Gene-Predictions/AllClassifiers.Merged.drugnome_ai_predictions.csv')
novel = df[df['known_gene'] == 0].sort_values('drugnome_ai_proba', ascending=False)
print(novel[['Gene_Name', 'drugnome_ai_proba', 'drugnome_ai_perc']].head(20).to_string(index=False))
"
```

---

## D: Inspect Master Table & Query Targets

This is the primary interactive module. The master table at `~/.claude/skills/drugnomeai/references/druggability_master_table.csv` contains all 20,080 genes with predictions from all 11 runs plus derived columns.

### Load master table
```bash
conda run -n drugnome python3 -c "
import pandas as pd
from pathlib import Path
master = pd.read_csv(Path('~/.claude/skills/drugnomeai/references/druggability_master_table.csv').expanduser())
print(f'Loaded: {len(master)} genes x {len(master.columns)} columns')
print('Derived columns:', [c for c in master.columns if c in [
    'pharos_mean_proba','tier_mean_proba','best_modality','best_modality_proba',
    'modality_specificity','composite_score','consensus_tier','n_runs_top_decile',
    'is_novel_everywhere','n_runs_above_75perc']])
"
```

### Top-N druggable targets
```python
import pandas as pd
from pathlib import Path
master = pd.read_csv(Path('~/.claude/skills/drugnomeai/references/druggability_master_table.csv').expanduser())

# Top N overall
top = master.sort_values('composite_score', ascending=False).head(N)

# Top N novel only
novel_top = master[master.is_novel_everywhere == 1].sort_values('composite_score', ascending=False).head(N)

# Top N by specific modality
sm_top = master[master.best_modality == 'Small Molecule'].sort_values('composite_score', ascending=False).head(N)
ab_top = master[master.best_modality == 'Antibody'].sort_values('composite_score', ascending=False).head(N)
protac_top = master[master.best_modality == 'PROTAC'].sort_values('composite_score', ascending=False).head(N)
```

### Filter by consensus tier
```python
high_tier = master[master.consensus_tier == 'High'].sort_values('composite_score', ascending=False)
moderate_tier = master[master.consensus_tier == 'Moderate'].sort_values('composite_score', ascending=False)
```

### Specific target query ("is GREM1 druggable?")

Use the built-in `assess_gene()` function for the richest output:
```bash
conda run -n drugnome python3 -c "
import sys; sys.path.insert(0, '.')
from drugnome_ai.post_analysis.druggability_assessment import *
from pathlib import Path
master = build_master_table(Path('output03182026'))
master = add_derived_columns(master)
print(assess_gene('GREM1', master))
"
```

Or query the pre-built CSV directly for quick lookups:
```bash
conda run -n drugnome python3 -c "
import pandas as pd
from pathlib import Path
master = pd.read_csv(Path('~/.claude/skills/drugnomeai/references/druggability_master_table.csv').expanduser())
gene = master[master.Gene_Name == 'GREM1'].iloc[0]
print(f\"Gene: {gene['Gene_Name']}\")
print(f\"Composite Score: {gene['composite_score']:.4f}\")
print(f\"Consensus Tier: {gene['consensus_tier']}\")
print(f\"Best Modality: {gene['best_modality']} (proba={gene['best_modality_proba']:.4f})\")
print(f\"Modality Specificity: {gene['modality_specificity']:.4f}\")
print(f\"PHAROS Mean: {gene['pharos_mean_proba']:.4f}\")
print(f\"Tier Mean: {gene['tier_mean_proba']:.4f}\")
print(f\"Runs in Top Decile: {int(gene['n_runs_top_decile'])} / 11\")
print(f\"Novel Everywhere: {'Yes' if gene['is_novel_everywhere'] else 'No'}\")
"
```

### Multi-target comparison
```python
genes = ['GREM1', 'IL11', 'TP53', 'EGFR']
subset = master[master.Gene_Name.isin(genes)][[
    'Gene_Name', 'composite_score', 'consensus_tier', 'best_modality',
    'best_modality_proba', 'modality_specificity', 'pharos_mean_proba',
    'tier_mean_proba', 'n_runs_top_decile', 'is_novel_everywhere'
]].sort_values('composite_score', ascending=False)
print(subset.to_string(index=False))
```

### Interpretation Guide

| composite_score | Interpretation |
|----------------|----------------|
| >= 0.8 | Strong druggability signal |
| 0.6 - 0.8 | Moderate signal, further validation needed |
| 0.3 - 0.6 | Weak signal |
| < 0.3 | Low/minimal signal |

| consensus_tier | Meaning |
|---------------|---------|
| High | >= 8 of 11 runs in top 10% (>= 90th percentile) |
| Moderate | >= 5 of 11 runs in top 25% (>= 75th percentile) |
| Low | Otherwise |

| modality_specificity | Meaning |
|---------------------|---------|
| > 0.15 | Strong modality recommendation |
| 0.05 - 0.15 | Moderate — best modality preferred but alternatives possible |
| < 0.05 | Multi-modal candidate (multiple approaches viable) |

For full biological definitions of PHAROS TDL, druggability tiers, modalities, and how to interpret per-run results, read `references/annotations.md` in this skill directory.

### Quick Biological Glossary (for interpreting results)

**PHAROS TDL:**
- Tclin = approved drug target (has >=1 approved drug with known MoA)
- Tchem = potent chemical tool target (bioactive <=30nM but no approved drug)
- Tbio = functionally annotated (OMIM/GO evidence but no chemical tools)
- Tdark = understudied (limited literature, no qualifying annotations)

**Druggability Tiers:**
- Tier 1 = clinical precedence (target of approved drug or clinical candidate)
- Tier 2 = discovery precedence (active compound but no clinical candidate)
- Tier 3A = predicted tractable by structure (druggable pocket/family)
- Tier 3B = predicted tractable by other evidence (druggable pathway/expression)

**Modalities:**
- Small Molecule = binds intracellular targets (enzymes, GPCRs, kinases, ion channels)
- Antibody = targets extracellular/secreted proteins (receptors, ligands, membrane proteins)
- PROTAC = degrades intracellular proteins via ubiquitin-proteasome (TFs, scaffolds, oncoproteins)

**PU Learning:**
- Known gene (known_gene=1) = in positive training set; established druggability evidence
- Novel gene (known_gene=0) = unlabelled; NOT assumed undruggable, just lacking prior evidence
- High-scoring novel = key prediction — model sees druggable features despite no prior evidence

### Recommendation Format

When a user asks "is X druggable?" or "druggability of X", present results as:

1. **One-line verdict:** e.g., "GREM1: Moderate druggability, Antibody-recommended"
2. **Key metrics table:** composite_score, consensus_tier, best_modality, modality_specificity
3. **Biological interpretation** (mandatory — translate every dimension into biological meaning):
   - **PHAROS interpretation:** Based on pharos_mean_proba and per-run known/novel status, explain what the PHAROS profile means. Example: "High score on tclin run = gene's features resemble approved drug targets with known MoA; novel in tclin = not currently an approved drug target but predicted to share characteristics with approved targets."
   - **Tier interpretation:** Based on tier_mean_proba and which tier runs score highest, explain tractability. Example: "High in tier_1 run = gene resembles targets with clinical precedence — approved drugs or clinical candidates exist for structurally similar targets."
   - **Modality interpretation:** Based on best_modality, explain the predicted therapeutic approach. Example: "Antibody = gene likely encodes extracellular/secreted protein amenable to monoclonal antibody targeting; Small Molecule = intracellular target with druggable pocket; PROTAC = degradable intracellular target that may lack traditional binding pocket."
   - **Known/Novel status:** For each dimension, explain whether the gene was in the positive training set (known = established evidence) or unlabelled (novel = model prediction without prior evidence). Highlight high-scoring novel genes as key discoveries.
4. **Per-run pattern analysis:** Use the pattern interpretation guidance from `references/annotations.md` Section 6 to identify and explain the gene's cross-run signature. Examples:
   - "High in Tclin but low in Tchem → resembles a clinical target but lacks chemical-tool features, suggesting biologic-only approach"
   - "High across all PHAROS runs → broad druggability features regardless of development level"
   - "High in Tier 1 + Tier 2 but low in Tier 3A/3B → established target class, not a novel tractability prediction"
   - "Novel everywhere + high composite → entirely model-predicted discovery candidate with no prior druggability evidence"
5. **Actionable guidance** (with biological rationale):
   - If composite >= 0.6: "Recommend further investigation as drug target"
   - If modality_specificity > 0.15: "Strong evidence for [modality] approach — [biological reason, e.g., 'extracellular localization supports antibody access']"
   - If modality_specificity < 0.05: "Multi-modal candidate — gene features overlap multiple target classes, consider screening across modalities"
   - If is_novel_everywhere: "Discovery candidate — entirely model-predicted, not in any training set; prioritize experimental validation of predicted target class"

---

## E: Single-Gene Druggability Profile

For deep profiling of a specific target, use `assess_gene()` from `drugnome_ai/post_analysis/druggability_assessment.py`. This outputs a comprehensive multi-line report covering all 11 runs grouped by dimension (PHAROS, Tier, Modality) with per-run probabilities, percentiles, known/novel status, and an interpretation section.

```bash
cd ~/claude_code/DrugnomeAI && conda run -n drugnome python3 -c "
import sys; sys.path.insert(0, '.')
from drugnome_ai.post_analysis.druggability_assessment import *
from pathlib import Path
master = build_master_table(Path('output03182026'))
master = add_derived_columns(master)
print(assess_gene('GENE_NAME_HERE', master))
"
```

Replace `GENE_NAME_HERE` with the target gene (case-insensitive).

---

## F: Cross-Run Assessment Generation

To regenerate the cross-run assessment (e.g., after new runs or with custom weights):

```bash
cd ~/claude_code/DrugnomeAI && conda run -n drugnome python3 \
  drugnome_ai/post_analysis/druggability_assessment.py \
  --output-dir output03182026 \
  --top-n 200

# With custom weights (must sum to 1.0)
conda run -n drugnome python3 \
  drugnome_ai/post_analysis/druggability_assessment.py \
  --output-dir output03182026 \
  --w-pharos 0.25 --w-tier 0.35 --w-modality 0.40

# Query a single gene directly
conda run -n drugnome python3 \
  drugnome_ai/post_analysis/druggability_assessment.py \
  --output-dir output03182026 --gene GREM1
```

**Default weights:** pharos=0.3, tier=0.3, modality=0.4

**Outputs generated:**
- `cross_run_assessment/druggability_master_table.csv` — Full 20,080-gene table
- `cross_run_assessment/top200_novel_targets.csv` — Top novel genes by composite score
- `cross_run_assessment/top_novel_small_molecule.csv` — Top novel small molecule targets
- `cross_run_assessment/top_novel_antibody.csv` — Top novel antibody targets
- `cross_run_assessment/top_novel_protac.csv` — Top novel PROTAC targets
- `cross_run_assessment/assessment_summary.txt` — Distribution statistics

---

## G: Data Source Updates

DrugnomeAI uses features from two primary source families:

| Source Family | Flag | Features |
|--------------|------|----------|
| PHAROS | `-d pharos` | Target development level annotations (Tclin, Tchem, Tbio, Tdark) |
| InterPro | `-d inter` | Protein domain/family features (`-x dom fam sup`) |

Additional feature sets from mantis-ml heritage:
- `-m` — Generic mantis-ml features (ExAC, Essential Mouse Genes, GnomAD, Genic Intolerance, GWAS, MGI)
- `-l` — Genic Intolerance Scores specifically

Data files live in `drugnome_ai/data/`. To refresh:
1. Download updated source files to appropriate subdirectories
2. Re-run with `-r pre` to recompile feature tables
3. Then `-r pu` and `-r post` for training and post-processing

---

## Parameter Reference

### CLI Flags (`drugnomeai`)

| Flag | Long | Description | Default |
|------|------|-------------|---------|
| `-o` | `--output-dir` | Output directory (required) | — |
| `-c` | `--config-file` | YAML config for disease-specific analysis | None |
| `-r` | `--run-tag` | Pipeline stage: `all`, `pre`, `boruta`, `pu`, `post`, `post_unsup`, `debug` | `all` |
| `-f` | `--fast` | Fast mode (4 classifiers only: ET, RF, SVC, GB) | off |
| `-s` | `--superv-models` | Specific classifiers: `et rf svc gb xgb dnn stack nb` | all 6 default |
| `-d` | `--data-source` | Data sources: `pharos`, `inter` | all |
| `-x` | `--inter-pro` | InterPro feature types: `dom`, `fam`, `sup` | all |
| `-t` | `--tier-tag` | Druggability tiers: `1`, `2`, `3A`, `3B` | all |
| `-p` | `--pharos-tag` | PHAROS categories: `tclin`, `tchem`, `tbio`, `tdark` | all |
| `-m` | `--mantis-ml` | Include mantis-ml generic features | off |
| `-l` | `--genic-intol` | Include Genic Intolerance Scores | off |
| `-k` | `--known-genes-file` | Custom known gene list (newline-separated) | None |
| `-n` | `--nthreads` | Parallel threads | 4 |
| `-i` | `--iterations` | Stochastic PU learning iterations | 10 |

Note: `-t`/`-p` and `-k` are mutually exclusive. `-t` and `-p` are mutually exclusive.

### Assessment Script Flags (`druggability_assessment.py`)

| Flag | Description | Default |
|------|-------------|---------|
| `--output-dir` | Directory containing the 11 run folders | `output03182026` |
| `--gene` | Query a single gene (prints profile) | None |
| `--top-n` | Number of top genes in ranked outputs | 200 |
| `--w-pharos` | PHAROS dimension weight | 0.3 |
| `--w-tier` | Tier dimension weight | 0.3 |
| `--w-modality` | Modality dimension weight | 0.4 |

---

## Output Format Reference

### Master Table Columns

**Per-run columns** (11 runs x 3 columns = 33):
- `{prefix}_proba` — Mean prediction probability across iterations
- `{prefix}_perc` — Percentile rank (0-100)
- `{prefix}_known` — 1 if gene was in training set, 0 if novel

**Prefixes:** `tclin`, `tchem`, `tclin_tchem`, `tier1`, `tier12`, `tier123A`, `tier123AB`, `tclin_tier1`, `small_mol`, `antibody`, `protac`

**Derived columns:**
- `pharos_mean_proba` — Mean across 3 PHAROS runs
- `tier_mean_proba` — Mean across 5 tier runs
- `best_modality` — Highest-scoring modality (Small Molecule, Antibody, or PROTAC)
- `best_modality_proba` — Probability of the best modality
- `modality_specificity` — Gap between best and second-best modality
- `composite_score` — Weighted sum: 0.3×pharos + 0.3×tier + 0.4×modality
- `consensus_tier` — High / Moderate / Low
- `n_runs_above_75perc` — Count of runs where gene is in top 25%
- `n_runs_top_decile` — Count of runs where gene is in top 10%
- `is_novel_everywhere` — 1 if gene is novel (unlabelled) in all 11 runs

### Directory Structure

```
output03182026/
├── run_all_iterations.log
├── pharos_tclin/          # Each of the 11 runs follows this structure:
│   ├── data/compiled_feature_tables/
│   ├── processed-feature-tables/
│   ├── supervised-learning/
│   │   ├── PU_{Classifier}.evaluation_metrics.tsv
│   │   ├── gene_predictions/
│   │   ├── gene_proba_predictions/
│   │   └── ranked-by-proba_predictions/
│   ├── Gene-Predictions/
│   │   └── AllClassifiers.Merged.drugnome_ai_predictions.csv
│   ├── Output-Figures/
│   └── unsupervised-learning/
├── ... (10 more run directories)
└── cross_run_assessment/
    ├── druggability_master_table.csv
    ├── top200_novel_targets.csv
    ├── top_novel_small_molecule.csv
    ├── top_novel_antibody.csv
    ├── top_novel_protac.csv
    └── assessment_summary.txt
```

---

## H: Scoring & Report Generation (Individual + Combo)

This section covers scoring arbitrary gene lists and optional combinations, with dynamic interpretation by Claude Code (or /codex-skill).

### Scoring Script

**Script:** `~/.claude/skills/drugnomeai/scripts/drugnome_combo_score.py`

The script produces clean numeric TSVs. Interpretation is NOT hardcoded — it is generated dynamically by the calling agent after TSV creation.

```bash
# Individual genes only (no combos)
cd ~/claude_code/DrugnomeAI && conda run -n drugnome python \
  ~/.claude/skills/drugnomeai/scripts/drugnome_combo_score.py \
  --genes "TYK2,JAK1,IL17A,IL23A,GREM1" \
  --output-dir .

# Individual genes + combos
cd ~/claude_code/DrugnomeAI && conda run -n drugnome python \
  ~/.claude/skills/drugnomeai/scripts/drugnome_combo_score.py \
  --genes "TYK2,JAK1,IL17A,IL17F,IL23A" \
  --combos "TYK2+JAK1" "IL17A+IL17F+IL23A" \
  --output-dir .

# With --report flag (signals that calling agent should generate interpretation)
cd ~/claude_code/DrugnomeAI && conda run -n drugnome python \
  ~/.claude/skills/drugnomeai/scripts/drugnome_combo_score.py \
  --genes "TYK2,JAK1" --combos "TYK2+JAK1" --report \
  --output-dir .
```

**Script flags:**

| Flag | Description | Required |
|------|-------------|----------|
| `--genes` | Comma-separated gene names | Yes |
| `--combos` | Space-separated combo strings (GENE1+GENE2+...) | No |
| `--master-table` | Path to master CSV (default: skill references table) | No |
| `--output-dir` | Output directory (default: current working directory) | No |
| `--report` | Flag placeholder — signals agent to generate report | No |

**Outputs:**
- `drugnome_individual_score.tsv` — always created
- `drugnome_combo_score.tsv` — only created when `--combos` is provided

**Append behavior:** If TSVs already exist, only missing genes/combos are computed and appended. Existing entries are preserved.

### Report Generation Workflow (--report flag)

When `--report` is requested (or user asks for "report", "interpretation", "scoring with analysis"):

**IMPORTANT:** If unsure whether the user wants a report, use AskUserQuestion to confirm before generating.

**Step 1: Run the scoring script** (produces TSVs with numeric scores)

**Step 2: Read the data** — read both TSVs plus the master table for full per-run details:
```bash
cd ~/claude_code/DrugnomeAI && conda run -n drugnome python3 -c "
import pandas as pd
from pathlib import Path
master = pd.read_csv(Path('~/.claude/skills/drugnomeai/references/druggability_master_table.csv').expanduser())
genes = ['GENE1', 'GENE2', ...]  # from user's list
subset = master[master.Gene_Name.isin(genes)]
# Print ALL columns for each gene for interpretation
for _, row in subset.iterrows():
    print(f'=== {row.Gene_Name} ===')
    for col in master.columns:
        print(f'  {col}: {row[col]}')
    print()
"
```

**Step 3: Dynamically interpret each gene** — For EACH gene, Claude Code analyzes:
- **Cross-run signature:** Which runs score high vs low? What pattern does this reveal?
- **PHAROS profile:** Known Tclin? Novel in Tchem? What does the gap between runs mean?
- **Tier pattern:** Uniformly high? Jump at tier123A? What does this say about tractability?
- **Modality profile:** How dominant is the best modality? Multi-modal or modality-specific?
- **Known/Novel status:** Known in which runs? Novel where? Is this a discovery candidate?
- **Actionable recommendation:** What therapeutic approach? What validation needed?

Use the interpretation patterns from `references/annotations.md` Section 6:
| Pattern | Interpretation |
|---------|----------------|
| High all PHAROS runs | Broad drug target features |
| High Tclin, low Tchem | Clinical target, biologic-only (no chemical tools) |
| High tiers uniformly | Strong structural tractability |
| Jump at tier123A/3B | Predicted tractable (not yet clinical) |
| High modality_specificity (>0.15) | Clear single modality recommendation |
| Low modality_specificity (<0.05) | Multi-modal candidate |
| Novel everywhere + high composite | Discovery candidate — model-predicted |

**Step 4: Dynamically interpret each combo** (if combos provided) — analyze:
- **Member strength distribution:** Which members are strong/weak anchors?
- **Modality concordance:** All same modality? Mixed? What does this mean for co-formulation?
- **Score heterogeneity:** Even combo or one member carrying it?
- **Biological rationale:** What pathways are being targeted? Synergistic or redundant?
- **Actionable recommendation:** Proceed? Validate weaker members? Replace a member?

**Step 5: Write outputs (INCREMENTAL — never overwrite existing interpretations)**

5a. **Add Interpretation column to TSVs:**
- Read the TSV into a DataFrame
- Check which rows already have a non-empty `Interpretation` (or `Combo_Interpretation`) value
- Only fill in interpretations for rows where the column is missing or NaN
- Write the full DataFrame back (preserves existing interpretations, adds new ones)

Example pattern:
```python
idf = pd.read_csv(individual_tsv, sep="\t")
# Only interpret genes that don't already have an interpretation
needs_interp = idf["Interpretation"].isna() | (idf["Interpretation"] == "")
for idx in idf[needs_interp].index:
    gene = idf.loc[idx, "Gene_Name"]
    idf.loc[idx, "Interpretation"] = NEW_INTERPRETATIONS[gene]  # from Claude Code
idf.to_csv(individual_tsv, sep="\t", index=False)
```

5b. **Append to Markdown report** `<output-dir>/drugnome_report.md`:
- If report file exists, read it and check which genes/combos already have sections
- Only append NEW sections for genes/combos not yet in the report
- If report file does not exist, create it with all sections
- Always regenerate the Summary Tables at the bottom (they should reflect ALL entries)

Example check:
```python
existing_report = REPORT_MD.read_text() if REPORT_MD.exists() else ""
for gene in new_genes:
    if f"### {gene} " not in existing_report:
        # Append this gene's section
```

### Important: Incremental Update Rules
| Scenario | TSV Behavior | Report Behavior |
|----------|-------------|-----------------|
| First run (no files exist) | Create with all entries + interpretations | Create full report |
| Re-run same genes/combos | Skip (already in TSV with interpretation) | Skip (already in report) |
| New genes/combos added | Append scores (scoring script) then fill interpretation for new rows only | Append new sections only, regenerate summary tables |
| Never | Overwrite existing interpretations | Remove existing report sections |

### When NOT to generate report
- User only asked for scores/TSVs without interpretation
- User explicitly says "no report" or "just scores"
- Default behavior (no --report flag) is scores only

### When to AskUserQuestion
- User's intent about report is ambiguous
- User asks for "druggability" of genes but doesn't specify report vs. quick lookup
- Example: "check druggability of X" — ask if they want a full report or just scores

---

## Troubleshooting

### GPU Issues
```bash
# Check GPU status
nvidia-smi

# Force CPU if GPU issues
export DRUGNOME_DNN_DEVICE=cpu

# Check PyTorch CUDA
conda run -n drugnome python3 -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}, Devices: {torch.cuda.device_count()}')"
```

### sklearn / Package Errors
```bash
# Verify environment
conda run -n drugnome python3 -c "import sklearn, xgboost, torch; print('All imports OK')"

# Reinstall if needed
conda activate drugnome && pip install -e .
```

### `set -e` Failures in Batch Script
The `run_all_iterations.sh` uses `set -e` — any non-zero exit kills the entire batch. If a single run fails:
1. Check the log: `grep "FAILED" output03182026/run_all_iterations.log`
2. Fix the issue for that specific run
3. Run just that one: `drugnomeai -o output03182026/<run_name> <flags> -s et rf gb svc xgb dnn stack nb -n 8 -i 10 -r all`

### Lock File Blocking
```bash
# If "Another instance is already running" error:
ls /tmp/drugnome_runs.lock
# Only remove if you're sure no other instance is running:
rm /tmp/drugnome_runs.lock
```

### Missing Merged Predictions
If `AllClassifiers.Merged.drugnome_ai_predictions.csv` is missing, run post-processing:
```bash
conda run -n drugnome drugnomeai -o output03182026/<run_name> -r post
```

### Master Table Not Found
Ensure all 11 runs completed, then regenerate:
```bash
cd ~/claude_code/DrugnomeAI && conda run -n drugnome python3 \
  drugnome_ai/post_analysis/druggability_assessment.py --output-dir output03182026
```
