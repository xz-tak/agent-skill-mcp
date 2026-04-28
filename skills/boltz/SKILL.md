---
name: boltz
description: "Predict biomolecular structures using Boltz 2.2.1 deep learning model. Use this skill whenever the user mentions protein structure prediction, protein folding, Boltz, structure prediction, 3D structure, predict structure from sequence, fold a protein, AlphaFold-like prediction, biomolecular structure, predict protein complex, structure of [gene name], or wants to predict/model/fold any protein, DNA, RNA, or ligand structure. Also triggers on gene names combined with 'structure', 'fold', or 'predict'. Supports single proteins, multi-chain complexes, protein-ligand, protein-DNA/RNA systems, and binding affinity estimation. Always use this skill even if the user just gives a gene name and says 'predict its structure' or 'what does it look like'."
user-invocable: true
---

# Boltz Protein Structure Prediction

Predict biomolecular structures from gene names, amino acid sequences, or custom YAML inputs using the Boltz 2.2.1 deep learning model. The pipeline checks RCSB PDB for existing experimental structures (with full entity metadata: sequence length, description, method, molecular weight), fetches canonical sequences from UniProt Swiss-Prot (with TrEMBL fallback), and runs GPU-accelerated predictions.

## Environment

- **Conda env**: `boltz` (Python 3.12, Boltz 2.2.1)
- **GPUs**: Auto-detected NVIDIA GPUs (A10G 23GB each). Pipeline assigns 1 job per GPU for parallel batch processing.
- **Skill scripts**: `~/.claude/skills/boltz/scripts/` — the canonical source of truth for all pipeline code
- **Default output**: `{workdir}/protein_structure/` for single-protein predictions. For other prediction types, output to a category-appropriate subdirectory under the workdir:
  - Single protein → `{workdir}/protein_structure/`
  - Protein-protein complex → `{workdir}/protein_complex/`
  - Protein-ligand complex → `{workdir}/protein_ligand/`
  - Protein-DNA/RNA complex → `{workdir}/protein_nucleic_acid/`
  - Batch predictions → `{workdir}/protein_structure/` (default) or user-specified
  - Within each category: PDB downloads to `{identifier}_pdb/`, Boltz predictions to `boltz_results_{identifier}/`, results log to `results_summary.json` — all in the same flat directory
  - If user requests a custom output directory, use that instead

## Step 0: Bootstrap Pipeline Scripts

**Before running any predictions**, ensure the pipeline scripts exist in the working directory. The skill bundles scripts at `~/.claude/skills/boltz/scripts/`. Copy them to the output directory if not already present:

```bash
SKILL_SCRIPTS="$HOME/.claude/skills/boltz/scripts"
WORKDIR="{workdir}/protein_structure"

# Create workdir and copy scripts (preserves existing outputs)
mkdir -p "$WORKDIR/modules"
cp -n "$SKILL_SCRIPTS/predict.py" "$WORKDIR/"
cp -n "$SKILL_SCRIPTS/config.yaml" "$WORKDIR/"
cp -n "$SKILL_SCRIPTS/modules/__init__.py" "$WORKDIR/modules/"
cp -n "$SKILL_SCRIPTS/modules/check_pdb.py" "$WORKDIR/modules/"
cp -n "$SKILL_SCRIPTS/modules/fetch_sequence.py" "$WORKDIR/modules/"
cp -n "$SKILL_SCRIPTS/modules/run_boltz.py" "$WORKDIR/modules/"
cp -n "$SKILL_SCRIPTS/modules/utils.py" "$WORKDIR/modules/"
```

Use `cp -n` (no-clobber) to avoid overwriting user customizations. If the user reports bugs or you need to force-update, use `cp` without `-n`.

### Skill directory structure

```
~/.claude/skills/boltz/
├── SKILL.md                              # This file (workflow guide)
├── references/
│   ├── pipeline-reference.md             # CLI reference, config, output format, module API, troubleshooting
│   └── boltz-yaml-format.md              # Boltz YAML schema for complex predictions
└── scripts/                              # Pipeline source code (canonical)
    ├── predict.py                        # Main entry point
    ├── config.yaml                       # Default configuration (single source of truth)
    └── modules/
        ├── __init__.py
        ├── check_pdb.py                  # PDB search + download + entity metadata
        ├── fetch_sequence.py              # UniProt canonical sequence lookup (Swiss-Prot first)
        ├── run_boltz.py                  # Boltz subprocess execution (no hardcoded defaults)
        └── utils.py                      # Config loading, GPU detection, shared utilities
```

**Reference docs** (read these when you need details beyond the workflow):
- [references/pipeline-reference.md](references/pipeline-reference.md) — Full CLI reference, config.yaml format, results_summary.json schema, module API signatures, troubleshooting guide, known issues
- [references/boltz-yaml-format.md](references/boltz-yaml-format.md) — Boltz YAML input schema for complex predictions (protein-protein, protein-ligand, protein-DNA/RNA, constraints, affinity)

## Step 1: Collect Parameters via AskUserQuestion

Before running anything, always present the user with choices. Parse their prompt for context (gene names, sequences, etc.) and pre-fill where possible, but always confirm.

### Question Set 1: Input Type

```
Question: "What would you like to predict?"
Options:
  - "Single protein by gene name" — e.g., BRCA1, EGFR, TP53
  - "Single protein from amino acid sequence" — paste or provide a raw AA sequence
  - "Multi-entity complex (protein-protein, protein-ligand, protein-DNA/RNA)" — Boltz multimer/complex mode
  - "Batch prediction from gene list" — provide a file with one gene per line
```

If the user's prompt already clearly specifies the input (e.g., "predict BRCA1 structure"), pre-select the matching option and ask for confirmation rather than forcing them through the full wizard.

### Question Set 2: Prediction Settings

```
Question: "Which prediction settings do you want?"
Options:
  - "Quick defaults (Recommended)" — boltz2 model, 3 recycling steps, 1 sample, MSA server, mmcif output
  - "High accuracy" — 5 recycling steps, 3 diffusion samples, potentials enabled
  - "Custom" — specify individual parameters
```

### Question Set 3: PDB Behavior (for gene name inputs)

```
Question: "How should we handle existing PDB structures?"
Options:
  - "Download PDB + run Boltz prediction (Recommended)" — get both experimental and predicted structures
  - "Download PDB only, skip Boltz if found" — save GPU time when experimental structure exists
  - "Skip PDB check, predict only" — go straight to Boltz prediction
```

### Question Set 4: Organism (for gene name inputs)

```
Question: "Which organism?"
Options:
  - "Human (Homo sapiens)" — default
  - "Mouse (Mus musculus)"
  - "Other" — specify organism name and NCBI taxonomy ID
```

### For Multi-Entity Complex Mode

Use a step-by-step wizard:
1. Ask for first entity type (protein/DNA/RNA/ligand) and its sequence/identifier
2. Ask "Add another entity?" — repeat until user says no
3. For ligands: accept CCD code or SMILES string
4. Ask about constraints (bond, pocket, contact) — optional
5. Ask about affinity prediction — if ligand present, offer it

Build the Boltz YAML automatically from collected entities.

**PDB complex search**: Before running Boltz, search RCSB PDB for existing experimental structures of the complex. This is important because PDB often has solved co-crystal structures. Search strategy:
- For protein-protein: search PDB for entries containing both gene names using the RCSB Search API with a group query (AND logic) on `rcsb_entity_source_organism.rcsb_gene_name.value` for each gene
- For protein-ligand: search PDB for entries with both the protein gene name and the ligand CCD code or similar compound
- For protein-DNA/RNA: search for entries with the protein gene name and polymer type "DNA"/"RNA"
- Download the best-resolution hit to `{complex_name}_pdb/` in the output directory
- Report the complex PDB structure alongside Boltz predictions so user can compare

```python
# Example: search PDB for complex with two proteins
import requests
query = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                           "operator": "exact_match", "value": "BRCA1"}},
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                           "operator": "exact_match", "value": "BARD1"}}
        ]
    },
    "return_type": "entry",
    "request_options": {
        "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        "paginate": {"start": 0, "rows": 5}
    }
}
resp = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query)
```

**Power user shortcut**: If user provides a path to an existing Boltz YAML file, validate it and run directly — skip the wizard.

## Step 2: Run the Pipeline

All commands run from the bootstrapped `{workdir}/protein_structure/` directory (set up in Step 0).

### For Single Protein (Gene Name)

```bash
cd {workdir}/protein_structure && CONDA_PREFIX="/home/sagemaker-user/.conda/envs/boltz" conda run -n boltz python predict.py {GENE_NAME} {flags}
```

Flags based on user choices:
- `--skip-boltz-if-pdb` — if user chose PDB-only mode
- `--skip-pdb` — if user chose predict-only mode
- `--recycling-steps N --diffusion-samples N --use-potentials` — for custom/high-accuracy
- `--output-format pdb|mmcif`
- `--organism "Name" --organism-id NNNN` — for non-human
- `--max-parallel N` — for batch mode

### For Raw Sequence

```bash
cd {workdir}/protein_structure && CONDA_PREFIX="/home/sagemaker-user/.conda/envs/boltz" conda run -n boltz python predict.py --raw-sequence {SEQUENCE}
```

### For Batch Mode

```bash
cd {workdir}/protein_structure && CONDA_PREFIX="/home/sagemaker-user/.conda/envs/boltz" conda run -n boltz python predict.py --input-file {FILE_PATH}
```

### For Multi-Entity Complex (Custom YAML)

When the wizard has collected multiple entities, build a Boltz YAML file:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: MVTPEG...
  - protein:
      id: B
      sequence: QLEDSE...
  - ligand:
      id: C
      smiles: 'CCO'
# Optional sections:
constraints:
  - pocket:
      binder: C
      contacts: [100, 200]
properties:
  affinity:
    binder: C
```

Write it to `{workdir}/protein_structure/{complex_name}.yaml`, then run:

```bash
cd {workdir}/protein_structure && CONDA_PREFIX="/home/sagemaker-user/.conda/envs/boltz" conda run -n boltz boltz predict {yaml_path} --out_dir . --devices 1 --use_msa_server --recycling_steps 3 --diffusion_samples 1 --output_format mmcif --model boltz2
```

Set the CUDA and library path environment:
```bash
export CUDA_VISIBLE_DEVICES=0
export LD_LIBRARY_PATH="/home/sagemaker-user/.conda/envs/boltz/lib/python3.12/site-packages/nvidia/cu13/lib:${LD_LIBRARY_PATH}"
```

## Step 3: Report Results

After prediction completes, present a structured summary:

```
## Prediction Results

| Gene | PDB ID | Resolution | PDB Seq Length | Method | Entity Description | Boltz Status | Confidence | Output Path |
|------|--------|-----------|----------------|--------|--------------------|-------------|------------|-------------|
| BRCA1 | 8RS8 | 1.31 Å | 1863 | X-RAY DIFFRACTION | BRCA1-associated... | completed | 0.742 | boltz_results_BRCA1/... |
```

Read `results_summary.json` for structured data. Report:
- PDB structure found (ID + resolution + entity sequence length + method + title) or "No experimental structure"
- PDB entity metadata: description, gene name, organism, molecular weight (kDa)
- Boltz prediction status + confidence score
- File paths for all outputs
- Any errors encountered

## Step 4: Offer Downstream Analysis

After reporting results, ask the user what they'd like to do next:

```
Question: "What would you like to do with the prediction results?"
Options:
  - "Visualize the structure" — generate a py3Dmol notebook cell or export for viewer
  - "Compare with PDB experimental structure" — RMSD calculation, structural alignment
  - "Analyze binding sites" — identify pockets, active sites
  - "Run pathway/interaction analysis" — use pathwaydb-query, interactdb-query, biobridge skills
  - "Predict binding affinity" — re-run with affinity properties (if ligand involved)
  - "I'm done for now"
```

For visualization, create a Python cell using py3Dmol:
```python
import py3Dmol
view = py3Dmol.view()
with open("{output_path}/{job_name}_model_0.cif") as f:
    view.addModel(f.read(), "cif")
view.setStyle({"cartoon": {"color": "spectrum"}})
view.zoomTo()
view.show()
```

For pathway/interaction analysis, invoke the relevant skills (biobridge, interactdb-query, pathwaydb-query) with the gene name.

## Step 5: Save User Preferences

After a successful run, if the user chose non-default settings, offer to save preferences:

1. **Memory**: Save organism preference, common Boltz params, frequently queried genes to the Claude memory system for future sessions
2. **Config**: Update `{workdir}/protein_structure/config.yaml` with their preferred defaults

Only save if the user explicitly agrees or asks to remember settings.

## Key Technical Details

- **Boltz YAML version must be 1**: The schema parser enforces `version: 1`
- **Output path convention**: Boltz creates `{out_dir}/boltz_results_{input_stem}/predictions/{input_stem}/`
- **Confidence JSON**: Located at `predictions/{name}/confidence_{name}_model_0.json`
- **CUDA fix**: Must set `LD_LIBRARY_PATH` to include the cu13 lib path for GPU inference to work
- **Timeout**: Default 3600s per job. Long sequences (>1000 aa) may need more. Timeouts are captured and logged — batch continues.
- **GPU parallelism**: Pipeline auto-detects GPUs with <10% utilization and assigns 1 job per GPU
- **UniProt canonical lookup**: The pipeline queries `reviewed:true` (Swiss-Prot) first to get canonical full-length sequences, falling back to TrEMBL if no reviewed entry exists. This prevents picking up short fragment isoforms.
- **PDB entity metadata**: For each PDB hit, the pipeline fetches entry-level metadata (title, method, deposition date) and entity-level metadata (sequence length, description, gene name, organism, molecular weight) via the RCSB Data API. The `pdb_entity_sequence_length` field in `results_summary.json` reflects the actual chain length in the crystal structure.
- **Default: Boltz runs for all genes** including those with PDB structures. Use `--skip-boltz-if-pdb` to skip Boltz for genes with experimental structures. Without this flag, you get both PDB download AND full-length Boltz prediction.
- **Auto re-run low confidence**: If a prediction scores below the threshold (default 0.5, configurable via `--rerun-threshold` or `boltz.rerun_confidence_threshold` in config.yaml), the pipeline automatically re-runs with boosted params (recycling=10, diffusion_samples=5, use_potentials=True) and keeps the better result.

## Reference Documentation

All detailed reference material lives inside the skill — no standalone files outside the skill directory:

- **[references/pipeline-reference.md](references/pipeline-reference.md)** — CLI flags, config.yaml format, results_summary.json schema, module API signatures, troubleshooting, known issues and past fixes
- **[references/boltz-yaml-format.md](references/boltz-yaml-format.md)** — Boltz YAML input schema for complex predictions (protein-protein, protein-ligand, protein-DNA/RNA, constraints, templates, modifications, affinity)

Consult `pipeline-reference.md` for error handling, troubleshooting, and config setup details.
