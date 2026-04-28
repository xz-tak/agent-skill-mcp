---
name: popv
description: >
  PopV SCDS2 consensus cell-type annotation pipeline for single-cell RNA-seq data.
  Use this skill whenever the user mentions popv, popV, cell type annotation with
  consensus voting, annotating cells with multiple classifiers, scVI+scanVI+CellTypist
  annotation, ontology-aware cell labeling, Tabula Sapiens/Muris models, or downloading
  popV reference models from HuggingFace. Also use when the user wants to create custom
  reference models from CellxGene Census for popV, generate model summaries, or run
  multi-classifier ensemble annotation on h5ad files. Triggers on phrases like
  "annotate cells", "popv annotation", "cell type prediction", "download popv models",
  "create popv reference model", "popv consensus vote", "run popv on my h5ad".
---

# PopV SCDS2 Cell-Type Annotation Pipeline

PopV uses 6-8 classifiers (scVI KNN, ScanVI, CellTypist, SVM, XGBoost, OnClass + BBKNN/Harmony KNN) with majority voting and Cell Ontology refinement to produce consensus cell-type annotations.

## Default Model Location (S3)

```
S3_MODEL_BASE = "s3://tec-rnd-sci-dev-gi2/gi2-xz/models/popv/data/huggingface_data"
```

Pre-trained Tabula Sapiens (human, 35 models) and Tabula Muris (mouse, 42 models) are stored here. Models are downloaded to `/tmp/popv_models/{MODEL}/` at runtime. Ontology files at `{S3_MODEL_BASE}/ontology/` are read directly from S3 where possible, downloaded to `/tmp/popv_ontology/` only when popv requires a local path.

**Source notebooks** at `/home/sagemaker-user/popv/`:
- `popv_query_annotation_scds2.ipynb` — Primary annotation workflow
- `popv_create_custom_model_scds2.ipynb` — Custom model training
- `popv_download_huggingface_collections_scds2.ipynb` — Model download
- `generate_model_summaries.py` — Metadata extraction

---

## Interactive Flow

When this skill triggers, Claude MUST follow this three-phase interactive flow before generating any code. Do NOT skip phases or assume parameter values.

### Phase 1: Workflow Selection

Use **AskUserQuestion** to determine which workflow the user needs:

```
Which popV workflow do you need?

1. **Download reference models** — Get pre-trained Tabula Muris/Sapiens from HuggingFace
2. **Create custom model** — Train a new reference model from CellxGene Census
3. **Query annotation** — Annotate cells in an h5ad file (most common)
4. **Model summaries / listing** — List available models or generate metadata summaries
```

Then proceed to the matching Phase 2 section below.

### Phase 2: Gather Parameters

#### Workflow 3 — Query Annotation (most common)

**Step 2a: Ask required inputs via AskUserQuestion**

Each of these MUST be explicitly confirmed by the user. Ask them together in one AskUserQuestion:

1. **Input h5ad** — File path or S3 URI to the query dataset (required)
2. **Organism** — `human` or `mouse` (required)
3. **Tissue or compartment** — e.g. `large intestine`, `liver`, `immune`, or `all` for pan-tissue (required)
4. **Output path** — Default: `{input_dir}/{shortname}_popv.h5ad`. Let user confirm or override
5. **Cell type / cluster column** — Which existing `obs` column to use for cluster-level comparison (e.g., `Cell_type`, `leiden`, `seurat_clusters`)

**Step 2a-confirm: Model selection confirmation (ALWAYS required)**

After the user answers Step 2a, use `select_model()` from `references/s3-io-patterns.md` to auto-select the model. Then ALWAYS present the selection via **AskUserQuestion** for confirmation:

```
Based on your data ({organism}, {tissue}), I recommend:
  Model: {MODEL}
  Collection: {COLLECTION}
  S3 path: s3://tec-rnd-sci-dev-gi2/gi2-xz/models/popv/data/huggingface_data/{COLLECTION}/{MODEL}/

1. **Proceed** with this model
2. **Use a different model** — specify which one (see references/available-models.md for full list)
```

The model directory defaults to `S3_MODEL_BASE` — no need to ask the user for a path. Models are downloaded to `/tmp/popv_models/{MODEL}/` at runtime.

**Step 2b: Inspect input h5ad — auto-detect remaining settings**

After the user confirms the 4 required inputs, run this via **Bash** to read h5ad metadata:

```python
python3 -c "
import scanpy as sc, json
adata = sc.read_h5ad('INPUT_PATH', backed='r')
info = {
    'n_obs': adata.n_obs, 'n_vars': adata.n_vars,
    'var_names_sample': list(adata.var_names[:5]),
    'gene_format': 'ensembl' if adata.var_names[0].startswith(('ENSG','ENSMUSG')) else 'symbol',
    'species_hint': 'mouse' if any(v.startswith('ENSMUSG') for v in adata.var_names[:100]) else 'human',
    'layers': list(adata.layers.keys()),
    'obs_columns': list(adata.obs.columns),
    'has_popv_cols': any(c.startswith('popv') for c in adata.obs.columns),
}
print(json.dumps(info, indent=2))
"
```

Replace `INPUT_PATH` with the user's actual input path. From the output, auto-detect:

| Setting | How to detect |
|---------|---------------|
| `counts_layer` | Pick the most likely raw counts layer: prefer `raw_counts` > `counts` > `spliced` > first available layer |
| `batch_list` | Find plausible batch columns from obs: `batch`, `donor_id`, `sample`, `protocol`, `library_id` |
| `gene_format` | `ensembl` if var_names start with ENSG/ENSMUSG, else `symbol` (determines if mygene mapping needed) |
| `species` | `mouse` if ENSMUSG prefixes found, else `human` (determines cross-species homolog mapping) |
| `SHORTNAME` | Derive from input filename (strip `.h5ad` extension) |
| `PREDICTION_MODEL` | Default `"inference"` |
| `N_SAMPLES_PER_LABEL` | Default `300` |
| `SOURCE` | Infer from model_dir path: `"huggingface_data"` if path contains it, else `"custom_scds2"` |
| `COLLECTION` | Infer from model_dir path structure (e.g., `tabula-sapiens`, `tabula-muris`, `Human`, `Mouse`) |

**Step 2c: Present validation summary and confirm**

Print a formatted summary of ALL settings — both user-confirmed and auto-detected:

```
## Configuration Summary

**User-confirmed:**
- Input: /path/to/sample.h5ad (45,231 cells x 33,694 genes)
- Model: popV_tabula_sapiens_Liver (from /data/models/huggingface_data/tabula-sapiens/)
- Output: /path/to/sample_popv.h5ad
- Cluster column: Cell_type

**Auto-detected:**
- Counts layer: raw_counts
- Batch columns: ["batch", "donor_id"]
- Gene format: symbol (mygene Ensembl mapping will be applied)
- Species: human
- Prediction model: inference
- N samples per label: 300
```

Then use **AskUserQuestion** with two choices:
1. **Proceed** — generate the annotation code
2. **I need to change something** — if selected, ask what to change and loop back

#### Workflow 1 — Download Reference Models

Use **AskUserQuestion** to ask:
1. **MODEL_DIR** — Where to save models (local path or S3 URI)
2. **Which collections** — Tabula Muris, Tabula Sapiens, or both (default: both)

Then confirm and proceed.

#### Workflow 2 — Create Custom Model

Use **AskUserQuestion** to ask all of:
1. **MODEL_DIR** — Where to save the model
2. **COLLECTION** — `"Human"`, `"Mouse"`, or `"NHP"`
3. **TISSUE** — Tissue group (e.g., `"Heart"`, `"Pancreas"`)
4. **MODEL** — Model name (e.g., `"popV_CellxGene_Human_Author_2022_Heart"`)
5. **REF_DATASET_ID** — CellxGene dataset UUID

Confirm summary and proceed. Defaults: `N_SAMPLES_PER_LABEL=300`, `HVG=4000`.

#### Workflow 4 — Model Summaries

Use **AskUserQuestion** to ask:
1. **MODEL_DIR** — Path to the model directory

Then list models or generate summaries as requested.

### Phase 3: Generate Code

After the user confirms in Phase 2, generate the notebook cells or script:
- Use **NotebookEdit** for notebook cells (default format, matching existing .ipynb style)
- Use **Write** for standalone .py scripts (only when user explicitly requests)
- Fill ALL parameters from the confirmed configuration — never leave placeholders
- Follow the code patterns in the reference files (`references/workflow-annotation.md`, etc.)

---

## Workflow Decision Tree

```
User request
├── "download models" / "get popv models" → Workflow 1: Download Reference Models
├── "create/train model" / "custom reference" → Workflow 2: Create Custom Models
├── "annotate" / "predict cell types" / "run popv" → Workflow 3: Query Annotation
└── "list models" / "model summary" / "what models" → Workflow 4: Model Summaries
```

---

## Workflow 1: Download Reference Models

Downloads pre-trained Tabula Muris (42 models) and Tabula Sapiens (35 models) from HuggingFace, plus Cell Ontology files.

**Parameters:** Gathered via Interactive Flow Phase 2 (Workflow 1) above. Read `references/workflow-download.md` for the full code pattern.

**Key code:**
```python
from huggingface_hub import get_collection, snapshot_download
from pathlib import Path

BASE_DIR = Path(MODEL_DIR) / "huggingface_data"
BASE_DIR.mkdir(parents=True, exist_ok=True)

# Collections to download
COLLECTIONS = {
    "tabula-muris": "popV/tabula-muris-6791cedaf0ecdeb1a8a4840e",
    "tabula-sapiens": "popV/tabula-sapiens-67627b2bd44ba09e9129589a",
}
# Ontology dataset
ONTOLOGY_REPO = "popV/ontology"
```

---

## Workflow 2: Create Custom Models

Train custom reference models from CellxGene Census datasets. Uses `prediction_mode="retrain"` to train all classifiers from scratch.

**Parameters:** Gathered via Interactive Flow Phase 2 (Workflow 2) above. Read `references/workflow-custom-model.md` for the full code pattern.

| Parameter | Description | Example |
|-----------|-------------|---------|
| `MODEL_DIR` | Where to save the model | `/data/popv_models/` |
| `COLLECTION` | Organism group | `"Human"`, `"Mouse"`, `"NHP"` |
| `TISSUE` | Tissue group | `"Heart"`, `"Pancreas"` |
| `MODEL` | Model name | `"popV_CellxGene_Human_Author_2022_Heart"` |
| `REF_DATASET_ID` | CellxGene dataset UUID | `"5500c673-1610-40a0-..."` |
| `N_SAMPLES_PER_LABEL` | Cells per label for balancing | `300` (default) |
| `HVG` | Highly variable genes | `4000` (default) |

---

## Workflow 3: Query Annotation (Primary)

Annotates a query h5ad dataset using pre-trained or custom popV models. This is the most common workflow.

**Parameters:** Gathered via Interactive Flow Phase 2 (Workflow 3) above. The 4 required inputs are collected via AskUserQuestion; remaining settings are auto-detected from the h5ad file.

| Parameter | Source | Default | Description |
|-----------|--------|---------|-------------|
| `INPUT_PATH` | **AskUser** | — | Path to query h5ad (file or S3 URI) |
| `organism` | **AskUser** | `"human"` | `"human"` or `"mouse"` |
| `tissue` | **AskUser** | — | Tissue name (e.g., `"large_intestine"`, `"liver"`) or `"all"` |
| `MODEL` | **Auto-select + confirm** | — | Auto-selected via `select_model()`, ALWAYS confirmed by user |
| `OUTPUT_PATH` | **AskUser** | `{input_dir}/{shortname}_popv.h5ad` | Output path (user confirms or overrides) |
| `cluster_group` | **AskUser** | — | Existing obs column for cluster comparison |
| `SHORTNAME` | Auto-detect | from filename | Derived from input filename |
| `SOURCE` | Fixed | `"huggingface_data"` | Models stored on S3 |
| `COLLECTION` | Auto-select | from organism | `"tabula-sapiens"` or `"tabula-muris"` |
| `counts_layer` | Auto-detect | `"raw_counts"` | Best-match raw counts layer |
| `batch_list` | Auto-detect | `["batch"]` | Plausible batch columns from obs |
| `gene_format` | Auto-detect | — | `symbol` or `ensembl` (drives mygene mapping) |
| `species` | Auto-detect | — | `human` or `mouse` (drives homolog mapping) |
| `PREDICTION_MODEL` | Auto-detect | `"inference"` | `"fast"` or `"inference"` |
| `N_SAMPLES_PER_LABEL` | Auto-detect | `300` | Cells per label for balancing |
| `TISSUE` | Auto-detect | `""` | Required only when `SOURCE == "custom_scds2"` |

### I/O Defaults

- **Input**: `{WORK_DIR}/{DATASET}/{SHORTNAME}.h5ad` (or S3 URI)
- **Output h5ad**: `{WORK_DIR}/{DATASET}/{SHORTNAME}_popv.h5ad` (or user-specified)
- **QC figures**: `{WORK_DIR}/{DATASET}/popV_qc/{SHORTNAME}/`
- If user specifies a custom output dir, both h5ad and QC figures go there

### Pipeline Steps (high level)

1. Read query h5ad, set X to raw counts layer, create batch_key
2. Map gene symbols to Ensembl IDs via `mygene` (if needed)
3. Map mouse genes to human homologs (if cross-species)
4. Load reference adata and ontology
5. Run `popv.preprocessing.Process_Query(...)` — harmonization + classifier training
6. Run `popv.annotation.annotate_data(adata)` — ensemble prediction
7. Filter to query cells (inference mode may change cell count)
8. Copy popv columns back to original adata, add ontology IDs
9. Compute cluster-level majority vote
10. Generate QC visualizations (UMAPs, agreement plots, score bars, marker genes)
11. Clean up obsoleted intermediate columns before saving:

```python
# --- Remove per-classifier raw predictions (redundant with final consensus) ---
# These are intermediate outputs consumed by majority vote. Only final consensus is needed.
drop_popv_cols = [c for c in adata.obs.columns
                  if any(c.startswith(f'popv_{clf}_prediction') for clf in
                         ['celltypist', 'knn_bbknn', 'knn_harmony', 'knn_on_scvi',
                          'onclass', 'scanvi', 'svm', 'xgboost'])
                  or c.startswith('popv_onclass_seen')
                  or c.startswith('popv_parent_ontology')
                  or c in ('popv_labels', 'popv_prediction_depth',
                           'popv_prediction_onclass_relative_depth')]
if drop_popv_cols:
    adata.obs.drop(columns=drop_popv_cols, inplace=True)
    print(f"Removed {len(drop_popv_cols)} intermediate popv columns")

# Kept: popv_prediction_ontology_name/id, popv_prediction_score,
#        popv_majority_vote_prediction_ontology_name/id, popv_majority_vote_score,
#        popv_model_split, cluster_popv, cluster_popv_id, cluster_popv_score
```

12. Save annotated h5ad with gzip compression

Read `references/workflow-annotation.md` for the complete code.

### S3 Support (Default)

Models default to `S3_MODEL_BASE = "s3://tec-rnd-sci-dev-gi2/gi2-xz/models/popv/data/huggingface_data"`. Read `references/s3-io-patterns.md` for all helpers. Key points:
- **Models**: downloaded to `/tmp/popv_models/{MODEL}/` (popv requires local file access)
- **Ontology for Process_Query**: downloaded to `/tmp/popv_ontology/` (popv needs local cl_obo_folder path)
- **Ontology cl.json for post-processing**: read directly from S3 via `read_json_from_s3()` (no download)
- **Input h5ad from S3**: download to `/tmp/`, read with scanpy
- **Output to S3**: upload annotated h5ad and QC plots after completion

---

## Workflow 4: Model Summaries & Listing

Generate metadata summaries for trained models or list available models. **Parameters:** Gathered via Interactive Flow Phase 2 (Workflow 4) above.

**To list available models:**
```python
import os
from pathlib import Path

model_base = Path(MODEL_DIR) / "huggingface_data"
for collection in sorted(model_base.iterdir()):
    if collection.is_dir() and collection.name != "ontology":
        print(f"\n{collection.name}:")
        for model in sorted(collection.iterdir()):
            if model.is_dir():
                print(f"  {model.name}")
```

**To generate model.txt summaries:**
```python
# Located at /home/sagemaker-user/popv/generate_model_summaries.py
from generate_model_summaries import summarize_models
summarize_models(MODEL_DIR)
```

Read `references/available-models.md` for the full inventory of pre-trained models.

---

## Batch Mode

Default is single-file annotation. For multiple files, guide the user with this loop pattern:

```python
# Each adata paired with its own model
jobs = [
    {"shortname": "sample_A", "model": "popV_tabula_sapiens_Liver", "collection": "tabula-sapiens"},
    {"shortname": "sample_B", "model": "popV_tabula_sapiens_Lung",  "collection": "tabula-sapiens"},
]

for job in jobs:
    SHORTNAME = job["shortname"]
    MODEL = job["model"]
    COLLECTION = job["collection"]
    # ... run full annotation pipeline for each ...
```

Each iteration is independent — different models can be used per sample.

---

## Output Format

Default to **notebook cells** (NotebookEdit) matching the existing .ipynb style. Generate **standalone .py scripts** only when the user explicitly asks. Both formats use the same code patterns from the reference files.

---

## Quick Reference: Key Functions

| Function | Purpose |
|----------|---------|
| `popv.preprocessing.Process_Query(...)` | Harmonize query with reference, train classifiers |
| `popv.annotation.annotate_data(adata)` | Run all classifiers, majority vote |
| `popv.visualization.make_agreement_plots(...)` | Classifier agreement heatmaps |
| `popv.visualization.prediction_score_bar_plot(...)` | Confidence distribution |
| `sc.pl.umap(adata, color=[...])` | UMAP visualizations |

## Troubleshooting

Read `references/troubleshooting.md` for common issues: gene ID mismatches, cell count changes in inference mode, OnClass seen != prediction, CUDA warnings, and more.
