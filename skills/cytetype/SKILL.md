---
name: cytetype
description: Automated cell type annotation for single-cell RNA-seq data using CyteType's multi-agent AI architecture with consensus validation. This skill should be used when users request cell type annotation of h5ad files, automated single-cell cluster labeling, CyteType analysis, or evidence-based cell annotation. Features GPU-aware compute (multi-GPU rsc / single-GPU rsc / scanpy CPU fallback), 3-agent consensus with GPT-5.2 reviewer harmonization, iterative cache-busted re-runs for convergence, and confidence-scored annotations. Use this skill whenever the user mentions cytetype, CyteType, cell type annotation with LLM, annotating cells with AI, cluster labeling, or cell annotation with confidence scores — even if they don't explicitly name CyteType.
---

# CyteType Consensus Cell Type Annotation

## Overview

CyteType annotates single-cell RNA-seq clusters using a multi-agent LLM architecture. This skill implements a **3-agent consensus + reviewer harmonization** strategy with GPU-aware compute for clustering.

**Pipeline:**
```
Input h5ad → Parameter Collection → GPU Detection → Cluster QC/Compute
  → 3 Independent CyteType Agents → GPT-5.2 Reviewer → Convergence Loop
  → Save Results → Marker Gene Scoring (CELLxGENE Cohen's d + bootstrap)
```

**Key outputs per cluster:**
- `cytetype_cluster` — harmonized annotation with gene markers
- `cytetype_ontologyTerm` — standard Cell Ontology term
- `cytetype_ontologyTermID` — CL ID (CL:XXXXXXX)
- `cytetype_cellState` — cell state/activation
- `cytetype_confidence` — numeric 0.0-1.0

## When to Use This Skill

Invoke when users request:
- "Annotate cell types in my h5ad file"
- "Run CyteType on my single-cell data"
- "Label clusters in my scRNA-seq data"
- "Annotate leiden clusters"
- "Cell type annotation with confidence scores"
- "CyteType consensus annotation"

## Interactive Flow

Follow this flow strictly. Use **AskUserQuestion** at each decision point.

### Phase 1: Input & Parameter Collection

**Step 1a: Get input h5ad path**

Ask the user for the h5ad file path. Then inspect it:

```python
import scanpy as sc, json
adata = sc.read_h5ad(INPUT_PATH, backed='r')
info = {
    'n_obs': adata.n_obs, 'n_vars': adata.n_vars,
    'obs_columns': list(adata.obs.columns),
    'uns_keys': list(adata.uns.keys()),
    'obsm_keys': list(adata.obsm.keys()),
    'layers': list(adata.layers.keys()),
}
```

**Step 1b: Ask key parameters via AskUserQuestion**

Collect all of these — present auto-detected defaults where possible:

1. **group_key** — Which obs column has cluster labels? (e.g., `leiden`, `louvain`, `cluster`). Auto-detect categorical columns as suggestions.
2. **rank_key** — Is `rank_genes_groups` in `adata.uns`? If not, offer to compute it.
3. **coordinates_key** — Which obsm key for visualization? Default `X_pca_harmony`, auto-detect available keys.
4. **study_context** — Biological context (organism, tissue, disease). Pre-fill from obs metadata columns (organism, tissue, disease_ontology, etc.), but always ask for confirmation.
5. **Compute option** — (A) Compute new leiden clustering + gene ranking, or (B) Use existing group_key/rank_key.

### Phase 2: Cluster Quality Check

**Step 3.1 — "Use Existing" path (option B from Step 1b)**

```python
cluster_sizes = adata.obs[group_key].value_counts()
min_size = cluster_sizes.min()
min_cluster = cluster_sizes.idxmin()
n_clusters = len(cluster_sizes)
```

If min cluster size <= 100 cells, use **AskUserQuestion**:
- "Cluster {min_cluster} has only {min_size} cells. Small clusters produce unreliable annotations."
  1. **Proceed anyway** — annotate all clusters including small ones
  2. **Recompute clusters** — run leiden to find a resolution where all clusters > 100 cells
  3. **Remove small clusters** — exclude clusters with < 100 cells

If min > 100: proceed directly to Phase 4.

**Step 3.2 — "Compute / Recompute" path (option A, or recompute from 3.1)**

Read `references/gpu_compute.md` for all code patterns used below.

**GPU detection and backend selection:**

1. Run `nvidia-smi` to detect GPUs with VRAM usage **strictly < 10%**
2. Select backend:
   - 2+ idle GPUs + `rapids_singlecell` + `dask_cuda` available → **multi-GPU**
   - 1 idle GPU + `rapids_singlecell` available → **single-GPU**
   - No idle GPU or rsc unavailable → **scanpy CPU**
3. If CPU fallback and dataset > 500K cells: **AskUserQuestion** to warn about slowness and offer to abort

**Multi-GPU setup** (if selected):
```python
from dask_cuda import LocalCUDACluster
from dask.distributed import Client
import rapids_singlecell as rsc
import scipy.sparse as sp

cluster = LocalCUDACluster(
    CUDA_VISIBLE_DEVICES=",".join(idle_gpu_ids),
    threads_per_worker=10, protocol="tcp",
    rmm_pool_size="10GB", rmm_maximum_pool_size="20GB",  # per-worker, ~87% of GPU VRAM
    rmm_allocator_external_lib_list="cupy",
)
client = Client(cluster)

# IMPORTANT: For large datasets where X > single-GPU VRAM, swap out X before
# anndata_to_GPU. Leiden only needs obsp connectivities, not X.
X_backup = adata.X
layers_backup = dict(adata.layers)
adata.X = sp.csr_matrix(adata.shape, dtype='float32')  # empty placeholder
adata.layers.clear()
rsc.get.anndata_to_GPU(adata)  # only transfers small obsp (~700 MB)
```

**Single-GPU setup** (if selected):
```python
import rapids_singlecell as rsc
# Same X-swap pattern if X > GPU VRAM
rsc.get.anndata_to_GPU(adata)
```

**Auto-compute neighbors if missing:**

Check `adata.obsp` for `connectivities`. If absent:
1. Look for `X_pca_harmony` in obsm → use as `use_rep`
2. Else look for `X_pca` → use as `use_rep`
3. Else **AskUserQuestion** which obsm key to use
4. Run `rsc.pp.neighbors()` (GPU backends) or `sc.pp.neighbors()` (CPU)

**Check for existing cytetype_leiden column:**

If `cytetype_leiden` already exists in `adata.obs`, use **AskUserQuestion**:
- "Column 'cytetype_leiden' already exists from a previous run."
  1. **Overwrite** — replace with new clustering
  2. **Use new name** — create `cytetype_leiden_2`

**Leiden resolution iteration:**

Start at resolution=1.0, decrease by 0.1 each round, floor at 0.2.

Each round:
1. Compute leiden → `cytetype_leiden`
   - **Multi-GPU**: `rsc.tl.leiden(adata, resolution=X, random_state=42, key_added='cytetype_leiden', use_dask=True)`
   - **Single GPU**: `rsc.tl.leiden(adata, resolution=X, random_state=42, key_added='cytetype_leiden')`
   - **Scanpy**: `sc.tl.leiden(adata, resolution=X, random_state=42, key_added='cytetype_leiden')`
2. Report summary stats: resolution, n_clusters, min/median/max cluster size
3. If min cluster > 100: **AskUserQuestion** to approve this resolution
4. If min cluster <= 100 and resolution > 0.2: decrease by 0.1, repeat
5. If resolution hits floor (0.2): **AskUserQuestion** — accept current or abort

Once approved, **restore X** (from backup) and compute `rank_genes_groups`:
- First restore X: `adata.X = X_backup` (and layers), then `rsc.get.anndata_to_CPU(adata)` and close dask cluster
- rank_genes_groups needs the full X matrix. For large datasets where X > single-GPU VRAM:
  - Try **single-GPU rsc** (`method='wilcoxon'`) — if X fits on one GPU
  - If CUDA OOM, auto-fallback to **scanpy CPU** (`method='wilcoxon'`)
- For datasets where X fits in distributed GPU memory, use **multi-GPU dask** with `method='wilcoxon_binned'`
- See `references/gpu_compute.md` Pattern B for the full fallback code

**Cleanup GPU resources after clustering:**
- Multi-GPU: `rsc.get.anndata_to_CPU(adata); client.close(); cluster.close()`
- Single GPU: `rsc.get.anndata_to_CPU(adata)`

Then go back to Step 1b to confirm parameters with `cytetype_leiden` as `group_key`.

### Phase 3: Compute rank_genes_groups (if needed)

If `rank_key` not in `adata.uns` and not already computed in Step 3.2:

**For multi-GPU (default when X > single-GPU VRAM):**
1. Convert h5ad to zarr: `adata.write_zarr(zarr_path)` — one-time, preserves full object
2. Load zarr with dask chunking: `read_elem_as_dask(X, (50_000, n_vars))` — X stays lazy
3. `anndata_to_GPU` + `persist()` + `rechunk()` — distributes X across all GPUs
4. `rsc.tl.rank_genes_groups(method='wilcoxon_binned')` — dask-compatible method

Read `references/gpu_compute.md` Sections 4-5 and 9 for complete code patterns.

```python
# Multi-GPU: zarr → dask → distribute → wilcoxon_binned
adata = load_zarr_for_multi_gpu(zarr_path, sparse_chunk_size=50_000)
transfer_zarr_to_multi_gpu(adata, n_gpus=8)
rsc.tl.rank_genes_groups(adata, groupby=group_key, method='wilcoxon_binned', pts=True)
```

**For single-GPU (if X fits in VRAM):**
```python
rsc.tl.rank_genes_groups(adata, groupby=group_key, method='wilcoxon', pts=True)
```

**For CPU fallback (no GPU, rsc unavailable, or GPU OOM on very large datasets):**
```python
# t-test_overestim_var is O(n) — completes in minutes even on 3M cells x 38K genes.
# Statistically near-identical to wilcoxon for n >> 1000 cells per cluster.
# Use when: >2M cells x >30K genes exhausts GPU VRAM even with zarr+dask+float32.
sc.tl.rank_genes_groups(adata, groupby=group_key, method='t-test_overestim_var', pts=True)
```

### Phase 4: CyteType Consensus Annotation

This is the core annotation pipeline. Generate a Python script and run it.

**Default LLM config** (can be modified by user):
```python
llm_configs = [{
    "provider": "openai",
    "name": "gpt-5.2",
    "apiKey": os.environ["OPENAI_API_KEY"],
    "modelSettings": {"reasoning_effort": "high"},
}]
```

Read `references/consensus_workflow.md` for the complete code pattern.

The pipeline has 4 steps:

#### Step 4.1: CyteType Initialization (once, ~25 min for large datasets)

```python
from cytetype import CyteType
annotator = CyteType(
    adata, group_key=GROUP_KEY, rank_key=RANK_KEY,
    n_top_genes=200, coordinates_key=COORDINATES_KEY,
    pcent_batch_size=5000, max_cells_per_group=1000,
    max_metadata_categories=500,
    vars_h5_path=f"{OUTPUT_DIR}/vars.h5",
    obs_duckdb_path=f"{OUTPUT_DIR}/obs.duckdb",
)
```

#### Step 4.2: 3 Independent CyteType Annotator Agents (~2h each, sequential)

Run 3 times **sequentially** with **cache-busting** (unique run_id in study_context):
```python
from uuid import uuid4
for run_idx in range(1, 4):
    run_id = uuid4().hex[:8]
    annotator.run(
        study_context=f"{STUDY_CONTEXT} [run_id: {run_id}]",
        llm_configs=llm_configs,
        results_prefix=f"round1_run{run_idx}",
        n_parallel_clusters=5,
        timeout_seconds=14400,
        override_existing_results=True,
    )
```

**Why cache-bust**: CyteType's API caches identical payloads. Without unique run_ids, subsequent runs return cached results (~66s) instead of fresh annotations (~2h). The `[run_id: ...]` suffix forces independent computation.

#### Step 4.3: GPT-5.2 Reviewer Agent (~5s per cluster)

For each cluster, send ALL annotations + top 200 marker genes + PopV reference (if present in adata.obs) to a GPT-5.2 reviewer that produces a harmonized annotation with numeric confidence.

**PopV reference**: Check `adata.obs` for columns matching `popv_*` or `popV_*`. If found, include the majority annotation for each cluster. If absent, skip silently.

The reviewer outputs JSON:
```json
{
  "annotation": "JCHAIN-high IGHA1+ plasma cell",
  "ontologyTerm": "IgA plasma cell",
  "ontologyTermID": "CL:0000987",
  "cellState": "IgA-secreting",
  "confidence": 0.95,
  "reasoning": "All 3 runs agree this is an IgA plasma cell...",
  "agreement_level": "semantic_agreement"
}
```

Accept at **>= 0.9 confidence**.

Read `references/reviewer_prompt.md` for the full system prompt and JSON schema.

#### Step 4.4: Cache-Busted Re-runs for Unresolved Clusters

For clusters with confidence < 0.9:
1. Run 3 more CyteType annotations with new cache-busted run_ids
2. Reviewer now has 6+ annotations (more evidence → higher confidence)
3. Repeat up to **max 10 batches** (up to 33 annotations per cluster)
4. Force-accept best confidence after max batches

### Phase 5: Report & Approval

**Step 5a: Report annotation status**

Present to user via **AskUserQuestion**:
```
## CyteType Annotation Results

Total clusters: 61
Passed (>= 0.9 confidence): 58
Failed (< 0.9 confidence): 3

Failed clusters:
  Cluster 27 (26,886 cells): conf=0.78 — mixed identity (enteroendocrine/cycling)
  Cluster 44 (818 cells): conf=0.78 — stress-response T/NK
  Cluster 60 (8 cells): conf=0.78 — donor-specific artifact

Options:
1. Accept all (including low-confidence)
2. Rerun Step 4 on failed clusters only
3. Flag failed clusters and proceed
```

**Step 5b: If user chooses rerun**, go back to Step 4.4 with only the failed clusters.

### Phase 6: Save Results

Save all outputs to `./cytetype/` (or user-specified directory).

**Step 6a: Remove obsoleted/intermediate obs columns and uns keys**

Before saving, drop all intermediate CyteType columns and uns artifacts:

```python
# --- obs: drop intermediate annotation columns ---
drop_cols = [c for c in adata.obs.columns
             if (c.startswith('round') and '_run' in c)
             or (c.startswith('batch') and '_run' in c)
             or c == 'cluster_original']
if drop_cols:
    adata.obs.drop(columns=drop_cols, inplace=True)
    logging.info(f"Removed {len(drop_cols)} intermediate obs columns")

# Verify only cytetype_ prefixed columns remain from this pipeline
cytetype_cols = [c for c in adata.obs.columns if 'cytetype' in c.lower()]
assert all(c.startswith('cytetype_') for c in cytetype_cols), \
    f"Non-standard CyteType columns found: {[c for c in cytetype_cols if not c.startswith('cytetype_')]}"

# --- uns: drop intermediate batch/round job details and __cytetype temp keys ---
drop_uns = [k for k in adata.uns.keys()
            if k.startswith('batch') or k.startswith('round') or k.startswith('__cytetype')]
for k in drop_uns:
    del adata.uns[k]
if drop_uns:
    logging.info(f"Removed {len(drop_uns)} intermediate uns keys: {drop_uns[:5]}...")

# --- var: drop __cytetype_gene_symbols* temp columns ---
# CyteType creates these during initialization for internal gene symbol mapping.
# Multiple suffixed copies (_1, _2, ...) accumulate across re-runs.
drop_var = [c for c in adata.var.columns if c.startswith('__cytetype')]
if drop_var:
    adata.var.drop(columns=drop_var, inplace=True)
    logging.info(f"Removed {len(drop_var)} intermediate var columns: {drop_var}")
```

**Step 6b: Save outputs**

```python
OUTPUT_DIR = Path("./cytetype")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Save annotated adata
adata.write(f"{OUTPUT_DIR}/integration_cytetype.h5ad", compression="gzip")

# Save harmonized annotations table (CSV + JSON)
# harmonized_annotations.csv: cluster_id, n_cells, annotation, ontologyTerm,
#   ontologyTermID, cellState, confidence, agreement_level, reasoning

# Save config and run log
# config.json — all parameters, LLM config, timing, backend used
# run.log — full execution log

# Cleanup temp artifacts
annotator.cleanup()
```

**Output columns added to adata.obs** (all with `cytetype_` prefix):
| Column | Description |
|--------|-------------|
| `cytetype_cluster` | Harmonized annotation with gene markers |
| `cytetype_ontologyTerm` | Standard Cell Ontology term |
| `cytetype_ontologyTermID` | CL ID (e.g., CL:0000987) |
| `cytetype_cellState` | Cell state/activation |
| `cytetype_confidence` | Numeric confidence 0.0-1.0 |

**Output files:**
| File | Description |
|------|-------------|
| `integration_cytetype.h5ad` | Full annotated AnnData (with `uns['cytetype_marker']`) |
| `harmonized_annotations.csv` | Per-cluster annotation table |
| `harmonized_annotations.json` | Full harmonization detail |
| `config.json` | Run configuration + timing |
| `run.log` | CyteType execution log |
| `cytetype_markers.csv` | Top markers per cell type (CELLxGENE Cohen's d) |
| `cytetype_markers.pkl` | Dict of marker DataFrames keyed by cell type |

### Phase 7: Marker Gene Scoring

After annotation and save, compute CELLxGENE-style marker scores (pairwise Cohen's d + bootstrap).

**Step 7a: Ask marker parameters via AskUserQuestion**

Present defaults — user can accept or customize:

```
## Marker Gene Scoring

CyteType annotation is complete. Computing marker genes using the CELLxGENE algorithm
(pairwise Cohen's d with bootstrapped 10th percentile).

Questions:
1. Top N marker genes per cell type? [200]
2. Bootstrap resamples (B)? [1000]
3. Export format? [CSV + pkl]
```

AskUserQuestion options:
1. **Top N genes**: 200 (recommended) / 500 / All genes
2. **Bootstrap B**: 1000 (recommended) / 500 / 100
3. **Export format**: CSV + pkl (recommended) / CSV only / pkl only

**Step 7b: Run marker scoring**

Use the `scripts/compute_markers.py` module:

```python
import sys
sys.path.insert(0, "SKILL_SCRIPTS_DIR")  # ~/.claude/skills/cytetype/scripts/
from compute_markers import compute_cellxgene_markers

markers = compute_cellxgene_markers(
    adata,
    output_dir=OUTPUT_DIR,      # same as Phase 6 output dir
    top_n=200,                  # from user choice
    bootstrap_b=1000,           # from user choice
    percentile=10,              # fixed: CELLxGENE standard
    seed=42,
    update_h5ad=True,           # patches h5ad via h5py (not adata.write_h5ad)
    h5ad_path=f"{OUTPUT_DIR}/integration_cytetype.h5ad",
)
```

**CRITICAL**: The h5ad update uses `h5py` to patch only `uns['cytetype_marker']` in-place. **Never** use `adata.write_h5ad()` for this — anndata's read/write round-trip promotes float32→float64 and strips gzip compression, bloating a 23GB file to 139GB.

**Algorithm**: For each cell type T and each gene:
1. Compute Cohen's d (pooled SD) vs each of the other cell types → 51 effect sizes
2. Bootstrap B=1000 resamples, take P10 of each resample, average → marker score
3. Specificity = fraction of pairwise d > 0 (simple fraction)
4. If pooled SD = 0, set d = 0

**Step 7c: Report top markers**

Show the user a summary of top markers for a few representative cell types (e.g., highest confidence, largest cluster). This serves as a sanity check.

**Output files** (saved to output_dir):
| File | Description |
|------|-------------|
| `cytetype_markers.csv` | Concatenated top markers for all cell types |
| `cytetype_markers.pkl` | Dict of DataFrames keyed by cell type name |

**Output DataFrame schema** (per cell type):
| Column | Description |
|--------|-------------|
| `cytetype_leiden` | Cluster ID |
| `cytetype_cluster` | Full cell type annotation name |
| `Symbol` | Gene name |
| `Effect Size` | Bootstrapped P10 of pairwise Cohen's d |
| `Specificity` | Fraction of pairwise comparisons where d > 0 |
| `Mean Expression` | Mean log-normalized expression in target |
| `% of Cells` | Percentage of cells expressing (> 0) |

Also updates `adata.uns['cytetype_marker']` in the CyteType output h5ad.

**IMPORTANT — h5ad storage format**: `uns['cytetype_marker']` stores a **single concatenated DataFrame** (not a dict-of-DataFrames). This is required because CyteType cell type names contain `/` characters (e.g., `ACTA2/TAGLN/MCAM/...`) which HDF5 interprets as path separators, corrupting dict keys and massively bloating file size. To get per-type markers from the h5ad:

```python
# From h5ad (concatenated DataFrame)
markers_df = adata.uns['cytetype_marker']
plasma = markers_df[markers_df['cytetype_cluster'].str.contains('plasma', case=False)]

# From pkl (dict of DataFrames — more convenient)
import pickle
with open('cytetype/cytetype_markers.pkl', 'rb') as f:
    markers = pickle.load(f)
plasma = markers['JCHAIN-high IGHA1+ ...']
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| LLM model | `gpt-5.2` | OpenAI model for CyteType and reviewer |
| reasoning_effort | `high` | `high` avoids timeouts; `xhigh` too slow |
| n_top_genes | 200 | Marker genes per cluster (CyteType API and reviewer) |
| n_parallel_clusters | 5 | API parallelism per run |
| timeout_seconds | 14400 | 4h per run |
| reviewer_threshold | 0.9 | Confidence to accept |
| max_batches | 10 | Cache-busted re-run batches |
| coordinates_key | `X_pca_harmony` | Harmony-corrected PCA preferred |
| leiden_start_resolution | 1.0 | Starting resolution for leiden iteration |
| leiden_step | -0.1 | Decrease per round |
| leiden_floor | 0.2 | Minimum resolution before forced ask |
| min_cluster_size | 100 | Minimum cells per cluster for reliable annotation |
| gpu_vram_threshold | 10% | Max VRAM usage to consider a GPU idle |
| marker_top_n | 200 | Top marker genes per cell type (Phase 7) |
| marker_bootstrap_b | 1000 | Bootstrap resamples for marker scoring |
| marker_percentile | 10 | Percentile of pairwise Cohen's d (P10) |
| marker_seed | 42 | Random seed for bootstrap reproducibility |

## Reference Files

- `references/gpu_compute.md` — GPU detection, 3-tier backend selection, leiden iteration, rank_genes_groups patterns
- `references/consensus_workflow.md` — Complete code patterns for the 3-agent consensus + reviewer pipeline
- `references/api_reference.md` — CyteType Python API documentation
- `references/reviewer_prompt.md` — GPT-5.2 reviewer system prompt and JSON schema
- `scripts/compute_markers.py` — CELLxGENE marker scoring module (Phase 7). Importable as `from compute_markers import compute_cellxgene_markers` or CLI via `python compute_markers.py --input h5ad`

## Troubleshooting

| Issue | Solution |
|-------|----------|
| API timeout (7200s) | Increase `timeout_seconds` to 14400; use `reasoning_effort="high"` not `xhigh` |
| Cached results (~66s per run) | Ensure cache-busting: unique `run_id` in `study_context` per run |
| `max_tokens` unsupported | Use `max_completion_tokens` for GPT-5.2 reviewer calls |
| `openai` module missing | `pip install openai` in the conda env |
| Low confidence on small clusters | Clusters < 100 cells are inherently ambiguous; flag as `potential_artifact` |
| OPENAI_API_KEY not set | Load from `~/.env` via `python-dotenv` |
| Ontology terms inconsistent | CyteType CL terms vary across runs; reviewer harmonization resolves this |
| rsc ImportError (numpy/cupy) | Fix: `pip install numpy>=2.0 cupy-cuda12x`. Or scanpy fallback activates automatically |
| CUDA OOM on leiden | Automatic fallback to scanpy CPU for very large datasets |
| rsc not installed | Skill auto-detects and falls back to `scanpy.tl`. No manual intervention needed |
| No GPUs with <10% VRAM | All GPUs busy. Wait or proceed with scanpy CPU (slow for >500K cells) |
| `wilcoxon` error with Dask | Multi-GPU path must use `method='wilcoxon_binned'`, not `wilcoxon` |
| nvidia-smi not found | No NVIDIA driver. Falls back to scanpy CPU automatically |
| Multi-GPU dask hangs | Check `protocol="tcp"` (not "ucx"). Ensure `client.close(); cluster.close()` on cleanup |
| Existing cytetype_leiden column | Skill warns and asks before overwriting. Choose overwrite or new key name |
| Marker scoring OOM | Bootstrap on (38844, 51) matrix is ~6MB. If OOM, reduce `top_n` or run on smaller gene set |
| Marker scoring slow | B=1000 on 52 types × 38.8K genes takes ~30 min. Use B=100 for quick exploratory runs |
| No `pts` in rank_genes_groups | Re-run `sc.tl.rank_genes_groups(..., pts=True)` before marker scoring |
| h5ad bloats after marker save | Two causes: (1) Cell type names with `/` corrupt HDF5 dict keys — script stores as single DataFrame to avoid. (2) `adata.write_h5ad()` promotes float32→float64 and strips gzip — script uses h5py direct patch instead. If bloated, use `h5py` to delete `uns/cytetype_marker` and re-patch |
