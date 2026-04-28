---
name: scimilarity
description: >
  Complete toolkit for scimilarity (v0.4.1) single-cell RNA-seq foundation model — cell embedding,
  annotation, search, and gene interpretation. Use this skill whenever users mention scimilarity,
  cell type annotation with kNN, cell search across reference atlases, single-cell embedding models,
  gene attribution with integrated gradients, cosine similarity for single cells, or reference atlas
  search. Also trigger when users want to embed cells into a 128-dim latent space, annotate cell types
  using a 7M+ cell reference, search for similar cells across public datasets, interpret which genes
  drive cell identity, or work with the SCimilarity pretrained model. Covers model download and setup,
  data preparation (log normalization, gene alignment), cell type annotation via kNN, six cell search
  modes (nearest, centroid, cluster, exhaustive variants), gene interpretation via Integrated Gradients,
  ontology utilities, visualization (circle packing, heatmaps), training/fine-tuning with metric
  learning, and all utility functions. Even if the user just asks about "cell similarity" or
  "finding similar cells" in scRNA-seq data, this skill likely applies.
---

# SCimilarity Toolkit

**SCimilarity** is a cell atlas foundation model that embeds single-cell RNA-seq expression data into a 128-dimensional L2-normalized latent space for cell type annotation and cell search across 7M+ reference cells.

**Publication**: Heimberg et al. "A cell atlas foundation model for scalable search of similar human cells." *Nature* (2024). DOI: 10.1038/s41586-024-08411-y

**Package**: `scimilarity` v0.4.1 (installed at `~/.conda/envs/scgpt/lib/python3.10/site-packages/scimilarity/`)

## Class Hierarchy

```
CellEmbedding              # Base: loads encoder, computes embeddings
  └── CellSearchKNN        # Adds kNN index loading and neighbor search
        ├── CellAnnotation  # Cell type prediction via kNN voting
        └── CellQuery       # Large-scale cell search (6 search modes)

Interpreter                 # Gene attribution via Integrated Gradients (standalone)
```

## Quick Start

### 1. Embed Cells
```python
from scimilarity import CellEmbedding, align_dataset, lognorm_counts

ce = CellEmbedding(model_path="/path/to/model")
adata = align_dataset(adata, ce.gene_order)
adata = lognorm_counts(adata)  # requires adata.layers["counts"]
embeddings = ce.get_embeddings(adata.X)  # shape: [n_cells x 128]
```

### 2. Annotate Cell Types
```python
from scimilarity import CellAnnotation

ca = CellAnnotation(model_path="/path/to/model")
adata = ca.annotate_dataset(adata)  # adds obs["celltype_hint"], obsm["X_scimilarity"]
# Or step by step:
embeddings = ca.get_embeddings(align_dataset(adata, ca.gene_order).X)
predictions, nn_idxs, nn_dists, stats = ca.get_predictions_knn(embeddings, k=50)
```

### 3. Search Similar Cells
```python
from scimilarity import CellQuery

cq = CellQuery(model_path="/path/to/model")
embeddings = cq.get_embeddings(align_dataset(adata, cq.gene_order).X)
nn_idxs, nn_dists, metadata = cq.search_nearest(embeddings, k=10000)
```

### 4. Interpret Genes
```python
from scimilarity import Interpreter

interp = Interpreter(ca.model, ca.gene_order)
attrs = interp.get_attributions(anchor_cells_X, negative_cells_X)
ranked = interp.get_ranked_genes(attrs)
interp.plot_ranked_genes(ranked, n_plot=15)
```

## Model Setup

See [references/model_setup.md](references/model_setup.md) for download instructions and directory structure.

**S3 location**: `s3://tec-rnd-sci-dev-gi2/gi2-xz/models/scimilarity/model_v1.1/`
**Zenodo source**: Record 10685499 (`model_v1.1.tar.gz`, ~28GB)

### S3 Cache (Recommended)

scimilarity requires local file paths. Use the S3 cache helper to auto-sync from S3:

```python
import sys
sys.path.insert(0, "/home/sagemaker-user/.claude/skills/scimilarity/scripts")
from s3_cache import get_model_path

model_path = get_model_path(tier="annotation")  # syncs ~9GB on first use, instant after
# Tiers: "embedding" (~250MB), "annotation" (~9GB), "full" (~28GB)
```

Or from CLI:
```bash
python /home/sagemaker-user/.claude/skills/scimilarity/scripts/s3_cache.py --tier annotation
# Cache location: /tmp/scimilarity/model_v1.1/
```

All pipeline scripts support `--s3` flag to auto-cache from S3:
```bash
python scripts/annotate_cells.py --input data.h5ad --s3 --output annotated.h5ad
```

**Expected directory structure after extraction**:
```
model_v1.1/
├── encoder.ckpt                    # Encoder weights
├── gene_order.tsv                  # Gene vocabulary (~48K genes)
├── layer_sizes.json                # Network architecture
├── label_ints.csv                  # Label-to-integer mapping
├── annotation/
│   ├── labelled_kNN.bin            # Annotation kNN index
│   └── reference_labels.tsv        # Reference cell type labels
└── cellsearch/
    ├── full_kNN.bin                # Search kNN index (~large)
    ├── full_kNN_meta.csv           # Cell metadata
    ├── cell_metadata/              # TileDB metadata store
    └── cell_embedding/             # TileDB embedding store
```

## Data Preparation

All scimilarity workflows require properly prepared data:

```python
from scimilarity.utils import filter_cells, align_dataset, lognorm_counts
from scimilarity.utils import convert_id2symbol, consolidate_duplicate_symbols

# 1. Raw counts must be in layers["counts"]
adata.layers["counts"] = adata.X.copy()

# 2. Optional: convert Ensembl IDs to gene symbols
adata = convert_id2symbol(adata, mapping_table="ensembl_mapping.tsv")

# 3. Optional: merge duplicate gene symbols
adata = consolidate_duplicate_symbols(adata)

# 4. Optional: QC filter cells
adata = filter_cells(adata, min_genes=400, mito_percent=30.0)

# 5. Align gene space (min 5000 gene overlap required)
adata = align_dataset(adata, model.gene_order)

# 6. Log normalize (CPM 1e4 + log1p)
adata = lognorm_counts(adata)
```

**Critical requirements**:
- `layers["counts"]` must contain raw integer counts (not normalized)
- Gene symbols (HGNC) in `var.index` (not Ensembl IDs)
- Minimum 5,000 gene overlap with model's gene_order
- No NaN values in expression matrix
- No negative values in counts

## Cell Annotation

`CellAnnotation` predicts cell types via kNN voting against a labeled reference.

```python
ca = CellAnnotation(model_path="/path/to/model")

# Full pipeline (handles alignment internally)
adata = ca.annotate_dataset(adata)
# Outputs in adata.obs: celltype_hint, min_dist, celltype_hits, celltype_hint_stat,
#   celltype_hits_weighted, celltype_hint_weighted_stat
# Outputs in adata.obsm: X_scimilarity (128-dim embeddings)

# Manual pipeline with more control
predictions, nn_idxs, nn_dists, stats = ca.get_predictions_knn(
    embeddings, k=50, ef=100, weighting=False
)
```

**Prediction stats interpretation**:
- `vsAll`: fraction of k neighbors belonging to predicted class (higher = more confident)
- `vs2nd`: ratio of best vs 2nd-best class (closer to 1.0 = ambiguous)
- `min_dist`: distance to nearest neighbor (lower = more similar to reference)
- `hits`: JSON dict of class counts among k neighbors

**Blocklist/Safelist** to control which cell types can be predicted:
```python
ca.blocklist_celltypes(["T cell", "B cell"])     # exclude these types
ca.safelist_celltypes(["macrophage", "monocyte"]) # only predict these
ca.reset_knn()                                    # clear all filters
```

## Cell Search

`CellQuery` searches for similar cells across the reference atlas. Six search modes:

### Decision Tree: Which Search Mode?

```
Do you have a precomputed embedding?
├── YES → search_nearest() or search_exhaustive()
└── NO → Do you want to search for a specific cell population?
    ├── YES (marked cells) → search_centroid_nearest() or search_centroid_exhaustive()
    └── YES (cluster labels) → search_cluster_centroids_nearest() or search_cluster_centroids_exhaustive()

Fast approximate (kNN index) vs Complete (exhaustive)?
├── Fast: *_nearest() methods — uses hnswlib, returns top-k
└── Complete: *_exhaustive() methods — brute-force cdist, returns all within max_dist
```

### search_nearest — Fast kNN search
```python
nn_idxs, nn_dists, metadata = cq.search_nearest(embeddings, k=10000, max_dist=0.05)
```

### search_centroid_nearest — Centroid of marked cells + kNN
```python
adata.obs["query"] = (adata.obs["celltype"] == "macrophage").astype(int)
centroid_emb, nn_idxs, nn_dists, metadata, qc = cq.search_centroid_nearest(
    adata, centroid_key="query", k=10000, qc=True
)
# qc["query_coherence"] — higher means query cells are internally consistent
```

### search_cluster_centroids_nearest — Per-cluster kNN
```python
centroid_embs, cluster_idx, nn_idxs, nn_dists, metadata = \
    cq.search_cluster_centroids_nearest(adata, cluster_key="leiden", k=10000)
# nn_idxs is a dict keyed by cluster label
```

### search_exhaustive — Brute-force within max_dist
```python
nn_idxs, nn_dists, metadata = cq.search_exhaustive(
    embeddings, max_dist=0.03, metadata_filter={"tissue": "lung"}, buffer_size=100000
)
```

### search_centroid_exhaustive / search_cluster_centroids_exhaustive
Same as above but with centroid/cluster preprocessing.

### Working with search results
```python
# Aggregate metadata at sample level
sample_stats = cq.compile_sample_metadata(nn_idxs[0],
    levels=["study", "sample", "tissue", "disease"])
# Columns: study, sample, tissue, disease, cells, fraction, total

# Add per-sample cell indices
metadata = cq.annotate_cell_index(metadata)
```

## Gene Interpretation

`Interpreter` uses Integrated Gradients (via Captum) to identify genes driving the difference between anchor and negative cell populations.

```python
from scimilarity import Interpreter

interp = Interpreter(encoder=ca.model, gene_order=ca.gene_order)

# anchors and negatives must have same shape [n_cells x n_genes]
# Use aligned, log-normalized expression data
attrs = interp.get_attributions(anchor_X, negative_X)  # [n_cells x n_genes]

ranked_genes = interp.get_ranked_genes(attrs)
# Columns: gene, gene_idx, attribution, attribution_std, cells

interp.plot_ranked_genes(ranked_genes, n_plot=15, filename="top_genes.pdf")
```

High attribution genes: (1) expressed more in anchors than negatives AND (2) strongly affect embedding distance. This identifies genes driving cell identity differences, not just differentially expressed genes.

## Visualization

See `scimilarity.visualizations` for circle packing and heatmap functions:

```python
from scimilarity.visualizations import hits_circles, hits_heatmap

# Circle packing of search results
hits_circles(metadata, levels=["tissue", "disease"],
    label_column="tissue", value_column="cells", filename="circles.pdf")

# Heatmap of search results
hits_heatmap(sample_metadata, rows="tissue", columns="disease",
    values="cells", filename="heatmap.pdf")
```

## Ontology Utilities

Import and traverse biomedical ontologies for cell type hierarchy analysis:

```python
from scimilarity.ontologies import (
    import_cell_ontology, import_uberon_ontology,
    get_all_ancestors, get_all_descendants, find_most_viable_parent,
    ontology_similarity, all_pair_similarities
)

cl = import_cell_ontology()  # Cell Ontology graph
ancestors = get_all_ancestors(cl, "CL:0000084")  # T cell ancestors
similarity = ontology_similarity(cl, "CL:0000084", "CL:0000236")  # path distance
```

Available ontologies: Cell Ontology (`import_cell_ontology`), Uberon tissues (`import_uberon_ontology`), DOID diseases (`import_doid_ontology`), MONDO diseases (`import_mondo_ontology`).

## Training & Fine-Tuning

For training custom models, see the full API in [references/api_reference.md](references/api_reference.md) under "Training Components".

```python
from scimilarity.training_models import MetricLearning
from scimilarity.anndata_data_models import MetricLearningDataModule

model = MetricLearning(
    n_genes=n_genes, latent_dim=128, hidden_dim=[1024, 1024, 1024],
    margin=0.05, negative_selection="semihard", sample_across_studies=True
)
# Train with PyTorch Lightning, then save:
model.save_all("/path/to/output")
```

## Utility Functions

Key utilities in `scimilarity.utils` — see [references/api_reference.md](references/api_reference.md) for complete signatures:

| Function | Purpose |
|----------|---------|
| `lognorm_counts(data)` | Normalize total 1e4 + log1p |
| `align_dataset(data, gene_order)` | Align gene space to model vocabulary |
| `filter_cells(data, min_genes, mito_percent)` | QC filter cells |
| `get_centroid(counts)` | Compute mean expression centroid |
| `get_cluster_centroids(data, gene_order, key)` | Per-cluster centroids |
| `get_dist2centroid(centroid_emb, X)` | Cosine distance to centroid |
| `pseudobulk_anndata(adata, groupby)` | Aggregate to pseudobulk |
| `convert_id2symbol(adata, mapping)` | Ensembl ID to gene symbol |
| `consolidate_duplicate_symbols(adata)` | Merge duplicate genes |

## Performance Tips

- **GPU**: Pass `use_gpu=True` to `CellEmbedding`/`CellAnnotation`/`CellQuery` for ~10x speedup
- **buffer_size**: Increase for GPU (e.g., 50000), decrease for limited memory
- **ef parameter**: Higher ef = more accurate but slower kNN search (default: 100 for annotation, k for query)
- **Exhaustive search**: Use `metadata_filter` to reduce search space
- **Embeddings**: L2-normalized, 128-dim float32 — cosine distance range [0, 2]

## Troubleshooting

See [references/troubleshooting.md](references/troubleshooting.md) for common errors and fixes.

## Scripts

Helper scripts in this skill's `scripts/` directory:
- `download_model.py` — Download model from Zenodo with progress
- `annotate_cells.py` — End-to-end cell annotation pipeline
- `search_cells.py` — Cell search pipeline (all 6 modes)
- `embed_cells.py` — Standalone cell embedding extraction
- `interpret_genes.py` — Gene attribution analysis
