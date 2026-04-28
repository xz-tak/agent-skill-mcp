# SCimilarity Model Setup

## Download from Zenodo

The pretrained model is available as a single archive from Zenodo record 10685499.

**DOI**: 10.5281/zenodo.10685499
**File**: `model_v1.1.tar.gz` (~28.2 GB)
**License**: CC BY-SA 4.0

### Automated Download

Use the download script from this skill:

```bash
conda activate scgpt
python /path/to/skills/scimilarity/scripts/download_model.py \
  --output-dir /path/to/destination
```

### Manual Download

```bash
MODEL_DIR=/home/sagemaker-user/claude_code/scimilarity
cd $MODEL_DIR

# Download (~28GB, use -c for resume support)
wget -c "https://zenodo.org/records/10685499/files/model_v1.1.tar.gz?download=1" \
  -O model_v1.1.tar.gz

# Extract
tar -xzf model_v1.1.tar.gz

# Verify
ls model_v1.1/
# Expected: encoder.ckpt, gene_order.tsv, layer_sizes.json, label_ints.csv,
#           annotation/, cellsearch/
```

### S3 Location (Persistent Storage)

The model is stored on S3 for persistent access across SageMaker sessions:

```
s3://tec-rnd-sci-dev-gi2/gi2-xz/models/scimilarity/model_v1.1/
```

### S3 Cache Helper (Recommended)

Since scimilarity requires local file paths, use the S3 cache helper:

```python
import sys
sys.path.insert(0, "/home/sagemaker-user/.claude/skills/scimilarity/scripts")
from s3_cache import get_model_path

# Auto-syncs from S3 to /tmp/scimilarity/model_v1.1/ on first use
model_path = get_model_path(tier="annotation")  # ~9GB, includes kNN index
model_path = get_model_path(tier="embedding")   # ~250MB, encoder only
model_path = get_model_path(tier="full")         # ~28GB, includes cell search

# Clear cache
from s3_cache import clear_cache
clear_cache()
```

Or from CLI:
```bash
python /home/sagemaker-user/.claude/skills/scimilarity/scripts/s3_cache.py --tier annotation
python /home/sagemaker-user/.claude/skills/scimilarity/scripts/s3_cache.py --clear
```

All pipeline scripts support `--s3` to auto-cache:
```bash
python annotate_cells.py --input data.h5ad --s3 --output result.h5ad
python embed_cells.py --input data.h5ad --s3 --output embedded.h5ad
```

## Directory Structure

After extraction, the model directory should contain:

```
model_v1.1/
├── encoder.ckpt                    # Encoder neural network weights (PyTorch state_dict)
├── gene_order.tsv                  # Gene vocabulary (one HGNC symbol per line, ~48K genes)
├── layer_sizes.json                # Network layer dimensions (infers architecture)
├── label_ints.csv                  # Cell type label to integer mapping
├── annotation/
│   ├── labelled_kNN.bin            # HNSWlib kNN index for cell type annotation
│   └── reference_labels.tsv        # Tab-separated: celltype_name\tstudy_id
└── cellsearch/
    ├── full_kNN.bin                # HNSWlib kNN index for cell search (large)
    ├── full_kNN_meta.csv           # Cell metadata (study, sample, tissue, disease, index)
    ├── cell_metadata/              # TileDB directory for cell metadata
    └── cell_embedding/             # TileDB directory for precomputed embeddings
```

## File Purposes

| File | Used By | Purpose |
|------|---------|---------|
| `encoder.ckpt` | All classes | Trained encoder weights for embedding generation |
| `gene_order.tsv` | All classes | Defines the gene vocabulary and ordering |
| `layer_sizes.json` | All classes | Infers network architecture (hidden dims, latent dim) |
| `label_ints.csv` | All classes | Maps cell type labels to integers |
| `annotation/labelled_kNN.bin` | CellAnnotation | kNN index over labeled reference cells |
| `annotation/reference_labels.tsv` | CellAnnotation | Cell type label for each reference cell |
| `cellsearch/full_kNN.bin` | CellQuery | kNN index over entire cell atlas |
| `cellsearch/full_kNN_meta.csv` | CellQuery | Metadata for each cell in the atlas |
| `cellsearch/cell_metadata/` | CellQuery | TileDB store for efficient metadata queries |
| `cellsearch/cell_embedding/` | CellQuery | TileDB store for precomputed embeddings |

## Usage Tiers

Depending on your use case, you may not need all files:

| Tier | Files Needed | Use Case | Size |
|------|-------------|----------|------|
| Embedding only | encoder.ckpt, gene_order.tsv, layer_sizes.json, label_ints.csv | Just compute embeddings | ~50MB |
| Annotation | Above + annotation/ | Cell type prediction | ~5GB |
| Full search | Above + cellsearch/ | Search across 7M+ cells | ~28GB |

## Verification

```python
from scimilarity import CellEmbedding, CellAnnotation

# Test embedding model loads
ce = CellEmbedding(model_path="/path/to/model_v1.1")
print(f"Genes: {ce.n_genes}, Latent dim: {ce.latent_dim}")
print(f"Labels: {len(ce.int2label)} cell types")

# Test annotation model loads
ca = CellAnnotation(model_path="/path/to/model_v1.1")
print(f"Reference labels: {len(ca.idx2label)}")
print(f"Available classes: {len(ca.classes)}")
```

## Neural Network Architecture

The encoder architecture (inferred from `layer_sizes.json`):

```
Input (n_genes) → InputDropout(0.4) → Linear → BatchNorm → PReLU
→ [Dropout(0.5) → Linear → BatchNorm → PReLU] × N hidden layers
→ Linear → L2-Normalize → Output (latent_dim=128)
```

Default hidden dimensions: `[1024, 1024]` (from model_v1.1)
Output: 128-dimensional L2-normalized vector (unit hypersphere)
Distance metric: Cosine distance = 1 - dot_product(emb1, emb2), range [0, 2]
