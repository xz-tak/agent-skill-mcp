---
name: scgpt
description: >
  Complete toolkit for scGPT (v0.2.4) — a generative pre-trained transformer for single-cell biology
  (Nature Methods 2024). Use this skill whenever users mention scGPT, single-cell transformers,
  gene expression prediction, in silico perturbation, perturbation prediction, cell type annotation
  with transformers, batch integration, multi-omic integration, GRN inference from transformers,
  cell embedding with scGPT, fine-tuning scGPT, masked language modeling for gene expression,
  or TransformerGenerator/TransformerModel from scGPT. Also trigger when users want to predict
  how genetic perturbations change cell state, do reverse phenotype perturbation analysis,
  generate cell embeddings from a pretrained transformer, or fine-tune a foundation model on
  scRNA-seq data. Covers model download, preprocessing, tokenization, all 4 task types (annotation,
  integration, perturbation, multiomic), the full training pipeline, and gene regulatory network
  inference. Even if the user just says "perturbation prediction" or "cell transformer" in a
  single-cell context, this skill likely applies.
---

# scGPT Toolkit

**scGPT** is a generative pre-trained transformer for single-cell biology that learns cell and gene representations from 33M+ human cells. It supports cell embedding, annotation, batch integration, perturbation prediction, multi-omic analysis, and GRN inference.

**Publication**: Cui et al. "scGPT: toward building a foundation model for single-cell multi-omics using generative AI." *Nature Methods* (2024). DOI: 10.1038/s41592-024-02201-0

**Package**: `scgpt` v0.2.4 (conda env: `scgpt`)

## Architecture Overview

```
TransformerModel               # Base: cell embedding, annotation, integration
  ├── GeneEncoder              # Gene token embeddings (nn.Embedding)
  ├── ContinuousValueEncoder   # Expression value → d_model MLP
  ├── TransformerEncoder       # Standard/Flash-attention transformer
  ├── ExprDecoder              # Expression prediction (MLM head)
  ├── ClsDecoder               # Cell type classification
  ├── MVCDecoder               # Masked value prediction from cell embedding
  └── AdversarialDiscriminator # Domain adversarial batch correction

TransformerGenerator           # Extends above for perturbation
  ├── pert_encoder             # nn.Embedding(3, d_model): 0=unpert, 1=pert, 2=pad
  └── AffineExprDecoder        # pred = coeff * input_values + bias

MultiOmicTransformerModel      # Extends TransformerModel with modality embeddings
```

## Quick Start

### 1. Load Model + Embed Cells
```python
from scgpt.tasks.cell_emb import embed_data

adata = embed_data(
    adata_or_file="data.h5ad",
    model_dir="/path/to/scgpt_model/",
    batch_size=64, device="cuda",
)
# adata.obsm["X_scGPT"] contains cell embeddings
```

### 2. Tokenize Data
```python
from scgpt.tokenizer import GeneVocab, tokenize_and_pad_batch

vocab = GeneVocab.from_file("/path/to/vocab.json")
tokenized = tokenize_and_pad_batch(
    data=adata.X, gene_ids=gene_ids,
    max_len=1200, vocab=vocab,
    pad_token="<pad>", pad_value=0, append_cls=True,
)
# Returns: {"genes": Tensor, "values": Tensor}
```

### 3. Perturbation Prediction
```python
from scgpt.model.generation_model import TransformerGenerator

model = TransformerGenerator(ntoken=len(vocab), d_model=512, nhead=8,
    d_hid=512, nlayers=12, nlayers_cls=3, n_cls=1, vocab=vocab)
# Load fine-tuned weights...
pred_expr = model.pred_perturb(batch_data, include_zero_gene="batch-wise",
    gene_ids=gene_ids, amp=True)
```

### 4. Fine-tune for Annotation
```python
from scgpt.trainer import prepare_data, prepare_dataloader, train
train_data, valid_data = prepare_data(tok_train, tok_valid, batch_labels_train,
    batch_labels_valid, config, epoch, celltype_labels_train, celltype_labels_valid)
loader = prepare_dataloader(train_data, batch_size=32)
train(model, loader, vocab, criterion_gep, criterion_dab, criterion_cls,
    scaler, optimizer, scheduler, device, config, logger, epoch)
```

## Model Setup

See [references/model_setup.md](references/model_setup.md) for download instructions.

**S3 location**: `s3://tec-rnd-sci-dev-gi2/gi2-xz/models/scgpt/`

### S3 Cache (Recommended)

scGPT requires local file paths. Use the S3 cache helper to auto-sync from S3:

```python
import sys; sys.path.insert(0, "/home/sagemaker-user/.claude/skills/scgpt/scripts")
from s3_cache import get_model_path

model_dir = get_model_path()  # syncs ~210MB on first use, instant after
# model_dir = "/tmp/scgpt"
```

Or from CLI:
```bash
python /home/sagemaker-user/.claude/skills/scgpt/scripts/s3_cache.py
# Cache location: /tmp/scgpt/
```

All pipeline scripts support `--s3` flag to auto-cache from S3:
```bash
python scripts/embed_cells.py --input data.h5ad --s3 --output embedded.h5ad
```

**Model files** (whole-human, ~208MB):
- `best_model.pt` — Pretrained weights
- `vocab.json` — Gene vocabulary (48,292 HGNC symbols)
- `args.json` — Model configuration/hyperparameters

**8 pretrained models available** (Google Drive): whole-human (recommended), continual-pretrained, brain, blood, heart, lung, kidney, pan-cancer.

## Preprocessing

```python
from scgpt.preprocess import Preprocessor

preprocessor = Preprocessor(
    use_key="X",                    # or specific layer
    filter_gene_by_counts=False,    # int or False
    filter_cell_by_counts=False,    # int or False
    normalize_total=1e4,            # CPM normalization (tp10k)
    log1p=True,                     # log1p transform
    subset_hvg=False,               # int for N HVGs, or False
    binning=None,                   # int for N bins, or None
)
preprocessor(adata, batch_key="batch")
# Results stored in adata.layers["X_normed"], ["X_log1p"], ["X_binned"]
```

**Normalization**: CPM 1e4 + log1p (same scale as scimilarity).

## Tokenization

The tokenizer converts expression matrices to gene-token sequences:

```python
from scgpt.tokenizer import (
    GeneVocab, tokenize_and_pad_batch, random_mask_value, get_default_gene_vocab
)

# Load vocabulary
vocab = GeneVocab.from_file("vocab.json")  # or get_default_gene_vocab()

# Map genes to vocab IDs
gene_ids = [vocab[g] for g in adata.var.index if g in vocab]

# Tokenize + pad
batch = tokenize_and_pad_batch(
    data=expression_matrix,  # (n_cells, n_genes)
    gene_ids=gene_ids,
    max_len=1200,            # max sequence length
    vocab=vocab,
    pad_token="<pad>", pad_value=0,
    append_cls=True,         # prepend <cls> token
    include_zero_gene=False, # only non-zero genes
)

# Apply random masking for training
masked_values = random_mask_value(batch["values"], mask_ratio=0.15, mask_value=-1)
```

**Special tokens**: `<pad>` (padding), `<cls>` (cell-level token, prepended).

## Four Task Types

scGPT supports 4 fine-tuning tasks, each using different loss combinations:

| Task | Model Class | Losses | Key Outputs |
|------|-------------|--------|-------------|
| `annotation` | TransformerModel | CLS | Cell type predictions |
| `integration` | TransformerModel | GEP + GEPC + ECS + DAB | Batch-corrected embeddings |
| `perturb` | TransformerGenerator | GEP + GEPC | Perturbed expression |
| `multiomic` | MultiOmicTransformerModel | GEP + GEPC | Joint embeddings |

**Loss acronyms**:
- **GEP**: Gene Expression Prediction (masked_mse_loss on masked genes)
- **GEPC**: Gene Expression Prediction from Cell embedding (MVC decoder)
- **CLS**: Classification (cross-entropy on cell type labels)
- **ECS**: Elastic Cell Similarity (cosine similarity regularization)
- **DAB**: Domain Adversarial Batch correction (adversarial discriminator)

## Cell Embedding

```python
from scgpt.tasks.cell_emb import get_batch_cell_embeddings, embed_data

# Quick path (handles everything)
adata = embed_data(adata, model_dir="/path/to/model/", batch_size=64, device="cuda")

# Manual path (more control)
embeddings = get_batch_cell_embeddings(
    adata, cell_embedding_mode="cls", model=model, vocab=vocab,
    max_length=1200, batch_size=64, gene_ids=gene_ids,
)
# Returns: (n_cells, d_model) numpy array, L2-normalized
```

**Embedding modes**: `"cls"` (CLS token), `"avg-pool"` (mean of all genes), `"w-pool"` (weighted by expression).

## Cell Type Annotation

Fine-tune with `config.task = "annotation"`, `config.CLS = True`:

```python
model = TransformerModel(ntoken=len(vocab), d_model=512, nhead=8,
    d_hid=512, nlayers=12, n_cls=num_celltypes, vocab=vocab,
    cell_emb_style="cls")

# Load pretrained, then fine-tune
load_pretrained(model, torch.load("best_model.pt"))
# ... training loop with criterion_cls = nn.CrossEntropyLoss()

# Predict
predictions = predict(model, test_loader, vocab, config, device)
```

## In Silico Perturbation

The perturbation system uses `TransformerGenerator` with `AffineExprDecoder`:

```python
model = TransformerGenerator(
    ntoken=len(vocab), d_model=512, nhead=8, d_hid=512,
    nlayers=12, nlayers_cls=3, n_cls=1, vocab=vocab,
    pert_pad_id=2,  # 3-state: 0=unpert, 1=pert, 2=pad
)

# After fine-tuning on perturbation data:
# batch_data.x[:, 0, :] = original expression values
# batch_data.x[:, 1, :] = perturbation flags (0/1/2)
pred_expr = model.pred_perturb(batch_data, include_zero_gene="batch-wise",
    gene_ids=gene_ids, amp=True)
# Returns: (batch_size, n_genes) predicted expression after perturbation
```

**How AffineExprDecoder works**: `pred = coeff * input_value + bias`
- For KO (value=0): `pred = bias` (learned basal expression)
- For KD (reduced value): `pred = reduced * coeff + bias` (proportional)

## Batch Integration

Fine-tune with `config.task = "integration"`:

```python
model = TransformerModel(
    ntoken=len(vocab), d_model=512, nhead=8, d_hid=512, nlayers=12,
    vocab=vocab, do_mvc=True, do_dab=True,
    use_batch_labels=True, num_batch_labels=n_batches,
    domain_spec_batchnorm="dsbn",  # Domain-specific batch norm
)
```

Uses GEP + GEPC + ECS + DAB losses for batch-corrected representations.

## GRN Inference

Extract gene-gene regulatory networks from transformer attention:

```python
from scgpt.tasks.grn import GeneEmbedding

gene_embs = GeneEmbedding(embeddings_dict)  # gene_name -> vector
similarities = gene_embs.compute_similarities("TP53")  # top similar genes
network = gene_embs.generate_network(threshold=0.5)     # NetworkX graph
metagenes = gene_embs.get_metagenes(gene_adata)          # Leiden clusters
```

## Training Pipeline

See [references/api_reference.md](references/api_reference.md) for complete trainer API.

```python
from scgpt.trainer import prepare_data, prepare_dataloader, train, evaluate

# Masking + batching
train_data, valid_data = prepare_data(
    tok_train, tok_valid, batch_labels_train, batch_labels_valid,
    config, epoch, celltype_labels_train, celltype_labels_valid
)

train_loader = prepare_dataloader(train_data, batch_size=32, shuffle=True)
valid_loader = prepare_dataloader(valid_data, batch_size=32)

# Train one epoch
train(model, train_loader, vocab, criterion_gep, criterion_dab,
    criterion_cls, scaler, optimizer, scheduler, device, config, logger, epoch)

# Evaluate
val_loss = evaluate(model, valid_loader, vocab, criterion_gep,
    criterion_dab, criterion_cls, device, config, epoch)
```

**Key config attributes**: `task`, `GEP`, `GEPC`, `CLS`, `ESC`, `DAR`, `mask_ratio`, `mask_value`, `dab_weight`, `explicit_zero_prob`

## Performance Tips

- **GPU**: Required for training; embedding works on CPU but slow
- **max_seq_len**: 1200 genes default; reduce for memory constraints
- **batch_size**: 32-64 typical; reduce for large models
- **AMP**: Always enabled via `torch.cuda.amp.GradScaler()`
- **Flash-attention**: Set `use_fast_transformer=True` for ~2x speedup
- **Gradient clipping**: 1.0 default

## Troubleshooting

See [references/troubleshooting.md](references/troubleshooting.md) for common errors.

## Scripts

- `scripts/download_model.py` — Download from HuggingFace or Google Drive
- `scripts/embed_cells.py` — Cell embedding extraction
- `scripts/annotate_cells.py` — Cell type annotation pipeline
- `scripts/perturbation.py` — In silico perturbation prediction
- `scripts/finetune.py` — Fine-tuning pipeline for all 4 tasks
