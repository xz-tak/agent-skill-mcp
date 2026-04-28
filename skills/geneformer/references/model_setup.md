# Geneformer Model Setup

Guide for downloading, configuring, and verifying Geneformer V2 models.

---

## 1. Model Variants Overview

| Property | V2-104M | V2-316M | V2-104M_CLcancer |
|----------|---------|---------|------------------|
| **Parameters** | 104M | 316M | 104M |
| **Transformer Layers** | 12 | 18 | 12 |
| **Attention Heads** | 12 | 18 | 12 |
| **Hidden Dimension** | 768 | 1152 | 768 |
| **Max Input Size** | 4096 | 4096 | 4096 |
| **Vocab Size** | 20275 | 20275 | 20275 |
| **Training Data** | ~104M cells | ~104M cells | ~104M + 14M cancer cells |
| **Approximate Size** | ~399MB | ~1.2GB | ~399MB |
| **Use Case** | General purpose (default) | Higher capacity tasks | Cancer-specific applications |

All V2 models use the gc104M gene dictionary set and share the same tokenizer vocabulary.

---

## 2. S3 Storage

Models are stored in the team S3 bucket:

```
s3://tec-rnd-sci-dev-gi2/gi2-xz/models/geneformer/
```

### Model Tiers

| Tier | S3 Path | Notes |
|------|---------|-------|
| `V2-104M` | `.../geneformer/Geneformer-V2-104M/` | Default. Recommended starting point. |
| `V2-316M` | `.../geneformer/Geneformer-V2-316M/` | Larger model. Use when 104M underfits. |
| `V2-104M-CLcancer` | `.../geneformer/Geneformer-V2-104M_CLcancer/` | Continued pretraining on cancer data. |

### Manual Download

```bash
# Download default model
aws s3 sync s3://tec-rnd-sci-dev-gi2/gi2-xz/models/geneformer/Geneformer-V2-104M/ \
    /tmp/geneformer/Geneformer-V2-104M/

# Download 316M model
aws s3 sync s3://tec-rnd-sci-dev-gi2/gi2-xz/models/geneformer/Geneformer-V2-316M/ \
    /tmp/geneformer/Geneformer-V2-316M/

# Download cancer model
aws s3 sync s3://tec-rnd-sci-dev-gi2/gi2-xz/models/geneformer/Geneformer-V2-104M_CLcancer/ \
    /tmp/geneformer/Geneformer-V2-104M_CLcancer/
```

---

## 3. S3 Cache (Recommended)

Use the S3 cache helper to automatically download and cache models locally.

### Python API

```python
from s3_cache import get_model_path

# Default: V2-104M
model_path = get_model_path(tier="V2-104M")
# Returns: /tmp/geneformer/Geneformer-V2-104M/

# 316M model
model_path = get_model_path(tier="V2-316M")
# Returns: /tmp/geneformer/Geneformer-V2-316M/

# Cancer model
model_path = get_model_path(tier="V2-104M-CLcancer")
# Returns: /tmp/geneformer/Geneformer-V2-104M_CLcancer/
```

The function downloads from S3 on first call, then serves from the local cache
at `/tmp/geneformer/` on subsequent calls.

### CLI

```bash
# Cache V2-104M (default)
python scripts/s3_cache.py --tier V2-104M

# Cache V2-316M
python scripts/s3_cache.py --tier V2-316M

# Cache cancer model
python scripts/s3_cache.py --tier V2-104M-CLcancer
```

---

## 4. Expected Directory Structure

Each model directory must contain these files:

```
Geneformer-V2-104M/
├── config.json              # Model architecture configuration
├── model.safetensors        # Model weights (safetensors format)
├── generation_config.json   # Generation configuration
└── training_args.bin        # Original training arguments
```

The `config.json` defines the model architecture (number of layers, hidden size,
attention heads, vocabulary size, max position embeddings). The model loads via
`transformers.AutoModel.from_pretrained(model_path)`.

---

## 5. Gene Dictionaries

Gene dictionaries are **bundled with the geneformer pip package** and are NOT
downloaded from S3. They are installed into the package's `site-packages` directory.

| File | Size | Description |
|------|------|-------------|
| `gene_median_dictionary_gc104M.pkl` | ~1.5MB | Per-gene median expression values across the training corpus. Used by the tokenizer for rank-value encoding normalization. |
| `token_dictionary_gc104M.pkl` | ~416KB | Maps Ensembl gene IDs to integer token IDs (vocabulary). Also contains special tokens (`<cls>`, `<pad>`, `<mask>`). |
| `gene_name_id_dict_gc104M.pkl` | ~1.6MB | Maps gene symbols to Ensembl IDs and vice versa. Used by `InSilicoPerturberStats` for human-readable output. |
| `ensembl_mapping_dict_gc104M.pkl` | ~3.8MB | Maps between Ensembl ID versions and handles gene ID cross-referencing. Used by the tokenizer when `collapse_gene_ids=True`. |

These are resolved automatically as default parameter values. To inspect paths:

```python
from geneformer.tokenizer import (
    GENE_MEDIAN_FILE,
    TOKEN_DICTIONARY_FILE,
    ENSEMBL_MAPPING_FILE,
)
from geneformer.in_silico_perturber_stats import ENSEMBL_DICTIONARY_FILE

print(f"Gene medians:    {GENE_MEDIAN_FILE}")
print(f"Token dict:      {TOKEN_DICTIONARY_FILE}")
print(f"Ensembl mapping: {ENSEMBL_MAPPING_FILE}")
print(f"Gene name dict:  {ENSEMBL_DICTIONARY_FILE}")
```

---

## 6. Quantization Setup

### Simple 8-bit Quantization (Inference)

Pass `quantize=True` to `Classifier` or use 8-bit loading directly:

```python
from geneformer import Classifier

clf = Classifier(
    classifier="cell",
    quantize=True,  # 8-bit inference quantization
    cell_state_dict={"state_key": "cell_type", "states": "all"},
    model_version="V2",
)
```

### Custom Quantization with BitsAndBytesConfig + LoRA

For 4-bit QLoRA training, pass a dict with `bitsandbytes_config` and `peft_config`:

```python
from geneformer import Classifier
from transformers import BitsAndBytesConfig
from peft import LoraConfig

quantize_config = {
    "bitsandbytes_config": BitsAndBytesConfig(
        load_in_4bit=True,
        bnb_4bit_use_double_quant=True,
        bnb_4bit_quant_type="nf4",
        bnb_4bit_compute_dtype="float16",
    ),
    "peft_config": LoraConfig(
        r=8,
        lora_alpha=32,
        target_modules=["query", "value"],
        lora_dropout=0.05,
        bias="none",
    ),
}

clf = Classifier(
    classifier="cell",
    quantize=quantize_config,
    cell_state_dict={"state_key": "cell_type", "states": "all"},
    model_version="V2",
)
```

### Requirements

```bash
pip install bitsandbytes peft
```

Both packages are included in the `geneformer` conda environment. If installing
manually, ensure compatible versions:
- `bitsandbytes >= 0.41.0`
- `peft >= 0.6.0`

---

## 7. Verification

Run this script to verify that the environment and model are correctly set up:

```python
import torch
from geneformer import TranscriptomeTokenizer, EmbExtractor, Classifier

# Check CUDA
print(f"CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"GPU: {torch.cuda.get_device_name(0)}")
    print(f"GPU memory: {torch.cuda.get_device_properties(0).total_mem / 1e9:.1f} GB")

# Check gene dictionaries are accessible
from geneformer.tokenizer import GENE_MEDIAN_FILE, TOKEN_DICTIONARY_FILE
import os
print(f"Gene median dict exists: {os.path.exists(GENE_MEDIAN_FILE)}")
print(f"Token dict exists:       {os.path.exists(TOKEN_DICTIONARY_FILE)}")

# Verify model loads
from s3_cache import get_model_path
model_path = get_model_path(tier="V2-104M")
print(f"Model path: {model_path}")
print(f"Model files: {os.listdir(model_path)}")

from transformers import AutoModel
model = AutoModel.from_pretrained(model_path)
print(f"Model loaded: {type(model).__name__}")
print(f"Parameters: {sum(p.numel() for p in model.parameters()) / 1e6:.1f}M")
print(f"Hidden size: {model.config.hidden_size}")
print(f"Num layers: {model.config.num_hidden_layers}")
print(f"Max position: {model.config.max_position_embeddings}")

# Quick tokenizer test
tk = TranscriptomeTokenizer(model_version="V2")
print(f"Tokenizer initialized with model_input_size={tk.model_input_size}")

print("\nVerification complete.")
```

Expected output for V2-104M:

```
CUDA available: True
GPU: NVIDIA A10G (or similar)
GPU memory: 22.5 GB
Gene median dict exists: True
Token dict exists:       True
Model path: /tmp/geneformer/Geneformer-V2-104M/
Model files: ['config.json', 'model.safetensors', 'generation_config.json', 'training_args.bin']
Model loaded: BertModel
Parameters: 104.0M
Hidden size: 768
Num layers: 12
Max position: 4096
Tokenizer initialized with model_input_size=4096

Verification complete.
```

---

## 8. Conda Environment

The recommended environment for Geneformer is the pre-configured conda environment:

```bash
conda activate geneformer
```

### Environment Details

| Component | Version |
|-----------|---------|
| Python | 3.12 |
| PyTorch | 2.11.0+cu130 |
| transformers | >= 4.40.0 |
| datasets | >= 2.18.0 |
| geneformer | 0.1.0 |
| scanpy | >= 1.10.0 |
| anndata | >= 0.10.0 |
| loompy | >= 3.0.0 |
| bitsandbytes | >= 0.41.0 |
| peft | >= 0.6.0 |
| ray[tune] | >= 2.9.0 |
| optuna | >= 3.5.0 |

### Installing from Scratch

If the conda environment is not available:

```bash
conda create -n geneformer python=3.12
conda activate geneformer

# PyTorch with CUDA
pip install torch==2.11.0 --index-url https://download.pytorch.org/whl/cu130

# Geneformer and dependencies
pip install geneformer
pip install scanpy anndata loompy
pip install bitsandbytes peft  # For quantization
pip install ray[tune] optuna   # For hyperparameter optimization
```

### GPU Requirements

| Model | Minimum GPU Memory (Inference) | Minimum GPU Memory (Fine-tuning) |
|-------|-------------------------------|----------------------------------|
| V2-104M | ~4GB | ~16GB (full), ~8GB (QLoRA) |
| V2-316M | ~8GB | ~40GB (full), ~16GB (QLoRA) |
| V2-104M_CLcancer | ~4GB | ~16GB (full), ~8GB (QLoRA) |

Recommended instances:
- `ml.g5.xlarge` (1x A10G, 24GB) -- sufficient for V2-104M inference and fine-tuning
- `ml.g5.2xlarge` (1x A10G, 24GB) -- more CPU/RAM for larger datasets
- `ml.p4d.24xlarge` (8x A100, 40GB each) -- V2-316M full fine-tuning or large-scale perturbation
