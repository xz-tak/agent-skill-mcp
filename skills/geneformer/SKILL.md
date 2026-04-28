---
name: geneformer
description: >
  Complete toolkit for Geneformer (v0.1.0) — a foundational transformer model pretrained on
  ~104M single-cell transcriptomes for context-aware predictions in network biology (Nature 2023,
  Nature Computational Science 2026). Use this skill whenever users mention geneformer, transcriptome
  tokenizer, rank value encoding, single-cell transformer perturbation, in silico perturbation with
  gene deletion or overexpression, cell state classification with transformers, gene classification,
  multi-task cell classification, cell embeddings from Geneformer, gene embeddings, reverse phenotype
  perturbation, GeneformerPretrainer, TranscriptomeTokenizer, EmbExtractor, InSilicoPerturber,
  InSilicoPerturberStats, MTLClassifier, Classifier, DataCollatorForCellClassification,
  DataCollatorForGeneClassification, or goal state shift analysis. Also trigger when users want to
  tokenize scRNA-seq data into rank value encodings, extract CLS or cell or gene embeddings,
  fine-tune a cell state or gene classifier, run in silico perturbation (delete, overexpress,
  inhibit, activate), compute perturbation statistics (goal_state_shift, vs_null, mixture_model),
  run multi-task learning with Optuna, or use quantized Geneformer models. Covers V2 model series:
  V2-104M (default, 104M params), V2-316M (316M params), and V2-104M_CLcancer (cancer continual
  learning). Even if the user just says "gene perturbation transformer" or "tokenize single cell
  data for Geneformer" or "cell state shift analysis", this skill applies.
---

# Geneformer Toolkit

**Geneformer** is a foundational transformer model pretrained on ~104M single-cell transcriptomes. It uses rank value encoding to convert gene expression into token sequences, enabling transfer learning for cell classification, gene classification, in silico perturbation, and multi-task learning.

**Publications**:
- Theodoris et al. "Transfer learning enables predictions in network biology." *Nature* (2023). DOI: 10.1038/s41586-023-06139-9
- Theodoris et al. "Scaling and quantization of Geneformer." *Nature Computational Science* (2026).

**Package**: `geneformer` v0.1.0 (installed in conda env `geneformer`, Python 3.12, torch 2.11.0+cu130)

## Class Hierarchy

```
TranscriptomeTokenizer          # Converts raw scRNA-seq to rank value encoding .dataset
EmbExtractor                    # Extracts cell/gene embeddings, generates state embeddings
  get_embs()                    # Low-level embedding extraction function
Classifier                      # Fine-tunes cell state or gene classifiers
InSilicoPerturber               # Runs in silico perturbation (delete/overexpress/inhibit/activate)
InSilicoPerturberStats          # Statistical analysis of perturbation results
MTLClassifier                   # Multi-task cell classification with Optuna
GeneformerPretrainer            # Continue pretraining (MLM, extends HF Trainer)
DataCollatorForCellClassification  # Batching for cell classification
DataCollatorForGeneClassification  # Batching for gene classification
```

## Quick Start

### 1. Tokenize Data
```python
from geneformer import TranscriptomeTokenizer

tk = TranscriptomeTokenizer({"cell_type": "cell_type", "organ_major": "organ"}, nproc=16)
tk.tokenize_data("data_directory", "output_directory", "output_prefix", file_format="h5ad")
```

### 2. Extract Embeddings
```python
import sys
sys.path.insert(0, "/home/sagemaker-user/.claude/skills/geneformer/scripts")
from s3_cache import get_model_path
from geneformer import EmbExtractor

model_path = get_model_path(tier="V2-104M")  # syncs ~399MB on first use, instant after

embex = EmbExtractor(model_type="Pretrained", num_classes=0, emb_mode="cls",
                     max_ncells=2000, emb_layer=-1, forward_batch_size=256, nproc=16)
embs = embex.extract_embs(model_path, "tokenized_data.dataset",
                          "output_dir", "output_prefix")
embex.plot_embs(embs, plot_style="umap", output_directory="output_dir",
                output_prefix="output_prefix")
```

### 3. Classify Cells
```python
from geneformer import Classifier

cc = Classifier(classifier="cell",
                cell_state_dict={"state_key": "disease", "states": ["healthy", "diseased"]},
                filter_data={"cell_type": ["Cardiomyocyte"]},
                training_args={"num_train_epochs": 3, "learning_rate": 5e-5,
                               "per_device_train_batch_size": 12, "warmup_steps": 500},
                freeze_layers=2, num_crossval_splits=1,
                forward_batch_size=200, nproc=16)
cc.prepare_data("tokenized.dataset", "output_dir", "my_classifier")
all_metrics = cc.validate(model_path, "output_dir/my_classifier_labeled.dataset",
                          "output_dir/my_classifier_id_class_dict.pkl",
                          "output_dir", "my_classifier", predict_eval=True)
cc.plot_conf_mat({"Geneformer": all_metrics["conf_matrix"]}, "output_dir", "my_classifier")
```

### 4. In Silico Perturbation (Reverse Phenotype)
```python
from geneformer import EmbExtractor, InSilicoPerturber, InSilicoPerturberStats

# Step 1: Get state embeddings
cell_states = {"state_key": "disease", "start_state": "dcm",
               "goal_state": "nf", "alt_states": ["hcm"]}
embex = EmbExtractor(model_type="CellClassifier", num_classes=3,
                     filter_data={"cell_type": ["Cardiomyocyte"]},
                     max_ncells=1000, emb_layer=0, summary_stat="exact_mean",
                     forward_batch_size=256, nproc=16)
state_embs = embex.get_state_embs(cell_states, model_path,
                                   "tokenized.dataset", "output_dir", "embs")

# Step 2: Run perturbation
isp = InSilicoPerturber(perturb_type="delete", genes_to_perturb="all",
                        model_type="CellClassifier", num_classes=3, emb_mode="cls",
                        filter_data={"cell_type": ["Cardiomyocyte"]},
                        cell_states_to_model=cell_states, state_embs_dict=state_embs,
                        max_ncells=2000, emb_layer=0, forward_batch_size=400, nproc=16)
isp.perturb_data(model_path, "tokenized.dataset", "isp_output_dir", "isp_prefix")

# Step 3: Get statistics
ispstats = InSilicoPerturberStats(mode="goal_state_shift", genes_perturbed="all",
                                   cell_states_to_model=cell_states)
ispstats.get_stats("isp_output_dir", None, "stats_output_dir", "stats_prefix")
```

### 5. Multi-Task Classification
```python
from geneformer import MTLClassifier

mc = MTLClassifier(task_columns=["cell_type", "disease"],
                   pretrained_path=model_path,
                   train_path="train.dataset", val_path="val.dataset",
                   test_path="test.dataset",
                   model_save_path="mtl_model", results_dir="mtl_results",
                   trials_result_path="mtl_results/trials.txt",
                   study_name="mtl_study", batch_size=8, n_trials=15, epochs=10,
                   use_attention_pooling=True, use_task_weights=True, seed=42)
mc.run_optuna_study()
mc.load_and_evaluate_test_model()
```

## Model Setup

See [references/model_setup.md](references/model_setup.md) for download instructions and directory structure.

**S3 location**: `s3://tec-rnd-sci-dev-gi2/gi2-xz/models/geneformer/`

### S3 Cache (Recommended)

```python
import sys
sys.path.insert(0, "/home/sagemaker-user/.claude/skills/geneformer/scripts")
from s3_cache import get_model_path

model_path = get_model_path(tier="V2-104M")   # default, ~399MB
# model_path = get_model_path(tier="V2-316M")  # larger, ~1.2GB
# model_path = get_model_path(tier="V2-104M-CLcancer")  # cancer variant
```

Or from CLI:
```bash
python /home/sagemaker-user/.claude/skills/geneformer/scripts/s3_cache.py
python /home/sagemaker-user/.claude/skills/geneformer/scripts/s3_cache.py --tier V2-316M
python /home/sagemaker-user/.claude/skills/geneformer/scripts/s3_cache.py --clear
```

All pipeline scripts support `--s3` flag (auto-cache default V2-104M) and `--s3-tier` to choose variant.

### Model Variants

| Model | Params | Layers | Hidden | Input | Size |
|-------|--------|--------|--------|-------|------|
| **V2-104M** (default) | 104M | 12 | 768 | 4096 | ~399MB |
| V2-316M | 316M | 18 | 1152 | 4096 | ~1.2GB |
| V2-104M_CLcancer | 104M | 12 | 768 | 4096 | ~399MB |

## Data Preparation

All Geneformer workflows require tokenized data. Input files must be `.loom`, `.h5ad`, or `.zarr` with:

- **Gene IDs**: `ensembl_id` attribute (Ensembl IDs, e.g. ENSG00000141510)
- **Cell counts**: `n_counts` attribute (total read counts per cell)
- **Raw counts**: Un-normalized expression values (no log-transform, no scaling)
- **Optional**: `filter_pass` attribute (binary, 1 = include cell)
- **Optional**: Custom cell metadata attributes passed via `custom_attr_name_dict`

```python
import scanpy as sc
adata = sc.read_h5ad("raw_data.h5ad")

# Ensure ensembl_id exists
if "ensembl_id" not in adata.var.columns:
    adata.var["ensembl_id"] = adata.var.index  # if index is already Ensembl IDs

# Ensure n_counts exists
if "n_counts" not in adata.obs.columns:
    adata.obs["n_counts"] = adata.X.sum(axis=1)

adata.write("prepared_data.h5ad")
```

## Tokenization

`TranscriptomeTokenizer` converts raw scRNA-seq data to rank value encodings.

```python
TranscriptomeTokenizer(
    custom_attr_name_dict=None,   # dict: {loom_col: dataset_col} for custom metadata
    nproc=1,                       # number of processes
    chunk_size=512,                # chunk size for h5ad/zarr
    model_input_size=4096,         # V2 default (DO NOT CHANGE for V2)
    special_token=True,            # V2 default: adds CLS/EOS tokens
    collapse_gene_ids=True,        # merge duplicate Ensembl IDs
    gene_median_file=GENE_MEDIAN_FILE,
    token_dictionary_file=TOKEN_DICTIONARY_FILE,
    gene_mapping_file=ENSEMBL_MAPPING_FILE,
)
```

**Key methods**:
- `tokenize_data(data_directory, output_directory, output_prefix, file_format="loom")` — tokenize all files in a directory
- `tokenize_anndata(adata)` — tokenize a single AnnData object, returns tokenized dict
- Output: HuggingFace `.dataset` format with `input_ids`, `length`, and custom metadata

## Embedding Extraction

`EmbExtractor` extracts cell or gene embeddings from pretrained or fine-tuned models.

```python
EmbExtractor(
    model_type="Pretrained",     # {"Pretrained", "CellClassifier", "GeneClassifier"}
    num_classes=0,               # number of classes (0 for pretrained)
    emb_mode="cls",              # {"cls", "cell", "gene"} — see decision tree below
    cell_emb_style="mean_pool",  # only option currently
    gene_emb_style="mean_pool",  # {"mean_pool", "all"}
    filter_data=None,            # dict to filter cells, e.g. {"cell_type": ["T cell"]}
    max_ncells=1000,             # max cells to process (None = all)
    emb_layer=-1,                # -1: 2nd-to-last layer (general), 0: last layer (task-specific)
    emb_label=None,              # list of column names to add as labels
    labels_to_plot=None,         # labels for coloring plots
    forward_batch_size=100,
    nproc=4,
    summary_stat=None,           # {None, "mean", "median", "exact_mean", "exact_median"}
    model_version="V2",
)
```

### Embedding Mode Decision Tree

```
What embedding do you need?
├── Cell-level embedding (one vector per cell)
│   ├── CLS token embedding → emb_mode="cls" (V2 recommended)
│   └── Mean-pooled gene embeddings → emb_mode="cell"
└── Gene-level embedding (one vector per gene per cell)
    └── emb_mode="gene"
```

**Key methods**:
- `extract_embs(model_directory, input_data_file, output_directory, output_prefix, output_torch_embs=False)` — extract and save embeddings
- `get_state_embs(cell_states_to_model, model_directory, input_data_file, output_directory, output_prefix)` — compute state embeddings for perturbation
- `plot_embs(embs, plot_style, output_directory, output_prefix)` — plot as "heatmap" or "umap"

## Cell & Gene Classification

`Classifier` fine-tunes Geneformer for cell state or gene classification.

```python
Classifier(
    classifier="cell",           # {"cell", "gene"}
    quantize=False,              # True for 8-bit, or dict with BitsAndBytesConfig + LoraConfig
    cell_state_dict=None,        # {"state_key": "disease", "states": ["healthy", "diseased"]}
    gene_class_dict=None,        # {"Label_A": ["ENSG...", ...], "Label_B": [...]}
    filter_data=None,            # dict to subset cells
    rare_threshold=0,            # remove cell states below this fraction
    max_ncells=None,             # max cells for fine-tuning
    max_ncells_per_class=None,   # max cells per class
    training_args=None,          # HuggingFace TrainingArguments as dict
    freeze_layers=0,             # number of transformer layers to freeze
    num_crossval_splits=1,       # {0: no split, 1: train/eval, 5: 5-fold CV}
    split_sizes={"train": 0.8, "valid": 0.1, "test": 0.1},
    no_eval=False,
    forward_batch_size=100,
    model_version="V2",
    nproc=4, ngpu=1,
)
```

**Workflow**: `prepare_data()` → `validate()` (or `train_all_data()`) → `plot_conf_mat()` / `plot_roc()` / `plot_predictions()`

**Key methods**:
- `prepare_data(input_data_file, output_directory, output_prefix)` — label and split data
- `validate(model_directory, prepared_input_data_file, id_class_dict_file, output_directory, output_prefix, predict_eval=True)` — train + cross-validate
- `train_all_data(model_directory, prepared_input_data_file, id_class_dict_file, output_directory, output_prefix)` — train on all data (no eval split)
- `evaluate_saved_model(model_directory, id_class_dict_file, test_data_file, output_directory, output_prefix)` — evaluate a saved model on test data
- `plot_conf_mat(conf_mat_dict, output_directory, output_prefix, custom_class_order=None)`
- `plot_roc(roc_metric_dict, model_style_dict, title, output_directory, output_prefix)`
- `plot_predictions(predictions_file, id_class_dict_file, title, output_directory, output_prefix)`

## In Silico Perturbation

`InSilicoPerturber` simulates gene perturbations and measures embedding shifts.

```python
InSilicoPerturber(
    perturb_type="delete",       # {"delete", "overexpress", "inhibit", "activate"}
    perturb_rank_shift=None,     # {None, 1, 2, 3} for inhibit/activate
    genes_to_perturb="all",      # "all" or list of Ensembl IDs
    combos=0,                    # 0: individual genes, 1: pairwise combinations
    anchor_gene=None,            # Ensembl ID for combination anchor
    model_type="Pretrained",     # see model_type options below
    num_classes=0,
    emb_mode="cls",              # {"cls", "cell", "cls_and_gene", "cell_and_gene"}
    filter_data=None,
    cell_states_to_model=None,   # for goal state shift analysis
    state_embs_dict=None,        # from EmbExtractor.get_state_embs()
    max_ncells=None,
    emb_layer=-1,
    forward_batch_size=100,
    nproc=4,
    model_version="V2",
)
```

### model_type Options

| model_type | Model Class | Quantized | Use Case |
|---|---|---|---|
| `"Pretrained"` | BertForMaskedLM | No | Zero-shot perturbation |
| `"CellClassifier"` | BertForSequenceClassification | No | Fine-tuned cell classifier |
| `"GeneClassifier"` | BertForTokenClassification | No | Fine-tuned gene classifier |
| `"MTLCellClassifier"` | BertForMaskedLM | No | Multi-task fine-tuned |
| `"Pretrained-Quantized"` | BertForMaskedLM | 8-bit | Memory-efficient zero-shot |
| `"MTLCellClassifier-Quantized"` | BertForMaskedLM | 8-bit | Memory-efficient MTL |

**Perturbation pipeline** (3 steps):
1. `EmbExtractor.get_state_embs()` → `state_embs_dict`
2. `InSilicoPerturber.perturb_data(model_dir, input_data, output_dir, prefix)` → raw pickle files
3. `InSilicoPerturberStats.get_stats()` → final CSV

## Perturbation Statistics

`InSilicoPerturberStats` aggregates and analyzes perturbation results.

```python
InSilicoPerturberStats(
    mode="mixture_model",        # see mode options below
    genes_perturbed="all",       # "all" or list of Ensembl IDs
    combos=0,                    # must match InSilicoPerturber setting
    anchor_gene=None,
    cell_states_to_model=None,   # must match InSilicoPerturber setting
    pickle_suffix="_raw.pickle",
    model_version="V2",
)
```

### Stats Mode Options

| Mode | Use Case | Requires cell_states_to_model |
|---|---|---|
| `"goal_state_shift"` | Genes whose deletion shifts cells toward goal state | Yes |
| `"vs_null"` | Perturbation vs. null distribution | No |
| `"mixture_model"` | Impact vs. no-impact mixture model (undirected) | No |
| `"aggregate_data"` | Aggregate cosine shifts for single perturbation | No |
| `"aggregate_gene_shifts"` | Aggregate gene-level shifts across perturbations | No |

**Method**: `get_stats(input_data_directory, null_dist_data_directory, output_directory, output_prefix)` — outputs CSV with gene rankings and significance.

## Multi-Task Learning

`MTLClassifier` fine-tunes on multiple classification tasks simultaneously with Optuna hyperparameter optimization.

```python
MTLClassifier(
    task_columns=None,           # list of .dataset column names for each task
    pretrained_path=None,        # path to pretrained model
    train_path=None, val_path=None, test_path=None,
    model_save_path=None,        # where to save the fine-tuned model
    results_dir=None,
    trials_result_path=None,     # path for Optuna trial results
    study_name=None,
    batch_size=4,
    n_trials=15,                 # number of Optuna trials
    epochs=None,
    use_attention_pooling=None,  # attention-based pooling (recommended)
    use_task_weights=None,       # learnable task weights
    distributed_training=None,   # DDP support
    hyperparameters=None,        # custom hyperparameter search space
    manual_hyperparameters=None, # fixed hyperparameters (skip Optuna)
    use_manual_hyperparameters=None,
    gradient_clipping=None,
    seed=None,
    gradient_accumulation_steps=None,
)
```

**Methods**:
- `run_optuna_study()` — hyperparameter optimization + training
- `load_and_evaluate_test_model()` — evaluate on test data

**Note**: Input dataset must contain a `unique_cell_id` column.

## Pretraining

`GeneformerPretrainer` extends HuggingFace's `Trainer` for masked language modeling (15% masking). Used for continued pretraining on new data — not for pretraining from scratch (requires Genecorpus dataset).

```python
from geneformer import GeneformerPretrainer

trainer = GeneformerPretrainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_dataset,
    token_dictionary=token_dict,
    example_lengths_file="lengths.pkl",
)
trainer.train()
```

## Data Collators

- `DataCollatorForCellClassification()` — dynamic padding + attention masks for cell classification
- `DataCollatorForGeneClassification()` — dynamic padding + label handling for gene-level tasks

Used internally by `Classifier`; pass to HuggingFace `Trainer` for custom training loops.

## Constants

| Constant | Description |
|----------|-------------|
| `GENE_MEDIAN_FILE` | Gene median dictionary (gc104M) |
| `TOKEN_DICTIONARY_FILE` | Token dictionary mapping Ensembl IDs to tokens (gc104M) |
| `ENSEMBL_DICTIONARY_FILE` | Gene name to Ensembl ID mapping (gc104M) |
| `ENSEMBL_MAPPING_FILE` | Ensembl ID mapping for gene symbol conversion (gc104M) |

All bundled with the pip package at `~/.conda/envs/geneformer/lib/python3.12/site-packages/geneformer/`.

## Quantization

Geneformer supports 8-bit and 4-bit quantization for memory-efficient inference and training.

### For Classification (Classifier)
```python
# Simple: 8-bit quantization (inference)
cc = Classifier(classifier="cell", quantize=True, ...)

# Custom: BitsAndBytesConfig + LoRA
from transformers import BitsAndBytesConfig
from peft import LoraConfig

cc = Classifier(classifier="cell",
    quantize={
        "bnb_config": BitsAndBytesConfig(load_in_8bit=True),
        "peft_config": LoraConfig(r=8, lora_alpha=32, target_modules=["query", "value"])
    }, ...)
```

### For Perturbation (InSilicoPerturber)
```python
# Use quantized model_type suffix
isp = InSilicoPerturber(model_type="Pretrained-Quantized", ...)    # 8-bit pretrained
isp = InSilicoPerturber(model_type="MTLCellClassifier-Quantized", ...)  # 8-bit MTL
```

**Requirements**: `pip install bitsandbytes peft`

## Performance Tips

| Parameter | Impact | Recommendation |
|-----------|--------|----------------|
| `forward_batch_size` | GPU memory / speed | Start at 100, increase if GPU allows |
| `max_ncells` | Processing time / memory | Use None for full data, 1000-5000 for exploration |
| `nproc` | CPU parallelization | Set to available CPU cores |
| `emb_layer` | Embedding specificity | -1 for general, 0 for task-specific |
| `freeze_layers` | Training speed / overfitting | 2-4 layers for small datasets |
| `quantize` | GPU memory (~50% reduction) | Use for V2-316M on <24GB GPU |

## Troubleshooting

See [references/troubleshooting.md](references/troubleshooting.md) for common errors and solutions.

## Scripts

| Script | Purpose |
|--------|---------|
| `scripts/s3_cache.py` | S3-to-local model cache |
| `scripts/ensure_env.py` | Validate conda environment |
| `scripts/tokenize_data.py` | CLI tokenization pipeline |
| `scripts/extract_embeddings.py` | CLI embedding extraction |
| `scripts/classify_cells.py` | CLI cell/gene classification |
| `scripts/perturb_cells.py` | CLI in silico perturbation |
| `scripts/perturb_stats.py` | CLI perturbation statistics |
| `scripts/mtl_classify.py` | CLI multi-task classification |
