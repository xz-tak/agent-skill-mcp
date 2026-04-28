# Geneformer API Reference (V2)

Complete API documentation for all Geneformer classes and utilities.
All classes default to `model_version="V2"` unless otherwise specified.

---

## Table of Contents

1. [TranscriptomeTokenizer](#1-transcriptometokenizer)
2. [EmbExtractor](#2-embextractor)
3. [get_embs (standalone function)](#3-get_embs-standalone-function)
4. [Classifier](#4-classifier)
5. [InSilicoPerturber](#5-insilicoperturer)
6. [InSilicoPerturberStats](#6-insilicoperturberstats)
7. [MTLClassifier](#7-mtlclassifier)
8. [GeneformerPretrainer](#8-geneformerpretrainer)
9. [DataCollators](#9-datacollators)
10. [Utility Functions](#10-utility-functions)
11. [Constants](#11-constants)

---

## 1. TranscriptomeTokenizer

Tokenizes single-cell transcriptomic data into rank-value encoded token sequences
for Geneformer input. Converts raw count matrices from `.loom`, `.h5ad`, or `.zarr`
files into HuggingFace Datasets with `input_ids` representing gene expression ranks.

### Module

`geneformer.tokenizer`

### Import

```python
from geneformer import TranscriptomeTokenizer
```

### Constructor

```python
TranscriptomeTokenizer(
    custom_attr_name_dict=None,
    nproc=1,
    chunk_size=512,
    model_input_size=4096,
    special_token=True,
    collapse_gene_ids=True,
    gene_median_file=GENE_MEDIAN_FILE,
    token_dictionary_file=TOKEN_DICTIONARY_FILE,
    gene_mapping_file=ENSEMBL_MAPPING_FILE,
    model_version="V2",
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `custom_attr_name_dict` | `dict` or `None` | `None` | Dictionary mapping custom column names to standard names. Keys are standard names (`"cell_type"`, `"disease"`, etc.), values are the column names in your data. These metadata columns are carried through to the output dataset. |
| `nproc` | `int` | `1` | Number of processes for parallel tokenization. Set higher for large datasets. |
| `chunk_size` | `int` | `512` | Number of cells processed per chunk during tokenization. Controls memory usage. |
| `model_input_size` | `int` | `4096` | Maximum number of tokens per cell. Gene tokens are truncated to this length after rank ordering. Must match the model's max input size. |
| `special_token` | `bool` | `True` | Whether to prepend `<cls>` token to each cell's token sequence. Required for V2 models. |
| `collapse_gene_ids` | `bool` | `True` | Whether to collapse Ensembl gene IDs that map to the same token. When `True`, expression values of duplicate gene IDs are summed. |
| `gene_median_file` | `str` or `Path` | `GENE_MEDIAN_FILE` | Path to the gene median dictionary pickle file. Used for rank-value encoding normalization. V2 default: `gene_median_dictionary_gc104M.pkl`. |
| `token_dictionary_file` | `str` or `Path` | `TOKEN_DICTIONARY_FILE` | Path to the token dictionary pickle file mapping Ensembl IDs to token IDs. V2 default: `token_dictionary_gc104M.pkl`. |
| `gene_mapping_file` | `str` or `Path` | `ENSEMBL_MAPPING_FILE` | Path to the Ensembl mapping dictionary pickle file for gene ID versioning. V2 default: `ensembl_mapping_dict_gc104M.pkl`. |
| `model_version` | `str` | `"V2"` | Model version. Determines which dictionary files to use. Use `"V2"` for all current models. |

### Methods

#### `tokenize_data`

```python
tokenize_data(
    data_directory,
    output_directory,
    output_prefix,
    file_format="loom",
)
```

Tokenize all files in a directory and save as a HuggingFace Dataset.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_directory` | `str` or `Path` | *required* | Directory containing input data files. |
| `output_directory` | `str` or `Path` | *required* | Directory where the output `.dataset` will be saved. |
| `output_prefix` | `str` | *required* | Prefix for output dataset directory name. Output is saved as `{output_prefix}.dataset`. |
| `file_format` | `str` | `"loom"` | Input file format. One of `"loom"`, `"h5ad"`, or `"zarr"`. |

**Returns:** `datasets.Dataset` -- HuggingFace Dataset with columns `input_ids`, `length`, and any custom metadata.

#### `tokenize_anndata`

```python
tokenize_anndata(adata)
```

Tokenize a single AnnData object in memory.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | `anndata.AnnData` | *required* | AnnData object with raw counts in `.X`, Ensembl IDs in `.var["ensembl_id"]`, and `n_counts` in `.obs`. |

**Returns:** `datasets.Dataset` -- HuggingFace Dataset with columns `input_ids`, `length`, and any custom metadata from `.obs`.

#### `tokenize_files`

```python
tokenize_files(data_directory, file_format)
```

Tokenize all files in a directory and return the tokenized cells and metadata without saving.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_directory` | `str` or `Path` | *required* | Directory containing input files. |
| `file_format` | `str` | *required* | File format: `"loom"`, `"h5ad"`, or `"zarr"`. |

**Returns:** `tuple[list, list]` -- `(tokenized_cells, cell_metadata)` where `tokenized_cells` is a list of token ID lists and `cell_metadata` is a list of metadata dicts.

#### `tokenize_loom`

```python
tokenize_loom(loom_data_directory)
```

Tokenize all `.loom` files in a directory. Convenience wrapper around `tokenize_files`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `loom_data_directory` | `str` or `Path` | *required* | Directory containing `.loom` files. |

**Returns:** `tuple[list, list]` -- Same as `tokenize_files`.

#### `create_dataset`

```python
create_dataset(
    tokenized_cells,
    cell_metadata,
    use_generator=False,
)
```

Create a HuggingFace Dataset from pre-tokenized cells and metadata.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tokenized_cells` | `list` | *required* | List of token ID lists from `tokenize_files` or `tokenize_loom`. |
| `cell_metadata` | `list` | *required* | List of metadata dicts corresponding to each cell. |
| `use_generator` | `bool` | `False` | Whether to use a generator for memory-efficient dataset creation. |

**Returns:** `datasets.Dataset`

### Input Requirements

- **File formats:** `.loom`, `.h5ad`, or `.zarr`
- **Expression data:** Raw counts (not normalized, not log-transformed) in the main matrix
- **Gene identifiers:** Ensembl IDs stored in:
  - `.loom`: row attribute `"ensembl_id"`
  - `.h5ad`: `adata.var["ensembl_id"]`
  - `.zarr`: var column `"ensembl_id"`
- **Total counts:** Pre-computed `n_counts` per cell in obs/col attributes
- Cells with zero expressed genes after filtering are dropped

### Output Format

HuggingFace `.dataset` directory containing:

| Column | Type | Description |
|--------|------|-------------|
| `input_ids` | `list[int]` | Rank-ordered gene token IDs, most expressed first. Prepended with `<cls>` token if `special_token=True`. |
| `length` | `int` | Number of tokens (including `<cls>` if present). |
| *custom columns* | varies | Any columns specified in `custom_attr_name_dict`. |

### Usage Example

```python
from geneformer import TranscriptomeTokenizer

tk = TranscriptomeTokenizer(
    custom_attr_name_dict={"cell_type": "cell_type", "disease": "disease"},
    nproc=8,
    model_input_size=4096,
    special_token=True,
    model_version="V2",
)

# Tokenize from directory
dataset = tk.tokenize_data(
    data_directory="/path/to/h5ad_files",
    output_directory="/path/to/output",
    output_prefix="my_dataset",
    file_format="h5ad",
)

# Tokenize a single AnnData in memory
import scanpy as sc
adata = sc.read_h5ad("/path/to/file.h5ad")
dataset = tk.tokenize_anndata(adata)
```

---

## 2. EmbExtractor

Extracts cell or gene embeddings from pretrained or fine-tuned Geneformer models.
Supports CLS token embeddings, mean-pooled cell embeddings, and per-gene embeddings.

### Module

`geneformer.emb_extractor`

### Import

```python
from geneformer import EmbExtractor
```

### Constructor

```python
EmbExtractor(
    model_type="Pretrained",
    num_classes=0,
    emb_mode="cls",
    cell_emb_style="mean_pool",
    gene_emb_style="mean_pool",
    filter_data=None,
    max_ncells=1000,
    emb_layer=-1,
    emb_label=None,
    labels_to_plot=None,
    forward_batch_size=100,
    nproc=4,
    summary_stat=None,
    save_tdigest=False,
    model_version="V2",
    token_dictionary_file=None,
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model_type` | `str` | `"Pretrained"` | Type of model. One of `"Pretrained"`, `"GeneClassifier"`, `"CellClassifier"`. Determines model loading behavior. |
| `num_classes` | `int` | `0` | Number of output classes. Set to 0 for pretrained models; set to the class count for fine-tuned classifiers. |
| `emb_mode` | `str` | `"cls"` | Embedding mode. `"cls"`: CLS token embedding only. `"cell"`: mean-pooled embedding across all gene tokens. `"gene"`: per-gene embeddings. |
| `cell_emb_style` | `str` | `"mean_pool"` | Pooling strategy for cell embeddings when `emb_mode="cell"`. Currently only `"mean_pool"` is supported. |
| `gene_emb_style` | `str` | `"mean_pool"` | Aggregation for gene embeddings across cells. `"mean_pool"`: average across cells. `"all"`: return all per-cell gene embeddings without aggregation. |
| `filter_data` | `dict` or `None` | `None` | Dictionary for filtering cells before extraction. Keys are column names, values are lists of values to keep. Example: `{"cell_type": ["T cell", "B cell"]}`. |
| `max_ncells` | `int` | `1000` | Maximum number of cells to process. Cells are randomly sampled if the dataset exceeds this limit. |
| `emb_layer` | `int` | `-1` | Transformer layer to extract embeddings from. `-1`: last layer. `0`: input embedding layer. Any integer in `[-num_layers, num_layers-1]`. |
| `emb_label` | `str` or `list[str]` or `None` | `None` | Column name(s) in the dataset to use as labels for the output embeddings DataFrame index or for grouping. |
| `labels_to_plot` | `list[str]` or `None` | `None` | List of label column names to use when generating UMAP/tSNE plots. |
| `forward_batch_size` | `int` | `100` | Batch size for forward passes through the model. |
| `nproc` | `int` | `4` | Number of processes for data loading. |
| `summary_stat` | `str` or `None` | `None` | Summary statistic for aggregating gene embeddings across cells. One of `None`, `"mean"`, `"median"`, `"exact_mean"`, `"exact_median"`. The `"exact_*"` variants compute full statistics (slower but precise); non-exact variants use t-digest approximation. |
| `save_tdigest` | `bool` | `False` | Whether to save the t-digest data structure for approximate summary statistics. Useful for incremental updates. |
| `model_version` | `str` | `"V2"` | Model version string. |
| `token_dictionary_file` | `str` or `Path` or `None` | `None` | Path to token dictionary. If `None`, uses the default V2 dictionary. |

### Methods

#### `extract_embs`

```python
extract_embs(
    model_directory,
    input_data_file,
    output_directory,
    output_prefix,
    output_torch_embs=False,
    cell_state=None,
)
```

Extract embeddings from a model and save to disk.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model_directory` | `str` or `Path` | *required* | Path to the pretrained or fine-tuned model directory. |
| `input_data_file` | `str` or `Path` | *required* | Path to the tokenized `.dataset` directory. |
| `output_directory` | `str` or `Path` | *required* | Directory for output files. |
| `output_prefix` | `str` | *required* | Prefix for output filenames. |
| `output_torch_embs` | `bool` | `False` | If `True`, save embeddings as a `.pt` PyTorch tensor file in addition to the pandas DataFrame. |
| `cell_state` | `dict` or `None` | `None` | Cell state specification for state-aware embedding extraction. Same format as `cell_states_to_model` in `InSilicoPerturber`. |

**Returns:** `pandas.DataFrame` -- DataFrame with embeddings. Shape is `(n_cells, emb_dim)` for cell/cls mode, or dict of DataFrames for gene mode.

#### `get_state_embs`

```python
get_state_embs(
    cell_states_to_model,
    model_directory,
    input_data_file,
    output_directory,
    output_prefix,
    output_torch_embs=True,
)
```

Extract state-specific embeddings for use with `InSilicoPerturber`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cell_states_to_model` | `dict` | *required* | Cell state specification dict. Format: `{"state_key": "column_name", "start_state": "state_A", "goal_state": "state_B", "alt_states": ["state_C"]}`. |
| `model_directory` | `str` or `Path` | *required* | Path to model directory. |
| `input_data_file` | `str` or `Path` | *required* | Path to tokenized dataset. |
| `output_directory` | `str` or `Path` | *required* | Output directory. |
| `output_prefix` | `str` | *required* | Output filename prefix. |
| `output_torch_embs` | `bool` | `True` | Whether to save as `.pt` files. |

**Returns:** `dict` -- Dictionary mapping state names to embedding tensors. Pass this as `state_embs_dict` to `InSilicoPerturber`.

#### `plot_embs`

```python
plot_embs(
    embs,
    plot_style,
    output_directory,
    output_prefix,
    max_ncells_to_plot=1000,
    kwargs_dict=None,
)
```

Generate dimensionality reduction plots of extracted embeddings.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `embs` | `pandas.DataFrame` | *required* | Embeddings DataFrame from `extract_embs`. |
| `plot_style` | `str` | *required* | Plot type. One of `"umap"`, `"tsne"`, `"pca"`. |
| `output_directory` | `str` or `Path` | *required* | Output directory for plot files. |
| `output_prefix` | `str` | *required* | Filename prefix for plots. |
| `max_ncells_to_plot` | `int` | `1000` | Maximum cells to include in the plot. |
| `kwargs_dict` | `dict` or `None` | `None` | Additional keyword arguments passed to the underlying plotting function (e.g., UMAP parameters). |

**Returns:** `None` -- Saves plot files to `output_directory`.

### Usage Example

```python
from geneformer import EmbExtractor

# Cell embeddings (CLS mode)
emb_extractor = EmbExtractor(
    model_type="Pretrained",
    emb_mode="cls",
    emb_layer=-1,
    max_ncells=5000,
    forward_batch_size=100,
    model_version="V2",
)

embs_df = emb_extractor.extract_embs(
    model_directory="/path/to/Geneformer-V2-104M",
    input_data_file="/path/to/tokenized.dataset",
    output_directory="/path/to/output",
    output_prefix="my_embs",
    output_torch_embs=True,
)

# Plot UMAP
emb_extractor.plot_embs(
    embs=embs_df,
    plot_style="umap",
    output_directory="/path/to/output",
    output_prefix="my_embs",
)

# State embeddings for perturbation analysis
state_embs = emb_extractor.get_state_embs(
    cell_states_to_model={
        "state_key": "disease",
        "start_state": "dcm",
        "goal_state": "nf",
        "alt_states": ["hcm"],
    },
    model_directory="/path/to/Geneformer-V2-104M",
    input_data_file="/path/to/tokenized.dataset",
    output_directory="/path/to/output",
    output_prefix="state_embs",
)
```

---

## 3. get_embs (standalone function)

Low-level function for extracting embeddings from a loaded model. Used internally
by `EmbExtractor` but exposed for advanced use cases requiring direct model access.

### Module

`geneformer.emb_extractor`

### Import

```python
from geneformer import get_embs
```

### Signature

```python
get_embs(
    model,
    filtered_input_data,
    emb_mode,
    layer_to_quant,
    pad_token_id,
    forward_batch_size,
    token_gene_dict,
    special_token=False,
    summary_stat=None,
    silent=False,
    save_tdigest=False,
    tdigest_path=None,
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model` | `transformers.PreTrainedModel` | *required* | Loaded Geneformer model (pretrained or fine-tuned). |
| `filtered_input_data` | `datasets.Dataset` | *required* | Tokenized dataset, optionally filtered. Must contain `input_ids` column. |
| `emb_mode` | `str` | *required* | Embedding mode: `"cls"`, `"cell"`, or `"gene"`. |
| `layer_to_quant` | `int` | *required* | Transformer layer index to extract from. Typically `-1` (last) or `0` (first). |
| `pad_token_id` | `int` | *required* | Token ID used for padding. Obtain from the token dictionary. |
| `forward_batch_size` | `int` | *required* | Batch size for model forward passes. |
| `token_gene_dict` | `dict` | *required* | Reverse mapping from token IDs to gene identifiers. |
| `special_token` | `bool` | `False` | Whether the input contains a `<cls>` special token at position 0. |
| `summary_stat` | `str` or `None` | `None` | For gene mode: aggregation statistic across cells. One of `None`, `"mean"`, `"median"`, `"exact_mean"`, `"exact_median"`. |
| `silent` | `bool` | `False` | If `True`, suppress progress bars. |
| `save_tdigest` | `bool` | `False` | Whether to save t-digest for incremental summary stats. |
| `tdigest_path` | `str` or `Path` or `None` | `None` | Path to save/load t-digest data. |

### Returns

- **Cell/CLS mode:** `torch.Tensor` of shape `(n_cells, emb_dim)`.
- **Gene mode with `summary_stat`:** `dict` mapping gene Ensembl IDs to aggregated embedding vectors.
- **Gene mode without `summary_stat`:** `dict` mapping gene Ensembl IDs to `torch.Tensor` of shape `(n_cells_expressing, emb_dim)`.

### Usage Example

```python
import pickle
import torch
from datasets import load_from_disk
from transformers import AutoModel
from geneformer import get_embs

# Load model and data
model = AutoModel.from_pretrained("/path/to/Geneformer-V2-104M")
model.eval()
if torch.cuda.is_available():
    model = model.to("cuda")

dataset = load_from_disk("/path/to/tokenized.dataset")

# Load token dictionary and build reverse mapping
with open("/path/to/token_dictionary_gc104M.pkl", "rb") as f:
    token_dict = pickle.load(f)
token_gene_dict = {v: k for k, v in token_dict.items()}

pad_token_id = token_dict.get("<pad>", 0)

embs = get_embs(
    model=model,
    filtered_input_data=dataset,
    emb_mode="cls",
    layer_to_quant=-1,
    pad_token_id=pad_token_id,
    forward_batch_size=64,
    token_gene_dict=token_gene_dict,
    special_token=True,
)
# embs: torch.Tensor of shape (n_cells, 768) for V2-104M
```

---

## 4. Classifier

Fine-tunes Geneformer for cell type classification or gene classification tasks.
Supports hyperparameter optimization via Ray Tune, cross-validation, and quantized training.

### Module

`geneformer.classifier`

### Import

```python
from geneformer import Classifier
```

### Constructor

```python
Classifier(
    classifier="cell",
    quantize=False,
    cell_state_dict=None,
    gene_class_dict=None,
    filter_data=None,
    rare_threshold=0,
    max_ncells=None,
    max_ncells_per_class=None,
    training_args=None,
    ray_config=None,
    freeze_layers=0,
    num_crossval_splits=1,
    split_sizes={"train": 0.8, "valid": 0.1, "test": 0.1},
    stratify_splits_col=None,
    no_eval=False,
    forward_batch_size=100,
    model_version="V2",
    token_dictionary_file=None,
    nproc=4,
    ngpu=1,
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `classifier` | `str` | `"cell"` | Classification type. `"cell"`: cell-level classification (e.g., cell type, disease state). `"gene"`: gene-level classification using per-gene embeddings. |
| `quantize` | `bool` or `dict` | `False` | Quantization config. `False`: full precision. `True`: 8-bit quantization for inference. `dict`: custom config with `"bitsandbytes_config"` (BitsAndBytesConfig) and optionally `"peft_config"` (LoraConfig) keys. |
| `cell_state_dict` | `dict` or `None` | `None` | For cell classification: specifies the target label column and classes. Format: `{"state_key": "column_name", "states": ["class_A", "class_B"]}` or `{"state_key": "column_name", "states": "all"}` to use all unique values. |
| `gene_class_dict` | `dict` or `None` | `None` | For gene classification: maps class labels to lists of Ensembl gene IDs. Format: `{"Label_A": ["ENSG00000123456", ...], "Label_B": ["ENSG00000654321", ...]}`. |
| `filter_data` | `dict` or `None` | `None` | Filter cells before training. Keys are column names, values are lists of values to keep. |
| `rare_threshold` | `int` | `0` | Minimum number of samples per class. Classes with fewer samples are dropped. Set to 0 to keep all classes. |
| `max_ncells` | `int` or `None` | `None` | Maximum total number of cells. `None` for no limit. |
| `max_ncells_per_class` | `int` or `None` | `None` | Maximum cells per class for balanced training. `None` for no limit. |
| `training_args` | `dict` or `None` | `None` | Custom HuggingFace `TrainingArguments` as a dict. If `None`, uses default training args (see `classifier_utils.get_default_train_args`). |
| `ray_config` | `dict` or `None` | `None` | Ray Tune hyperparameter search space configuration. If `None`, hyperparameter optimization is disabled. |
| `freeze_layers` | `int` | `0` | Number of transformer encoder layers to freeze from the bottom. `0`: no freezing (full fine-tuning). Freezing lower layers speeds training and can prevent overfitting. |
| `num_crossval_splits` | `int` | `1` | Number of cross-validation splits. `0` or `1`: no cross-validation. `5`: 5-fold cross-validation. |
| `split_sizes` | `dict` | `{"train": 0.8, "valid": 0.1, "test": 0.1}` | Proportions for train/validation/test splits. Must sum to 1.0. |
| `stratify_splits_col` | `str` or `None` | `None` | Column name to stratify splits on. Ensures proportional class representation in each split. |
| `no_eval` | `bool` | `False` | If `True`, skip evaluation during training. Useful for final training on all data. |
| `forward_batch_size` | `int` | `100` | Batch size for forward passes during evaluation. |
| `model_version` | `str` | `"V2"` | Model version. |
| `token_dictionary_file` | `str` or `Path` or `None` | `None` | Path to token dictionary. `None` uses V2 default. |
| `nproc` | `int` | `4` | Number of data loading processes. |
| `ngpu` | `int` | `1` | Number of GPUs for training. |

### Methods

#### `prepare_data`

```python
prepare_data(
    input_data_file,
    output_directory,
    output_prefix,
    split_id_dict=None,
)
```

Prepare and split the tokenized dataset for training.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_data_file` | `str` or `Path` | *required* | Path to the tokenized `.dataset` directory. |
| `output_directory` | `str` or `Path` | *required* | Output directory for prepared data splits. |
| `output_prefix` | `str` | *required* | Prefix for output files. |
| `split_id_dict` | `dict` or `None` | `None` | Pre-defined split assignments. Keys: `"train"`, `"valid"`, `"test"`. Values: lists of cell indices. If `None`, splits are created according to `split_sizes`. |

**Returns:** `dict` -- Dictionary with paths to prepared data files and the `id_class_dict`.

#### `validate`

```python
validate(
    model_directory,
    prepared_input_data_file,
    id_class_dict_file,
    output_directory,
    output_prefix,
    predict_eval=True,
)
```

Train on the train split and evaluate on the validation/test splits.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model_directory` | `str` or `Path` | *required* | Path to pretrained model. |
| `prepared_input_data_file` | `str` or `Path` | *required* | Path to prepared dataset from `prepare_data`. |
| `id_class_dict_file` | `str` or `Path` | *required* | Path to the ID-to-class mapping pickle file from `prepare_data`. |
| `output_directory` | `str` or `Path` | *required* | Output directory. |
| `output_prefix` | `str` | *required* | Output prefix. |
| `predict_eval` | `bool` | `True` | Whether to run prediction on the eval set and output metrics. |

**Returns:** `dict` -- Evaluation metrics including accuracy, F1, confusion matrix, and ROC data.

#### `train_all_data`

```python
train_all_data(
    model_directory,
    prepared_input_data_file,
    id_class_dict_file,
    output_directory,
    output_prefix,
)
```

Train on the entire dataset (no eval split). Use after hyperparameter tuning for final model.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model_directory` | `str` or `Path` | *required* | Path to pretrained model. |
| `prepared_input_data_file` | `str` or `Path` | *required* | Prepared dataset path. |
| `id_class_dict_file` | `str` or `Path` | *required* | ID-to-class mapping path. |
| `output_directory` | `str` or `Path` | *required* | Output directory for the fine-tuned model. |
| `output_prefix` | `str` | *required* | Output prefix. |

**Returns:** `None` -- Saves fine-tuned model to `output_directory`.

#### `hyperopt_classifier`

```python
hyperopt_classifier(
    model_directory,
    num_classes,
    train_data,
    eval_data,
    output_directory,
    output_prefix,
)
```

Run hyperparameter optimization using Ray Tune.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model_directory` | `str` or `Path` | *required* | Pretrained model path. |
| `num_classes` | `int` | *required* | Number of classification classes. |
| `train_data` | `datasets.Dataset` | *required* | Training split. |
| `eval_data` | `datasets.Dataset` | *required* | Evaluation split. |
| `output_directory` | `str` or `Path` | *required* | Output directory. |
| `output_prefix` | `str` | *required* | Output prefix. |

**Returns:** `dict` -- Best hyperparameters found.

#### `train_classifier`

```python
train_classifier(
    model_directory,
    num_classes,
    train_data,
    eval_data,
    output_directory,
    output_prefix,
)
```

Train a classifier with fixed hyperparameters (no Ray Tune).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model_directory` | `str` or `Path` | *required* | Pretrained model path. |
| `num_classes` | `int` | *required* | Number of classes. |
| `train_data` | `datasets.Dataset` | *required* | Training split. |
| `eval_data` | `datasets.Dataset` | *required* | Evaluation split. |
| `output_directory` | `str` or `Path` | *required* | Output directory. |
| `output_prefix` | `str` | *required* | Output prefix. |

**Returns:** `transformers.Trainer` -- The trained Trainer object.

#### `evaluate_model`

```python
evaluate_model(
    model,
    num_classes,
    id_class_dict,
    eval_data,
    predict_trainer,
    output_directory,
    output_prefix,
)
```

Evaluate a loaded model on evaluation data.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model` | `transformers.PreTrainedModel` | *required* | Fine-tuned model. |
| `num_classes` | `int` | *required* | Number of classes. |
| `id_class_dict` | `dict` | *required* | Mapping from integer IDs to class names. |
| `eval_data` | `datasets.Dataset` | *required* | Evaluation dataset. |
| `predict_trainer` | `transformers.Trainer` | *required* | Trainer object configured for prediction. |
| `output_directory` | `str` or `Path` | *required* | Output directory. |
| `output_prefix` | `str` | *required* | Output prefix. |

**Returns:** `dict` -- Evaluation metrics.

#### `evaluate_saved_model`

```python
evaluate_saved_model(
    model_directory,
    id_class_dict_file,
    test_data_file,
    output_directory,
    output_prefix,
)
```

Load a saved fine-tuned model and evaluate on test data.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model_directory` | `str` or `Path` | *required* | Path to fine-tuned model. |
| `id_class_dict_file` | `str` or `Path` | *required* | Path to ID-to-class mapping pickle. |
| `test_data_file` | `str` or `Path` | *required* | Path to test dataset. |
| `output_directory` | `str` or `Path` | *required* | Output directory. |
| `output_prefix` | `str` | *required* | Output prefix. |

**Returns:** `dict` -- Test evaluation metrics.

#### `plot_conf_mat`

```python
plot_conf_mat(
    conf_mat_dict,
    output_directory,
    output_prefix,
    custom_class_order=None,
)
```

Plot a confusion matrix heatmap.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `conf_mat_dict` | `dict` | *required* | Confusion matrix data from evaluation. |
| `output_directory` | `str` or `Path` | *required* | Output directory. |
| `output_prefix` | `str` | *required* | Output prefix. |
| `custom_class_order` | `list` or `None` | `None` | Custom ordering of class labels for the plot axes. |

**Returns:** `None` -- Saves plot to file.

#### `plot_roc`

```python
plot_roc(
    roc_metric_dict,
    model_style_dict,
    title,
    output_directory,
    output_prefix,
)
```

Plot ROC curves.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `roc_metric_dict` | `dict` | *required* | ROC metrics from evaluation (FPR, TPR, AUC per class). |
| `model_style_dict` | `dict` | *required* | Style configuration for each model curve (color, linestyle). |
| `title` | `str` | *required* | Plot title. |
| `output_directory` | `str` or `Path` | *required* | Output directory. |
| `output_prefix` | `str` | *required* | Output prefix. |

**Returns:** `None` -- Saves ROC plot to file.

#### `plot_predictions`

```python
plot_predictions(
    predictions_file,
    id_class_dict_file,
    title,
    output_directory,
    output_prefix,
    custom_class_order=None,
)
```

Plot prediction distributions from a saved predictions file.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `predictions_file` | `str` or `Path` | *required* | Path to predictions pickle file. |
| `id_class_dict_file` | `str` or `Path` | *required* | Path to ID-to-class dict pickle. |
| `title` | `str` | *required* | Plot title. |
| `output_directory` | `str` or `Path` | *required* | Output directory. |
| `output_prefix` | `str` | *required* | Output prefix. |
| `custom_class_order` | `list` or `None` | `None` | Custom ordering of classes. |

**Returns:** `None` -- Saves prediction plot to file.

### cell_state_dict Format

```python
# Specify exact states
cell_state_dict = {
    "state_key": "disease",          # Column name in dataset
    "states": ["healthy", "diseased"]  # Classes to classify
}

# Use all unique values in the column
cell_state_dict = {
    "state_key": "cell_type",
    "states": "all"
}
```

### gene_class_dict Format

```python
gene_class_dict = {
    "Transcription_Factor": [
        "ENSG00000136997",  # MYC
        "ENSG00000136826",  # KLF4
        "ENSG00000181449",  # SOX2
    ],
    "Kinase": [
        "ENSG00000198793",  # MTOR
        "ENSG00000142208",  # AKT1
    ],
}
```

### Usage Example

```python
from geneformer import Classifier

# Cell type classifier
cc = Classifier(
    classifier="cell",
    cell_state_dict={
        "state_key": "cell_type",
        "states": "all",
    },
    filter_data={"tissue": ["heart"]},
    max_ncells=50000,
    max_ncells_per_class=5000,
    rare_threshold=100,
    freeze_layers=4,
    num_crossval_splits=5,
    forward_batch_size=100,
    model_version="V2",
    nproc=8,
    ngpu=1,
)

# Prepare data
prep_output = cc.prepare_data(
    input_data_file="/path/to/tokenized.dataset",
    output_directory="/path/to/prepared",
    output_prefix="heart_celltype",
)

# Train and evaluate with cross-validation
results = cc.validate(
    model_directory="/path/to/Geneformer-V2-104M",
    prepared_input_data_file=prep_output["prepared_data"],
    id_class_dict_file=prep_output["id_class_dict"],
    output_directory="/path/to/results",
    output_prefix="heart_celltype",
)

# Plot confusion matrix
cc.plot_conf_mat(
    conf_mat_dict=results["conf_mat"],
    output_directory="/path/to/results",
    output_prefix="heart_celltype",
)
```

---

## 5. InSilicoPerturber

Performs in silico perturbation experiments by deleting, overexpressing, inhibiting,
or activating genes and measuring the effect on cell embeddings or cell state shifts.

### Module

`geneformer.in_silico_perturber`

### Import

```python
from geneformer import InSilicoPerturber
```

### Constructor

```python
InSilicoPerturber(
    perturb_type="delete",
    perturb_rank_shift=None,
    genes_to_perturb="all",
    combos=0,
    anchor_gene=None,
    model_type="Pretrained",
    num_classes=0,
    emb_mode="cell",
    cell_emb_style="mean_pool",
    filter_data=None,
    cell_states_to_model=None,
    state_embs_dict=None,
    max_ncells=None,
    cell_inds_to_perturb="all",
    emb_layer=-1,
    model_version="V2",
    token_dictionary_file=None,
    forward_batch_size=100,
    nproc=4,
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `perturb_type` | `str` | `"delete"` | Type of perturbation. `"delete"`: remove gene tokens entirely. `"overexpress"`: move gene to rank 1 (highest expression). `"inhibit"`: shift gene to lower rank (lower expression). `"activate"`: shift gene to higher rank (higher expression). |
| `perturb_rank_shift` | `int` or `None` | `None` | Number of rank positions to shift for `"inhibit"` or `"activate"` perturbations. Positive integer. Required when `perturb_type` is `"inhibit"` or `"activate"`. Typical values: `1`, `2`, or `3`. |
| `genes_to_perturb` | `str` or `list[str]` | `"all"` | Genes to perturb. `"all"`: perturb every gene in each cell individually. `list`: list of Ensembl IDs to perturb. |
| `combos` | `int` | `0` | Combination perturbation mode. `0`: single gene perturbations only. `1`: pairwise combinations of `genes_to_perturb` (or each gene paired with `anchor_gene`). |
| `anchor_gene` | `str` or `None` | `None` | Ensembl ID of an anchor gene. When `combos=1`, each gene in `genes_to_perturb` is perturbed in combination with this anchor gene. |
| `model_type` | `str` | `"Pretrained"` | Model type. One of `"Pretrained"`, `"GeneClassifier"`, `"CellClassifier"`, `"MTLCellClassifier"`, `"Pretrained-Quantized"`, `"MTLCellClassifier-Quantized"`. |
| `num_classes` | `int` | `0` | Number of classes for classifier models. |
| `emb_mode` | `str` | `"cell"` | Embedding mode for measuring perturbation effects. `"cls"`: CLS token only. `"cell"`: mean-pooled cell embedding. `"cls_and_gene"`: CLS plus per-gene shift. `"cell_and_gene"`: cell embedding plus per-gene shift. |
| `cell_emb_style` | `str` | `"mean_pool"` | Cell embedding pooling strategy. |
| `filter_data` | `dict` or `None` | `None` | Cell filtering dictionary. |
| `cell_states_to_model` | `dict` or `None` | `None` | Cell state transition specification. See format below. |
| `state_embs_dict` | `dict` or `None` | `None` | Pre-computed state embeddings from `EmbExtractor.get_state_embs()`. Required when `cell_states_to_model` is set. |
| `max_ncells` | `int` or `None` | `None` | Maximum number of cells to perturb. |
| `cell_inds_to_perturb` | `str` or `list[int]` | `"all"` | Cell indices to perturb. `"all"` or a list of integer indices. |
| `emb_layer` | `int` | `-1` | Layer for embedding extraction. |
| `model_version` | `str` | `"V2"` | Model version. |
| `token_dictionary_file` | `str` or `Path` or `None` | `None` | Token dictionary path. |
| `forward_batch_size` | `int` | `100` | Batch size for forward passes. |
| `nproc` | `int` | `4` | Number of data loading processes. |

### Methods

#### `perturb_data`

```python
perturb_data(
    model_directory,
    input_data_file,
    output_directory,
    output_prefix,
)
```

Run in silico perturbation on all specified genes and cells.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model_directory` | `str` or `Path` | *required* | Path to model directory. |
| `input_data_file` | `str` or `Path` | *required* | Path to tokenized dataset. |
| `output_directory` | `str` or `Path` | *required* | Output directory for raw perturbation results (pickle files). |
| `output_prefix` | `str` | *required* | Output prefix. |

**Returns:** `None` -- Saves per-gene perturbation results as pickle files in `output_directory`.

### cell_states_to_model Format

Defines the cell state transition to model. The perturbation analysis measures
how each gene perturbation shifts cell embeddings toward or away from defined states.

```python
cell_states_to_model = {
    "state_key": "disease",       # Column name in the dataset
    "start_state": "dcm",         # Starting cell state (cells to perturb)
    "goal_state": "nf",           # Target state (desired direction)
    "alt_states": ["hcm"],        # Alternative states for comparison
}
```

- `state_key`: The obs column containing cell state labels.
- `start_state`: Filter cells to only those in this state for perturbation.
- `goal_state`: The desired target state. Perturbation shifts toward this state are positive.
- `alt_states`: Additional states to compute shifts toward/away from.

### state_embs_dict Format

Pre-computed mean embeddings per state, obtained from `EmbExtractor.get_state_embs()`:

```python
# Typically obtained as:
state_embs_dict = emb_extractor.get_state_embs(
    cell_states_to_model=cell_states_to_model,
    model_directory="/path/to/model",
    input_data_file="/path/to/data.dataset",
    output_directory="/path/to/output",
    output_prefix="state_embs",
)
# Returns: {"dcm": tensor, "nf": tensor, "hcm": tensor}
```

### Usage Example

```python
from geneformer import EmbExtractor, InSilicoPerturber

# Step 1: Get state embeddings
emb_extractor = EmbExtractor(
    model_type="Pretrained",
    emb_mode="cell",
    max_ncells=5000,
    model_version="V2",
)

cell_states = {
    "state_key": "disease",
    "start_state": "dcm",
    "goal_state": "nf",
    "alt_states": ["hcm"],
}

state_embs = emb_extractor.get_state_embs(
    cell_states_to_model=cell_states,
    model_directory="/path/to/Geneformer-V2-104M",
    input_data_file="/path/to/tokenized.dataset",
    output_directory="/path/to/output",
    output_prefix="state_embs",
)

# Step 2: Run perturbation
isp = InSilicoPerturber(
    perturb_type="delete",
    genes_to_perturb="all",
    combos=0,
    model_type="Pretrained",
    emb_mode="cell",
    cell_states_to_model=cell_states,
    state_embs_dict=state_embs,
    max_ncells=2000,
    emb_layer=-1,
    forward_batch_size=100,
    model_version="V2",
)

isp.perturb_data(
    model_directory="/path/to/Geneformer-V2-104M",
    input_data_file="/path/to/tokenized.dataset",
    output_directory="/path/to/perturbation_output",
    output_prefix="dcm_to_nf",
)
```

---

## 6. InSilicoPerturberStats

Computes summary statistics and rankings from raw in silico perturbation results.
Aggregates per-cell perturbation effects into gene-level scores.

### Module

`geneformer.in_silico_perturber_stats`

### Import

```python
from geneformer import InSilicoPerturberStats
```

### Constructor

```python
InSilicoPerturberStats(
    mode="goal_state_shift",
    genes_perturbed="all",
    combos=0,
    anchor_gene=None,
    cell_states_to_model=None,
    pickle_suffix="_raw.pickle",
    model_version="V2",
    token_dictionary_file=TOKEN_DICTIONARY_FILE,
    gene_name_id_dictionary_file=ENSEMBL_DICTIONARY_FILE,
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mode` | `str` | `"goal_state_shift"` | Analysis mode. `"goal_state_shift"`: rank genes by shift toward goal state. `"vs_null"`: compare perturbation effect to null distribution. `"mixture_model"`: fit mixture model to identify significant perturbations. `"aggregate_data"`: aggregate raw perturbation data across cells. `"aggregate_gene_shifts"`: aggregate per-gene embedding shifts. |
| `genes_perturbed` | `str` or `list[str]` | `"all"` | Which genes were perturbed. Must match the `genes_to_perturb` setting used in `InSilicoPerturber`. |
| `combos` | `int` | `0` | Combination mode. Must match the `combos` setting from `InSilicoPerturber`. |
| `anchor_gene` | `str` or `None` | `None` | Anchor gene Ensembl ID. Must match `InSilicoPerturber` if used. |
| `cell_states_to_model` | `dict` or `None` | `None` | Cell state specification. Must match the `InSilicoPerturber` configuration. |
| `pickle_suffix` | `str` | `"_raw.pickle"` | Suffix of the raw perturbation pickle files to load. |
| `model_version` | `str` | `"V2"` | Model version. |
| `token_dictionary_file` | `str` or `Path` | `TOKEN_DICTIONARY_FILE` | Path to token dictionary for mapping token IDs back to gene names. |
| `gene_name_id_dictionary_file` | `str` or `Path` | `ENSEMBL_DICTIONARY_FILE` | Path to gene name-to-ID dictionary for converting Ensembl IDs to gene symbols. |

### Methods

#### `get_stats`

```python
get_stats(
    input_data_directory,
    null_dist_data_directory,
    output_directory,
    output_prefix,
)
```

Compute statistics from raw perturbation results.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_data_directory` | `str` or `Path` | *required* | Directory containing raw perturbation pickle files from `InSilicoPerturber.perturb_data`. |
| `null_dist_data_directory` | `str` or `Path` or `None` | *required* | Directory containing null distribution data (from a control perturbation run). Required for `"vs_null"` mode; can be `None` for other modes. |
| `output_directory` | `str` or `Path` | *required* | Output directory for statistics results. |
| `output_prefix` | `str` | *required* | Output prefix. |

**Returns:** `pandas.DataFrame` -- DataFrame with gene-level perturbation statistics, ranked by effect size. Columns depend on the analysis mode.

### Usage Example

```python
from geneformer import InSilicoPerturberStats

cell_states = {
    "state_key": "disease",
    "start_state": "dcm",
    "goal_state": "nf",
    "alt_states": ["hcm"],
}

# Compute goal state shift rankings
ispstats = InSilicoPerturberStats(
    mode="goal_state_shift",
    genes_perturbed="all",
    combos=0,
    cell_states_to_model=cell_states,
    model_version="V2",
)

stats_df = ispstats.get_stats(
    input_data_directory="/path/to/perturbation_output",
    null_dist_data_directory=None,
    output_directory="/path/to/stats_output",
    output_prefix="dcm_to_nf_stats",
)

# Top genes shifting cells toward goal state
print(stats_df.head(20))
```

---

## 7. MTLClassifier

Multi-task learning classifier that simultaneously predicts multiple cell state
labels. Uses Optuna for hyperparameter optimization and supports distributed training.

### Module

`geneformer.mtl_classifier`

### Import

```python
from geneformer import MTLClassifier
```

### Constructor

```python
MTLClassifier(
    task_columns,
    train_path,
    val_path,
    test_path,
    pretrained_path,
    model_save_path,
    results_dir,
    trials_result_path,
    batch_size=4,
    n_trials=15,
    study_name="mtl_study",
    max_layers_to_freeze=4,
    epochs=10,
    tensorboard_log_dir=None,
    distributed_training=False,
    master_addr="localhost",
    master_port="12355",
    use_attention_pooling=False,
    use_task_weights=False,
    hyperparameters=None,
    manual_hyperparameters=None,
    use_manual_hyperparameters=False,
    use_wandb=False,
    wandb_project=None,
    gradient_clipping=False,
    max_grad_norm=1.0,
    seed=42,
    gradient_accumulation_steps=1,
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `task_columns` | `list[str]` | *required* | List of column names in the dataset to use as classification targets. Each column becomes a separate classification head. |
| `train_path` | `str` or `Path` | *required* | Path to the training dataset (HuggingFace `.dataset` directory). |
| `val_path` | `str` or `Path` | *required* | Path to the validation dataset. |
| `test_path` | `str` or `Path` | *required* | Path to the test dataset. |
| `pretrained_path` | `str` or `Path` | *required* | Path to the pretrained Geneformer model directory. |
| `model_save_path` | `str` or `Path` | *required* | Directory to save the fine-tuned MTL model. |
| `results_dir` | `str` or `Path` | *required* | Directory for evaluation results and metrics. |
| `trials_result_path` | `str` or `Path` | *required* | Path to save Optuna trial results. |
| `batch_size` | `int` | `4` | Training batch size per device. |
| `n_trials` | `int` | `15` | Number of Optuna hyperparameter optimization trials. |
| `study_name` | `str` | `"mtl_study"` | Name for the Optuna study. |
| `max_layers_to_freeze` | `int` | `4` | Maximum number of layers Optuna can freeze during search. |
| `epochs` | `int` | `10` | Number of training epochs. |
| `tensorboard_log_dir` | `str` or `Path` or `None` | `None` | Directory for TensorBoard logs. `None` disables TensorBoard logging. |
| `distributed_training` | `bool` | `False` | Whether to use distributed data parallel training across multiple GPUs. |
| `master_addr` | `str` | `"localhost"` | Master address for distributed training. |
| `master_port` | `str` | `"12355"` | Master port for distributed training. |
| `use_attention_pooling` | `bool` | `False` | Whether to use attention-based pooling instead of mean pooling for the CLS representation. |
| `use_task_weights` | `bool` | `False` | Whether to learn per-task loss weights during training. |
| `hyperparameters` | `dict` or `None` | `None` | Custom hyperparameter search space for Optuna. |
| `manual_hyperparameters` | `dict` or `None` | `None` | Fixed hyperparameters to use instead of Optuna search. |
| `use_manual_hyperparameters` | `bool` | `False` | If `True`, use `manual_hyperparameters` instead of running Optuna. |
| `use_wandb` | `bool` | `False` | Whether to log to Weights & Biases. |
| `wandb_project` | `str` or `None` | `None` | W&B project name. |
| `gradient_clipping` | `bool` | `False` | Whether to apply gradient clipping. |
| `max_grad_norm` | `float` | `1.0` | Maximum gradient norm when clipping is enabled. |
| `seed` | `int` | `42` | Random seed for reproducibility. |
| `gradient_accumulation_steps` | `int` | `1` | Number of gradient accumulation steps before weight update. Effective batch size = `batch_size * gradient_accumulation_steps`. |

### Methods

#### `run_optuna_study`

```python
run_optuna_study()
```

Run the full Optuna hyperparameter search, train with the best parameters, and save the model.

**Returns:** `None` -- Saves the best model to `model_save_path` and trial results to `trials_result_path`.

#### `load_and_evaluate_test_model`

```python
load_and_evaluate_test_model()
```

Load the saved best model and evaluate on the test set.

**Returns:** `dict` -- Per-task evaluation metrics including accuracy, F1, and confusion matrices.

### Usage Example

```python
from geneformer import MTLClassifier

mtl = MTLClassifier(
    task_columns=["cell_type", "disease", "tissue"],
    train_path="/path/to/train.dataset",
    val_path="/path/to/val.dataset",
    test_path="/path/to/test.dataset",
    pretrained_path="/path/to/Geneformer-V2-104M",
    model_save_path="/path/to/mtl_model",
    results_dir="/path/to/mtl_results",
    trials_result_path="/path/to/mtl_results/trials.pkl",
    batch_size=8,
    n_trials=20,
    study_name="cell_mtl",
    max_layers_to_freeze=6,
    epochs=15,
    seed=42,
)

# Run hyperparameter optimization and train
mtl.run_optuna_study()

# Evaluate on test set
test_results = mtl.load_and_evaluate_test_model()
for task, metrics in test_results.items():
    print(f"{task}: accuracy={metrics['accuracy']:.3f}, f1={metrics['f1']:.3f}")
```

---

## 8. GeneformerPretrainer

Custom Trainer for pretraining Geneformer with masked language modeling (MLM)
on tokenized single-cell data. Extends HuggingFace `Trainer` with custom
data collation for variable-length gene token sequences.

### Module

`geneformer.pretrainer`

### Import

```python
from geneformer.pretrainer import GeneformerPretrainer
```

### Constructor

Extends `transformers.Trainer`. Accepts all standard `Trainer` arguments plus:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `token_dictionary` | `dict` | *required* | Token dictionary mapping gene Ensembl IDs to token IDs. Loaded from the pickle file. |
| `mlm` | `bool` | `True` | Whether to perform masked language modeling. |
| `mlm_probability` | `float` | `0.15` | Probability of masking each gene token during training. |
| `example_lengths_file` | `str` or `Path` | *required* | Path to a pickle file containing the length of each example in the training dataset. Used for efficient batching by length. |

### Internal Components

- **GeneformerPreCollator**: Pads variable-length token sequences within a batch to the same length.
- **DataCollatorForLanguageModeling**: Applies MLM masking to padded batches.

### Usage Example

```python
from transformers import AutoModelForMaskedLM, TrainingArguments
from datasets import load_from_disk
from geneformer.pretrainer import GeneformerPretrainer
import pickle

# Load token dictionary
with open("/path/to/token_dictionary_gc104M.pkl", "rb") as f:
    token_dict = pickle.load(f)

# Load tokenized training data
train_dataset = load_from_disk("/path/to/tokenized_train.dataset")

# Configure training
training_args = TrainingArguments(
    output_dir="/path/to/pretrain_output",
    per_device_train_batch_size=32,
    num_train_epochs=3,
    learning_rate=1e-4,
    warmup_steps=10000,
    weight_decay=0.01,
    save_steps=5000,
    logging_steps=100,
    fp16=True,
)

# Load model
model = AutoModelForMaskedLM.from_pretrained("/path/to/Geneformer-V2-104M")

# Create trainer
trainer = GeneformerPretrainer(
    model=model,
    args=training_args,
    train_dataset=train_dataset,
    token_dictionary=token_dict,
    mlm=True,
    mlm_probability=0.15,
    example_lengths_file="/path/to/example_lengths.pkl",
)

trainer.train()
trainer.save_model("/path/to/continued_pretrain_model")
```

---

## 9. DataCollators

Custom data collators for fine-tuning Geneformer on classification tasks.
Handle padding and label preparation for variable-length gene token sequences.

### Module

`geneformer.collator_for_classification`

### DataCollatorForCellClassification

```python
from geneformer import DataCollatorForCellClassification
```

Data collator for cell-level classification tasks. Pads `input_ids` to batch max
length and prepares `labels` for the classification head.

| Attribute | Type | Description |
|-----------|------|-------------|
| `padding` | `bool` or `str` | Padding strategy. Default: `True` (pad to longest in batch). |
| `max_length` | `int` or `None` | Maximum sequence length. `None` uses model's max. |
| `pad_to_multiple_of` | `int` or `None` | Pad to multiple of this value. |
| `return_tensors` | `str` | Return format. Default: `"pt"` (PyTorch). |

### DataCollatorForGeneClassification

```python
from geneformer import DataCollatorForGeneClassification
```

Data collator for gene-level classification tasks. Pads sequences and prepares
per-gene labels.

| Attribute | Type | Description |
|-----------|------|-------------|
| `padding` | `bool` or `str` | Padding strategy. Default: `True`. |
| `max_length` | `int` or `None` | Maximum sequence length. |
| `pad_to_multiple_of` | `int` or `None` | Pad to multiple of this value. |
| `return_tensors` | `str` | Return format. Default: `"pt"`. |

### Usage Example

```python
from geneformer import DataCollatorForCellClassification
from transformers import Trainer, TrainingArguments

collator = DataCollatorForCellClassification()

training_args = TrainingArguments(
    output_dir="/path/to/output",
    per_device_train_batch_size=16,
    num_train_epochs=5,
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=train_dataset,
    eval_dataset=eval_dataset,
    data_collator=collator,
)

trainer.train()
```

---

## 10. Utility Functions

### perturber_utils

Module: `geneformer.perturber_utils`

| Function | Signature | Description |
|----------|-----------|-------------|
| `load_model` | `load_model(model_type, num_classes, model_directory, mode="eval", quantize=False)` | Load a Geneformer model for perturbation or embedding extraction. Returns the model on the appropriate device. |
| `load_and_filter` | `load_and_filter(filter_data, nproc, input_data_file)` | Load a tokenized dataset and apply filtering. Returns filtered `datasets.Dataset`. |
| `filter_by_dict` | `filter_by_dict(dataset, filter_dict, nproc)` | Filter a dataset by a dictionary of column-value pairs. |
| `pad_tensor_list` | `pad_tensor_list(tensor_list, dynamic_or_constant, pad_token_id, model_input_size)` | Pad a list of tensors to uniform length. `"dynamic"` pads to batch max; `"constant"` pads to `model_input_size`. |
| `gen_attention_mask` | `gen_attention_mask(padded_input, pad_token_id)` | Generate attention mask from padded input IDs. Returns tensor of 1s (real tokens) and 0s (padding). |
| `get_model_input_size` | `get_model_input_size(model)` | Return the model's maximum input size from its config. |
| `get_model_emb_dims` | `get_model_emb_dims(model)` | Return the model's hidden dimension size. |
| `mean_nonpadding_embs` | `mean_nonpadding_embs(embs, original_lens, dim=1)` | Compute mean embeddings excluding padding positions. |
| `quant_layers` | `quant_layers(model)` | Return the number of quantizable transformer layers. |
| `move_to_cuda` | `move_to_cuda(batch)` | Move a batch dict of tensors to CUDA device. |
| `validate_cell_states_to_model` | `validate_cell_states_to_model(cell_states_to_model)` | Validate the format of a `cell_states_to_model` dict. Raises `ValueError` on invalid format. |
| `flatten_list` | `flatten_list(nested_list)` | Flatten a nested list of lists into a single list. |
| `downsample_and_sort` | `downsample_and_sort(dataset, max_ncells)` | Downsample to `max_ncells` and sort by sequence length for efficient batching. |

### classifier_utils

Module: `geneformer.classifier_utils`

| Function | Signature | Description |
|----------|-----------|-------------|
| `downsample_and_shuffle` | `downsample_and_shuffle(dataset, max_ncells, seed=42)` | Downsample a dataset and shuffle. |
| `subsample_by_class` | `subsample_by_class(dataset, label_col, max_per_class, seed=42)` | Subsample each class to at most `max_per_class` examples. |
| `remove_rare` | `remove_rare(dataset, label_col, threshold)` | Remove classes with fewer than `threshold` examples. |
| `label_classes` | `label_classes(dataset, cell_state_dict)` | Add integer class labels to dataset based on `cell_state_dict`. |
| `label_gene_classes` | `label_gene_classes(dataset, gene_class_dict, token_dictionary)` | Add per-gene class labels for gene classification. |
| `prep_gene_classifier_train_eval_split` | `prep_gene_classifier_train_eval_split(data, targets, train_index, eval_index)` | Split gene classification data into train/eval sets. |
| `compute_metrics` | `compute_metrics(eval_pred)` | Compute accuracy, macro F1, and weighted F1 from evaluation predictions. |
| `get_default_train_args` | `get_default_train_args(output_dir, epochs=10)` | Return default `TrainingArguments` dict for cell classification. |

### evaluation_utils

Module: `geneformer.evaluation_utils`

| Function | Signature | Description |
|----------|-----------|-------------|
| `classifier_predict` | `classifier_predict(model, eval_dataset, collator, forward_batch_size)` | Run prediction on evaluation data and return logits and labels. |
| `preprocess_classifier_batch` | `preprocess_classifier_batch(batch, max_len)` | Preprocess a batch for classifier prediction (pad, create attention mask). |
| `compute_metrics` | `compute_metrics(eval_pred)` | Same as `classifier_utils.compute_metrics`. Compute accuracy and F1. |

---

## 11. Constants

Module: `geneformer.tokenizer` (V2 defaults, gc104M)

| Constant | Description | Filename |
|----------|-------------|----------|
| `GENE_MEDIAN_FILE` | Path to gene median dictionary for rank-value encoding normalization. | `gene_median_dictionary_gc104M.pkl` |
| `TOKEN_DICTIONARY_FILE` | Path to token dictionary mapping Ensembl IDs to integer token IDs. | `token_dictionary_gc104M.pkl` |
| `ENSEMBL_DICTIONARY_FILE` | Path to gene name-to-Ensembl-ID dictionary for symbol lookup. | `gene_name_id_dict_gc104M.pkl` |
| `ENSEMBL_MAPPING_FILE` | Path to Ensembl mapping dictionary for gene ID versioning and cross-referencing. | `ensembl_mapping_dict_gc104M.pkl` |

These files are bundled with the `geneformer` pip package and are resolved automatically
when using default parameter values. Override paths only when using custom dictionaries.

```python
# Accessing default paths
from geneformer.tokenizer import (
    GENE_MEDIAN_FILE,
    TOKEN_DICTIONARY_FILE,
    ENSEMBL_MAPPING_FILE,
)
print(GENE_MEDIAN_FILE)
# /path/to/site-packages/geneformer/gene_median_dictionary_gc104M.pkl
```

---

## Appendix: Common Patterns

### End-to-End Workflow: Tokenize, Embed, Classify

```python
from geneformer import TranscriptomeTokenizer, EmbExtractor, Classifier

# 1. Tokenize
tk = TranscriptomeTokenizer(
    custom_attr_name_dict={"cell_type": "cell_type"},
    nproc=8,
    model_version="V2",
)
dataset = tk.tokenize_data(
    data_directory="/data/h5ad",
    output_directory="/data/tokenized",
    output_prefix="my_data",
    file_format="h5ad",
)

# 2. Extract embeddings (optional, for exploration)
ee = EmbExtractor(
    emb_mode="cls",
    max_ncells=10000,
    model_version="V2",
)
embs = ee.extract_embs(
    model_directory="/models/Geneformer-V2-104M",
    input_data_file="/data/tokenized/my_data.dataset",
    output_directory="/results/embs",
    output_prefix="my_embs",
)
ee.plot_embs(embs, "umap", "/results/embs", "my_embs")

# 3. Fine-tune classifier
clf = Classifier(
    classifier="cell",
    cell_state_dict={"state_key": "cell_type", "states": "all"},
    freeze_layers=4,
    num_crossval_splits=5,
    model_version="V2",
)
prep = clf.prepare_data(
    input_data_file="/data/tokenized/my_data.dataset",
    output_directory="/data/prepared",
    output_prefix="my_clf",
)
results = clf.validate(
    model_directory="/models/Geneformer-V2-104M",
    prepared_input_data_file=prep["prepared_data"],
    id_class_dict_file=prep["id_class_dict"],
    output_directory="/results/clf",
    output_prefix="my_clf",
)
```

### End-to-End Workflow: In Silico Perturbation

```python
from geneformer import EmbExtractor, InSilicoPerturber, InSilicoPerturberStats

# 1. Compute state embeddings
ee = EmbExtractor(
    emb_mode="cell",
    max_ncells=5000,
    model_version="V2",
)
cell_states = {
    "state_key": "disease",
    "start_state": "dcm",
    "goal_state": "nf",
    "alt_states": ["hcm"],
}
state_embs = ee.get_state_embs(
    cell_states_to_model=cell_states,
    model_directory="/models/Geneformer-V2-104M",
    input_data_file="/data/tokenized/cardiac.dataset",
    output_directory="/results/state_embs",
    output_prefix="cardiac",
)

# 2. Run perturbations
isp = InSilicoPerturber(
    perturb_type="delete",
    genes_to_perturb="all",
    model_type="Pretrained",
    emb_mode="cell",
    cell_states_to_model=cell_states,
    state_embs_dict=state_embs,
    max_ncells=2000,
    model_version="V2",
)
isp.perturb_data(
    model_directory="/models/Geneformer-V2-104M",
    input_data_file="/data/tokenized/cardiac.dataset",
    output_directory="/results/perturbations",
    output_prefix="cardiac_perturb",
)

# 3. Compute statistics
stats = InSilicoPerturberStats(
    mode="goal_state_shift",
    genes_perturbed="all",
    cell_states_to_model=cell_states,
    model_version="V2",
)
stats_df = stats.get_stats(
    input_data_directory="/results/perturbations",
    null_dist_data_directory=None,
    output_directory="/results/stats",
    output_prefix="cardiac_stats",
)
print(stats_df.head(20))
```

### Quantized Training (4-bit QLoRA)

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
    freeze_layers=0,
    model_version="V2",
)
```
