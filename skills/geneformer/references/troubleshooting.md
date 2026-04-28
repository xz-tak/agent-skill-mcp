# Geneformer Troubleshooting Guide

## Common Errors

### Error: Missing ensembl_id Column
**Symptom**: KeyError on `ensembl_id` or empty tokenization output with no genes detected.
**Cause**: The gene attribute `ensembl_id` is not found in the `.loom` or `.h5ad` file. Geneformer requires Ensembl gene IDs to map genes to its token dictionary.
**Fix**:
```python
# If your .var index already contains Ensembl IDs (e.g., ENSG00000141510)
adata.var["ensembl_id"] = adata.var.index

# If Ensembl IDs are in a different column
adata.var["ensembl_id"] = adata.var["gene_ids"]

# For .loom files, set the row attribute during creation
loompy.create("output.loom", matrix, row_attrs={"ensembl_id": ensembl_ids}, col_attrs=col_attrs)
```

### Error: Missing n_counts Attribute
**Symptom**: KeyError on `n_counts` during tokenization.
**Cause**: The cell attribute `n_counts` (total UMI counts per cell) is not present in the data. The tokenizer requires this to rank-order genes by expression.
**Fix**:
```python
import numpy as np

# For dense matrices
adata.obs["n_counts"] = np.array(adata.X.sum(axis=1)).flatten()

# For sparse matrices
adata.obs["n_counts"] = np.array(adata.X.sum(axis=1)).flatten()

# Verify
assert "n_counts" in adata.obs.columns
```

### Error: CLS Token Assertion Error
**Symptom**: `AssertionError: "First token is not <cls> token value"` when extracting embeddings.
**Cause**: Using `emb_mode="cls"` with data that was tokenized without `special_token=True`. CLS embedding requires the `<cls>` token to be prepended during tokenization.
**Fix**:
```python
# Re-tokenize with special_token=True (this is the V2 default)
tokenizer = TranscriptomeTokenizer(
    custom_attr_name_dict={"cell_type": "cell_type"},
    special_token=True,  # Required for CLS mode
)
tokenizer.tokenize_data("input_dir", "output_dir", "tokenized_data")

# Alternatively, use emb_mode="cell" which does not require CLS token
embex = EmbExtractor(
    model_type="Pretrained",
    emb_mode="cell",  # Mean pooling, no CLS needed
)
```

### Error: CUDA Out of Memory
**Symptom**: `RuntimeError: CUDA out of memory. Tried to allocate X MiB`.
**Cause**: The model or batch size exceeds available GPU memory. V2-316M requires significantly more VRAM than V2-104M.
**Fix**:
```python
# Option 1: Reduce forward batch size
embex = EmbExtractor(
    model_type="Pretrained",
    forward_batch_size=50,  # Default is 200; try 50-100
)

# Option 2: Use the smaller model
model_directory = "path/to/geneformer-v2-104M"

# Option 3: Enable quantization (requires bitsandbytes + peft)
embex = EmbExtractor(
    model_type="Pretrained",
    quantize=True,  # 8-bit quantization
)

# Option 4: Limit cell count
embex = EmbExtractor(
    model_type="Pretrained",
    max_ncells=1000,
)
```

### Error: Empty Tokenized Dataset
**Symptom**: Tokenized dataset has 0 rows or 0 genes after tokenization.
**Cause**: No genes in the input data overlap with Geneformer's token dictionary, or all cells were filtered out by quality criteria (e.g., `filter_pass` column).
**Fix**:
```python
# Check that ensembl_id values are valid Ensembl IDs
print(adata.var["ensembl_id"].head())
# Should look like: ENSG00000141510, ENSG00000012048, etc.

# Check overlap with token dictionary
from geneformer import TranscriptomeTokenizer
tk = TranscriptomeTokenizer()
token_dict = tk.gene_token_dict
overlap = set(adata.var["ensembl_id"]) & set(token_dict.keys())
print(f"Overlap: {len(overlap)} / {len(adata.var)} genes")

# If filter_pass column exists, check how many cells pass
if "filter_pass" in adata.obs.columns:
    print(f"Cells passing filter: {adata.obs['filter_pass'].sum()} / {len(adata.obs)}")
```

### Error: Duplicate Ensembl IDs
**Symptom**: Warning about duplicate gene IDs during tokenization; unexpected gene count in output.
**Cause**: Multiple gene entries in the input map to the same Ensembl ID (common after gene name to ID conversion).
**Fix**:
```python
# Use collapse_gene_ids=True (this is the default behavior)
tokenizer = TranscriptomeTokenizer(
    custom_attr_name_dict={"cell_type": "cell_type"},
    collapse_gene_ids=True,  # Sums expression for duplicate IDs
)

# Or manually deduplicate before tokenization
adata.var_names_make_unique()
adata = adata[:, ~adata.var["ensembl_id"].duplicated(keep="first")]
```

### Error: Dataset Format Errors
**Symptom**: `FileNotFoundError` or `"Not a valid dataset"` error when loading tokenized data.
**Cause**: Attempting to load a `.dataset` directory with `scanpy.read_h5ad()`, or loading an `.h5ad` file with `datasets.load_from_disk()`. Tokenized output is in HuggingFace Dataset format, not AnnData.
**Fix**:
```python
# WRONG: Tokenized data is NOT an h5ad file
# adata = sc.read_h5ad("tokenized.dataset")

# CORRECT: Use HuggingFace datasets
from datasets import load_from_disk
tokenized_dataset = load_from_disk("tokenized.dataset")

# To inspect the dataset
print(tokenized_dataset)
print(tokenized_dataset[0])  # First cell's token sequence
```

### Error: Pickle Suffix Mismatch
**Symptom**: `InSilicoPerturberStats.get_stats()` finds no files or returns empty results.
**Cause**: The `pickle_suffix` parameter does not match the actual suffix of output files generated by `InSilicoPerturber`.
**Fix**:
```python
import os

# Check actual file names in the ISP output directory
isp_output_dir = "path/to/isp_output/"
files = os.listdir(isp_output_dir)
print([f for f in files if f.endswith(".pickle")])
# Example output: ['cell_type_raw.pickle', 'cell_type_raw.pickle']

# Match the suffix exactly (default is "_raw.pickle")
ispstats = InSilicoPerturberStats(
    pickle_suffix="_raw.pickle",  # Must match actual file suffix
)
```

### Error: Quantization Dependencies Missing
**Symptom**: `ImportError: No module named 'bitsandbytes'` or `ImportError: No module named 'peft'`.
**Cause**: Quantization requires additional packages that are not part of the base Geneformer installation.
**Fix**:
```bash
pip install bitsandbytes peft

# For CUDA compatibility issues with bitsandbytes
pip install bitsandbytes --prefer-binary

# Verify installation
python -c "import bitsandbytes; import peft; print('OK')"
```

### Error: MTL Missing unique_cell_id
**Symptom**: `KeyError: "unique_cell_id"` during `MTLClassifier` training or evaluation.
**Cause**: The tokenized dataset lacks the required `unique_cell_id` column that MTLClassifier uses for cell-level tracking.
**Fix**:
```python
from datasets import load_from_disk

dataset = load_from_disk("tokenized.dataset")

# Add unique_cell_id column
dataset = dataset.map(
    lambda x, idx: {"unique_cell_id": str(idx)},
    with_indices=True,
)

# Save updated dataset
dataset.save_to_disk("tokenized_with_ids.dataset")
```

### Error: S3 Access Errors
**Symptom**: `"An error occurred (AccessDenied)"` or `NoCredentialsError` when accessing model files on S3.
**Cause**: AWS credentials are not configured, have expired, or lack permissions for the target S3 bucket.
**Fix**:
```bash
# Configure credentials
aws configure
# Enter Access Key ID, Secret Access Key, region (us-east-1), output format (json)

# Verify access
aws s3 ls s3://tec-rnd-sci-dev-gi2/ --region us-east-1

# If using temporary credentials, check expiration
aws sts get-caller-identity

# For SageMaker notebooks, ensure the execution role has S3 access
# Check IAM policy includes: s3:GetObject, s3:ListBucket for the bucket
```

---

## Performance Tuning

| Parameter | Low Memory | Balanced | High Performance |
|-----------|-----------|----------|-----------------|
| forward_batch_size | 50 | 200 | 500 |
| max_ncells | 1000 | 5000 | None (all) |
| model tier | V2-104M | V2-104M | V2-316M |
| quantize | True | False | False |
| nproc | 4 | 8 | 16 |
