# SCimilarity Troubleshooting

## Common Errors

### Gene Overlap Too Low

**Error**: `align_dataset` raises error about insufficient gene overlap.

**Cause**: Your data has fewer than 5,000 genes in common with the model's `gene_order.tsv`.

**Fix**:
1. Ensure gene names are HGNC symbols (not Ensembl IDs): `convert_id2symbol(adata, mapping_table)`
2. Check for duplicate gene symbols: `consolidate_duplicate_symbols(adata)`
3. Verify: `len(set(adata.var.index) & set(model.gene_order))` should be >5000
4. If using mouse data, scimilarity is trained on human genes — map orthologs first

### NaN in Embeddings

**Error**: `RuntimeError: NaN detected in embeddings`

**Cause**: Input expression matrix contains NaN values.

**Fix**:
```python
import numpy as np
from scipy.sparse import issparse
if issparse(adata.X):
    assert not np.isnan(adata.X.data).any()
else:
    assert not np.isnan(adata.X).any()
```

### Missing layers["counts"]

**Error**: KeyError when calling `lognorm_counts()` or `filter_cells()`.

**Cause**: Raw counts not stored in `adata.layers["counts"]`.

**Fix**:
```python
adata.layers["counts"] = adata.X.copy()  # Before any normalization
```

### kNN Index Not Found

**Error**: `Warning: No KNN index found at /path/to/knn.bin`

**Cause**: kNN index file missing or wrong path.

**Fix**:
1. Verify the file exists: `ls model_path/annotation/labelled_kNN.bin`
2. For CellQuery: `ls model_path/cellsearch/full_kNN.bin`
3. If using custom filenames: pass `filenames={"knn": "your_file.bin"}`

### Negative Values in Expression

**Cause**: Data was already processed (e.g., scaled/centered).

**Fix**: Use raw counts from `layers["counts"]`, not transformed `.X`.

### Out of Memory

**Cause**: Too many cells being processed at once.

**Fix**:
- Reduce `buffer_size` in `get_embeddings()` (default 10000)
- For exhaustive search, reduce `buffer_size` parameter (default 100000)
- Process in chunks and concatenate embeddings

### TileDB Errors with CellQuery

**Cause**: TileDB arrays not present or corrupted.

**Fix**:
1. Verify TileDB directories exist: `cell_metadata/` and `cell_embedding/`
2. Check TileDB version compatibility: requires `tiledb >= 0.18.2`
3. If only using kNN search (not exhaustive), set `load_knn=True` and avoid TileDB-dependent methods

### GPU Memory Errors

**Fix**:
```python
# Use CPU for embedding, GPU only when needed
ce = CellEmbedding(model_path="/path/to/model", use_gpu=False)
# Reduce buffer_size for GPU
embeddings = ce.get_embeddings(adata.X, buffer_size=5000)
```

## Performance Tuning

| Parameter | Default | Increase For | Decrease For |
|-----------|---------|-------------|--------------|
| `buffer_size` (embeddings) | 10000 | GPU with large VRAM | Limited memory |
| `k` (neighbors) | 50 (annotation), 10000 (search) | More comprehensive results | Speed |
| `ef` (hnswlib) | 100 (annotation), k (search) | Better recall accuracy | Speed |
| `max_dist` | None | Filtering to close matches | N/A |
| `buffer_size` (exhaustive) | 100000 | More memory available | Limited memory |

## Compatibility Notes

- **Python**: >= 3.10
- **PyTorch**: >= 1.10.1
- **Normalization**: CPM 1e4 + log1p (same as scGPT)
- **Expression scale**: Both scimilarity and scGPT use identical normalization, so outputs are directly compatible
- **Species**: Trained on human data; for mouse, map orthologs to human genes first
