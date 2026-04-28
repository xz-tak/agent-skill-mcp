# PopV Troubleshooting Guide

## Common Issues and Solutions

### 1. Gene ID Mismatch

**Symptom:** Low percentage of reference vars found in query data (e.g., "Found 10% reference vars in query data").

**Cause:** Query data uses gene symbols (e.g., `TP53`) while the reference expects Ensembl IDs (e.g., `ENSG00000141510`).

**Solution:** Map gene symbols to Ensembl IDs before running popV:
```python
import mygene
mg = mygene.MyGeneInfo()
geneid_query = mg.querymany(query_adata.var_names, scopes='symbol', fields='ensembl.gene', species='human')
geneidmap_dict = {
    res["query"]: res["ensembl"]["gene"]
    if "ensembl" in res and "gene" in res["ensembl"]
    else "not_mapped"
    for res in geneid_query
}
query_adata.var["gene_ids"] = query_adata.var.index.map(geneidmap_dict)
query_adata.var_names = query_adata.var["gene_ids"]
query_adata = query_adata[:, ~query_adata.var_names.duplicated()]
```

**Expected match rate:** Typically 70-90% of reference vars should be found in query data. Below 50% suggests a serious mismatch.

---

### 2. Cell Count Changes in Inference Mode

**Symptom:** `AssertionError: Number of cells changed during preprocessing` when using `prediction_mode="inference"`.

**Cause:** In inference mode, popV creates a joint embedding of query + reference cells. The output `adata` contains both. After annotation, you need to filter back to query cells only.

**Solution:**
```python
if PREDICTION_MODEL != "fast":
    adata = adata[adata.obs_names.isin(query_adata.obs_names), :]
```

Then when copying popv columns back to the original adata, handle the size mismatch:
```python
if query_adata2.n_obs == adata.n_obs:
    query_adata2.obs[popv_cols] = adata.obs[popv_cols].loc[query_adata2.obs.index]
else:
    for col in popv_cols:
        query_adata2.obs[col] = query_adata2.obs.index.map(adata.obs[col].to_dict())
        if col.endswith(('score', 'probabilities')):
            query_adata2.obs[col] = query_adata2.obs[col].fillna(0)
        elif not col.endswith('depth'):
            query_adata2.obs[col] = query_adata2.obs[col].fillna('unknown_celltype_label')
```

---

### 3. OnClass Seen != Prediction

**Symptom:** `popv_onclass_seen` and `popv_onclass_prediction` columns differ for many cells.

**Explanation:** This is expected behavior in inference mode. OnClass can predict "unseen" cell types — types not in the training data but present in the Cell Ontology. `popv_onclass_seen` is the closest training label; `popv_onclass_prediction` may be a more specific ontology term.

In **fast mode**, these should be 100% identical. If they differ in fast mode, it indicates a configuration issue.

**Action:** If they match 100%, clean up the redundant column:
```python
if count_same == total:
    adata.obs.drop(columns=["popv_onclass_seen", "popv_onclass_seen_probabilities"],
                   inplace=True, errors="ignore")
```

---

### 4. CUDA / GPU Warnings

**Symptom:** Messages like `Unable to register cuFFT factory`, `Can't find libdevice directory`.

**Impact:** These are warnings, not errors. PopV will fall back to CPU for operations that fail on GPU.

**Solution:** Ignore unless performance is critically slow. To suppress:
```python
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # Force CPU-only
```

Or to use specific GPUs:
```python
os.environ["CUDA_VISIBLE_DEVICES"] = "0,1"  # Use GPUs 0 and 1
```

---

### 5. Raw Counts Not in X

**Symptom:** Errors during preprocessing, or `popV` produces nonsensical results.

**Cause:** `query_adata.X` contains normalized/log-transformed data instead of raw integer counts.

**Solution:** Copy raw counts from the appropriate layer:
```python
# Common layer names for raw counts:
query_adata.X = query_adata.layers["raw_counts"].copy()  # or "counts", "spliced"

# For CellxGene format:
query_adata.X = query_adata.raw.X
```

Verify with:
```python
import numpy as np
print(query_adata.X[:5, :5].toarray())  # Should be non-negative integers
```

---

### 6. Ontology Mapping Failures

**Symptom:** `Unmatched values in column 'popv_*_ontology_name': ['unknown_celltype_label']`

**Impact:** This is expected — `unknown_celltype_label` is a placeholder, not a real Cell Ontology term.

**Action:** No action needed. These cells couldn't be confidently assigned. The `_ontology_id` column will be `NaN` for these cells.

---

### 7. Memory Issues with Large Datasets

**Symptom:** Out of memory errors during `Process_Query` or `annotate_data`.

**Solutions:**
1. Use `prediction_mode="fast"` instead of `"inference"` — processes only query cells
2. Reduce `n_samples_per_label` (e.g., from 300 to 100)
3. Subset the query dataset into batches and annotate separately
4. Use a tissue-specific model instead of `All_Cells`

---

### 8. Precomputed PCs Conflict

**Symptom:** Error during preprocessing about PCA dimensions mismatch.

**Cause:** Precomputed principal components in `query_adata.varm["PCs"]` conflict with popV's recomputation.

**Solution:** Remove before running popV:
```python
if "PCs" in query_adata.varm.keys():
    del query_adata.varm["PCs"]
```

---

### 9. Batch Key Issues

**Symptom:** Poor integration or errors about missing `batch_key`.

**Cause:** The batch key column doesn't exist or has too many unique values.

**Solution:** Create a composite batch key from available metadata:
```python
batch_list = ["donor_id", "assay", "tissue"]  # Adjust to available columns
query_adata.obs['batch_key'] = query_adata.obs.apply(
    lambda row: "_".join([str(row[ele]) for ele in batch_list if ele in row.index]),
    axis=1
)
```

If no batch information exists, use a constant:
```python
query_adata.obs['batch_key'] = "single_batch"
```

---

### 10. Model Not Found

**Symptom:** `FileNotFoundError` when loading reference adata or scVI model.

**Checklist:**
1. Verify `MODEL_DIR` points to the correct location
2. Check that the specific model was downloaded: `ls {MODEL_DIR}/{SOURCE}/{COLLECTION}/{MODEL}/`
3. For HuggingFace models, the reference file is `minified_ref_adata.h5ad`
4. For custom models, the reference file is `{MODEL}.h5ad`
5. Verify ontology files exist: `ls {MODEL_DIR}/{SOURCE}/ontology/`

---

### 11. Prediction Mode Selection

| Mode | Behavior | Speed | Use When |
|------|----------|-------|----------|
| `"fast"` | Uses only pretrained classifiers on query cells | Fast | Quick annotation, limited compute |
| `"inference"` | Retrains on joint embedding of query + reference | Slow | Best accuracy, deeper integration |

The `"inference"` mode produces a joint embedding which is useful for visualization but takes significantly longer. For most use cases, `"fast"` mode provides good results.

---

### 12. Cross-Species Annotation Issues

**Symptom:** Very few genes mapped when using a human model on mouse data.

**Checklist:**
1. Ensure the homologs file exists at `{MODEL_DIR}/aux_files/mouse_human_homologs_one_to_one.txt`
2. Verify gene IDs are in Ensembl format before mapping
3. One-to-one homologs only cover ~15-16k genes — some loss is expected
4. Use species-matched models when possible for best results
