# Geneformer Quick Reference

## 1. Decision Tree: Which Class to Use?

```
What do you want to do?
|
+-- Convert raw scRNA-seq to tokenized format
|   --> TranscriptomeTokenizer
|
+-- Extract embeddings from cells/genes
|   --> EmbExtractor
|
+-- Classify cells by state/type
|   --> Classifier (classifier="cell")
|
+-- Classify genes by function
|   --> Classifier (classifier="gene")
|
+-- Simulate gene perturbation effects
|   +-- Run perturbation --> InSilicoPerturber
|   +-- Analyze results  --> InSilicoPerturberStats
|
+-- Multi-task classification
|   --> MTLClassifier
|
+-- Continue pretraining on custom corpus
    --> GeneformerPretrainer
```

## 2. Decision Tree: Which Model Tier?

```
Choose your model:
|
+-- Standard analysis (default, fast)
|   --> V2-104M (104M params, ~399MB)
|
+-- Maximum accuracy (large, slower)
|   --> V2-316M (316M params, ~1.2GB)
|
+-- Cancer-specific analysis
    --> V2-104M-CLcancer (cancer continual learning)
```

## 3. Decision Tree: Perturbation Type?

```
What perturbation to simulate?
|
+-- Gene knockout (remove gene expression)
|   --> perturb_type="delete"
|
+-- Gene overexpression (max rank)
|   --> perturb_type="overexpress"
|
+-- Partial inhibition (shift rank down)
|   --> perturb_type="inhibit", perturb_rank_shift={1,2,3}
|
+-- Partial activation (shift rank up)
    --> perturb_type="activate", perturb_rank_shift={1,2,3}
```

## 4. Decision Tree: Stats Mode?

```
What analysis do you need?
|
+-- Genes that shift cells toward a desired state
|   --> mode="goal_state_shift" (requires cell_states_to_model)
|
+-- Compare perturbation vs null distribution
|   --> mode="vs_null" (requires null_dist_data)
|
+-- Find impactful perturbations (undirected)
|   --> mode="mixture_model"
|
+-- Aggregate shifts for a single perturbation across cells
|   --> mode="aggregate_data"
|
+-- Aggregate gene-level shifts across perturbations
    --> mode="aggregate_gene_shifts"
```

## 5. Decision Tree: Embedding Mode?

```
What embedding do you need?
|
+-- Cell-level (one vector per cell)
|   +-- CLS token (V2 recommended)       --> emb_mode="cls"
|   +-- Mean pooling of gene embeddings   --> emb_mode="cell"
|
+-- Gene-level (one vector per gene)
|   --> emb_mode="gene"
|
+-- Both cell and gene
    +-- CLS + gene  --> emb_mode="cls_and_gene"
    +-- Cell + gene --> emb_mode="cell_and_gene"
```

## 6. Common Workflows Cheat Sheet

### Tokenize raw data
1. Prepare `.h5ad` or `.loom` with `ensembl_id` and `n_counts`
2. Create `TranscriptomeTokenizer(custom_attr_name_dict={...})`
3. Call `tokenizer.tokenize_data(input_dir, output_dir, output_prefix)`
4. Load result: `datasets.load_from_disk("output.dataset")`

### Extract embeddings
1. Load tokenized dataset with `datasets.load_from_disk()`
2. Create `EmbExtractor(model_type="Pretrained", emb_mode="cls")`
3. Call `embex.extract_embs(model_dir, input_data, output_dir, output_prefix)`
4. Result is a pandas DataFrame with cell metadata + embedding columns

### Classify cells
1. Prepare tokenized + labeled dataset (label column in dataset)
2. Create `Classifier(classifier="cell", cell_state_dict={...})`
3. Call `cc.prepare_data(input_data, ...)` to split train/test
4. Call `cc.train(model_dir, train_data, ...)` to fine-tune
5. Call `cc.evaluate(model_dir, test_data, ...)` for metrics

### Perturb (goal state shift)
1. Define `cell_states_to_model = {"state_key": {"start_state": ..., "goal_state": ...}}`
2. Create `InSilicoPerturber(perturb_type="delete", cell_states_to_model=...)`
3. Call `isp.perturb_data(model_dir, input_data, output_dir)`
4. Create `InSilicoPerturberStats(mode="goal_state_shift", cell_states_to_model=...)`
5. Call `ispstats.get_stats(isp_output_dir, ...)` for ranked gene list

### Perturb (undirected / mixture model)
1. Create `InSilicoPerturber(perturb_type="delete")`
2. Call `isp.perturb_data(model_dir, input_data, output_dir)`
3. Create `InSilicoPerturberStats(mode="mixture_model")`
4. Call `ispstats.get_stats(isp_output_dir, ...)` for impact scores

### Multi-task classification
1. Ensure dataset has `unique_cell_id` column
2. Create `MTLClassifier(task_columns=["task1", "task2"])`
3. Call `mtl.prepare_data(input_data, ...)` to split
4. Call `mtl.train(model_dir, train_data, ...)` to fine-tune
5. Call `mtl.evaluate(model_dir, test_data, ...)` for per-task metrics

## 7. Parameter Quick Reference Tables

### TranscriptomeTokenizer

| Parameter | Type | Default | Common Values |
|-----------|------|---------|---------------|
| custom_attr_name_dict | dict | None | `{"cell_type": "cell_type"}` |
| nproc | int | 1 | 4, 8, 16 |
| gene_median_file | str | built-in | path to custom median file |
| token_dictionary_file | str | built-in | path to custom token dict |
| special_token | bool | True (V2) | True for CLS mode |
| collapse_gene_ids | bool | True | True |
| gene_mapping_file | str | None | path to ID mapping file |

### EmbExtractor

| Parameter | Type | Default | Common Values |
|-----------|------|---------|---------------|
| model_type | str | "Pretrained" | "Pretrained", "GeneClassifier", "CellClassifier" |
| emb_mode | str | "cls" | "cls", "cell", "gene", "cls_and_gene", "cell_and_gene" |
| emb_layer | int | -1 | -1 (last), 0 (first), specific layer index |
| forward_batch_size | int | 200 | 50, 100, 200, 500 |
| max_ncells | int | None | 1000, 5000, None |
| summary_stat | str | "exact_mean" | "exact_mean", "exact_median" |
| quantize | bool | False | True for low memory |
| nproc | int | 1 | 4, 8 |

### Classifier

| Parameter | Type | Default | Common Values |
|-----------|------|---------|---------------|
| classifier | str | "cell" | "cell", "gene" |
| cell_state_dict | dict | None | `{"state_key": {"start_state": ..., "goal_state": ...}}` |
| forward_batch_size | int | 200 | 50, 100, 200 |
| epochs | int | 10 | 5, 10, 20 |
| learning_rate | float | 5e-5 | 1e-5, 5e-5, 1e-4 |
| freeze_layers | int | 0 | 0, 2, 4 |
| num_crossval_splits | int | 1 | 1, 5, 10 |
| quantize | bool | False | True for low memory |

### InSilicoPerturber

| Parameter | Type | Default | Common Values |
|-----------|------|---------|---------------|
| perturb_type | str | "delete" | "delete", "overexpress", "inhibit", "activate" |
| perturb_rank_shift | set | None | {1}, {1,2,3} (for inhibit/activate) |
| genes_to_perturb | list | "all" | "all", ["ENSG00000141510"] |
| cell_states_to_model | dict | None | `{"state_key": {"start_state": ..., "goal_state": ...}}` |
| emb_mode | str | "cls" | "cls", "cell" |
| forward_batch_size | int | 200 | 50, 100, 200 |
| max_ncells | int | None | 1000, 5000 |
| nproc | int | 1 | 4, 8 |
| quantize | bool | False | True for low memory |

### InSilicoPerturberStats

| Parameter | Type | Default | Common Values |
|-----------|------|---------|---------------|
| mode | str | required | "goal_state_shift", "vs_null", "mixture_model", "aggregate_data", "aggregate_gene_shifts" |
| cell_states_to_model | dict | None | required for "goal_state_shift" |
| pickle_suffix | str | "_raw.pickle" | must match ISP output files |
| null_dist_data | str | None | path to null distribution (for "vs_null") |

### MTLClassifier

| Parameter | Type | Default | Common Values |
|-----------|------|---------|---------------|
| task_columns | list | required | ["cell_type", "disease"] |
| forward_batch_size | int | 200 | 50, 100, 200 |
| epochs | int | 10 | 5, 10, 20 |
| learning_rate | float | 5e-5 | 1e-5, 5e-5 |
| freeze_layers | int | 0 | 0, 2, 4 |
| quantize | bool | False | True for low memory |

## 8. File Format Requirements

| Format | Gene Attribute | Cell Attribute | Expression Data |
|--------|---------------|----------------|-----------------|
| `.h5ad` (AnnData) | `adata.var["ensembl_id"]` -- Ensembl gene IDs | `adata.obs["n_counts"]` -- total UMI counts per cell | `adata.X` -- raw counts (not normalized) |
| `.loom` | Row attribute `ensembl_id` | Column attribute `n_counts` | Main matrix -- raw counts |
| `.zarr` | `var/ensembl_id` array | `obs/n_counts` array | `X` array -- raw counts |
| Tokenized `.dataset` | N/A (genes encoded as tokens) | Metadata columns preserved from input | `input_ids` -- rank-ordered token list per cell |
