# SCimilarity Complete API Reference

## Table of Contents
1. [CellEmbedding](#cellembedding)
2. [CellSearchKNN](#cellsearchknn)
3. [CellAnnotation](#cellannotation)
4. [CellQuery](#cellquery)
5. [Interpreter](#interpreter)
6. [Neural Network Models](#neural-network-models)
7. [Training Components](#training-components)
8. [Utility Functions](#utility-functions)
9. [Ontology Functions](#ontology-functions)
10. [Visualization Functions](#visualization-functions)
11. [Data Model Classes](#data-model-classes)

---

## CellEmbedding

**Module**: `scimilarity.cell_embedding`
**Import**: `from scimilarity import CellEmbedding`

### Constructor

```python
CellEmbedding(model_path: str, use_gpu: bool = False)
```

**Parameters**:
- `model_path`: Directory containing encoder.ckpt, gene_order.tsv, layer_sizes.json, label_ints.csv
- `use_gpu`: Use CUDA GPU for inference

**Attributes**:
- `model_path` (str): Path to model directory
- `gene_order` (list): Ordered gene symbols from gene_order.tsv
- `n_genes` (int): Number of genes in vocabulary
- `latent_dim` (int): Embedding dimensionality (typically 128)
- `model` (Encoder): PyTorch encoder network
- `int2label` (dict): Integer to cell type label mapping
- `label2int` (dict): Cell type label to integer mapping
- `filenames` (dict): Paths to model files

### get_embeddings

```python
get_embeddings(
    X: Union[scipy.sparse.csr_matrix, scipy.sparse.csc_matrix, numpy.ndarray],
    num_cells: int = -1,
    buffer_size: int = 10000,
) -> numpy.ndarray
```

**Parameters**:
- `X`: Gene-aligned, log-normalized (tp10k) expression matrix
- `num_cells`: Number of cells to embed (-1 = all)
- `buffer_size`: Batch size for processing

**Returns**: 2D numpy array [num_cells x latent_dim], L2-normalized

**Raises**: RuntimeError if NaN detected or unknown data type

---

## CellSearchKNN

**Module**: `scimilarity.cell_search_knn`
**Inherits**: CellEmbedding

### Constructor

```python
CellSearchKNN(model_path: str, knn_type: str, use_gpu: bool = False)
```

**Parameters**:
- `knn_type`: "hnswlib" or "tiledb_vector_search"

### load_knn_index

```python
load_knn_index(knn_file: str, memory_budget: int = 50000000)
```

Loads hnswlib `.bin` file or tiledb_vector_search directory.

### get_nearest_neighbors

```python
get_nearest_neighbors(
    embeddings: numpy.ndarray,
    k: int = 50,
    ef: int = 100,
) -> Tuple[numpy.ndarray, numpy.ndarray]
```

**Returns**: (nn_idxs, nn_dists) both shape [num_embeddings x k]

---

## CellAnnotation

**Module**: `scimilarity.cell_annotation`
**Import**: `from scimilarity import CellAnnotation`
**Inherits**: CellSearchKNN

### Constructor

```python
CellAnnotation(
    model_path: str,
    use_gpu: bool = False,
    filenames: Optional[dict] = None,
)
```

Loads: encoder.ckpt + annotation/labelled_kNN.bin + annotation/reference_labels.tsv

**Additional Attributes**:
- `annotation_path` (str): Path to annotation/ subdirectory
- `idx2label` (dict): kNN index to cell type label mapping
- `idx2study` (dict): kNN index to study ID mapping
- `safelist` (set|None): Allowed cell types
- `blocklist` (set|None): Blocked cell types
- `classes` (property, set): All viable prediction classes

### get_predictions_knn

```python
get_predictions_knn(
    embeddings: numpy.ndarray,
    k: int = 50,
    ef: int = 100,
    weighting: bool = False,
    disable_progress: bool = False,
) -> Tuple[pandas.Series, numpy.ndarray, numpy.ndarray, pandas.DataFrame]
```

**Returns**: (predictions, nn_idxs, nn_dists, stats)

**Stats DataFrame columns**:
| Column | Description |
|--------|-------------|
| `hits` | JSON dict of class counts among k neighbors |
| `hits_weighted` | JSON dict of inverse-distance-weighted counts |
| `min_dist` | Minimum distance to any neighbor |
| `max_dist` | Maximum distance to any neighbor |
| `vs2nd` | sum(best) / sum(best + 2nd_best) |
| `vsAll` | sum(best) / sum(all_hits) |
| `vs2nd_weighted` | Weighted version of vs2nd |
| `vsAll_weighted` | Weighted version of vsAll |

### annotate_dataset

```python
annotate_dataset(data: anndata.AnnData) -> anndata.AnnData
```

Adds to `obs`: celltype_hint, min_dist, celltype_hits, celltype_hits_weighted, celltype_hint_stat, celltype_hint_weighted_stat
Adds to `obsm`: X_scimilarity (128-dim embeddings)

Internally calls `align_dataset` — input data does NOT need to be pre-aligned.

### blocklist_celltypes

```python
blocklist_celltypes(labels: Union[List[str], Set[str]])
```

Exclude cell types from predictions. Persists across calls. Clears safelist.

### safelist_celltypes

```python
safelist_celltypes(labels: Union[List[str], Set[str]])
```

Only predict these cell types. Persists across calls. Clears blocklist.

### reset_knn

```python
reset_knn()
```

Clear all blocklist/safelist filters.

---

## CellQuery

**Module**: `scimilarity.cell_query`
**Import**: `from scimilarity import CellQuery`
**Inherits**: CellSearchKNN

### Constructor

```python
CellQuery(
    model_path: str,
    use_gpu: bool = False,
    filenames: Optional[dict] = None,
    metadata_tiledb_uri: str = "cell_metadata",
    embedding_tiledb_uri: str = "cell_embedding",
    knn_type: str = "hnswlib",
    load_knn: bool = True,
)
```

**Additional Attributes**:
- `cellsearch_path` (str): Path to cellsearch/ subdirectory
- `cell_metadata` (pd.DataFrame): Metadata for all reference cells
- `study_sample_index` (pd.Series): Indexed lookup for study/sample boundaries
- `embedding_tiledb_uri` (str): Path to TileDB embedding store

### search_nearest

```python
search_nearest(
    embeddings: numpy.ndarray,
    k: int = 10000,
    ef: int = None,           # defaults to k
    max_dist: Optional[float] = None,
) -> Tuple[List[numpy.ndarray], List[numpy.ndarray], pandas.DataFrame]
```

**Returns**: (nn_idxs, nn_dists, metadata)
- If `max_dist` set: k overridden to 1M, results filtered by distance

**Metadata columns**: study, sample, tissue, disease, data_type, index, embedding_idx, query_nn_dist

### search_centroid_nearest

```python
search_centroid_nearest(
    adata: anndata.AnnData,
    centroid_key: str,
    k: int = 10000,
    ef: int = None,
    max_dist: Optional[float] = None,
    qc: bool = True,
    qc_params: dict = {"k_clusters": 10},
    random_seed: int = 4,
) -> Tuple[numpy.ndarray, List[numpy.ndarray], List[numpy.ndarray],
          pandas.DataFrame, dict]
```

**Requirements**: `adata.layers["counts"]`, `adata.obs[centroid_key]` contains 0/1 values

**Returns**: (centroid_embedding, nn_idxs, nn_dists, metadata, qc_stats)
- `qc_stats["query_coherence"]`: Higher = query cells are more internally consistent

### search_cluster_centroids_nearest

```python
search_cluster_centroids_nearest(
    adata: anndata.AnnData,
    cluster_key: str,
    cluster_label: Optional[str] = None,
    k: int = 10000,
    ef: int = None,
    skip_null: bool = True,
    max_dist: Optional[float] = None,
) -> Tuple[numpy.ndarray, list, Dict[str, numpy.ndarray],
          Dict[str, numpy.ndarray], pandas.DataFrame]
```

**Returns**: (centroid_embeddings, cluster_idx, nn_idxs_dict, nn_dists_dict, metadata)
- Dicts keyed by cluster label

### search_exhaustive

```python
search_exhaustive(
    embeddings: numpy.ndarray,
    max_dist: float = 0.03,
    metadata_filter: Optional[dict] = None,
    buffer_size: int = 100000,
) -> Tuple[List[numpy.ndarray], List[numpy.ndarray], pandas.DataFrame]
```

**Parameters**:
- `metadata_filter`: Dict of {column: value} to filter reference cells
- `buffer_size`: Cells processed per batch (controls memory usage)

Uses `scipy.spatial.distance.cdist` with cosine metric. Sorted by lowest distance.

### search_centroid_exhaustive

```python
search_centroid_exhaustive(
    adata: anndata.AnnData,
    centroid_key: str,
    max_dist: float = 0.03,
    metadata_filter: Optional[dict] = None,
    qc: bool = True,
    qc_params: dict = {"k_clusters": 10},
    buffer_size: int = 100000,
    random_seed: int = 4,
) -> Tuple[numpy.ndarray, List[numpy.ndarray], List[numpy.ndarray],
          pandas.DataFrame, dict]
```

### search_cluster_centroids_exhaustive

```python
search_cluster_centroids_exhaustive(
    adata: anndata.AnnData,
    cluster_key: str,
    cluster_label: Optional[str] = None,
    max_dist: float = 0.03,
    metadata_filter: Optional[dict] = None,
    buffer_size: int = 100000,
    skip_null: bool = True,
) -> Tuple[numpy.ndarray, list, Dict[str, numpy.ndarray],
          Dict[str, numpy.ndarray], pandas.DataFrame]
```

### compile_sample_metadata

```python
compile_sample_metadata(
    nn_idxs: numpy.ndarray,
    levels: list = ["study", "sample", "tissue", "disease"],
) -> pandas.DataFrame
```

**Returns**: DataFrame with columns matching `levels` + cells, fraction, total

### annotate_cell_index

```python
annotate_cell_index(metadata: pandas.DataFrame) -> pandas.DataFrame
```

Adds `cell_index` column — the per-sample cell position (not global index).
Requires columns: study, sample, data_type, index.

### get_precomputed_embeddings

```python
get_precomputed_embeddings(idx: Union[slice, List[int]]) -> numpy.ndarray
```

Fast retrieval from TileDB embedding store.

---

## Interpreter

**Module**: `scimilarity.interpreter`
**Import**: `from scimilarity import Interpreter`

### Constructor

```python
Interpreter(encoder: torch.nn.Module, gene_order: list)
```

Uses Captum's `IntegratedGradients` for attribution.

### get_attributions

```python
get_attributions(
    anchors: Union[torch.Tensor, numpy.ndarray, scipy.sparse.csr_matrix],
    negatives: Union[torch.Tensor, numpy.ndarray, scipy.sparse.csr_matrix],
) -> numpy.ndarray
```

**Requirements**: `anchors.shape == negatives.shape`

**Returns**: [num_cells x num_genes] attribution matrix

Attribution logic: Integrates from negatives to anchors, masks by `anchor > negative`, takes absolute values. High attribution = gene is more expressed in anchor AND strongly affects embedding distance.

### get_ranked_genes

```python
get_ranked_genes(attrs: numpy.ndarray) -> pandas.DataFrame
```

**Returns**: DataFrame with columns: gene, gene_idx, attribution, attribution_std, cells. Sorted by mean attribution descending.

### plot_ranked_genes

```python
plot_ranked_genes(
    attrs_df: pandas.DataFrame,
    n_plot: int = 15,
    filename: Optional[str] = None,
)
```

Barplot with 95% CI error bars. Saves to filename if provided.

---

## Neural Network Models

**Module**: `scimilarity.nn_models`

### Encoder

```python
Encoder(
    n_genes: int,
    latent_dim: int = 128,
    hidden_dim: list = [1024, 1024],
    dropout: float = 0.5,
    input_dropout: float = 0.4,
)
```

**Architecture**: Input → InputDropout → [Linear, BatchNorm, PReLU] → [Dropout, Linear, BatchNorm, PReLU] × N → Linear → L2-Normalize

**Methods**:
- `forward(x) -> torch.Tensor` — L2-normalized embeddings
- `save_state(filename: str)` — Save state_dict
- `load_state(filename: str, use_gpu: bool = False)` — Load state_dict

### Decoder

```python
Decoder(
    n_genes: int,
    latent_dim: int = 128,
    hidden_dim: list = [1024, 1024],
    dropout: float = 0.5,
)
```

Mirror architecture of Encoder (reversed layers). Used for reconstruction loss during training.

---

## Training Components

**Module**: `scimilarity.training_models`

### MetricLearning (PyTorch Lightning)

```python
MetricLearning(
    n_genes: int,
    latent_dim: int = 128,
    hidden_dim: list = [1024, 1024, 1024],
    dropout: float = 0.5,
    input_dropout: float = 0.4,
    triplet_loss_weight: float = 0.001,
    margin: float = 0.05,
    negative_selection: str = "semihard",    # "semihard", "hardest", "random"
    sample_across_studies: bool = False,
    perturb_labels: bool = False,
    perturb_labels_fraction: float = 0.5,
    lr: float = 5e-3,
    l1_lambda: float = 0.0001,
    l2_lambda: float = 0.01,
    max_epochs: int = 500,
    cosine_annealing_tmax: Optional[int] = None,
    track_triplets: Optional[str] = None,
)
```

**Loss**: `triplet_loss_weight * triplet_loss + reconstruction_mse_loss`
**Optimizer**: AdamW with CosineAnnealingLR

**Key Methods**:
- `forward(x) -> Tuple[embeddings, reconstruction]`
- `training_step(batch, batch_idx) -> loss`
- `validation_step(batch, batch_idx) -> loss`
- `save_all(model_path)` — Saves encoder.ckpt, decoder.ckpt, layer_sizes.json, gene_order.tsv, label_ints.csv, hyperparameters.json, metadata.json
- `load_state(encoder_filename, decoder_filename, use_gpu=False, freeze=False)`

### TripletSelector

**Module**: `scimilarity.triplet_selector`

```python
TripletSelector(
    margin: float,
    negative_selection: str,          # "semihard", "hardest", "random"
    perturb_labels: bool = False,
    perturb_labels_fraction: float = 0.5,
)
```

**Key Methods**:
- `get_triplets_idx(embeddings, labels, int2label, studies=None) -> Tuple[triplets, num_hard, num_viable]`

### TripletLoss

```python
TripletLoss(margin: float, triplet_selector: TripletSelector)
```

---

## Utility Functions

**Module**: `scimilarity.utils`
**Import**: `from scimilarity.utils import <function>` or `from scimilarity import align_dataset, lognorm_counts`

### Data Preparation

```python
lognorm_counts(data: anndata.AnnData) -> anndata.AnnData
```
Log normalize: `normalize_total(target_sum=1e4)` + `log1p()`. Requires `layers["counts"]`.

```python
align_dataset(
    data: anndata.AnnData,
    target_gene_order: list,
    keep_obsm: bool = True,
    gene_overlap_threshold: int = 5000,
) -> anndata.AnnData
```
Aligns gene space to target order. Missing genes filled with zeros (sparse). Minimum overlap: 5000.

```python
filter_cells(
    data: anndata.AnnData,
    min_genes: int = 400,
    mito_prefix: Optional[str] = None,  # auto-detects "MT-" or "mt-"
    mito_percent: float = 30.0,
) -> anndata.AnnData
```
QC filter on gene count and mitochondrial percentage. Requires `layers["counts"]`.

```python
convert_id2symbol(adata: anndata.AnnData, mapping_table: str) -> anndata.AnnData
```
Convert Ensembl IDs to gene symbols using TSV mapping table.

```python
consolidate_duplicate_symbols(adata: anndata.AnnData) -> anndata.AnnData
```
Merge duplicate gene symbols by summing counts.

### Centroid Operations

```python
get_centroid(counts: Union[scipy.sparse.csr_matrix, numpy.ndarray]) -> numpy.ndarray
```
Compute mean expression centroid with log normalization (1e4).

```python
get_cluster_centroids(
    data: anndata.AnnData,
    target_gene_order: numpy.ndarray,
    cluster_key: str,
    cluster_label: Optional[str] = None,
    skip_null: bool = True,
) -> Tuple[numpy.ndarray, list]
```
Returns: (centroids, cluster_idx). Centroids are log-normalized.

```python
get_dist2centroid(
    centroid_embedding: numpy.ndarray,
    X: Union[scipy.sparse.csr_matrix, numpy.ndarray],
) -> numpy.ndarray
```
Cosine distances from cells to centroid embedding.

### Pseudobulk

```python
pseudobulk_anndata(
    adata: anndata.AnnData,
    groupby_labels: Union[str, list],
    qc_filters: Optional[dict] = None,  # mito_percent, min_counts, min_genes
    min_num_cells: int = 1,
    only_orig_genes: bool = False,
) -> anndata.AnnData
```
Returns AnnData with `layers["counts"]` and `layers["detection"]`.

### TileDB Operations

```python
write_array_to_tiledb(tdb, arr, value_type, row_start=0, batch_size=100000)
write_csr_to_tiledb(tdb, matrix, value_type, row_start=0, batch_size=25000)
optimize_tiledb_array(tiledb_array_uri, config=None, verbose=True)
query_tiledb_df(tdb, query_condition, attrs=None) -> pandas.DataFrame
embedding_from_tiledb(cell_idx, embedding_tdb_uri, config=None) -> numpy.ndarray
```

### DataFrame Helpers

```python
subset_by_unique_values(adata, groupby_column, subset_column, unique_threshold) -> anndata.AnnData
subset_by_frequency(adata, groupby_column, min_frequency) -> anndata.AnnData
categorize_and_sort_by_score(metadata, score_column, category_column,
    max_categories=None, ascending=False) -> pandas.DataFrame
clean_tissues(tissues: pandas.Series) -> pandas.Series
clean_diseases(diseases: pandas.Series) -> pandas.Series
```

---

## Ontology Functions

**Module**: `scimilarity.ontologies`

### Import Ontologies

```python
import_cell_ontology(url="http://purl.obolibrary.org/obo/cl/cl-basic.obo") -> networkx.DiGraph
import_uberon_ontology(url="http://purl.obolibrary.org/obo/uberon/basic.obo") -> networkx.DiGraph
import_doid_ontology(url="http://purl.obolibrary.org/obo/doid.obo") -> networkx.DiGraph
import_mondo_ontology(url="http://purl.obolibrary.org/obo/mondo.obo") -> networkx.DiGraph
```

### Graph Traversal

```python
get_id_mapper(graph) -> dict                    # term ID -> name
get_children(graph, node, node_list=None) -> set
get_parents(graph, node, node_list=None) -> set
get_siblings(graph, node, node_list=None) -> set
get_all_ancestors(graph, node, node_list=None, inclusive=False) -> set
get_all_descendants(graph, nodes, node_list=None, inclusive=False) -> set
get_lowest_common_ancestor(graph, node1, node2) -> str
find_most_viable_parent(graph, node, node_list) -> str
```

### Similarity

```python
ontology_similarity(graph, node1, node2, restricted_set=None) -> int
```
Minimum path distance between nodes.

```python
all_pair_similarities(graph, nodes, restricted_set=None) -> pandas.DataFrame
```
Pairwise distance matrix.

```python
ontology_silhouette_width(graph, nodes, restricted_set=None) -> float
```

---

## Visualization Functions

**Module**: `scimilarity.visualizations`

### Circle Packing

```python
aggregate_counts(data: pandas.DataFrame, levels: List[str]) -> dict
assign_size(data_dict, data, levels, size_column, name_column) -> dict
assign_suffix(data_dict, suffix) -> dict
assign_colors(data_dict, colors_dict) -> dict
get_children_data(data_dict) -> List[dict]
circ_dict2data(circ_dict) -> List[dict]

draw_circles(
    circ_dict: dict,
    show_value: bool = True,
    show_label: bool = True,
    figsize: Tuple = (12, 12),
    legend_loc: str = "upper left",
    filename: Optional[str] = None,
)
```

### Convenience Functions

```python
hits_circles(
    data: pandas.DataFrame,
    levels: List[str],
    label_column: str,
    value_column: str,
    figsize: Tuple = (20, 20),
    filename: Optional[str] = None,
)
```
Visualize search hits as hierarchical circle packing.

```python
hits_heatmap(
    data: pandas.DataFrame,
    rows: str,
    columns: str,
    values: str,
    figsize: Tuple = (12, 12),
    cmap: str = "YlOrRd",
    filename: Optional[str] = None,
)
```
Visualize search hits as clustered heatmap.

---

## Data Model Classes

### AnnData-based (`scimilarity.anndata_data_models`)

```python
class scDataset(torch.utils.data.Dataset)
    # __getitem__ returns: (expression_vector, label, study)

class scCollator
    # Collates batch data into tensors

class MetricLearningDataModule(pl.LightningDataModule)
    # DataModule for h5ad-based training
```

### Zarr-based (`scimilarity.zarr_data_models`)

```python
class ZarrDataset
    # Properties: dataset_info, shape
    # Methods: get_cell(idx), get_obs(field), var_index, etc.

class scDataset(Dataset)
    # Multi-zarr dataset wrapper

class MetricLearningDataModule(pl.LightningDataModule)
    # Params: batch_size, num_workers, obs_field, gene_order path
```

### TileDB-based (`scimilarity.tiledb_data_models`)

```python
class scDataset(Dataset)     # TileDB-backed dataset
class scSampler(Sampler)     # Weighted random sampler with dynamic weights
class scCollator              # TileDB data collation
class CellMultisetDataModule(pl.LightningDataModule)  # TileDB training
```
