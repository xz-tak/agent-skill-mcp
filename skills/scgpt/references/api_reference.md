# scGPT Complete API Reference

## Table of Contents
1. [TransformerModel](#transformermodel)
2. [TransformerGenerator](#transformergenerator)
3. [MultiOmicTransformerModel](#multiomictransformermodel)
4. [Encoder Classes](#encoder-classes)
5. [Decoder Classes](#decoder-classes)
6. [Tokenizer](#tokenizer)
7. [Preprocessor](#preprocessor)
8. [DataCollator](#datacollator)
9. [Trainer Functions](#trainer-functions)
10. [Cell Embedding Tasks](#cell-embedding-tasks)
11. [GRN Inference](#grn-inference)
12. [Loss Functions](#loss-functions)
13. [Utility Functions](#utility-functions)
14. [Domain Adaptation](#domain-adaptation)

---

## TransformerModel

**Module**: `scgpt.model.model`

```python
TransformerModel(
    ntoken: int,                          # vocab size
    d_model: int,                         # embedding dim (256-768)
    nhead: int,                           # attention heads
    d_hid: int,                           # FFN hidden dim
    nlayers: int,                         # transformer layers
    nlayers_cls: int = 3,                 # classification head layers
    n_cls: int = 1,                       # number of classes
    vocab: Any = None,                    # GeneVocab object
    dropout: float = 0.5,
    pad_token: str = "<pad>",
    pad_value: int = 0,
    do_mvc: bool = False,                 # masked value cell prediction
    do_dab: bool = False,                 # domain adversarial batch
    use_batch_labels: bool = False,
    num_batch_labels: Optional[int] = None,
    domain_spec_batchnorm: Union[bool, str] = False,  # "dsbn" or False
    input_emb_style: str = "continuous",  # "continuous"|"category"|"scaling"
    n_input_bins: Optional[int] = None,
    cell_emb_style: str = "cls",          # "cls"|"avg-pool"|"w-pool"
    mvc_decoder_style: str = "inner product",
    ecs_threshold: float = 0.3,
    explicit_zero_prob: bool = False,
    use_fast_transformer: bool = False,
    fast_transformer_backend: str = "flash",
    pre_norm: bool = False,
)
```

### Methods

**forward**:
```python
forward(
    src: Tensor,                    # gene IDs (batch, seq_len)
    values: Tensor,                 # expression values (batch, seq_len)
    src_key_padding_mask: Tensor,   # bool mask (batch, seq_len)
    batch_labels: Optional[Tensor] = None,
    CLS: bool = False,
    CCE: bool = False,
    MVC: bool = False,
    ECS: bool = False,
    do_sample: bool = False,
) -> Mapping[str, Tensor]
```
Returns dict: `mlm_output`, `cell_emb`, optionally `cls_output`, `mvc_output`, `loss_cce`, `loss_ecs`, `dab_output`

**encode_batch**:
```python
encode_batch(src, values, src_key_padding_mask, batch_size,
    batch_labels=None, output_to_cpu=True, time_step=None, return_np=False) -> Tensor
```

**generate**:
```python
generate(cell_emb, src, values=None, src_key_padding_mask=None,
    gen_iters=1, batch_labels=None) -> Tensor
```
Generate expression values given cell embedding.

---

## TransformerGenerator

**Module**: `scgpt.model.generation_model`

Extends TransformerModel with perturbation support.

```python
TransformerGenerator(
    # Same as TransformerModel plus:
    pert_pad_id: int = 2,
    decoder_activation: Optional[str] = None,
    decoder_adaptive_bias: bool = False,
)
```

### Additional Methods

**forward**: Same as TransformerModel but takes `input_pert_flags: Tensor` instead of `batch_labels`.

**pred_perturb**:
```python
pred_perturb(
    batch_data,                          # batch with x[:, 0]=values, x[:, 1]=pert_flags
    include_zero_gene: str = "batch-wise",  # "all"|"batch-wise"|None
    gene_ids = None,                     # vocab-mapped gene IDs
    amp: bool = True,                    # mixed precision
) -> Tensor                              # (batch_size, n_genes)
```

---

## MultiOmicTransformerModel

**Module**: `scgpt.model.multiomic_model`

```python
MultiOmicTransformerModel(
    # Same as TransformerModel plus:
    use_mod: bool = False,
    ntokens_mod: Optional[int] = None,
    vocab_mod: Optional[Any] = None,
)
```

Forward takes additional `mod_types: Optional[Tensor]`.

---

## Encoder Classes

**Module**: `scgpt.model.model`

### GeneEncoder
```python
GeneEncoder(num_embeddings: int, embedding_dim: int, padding_idx: Optional[int] = None)
forward(x: Tensor) -> Tensor  # (batch, seq_len, embsize)
```

### ContinuousValueEncoder
```python
ContinuousValueEncoder(d_model: int, dropout: float = 0.1, max_value: int = 512)
forward(x: Tensor) -> Tensor  # (batch, seq_len, d_model)
```
MLP: Linear(1, d_model) -> ReLU -> Linear(d_model, d_model) -> LayerNorm. Clamps to max_value.

### CategoryValueEncoder
```python
CategoryValueEncoder(num_embeddings: int, embedding_dim: int, padding_idx: Optional[int] = None)
forward(x: Tensor) -> Tensor  # (batch, seq_len, embsize)
```

### BatchLabelEncoder
```python
BatchLabelEncoder(num_embeddings: int, embedding_dim: int, padding_idx: Optional[int] = None)
forward(x: Tensor) -> Tensor  # (batch, embsize)
```

---

## Decoder Classes

**Module**: `scgpt.model.model` and `scgpt.model.generation_model`

### ExprDecoder
```python
ExprDecoder(d_model: int, explicit_zero_prob: bool = False, use_batch_labels: bool = False)
forward(x: Tensor) -> Dict[str, Tensor]  # {"pred": (batch, seq_len), "zero_probs": ...}
```

### ClsDecoder
```python
ClsDecoder(d_model: int, n_cls: int, nlayers: int = 3, activation: callable = nn.ReLU)
forward(x: Tensor) -> Tensor  # (batch, n_cls)
```
Residual MLP with layer norms.

### MVCDecoder
```python
MVCDecoder(
    d_model: int,
    arch_style: str = "inner product",  # "inner product"|"concat query"|"sum query"
    query_activation: nn.Module = nn.Sigmoid,
    hidden_activation: nn.Module = nn.PReLU,
    explicit_zero_prob: bool = False,
    use_batch_labels: bool = False,
)
forward(cell_emb: Tensor, gene_embs: Tensor) -> Union[Tensor, Dict[str, Tensor]]
```

### AffineExprDecoder
```python
AffineExprDecoder(
    d_model: int,
    explicit_zero_prob: bool = False,
    activation: Optional[str] = None,
    tanh_coeff: bool = False,
    adaptive_bias: bool = False,
)
forward(x: Tensor, values: Tensor) -> Dict[str, Tensor]
# Returns {"pred": coeff * values + bias}
```

### AdversarialDiscriminator
```python
AdversarialDiscriminator(d_model: int, n_cls: int, nlayers: int = 3,
    activation: callable = nn.LeakyReLU, reverse_grad: bool = False)
forward(x: Tensor) -> Tensor  # (batch, n_cls)
```

---

## Tokenizer

**Module**: `scgpt.tokenizer.gene_tokenizer`

### GeneVocab
```python
GeneVocab(gene_list_or_vocab: Union[List[str], Vocab],
    specials: Optional[List[str]] = None, special_first: bool = True,
    default_token: Optional[str] = "<pad>")

# Class methods
GeneVocab.from_file(file_path: Union[Path, str]) -> GeneVocab  # .json or .pkl
GeneVocab.from_dict(token2idx: Dict[str, int], default_token="<pad>") -> GeneVocab

# Instance methods
save_json(file_path)
set_default_token(default_token: str)
```

### Functions
```python
tokenize_batch(data, gene_ids, return_pt=True, append_cls=True,
    include_zero_gene=False, cls_id="<cls>", mod_type=None) -> List[Tuple]

pad_batch(batch, max_len, vocab, pad_token="<pad>", pad_value=0,
    cls_appended=True, vocab_mod=None) -> Dict[str, Tensor]

tokenize_and_pad_batch(data, gene_ids, max_len, vocab, pad_token, pad_value,
    append_cls=True, include_zero_gene=False, cls_token="<cls>",
    return_pt=True, mod_type=None, vocab_mod=None) -> Dict[str, Tensor]
# Returns: {"genes": Tensor, "values": Tensor, "mod_types": Tensor (optional)}

random_mask_value(values: Tensor, mask_ratio=0.15,
    mask_value=-1, pad_value=0) -> Tensor

get_default_gene_vocab() -> GeneVocab  # 48,292 HGNC symbols
```

---

## Preprocessor

**Module**: `scgpt.preprocess`

```python
Preprocessor(
    use_key: Optional[str] = None,
    filter_gene_by_counts: Union[int, bool] = False,
    filter_cell_by_counts: Union[int, bool] = False,
    normalize_total: Union[float, bool] = 1e4,
    result_normed_key: Optional[str] = "X_normed",
    log1p: bool = False,
    result_log1p_key: str = "X_log1p",
    subset_hvg: Union[int, bool] = False,
    hvg_use_key: Optional[str] = None,
    hvg_flavor: str = "seurat_v3",
    binning: Optional[int] = None,
    result_binned_key: str = "X_binned",
)

__call__(adata: AnnData, batch_key: Optional[str] = None) -> Dict
check_logged(adata: AnnData, obs_key: Optional[str] = None) -> bool
```

### Helper Functions
```python
binning(row: Union[np.ndarray, Tensor], n_bins: int) -> Union[np.ndarray, Tensor]
```

---

## DataCollator

**Module**: `scgpt.data_collator`

```python
@dataclass
class DataCollator:
    do_padding: bool = True
    pad_token_id: Optional[int] = None
    pad_value: int = 0
    do_mlm: bool = True
    do_binning: bool = True
    mlm_probability: float = 0.15
    mask_value: int = -1
    max_length: Optional[int] = None
    sampling: bool = True
    keep_first_n_tokens: int = 1  # preserve CLS

    __call__(examples: List[Dict[str, Tensor]]) -> Dict[str, Tensor]
    # Returns: {"gene": Tensor, "expr": Tensor, "masked_expr": Tensor}
```

### SubsetsBatchSampler

**Module**: `scgpt.data_sampler`

```python
SubsetsBatchSampler(
    subsets: List[Sequence[int]],
    batch_size: int,
    intra_subset_shuffle: bool = True,
    inter_subset_shuffle: bool = True,
    drop_last: bool = False,
)
```
Ensures each batch comes from a single subset/domain.

---

## Trainer Functions

**Module**: `scgpt.trainer`

```python
prepare_data(tokenized_train, tokenized_valid, train_batch_labels,
    valid_batch_labels, config, epoch, train_celltype_labels=None,
    valid_celltype_labels=None, sort_seq_batch=False) -> Tuple[Dict, Dict]
# Returns dicts with: gene_ids, values, target_values, batch_labels,
#   celltype_labels (annotation), mod_types (multiomic)

prepare_dataloader(data_pt, batch_size, shuffle=False,
    intra_domain_shuffle=False, drop_last=False, num_workers=0,
    per_seq_batch_sample=False) -> DataLoader

train(model, loader, vocab, criterion_gep_gepc, criterion_dab,
    criterion_cls, scaler, optimizer, scheduler, device, config,
    logger, epoch) -> None

evaluate(model, loader, vocab, criterion_gep_gepc, criterion_dab,
    criterion_cls, device, config, epoch) -> float

predict(model, loader, vocab, config, device) -> np.ndarray

test(model, adata, gene_ids, vocab, config, device, logger)
    -> Tuple[np.ndarray, np.ndarray, Dict]
# Returns: (predictions, labels, metrics_dict)

eval_testdata(model, adata_t, gene_ids, vocab, config, logger,
    include_types=["cls"]) -> Optional[Dict]
# Stores embeddings in adata.obsm["X_scGPT"]

define_wandb_metrcis() -> None
```

**Config attributes used by trainer**: `task` ("annotation"|"integration"|"perturb"|"multiomic"), `GEP`, `GEPC`, `CLS`, `ESC`, `DAR`, `mask_ratio`, `mask_value`, `pad_value`, `pad_token`, `dab_weight`, `explicit_zero_prob`, `amp`

---

## Cell Embedding Tasks

**Module**: `scgpt.tasks.cell_emb`

```python
get_batch_cell_embeddings(
    adata,
    cell_embedding_mode: str = "cls",
    model = None,
    vocab = None,
    max_length: int = 1200,
    batch_size: int = 64,
    model_configs = None,  # {"pad_value": 0, "pad_token": "<pad>"}
    gene_ids = None,
    use_batch_labels: bool = False,
) -> np.ndarray  # (n_cells, d_model), L2-normalized

embed_data(
    adata_or_file: Union[AnnData, PathLike],
    model_dir: PathLike,
    gene_col: str = "feature_name",
    max_length: int = 1200,
    batch_size: int = 64,
    obs_to_save: Optional[list] = None,
    device: Union[str, torch.device] = "cuda",
    use_fast_transformer: bool = True,
    return_new_adata: bool = False,
) -> AnnData  # with obsm["X_scGPT"]
```

---

## GRN Inference

**Module**: `scgpt.tasks.grn`

```python
class GeneEmbedding:
    __init__(embeddings: Mapping)  # gene_name -> vector dict

    # Class methods
    read_embedding(filename) -> Dict[str, np.ndarray]

    # Analysis
    compute_similarities(gene, subset=None, feature_type=None) -> pd.DataFrame
    generate_vector(genes) -> List
    generate_weighted_vector(genes, weights) -> List
    generate_network(threshold: float = 0.5) -> nx.Graph
    get_similar_genes(vector) -> pd.DataFrame

    # Clustering
    get_adata(resolution=20) -> AnnData
    get_metagenes(gdata) -> Dict
    cluster_definitions_as_df(top_n=20) -> pd.DataFrame

    # Visualization
    plot_similarities(gene, n_genes=10, save=None)
    plot_metagene(gdata, mg=None, title="Gene Embedding")
    plot_metagenes_scores(adata, metagenes, column, plot=None)
    score_metagenes(adata, metagenes)
```

---

## Loss Functions

**Module**: `scgpt.loss`

```python
masked_mse_loss(input: Tensor, target: Tensor, mask: Tensor) -> Tensor
# MSE only on masked positions, normalized by mask count

criterion_neg_log_bernoulli(input: Tensor, target: Tensor, mask: Tensor) -> Tensor
# Negative log-likelihood for zero/non-zero prediction

masked_relative_error(input: Tensor, target: Tensor, mask: Tensor) -> Tensor
# Mean |input - target| / (target + 1e-6) on masked positions
```

---

## Utility Functions

**Module**: `scgpt.utils.util`

```python
set_seed(seed: int) -> None
add_file_handler(logger, log_file_path) -> None
category_str2int(category_strs: List[str]) -> List[int]
isnotebook() -> bool
get_free_gpu() -> int
get_git_commit() -> str
histogram(*data, label=["train", "valid"], ...) -> axes.Axes
tensorlist2tensor(tensorlist, pad_value) -> Tensor
map_raw_id_to_vocab_id(raw_ids, gene_ids) -> Union[np.ndarray, Tensor]

load_pretrained(
    model: nn.Module,
    pretrained_params: dict,
    strict: bool = False,
    prefix: Optional[List[str]] = None,
    verbose: bool = True,
) -> nn.Module
# Loads pretrained weights, handles FlashAttention key conversion

find_required_colums(adata, id, configs_dir, update=False) -> List[Optional[str]]
```

---

## Domain Adaptation

**Module**: `scgpt.model.dsbn`

```python
DomainSpecificBatchNorm1d(num_features, num_domains, eps=1e-5,
    momentum=0.1, affine=True, track_running_stats=True)
forward(x: Tensor, domain_label: int) -> Tensor
```

**Module**: `scgpt.model.grad_reverse`

```python
grad_reverse(x: Tensor, lambd: float = 1.0) -> Tensor
# Reverses gradients with scaling factor
```

### FlashTransformerEncoderLayer

```python
FlashTransformerEncoderLayer(d_model, nhead, dim_feedforward=2048,
    dropout=0.1, activation="relu", layer_norm_eps=1e-5,
    batch_first=True, norm_scheme="post")
forward(src, src_mask=None, src_key_padding_mask=None) -> Tensor
```
