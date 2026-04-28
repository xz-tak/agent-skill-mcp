# scGPT Troubleshooting

## Common Errors

### CUDA Out of Memory

**Fix**:
- Reduce `batch_size` (try 16 or 8)
- Reduce `max_seq_len` (try 600 instead of 1200)
- Use AMP: `scaler = torch.cuda.amp.GradScaler()`
- Use `use_fast_transformer=True` for flash-attention (lower memory)

### FlashAttention Key Mismatch

**Error**: `Missing keys in state_dict` when loading pretrained model with flash-attention

**Cause**: Pretrained model uses standard attention, but you set `use_fast_transformer=True`

**Fix**: `load_pretrained()` handles this automatically. Make sure to use it instead of `model.load_state_dict()` directly.

### Gene Not in Vocabulary

**Error**: KeyError when looking up gene in vocab

**Fix**:
```python
# Filter to genes in vocab
gene_ids = []
valid_genes = []
for gene in adata.var.index:
    if gene in vocab:
        gene_ids.append(vocab[gene])
        valid_genes.append(gene)
adata = adata[:, valid_genes]
```

### Sequence Too Long

**Error**: `max_len exceeded` or memory issues

**Cause**: Cell has more non-zero genes than `max_length`

**Fix**: The tokenizer automatically samples/truncates. Ensure `max_length` is set:
```python
tokenize_and_pad_batch(data, gene_ids, max_len=1200, ...)
```

### Perturbation Prediction Fails

**Error**: Shapes don't match in `pred_perturb()`

**Cause**: `batch_data.x` must be shape `(n_cells, 2, n_genes)` where:
- `x[:, 0, :]` = original expression values
- `x[:, 1, :]` = perturbation flags (0=unpert, 1=pert, 2=pad)

**Fix**: Ensure perturbation data format matches expected structure.

### NaN Loss During Training

**Fix**:
- Check for NaN values in input data
- Reduce learning rate (try 1e-5)
- Enable gradient clipping (already 1.0 by default)
- Check that mask values are correct (-1 for masked, 0 for padding)

### Mask Value Confusion

**Important distinctions**:
- `pad_value = 0`: Value for padding tokens (ignored in attention)
- `mask_value = -1`: Value for masked genes (training target)
- These are NOT the same thing

### Batch Label Mismatch

**Error**: Domain-specific batch norm crashes

**Fix**: Ensure `num_batch_labels` matches actual number of unique batches:
```python
n_batches = adata.obs["batch"].nunique()
model = TransformerModel(..., num_batch_labels=n_batches)
```

### DataCollator Binning Issues

**Error**: Values out of range after binning

**Fix**: Set `do_binning=False` in DataCollator if using continuous values, or ensure `n_input_bins` is set in the model.

## Performance Tuning

| Parameter | Default | Increase For | Decrease For |
|-----------|---------|-------------|--------------|
| `batch_size` | 32-64 | GPU with large VRAM | OOM errors |
| `max_seq_len` | 1200 | More genes per cell | Speed/memory |
| `d_model` | 512 | Better representation | Speed |
| `nlayers` | 12 | Deeper model | Speed/memory |
| `mask_ratio` | 0.15 | Harder training | Easier convergence |
| `learning_rate` | 1e-4 | Faster fine-tuning | Stability |

## Task-Specific Tips

### Annotation
- Use `cell_emb_style="cls"` for best results
- Classification layers (`nlayers_cls=3`) are sufficient for most cases
- Fine-tune with `config.CLS=True`, `config.GEP=False`

### Integration
- Enable `do_dab=True` and `domain_spec_batchnorm="dsbn"`
- Use `SubsetsBatchSampler` to ensure per-batch sampling
- Weight DAB loss: `config.dab_weight=1.0`

### Perturbation
- Use `TransformerGenerator` (not `TransformerModel`)
- Fine-tune `pert_encoder` and `AffineExprDecoder` from scratch
- Freeze pretrained layers initially, then unfreeze

### Multi-omic
- Use `MultiOmicTransformerModel`
- Set `use_mod=True` with modality vocabulary

## Compatibility

- **Python**: >= 3.10
- **PyTorch**: >= 1.13.0
- **CUDA**: 11.7+ recommended
- **scanpy**: >= 1.9.1
- **Normalization**: CPM 1e4 + log1p (identical to scimilarity)
