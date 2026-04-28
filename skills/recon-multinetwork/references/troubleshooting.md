# ReCoN Multinetwork Troubleshooting

Known issues, bugs, and workarounds discovered during development.

---

## 1. multixrank_patch: Division by Zero in Sparse CCC Networks

**Symptom**: `NaN` propagation in ReCoN results, especially for conditions with fewer CCC interactions (e.g., IPF with 14,513 vs SSc with 20,051 interactions).

**Root cause**: In `multixrank.TransitionMatrix.get_normalization_bipartite_alpha_beta`, when a node has no inter-layer connections, `norm` becomes 0. The original code divides by `norm`, producing `NaN` that propagates through the random walk.

**Fix**: Apply the monkey-patch before any ReCoN analysis:
```python
from multixrank_patch import apply_transition_matrix_patch
apply_transition_matrix_patch()
```

The patch uses `epsilon = 1e-10` when `norm == 0`, allowing minimal signal flow for isolated nodes while preventing `NaN` propagation.

**When it occurs**: 10+ cell types with sparse CCC network where some cell-type pairs have no interactions.

---

## 2. receptor_grn Mutation Bug

**Symptom**: After processing the first condition, subsequent conditions produce garbage results with double-suffixed node names like `"gene_receptor_receptor"`.

**Root cause**: `multicell_targets()` modifies the `receptor_grn` DataFrame in-place, adding `"_receptor"` suffix to source column values. If the same DataFrame is reused for the next condition, the suffix gets applied again.

**Fix**: Always `.copy()` receptor_grn before passing to multicell_targets:
```python
receptor_grn_copy = receptor_grn.copy()
direct, indirect, stats = run_recon_analysis(
    condition, grns, ccc, receptor_grn_copy, seeds
)
```

---

## 3. direct/indirect .copy() Before combine_effects

**Symptom**: Saved direct/indirect CSVs contain `NaN` columns, or condition comparisons show unexpected zeros.

**Root cause**: `combine_effects()` normalizes the `direct_effect` DataFrame in-place. It divides each column by its sum. If a column sum is ~0, this produces `NaN` that overwrites the original data.

**Fix**: Save copies of direct and indirect BEFORE calling combine_effects:
```python
direct_copy = direct.copy()
indirect_copy = indirect.copy()
combined = combine_effects(direct, indirect, alpha=0.8)

# Save the preserved copies, NOT the originals (which are now modified)
direct_copy.to_csv("direct_effects.csv")
indirect_copy.to_csv("indirect_effects.csv")
```

---

## 4. CellChat Column Name Standardization

**Symptom**: KeyError on `celltype_source` or `celltype_target` when using CellChat CCC data.

**Root cause**: CellChat uses different column naming than CellPhoneDB:
- CellChat: `source` = cell type, `target` = cell type, `ligand` = gene, `receptor` = gene
- CellPhoneDB: `celltype_source`, `celltype_target`, `source` = gene, `target` = gene

**Fix**: Rename CellChat columns to CellPhoneDB convention before use:
```python
ccc = ccc.rename(columns={
    'source': 'celltype_source',
    'target': 'celltype_target',
    'ligand': 'source',
    'receptor': 'target',
})
ccc['lr_means'] = ccc['prob']
```

This renaming is already built into M3 and M4 when `ccc_source == "cellchat"`.

---

## 5. Fill NaN lr_means with 0 (NOT 1.0)

**Symptom**: False strong CCC interactions appear in results, inflating scores for non-interacting cell-type pairs.

**Root cause**: Some CCC results have `NaN` for lr_means (e.g., when a gene pair wasn't tested for a cell-type pair). Filling with 1.0 creates false positives because 1.0 implies a strong interaction.

**Fix**: Always fill with 0:
```python
ccc['lr_means'] = ccc['lr_means'].fillna(0)
```

**Why 0**: Zero means "no detected interaction," which is the correct biological interpretation of a missing score. Using 1.0 would mean "maximum interaction strength," which is incorrect.

---

## 6. Memory Issues with Large Datasets

### Symptoms
- OOM (Out of Memory) during GRN building
- Kernel killed during multicell_targets()
- Slow progress with high swap usage

### Tips

**GRN building (M2)**:
- Process cell types one at a time with `--celltype` flag
- Reduce `nb_features_selected` from 100k to 50k
- Reduce `n_cpus` to limit per-process memory

**ReCoN analysis (M4)**:
- Reduce `n_jobs` (each job holds a copy of the network in memory)
- Increase `min_grn_weight` to reduce edge count
- On 373Gi instance: `n_jobs=16` keeps memory ~50%

**Multinetwork (M5)**:
- Use Parquet format (built-in) instead of CSV for large networks
- Increase `grn_score_threshold` to reduce edges

**Differential cascades (M6)**:
- `max_cascades_per_cellpair=50000` prevents combinatorial explosion
- Increase `edge_weight_threshold` to reduce cascade count
- Use `--disease ssc` or `--disease ipf` to process one disease at a time

### General memory tips:
- Use per-condition h5ad files (M1 output) instead of the full dataset
- Call `gc.collect()` between conditions
- Monitor with `htop` or `nvidia-smi`

---

## 7. Sankey Diagrams Empty or Missing Layers

**Symptom**: Sankey HTML files are created but show no links, or some layers are missing.

**Possible causes**:
1. No valid seeds found in the network -- check seed gene names match h5ad var_names
2. GRN weight filter too strict -- reduce `min_grn_weight`
3. CCC network too sparse for the focal cell type -- try different focal cell types
4. Receptor-gene network has no overlap with CCC receptors -- check species match (human only)

**Debug**: Check the data CSVs in `sankey/data/` directory. Each Sankey has companion tables showing the underlying network layers.

---

## 8. CIRCE Fails for Specific Cell Types

**Symptom**: `ValueError` or empty output from CIRCE co-accessibility computation.

**Root cause**: Cell type has too few cells in scATAC data, or the scATAC cell type mapping doesn't match.

**Fix**:
- Check `scatac_celltype_mapping` maps RNA cell types to correct ATAC peak files
- Increase `min_cells_scatac` filter
- Skip CIRCE for rare cell types (GRN will use RNA-only layers)

---

## 9. Merged CCC Has Duplicate Interactions

**Symptom**: Some ligand-receptor pairs appear multiple times in merged CCC with different scores.

**Root cause**: Both CellPhoneDB and CellChat detected the same interaction. The merge script uses percentile-rank to combine scores.

**Expected behavior**: The merged file has a `source_method` column indicating origin. Duplicates are resolved by taking the maximum percentile rank across methods. This is correct behavior -- not a bug.

---

## 10. 5-Layer GRN Score is 0 for Known TF-Gene Pairs

**Symptom**: A known TF-gene regulatory pair has score 0 in the 5-layer GRN output.

**Root cause**: The 5-layer score comes from Random Walk with Restart across all 5 layers. If the TF has low expression in that cell type/condition, the RNA layer edge is weak and RWR may not propagate sufficient signal.

**Fix**: This is expected behavior. The 5-layer approach captures regulatory potential conditioned on the current cellular state. Check the RNA-only GRN (`{ct}_{cond}_rna_network.csv`) for the raw GRNBoost2 weight.
