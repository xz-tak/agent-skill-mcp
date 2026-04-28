# ReCoN Multinetwork Parameter Guide

All parameters are fields on `ReconConfig` (defined in `scripts/config.py`).
Default values shown in parentheses. Parameters grouped by the module that primarily uses them.

---

## Paths (All Modules)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `h5ad_path` | `""` (required) | Path to input scRNA-seq h5ad file. Must contain raw counts, cell type annotations, and condition labels. |
| `output_dir` | `"results/"` | Root output directory. Each module creates subdirectories (grn/, ccc/, recon/, etc.). |
| `data_dir` | `"data/"` | Auxiliary data directory for peaks, mappings, etc. |
| `scatac_path` | `None` | Path to scATAC-seq h5ad file. Optional; enables ATAC layers in 5-layer GRN. |
| `scatac_metadata_path` | `None` | Path to scATAC metadata CSV. Required if scatac_path is set. |
| `cellchat_data_path` | `None` | Path to pre-computed CellChat RDS/CSV data. Optional; used by M3 if CellChat is in ccc_compute_methods. |

---

## Data Columns (M1)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `condition_col` | `"condition"` | Column in `adata.obs` containing condition labels (e.g., "ssc", "ipf", "normal"). |
| `celltype_col` | `"cluster"` | Column in `adata.obs` containing cell type annotations. |
| `disease_conditions` | `[]` (required) | List of disease condition names, e.g., `["ssc", "ipf"]`. |
| `normal_condition` | `"normal"` | Name of the normal/control condition. |

---

## GRN Parameters (M2)

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| `ref_genome` | `"hg38"` | hg38, hg19, mm10 | Reference genome for motif scanning. Human only recommended. |
| `circe_window` | `500000` | 100k-1M | CIRCE co-accessibility window in bp. Larger = more distal regulatory links but slower. |
| `min_cells_grn` | `20` | 10-100 | Minimum cells per cell-type x condition to build a GRN. Lower captures rare types but noisier. |
| `n_cpus` | `8` | 1-32 | CPUs per GRN build. For parallel cell-type processing, keep moderate (8) to avoid memory issues. |
| `motif_fpr` | `0.01` | 0.001-0.05 | False positive rate for TF motif scanning. Lower = stricter motif matches. |
| `tss_distance` | `10000` | 1k-100k | TSS proximity distance for ATAC-RNA links (bp). 10kb is standard. |
| `min_features_scatac` | `300` | 100-1000 | Minimum features for scATAC QC filtering. |
| `min_cells_scatac` | `50` | 20-200 | Minimum cells for scATAC QC filtering. |
| `nb_features_selected` | `100000` | 50k-500k | Number of features to select from scATAC for CIRCE. Higher = more comprehensive but slower/more memory. |
| `scatac_celltype_mapping` | `{}` | dict | Maps RNA cell types to lists of ATAC peak files, e.g., `{"Fibroblast": ["Fibro_General.bed.gz"]}`. |

---

## CCC Compute Parameters (M3)

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| `ccc_compute_methods` | `["cellphonedb", "cellchat"]` | List | Methods to compute. Options: "cellphonedb" (via LIANA+), "cellchat" (external R or pre-computed). |
| `resource_name` | `"consensus"` | consensus, cellphonedb, cellchat | LIANA+ LR resource database. "consensus" merges multiple databases. |
| `expr_prop` | `0.1` | 0.05-0.3 | Minimum expression proportion. A gene must be expressed in at least this fraction of cells in a cell type. |
| `min_lr_means` | `0.5` | 0-2.0 | Minimum LR mean score threshold for filtering CCC results. Lower = more interactions but more noise. |

---

## CCC Source Selection (M4-M7)

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `ccc_source` | `"merged"` | merged, cellphonedb, cellchat | Which CCC results to use for downstream modules. "merged" combines methods with percentile-rank normalization. Affects output directory naming (`{ccc_source}_ccc/`). |

---

## ReCoN Parameters (M4)

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| `restart_proba` | `0.6` | (0, 1) | Random walk restart probability. Lower = deeper network exploration (signal propagates further). Higher = seed-biased results. 0.6 is empirically validated. |
| `alpha` | `0.8` | [0, 1] | Weight for indirect effects in `combine_effects()`. Higher = more weight on indirect (cross-cell) effects vs direct intracellular effects. |
| `min_grn_weight` | `1.0` | 0-10 | Minimum GRNBoost2 edge weight. Filters weak regulatory edges. 1.0 keeps ~50% of edges in typical datasets. |
| `n_jobs` | `16` | 1-32 | Parallel jobs for multicell_targets(). Keep at 16 for ~50% memory usage on typical instances. |
| `extend_seeds` | `True` | bool | If True, ReCoN extends seed gene names to "gene-CellType" format across all cell types. Always True unless seeds are already formatted. |

---

## Seed Genes (M4)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `seeds` | `[]` | Inline list of seed gene names (e.g., `["TGFB1", "IL6", "TNF"]`). Mutually exclusive with seeds_file. |
| `seeds_file` | `None` | Path to text file with one gene per line. Typically DEG + CCC genes from preliminary analysis. |

---

## Sankey Visualization Seeds (M7)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `seed_categories` | `{}` | Dict mapping category names to gene lists for Sankey diagrams, e.g., `{"fibrosis": ["TGFB1", "IL6"], "inflammatory": ["TNF", "IL1B"]}`. |
| `focal_celltypes` | `[]` | Cell types to generate Sankey diagrams for (receiver cells). |
| `ligand_source_cells` | `[]` | Cell types that act as ligand sources in Sankey visualization. |

---

## Target Prediction Parameters (M8)

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| `target_genes` | `[]` | list | Genes to highlight in plots and rankings. If empty, auto-selects top 20 by mean combined score after first condition's RWR. |
| `focus_cell_types` | `[]` | list | Cell types to focus analysis on. If empty, uses all discovered cell types. |
| `seed_type` | `"receptor_activation"` | string | Label for the seed perturbation type. Informational only. |
| `prediction_output_dir` | `None` | path | Override output directory for M8. If None, uses `results/target_prediction/{ccc_source}_ccc/`. |
| `top_tfs_sankey` | `20` | >= 1 | Maximum TFs shown per Sankey diagram. |
| `top_grn_genes_sankey` | `20` | >= 1 | Maximum downstream genes per TF in Sankey (when no target_genes specified). |
| `min_rtf_weight` | `0.01` | >= 0 | Minimum NicheNet receptor-TF weight for Sankey edges. |
| `ppi_min_score` | `400` | [0, 1000] | Minimum STRING combined score for PPI partners. 0 disables PPI layer. Used by both M7 and M8. |
| `min_sankey_grn_weight` | `0.5` | >= 0 | Minimum GRN weight for edges in M8 Sankey diagrams. |

### Parameter Interactions (M8)
- `target_genes` empty + first condition completes → auto-selects top 20 genes by mean combined score
- `focus_cell_types` empty → uses all cell types discovered from GRN files
- `ppi_min_score = 0` → disables PPI 5-layer Sankey in both M7 and M8

---

## Multinetwork Parameters (M5)

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| `grn_score_threshold` | `0.001` | 0-0.1 | Minimum 5-layer GRN score to include in multinetwork. 0 = all edges, 0.001 = ~1.6k edges per cell type. |
| `module_file` | `None` | path | Optional module weights file for weighted network integration. |

---

## Differential Cascade Parameters (M6)

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| `n_permutations` | `1000` | 100-10000 | Permutations for empirical p-value computation. Higher = more precise p-values but slower. 1000 is standard. |
| `edge_weight_threshold` | `0.5` | 0-1.0 | Legacy parameter (currently unused). After percentile-rank normalization, all CCC edges are included in cascade enumeration regardless of this threshold. |
| `min_cascade_grn_weight` | `0.01` | 0-1.0 | Minimum GRN weight for edges in cascade enumeration. Lower = more comprehensive cascades. |
| `max_cascades_per_cellpair` | `50000` | 1k-500k | Maximum cascades enumerated per cell-type pair. Prevents combinatorial explosion for highly connected pairs. |
| `fdr_threshold` | `0.05` | 0.01-0.1 | False discovery rate threshold for significant cascades/edges. |

---

## Parameter Interactions and Tips

### Memory-sensitive combinations
- `n_jobs` x `n_cpus`: Total CPU load. On a 373Gi instance, `n_jobs=16` + `n_cpus=8` keeps memory ~50%.
- `max_cascades_per_cellpair` x number of cell-type pairs: With 10 cell types (100 pairs), 50k max = up to 5M cascades.
- `nb_features_selected`: 100k features with 10 cell types can use ~30GB for CIRCE.

### Quality vs speed tradeoffs
- `restart_proba`: 0.3 = deep exploration (slow), 0.9 = seed-local (fast). Default 0.6 balances both.
- `n_permutations`: 100 for quick exploratory runs, 1000+ for publication.
- `edge_weight_threshold`: 0.3 = many cascades (hours), 0.7 = few cascades (minutes).

### Score normalization chain
1. M3 computes raw CCC scores (lr_means for CellPhoneDB, prob for CellChat)
2. Merged CCC uses percentile-rank normalization across methods (0-1 scale)
3. M4 fills NaN lr_means with 0 (NOT 1.0 -- filling with 1.0 would create false strong interactions)
4. M5 normalizes GRN scores to 0-1 using max score
5. M6 uses raw edge weights for permutation testing
