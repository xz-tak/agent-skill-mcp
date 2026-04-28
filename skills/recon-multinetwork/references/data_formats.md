# ReCoN Multinetwork Data Formats

Input/output schemas per module. Column names are exact.

---

## M1: Data Preparation

### Input
- **h5ad file**: AnnData with `.X` (raw counts), `.obs[condition_col]`, `.obs[celltype_col]`
- **scATAC h5ad** (optional): AnnData with peaks as features
- **scATAC metadata** (optional): CSV with cell barcodes and annotations

### Output
| File | Columns / Format | Description |
|------|------------------|-------------|
| `results/adata_{condition}.h5ad` | Per-condition AnnData subsets | One file per condition |
| `data/celltype_atac_mapping.csv` | celltype, peak_files | Maps RNA cell types to ATAC peak files |

---

## M2: GRN Pipeline

### Input
- Per-condition h5ad files from M1
- scATAC h5ad (optional)
- CIRCE co-accessibility (optional)

### Output
| File | Columns | Description |
|------|---------|-------------|
| `results/grn/{ct}_{cond}_rna_network.csv` | source, target, weight | GRNBoost2 TF->gene edges. Weight = importance score. |
| `results/grn/{ct}_{cond}_tf_network.csv` | source, target, weight | TF-TF correlation network |
| `results/grn/{ct}_tf_atac_links.csv` | source, target, weight | TF motif -> ATAC peak links |
| `results/grn/{ct}_atac_rna_links.csv` | source, target, weight | ATAC peak -> gene (TSS proximity) |
| `results/grn/{ct}_{cond}_5layer_grn.csv` | seed, target, score, path | Integrated 5-layer GRN via RWR propagation |
| `results/grn/lung_scatac_prep.h5ad` | AnnData | Preprocessed scATAC reference |

Key: `{ct}` = cell type (lowercase, spaces as underscores), `{cond}` = condition name.

---

## M3: CCC Analysis

### Input
- Per-condition h5ad files from M1

### CellPhoneDB Output
| File | Columns | Description |
|------|---------|-------------|
| `results/ccc/cellphonedb/{cond}_ccc.csv` | source, target, celltype_source, celltype_target, lr_means, pval, ligand_complex, receptor_complex | Raw CellPhoneDB results via LIANA+ |

### CellChat Output
| File | Columns | Description |
|------|---------|-------------|
| `results/ccc/cellchat/{cond}_ccc.csv` | source, target, ligand, receptor, prob, pval, pathway_name, annotation | CellChat results. NOTE: column names differ from CellPhoneDB (source/target = cell types, ligand/receptor = genes) |

### Merged Output
| File | Columns | Description |
|------|---------|-------------|
| `results/ccc/merged/{cond}_ccc.csv` | source, target, celltype_source, celltype_target, lr_means, pct_rank_overall, source_method, lr_means_cellphonedb, prob_cellchat | Merged CCC with percentile-rank normalized lr_means |

### CellChat Column Standardization
When loading CellChat data downstream, columns are renamed:
```
CellChat original   ->  Standardized
source              ->  celltype_source
target              ->  celltype_target
ligand              ->  source
receptor            ->  target
prob                ->  lr_means
```

---

## M4: ReCoN Analysis

### Input
- GRN files from M2 (`{ct}_{cond}_rna_network.csv`)
- CCC files from M3 (based on `ccc_source` setting)
- NicheNet receptor-gene network (built-in to ReCoN package)
- Seeds file or seed gene list

### Output
| File | Columns | Description |
|------|---------|-------------|
| `results/recon/{ccc}_ccc/{cond}_direct_effects.csv` | Index: genes, Columns: cell types | Direct intracellular effects from multicell_targets() |
| `results/recon/{ccc}_ccc/{cond}_indirect_effects.csv` | Index: genes, Columns: cell types | Indirect cross-cell effects |
| `results/recon/{ccc}_ccc/{cond}_combined_effects.csv` | Index: genes, Columns: cell types | alpha * indirect + (1-alpha) * direct |
| `results/recon/{ccc}_ccc/comparison_{disease}_vs_normal.csv` | celltype, mean_fibrotic, mean_normal, mean_log2_fc, correlation, n_upregulated, n_downregulated | Per-cell-type condition comparison |
| `results/recon/{ccc}_ccc/recon_stats.json` | JSON | Run statistics per condition |

**Critical**: Direct/indirect DataFrames must be `.copy()`'d before `combine_effects()` -- that function modifies inputs in-place.

---

## M5: Multinetwork

### Input
- 5-layer GRN files from M2 (`{ct}_{cond}_5layer_grn.csv`)
- CCC files from M3
- NicheNet receptor-gene network
- Per-condition h5ad (optional, for expression overlay)
- Module weights file (optional)

### Output
| File | Columns | Description |
|------|---------|-------------|
| `results/multinetwork/{ccc}_ccc/{cond}_multinetwork.parquet` | source_node, target_node, weight, edge_type, interaction | Unified multinetwork with standardized node format `GENE:CellType` |
| `results/multinetwork/{ccc}_ccc/{cond}_node_metadata.parquet` | node, gene, celltype, layers, degree | Node-level metadata |
| `results/multinetwork/{ccc}_ccc/network_stats.json` | JSON | Edge counts by type and condition |

Node format in multinetwork: `GENE:CellType` (colon-separated).

---

## M6: Differential Cascades

### Input
- CCC files from M3
- GRN files from M2 (`{ct}_{cond}_rna_network.csv`)
- NicheNet receptor-gene network

### Cascade Structure
```
Ligand:Cell_A -> Receptor:Cell_B -> TF:Cell_B -> Gene:Cell_B
    (CCC)           (Receptor-TF)        (GRN)
```

### Output
| File | Columns | Description |
|------|---------|-------------|
| `results/differential_cascades/{ccc}_ccc/{disease}/cascade_results.csv` | cascade_id, ligand, receptor, tf, gene, cell_source, cell_target, score_disease, score_normal, diff, pval, padj, ci_low, ci_high | Full cascade statistics |
| `results/differential_cascades/{ccc}_ccc/{disease}/cellpair_results.csv` | cell_source, cell_target, n_cascades, n_significant, fisher_pval, mean_diff | Cell-pair aggregation with Fisher's combined p-value |
| `results/differential_cascades/{ccc}_ccc/{disease}/edge_results.csv` | edge_key, edge_type, weight_disease, weight_normal, diff, pval, padj | Edge-level differential statistics |
| `results/differential_cascades/{ccc}_ccc/summary_stats.json` | JSON | Combined summary across disease comparisons |

Edge key format: `SRC_GENE::SRC_CELLTYPE|TGT_GENE::TGT_CELLTYPE|LAYER_TYPE`
where LAYER_TYPE is `ccc`, `rtf`, or `grn`.

### P-value Computation

Two-tailed t-test with adaptive degrees of freedom:
- **Formula**: `df = 6 / kurtosis + 4` (for excess kurtosis > 0, bounded to [4.5, 100])
- **Rationale**: Cascade scores often have heavy-tailed distributions (high kurtosis). Standard normal distribution assumptions break down. Estimating df from kurtosis adapts the t-distribution to match the actual data distribution.
- **When excess kurtosis ≤ 0**: Use `df = 4.5` (minimum)
- **When excess kurtosis > 25**: Cap at `df = 100` (maximum reasonable freedom)

---

## M7: Visualization

### Input
- ReCoN results from M4 (combined effects CSVs)
- CCC files from M3
- GRN files from M2
- Differential cascade results from M6 (for Sankey hover tooltips)

### Output
| File | Format | Description |
|------|--------|-------------|
| `results/figures/{ccc}_ccc/coordination_heatmap_{cond}.png` | PNG | Cell type correlation matrix |
| `results/figures/{ccc}_ccc/top_genes_heatmap_{cond}.png` | PNG | Top 30 affected genes |
| `results/figures/{ccc}_ccc/condition_comparison_{disease}_vs_normal.png` | PNG | Scatter plots per cell type |
| `results/figures/{ccc}_ccc/differential_effects_{disease}_vs_normal.png` | PNG | Volcano-style MA plots |
| `results/figures/{ccc}_ccc/ccc_network_{cond}.png` | PNG | CCC strength heatmap |
| `results/figures/{ccc}_ccc/analysis_summary_{disease}.png` | PNG | 6-panel summary figure |
| `results/figures/{ccc}_ccc/fibrosis_markers_{disease}.png` | PNG | Fibrosis gene bar plot |
| `results/figures/{ccc}_ccc/sankey/*.html` | HTML (Plotly) | Interactive Sankey diagrams |
| `results/figures/{ccc}_ccc/sankey/data/*.csv` | CSV | Underlying data for each Sankey |

### Sankey Diagram Types (per focal cell type x condition)
- `sankey_intracell_{cond}_{ct}.html` -- 3-layer: Receptor -> TF -> Gene
- `sankey_ligand_{cond}_{ct}.html` -- 4-layer: Ligand -> Receptor -> TF -> Gene
- `sankey_intercell_{cond}_{ct}.html` -- 6-layer: Upstream Receptor -> Upstream TF -> Ligand -> Receptor -> TF -> Gene

---

## Seeds File Format

Plain text, one gene per line:
```
TGFB1
IL6
TNF
PDGFA
CCL2
```

No headers. Gene names must match the h5ad `.var_names` (typically HGNC symbols for human).

---

## M8: Target Prediction

### Input
- GRN files from M2 (`{ct}_{cond}_rna_network.csv`)
- CCC files from M3 (based on `ccc_source` setting)
- NicheNet receptor-gene network (built-in to ReCoN package)
- Seeds from config (`seeds` or `seeds_file`)

### Per-Condition Output
| File | Columns | Description |
|------|---------|-------------|
| `results/target_prediction/{ccc}_ccc/{cond}/direct_effects.csv` | Index: genes, Columns: cell types | Direct intracellular effects |
| `results/target_prediction/{ccc}_ccc/{cond}/indirect_effects.csv` | Index: genes, Columns: cell types | Indirect cross-cell effects |
| `results/target_prediction/{ccc}_ccc/{cond}/combined_effects.csv` | Index: genes, Columns: cell types | alpha * indirect + (1-alpha) * direct |
| `results/target_prediction/{ccc}_ccc/{cond}/gene_rankings.csv` | gene, {ct}_rank..., mean_rank, is_target | Rankings across focus cell types |
| `results/target_prediction/{ccc}_ccc/{cond}/top_genes_per_celltype.csv` | gene, celltype, rank, score, is_target | All genes ranked per cell type |

### Cross-Condition Output
| File | Format | Description |
|------|--------|-------------|
| `results/target_prediction/{ccc}_ccc/differential/{disease}_vs_normal_all.csv` | CSV | All genes x cell types with log2FC |
| `results/target_prediction/{ccc}_ccc/differential/{disease}_vs_normal_pivot.csv` | CSV | Pivoted log2FC matrix |
| `results/target_prediction/{ccc}_ccc/differential/heatmap_target_genes_{disease}_vs_normal.html` | HTML (Plotly) | Target gene log2FC heatmap |
| `results/target_prediction/{ccc}_ccc/{cond}/scatter_*.html` | HTML (Plotly) | Pairwise scatter plots |
| `results/target_prediction/{ccc}_ccc/sankey/sankey_{cond}_{ct}.html` | HTML (Plotly) | 4-layer treatment Sankey |
| `results/target_prediction/{ccc}_ccc/setup.json` | JSON | Run configuration and statistics |

---

## M7: PPI 5-Layer Sankey Output (when ppi_min_score > 0)

| File | Format | Description |
|------|--------|-------------|
| `results/figures/{ccc}_ccc/sankey/sankey_ppi_{cond}_{ct}.html` | HTML (Plotly) | 5-layer PPI Sankey: Ligand -> Seed -> PPI Partner -> TF -> Gene |
| `results/figures/{ccc}_ccc/sankey/ppi_partners.csv` | CSV | Cached STRING PPI partners with combined_score and is_tf flag |
