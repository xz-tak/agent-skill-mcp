---
name: recon-multinetwork
description: >
  ReCoN multicellular coordination network analysis - end-to-end pipeline for
  integrating GRNs with cell-cell communication from scRNA-seq and optional
  scATAC-seq data. Features: gene regulatory network inference (5-layer HuMMuS with
  RNA, ATAC, TF-ATAC motif, ATAC-RNA TSS proximity layers), cell-cell communication
  (CellPhoneDB/CellChat), differential cascade analysis (L→R→TF→Gene), and interactive
  Sankey visualization with cascade hover tooltips. Multi-condition support (SSc/IPF
  vs normal, or any custom diseases). Use when analyzing multicellular signaling,
  CCC-GRN integration, Sankey cascade visualization, or comparing disease vs normal
  regulatory cascades. Supports CellPhoneDB, CellChat, and merged CCC methods. Human only.
---

# ReCoN Multinetwork Analysis Skill

## Overview

This skill runs the ReCoN (Regulatory Communication Network) pipeline for
multicellular coordination analysis from single-cell RNA-seq (and optional
scATAC-seq) data. It integrates gene regulatory networks (GRNs) with
cell-cell communication (CCC) to discover signaling cascades of the form:

```
Ligand:CellA -> Receptor:CellB -> TF:CellB -> Gene:CellB
```

## When to Use

Activate this skill when the user asks about:
- Multicellular signaling or coordination analysis
- Gene regulatory network (GRN) inference from scRNA-seq
- Cell-cell communication (CCC) analysis (CellPhoneDB, CellChat)
- CCC-GRN integration or multicellular target prediction
- Differential cascade analysis between conditions (disease vs normal)
- Sankey visualization of ligand-receptor-TF-gene cascades
- ReCoN, HuMMuS, or 5-layer GRN analysis
- Target prediction or treatment simulation

## Quick Start: Parameter Collection

Before running ANY module, collect required parameters from the user via
AskUserQuestion. The pipeline is configured through a single `ReconConfig`
JSON file.

### Initial Questions (ask before M1)

1. **h5ad path**: "What is the path to your scRNA-seq h5ad file?"
2. **scATAC**: "Do you have scATAC-seq data? If yes, what is the h5ad path?"
3. **Condition column**: "What is the column name in adata.obs for conditions? (e.g., 'condition', 'disease')"
4. **Cell type column**: "What is the column name for cell type annotations? (e.g., 'cluster', 'celltype')"
5. **Disease/Normal names**: "What are the disease condition names and the normal/control name? (e.g., disease=['ssc', 'ipf'], normal='normal')"
6. **Execution mode**: "Run the full pipeline (M1-M8), or specific modules? If resuming, which module to start from?"

### Required AskUserQuestion Gates (MUST ask if config field is empty)

These parameters have no sensible defaults and the modules will **fail or produce
empty output** without them. You MUST use AskUserQuestion to collect them before
running the corresponding module.

| Before Module | Config Field | Question to Ask |
|---------------|-------------|-----------------|
| **M4** (ReCoN) | `seeds` | "What seed genes should ReCoN propagate? Provide a gene list (e.g., CDH11, TGFB1) or a file path with one gene per line." |
| **M7** (Sankey) | `seed_categories` | "What seed gene categories for Sankey diagrams? Provide a dict of category name -> gene list (e.g., fibrosis: [TGFB1, COL1A1], inflammatory: [IL6, TNF])." |
| **M7** (Sankey) | `focal_celltypes` | "Which cell types should receive Sankey diagrams? (e.g., Fibroblast, Myeloid, Endothelial)" |
| **M7** (Sankey) | `ligand_source_cells` | "Which cell types act as ligand sources? (e.g., Myeloid, Epithelial, T)" |
| **M8** (Target Pred) | `seeds` | Same as M4 — reuses `config.seeds`. |
| **M8** (Target Pred) | `target_genes` (optional) | "Any specific target/readout genes to highlight? (leave empty for auto top-20, or provide e.g., IL6, MMP1, COL1A1)" |
| **M8** (Target Pred) | `focus_cell_types` (optional) | "Focus cell types for prediction analysis? (leave empty for all, or e.g., Fibroblast, Myeloid)" |
| **M8** (GSEA) | `gsea_gene_sets` | "Gene set libraries for GSEA? (default: ['MSigDB_Hallmark_2020']). Confirm or provide additional libraries." |

**Rule**: Before running M4, M7 Sankey, or M8, check if the required field is
populated in the config. If empty, use AskUserQuestion. Do NOT silently run with
empty seeds/categories — the module will either error out or skip entirely.

### Build Config JSON

After collecting answers, create a `ReconConfig` JSON:

```python
from scripts.config import ReconConfig

config = ReconConfig(
    h5ad_path="/path/to/data.h5ad",
    output_dir="results/",
    condition_col="condition",
    celltype_col="cluster",
    disease_conditions=["ssc", "ipf"],
    normal_condition="normal",
    # ... other params with defaults
)
config.to_json(Path("config.json"))
```

See `references/parameter_guide.md` for all ~40 parameters with defaults and ranges.

---

## Pipeline Modules

### Module Dependencies

```
M1 (Data Prep) -> M2 (GRN) -> M3 (CCC) -> M4 (ReCoN) -> M5 (Multinetwork)
                                    |                          |
                                    +-------> M6 (Differential) <---+
                                    |                |
                                    |          M7 (Visualization)
                                    |
                                    +-------> M8 (Target Prediction)
```

M1-M3 can run independently. M4 requires M2+M3. M5 requires M2+M3. M6 requires M2+M3. M7 requires M4+M6. M8 requires M2+M3.

### M1: Data Preparation

**Purpose**: Load h5ad, subset by condition, create per-condition files, download ATAC peaks.

**JIT Questions**:
- "Confirm the cell type and condition column names look correct?" (show unique values)
- If scATAC: "Provide the cell-type to ATAC peak file mapping"

**CLI**: `python scripts/m1_data_prep.py --config config.json`

**Outputs**: `results/adata_{condition}.h5ad`, cell type mapping

### M2: GRN Pipeline

**Purpose**: Build 5-layer HuMMuS GRNs per cell type per condition.

**JIT Questions**:
- "How many CPUs to use per cell type? (default: 8)"
- If scATAC available: "Confirm scATAC cell type mapping?"
- "Process all cell types or a specific subset?"

**CLI**: `python scripts/m2_grn_pipeline.py --config config.json`
Or per cell type: `python scripts/m2_grn_pipeline.py --config config.json --celltype fibroblast`

**Outputs**: `results/grn/{ct}_{cond}_rna_network.csv`, `*_5layer_grn.csv`

**Tip**: Run cell types in parallel via tmux for 10x speedup.

### M3: CCC Analysis

**Purpose**: Run CellPhoneDB (via LIANA+) and/or CellChat, merge results.

**JIT Questions**:
- "Which CCC methods: cellphonedb, cellchat, or both? (default: both)"
- "Minimum expression proportion? (default: 0.1)"
- If CellChat: "Path to pre-computed CellChat data?"

**CLI**: `python scripts/m3_ccc_analysis.py --config config.json`

**Outputs**: `results/ccc/{method}/{cond}_ccc.csv`, `results/ccc/merged/{cond}_ccc.csv`

### M4: ReCoN Analysis

**Purpose**: Run multicell_targets() to compute direct/indirect effects across cell types.

**JIT Questions**:
- "Which CCC source for downstream analysis: merged, cellphonedb, or cellchat? (default: merged)"
- "Seed genes: provide a list, a file path, or use DEG-derived seeds?"
- "Restart probability? (default: 0.6, lower = deeper exploration)"

**CLI**: `python scripts/m4_recon_analysis.py --config config.json`

**Outputs**: `results/recon/{ccc}_ccc/{cond}_{direct,indirect,combined}_effects.csv`

**Critical bugs to watch for**: See `references/troubleshooting.md` items 1-3.

### M5: Multinetwork

**Purpose**: Merge CCC + GRN + Receptor-Gene into unified network with `GENE:CellType` nodes.

**JIT Questions**:
- "GRN score threshold? (default: 0.001, 0 = all edges)"
- "Module weights file available?"

**CLI**: `python scripts/m5_multinetwork.py --config config.json`

**Outputs**: `results/multinetwork/{ccc}_ccc/{cond}_multinetwork.parquet`

### M6: Differential Cascades

**Purpose**: Enumerate L->R->TF->Gene cascades, permutation-test disease vs normal.

**JIT Questions**:
- "Number of permutations? (default: 1000, use 100 for quick test)"
- "Edge weight threshold? (default: 0.5, lower = more cascades but slower)"
- "Which disease to compare: ssc, ipf, or both?"

**CLI**: `python scripts/m6_differential.py --config config.json`

**Outputs**: `results/differential_cascades/{ccc}_ccc/{disease}/cascade_results.csv`

### M7: Visualization

**Purpose**: Generate heatmaps, comparison plots, and interactive Sankey diagrams.

**JIT Questions**:
- "Which cell types to focus Sankey diagrams on? (e.g., Fibroblast, Pericytes)"
- "Which cell types are ligand sources? (e.g., Myeloid, Epithelial)"
- "Seed gene categories for Sankey? (provide dict of category->gene lists)"

**CLI**: `python scripts/m7_visualization.py --config config.json`

**Outputs**: `results/figures/{ccc}_ccc/*.png`, `results/figures/{ccc}_ccc/sankey/*.html`
- PPI 5-layer Sankey diagrams (if ppi_min_score > 0 and seed_categories configured)

### M8: Target Prediction

**Purpose**: Predict cell-type-specific effects of seed gene perturbation across conditions. Generalizes CDH11 treatment pipeline into reusable module.

**JIT Questions**:
- "Target genes to highlight? (provide list, or leave empty for auto top-20)"
- "Focus cell types for analysis? (provide list, or all cell types)"
- "Seed type: receptor_activation or other?"

**CLI**: `python scripts/m8_target_prediction.py --config config.json`
Or standalone: `python scripts/m8_target_prediction.py --config config.json --conditions ssc ipf --sankey-only --skip-plots`

**Outputs**: `results/target_prediction/{ccc}_ccc/{cond}/combined_effects.csv`, scatter plots, Sankey plots, differential heatmaps

**Standalone flags**: `--conditions`, `--sankey-only`, `--skip-plots`, `--output-dir`

#### M8 GSEA Enrichment

**Purpose**: Run Gene Set Enrichment Analysis (prerank) on gene rankings per cell type per condition. Auto-detects cell types and custom target gene sets from gene_rankings.csv.

**JIT Questions**:
- "Gene set libraries for GSEA? (default: ['MSigDB_Hallmark_2020'])"
  Show auto-detected custom set name and gene count for confirmation.

**Outputs per condition**: `target_gsea/` subdirectory containing:
- `{celltype}/` — enrichment plots for FDR<0.05 pathways + gseapy reports
- `gsea_all_celltypes.csv` — concatenated NES/FDR results
- `gsea_heatmap.png` — clustered NES heatmap (significant pathways only)

**Cross-condition output**:
- `gsea_cross_condition_heatmap.png` — NES comparison across all conditions

**Heatmap colormap logic**:
- Both positive & negative NES -> `bwr` centered at 0
- All positive NES -> `YlOrRd` (data range)
- All negative NES -> `YlGnBu` (data range)

**Config parameters**: `gsea_gene_sets`, `gsea_min_size`, `gsea_max_size`, `gsea_permutations`, `gsea_fdr_threshold`

---

## Full Pipeline Run

```bash
# Run everything M1-M8
python scripts/run_pipeline.py --config config.json

# Resume from M4
python scripts/run_pipeline.py --config config.json --start-from 4

# Run only M4-M8
python scripts/run_pipeline.py --config config.json --start-from 4 --end-at 8
```

---

## CLI Reference

All modules accept `--config config.json`. Additional per-module flags:

| Module | Extra Flags |
|--------|-------------|
| M2 | `--celltype <name>` (process single cell type) |
| M4 | `--ccc-source {merged,cellphonedb,cellchat}` |
| M5 | `--ccc-source {merged,cellphonedb,cellchat}` |
| M6 | `--ccc-source`, `--disease {ssc,ipf,both}`, `--n-permutations`, `--n-jobs` |
| M7 | `--ccc-source {merged,cellphonedb,cellchat}` |
| M8 | `--conditions`, `--sankey-only`, `--skip-plots`, `--output-dir`, `--gsea-gene-sets`, `--gsea-min-size`, `--gsea-max-size`, `--gsea-permutations`, `--gsea-fdr-threshold` |
| run_pipeline | `--start-from N`, `--end-at N` (N=1-8) |

---

## Reference Documentation

- **Parameter Guide**: `references/parameter_guide.md` -- All ~40 parameters with defaults, ranges, impact
- **Data Formats**: `references/data_formats.md` -- Input/output schemas per module
- **Troubleshooting**: `references/troubleshooting.md` -- Known bugs and workarounds
- **SSc Lung Atlas Example**: `references/ssc_lung_atlas_example.md` -- Complete walkthrough with config, commands, expected outputs
