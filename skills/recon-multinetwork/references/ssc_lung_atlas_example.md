# SSc Lung Atlas: Complete Walkthrough

End-to-end example using the SSc (Systemic Sclerosis) Lung Atlas dataset with
SSc, IPF, and Normal conditions across 10 cell types.

---

## Environment Setup

```bash
# Create conda environment
conda create -n recon python=3.10 -y
conda activate recon

# Core packages
pip install scanpy anndata pandas numpy scipy
pip install liana  # CellPhoneDB via LIANA+
pip install recon  # ReCoN multicellular analysis (includes HuMMuS, multixrank)
pip install plotly adjustText seaborn matplotlib
pip install joblib tqdm pyarrow  # For parallel processing and parquet

# Optional for scATAC
pip install episcanpy pybiomart
```

---

## Data Paths

```
# scRNA-seq (integrated atlas, ~439K cells)
/home/sagemaker-user/data_process/ssc/ssc_lung_atlas/integration_scanpy_11202025.h5ad

# scATAC-seq peaks (Zhang 2021)
s3://tec-rnd-sci-dev-gi2/gi2-xz/insights/project/Zhang_2021_scatac/peaks/

# CellChat pre-computed (from R)
/home/sagemaker-user/data_process/ssc/ssc_lung_atlas/cellchat/

# Module weights (for multinetwork)
/home/sagemaker-user/data_process/ssc/ssc_lung_atlas/lung_moduleweights.txt

# Output base
/home/sagemaker-user/claude_code/recon/ssc_lung_atlas/results/
```

---

## ATAC Cell Type Mapping Guide

Before setting `scatac_celltype_mapping`, understand the mapping structure:

**Format**: `{"RNA_CellType": ["peak_file1.bed.gz", "peak_file2.bed.gz", ...]}`

Each RNA cell type is mapped to one or more ATAC peak files. For the SSc dataset, this mapping comes from Zhang et al. 2021 scATAC data:

| RNA Cell Type | ATAC Peak Files | Source |
|---|---|---|
| Myeloid | Macrophage_Gen_or_Alv.bed.gz, Macrophage_General.bed.gz, Mast.bed.gz | Zhang 2021 |
| Epithelial | Alveolar_Type_1.bed.gz, Alveolar_Type_2.bed.gz, Club.bed.gz, Ciliated.bed.gz, Airway_Goblet.bed.gz | Zhang 2021 |
| Fibroblast | Fibro_General.bed.gz, Fibro_Muscle.bed.gz | Zhang 2021 |
| Endothelial | Endothelial_General_1.bed.gz, Endothelial_General_2.bed.gz, Alveolar_Cap_Endo.bed.gz, Lymphatic.bed.gz | Zhang 2021 |

**For other datasets:**
1. Identify cell type-specific peak files from your ATAC source
2. Create a mapping that links RNA cell types (from your h5ad) to corresponding ATAC peak BED files
3. Set `scatac_path` to the directory containing these .bed.gz files
4. If no scATAC data available, set `scatac_path: null` to build RNA-only GRNs

---

## Full Config JSON

```json
{
  "h5ad_path": "/home/sagemaker-user/data_process/ssc/ssc_lung_atlas/integration_scanpy_11202025.h5ad",
  "output_dir": "/home/sagemaker-user/claude_code/recon/ssc_lung_atlas/results",
  "data_dir": "/home/sagemaker-user/claude_code/recon/ssc_lung_atlas/data",
  "scatac_path": "/home/sagemaker-user/claude_code/recon/ssc_lung_atlas/results/grn/lung_scatac_prep.h5ad",
  "scatac_metadata_path": null,
  "cellchat_data_path": "/home/sagemaker-user/data_process/ssc/ssc_lung_atlas/cellchat",

  "condition_col": "condition",
  "celltype_col": "cluster_l2",
  "disease_conditions": ["ssc", "ipf"],
  "normal_condition": "normal",

  "ref_genome": "hg38",
  "circe_window": 500000,
  "min_cells_grn": 20,
  "n_cpus": 8,
  "motif_fpr": 0.01,
  "tss_distance": 10000,
  "min_features_scatac": 300,
  "min_cells_scatac": 50,
  "nb_features_selected": 100000,
  "scatac_celltype_mapping": {
    "Myeloid": ["Macrophage_Gen_or_Alv.bed.gz", "Macrophage_General.bed.gz", "Mast.bed.gz"],
    "Epithelial": ["Alveolar_Type_1.bed.gz", "Alveolar_Type_2.bed.gz", "Club.bed.gz", "Cilliated.bed.gz", "Airway_Goblet.bed.gz"],
    "Fibroblast": ["Fibro_General.bed.gz", "Fibro_Muscle.bed.gz"],
    "Endothelial": ["Endothelial_General_1.bed.gz", "Endothelial_General_2.bed.gz", "Alveolar_Cap_Endo.bed.gz", "Lymphatic.bed.gz"],
    "T": ["CD4_T.bed.gz", "CD8_T.bed.gz"],
    "B": ["B_Cell.bed.gz", "Plasma.bed.gz"],
    "NK": ["NK.bed.gz"],
    "Smooth muscle": ["Smooth_Muscle.bed.gz"],
    "Pericytes": ["Pericyte.bed.gz"],
    "Mesothelial": ["Mesothelial.bed.gz"]
  },

  "ccc_compute_methods": ["cellphonedb", "cellchat"],
  "resource_name": "consensus",
  "expr_prop": 0.1,
  "min_lr_means": 0.5,

  "ccc_source": "merged",

  "restart_proba": 0.6,
  "alpha": 0.8,
  "min_grn_weight": 1.0,
  "n_jobs": 16,
  "extend_seeds": true,

  "seeds": [],
  "seeds_file": null,

  "seed_categories": {
    "fibrosis": ["TGFB1", "TGFB2", "TGFB3", "IL6", "IL1B", "TNF", "PDGFA", "PDGFB", "CTGF", "CCL2", "SPP1", "IL11", "COL1A1", "COL3A1", "FN1"],
    "cgas_sting": ["CGAS", "MB21D1", "STING1", "TMEM173", "TBK1", "IRF3", "IRF7", "IFNB1", "CXCL10", "ISG15"],
    "nlrp3": ["NLRP3", "PYCARD", "CASP1", "IL1B", "IL18", "GSDMD", "TXNIP", "NEK7"],
    "taci_baff_april": ["TNFRSF13B", "TNFRSF13C", "TNFRSF17", "TNFSF13B", "TNFSF13", "NFKB1", "RELA", "NFKB2"],
    "irf5_associated": ["IRF5", "IRF8", "STAT2", "STAT3", "NFKB1", "RELA", "STAT1", "IRF3", "IRF4", "IRF9", "IFNB1", "BIRC2", "BIRC3"],
    "grem1_associated": ["GREM1", "GREM2", "SHH", "IHH", "DHH", "BMPR1A", "TGFBR3", "BMP2", "BMP4", "BMP7", "NOTCH1", "WNT5A"],
    "pcolce_associated": ["PCOLCE", "PCOLCE2", "BMP1", "TLL1", "TLL2", "COL1A1", "COL1A2", "COL3A1", "LOX", "LOXL2"],
    "il23_associated": ["IL23A", "IL12B", "IL23R", "IL17A", "IL17F", "RORC", "STAT3", "FOXP3", "JAK2", "TYK2"],
    "glp2r_associated": ["GLP2R", "GCG", "TCF4", "TCF12", "PPARD", "PPARA", "PPARG", "ATF3", "FOS", "CREB1", "HIF1A", "EGR1"]
  },
  "focal_celltypes": ["Fibroblast", "Pericytes", "Smooth muscle", "B", "Myeloid"],
  "ligand_source_cells": ["Myeloid", "Epithelial", "Endothelial", "Fibroblast", "B"],

  "n_permutations": 1000,
  "edge_weight_threshold": 0.5,
  "min_cascade_grn_weight": 0.01,
  "max_cascades_per_cellpair": 50000,
  "fdr_threshold": 0.05,

  "grn_score_threshold": 0.001,
  "module_file": "/home/sagemaker-user/data_process/ssc/ssc_lung_atlas/lung_moduleweights.txt"
}
```

---

## Module-by-Module Run Commands

### M1: Data Preparation
```bash
cd /home/sagemaker-user/claude_code/recon/ssc_lung_atlas/scripts
python 01_data_prep.py
```
Subsets atlas to 3 conditions, creates per-condition h5ad files, downloads ATAC peaks from S3.

### M2: GRN Pipeline (parallel by cell type)
```bash
# Run all 10 cell types in parallel (tmux recommended)
for ct in myeloid epithelial fibroblast endothelial t b nk smooth_muscle pericytes mesothelial; do
    python 02c_5layer_grn.py --celltype $ct &
done
wait
```
Builds 30 GRNs (10 cell types x 3 conditions). scATAC processing (02a, 02b) runs first if needed.

### M3: CCC Analysis
```bash
python 03_ccc_analysis.py
```
Runs CellPhoneDB via LIANA+. CellChat data loaded from pre-computed path. Merges results.

### M4: ReCoN Analysis
```bash
python 04_recon_analysis.py --ccc-source merged
```
Runs `multicell_targets()` for each condition using merged CCC + merged GRNs.

### M5: Multinetwork
```bash
python generate_multinetwork.py --ccc-source merged
```
Merges CCC + GRN + Receptor-Gene into unified parquet multinetwork.

### M6: Differential Cascades
```bash
python 05_differential_cascades.py --ccc-source merged --disease both
```
Enumerates L->R->TF->Gene cascades, permutation-tests disease vs normal.

### M7: Visualization
```bash
python 06_visualization.py --ccc-source merged
```
Generates heatmaps, comparison plots, and interactive Sankey diagrams.

---

## Expected Outputs

### Cell Counts (approximate)
- Total: ~439,000 cells
- SSc: ~180,000
- IPF: ~100,000
- Normal: ~159,000
- 10 cell types: Myeloid, Epithelial, Fibroblast, Endothelial, T, B, NK, Smooth muscle, Pericytes, Mesothelial

### GRN Counts
- 30 RNA networks (10 cell types x 3 conditions)
- ~5,000-50,000 edges per cell-type/condition after weight > 1.0 filter
- 30 five-layer GRNs with integrated scores

### CCC Counts
- CellPhoneDB: ~15,000-25,000 interactions per condition
- CellChat: ~10,000-20,000 interactions per condition
- Merged: ~20,000-30,000 unique interactions per condition
- SSc: ~20,051 merged interactions
- IPF: ~14,513 merged interactions (sparser -- triggers multixrank_patch)
- Normal: ~18,000 merged interactions

### Sankey Count
- 5 focal cell types x 3 conditions x 3 diagram types = 45 HTML files
- Plus companion CSV data tables in sankey/data/

---

## Key Findings (SSc Lung Atlas)

### Most Disrupted Cell Types
**Pericytes** showed the most disrupted multicellular coordination in SSc, with:
- Lowest correlation between SSc and Normal combined effects
- Highest number of upregulated genes (log2FC > 1)
- Key signaling through PDGFB/PDGFRB axis

### Top Cascade Pathways
1. **Fibrosis cascades**: TGFB1->TGFBR1->SMAD3->COL1A1 (Myeloid->Fibroblast)
2. **Inflammatory cascades**: IL6->IL6R->STAT3->acute phase genes (Myeloid->various)
3. **Vascular remodeling**: PDGFB->PDGFRB->downstream (Endothelial->Pericytes)

### SSc vs IPF Comparison
- SSc shows stronger inflammatory component (higher IL6, TNF cascades)
- IPF shows stronger fibrotic component (higher TGFB pathway)
- Both share disrupted Pericyte and Fibroblast coordination

---

## SANKEY_SEEDS Categories

Nine seed categories used for Sankey visualization, covering major signaling pathways:

| Category | Genes | Rationale |
|----------|-------|-----------|
| **fibrosis** | TGFB1, TGFB2, TGFB3, IL6, IL1B, TNF, PDGFA, PDGFB, CTGF, CCL2, SPP1, IL11, COL1A1, COL3A1, FN1 | Core fibrosis signaling + ECM components |
| **cgas_sting** | CGAS, MB21D1, STING1, TMEM173, TBK1, IRF3, IRF7, IFNB1, CXCL10, ISG15 | cGAS-STING innate immune pathway |
| **nlrp3** | NLRP3, PYCARD, CASP1, IL1B, IL18, GSDMD, TXNIP, NEK7 | NLRP3 inflammasome |
| **taci_baff_april** | TNFRSF13B, TNFRSF13C, TNFRSF17, TNFSF13B, TNFSF13, NFKB1, RELA, NFKB2 | TACI/BAFF/APRIL B-cell signaling |
| **irf5_associated** | IRF5, IRF8, STAT2, STAT3, NFKB1, RELA, STAT1, IRF3, IRF4, IRF9, IFNB1, BIRC2, BIRC3 | IRF5 transcription factor network |
| **grem1_associated** | GREM1, GREM2, SHH, IHH, DHH, BMPR1A, TGFBR3, BMP2, BMP4, BMP7, NOTCH1, WNT5A | GREM1 BMP antagonism + Hedgehog/Notch |
| **pcolce_associated** | PCOLCE, PCOLCE2, BMP1, TLL1, TLL2, COL1A1, COL1A2, COL3A1, LOX, LOXL2 | Procollagen processing + crosslinking |
| **il23_associated** | IL23A, IL12B, IL23R, IL17A, IL17F, RORC, STAT3, FOXP3, JAK2, TYK2 | IL-23/IL-17 axis + Th17 differentiation |
| **glp2r_associated** | GLP2R, GCG, TCF4, TCF12, PPARD, PPARA, PPARG, ATF3, FOS, CREB1, HIF1A, EGR1 | GLP2R metabolic/epithelial repair signaling |

Total unique seed genes across all categories: ~90.
