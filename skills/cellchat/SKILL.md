---
name: cellchat
description: CellChat cell-cell communication analysis - end-to-end workflow for single-cell RNA-seq data. Use when users request cell-cell communication analysis, ligand-receptor interaction inference, or comparative analysis across conditions.
---

# CellChat Cell-Cell Communication Analysis

## Overview

End-to-end workflow for inferring and analyzing cell-cell communication from single-cell RNA-seq data using CellChat v2.

**Citation:** Jin et al., "CellChat for systematic analysis of cell-cell communication from single-cell transcriptomics", Nature Protocols 2024

## When to Use This Skill

Invoke this skill when users request:
- "Analyze cell-cell communication in my single-cell data"
- "Find ligand-receptor interactions between cell types"
- "Compare signaling between conditions"
- "Run CellChat analysis"
- User has h5ad or Seurat RDS file with cell type annotations

## End-to-End Workflow

### Step 1: Confirm Running Environment

**Purpose:** Verify CellChat and dependencies are installed

```bash
Rscript /home/sagemaker-user/.claude/skills/cellchat/scripts/check_environment.R
```

**Expected output:**
- R version (4.0+)
- CellChat version (v2+)
- Seurat, Matrix, reticulate status
- Python/UMAP availability

**If environment issues:** The script provides installation instructions. Key commands:
```r
# Install CellChat
devtools::install_github("jinworks/CellChat")

# Install Python UMAP (for similarity plots)
# pip install umap-learn
```

### Step 2: Validate and Convert Input File

**Purpose:** Confirm valid single-cell object, convert if needed, inspect metadata

#### For H5AD Files (MUST CONVERT FIRST)

**Convert to Seurat RDS:**
```bash
Rscript /home/sagemaker-user/.claude/skills/cellchat/scripts/convert_h5ad_to_seurat.R <file.h5ad> [output.rds] [output.log]
```

**Arguments:**
- `input.h5ad` - Path to input h5ad file (required)
- `output.rds` - Path to output RDS file (optional, defaults to same directory with .rds extension)
- `output.log` - Path to log file (optional, defaults to output path with .log extension)

**Examples:**
```bash
# Convert in same directory (creates data.rds and data.log)
Rscript convert_h5ad_to_seurat.R /path/to/data.h5ad

# Convert to specific output directory
Rscript convert_h5ad_to_seurat.R /path/to/data.h5ad /output/dir/converted.rds

# With custom log file
Rscript convert_h5ad_to_seurat.R /path/to/data.h5ad /output/dir/converted.rds /output/dir/conversion.log
```

**What the script does:**
1. Loads h5ad using anndata R package
2. Extracts counts from `raw.X` (preferred) or counts-like layers
3. Extracts normalized data from `X` (assumed log1p-normalized)
4. Copies all layers into Seurat assay layers
5. Converts obsm (PCA, UMAP, etc.) to Seurat reductions
6. Saves Seurat v5 RDS with comprehensive log

**IMPORTANT:** CellChat analysis script expects Seurat RDS input. Always convert h5ad → RDS first.

#### For RDS Files

Inspect the Seurat object:
```bash
Rscript /home/sagemaker-user/.claude/skills/cellchat/scripts/inspect_data.R <file.rds>
```

**Output shows:**
- Cell and gene counts
- Available metadata columns
- Suggested cell_type and condition columns
- Assays and reductions

**Present to user:** Show metadata columns with example values so they can choose parameters.

### Step 3: Align Parameters with User

**Purpose:** Confirm analysis configuration before running

#### 3a. Present Available Metadata Columns

The inspect script outputs metadata organized by category. Present these to the user in a **structured table format**:

| Column | Values (n) | Sample Values | Rationale |
|--------|------------|---------------|-----------|
| **Cell Type Candidates** |
| cluster | 3 | FibroblastTelo, FibroblastTropho, PeriMyofib | Best for stromal subtype analysis |
| grem1_subcluster | 5 | FibroblastTelo_CXCL14+BMP5+F3, ... | Most granular GREM1+ subtypes |
| cluster_popv | 2 | fibroblast, myofibroblast cell | Coarse ontology-based |
| **Condition Candidates** |
| comb_condition | 4 | CRC_Normal, CD_InterFibrosis, CD_Normal, CD_MatureFibrosis | Disease + fibrosis stage |
| Health | 2 | CRC, CD | Simple disease comparison |
| condition | 5 | CD_nostricture_InterFibrosis, ... | Detailed with stricture status |

**Ask user to select:**
1. Cell type column
2. Condition column
3. Which conditions to include (if not all)

#### 3b. Confirm Comparison Pairs

**IMPORTANT:** When multiple conditions are selected, ask user to confirm which **pairwise comparisons** to run:

Example prompt:
```
Selected conditions: CD_Normal, CD_InterFibrosis, CD_MatureFibrosis

Possible comparison pairs:
  1. CD_InterFibrosis vs CD_Normal (reference)
  2. CD_MatureFibrosis vs CD_Normal (reference)
  3. CD_InterFibrosis vs CD_MatureFibrosis

Which comparisons should be included? [1,2 recommended for fibrosis progression]
```

The `--reference` parameter can specify which condition is the baseline for differential analysis.

#### 3c. Script Parameters Reference

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--input` / `-i` | Yes | - | Path to Seurat RDS file |
| `--cell-type` / `-t` | Yes | - | Column with cell type labels |
| `--condition` / `-c` | No | None | Column for condition comparison |
| `--conditions` | No | All | Specific conditions (comma-separated) |
| `--reference` | No | First condition | Reference condition for differential analysis |
| `--output` / `-o` | No | cellchat_results | Output directory |
| `--species` / `-s` | No | human | human or mouse |
| `--min-cells` | No | 10 | Min cells per cell type |
| `--trim` | No | 0.1 | Trim value for truncatedMean |
| `--dry-run` | No | false | Print config without running |

### Step 4: Generate and Confirm Plan

Present analysis plan based on parameters:

```
CellChat Analysis Plan
======================
Input:        /path/to/data.rds
Cell Type:    cluster (3 types: FibroblastTelo, FibroblastTropho, PeriMyofib)
Condition:    comb_condition
Conditions:   CD_Normal, CD_InterFibrosis, CD_MatureFibrosis
Reference:    CD_Normal
Output:       ./cellchat_results/
Species:      human
Database:     CellChatDB v2 (excl. Non-protein Signaling)

Comparison Pairs (differential analysis):
  - CD_InterFibrosis vs CD_Normal
  - CD_MatureFibrosis vs CD_Normal

Cell Distribution by Condition:
  - CD_Normal: 406 cells
  - CD_InterFibrosis: 361 cells
  - CD_MatureFibrosis: 1398 cells

Analysis Steps:
1. Load RDS and extract expression matrix + metadata
2. Setup CellChatDB (human/mouse, exclude non-protein)
3. Create CellChat objects for each condition
4. For each condition:
   - Identify over-expressed genes/interactions
   - Compute communication probability (truncatedMean)
   - Filter communications (min_cells threshold)
   - Compute pathway-level signaling
   - Calculate network centrality
5. Merge conditions for comparative analysis
6. Run differential analysis for each comparison pair
7. Generate visualizations (heatmaps, circles, bubbles, etc.)
8. Generate comprehensive markdown report

Proceed? [Confirm/Modify]
```

### Step 5: Execute Analysis

1. **Copy script to output directory (for reproducibility):**
```bash
cp /home/sagemaker-user/.claude/skills/cellchat/scripts/cellchat_analysis.R <output_dir>/
```

2. **Run analysis with parameters:**
```bash
cd <output_dir>
Rscript cellchat_analysis.R \
  -i /path/to/data.rds \
  -t cluster \
  -c condition \
  -o ./cellchat_results \
  -s human
```

**Dry run (print config only):**
```bash
Rscript cellchat_analysis.R -i data.rds -t cluster -c condition --dry-run
```

**IMPORTANT - Pipeline Execution Guidelines:**
- The analysis takes 10-30 minutes depending on dataset size
- **DO NOT stop the pipeline prematurely** - it may appear idle during computation-heavy steps (e.g., computing communication probability, similarity analysis)
- Steps 4 (CellChat object creation), 13 (functional similarity), and 14 (structural similarity) are particularly time-intensive
- Only stop the pipeline if explicitly requested by the user
- Progress can be monitored by checking the output directory for new files

### Step 6: Review and Report Results

**IMPORTANT - ULTRATHINK Deep Biological Interpretation Required**

After analysis completes, Claude Code MUST apply **ULTRATHINK** reasoning to generate deep biological interpretations. This is NOT optional - all reports must include expert-level biological analysis.

#### ULTRATHINK Protocol for Result Interpretation

When interpreting CellChat outputs, Claude Code must act as a **professional biologist** with expertise in:
- **Molecular Biology**: Ligand-receptor signaling mechanisms, pathway crosstalk, protein interactions
- **Disease Biology**: Pathophysiology relevant to the SPECIFIC study context (e.g., fibrosis, cancer, inflammation, development, regeneration)
- **Computational Biology**: Network analysis interpretation, statistical significance, data integration
- **Context Adaptation**: Automatically focus on the disease, tissue, and biological processes SPECIFIC to the current analysis

**For EACH visualization and table, apply this ULTRATHINK framework:**

1. **MOLECULAR MECHANISM ANALYSIS**
   - What specific molecular interactions drive observed patterns?
   - Which ligands/receptors are responsible for altered signaling?
   - How do these interact with canonical pathways relevant to the biological context (e.g., TGF-β, WNT, NF-κB, NOTCH, Hedgehog)?
   - What is the biochemical basis for observed communication changes?

2. **BIOLOGICAL PROCESS IMPLICATIONS**
   - How do changes relate to the SPECIFIC pathophysiology being studied?
   - What stage of the biological process do patterns suggest (initiation, progression, resolution, homeostasis)?
   - Are patterns consistent with known biology from literature and databases?
   - What cellular/tissue-level consequences would these changes produce?

3. **CELL TYPE FUNCTIONAL INTERPRETATION**
   - What is the biological identity and function of each cell type IN THIS CONTEXT?
   - Why would specific cell pairs show altered communication?
   - What does autocrine vs paracrine balance indicate about tissue state?
   - How do cell type roles change across conditions/timepoints?

4. **THERAPEUTIC/TRANSLATIONAL OPPORTUNITY ASSESSMENT**
   - Which pathways represent potential intervention points?
   - Are there existing therapeutics or tool compounds that modulate identified pathways?
   - What would be predicted outcomes of pathway modulation?
   - Are identified targets relevant to the disease/biological context?

5. **EXPERIMENTAL VALIDATION PRIORITIES**
   - Which findings have highest confidence and biological plausibility?
   - What orthogonal experiments would validate predictions (in vitro, in vivo, clinical)?
   - Are there publicly available datasets for cross-validation?
   - What are the key uncertainties requiring further investigation?

**Generate the comprehensive report using:**
```bash
python /home/sagemaker-user/.claude/skills/cellchat/scripts/generate_comprehensive_report.py -i <results_dir>
```

Then apply ULTRATHINK interpretation to ALL sections of the generated report.

#### Section-by-Section ULTRATHINK Requirements

When interpreting the comprehensive report, Claude Code MUST provide deep biological insights for each section:

| Section | ULTRATHINK Requirements |
|---------|------------------------|
| **4.1 Heatmaps** | Interpret autocrine/paracrine balance; identify dominant signaling axes; explain WHY specific cell pairs communicate intensely; note asymmetric patterns (A→B vs B→A); use WebFetch for cell type communication roles |
| **4.3 Signaling Roles** | Identify communication HUBS; classify cells as senders/receivers/hubs with biological rationale; explain significance in disease context; identify therapeutic targets among hub populations |
| **4.4 Comparison Summary** | **LIST SPECIFIC PATHWAYS UP/DOWN** (not just counts); follow significance hierarchy: FDR<0.05 (high confidence) → p<0.05 (nominally significant) → CAUTION if non-significant; provide biological meaning for each altered pathway using PATHWAY_KNOWLEDGE |
| **4.5 Differential** | Explain WHY specific pairs show increased/decreased communication; identify NOVEL vs LOST communication axes; use WebFetch for known cell-cell interactions; distinguish true signaling changes from cell expansion |
| **4.6 Information Flow** | Explain WHY dominant pathways are most active; identify condition-specific pathway signatures; use WebFetch for TOP 3 pathways in disease context; note HIGH VARIABILITY pathways as key drivers |
| **4.7 Pathway Stats** | List FDR-significant pathways first, then p<0.05 if none FDR-significant; provide biological interpretation from PATHWAY_KNOWLEDGE; identify therapeutic opportunities |
| **6. Discussion** | SYNTHESIZE all findings into coherent biological narrative; use WebFetch for literature validation; generate TESTABLE HYPOTHESES; propose SPECIFIC EXPERIMENTS; identify THERAPEUTIC TARGETS and BIOMARKERS |

**For Section 4.4 specifically:**
- If FDR<0.05 pathways exist → List them with **HIGH CONFIDENCE** interpretation
- If only p<0.05 pathways → List them with **⚠️ CAUTION** note requiring validation
- If nothing significant → NOTE that changes are **subtle/distributed** and require further investigation
- ALWAYS use WebFetch to retrieve study-specific context for top altered pathways

#### Web Search Integration for Evidence-Based Interpretation

Claude Code SHOULD use **WebFetch** to retrieve relevant biological context from the MOST APPROPRIATE sources based on task context:

1. **Automatically select relevant databases/resources** which may include but are not limited to:
   - **Literature**: PubMed, Google Scholar, preprint servers for recent findings
   - **Gene/Protein**: UniProt, GeneCards, NCBI Gene, Ensembl
   - **Pathways**: KEGG, Reactome, WikiPathways, STRING
   - **Disease**: OMIM, DisGeNET, ClinVar, GWAS Catalog
   - **Drug/Target**: DrugBank, ChEMBL, OpenTargets
   - **Cell/Tissue**: Human Protein Atlas, GTEx, CellMarker
   - **Any other relevant scientific databases** based on the specific study context

2. **When to Search**:
   - Novel or unexpected pathway findings (high log2FC, significant p-value)
   - Cell type combinations not well-characterized
   - Therapeutic target validation
   - Cross-referencing with known biological mechanisms
   - Identifying existing literature support for observed patterns

3. **Context-Adaptive Search Strategy**:
   - Focus retrieval on the SPECIFIC targets, disease, and mechanisms relevant to THIS study
   - Prioritize recent publications and authoritative sources
   - Adapt search queries based on the biological context (tissue type, disease, cell populations)
   - Search for both supporting and contradicting evidence

4. **Integrate Search Findings**:
   - Cite relevant literature and database entries in interpretations
   - Compare CellChat findings with published/curated knowledge
   - Identify concordance or discordance with known biology
   - Use retrieved knowledge to strengthen biological conclusions

---

**Output Structure:**
- Root directory: `01_*`, `02_*` files (flat)
- Subfolders: `03_heatmaps/` through `18_pathway_stats/` and `objects/`

**Root Directory Files:**

| File | Description |
|------|-------------|
| `01_data_summary.csv` | Input data statistics (cells, genes, cell types, conditions) |
| `01_celltype_by_condition.csv` | Cross-tabulation of cell types by condition |
| `02_database_categories.csv` | CellChatDB signaling category counts |
| `02_database_interactions.csv` | Full L-R interaction database used |
| `cellchat_report.md` | Detailed markdown report with interpretation guide |
| `cellchat_comprehensive_report.md` | Comprehensive report with embedded images and data-driven interpretations (generated by Python script) |

**Subfolder Contents:**

| Subfolder | Plots | Unified Data Table |
|-----------|-------|-------------------|
| `03_heatmaps/` | `{cond}_count_heatmap.png`, `{cond}_weight_heatmap.png` | `heatmap_data.csv` |
| `04_circle_plots/` | `{cond}_circle_plot.png` | `circle_plot_data.csv` |
| `05_bubble_plots/` | `{cond}_{category}_bubble.png` | `bubble_plot_data.csv` |
| `06_signaling_roles/` | `{cond}_signaling_roles.png` | `signaling_roles_data.csv` |
| `07_LR_pairs/` | - | `LR_pairs_data.csv` |
| `08_pathways/` | - | `pathways_data.csv` |
| `09_compare_total/` | `compare_interactions.png` | `compare_interactions.csv` |
| `10_differential/` | `{cmp}_vs_{ref}_count.png`, `{cmp}_vs_{ref}_weight.png` | `differential_data.csv` |
| `11_info_flow/` | `information_flow.png` | `info_flow_data.csv` |
| `12_pathway_networks/` | `{cond}_{pathway}_network.png` | `pathway_networks_data.csv` |
| `13_functional_sim/` | `functional_similarity.png`, `estimateNumCluster_functional.pdf` | `functional_clusters.csv` |
| `14_structural_sim/` | `structural_similarity.png`, `estimateNumCluster_structural.pdf` | `structural_clusters.csv` |
| `15_summary/` | - | `interaction_summary.csv` |
| `16_pathway_matrix/` | - | `pathway_matrix.csv` |
| `17_pathway_flow/` | - | `pathway_info_flow.csv` |
| `18_pathway_stats/` | `{cmp}_vs_{ref}_pathway_logfc.png` | `pathway_stats_data.csv` |
| `objects/` | - | `cellchat_list.rds`, `cellchat_merged.rds` |

**Data Table Key Columns:**
- All tables include `Condition` column to filter by condition
- Differential tables include `Comparison_Pair` to filter by comparison (e.g., `CD_InterFibrosis_vs_CD_Normal`)

**Detailed Column Descriptions:**

| Data Table | Key Columns | Purpose |
|------------|-------------|---------|
| `heatmap_data.csv` | Condition, Sender, Receiver, Count, Weight | Count=interactions, Weight=probability sum |
| `circle_plot_data.csv` | Condition, CellType, CellCount | Node sizes in circle plots |
| `bubble_plot_data.csv` | Condition, annotation, source, target, ligand, receptor, prob, pval, pathway_name | L-R pairs for bubble plots; filter by annotation (Secreted/ECM/Contact) |
| `signaling_roles_data.csv` | Condition, CellType, Outgoing_Strength, Incoming_Strength | X/Y axes in scatter plots |
| `LR_pairs_data.csv` | Condition, source, target, ligand, receptor, prob, pval, pathway_name | All significant L-R pairs |
| `pathways_data.csv` | Condition, Pathway | Active pathways per condition |
| `compare_interactions.csv` | Condition, Total_Count, Total_Weight, Pct_Change_*, Is_Reference | Bar plot comparison data |
| `differential_data.csv` | Comparison_Pair, Sender, Receiver, Diff_Count, Diff_Weight, Count_Direction | UP/DOWN changes |
| `info_flow_data.csv` | Condition, Pathway, Info_Flow, Pathway_Rank | Pathway activity for information flow plots |
| `pathway_networks_data.csv` | Condition, Pathway, Sender, Receiver, Prob, Filename | Edge data for pathway network circle plots |
| `functional_clusters.csv` | Pathway, Condition, Cluster, UMAP1, UMAP2, Pathway_Condition | Functional clustering with parsed pathway/condition |
| `structural_clusters.csv` | Pathway, Condition, Cluster, UMAP1, UMAP2, Pathway_Condition | Structural clustering with parsed pathway/condition |
| `interaction_summary.csv` | Condition, Total_Interactions, Total_Strength, Num_LR_Pairs, Num_Pathways | Per-condition summary |
| `pathway_matrix.csv` | Pathway, {Condition} columns (1/0) | Pathway presence matrix |
| `pathway_info_flow.csv` | Pathway, {Condition} columns (probability sums) | Information flow values |
| `pathway_stats_data.csv` | Comparison_Pair, Pathway, Log2FC, Direction, PValue, FDR, Significant | Differential pathway stats |

---

## Step-by-Step Interpretation Guide

### 01-02: Input Data and Database (Root Directory)
**Purpose:** Validate input data quality and understand the L-R interaction database.
- Check cell counts per type (>10 recommended)
- Imbalanced composition may bias results
- Use `02_database_interactions.csv` to find specific L-R pairs

### 03: Interaction Heatmaps
**Purpose:** Identify which cell type pairs communicate most actively.
- **Rows = Senders**, **Columns = Receivers**
- **Count heatmap**: Number of L-R pairs (diversity)
- **Weight heatmap**: Sum of probabilities (strength)
- Diagonal = autocrine signaling
- **Questions:** Which cell types are signaling hubs? Autocrine vs paracrine?

### 04: Circle Plots
**Purpose:** Visualize overall communication network topology.
- **Node size** = Cell count
- **Edge width** = Interaction count
- Hub nodes = central to signaling network
- **Questions:** Overall network structure? Communication hubs?

### 05: Bubble Plots
**Purpose:** Examine specific L-R pairs driving communication.
- **Bubble size** = Communication probability
- **Bubble color** = P-value significance
- Split by category (Secreted/ECM/Contact)
- **Questions:** Which L-R pairs mediate communication? Active pathways per cell pair?

### 06: Signaling Roles
**Purpose:** Classify cell types as senders, receivers, or both.
- **X-axis** = Outgoing signal strength
- **Y-axis** = Incoming signal strength
- Upper-left = Receivers, Lower-right = Senders, Upper-right = Hubs
- **Questions:** Which cells drive vs respond to signals?

### 07: L-R Pairs Data
**Purpose:** Complete list of significant interactions for custom analysis.
- Filter by `pathway_name`, `source`, `target`
- Sort by `prob` for strongest interactions
- Compare across conditions using `Condition` column

### 08: Pathways Data
**Purpose:** List active signaling pathways per condition.
- Compare lists to find condition-specific signaling

### 09: Total Interactions Comparison
**Purpose:** Compare overall communication activity across conditions.
- Left: Interaction COUNT, Right: Interaction WEIGHT
- % change relative to reference condition
- **Questions:** Does disease increase/decrease communication?

### 10: Differential Interactions
**Purpose:** Identify altered cell type pair communication.
- **RED** = Increased in test vs reference
- **BLUE** = Decreased in test vs reference
- Filter CSV by `Comparison_Pair` for specific comparisons
- **Questions:** Which pairs show altered communication?

### 11: Information Flow
**Purpose:** Rank signaling pathways by activity.
- Stacked: Relative contribution per condition
- Grouped: Direct comparison
- **Questions:** Most active pathways? Condition-specific?

### 12: Pathway Networks
**Purpose:** Visualize individual pathway communication patterns.
- Circle plots per pathway per condition
- Compare same pathway across conditions

### 13: Functional Similarity
**Purpose:** Group pathways by sender-receiver patterns.
- UMAP embedding - close pathways = similar patterns
- Clusters = functionally related signaling programs
- `estimateNumCluster_functional.pdf` - silhouette/gap plots for optimal cluster number

### 14: Structural Similarity
**Purpose:** Group pathways by network topology.
- Different from functional: focuses on structure, not strength
- Identifies potentially redundant pathways
- `estimateNumCluster_structural.pdf` - silhouette/gap plots for optimal cluster number

### 15: Interaction Summary
**Purpose:** Quick reference statistics per condition.
- Total interactions, strength, L-R pairs, pathways

### 16: Pathway Matrix
**Purpose:** Binary presence/absence across conditions.
- 1 = active, 0 = not active
- Find condition-specific vs ubiquitous pathways

### 17: Pathway Information Flow
**Purpose:** Quantitative pathway activity values.
- Sum of communication probabilities per pathway
- Underlies 11_info_flow plots

### 18: Pathway Statistics
**Purpose:** Statistical comparison between conditions.
- **Log2FC**: Positive = UP, Negative = DOWN
- **FDR < 0.05** = Significant (marked with *)
- **Questions:** Which pathways significantly altered? Priority for validation?

---

## API Reference

The `references/` directory contains detailed CellChat API documentation. Use these when you need to customize the analysis or troubleshoot issues.

### Key Reference Files

| Reference File | Content | When to Consult |
|----------------|---------|-----------------|
| `getting_started.md` | Object creation, basic workflow, data formats | Setup issues, data format problems, understanding CellChat objects |
| `communication_inference.md` | `computeCommunProb` parameters, probability models | Adjusting sensitivity (triMean vs truncatedMean), tuning trim values |
| `comparison_analysis.md` | `mergeCellChat`, `liftCellChat`, differential analysis | Different cell compositions between conditions, multi-sample comparison |
| `visualization_tools.md` | Plot functions, customization options | Custom figures, styling, specific visualization needs |
| `database_management.md` | Custom L-R pairs, `updateCellChatDB` | Adding custom interactions, modifying database |
| `centrality_analysis.md` | Network centrality metrics | Identifying signaling hubs, flow analysis |
| `network_analysis.md` | Network topology, similarity analysis | Understanding pathway clustering, similarity metrics |
| `spatial_tutorials.md` | Spatial transcriptomics workflow | Visium, MERFISH, seqFISH data analysis |
| `r_source_code.md` | Full R source code reference | Deep debugging, understanding implementation details |

### Key Function Patterns (Quick Reference)

**Creating CellChat Object:**
```r
cellchat <- createCellChat(object = expression_matrix,
                            meta = metadata,
                            group.by = "cell_type_column")
```
- `object`: genes x cells expression matrix (log-normalized)
- `meta`: data.frame with cell annotations
- `group.by`: column name for cell type labels

**Setting Database:**
```r
CellChatDB <- CellChatDB.human  # or CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB)  # removes non-protein signaling
cellchat@DB <- CellChatDB.use
```

**Communication Inference:**
```r
cellchat <- computeCommunProb(cellchat,
                               type = "truncatedMean",  # or "triMean"
                               trim = 0.1,              # expression cutoff
                               raw.use = TRUE,          # use expression directly
                               population.size = FALSE) # don't normalize by cell count
```
- `type = "truncatedMean"`: More sensitive, better for heterogeneous data
- `type = "triMean"`: More stringent, reduces false positives
- `trim`: Higher values = more stringent filtering

**Merging for Comparison:**
```r
cellchat.merged <- mergeCellChat(list(cond1 = cc1, cond2 = cc2),
                                  add.names = c("Condition1", "Condition2"))
```

**Key Visualizations:**
```r
# Heatmap of interaction counts/weights
netVisual_heatmap(cellchat, measure = "count")  # or "weight"

# Circle plot of aggregated network
netVisual_circle(cellchat@net$count, vertex.weight = table(cellchat@idents))

# Bubble plot of L-R pairs
netVisual_bubble(cellchat, remove.isolate = TRUE)

# Differential interactions
netVisual_diffInteraction(cellchat.merged, comparison = c(1, 2))

# Information flow comparison
rankNet(cellchat.merged, mode = "comparison")
```

For detailed function signatures, parameters, and examples, see the corresponding reference file in `references/`.

---

## Common Issues and Solutions

| Issue | Cause | Solution |
|-------|-------|----------|
| No pathways detected | Expression too sparse | Use `type = "truncatedMean"` with lower `trim` (0.05) |
| Memory errors | Large dataset | Use `subsetData()` to filter to signaling genes first; reduce parallel workers |
| Condition has few cells | Sparse condition | Increase `min_cells` threshold or exclude condition |
| UMAP plots fail | Python UMAP missing | Install: `pip install umap-learn` |
| Old CellChat object | v1 → v2 migration | Run `updateCellChat(cellchat)` |
| Different cell types across conditions | Cell composition varies | Use `liftCellChat()` to align cell types |
| "No significant interactions" | Too stringent filtering | Lower `trim`, reduce `min.cells` in `filterCommunication()` |

---

## Resources

- **GitHub:** https://github.com/jinworks/CellChat
- **Nature Protocols 2024:** CellChat v2 publication
- **Nature Communications 2021:** Original CellChat publication
- **CellChat Explorer:** http://www.cellchat.org/
- **Tutorials:** https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
