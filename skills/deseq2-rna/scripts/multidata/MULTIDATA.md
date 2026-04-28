# Multi-Dataset Comparison Module

Independent module for comparing DESeq2 results with external datasets. This module is **detached from the main 10-step pipeline** and invoked on-demand when users explicitly request multi-dataset comparison analysis.

## When to Use

Use this module when users:
- Request cross-dataset comparison (e.g., "Compare our results with published data")
- Want to identify reversal patterns across experiments (genes UP by stimulation, DOWN by treatment)
- Need heatmaps comparing in-house DESeq2 results with external datasets
- Ask for biomarker tables with significance annotations across multiple sources

## Default Output Directory

All results are saved to `deg_multi_comparison/` in the working directory:

```
deg_multi_comparison/
├── {prefix}_UP_by_stim_DOWN_by_treatment.tsv
├── {prefix}_DOWN_by_stim_UP_by_treatment.tsv
├── {prefix}_UP_stim_DOWN_treatment_heatmap.png        # Static heatmap
├── {prefix}_UP_stim_DOWN_treatment_heatmap_table.tsv  # Table matching PNG row order
├── {prefix}_UP_stim_DOWN_treatment_heatmap.html       # Interactive HTML heatmap
├── {prefix}_DOWN_stim_UP_treatment_heatmap.png
├── {prefix}_DOWN_stim_UP_treatment_heatmap_table.tsv
├── {prefix}_DOWN_stim_UP_treatment_heatmap.html
├── config.json
└── gsea/                              (pathway mode only)
    └── {prefix}_gsea_all.txt
```

### Heatmap Outputs (Always Generated Together)

| File | Description |
|------|-------------|
| `*_heatmap.png` | Static PNG heatmap (ComplexHeatmap) |
| `*_heatmap_table.tsv` | TSV table with rows in same order as PNG |
| `*_heatmap.html` | Interactive HTML heatmap (heatmaply/plotly) |

All three outputs have **consistent row ordering** (score-sorted by default, or clustered if configured).

## Analysis Modes

| Mode | Description | Use Case |
|------|-------------|----------|
| **Discovery** | Find top N genes matching reversal criteria | Identify candidate biomarkers |
| **Gene List** | Analyze specific genes regardless of significance | Validate known targets |
| **Pathway Analysis** | Compare pathway enrichments (GSEA) across datasets | Understand pathway-level effects |

## Interactive Workflow

### Step 1: Confirm Multi-Dataset Analysis Request

When user requests cross-dataset comparison, confirm the intent:

```
header: "Analysis type"
question: "You want to compare DESeq2 results with external datasets?"
options:
  - label: "Yes - Discovery mode (Recommended)"
    description: "Find top genes with reversal patterns across datasets"
  - label: "Yes - Gene list mode"
    description: "Analyze specific genes across datasets"
  - label: "Yes - Pathway mode"
    description: "Compare pathway enrichments across datasets"
```

### Step 2: Select DESeq2 Comparisons

```
header: "DESeq2 comparisons"
question: "Which comparisons from your DESeq2 analysis should be included?"
options: [List available comparisons from summstats_all.txt]
multiSelect: true
```

### Step 3: Collect External Data Sources

For each external dataset, gather:
- **File path**: Path to data file (Excel, TSV, CSV)
- **Column mappings**: gene_col, log2fc_col, padj_col
- **Label**: Short identifier for this source
- **Group**: Which column group it belongs to
- **Significance cutoffs**: padj threshold, optional log2FC cutoff

### Step 4: Configure Column Groups

Organize columns into logical groups for heatmap visualization:

```
header: "Column groups"
question: "How should columns be grouped?"
options:
  - label: "Auto-group (Recommended)"
    description: "Group by source type (External, InHouse, etc.)"
  - label: "Custom grouping"
    description: "Define custom column groups"
```

### Step 5: Configure Score Logic

Define which direction contributes to the score for each group:

**For UP table** (genes UP by stimulation, DOWN by treatment):
- External disease data: usually "up" (disease = elevated)
- Stimulation data: "up" (stimulation induces)
- Treatment data: "down" (treatment reverses)

**For DOWN table**: opposite directions apply.

### Step 6: Configure Heatmap Options

Configure heatmap visualization settings. Heatmaps are generated using R (ComplexHeatmap package).

```
header: "Heatmap options"
question: "What color scale should be used for the heatmap?"
options:
  - label: "Blue-White-Red (Recommended)"
    description: "Blue for negative, white for zero, red for positive log2FC/NES"
  - label: "Purple-White-Orange"
    description: "Alternative color scheme for colorblind accessibility"
  - label: "Custom"
    description: "Specify custom colors"
```

**Additional heatmap options** (collect as needed):
- **Figure width**: Default 14 inches
- **Figure height**: Default 10 inches (auto-scales with gene count)
- **Font size**: Default 9pt for labels
- **Star font size**: Default 6pt for significance annotations
- **Column gap**: Default 3mm between column groups
- **Cluster rows**: Default FALSE (score-sorted); set TRUE for hierarchical clustering
- **Clustering method**: Default "ward.D2" (when clustering enabled)
- **HTML column font size**: Default 12pt for interactive heatmap labels

**Column label mapping** (optional): Create short display labels for columns:
```json
{
  "column_mapping": {
    "InHouse_TGFb_vs_Ctrl": "TGFb\nvs Ctrl",
    "InHouse_Drug_vs_TGFb": "Drug\nvs TGFb"
  }
}
```

### Step 7: Run Analysis

Execute the adapter script (uses R `generate_heatmaps.R` for visualization):

```bash
python ~/.claude/skills/deseq2-rna/scripts/multidata/multidata_adapter.py \
  --deseq2_output /path/to/deseq2_output \
  --external_config external_sources.json \
  --mode discovery \
  --top_n 50
```

**With custom heatmap settings:**
```bash
python ~/.claude/skills/deseq2-rna/scripts/multidata/multidata_adapter.py \
  --deseq2_output /path/to/deseq2_output \
  --mode discovery \
  --heatmap_colors "blue,white,red" \
  --heatmap_width 16 \
  --heatmap_height 12 \
  --heatmap_fontsize 10
```

## Column Mapping (DESeq2 -> multidata)

| DESeq2 Column | Multidata Parameter | Description |
|---------------|---------------------|-------------|
| `symbol` | `gene_col` | Gene symbol column |
| `log2FoldChange` | `log2fc_col` | Log2 fold change |
| `padj` | `padj_col` | Adjusted p-value |
| `comparison` | `contrast_filter.column` | Comparison name for filtering |

## External Sources Configuration File

Create a JSON file with external dataset configurations:

```json
{
  "sources": [
    {
      "file": "/path/to/external_data.xlsx",
      "sheet": "DEG_results",
      "gene_col": "Gene",
      "log2fc_col": "log2FC",
      "padj_col": "FDR",
      "label": "External_Study1",
      "group": "External",
      "log2fc_cutoff": 0.5
    }
  ],
  "column_groups": {
    "External": ["External_Study1", "External_Study2"],
    "InHouse": ["InHouse_Stim_vs_Ctrl", "InHouse_Drug_vs_Stim"]
  },
  "score_logic": {
    "UP": {
      "External": "up",
      "InHouse": "up"
    },
    "DOWN": {
      "External": "down",
      "InHouse": "down"
    }
  },
  "heatmap": {
    "color_scale": ["blue", "white", "red"],
    "width": 14,
    "height": 10,
    "fontsize": 9,
    "star_fontsize": 6,
    "column_gap": 3,
    "cluster_rows": false,
    "clustering_method": "ward.D2",
    "html_fontsize_col": 12
  },
  "column_mapping": {
    "InHouse_Stim_vs_Ctrl": "Stim\nvs Ctrl",
    "InHouse_Drug_vs_Stim": "Drug\nvs Stim"
  }
}
```

## Heatmap Generation (R-based)

Heatmaps are generated using the R script `generate_heatmaps.R` which produces three synchronized outputs:

### Static PNG Heatmap
- **ComplexHeatmap** package for publication-quality heatmaps
- **circlize** for color mapping
- Row annotation showing gene scores
- Column grouping with labeled separators
- Significance stars overlaid on cells

### Table TSV Export
- Rows in exact same order as PNG heatmap
- Includes all original columns (Gene/Pathway, Score, data columns)
- Useful for downstream analysis or sharing

### Interactive HTML Heatmap
- **heatmaply** package (plotly-based)
- Same row order and BWR color scale as PNG
- Hover tooltips showing exact values
- Column and row side color annotations
- Self-contained HTML (no external dependencies)
- Row labels hidden if >100 rows for readability

### Row Ordering Options

| Setting | Behavior |
|---------|----------|
| `cluster_rows: false` (default) | Rows sorted by Score (descending) |
| `cluster_rows: true` | Hierarchical clustering with dendrogram |

When clustering is enabled, the clustering method can be configured via `clustering_method` (default: "ward.D2").

The R script is invoked automatically by the adapter via `conda run -n r_env`.

## Scripts (Local)

All multidata scripts are included in this skill:

```
~/.claude/skills/deseq2-rna/scripts/multidata/
├── multidata_adapter.py      # Main adapter (orchestrates analysis)
├── generate_tables.py        # Gene-level tables with significance stars
├── generate_heatmaps.R       # ComplexHeatmap visualization (R)
├── generate_pathway_tables.py # Pathway-level tables (NES values)
└── run_gsea.R                # GSEA analysis for pathway mode
```

### Script Functions

**generate_tables.py**:
- `format_with_stars()` - Format log2FC with significance stars
- `is_significant()` - Check if value is significant in a direction
- `load_data_source()` - Load and process a data source
- `calculate_score()` - Calculate gene score based on logic
- `create_table()` - Generate the full table

**generate_heatmaps.R**:
- `extract_lfc()` - Parse numeric log2FC from formatted string
- `extract_sig()` - Extract significance stars
- `create_heatmap()` - Generate PNG, TSV table, and interactive HTML (all with consistent row order)

**generate_pathway_tables.py**:
- `format_nes_with_stars()` - Format NES with significance stars
- `load_gsea_results()` - Load combined GSEA results
- `discover_top_pathways()` - Find top pathways by score
- `create_pathway_table()` - Generate pathway table

## Example Commands

```bash
# Discovery mode (default output: deg_multi_comparison/)
python ~/.claude/skills/deseq2-rna/scripts/multidata/multidata_adapter.py \
  --deseq2_output /path/to/deseq2_output \
  --mode discovery \
  --top_n 50

# Gene list mode
python ~/.claude/skills/deseq2-rna/scripts/multidata/multidata_adapter.py \
  --deseq2_output /path/to/deseq2_output \
  --mode gene_list \
  --genes "GREM1,IL11,NOG,CHRD"

# Pathway analysis mode
python ~/.claude/skills/deseq2-rna/scripts/multidata/multidata_adapter.py \
  --deseq2_output /path/to/deseq2_output \
  --mode pathway_analysis \
  --top_n 100

# With external datasets
python ~/.claude/skills/deseq2-rna/scripts/multidata/multidata_adapter.py \
  --deseq2_output /path/to/deseq2_output \
  --external_config external_sources.json \
  --mode discovery

# Custom output directory
python ~/.claude/skills/deseq2-rna/scripts/multidata/multidata_adapter.py \
  --deseq2_output /path/to/deseq2_output \
  --output_dir /custom/path
```

## Output Format

### TSV Table Structure

| Column | Description |
|--------|-------------|
| Gene/Pathway | Gene symbol or pathway name |
| Score | Total score (sum of direction-matching significant hits) |
| [Source columns] | log2FC/NES values with significance stars |
| n_group_direction | Count of significant values per group in expected direction |

### Significance Annotation

Values are formatted as `{value}{stars}`:
- `****` = padj < 0.0001
- `***` = padj < 0.001
- `**` = padj < 0.01
- `*` = padj < 0.05
- `.` = padj < 0.05 but |log2FC| below cutoff
- (no annotation) = not significant
- `0` = gene not found in data source

Example: `-1.45***` means log2FC = -1.45 with padj < 0.001

### Heatmap Organization

- **Rows**: Genes/pathways sorted by Score (descending)
- **Columns**: Grouped by data source type with labeled separators
- **Cell colors**: Blue-white-red scale for log2FC/NES
- **Cell text**: Significance stars overlaid
- **Row annotation**: Score bar (grey to green)

## Reference Documentation

- `references/multidata_integration.md` - Adapter documentation, column mapping, score logic patterns
- `references/multidata_parameter_guide.md` - Detailed parameter documentation for all scripts
