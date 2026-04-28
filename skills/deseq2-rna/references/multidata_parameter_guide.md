# DEG Multi-Dataset Analysis - Parameter Guide

Comprehensive documentation for all configuration parameters.

## Configuration File Structure

```json
{
  "mode": "...",
  "genes": [...],
  "top_n": N,
  "data_sources": [...],
  "column_groups": {...},
  "score_logic": {...},
  "column_mapping": {...},
  "heatmap": {...},
  "output": {...}
}
```

## Mode Parameters

### mode
- **Type**: String
- **Values**: `"discovery"` | `"gene_list"` | `"pathway_analysis"`
- **Default**: `"gene_list"`
- **Description**: Analysis mode selection

**Discovery Mode**: Finds top N genes that best match the reversal pattern. Useful when exploring data to identify candidate biomarkers.

**Gene List Mode**: Analyzes a specific set of genes regardless of their scores. Useful for hypothesis-driven analysis of known gene sets.

**Pathway Analysis Mode**: Runs GSEA on each data source and identifies pathways with reversal patterns. Useful for understanding biological processes rather than individual genes.

### genes
- **Type**: Array of strings
- **Required for**: Gene List mode
- **Example**: `["GREM1", "GREM2", "NOG", "CHRD"]`
- **Description**: Gene symbols to analyze. Case-sensitive, should match gene symbols in data sources.

### top_n
- **Type**: Integer
- **Default**: `50` (Discovery mode), `100` (Pathway mode)
- **Used in**: Discovery mode, Pathway Analysis mode
- **Description**: Number of top-scoring genes/pathways to return. Items are ranked by score and the top N are selected.

## Gene Set Parameters (Pathway Analysis Mode)

### gene_sets
Configuration for which gene sets to use in GSEA.

```json
{
  "hallmark": true,
  "gobp": true,
  "custom_gmt": "/path/to/custom.gmt"
}
```

#### hallmark
- **Type**: Boolean
- **Default**: `true`
- **Description**: Include MSigDB Hallmark gene sets (50 curated gene sets).

#### gobp
- **Type**: Boolean
- **Default**: `true`
- **Description**: Include GO Biological Process gene sets (~7500 gene sets).

#### custom_gmt
- **Type**: String (path) or null
- **Default**: `null`
- **Description**: Path to custom GMT file with user-defined gene sets.

## GSEA Settings (Pathway Analysis Mode)

### gsea_settings
Configuration for GSEA algorithm parameters.

```json
{
  "seed": 42,
  "pvalueCutoff": 1,
  "minGSSize": 3,
  "maxGSSize": "Inf",
  "nPermSimple": 10000
}
```

#### seed
- **Type**: Integer
- **Default**: `42`
- **Description**: Random seed for reproducibility.

#### pvalueCutoff
- **Type**: Number
- **Default**: `1`
- **Description**: P-value cutoff for GSEA results. Set to 1 to return all results for downstream filtering.

#### minGSSize
- **Type**: Integer
- **Default**: `3`
- **Description**: Minimum gene set size to include in analysis.

#### maxGSSize
- **Type**: Integer or `"Inf"`
- **Default**: `"Inf"`
- **Description**: Maximum gene set size to include in analysis. Use `"Inf"` for no upper limit.

#### nPermSimple
- **Type**: Integer
- **Default**: `10000`
- **Description**: Number of permutations for p-value estimation.

## Data Source Parameters

### data_sources
Array of data source configurations. Each source represents one dataset or comparison.

```json
{
  "file": "/path/to/data.xlsx",
  "sheet": "SheetName",
  "gene_col": "gene_symbol",
  "log2fc_col": "log2FoldChange",
  "padj_col": "padj",
  "label": "Source_Label",
  "group": "Group_Name",
  "log2fc_cutoff": 1.0,
  "padj_cutoff": 0.05,
  "contrast_filter": {
    "column": "Contrast",
    "value": "treatment_vs_control"
  }
}
```

#### file
- **Type**: String (path)
- **Required**: Yes
- **Formats**: Excel (.xlsx, .xls), TSV (.tsv, .txt), CSV (.csv)
- **Description**: Path to the data file containing DEG results.

#### sheet
- **Type**: String
- **Required for**: Excel files
- **Description**: Sheet name containing the data. Ignored for TSV/CSV files.

#### gene_col
- **Type**: String
- **Required**: Yes
- **Description**: Column name containing gene symbols. Common values: `"gene_symbol"`, `"HGNC"`, `"symbol"`, `"gene_name"`

#### log2fc_col
- **Type**: String
- **Required**: Yes
- **Description**: Column name containing log2 fold change values. Common values: `"log2FoldChange"`, `"logFC"`, `"Log<sub>2</sub>(Fold Change)"`

#### padj_col
- **Type**: String
- **Required**: Yes
- **Description**: Column name containing adjusted p-values. Common values: `"padj"`, `"FDR"`, `"Adjusted p-value"`, `"adj.P.Val"`

#### label
- **Type**: String
- **Required**: Yes
- **Description**: Short label for this data source. Used as column header in output tables. Examples: `"TGFb_stim"`, `"GREM1_KD"`, `"Matri_IvH"`

#### group
- **Type**: String
- **Required**: Yes
- **Description**: Logical group this source belongs to. Used for heatmap organization and score calculation. Examples: `"Stimulation"`, `"Treatment"`, `"External"`

#### log2fc_cutoff
- **Type**: Number or null
- **Default**: `null`
- **Description**: Minimum |log2FC| required for significance. If `null`, any magnitude is accepted. Use `0.5` or `1.0` for stricter criteria.

#### padj_cutoff
- **Type**: Number
- **Default**: `0.05`
- **Description**: Maximum adjusted p-value for significance.

#### contrast_filter
- **Type**: Object or null
- **Default**: `null`
- **Description**: Filter rows by a specific column value. Useful when a single file contains multiple contrasts.

**Example**: To select only rows where `Contrast` column contains "grem1":
```json
"contrast_filter": {
  "column": "Contrast",
  "value": "grem1"
}
```

## Column Organization

### column_groups
- **Type**: Object
- **Description**: Groups data source labels into logical categories for heatmap visualization and score calculation.

```json
{
  "External_CD": ["Matri_IvH", "Matri_FvH", "FFPE_IvH", "FFPE_FvH"],
  "External_3D": ["3D_GREM1_KD"],
  "Stimulation": ["TGFb_stim", "TGFb_BMP4_stim"],
  "Antibody": ["TGFb_Gini", "TGFb_BMP4_Gini"]
}
```

Groups appear in the heatmap in the order defined. Columns within each group maintain their specified order.

### column_mapping
- **Type**: Object
- **Description**: Maps column labels to display names for heatmap column headers.

```json
{
  "Matri_IvH": "Inf vs H\n(Matri)",
  "TGFb_stim": "TGFb\nvs NT",
  "TGFb_Gini": "TGFb\nGini"
}
```

Use `\n` for line breaks in display labels.

## Score Logic

### score_logic
- **Type**: Object
- **Description**: Defines which direction contributes to the score for each table type and group.

```json
{
  "UP": {
    "External_CD": "up",
    "External_3D": "down",
    "Stimulation": "up",
    "Antibody": "down"
  },
  "DOWN": {
    "External_CD": "down",
    "External_3D": "up",
    "Stimulation": "down",
    "Antibody": "up"
  }
}
```

**UP Table Logic** (genes upregulated by stimulus, reversed by treatment):
- `"up"`: Significant positive log2FC contributes to score
- `"down"`: Significant negative log2FC contributes to score

**DOWN Table Logic** (genes downregulated by stimulus, reversed by treatment):
- Typically opposite of UP table logic

**Score Calculation**: For each gene, score = count of significant values in the expected direction across all groups.

## Heatmap Parameters

### heatmap
Configuration for heatmap visualization.

```json
{
  "color_scale": ["blue", "white", "red"],
  "row_annotation": "Score",
  "width": 14,
  "height": 10,
  "fontsize": 9,
  "star_fontsize": 6,
  "column_gap": 3
}
```

#### color_scale
- **Type**: Array of 3 colors
- **Default**: `["blue", "white", "red"]`
- **Description**: Colors for negative, zero, positive log2FC values.
- **Alternatives**: `["purple", "white", "orange"]`, `["green", "white", "magenta"]`

#### row_annotation
- **Type**: String
- **Default**: `"Score"`
- **Description**: Column to use for row annotation bar.

#### width
- **Type**: Number (inches)
- **Default**: `14`
- **Description**: Plot width. Adjust based on number of columns.

#### height
- **Type**: Number (inches)
- **Default**: `10`
- **Description**: Minimum plot height. Actual height scales with gene count.

#### fontsize
- **Type**: Number (points)
- **Default**: `9`
- **Description**: Font size for gene names and column labels.

#### star_fontsize
- **Type**: Number (points)
- **Default**: `6`
- **Description**: Font size for significance stars overlaid on cells.

#### column_gap
- **Type**: Number (mm)
- **Default**: `3`
- **Description**: Gap between column groups.

## Output Parameters

### output
```json
{
  "prefix": "BMP_antagonist",
  "directory": "deg_multidata_output"
}
```

#### prefix
- **Type**: String
- **Default**: `"analysis"`
- **Description**: Prefix for all output files.

**Generated files**:
- `{prefix}_UP_by_stim_DOWN_by_treatment.tsv`
- `{prefix}_DOWN_by_stim_UP_by_treatment.tsv`
- `{prefix}_UP_stim_DOWN_treatment_heatmap.png`
- `{prefix}_DOWN_stim_UP_treatment_heatmap.png`

#### directory
- **Type**: String (path)
- **Default**: `deg_multidata_output` (in working directory)
- **Description**: Output directory for all files. Created if it doesn't exist.

## Complete Example Configuration

```json
{
  "mode": "gene_list",
  "genes": ["NOG", "CHRD", "GREM1", "GREM2", "FST", "FSTL1"],

  "data_sources": [
    {
      "file": "/data/engitix_dataset_de.xlsx",
      "sheet": "Engitix.Matrisome.CD.InfvsH",
      "gene_col": "HGNC",
      "log2fc_col": "Log<sub>2</sub>(Fold Change)",
      "padj_col": "Adjusted p-value",
      "label": "Matri_IvH",
      "group": "External_CD",
      "log2fc_cutoff": null
    },
    {
      "file": "/data/GREM1_IAF_summstats_all.txt",
      "gene_col": "symbol",
      "log2fc_col": "log2FoldChange",
      "padj_col": "padj",
      "label": "TGFb_stim",
      "group": "Stimulation",
      "log2fc_cutoff": 1.0,
      "contrast_filter": {
        "column": "comparison",
        "value": "TGFb_RSV_IgG1_vs_Non_treatment"
      }
    }
  ],

  "column_groups": {
    "External CD": ["Matri_IvH"],
    "Stimulation": ["TGFb_stim"],
    "Antibody": ["TGFb_Gini"]
  },

  "score_logic": {
    "UP": {
      "External CD": "up",
      "Stimulation": "up",
      "Antibody": "down"
    },
    "DOWN": {
      "External CD": "down",
      "Stimulation": "down",
      "Antibody": "up"
    }
  },

  "column_mapping": {
    "Matri_IvH": "Inf vs H\n(Matri)",
    "TGFb_stim": "TGFb\nvs NT",
    "TGFb_Gini": "TGFb\nGini"
  },

  "heatmap": {
    "color_scale": ["blue", "white", "red"],
    "width": 12,
    "height": 8
  },

  "output": {
    "prefix": "BMP_antagonist",
    "directory": "deg_multidata_output"
  }
}
```

## Pathway Analysis Mode Example Configuration

```json
{
  "mode": "pathway_analysis",

  "gene_sets": {
    "hallmark": true,
    "gobp": true,
    "custom_gmt": null
  },

  "gsea_settings": {
    "seed": 42,
    "pvalueCutoff": 1,
    "minGSSize": 3,
    "maxGSSize": "Inf"
  },

  "top_n": 100,

  "data_sources": [
    {
      "file": "/data/deg_results.xlsx",
      "sheet": "TGFb_vs_Control",
      "gene_col": "gene_symbol",
      "log2fc_col": "log2FoldChange",
      "padj_col": "padj",
      "label": "TGFb_stim",
      "group": "Stimulation"
    },
    {
      "file": "/data/deg_results.xlsx",
      "sheet": "Antibody_vs_TGFb",
      "gene_col": "gene_symbol",
      "log2fc_col": "log2FoldChange",
      "padj_col": "padj",
      "label": "TGFb_Ab",
      "group": "Treatment"
    }
  ],

  "column_groups": {
    "Stimulation": ["TGFb_stim"],
    "Treatment": ["TGFb_Ab"]
  },

  "score_logic": {
    "UP": {
      "Stimulation": "up",
      "Treatment": "down"
    },
    "DOWN": {
      "Stimulation": "down",
      "Treatment": "up"
    }
  },

  "heatmap": {
    "color_scale": ["blue", "white", "red"],
    "width": 14,
    "height": 12
  },

  "output": {
    "prefix": "pathway_reversal",
    "directory": "deg_multidata_output"
  }
}
```

## Pathway Mode Output Format

### NES Significance Annotation

Values are formatted as `{NES}{stars}`:
- `****` = padj < 0.0001
- `***` = padj < 0.001
- `**` = padj < 0.01
- `*` = padj < 0.05
- (no annotation) = not significant

Example: `1.85***` means NES = 1.85 with padj < 0.001

### Pathway Table Columns

| Column | Description |
|--------|-------------|
| Pathway | Pathway/gene set name (e.g., HALLMARK_EMT) |
| Score | Count of comparisons with NES in expected direction |
| [Source columns] | NES values with significance stars |
| n_group_direction | Count per group in expected direction |

## Troubleshooting

### Common Issues

**"No data columns found"**
- Check that `column_groups` labels match `data_sources` labels exactly
- Verify case sensitivity in label names

**"Gene not found in source"**
- Gene symbols are case-sensitive
- Check source file for exact gene symbol format

**Empty heatmap cells**
- Gene may not be present in that data source
- Displayed as `0` in TSV and white in heatmap

**Stars not showing**
- padj >= 0.05 (not significant)
- Check padj column name and data format

### Validating Configuration

Before running full analysis:
1. Check file paths exist
2. Verify column names in source files
3. Ensure labels are unique across sources
4. Confirm groups contain valid labels
