# Multi-Dataset Integration Reference

This document provides detailed documentation for the Multi-Dataset Comparison module, which enables comparison of DESeq2 results with external datasets.

## Overview

The multi-dataset module bridges DESeq2 output with the `deg-multidata-analysis` skill to generate:
- TSV tables with log2FC/NES values and significance annotations
- Heatmaps organized by data source groups
- Score rankings based on reversal patterns

## Architecture

```
DESeq2 Output                    External Data
     │                               │
     ▼                               ▼
{prefix}_summstats_all.txt     user-provided files
     │                               │
     └──────────────┬────────────────┘
                    │
                    ▼
          multidata_adapter.py
                    │
     ┌──────────────┼──────────────┐
     │              │              │
     ▼              ▼              ▼
generate_    generate_      run_gsea.R
tables.py   heatmaps.R    (pathway mode)
     │              │              │
     └──────────────┼──────────────┘
                    ▼
          deg_multi_comparison/
```

## Column Mapping Reference

### DESeq2 to Multidata

| DESeq2 Output Column | Multidata Config Key | Description |
|---------------------|---------------------|-------------|
| `symbol` | `gene_col` | Gene symbol |
| `log2FoldChange` | `log2fc_col` | Log2 fold change value |
| `padj` | `padj_col` | Adjusted p-value (BH) |
| `comparison` | `contrast_filter.column` | Comparison name for row filtering |
| `baseMean` | - | Not used by multidata |
| `pvalue` | - | Not used by multidata |

### External Data Configuration

For each external data source, specify:

```json
{
  "file": "/path/to/data.xlsx",     // Required: path to data file
  "sheet": "Sheet1",                 // Excel only: sheet name
  "gene_col": "Gene",                // Required: gene symbol column
  "log2fc_col": "log2FC",            // Required: log2 fold change column
  "padj_col": "FDR",                 // Required: adjusted p-value column
  "label": "Study_Name",             // Required: unique identifier
  "group": "External",               // Required: column group for heatmap
  "log2fc_cutoff": 0.5,              // Optional: |log2FC| threshold for stars
  "padj_cutoff": 0.05,               // Optional: padj threshold (default 0.05)
  "contrast_filter": {               // Optional: filter rows by value
    "column": "Contrast",
    "value": "Disease_vs_Normal"
  }
}
```

## Score Logic Configuration

Score logic defines which direction (up/down) contributes to the score for each column group.

### Reversal Pattern Logic

For identifying genes that are **UP in disease/stimulation** and **DOWN with treatment**:

```json
{
  "score_logic": {
    "UP": {
      "External": "up",      // Disease data: upregulated = good
      "Stimulation": "up",   // Stim effect: upregulated = good
      "Treatment": "down"    // Treatment: downregulated = reversal
    },
    "DOWN": {
      "External": "down",    // Disease data: downregulated = good
      "Stimulation": "down", // Stim effect: downregulated = good
      "Treatment": "up"      // Treatment: upregulated = reversal
    }
  }
}
```

### Common Patterns

#### Pattern 1: Disease Reversal
Finding genes UP in disease, DOWN with drug:
```json
{
  "UP": {"Disease": "up", "Drug": "down"},
  "DOWN": {"Disease": "down", "Drug": "up"}
}
```

#### Pattern 2: Knockdown Validation
Comparing knockdown with overexpression:
```json
{
  "UP": {"KD": "down", "OE": "up"},
  "DOWN": {"KD": "up", "OE": "down"}
}
```

#### Pattern 3: Multi-source Consensus
Finding genes consistent across studies:
```json
{
  "UP": {"Study1": "up", "Study2": "up", "Study3": "up"},
  "DOWN": {"Study1": "down", "Study2": "down", "Study3": "down"}
}
```

## Column Groups Configuration

Column groups organize data columns in the heatmap with visual separators.

```json
{
  "column_groups": {
    "External": ["PublishedStudy1", "PublishedStudy2"],
    "Stimulation": ["TGFb_vs_Ctrl", "IL6_vs_Ctrl"],
    "Treatment": ["DrugA_vs_Stim", "DrugB_vs_Stim"]
  }
}
```

**Key rules:**
- Group names appear as column headers in heatmap
- Order of groups determines left-to-right arrangement
- Each column label must be unique across all groups

## Significance Annotation

Values in output tables are formatted with significance stars:

| Stars | Meaning |
|-------|---------|
| `****` | padj < 0.0001 |
| `***` | padj < 0.001 |
| `**` | padj < 0.01 |
| `*` | padj < 0.05 |
| `.` | padj < 0.05 but \|log2FC\| below cutoff |
| (none) | Not significant |
| `0` | Gene not found in data source |

Example: `-1.45***` means log2FC = -1.45 with padj < 0.001

## Mode-Specific Behavior

### Discovery Mode

1. Loads all genes from all data sources
2. Calculates score for each gene based on score_logic
3. Ranks genes by total score
4. Returns top N genes
5. Generates UP and DOWN tables

**Score calculation:**
- For each column group, check if gene is significant in expected direction
- Score = count of columns matching expected pattern

### Gene List Mode

1. Uses user-provided gene list
2. Extracts values for those genes from all sources
3. Calculates scores (same as discovery)
4. Generates tables with all specified genes

### Pathway Analysis Mode

1. Reads GSEA results (gsea_all.txt) from DESeq2 output
2. Uses NES values instead of log2FC
3. Applies same score logic to pathway enrichments
4. Generates pathway-level tables and heatmaps

## Example Configurations

### Example 1: Comparing with Published Fibrosis Data

```json
{
  "sources": [
    {
      "file": "/data/published/fibrosis_deg.xlsx",
      "sheet": "IPF_vs_Normal",
      "gene_col": "Gene.symbol",
      "log2fc_col": "logFC",
      "padj_col": "adj.P.Val",
      "label": "IPF_Published",
      "group": "External",
      "log2fc_cutoff": 0.5
    },
    {
      "file": "/data/published/nash_deg.csv",
      "gene_col": "Symbol",
      "log2fc_col": "log2FoldChange",
      "padj_col": "padj",
      "label": "NASH_Published",
      "group": "External"
    }
  ],
  "column_groups": {
    "External": ["IPF_Published", "NASH_Published"],
    "InHouse": ["InHouse_TGFb_vs_Ctrl", "InHouse_Drug_vs_TGFb"]
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
  }
}
```

### Example 2: BMP Antagonist Analysis

```json
{
  "sources": [
    {
      "file": "/data/external/matrisome_deg.xlsx",
      "sheet": "CD_Inflamed",
      "gene_col": "Gene",
      "log2fc_col": "log2FC",
      "padj_col": "FDR",
      "label": "Matrisome_CD",
      "group": "External"
    }
  ],
  "column_groups": {
    "External": ["Matrisome_CD"],
    "KD": ["InHouse_GREM1_KD"],
    "Stimulation": ["InHouse_TGFb_vs_Ctrl"],
    "Antibody": ["InHouse_AntiGREM1_vs_TGFb"]
  },
  "score_logic": {
    "UP": {
      "External": "up",
      "KD": "down",
      "Stimulation": "up",
      "Antibody": "down"
    },
    "DOWN": {
      "External": "down",
      "KD": "up",
      "Stimulation": "down",
      "Antibody": "up"
    }
  }
}
```

## Troubleshooting

### Common Issues

#### "Column not found" error
- Check that `gene_col`, `log2fc_col`, `padj_col` match exact column names
- Column names are case-sensitive
- Use quotes around column names with spaces

#### "No genes found" in output
- Check that gene symbols match between sources (HGNC symbols recommended)
- Verify contrast_filter is correctly filtering rows
- Ensure data file has data (not just headers)

#### Empty heatmap cells
- Gene not found in that data source (shows as "0")
- This is normal when comparing across species or platforms

#### Low scores
- Check score_logic directions match your hypothesis
- Verify significance thresholds aren't too strict
- Ensure column_groups include all relevant columns

### Debugging Commands

```bash
# Check DESeq2 output structure
head -5 deg/*_summstats_all.txt

# List available comparisons
cut -f9 deg/*_summstats_all.txt | sort -u

# Validate external file columns
head -1 /path/to/external_data.xlsx

# Test adapter with verbose output
python multidata_adapter.py --deseq2_output ./output --mode discovery 2>&1 | tee debug.log
```

## Best Practices

1. **Consistent gene symbols**: Use HGNC symbols across all sources
2. **Meaningful labels**: Use descriptive labels that fit in heatmap columns
3. **Logical grouping**: Group related data sources together
4. **Clear score logic**: Document the biological rationale for each direction
5. **Verify before running**: Check a few genes manually in source files
6. **Start with discovery**: Find top genes first, then switch to gene_list for targeted analysis
