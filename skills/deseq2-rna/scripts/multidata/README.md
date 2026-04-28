# Multi-Dataset DEG Analysis Adapter

Bridge script that adapts DESeq2 output for use with `deg-multidata-analysis` tools.

## Quick Start

```bash
# Discovery mode - find top genes with reversal patterns
python multidata_adapter.py \
  --deseq2_output /path/to/deseq2_output \
  --mode discovery \
  --top_n 50

# Gene list mode - analyze specific genes
python multidata_adapter.py \
  --deseq2_output /path/to/deseq2_output \
  --mode gene_list \
  --genes "GREM1,IL11,NOG,CHRD"

# Pathway mode - compare pathway enrichments
python multidata_adapter.py \
  --deseq2_output /path/to/deseq2_output \
  --mode pathway_analysis \
  --top_n 100
```

## Output

Results are saved to `deg_multi_comparison/` (default) in the working directory:

```
deg_multi_comparison/
├── {prefix}_UP_by_stim_DOWN_by_treatment.tsv
├── {prefix}_DOWN_by_stim_UP_by_treatment.tsv
├── {prefix}_UP_stim_DOWN_treatment_heatmap.png
├── {prefix}_DOWN_stim_UP_treatment_heatmap.png
└── config.json
```

## Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--deseq2_output` | Path to DESeq2 output directory | Required |
| `--external_config` | Path to external sources JSON config | None |
| `--mode` | Analysis mode: discovery, gene_list, pathway_analysis | discovery |
| `--genes` | Comma-separated gene list (gene_list mode) | None |
| `--top_n` | Number of top genes/pathways | 50 |
| `--comparisons` | Comma-separated comparison names to include | All |
| `--output_dir` | Output directory | deg_multi_comparison |
| `--prefix` | Output file prefix | multidata |
| `--conda_env` | Conda environment for R scripts | r_env |

## External Data Configuration

Create a JSON file with external dataset configurations:

```json
{
  "sources": [
    {
      "file": "/path/to/data.xlsx",
      "sheet": "DEG_results",
      "gene_col": "Gene",
      "log2fc_col": "log2FC",
      "padj_col": "FDR",
      "label": "External_Study",
      "group": "External",
      "log2fc_cutoff": 0.5
    }
  ],
  "column_groups": {
    "External": ["External_Study"],
    "InHouse": ["InHouse_Comparison1"]
  },
  "score_logic": {
    "UP": {"External": "up", "InHouse": "up"},
    "DOWN": {"External": "down", "InHouse": "down"}
  }
}
```

Then run:

```bash
python multidata_adapter.py \
  --deseq2_output /path/to/deseq2_output \
  --external_config external_sources.json
```

## Dependencies

This script invokes tools from `deg-multidata-analysis`:
- `generate_tables.py` - Gene-level TSV tables
- `generate_pathway_tables.py` - Pathway-level tables
- `generate_heatmaps.R` - Heatmap visualization

## Reference

See `references/multidata_integration.md` for detailed documentation.
