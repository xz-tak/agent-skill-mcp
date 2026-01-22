# CELLxGENE Census Query Script

A flexible Python script to query single-cell data from CELLxGENE Census with comprehensive filtering options.

## Features

- **Flexible filtering** by:
  - Species (default: human)
  - Tissue (partial or exact match)
  - Cell type (partial or exact match)
  - Disease (partial or exact match)
  - Sex
  - Development stage (default: adult, with intelligent interpretation)
  - Drug treatment
  - Gene sets

- **Case-insensitive matching** for all text filters
- **Partial matching** for string inputs, **any matching** for list/set inputs
- **Automatic summarization** of query results
- **Interactive download confirmation**
- **Comprehensive logging** of filters and results
- **Saves**:
  - AnnData object (.h5ad)
  - Metadata CSV
  - Summary JSON
  - Query log JSON

## Installation

Ensure you have the required dependencies:

```bash
conda activate claude_test
pip install cellxgene-census anndata pandas numpy
```

## Command Line Usage

### Basic Examples

```bash
# Query human lung tissue (adult by default)
python scripts/query_cellxgene.py --tissue lung --output lung_adult

# Query specific cell type
python scripts/query_cellxgene.py --tissue lung --cell-type "epithelial cell" --output lung_epithelial

# Query with multiple tissues
python scripts/query_cellxgene.py --tissue "lung,intestine" --cell-type "T cell" --output tcells_multi_tissue

# Query with disease filter
python scripts/query_cellxgene.py --tissue lung --disease "COVID-19" --output covid_lung

# Query with specific genes
python scripts/query_cellxgene.py --tissue liver --genes "APOE,APOC1,APOC2,APOC3" --output liver_apoe

# Query multiple species
python scripts/query_cellxgene.py --species "human,mouse" --tissue brain --cell-type neuron --output brain_neurons

# Query by sex
python scripts/query_cellxgene.py --tissue "prostate" --sex male --output male_prostate

# Query by development stage
python scripts/query_cellxgene.py --tissue brain --development-stage embryonic --output embryonic_brain

# Non-interactive mode (no confirmation prompt)
python scripts/query_cellxgene.py --tissue lung --output lung_data --no-interactive
```

### Command Line Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--species` | Species name(s), comma-separated | human |
| `--tissue` | Tissue name(s), comma-separated | None |
| `--cell-type` | Cell type name(s), comma-separated | None |
| `--disease` | Disease name(s), comma-separated | None |
| `--sex` | Sex, comma-separated | None |
| `--development-stage` | Development stage(s), comma-separated | adult |
| `--drug-treatment` | Drug treatment(s), comma-separated | None |
| `--genes` | Gene symbols, comma-separated | None (all genes) |
| `--output` | Output file prefix | cellxgene_query |
| `--output-dir` | Output directory | . (current dir) |
| `--no-interactive` | Skip confirmation prompt | False |

## Python API Usage

You can also use the script as a Python library:

```python
from query_cellxgene import CellxGeneQuery
from pathlib import Path

# Create querier (use context manager for automatic cleanup)
with CellxGeneQuery() as querier:
    # Query data
    result = querier.query_data(
        tissue="lung",
        cell_type="T cell",
        development_stage="adult"
    )

    if result:
        experiment, value_filter, var_filter, obs_df = result

        # Execute query
        with experiment.axis_query(
            measurement_name="RNA",
            obs_query=None if value_filter is None else {"value_filter": value_filter},
            var_query=None if var_filter is None else {"value_filter": var_filter},
        ) as query:
            # Get AnnData
            adata = query.to_anndata(X_name="raw")
            obs_df_filtered = query.obs().concat().to_pandas()

        # Generate summary
        summary = querier.generate_summary(adata, obs_df_filtered, value_filter or "None")
        querier.print_summary(summary)

        # Save results
        querier.save_results(
            adata,
            obs_df_filtered,
            summary,
            output_dir=Path("./output"),
            prefix="my_query"
        )
```

### More Python Examples

```python
# Example 1: Query with gene set
with CellxGeneQuery() as querier:
    result = querier.query_data(
        tissue="liver",
        gene_set=["APOE", "APOC1", "APOC2"]
    )

# Example 2: Query multiple conditions
with CellxGeneQuery() as querier:
    result = querier.query_data(
        species=["human", "mouse"],
        tissue=["lung", "intestine"],
        cell_type="epithelial",
        disease="normal",
        development_stage="adult"
    )

# Example 3: Query by sex and disease
with CellxGeneQuery() as querier:
    result = querier.query_data(
        tissue="lung",
        disease="COVID-19",
        sex="male",
        development_stage="adult"
    )
```

## Filter Matching Logic

### String Input (Single Value)
- **Partial match** for tissue, cell_type, disease
- Example: `tissue="lung"` matches "lung", "left lung", "right lung"

### List/Set Input (Multiple Values)
- **Any match** (checks if any element matches)
- Example: `tissue=["lung", "intestine"]` matches any tissue containing "lung" OR "intestine"

### Case Insensitivity
All text matching is case-insensitive.

### Development Stage Interpretation
The script intelligently interprets development stages:
- `"adult"` → matches "adult", "mature"
- `"embryonic"` → matches "embryo", "embryonic"
- `"fetal"` → matches "fetus", "fetal"
- `"postnatal"` → matches "postnatal", "newborn", "infant", "child"

## Output Files

For each query, the following files are generated:

1. **`{prefix}_data.h5ad`**: AnnData object with expression data and metadata
2. **`{prefix}_metadata.csv`**: Observation metadata in CSV format
3. **`{prefix}_summary.json`**: Summary statistics in JSON format
4. **`{prefix}_log.json`**: Query log with filters and timestamps

## Summary Statistics

The summary includes:
- Total number of cells
- Total number of genes
- Number of unique donors
- Breakdown by:
  - Organism
  - Tissue
  - Cell type
  - Disease
  - Sex
  - Development stage
  - Assay

## Example Workflow

```bash
# 1. Query lung T cells
python scripts/query_cellxgene.py \
    --tissue lung \
    --cell-type "T cell" \
    --development-stage adult \
    --output lung_tcells

# Output:
# ==================================================
# QUERY RESULTS SUMMARY
# ==================================================
# Timestamp: 2025-12-05T10:30:45.123456
# Filters Applied: (tissue_general == "lung") and (cell_type == "T cell") ...
#
# Total Cells: 125,432
# Total Genes: 36,601
# Unique Donors: 45
#
# Organisms:
#   Homo sapiens: 125,432
#
# Tissues:
#   lung: 125,432
#
# Cell Types:
#   T cell: 98,234
#   CD4+ T cell: 15,678
#   CD8+ T cell: 11,520
# ...
#
# Do you want to download and save this data? (yes/no): yes
#
# Saving AnnData object to: ./lung_tcells_data.h5ad
# Saving metadata to: ./lung_tcells_metadata.csv
# Saving summary to: ./lung_tcells_summary.json
# Saving query log to: ./lung_tcells_log.json
```

## Running Examples

See pre-built examples in `example_query.py`:

```bash
# Run specific example
python scripts/example_query.py --example 1

# Run all examples (warning: downloads a lot of data!)
python scripts/example_query.py
```

Available examples:
1. Simple tissue query (lung)
2. Cell type across multiple tissues (T cells in lung and intestine)
3. Gene set query (APOE family in liver)
4. Disease query (COVID-19 lung samples)
5. Multi-species query (brain neurons from human and mouse)

## Tips

1. **Start small**: Begin with restrictive filters (specific tissue + cell type) to get manageable data sizes
2. **Check summary first**: Review the summary before downloading to ensure it matches expectations
3. **Use gene sets**: If analyzing specific genes, use `--genes` to reduce data size
4. **Save logs**: Query logs help reproduce results and track analyses
5. **Species default**: If not specified, queries default to human
6. **Development stage default**: If not specified, queries default to adult

## Troubleshooting

### Issue: Query returns no results
- Check filter spelling and capitalization (though matching is case-insensitive)
- Try broader filters (e.g., partial tissue name instead of exact)
- Remove some filters to see what data is available

### Issue: Query is too slow
- Add more restrictive filters (tissue, cell type, development stage)
- Use `--genes` to limit to specific genes
- Consider querying smaller tissue regions

### Issue: Out of memory
- Use more restrictive filters to reduce data size
- Query specific genes instead of all genes
- Process data in batches (query different tissues separately)

## Notes

- **Drug treatment filtering**: Currently not directly supported in CELLxGENE Census schema (warning will be shown)
- **Data freshness**: Queries use the latest CELLxGENE Census release
- **Performance**: First query may be slow due to Census initialization; subsequent queries are faster

## Support

For issues with:
- **CELLxGENE Census**: https://chanzuckerberg.github.io/cellxgene-census/
- **This script**: Check logs in `{prefix}_log.json` and summary in `{prefix}_summary.json`
