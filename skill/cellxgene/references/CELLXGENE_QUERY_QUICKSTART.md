# CELLxGENE Query Quick Start

## Installation

```bash
conda activate claude_test
pip install cellxgene-census  # Already installed
```

## Quick Examples

### 1. Simple Tissue Query
```bash
cd /home/sagemaker-user/claude_code/gene-expression-specificity-cellxgene
conda run -n claude_test python scripts/query_cellxgene.py \
    --tissue lung \
    --output lung_adult
```

### 2. Query Specific Cell Type
```bash
conda run -n claude_test python scripts/query_cellxgene.py \
    --tissue intestine \
    --cell-type "T cell" \
    --output intestine_tcells
```

### 3. Query with Gene Set
```bash
conda run -n claude_test python scripts/query_cellxgene.py \
    --tissue liver \
    --genes "APOE,APOC1,APOC2,APOC3" \
    --output liver_apoe_genes
```

### 4. Query Disease Samples
```bash
conda run -n claude_test python scripts/query_cellxgene.py \
    --tissue lung \
    --disease "COVID-19" \
    --output covid_lung
```

### 5. Multi-Tissue Query
```bash
conda run -n claude_test python scripts/query_cellxgene.py \
    --tissue "lung,intestine" \
    --cell-type "epithelial" \
    --output multi_tissue_epithelial
```

### 6. Query by Development Stage
```bash
conda run -n claude_test python scripts/query_cellxgene.py \
    --tissue brain \
    --development-stage embryonic \
    --output embryonic_brain
```

### 7. Query Multiple Species
```bash
conda run -n claude_test python scripts/query_cellxgene.py \
    --species "human,mouse" \
    --tissue brain \
    --cell-type neuron \
    --output brain_neurons_multispecies
```

### 8. Non-Interactive Mode (No Confirmation)
```bash
conda run -n claude_test python scripts/query_cellxgene.py \
    --tissue lung \
    --output lung_data \
    --no-interactive
```

## Common Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--species` | string/list | human | Species name(s) |
| `--tissue` | string/list | None | Tissue name(s) |
| `--cell-type` | string/list | None | Cell type name(s) |
| `--disease` | string/list | None | Disease name(s) |
| `--sex` | string/list | None | Sex |
| `--development-stage` | string/list | adult | Development stage(s) |
| `--genes` | list | None | Gene symbols |
| `--output` | string | cellxgene_query | Output file prefix |
| `--output-dir` | string | . | Output directory |
| `--no-interactive` | flag | False | Skip confirmation |

## Output Files

Each query creates 4 files:
1. `{prefix}_data.h5ad` - AnnData object
2. `{prefix}_metadata.csv` - Cell metadata
3. `{prefix}_summary.json` - Summary statistics
4. `{prefix}_log.json` - Query log

## Python API Usage

```python
from scripts.query_cellxgene import CellxGeneQuery
from pathlib import Path

with CellxGeneQuery() as querier:
    result = querier.query_data(
        tissue="lung",
        cell_type="T cell",
    )

    if result:
        experiment, value_filter, var_filter, obs_df = result

        with experiment.axis_query(
            measurement_name="RNA",
            obs_query=None if value_filter is None else {"value_filter": value_filter},
        ) as query:
            adata = query.to_anndata(X_name="raw")
            obs_df_filtered = query.obs().concat().to_pandas()

        summary = querier.generate_summary(adata, obs_df_filtered, value_filter or "None")
        querier.print_summary(summary)

        querier.save_results(
            adata, obs_df_filtered, summary,
            output_dir=Path("./output"),
            prefix="my_query"
        )
```

## Tips

1. **Start Small**: Use restrictive filters first (tissue + cell type)
2. **Check Summary**: Review before downloading
3. **Use Gene Sets**: Reduce data size with `--genes`
4. **Save Logs**: Track analyses with log files
5. **Case-Insensitive**: All matching is case-insensitive
6. **Partial Match**: String inputs use partial matching
7. **Any Match**: List inputs match any element

## Documentation

Full documentation: `scripts/QUERY_README.md`

## Examples

Run pre-built examples:
```bash
# Run example 1 (lung tissue)
conda run -n claude_test python scripts/example_query.py --example 1

# Run example 3 (gene set query)
conda run -n claude_test python scripts/example_query.py --example 3
```

Available examples (1-5):
1. Simple tissue query (lung)
2. Cell type across multiple tissues
3. Gene set query (APOE family)
4. Disease query (COVID-19)
5. Multi-species query

## Workflow Example

```bash
# 1. Query data
conda run -n claude_test python scripts/query_cellxgene.py \
    --tissue lung \
    --cell-type "epithelial cell" \
    --output lung_epithelial

# 2. Review summary (printed to console)
# Total Cells: 50,000
# Total Genes: 36,000
# ...

# 3. Confirm download (y/n)
# yes

# 4. Files saved:
# - lung_epithelial_data.h5ad
# - lung_epithelial_metadata.csv
# - lung_epithelial_summary.json
# - lung_epithelial_log.json

# 5. Load and analyze
python
>>> import scanpy as sc
>>> adata = sc.read_h5ad("lung_epithelial_data.h5ad")
>>> adata
AnnData object with n_obs × n_vars = 50000 × 36000
```

## Common Use Cases

### Marker Gene Discovery
```bash
conda run -n claude_test python scripts/query_cellxgene.py \
    --tissue intestine \
    --cell-type "epithelial" \
    --genes "CDH1,EPCAM,KRT8,KRT18,VIL1" \
    --output intestine_epithelial_markers
```

### Disease vs Normal Comparison
```bash
# Query normal samples
conda run -n claude_test python scripts/query_cellxgene.py \
    --tissue lung \
    --disease "normal" \
    --output lung_normal

# Query disease samples
conda run -n claude_test python scripts/query_cellxgene.py \
    --tissue lung \
    --disease "COVID-19" \
    --output lung_covid
```

### Cell Type Enrichment
```bash
conda run -n claude_test python scripts/query_cellxgene.py \
    --tissue "lung,intestine,liver" \
    --cell-type "macrophage" \
    --output macrophage_multi_tissue
```

## Troubleshooting

### No Results
- Check filter spelling
- Try broader filters
- Remove some filters

### Too Slow
- Add more restrictive filters
- Use `--genes` for specific genes
- Query smaller tissues

### Out of Memory
- More restrictive filters
- Query specific genes only
- Process in batches

## Contact

For issues or questions, see the full documentation in `scripts/QUERY_README.md`
