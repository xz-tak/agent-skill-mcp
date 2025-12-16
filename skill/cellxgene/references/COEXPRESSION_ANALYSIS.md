# Gene Coexpression Analysis Guide

Comprehensive guide for analyzing pairwise gene coexpression from CELLxGENE Census data.

## Overview

The `analyze_coexpression.py` script computes pairwise gene correlations from AnnData objects and generates publication-quality visualizations with hierarchical clustering.

### Key Features
- **Flexible Input**: Accepts AnnData files from CELLxGENE queries or local files
- **Gene Lists**: Supports file-based or comma-separated gene lists
- **Correlation Methods**: Pearson (linear) or Spearman (rank-based)
- **Metadata Filtering**: Filter by cell type, tissue, disease, or custom expressions
- **Hierarchical Clustering**: Optional clustering for pattern discovery
- **Comprehensive Output**: Correlation matrix CSV, heatmap PNG, and summary report

## Quick Start

```bash
# Basic usage with comma-separated genes
python scripts/analyze_coexpression.py \
  --input data.h5ad \
  --genes "GENE1,GENE2,GENE3,GENE4" \
  --output my_coexp

# With gene list file
python scripts/analyze_coexpression.py \
  --input data.h5ad \
  --genes my_genes.txt \
  --output my_coexp

# With filtering and clustering
python scripts/analyze_coexpression.py \
  --input lung_data.h5ad \
  --genes immune_genes.txt \
  --cell-type "T cell" \
  --method spearman \
  --cluster \
  --output tcell_coexp
```

## Input Requirements

### AnnData File
- Must be in `.h5ad` format
- Can be from CELLxGENE query or any standard AnnData object
- Should contain normalized or raw expression data
- Metadata fields optional but enable filtering

### Gene List
**Format 1: File** (one gene per line):
```
GENE1
GENE2
GENE3
GENE4
```

**Format 2: Command Line** (comma-separated):
```bash
--genes "GENE1,GENE2,GENE3,GENE4"
```

**Minimum**: 2 genes required for correlation analysis

## Command Line Arguments

### Required
- `--input`, `-i`: Input AnnData file (.h5ad)
- `--genes`, `-g`: Gene list file or comma-separated names

### Output
- `--output`, `-o`: Output prefix (default: coexpression)
- `--output-dir`: Output directory (default: current directory)

### Analysis Parameters
- `--method`: Correlation method
  - `pearson` (default): Linear correlation, good for normally distributed data
  - `spearman`: Rank-based correlation, robust to outliers
- `--use-raw`: Use `.raw.X` instead of `.X` layer
- `--min-cells`: Minimum cells expressing each gene (default: 10)

### Filtering Parameters
- `--cell-type`: Filter by cell type (partial, case-insensitive match)
- `--tissue`: Filter by tissue (partial, case-insensitive match)
- `--disease`: Filter by disease (partial, case-insensitive match)
- `--metadata-filter`: Custom filter expression (e.g., `"sex == 'female' and age > 30"`)

### Visualization Parameters
- `--figsize`: Figure size in inches (width height) (default: 10 8)
- `--cmap`: Colormap name (default: RdBu_r)
- `--cluster`: Enable hierarchical clustering (groups similar genes)
- `--no-annot`: Disable correlation value annotations on heatmap
- `--title`: Custom plot title

## Output Files

Each analysis generates 3 files:

### 1. Correlation Matrix CSV
**Filename**: `{prefix}_correlation_matrix.csv`

Symmetric matrix with pairwise correlations:
```csv
,GENE1,GENE2,GENE3,GENE4
GENE1,1.000,0.856,0.234,-0.123
GENE2,0.856,1.000,0.567,0.012
GENE3,0.234,0.567,1.000,0.789
GENE4,-0.123,0.012,0.789,1.000
```

### 2. Heatmap PNG
**Filename**: `{prefix}_heatmap.png`

High-resolution (300 DPI) heatmap with:
- Color indicating correlation strength (-1 to 1)
- Optional hierarchical clustering (dendrograms)
- Optional correlation values in cells
- Colorbar with correlation scale

### 3. Summary Report TXT
**Filename**: `{prefix}_summary.txt`

Includes:
- Analysis parameters (method, genes, filters)
- Correlation statistics (mean, median, min, max)
- List of genes analyzed
- List of genes not found
- Top 20 strongest correlations (ranked by absolute value)

## Usage Examples

### Example 1: T Cell Markers

```bash
# Query T cells from lung
python scripts/query_cellxgene.py \
  --tissue lung \
  --cell-type "T cell" \
  --output lung_tcells

# Analyze T cell marker coexpression
python scripts/analyze_coexpression.py \
  --input lung_tcells_data.h5ad \
  --genes "CD3D,CD3E,CD4,CD8A,CD8B,GZMB,PRF1,CCR7,SELL" \
  --method pearson \
  --cluster \
  --title "T Cell Marker Coexpression" \
  --output tcell_markers_coexp
```

**Interpretation**:
- CD3D, CD3E high correlation → co-regulated TCR components
- CD4 vs CD8A/B negative/low correlation → distinct T cell subsets
- GZMB, PRF1 high correlation → cytotoxic program
- CCR7, SELL high correlation → naive/central memory markers

### Example 2: APOE Family in Liver

```bash
# Query liver data
python scripts/query_cellxgene.py \
  --tissue liver \
  --disease normal \
  --output liver_normal

# Analyze APOE family coexpression
python scripts/analyze_coexpression.py \
  --input liver_normal_data.h5ad \
  --genes "APOE,APOC1,APOC2,APOC3,APOC4,APOA1,APOA2" \
  --method pearson \
  --cluster \
  --title "APOE Family Coexpression in Liver" \
  --output apoe_liver_coexp

# Alternative: Only in hepatocytes
python scripts/analyze_coexpression.py \
  --input liver_normal_data.h5ad \
  --genes "APOE,APOC1,APOC2,APOC3,APOC4,APOA1,APOA2" \
  --cell-type hepatocyte \
  --cluster \
  --output apoe_hepatocyte_coexp
```

### Example 3: Interferon Response Genes

```bash
# Query immune cells
python scripts/query_cellxgene.py \
  --tissue "lung,blood" \
  --cell-type immune \
  --output immune_cells

# Analyze IFN response gene coexpression
python scripts/analyze_coexpression.py \
  --input immune_cells_data.h5ad \
  --genes "IFI27,IFI44L,IFI6,IFIT1,IFIT2,IFIT3,ISG15,ISG20,MX1,MX2,OAS1,OAS2,OAS3" \
  --method spearman \
  --cluster \
  --title "Interferon Response Gene Coexpression" \
  --output ifn_response_coexp
```

### Example 4: Disease-Specific Analysis

```bash
# Compare COVID-19 vs normal
python scripts/query_cellxgene.py \
  --tissue lung \
  --output lung_all

# Normal samples
python scripts/analyze_coexpression.py \
  --input lung_all_data.h5ad \
  --genes inflammatory_genes.txt \
  --disease "normal" \
  --cluster \
  --output normal_inflam_coexp

# COVID-19 samples
python scripts/analyze_coexpression.py \
  --input lung_all_data.h5ad \
  --genes inflammatory_genes.txt \
  --disease "COVID-19" \
  --cluster \
  --output covid_inflam_coexp

# Compare correlation matrices to identify rewiring
```

### Example 5: Sex-Specific Differences

```bash
# Query data
python scripts/query_cellxgene.py \
  --tissue heart \
  --output heart_data

# Male coexpression
python scripts/analyze_coexpression.py \
  --input heart_data.h5ad \
  --genes cardiac_genes.txt \
  --metadata-filter "sex == 'male'" \
  --output male_cardiac_coexp

# Female coexpression
python scripts/analyze_coexpression.py \
  --input heart_data.h5ad \
  --genes cardiac_genes.txt \
  --metadata-filter "sex == 'female'" \
  --output female_cardiac_coexp
```

## Correlation Methods

### Pearson Correlation
**Use when**:
- Linear relationships expected
- Data is approximately normally distributed
- Interested in linear co-regulation

**Formula**: Measures linear correlation between two variables

**Range**: -1 (perfect negative) to +1 (perfect positive)

**Interpretation**:
- > 0.7: Strong positive correlation
- 0.3-0.7: Moderate positive correlation
- -0.3-0.3: Weak/no correlation
- < -0.7: Strong negative correlation

### Spearman Correlation
**Use when**:
- Non-linear monotonic relationships
- Data has outliers or is non-normal
- Interested in rank-order relationships

**Formula**: Pearson correlation on rank-transformed data

**Advantages**:
- Robust to outliers
- Detects monotonic relationships
- No distribution assumptions

## Hierarchical Clustering

Enable with `--cluster` flag.

**Method**: Average linkage clustering

**Purpose**:
- Groups genes with similar correlation patterns
- Reveals modules of co-regulated genes
- Makes heatmap easier to interpret

**Interpretation**:
- Genes close together in clustered heatmap have similar correlation profiles
- Distinct blocks indicate gene modules
- Dendrogram height indicates dissimilarity

## Filtering Strategies

### Filter by Cell Type
```bash
--cell-type "T cell"  # Partial match
```
Analyzes coexpression specifically in T cells, removing other cell types.

**Use when**:
- Cell type-specific expression programs expected
- Mixed cell type data but interest in one type
- Comparing coexpression across cell types

### Filter by Tissue
```bash
--tissue lung  # Partial match
```
Restricts analysis to specific tissue context.

**Use when**:
- Multi-tissue dataset
- Tissue-specific regulation expected
- Comparing tissue-specific patterns

### Filter by Disease
```bash
--disease "COVID-19"  # Partial match
```
Analyzes coexpression in disease vs normal.

**Use when**:
- Studying disease-induced rewiring
- Comparing healthy vs disease states
- Identifying disease biomarkers

### Custom Metadata Filter
```bash
--metadata-filter "sex == 'female' and development_stage == 'adult'"
```
Complex filtering using pandas query syntax.

**Use when**:
- Multiple criteria needed
- Specific cohort definition required
- Advanced filtering beyond standard options

## Interpretation Guidelines

### High Positive Correlation (> 0.7)
- Genes likely co-regulated
- May share transcription factors
- Often in same pathway/complex
- Could be functionally related

### Moderate Correlation (0.3-0.7)
- Some shared regulation
- May be indirectly related
- Context-dependent co-expression

### Low/No Correlation (-0.3 to 0.3)
- Independent regulation
- Different pathways/processes
- Not co-expressed in this context

### Negative Correlation (< -0.3)
- Mutually exclusive expression
- Antagonistic pathways
- Different cell states/subtypes
- Temporal dynamics

### Clustered Modules
- Genes clustering together form co-expression modules
- Modules often represent:
  - Shared biological pathways
  - Common regulatory programs
  - Coordinated cellular responses
  - Functional complexes

## Troubleshooting

### Error: "Need at least 2 genes in dataset"
**Cause**: Genes not found in AnnData var_names

**Solution**:
- Check gene symbol spelling
- Verify genes are in dataset (check adata.var_names)
- Use correct gene symbol format (official HUGO symbols)

### Warning: "Genes expressed in fewer than X cells"
**Cause**: Low gene expression

**Solution**:
- Lower `--min-cells` threshold
- Check if genes are truly expressed in this context
- Consider removing low-expression genes from analysis

### Issue: Correlation matrix all zeros/NaN
**Cause**: No variance in gene expression

**Solution**:
- Check if data is normalized
- Verify genes have non-zero expression
- Use `--use-raw` if normalized data in .X

### Issue: Heatmap too crowded
**Cause**: Too many genes

**Solution**:
- Reduce number of genes analyzed
- Increase `--figsize` (e.g., `--figsize 15 12`)
- Use `--no-annot` to remove value labels
- Split into multiple analyses

### Issue: Clustering doesn't make sense
**Cause**: Noise, batch effects, or inappropriate method

**Solution**:
- Use Spearman instead of Pearson
- Filter to specific cell type/condition
- Check for batch effects in data
- Verify data is properly normalized

## Best Practices

1. **Start with established gene sets**
   - Use known pathway members
   - Include positive/negative controls
   - Test well-characterized gene relationships

2. **Filter appropriately**
   - Use cell type filtering for cell type-specific analyses
   - Remove low-quality cells if available
   - Consider disease/condition context

3. **Choose correlation method wisely**
   - Pearson for linear, normal data
   - Spearman for robust, rank-based analysis
   - Try both and compare

4. **Use clustering for discovery**
   - Enable clustering for exploratory analysis
   - Look for modules of co-expressed genes
   - Validate findings with literature

5. **Validate findings**
   - Check top correlations make biological sense
   - Cross-reference with known pathways
   - Validate in independent datasets

6. **Compare conditions**
   - Run separate analyses for disease vs normal
   - Compare correlation matrices
   - Look for rewired relationships

## Python API Usage

For integration into analysis pipelines:

```python
import scanpy as sc
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
adata = sc.read_h5ad("data.h5ad")

# Get genes
genes = ["GENE1", "GENE2", "GENE3"]
gene_indices = [adata.var_names.get_loc(g) for g in genes]

# Extract expression
if sparse.issparse(adata.X):
    expr = adata.X[:, gene_indices].toarray()
else:
    expr = adata.X[:, gene_indices]

# Compute correlations
corr_matrix = np.corrcoef(expr.T)  # Pearson
corr_df = pd.DataFrame(corr_matrix, index=genes, columns=genes)

# Plot
sns.clustermap(corr_df, cmap='RdBu_r', center=0, vmin=-1, vmax=1)
plt.savefig("heatmap.png", dpi=300)
```

## Additional Resources

- **CELLxGENE Census**: https://chanzuckerberg.github.io/cellxgene-census/
- **Scanpy**: https://scanpy.readthedocs.io/
- **Gene Ontology**: http://geneontology.org/
- **KEGG Pathways**: https://www.genome.jp/kegg/pathway.html
- **Reactome**: https://reactome.org/
