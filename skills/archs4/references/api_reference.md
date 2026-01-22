# ARCHS4 API Reference

**Version:** 1.1
**Last Updated:** 2025-12-17

## Overview

Complete API reference for the ARCHS4 gene expression database wrapper. The ARCHS4 API provides access to processed RNA-seq data from GEO, enabling gene expression analysis across tissues, cell types, diseases, and treatments.

**Base URL:** https://maayanlab.cloud
**Species Support:** Human, Mouse

---

## API Functions

### 1. get_service_status()

Check if the ARCHS4 API service is online.

**Signature:**
```python
get_service_status() -> dict
```

**Returns:** Service status information

**Example:**
```python
status = get_service_status()
print(status)
```

---

### 2. quicksearch_metadata(query, species='human')

Search for samples by any metadata term (disease, treatment, GEO series, etc.).

**Signature:**
```python
quicksearch_metadata(query: str, species: str = 'human') -> pd.DataFrame
```

**Parameters:**
- `query` (str): Search term (disease, treatment, GSE series, etc.)
- `species` (str): 'human' or 'mouse' (default: 'human')

**Returns:** DataFrame of matching samples with metadata

**Example:**
```python
# Find cancer samples
samples = quicksearch_metadata('cancer', species='human')
sample_ids = samples.iloc[:, 0].tolist()
```

**Use Cases:**
- Find disease-specific samples
- Find treatment samples
- Query specific GEO series
- Build custom sample sets

---

### 3. gene_expression(gene, species='human')

Get tissue expression atlas for a single gene.

**Signature:**
```python
gene_expression(gene: str, species: str = 'human') -> pd.DataFrame
```

**Parameters:**
- `gene` (str): Gene symbol (e.g., 'TP53')
- `species` (str): 'human' or 'mouse' (default: 'human')

**Returns:** DataFrame with tissue expression statistics (min, q1, median, q3, max, mean)

**Example:**
```python
# Get TP53 expression across tissues
df = gene_expression('TP53')
print(df[['tissue', 'median', 'mean']].head())
```

**Output Columns:** tissue, min, q1, median, q3, max, mean

---

### 4. gene_expression_analysis(genes, species='human', tissue_filter=None, n_tissues=10, output_prefix=None, output_dir=None, generate_plot=True)

**NEW v1.1:** Multi-gene expression analysis with smart filter suggestions.

**Signature:**
```python
gene_expression_analysis(
    genes: Union[str, List[str]],
    species: str = 'human',
    tissue_filter: Union[str, List[str], None] = None,
    n_tissues: int = 10,
    output_prefix: Optional[str] = None,
    output_dir: Optional[str] = None,
    generate_plot: bool = True
) -> dict
```

**Parameters:**
- `genes` (str or List[str]): Gene symbol(s) to analyze
- `species` (str): 'human' or 'mouse' (default: 'human')
- `tissue_filter` (str, List[str], or None): Filter tissues by name(s) - **supports multiple filters with OR logic**
- `n_tissues` (int): Number of top tissues to show (default: 10)
- `output_prefix` (str or None): Prefix for output files
- `output_dir` (str or None): Directory where output files will be saved. If None, uses current working directory (where Python was launched). Can be absolute or relative path.
- `generate_plot` (bool): Generate interactive HTML boxplot (default: True)

**Returns:** Dictionary with:
- `combined_df`: Integrated expression table
- `top_tissues`: List of selected tissues
- `n_tissues`: Number of tissues analyzed

**Filtering Features (v1.1):**
- **Multiple filters:** `tissue_filter=['liver', 'hepatocyte', 'bile']`
- **Substring matching:** Case-insensitive, primary method
- **Fuzzy matching:** 40% similarity threshold (fallback)
- **Smart suggestions:** Shows "Did you mean...?" when no matches found

**Example 1: Single filter**
```python
results = gene_expression_analysis(
    genes=['TP53'],
    tissue_filter='liver'
)
```

**Example 2: Multiple filters (NEW)**
```python
results = gene_expression_analysis(
    genes=['ALB', 'AFP', 'TTR'],
    tissue_filter=['liver', 'hepatocyte', 'bile'],  # OR logic
    output_prefix='liver_markers'
)
```

**Example 3: Typo tolerance**
```python
# Typo 'livar' triggers smart suggestions
results = gene_expression_analysis(
    genes=['ALB'],
    tissue_filter='livar'  # System suggests: "Did you mean: liver?"
)
```

**Example 4: Custom output directory**
```python
results = gene_expression_analysis(
    genes=['ALB', 'AFP', 'TTR'],
    tissue_filter='liver',
    output_prefix='liver_markers',
    output_dir='/path/to/output'  # Files saved to specified directory
)
```

**Outputs:**
- `{prefix}_expression_table.csv`: Long-format table with all data
- `{prefix}_boxplot.html`: Interactive Plotly visualization

---

### 5. gene_correlation(gene, meta, species='human', k=200, samples=None, filter_genes=None)

Find genes correlated with a query gene in specific samples.

**Signature:**
```python
gene_correlation(
    gene: str,
    meta: str,
    species: str = 'human',
    k: int = 200,
    samples: Optional[List[str]] = None,
    filter_genes: Optional[List[str]] = None
) -> pd.DataFrame
```

**Parameters:**
- `gene` (str): Gene symbol to examine
- `meta` (str): Metadata filter/search term
- `species` (str): 'human' or 'mouse' (default: 'human')
- `k` (int): Number of top correlated genes (default: 200)
- `samples` (List[str] or None): Specific GSM sample IDs
- `filter_genes` (List[str] or None): Specific genes to filter for

**Returns:** DataFrame of correlated genes with correlation coefficients

**Example:**
```python
# Find genes correlated with TP53 in cancer
results = gene_correlation('TP53', 'cancer', k=100)
```

**Note:** API limits to ~200 samples internally (server-side limitation)

---

### 6. gene_correlation_pairwise(genes, meta, species='human', samples=None, n_samples=None, output_prefix=None, output_dir=None, generate_heatmap=True)

Calculate pairwise correlations between multiple genes with bootstrap aggregation (v1.2) and bidirectional query support (v1.4).

**Signature:**
```python
gene_correlation_pairwise(
    genes: List[str],
    meta: str,
    species: str = 'human',
    samples: Optional[List[str]] = None,
    n_samples: Optional[int] = None,
    output_prefix: Optional[str] = None,
    output_dir: Optional[str] = None,
    generate_heatmap: bool = True
) -> Dict[str, Any]
```

**Parameters:**
- `genes` (List[str]): List of gene symbols (minimum 2)
- `meta` (str): Metadata filter/search term
- `species` (str): 'human' or 'mouse' (default: 'human')
- `samples` (List[str] or None): Specific GSM sample IDs
- `n_samples` (int or None): Number of samples. If > 200, uses bootstrap aggregation
- `output_prefix` (str or None): Prefix for output files
- `output_dir` (str or None): Directory where output files will be saved. If None, uses current working directory (where Python was launched). Can be absolute or relative path.
- `generate_heatmap` (bool): Generate interactive heatmap (default: True)

**Returns:** Dictionary with keys:
- `correlation_matrix`: Square correlation matrix DataFrame
- `pairwise_results`: Long-format DataFrame with p-values and n_batches
- `p_value_matrix`: Square p-value matrix DataFrame
- `heatmap_file`: Path to HTML heatmap file
- `matrix_file`: Path to CSV matrix file
- `n_batches`: Number of bootstrap batches used (v1.2)
- `total_samples_used`: Total unique samples analyzed (v1.2)

**Bidirectional Query Logic (v1.4):**
Automatically handles asymmetric ARCHS4 correlation matrices by trying both query directions:
1. First tries: `gene_correlation(gene1, filter_genes=[gene2])`
2. If empty, tries: `gene_correlation(gene2, filter_genes=[gene1])`
3. Uses result from whichever direction succeeds
4. Logs which direction worked for transparency

This increases success rate by ~3x for asymmetric gene pairs.

**Bootstrap Aggregation (v1.2):**
When `n_samples > 200`, samples are automatically divided into batches of 200 and processed separately. Correlations are aggregated using Fisher's z-transformation for robust estimates.

**Example:**
```python
# Default: ~200 samples
results = gene_correlation_pairwise(
    genes=['KRT14', 'KRT5', 'KRT6A', 'KRT17'],
    meta='keratinocyte',
    output_prefix='keratin_pairwise'
)
matrix = results['correlation_matrix']
pairwise = results['pairwise_results']

# Bootstrap aggregation with 1600 samples (8 batches)
results = gene_correlation_pairwise(
    genes=['TYK2', 'JAK1'],
    meta='intestine',
    n_samples=1600,
    output_prefix='intestine_analysis'
)
print(f"Used {results['n_batches']} batches with {results['total_samples_used']} samples")

# Save files to specific directory
results = gene_correlation_pairwise(
    genes=['TYK2', 'JAK1'],
    meta='intestine',
    output_prefix='intestine_analysis',
    output_dir='/path/to/output'  # Files saved to specified directory
)

# Bidirectional query example (v1.4)
# Function automatically tries both directions
results = gene_correlation_pairwise(
    genes=['CDKN2D', 'ITGA4'],
    meta='intestine',
    output_prefix='asymmetric_pair'
)

# Output shows diagnostic logging:
# Computing correlation: CDKN2D - ITGA4
#   Forward query empty, trying reverse direction...
#   ✓ Reverse direction succeeded: r=-0.1153
```

**Outputs:**
- `{prefix}_matrix.csv`: Correlation matrix
- `{prefix}_pairwise.csv`: Pairwise results with p-values, n_samples, n_batches
- `{prefix}_heatmap.html`: Interactive clustered heatmap

**Diagnostic Output (v1.4):**
The function logs which query direction succeeded:
- `r=0.4620 (n=200)` - Forward direction worked
- `r=-0.1153 (n=200) [reversed]` - Reverse direction worked
- `Both directions failed - no correlation data available` - Neither direction worked

---

### 7. diffexp(gene, meta, species='human', fdr_cutoff=0.1)

Differential expression analysis for samples matching metadata.

**Signature:**
```python
diffexp(
    gene: str,
    meta: str,
    species: str = 'human',
    fdr_cutoff: float = 0.1
) -> pd.DataFrame
```

**Parameters:**
- `gene` (str): Gene symbol
- `meta` (str): Metadata filter/search term
- `species` (str): 'human' or 'mouse' (default: 'human')
- `fdr_cutoff` (float): FDR significance cutoff (default: 0.1)

**Returns:** DataFrame of differentially expressed genes

**Example:**
```python
# Differential expression in keratinocytes
results = diffexp('KRT14', 'keratinocyte', fdr_cutoff=0.05)
top_genes = results.head(20)
```

**Output Columns:** gene, fdr, t-statistic, log2_fold_change

---

### 8. knn_gene_sets(query_gene_set, species='human', top_n=20)

Find similar gene sets using k-NN search.

**Signature:**
```python
knn_gene_sets(
    query_gene_set: List[str],
    species: str = 'human',
    top_n: int = 20
) -> pd.DataFrame
```

**Parameters:**
- `query_gene_set` (List[str]): List of gene symbols
- `species` (str): 'human' or 'mouse' (default: 'human')
- `top_n` (int): Number of similar sets to return (default: 20)

**Returns:** DataFrame of similar gene sets with similarity scores

**Example:**
```python
# Find similar gene sets for T cell markers
gene_set = ['CD3D', 'CD3E', 'CD4', 'CD8A']
similar = knn_gene_sets(gene_set, top_n=10)
```

---

## Tissue & Cell Type Filtering

### Hierarchical Structure

Tissues follow format: `System.Organ.Cell_Type`

Example: `System.Immune System.Lymphoid.TLYMPHOCYTE`

### Filtering Strategies

#### 1. Substring Matching (Primary)
Case-insensitive substring search.

**Examples:**
- `'hepat'` → matches HEPATOCYTE, HEPATIC STELLATE CELL
- `'lymph'` → matches TLYMPHOCYTE, BLYMPHOCYTE
- `'neuro'` → matches NEURON, NEURONAL CELL

#### 2. Fuzzy Matching (Fallback)
40% similarity threshold for typo tolerance.

**Examples:**
- `'hepatocyt'` → matches HEPATOCYTE
- `'keritinocyte'` → matches KERATINOCYTE

#### 3. Smart Suggestions (v1.1 - Help)
30% similarity threshold, shows closest terms.

**Example:**
- Input: `'livar'`
- Output: "💡 Did you mean: liver, ovary, muscular?"

### Multiple Filters (v1.1)

Combine multiple terms with OR logic:

```python
tissue_filter=['liver', 'hepatocyte', 'bile']
```

Matches: All tissues containing 'liver' OR 'hepatocyte' OR 'bile'

---

## Disease & Treatment Filtering

Disease and treatment filters require a two-step process:

**Step 1:** Find samples using `quicksearch_metadata()`
```python
samples = quicksearch_metadata('cancer', species='human')
sample_ids = samples.iloc[:, 0].tolist()
```

**Step 2:** Use samples in analysis
```python
results = gene_correlation('TP53', 'cancer', samples=sample_ids[:500])
```

**Available Disease Terms:**
- cancer, tumor, carcinoma, adenocarcinoma
- lymphoma, leukemia
- diabetes, alzheimer, parkinson
- covid, inflammation, infection, fibrosis
- normal, healthy (controls)

**Available Treatment Terms:**
- treatment, treated, untreated, control
- drug, inhibitor, stimulation
- knockout, knockdown, overexpression

---

## Best Practices

### 1. Start with Defaults
```python
gene_expression_analysis(genes=['TP53'])
```
See what's highly expressed before filtering.

### 2. Use Broad Filters First
```python
tissue_filter='immune'
```
Explore organ system before narrowing.

### 3. Combine Related Terms
```python
tissue_filter=['liver', 'hepat', 'bile']
```
Comprehensive tissue coverage with OR logic.

### 4. Use Short Substrings
```python
'hepat' instead of 'hepatocyte'
```
Catches more variations.

### 5. Check Matched Tissues
```python
print(results['top_tissues'])
```
Verify filters matched expected tissues.

### 6. Iterate Based on Results
Adjust filters after seeing initial matches.

---

## Limitations

### API Sample Limit (MITIGATED in v1.2)
ARCHS4 API limits correlation to ~200 samples per query (server-side).

**Impact:** Each individual API call can only use 200 samples.

**Mitigation (v1.2):** Bootstrap aggregation automatically divides samples into batches of 200, processes each independently, and aggregates using Fisher's z-transformation. This allows analyzing full datasets (e.g., 1600+ samples).

### Asymmetric Correlation Matrix (FIXED in v1.4)
Some gene pairs return correlations in one direction but not the other due to ARCHS4's precomputed correlation matrix structure.

**Impact:** `gene_correlation('A', filter_genes=['B'])` may return empty while `gene_correlation('B', filter_genes=['A'])` returns data.

**Fix (v1.4):** `gene_correlation_pairwise()` automatically tries both directions and uses whichever succeeds. Increases success rate by ~3x.

### Zero/Weak Correlations (Biological)
Gene pairs with correlation below threshold (|r| < ~0.06) are not stored in ARCHS4's precomputed matrices.

**Impact:** Some gene pairs return empty in both directions (e.g., genes in different pathways).

**This is expected:** Not all genes are correlated. Use alternative databases (STRING, BioGRID) or compute from raw data if needed.

### Disease/Treatment Filtering
Not directly available in `gene_expression_analysis()`.

**Workaround:** Use `quicksearch_metadata()` + `samples` parameter.

### Hierarchical Tissue Names
Tissue names are hierarchical paths.

**Workaround:** Use substring matching to catch any part of path.

---

## Performance

- Single gene query: ~2-3 seconds
- Multiple genes: ~5-20 seconds (depends on count)
- Suggestion overhead: < 0.1 seconds (negligible)

---

## Version History

### v1.4 (2025-12-17)
- ✅ Bidirectional correlation queries in gene_correlation_pairwise()
- ✅ Automatically tries reverse direction when forward query returns empty
- ✅ Handles asymmetric ARCHS4 correlation matrix indexing
- ✅ 3x improvement in correlation retrieval success rate
- ✅ Diagnostic logging shows which direction succeeded

### v1.3 (2025-12-17)
- ✅ Added output_dir parameter to gene_expression_analysis()
- ✅ Added output_dir parameter to gene_correlation_pairwise()
- ✅ Files now saved to working directory by default
- ✅ Support for custom output directories

### v1.2 (2025-12-17)
- ✅ Bootstrap aggregation for n_samples > 200
- ✅ Overcomes 200-sample API limit via batch processing
- ✅ Fisher's z-transformation for robust correlation aggregation
- ✅ Support for analyzing full datasets (e.g., 1600+ intestinal samples)
- ✅ n_batches and total_samples_used in results

### v1.1 (2025-12-17)
- ✅ Smart filter suggestions
- ✅ Multiple tissue filters with OR logic
- ✅ Schema-aware vocabulary
- ✅ "Did you mean?" messages

### v1.0 (2025-12-16)
- ✅ 15 API functions
- ✅ Gene expression analysis
- ✅ Pairwise correlations
- ✅ Differential expression
- ✅ Interactive visualizations
