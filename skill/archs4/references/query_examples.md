# ARCHS4 Query Examples

**Version:** 1.1
**Last Updated:** 2025-12-17

Complete collection of working examples for common ARCHS4 queries.

---

## Basic Expression Queries

### Example 1: Single Gene Expression Atlas

Query tissue expression for a single gene across all tissues.

```python
from archs4_api import gene_expression

# Get TP53 expression across all tissues
df = gene_expression('TP53', species='human')

# View top expressing tissues
print(df[['tissue', 'median', 'mean']].head(10))

# Save results
df.to_csv('tp53_expression.csv', index=False)
```

**Output:**
- DataFrame with columns: tissue, min, q1, median, q3, max, mean
- 72 rows (one per tissue)

**Use case:** Understand where a gene is normally expressed.

---

### Example 2: Multi-Gene Expression Analysis

Analyze multiple genes with integrated visualization.

```python
from archs4_api import gene_expression_analysis

# Analyze keratin family expression
results = gene_expression_analysis(
    genes=['KRT14', 'KRT5', 'KRT6A', 'KRT17'],
    output_prefix='keratin_family',
    generate_plot=True
)

# Access results
print(f"Analyzed {results['n_tissues']} tissues")
print(f"Top tissue: {results['top_tissues'][0]}")
```

**Output:**
- `keratin_family_expression_table.csv` - Long-format data
- `keratin_family_boxplot.html` - Interactive visualization

**Use case:** Compare expression patterns of related genes.

---

### Example 3: Single Tissue Filter

Filter to specific tissue type.

```python
from archs4_api import gene_expression_analysis

# ALB expression in liver tissues only
results = gene_expression_analysis(
    genes=['ALB'],
    tissue_filter='liver',
    output_prefix='alb_liver'
)

# Check which tissues matched
print("Matched tissues:")
for tissue in results['top_tissues']:
    print(f"  - {tissue}")
```

**Expected matches:**
- HEPATOCYTE
- LIVER
- HEPATIC STELLATE CELL

**Use case:** Focus analysis on relevant tissue context.

---

### Example 4: Multiple Tissue Filters (v1.1)

**NEW:** Combine multiple filter terms with OR logic.

```python
from archs4_api import gene_expression_analysis

# Comprehensive liver analysis
results = gene_expression_analysis(
    genes=['ALB', 'AFP', 'TTR'],
    tissue_filter=['liver', 'hepatocyte', 'bile'],
    output_prefix='liver_comprehensive'
)
```

**Matches:** All tissues containing 'liver' OR 'hepatocyte' OR 'bile'

**Use case:** Ensure comprehensive coverage of related tissues.

---

### Example 5: System-Level Filter

Filter by entire organ system.

```python
from archs4_api import gene_expression_analysis

# Analyze immune system expression
results = gene_expression_analysis(
    genes=['CD3D', 'CD4', 'CD8A'],
    tissue_filter='immune',
    n_tissues=15,
    output_prefix='tcell_markers'
)
```

**Matches:** All 13 immune system tissues

**Use case:** Survey expression across an entire physiological system.

---

## Correlation Queries

### Example 6: Gene Correlation

Find genes correlated with a query gene in specific context.

```python
from archs4_api import gene_correlation

# Find genes correlated with KRT14 in keratinocytes
results = gene_correlation(
    gene='KRT14',
    meta='keratinocyte',
    species='human',
    k=100  # Top 100 correlated genes
)

# View top correlations
print(results[['gene', 'correlation']].head(20))

# Filter for specific genes of interest
keratins = results[results['gene'].str.startswith('KRT')]
print(f"\nFound {len(keratins)} correlated keratins")
```

**Output:** DataFrame with correlated genes and correlation coefficients

**Use case:** Discover co-expressed genes in specific cell types.

---

### Example 7: Pairwise Gene Correlations

Calculate correlations between multiple genes.

```python
from archs4_api import gene_correlation_pairwise

# Pairwise correlations for keratin family
matrix, pairwise = gene_correlation_pairwise(
    genes=['KRT14', 'KRT5', 'KRT6A', 'KRT17'],
    meta='keratinocyte',
    output_prefix='keratin_pairwise',
    generate_heatmap=True
)

# View correlation matrix
print("\nCorrelation Matrix:")
print(matrix.round(3))

# Find highest correlations
top_pairs = pairwise.nlargest(5, 'correlation')
print("\nTop 5 Correlations:")
print(top_pairs[['gene1', 'gene2', 'correlation', 'p_value']])
```

**Output:**
- Correlation matrix (DataFrame)
- Pairwise results with p-values (DataFrame)
- Interactive clustered heatmap (HTML)

**Use case:** Understand co-regulation patterns within gene families.

---

### Example 8: Correlation with Custom Samples

Use disease-specific samples for correlation analysis.

```python
from archs4_api import quicksearch_metadata, gene_correlation

# Step 1: Find cancer samples
cancer_samples = quicksearch_metadata('liver cancer', species='human')
sample_ids = cancer_samples.iloc[:, 0].tolist()

print(f"Found {len(sample_ids)} cancer samples")

# Step 2: Correlation in cancer context
results = gene_correlation(
    gene='TP53',
    meta='liver cancer',
    samples=sample_ids[:500],  # Use first 500
    k=200
)

# View results
print(f"\nTop genes correlated with TP53 in liver cancer:")
print(results[['gene', 'correlation']].head(10))
```

**Use case:** Disease-specific co-expression analysis.

---

## Differential Expression

### Example 9: Differential Expression Analysis

Find genes differentially expressed in specific cell types.

```python
from archs4_api import diffexp

# Differential expression in keratinocytes
results = diffexp(
    gene='KRT14',
    meta='keratinocyte',
    species='human',
    fdr_cutoff=0.05
)

# View top results
print(f"Found {len(results)} significant genes")
print("\nTop 20 differentially expressed:")
print(results[['gene', 'fdr', 'log2_fold_change']].head(20))

# Filter by fold change
upregulated = results[results['log2_fold_change'] > 2]
print(f"\n{len(upregulated)} genes with >2-fold upregulation")
```

**Output:** DataFrame with gene, FDR, t-statistic, log2 fold change

**Use case:** Identify cell type-specific markers.

---

## Disease-Specific Queries

### Example 10: Cancer Sample Analysis

Two-step workflow for disease-specific analysis.

```python
from archs4_api import quicksearch_metadata, gene_expression_analysis, gene_correlation

# Step 1: Find cancer samples
samples = quicksearch_metadata('breast cancer', species='human')
sample_ids = samples.iloc[:, 0].tolist()

print(f"Found {len(sample_ids)} breast cancer samples")

# Step 2a: Gene correlation in cancer
corr_results = gene_correlation(
    gene='BRCA1',
    meta='breast cancer',
    samples=sample_ids[:500],
    k=100
)

# Step 2b: Differential expression
diff_results = diffexp(
    gene='BRCA1',
    meta='breast cancer',
    fdr_cutoff=0.1
)

print(f"\nTop genes correlated with BRCA1:")
print(corr_results[['gene', 'correlation']].head(10))

print(f"\nTop differentially expressed:")
print(diff_results[['gene', 'fdr', 'log2_fold_change']].head(10))
```

**Use case:** Comprehensive disease-specific gene analysis.

---

### Example 11: COVID-19 Sample Analysis

```python
from archs4_api import quicksearch_metadata, gene_correlation

# Find COVID-19 samples
covid_samples = quicksearch_metadata('covid', species='human')
sample_ids = covid_samples.iloc[:, 0].tolist()

print(f"Found {len(sample_ids)} COVID-19 samples")

# Analyze ACE2 receptor correlations
results = gene_correlation(
    gene='ACE2',
    meta='covid',
    samples=sample_ids,
    k=200
)

# Filter for immune genes
immune_genes = results[results['gene'].str.contains('IL|TNF|IFN|CCL|CXCL', case=False)]
print(f"\nFound {len(immune_genes)} immune-related genes")
print(immune_genes[['gene', 'correlation']].head(20))
```

**Use case:** Infection-specific immune response analysis.

---

## Gene Set Queries

### Example 12: k-NN Gene Set Search

Find similar gene sets for your query set.

```python
from archs4_api import knn_gene_sets

# Define T cell marker gene set
tcell_markers = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B']

# Find similar gene sets
similar_sets = knn_gene_sets(
    query_gene_set=tcell_markers,
    species='human',
    top_n=20
)

# View results
print("Similar gene sets:")
print(similar_sets[['name', 'similarity', 'size']].head(10))
```

**Output:** DataFrame of similar gene sets with similarity scores

**Use case:** Discover related gene sets and biological processes.

---

## Multi-Step Workflows

### Workflow 1: Complete Cell Type Analysis

Comprehensive analysis of a cell type marker.

```python
from archs4_api import (
    gene_expression_analysis,
    gene_correlation,
    gene_correlation_pairwise,
    diffexp
)

MARKER = 'KRT14'
CELL_TYPE = 'keratinocyte'
OUTPUT_PREFIX = 'keratinocyte_analysis'

# Step 1: Expression atlas
print("Step 1: Expression atlas")
expr_results = gene_expression_analysis(
    genes=[MARKER],
    tissue_filter='skin',
    output_prefix=f'{OUTPUT_PREFIX}_expression'
)

# Step 2: Find correlated genes
print("\nStep 2: Correlated genes")
corr_results = gene_correlation(
    gene=MARKER,
    meta=CELL_TYPE,
    k=50
)

# Step 3: Top correlated keratins
top_keratins = corr_results[corr_results['gene'].str.startswith('KRT')]['gene'].head(5).tolist()
top_keratins.insert(0, MARKER)

print(f"\nTop correlated keratins: {top_keratins}")

# Step 4: Pairwise correlations
print("\nStep 4: Pairwise correlations")
matrix, pairwise = gene_correlation_pairwise(
    genes=top_keratins,
    meta=CELL_TYPE,
    output_prefix=f'{OUTPUT_PREFIX}_pairwise'
)

# Step 5: Differential expression
print("\nStep 5: Differential expression")
diff_results = diffexp(
    gene=MARKER,
    meta=CELL_TYPE,
    fdr_cutoff=0.05
)

print("\n✓ Analysis complete!")
print(f"  Expression: {OUTPUT_PREFIX}_expression_table.csv")
print(f"  Pairwise: {OUTPUT_PREFIX}_pairwise_matrix.csv")
print(f"  Top DE genes: {len(diff_results)}")
```

**Use case:** Complete characterization of cell type-specific markers.

---

### Workflow 2: Gene Family Analysis

Systematic analysis of a gene family.

```python
from archs4_api import gene_expression_analysis, gene_correlation_pairwise

# Define gene family
gene_family = ['ALB', 'AFP', 'TTR', 'SERPINA1', 'FGA']
tissue_context = ['liver', 'hepatocyte']

# Step 1: Expression patterns
print("Step 1: Expression analysis")
expr_results = gene_expression_analysis(
    genes=gene_family,
    tissue_filter=tissue_context,
    output_prefix='liver_secreted_proteins_expr'
)

# Step 2: Co-expression analysis
print("\nStep 2: Co-expression analysis")
matrix, pairwise = gene_correlation_pairwise(
    genes=gene_family,
    meta='hepatocyte',
    output_prefix='liver_secreted_proteins_corr'
)

# Step 3: Analyze results
print("\nCorrelation Matrix:")
print(matrix.round(3))

highly_correlated = pairwise[pairwise['correlation'] > 0.7]
print(f"\n{len(highly_correlated)} highly correlated pairs (r > 0.7)")
```

**Use case:** Understand co-regulation within gene families.

---

### Workflow 3: Treatment Response Analysis

Compare treated vs control samples.

```python
from archs4_api import quicksearch_metadata, gene_correlation, diffexp

GENE = 'NR3C1'
TREATMENT = 'dexamethasone'

# Step 1: Find treated samples
print("Step 1: Finding samples")
treated_samples = quicksearch_metadata(f'{TREATMENT} treated', species='human')
control_samples = quicksearch_metadata(f'{TREATMENT} control', species='human')

treated_ids = treated_samples.iloc[:, 0].tolist()
control_ids = control_samples.iloc[:, 0].tolist()

print(f"  Treated: {len(treated_ids)} samples")
print(f"  Control: {len(control_ids)} samples")

# Step 2: Correlation in treated
print("\nStep 2: Treated sample correlations")
treated_corr = gene_correlation(
    gene=GENE,
    meta=TREATMENT,
    samples=treated_ids[:500],
    k=100
)

# Step 3: Differential expression
print("\nStep 3: Differential expression")
diff_expr = diffexp(
    gene=GENE,
    meta=TREATMENT,
    fdr_cutoff=0.05
)

print(f"\nTop genes correlated in treated samples:")
print(treated_corr[['gene', 'correlation']].head(10))

print(f"\nTop differentially expressed genes:")
print(diff_expr[['gene', 'fdr', 'log2_fold_change']].head(10))
```

**Use case:** Drug/treatment response profiling.

---

## Common Use Cases

### Use Case 1: Find Tissue-Specific Markers

```python
from archs4_api import gene_expression

# Check if gene is tissue-specific
df = gene_expression('ALB')

# Find top expressing tissue
top_tissue = df.iloc[0]
print(f"Top expression: {top_tissue['tissue']}")
print(f"  Median: {top_tissue['median']:.2f}")
print(f"  Mean: {top_tissue['mean']:.2f}")

# Check specificity
second_tissue = df.iloc[1]
fold_difference = top_tissue['median'] / second_tissue['median']
print(f"\nFold difference vs 2nd: {fold_difference:.2f}x")

if fold_difference > 10:
    print("✓ Highly tissue-specific marker!")
```

---

### Use Case 2: Validate Known Co-Expression

```python
from archs4_api import gene_correlation

# Check if two genes are correlated in expected context
results = gene_correlation(
    gene='CD3D',
    meta='T cell',
    k=200,
    filter_genes=['CD3E', 'CD4', 'CD8A']  # Filter for specific genes
)

print("Co-expression with known T cell markers:")
print(results[['gene', 'correlation']])
```

---

### Use Case 3: Explore Unknown Gene

```python
from archs4_api import gene_expression, gene_correlation, knn_gene_sets

UNKNOWN_GENE = 'ENSG00000...'  # Replace with your gene

# Step 1: Where is it expressed?
expr = gene_expression(UNKNOWN_GENE)
print(f"Top 3 expressing tissues:")
print(expr[['tissue', 'median']].head(3))

# Step 2: What is it correlated with?
# Use top expressing tissue as context
top_tissue = expr.iloc[0]['tissue'].split('.')[-1].lower()
corr = gene_correlation(UNKNOWN_GENE, top_tissue, k=20)
print(f"\nTop correlated genes in {top_tissue}:")
print(corr[['gene', 'correlation']].head(10))

# Step 3: What gene sets is it similar to?
similar = knn_gene_sets([UNKNOWN_GENE], top_n=10)
print(f"\nSimilar gene sets:")
print(similar[['name', 'similarity']].head(5))
```

---

## Troubleshooting Examples

### Handle Typos with Smart Suggestions

```python
from archs4_api import gene_expression_analysis

# Intentional typo - system will suggest correct term
results = gene_expression_analysis(
    genes=['ALB'],
    tissue_filter='livar'  # Typo: should be 'liver'
)

# Output shows:
# 'livar': No matches found
#   💡 Did you mean: liver, ovary, muscular?

# Correct and retry
results = gene_expression_analysis(
    genes=['ALB'],
    tissue_filter='liver'  # Fixed
)
```

---

### Check API Status

```python
from archs4_api import get_service_status

# Check if API is accessible
status = get_service_status()
print(f"Service status: {status}")
```

---

### Verify Sample Count

```python
from archs4_api import quicksearch_metadata

# Check how many samples are available
samples = quicksearch_metadata('your_term_here')
print(f"Found {len(samples)} samples")

# View sample metadata
print(samples.head())
```

---

## Performance Tips

### Tip 1: Batch Gene Queries
```python
# Instead of multiple separate queries
genes = ['TP53', 'MYC', 'BRCA1']

# Use gene_expression_analysis for efficiency
results = gene_expression_analysis(genes=genes)
```

### Tip 2: Limit Tissues Early
```python
# Use tissue_filter to reduce data
results = gene_expression_analysis(
    genes=['TP53'],
    tissue_filter='liver',  # Much faster than processing all 72 tissues
    n_tissues=5
)
```

### Tip 3: Reuse Sample Lists
```python
# Query once, use many times
cancer_samples = quicksearch_metadata('cancer')
sample_ids = cancer_samples.iloc[:, 0].tolist()

# Reuse for multiple genes
for gene in ['TP53', 'MYC', 'BRCA1']:
    results = gene_correlation(gene, 'cancer', samples=sample_ids)
    # Process results...
```

---

## Version History

### v1.1 (2025-12-17)
- ✅ Added multiple tissue filter examples
- ✅ Added smart suggestion examples
- ✅ Added troubleshooting section

### v1.0 (2025-12-16)
- ✅ 14+ complete examples
- ✅ 3 multi-step workflows
- ✅ 6 common use cases
