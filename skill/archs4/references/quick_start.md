# ARCHS4 Quick Start Guide

**Version:** 1.1
**Get started with ARCHS4 in 5 minutes**

---

## Installation

The `archs4_api.py` script is included in this skill under `scripts/`.

**Requirements:**
```bash
pip install requests pandas plotly numpy scipy difflib
```

**Import:**
```python
from archs4_api import (
    gene_expression,
    gene_expression_analysis,
    gene_correlation,
    gene_correlation_pairwise,
    diffexp,
    quicksearch_metadata,
    knn_gene_sets
)
```

---

## 5-Minute Tutorial

### Step 1: Single Gene Expression (30 seconds)

```python
from archs4_api import gene_expression

# Where is TP53 expressed?
df = gene_expression('TP53')

# View top 5 tissues
print(df[['tissue', 'median', 'mean']].head(5))
```

**Output:**
```
                  tissue  median   mean
HEPATOCYTE             15.23   15.45
LIVER                  15.01   15.12
KIDNEY                 12.34   12.50
...
```

---

### Step 2: Multi-Gene Analysis (1 minute)

```python
from archs4_api import gene_expression_analysis

# Analyze multiple genes with visualization
results = gene_expression_analysis(
    genes=['ALB', 'AFP', 'TTR'],
    output_prefix='liver_markers'
)

# Check results
print(f"Analyzed {results['n_tissues']} tissues")
print(f"Top tissue: {results['top_tissues'][0]}")
```

**Output:**
- `liver_markers_expression_table.csv`
- `liver_markers_boxplot.html` (open in browser)

---

### Step 3: Tissue Filtering (1 minute)

```python
from archs4_api import gene_expression_analysis

# Focus on liver tissues only
results = gene_expression_analysis(
    genes=['ALB'],
    tissue_filter='liver',  # NEW: substring matching
    output_prefix='alb_liver'
)

# NEW v1.1: Multiple filters with OR logic
results = gene_expression_analysis(
    genes=['ALB'],
    tissue_filter=['liver', 'hepatocyte', 'bile'],
    output_prefix='alb_comprehensive'
)
```

**Matches:** All tissues containing 'liver' OR 'hepatocyte' OR 'bile'

---

### Step 4: Find Correlated Genes (1 minute)

```python
from archs4_api import gene_correlation

# What genes are correlated with KRT14 in keratinocytes?
results = gene_correlation(
    gene='KRT14',
    meta='keratinocyte',
    k=50  # Top 50 genes
)

# View results
print(results[['gene', 'correlation']].head(10))
```

**Output:**
```
    gene  correlation
KRT5       0.95
KRT6A      0.92
KRT17      0.89
...
```

---

### Step 5: Pairwise Correlations (1 minute)

```python
from archs4_api import gene_correlation_pairwise

# Correlations between gene family members
matrix, pairwise = gene_correlation_pairwise(
    genes=['KRT14', 'KRT5', 'KRT6A', 'KRT17'],
    meta='keratinocyte',
    output_prefix='keratin_pairwise'
)

# View matrix
print(matrix.round(3))
```

**Output:**
- Correlation matrix
- Pairwise results with p-values
- Interactive heatmap (HTML)

---

### Step 6: Differential Expression (30 seconds)

```python
from archs4_api import diffexp

# What genes are differentially expressed in keratinocytes?
results = diffexp(
    gene='KRT14',
    meta='keratinocyte',
    fdr_cutoff=0.05
)

# View top results
print(results[['gene', 'fdr', 'log2_fold_change']].head(10))
```

---

## Common Patterns

### Pattern 1: Cell Type Marker Discovery

```python
# Step 1: Check if gene is cell-specific
df = gene_expression('GENE_NAME')
top_tissue = df.iloc[0]['tissue']

# Step 2: Find correlated genes in that tissue
corr = gene_correlation('GENE_NAME', top_tissue.split('.')[-1].lower())

# Step 3: Check differential expression
diff = diffexp('GENE_NAME', top_tissue.split('.')[-1].lower())
```

---

### Pattern 2: Disease-Specific Analysis

```python
from archs4_api import quicksearch_metadata, gene_correlation

# Step 1: Find disease samples
samples = quicksearch_metadata('cancer', species='human')
sample_ids = samples.iloc[:, 0].tolist()

# Step 2: Analyze in disease context
results = gene_correlation(
    gene='TP53',
    meta='cancer',
    samples=sample_ids[:500]
)
```

---

### Pattern 3: Gene Family Co-Expression

```python
from archs4_api import gene_correlation_pairwise

# Define family
genes = ['GENE1', 'GENE2', 'GENE3', 'GENE4']

# Analyze co-expression
matrix, pairwise = gene_correlation_pairwise(
    genes=genes,
    meta='cell_type',
    generate_heatmap=True
)

# Find highly correlated pairs
high_corr = pairwise[pairwise['correlation'] > 0.8]
```

---

## Tissue Filtering Guide

### Simple Filters

```python
# Organ-level
tissue_filter='liver'

# Cell type
tissue_filter='hepatocyte'

# System-level
tissue_filter='immune'  # All 13 immune tissues
```

### Multiple Filters (v1.1)

```python
# Comprehensive coverage with OR logic
tissue_filter=['liver', 'hepatocyte', 'bile']

# Multiple cell types
tissue_filter=['T cell', 'B cell', 'NK cell']

# Brain tissues
tissue_filter=['brain', 'neuron', 'astrocyte']
```

### Smart Suggestions (v1.1)

If you make a typo, the system helps:

```python
tissue_filter='livar'  # Typo
```

**Output:**
```
'livar': No matches found
  💡 Did you mean: liver, ovary, muscular?
```

---

## Key Functions Reference

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `gene_expression()` | Single gene atlas | gene, species |
| `gene_expression_analysis()` | Multi-gene analysis | genes, tissue_filter, n_tissues |
| `gene_correlation()` | Find correlated genes | gene, meta, k, samples |
| `gene_correlation_pairwise()` | Pairwise correlations | genes, meta, generate_heatmap |
| `diffexp()` | Differential expression | gene, meta, fdr_cutoff |
| `quicksearch_metadata()` | Find samples | query, species |
| `knn_gene_sets()` | Similar gene sets | query_gene_set, top_n |

---

## Available Tissues (72 total)

### Major Organ Systems

- **Nervous System (14):** brain, neuron, astrocyte, oligodendrocyte, microglia
- **Immune System (13):** T cell, B cell, NK cell, monocyte, macrophage, dendritic
- **Digestive System (12):** liver, hepatocyte, pancreas, stomach, intestine, colon
- **Urogenital (9):** kidney, breast, ovary, testis, bladder
- **Connective (6):** fibroblast, adipocyte, stromal, chondrocyte, bone
- **Skin (6):** keratinocyte, melanocyte, basal cell
- **Muscular (5):** skeletal muscle, smooth muscle, myoblast
- **Cardiovascular (4):** heart, endothelial, vascular
- **Respiratory (3):** lung, trachea

**Full list:** See `references/filter_reference.md`

---

## Disease & Treatment Filters

**Two-step process:**

```python
# Step 1: Find samples
samples = quicksearch_metadata('disease_or_treatment')
sample_ids = samples.iloc[:, 0].tolist()

# Step 2: Use in analysis
results = gene_correlation('GENE', 'context', samples=sample_ids)
```

**Common terms:**
- **Disease:** cancer, tumor, diabetes, alzheimer, covid, inflammation
- **Treatment:** treatment, drug, inhibitor, knockout, knockdown

---

## Best Practices

### 1. Start with Defaults
```python
gene_expression_analysis(genes=['GENE'])
```
See top 10 tissues before filtering.

### 2. Use Substrings
```python
tissue_filter='hepat'  # Matches HEPATOCYTE, HEPATIC STELLATE CELL
```
Short substrings catch more variations.

### 3. Combine Related Terms
```python
tissue_filter=['liver', 'hepat', 'bile']
```
Comprehensive coverage with OR logic.

### 4. Check Matched Tissues
```python
print(results['top_tissues'])
```
Verify your filter worked as expected.

### 5. Iterate
Refine filters based on results.

---

## Common Mistakes

### ❌ Not Using Multiple Filters

```python
# Inefficient: separate queries
result1 = gene_expression_analysis(genes=['ALB'], tissue_filter='liver')
result2 = gene_expression_analysis(genes=['ALB'], tissue_filter='hepatocyte')
```

### ✅ Use Multiple Filters

```python
# Efficient: single query
results = gene_expression_analysis(
    genes=['ALB'],
    tissue_filter=['liver', 'hepatocyte', 'bile']
)
```

---

### ❌ Full Hierarchical Paths

```python
tissue_filter='System.Digestive System.Liver.HEPATOCYTE'
```

### ✅ Simple Substrings

```python
tissue_filter='hepatocyte'  # Much easier!
```

---

### ❌ Exact Case

```python
tissue_filter='HEPATOCYTE'  # Works but unnecessary
```

### ✅ Lowercase

```python
tissue_filter='hepatocyte'  # Case-insensitive matching
```

---

## Troubleshooting

### Issue: "No matches found"

**Solution:** Check smart suggestions in output
```python
# System automatically suggests alternatives
tissue_filter='livar'  # Typo

# Output shows:
# 💡 Did you mean: liver, ovary, muscular?
```

### Issue: "API timeout"

**Solution:** Check service status
```python
from archs4_api import get_service_status
status = get_service_status()
```

### Issue: "Not enough samples"

**Solution:** Check sample count
```python
samples = quicksearch_metadata('your_term')
print(f"Found {len(samples)} samples")
```

---

## Performance Tips

### 1. Filter Early
```python
# Faster with tissue_filter
results = gene_expression_analysis(
    genes=['TP53'],
    tissue_filter='liver',  # Process fewer tissues
    n_tissues=5
)
```

### 2. Reuse Sample Lists
```python
# Query once
samples = quicksearch_metadata('cancer')
sample_ids = samples.iloc[:, 0].tolist()

# Use multiple times
for gene in genes:
    results = gene_correlation(gene, 'cancer', samples=sample_ids)
```

### 3. Batch Genes
```python
# Instead of multiple calls
results = gene_expression_analysis(genes=['GENE1', 'GENE2', 'GENE3'])
```

---

## Next Steps

### Learn More
- **Complete API:** `references/api_reference.md`
- **All Filters:** `references/filter_reference.md`
- **More Examples:** `references/query_examples.md`

### Try These
1. Explore expression of your favorite gene
2. Find tissue-specific markers
3. Discover correlated gene sets
4. Analyze disease-specific patterns

---

## Quick Reference Card

```python
# Single gene expression
gene_expression('TP53')

# Multi-gene with filters
gene_expression_analysis(
    genes=['ALB', 'AFP'],
    tissue_filter=['liver', 'hepatocyte']
)

# Find correlated genes
gene_correlation('KRT14', 'keratinocyte', k=50)

# Pairwise correlations
gene_correlation_pairwise(
    genes=['GENE1', 'GENE2', 'GENE3'],
    meta='cell_type'
)

# Differential expression
diffexp('MARKER', 'cell_type', fdr_cutoff=0.05)

# Find disease samples
samples = quicksearch_metadata('disease')
```

---

## Support

For detailed information:
- **API Reference:** `references/api_reference.md`
- **Filter Catalog:** `references/filter_reference.md`
- **14+ Examples:** `references/query_examples.md`

**Version:** 1.1 (2025-12-17)
**Features:** Multiple filters, smart suggestions, 72 tissues, 15 API functions
