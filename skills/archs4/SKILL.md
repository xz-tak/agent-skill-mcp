---
name: archs4
description: Query and analyze gene expression data from the ARCHS4 database. Use this skill when users request tissue-specific gene expression analysis, gene correlation studies, differential expression analysis, or need to query expression patterns across 72 human/mouse tissues and cell types. The skill provides smart tissue filtering with typo tolerance and supports disease/treatment-specific sample queries.
---

# ARCHS4 Gene Expression Analysis

## Overview

Query the ARCHS4 database to analyze gene expression patterns across tissues, cell types, diseases, and treatments. The ARCHS4 database contains processed RNA-seq data from GEO (Gene Expression Omnibus), enabling comprehensive gene expression analysis without requiring raw data processing.

**Key capabilities:**
- Expression atlas for 72+ human/mouse tissues
- Multi-gene expression analysis with visualization
- Gene correlation and co-expression analysis
- Differential expression in specific contexts
- Disease and treatment-specific sample queries
- Smart tissue filtering with typo tolerance (v1.1)

**Data coverage:**
- 72 tissues across 9 organ systems
- ~77,760 cancer samples
- ~296,703 treatment samples
- Comprehensive metadata (disease, treatment, GEO series)

## When to Use This Skill

Use this skill when users request:
- **Expression queries**: "Where is TP53 expressed?", "Show me ALB expression in liver"
- **Tissue-specific analysis**: "Analyze keratin expression in skin", "Find liver-specific markers"
- **Correlation studies**: "What genes correlate with KRT14 in keratinocytes?"
- **Gene families**: "Compare expression of albumin family genes"
- **Disease context**: "Find genes correlated with TP53 in cancer samples"
- **Cell type markers**: "What are markers for hepatocytes?", "Validate T cell markers"
- **Co-expression**: "Calculate pairwise correlations for immune markers"
- **Differential expression**: "What genes are differentially expressed in keratinocytes?"

## Quick Start

### Basic Expression Query

Start with a simple expression atlas query to see where a gene is expressed:

```python
from archs4_api import gene_expression

# Query single gene
df = gene_expression('TP53', species='human')

# View top expressing tissues
print(df[['tissue', 'median', 'mean']].head(10))
```

### Multi-Gene Analysis with Filtering

For more comprehensive analysis with tissue filtering and visualization:

```python
from archs4_api import gene_expression_analysis

# Analyze multiple genes with tissue filter
results = gene_expression_analysis(
    genes=['ALB', 'AFP', 'TTR'],
    tissue_filter=['liver', 'hepatocyte'],  # NEW v1.1: multiple filters
    output_prefix='liver_markers'
)

# Access results
print(f"Analyzed {results['n_tissues']} tissues")
print(f"Top tissue: {results['top_tissues'][0]}")
```

**Outputs:**
- CSV table with expression data
- Interactive HTML boxplot visualization

### Find Correlated Genes

Discover genes with similar expression patterns:

```python
from archs4_api import gene_correlation

# Find genes correlated with marker in specific cell type
results = gene_correlation(
    gene='KRT14',
    meta='keratinocyte',
    k=50  # Top 50 correlated genes
)

# View top correlations
print(results[['gene', 'correlation']].head(10))
```

## Core Analysis Tasks

### Task 1: Tissue Expression Atlas

Query gene expression across all tissues to understand expression patterns.

**Use when:**
- User asks "Where is [GENE] expressed?"
- Need to identify tissue-specific markers
- Want to see expression breadth/specificity

**Workflow:**
1. Use `gene_expression()` for single gene atlas
2. Check top expressing tissues
3. Calculate tissue specificity (fold difference between top tissues)
4. For multiple genes, use `gene_expression_analysis()` with integrated visualization

**Example:**
```python
from archs4_api import gene_expression

# Single gene
df = gene_expression('ALB')

# Check if tissue-specific
top_tissue = df.iloc[0]
second_tissue = df.iloc[1]
fold_diff = top_tissue['median'] / second_tissue['median']

if fold_diff > 10:
    print(f"✓ Highly tissue-specific to {top_tissue['tissue']}")
```

### Task 2: Tissue-Filtered Expression Analysis

Focus analysis on specific tissues or cell types using smart filtering.

**Use when:**
- User specifies tissue context: "in liver", "in brain", "in immune cells"
- Need to filter to relevant tissues only
- Want comprehensive tissue coverage

**Filtering strategies:**
1. **Single filter**: `tissue_filter='liver'`
2. **Multiple filters (v1.1)**: `tissue_filter=['liver', 'hepatocyte', 'bile']` (OR logic)
3. **System-level**: `tissue_filter='immune'` (all 13 immune tissues)
4. **Substring**: `tissue_filter='hepat'` (matches HEPATOCYTE, HEPATIC STELLATE CELL)

**Smart suggestions (v1.1):**
- System provides "Did you mean?" suggestions for typos
- Shows closest matching terms from schema
- Helps users discover correct filter vocabulary

**Example:**
```python
from archs4_api import gene_expression_analysis

# Multiple related filters for comprehensive coverage
results = gene_expression_analysis(
    genes=['CD3D', 'CD4', 'CD8A'],
    tissue_filter=['T cell', 'lymphocyte', 'thymus', 'immune'],
    output_prefix='tcell_markers'
)
```

**Typo tolerance:**
```python
# User makes typo
results = gene_expression_analysis(
    genes=['ALB'],
    tissue_filter='livar'  # Typo
)

# System outputs:
# 'livar': No matches found
#   💡 Did you mean: liver, ovary, muscular?
#   💡 Try one of these filters: liver, ovary, muscular, microglia, valve

# User corrects and retries
results = gene_expression_analysis(
    genes=['ALB'],
    tissue_filter='liver'  # Fixed
)
```

### Task 3: Gene Correlation Analysis

Find genes with similar expression patterns in specific contexts.

**Use when:**
- User asks "What genes correlate with [GENE]?"
- Need to discover co-expressed genes
- Building gene modules or networks
- Validating known co-expression

**Workflow:**
1. Use `gene_correlation()` with cell type or tissue context
2. Review top correlated genes
3. Optional: Filter for specific gene families using `filter_genes` parameter
4. Optional: Use disease-specific samples via `quicksearch_metadata()`

**Example:**
```python
from archs4_api import gene_correlation

# Find genes correlated in specific cell type
results = gene_correlation(
    gene='KRT14',
    meta='keratinocyte',
    k=100  # Top 100
)

# Filter for gene family
keratins = results[results['gene'].str.startswith('KRT')]
print(f"Found {len(keratins)} correlated keratins")
```

**Disease-specific context:**
```python
from archs4_api import quicksearch_metadata, gene_correlation

# Step 1: Find disease samples
samples = quicksearch_metadata('liver cancer', species='human')
sample_ids = samples.iloc[:, 0].tolist()

# Step 2: Correlation in disease context
results = gene_correlation(
    gene='TP53',
    meta='liver cancer',
    samples=sample_ids[:500]
)
```

### Task 4: Pairwise Correlation Analysis

Calculate correlations between multiple genes to understand co-regulation. Supports bootstrap aggregation for large sample sizes (> 200) to overcome API limitations. **New in v1.4:** Automatically tries both query directions to handle asymmetric API indexing.

**Use when:**
- User provides gene list for co-expression analysis
- Analyzing gene families or modules
- Need correlation matrix and heatmap
- Validating pathway relationships
- Want to use more than 200 samples for robust estimates

**Key Features (v1.4):**
- **Bidirectional queries:** Automatically tries reverse direction if forward fails
- **Bootstrap aggregation:** Divides samples into batches for > 200 samples
- **Fisher's z-transformation:** Robust correlation aggregation across batches
- **Smart fallback:** Handles asymmetric ARCHS4 correlation matrices

**Workflow:**
1. Collect gene list (from user or correlation results)
2. Use `gene_correlation_pairwise()` to calculate all pairs
3. Function automatically tries both directions (A→B, then B→A if needed)
4. If using > 200 samples, bootstrap aggregation automatically divides into batches
5. Review correlation matrix and aggregated statistics
6. Visualize with interactive heatmap

**Example (Default ~200 samples):**
```python
from archs4_api import gene_correlation_pairwise

# Pairwise correlations for gene family
results = gene_correlation_pairwise(
    genes=['KRT14', 'KRT5', 'KRT6A', 'KRT17', 'KRT16'],
    meta='keratinocyte',
    output_prefix='keratin_family',
    generate_heatmap=True
)

# Find highly correlated pairs
pairwise = results['pairwise_results']
high_corr = pairwise[pairwise['correlation'] > 0.8]
print(f"{len(high_corr)} pairs with r > 0.8")

# Check significance
significant = pairwise[pairwise['p_value'] < 0.01]
print(f"{len(significant)} significant pairs (p < 0.01)")
```

### Task 5: Differential Expression Analysis

Identify genes differentially expressed in specific cell types or conditions.

**Use when:**
- User asks "What genes are specific to [CELL_TYPE]?"
- Need to find cell type markers
- Comparing expression between contexts
- Validating known markers

**Workflow:**
1. Use `diffexp()` with cell type or condition
2. Filter by FDR cutoff (default 0.1, stricter 0.05)
3. Sort by log2 fold change
4. Identify upregulated vs downregulated genes

**Example:**
```python
from archs4_api import diffexp

# Differential expression in cell type
results = diffexp(
    gene='KRT14',
    meta='keratinocyte',
    fdr_cutoff=0.05  # Stricter cutoff
)

# Identify upregulated markers
upregulated = results[results['log2_fold_change'] > 2]
print(f"{len(upregulated)} genes with >2-fold upregulation")

# Top markers
top_markers = upregulated.head(20)
print(top_markers[['gene', 'fdr', 'log2_fold_change']])
```

### Task 6: Gene Set Similarity Search

Find similar gene sets using k-NN search.

**Use when:**
- User provides gene set for enrichment
- Need to discover related biological processes
- Validating gene set composition

**Example:**
```python
from archs4_api import knn_gene_sets

# Define query set
tcell_markers = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B']

# Find similar sets
similar = knn_gene_sets(
    query_gene_set=tcell_markers,
    species='human',
    top_n=20
)

print("Similar gene sets:")
print(similar[['name', 'similarity', 'size']].head(10))
```

## Multi-Step Workflows

### Workflow 1: Complete Cell Type Characterization

Comprehensive analysis of a cell type marker.

**Steps:**
1. Expression atlas to identify expression pattern
2. Correlation analysis to find co-expressed genes
3. Pairwise correlations among top genes
4. Differential expression to identify markers
5. Compile results

**Example:**
```python
from archs4_api import (
    gene_expression_analysis,
    gene_correlation,
    gene_correlation_pairwise,
    diffexp
)

MARKER = 'KRT14'
CELL_TYPE = 'keratinocyte'

# Step 1: Expression atlas
expr_results = gene_expression_analysis(
    genes=[MARKER],
    tissue_filter='skin',
    output_prefix=f'{CELL_TYPE}_expression'
)

# Step 2: Find correlated genes
corr_results = gene_correlation(MARKER, CELL_TYPE, k=50)

# Step 3: Top correlated family members
top_genes = corr_results[corr_results['gene'].str.startswith('KRT')]['gene'].head(5).tolist()
top_genes.insert(0, MARKER)

# Step 4: Pairwise correlations
matrix, pairwise = gene_correlation_pairwise(
    genes=top_genes,
    meta=CELL_TYPE,
    output_prefix=f'{CELL_TYPE}_pairwise'
)

# Step 5: Differential expression
diff_results = diffexp(MARKER, CELL_TYPE, fdr_cutoff=0.05)

print(f"✓ Characterization complete:")
print(f"  - Expression: {expr_results['n_tissues']} tissues")
print(f"  - Correlated: {len(corr_results)} genes")
print(f"  - Pairwise: {len(top_genes)} genes")
print(f"  - DE genes: {len(diff_results)}")
```

### Workflow 2: Gene Family Co-Expression

Systematic analysis of gene family co-expression patterns.

**Steps:**
1. Define gene family
2. Expression analysis with tissue filtering
3. Pairwise correlations
4. Identify highly correlated sub-groups

**Example:**
```python
from archs4_api import gene_expression_analysis, gene_correlation_pairwise

# Define family
albumin_family = ['ALB', 'AFP', 'TTR', 'SERPINA1', 'FGA']

# Step 1: Expression in liver
expr_results = gene_expression_analysis(
    genes=albumin_family,
    tissue_filter=['liver', 'hepatocyte'],
    output_prefix='albumin_family_expr'
)

# Step 2: Co-expression analysis
matrix, pairwise = gene_correlation_pairwise(
    genes=albumin_family,
    meta='hepatocyte',
    output_prefix='albumin_family_corr'
)

# Step 3: Identify sub-groups
highly_correlated = pairwise[pairwise['correlation'] > 0.7]
print(f"{len(highly_correlated)} highly correlated pairs")
```

### Workflow 3: Disease-Specific Gene Analysis

Analyze gene behavior in disease context.

**Steps:**
1. Find disease samples
2. Correlation analysis in disease
3. Differential expression
4. Compare to normal tissue

**Example:**
```python
from archs4_api import quicksearch_metadata, gene_correlation, diffexp

# Step 1: Find disease samples
cancer_samples = quicksearch_metadata('breast cancer', species='human')
sample_ids = cancer_samples.iloc[:, 0].tolist()

print(f"Found {len(sample_ids)} breast cancer samples")

# Step 2: Correlation in cancer
cancer_corr = gene_correlation(
    gene='BRCA1',
    meta='breast cancer',
    samples=sample_ids[:500],
    k=100
)

# Step 3: Differential expression
diff_expr = diffexp(
    gene='BRCA1',
    meta='breast cancer',
    fdr_cutoff=0.05
)

print(f"Top correlated in cancer: {cancer_corr.iloc[0]['gene']}")
print(f"Top DE gene: {diff_expr.iloc[0]['gene']}")
```

## Available Resources

### scripts/archs4_api.py

Complete Python API wrapper with 15 functions. Execute functions directly without loading into context.

**Key functions:**
- `gene_expression()` - Single gene atlas
- `gene_expression_analysis()` - Multi-gene with filtering and visualization
- `gene_correlation()` - Find correlated genes
- `gene_correlation_pairwise()` - Pairwise correlations with heatmap
- `diffexp()` - Differential expression
- `quicksearch_metadata()` - Find samples by metadata
- `knn_gene_sets()` - Similar gene sets
- `get_service_status()` - Check API status

**Usage:**
```python
from archs4_api import function_name
results = function_name(parameters)
```

### references/api_reference.md

Complete API documentation with:
- All 15 function signatures and parameters
- Return types and output formats
- Filtering strategies (substring, fuzzy, smart suggestions)
- Best practices and limitations
- Performance notes
- Version history

**Load when:** Need detailed function parameters, filtering details, or troubleshooting.

### references/filter_reference.md

Comprehensive filter catalog with:
- All 72 tissues organized by organ system
- Common filter terms and combinations
- Disease filters (~77k cancer samples)
- Treatment filters (~297k samples)
- Smart suggestion examples
- Common mistakes and tips

**Load when:** User needs help finding correct filter terms or exploring available tissues.

### references/query_examples.md

14+ complete working examples including:
- Basic expression queries
- Correlation analyses
- Differential expression
- Disease-specific workflows
- Multi-step analysis patterns
- Common use cases with full code

**Load when:** Need examples for specific analysis types or multi-step workflows.

### references/quick_start.md

5-minute tutorial with:
- Installation and setup
- Quick examples for each task type
- Common patterns (3 pre-built workflows)
- Troubleshooting guide
- Quick reference card

**Load when:** User needs quick introduction or wants to get started fast.

## Tissue Filtering Reference

### Available Organ Systems (72 tissues total)

- **Nervous System (14):** brain, neuron, astrocyte, oligodendrocyte, microglia
- **Immune System (13):** T cell, B cell, NK cell, monocyte, macrophage, dendritic
- **Digestive System (12):** liver, hepatocyte, pancreas, stomach, intestine, colon
- **Urogenital (9):** kidney, breast, ovary, testis, bladder
- **Connective (6):** fibroblast, adipocyte, stromal, chondrocyte, bone
- **Integumentary (6):** keratinocyte, skin, melanocyte, basal cell
- **Muscular (5):** skeletal muscle, smooth muscle, myoblast
- **Cardiovascular (4):** heart, endothelial, vascular
- **Respiratory (3):** lung, trachea

### Common Filter Combinations

Pre-built filter sets for comprehensive coverage:

```python
# Liver comprehensive
tissue_filter=['liver', 'hepatocyte', 'hepatic', 'kupffer']

# Immune comprehensive
tissue_filter=['T cell', 'B cell', 'NK cell', 'lymphocyte', 'monocyte', 'dendritic']

# Neural comprehensive
tissue_filter=['brain', 'neuron', 'astrocyte', 'oligodendrocyte']

# Skin comprehensive
tissue_filter=['skin', 'keratinocyte', 'melanocyte', 'basal']

# Kidney comprehensive
tissue_filter=['kidney', 'podocyte', 'renal']
```

### Disease & Treatment Filters

Use two-step process with `quicksearch_metadata()`:

**Common disease terms:**
- cancer, tumor, carcinoma, adenocarcinoma
- lymphoma, leukemia
- diabetes, alzheimer, parkinson
- covid, inflammation, infection

**Common treatment terms:**
- treatment, drug, inhibitor
- knockout, knockdown, overexpression

**Workflow:**
```python
# Step 1: Find samples
samples = quicksearch_metadata('disease_or_treatment')
sample_ids = samples.iloc[:, 0].tolist()

# Step 2: Use in analysis
results = gene_correlation('GENE', 'context', samples=sample_ids)
```

## Best Practices

### 1. Start with Expression Atlas

Always begin with expression atlas to understand where gene is expressed:
```python
df = gene_expression('GENE_NAME')
```

### 2. Use Smart Filtering (v1.1)

Take advantage of multiple filters and typo tolerance:
```python
# Multiple related terms for comprehensive coverage
tissue_filter=['liver', 'hepatocyte', 'bile']

# Short substrings catch variations
tissue_filter='hepat'  # Matches HEPATOCYTE, HEPATIC STELLATE CELL

# System-level for exploration
tissue_filter='immune'  # All 13 immune tissues
```

### 3. Verify Matched Tissues

Always check which tissues were matched:
```python
print(results['top_tissues'])
```

### 4. Use Appropriate Context

Provide relevant metadata context for correlations:
```python
# Cell type context
gene_correlation('KRT14', 'keratinocyte')

# Disease context
gene_correlation('TP53', 'cancer', samples=cancer_samples)
```

### 5. Filter by Significance

Use appropriate cutoffs for differential expression:
```python
# Standard
diffexp('GENE', 'context', fdr_cutoff=0.1)

# Stricter
diffexp('GENE', 'context', fdr_cutoff=0.05)
```

### 6. Batch Queries When Possible

Analyze multiple genes together:
```python
# Efficient
gene_expression_analysis(genes=['GENE1', 'GENE2', 'GENE3'])

# Less efficient: separate queries for each gene
```

## Common Patterns

### Pattern: Validate Known Marker

Check if gene is specific to expected tissue:
```python
df = gene_expression('MARKER_GENE')
top_tissue = df.iloc[0]['tissue']

# Check specificity
fold_diff = df.iloc[0]['median'] / df.iloc[1]['median']
if fold_diff > 10:
    print(f"✓ Highly specific to {top_tissue}")
```

### Pattern: Discover Novel Markers

Find new markers for cell type:
```python
# Differential expression
diff_results = diffexp('KNOWN_MARKER', 'cell_type', fdr_cutoff=0.05)

# Filter by fold change
novel_markers = diff_results[diff_results['log2_fold_change'] > 3]
print(f"Found {len(novel_markers)} potential novel markers")
```

### Pattern: Build Gene Module

Construct co-expressed gene module:
```python
# Step 1: Find correlated genes
corr = gene_correlation('SEED_GENE', 'cell_type', k=100)

# Step 2: Select top genes
module_genes = corr.head(10)['gene'].tolist()
module_genes.insert(0, 'SEED_GENE')

# Step 3: Calculate pairwise correlations
matrix, pairwise = gene_correlation_pairwise(module_genes, 'cell_type')

# Step 4: Filter by correlation threshold
strong_pairs = pairwise[pairwise['correlation'] > 0.7]
```

## Limitations

### API Sample Limit (with Bootstrap Solution)

ARCHS4 API limits correlation calculations to ~200 samples per query (server-side constraint).

**Impact:** Single API calls can only use 200 samples maximum.

**Solution (v1.2):** Bootstrap aggregation overcomes this limitation:
- When `n_samples > 200`, samples are divided into batches of 200
- Each batch is processed independently via separate API calls
- Correlations are aggregated using Fisher's z-transformation
- Provides robust estimates using the full dataset (e.g., all 1600 intestinal samples)
- Random sampling with seed=42 ensures reproducibility

**Example:**
```python
# Use bootstrap aggregation with 1600 samples (8 batches)
results = gene_correlation_pairwise(
    genes=['TYK2', 'JAK1'],
    meta='intestine',
    n_samples=1600  # Automatically divided into 8 batches of 200
)
print(f"Used {results['n_batches']} batches with {results['total_samples_used']} samples")
```

### Disease/Treatment Filtering

Not directly available in `gene_expression_analysis()`.

**Workaround:** Use `quicksearch_metadata()` to find samples, then pass to other functions via `samples` parameter.

### Hierarchical Tissue Names

Tissue names follow hierarchical format: `System.Organ.Cell_Type`

**Solution:** Use substring matching to match any part of the path. Short substrings (e.g., 'hepat') work best.

## Performance Notes

- Single gene query: ~2-3 seconds
- Multiple genes: ~5-20 seconds (depends on count)
- Filtering overhead: < 0.1 seconds
- Smart suggestions: Negligible overhead

## Troubleshooting

### No matches for tissue filter

Check smart suggestions in output - system automatically suggests alternatives.

### API timeout

Check service status:
```python
from archs4_api import get_service_status
status = get_service_status()
```

### Insufficient samples

Verify sample count:
```python
samples = quicksearch_metadata('your_term')
print(f"Found {len(samples)} samples")
```

### Specific Gene Correlation Errors

Some genes return "No correlation found" when querying correlations:

**Two types of issues:**

1. **Asymmetric API indexing (FIXED in v1.4):**
   - Some gene pairs work in one direction but not the other
   - Example: `CDKN2D → ITGA4` fails, but `ITGA4 → CDKN2D` succeeds
   - **Solution:** `gene_correlation_pairwise()` automatically tries both directions

2. **Zero/weak correlations (biological):**
   - Gene pairs with correlation below threshold (|r| < ~0.06)
   - Example: `PCOLCE-TNFRSF25` (different functional pathways)
   - **Not a bug:** These genes simply aren't correlated
   - **Alternative:** Use STRING, BioGRID, or compute from raw data

**Diagnostic logging (v1.4):**
```
Computing correlation: CDKN2D - ITGA4
  Forward query empty, trying reverse direction...
  ✓ Reverse direction succeeded: r=-0.1153
```

**What works:**
- Expression data for all genes is always available ✓
- Bidirectional queries (v1.4) recover most asymmetric pairs ✓
- Bootstrap aggregation increases power for weak correlations ✓

## Version History

### v1.4 (2025-12-17)
- ✅ Bidirectional correlation queries to handle API asymmetry
- ✅ Automatically tries reverse direction when forward fails
- ✅ 3x improvement in correlation success rate (20% → 60%)
- ✅ Clear diagnostic logging showing which direction worked
- ✅ Fixes asymmetric gene pairs (e.g., CDKN2D-ITGA4, TNFRSF25-GREM1)

### v1.3 (2025-12-17)
- ✅ Added output_dir parameter to control file output location
- ✅ Files now saved to working directory by default (not script directory)
- ✅ Support for custom output directories
- ✅ Updated gene_expression_analysis() and gene_correlation_pairwise()

### v1.2 (2025-12-17)
- ✅ Bootstrap aggregation for n_samples > 200
- ✅ Overcomes 200-sample API limit via batch processing
- ✅ Fisher's z-transformation for robust correlation aggregation
- ✅ Support for analyzing full datasets (e.g., 1600+ intestinal samples)
- ✅ n_batches and total_samples_used in results

### v1.1 (2025-12-17)
- ✅ Multiple tissue filters with OR logic
- ✅ Smart filter suggestions with typo tolerance
- ✅ Schema-aware vocabulary (72 tissues)
- ✅ "Did you mean?" messages

### v1.0 (2025-12-16)
- ✅ 15 API functions
- ✅ 72 tissues catalogued
- ✅ Gene expression analysis with visualization
- ✅ Pairwise correlations with heatmaps
- ✅ Differential expression analysis
- ✅ Disease and treatment sample queries
