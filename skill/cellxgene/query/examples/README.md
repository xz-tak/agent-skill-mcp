# CELLxGENE Census Query Examples

This directory contains complete example workflows demonstrating advanced use cases of the CELLxGENE Census query tools.

## Examples

### IBD Intestinal Tissue Gene Coexpression Analysis

**File:** `ibd_coexpression_comprehensive.py`

**Description:** Comprehensive workflow for analyzing gene coexpression patterns in inflammatory bowel disease (IBD) intestinal tissues across multiple cell types.

**Features:**
- Multi-cell type querying with regex patterns (fibroblast|myo|smooth muscle|pericyte)
- Analysis of 5 gene lists representing key IBD pathways
- Automated correlation analysis with statistical significance testing
- Interactive HTML heatmaps with hierarchical clustering
- Interpretive biological reporting with clinical context
- Handles ~190,000 cells across 3 cell populations

**Key Components:**

1. **Data Acquisition**
   - Queries CELLxGENE Census with tissue and disease filters
   - Uses regex patterns to capture related cell types
   - Filters for specific genes to reduce data transfer

2. **Statistical Analysis**
   - Pearson correlation coefficients for gene pairs
   - P-value calculation with significance testing
   - Expression frequency analysis

3. **Visualization**
   - Interactive Plotly heatmaps with hover details
   - Hierarchical clustering for pattern discovery
   - Color scale centered at zero for correlation interpretation
   - P-value significance annotations (* p<0.05, ** p<0.01, *** p<0.001)

4. **Reporting**
   - Automated interpretive markdown report generation
   - Biological context for each gene list
   - Cross-cell type comparisons
   - Key findings summary with clinical implications

**Usage:**

```python
# The script is designed to be run directly
python ibd_coexpression_comprehensive.py

# Key configuration at top of file:
TISSUES = ["intestine", "colon", "rectum", "ileum", "sigmoid"]
DISEASES = ["crohn's disease", "ulcerative colitis", "inflammatory bowel disease"]
CELL_TYPES = {
    "fibroblast": "fibroblast|myo|smooth muscle|pericyte",
    "immune": "T cell|B cell|plasma cell|macrophage|monocyte|dendritic",
    "endothelial": "endothelial"
}
GENE_LISTS = {
    "list1": ["TYK2", "JAK1"],           # JAK-STAT signaling
    "list2": ["TNFRSF25", "GREM1"],      # TNFR/BMP
    "list3": ["TNFRSF25", "PCOLCE"],     # TNFR/ECM
    "list4": ["CDKN2D", "ITGA4", "ITGB1"], # Cell cycle/Integrin
    "list5": ["CDKN2D", "PCOLCE"]        # Cell cycle/ECM
}
```

**Outputs:**

```
ibd_coexpression_comprehensive_results/
├── fibroblast/
│   ├── fibroblast_data.h5ad (AnnData with 12,612 cells)
│   ├── fibroblast_results.json (metadata)
│   ├── fibroblast_list1_correlation.csv
│   ├── fibroblast_list1_pvalues.csv
│   ├── fibroblast_list1_heatmap.html (interactive)
│   └── ... (for all 5 gene lists)
├── immune/
│   └── ... (same structure, 170,069 cells)
├── endothelial/
│   └── ... (same structure, 6,747 cells)
└── comprehensive_report.md (interpretive summary)
```

**Adaptation for Your Data:**

To adapt this workflow for your research:

1. **Modify tissue filters:**
   ```python
   TISSUES = ["liver", "kidney"]  # Change to your tissues of interest
   ```

2. **Update disease conditions:**
   ```python
   DISEASES = ["diabetes mellitus", "normal"]  # Your conditions
   ```

3. **Customize cell types:**
   ```python
   CELL_TYPES = {
       "hepatocyte": "hepatocyte",
       "kupffer": "kupffer|macrophage"  # Use regex for related types
   }
   ```

4. **Define your gene lists:**
   ```python
   GENE_LISTS = {
       "list1": ["GENE1", "GENE2"],
       "list2": ["GENE3", "GENE4", "GENE5"]
   }
   ```

5. **Customize interpretations:**
   - Edit the `generate_interpretive_report()` function
   - Add biological context specific to your pathways
   - Modify correlation strength thresholds if needed

**Key Functions:**

- `query_cellxgene_data()`: Fetches data from Census with filters
- `compute_correlation_matrix()`: Calculates pairwise correlations
- `plot_correlation_heatmap()`: Creates interactive visualizations
- `analyze_cell_type()`: Orchestrates analysis for one cell type
- `generate_interpretive_report()`: Creates comprehensive markdown report

**Dependencies:**

```bash
pip install cellxgene-census scanpy plotly pandas numpy scipy matplotlib seaborn
```

**Performance Notes:**

- Queries ~190,000 cells in ~2-3 minutes
- Memory usage scales with cell count
- Consider filtering by specific genes to reduce data transfer
- Interactive heatmaps work best with <20 genes per visualization

**Citation:**

If you use this workflow in your research, please cite:
- CELLxGENE Census: https://chanzuckerberg.github.io/cellxgene-census/
- Your analysis should acknowledge the data sources from Census

## Creating Your Own Workflows

This example demonstrates the recommended pattern for complex analyses:

1. **Configuration at top** - Makes it easy to modify parameters
2. **Modular functions** - Each function has a single responsibility
3. **Progress logging** - Informative print statements for tracking
4. **Error handling** - Graceful handling of missing data
5. **Multiple output formats** - CSV for data, HTML for visualization, MD for reports
6. **Interpretive context** - Don't just report numbers, explain biological significance

## Additional Resources

- **Query documentation:** `../SKILL.md`
- **Coexpression analysis:** `../scripts/analyze_coexpression.py`
- **Metadata inspection:** `../scripts/inspect_metadata_fields.py`
- **Census documentation:** https://chanzuckerberg.github.io/cellxgene-census/
