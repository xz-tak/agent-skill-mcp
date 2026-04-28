# Output Schema

## Directory Structure

```
<output_dir>/
├── pseudobulk_<input_name>_DE.tsv       # Main DEG results
├── pseudobulk_aggregation_qc.tsv        # Aggregation QC table
├── combined_markers_clusters.tsv         # Cluster markers
├── fgsea_custom_signatures.tsv           # Custom signature results (if provided)
│
├── fgsea_groups_clusters/                # Pathway enrichment per comparison
│   ├── condition_clusters_reactome_cluster.tsv
│   ├── condition_clusters_msigdb_cluster.tsv
│   ├── condition_clusters_custom_cluster.tsv  # If custom signatures
│   └── *_pathways.pdf                    # Pathway enrichment plots
│
├── fgsea_clusters/                       # Cluster-level pathway analysis
│   ├── clusters_reactome_cluster.tsv
│   ├── clusters_reactome_pathways.pdf
│   ├── clusters_msigdbC8_clusters.tsv
│   ├── clusters_msigdbC8_pathways.pdf
│   ├── clusters_custom_clusters.tsv      # If custom signatures
│   └── clusters_custom_pathways.pdf
│
├── speckle_diffprop/                     # Cell composition results
│   └── comb_condition_<group2>_vs_<group1>.tsv
│
├── volcano_plots/                        # Volcano plots per cluster/comparison
│   └── Volcano_plot_for_DEG_<group1>_vs_<group2>_<cluster>.png
│
└── report_data/                          # Filtered data for interpretation
    ├── <comparison>_top_deg.tsv
    ├── <comparison>_top_pathways_*.tsv
    ├── <comparison>_speckle.tsv
    ├── <comparison>_top_custom_signatures.tsv  # If custom signatures
    ├── CLAUDE_CODE_INSTRUCTIONS.md
    └── report_summary.txt
```

---

## File Formats

### pseudobulk_*_DE.tsv

Main differential expression results across all comparisons and clusters.

| Column | Description |
|--------|-------------|
| `gene_name` | Gene symbol |
| `logFC` | Log2 fold change (group2 vs group1) |
| `PValue` | Raw p-value |
| `FDR` | BH-adjusted p-value |
| `group1` | Reference group |
| `group2` | Comparison group |
| `cluster` | Cell type/cluster |
| `comparison` | Comparison string (group2vsgroup1) |
| `clust_res` | Cluster resolution level |

**Interpretation:**
- Positive logFC: Gene upregulated in group2 relative to group1
- Negative logFC: Gene downregulated in group2 relative to group1
- FDR < 0.05: Statistically significant

---

### fgsea_*_cluster.tsv

Pathway enrichment results from fGSEA.

| Column | Description |
|--------|-------------|
| `pathway` | Pathway name |
| `pval` | Raw p-value |
| `padj` | BH-adjusted p-value |
| `NES` | Normalized Enrichment Score |
| `size` | Number of genes in pathway |
| `leadingEdge` | Comma-separated leading edge genes |
| `cluster` | Cell type/cluster |
| `database` | Source database (reactome/msigdb/custom) |
| `comparison` | Comparison string |

**Interpretation:**
- Positive NES: Pathway enriched in upregulated genes
- Negative NES: Pathway enriched in downregulated genes
- Leading edge: Core genes driving the enrichment

---

### combined_markers_clusters.tsv

Cluster marker genes identified by scran's findMarkers.

| Column | Description |
|--------|-------------|
| `feature` | Gene symbol |
| `cluster` | Cell type/cluster |
| `p.adjusted` | Bonferroni-adjusted p-value |
| `p.value` | Raw p-value |
| `summary.logFC` | Summary log fold change |
| `Top` | Rank within cluster |

---

### speckle_diffprop/*.tsv

Cell composition analysis results from speckle.

| Column | Description |
|--------|-------------|
| `BaselineProp.group1` | Proportion in reference group |
| `BaselineProp.group2` | Proportion in comparison group |
| `PropRatio` | Ratio of proportions |
| `PropMean` | Mean proportion |
| `Fstatistic` | F-statistic |
| `P.Value` | Raw p-value |
| `FDR` | BH-adjusted p-value |

**Row names:** Cell type names

**Interpretation:**
- PropRatio > 1: Cell type proportion higher in group2
- PropRatio < 1: Cell type proportion lower in group2
- FDR < 0.05: Significant composition change

---

### report_data/*_top_deg.tsv

Filtered top DEG results for Claude Code interpretation.

Same columns as main DEG file, but filtered to:
- Top k genes by |logFC| within FDR cutoff
- Split equally between up/down regulated

---

### report_data/*_top_pathways_*.tsv

Filtered top pathway results for Claude Code interpretation.

Same columns as fgsea files, but filtered to:
- Top k pathways by |NES| within padj cutoff
- Split equally between positive/negative NES

---

### report_data/*_speckle.tsv

Filtered cell composition results.

Same format as speckle output, filtered to FDR significant results.

---

### report_data/CODEX_INSTRUCTIONS.md

Instructions for Codex XHIGH interpretation.

Contains:
- List of files to process
- Analysis guidelines for each data type
- Filter parameters used
- Questions to address in interpretation

---

## Custom Signatures Output

When custom signatures are provided (`--custom_signatures`):

### fgsea_custom_signatures.tsv
Root-level file with custom signature enrichment results.

### fgsea_groups_clusters/condition_clusters_custom_cluster.tsv
Comparison-level custom signature enrichment.

### fgsea_clusters/clusters_custom_clusters.tsv
Cluster-level custom signature enrichment.

### report_data/<comparison>_top_custom_signatures.tsv
Filtered custom signatures for interpretation.

---

## QC Files

### pseudobulk_aggregation_qc.tsv

Matrix showing number of pseudobulk samples per condition × cell type.

- Rows: Conditions
- Columns: Cell types
- Values: Number of aggregated samples

Use this to identify:
- Cell types missing in certain conditions
- Unbalanced designs
- Potential confounding
