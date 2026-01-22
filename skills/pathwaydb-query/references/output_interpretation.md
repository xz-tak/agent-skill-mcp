# Output Formats and Interpretation Guide

This reference explains how to interpret and present pathway query results to users. Use this when processing query outputs.

## Single Gene Query Outputs

### pathway_query.py Output Formats

#### Text Format (Default)
```
================================================================================
Pathway Database Query Results for: TP53
================================================================================

Summary:
  Total databases queried: 3
  Databases with results: 3
  Total pathways/gene sets found: 68
    ✓ KEGG: 8 results
    ✓ MSIGDB: 15 results
    ✓ REACTOME: 45 results

--------------------------------------------------------------------------------
KEGG Results
--------------------------------------------------------------------------------
Organism: hsa
Total pathways: 8

Gene: hsa:7157
  path:hsa04010: MAPK signaling pathway - Homo sapiens (human)
  path:hsa04115: p53 signaling pathway - Homo sapiens (human)
  path:hsa04210: Apoptosis - Homo sapiens (human)
  ...
```

**How to present:** Summarize key statistics, highlight most relevant pathways

#### JSON Format (--output json)
```json
{
  "gene_symbol": "TP53",
  "databases": {
    "kegg": {
      "query": {"symbol": "TP53", "organism": "hsa"},
      "genes_found": [...],
      "total_pathways": 8
    },
    "msigdb": {
      "query": {"gene_symbol": "TP53", "collection": "H"},
      "gene_sets": [...],
      "total_gene_sets": 15
    },
    "reactome": {
      "query": {"gene_id": "TP53", "species": 9606},
      "terms": [...],
      "total_terms": 45
    }
  }
}
```

**How to use:** Parse with Python/R for programmatic analysis, extract specific fields

---

## Multi-Gene Query Outputs

The `multi_gene_analysis.py` script generates multiple file types. Always mention all generated files to the user.

### File Naming Convention
All files use the format: `{output_prefix}_{type}.{extension}`

Example: If output prefix is "cancer_genes":
- `cancer_genes_wide.csv`
- `cancer_genes_long.csv`
- `cancer_genes_summary.csv`
- `cancer_genes_wide.xlsx`
- `cancer_genes_kegg.png`
- etc.

---

### 1. Wide Format CSV (`{prefix}_wide.csv`)

**Structure:**
```csv
Gene,KEGG_Pathways,KEGG_Count,Reactome_Pathways,Reactome_Count,MSigDB_GeneSets,MSigDB_Count,Total_Count
TP53,path:hsa04115;path:hsa04210;...,68,Pyroptosis;Activation of NOXA;...,45,HALLMARK_P53_PATHWAY;...,4,117
BRCA1,path:hsa03440;path:hsa04120;...,12,HDR through HRR;...,27,HALLMARK_E2F_TARGETS;...,4,43
```

**Characteristics:**
- One row per gene
- Pathways in semicolon-separated lists
- Count columns for quick reference
- Total count column sums all databases

**Best for:**
- Quick overview
- Comparing total pathway counts per gene
- Sharing with non-technical users

**How to present:**
"I've created a summary table with all pathways. Each gene has X pathways total, with Y from KEGG, Z from Reactome, and W from MSigDB."

**Reading with pandas:**
```python
import pandas as pd
wide_df = pd.read_csv("results_wide.csv")

# Access pathways for a gene
tp53_kegg = wide_df[wide_df['Gene'] == 'TP53']['KEGG_Pathways'].values[0].split(';')
```

---

### 2. Long Format CSV (`{prefix}_long.csv`)

**Structure:**
```csv
Gene,Database,Pathway
TP53,KEGG,p53 signaling pathway - Homo sapiens (human)
TP53,KEGG,Apoptosis - Homo sapiens (human)
TP53,REACTOME,Pyroptosis
TP53,MSIGDB,HALLMARK_P53_PATHWAY
BRCA1,KEGG,Homologous recombination - Homo sapiens (human)
BRCA1,REACTOME,HDR through HRR
```

**Characteristics:**
- One row per gene-pathway pair
- Normalized format (tidy data)
- Easy to filter, group, and analyze

**Best for:**
- Downstream statistical analysis
- Filtering specific pathways
- Counting pathway occurrences
- Merging with other datasets

**How to present:**
"I've created a detailed table with one row per gene-pathway pair. This format is ideal for further analysis in R or Python."

**Reading with pandas:**
```python
import pandas as pd
long_df = pd.read_csv("results_long.csv")

# Find pathways shared by all genes
pathway_counts = long_df.groupby('Pathway')['Gene'].nunique()
shared_pathways = pathway_counts[pathway_counts == 3].index.tolist()

# Filter by database
kegg_only = long_df[long_df['Database'] == 'KEGG']

# Count pathways per gene
counts = long_df.groupby('Gene').size()
```

**Reading with R:**
```r
library(tidyverse)
long_df <- read_csv("results_long.csv")

# Find shared pathways
shared <- long_df %>%
  group_by(Pathway) %>%
  summarize(n_genes = n_distinct(Gene)) %>%
  filter(n_genes == 3)

# Pathways unique to TP53
tp53_unique <- long_df %>%
  group_by(Pathway) %>%
  filter(n_distinct(Gene) == 1, Gene == "TP53")
```

---

### 3. Summary CSV (`{prefix}_summary.csv`)

**Structure:**
```csv
Gene,Database,Pathway_Count
TP53,KEGG,68
TP53,REACTOME,45
TP53,MSIGDB,4
BRCA1,KEGG,12
BRCA1,REACTOME,27
BRCA1,MSIGDB,4
```

**Characteristics:**
- Pathway counts only (no pathway names)
- One row per gene-database combination
- Compact summary

**Best for:**
- Quick comparison of pathway enrichment
- Visualizing count differences
- Summary statistics

**How to present:**
"Here's a summary of pathway counts: TP53 has 68 KEGG pathways, 45 Reactome pathways, and 4 MSigDB gene sets."

**Visualization example:**
```python
import pandas as pd
import matplotlib.pyplot as plt

summary_df = pd.read_csv("results_summary.csv")
summary_pivot = summary_df.pivot(index='Gene', columns='Database', values='Pathway_Count')

summary_pivot.plot(kind='bar', figsize=(10, 6))
plt.title("Pathway Counts by Gene and Database")
plt.ylabel("Number of Pathways")
plt.show()
```

---

### 4. Excel File (`{prefix}_wide.xlsx`)

**Structure:**
Multi-sheet workbook with 4 sheets:
1. **Summary** - Wide format overview (same as wide CSV)
2. **KEGG** - All KEGG pathways (Gene | Pathway columns)
3. **REACTOME** - All Reactome pathways (Gene | Pathway columns)
4. **MSIGDB** - All MSigDB gene sets (Gene | Pathway columns)

**Characteristics:**
- All data in one file
- Separate sheets for each database
- Easy to share and browse

**Best for:**
- Sharing with collaborators
- Excel-based analysis
- Non-programmatic review

**How to present:**
"I've created an Excel file with multiple sheets. The Summary sheet has an overview, and each database has its own sheet with detailed results."

**Reading with pandas:**
```python
import pandas as pd

# Read specific sheet
summary_df = pd.read_excel("results_wide.xlsx", sheet_name="Summary")
kegg_df = pd.read_excel("results_wide.xlsx", sheet_name="KEGG")

# Read all sheets
all_sheets = pd.read_excel("results_wide.xlsx", sheet_name=None)
```

---

## UpSet Plot Outputs

### Generated Plot Files

For output prefix "results", the following PNG files are generated (300 DPI):
1. **`results_kegg.png`** - KEGG pathway overlap across genes
2. **`results_reactome.png`** - Reactome pathway overlap across genes
3. **`results_msigdb.png`** - MSigDB gene set overlap across genes
4. **`results_combined.png`** - All databases combined

### Understanding UpSet Plots

**Visual Structure:**
```
Gene Matrix (bottom):
          TP53  BRCA1  EGFR
Bar 1:     ●     ●      ●      Height: 25 pathways
Bar 2:     ●     ●              Height: 12 pathways
Bar 3:     ●            ●      Height: 8 pathways
Bar 4:     ●                   Height: 45 pathways
```

**Interpretation:**
- **Filled dots (●)** indicate which genes share pathways
- **Bar height** = number of pathways in that intersection
- **Bars sorted by size** (most pathways first)

**Key Insights:**

1. **Shared by all genes** (all dots filled):
   - Core pathways common to all queried genes
   - May indicate related biological functions

2. **Pairwise overlaps** (two dots filled):
   - Pathways shared by exactly two genes
   - Shows specific functional relationships

3. **Gene-specific** (one dot filled):
   - Pathways unique to one gene
   - Indicates specialized functions

### How to Present UpSet Plots

**Example presentation:**
```
I've generated UpSet plots showing pathway overlap:

KEGG Analysis:
- 25 pathways are shared by all 3 genes (TP53, BRCA1, EGFR)
- 12 pathways are shared by TP53 and BRCA1 only
- 45 pathways are unique to TP53
- See results_kegg.png for full visualization

The combined plot (results_combined.png) shows pathway overlap across all databases.
```

### Common Plot Parameters

**`--max-bars N`** (default: 20)
- Shows top N intersections only
- Reduce if plot is too cluttered
- Typical values: 10-25

**`--min-intersection N`** (default: 1)
- Only show intersections with ≥N pathways
- Increase to filter out small overlaps
- Typical values: 2-5

**Example:**
```bash
--max-bars 15 --min-intersection 3
```
Shows top 15 intersections with at least 3 pathways each

---

## Interpreting Results for Users

### High Pathway Count (>50)
"[GENE] is found in X pathways, indicating it's involved in many biological processes. This is common for tumor suppressor genes and signaling hubs."

### Low Pathway Count (<10)
"[GENE] appears in X pathways, suggesting more specialized or less well-annotated functions in these databases."

### High Overlap Between Genes
"These genes share X pathways, suggesting they may be functionally related or involved in similar biological processes."

### Low Overlap Between Genes
"These genes have limited pathway overlap (X shared pathways), indicating distinct biological functions despite [user's context]."

### Database-Specific Insights

**Many KEGG results:**
"KEGG found X pathways, indicating strong metabolic and signaling pathway associations."

**Many Reactome results:**
"Reactome found X pathways, indicating detailed molecular-level annotations are available."

**Many MSigDB results:**
"MSigDB collection [X] found Y gene sets, indicating this gene is well-represented in [cancer/immune/experimental] signatures."

---

## Presenting Results Based on User Request

### User wants "summary"
Show:
- Total pathway counts per database
- Top 3-5 most relevant pathways
- Mention file locations for details

### User wants "comprehensive results"
Show:
- Full pathway counts
- List all pathways (or top 20 if too many)
- Mention all generated files
- Offer to filter or analyze further

### User wants "comparison" (multi-gene)
Show:
- Pathway overlap statistics from UpSet analysis
- Shared vs unique pathways
- Highlight interesting intersections
- Mention plot files

### User wants "export"
List all generated files:
```
Generated files:
- results_wide.csv - Overview table with semicolon-separated pathways
- results_long.csv - Detailed table (one row per gene-pathway pair)
- results_summary.csv - Pathway counts by gene and database
- results_wide.xlsx - Multi-sheet Excel file
- results_kegg.png - KEGG pathway overlap visualization
- results_reactome.png - Reactome pathway overlap visualization
- results_msigdb.png - MSigDB gene set overlap visualization
- results_combined.png - Combined pathway overlap across all databases
```

---

## Common Analysis Patterns

### Finding Shared Pathways

**From long format:**
```python
import pandas as pd

long_df = pd.read_csv("results_long.csv")
n_genes = long_df['Gene'].nunique()

pathway_counts = long_df.groupby('Pathway')['Gene'].nunique()
shared_pathways = pathway_counts[pathway_counts == n_genes].index.tolist()

print(f"Pathways shared by all {n_genes} genes:")
for pathway in shared_pathways:
    print(f"  - {pathway}")
```

### Finding Gene-Specific Pathways

**From long format:**
```python
import pandas as pd

long_df = pd.read_csv("results_long.csv")

# Pathways unique to TP53
tp53_pathways = set(long_df[long_df['Gene'] == 'TP53']['Pathway'])
other_pathways = set(long_df[long_df['Gene'] != 'TP53']['Pathway'])
tp53_specific = tp53_pathways - other_pathways

print(f"Pathways unique to TP53: {len(tp53_specific)}")
```

### Comparing Pathway Enrichment

**From summary format:**
```python
import pandas as pd

summary_df = pd.read_csv("results_summary.csv")

# Pivot to compare
pivot = summary_df.pivot(index='Gene', columns='Database', values='Pathway_Count')
print(pivot)

# Find gene with most KEGG pathways
max_kegg = pivot['KEGG'].idxmax()
print(f"Gene with most KEGG pathways: {max_kegg}")
```

---

## Error Messages and Handling

### "No results found for [GENE] in [DATABASE]"
**Interpretation:** Gene not annotated in that database
**User message:** "[GENE] was not found in [DATABASE]. This could mean the gene symbol is not recognized, or the gene is not annotated in this database. You might try alternative gene names or a different database."

### "Found 0 pathways"
**Interpretation:** Query succeeded but gene has no pathway annotations
**User message:** "The query completed successfully, but [GENE] doesn't appear in any pathways in [DATABASE]. This is unusual for well-studied genes and may indicate the gene symbol needs verification."

### Empty UpSet plot or very few bars
**Interpretation:** Genes have little to no pathway overlap
**User message:** "The genes have very little pathway overlap. This suggests they have distinct biological functions. See individual pathway lists in the CSV files for details."

---

## File Size Expectations

### Small Query (1 gene, 1 database)
- CSV files: ~1-10 KB
- No plots

### Medium Query (3 genes, all databases, H+C2)
- CSV files: ~50-200 KB each
- Excel file: ~100-300 KB
- PNG plots: ~50-100 KB each

### Large Query (10 genes, all databases, many collections)
- CSV files: ~500 KB - 2 MB each
- Excel file: ~1-5 MB
- PNG plots: ~100-200 KB each

---

## Best Practices for Presenting Results

1. **Start with summary statistics** - give user the big picture first
2. **Highlight interesting findings** - shared pathways, unusually high/low counts
3. **Always mention file locations** - tell user where to find detailed results
4. **Offer follow-up analysis** - "Would you like me to filter for specific pathways?" or "Should I analyze the overlap in more detail?"
5. **Explain plot files** - briefly describe what each UpSet plot shows
6. **Suggest next steps** - "You can open the Excel file to browse pathways by database" or "The long CSV format is best for statistical analysis in R or Python"

---

## Quick Reference: Which File to Use?

| User Need | Recommended File | Why |
|-----------|------------------|-----|
| Quick overview | `_wide.csv` or `_summary.csv` | Easy to scan |
| Detailed analysis | `_long.csv` | Best for filtering, grouping |
| Sharing with collaborators | `_wide.xlsx` | Multi-sheet, no code needed |
| Pathway overlap visualization | `_kegg.png`, `_reactome.png`, `_msigdb.png` | Visual intersection |
| R/Python analysis | `_long.csv` | Tidy format |
| Counting pathways | `_summary.csv` | Quick stats |
| Finding shared pathways | `_long.csv` or plots | Group by pathway |
| Finding gene-specific pathways | `_long.csv` | Filter by gene |
