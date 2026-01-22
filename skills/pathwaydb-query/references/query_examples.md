# Query Examples and Common Workflows

This reference provides practical examples of common pathway query workflows and patterns. Use these examples as templates for similar user requests.

## Basic Query Patterns

### Pattern 1: Single Gene, All Databases

**User Request Examples:**
- "Find all pathways for TP53"
- "What pathways is BRCA1 involved in?"
- "Query pathway databases for EGFR"

**Implementation:**
```bash
conda run -n claude_test python scripts/pathway_query.py TP53
```

**When to use:** User wants comprehensive pathway information for one gene

---

### Pattern 2: Single Gene, Specific Database

**User Request Examples:**
- "Get KEGG pathways for TP53"
- "What Reactome pathways contain BRCA1?"
- "Find MSigDB Hallmark gene sets with EGFR"

**Implementation:**
```bash
# KEGG only
conda run -n claude_test python scripts/kegg_api.py TP53

# Reactome only
conda run -n claude_test python scripts/reactome_api.py TP53

# MSigDB only (Hallmark)
conda run -n claude_test python scripts/msigdb_api.py TP53 --collection H
```

**When to use:** User specifies a particular database or you need results from just one source

---

### Pattern 3: Single Gene, Custom Parameters

**User Request Examples:**
- "Find mouse pathways for Trp53"
- "Get MSigDB C2 curated pathways for TP53"
- "Query Reactome reactions (not pathways) for BRCA1"

**Implementation:**
```bash
# Mouse KEGG pathways
conda run -n claude_test python scripts/kegg_api.py Trp53 mmu

# MSigDB C2 collection
conda run -n claude_test python scripts/msigdb_api.py TP53 --collection C2

# Reactome reactions instead of pathways
conda run -n claude_test python scripts/reactome_api.py BRCA1 --map-to reactions
```

**When to use:** User specifies organism, species, collection, or other parameters

---

### Pattern 4: Multiple Genes with Visualization

**User Request Examples:**
- "Compare pathways for TP53, BRCA1, and EGFR"
- "Find shared pathways between these cancer genes"
- "Analyze pathway overlap for TP53, BRCA1, EGFR with visualization"

**Implementation:**
```bash
conda run -n claude_test python scripts/multi_gene_analysis.py TP53 BRCA1 EGFR \
  --output comparison \
  --msigdb-collections H C2 C6
```

**Outputs:**
- `comparison_wide.csv` - Overview table
- `comparison_long.csv` - Detailed gene-pathway pairs
- `comparison_summary.csv` - Pathway counts
- `comparison_wide.xlsx` - Multi-sheet Excel file
- `comparison_kegg.png` - KEGG UpSet plot
- `comparison_reactome.png` - Reactome UpSet plot
- `comparison_msigdb.png` - MSigDB UpSet plot
- `comparison_combined.png` - Combined UpSet plot

**When to use:** User has multiple genes and wants to see overlap/intersection

---

### Pattern 5: Multiple Genes, Tables Only (No Plots)

**User Request Examples:**
- "Get pathways for 10 genes and export to Excel"
- "Query pathways for these genes without visualization"
- "I need a table of pathways for TP53, BRCA1, EGFR, KRAS, MYC"

**Implementation:**
```bash
conda run -n claude_test python scripts/multi_gene_analysis.py \
  TP53 BRCA1 EGFR KRAS MYC \
  --output gene_panel \
  --no-plot
```

**When to use:**
- User has many genes (>5)
- User only needs tabular data
- Visualization is not requested

---

## Advanced Workflows

### Workflow 1: Comprehensive Cancer Gene Panel Analysis

**Scenario:** User wants to analyze a panel of cancer-related genes across multiple MSigDB collections with visualization

**Command:**
```bash
conda run -n claude_test python scripts/multi_gene_analysis.py \
  TP53 BRCA1 BRCA2 EGFR KRAS MYC PTEN ATM \
  --output cancer_panel_comprehensive \
  --msigdb-collections H C2 C6 \
  --max-bars 20 \
  --min-intersection 2
```

**What it does:**
1. Queries KEGG, Reactome, and MSigDB (H, C2, C6) for 8 genes
2. Creates UpSet plots showing pathway overlap (top 20 intersections, minimum 2 pathways)
3. Exports CSV/Excel tables with all results

**Use when:** Comprehensive cancer pathway analysis is needed

---

### Workflow 2: Cross-Species Pathway Comparison

**Scenario:** User wants to compare pathways between human and mouse orthologs

**Commands:**
```bash
# Query human gene
conda run -n claude_test python scripts/pathway_query.py TP53 \
  --kegg-organism hsa \
  --reactome-species 9606 \
  --export human_tp53.json

# Query mouse ortholog
conda run -n claude_test python scripts/pathway_query.py Trp53 \
  --kegg-organism mmu \
  --reactome-species 10090 \
  --export mouse_trp53.json
```

**Follow-up:** Compare the JSON outputs programmatically or manually review

**Use when:** User needs to compare pathways across species

---

### Workflow 3: MSigDB Collection Exploration

**Scenario:** User wants to explore which MSigDB collections contain their gene

**Commands:**
```bash
# Query Hallmark
conda run -n claude_test python scripts/msigdb_api.py TP53 --collection H

# Query Curated Pathways
conda run -n claude_test python scripts/msigdb_api.py TP53 --collection C2

# Query Oncogenic Signatures
conda run -n claude_test python scripts/msigdb_api.py TP53 --collection C6

# Query ALL collections at once
conda run -n claude_test python scripts/msigdb_api.py TP53 --all
```

**Use when:** User wants to understand which MSigDB collections are relevant for their gene

---

### Workflow 4: Detailed Reactome Reaction Analysis

**Scenario:** User wants molecular-level reaction details from Reactome (not just pathways)

**Command:**
```bash
conda run -n claude_test python scripts/reactome_api.py TP53 \
  --map-to reactions \
  --output json \
  > tp53_reactions.json
```

**Use when:** User needs detailed molecular mechanisms or reaction-level information

---

### Workflow 5: Batch Gene Query with Programmatic Analysis

**Scenario:** User has a list of genes and wants to query all of them

**Python Script Example:**
```python
from pathway_query import query_all_databases
import pandas as pd

genes = ["TP53", "BRCA1", "EGFR", "KRAS", "MYC"]
results = []

for gene in genes:
    print(f"Querying {gene}...")
    result = query_all_databases(gene, parallel=True)

    # Extract counts
    results.append({
        "Gene": gene,
        "KEGG_Pathways": result["databases"]["kegg"]["total_pathways"],
        "Reactome_Pathways": result["databases"]["reactome"]["total_terms"],
        "MSigDB_GeneSets": result["databases"]["msigdb"]["total_gene_sets"]
    })

# Create summary table
df = pd.DataFrame(results)
df.to_csv("batch_query_summary.csv", index=False)
print(df)
```

**When to use:** User has many genes (10+) and needs programmatic control

---

## Common User Request Patterns

### "Find pathways for [GENE]"
→ Use Pattern 1 (Single Gene, All Databases)

### "Compare pathways between [GENE1], [GENE2], [GENE3]"
→ Use Pattern 4 (Multiple Genes with Visualization)

### "Get [DATABASE] pathways for [GENE]"
→ Use Pattern 2 (Single Gene, Specific Database)

### "Analyze [GENE] in [ORGANISM/SPECIES]"
→ Use Pattern 3 (Single Gene, Custom Parameters) with organism/species codes

### "Find shared pathways for these genes: [LIST]"
→ Use Pattern 4 and highlight shared pathways in results

### "Query MSigDB [COLLECTION] for [GENE]"
→ Use Pattern 2 or 3 with specific MSigDB collection

### "Get pathway data for [GENES] and export to Excel"
→ Use Pattern 5 or execute query and show user the Excel file location

### "Visualize pathway overlap for [GENES]"
→ Use Pattern 4, emphasize the UpSet plots generated

---

## Output Handling Patterns

### When User Wants Summary
Show key statistics:
```
Found X pathways total for [GENE]:
- KEGG: X pathways
- Reactome: X pathways
- MSigDB: X gene sets
```

### When User Wants Full Results
Read and display the generated files:
- Use `pandas` to read CSV files
- Show top N results or filtered results
- Mention file locations for download

### When User Wants Visualization
Always mention:
- Where the PNG files are saved
- What each UpSet plot shows
- How to interpret the plots (shared vs unique pathways)

### When User Wants Export
Mention all generated files:
- CSV files (wide, long, summary)
- Excel file (multi-sheet)
- PNG plots (if applicable)
- JSON exports (if requested)

---

## Parameter Selection Guidelines

### Selecting KEGG Organism
- Default to `hsa` (human) unless specified
- Common organisms: hsa (human), mmu (mouse), rno (rat)
- Reference: `database_guide.md` for full list

### Selecting Reactome Species
- Default to `9606` (human) unless specified
- Use taxonomy IDs or common names
- Reference: `database_guide.md` for full list

### Selecting MSigDB Collections
**Default:** H (Hallmark) - quick overview
**Comprehensive:** H, C2, C6 - pathways + oncogenic signatures
**Immune-focused:** H, C7 - immunologic signatures
**All collections:** Use `--all` flag (slow but complete)

### Setting UpSet Plot Parameters
- `--max-bars`: Default 20, reduce to 10-15 if plot is cluttered
- `--min-intersection`: Default 1, increase to 2-3 to show only significant overlaps
- Use both together for cleaner plots with many genes

---

## Error Handling Patterns

### Gene Not Found in Database
**Response:**
"I searched [DATABASE] but didn't find results for [GENE]. This could mean:
1. Gene symbol spelling may differ (try alternative names)
2. Gene may not be annotated in this database
3. For KEGG: try different organism code
4. For Reactome: try different identifier type (UniProt/NCBI/ENSEMBL)"

### Slow Query Warning
**When:** User queries MSigDB with `--all` or many collections
**Response:**
"Querying all MSigDB collections may take 2-3 minutes as it downloads gene set files. I'll proceed, but you can also query specific collections (H, C2, C6) for faster results."

### Multiple Genes Without Clear Intent
**When:** User provides multiple genes but unclear if they want comparison or separate queries
**Ask:**
"Would you like me to:
1. Compare pathways across these genes (with UpSet visualization)
2. Query each gene separately and provide individual results
3. Generate a summary table only"

---

## Integration with Other Analysis

### Exporting for R
Recommend using the `_long.csv` format:
```bash
# Long format is best for tidyverse
conda run -n claude_test python scripts/multi_gene_analysis.py TP53 BRCA1 EGFR \
  --output for_R
```

Then in R:
```r
library(tidyverse)
df <- read_csv("for_R_long.csv")
```

### Exporting for Python/Pandas
Use `_long.csv` or `_wide.csv`:
```python
import pandas as pd
long_df = pd.read_csv("results_long.csv")
wide_df = pd.read_csv("results_wide.csv")
```

### Exporting for GSEA
MSigDB queries are directly compatible with GSEA software. Suggest:
1. Export gene lists from query
2. Use MSigDB collections in GSEA
3. Reference the gene set names from query results

---

## Performance Tips

1. **Single gene query:** ~5-10 seconds across all databases
2. **Multi-gene query (3 genes, H collection):** ~30 seconds
3. **Multi-gene query (5 genes, H+C2+C6):** ~60-90 seconds
4. **First MSigDB query per collection:** Slow (downloads GMT file)
5. **Subsequent MSigDB queries:** Fast (uses cache)

**Recommendation:** For large gene lists (>10 genes), use `--no-plot` and focus on table exports.

---

## Quick Reference Command Templates

```bash
# Single gene, all databases
conda run -n claude_test python scripts/pathway_query.py [GENE]

# Single gene, specific database
conda run -n claude_test python scripts/kegg_api.py [GENE]
conda run -n claude_test python scripts/reactome_api.py [GENE]
conda run -n claude_test python scripts/msigdb_api.py [GENE] --collection [H|C2|C6]

# Multiple genes with viz
conda run -n claude_test python scripts/multi_gene_analysis.py [GENE1] [GENE2] [GENE3] \
  --output [PREFIX] \
  --msigdb-collections H C2 C6

# Multiple genes, tables only
conda run -n claude_test python scripts/multi_gene_analysis.py [GENE1] [GENE2] [GENE3] \
  --output [PREFIX] \
  --no-plot

# Custom organism/species
conda run -n claude_test python scripts/pathway_query.py [GENE] \
  --kegg-organism [ORG] \
  --reactome-species [TAXID]

# Export to JSON
conda run -n claude_test python scripts/pathway_query.py [GENE] \
  --output json \
  --export [FILE].json
```
