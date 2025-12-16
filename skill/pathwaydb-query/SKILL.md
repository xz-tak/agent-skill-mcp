---
name: pathwaydb-query
description: Query pathway databases (KEGG, Reactome, MSigDB) to find gene-associated pathways and terms, with support for multi-gene analysis and visualization. Build pathway similarity networks using Jaccard indices. Extract gene-specific subnetworks with comprehensive centrality analysis. Use this skill when users request pathway information for genes, ask to query KEGG/Reactome/MSigDB databases, want to compare pathways across multiple genes, need to visualize pathway overlap with UpSet plots, want to construct pathway similarity networks, or need to analyze pathway centrality and importance for specific genes.
---

# Pathway Database Query

## Overview

This skill enables comprehensive querying of three major pathway databases—KEGG, Reactome, and MSigDB—to identify all pathways, gene sets, and biological terms associated with target genes. The skill supports both single-gene queries across all databases, multi-gene comparative analysis with UpSet plot visualization, and large-scale pathway similarity network construction.

**Key Capabilities:**
- Query KEGG (Kyoto Encyclopedia of Genes and Genomes) for pathway annotations
- Query Reactome for detailed pathway and reaction information
- Query MSigDB (Molecular Signatures Database) across multiple gene set collections
- Perform multi-gene analysis with pathway overlap visualization (UpSet plots)
- Build pathway similarity networks using Jaccard index across all databases
- Access pre-computed comprehensive pathway network (630M+ edges, 53K+ pathways)
- Extract gene-specific subnetworks with comprehensive centrality metrics
- Compute degree, betweenness, closeness, eigenvector centrality, PageRank, HITS, clustering
- Classify pathways as hubs, bridges, leaves, or regular nodes
- Export results in multiple formats (CSV, Excel, JSON, PNG, Parquet)
- Support organism/species-specific queries and custom parameters

## When to Use This Skill

Invoke this skill when users request:

**Explicit Database Queries:**
- "Query KEGG/Reactome/MSigDB for [GENE]"
- "Find pathways for [GENE] in KEGG"
- "What pathways is [GENE] involved in?"
- "Get MSigDB Hallmark gene sets containing [GENE]"

**General Pathway Information:**
- "What pathways is [GENE] associated with?"
- "Find all biological pathways for [GENE]"
- "Which pathways contain [GENE]?"

**Multi-Gene Comparisons:**
- "Compare pathways for [GENE1], [GENE2], [GENE3]"
- "Find shared pathways between these genes"
- "Analyze pathway overlap for [GENE LIST]"
- "Visualize pathway intersections for these genes"

**Specialized Queries:**
- "Get mouse pathways for Trp53" (organism-specific)
- "Query MSigDB C2 curated pathways for [GENE]" (collection-specific)
- "Find Reactome reactions for [GENE]" (reaction vs pathway)

**Pathway Network Queries:**
- "Build a pathway similarity network"
- "Find similar pathways across databases"
- "Create Jaccard index network for pathways"
- "Which pathways are most similar to [PATHWAY]?"
- "Analyze pathway relationships using the pre-computed network"

**Gene-Specific Network Analysis:**
- "What are the most important pathways for [GENE]?"
- "Analyze pathway centrality for [GENE]"
- "Which [GENE] pathways are most connected/central?"
- "Find hub pathways related to [GENE]"
- "Compute pathway network metrics for [GENE]"
- "Identify bridging pathways for [GENE]"

## Quick Start

### Environment Setup

Always activate the conda environment before executing scripts:

```bash
conda activate claude_test
```

### Basic Query Patterns

**Single gene, all databases:**
```bash
conda run -n claude_test python scripts/pathway_query.py TP53
```

**Multiple genes with visualization:**
```bash
conda run -n claude_test python scripts/multi_gene_analysis.py TP53 BRCA1 EGFR --output results
```

**Specific database:**
```bash
conda run -n claude_test python scripts/kegg_api.py TP53
conda run -n claude_test python scripts/reactome_api.py TP53
conda run -n claude_test python scripts/msigdb_api.py TP53 --collection H
```

## Single Gene Queries

### Query All Databases Simultaneously

Use `scripts/pathway_query.py` to query KEGG, Reactome, and MSigDB in parallel for a single gene.

**Basic usage:**
```bash
conda run -n claude_test python scripts/pathway_query.py [GENE_SYMBOL]
```

**Common parameters:**
- `--kegg-organism [ORG]` - KEGG organism code (default: "hsa" for human)
- `--reactome-species [TAXID]` - Reactome species taxonomy ID (default: 9606 for human)
- `--msigdb-collection [COLL]` - MSigDB collection (default: "H" for Hallmark)
- `--output json` - Output results in JSON format
- `--export [FILE]` - Export results to JSON file

**Example:**
```bash
conda run -n claude_test python scripts/pathway_query.py TP53 \
  --kegg-organism hsa \
  --msigdb-collection C2 \
  --export tp53_results.json
```

**Output interpretation:**
- Displays summary statistics (total pathways found per database)
- Lists all pathways from KEGG, Reactome, and MSigDB
- Can export to JSON for programmatic analysis

**When to use:**
- User wants comprehensive pathway information for one gene
- User doesn't specify a particular database
- Need results from all three databases for comparison

### Query Specific Databases

Use individual API scripts for database-specific queries with more control.

#### KEGG Queries

**Script:** `scripts/kegg_api.py`

**Usage:**
```bash
conda run -n claude_test python scripts/kegg_api.py [GENE] [ORGANISM]
```

**Parameters:**
- `[GENE]` - Gene symbol (e.g., "TP53", "BRCA1")
- `[ORGANISM]` - 3-letter organism code (default: "hsa")
  - Common: hsa (human), mmu (mouse), rno (rat), dre (zebrafish)
  - See `references/database_guide.md` for full list

**Options:**
- `--output json` - JSON output format
- `--verbose` - Detailed logging

**Example:**
```bash
# Human TP53
conda run -n claude_test python scripts/kegg_api.py TP53 hsa

# Mouse Trp53
conda run -n claude_test python scripts/kegg_api.py Trp53 mmu
```

**When to use:**
- User specifically requests KEGG pathways
- Need organism-specific pathway annotations
- Cross-species pathway comparison

#### Reactome Queries

**Script:** `scripts/reactome_api.py`

**Usage:**
```bash
conda run -n claude_test python scripts/reactome_api.py [GENE]
```

**Parameters:**
- `--species [ID]` - Species taxonomy ID or name (default: 9606)
  - Examples: 9606 or "human", 10090 or "mouse"
  - See `references/database_guide.md` for full list
- `--resource [TYPE]` - Identifier type (default: "UniProt")
  - Options: "UniProt", "NCBI", "ENSEMBL"
- `--map-to [TARGET]` - Query target (default: "pathways")
  - Options: "pathways", "reactions"

**Options:**
- `--output json` - JSON output format
- `--verbose` - Detailed logging

**Examples:**
```bash
# Basic human query
conda run -n claude_test python scripts/reactome_api.py TP53

# Mouse query
conda run -n claude_test python scripts/reactome_api.py Trp53 --species 10090

# Reactions instead of pathways
conda run -n claude_test python scripts/reactome_api.py TP53 --map-to reactions
```

**When to use:**
- User specifically requests Reactome pathways
- Need detailed molecular-level reaction information
- User wants hierarchical pathway organization

#### MSigDB Queries

**Script:** `scripts/msigdb_api.py`

**Usage:**
```bash
conda run -n claude_test python scripts/msigdb_api.py [GENE]
```

**Parameters:**
- `--collection [COLL]` or `-c [COLL]` - Single collection code (default: "H")
- `--collections [COLL1 COLL2...]` - Multiple collections
- `--all` - Query all collections (slow but comprehensive)
- `--version [VER]` or `-v [VER]` - MSigDB version (default: "2025.1.Hs")

**MSigDB Collections:**
- **H** - Hallmark gene sets (50 sets, well-characterized)
- **C2** - Curated pathways (7,561 sets, comprehensive)
- **C5** - GO terms (16,228 sets, functional annotation)
- **C6** - Oncogenic signatures (189 sets, cancer-relevant)
- **C7** - Immunologic signatures (5,219 sets, immune-related)
- **C8** - Cell type signatures (830 sets, single-cell)

*See `references/database_guide.md` for detailed collection descriptions*

**Options:**
- `--output json` - JSON output format
- `--verbose` - Detailed logging

**Examples:**
```bash
# Hallmark gene sets (default, fast)
conda run -n claude_test python scripts/msigdb_api.py TP53

# Specific collection
conda run -n claude_test python scripts/msigdb_api.py TP53 --collection C2

# Multiple collections
conda run -n claude_test python scripts/msigdb_api.py TP53 --collections H C2 C6

# All collections (comprehensive but slow)
conda run -n claude_test python scripts/msigdb_api.py TP53 --all
```

**When to use:**
- User specifically requests MSigDB gene sets
- User mentions specific collections (Hallmark, curated pathways, etc.)
- Need gene set enrichment signatures
- Cancer or immunology research context

**Collection selection guidelines:**
- **Default:** H (quick overview)
- **Comprehensive pathway analysis:** H, C2, C6
- **Cancer research:** H, C2, C6
- **Immunology:** H, C7
- **Single-cell studies:** H, C8
- **Functional annotation:** C5 (GO terms)

## Multi-Gene Analysis with Visualization

Use `scripts/multi_gene_analysis.py` for comparing pathways across multiple genes with UpSet plot visualization.

### Basic Multi-Gene Query

**Usage:**
```bash
conda run -n claude_test python scripts/multi_gene_analysis.py [GENE1] [GENE2] [GENE3] --output [PREFIX]
```

**Parameters:**
- Genes: List 2+ gene symbols as positional arguments
- `--output [PREFIX]` - Output file prefix (default: "multi_gene_analysis")
- `--msigdb-collections [COLL1 COLL2...]` - MSigDB collections to query (default: ["H"])
- `--kegg-organism [ORG]` - KEGG organism code (default: "hsa")
- `--reactome-species [ID]` - Reactome species ID (default: 9606)
- `--max-bars [N]` - Maximum bars to show in UpSet plots (default: 20)
- `--min-intersection [N]` - Minimum pathways per intersection to display (default: 1)
- `--no-plot` - Skip plot generation (tables only)

**Example:**
```bash
conda run -n claude_test python scripts/multi_gene_analysis.py TP53 BRCA1 EGFR \
  --output cancer_genes \
  --msigdb-collections H C2 C6 \
  --max-bars 15
```

### Output Files Generated

For output prefix "results", the script generates:

**CSV Tables:**
- `results_wide.csv` - One row per gene, pathways in semicolon-separated lists
- `results_long.csv` - One row per gene-pathway pair (best for analysis)
- `results_summary.csv` - Pathway counts by gene and database

**Excel File:**
- `results_wide.xlsx` - Multi-sheet workbook (Summary, KEGG, REACTOME, MSIGDB sheets)

**UpSet Plots (PNG, 300 DPI):**
- `results_kegg.png` - KEGG pathway overlap
- `results_reactome.png` - Reactome pathway overlap
- `results_msigdb.png` - MSigDB gene set overlap
- `results_combined.png` - All databases combined

*See `references/output_interpretation.md` for detailed file format descriptions*

### UpSet Plot Interpretation

UpSet plots visualize pathway overlap patterns across genes:

**Visual components:**
- **Dots (●)** indicate which genes share pathways
- **Bar height** shows number of pathways in that intersection
- **Bars sorted by size** (most pathways first)

**Key insights:**
- All dots filled → pathways shared by all genes (core biological functions)
- Two dots filled → pairwise overlap (specific functional relationships)
- One dot filled → gene-specific pathways (unique functions)

**Example interpretation:**
```
If analyzing TP53, BRCA1, EGFR:
- Bar with all 3 dots: Pathways common to all (e.g., DNA damage response)
- Bar with TP53 + BRCA1 dots: Shared tumor suppressor pathways
- Bar with only TP53 dot: TP53-specific pathways (e.g., p53 signaling)
```

### When to Use Multi-Gene Analysis

**Use this approach when:**
- User provides 2+ genes to compare
- User asks about "shared pathways" or "pathway overlap"
- User wants to "compare" or "analyze" multiple genes
- User requests visualization of pathway intersections

**Recommended gene counts:**
- **2-5 genes:** Ideal for UpSet plots and detailed comparison
- **5-10 genes:** Good for pathway enrichment, consider `--max-bars` reduction
- **10+ genes:** Use `--no-plot` and focus on table exports

**Visualization parameters:**
- Default (`--max-bars 20`) works well for 2-5 genes
- Reduce to `--max-bars 10-15` if plot is cluttered
- Increase `--min-intersection 2-3` to show only significant overlaps

## Understanding and Presenting Results

### Reading Query Outputs

Refer to `references/output_interpretation.md` for comprehensive output format documentation.

**Quick reference:**

1. **Text output** - Parse summary statistics and pathway lists from terminal output
2. **JSON output** - Load with pandas/json for programmatic analysis
3. **Wide CSV** - Quick overview, one row per gene
4. **Long CSV** - Best for filtering, grouping, and statistical analysis
5. **Summary CSV** - Pathway count comparison
6. **Excel file** - Multi-sheet format for sharing

### Common Analysis Tasks

**Find shared pathways:**
```python
import pandas as pd
long_df = pd.read_csv("results_long.csv")
n_genes = long_df['Gene'].nunique()
pathway_counts = long_df.groupby('Pathway')['Gene'].nunique()
shared = pathway_counts[pathway_counts == n_genes].index.tolist()
```

**Find gene-specific pathways:**
```python
import pandas as pd
long_df = pd.read_csv("results_long.csv")
tp53_pathways = set(long_df[long_df['Gene'] == 'TP53']['Pathway'])
other_pathways = set(long_df[long_df['Gene'] != 'TP53']['Pathway'])
tp53_specific = tp53_pathways - other_pathways
```

*See `references/output_interpretation.md` for more examples*

### Presenting Results to Users

**For summary requests:**
- Show total pathway counts per database
- Highlight top 3-5 most relevant pathways
- Mention file locations for detailed results

**For comprehensive results:**
- Display full pathway counts
- List key pathways (top 20 if many results)
- Enumerate all generated files
- Offer to filter or analyze further

**For multi-gene comparisons:**
- Report pathway overlap statistics
- Highlight shared vs unique pathways
- Point to UpSet plot files for visualization
- Interpret interesting intersections

**For exports:**
- List all generated files with brief descriptions
- Recommend which file to use for their specific needs
- Suggest next steps (R/Python analysis, Excel browsing, etc.)

*See `references/output_interpretation.md` for detailed presentation guidelines*

## Pathway Similarity Network Analysis

### Overview

The skill includes tools to build and analyze pathway similarity networks based on gene set overlap. Each pathway from KEGG, Reactome, and MSigDB is represented as a node, and edges connect pathways that share genes, weighted by their Jaccard similarity index (intersection over union).

**Use cases:**
- Identify functionally related pathways across databases
- Find pathway clusters and communities
- Discover pathway redundancy and complementarity
- Build network visualizations for systems biology analysis
- Integrate with graph analysis tools (NetworkX, Cytoscape, igraph)

### Pre-Computed Comprehensive Network

A pre-computed pathway similarity network is available in the skill:

**Location:** `data/all_pathway_network_12052025.parquet`

**Network Statistics:**
- **Nodes:** 52,870 pathways
  - KEGG human: ~370 pathways
  - MSigDB all collections: ~50,000 gene sets
  - Reactome human: ~2,500 pathways
- **Edges:** 630,213,753 pairwise comparisons
- **Format:** Parquet (11.53 GB compressed from 62 GB CSV)
- **Columns:** `Pathway1`, `Pathway2`, `Jaccard_Index`

**Node Format:**
- KEGG: `pathway_name (KEGG:pathway_id)` - e.g., `p53 signaling pathway (KEGG:hsa04115)`
- MSigDB: `gene_set_name (MSigDB:collection)` - e.g., `HALLMARK_APOPTOSIS (MSigDB:H)`
- Reactome: `pathway_name (Reactome:stId)` - e.g., `TP53 Regulates Transcription (Reactome:R-HSA-6796648)`

**Citation for network:**
```
Pre-computed pathway similarity network (KEGG, MSigDB, Reactome)
Generated: December 5, 2025
Databases: KEGG hsa (human), MSigDB v2025.1.Hs (all collections H, C1-C8), Reactome v84 (human, species 9606)
Method: Jaccard similarity index (gene set overlap)
Total edges: 630,213,753 pairwise comparisons
Available at: .claude/skills/pathwaydb-query/data/all_pathway_network_12052025.parquet
```

### Building Custom Pathway Networks

Use `network_builder/pathway_network_builder.py` to build custom networks with different parameters.

**Basic usage (builds full network with defaults):**
```bash
conda run -n claude_test python network_builder/pathway_network_builder.py \
  --output custom_network.parquet
```

**Default behavior:**
- All three databases: KEGG, MSigDB (all collections), Reactome
- Human organism: KEGG (hsa), Reactome (9606)
- No Jaccard filtering (min-jaccard = 0.0)
- Output format: Parquet (compressed)

**Common parameters:**
- `--databases [kegg|msigdb|reactome]` - Select specific databases
- `--msigdb-collections [H C2 C6...]` - Select MSigDB collections (default: all)
- `--min-jaccard [THRESHOLD]` - Filter edges by minimum similarity (default: 0.0)
- `--kegg-org [ORG]` - KEGG organism code (default: hsa)
- `--species [TAXID]` - Reactome species (default: 9606)
- `--max-pathways-per-db [N]` - Limit pathways for testing
- `--output [FILE]` - Output filename (.parquet or .csv)

**Examples:**

```bash
# Build network with only MSigDB Hallmark gene sets
conda run -n claude_test python network_builder/pathway_network_builder.py \
  --databases msigdb \
  --msigdb-collections H \
  --min-jaccard 0.1 \
  --output hallmark_network.parquet

# Build network with KEGG and MSigDB curated pathways
conda run -n claude_test python network_builder/pathway_network_builder.py \
  --databases kegg msigdb \
  --msigdb-collections C2 \
  --min-jaccard 0.05 \
  --output curated_network.parquet

# Build test network with limited pathways
conda run -n claude_test python network_builder/pathway_network_builder.py \
  --max-pathways-per-db 50 \
  --min-jaccard 0.05 \
  --output test_network.parquet

# Build mouse pathway network
conda run -n claude_test python network_builder/pathway_network_builder.py \
  --kegg-org mmu \
  --species 10090 \
  --databases kegg reactome \
  --output mouse_network.parquet
```

**Performance notes:**
- Full network with all databases: 15-40 minutes
- KEGG collection: ~40 seconds (~370 pathways)
- MSigDB single collection: <1 second to 2 minutes depending on size
- Reactome collection: ~5 minutes (~2,500 pathways)
- Parquet output is 5-10x smaller than CSV

**Output files:**
- `.parquet` format (default, recommended) - Compressed, fast to read
- `.csv` format (optional) - For tools without parquet support

### Analyzing Pathway Networks

Use `network_builder/analyze_network_example.py` to compute network statistics.

**Basic usage:**
```bash
conda run -n claude_test python network_builder/analyze_network_example.py \
  network.parquet
```

**Outputs:**
- Network statistics (nodes, edges, Jaccard distribution)
- Edge counts by database combination
- Top 10 most similar pathway pairs
- Degree distribution (most connected pathways)

**Example with pre-computed network:**
```bash
conda run -n claude_test python network_builder/analyze_network_example.py \
  data/all_pathway_network_12052025.parquet
```

**Programmatic analysis:**
```python
import pandas as pd
import pyarrow.parquet as pq

# Load network (fast, columnar format)
df = pd.read_parquet('data/all_pathway_network_12052025.parquet')

# Filter by Jaccard threshold
high_similarity = df[df['Jaccard_Index'] > 0.1]

# Find pathways similar to a specific pathway
target = "HALLMARK_APOPTOSIS (MSigDB:H)"
similar = df[
    (df['Pathway1'] == target) | (df['Pathway2'] == target)
].sort_values('Jaccard_Index', ascending=False).head(20)

# Extract database from pathway name
df['DB1'] = df['Pathway1'].str.extract(r'\(([^:]+):')[0]
df['DB2'] = df['Pathway2'].str.extract(r'\(([^:]+):')[0]

# Cross-database edges
cross_db = df[df['DB1'] != df['DB2']]
```

**Integration with graph analysis tools:**

```python
# NetworkX
import networkx as nx
G = nx.from_pandas_edgelist(df, 'Pathway1', 'Pathway2', 'Jaccard_Index')
communities = nx.community.louvain_communities(G, weight='Jaccard_Index')

# igraph (R)
library(igraph)
edges <- read_parquet('network.parquet')
g <- graph_from_data_frame(edges, directed=FALSE)
communities <- cluster_louvain(g, weights=E(g)$Jaccard_Index)

# Cytoscape
# Import parquet/CSV directly as network
# Style nodes by database, edge width by Jaccard_Index
```

### Converting Between Formats

Use `network_builder/convert_csv_to_parquet.py` to convert existing CSV networks to Parquet:

```bash
conda run -n claude_test python network_builder/convert_csv_to_parquet.py \
  input.csv output.parquet
```

**Benefits of Parquet:**
- 5-10x smaller file size
- 3-5x faster to read
- Better type preservation
- Efficient columnar filtering

### When to Use Network Analysis

**Use pre-computed network when:**
- User wants comprehensive pathway relationships
- Need to identify similar pathways across databases
- Building pathway visualization or clustering analysis
- Comparative pathway studies

**Build custom network when:**
- Need specific organism (mouse, rat, etc.)
- Want subset of databases or collections
- Require different Jaccard threshold
- Need fresh data with updated databases

**Analysis script when:**
- Need quick statistics about network
- Want to identify top similar pathways
- Checking network connectivity patterns
- Generating reports about pathway relationships

### Network Documentation

Comprehensive documentation available in `network_builder/`:
- `README_pathway_network.md` - Full usage guide with examples
- `PARQUET_FORMAT.md` - Format specifications and benefits
- `DEFAULT_BEHAVIOR.md` - Detailed explanation of defaults

*See `references/output_interpretation.md` for detailed presentation guidelines*

## Gene-Specific Subnetwork Analysis

### Overview

For each gene query, extract the relevant subnetwork from the full pathway network and compute comprehensive centrality metrics. This identifies which pathways are most central and interconnected within the gene's functional context.

**Use case:** After querying pathways for a gene (e.g., TP53), extract those pathways from the network to analyze their relationships, similarity, and importance within the broader pathway landscape.

### Quick Start

```bash
# Step 1: Query pathways for gene (creates pathway list)
conda run -n claude_test python scripts/pathway_query.py TP53 --export tp53_pathways.json

# Step 2: Create pathway list CSV (extract pathway names)
# (pathway names must match network format: "pathway_name (database:id)")

# Step 3: Run subnetwork analysis
conda run -n claude_test python network_builder/gene_subnetwork_analysis.py \
  --gene TP53 \
  --network data/all_pathway_network_12052025.parquet \
  --pathways tp53_pathway_list.csv \
  --min-jaccard 0.05 \
  --output tp53_subnetwork
```

### Usage

**Basic usage (seed pathways only, default):**
```bash
conda run -n claude_test python network_builder/gene_subnetwork_analysis.py \
  --gene GENE_SYMBOL \
  --network data/all_pathway_network_12052025.parquet \
  --pathways pathway_results.csv \
  --output gene_subnetwork
```

**With neighbor expansion:**
```bash
conda run -n claude_test python network_builder/gene_subnetwork_analysis.py \
  --gene TP53 \
  --network data/all_pathway_network_12052025.parquet \
  --pathways tp53_pathways.csv \
  --min-jaccard 0.1 \
  --include-neighbors \
  --max-distance 1 \
  --output tp53_subnetwork
```

**With manual pathway list:**
```bash
conda run -n claude_test python network_builder/gene_subnetwork_analysis.py \
  --gene BRCA1 \
  --network data/all_pathway_network_12052025.parquet \
  --pathway-list "HALLMARK_DNA_REPAIR (MSigDB:H),p53 signaling pathway (KEGG:hsa04115)" \
  --output brca1_subnetwork
```

### Parameters

- `--gene` (required): Gene symbol for labeling output
- `--network` (required): Path to full pathway network (parquet or CSV)
- `--pathways`: Pathway query result files (CSV/Excel) - accepts multiple files
- `--pathway-list`: Comma-separated pathway names (alternative to --pathways)
- `--min-jaccard`: Minimum Jaccard index threshold (default: 0.0, recommend: 0.05-0.1)
- `--include-neighbors`: Expand subnetwork to include neighboring pathways (default: False)
- `--max-distance`: Maximum distance for neighbors (default: 1, only applies with --include-neighbors)
- `--output`: Output file prefix (default: gene_subnetwork)

### Output Files

**1. Centrality Table** (`{output}_centrality.csv`)
- One row per pathway in subnetwork
- Columns: Pathway, Database, Is_Seed, Node_Type, Degree, Degree_Centrality,
  Closeness_Centrality, Betweenness_Centrality, Eigenvector_Centrality,
  PageRank, Hub_Score, Authority_Score, Clustering_Coefficient
- Sorted by Degree_Centrality (highest first)

**2. Comprehensive Analysis** (`{output}_analysis.xlsx`)
Multi-sheet Excel workbook:
- **Centrality**: Full centrality metrics table
- **Network_Stats**: Overall subnetwork statistics
- **Top_Hubs**: Top 20 most connected pathways
- **Top_Bridges**: Top 20 bridging pathways (betweenness)
- **Seed_Pathways**: Metrics for input gene's pathways

**3. Edge List** (`{output}_edgelist.csv`)
- Network edges with Source, Target columns
- Import into Cytoscape, Gephi, or other visualization tools

### Centrality Metrics Explained

- **Degree**: Number of direct connections to other pathways
- **Degree Centrality**: Normalized degree (0-1 scale)
- **Closeness Centrality**: Average distance to all other pathways (based on distance = 1/Jaccard)
- **Betweenness Centrality**: Frequency on shortest paths between pathways (bridging role)
- **Eigenvector Centrality**: Influence based on connections to influential pathways (uses Jaccard as weight)
- **PageRank**: Google's ranking algorithm (uses Jaccard as weight)
- **Hub Score**: Authority as information source (many connections)
- **Authority Score**: Authority as information destination (quality connections)
- **Clustering Coefficient**: How interconnected neighbors are (weighted by Jaccard)

**Weight handling:**
- Similarity-based metrics (Eigenvector, PageRank, Clustering): Use Jaccard Index directly (higher = stronger connection)
- Distance-based metrics (Closeness, Betweenness): Use Distance = 1/Jaccard (higher Jaccard = lower distance)

### Node Classification

Pathways automatically classified as:
- **Hub**: Highly connected (degree ≥ 50% of average × network size)
- **Bridge**: High betweenness (> 0.1), connects different pathway modules
- **Leaf**: Single connection only (degree = 1)
- **Regular**: Standard connectivity

### Neighbor Expansion

**Without `--include-neighbors` (default):**
- Subnetwork contains only seed pathways directly associated with the gene
- Example: TP53 in 3 pathways → 3 node subnetwork
- Best for analyzing gene's direct pathway involvement

**With `--include-neighbors`:**
- Subnetwork includes seed pathways + all connected pathways
- `--max-distance 1`: Direct neighbors only
- `--max-distance 2`: Neighbors of neighbors (2-hop)
- Example: TP53 in 3 pathways → 763 node subnetwork with neighbors
- Best for understanding broader biological context

### Example Results Interpretation

```
================================================================================
Gene-Specific Pathway Subnetwork Analysis: TP53
================================================================================
Seed pathways: 3
Subnetwork nodes: 3
Subnetwork edges: 0
Network density: 0.0000
Average degree: 0.00

Node classification:
  Hubs: 3
  Bridges: 0
  Leaves: 0
  Regular: 0
```

**Interpretation (without neighbors):**
- 3 seed pathways for TP53 (from pathway query)
- 0 edges between them at this Jaccard threshold
- Suggests the 3 TP53 pathways have distinct gene sets
- May need lower Jaccard threshold or neighbor expansion

```
================================================================================
Gene-Specific Pathway Subnetwork Analysis: TP53
================================================================================
Seed pathways: 3
Subnetwork nodes: 763
Subnetwork edges: 58854
Network density: 0.2025
Average degree: 154.27

Node classification:
  Hubs: 762
  Bridges: 0
  Leaves: 1
  Regular: 0

Top 5 most central pathways (by degree):
  HALLMARK_E2F_TARGETS (MSigDB:H)
    Degree: 627, Centrality: 0.8228
```

**Interpretation (with neighbors):**
- Started with 3 seed pathways, expanded to 763 related pathways
- Dense network (20% of possible connections exist)
- Most pathways are hubs (highly interconnected module)
- HALLMARK_E2F_TARGETS most central (82% centrality)
- Suggests TP53 pathways are part of tightly connected cell cycle/DNA repair module

### Integration with Pathway Queries

**Typical workflow:**

```bash
# 1. Query pathways for gene across all databases
conda run -n claude_test python scripts/pathway_query.py TP53 \
  --kegg-organism hsa \
  --msigdb-collection H \
  --export tp53_query.json

# 2. Extract pathway names to CSV (manual or programmatic)
# Create file with column "Pathway" containing pathway names in network format:
# HALLMARK_P53_PATHWAY (MSigDB:H)
# p53 signaling pathway (KEGG:hsa04115)
# TP53 Regulates Transcription (Reactome:R-HSA-6796648)

# 3. Run subnetwork analysis
conda run -n claude_test python network_builder/gene_subnetwork_analysis.py \
  --gene TP53 \
  --network data/all_pathway_network_12052025.parquet \
  --pathways tp53_pathway_list.csv \
  --min-jaccard 0.05 \
  --output tp53_subnetwork

# 4. Summarize and interpret results
conda run -n claude_test python network_builder/summarize_centrality_analysis.py \
  tp53_subnetwork_centrality.csv \
  --top 10

# 5. Detailed analysis
# - Open tp53_subnetwork_analysis.xlsx for interactive exploration
# - Check tp53_subnetwork_centrality.csv for detailed metrics
# - Import tp53_subnetwork_edgelist.csv into Cytoscape for visualization
```

### Performance Notes

- Loading full network (11GB parquet): 30 seconds - 2 minutes
- Extracting subnetwork: <1 second
- Computing centrality metrics: 30 seconds - 2 minutes (depends on size)
- Memory usage: 5-10 GB (with min-jaccard filtering)

### Tips

1. **Use Jaccard filtering**: `--min-jaccard 0.05` or `0.1` reduces network size and computation time

2. **Choose neighbor expansion based on goal:**
   - No neighbors: Analyze gene's direct pathway roles
   - With neighbors: Understand broader functional context

3. **Pathway name format matters:**
   - Must match network format: `"pathway_name (database:id)"`
   - Examples:
     - `HALLMARK_APOPTOSIS (MSigDB:H)`
     - `p53 signaling pathway (KEGG:hsa04115)`
     - `TP53 Regulates Transcription (Reactome:R-HSA-6796648)`

4. **For visualization:**
   - Import edge list CSV into Cytoscape
   - Style nodes by Database column
   - Size nodes by Degree or Centrality
   - Color by Node_Type (Hub/Bridge/Leaf/Regular)

## Advanced Usage Patterns

### Cross-Species Pathway Comparison

Compare pathways between human and mouse orthologs:

```bash
# Human TP53
conda run -n claude_test python scripts/pathway_query.py TP53 \
  --kegg-organism hsa \
  --reactome-species 9606 \
  --export human_tp53.json

# Mouse Trp53
conda run -n claude_test python scripts/pathway_query.py Trp53 \
  --kegg-organism mmu \
  --reactome-species 10090 \
  --export mouse_trp53.json
```

Then compare JSON outputs programmatically or review manually.

### Batch Gene Queries

For querying many genes programmatically:

```python
from pathway_query import query_all_databases
import pandas as pd

genes = ["TP53", "BRCA1", "EGFR", "KRAS", "MYC"]
results = []

for gene in genes:
    print(f"Querying {gene}...")
    result = query_all_databases(gene, parallel=True)
    results.append({
        "Gene": gene,
        "KEGG_Pathways": result["databases"]["kegg"]["total_pathways"],
        "Reactome_Pathways": result["databases"]["reactome"]["total_terms"],
        "MSigDB_GeneSets": result["databases"]["msigdb"]["total_gene_sets"]
    })

df = pd.DataFrame(results)
df.to_csv("batch_query_summary.csv", index=False)
```

### MSigDB Collection Exploration

Explore which collections contain a gene:

```bash
# Query each collection separately
conda run -n claude_test python scripts/msigdb_api.py TP53 --collection H
conda run -n claude_test python scripts/msigdb_api.py TP53 --collection C2
conda run -n claude_test python scripts/msigdb_api.py TP53 --collection C6

# Or query all at once (slower)
conda run -n claude_test python scripts/msigdb_api.py TP53 --all
```

## Troubleshooting Common Issues

### Gene Not Found

**Symptom:** No results for a gene in specific database

**Solutions:**
1. Verify gene symbol spelling (case-sensitive for some databases)
2. Try alternative gene names (e.g., "CDKN1A" vs "P21")
3. For KEGG: check organism code is correct
4. For Reactome: try different identifier types (UniProt, NCBI, ENSEMBL)
5. Some genes may not be annotated in certain databases

### Slow MSigDB Queries

**Symptom:** First query to MSigDB collection takes 1-2 minutes

**Explanation:** GMT files are downloaded on first use; subsequent queries use cache

**Solutions:**
- Limit to specific collections (H, C2, C6) instead of `--all`
- First download is slow; following queries are fast
- Collections C5 and C7 are very large

### Empty or Sparse UpSet Plots

**Symptom:** Very few bars or no visualization

**Explanation:** Genes have little pathway overlap

**Solutions:**
- Reduce `--min-intersection` to 1 (show all overlaps)
- Check individual pathway lists in CSV files
- Inform user that genes have distinct biological functions

### No Visualization Generated

**Symptom:** CSV files created but no PNG plots

**Check:**
1. Verify `upsetplot` is installed: `pip list | grep upsetplot`
2. If using `--no-plot` flag, plots are intentionally skipped
3. Check for error messages in terminal output

## Reference Documentation

This skill includes comprehensive reference documentation in the `references/` directory. Load these files when detailed information is needed:

### references/database_guide.md

**Contains:**
- Detailed descriptions of KEGG, Reactome, and MSigDB
- Complete organism code lists (KEGG)
- Species taxonomy ID lists (Reactome)
- MSigDB collection descriptions and use cases
- Database comparison and selection guidelines
- API rate limits and best practices

**When to reference:**
- User asks about organism/species codes
- Need to explain MSigDB collections
- Choosing which database is best for user's needs
- Understanding database-specific limitations

### references/query_examples.md

**Contains:**
- Common query patterns with examples
- Advanced workflow templates
- User request → command mapping
- Batch processing examples
- Performance optimization tips

**When to reference:**
- Uncertain which command to use for user's request
- Need workflow examples for complex queries
- User wants batch processing or programmatic usage
- Optimizing queries for performance

### references/output_interpretation.md

**Contains:**
- Detailed output format descriptions
- File format comparisons (wide vs long vs summary)
- UpSet plot interpretation guidelines
- Result presentation strategies
- Analysis code examples (pandas, R)

**When to reference:**
- Need to explain output files to user
- User asks how to analyze results
- Presenting multi-gene comparison findings
- Recommending which file format to use

## Resources

This skill includes:

### scripts/

Executable Python tools for querying pathway databases:

- **pathway_query.py** - Unified interface to query all databases for single gene
- **multi_gene_analysis.py** - Multi-gene analysis with UpSet plot visualization
- **kegg_api.py** - KEGG pathway database API
- **reactome_api.py** - Reactome pathway database API
- **msigdb_api.py** - MSigDB gene set database API

Execute these scripts directly via conda environment. Scripts may also be read for understanding implementation details.

### network_builder/

Tools for building and analyzing pathway similarity networks:

- **pathway_network_builder.py** - Build custom pathway networks with Jaccard similarity
- **analyze_network_example.py** - Compute network statistics and identify top connections
- **gene_subnetwork_analysis.py** - Extract gene-specific subnetworks with centrality metrics
- **summarize_centrality_analysis.py** - Summarize and interpret centrality results with biological insights
- **convert_csv_to_parquet.py** - Convert CSV networks to efficient Parquet format
- **README_pathway_network.md** - Comprehensive usage guide and examples
- **PARQUET_FORMAT.md** - Format specifications and benefits
- **DEFAULT_BEHAVIOR.md** - Detailed explanation of default parameters

### data/

Pre-computed pathway networks:

- **all_pathway_network_12052025.parquet** - Comprehensive pathway similarity network
  - 52,870 pathways (KEGG, MSigDB all collections, Reactome)
  - 630,213,753 edges with Jaccard similarity weights
  - 11.53 GB Parquet format (compressed from 62 GB CSV)
  - Generated December 5, 2025

### references/

Documentation loaded into context when needed:

- **database_guide.md** - Comprehensive database reference (organisms, collections, parameters)
- **query_examples.md** - Query patterns and workflow examples
- **output_interpretation.md** - Output format guide and analysis templates

Reference these files for detailed information beyond what's in this SKILL.md.

---

## Key Principles

1. **Always activate conda environment** before script execution
2. **Default to comprehensive queries** (all databases) unless user specifies otherwise
3. **Use multi-gene analysis** when user provides 2+ genes
4. **Use pre-computed network** for pathway similarity queries when possible
5. **Reference documentation files** for detailed organism/collection information
6. **Present results with context** - explain what the numbers mean biologically
7. **List all output files** - users need to know what was generated
8. **Cite the pre-computed network** when using data/all_pathway_network_12052025.parquet
9. **Offer follow-up analysis** - filtering, visualization, statistical analysis
10. **Handle errors gracefully** - gene not found is common, provide alternatives
