# Pathway Similarity Network Builder

This script builds a pathway similarity network by collecting pathways from KEGG, MSigDB, and Reactome databases, then computing pairwise Jaccard similarity indices between all pathways.

## Overview

The network construction follows these steps:

1. **Collect pathways** from selected databases with their associated gene sets
2. **Compute Jaccard similarity** between all pathway pairs: `Jaccard = |intersection| / |union|`
3. **Output edge list** where each row represents a connection between two pathways

## Node Format

Each pathway node is formatted as: `pathway_name (database:access_id)`

Examples:
- `p53 signaling pathway - Homo sapiens (human) (KEGG:hsa04115)`
- `HALLMARK_APOPTOSIS (MSigDB:H)`
- `TP53 Regulates Transcription of DNA Repair Genes (Reactome:R-HSA-6796648)`

## Usage

### Basic Usage

```bash
# Build network with defaults (all databases, all collections, Jaccard >= 0)
# KEGG human (~370 pathways) + MSigDB all collections (~50,000) + Reactome human (~2,500)
# Warning: This will take 10+ minutes and generate millions of edges
python pathway_network_builder.py --output pathway_network.csv

# Use only specific databases
python pathway_network_builder.py --databases kegg msigdb --output network.csv

# Use specific MSigDB collections (instead of all)
python pathway_network_builder.py --msigdb-collections H C2 --output network.csv

# Filter edges by minimum Jaccard index
python pathway_network_builder.py --min-jaccard 0.1 --output network.csv
```

### Default Behavior

When run without arguments, the script:
- ✅ Collects from **all three databases** (KEGG, MSigDB, Reactome)
- ✅ Uses **all MSigDB collections** (H, C1-C8 = ~50,000 gene sets)
- ✅ Includes edges with **any Jaccard index** (min-jaccard = 0.0)
- ✅ Uses **human** organism (KEGG: hsa, Reactome: 9606)

### Testing with Limited Pathways

For testing or faster execution, limit the number of pathways:

```bash
# Collect only first 20 pathways from each database
python pathway_network_builder.py --max-pathways-per-db 20 --output test_network.csv
```

### Database-Specific Options

#### KEGG
```bash
# Specify organism (default: hsa for human)
python pathway_network_builder.py --databases kegg --kegg-org mmu --output mouse_kegg.csv
```

#### MSigDB
```bash
# Use specific collections (H, C1-C8)
python pathway_network_builder.py --databases msigdb --msigdb-collections H C2 --output msigdb_selected.csv

# Use all collections (default behavior)
python pathway_network_builder.py --databases msigdb --output msigdb_all.csv

# Specify MSigDB version
python pathway_network_builder.py --databases msigdb --msigdb-version 2025.1.Hs
```

#### Reactome
```bash
# Specify species by taxonomy ID (default: 9606 for human)
python pathway_network_builder.py --databases reactome --species 10090 --output mouse_reactome.csv
```

## Output Format

The script generates a CSV file with three columns:

| Pathway1 | Pathway2 | Jaccard_Index |
|----------|----------|---------------|
| pathway_name (db:id) | pathway_name (db:id) | 0.0 to 1.0 |

The output is sorted by Jaccard index in descending order (highest similarity first).

### Example Output

```csv
Pathway1,Pathway2,Jaccard_Index
Pentose and glucuronate interconversions (KEGG:hsa00040),Ascorbate and aldarate metabolism (KEGG:hsa00053),0.5227
Carbon metabolism (KEGG:hsa01200),Biosynthesis of amino acids (KEGG:hsa01230),0.3617
HALLMARK_ESTROGEN_RESPONSE_EARLY (MSigDB:H),HALLMARK_ESTROGEN_RESPONSE_LATE (MSigDB:H),0.3378
```

## Important Notes

### Gene ID Formats

Different databases use different gene identifier formats:

- **KEGG**: Uses gene IDs like `hsa:7157` (organism:entrez_id)
- **MSigDB**: Uses gene symbols like `TP53`
- **Reactome**: Uses gene names/symbols like `TP53`

**Impact**: Cross-database comparisons (e.g., KEGG pathway vs MSigDB gene set) will have low or zero Jaccard scores because the gene identifiers don't match. Within-database comparisons work correctly.

To enable accurate cross-database comparisons, you would need to implement gene ID normalization (e.g., converting all identifiers to gene symbols).

### Performance Considerations

#### API Rate Limiting

- **KEGG**: Rate limited to ~10 requests/second (0.1s delay between pathways)
- **Reactome**: Rate limited to ~10 requests/second (0.1s delay between pathways)
- **MSigDB**: No rate limiting (uses cached local GMT files)

#### Time Estimates

For human databases:
- **KEGG**: ~370 pathways → ~40 seconds to collect
- **MSigDB Hallmark (H)**: 50 gene sets → ~1 second to collect
- **MSigDB C2**: ~6000 gene sets → ~10 seconds to collect
- **Reactome**: ~2500 top-level pathways → ~5 minutes to collect

Network building time depends on total pathway count:
- 50 pathways = 1,225 pairs → <1 second
- 500 pathways = 124,750 pairs → ~5 seconds
- 5,000 pathways = 12,497,500 pairs → ~5 minutes

### Memory Usage

The script loads all pathway-gene mappings into memory. For large collections:
- KEGG (hsa, ~370 pathways): ~10 MB
- MSigDB (all collections, ~50,000 gene sets): ~500 MB
- Reactome (human, ~2,500 pathways): ~50 MB

## Complete Examples

### Example 1: Build KEGG-MSigDB network with filtering

```bash
python pathway_network_builder.py \
    --databases kegg msigdb \
    --kegg-org hsa \
    --msigdb-collections H \
    --min-jaccard 0.05 \
    --output kegg_msigdb_network.csv
```

### Example 2: Test with limited pathways from all databases

```bash
python pathway_network_builder.py \
    --databases kegg msigdb reactome \
    --max-pathways-per-db 20 \
    --min-jaccard 0.1 \
    --output test_all_db.csv
```

### Example 3: Build large MSigDB network with selected collections

```bash
python pathway_network_builder.py \
    --databases msigdb \
    --msigdb-collections C2 C5 \
    --min-jaccard 0.05 \
    --output msigdb_c2_c5_network.csv
```

This will create a large network from curated (C2) and ontology (C5) gene sets.

### Example 4: Build comprehensive network with ALL MSigDB collections (default)

```bash
python pathway_network_builder.py \
    --databases msigdb \
    --min-jaccard 0.05 \
    --output msigdb_all_collections.csv
```

This uses all MSigDB collections (H, C1-C8) = ~50,000 gene sets.

### Example 4: Mouse KEGG network

```bash
python pathway_network_builder.py \
    --databases kegg \
    --kegg-org mmu \
    --min-jaccard 0.1 \
    --output mouse_kegg_network.csv
```

## Downstream Analysis

The edge list CSV can be imported into network analysis tools:

### Python (NetworkX)
```python
import pandas as pd
import networkx as nx

# Load edge list
df = pd.read_csv('pathway_network.csv')

# Create graph
G = nx.from_pandas_edgelist(
    df,
    source='Pathway1',
    target='Pathway2',
    edge_attr='Jaccard_Index'
)

# Analyze
print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")
print(f"Density: {nx.density(G):.4f}")

# Find communities
from networkx.algorithms import community
communities = community.louvain_communities(G, weight='Jaccard_Index')
```

### R (igraph)
```R
library(igraph)

# Load edge list
edges <- read.csv('pathway_network.csv')

# Create graph
g <- graph_from_data_frame(edges, directed=FALSE)

# Analyze
print(paste("Nodes:", vcount(g)))
print(paste("Edges:", ecount(g)))

# Detect communities
communities <- cluster_louvain(g, weights=E(g)$Jaccard_Index)
```

### Cytoscape
1. Import CSV as network (File → Import → Network from File)
2. Set "Pathway1" as Source, "Pathway2" as Target
3. Set "Jaccard_Index" as edge weight attribute
4. Apply layout algorithms (e.g., Prefuse Force Directed)
5. Style nodes by database (color by prefix: KEGG, MSigDB, Reactome)

## Troubleshooting

### No edges found
- Lower `--min-jaccard` threshold (try 0.01 or 0.0)
- Check that pathways were collected successfully (look at log output)

### Script is slow
- Use `--max-pathways-per-db` for testing
- Consider using only one database
- For KEGG/Reactome, be patient with API rate limiting

### Memory errors
- Reduce number of pathways with `--max-pathways-per-db`
- Process one database at a time
- For MSigDB, use smaller collections (H instead of C2)

## Dependencies

Required Python packages:
```
pandas
requests
gseapy (for MSigDB)
```

Install from the existing skill environment:
```bash
conda activate claude_test
```

## Credits

This script uses the pathwaydb-query skill infrastructure for database access:
- KEGG API: https://rest.kegg.jp
- Reactome API: https://reactome.org/ContentService
- MSigDB via gseapy: https://www.gsea-msigdb.org
