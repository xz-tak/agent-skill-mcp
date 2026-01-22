# Pathway Network Builder - Skill Addition Summary

**Date Added:** December 11, 2025
**Added by:** Claude Code
**Version:** 1.0

## Overview

The pathwaydb-query skill has been enhanced with comprehensive pathway similarity network analysis capabilities. This addition enables building and analyzing large-scale pathway networks using Jaccard similarity indices across KEGG, Reactome, and MSigDB databases.

## What Was Added

### 1. Network Builder Scripts (`network_builder/`)

Three Python scripts for building and analyzing pathway networks:

- **`pathway_network_builder.py`** - Main script to build custom pathway networks
  - Collects all pathways from KEGG, MSigDB, Reactome
  - Computes pairwise Jaccard similarity (gene set overlap)
  - Outputs edge list with pathway pairs and similarity weights
  - Default: Parquet format (5-10x smaller than CSV)
  - Configurable: databases, collections, organisms, thresholds

- **`analyze_network_example.py`** - Network analysis and statistics
  - Computes network statistics (nodes, edges, distributions)
  - Identifies top similar pathway pairs
  - Analyzes degree distribution and connectivity
  - Supports both CSV and Parquet formats

- **`convert_csv_to_parquet.py`** - Format conversion utility
  - Converts large CSV networks to efficient Parquet format
  - Uses PyArrow streaming for 60+ GB files
  - Compression ratio: ~5-10x smaller file size

### 2. Pre-Computed Comprehensive Network (`data/`)

A complete pathway similarity network covering all major databases:

**File:** `data/all_pathway_network_12052025.parquet`

**Statistics:**
- **Size:** 11.53 GB (compressed from 62 GB CSV)
- **Nodes:** 52,870 pathways
  - KEGG human (hsa): ~370 pathways
  - MSigDB v2025.1.Hs (H, C1-C8): ~50,000 gene sets
  - Reactome v84 (human, 9606): ~2,500 pathways
- **Edges:** 630,213,753 pairwise comparisons
- **Metric:** Jaccard similarity index (gene set intersection/union)
- **Generated:** December 5, 2025

**Network Format:**
```
Pathway1,Pathway2,Jaccard_Index
p53 signaling pathway (KEGG:hsa04115),Apoptosis (KEGG:hsa04210),0.3456
HALLMARK_APOPTOSIS (MSigDB:H),HALLMARK_P53_PATHWAY (MSigDB:H),0.2891
TP53 Regulates Transcription (Reactome:R-HSA-6796648),DNA Damage Response (Reactome:R-HSA-69473),0.4123
```

### 3. Documentation

Three comprehensive documentation files in `network_builder/`:

- **`README_pathway_network.md`** - Complete usage guide
  - Installation and setup
  - Command-line usage examples
  - Parameter explanations
  - Output format descriptions
  - Integration with graph analysis tools (NetworkX, igraph, Cytoscape)
  - Performance benchmarks
  - Troubleshooting guide

- **`PARQUET_FORMAT.md`** - Format specifications
  - Why Parquet vs CSV
  - Compression benefits and benchmarks
  - Reading/writing Parquet files
  - Tool compatibility
  - Migration guide

- **`DEFAULT_BEHAVIOR.md`** - Default parameters reference
  - Detailed explanation of all defaults
  - Expected output sizes and runtimes
  - Recommendations for different use cases
  - Override examples

### 4. Updated SKILL.md

The main skill documentation now includes:

- New pathway network capabilities in description
- "Pathway Similarity Network Analysis" section
- Citation format for the pre-computed network
- Usage examples for network builder
- Integration examples (NetworkX, igraph, Cytoscape)
- Updated Resources section listing network_builder/ and data/
- Updated Key Principles to include network usage

## Citation Format

When using the pre-computed network, cite as:

```
Pre-computed pathway similarity network (KEGG, MSigDB, Reactome)
Generated: December 5, 2025
Databases: KEGG hsa (human), MSigDB v2025.1.Hs (all collections H, C1-C8), Reactome v84 (human, species 9606)
Method: Jaccard similarity index (gene set overlap)
Total edges: 630,213,753 pairwise comparisons
Available at: .claude/skills/pathwaydb-query/data/all_pathway_network_12052025.parquet
```

## Use Cases

The network addition enables:

1. **Pathway Similarity Analysis**
   - Find pathways similar to a target pathway
   - Identify functionally related pathways across databases
   - Discover pathway redundancy and complementarity

2. **Network Clustering**
   - Detect pathway communities
   - Identify functional modules
   - Find pathway hubs and connectors

3. **Cross-Database Integration**
   - Link equivalent pathways across KEGG, MSigDB, Reactome
   - Build unified pathway representations
   - Integrate multiple pathway annotations

4. **Systems Biology**
   - Pathway enrichment with network context
   - Multi-scale pathway analysis
   - Integration with other omics networks

5. **Visualization**
   - Create pathway network visualizations in Cytoscape
   - Build interactive network apps
   - Generate publication-quality figures

## Technical Details

### Default Behavior

Running with no arguments:
```bash
python network_builder/pathway_network_builder.py
```

Collects:
- All KEGG human pathways (~370)
- All MSigDB collections H, C1-C8 (~50,000 gene sets)
- All Reactome human pathways (~2,500)
- No Jaccard filtering (includes all edges)
- Output: `pathway_network_edges.parquet`

### Performance

| Operation | Time | Memory |
|-----------|------|--------|
| Full network build | 15-40 min | ~2-4 GB |
| KEGG collection | ~40 sec | ~100 MB |
| MSigDB single collection | <1-2 min | ~200 MB |
| Reactome collection | ~5 min | ~500 MB |
| Load network (Parquet) | ~10 sec | ~25 GB |
| Load network (filtered) | ~5 sec | ~5 GB |

### File Sizes

| Network Type | CSV Size | Parquet Size | Compression |
|--------------|----------|--------------|-------------|
| Full (630M edges) | 62 GB | 11.53 GB | 5.38x |
| Filtered (J>0.05) | ~2-5 GB | ~0.5-1 GB | ~5x |
| Hallmark only | ~50 MB | ~10 MB | ~5x |

## Example Workflows

### Workflow 1: Find Similar Pathways

```python
import pandas as pd

# Load network
df = pd.read_parquet('data/all_pathway_network_12052025.parquet')

# Find pathways similar to target
target = "HALLMARK_APOPTOSIS (MSigDB:H)"
similar = df[
    ((df['Pathway1'] == target) | (df['Pathway2'] == target)) &
    (df['Jaccard_Index'] > 0.1)
].sort_values('Jaccard_Index', ascending=False)

print(f"Found {len(similar)} pathways similar to {target}")
print(similar.head(10))
```

### Workflow 2: Build Custom Network

```bash
# Build network with only curated pathways
conda run -n claude_test python network_builder/pathway_network_builder.py \
  --databases kegg msigdb \
  --msigdb-collections C2 \
  --min-jaccard 0.05 \
  --output curated_network.parquet

# Analyze the network
conda run -n claude_test python network_builder/analyze_network_example.py \
  curated_network.parquet
```

### Workflow 3: Network Clustering

```python
import pandas as pd
import networkx as nx

# Load and filter network
df = pd.read_parquet('data/all_pathway_network_12052025.parquet')
df_filtered = df[df['Jaccard_Index'] > 0.1]

# Create graph
G = nx.from_pandas_edgelist(
    df_filtered,
    source='Pathway1',
    target='Pathway2',
    edge_attr='Jaccard_Index'
)

# Detect communities
communities = nx.community.louvain_communities(G, weight='Jaccard_Index')

print(f"Found {len(communities)} pathway communities")
for i, comm in enumerate(communities[:5], 1):
    print(f"\nCommunity {i}: {len(comm)} pathways")
    for pathway in list(comm)[:3]:
        print(f"  - {pathway}")
```

## Integration Points

The network builder integrates with existing skill capabilities:

1. **Gene Queries** → Find pathways → Analyze pathway network
2. **Multi-gene Analysis** → Identify shared pathways → Check similarity
3. **Database Queries** → Build custom network → Cluster pathways

## Future Enhancements

Potential additions for future versions:

1. Gene ID normalization for accurate cross-database comparison
2. Statistical significance testing (hypergeometric, permutation)
3. Pathway enrichment with network context
4. Interactive network visualization tools
5. Additional similarity metrics (overlap coefficient, cosine)
6. Temporal network updates tracking
7. Integration with other pathway databases (WikiPathways, BioCarta)

## Files Added

### Scripts (executable)
```
network_builder/
├── pathway_network_builder.py       (23 KB)
├── analyze_network_example.py       (4 KB)
└── convert_csv_to_parquet.py        (3 KB)
```

### Documentation
```
network_builder/
├── README_pathway_network.md        (35 KB)
├── PARQUET_FORMAT.md                (18 KB)
└── DEFAULT_BEHAVIOR.md              (15 KB)
```

### Data
```
data/
└── all_pathway_network_12052025.parquet  (11.53 GB)
```

### Updated
```
SKILL.md                             (Updated with network section)
```

**Total Addition Size:** ~11.6 GB (11.53 GB network + ~100 KB scripts/docs)

## Maintenance

### Updating the Pre-Computed Network

To regenerate the comprehensive network with updated databases:

```bash
cd /home/sagemaker-user/.claude/skills/pathwaydb-query

conda run -n claude_test python network_builder/pathway_network_builder.py \
  --output data/all_pathway_network_YYYYMMDD.parquet

# Update SKILL.md citation with new date and statistics
```

Recommended update frequency: Every 6-12 months or when databases release major updates.

### Version History

- **v1.0 (December 11, 2025)** - Initial addition
  - Network builder scripts
  - Pre-computed network (December 5, 2025)
  - Comprehensive documentation
  - SKILL.md updates

## Contact & Support

For questions or issues with the network functionality:

1. Check `network_builder/README_pathway_network.md` for usage
2. Review `PARQUET_FORMAT.md` for format issues
3. See `DEFAULT_BEHAVIOR.md` for parameter questions
4. Refer to `SKILL.md` for integration examples

---

**End of Addition Summary**
