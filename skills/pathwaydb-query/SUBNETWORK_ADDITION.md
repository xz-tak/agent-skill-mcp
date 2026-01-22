# Gene-Specific Subnetwork Analysis - Skill Addition Summary

**Date Added:** December 12, 2025
**Version:** 2.0
**Feature:** Gene-specific pathway subnetwork extraction and centrality analysis

---

## Overview

This update adds comprehensive **gene-specific subnetwork analysis capabilities** to the pathwaydb-query skill. After querying pathways for a gene, users can now extract the relevant subnetwork from the full pathway network and compute detailed centrality metrics to understand pathway importance, interconnectivity, and biological significance.

## What Was Added

### 1. Core Analysis Script

**File:** `network_builder/gene_subnetwork_analysis.py`

Extracts gene-specific pathway subnetworks and computes comprehensive centrality metrics:

**Key Features:**
- Loads full pathway network (630M+ edges)
- Extracts seed pathways from query results
- Optional neighbor expansion (--include-neighbors)
- Computes 9 centrality metrics:
  - Degree & Degree Centrality
  - Closeness Centrality (distance-based, using Distance = 1/Jaccard)
  - Betweenness Centrality (distance-based)
  - Eigenvector Centrality (similarity-based, using Jaccard)
  - PageRank (similarity-based)
  - Hub & Authority Scores (HITS algorithm)
  - Clustering Coefficient (similarity-based)
- Classifies nodes: Hub, Bridge, Leaf, Regular
- Exports to CSV and Excel with multiple sheets

**Usage:**
```bash
conda run -n claude_test python network_builder/gene_subnetwork_analysis.py \
  --gene TP53 \
  --network data/all_pathway_network_12052025.parquet \
  --pathways tp53_pathways.csv \
  --min-jaccard 0.05 \
  --output tp53_subnetwork
```

**Output Files:**
- `{output}_centrality.csv` - Full centrality metrics table
- `{output}_analysis.xlsx` - Multi-sheet Excel workbook (Centrality, Network_Stats, Top_Hubs, Top_Bridges, Seed_Pathways)
- `{output}_edgelist.csv` - Network edges for visualization (Cytoscape, Gephi)

### 2. Summarization and Interpretation Script

**File:** `network_builder/summarize_centrality_analysis.py`

Provides human-readable summaries with biological interpretations:

**Key Features:**
- Network structure summary (node counts, connectivity, database distribution)
- Top pathways by multiple centrality metrics
- Automated biological insights:
  - Network structure interpretation (highly interconnected vs sparse)
  - Bridge pathway identification
  - Seed pathway importance
  - Database coverage analysis
- Actionable recommendations for further analysis
- Formatted terminal output with clear sections

**Usage:**
```bash
conda run -n claude_test python network_builder/summarize_centrality_analysis.py \
  tp53_subnetwork_centrality.csv \
  --top 10
```

**Output:**
- Terminal-formatted report with:
  - Network structure statistics
  - Top pathways by degree, betweenness, PageRank, closeness
  - Biological insights with interpretation
  - Recommendations for next steps

### 3. Updated Documentation

**SKILL.md Updates:**
- Added "Gene-Specific Subnetwork Analysis" section with:
  - Overview and use case
  - Quick start guide
  - Detailed usage examples
  - Parameter descriptions
  - Output file explanations
  - Centrality metrics explained (with weight handling clarification)
  - Node classification definitions
  - Neighbor expansion details
  - Example results interpretation (with and without neighbors)
  - Integration workflow with pathway queries
  - Performance notes and tips
- Updated skill description to include subnetwork analysis
- Added "Gene-Specific Network Analysis" trigger phrases to "When to Use This Skill"
- Updated Key Capabilities to include centrality metrics
- Updated Resources section with new scripts

### 4. Example Files

**Location:** `examples/`

Test data and results demonstrating the workflow:
- `test_tp53_pathways.csv` - Example pathway list for TP53
- `tp53_no_neighbors_*` - Results without neighbor expansion
- Other test outputs

## Technical Implementation

### Weight Handling

**Critical Feature:** Proper handling of similarity vs distance metrics

- **Similarity-based metrics** (higher Jaccard = stronger connection):
  - Eigenvector Centrality: Uses `Jaccard_Index` as weight
  - PageRank: Uses `Jaccard_Index` as weight
  - Clustering Coefficient: Uses `Jaccard_Index` as weight

- **Distance-based metrics** (higher Jaccard = lower distance):
  - Closeness Centrality: Uses `Distance = 1/Jaccard_Index`
  - Betweenness Centrality: Uses `Distance = 1/Jaccard_Index`

This ensures metrics correctly interpret pathway similarity as connection strength.

### Node Classification

Automatic classification based on network position:

- **Hub**: Degree ≥ 50% of (average degree × network size)
  - Highly connected pathways, core biological functions

- **Bridge**: Betweenness Centrality > 0.1
  - Connects different pathway modules, critical for information flow

- **Leaf**: Degree = 1
  - Single connection, peripheral pathways

- **Regular**: Standard connectivity
  - Typical interconnection level

### Neighbor Expansion

**Default Behavior:** Neighbors NOT included (`--include-neighbors` disabled by default)

- Without neighbors: Analyze only gene's direct pathways
  - Example: TP53 in 3 pathways → 3 node subnetwork
  - Best for: Direct gene pathway involvement

- With neighbors: Include connected pathways
  - Example: TP53 in 3 pathways → 763 node subnetwork (with neighbors)
  - Best for: Broader functional context
  - Control distance: `--max-distance 1` (direct) or `2` (2-hop)

## Integration Points

### Workflow Integration

```
1. Query pathways for gene
   └─> scripts/pathway_query.py

2. Extract pathway names
   └─> Manual or programmatic extraction
   └─> Format: "pathway_name (database:id)"

3. Run subnetwork analysis
   └─> network_builder/gene_subnetwork_analysis.py

4. Summarize results
   └─> network_builder/summarize_centrality_analysis.py

5. Visualize
   └─> Import edgelist to Cytoscape/Gephi
   └─> Open Excel workbook for interactive exploration
```

### Existing Skill Integration

- Uses pre-computed network: `data/all_pathway_network_12052025.parquet`
- Accepts pathway query outputs from existing scripts
- Complements multi-gene analysis capabilities
- Extends pathway network analysis functionality

## Use Cases

### 1. Identify Most Important Pathways

**Question:** "What are the most central pathways for TP53?"

**Approach:**
```bash
# Run subnetwork analysis
python network_builder/gene_subnetwork_analysis.py \
  --gene TP53 --network data/all_pathway_network_12052025.parquet \
  --pathways tp53_pathways.csv --output tp53_analysis

# Summarize
python network_builder/summarize_centrality_analysis.py \
  tp53_analysis_centrality.csv --top 10
```

**Output:** Top 10 pathways by degree centrality, betweenness, PageRank

### 2. Find Bridge Pathways

**Question:** "Which pathways connect different TP53 functional modules?"

**Approach:** Look for pathways with high betweenness centrality (Node_Type = Bridge)

**Insight:** Bridge pathways are critical for coordinating multiple biological processes

### 3. Compare Direct vs Contextual Involvement

**Question:** "How do TP53's direct pathways fit into the broader pathway landscape?"

**Approach:**
```bash
# Without neighbors (direct involvement)
python network_builder/gene_subnetwork_analysis.py \
  --gene TP53 ... --output tp53_direct

# With neighbors (broader context)
python network_builder/gene_subnetwork_analysis.py \
  --gene TP53 ... --include-neighbors --output tp53_context

# Compare results
```

### 4. Visualize Pathway Relationships

**Question:** "How are TP53 pathways interconnected?"

**Approach:**
1. Generate edgelist: `tp53_subnetwork_edgelist.csv`
2. Import to Cytoscape
3. Style:
   - Node color by Database
   - Node size by Degree_Centrality
   - Node shape by Node_Type
   - Edge width by Jaccard_Index

## Performance Characteristics

| Operation | Time | Memory |
|-----------|------|--------|
| Load network (11GB parquet) | 30-120 sec | 10-15 GB |
| Load filtered (J>0.05) | 10-30 sec | 5-8 GB |
| Extract subnetwork | <1 sec | minimal |
| Compute centrality (small, <10 nodes) | <5 sec | minimal |
| Compute centrality (medium, 100-500 nodes) | 30-60 sec | ~1 GB |
| Compute centrality (large, 500+ nodes) | 1-3 min | ~2-4 GB |
| Summarization | <1 sec | minimal |

**Optimization Tips:**
- Use `--min-jaccard 0.05` or `0.1` to reduce network size
- Start without neighbors, add if needed
- Use parquet format (5-10x faster than CSV)

## Examples

### Example 1: Small Seed-Only Network

```bash
conda run -n claude_test python network_builder/gene_subnetwork_analysis.py \
  --gene TP53 \
  --network data/all_pathway_network_12052025.parquet \
  --pathways tp53_pathways.csv \
  --min-jaccard 0.1 \
  --output tp53_seed
```

**Result:**
```
Seed pathways: 3
Subnetwork nodes: 3
Subnetwork edges: 0
Network density: 0.0000
```

**Interpretation:** 3 TP53 pathways have no connections at J>0.1 threshold → Distinct pathway roles

### Example 2: Large Network with Neighbors

```bash
conda run -n claude_test python network_builder/gene_subnetwork_analysis.py \
  --gene TP53 \
  --network data/all_pathway_network_12052025.parquet \
  --pathways tp53_pathways.csv \
  --min-jaccard 0.05 \
  --include-neighbors \
  --output tp53_neighborhood
```

**Result:**
```
Seed pathways: 3
Subnetwork nodes: 763
Subnetwork edges: 58,854
Network density: 0.2025
Average degree: 154.27
Hubs: 762
```

**Interpretation:** TP53 pathways part of highly interconnected module (cell cycle/DNA repair)

### Example 3: Summarization Output

```bash
conda run -n claude_test python network_builder/summarize_centrality_analysis.py \
  tp53_neighborhood_centrality.csv --top 5
```

**Output Sections:**
1. **Network Structure:** Total pathways, node types, connectivity, database distribution
2. **Top Pathways:** By degree, betweenness, PageRank, closeness
3. **Seed Pathways:** Metrics for gene's direct pathways
4. **Biological Insights:** Automated interpretations
5. **Recommendations:** Next steps for analysis

## Biological Insights Generated

The summarization script automatically provides:

### Network Structure Insights
- **Highly interconnected**: >70% hubs → Integrated biological process, functional redundancy
- **Sparse network**: >50% leaves → Diverse, independent processes, pleiotropic effects

### Bridge Pathway Insights
- Identifies pathways connecting different modules
- Critical for information flow between processes
- Disruption affects multiple downstream pathways

### Seed Pathway Centrality
- High centrality seed pathways: Core functions with broad effects
- Low centrality: Specialized, specific roles

### Database Coverage
- Dominant database indicates curation focus
- Suggests complementary database queries

## Updated Files Summary

### New Files Created
```
network_builder/
├── gene_subnetwork_analysis.py           (16 KB)
└── summarize_centrality_analysis.py      (14 KB)

examples/
├── test_tp53_pathways.csv                (< 1 KB)
├── tp53_no_neighbors_centrality.csv      (< 1 KB)
├── tp53_no_neighbors_analysis.xlsx       (15 KB)
└── tp53_no_neighbors_edgelist.csv        (< 1 KB)

SUBNETWORK_ADDITION.md                    (This file, 15 KB)
```

### Modified Files
```
SKILL.md                                   (Updated, +240 lines)
  - Added "Gene-Specific Subnetwork Analysis" section
  - Updated description, capabilities, when to use
  - Updated Resources section
```

**Total Addition Size:** ~60 KB scripts + documentation + example files

## Future Enhancements

Potential improvements for future versions:

1. **Automated pathway name extraction** from query JSON outputs
2. **Comparative subnetwork analysis** (compare multiple genes)
3. **Pathway enrichment with network context** (enrichment weighted by centrality)
4. **Community detection** in subnetworks (identify functional modules)
5. **Time-series analysis** (how subnetwork changes with different thresholds)
6. **Interactive visualization** (web-based network browser)
7. **Statistical significance** testing for centrality differences
8. **Integration with gene expression data** (overlay expression on network)

## Citation

When using gene-specific subnetwork analysis, cite:

```
Gene-specific pathway subnetwork analysis using the pathwaydb-query skill
Date: December 12, 2025
Method: Jaccard similarity-based subnetwork extraction with comprehensive centrality metrics
Network: Pre-computed pathway similarity network (KEGG, MSigDB, Reactome)
  - 52,870 pathways, 630M+ edges
  - Generated: December 5, 2025
Tools: gene_subnetwork_analysis.py, summarize_centrality_analysis.py
Available at: .claude/skills/pathwaydb-query/
```

## Contact & Support

For questions or issues with subnetwork analysis:

1. Check `SKILL.md` Gene-Specific Subnetwork Analysis section
2. Review script help: `python gene_subnetwork_analysis.py --help`
3. See examples in `examples/` directory
4. Refer to `network_builder/README_pathway_network.md` for network details

---

**End of Subnetwork Analysis Addition Summary**
