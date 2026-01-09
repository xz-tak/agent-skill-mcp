# PrimeKG Shortest Path Analysis

## Overview

PrimeKG analysis uses graph traversal (BFS) to find shortest paths between genes and diseases in the PrimeKG knowledge graph. Unlike the neural approaches (BioBridge, ULTRA), PrimeKG provides explicit path information showing how entities connect through intermediate nodes.

**Key Characteristics:**
- Uses breadth-first search (BFS) for shortest path discovery
- Provides explicit paths with intermediate nodes (mechanistic insight)
- Score based on path length (exponential decay)
- Combo scores are averages of components (no synergy detection)
- Includes interactive network visualization

## When to Use PrimeKG

| Use Case | Recommendation |
|----------|----------------|
| Network topology analysis | ✓ Best choice |
| Finding mechanistic connections | ✓ Excellent - shows intermediate nodes |
| Direct vs indirect association | ✓ Path length indicates directness |
| Visualizing gene-disease networks | ✓ Built-in pyvis visualization |
| Synergy/dilution detection | Use BioBridge or ULTRA instead |
| Embedding similarity queries | Use BioBridge instead |
| Complex logical queries (e.g. intersection, projection) | Use ULTRA instead |

---

## Scoring Methodology

### Path Length Score Formula

PrimeKG uses exponential decay based on path length:

```
score = 0.9^(path_length - 1)
```

| Path Length | Hops | Score | Assessment | Interpretation |
|-------------|------|-------|------------|----------------|
| 1 | Direct | 1.00 | **Direct** | Curated direct association in KG |
| 2 | 1 intermediate | 0.90 | **Strong** | One intermediate node |
| 3 | 2 intermediates | 0.81 | **Moderate-Strong** | Two intermediate nodes |
| 4 | 3 intermediates | 0.729 | **Moderate** | Three intermediate nodes |
| 5 | 4 intermediates | 0.656 | Weak | Four intermediate nodes |
| 6 | 5 intermediates | 0.590 | Weak | Five intermediate nodes |
| 7+ | 6+ intermediates | < 0.53 | Very Weak | Distant connection |
| No path | - | 0.00 | **None** | No connection in graph |

### Why Exponential Decay?

- **Direct connections (1 hop)** represent curated associations and should score highest
- **Each additional hop** introduces uncertainty and potential noise
- **0.9 decay factor** balances penalizing distance while preserving signal for 2-4 hop paths
- **Score of 0** for disconnected entities provides clear signal

### Combo Score Calculation

For combos, PrimeKG uses **simple averaging** of component scores:

```
combo_score = mean(component_scores)
```

**Important**: Because averaging is used, PrimeKG combo analysis:
- Cannot detect synergy/dilution (delta is always ~0)
- Shows which component contributes most/least to combo score
- Provides path structure comparison across components

---

## Path Interpretation

### What Paths Tell You

The intermediate nodes in a path reveal **mechanistic connections**:

```
TYK2 → STAT3 → inflammatory bowel disease
```

This path suggests:
- TYK2 signals through STAT3
- STAT3 is directly associated with IBD
- TYK2's IBD link may be mediated by STAT3 signaling

### Common Intermediate Node Types

| Node Type | Interpretation |
|-----------|----------------|
| Gene/protein | Signaling pathway connection |
| Biological process | Shared functional involvement |
| Molecular function | Shared biochemical activity |
| Anatomy/tissue | Shared tissue expression or localization |
| Phenotype | Shared phenotypic manifestation |
| Drug | Shared therapeutic targeting |

### Path Quality Indicators

| Indicator | Good Sign | Concern |
|-----------|-----------|---------|
| Path length | 1-3 hops | 5+ hops |
| Intermediate types | Proteins, pathways | Anatomy only |
| Multiple paths | Converging evidence | Single path dependency |
| Consistent intermediates | Same hub across genes | Unrelated intermediates |

---

## Running PrimeKG Analysis

### Script Location

```
scripts/kg_association_shortest_path.py  (if implemented)
```

Or use the PrimeKG MCP tools directly.

### Using PrimeKG MCP Tools

PrimeKG provides several MCP tools for graph queries:

#### Search for Entities

```python
# Find entity IDs
result = mcp__primekg__search_entities(
    search_term="TYK2",
    entity_type="gene/protein",
    limit=5
)
# Returns: entity_id like "NCBI:7297"
```

#### Get Entity Connections

```python
# Get direct connections for an entity
result = mcp__primekg__get_entity_connections(
    entity_id="NCBI:7297",  # TYK2
    relation_type="associated with",
    limit=null  # IMPORTANT: Always use null to get ALL connections
)
```

#### Get Entity Neighborhood (Multi-hop)

```python
# Get 2-hop neighborhood
result = mcp__primekg__get_entity_neighborhoods(
    entity_id="NCBI:7297",  # TYK2
    max_depth=2,
    relation_types=["associated with", "ppi"],
    max_neighbors_per_level=10
)
```

---

## MCP Query Best Practices

### Always Use `limit=null` for Complete Data

**CRITICAL:** When querying PrimeKG for entity connections, **always set `limit=null`** to retrieve all connections. Limited queries may miss important associations.

```python
# CORRECT: Get all connections
result = mcp__primekg__get_entity_connections(
    entity_id="NCBI:7297",  # TYK2
    relation_type="associated with",
    limit=null  # No limit - get ALL connections
)
# Returns: total_connections and returned_connections will match

# INCORRECT: Limited results may miss important associations
result = mcp__primekg__get_entity_connections(
    entity_id="NCBI:7297",
    relation_type="associated with",
    limit=100  # May miss associations!
)
```

### Query All Entities Comprehensively

For complete analysis, query connections for:

1. **All target genes** - get complete disease associations
2. **All diseases** - get complete gene associations (for shared gene detection)
3. **Target PPI partners** - for shared PPI detection between combo targets

```python
# Query workflow for TYK2+JAK1 combo analysis on IBD diseases

# Step 1: Get ALL disease associations for each target
tyk2_diseases = mcp__primekg__get_entity_connections(
    entity_id="NCBI:7297", relation_type="associated with", limit=null
)
jak1_diseases = mcp__primekg__get_entity_connections(
    entity_id="NCBI:3716", relation_type="associated with", limit=null
)

# Step 2: Get ALL gene associations for each disease
ibd_genes = mcp__primekg__get_entity_connections(
    entity_id="MONDO_grouped:...", relation_type="associated with", limit=null
)
crohn_genes = mcp__primekg__get_entity_connections(
    entity_id="MONDO_grouped:5011_5535", relation_type="associated with", limit=null
)
uc_genes = mcp__primekg__get_entity_connections(
    entity_id="MONDO:5101", relation_type="associated with", limit=null
)

# Step 3: Get PPI partners for each target
tyk2_ppi = mcp__primekg__get_entity_connections(
    entity_id="NCBI:7297", relation_type="ppi", limit=null
)
jak1_ppi = mcp__primekg__get_entity_connections(
    entity_id="NCBI:3716", relation_type="ppi", limit=null
)
```

### Detecting Shared Associations

Identifying shared genes across diseases is critical for understanding common mechanisms:

```python
# Extract gene names from each disease query
ibd_gene_set = {c["target_name"] for c in ibd_genes["connections"] if c["direction"] == "outgoing"}
crohn_gene_set = {c["target_name"] for c in crohn_genes["connections"] if c["direction"] == "outgoing"}
uc_gene_set = {c["target_name"] for c in uc_genes["connections"] if c["direction"] == "outgoing"}

# Find shared genes
shared_all_three = ibd_gene_set & crohn_gene_set & uc_gene_set
shared_ibd_crohn = ibd_gene_set & crohn_gene_set - uc_gene_set
shared_ibd_uc = ibd_gene_set & uc_gene_set - crohn_gene_set
shared_crohn_uc = crohn_gene_set & uc_gene_set - ibd_gene_set

# Prioritize shared genes in visualization (they appear first)
print(f"Shared across all 3 diseases: {shared_all_three}")
# Example: {'IL10', 'IL12B', 'HLA-DRB1'}
```

### Detecting Shared PPI Partners

For combo analysis, find genes that interact with ALL targets:

```python
# Extract PPI partners
tyk2_partners = {c["target_name"] for c in tyk2_ppi["connections"]}
jak1_partners = {c["target_name"] for c in jak1_ppi["connections"]}

# Shared PPI = interacts with BOTH targets
shared_ppi = tyk2_partners & jak1_partners
# Example: {'HSP90AB1', 'PDGFRB', 'MAS1', 'PRMT5', 'GHR'}

# Unique PPI = interacts with only ONE target
tyk2_unique_ppi = tyk2_partners - jak1_partners
jak1_unique_ppi = jak1_partners - tyk2_partners
```

### Visualization Priority System

When building visualizations, nodes should be prioritized in this order:

| Priority | Category | Description | Color |
|----------|----------|-------------|-------|
| 100 | Target genes | Your input targets (e.g., TYK2, JAK1) | Blue |
| 95 | Diseases | Target diseases (IBD, Crohn, UC) | Red |
| 90 | Target-target PPI | Direct interactions between targets | Red edge (thick) |
| 85 | Shared intermediary | Nodes in paths for multiple targets | Gray |
| 80-85 | Shared disease genes | Genes associated with 2+ diseases | Orange (+5 boost) |
| 70 | Shared PPI partners | Genes interacting with all targets | Purple |
| 60 | Unique PPI partners | Target-specific interaction partners | Green |

### JSON Structure for Complete Visualization

Include complete disease gene data for shared gene detection:

```json
{
  "genes": ["TYK2", "JAK1"],
  "diseases": ["inflammatory bowel disease", "Crohn disease", "ulcerative colitis"],
  "individual": { ... },
  "combo": { ... },

  "ppi_data": {
    "shared_ppi": ["HSP90AB1", "PDGFRB", "MAS1", "PRMT5", "GHR"],
    "target_ppi": {
      "TYK2": ["STAT1", "IRF3", "NFKB1", "LYN", "JAK3"],
      "JAK1": ["SOCS3", "EGFR", "IL2RG", "IFNAR1", "IFNAR2"]
    }
  },

  "disease_genes": {
    "inflammatory bowel disease": ["IL10", "IL12B", "HLA-DRB1", "IL1B", "CCR6", ...],
    "Crohn disease": ["IL10", "IL12B", "HLA-DRB1", "JAK2", "IFNG", ...],
    "ulcerative colitis": ["IL10", "IL12B", "HLA-DRB1", "STAT3", "JAK2", ...]
  }
}
```

The visualization script will automatically:
- Detect genes shared across 2+ diseases (IL10, IL12B, HLA-DRB1)
- Boost their priority (+5)
- Connect them to ALL associated disease nodes
- Display "(Shared across diseases)" in node tooltips

### Manual BFS Shortest Path

For shortest path analysis, implement BFS traversal:

```python
from collections import deque

def bfs_shortest_path(graph, source_id, target_id, max_depth=10):
    """
    Find shortest path using BFS.

    Args:
        graph: Adjacency list {node_id: [(neighbor_id, relation), ...]}
        source_id: Starting node (e.g., gene ID)
        target_id: Target node (e.g., disease ID)
        max_depth: Maximum path length to search

    Returns:
        (path, length) or None if no path found
    """
    if source_id not in graph or target_id not in graph:
        return None

    queue = deque([(source_id, [source_id])])
    visited = {source_id}

    while queue:
        current, path = queue.popleft()

        if len(path) > max_depth:
            continue

        if current == target_id:
            return (path, len(path) - 1)

        for neighbor, relation in graph.get(current, []):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))

    return None

def calculate_score(path_length):
    """Calculate score from path length using exponential decay."""
    if path_length is None or path_length <= 0:
        return 0.0
    return 0.9 ** (path_length - 1)
```

---

## Network Visualization

### Visualization Script

**Script**: `scripts/primekg_visualization.py`

Creates interactive network visualizations using pyvis with:
- Color-coded node types
- Interactive drag and zoom
- Path highlighting
- Legend with explanations

**IMPORTANT:** Visualization generation is REQUIRED for all PrimeKG analyses.

### Required Visualizations

For every kg-association analysis, generate:

1. **Full Network Visualization** - Shows all genes and diseases with paths
2. **Combo Subgraph(s)** - One per combo, showing only that combo's paths

### Usage Examples

**STEP 1: Generate Full Network (REQUIRED)**

```bash
# Format: {genes}_{disease}_network.html
pixi run python scripts/primekg_visualization.py \
    --json-results ./kgpred_<context>/primekg/primekg_shortest_paths.json \
    --output ./kgpred_<context>/primekg/{genes}_{disease}_network.html \
    --title "<Genes> <Disease> Network (PrimeKG)"
```

Example:
```bash
pixi run python scripts/primekg_visualization.py \
    --json-results ./kgpred_IBD_bispecific/primekg/primekg_shortest_paths.json \
    --output ./kgpred_IBD_bispecific/primekg/TYK2_JAK1_IBD_network.html \
    --title "TYK2+JAK1 IBD Network (PrimeKG)"
```

**STEP 2: Generate Combo Subgraph(s) (REQUIRED for each combo)**

```bash
# Format: {genes}_{disease}_subgraph.html
pixi run python scripts/primekg_visualization.py \
    --json-results ./kgpred_<context>/primekg/primekg_shortest_paths.json \
    --combo "TYK2+JAK1" \
    --output ./kgpred_<context>/primekg/TYK2_JAK1_IBD_subgraph.html
```

**Direct parameters (without JSON file):**

```bash
pixi run python scripts/primekg_visualization.py \
    --genes "TYK2,JAK1" \
    --diseases "inflammatory bowel disease" \
    --paths '{"TYK2": {"IBD": {"path": ["TYK2", "STAT3", "IBD"], "score": 0.9, "path_found": true}}}' \
    --output ./visualization.html
```

### Visualization Workflow

```
1. Complete PrimeKG MCP queries
      ↓
2. Save results to primekg/primekg_shortest_paths.json
      ↓
3. Generate full network visualization
      ↓
4. Generate combo subgraph for each combo
      ↓
5. Reference visualizations in report
```

### Python API

```python
import sys
sys.path.append('/path/to/kg-association/scripts')
from primekg_visualization import create_network_visualization, create_combo_visualization

# Create main visualization
create_network_visualization(
    genes=["TYK2", "JAK1", "ITGA4"],
    diseases=["inflammatory bowel disease", "Crohn disease"],
    paths_data=individual_results,
    title="IBD Gene-Disease Network",
    output_path="./network.html"
)

# Create combo-specific visualization
create_combo_visualization(
    combo=["TYK2", "JAK1"],
    diseases=["inflammatory bowel disease"],
    paths_data=individual_results,
    output_path="./tyk2_jak1_subgraph.html"
)
```

### Visualization Color Scheme

**Node Colors:**

| Node Type | Color | Shape |
|-----------|-------|-------|
| Disease | Red (#e74c3c) | Diamond |
| Target Gene | Blue (#3498db) | Circle |
| Intermediate Node | Green (#2ecc71) | Circle |
| GO Term / Pathway | Purple (#9b59b6) | Circle |

**Edge Colors:**

| Edge Type | Color |
|-----------|-------|
| From Gene | Blue (#3498db) |
| Path Connection | Green (#2ecc71) |
| To Disease | Red (#e74c3c) |

### Required Input Structure for Visualization

```json
{
  "genes": ["TYK2", "JAK1", "ITGA4"],
  "diseases": ["inflammatory bowel disease", "Crohn disease"],
  "individual": {
    "TYK2": {
      "inflammatory bowel disease": {
        "path_found": true,
        "path": ["TYK2", "STAT3", "inflammatory bowel disease"],
        "path_details": [
          {"id": "NCBI:7297", "name": "TYK2", "type": "gene/protein"},
          {"id": "NCBI:6774", "name": "STAT3", "type": "gene/protein"},
          {"id": "MONDO:0005265", "name": "inflammatory bowel disease", "type": "disease"}
        ],
        "score": 0.9,
        "path_length": 2
      }
    }
  }
}
```

**Required fields:**
- `genes`: List of gene names
- `diseases`: List of disease names
- `individual`: Dict mapping gene → disease → path_result

**Required per path_result:**
- `path_found`: Boolean
- `path`: List of node names (strings)
- `score`: Float (0-1)

**Optional per path_result:**
- `path_details`: List of {id, name, type} dicts (enables type-based coloring)
- `path_length`: Integer

---

## Report Formats

### Part 1: Individual Analysis

```markdown
## PrimeKG Individual Gene-Disease Analysis

### Inflammatory Bowel Disease

| Gene | Path Length | Score | Path | Assessment |
|------|-------------|-------|------|------------|
| TYK2 | 2 | 0.90 | TYK2 → STAT3 → IBD | Strong |
| JAK1 | 2 | 0.90 | JAK1 → STAT1 → IBD | Strong |
| ITGA4 | 1 | 1.00 | ITGA4 → IBD | Direct |
| IL17A | 3 | 0.81 | IL17A → IL17R → NF-κB → IBD | Moderate-Strong |

**Key Findings:**
- ITGA4 has direct association (1 hop) - curated link in PrimeKG
- TYK2 and JAK1 connect through STAT proteins (shared JAK-STAT pathway)
- IL17A requires 3 hops, suggesting more indirect mechanism
```

### Part 2: Combo Analysis

```markdown
## PrimeKG Combo Analysis

### Inflammatory Bowel Disease

| Combo | Combo Score | Component Scores | Max | Min | Classification |
|-------|-------------|------------------|-----|-----|----------------|
| TYK2+JAK1 | 0.90 | TYK2: 0.90, JAK1: 0.90 | 0.90 | 0.90 | NEAR-ADDITIVE |
| ITGA4+ITGB7 | 0.95 | ITGA4: 1.00, ITGB7: 0.90 | 1.00 | 0.90 | NEAR-ADDITIVE |

**Notes:**
- PrimeKG combo scores = average of component scores
- Delta is always ~0 (no synergy detection with averaging)
- Max/Min show score range across components
```

### Part 3: Comparative Analysis

```markdown
## PrimeKG Part 3: Comparative Analysis

### TYK2+JAK1 Path Structure Analysis

#### Individual Paths
| Gene | Path | Length | Score |
|------|------|--------|-------|
| TYK2 | TYK2 → STAT3 → IBD | 2 | 0.90 |
| JAK1 | JAK1 → STAT1 → IBD | 2 | 0.90 |

#### Shared Intermediate Analysis
- Both genes connect to IBD through STAT proteins
- STAT3 and STAT1 are both transcription factors in JAK-STAT pathway
- Convergent biology: both genes activate similar downstream signaling

#### Path Structure Insights
1. **Equal path lengths**: Both genes are similarly "distant" from IBD in the graph
2. **Related intermediates**: STAT proteins suggest mechanistic similarity
3. **No direct connections**: Neither gene has curated direct IBD association

#### Mechanistic Interpretation
TYK2 and JAK1 appear to influence IBD through overlapping JAK-STAT signaling:
- TYK2 → STAT3 activation → inflammatory gene expression
- JAK1 → STAT1 activation → interferon response
- Both pathways converge on inflammatory processes in IBD
```

---

## Output Files

### Data Files (in `primekg/`)

| File | Contents |
|------|----------|
| `primekg_shortest_paths.json` | Complete results with paths, scores, metadata |

### Markdown Reports (in `primekg/`)

| File | Generated When | Contents |
|------|----------------|----------|
| `Part1_Individual_Analysis.md` | Always | Individual gene-disease paths and scores |
| `Part2_Combo_Analysis.md` | If combos provided | Combo average scores |
| `Part3_Comparative_Analysis.md` | If combos provided | Path structure comparison |

### Visualizations (in `primekg/`)

| File | Generated When | Contents |
|------|----------------|----------|
| `{genes}_{disease}_network.html` | Always (REQUIRED) | Full gene-disease network |
| `{genes}_{disease}_subgraph.html` | For each combo (REQUIRED) | Combo-specific subgraph |

**Naming Convention:**
- `{genes}`: Target gene names joined by underscore (e.g., `TYK2_JAK1`)
- `{disease}`: Disease context abbreviation (e.g., `IBD`, `Crohn`)

**Examples:**
- `TYK2_JAK1_IBD_network.html` - Full network for TYK2+JAK1 IBD analysis
- `TYK2_JAK1_IBD_subgraph.html` - Combo subgraph for TYK2+JAK1

---

## JSON Output Structure

**IMPORTANT:** Always include `pct_rank` field for consistency with BioBridge and ULTRA.
For PrimeKG, `pct_rank = score` (the exponential decay score serves as the percentile rank equivalent).

```json
{
  "genes": ["TYK2", "JAK1"],
  "diseases": ["inflammatory bowel disease", "Crohn disease"],
  "combos": ["TYK2+JAK1"],
  "individual": {
    "TYK2": {
      "inflammatory bowel disease": {
        "path_found": true,
        "path": ["TYK2", "inflammatory bowel disease"],
        "path_details": [
          {"id": "NCBI:7297", "name": "TYK2", "type": "gene/protein"},
          {"id": "MONDO:0005265", "name": "inflammatory bowel disease", "type": "disease"}
        ],
        "path_length": 1,
        "score": 1.0,
        "pct_rank": 1.0,
        "source_id": "NCBI:7297",
        "target_id": "MONDO:0005265"
      }
    },
    "JAK1": {
      "inflammatory bowel disease": {
        "path_found": true,
        "path": ["JAK1", "TYK2", "inflammatory bowel disease"],
        "path_details": [
          {"id": "NCBI:3716", "name": "JAK1", "type": "gene/protein"},
          {"id": "NCBI:7297", "name": "TYK2", "type": "gene/protein"},
          {"id": "MONDO:0005265", "name": "inflammatory bowel disease", "type": "disease"}
        ],
        "path_length": 2,
        "score": 0.9,
        "pct_rank": 0.9,
        "source_id": "NCBI:3716",
        "target_id": "MONDO:0005265",
        "intermediate": "TYK2"
      }
    }
  },
  "combo": {
    "TYK2+JAK1": {
      "inflammatory bowel disease": {
        "combo_score": 0.95,
        "combo_pct_rank": 0.95,
        "individual_mean": 0.95,
        "delta": 0.0,
        "classification": "NEAR-ADDITIVE",
        "component_scores": {
          "TYK2": {"score": 1.0, "pct_rank": 1.0, "path_length": 1, "path": ["TYK2", "IBD"]},
          "JAK1": {"score": 0.9, "pct_rank": 0.9, "path_length": 2, "path": ["JAK1", "TYK2", "IBD"]}
        }
      }
    }
  },
  "summary": {
    "total_pairs": 4,
    "paths_found": 4,
    "avg_score": 0.95,
    "avg_path_length": 1.5
  }
}
```

**Required fields for cross-method reporting:**
- `pct_rank` - Always include (= score for PrimeKG)
- `combo_pct_rank` - For combo results

### Extended JSON Structure for Enhanced Visualization

To generate rich visualizations with PPI partners and disease-associated genes, include these additional fields:

```json
{
  "genes": ["TYK2", "JAK1"],
  "diseases": ["inflammatory bowel disease"],
  "individual": { ... },
  "combo": { ... },

  "ppi_data": {
    "shared_ppi": ["STAT1", "STAT3", "HSP90AB1", "SOCS1"],
    "target_ppi": {
      "TYK2": ["IFNAR1", "IFNAR2", "IL10RA", "IRF3"],
      "JAK1": ["IL2RG", "IL6R", "EGFR", "BRCA1"]
    }
  },

  "disease_genes": {
    "inflammatory bowel disease": ["NOD2", "IL23R", "ATG16L1", "CARD9"],
    "Crohn disease": ["NOD2", "IL23R", "IRGM", "ATG16L1"]
  }
}
```

**Field Descriptions:**

| Field | Description | Source |
|-------|-------------|--------|
| `ppi_data.shared_ppi` | Genes with PPI to ALL combo targets | PrimeKG `ppi` relation intersection |
| `ppi_data.target_ppi` | Genes with PPI to specific target only | PrimeKG `ppi` relation minus shared |
| `disease_genes` | Genes associated with each disease | PrimeKG `associated with` relation |

**How to Populate from MCP Queries:**

```python
# 1. Get PPI partners for each target
tyk2_ppi = mcp__primekg__get_entity_connections(
    entity_id="NCBI:7297", relation_type="ppi", limit=100
)
jak1_ppi = mcp__primekg__get_entity_connections(
    entity_id="NCBI:3716", relation_type="ppi", limit=100
)

# 2. Find shared PPI (intersection)
tyk2_partners = {c["target_name"] for c in tyk2_ppi["connections"]}
jak1_partners = {c["target_name"] for c in jak1_ppi["connections"]}
shared_ppi = list(tyk2_partners & jak1_partners)

# 3. Find unique PPI per target
tyk2_unique = list(tyk2_partners - jak1_partners)
jak1_unique = list(jak1_partners - tyk2_partners)

# 4. Get disease-associated genes
ibd_genes = mcp__primekg__get_entity_connections(
    entity_id="MONDO:5265", relation_type="associated with", limit=50
)

# 5. Build ppi_data structure
ppi_data = {
    "shared_ppi": shared_ppi[:10],  # Top 10
    "target_ppi": {
        "TYK2": tyk2_unique[:8],
        "JAK1": jak1_unique[:8]
    }
}
```

---

## Comparison with BioBridge and ULTRA

| Aspect | PrimeKG | BioBridge | ULTRA |
|--------|---------|-----------|-------|
| **Method** | Graph traversal (BFS) | Embedding similarity | Foundation model inference |
| **Score basis** | Path length | Cosine similarity | Model prediction |
| **Combo method** | Average of components | Mean embedding | Logical intersection |
| **Synergy detection** | No (delta always ~0) | Yes | Yes |
| **Mechanistic insight** | Yes (intermediate nodes) | Limited | Limited |
| **Path information** | Full paths available | No paths | No paths |
| **Visualization** | Built-in pyvis | Not included | Not included |

### When to Combine Methods

For highest confidence, run all three methods and look for **agreement**:

| Agreement | Confidence |
|-----------|------------|
| All 3 methods show strong association | **HIGH** |
| 2 of 3 methods agree | **MODERATE** |
| All 3 disagree | **LOW** - investigate further |

PrimeKG's unique value:
- Provides **mechanistic interpretation** through intermediate nodes
- Identifies **direct associations** (1 hop = curated link)
- Shows **network topology** through visualization

---

## Troubleshooting

### Common Issues

**Issue: No path found for expected association**
- Gene/disease name may not match PrimeKG naming (try aliases)
- Association may not be in PrimeKG (use BioBridge/ULTRA as fallback)
- Check max_depth parameter (increase if needed)

**Issue: Very long paths (5+ hops)**
- Association is likely indirect/weak
- Consider if the intermediate nodes make biological sense
- Score will be low (<0.66) - treat with caution

**Issue: Visualization too crowded**
- Use combo-specific visualizations for focused view
- Reduce number of genes/diseases displayed
- Filter to specific diseases

**Issue: Unexpected intermediate nodes**
- PrimeKG includes many entity types (anatomy, phenotypes, etc.)
- Some paths may go through non-intuitive intermediates
- Cross-reference with BioBridge/ULTRA for validation

### Requirements

- **pyvis**: `pip install pyvis` (for visualization)
- **PrimeKG data**: `primekg.csv` and node files in data directory
- **polars or pandas**: For data processing

### Validation Checklist

- [ ] Gene and disease names resolve correctly in PrimeKG
- [ ] Path length scoring uses correct formula (0.9^(n-1))
- [ ] Combo scores are averages (not attempting synergy detection)
- [ ] Visualization includes all path nodes
- [ ] JSON output includes both path and path_details