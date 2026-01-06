---
name: interactdb-query
description: Query protein-protein interaction databases (STRING, IntAct, BioGRID) for single-gene neighbor analysis and multi-gene shortest path finding with comprehensive filtering options. This skill should be used when analyzing protein interactions, finding network connections between genes, conducting systems biology research requiring interaction data with entity and edge annotations, validating pathway connections, performing drug target discovery, or preparing interaction networks for visualization.
---

# Interaction Database Query

Query three major protein-protein interaction databases (STRING, IntAct, BioGRID) with support for single-gene neighbor queries and multi-gene shortest path analysis. All queries return comprehensive annotations including confidence scores, experimental methods, publications, and path tracking.

## Core Capabilities

1. **Single-Gene Queries**: Find all interaction neighbors of a target gene with confidence scores and evidence types
2. **Multi-Gene Shortest Paths**: Discover how multiple genes connect through the network using Dijkstra's algorithm
3. **Multi-Hop Expansion**: Automatically expand to 2-hop and 3-hop neighbors to reach desired result counts (NEW in v1.1)
4. **Unified Interface**: Query all three databases with a single function call

## Important Notes

**BioGRID API Key Requirement:**
- BioGRID queries require an API key set in environment variable `BIOGRID_API_KEY`
- Register at https://thebiogrid.org/ to obtain a free API key
- STRING and IntAct do not require API keys

**Default Output Location:**
- Results are saved to the current working directory by default
- Specify `output_dir` parameter to save to a different location

## Primary Workflow: Unified Query Interface

To query all three databases at once, use `scripts/unified_query.py`:

```python
from scripts.unified_query import query_single_gene_all_databases

# Query single gene across all databases
results = query_single_gene_all_databases(
    gene="TP53",
    species=9606,           # 9606=human, 10090=mouse
    top_n=100,
    export_results=True,
    output_dir=None,        # Optional: None uses current working directory
    min_combined_score=400  # Optional: filter threshold
)

# Results structure:
# results['string']  -> DataFrame of STRING neighbors
# results['intact']  -> DataFrame of IntAct neighbors
# results['biogrid'] -> List of BioGRID NeighborRecord objects
# results['exports'] -> Dict of exported file paths
```

To find shortest paths between multiple genes:

```python
from scripts.unified_query import query_shortest_paths_all_databases

# Find paths between genes across all databases
paths = query_shortest_paths_all_databases(
    gene_list=["TP53", "MDM2", "ATM"],
    species=9606,
    max_distance=3,
    export_results=True,
    output_dir=None,        # Optional: None uses current working directory
    min_combined_score=400
)

# Results structure:
# paths['string']  -> Dict mapping (gene_a, gene_b) tuples to path info
# paths['intact']  -> Dict mapping (gene_a, gene_b) tuples to path info
# paths['biogrid'] -> Dict mapping (gene_a, gene_b) tuples to path info
# paths['exports'] -> Dict of exported JSON file paths
```

**When to use unified interface**:
- Comparative analysis across databases
- Cross-validation of findings
- Comprehensive network analysis
- When BioGRID API key is configured

## Database-Specific Workflows

### STRING Database (Fast, No API Key Required)

**Use when**: Need fast, comprehensive results with multiple evidence types

```python
from scripts.string_api import get_string_neighbors, StringClient

# Single-gene query with evidence filtering
df = get_string_neighbors(
    gene="TP53",
    species=9606,
    top_n=50,
    min_score=700,              # Combined confidence (0-999)
    min_experimental=800,       # Experimental evidence threshold
    network_type="functional"   # or "physical"
)
df.to_csv("tp53_neighbors.csv", index=False)

# Shortest paths between multiple genes
client = StringClient()
paths = client.find_shortest_paths(
    gene_list=["TP53", "MDM2", "ATM"],
    species=9606,
    max_distance=3,
    min_combined_score=700,
    min_experimental_score=500  # Optional: require experimental evidence
)

# Display results
for (gene_a, gene_b), info in paths.items():
    path_str = " → ".join(info['path'])
    print(f"{gene_a} ↔ {gene_b}: {path_str}")
    print(f"  Hops: {info['hops']}, Distance: {info['distance']:.2f}")
    print(f"  Edge scores: {info['scores']}")
```

**STRING Evidence Channels** (all 0-999):
- `combined_score`: Overall confidence (≥700 recommended)
- `experimental_score`: Lab experiments (≥800 for high confidence)
- `database_score`: Curated databases (≥500 recommended)
- `textmining_score`: Literature co-mentions
- `coexpression_score`: Gene expression patterns (≥600 recommended)

See `references/database_comparison.md` for complete evidence channel descriptions and threshold recommendations.

### IntAct Database (Detailed Curation, REST API)

**Use when**: Need detailed experimental methods, PSI-MI ontology annotations, or cross-species filtering

**NEW in v1.2**: IntAct now uses REST API instead of PSICQUIC for improved reliability and performance:
- **25% faster** query times (~3s vs ~4s)
- **95% success rate** (vs 50% with PSICQUIC)
- **Simpler JSON parsing** (vs MITAB format)
- **Automatic UniProt AC mapping** for precise gene searches

```python
from scripts.intact_api import get_direct_neighbors, find_shortest_paths_intact

# Single-gene query with organism filtering
df = get_direct_neighbors(
    gene="TP53",
    species="human",  # Automatically converts to taxid "9606"
    top_n=100,
    organism_filter="homo sapiens,mus musculus"  # Optional: multi-species
)
df.to_csv("tp53_intact.csv", index=False)

# Shortest paths with confidence filtering
paths = find_shortest_paths_intact(
    gene_list=["TP53", "MDM2", "EP300"],
    species="human",  # Automatically converts to taxid "9606"
    max_distance=3,
    min_miscore=0.7,           # MI-score threshold (0.0-1.0)
    organism_filter="homo sapiens"
)

for (gene_a, gene_b), info in paths.items():
    print(f"{gene_a} → {gene_b}: {' → '.join(info['path'])}")
    print(f"  MI-scores: {[f'{s:.3f}' for s in info['scores']]}")
    print(f"  Algorithm: {info['algorithm']}")  # Shows "Dijkstra"
    print(f"  Weight formula: {info['weight_formula']}")  # Shows "weight = 1.0 - miscore"
```

**IntAct Confidence Levels** (MI-score):
- 0.4-0.69: Moderate confidence (default threshold)
- 0.7-0.89: High confidence
- 0.9-1.0: Very high confidence

**IntAct REST API Features**:
- **UniProt AC mapping**: Gene symbols automatically converted to UniProt accessions (e.g., TP53 → P04637)
- **Precise searches**: UniProt AC queries return exact gene interactions (not substring matches)
- **Prioritized entries**: Swiss-Prot (reviewed) entries prioritized over TrEMBL (unreviewed)
- **Weighted shortest paths**: Edge weight = 1.0 - miscore (higher MI-score = shorter distance)
- **JSON responses**: Direct field access (`intactMiscore`) instead of MITAB parsing

**Backward Compatibility**:
- Function signatures unchanged (same parameters as PSICQUIC version)
- Automatic species conversion (e.g., "human" → "9606", "mouse" → "10090")
- Old PSICQUIC version backed up as `intact_api_psicquic_backup.py`

### BioGRID Database (Genetic Interactions, Multi-Hop)

**Use when**: Need genetic interactions, manual curation focus, or multi-hop expansion for sparse genes

**NEW in v1.1**: BioGRID now supports multi-hop expansion (1→2→3 hops) with 5-minute default timeout.

```python
from scripts.biogrid_api import BioGRIDClient
import os

# Initialize with API key (register at https://thebiogrid.org/)
api_key = "d3367ed24eeea8fe8718f4993ed63ec9"
client = BioGRIDClient(api_key)  # Default timeout: 5 minutes

# Multi-hop expansion for sparse genes
# Example: GREM1 has only 4 direct neighbors
# With max_hops=3, automatically expands to find 30 total
neighbors = client.get_neighbors(
    seed_gene="GREM1",
    tax_id="9606",
    max_hops=3,              # Expand to 2-hop, 3-hop automatically
    max_neighbors=30,
    experimental_system_types=["physical"]  # or ["genetic"] or both
)

# Convert to DataFrame for export
import pandas as pd
from dataclasses import asdict
df = pd.DataFrame([asdict(n) for n in neighbors])
df.to_csv("grem1_biogrid.csv", index=False)

# Shortest paths
paths = client.find_shortest_paths(
    gene_list=["TP53", "MDM2", "ATM"],
    tax_id="9606",
    max_distance=3,
    min_score=0.5
)
```

**BioGRID Multi-Hop Performance** (NEW in v1.1):
- 1-hop queries: ~5-10 seconds (unchanged)
- Multi-hop queries: ~1-5 minutes depending on gene connectivity
- Default timeout: 5 minutes (300 seconds)
- Best for: Genes with sparse direct interactions (<20 neighbors)

See `references/CHANGELOG.md` for Version 1.1.0 multi-hop enhancement details.

## Common Use Cases

### Use Case 1: Drug Target Discovery

Find how a drug target connects to disease genes and identify bridge proteins as potential co-targets:

```python
from scripts.string_api import StringClient

target = "EGFR"
disease_genes = ["BRCA1", "BRCA2", "TP53", "ATM"]

client = StringClient()
paths = client.find_shortest_paths(
    gene_list=[target] + disease_genes,
    species=9606,
    max_distance=3,
    min_combined_score=700,
    min_experimental_score=800  # Require experimental evidence
)

# Identify bridge proteins
bridge_proteins = set()
for (gene_a, gene_b), info in paths.items():
    if gene_a == target and gene_b in disease_genes:
        bridges = info['path'][1:-1]  # Intermediate nodes
        if bridges:
            bridge_proteins.update(bridges)
            print(f"{target} → {gene_b} via {', '.join(bridges)}")

print(f"\nPotential co-targets: {', '.join(bridge_proteins)}")
```

### Use Case 2: Pathway Validation

Verify expected pathway connections and identify missing links:

```python
from scripts.string_api import get_string_neighbors, StringClient
from itertools import combinations

# Step 1: Get neighbors of pathway seed
pathway_seed = "TP53"
neighbors_df = get_string_neighbors(
    gene=pathway_seed,
    species=9606,
    top_n=50,
    min_score=700,
    min_experimental=800
)

# Step 2: Check connectivity among pathway members
pathway_genes = ["TP53", "MDM2", "ATM", "CHEK2", "BRCA1"]
client = StringClient()
paths = client.find_shortest_paths(
    gene_list=pathway_genes,
    species=9606,
    max_distance=2,
    min_combined_score=700
)

# Step 3: Categorize connections
direct = []
indirect = []
missing = []

for pair in combinations(pathway_genes, 2):
    key = tuple(sorted(pair))
    if key in paths:
        if paths[key]['hops'] == 1:
            direct.append(pair)
        else:
            indirect.append((pair, paths[key]['path'][1:-1]))
    else:
        missing.append(pair)

print(f"Direct connections: {len(direct)}")
print(f"Indirect connections: {len(indirect)} (with bridge proteins)")
print(f"Missing connections: {len(missing)}")
```

### Use Case 3: Cross-Database Consensus

Find paths supported by multiple databases for high-confidence results:

```python
from scripts.unified_query import query_shortest_paths_all_databases

genes = ["TP53", "MDM2", "ATM", "CHEK2"]

# Query all databases
results = query_shortest_paths_all_databases(
    gene_list=genes,
    species=9606,
    max_distance=3,
    min_combined_score=700,  # STRING threshold
    min_miscore=0.7,         # IntAct threshold
    min_score=0.7            # BioGRID threshold
)

# Find consensus paths (present in 2+ databases)
all_pairs = set(list(results['string'].keys()) +
                list(results['intact'].keys()) +
                list(results['biogrid'].keys()))

consensus = {}
for pair in all_pairs:
    db_count = sum([
        pair in results['string'],
        pair in results['intact'],
        pair in results['biogrid']
    ])
    if db_count >= 2:
        consensus[pair] = {
            'string': results['string'].get(pair),
            'intact': results['intact'].get(pair),
            'biogrid': results['biogrid'].get(pair)
        }

print(f"Consensus paths (≥2 databases): {len(consensus)}")
```

### Use Case 4: Network Visualization Export

Prepare data for Cytoscape, NetworkX, or other visualization tools:

```python
from scripts.string_api import StringClient
import pandas as pd

genes = ["TP53", "MDM2", "ATM", "CHEK2", "BRCA1"]
client = StringClient()
paths = client.find_shortest_paths(
    gene_list=genes,
    species=9606,
    max_distance=3,
    min_combined_score=700
)

# Extract unique edges
edges = []
for (gene_a, gene_b), info in paths.items():
    for i in range(len(info['path']) - 1):
        edges.append({
            'source': info['path'][i],
            'target': info['path'][i+1],
            'score': info['scores'][i],
            'weight': 1000 - info['scores'][i],  # For layout algorithms
            'query_pair': f"{gene_a}-{gene_b}"
        })

# Remove duplicates (keep highest score)
edge_df = pd.DataFrame(edges)
edge_df = edge_df.sort_values('score', ascending=False).drop_duplicates(
    subset=['source', 'target'], keep='first'
)
edge_df.to_csv("network_edges.csv", index=False)

# Create node list with attributes
nodes = set()
for edge in edges:
    nodes.add(edge['source'])
    nodes.add(edge['target'])

node_df = pd.DataFrame({
    'node': list(nodes),
    'is_query': [n in genes for n in nodes]
})
node_df.to_csv("network_nodes.csv", index=False)

print(f"Exported {len(edge_df)} edges and {len(node_df)} nodes")
print(f"Query genes: {sum(node_df['is_query'])}, Bridge proteins: {sum(~node_df['is_query'])}")
```

## Cross-Database Comparative Analysis

**Recommendation**: Always query multiple databases to maximize confidence and discover complementary information. Different databases excel at different interaction types and experimental methods.

### Comprehensive Analysis Workflow

Query all three databases and compare results to identify:
1. **Consensus interactions** - present in 2+ databases (high confidence)
2. **Database-specific interactions** - unique to one database (validate carefully)
3. **Complementary evidence** - different experimental methods supporting same interaction

```python
from scripts.unified_query import query_single_gene_all_databases
import pandas as pd

# Query all databases
results = query_single_gene_all_databases(
    gene="TP53",
    species=9606,
    top_n=100,
    export_results=True,
    output_dir="./results"
)

# Generate cross-database summary
summary = {
    'STRING': len(results['string']) if isinstance(results['string'], pd.DataFrame) else 0,
    'IntAct': len(results['intact']) if isinstance(results['intact'], pd.DataFrame) else 0,
    'BioGRID': len(results['biogrid']) if isinstance(results['biogrid'], list) else 0
}

print(f"\n{'='*70}")
print(f"CROSS-DATABASE SUMMARY: TP53 Interaction Network")
print(f"{'='*70}")
print(f"STRING:  {summary['STRING']:>4} neighbors")
print(f"IntAct:  {summary['IntAct']:>4} neighbors")
print(f"BioGRID: {summary['BioGRID']:>4} neighbors")
print(f"{'='*70}\n")

# Find consensus neighbors (present in multiple databases)
if isinstance(results['string'], pd.DataFrame) and isinstance(results['intact'], pd.DataFrame):
    string_genes = set(results['string']['preferred_name'].str.upper())
    intact_genes = set(results['intact']['neighbor_name'].str.upper())

    consensus = string_genes & intact_genes
    string_only = string_genes - intact_genes
    intact_only = intact_genes - string_genes

    print(f"Consensus (STRING ∩ IntAct): {len(consensus)} genes")
    print(f"STRING-specific: {len(string_only)} genes")
    print(f"IntAct-specific: {len(intact_only)} genes")
    print(f"\nConsensus genes (high confidence):")
    print(f"  {list(consensus)[:10]}")
```

### Comparative Metrics

| Metric | STRING | IntAct | BioGRID | Interpretation |
|--------|--------|--------|---------|----------------|
| **Coverage** | Broad | Moderate | Focused | STRING finds most interactions |
| **Curation** | Automated | Manual | Manual | IntAct/BioGRID highest quality |
| **Methods** | All evidence | Experimental | Physical+Genetic | Each provides unique view |
| **Speed** | Fast (~1s) | Moderate (~4s) | Slow (~5-10s) | STRING for rapid exploration |

**Best Practice**: Start with STRING for comprehensive coverage, then validate critical findings with IntAct and BioGRID for experimental details.

### Cross-Database Shortest Paths

Compare path structures across databases to identify the most supported routes:

```python
from scripts.unified_query import query_shortest_paths_all_databases

genes = ["TP53", "MDM2", "ATM", "CHEK2"]

paths = query_shortest_paths_all_databases(
    gene_list=genes,
    species=9606,
    max_distance=3,
    export_results=True,
    output_dir="./results"
)

# Analyze path agreement
for (gene_a, gene_b) in [(genes[i], genes[j]) for i in range(len(genes)) for j in range(i+1, len(genes))]:
    key = tuple(sorted([gene_a, gene_b]))

    found_in = []
    if isinstance(paths['string'], dict) and key in paths['string']:
        found_in.append('STRING')
    if isinstance(paths['intact'], dict) and key in paths['intact']:
        found_in.append('IntAct')
    if isinstance(paths['biogrid'], dict) and key in paths['biogrid']:
        found_in.append('BioGRID')

    if len(found_in) >= 2:
        print(f"✓ {gene_a} ↔ {gene_b}: Confirmed in {', '.join(found_in)}")

        # Compare path lengths
        if 'STRING' in found_in:
            string_path = paths['string'][key]
            print(f"  STRING:  {' → '.join(string_path['path'])} ({string_path['hops']} hops)")
        if 'IntAct' in found_in:
            intact_path = paths['intact'][key]
            print(f"  IntAct:  {' → '.join(intact_path['path'])} ({intact_path['hops']} hops)")
```

**When Paths Disagree**:
- Different path lengths → Databases have different interaction coverage
- Different intermediates → Alternative biological routes may exist
- Missing in one DB → May indicate database-specific bias or incomplete data

**Action**: Report consensus paths in publications; investigate discrepancies for novel insights.

## Database Selection Quick Guide

| Database | Best For | API Key | Query Time | Coverage |
|----------|----------|---------|------------|----------|
| **STRING** | High-throughput, comprehensive coverage | No | ~1s | 24B+ interactions |
| **IntAct** | Detailed curation, experimental methods | No | ~4s | 1M+ interactions |
| **BioGRID** | Genetic interactions, manual curation, multi-hop | Required | ~5-10s (1-hop) | 2M+ interactions |

**Decision Tree**:
1. **Need genetic interactions?** → Use BioGRID
2. **Need detailed experimental methods?** → Use IntAct
3. **Need fast, comprehensive results?** → Use STRING (recommended default)
4. **Need cross-validation?** → Use unified interface to query all three

## Filter Threshold Recommendations

**STRING**:
- Exploratory: `min_score=400`
- Standard: `min_score=700, min_experimental=500`
- High-confidence: `min_score=900, min_experimental=900`

**IntAct**:
- Exploratory: `min_miscore=0.4` (default)
- Standard: `min_miscore=0.6`
- High-confidence: `min_miscore=0.9`

**BioGRID**:
- Exploratory: `min_score=0.3`
- Standard: `min_score=0.5`
- High-confidence: `min_score=0.7, throughput_tag="low"`

See `references/database_comparison.md` for comprehensive threshold guidance and filter combination strategies.

## Parameter Selection Guide by Research Scenario

Choose parameters based on research goals, publication requirements, and biological context.

### Scenario 1: Exploratory Analysis - Discover Novel Connections

**Goal**: Cast a wide net to identify all potential interactions for hypothesis generation.

**Recommended Parameters**:
```python
# STRING - Broad coverage
results = query_single_gene_all_databases(
    gene="NOVEL_GENE",
    species=9606,
    top_n=200,                    # Request many neighbors
    min_combined_score=400,        # Permissive threshold
    max_hops=3                     # Include indirect neighbors
)
```

**Why These Parameters**:
- `min_score=400`: Includes moderate-confidence interactions (not just high-confidence)
- `top_n=200`: Captures broader network context
- `max_hops=3`: Discovers indirect relationships that may be biologically relevant
- No evidence-type filters: Maximizes coverage across all evidence types

**Analysis Approach**:
1. Query all three databases with permissive filters
2. Identify consensus interactions (present in 2+ databases)
3. Prioritize interactions with high combined scores for follow-up
4. Use database-specific interactions as exploratory leads

### Scenario 2: Publication-Quality Validation

**Goal**: Report only high-confidence, experimentally validated interactions for publication.

**Recommended Parameters**:
```python
# STRING - High experimental evidence
results = query_single_gene_all_databases(
    gene="TARGET_GENE",
    species=9606,
    top_n=50,
    min_combined_score=700,         # High confidence
    min_experimental_score=800,     # Strong experimental support
    min_database_score=500,         # Curated database evidence
    organism_filter="homo sapiens"  # IntAct: species-specific
)
```

**Why These Parameters**:
- `min_experimental=800`: Prioritizes direct experimental evidence (Y2H, Co-IP, etc.)
- `min_database=500`: Includes curated pathway databases (KEGG, Reactome)
- `organism_filter`: Ensures species-specific interactions (avoid orthology inference)
- Multiple evidence requirements: Reduces false positives

**Analysis Approach**:
1. Query with strict filters to minimize false positives
2. Require consensus across 2+ databases
3. Manually review experimental methods in IntAct for critical interactions
4. Report detection methods and PubMed IDs in supplementary materials

### Scenario 3: Cross-Species Comparative Study

**Goal**: Compare interaction networks across human and model organisms.

**Recommended Parameters**:
```python
# IntAct - Cross-species support
human_results = get_direct_neighbors(
    gene="GENE_SYMBOL",
    species="human",
    top_n=100,
    min_miscore=0.6,
    organism_filter="homo sapiens"   # Human-only
)

mouse_results = get_direct_neighbors(
    gene="Gene_symbol",              # Note: mouse uses different capitalization
    species="mouse",
    top_n=100,
    min_miscore=0.6,
    organism_filter="mus musculus"   # Mouse-only
)

# Find conserved interactions
human_genes = set(human_results['neighbor_name'].str.upper())
mouse_genes = set(mouse_results['neighbor_name'].str.upper())
conserved = human_genes & mouse_genes
```

**Why These Parameters**:
- `organism_filter`: Critical for species-specific filtering
- `min_miscore=0.6`: Moderate confidence to balance coverage and quality
- IntAct database: Best organism annotation and cross-species data

**Analysis Approach**:
1. Query each species independently with organism filter
2. Identify conserved interactions (present in both species)
3. Species-specific interactions may indicate lineage-specific evolution
4. Use STRING for broader evolutionary comparison across more species

### Scenario 4: Genetic Interaction Analysis

**Goal**: Study synthetic lethality and genetic dependencies for drug target discovery.

**Recommended Parameters**:
```python
# BioGRID - Genetic interactions
from scripts.biogrid_api import BioGRIDClient

client = BioGRIDClient(api_key)
neighbors = client.get_neighbors(
    seed_gene="TARGET_GENE",
    tax_id="9606",
    max_hops=2,                           # Include indirect genetic interactions
    max_neighbors=100,
    experimental_system_types=["genetic"], # Genetic interactions only
    throughput_tag="low",                 # Manual curation preferred
    min_score=0.5
)
```

**Why These Parameters**:
- `experimental_system_types=["genetic"]`: Filters for synthetic lethality, suppression, enhancement
- `throughput_tag="low"`: Prioritizes manually curated small-scale studies
- `max_hops=2`: Captures genetic interaction cascades
- BioGRID: Best genetic interaction coverage

**Analysis Approach**:
1. Identify direct genetic interactors (synthetic lethal partners)
2. Explore 2-hop genetic interactions for network effects
3. Cross-reference with physical interactions to understand mechanism
4. Validate critical findings with STRING/IntAct for supporting evidence

### Scenario 5: Drug Target Network Analysis

**Goal**: Map drug target neighborhood to identify on-target effects and potential off-target interactions.

**Recommended Parameters**:
```python
# STRING - Physical interactions with high experimental evidence
target_network = query_single_gene_all_databases(
    gene="DRUG_TARGET",
    species=9606,
    top_n=50,
    min_combined_score=700,
    min_experimental_score=800,      # Direct experimental evidence
    network_type="physical"          # Physical binding only
)

# Find paths to known disease genes
from scripts.unified_query import query_shortest_paths_all_databases

disease_genes = ["DISEASE_GENE1", "DISEASE_GENE2", "DISEASE_GENE3"]
paths = query_shortest_paths_all_databases(
    gene_list=[DRUG_TARGET] + disease_genes,
    species=9606,
    max_distance=3,
    min_combined_score=700,
    min_experimental_score=800
)
```

**Why These Parameters**:
- `network_type="physical"`: Direct binding interactions relevant for drug mechanism
- `min_experimental=800`: High confidence for clinical relevance
- `max_distance=3`: Captures direct and indirect pathway effects
- Unified interface: Cross-validates findings across databases

**Analysis Approach**:
1. Map direct physical interactors (potential on-target effects)
2. Identify bridge proteins connecting target to disease pathways
3. Cross-reference with BioGRID genetic interactions for functional validation
4. Use IntAct to examine experimental methods and tissue specificity

### Scenario 6: Pathway Reconstruction

**Goal**: Reconstruct signaling or metabolic pathway from seed genes.

**Recommended Parameters**:
```python
# STRING - Functional network with coexpression evidence
pathway_seeds = ["GENE1", "GENE2", "GENE3"]

# Find connections between pathway members
paths = query_shortest_paths_all_databases(
    gene_list=pathway_seeds,
    species=9606,
    max_distance=2,                   # Pathways are tightly connected
    min_combined_score=700,
    min_coexpression_score=600,       # Co-regulation evidence
    min_database_score=500,           # Known pathway databases
    network_type="functional"
)
```

**Why These Parameters**:
- `max_distance=2`: Pathway genes are typically directly or 1-hop connected
- `min_coexpression=600`: Pathway members often co-regulated
- `min_database=500`: Known pathway databases (KEGG, Reactome) support
- `network_type="functional"`: Includes functional relationships beyond physical

**Analysis Approach**:
1. Start with known pathway members as seeds
2. Identify direct connections (validated pathway edges)
3. Bridge proteins suggest missing pathway components
4. Validate pathway structure with IntAct experimental methods

### Parameter Quick Reference

| Research Goal | Database | Key Parameters | Rationale |
|--------------|----------|----------------|-----------|
| **Exploratory** | All | `min_score=400, top_n=200` | Broad coverage |
| **Publication** | All | `min_experimental=800, consensus≥2` | High confidence |
| **Cross-species** | IntAct | `organism_filter="species"` | Species-specific |
| **Genetic interactions** | BioGRID | `experimental_system_types=["genetic"]` | Synthetic lethality |
| **Drug target** | STRING | `network_type="physical", min_experimental=800` | Direct binding |
| **Pathway** | STRING | `min_coexpression=600, min_database=500` | Co-regulation |

## Output Annotations

All queries return rich annotations:

**Entity Annotations**:
- Gene symbols and database-specific IDs (STRING ID, UniProt AC, Entrez ID)
- Organism taxonomy ID and name
- Hop distance from seed gene (1-hop, 2-hop, 3-hop)
- Path from seed to neighbor (e.g., "TP53-MDM2-ATM")

**Edge Annotations**:
- Confidence scores (database-specific)
- Experimental methods and detection systems
- Evidence types (experimental, database, textmining, etc.)
- Supporting publications (PubMed IDs)
- Interaction count across experiments

**Path Metadata** (shortest paths):
- Complete node path (e.g., `['TP53', 'MDM2', 'ATM']`)
- Hop count (number of edges)
- Total distance (sum of edge weights)
- Confidence scores per edge
- Algorithm name and weight formula

## Testing and Validation

To verify installation and test all functionality:

```bash
cd /home/sagemaker-user/.claude/skills/interactdb-query/scripts

# Quick validation test
python test_shortest_paths.py

# Comprehensive test with all filters
python comprehensive_test.py
```

Tests validate:
- STRING, IntAct, BioGRID single-gene queries
- Shortest path algorithms with various filters
- Multi-hop expansion (BioGRID)
- Evidence type filtering
- Export functionality

## Dependencies

**Required Python Packages**:
```bash
pip install requests pandas
```

**BioGRID API Key** (Required for BioGRID queries):
1. Register at https://thebiogrid.org/ to obtain a free API key
2. Set the environment variable:
```bash
export BIOGRID_API_KEY="your_key_here"
```

**Note:** STRING and IntAct work without API keys. Only BioGRID requires authentication.

## Performance and Timeouts

| Query Type | Database | Typical Time | Notes |
|------------|----------|--------------|-------|
| Single-gene (1-hop) | STRING | ~1s | Fastest option |
| Single-gene (1-hop) | IntAct | ~4s | Detailed annotations |
| Single-gene (1-hop) | BioGRID | ~5-10s | Manual curation focus |
| Single-gene (multi-hop) | BioGRID | ~1-5 min | NEW in v1.1, auto-expands |
| Shortest paths (3 genes) | STRING | ~1-2s | Fast, comprehensive |
| Shortest paths (5 genes) | STRING | ~2-5s | Most pairs found |
| Unified query (all DBs) | All | ~15-30s | Includes BioGRID multi-hop |

**Tips for Optimal Performance**:
- Start with STRING for initial exploration (fastest)
- Use permissive filters first, then increase thresholds
- For BioGRID multi-hop, expect 1-5 minute queries for sparse genes
- Export intermediate results to avoid re-querying
- For large gene lists (>10), consider chunking

## No Path Validation Workflow

When shortest path queries return no results, validate whether this is due to filtering thresholds or genuine biological disconnection using the following workflow:

### Automatic Validation Protocol

**Step 1: Initial Query**
- Query with standard or user-specified filters (e.g., `min_combined_score=400`)
- If paths found, proceed with results
- If no paths found, proceed to validation

**Step 2: Zero-Threshold Validation**
- Automatically re-query with least restrictive criteria for each database:

  **STRING:**
  ```python
  from scripts.string_api import StringClient

  client = StringClient()
  validation_paths = client.find_shortest_paths(
      gene_list=genes,
      species=9606,
      min_combined_score=0,          # Accept any score
      min_experimental_score=0,
      min_database_score=0,
      min_textmining_score=0,
      min_coexpression_score=0,
      max_network_expansion=20,      # Extended search
      max_distance=100               # Extended path length
  )
  ```

  **BioGRID:**
  ```python
  from scripts.biogrid_api import BioGRIDClient

  client = BioGRIDClient(api_key)
  validation_paths = client.find_shortest_paths(
      gene_list=genes,
      tax_id="9606",
      max_distance=100,              # Extended path length
      min_score=0.0                  # Accept any score
  )
  ```

  **IntAct:**
  ```python
  from scripts.intact_api import find_shortest_paths_intact

  validation_paths = find_shortest_paths_intact(
      gene_list=genes,
      species="human",
      max_distance=100,              # Extended path length
      min_miscore=0.0                # Accept any score
  )
  ```

**Step 3: Interpret Results**

If validation finds paths:
- **Conclusion**: Path exists but was filtered by threshold
- **Interpretation**: Genes are connected via low-confidence interactions
- **Recommendation**: Review validation path scores; may include speculative or weak evidence
- **Action**: Decide if lower threshold appropriate for analysis

If validation still finds no paths:
- **Conclusion**: Genes are genuinely disconnected
- **Interpretation**: Operate in separate biological pathways/network communities
- **Evidence**: No shared neighbors, no connecting proteins, distinct functions
- **Action**: Accept as valid biological finding

### Recording Validation Results

Document validation in output:

```python
{
  "query": "GENE_A ↔ GENE_B",
  "initial_result": "no_path",
  "initial_parameters": {"min_combined_score": 400},
  "validation_performed": true,
  "validation_parameters": {
    "min_combined_score": 0,
    "max_network_expansion": 20
  },
  "validation_result": "no_path",
  "interpretation": "Genuinely disconnected - genes operate in separate network communities",
  "evidence": {
    "network_communities": {
      "GENE_A": ["neighbor1", "neighbor2", "..."],
      "GENE_B": ["neighbor1", "neighbor2", "..."]
    },
    "shared_neighbors": 0,
    "biological_context": "GENE_A: [pathway], GENE_B: [pathway]"
  }
}
```

### Example: Validation in Practice

```python
from scripts.string_api import StringClient, get_string_neighbors

client = StringClient()

# Initial query
print("Initial query with min_score=400:")
paths = client.find_shortest_paths(
    gene_list=['TNFRSF25', 'GREM1'],
    species=9606,
    min_combined_score=400
)

if not paths:
    print("✗ No paths found - performing validation...")

    # Validation with zero threshold
    validation_paths = client.find_shortest_paths(
        gene_list=['TNFRSF25', 'GREM1'],
        species=9606,
        min_combined_score=0,
        max_network_expansion=20,
        max_distance=100
    )

    if validation_paths:
        print("✓ Validation: Path exists (filtered by threshold)")
        for (gene_a, gene_b), info in validation_paths.items():
            print(f"  Path: {' → '.join(info['path'])}")
            print(f"  Edge scores: {info['scores']}")
            print(f"  Interpretation: Low-confidence connection")
    else:
        print("✗ Validation: Genuinely disconnected")

        # Check network communities
        df_a = get_string_neighbors('TNFRSF25', species=9606, top_n=10, min_score=0)
        df_b = get_string_neighbors('GREM1', species=9606, top_n=10, min_score=0)

        print(f"\n  TNFRSF25 neighbors: {', '.join(df_a['preferred_name'].head(5).tolist())}")
        print(f"  GREM1 neighbors: {', '.join(df_b['preferred_name'].head(5).tolist())}")
        print(f"\n  Interpretation: Separate biological pathways")
        print(f"    TNFRSF25: TNF/death receptor signaling")
        print(f"    GREM1: BMP/morphogen signaling")
```

## Troubleshooting

**No paths found**:
1. **First**: Run validation with zero threshold (see "No Path Validation Workflow" above)
2. Lower score thresholds (e.g., 700 → 400 → 0)
3. Increase max_distance (e.g., 2 → 3 → 100)
4. Remove evidence type filters temporarily
5. Verify gene name spelling and species ID
6. Check if genes exist: `get_string_neighbors(gene, species=9606, top_n=1)`

**BioGRID timeout**:
- Reduce max_hops (use 1 or 2 instead of 3)
- Reduce max_neighbors (use 50 instead of 100)
- For well-connected genes, use 1-hop only
- Consider using STRING instead for initial exploration

**Too many results**:
- Increase score thresholds (e.g., 400 → 700)
- Decrease max_distance (e.g., 3 → 2)
- Add evidence type filters (e.g., `min_experimental=800`)
- Use `throughput_tag="low"` for BioGRID (manual curation only)

## Version Information

**Current Version**: 1.2.1 (January 2, 2026)

**Recent Updates**:
- **v1.2.1 (Jan 2026)**: Critical bug fix for shortest path queries - individual gene queries with parallel execution
  - **CRITICAL FIX**: Changed from batched gene queries to individual gene queries
  - **Problem**: Querying multiple genes together (e.g., `"GENE1|GENE2|GENE3"`) caused STRING API to return limited/different results
  - **Solution**: Query each gene separately with parallel threading (ThreadPoolExecutor, max 10 workers)
  - **STRING changes**:
    - Increased `add_nodes` from 10 to 1000 for better edge coverage
    - Increased node limit from 1000 to 5000 nodes
    - Query each gene in frontier separately with concurrent execution
  - **IntAct changes**: Query each gene separately with parallel threading for graph building
  - **BioGRID changes**: Query each gene separately with parallel threading
  - **Impact**: Revealed 11 additional paths previously missed (73% increase in path discovery)
  - **Performance**: ~5x speedup with parallel queries (e.g., 5 genes in 1.93s vs ~10s sequential)
  - **Backward compatibility**: Function signatures unchanged, all parameters preserved
- **v1.2.0 (Dec 2025)**: IntAct REST API migration (replaced PSICQUIC)
- **v1.1.2 (Dec 2025)**: Parameter standardization across databases
  - **STRING `network_type` now optional (defaults to None → "functional")**:
    - Changed from required `str = "functional"` to `Optional[str] = None`
    - When None, internally defaults to "functional" which captures **ALL interaction types**
    - Clarified: "functional" = physical + predicted + functional (broadest category)
    - Clarified: "physical" = only direct physical binding (subset)
  - **Default parameters standardized across all databases**:
    - `max_distance`: 50 (STRING, IntAct, BioGRID)
    - `max_network_expansion`: 20 (STRING)
    - `min_score`: 0.4 (BioGRID)
    - All databases now use consistent timeout (5 minutes for BioGRID)
  - **Aligned with IntAct/BioGRID** which don't have network_type parameter
- **v1.1.1 (Dec 2025)**: Critical bug fix and validation workflow
  - Fixed network expansion in `find_shortest_paths` (changed `add_nodes=0` to `add_nodes=10`)
  - Paths now correctly found between connected genes (e.g., TYK2 ↔ JAK1)
  - Added "No Path Validation Workflow" for distinguishing filtered vs. genuinely disconnected genes
  - Documented zero-threshold validation protocol for confirming biological disconnections
  - No performance impact, queries complete in 1-5 seconds as before
- **v1.1.0 (Nov 2025)**: BioGRID multi-hop expansion, extended timeout, unified interface

See `references/CHANGELOG.md` for complete version history and migration guide.

## References

- `references/database_comparison.md`: Complete database specifications, score systems, threshold recommendations
- `references/CHANGELOG.md`: Version history, migration guide, known issues
- `references/SKILL_UPDATE_SUMMARY.md`: Version 1.1.0 implementation details

## Bundled Scripts

- `scripts/string_api.py`: STRING database client with multi-hop BFS
- `scripts/intact_api.py`: IntAct database client with REST API (v1.2.0+, replaced PSICQUIC)
- `scripts/intact_api_psicquic_backup.py`: Original PSICQUIC implementation (backup)
- `scripts/biogrid_api.py`: BioGRID database client with multi-hop expansion
- `scripts/unified_query.py`: Unified interface for all three databases
- `scripts/test_shortest_paths.py`: Validation test suite
- `scripts/comprehensive_test.py`: Extended test suite with filter combinations
