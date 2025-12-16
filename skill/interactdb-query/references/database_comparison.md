# Database Comparison and Parameter Reference

Comprehensive reference for STRING, IntAct, and BioGRID interaction databases including score systems, parameter specifications, and filtering options.

## Database Overview

| Feature | STRING | IntAct | BioGRID |
|---------|--------|--------|---------|
| **Data Source** | Multi-source aggregator | Curated literature | Curated literature |
| **API Access** | Free, no key | Free, PSICQUIC | Requires API key |
| **Total Interactions** | 24B+ | 1M+ | 2M+ |
| **Organisms** | 5000+ | 500+ | 70+ |
| **Update Frequency** | Regular | Regular | Regular |
| **Score Range** | 0-999 (per channel) | 0.0-1.0 (MI-score) | Varies (QUANTITATION) |
| **Evidence Types** | 8 channels | PSI-MI ontology | Physical/Genetic |
| **Best For** | High-throughput, confident interactions | Detailed curation, publications | Genetic interactions |

## STRING Database

### Evidence Channels

All scores range from 0-999:

1. **combined_score**: Overall confidence combining all evidence
   - <400: Low confidence
   - 400-699: Medium confidence
   - 700-899: High confidence
   - 900-999: Highest confidence

2. **experimental_score**: From experimental data repositories
   - High values indicate direct experimental evidence
   - Recommended threshold: 800+ for high-confidence experiments

3. **database_score**: From curated interaction databases
   - Includes KEGG, Reactome, etc.
   - Recommended threshold: 500+ for curated evidence

4. **textmining_score**: From literature co-mentions
   - Text mining of PubMed abstracts
   - Often high but may indicate co-occurrence rather than direct interaction
   - Use with caution or combine with other evidence

5. **coexpression_score**: From gene expression patterns
   - Co-expression across conditions
   - Recommended threshold: 600+ for functional associations

6. **neighborhood_score**: From genomic neighborhood
   - Gene proximity in bacterial genomes
   - Most relevant for prokaryotes

7. **fusion_score**: From gene fusion events
   - Gene fusion across species
   - Indicates functional coupling

8. **cooccurrence_score**: From phylogenetic co-occurrence
   - Co-occurrence across species
   - Indicates functional association

### Network Types

- **functional**: All evidence types (default, recommended)
- **physical**: Only direct physical interactions (more restrictive)

### Common Parameter Combinations

```python
# High-confidence experimental interactions
min_score=700, min_experimental=800

# Curated database interactions
min_score=700, min_database=700

# Multi-evidence validation
min_score=700, min_experimental=500, min_database=500

# Conservative (highest confidence)
min_score=900, min_experimental=900

# Permissive (exploratory)
min_score=400
```

## IntAct Database

### MI-Score (Molecular Interaction Score)

Range: 0.0 (low) to 1.0 (high confidence)

**Confidence Levels**:
- 0.0-0.39: Low confidence
- 0.4-0.69: Moderate confidence (default threshold: 0.4)
- 0.7-0.89: High confidence
- 0.9-1.0: Very high confidence

**Recommended Thresholds**:
- Exploratory analysis: 0.4
- Standard analysis: 0.6
- High-confidence analysis: 0.7
- Publication-quality: 0.9

### PSI-MI Ontology Annotations

IntAct provides rich PSI-MI (Proteomics Standards Initiative - Molecular Interactions) ontology terms:

**Detection Methods** (examples):
- "two hybrid" - Yeast two-hybrid
- "anti bait coimmunoprecipitation" - Co-IP
- "pull down" - Pull-down assay
- "fluorescence microscopy" - Imaging-based
- "surface plasmon resonance" - Biophysical

**Interaction Types** (examples):
- "physical association" - Direct or indirect physical interaction
- "direct interaction" - Confirmed direct binding
- "association" - General association
- "colocalization" - Same cellular location

### Organism Filtering

IntAct supports cross-species filtering:

**Common Organisms**:
- "homo sapiens" (Human)
- "mus musculus" (Mouse)
- "rattus norvegicus" (Rat)
- "drosophila melanogaster" (Fruit fly)
- "saccharomyces cerevisiae" (Yeast)
- "caenorhabditis elegans" (C. elegans)

**Usage**: Comma-separated list for multiple organisms

```python
organism_filter="homo sapiens,mus musculus"
```

## BioGRID Database

### QUANTITATION Score

BioGRID uses QUANTITATION field for confidence scores:
- Range: Varies by experimental system
- Often represents p-values, confidence scores, or numeric measurements
- Not all interactions have QUANTITATION values

**Recommended Thresholds**:
- Exploratory: 0.3-0.5
- Standard: 0.5-0.7
- High-confidence: 0.7-0.9

### Evidence Types

**experimental_system_types** parameter:

1. **"physical"**: Direct physical interactions
   - Two-hybrid assays
   - Co-immunoprecipitation
   - Protein-fragment complementation
   - Biochemical assays

2. **"genetic"**: Genetic interactions
   - Synthetic lethality
   - Dosage rescue
   - Phenotypic enhancement
   - Suppression

**Usage**: List of desired types

```python
experimental_system_types=["physical"]  # Physical only
experimental_system_types=["genetic"]   # Genetic only
experimental_system_types=["physical", "genetic"]  # Both
experimental_system_types=None  # All types (default)
```

### Throughput Filtering

**throughput_tag** parameter:

- **"any"**: All throughput levels (default)
- **"low"**: Low-throughput (manual curation, small-scale studies)
- **"high"**: High-throughput (large-scale screens)

**Recommendation**: Use "low" for highest confidence manually curated interactions

## Multi-Hop Expansion

### BioGRID Multi-Hop Capability

BioGRID supports multi-hop network expansion via the `max_hops` parameter to automatically fill `top_n` neighbor requests:

- **max_hops=1**: Direct neighbors only
- **max_hops=2**: Expand to 2-hop neighbors
- **max_hops=3**: Expand to 3-hop neighbors (default in unified query)

**Performance**: Multi-hop queries may take 1-5 minutes depending on gene connectivity. Default timeout is 5 minutes (300 seconds).

**Use Case**: Particularly useful for genes with sparse direct interactions. Example: GREM1 has only 4 direct neighbors, but with max_hops=3 and top_n=30, finds 4 (1-hop) + 26 (2-hop) = 30 total neighbors.

**Note**: STRING and IntAct only support direct neighbor queries via their simple APIs. For multi-hop paths between specific genes, use the shortest path functions available in all three databases.

## Shortest Path Algorithm Details

### Edge Weight Conversions

Different databases require different weight conversion formulas to ensure higher scores = shorter distances:

1. **STRING**:
   ```
   weight = 1000 - combined_score
   ```
   - Score 999 → weight 1 (very close)
   - Score 400 → weight 600 (distant)

2. **IntAct**:
   ```
   weight = 1.0 - miscore
   ```
   - MI-score 0.9 → weight 0.1 (very close)
   - MI-score 0.4 → weight 0.6 (distant)

3. **BioGRID**:
   ```
   weight = 1.0 / (score + 0.01)
   ```
   - Score 1.0 → weight ~0.99 (very close)
   - Score 0.1 → weight ~9.1 (distant)
   - +0.01 prevents division by zero

### Output Metadata

All shortest path functions return:

```python
{
    'path': ['GENE_A', 'INTERMEDIATE', 'GENE_B'],
    'distance': 1.234,  # Sum of edge weights
    'hops': 2,  # Number of edges
    'scores': [999, 850],  # Original confidence scores per edge
    'algorithm': 'Dijkstra',
    'weight_formula': 'weight = 1000 - combined_score'
}
```

## Performance Optimization

### Score Threshold Impact

| Threshold | Results | Query Speed | Recommendation |
|-----------|---------|-------------|----------------|
| Very permissive (<400) | Many results | Slower | Initial exploration |
| Moderate (400-700) | Balanced | Normal | Standard analysis |
| Strict (700-900) | Fewer results | Faster | High-confidence |
| Very strict (>900) | Few results | Fastest | Publication |

### Max Distance Impact

| max_distance | Pairs Checked | Graph Size | Query Time |
|--------------|---------------|------------|------------|
| 1 | Direct only | Small | Fast (~1s) |
| 2 | 1-2 hops | Medium | Medium (~2-5s) |
| 3 | 1-3 hops | Large | Slower (~5-15s) |
| >3 | Many hops | Very large | Very slow (>15s) |

**Recommendation**: Start with max_distance=2, increase to 3 only if needed.

## Common Taxon IDs

| Organism | Taxon ID | Common Name |
|----------|----------|-------------|
| Homo sapiens | 9606 | Human |
| Mus musculus | 10090 | Mouse |
| Rattus norvegicus | 10116 | Rat |
| Drosophila melanogaster | 7227 | Fruit fly |
| Caenorhabditis elegans | 6239 | C. elegans |
| Saccharomyces cerevisiae | 4932 | Yeast |
| Escherichia coli K-12 | 511145 | E. coli K-12 |
| Arabidopsis thaliana | 3702 | Thale cress |
| Danio rerio | 7955 | Zebrafish |

## API Rate Limits and Constraints

### STRING
- No API key required
- No explicit rate limits
- Recommended: <100 queries/minute for courtesy

### IntAct (PSICQUIC)
- No API key required
- Max results per query: Configurable (default 50000)
- Recommended: <50 queries/minute

### BioGRID
- **Requires API key** (register at https://thebiogrid.org/)
- **Default timeout**: 5 minutes (300 seconds) to support multi-hop queries
- Rate limit: Varies by key tier
- Max results per query: Configurable (default 10000)
- Recommended: <20 queries/minute
- **Note**: Multi-hop queries (max_hops > 1) may take 1-5 minutes

## Filter Combination Strategies

### Conservative (High Confidence)

**STRING**:
```python
min_combined_score=900
min_experimental_score=900
```

**IntAct**:
```python
min_miscore=0.9
```

**BioGRID**:
```python
min_score=0.8
throughput_tag="low"
experimental_system_types=["physical"]
```

### Balanced (Standard Analysis)

**STRING**:
```python
min_combined_score=700
min_experimental_score=500
```

**IntAct**:
```python
min_miscore=0.6
```

**BioGRID**:
```python
min_score=0.5
throughput_tag="any"
```

### Exploratory (Permissive)

**STRING**:
```python
min_combined_score=400
```

**IntAct**:
```python
min_miscore=0.4
```

**BioGRID**:
```python
min_score=0.3
throughput_tag="any"
```

## Cross-Database Validation

To validate findings across databases:

1. **Query all three databases** with equivalent thresholds:
   - STRING: min_score=700
   - IntAct: min_miscore=0.7
   - BioGRID: min_score=0.7

2. **Find consensus interactions**: Present in 2+ databases

3. **Examine database-specific interactions**: Unique to one database

4. **Compare path structures**: Same endpoints, different intermediates

## Troubleshooting Common Issues

### No Results Returned

**Problem**: Query returns zero results

**Solutions**:
1. Lower score thresholds (e.g., 700 → 400)
2. Increase max_distance (e.g., 2 → 3)
3. Remove evidence type filters temporarily
4. Check gene name spelling and case
5. Try alternative gene names/aliases

### Too Many Results

**Problem**: Query returns excessive results (>1000)

**Solutions**:
1. Increase score thresholds (e.g., 400 → 700)
2. Decrease max_distance (e.g., 3 → 2)
3. Add evidence type filters (e.g., min_experimental=800)
4. Reduce top_n parameter
5. Focus on specific network type (e.g., "physical")

### Slow Queries

**Problem**: Queries take >30 seconds

**Solutions**:
1. Use STRING instead of IntAct/BioGRID
2. Reduce gene list size (<5 genes)
3. Decrease max_distance (use 2 instead of 3)
4. Increase score thresholds to reduce graph size
5. Cache intermediate results

### Disconnected Graphs

**Problem**: No paths found between query genes

**Solutions**:
1. Lower score thresholds to allow weaker connections
2. Increase max_distance to allow longer paths
3. Check if genes are in the same organism
4. Try different database (STRING has broader coverage)
5. Verify gene names are correct

## Best Practice Workflows

### Publication-Quality Analysis

1. Query all three databases with strict thresholds
2. Find consensus interactions (present in 2+ databases)
3. Verify paths have high experimental evidence scores
4. Export with all annotations for supplementary material
5. Document exact parameters used

### Exploratory Analysis

1. Start with STRING (fastest, broadest coverage)
2. Use permissive thresholds (min_score=400)
3. Identify interesting patterns
4. Follow up with IntAct for detailed annotations
5. Cross-validate critical findings with BioGRID

### Network Visualization

1. Query with balanced thresholds
2. Export edge and node lists
3. Include confidence scores as edge weights
4. Mark query genes vs bridge proteins
5. Use layout algorithms weighted by confidence
