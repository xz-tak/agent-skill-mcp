# MCP Tool Usage Reference

## Available MCP Tools

### 1. BioBridge MCP (`mcp__biobridge__predict_associations`)

**Purpose**: Neural knowledge graph link prediction using BioBridge multimodal embeddings.

**Individual Entity Query**:
```python
mcp__biobridge__predict_associations(
    context="Find diseases associated with GREM1 gene",
    override_head_name="GREM1",
    override_head_type="gene/protein",
    override_tail_name="Crohn disease",  # Optional: specific target
    override_tail_type="disease",
    relation_hint="associated with",
    topk=1,  # For specific pair validation
    include_relation_catalog=False
)
```

**Combo/Signature Query (Mean Embeddings)**:
```python
mcp__biobridge__predict_associations(
    context="Find diseases associated with TYK2+JAK1 gene signature",
    override_head_names=["TYK2", "JAK1"],  # List of entities
    override_head_type="gene/protein",     # Required for multi-head
    override_tail_name="Crohn disease",
    override_tail_type="disease",
    relation_hint="associated with",
    topk=1
)
```

**Key Parameters**:
- `context` (required): Natural language description
- `override_head_name`: Single entity name
- `override_head_names`: List of entities (for signatures/combos)
- `override_head_type`: Entity type (required for multi-head)
- `override_tail_name`: Target entity name (optional)
- `override_tail_type`: Target entity type
- `relation_hint`: Relation family (e.g., "associated with")
- `topk`: Number of results (default: 25)

**Output Fields**:
- `cos_sim`: Cosine similarity score (0-1)
- `pct_rank`: Percentile rank (0-1, higher = stronger)
- `node_index`: KG node identifier
- `node_name`/`mondo_name`: Entity name

---

### 2. ULTRA MCP (`mcp__ultra-inference__predict_tail_entities`)

**Purpose**: Foundation model for single entity link prediction on PrimeKG.

**Individual Entity Query**:
```python
mcp__ultra-inference__predict_tail_entities(
    head_entity="GREM1",  # Entity name or ID
    relation="associated with",
    top_k=100
)
```

**Key Parameters**:
- `head_entity`: Entity name or ID (e.g., "GREM1", "MONDO:5301")
- `relation`: Relation label (e.g., "associated with", "ppi")
- `top_k`: Number of predictions (default: None for all)

**Output Fields**:
- `entity_id`: Entity identifier
- `entity_name`: Entity name
- `entity_type`: Entity type
- `score`: Model prediction score
- `percentile_rank`: Rank percentile
- `output_file`: Path to full results parquet

**Note**: ULTRA single-hop queries only. For combo/intersection queries, use UltraQuery.

---

### 3. UltraQuery MCP (`mcp__ultraquery-inference__answer_complex_query`)

**Purpose**: Complex logical queries with intersection, union, and multi-hop reasoning.

**2-Gene Intersection Query (2i)**:
```python
mcp__ultraquery-inference__answer_complex_query(
    query_structure=[
        ["TYK2", ["associated with"]],
        ["JAK1", ["associated with"]]
    ],
    top_k=25
)
```

**3-Gene Intersection Query (3i)**:
```python
mcp__ultraquery-inference__answer_complex_query(
    query_structure=[
        ["ITGA4", ["associated with"]],
        ["ITGB7", ["associated with"]],
        ["CDKN2D", ["associated with"]]
    ],
    top_k=25
)
```

**Multi-hop Query (2p)**:
```python
mcp__ultraquery-inference__answer_complex_query(
    query_structure=["GREM1", ["ppi", "associated with"]],
    top_k=25
)
```

**Key Parameters**:
- `query_structure`: Nested list in BetaE format
- `top_k`: Number of filtered predictions (default: 25)

**Supported Query Types**:
- `1p`: Single projection `["entity", ["relation"]]`
- `2p`: Two projections `["entity", ["r1", "r2"]]`
- `2i`: Two intersections `[["e1", ["r1"]], ["e2", ["r2"]]]`
- `3i`: Three intersections
- `ip`: Intersection then projection
- `pi`: Projection then intersection

**Output Fields**:
- `predictions`: List with rank, entity_id, entity_name, score, percentile_rank
- `query_type`: Detected query pattern
- `filtered_predictions`: Count after schema filtering
- `output_file_filtered`: Path to filtered results

---

## Entity Type Reference

Valid entity types for PrimeKG:
- `gene/protein` - Genes and proteins
- `disease` - Diseases (MONDO IDs)
- `drug` - Chemical compounds
- `effect/phenotype` - Effects and phenotypes
- `pathway` - Biological pathways
- `biological_process` - GO biological processes
- `molecular_function` - GO molecular functions
- `cellular_component` - GO cellular components
- `anatomy` - Anatomical terms
- `exposure` - Environmental exposures

---

## Relation Reference

Common relations:
- `associated with` - Gene-disease, protein-phenotype associations
- `ppi` - Protein-protein interactions
- `treats` - Drug-disease therapeutic relationships
- `side effect` - Drug adverse effects
- `interacts with` - General interactions
- `regulates` - Regulatory relationships
- `participates in` - Pathway/process participation

---

## Error Handling

**Entity Not Found**:
- BioBridge: Returns empty results or error in response
- ULTRA: Returns error message with suggestions

**Invalid Relation**:
- UltraQuery validates against PrimeKG schema
- Returns clear error explaining why query is invalid

**Recommended Approach**:
1. Validate entity existence with single query first
2. If not found, try alternative names/aliases
3. Check entity type is correct
