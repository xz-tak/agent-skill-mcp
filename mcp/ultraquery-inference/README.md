# UltraQuery Inference MCP Server

Inference-only MCP server for answering complex logical queries on PrimeKG using the UltraQuery model.

## Overview

This server provides on-demand inference for complex logical queries that involve:
- **Multi-hop reasoning**: Follow paths of multiple relations
- **Intersection**: Find entities satisfying multiple conditions simultaneously
- **Union**: Find entities satisfying any of several conditions
- **Negation**: Exclude entities matching certain patterns

Unlike simple link prediction `(head, relation, ?)`, UltraQuery can answer questions like:
- "What diseases are associated with proteins that interact with GREM1?" (2-hop)
- "What proteins interact with both GREM1 and IL11?" (intersection)
- "What diseases are associated with proteins that interact with GREM1 but NOT IL11?" (negation)

## Installation

### Prerequisites

```bash
# Install dependencies (from ULTRA directory)
cd ../
pixi install
# Or use pip:
# pip install torch torch-geometric torch-scatter fastmcp polars easydict
```

### Model Checkpoints

The server looks for two checkpoints:

1. **Ultra checkpoint** (base model for link prediction):
   - Default: `../ckpts/ultra_primekg_50g_ft_epoch_1.pth`
   - Environment variable: `ULTRA_CHECKPOINT_PATH`

2. **UltraQuery checkpoint** (query-specific layers):
   - Default: `../ckpts/ultraquery_primekg_ft_epoch_1.pth`
   - Environment variable: `ULTRAQUERY_CHECKPOINT_PATH`

**Note**: If the UltraQuery checkpoint is not available, the server will use only the base Ultra model, which may have reduced performance on complex queries.

### Graph Schema Setup

On first use, generate the PrimeKG schema file:

```bash
pixi run python build_schema.py
```

This creates `primekg_schema.json` which contains the complete graph schema for query validation and result filtering.

## Usage

### Starting the Server

```bash
# From ultraquery-inference-mcp directory
python ultraquery_inference_server.py
```

The server runs on HTTP by default (stateless mode).

### Query Format

Queries use the **BetaE nested list format**:

```python
# General structure:
# - Entities: "entity_id" or "entity_name"
# - Relations: "relation_label"
# - Operations:
#   - Projection: ["entity", ["relation"]]
#   - Multi-hop: ["entity", ["r1", "r2", ...]]
#   - Intersection: [["e1", ["r1"]], ["e2", ["r2"]]]
#   - Union: [["e1", ["r1"]], ["e2", ["r2"]], ["u"]]
#   - Negation: ["entity", ["relation", "n"]]
```

## Supported Query Types

| Type | Description | Example |
|------|-------------|---------|
| **1p** | One-hop projection | `["GREM1", ["ppi"]]` |
| **2p** | Two-hop projection | `["GREM1", ["ppi", "associated with"]]` |
| **3p** | Three-hop projection | `["entity", ["r1", "r2", "r3"]]` |
| **2i** | Two-way intersection | `[["GREM1", ["ppi"]], ["IL11", ["ppi"]]]` |
| **3i** | Three-way intersection | `[["e1", ["r1"]], ["e2", ["r2"]], ["e3", ["r3"]]]` |
| **ip** | Intersection → projection | `[[["GREM1", ["ppi"]], ["IL11", ["ppi"]]], ["associated with"]]` |
| **pi** | Projection → intersection | `[["entity", ["r1", "r2"]], ["entity2", ["r3"]]]` |
| **2in** | Intersection with negation | `[["e1", ["r1"]], ["e2", ["r2", "n"]]]` |
| **3in, inp, pin, pni** | Complex negation patterns | Various combinations |
| **2u-DNF** | Union (DNF) | `[["e1", ["r1"]], ["e2", ["r2"]], ["u"]]` |
| **up-DNF** | Union → projection | `[[["e1", ["r1"]], ["e2", ["r2"]], ["u"]], ["r3"]]` |

## API Tools

### 1. `answer_complex_query`

Answer a complex logical query.

**Parameters:**
- `query_structure` (required): Nested list in BetaE format
- `top_k` (optional): Number of top predictions to return (default: 25)

**Returns:**
```json
{
  "success": true,
  "query_type": "2i",
  "query_readable": "A <- projection_5(100)\nB <- projection_5(200)\nC <- intersection(A, B)",
  "entity_map": {
    "100": {"entity_id": "NCBI:9796", "entity_name": "GREM1"},
    "200": {"entity_id": "NCBI:3589", "entity_name": "IL11"}
  },
  "relation_map": {"5": "ppi"},
  "predictions": [
    {
      "rank": 1,
      "entity_id": "NCBI:5340",
      "entity_name": "PLG",
      "entity_type": "gene/protein",
      "score": 0.89,
      "schema_match": true
    }
  ],
  "total_predictions": 129312,
  "filtered_predictions": 32391,
  "expected_tail_types": ["disease", "effect/phenotype"],
  "inference_time_seconds": 1.23
}
```

### 2. `list_query_types`

List all supported query types with examples.

**Returns:**
```json
{
  "success": true,
  "query_types": ["1p", "2p", "3p", "2i", "3i", ...],
  "structures": {...},
  "descriptions": {...}
}
```

### 3. `setup_primekg_data`

Check if PrimeKG data is available and optionally install it if missing.

This tool ensures that PrimeKG data required for UltraQuery inference is properly set up. It can check the status of existing data or download and process the complete PrimeKG dataset from Harvard Dataverse.

**Parameters:**
- `dataset_path` (optional): Path to dataset directory. If None, uses default from server config.
- `force_redownload` (optional): If True, re-download even if data exists (default: False)
- `check_only` (optional): If True, only check status without installing (default: False)

**Returns:**
```json
{
  "success": true,
  "status": "available",
  "data_path": "/path/to/data/primekg1",
  "exists": {
    "primekg.csv": true,
    "train.txt": true,
    "test.txt": true,
    "valid.txt": true,
    "nodes.txt": true
  },
  "message": "PrimeKG data is available and ready to use"
}
```

Or when installing:
```json
{
  "success": true,
  "status": "installed",
  "data_path": "/path/to/data/primekg1",
  "nodes_count": 129375,
  "train_edges": 6480398,
  "test_edges": 810050,
  "valid_edges": 810050,
  "message": "PrimeKG data successfully installed with 129,375 nodes and 6,480,398 training edges"
}
```

**Usage Examples:**

Check if data is available:
```python
# Check without installing
result = setup_primekg_data(check_only=True)
if result["status"] == "missing":
    print("Data needs to be installed")
```

Install PrimeKG data:
```python
# Install if missing
result = setup_primekg_data()
if result["status"] == "installed":
    print(f"Installed {result['nodes_count']:,} nodes")
```

Force re-download:
```python
# Re-download and reprocess
result = setup_primekg_data(force_redownload=True)
```

**Notes:**
- Download size: ~500 MB (primekg.csv)
- Processing time: ~2-5 minutes
- Disk space required: ~2 GB total
- Data source: Harvard Dataverse

## Example Queries

### Example 1: Two-hop Projection (2p)

**Question**: What diseases are associated with proteins that interact with GREM1?

```python
query = ["GREM1", ["ppi", "associated with"]]
```

This query:
1. Finds proteins that interact with GREM1 (via `ppi` relation)
2. Finds diseases associated with those proteins (via `associated with` relation)

**Expected tail types**: `{disease, effect/phenotype}`

### Example 2: Intersection (2i)

**Question**: What proteins interact with both GREM1 and IL11?

```python
query = [["GREM1", ["ppi"]], ["IL11", ["ppi"]]]
```

This query finds the intersection of:
- Proteins that interact with GREM1
- Proteins that interact with IL11

**Expected tail types**: `{gene/protein}`

### Example 3: Intersection then Projection (ip)

**Question**: What diseases are associated with proteins that interact with both GREM1 and IL11?

```python
query = [[["GREM1", ["ppi"]], ["IL11", ["ppi"]]], ["associated with"]]
```

This query:
1. Finds proteins that interact with both GREM1 and IL11 (intersection)
2. Finds diseases associated with those proteins (projection)

**Expected tail types**: `{disease, effect/phenotype}`

### Example 4: Intersection with Negation (2in)

**Question**: What proteins interact with GREM1 but NOT IL11?

```python
query = [["GREM1", ["ppi"]], ["IL11", ["ppi", "n"]]]
```

This query finds proteins in:
- (Proteins that interact with GREM1) AND NOT (Proteins that interact with IL11)

**Note**: Negation queries do not use schema filtering due to complex set operations.

### Example 5: Using Entity IDs

You can also use entity IDs directly:

```python
# Using entity IDs instead of names
query = [["NCBI:9796", ["ppi"]], ["NCBI:3589", ["ppi"]]]
# Equivalent to: [["GREM1", ["ppi"]], ["IL11", ["ppi"]]]
```

## Python Client Example

```python
import requests

# Server URL
url = "http://localhost:8000/answer_complex_query"

# Example: Intersection query
query_structure = [["GREM1", ["ppi"]], ["IL11", ["ppi"]]]

response = requests.post(url, json={
    "query_structure": query_structure,
    "top_k": 25
})

result = response.json()

if result["success"]:
    print(f"Query type: {result['query_type']}")
    print(f"\nTop predictions:")
    for pred in result["predictions"][:10]:
        print(f"  {pred['rank']}. {pred['entity_name']} ({pred['entity_id']}) - Score: {pred['score']:.4f}")
else:
    print(f"Error: {result['error']}")
```

## Configuration

Environment variables:
- `ULTRA_CHECKPOINT_PATH`: Path to Ultra checkpoint (default: `../ckpts/ultra_primekg_50g_ft_epoch_1.pth`)
- `ULTRAQUERY_CHECKPOINT_PATH`: Path to UltraQuery checkpoint (default: `../ckpts/ultraquery_primekg_ft_epoch_1.pth`)
- `ULTRA_OUTPUT_DIR`: Output directory for results (default: `./output`)

## Architecture

### Model Components

1. **Base Ultra Model**: Link prediction foundation model
   - `RelNBFNet`: Learns relation representations
   - `EntityNBFNet`: Performs message passing on entity graph

2. **UltraQuery Wrapper**: Query execution layer
   - **Logic System**: Product fuzzy logic (configurable: product/godel/lukasiewicz)
   - **Query Execution**: Interprets postfix notation and executes operations
   - **Operations**:
     - Projection: Follow a relation
     - Intersection: Fuzzy AND operation
     - Union: Fuzzy OR operation
     - Negation: Fuzzy NOT operation

### Inference Pipeline

```
User Query (nested lists)
    ↓
Parse & Validate (entities/relations → indices)
    ↓
Validate against graph schema
    ↓
Convert to Query object (postfix notation)
    ↓
UltraQuery.forward(graph, query)
    ↓
Predictions (scores for all entities)
    ↓
Filter by expected tail types
    ↓
Top-k results with metadata
```

## Graph Schema System

The server uses a comprehensive graph schema for query validation and result filtering.

### Entity Types (10 types)

| Entity Type          | Occurrence Count | Description                                      |
| -------------------- | ---------------: | ------------------------------------------------ |
| drug                 |        5,611,392 | Pharmaceutical compounds from DrugBank           |
| gene/protein         |        5,262,458 | Genes and proteins from NCBI                     |
| anatomy              |        3,132,308 | Anatomical structures from UBERON                |
| disease              |          682,488 | Diseases and disorders from MONDO                |
| effect/phenotype     |          514,192 | Phenotypes and clinical effects from HPO         |
| biological_process   |          504,404 | Biological processes from GO                     |
| molecular_function   |          193,446 | Molecular functions from GO                      |
| cellular_component   |          186,204 | Cellular components from GO                      |
| pathway              |           95,432 | Biological pathways from Reactome                |
| exposure             |           18,672 | Environmental exposures                          |

### Relations (18 types)

| Relation                    | Edge Count | Description                                      |
| --------------------------- | ---------: | ------------------------------------------------ |
| expression present          |  3,036,406 | Gene/protein expressed in anatomy                |
| synergistic interaction     |  2,672,628 | Synergistic drug-drug interactions               |
| interacts with              |    686,550 | Gene/protein interactions with GO terms/pathways |
| ppi                         |    642,150 | Protein-protein interactions                     |
| phenotype present           |    300,634 | Phenotypes present in diseases                   |
| parent-child                |    281,744 | Hierarchical relationships (ontology structure)  |
| associated with             |    167,482 | Disease-gene associations                        |
| side effect                 |    129,568 | Drug side effects (adverse reactions)            |
| contraindication            |     61,350 | Drugs contraindicated for diseases               |
| expression absent           |     39,774 | Gene/protein not expressed in anatomy            |
| target                      |     32,760 | Drug-target relationships                        |
| indication                  |     18,776 | Drug indications for diseases                    |
| enzyme                      |     10,634 | Gene/protein acts as enzyme for drug             |
| transporter                 |      6,184 | Gene/protein transports drug                     |
| off-label use               |      5,136 | Off-label drug uses for diseases                 |
| linked to                   |      4,608 | Exposure linked to disease                       |
| phenotype absent            |      2,386 | Phenotypes absent in diseases                    |
| carrier                     |      1,728 | Gene/protein carries drug                        |

### Schema Validation

All queries are validated against the PrimeKG schema before inference:

**Validation Rules:**
1. **Single hop (1p)**: `(head_type, relation)` must exist in schema
2. **Multi-hop (2p, 3p)**: Each step in path must be valid
3. **Intersection (2i, 3i)**: All branches must be valid AND have non-empty intersection
4. **Compound (ip, pi)**: All sub-queries must be valid AND final intersection non-empty

**Example Valid Combinations:**
- `gene/protein` → `associated with` → `{disease, effect/phenotype}`
- `gene/protein` → `ppi` → `gene/protein`
- `drug` → `target` → `gene/protein`
- `disease` → `phenotype present` → `effect/phenotype`

**Example Invalid Combinations:**
- `drug` → `expression present` → ❌ (drugs are not expressed)
- `disease` → `ppi` → ❌ (diseases don't have protein interactions)

### Expected Tail Type Computation

For each query type, the schema defines what entity types are valid as answers:

**1p (Single hop):**
```python
query = ["GREM1", ["associated with"]]
# GREM1 = gene/protein
# schema[gene/protein, associated with] = {disease, effect/phenotype}
# Expected: disease or effect/phenotype entities only
```

**2i (Intersection):**
```python
query = [["TYK2", ["associated with"]], ["JAK1", ["associated with"]]]
# Both are gene/protein
# schema[gene/protein, associated with] = {disease, effect/phenotype}
# Intersection: {disease, effect/phenotype} ∩ {disease, effect/phenotype}
# Expected: disease or effect/phenotype entities only
```

**2p (Two-hop path):**
```python
query = ["GREM1", ["ppi", "associated with"]]
# GREM1 = gene/protein
# Step 1: schema[gene/protein, ppi] = {gene/protein}
# Step 2: schema[gene/protein, associated with] = {disease, effect/phenotype}
# Expected: disease or effect/phenotype entities only
```

### Result Filtering

After model inference:
1. **Load expected tail types** using schema
2. **Filter predictions** by entity type:
   ```python
   if entity_type in expected_tail_types:
       entity["schema_match"] = True
   else:
       entity["schema_match"] = False
   ```
3. **Save both versions**:
   - `predictions_all.parquet`: All predictions with `schema_match` flag
   - `predictions_filtered.parquet`: Only schema-valid predictions
4. **Return filtered results** in API response

### Schema File Management

**Location:** `primekg_schema.json` (auto-generated from PrimeKG data)

**Contents:**
- 45 valid `(head_type, relation)` → `tail_type(s)` mappings
- Complete coverage of all 8.1M edges in PrimeKG
- Entity type and relation statistics

**Regenerate schema:**
```bash
pixi run python build_schema.py
```

**When to regenerate:**
- After updating PrimeKG data
- When adding new entity types or relations
- When setting up a new environment (first time)

**Performance:** Schema loading is ~100x faster than dynamic computation (milliseconds vs seconds)

## Troubleshooting

### Error: "Schema file not found"

**Solution:**
```bash
pixi run python build_schema.py
```

### Error: "Checkpoint not found"

Ensure checkpoint files exist at the specified paths. You can:
1. Train your own checkpoints using `script/run_query.py`
2. Download pre-trained checkpoints (if available)
3. Use only the Ultra checkpoint (UltraQuery checkpoint is optional)

### Error: "Entity/Relation not found"

- Check entity ID format: Should be `SOURCE:ID` (e.g., `MONDO:5301`, `NCBI:9796`)
- Check relation labels: Use `ppi`, `associated with`, etc.
- The server provides suggestions for similar entities/relations

### Error: "Invalid query" with schema validation

**Example error:**
```
Invalid query: No edges in graph schema for (drug, expression present) (entity: Aspirin).
This combination does not exist in PrimeKG.
```

**Solution:** Check the schema reference table above for valid combinations. Common issues:
- Wrong entity type for relation (e.g., drugs don't have `expression present`)
- Empty intersection (e.g., trying to intersect incompatible entity types)
- Invalid multi-hop path (intermediate type doesn't support next relation)

### GPU Out of Memory

The server automatically falls back to CPU if GPU OOM occurs. You can also:
- Set `device: "cpu"` in config
- Reduce model size by using smaller checkpoints

## Performance Notes

- **First request**: Model loading takes 10-30 seconds
- **Subsequent requests**: < 2 seconds per query
- **Model caching**: Model is loaded once and reused (singleton pattern)
- **Schema loading**: Cached after first load (~10-50 ms initial load)
- **Batch processing**: Currently single-query mode (future: batch support)

## Differences from `run_query.py`

| Feature | run_query.py | ultraquery_inference_server.py |
|---------|--------------|-------------------------------|
| **Purpose** | Training + evaluation | Inference only |
| **Input** | Pre-built query datasets | Custom queries from users |
| **Output** | Metrics (MRR, AUROC) | Ranked entity predictions |
| **Mode** | Batch processing | On-demand API |
| **Interface** | CLI script | MCP server (HTTP) |
| **Training** | ✅ Supports training | ❌ Inference only |
| **Schema Validation** | ❌ No validation | ✅ Pre-inference validation |

## Advanced Usage

### Inspecting the Schema

```python
import json

# Load schema
with open("primekg_schema.json", "r") as f:
    schema_data = json.load(f)

# Check valid tail types for a head-relation pair
key = "gene/protein|associated with"
tail_types = schema_data["schema"][key]
print(f"Valid tail types: {tail_types}")
# Output: ['disease', 'effect/phenotype']

# List all relations
for relation, count in schema_data["relation_counts"].items():
    print(f"{relation:30s} {count:>10,} edges")
```

### Custom Schema Generation

If you have a modified PrimeKG dataset, regenerate the schema:

```bash
# Ensure your dataset is in data/primekg1/
pixi run python build_schema.py

# The script will:
# 1. Load all edges from train + valid + test
# 2. Build complete schema mappings
# 3. Save to primekg_schema.json
# 4. Print comprehensive statistics
```

## Citation

If you use UltraQuery, please cite:

```bibtex
@article{galkin2023ultra,
  title={Towards Foundation Models for Knowledge Graph Reasoning},
  author={Galkin, Mikhail and Trivedi, Rakshit and Maheshwari, Gaurav and Zhu, Renjie and Rasul, Kashif and Xu, Yujia and Jiang, Yali and Yuan, Jinwen and Wu, Linhao and Wu, Yifan and others},
  journal={arXiv preprint arXiv:2310.04562},
  year={2023}
}
```

## References

- **UltraQuery Paper**: [Complex Logical Query Answering on Knowledge Graphs](https://arxiv.org/abs/2301.02334)
- **PrimeKG**: [Precision Medicine Knowledge Graph](https://github.com/mims-harvard/PrimeKG)
- **ULTRA**: [Towards Foundation Models for Knowledge Graph Reasoning](https://arxiv.org/abs/2310.04562)

## License

Same license as the parent ULTRA project.
