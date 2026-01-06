# AI-Sci MCP Services

Centralized repository of Model Context Protocol (MCP) servers for the AI-Sci team, providing AI agents with access to biomedical knowledge graphs and machine learning models.

## Servers

### BioBridge (`biobridge/`)

A multimodal knowledge graph MCP server for predicting biomedical associations using neural retrieval over BioBridge embeddings. Predicts associations between biomedical entities (gene-disease, drug-phenotype, etc.) using pre-trained embeddings and LLM-assisted entity selection. Supports multiple entity types including genes/proteins, drugs, diseases, pathways, phenotypes, anatomy, and exposures.

### PrimeKG (`primekg/`)

A comprehensive query interface for the PrimeKG biomedical knowledge graph containing 129,000+ entities and 8M+ relationships. Enables querying and exploration of the knowledge graph through entity search, connection retrieval, and multi-hop neighborhood queries. See [primekg/README.md](primekg/README.md) for detailed documentation.

### ULTRA (`ultra/`)

(Deprecated and will be removed in the future). Query interface for ULTRA knowledge graph foundation model predictions on PrimeKG data stored in S3. Retrieves pre-computed model predictions for concept-relation pairs and ranks specific triples using DuckDB-powered S3 queries over parquet files. 

### ULTRA Inference (`ultra-inference/`)

On-demand inference MCP server for the ULTRA foundation model on PrimeKG biomedical knowledge graph. Provides zero-shot link prediction to predict tail entities for any (head_entity, relation) query without requiring pre-computed results. Includes data setup tools, input validation, GPU support with CPU fallback, and full data transformation pipeline with novelty detection and schema validation. Model loads lazily on first request and stays in memory for fast subsequent queries.

### UltraQuery Inference (`ultraquery-inference/`)

Inference-only MCP server for answering complex logical queries on PrimeKG using the UltraQuery model. Supports multi-hop reasoning, intersection, union, and negation operations to answer queries like "What diseases are associated with proteins that interact with both GREM1 and IL11?" Uses graph schema validation for query validation and result filtering across 10 entity types and 18 relation types.

## Repository Structure

```
ai-sci-mcp-services/
├── biobridge/                 # BioBridge MCP server
│   ├── biobridge_mcp_server.py
│   └── src/                   # Model components
├── primekg/                   # PrimeKG MCP server
│   ├── primekg_server.py
│   └── README.md              # Detailed documentation
├── ultra/                     # ULTRA MCP server (pre-computed)
│   ├── ultra_results_server.py
│   └── similarity_search.py
├── ultra-inference/           # ULTRA Inference MCP server (on-demand)
│   ├── ultra_inference_server.py
│   ├── primekg_setup.py
│   ├── src/ultra/             # ULTRA model implementation
│   └── README.md              # Detailed documentation
├── ultraquery-inference/      # UltraQuery Inference MCP server
│   ├── ultraquery_inference_server.py
│   ├── build_schema.py
│   ├── src/ultra/             # ULTRA model implementation
│   └── README.md              # Detailed documentation
└── README.md                  # This file
```

Each server directory contains its own `pixi.toml` for dependency management and can be run independently. See individual server directories for installation and usage instructions.
