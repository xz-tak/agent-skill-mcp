# PrimeKG MCP Server

A Model Context Protocol (MCP) server for querying and exploring the PrimeKG biomedical knowledge graph.

## Overview

This MCP server provides a set of tools for AI agents to interact with PrimeKG, a comprehensive biomedical knowledge graph containing over 129,000 entities and 8 million relationships across multiple biomedical domains.

## Features

The server provides the following tools:

### 1. Metadata Queries
- `get_entity_types()` - Get all unique entity types and their counts
- `get_relation_types()` - Get all relation types and their counts
- `get_data_sources()` - Get all data sources/resources included in PrimeKG
- `get_statistics()` - Get comprehensive dataset statistics

### 2. Entity Operations
- `check_entity_exists(entity_id)` - Check if an entity exists and get its details
- `search_entities(search_term, entity_type?, source?, limit?)` - Search for entities by name
- `get_entity_connections(entity_id, relation_type?, limit?)` - Get all connections for an entity
- `get_entity_neighborhoods(entity_id, max_depth?, relation_types?, max_neighbors_per_level?)` - Get multi-hop neighborhood

### 3. Relation Operations
- `check_relation_exists(relation_type)` - Check if a relation type exists

## Installation

1. Ensure you have the required dependencies:
```bash
cd /path/to/PrimeKG-Ultra/mcp
pixi install
# Or manually:
pip install fastmcp polars
```

2. Set the data directory (optional, defaults to ../data):
```bash
export PRIMEKG_DATA_DIR=/path/to/PrimeKG-Ultra/data
```

## Usage

### Running the Server

```bash
cd /path/to/PrimeKG-Ultra/mcp
python primekg_server.py
```

Or use it as a uvx tool:
```bash
uvx --from . primekg-server
```

### Example Queries

#### Get Entity Types
```python
# Returns all entity types with counts
get_entity_types()
# Example output:
# {
#   "entity_types": ["gene/protein", "drug", "disease", ...],
#   "count_by_type": {"gene/protein": 34567, "drug": 8765, ...},
#   "total_entities": 129375
# }
```

#### Check Entity Exists
```python
check_entity_exists("NCBI:9796")
# Returns entity details including name, type, source
```

#### Search for Entities
```python
search_entities("insulin", entity_type="gene/protein", limit=10)
# Returns up to 10 gene/protein entities matching "insulin"
```

#### Get Entity Connections
```python
get_entity_connections("NCBI:9796", relation_type="protein_protein", limit=50)
# Returns up to 50 protein-protein interactions for this gene
```

#### Get Multi-hop Neighborhood
```python
get_entity_neighborhoods("DRUGBANK:DB00001", max_depth=2, max_neighbors_per_level=10)
# Returns 1-hop and 2-hop neighbors of the drug entity
```

## Data Format

### Entity ID Format
Entities are identified by: `SOURCE:ID`
- Examples: `NCBI:9796`, `DRUGBANK:DB00001`, `MONDO:0005015`

### Entity Types
Common entity types include:
- `gene/protein` - Genes and proteins
- `drug` - Drugs and compounds
- `disease` - Diseases and disorders
- `biological_process` - GO biological processes
- `molecular_function` - GO molecular functions
- `cellular_component` - GO cellular components
- `pathway` - Biological pathways
- `anatomy` - Anatomical structures
- `phenotype` - Phenotypic traits
- `exposure` - Environmental exposures

### Relation Types
Common relation types include:
- `protein_protein` (ppi) - Protein-protein interactions
- `drug_drug` - Drug-drug interactions
- `disease_protein` - Disease-protein associations
- `drug_protein` - Drug-protein targets
- `drug_disease` - Drug-disease indications
- And many more...

## Configuration

### Environment Variables
- `PRIMEKG_DATA_DIR` - Path to PrimeKG data directory (default: ../data)

### Data Files Required
The server expects the following files in the data directory:
- `nodes.txt` - Node information (TSV format)
- `primekg.csv` - Complete edge list with metadata
- `ent2name_dict.pkl` (optional) - Entity to human-readable name mappings

## Performance Notes

- **First Query**: The first query will load the entire dataset into memory (~8M edges, ~130K nodes). This takes 10-30 seconds depending on your system.
- **Subsequent Queries**: All subsequent queries are served from memory and are very fast.
- **Memory Usage**: The server requires approximately 2-3 GB of RAM to hold the dataset in memory.
- **Neighborhood Queries**: Multi-hop queries are limited to depth 2 to avoid exponential explosion. Use `max_neighbors_per_level` to control result size.

## Integration with Claude Code

To use this MCP server with Claude Code:

1. Add to your MCP configuration file (e.g., `~/.config/claude-code/mcp.json`):
```json
{
  "mcpServers": {
    "primekg": {
      "command": "python",
      "args": ["/path/to/PrimeKG-Ultra/mcp/primekg_server.py"],
      "env": {
        "PRIMEKG_DATA_DIR": "/path/to/PrimeKG-Ultra/data"
      }
    }
  }
}
```

2. Restart Claude Code

3. The PrimeKG tools will now be available in your Claude Code sessions

## Example Use Cases

### 1. Explore Drug Targets
```python
# Find a drug
drug = search_entities("aspirin", entity_type="drug", limit=1)
drug_id = drug["results"][0]["entity_id"]

# Get its protein targets
targets = get_entity_connections(drug_id, relation_type="drug_protein", limit=20)
```

### 2. Investigate Disease-Gene Associations
```python
# Find a disease
disease = search_entities("diabetes", entity_type="disease", limit=1)
disease_id = disease["results"][0]["entity_id"]

# Get associated genes/proteins
genes = get_entity_connections(disease_id, relation_type="disease_protein", limit=50)
```

### 3. Analyze Protein Interaction Networks
```python
# Get a protein and its interaction partners
protein_network = get_entity_neighborhoods(
    "NCBI:9796",
    max_depth=2,
    relation_types=["protein_protein"],
    max_neighbors_per_level=10
)
```

## Troubleshooting

### "File not found" errors
- Check that `PRIMEKG_DATA_DIR` points to the correct directory
- Ensure `nodes.txt` and `primekg.csv` exist in the data directory

### Slow performance
- First query will always be slow (data loading)
- Reduce `limit` parameters in queries
- For neighborhood queries, reduce `max_depth` or `max_neighbors_per_level`

### Out of memory errors
- The dataset is large. Ensure you have at least 4 GB of free RAM
- Consider processing queries in smaller batches

## License

This MCP server is part of the PrimeKG-Ultra project. See the main repository for license information.

## Citation

If you use PrimeKG in your research, please cite:

```bibtex
@article{chandak2023building,
  title={Building a knowledge graph to enable precision medicine},
  author={Chandak, Payal and Huang, Kexin and Zitnik, Marinka},
  journal={Scientific Data},
  volume={10},
  number={1},
  pages={67},
  year={2023},
  publisher={Nature Publishing Group UK London}
}
```
