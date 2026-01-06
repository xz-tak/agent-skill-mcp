# ULTRA Inference MCP Server

On-demand inference for the ULTRA foundation model on PrimeKG biomedical knowledge graph via Model Context Protocol (MCP).

## Overview

This is a **standalone** MCP server that provides zero-shot link prediction on PrimeKG biomedical knowledge graphs using the ULTRA foundation model. All ULTRA dependencies are self-contained in the `src/` directory, with data and outputs stored in the local `data/` and `output/` directories, making this server fully independent and deployable separately from the main ULTRA repository.

## Key Features

- **Standalone**: All ULTRA dependencies contained in `src/ultra/` - no parent repository required
- **Data Setup**: Download and process PrimeKG data with train/test/valid splits
- **Zero-shot Inference**: Predict tail entities for any (head, relation) query without pre-training
- **On-demand Inference**: Query any (head_entity, relation) pair without pre-computed results
- **Lazy Loading**: Model loads on first request and stays in memory for fast subsequent queries
- **Input Validation**: Validates entities and relations against PrimeKG with helpful error messages
- **Data Transformation**: Full pipeline including ID translation, novelty detection, and schema validation
- **GPU Support**: Automatic GPU detection with CPU fallback on OOM
- **PrimeKG Integration**: Works seamlessly with PrimeKG biomedical knowledge graph data

## Directory Structure

```
ultra-inference/
├── ultra_inference_server.py   # Main MCP server implementation
├── primekg_setup.py             # PrimeKG data download and processing
├── .gitignore                   # Git ignore file (excludes data/ and output/)
├── src/                         # Standalone ultra package dependencies
│   └── ultra/                   # ULTRA model implementation
│       ├── __init__.py
│       ├── models.py            # Ultra, RelNBFNet, EntityNBFNet classes
│       ├── util.py              # Utility functions
│       ├── tasks.py             # Training/inference tasks
│       ├── data_util.py         # PrimeKG data processing
│       ├── datasets.py          # Dataset classes (PrimeKG1-5)
│       ├── base_nbfnet.py       # Base NBFNet architecture
│       ├── layers.py            # GNN layers
│       └── rspmm/               # Custom CUDA kernel for efficient message passing
│           ├── __init__.py
│           ├── rspmm.py
│           └── source/
│               ├── rspmm.cpp
│               ├── rspmm.cu
│               └── rspmm.h
├── ckpts/                       # Model checkpoints directory
│   ├── README.md                # Instructions for obtaining checkpoints
│   └── *.pth                    # Checkpoint files (gitignored)
├── data/                        # Self-contained data directory (gitignored)
│   ├── primekg.csv              # Downloaded from Harvard Dataverse
│   ├── *.pkl                    # Auto-generated dictionary files
│   └── primekg1/                # Dataset-specific files
│       ├── raw/                 # Train/test/valid splits + nodes
│       └── processed/           # PyTorch Geometric processed data
├── output/                      # Inference results (gitignored)
│   └── {entity_id}/
│       └── {relation}/
│           └── predictions.parquet
├── pixi.toml                    # Pixi environment configuration
└── README.md                    # This file
```

## Installation

### 1. Install Dependencies

```bash
cd ultra-inference
pixi install
```

The environment includes:
- PyTorch and PyTorch Geometric
- Polars for data processing
- FastMCP for MCP server implementation
- All ULTRA model dependencies

### 2. Model Checkpoint

**The checkpoint downloads automatically** on first inference request. No manual setup required!

## Usage

### Start the MCP Server

```bash
cd ultra-inference
pixi run python ultra_inference_server.py
```

## MCP Tools

The server exposes two main tools:

### 1. setup_primekg_data

Download and setup PrimeKG data with train/test/valid splits.

**Function Signature:**

```python
setup_primekg_data(
    dataset_name: str = "PrimeKG1",  # Name of dataset (e.g., "PrimeKG1", "PrimeKG2")
    force_redownload: bool = False,  # If True, re-download even if files exist
    train_frac: float = 0.8,         # Training set fraction
    test_frac: float = 0.1,          # Test set fraction
    valid_frac: float = 0.1,         # Validation set fraction
    seed: int = 42                   # Random seed for reproducibility
)
```

**What it does:**
1. Checks if PrimeKG data already exists
2. Downloads primekg.csv from Harvard Dataverse if needed
3. Processes the data into train/test/valid splits
4. Creates nodes file with entity metadata

**Example via MCP client:**

```python
result = mcp_client.call_tool(
    "setup_primekg_data",
    dataset_name="PrimeKG1",
    force_redownload=False
)

if result["success"]:
    print(f"Status: {result['status']}")
    print(f"Nodes: {result['nodes_count']:,}")
    print(f"Train edges: {result['train_edges']:,}")
    print(f"Test edges: {result['test_edges']:,}")
    print(f"Valid edges: {result['valid_edges']:,}")
else:
    print(f"Error: {result['error']}")
```

**Return Format:**

```python
{
    "success": True,
    "status": "completed",  # or "already_exists"
    "nodes_count": 129375,
    "train_edges": 6553281,
    "test_edges": 819160,
    "valid_edges": 819160,
    "dataset_path": "/path/to/data/primekg1"
}
```

### 2. predict_tail_entities

Predict tail entities for a given head entity and relation.

**Function Signature:**

```python
predict_tail_entities(
    head_entity: str,      # Entity ID (e.g., "NCBI:7297") or name (e.g., "TYK2")
    relation: str,         # Relation label (e.g., "associated_with", "associated with")
    top_k: int = None      # Number of predictions to return (None = all predictions)
)
```

**Parameters:**
- `head_entity` (required): Entity ID (e.g., "MONDO:5301") or name (e.g., "Crohn disease")
- `relation` (required): Relation label (e.g., "associated_with", "associated with")
- `top_k` (optional, default=None): Number of top predictions to return (None returns all predictions)

**Example via MCP client:**

```python
result = mcp_client.call_tool(
    "predict_tail_entities",
    head_entity="MONDO:5301",
    relation="associated_with",
    top_k=50
)

if result["success"]:
    print(f"Predictions saved to: {result['output_file']}")
    print(f"Total predictions: {result['total_predictions']}")
    print(f"Inference time: {result['inference_time_seconds']}s")
    print(f"\nTop predictions:\n{result['preview']}")
else:
    print(f"Error: {result['error']}")
```

**Return Format:**

```python
{
    "success": True,
    "head_entity": "MONDO:5301",
    "head_name": "Crohn disease",
    "relation": "associated_with",
    "output_file": "/path/to/output/MONDO_5301/associated_with/predictions.parquet",
    "total_predictions": 100,
    "inference_time_seconds": 2.34,
    "preview": "shape: (10, 11)\n┌───────────┬────────┬─────────────┬─────────────┬───┬────────────┬────────────┬──────┬────────────┐\n│ h_label ┆ h_name ┆ h_type ┆ r_label ┆ … ┆ t_pred_typ ┆ edge_in_pr ┆ rank ┆ percentile │\n│ --- ┆ --- ┆ --- ┆ --- ┆ ┆ e ┆ imekg ┆ --- ┆ _rank │\n│ str ┆ str ┆ str ┆ str ┆ ┆ --- ┆ --- ┆ i64 ┆ --- │\n│ ┆ ┆ ┆ ┆ ┆ str ┆ bool ┆ ┆ f64 │\n╞═══════════╪════════╪═════════════╪═════════════╪═══╪════════════╪════════════╪══════╪════════════╡\n..."
}
```

## Output Format

Results are saved as Parquet files with the following columns:

- `h_label`: Head entity ID
- `h_name`: Head entity name
- `h_type`: Head entity type
- `r_label`: Relation label
- `t_pred_label`: Predicted tail entity ID
- `t_pred_name`: Predicted tail entity name
- `t_pred_score`: Prediction score (cosine similarity)
- `t_pred_type`: Predicted tail entity type
- `edge_in_primekg`: Boolean - whether this edge already exists in PrimeKG (novelty detection)
- `rank`: Rank of this prediction (1 = best)
- `percentile_rank`: Percentile rank (0.99 = 99th percentile)

## Configuration

### Environment Variables

- `ULTRA_OUTPUT_DIR`: Override output directory (default: `./output`)

**Example:**
```bash
export ULTRA_OUTPUT_DIR=/path/to/output
pixi run python ultra_inference_server.py
```

**Note:** All paths are self-contained within the MCP directory. Checkpoints download automatically on first use.

### Model Configuration

The server uses the following default model configuration:

- **Dataset**: PrimeKG1
- **Dataset Root**: `./data/`
- **Model**: ULTRA with RelNBFNet + EntityNBFNet
- **Checkpoint**: `./ckpts/ultra_primekg_50g_ft_epoch_1.pth` (fine-tuned on PrimeKG)
- **Output Directory**: `./output/`
- **Message Function**: DistMult
- **Hidden Dimensions**: [64, 64, 64, 64, 64, 64]
- **Device**: Auto-detected (CUDA if available, else CPU)
- **Batch Size**: 4 (configurable in code)

All paths are self-contained within the MCP directory for easy deployment.

## Architecture

The server consists of four main components:

1. **ModelManager**: Lazy-loaded singleton for model, dataset, and ID dictionaries
2. **Validation Layer**: Validates entities and relations with fuzzy matching and helpful suggestions
3. **Inference Engine**: Executes ULTRA model on single queries with strict negative masking
4. **Transformation Pipeline**: Applies `translate_hrt`, `filter_process_results`, `structure_results`, and adds percentile ranks

### ULTRA Model

ULTRA uses a dual-graph architecture:
- **RelNBFNet**: Operates on the relation graph to learn relation representations
- **EntityNBFNet**: Performs message passing on the entity graph using relation representations
- **rspmm Kernel**: Custom CUDA kernel for efficient relational message passing with O(V) instead of O(E) complexity

## Performance

- **First query**: ~5-10 seconds (includes model loading)
- **Subsequent queries**: ~1-3 seconds (model stays in memory)
- **GPU**: Significantly faster than CPU for inference
- **CUDA kernel**: Compiles automatically on first use (requires CUDA_HOME)

## Error Handling

The server provides helpful error messages for common issues:

### Invalid Entity

```
Entity 'XYZ' not found in PrimeKG. Did you mean:
  MONDO:5301 (Crohn disease),
  MONDO:5011 (Crohn's disease of ileum)
```

### Invalid Relation

```
Relation 'treats' not found. Available relations include:
  associated with, drug_protein, contraindication, indication, ...
```

### GPU Out of Memory

The server automatically falls back to CPU if GPU memory is insufficient:
```
GPU out of memory, falling back to CPU
```

### Empty Predictions

If no valid predictions remain after schema filtering, the server returns a warning message with details about filtering.

## Troubleshooting

### Checkpoint download fails

If automatic checkpoint download fails:
- Check your internet connection
- Verify GitHub access (https://github.com/roger-tu/ULTRA)
- Check the logs for specific error messages

The checkpoint is downloaded from:
https://raw.githubusercontent.com/roger-tu/ULTRA/b149f9d42921047b58475d2d18929d864a2321a7/ckpts/ultra_primekg_50g_ft_epoch_1.pth

If download continues to fail, check firewall/proxy settings.

### Import errors

Make sure you're running from the mcp directory so the path manipulation works:
```bash
cd ultra-inference
pixi run python ultra_inference_server.py
```

The server uses `sys.path.insert(0, "src")` to import from the standalone `src/ultra/` directory.

### CUDA errors

Set environment variable to force CPU:
```bash
export CUDA_VISIBLE_DEVICES=""
pixi run python ultra_inference_server.py
```

### rspmm compilation fails

The rspmm CUDA kernel requires:
- CUDA toolkit installed
- CUDA_HOME environment variable set
- Ninja build system

If compilation fails, the model will fall back to slower PyTorch operations.

## Integration

This standalone MCP server can be integrated into:
- **Claude Code**: Add to MCP configuration
- **Other MCP clients**: Point to this directory in MCP settings

The server is fully self-contained and can be deployed independently of the main ULTRA repository.

## License

See parent ULTRA project license.
