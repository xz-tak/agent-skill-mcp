# CyteType API Reference

## Overview

CyteType is an automated cell type annotation tool for single-cell RNA-seq data using multi-agent AI architecture.

**Installation:**
```bash
pip install cytetype
```

**GitHub:** https://github.com/NygenAnalytics/CyteType

## CyteType Class

### Initialization

```python
from cytetype import CyteType

annotator = CyteType(
    adata,
    group_key="leiden",
    rank_key="rank_genes_groups",
    gene_symbols_column="gene_symbols",
    n_top_genes=50,
    aggregate_metadata=True,
    min_percentage=10,
    pcent_batch_size=2000,
    coordinates_key="X_umap",
    max_cells_per_group=1000,
)
```

### Initialization Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | `anndata.AnnData` | Required | AnnData object with log1p-normalized data |
| `group_key` | `str` | Required | Column in `adata.obs` with cluster labels |
| `rank_key` | `str` | `"rank_genes_groups"` | Key in `adata.uns` with DE results |
| `gene_symbols_column` | `str` | `"gene_symbols"` | Column in `adata.var` with gene symbols |
| `n_top_genes` | `int` | `50` | Top marker genes per cluster |
| `aggregate_metadata` | `bool` | `True` | Aggregate metadata from AnnData |
| `min_percentage` | `int` | `10` | Min percentage for cluster context |
| `pcent_batch_size` | `int` | `2000` | Batch size for expression % calculation |
| `coordinates_key` | `str` | `"X_umap"` | Coordinates key in `adata.obsm` |
| `max_cells_per_group` | `int` | `1000` | Max cells per group for visualization |

### run() Method

```python
adata = annotator.run(
    study_context="Human PBMC from healthy donor",
    llm_configs=[{
        "provider": "openai",
        "name": "gpt-5.2",
        "apiKey": "your-api-key",
        "baseUrl": "https://api.openai.com/v1",
        "modelSettings": {
            "temperature": 0.0,
            "max_tokens": 4096
        }
    }],
    metadata={"experiment": "PBMC_study"},
    n_parallel_clusters=2,
    results_prefix="cytetype",
    poll_interval_seconds=10,
    timeout_seconds=7200,
    save_query=True,
    query_filename="query.json",
    auth_token=None,
    show_progress=True,
)
```

### run() Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `study_context` | `str` | Required | Biological context (organism, tissue, disease) |
| `llm_configs` | `list[dict]` | `None` | LLM configuration (see below) |
| `metadata` | `dict` | `None` | Custom metadata for report header |
| `n_parallel_clusters` | `int` | `2` | Parallel requests (max 50) |
| `results_prefix` | `str` | `"cytetype"` | Prefix for result columns |
| `poll_interval_seconds` | `int` | `10` | Polling interval |
| `timeout_seconds` | `int` | `7200` | Max wait time (2 hours) |
| `api_url` | `str` | Default API | Custom API endpoint |
| `save_query` | `bool` | `True` | Save query JSON |
| `query_filename` | `str` | `"query.json"` | Query filename |
| `auth_token` | `str` | `None` | Bearer auth token |
| `show_progress` | `bool` | `True` | Show progress updates |

## LLM Configuration

### Supported Providers

- `anthropic` - Anthropic Claude models
- `bedrock` - AWS Bedrock
- `google` - Google AI (Gemini)
- `groq` - Groq
- `mistral` - Mistral AI
- `openai` - OpenAI (default)
- `openrouter` - OpenRouter

### Configuration Schema

```python
llm_config = {
    "provider": "openai",           # Required: provider name
    "name": "gpt-5.2",              # Required: model name
    "apiKey": "sk-...",             # Required for non-Bedrock
    "baseUrl": "https://...",       # Optional: custom endpoint
    "modelSettings": {              # Optional: model parameters
        "temperature": 0.0,
        "max_tokens": 4096,
    }
}
```

### AWS Bedrock Configuration

```python
llm_config = {
    "provider": "bedrock",
    "name": "anthropic.claude-3-sonnet",
    "awsAccessKeyId": "AKIA...",
    "awsSecretAccessKey": "...",
    "awsDefaultRegion": "us-east-1",
}
```

## Output Schema

### AnnData Columns Added

After annotation, these columns are added to `adata.obs`:

| Column | Description |
|--------|-------------|
| `{prefix}_annotation_{group_key}` | Primary cell type annotation |
| `{prefix}_cellOntologyTerm_{group_key}` | Cell Ontology term |
| `{prefix}_cellOntologyTermID_{group_key}` | Cell Ontology ID (CL:XXXXXXX) |
| `{prefix}_cellState_{group_key}` | Cell state/activation |

### Results Structure

Results stored in `adata.uns["{prefix}_results"]`:

```python
{
    "job_id": "uuid-string",
    "result": {
        "annotations": [
            {
                "clusterId": "1",
                "annotation": "CD4+ T cell",
                "granularAnnotation": "Naive CD4+ T cell",
                "cellState": "Resting",
                "ontologyTerm": "CD4-positive, alpha-beta T cell",
                "ontologyTermID": "CL:0000624",
                "confidence": 0.92,
                "supportingMarkers": ["CD4", "TCF7", "CCR7"],
                "conflictingMarkers": [],
                "missingExpression": [],
                "unexpectedExpression": [],
                "corroboratingPapers": [
                    {"title": "...", "pmid": "..."}
                ]
            },
            # ... more annotations
        ]
    }
}
```

### Annotation Fields

| Field | Type | Description |
|-------|------|-------------|
| `clusterId` | `str` | Cluster identifier (1-indexed) |
| `annotation` | `str` | Primary cell type |
| `granularAnnotation` | `str` | More specific phenotype |
| `cellState` | `str` | Activation/state |
| `ontologyTerm` | `str` | Cell Ontology term name |
| `ontologyTermID` | `str` | Cell Ontology ID |
| `confidence` | `float` | Confidence score (0-1) |
| `supportingMarkers` | `list[str]` | Markers supporting annotation |
| `conflictingMarkers` | `list[str]` | Markers conflicting with annotation |
| `missingExpression` | `list[str]` | Expected markers not expressed |
| `unexpectedExpression` | `list[str]` | Unexpected markers expressed |
| `corroboratingPapers` | `list[dict]` | Supporting literature |

## get_results() Method

Retrieve results with automatic API fallback:

```python
results = annotator.get_results()

# Or from adata.uns
import json
results = json.loads(adata.uns["cytetype_results"]["result"])

# Access annotations
for ann in results["annotations"]:
    print(f"Cluster {ann['clusterId']}: {ann['annotation']} ({ann['confidence']:.2f})")
```

## Prerequisites

### AnnData Requirements

1. **Normalized expression**: `adata.X` should contain log1p-normalized values
2. **Cluster labels**: Column in `adata.obs` with cluster assignments
3. **DE results**: Run `sc.tl.rank_genes_groups()` first

```python
import scanpy as sc

# Ensure data is preprocessed
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, key_added='cluster')

# REQUIRED: Run DE analysis
sc.tl.rank_genes_groups(adata, groupby='cluster')

# Optional: UMAP for visualization
sc.tl.umap(adata)
```

### Gene Symbol Handling

If `adata.var_names` contains gene symbols (not Ensembl IDs):

```python
# Copy var_names to gene_symbols column
adata.var['gene_symbols'] = adata.var_names

# Or specify different column
annotator = CyteType(adata, group_key='cluster', gene_symbols_column='gene_name')
```

## Example Workflow

```python
import scanpy as sc
from cytetype import CyteType
from dotenv import load_dotenv
import os

# Load API key from .env
load_dotenv()
api_key = os.getenv("OPENAI_API_KEY")
base_url = os.getenv("OPENAI_API_BASE")

# Load and preprocess data
adata = sc.read_h5ad("pbmc.h5ad")

# Initialize annotator
annotator = CyteType(
    adata,
    group_key='cluster',
    n_top_genes=50,
)

# Configure LLM
llm_config = [{
    "provider": "openai",
    "name": "gpt-5.2",
    "apiKey": api_key,
    "baseUrl": base_url,
    "modelSettings": {"temperature": 0.0}
}]

# Run annotation
adata = annotator.run(
    study_context="Human PBMC from healthy donor",
    llm_configs=llm_config,
)

# View results
print(adata.obs['cytetype_annotation_cluster'].value_counts())

# Plot
sc.pl.umap(adata, color='cytetype_annotation_cluster')

# Save
adata.write_h5ad("pbmc_annotated.h5ad")
```

## Error Handling

### Common Errors

**"rank_genes_groups not found"**
```python
# Solution: Run DE analysis first
sc.tl.rank_genes_groups(adata, groupby='cluster')
```

**"API key not found"**
```python
# Solution: Set environment variable or pass directly
llm_config = [{"provider": "openai", "name": "gpt-5.2", "apiKey": "sk-..."}]
```

**"Timeout exceeded"**
```python
# Solution: Increase timeout
adata = annotator.run(study_context="...", timeout_seconds=14400)  # 4 hours
```

## Resources

- **Documentation**: https://github.com/NygenAnalytics/CyteType/tree/master/docs
- **Colab Tutorial**: https://colab.research.google.com/drive/1aRLsI3mx8JR8u5BKHs48YUbLsqRsh2N7
- **Discord Support**: https://discord.gg/V6QFM4AN
- **R Wrapper (CyteTypeR)**: https://github.com/NygenAnalytics/CyteTypeR
