# Agent Skills & MCP Servers

A collection of self-developed Claude Code skills and Model Context Protocol (MCP) servers for biomedical data analysis and computational biology workflows.

## Overview

This repository contains custom-built tools that extend Claude's capabilities for biomedical research, enabling seamless integration with biological databases, knowledge graphs, and analysis pipelines.

## Repository Structure

```
agent-skill-mcp/
├── skill/           # Claude Code skills for specialized workflows
└── mcp/            # MCP servers for external service integration
```

## Skills

Claude Code skills provide specialized workflows and domain knowledge for specific tasks:

### 1. **BioBridge** (`skill/biobridge/`)
Predict biomedical entity associations using the BioBridge multimodal knowledge graph. Supports:
- Link prediction between genes, diseases, drugs, and phenotypes
- Entity matching with user confirmation
- Association scoring using neural retrieval
- Results export and summarization

### 2. **CELLxGENE Census** (`skill/cellxgene/`)
Comprehensive toolkit for working with single-cell RNA-seq data from CELLxGENE Census:
- **Query subskill**: Download and filter expression data with flexible metadata filtering
- **Specificity subskill**: Extract cell type marker genes and analyze expression specificity
- Gene coexpression analysis with visualization

### 3. **Disease Module Analysis** (`skill/disease-module-analysis/`)
Specialized workflows for disease module analysis using single-cell RNA-seq data:
- Disease module training and transformation
- Gene set enrichment analysis (pathway/ontology)
- Perturbation enrichment (TF/Gene/LINCS)
- Regulator identification
- NCBI metadata extraction

### 4. **Protein Interaction Query** (`skill/interactdb-query/`)
Query protein-protein interaction databases (STRING, IntAct, BioGRID):
- Single-gene neighbor analysis
- Multi-gene shortest path finding
- Comprehensive filtering with entity and edge annotations
- Network visualization preparation

### 5. **Pathway Database Query** (`skill/pathwaydb-query/`)
Query pathway databases (KEGG, Reactome, MSigDB):
- Find gene-associated pathways and terms
- Multi-gene comparative analysis with UpSet plots
- Pathway similarity network construction using Jaccard indices
- Gene-specific subnetwork extraction with centrality analysis

### 6. **Omicsoft DEG Analysis** (`skill/omicsoft-deg-analysis/`)
Analyze pre-computed differential gene expression statistics:
- Query bulk DEG data from h5ad files
- Filter by disease, study, or custom criteria
- Signature-based expression pattern analysis
- GSEA enrichment with leading edge annotation

## MCP Servers

MCP servers enable Claude to interact with external services and data sources:

### **BioBridge MCP** (`mcp/biobridge/`)
Model Context Protocol server for the BioBridge multimodal knowledge graph:
- `predict_associations` tool for inferring biomedical links
- Neural retrieval over pre-trained embeddings
- Support for gene-disease, drug-phenotype, and other entity associations
- Flexible entity/relation inputs with cosine similarity scoring

## Getting Started

Each skill and MCP server contains its own documentation:
- Skills: See `SKILL.md` in each skill directory
- MCP servers: See implementation files and inline documentation

## Contact

**Maintainer:** Xinghao Zhang
**Email:** xinghao.zhang@takeda.com

## Contributing

This repository is maintained for internal use. For questions, issues, or contributions, please contact the maintainer.

## License

Internal use only. Please contact the maintainer for licensing information.
