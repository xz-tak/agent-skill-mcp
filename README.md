# Agent Skills & MCP Servers

A collection of self-developed Claude Code skills and Model Context Protocol (MCP) servers for biomedical data analysis and computational biology workflows.

## Overview

This repository contains custom-built tools that extend Claude's capabilities for biomedical research, enabling seamless integration with biological databases, knowledge graphs, and analysis pipelines.

## Repository Structure

```
agent-skill-mcp/
├── skills/          # Claude Code skills for specialized workflows
└── mcp/             # MCP servers for external service integration
```

## Skills

Claude Code skills provide specialized workflows and domain knowledge for specific tasks:

### 1. **AnnData to Seurat Conversion** (`skills/anndata-seurat-conversion/`)
Convert single-cell RNA-seq data from AnnData h5ad format to Seurat RDS format:
- Full preservation of expression data (counts and normalized)
- Cell metadata and dimensional reductions (PCA, UMAP, t-SNE)
- Automatic counts source selection (raw.X, layers, or X fallback)
- Comprehensive logging and validation

### 2. **ARCHS4** (`skills/archs4/`)
Query and analyze gene expression data from the ARCHS4 database:
- Expression atlas for 72+ human/mouse tissues
- Multi-gene expression analysis with visualization
- Gene correlation and co-expression analysis
- Differential expression in specific contexts
- Disease and treatment-specific sample queries
- Smart tissue filtering with typo tolerance

### 3. **CellChat** (`skills/cellchat/`)
CellChat cell-cell communication analysis for single-cell RNA-seq data:
- End-to-end workflow for inferring ligand-receptor interactions
- Multi-condition comparative analysis
- CellChatDB v2 with human/mouse support
- Comprehensive visualizations (heatmaps, circle plots, bubble plots)
- Differential signaling and pathway analysis

### 4. **CELLxGENE Census** (`skills/cellxgene/`)
Comprehensive toolkit for working with single-cell RNA-seq data from CELLxGENE Census:
- **Query subskill**: Download and filter expression data with flexible metadata filtering
- **Specificity subskill**: Extract cell type marker genes and analyze expression specificity
- Gene coexpression analysis with visualization

### 5. **Cortellis** (`skills/cortellis/`)
Unified toolkit for Cortellis Drug Discovery Intelligence and OFF-X:
- **CI Web**: Build Excel-to-interactive CI dashboards with competition/opportunity scoring
- **API Access**: Query targets/drugs via API for structured analysis
- **Target-Drug Web**: Automate web exports across Cortellis categories
- **OFF-X Web**: Export safety evidence and adverse events

### 6. **DESeq2 RNA-seq Analysis** (`skills/deseq2-rna/`)
Interactive bulk RNA-seq differential expression analysis using DESeq2:
- Raw count data processing with sample metadata
- Multi-factor design formula support with collinearity testing
- Pathway enrichment (GSEA with MSigDB, g:Profiler ORA)
- Volcano plots, PCA, and comprehensive visualizations
- Interactive HTML plots with hover tooltips

### 7. **Disease Module Analysis** (`skills/disease-module-analysis/`)
Specialized workflows for disease module analysis using single-cell RNA-seq data:
- Disease module training and transformation
- Gene set enrichment analysis (pathway/ontology)
- Perturbation enrichment (TF/Gene/LINCS)
- Regulator identification
- NCBI metadata extraction

### 8. **Protein Interaction Query** (`skills/interactdb-query/`)
Query protein-protein interaction databases (STRING, IntAct, BioGRID):
- Single-gene neighbor analysis
- Multi-gene shortest path finding
- Comprehensive filtering with entity and edge annotations
- Network visualization preparation

### 9. **KG Association** (`skills/kg-association/`)
Analyze biomedical entity associations using multiple knowledge graph methods:
- **BioBridge**: Neural KG link prediction with mean embeddings for combos
- **ULTRA/UltraQuery**: Foundation model predictions with intersection queries
- **PrimeKG**: Shortest path analysis using graph traversal
- Structured markdown reports with biological interpretations
- Cross-method comparison and confidence assessment

### 10. **NCBI Paper Query** (`skills/ncbi-paper-query/`)
Query PubMed for publications and extract omics dataset accessions:
- Disease/tissue/organism-specific literature search
- Omics accession extraction (GEO, SRA, ArrayExpress, PRIDE, etc.)
- Multiple extraction modes (abstract-only, web scrape, PDF download)
- Journal impact factor filtering and year cutoffs
- Institutional authentication support for subscription publishers

### 11. **Omicsoft DEG Analysis** (`skills/omicsoft-deg-analysis/`)
Analyze pre-computed differential gene expression statistics:
- Query bulk DEG data from h5ad files
- Filter by disease, study, or custom criteria
- Signature-based expression pattern analysis
- GSEA enrichment with leading edge annotation

### 12. **Pathway Database Query** (`skills/pathwaydb-query/`)
Query pathway databases (KEGG, Reactome, MSigDB):
- Find gene-associated pathways and terms
- Multi-gene comparative analysis with UpSet plots
- Pathway similarity network construction using Jaccard indices
- Gene-specific subnetwork extraction with centrality analysis

### 13. **Pseudobulk DEG with Speckle** (`skills/pseudobulkdge-speckle/`)
Single-cell RNA-seq pseudobulk differential expression analysis:
- Cell type-specific DEG analysis using scran/edgeR
- Pathway enrichment with fgsea (Reactome, MSigDB)
- Differential cell composition using speckle
- Support for custom gene signatures
- Interactive volcano plots with hover tooltips

### 14. **Target Prioritization Report** (`skills/target-prioritization-report/`)
Create traceable target-prioritization pipelines and interactive HTML reports:
- Data loaders for Clinical/Disease/Safety/Opportunity/Novelty scoring
- Single-file offline interactive HTML reports with Plotly
- Support for both single genes and combinations
- Configurable scoring weights and normalizers

## MCP Servers

MCP servers enable Claude to interact with external services and data sources:

### 1. **BioBridge MCP** (`mcp/biobridge/`)
Multimodal knowledge graph MCP server for predicting biomedical associations:
- `predict_associations` tool for inferring biomedical links
- Neural retrieval over pre-trained embeddings
- Support for gene-disease, drug-phenotype, and other entity associations
- Flexible entity/relation inputs with cosine similarity scoring

### 2. **PrimeKG MCP** (`mcp/primekg/`)
Comprehensive query interface for the PrimeKG biomedical knowledge graph:
- 129,000+ entities and 8M+ relationships
- Entity search and connection retrieval
- Multi-hop neighborhood queries
- Graph exploration and traversal

*Source: [ai-sci-mcp-services](https://github.com/oneTakeda/ai-sci-mcp-services) by Roger Tu*

### 3. **ULTRA Inference MCP** (`mcp/ultra-inference/`)
On-demand inference server for the ULTRA foundation model on PrimeKG:
- Zero-shot link prediction for any (head_entity, relation) query
- No pre-computed results required
- GPU support with CPU fallback
- Full data transformation pipeline with novelty detection

*Source: [ai-sci-mcp-services](https://github.com/oneTakeda/ai-sci-mcp-services) by Roger Tu*

### 4. **UltraQuery Inference MCP** (`mcp/ultraquery-inference/`)
Complex logical query answering on PrimeKG using UltraQuery model:
- Multi-hop reasoning with intersection, union, and negation operations
- Supports queries like "What diseases are associated with proteins that interact with both GREM1 and IL11?"
- Graph schema validation for query validation and result filtering
- 10 entity types and 18 relation types

*Source: [ai-sci-mcp-services](https://github.com/oneTakeda/ai-sci-mcp-services) by Roger Tu*

## Getting Started

Each skill and MCP server contains its own documentation:
- Skills: See `SKILL.md` in each skill directory
- MCP servers: See `README.md` in each MCP directory

## Contact

**Maintainer:** Xinghao Zhang
**Email:** xinghao.zhang@takeda.com

## Contributing

This repository is maintained for internal use. For questions, issues, or contributions, please contact the maintainer.

## License

Internal use only. Please contact the maintainer for licensing information.
