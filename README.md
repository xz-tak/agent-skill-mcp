# Agent Skills & MCP Servers

A collection of self-developed Claude Code skills and Model Context Protocol (MCP) servers for biomedical data analysis and computational biology workflows.

## Installation

### Via Claude Code Plugin System

```
/plugin marketplace add tak-xz/agent-skill-mcp
/plugin install biomedical-skills@agent-skill-mcp
```

### Manual Installation

```bash
# Clone to plugins directory
git clone https://github.com/xz-tak/agent-skill-mcp.git ~/.claude/plugins/marketplaces/agent-skill-mcp

# Or clone to project-local plugins
git clone https://github.com/xz-tak/agent-skill-mcp.git ./.claude/plugins/agent-skill-mcp
```

### Standalone Skills (No Plugin System)

```bash
# Copy specific skills to user's skills directory
git clone https://github.com/xz-tak/agent-skill-mcp.git /tmp/agent-skill-mcp
cp -r /tmp/agent-skill-mcp/skills/archs4 ~/.claude/skills/
cp -r /tmp/agent-skill-mcp/skills/cellxgene ~/.claude/skills/
```

## Overview

This repository contains custom-built tools that extend Claude's capabilities for biomedical research, enabling seamless integration with biological databases, knowledge graphs, and analysis pipelines.

## Repository Structure

```
agent-skill-mcp/
├── skills/          # Claude Code skills for specialized workflows
└── mcp/             # MCP servers for external service integration
```

## Skills

Claude Code skills provide specialized workflows and domain knowledge for specific tasks.

### Single-Cell RNA-seq

| Skill | Description |
|-------|-------------|
| [CELLxGENE Census](skills/cellxgene/) | Query and download expression data from CELLxGENE Census with flexible metadata filtering, cell type marker extraction, and specificity analysis |
| [CellChat](skills/cellchat/) | Cell-cell communication analysis: ligand-receptor interaction inference, multi-condition comparison, CellChatDB v2 |
| [CyteType](skills/cytetype/) | Automated cell type annotation using multi-agent AI architecture with GPU-aware compute, 3-agent consensus, and confidence scoring |
| [PopV](skills/popv/) | PopV SCDS2 consensus cell-type annotation pipeline with multiple classifiers (scVI, scanVI, CellTypist) |
| [scGPT](skills/scgpt/) | Generative pre-trained transformer for single-cell biology: gene expression prediction, in silico perturbation, cell embedding |
| [scimilarity](skills/scimilarity/) | Single-cell foundation model for cell embedding, annotation, search across reference atlases, and gene interpretation |
| [Geneformer](skills/geneformer/) | Foundational transformer pretrained on ~104M single-cell transcriptomes for context-aware predictions in network biology |
| [Pseudobulk DEG + Speckle](skills/pseudobulkdge-speckle/) | Pseudobulk differential expression (scran/edgeR), pathway enrichment (fgsea), and differential cell composition (speckle) |
| [Disease Module Analysis](skills/disease-module-analysis/) | Disease module training/transformation, gene set enrichment, perturbation enrichment, and regulator identification |
| [ReCoN Multinetwork](skills/recon-multinetwork/) | Multicellular coordination network analysis integrating GRNs with cell-cell communication from scRNA-seq and optional scATAC-seq |

### Bulk RNA-seq

| Skill | Description |
|-------|-------------|
| [DESeq2 RNA-seq](skills/deseq2-rna/) | Interactive bulk RNA-seq differential expression analysis with multi-factor designs, GSEA pathway enrichment, and interactive visualizations |
| [AnnData-Seurat Conversion](skills/anndata-seurat-conversion/) | Bidirectional h5ad/RDS conversion preserving expression data, metadata, reductions, and layers |
| [Omicsoft Analysis](skills/omicsoft-analysis/) | Analyze pre-computed differential expression, expression data, and enrichment from Omicsoft h5ad files or TileDB-SOMA on S3 |

### Pathway & Network Analysis

| Skill | Description |
|-------|-------------|
| [Pathway Database Query](skills/pathwaydb-query/) | Query KEGG, Reactome, MSigDB for gene-pathway associations; build Jaccard similarity networks; pathway module clustering (127 modules, 26,881 pathways) |
| [Protein Interaction Query](skills/interactdb-query/) | Query STRING, IntAct, BioGRID for PPI neighbors and multi-gene shortest paths |
| [KG Association](skills/kg-association/) | Multi-method biomedical entity association analysis using BioBridge, ULTRA/UltraQuery, and PrimeKG |

### Gene Expression & Databases

| Skill | Description |
|-------|-------------|
| [ARCHS4](skills/archs4/) | Gene expression atlas for 72+ human/mouse tissues with correlation, differential expression, and tissue-specific analysis |
| [NCBI Paper Query](skills/ncbi-paper-query/) | PubMed literature search with omics accession extraction (GEO, SRA, ArrayExpress) and impact factor filtering |

### Genetics & Proteomics

| Skill | Description |
|-------|-------------|
| [Genetics GSP](skills/genetics-gsp/) | Aggregated human genetics data from GSP: GWAS, Mendelian genetics (OMIM), gene burden for target safety and disease associations |
| [Genetics UKB-PPP](skills/genetics-ukbppp/) | UKB-PPP plasma proteomics disease associations with logistic regression, ARD, and CoxPH models |

### Drug Discovery & Target Assessment

| Skill | Description |
|-------|-------------|
| [Cortellis](skills/cortellis/) | Cortellis Drug Discovery Intelligence + OFF-X: CI dashboards, API queries, web export automation |
| [DrugBank](skills/drugbank/) | Query DrugBank for drug-target interactions, indications, pharmacology, and adverse effects |
| [GOSTAR](skills/gostar/) | Query GOSTAR medicinal chemistry database for SAR data, bioassay/bioactivity, and pharmacological data |
| [DrugnomeAI](skills/drugnomeai/) | Gene druggability predictions using PU learning with consensus tiers, composite scores, and modality prediction |
| [Target Prioritization Report](skills/target-prioritization-report/) | Traceable target-prioritization pipelines with interactive HTML reports and configurable scoring |
| [AgenticBoost Prompt](skills/agenticboost-prompt/) | Generate target evaluation documents from template, auto-populating Cortellis, pathway, and PPI data |

### Structure Prediction

| Skill | Description |
|-------|-------------|
| [Boltz](skills/boltz/) | Biomolecular structure prediction using Boltz 2.2.1: proteins, complexes, protein-ligand, protein-DNA/RNA systems |

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
