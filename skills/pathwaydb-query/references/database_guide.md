# Pathway Database Reference Guide

## Overview

This guide provides detailed information about the three pathway databases supported by this skill: KEGG, Reactome, and MSigDB. Reference this document when selecting database-specific parameters or understanding query options.

## KEGG (Kyoto Encyclopedia of Genes and Genomes)

### About KEGG
KEGG is a comprehensive database of biological pathways, diseases, drugs, and genomes. It integrates genomic, chemical, and systemic functional information.

### Organism Codes
KEGG uses 3-letter organism codes:

**Common Organisms:**
- `hsa` - Homo sapiens (human)
- `mmu` - Mus musculus (mouse)
- `rno` - Rattus norvegicus (rat)
- `dre` - Danio rerio (zebrafish)
- `dme` - Drosophila melanogaster (fruit fly)
- `cel` - Caenorhabditis elegans (roundworm)
- `sce` - Saccharomyces cerevisiae (yeast)
- `eco` - Escherichia coli K-12
- `ath` - Arabidopsis thaliana (thale cress)

**Other Model Organisms:**
- `bta` - Bos taurus (cow)
- `cfa` - Canis familiaris (dog)
- `gga` - Gallus gallus (chicken)
- `ptr` - Pan troglodytes (chimpanzee)
- `ssc` - Sus scrofa (pig)
- `xla` - Xenopus laevis (African clawed frog)

### KEGG Query Parameters
- `org` (organism): 3-letter organism code (default: "hsa")
- Gene input: Gene symbol (e.g., "TP53", "BRCA1")

### KEGG Output Format
Returns pathway IDs and names:
```
path:hsa04115 - p53 signaling pathway - Homo sapiens (human)
path:hsa04210 - Apoptosis - Homo sapiens (human)
```

### KEGG Limitations
- Gene symbol matching may vary by organism
- Some organisms have limited pathway annotations
- Pathway IDs are organism-specific (e.g., hsa:04115 vs mmu:04115)

---

## Reactome

### About Reactome
Reactome is a manually curated, peer-reviewed database of human biological pathways and reactions. It provides detailed molecular-level information about cellular processes.

### Species Taxonomy IDs
Reactome uses NCBI Taxonomy IDs:

**Common Species:**
- `9606` - Homo sapiens (human)
- `10090` - Mus musculus (mouse)
- `10116` - Rattus norvegicus (rat)
- `7955` - Danio rerio (zebrafish)
- `6239` - Caenorhabditis elegans (roundworm)
- `7227` - Drosophila melanogaster (fruit fly)
- `9615` - Canis familiaris (dog)
- `9913` - Bos taurus (cow)
- `9031` - Gallus gallus (chicken)

**Alternative Names:**
The API also accepts common names:
- "human" → 9606
- "mouse" → 10090
- "rat" → 10116
- "zebrafish" → 7955

### Reactome Query Parameters
- `species`: Taxonomy ID or common name (default: 9606)
- `resource`: Identifier type - "UniProt", "NCBI", "ENSEMBL" (default: "UniProt")
- `map_to`: Query type - "pathways" or "reactions" (default: "pathways")

### Reactome Output Format
Returns pathway/reaction IDs and names:
```
R-HSA-6804116 - TP53 Regulates Transcription of Genes Involved in G1 Cell Cycle Arrest (Homo sapiens)
R-HSA-6803207 - TP53 Regulates Transcription of Caspase Activators and Caspases (Homo sapiens)
```

### Reactome Features
- Hierarchical pathway organization
- Detailed reaction-level information available
- Cross-references to other databases (UniProt, NCBI, etc.)
- Focus on human pathways with inferred orthologous pathways for other species

---

## MSigDB (Molecular Signatures Database)

### About MSigDB
MSigDB is a collection of annotated gene sets from various sources including published studies, pathway databases, and computational analyses. It's widely used for gene set enrichment analysis (GSEA).

### MSigDB Collections

#### H - Hallmark Gene Sets (50 gene sets)
**Description:** Well-defined biological states and processes with coherent expression
**Best for:** High-level biological process overview
**Examples:** HALLMARK_APOPTOSIS, HALLMARK_P53_PATHWAY, HALLMARK_DNA_REPAIR

#### C1 - Positional Gene Sets (299 gene sets)
**Description:** Gene sets corresponding to human chromosome cytogenetic bands
**Best for:** Chromosomal location analysis
**Examples:** CHR1P36, CHR3Q21

#### C2 - Curated Gene Sets (7,561 gene sets)
**Description:** Pathways from online databases, publications, and knowledge of domain experts
**Subcollections:**
- **CP (Canonical Pathways)**: From pathway databases
  - **CP:KEGG**: KEGG pathways
  - **CP:REACTOME**: Reactome pathways
  - **CP:BIOCARTA**: BioCarta pathways
  - **CP:PID**: Pathway Interaction Database
  - **CP:WIKIPATHWAYS**: WikiPathways
- **CGP (Chemical and Genetic Perturbations)**: Gene expression signatures from perturbational experiments

**Best for:** Comprehensive pathway analysis
**Examples:** KEGG_P53_SIGNALING_PATHWAY, REACTOME_DNA_REPAIR

#### C3 - Regulatory Target Gene Sets (3,735 gene sets)
**Description:** Gene sets based on regulatory motifs in DNA
**Subcollections:**
- **MIR (microRNA targets)**: Predicted microRNA targets
- **TFT (transcription factor targets)**: Predicted transcription factor binding sites

**Best for:** Regulatory network analysis
**Examples:** MIR-124, FOXO3_TARGET_GENES

#### C4 - Computational Gene Sets (858 gene sets)
**Description:** Gene sets defined by mining large collections of cancer-oriented microarray data
**Best for:** Cancer-specific computational signatures
**Examples:** MODULE_123, CLUSTER_456

#### C5 - Ontology Gene Sets (16,228 gene sets)
**Description:** Gene sets from Gene Ontology (GO) annotations
**Subcollections:**
- **BP (Biological Process)**
- **CC (Cellular Component)**
- **MF (Molecular Function)**
- **HPO (Human Phenotype Ontology)**

**Best for:** Functional annotation, GO term analysis
**Examples:** GO_DNA_REPAIR, GO_APOPTOTIC_PROCESS

#### C6 - Oncogenic Signature Gene Sets (189 gene sets)
**Description:** Gene sets defined directly from microarray data from cancer gene perturbations
**Best for:** Cancer pathway analysis, oncogene/tumor suppressor studies
**Examples:** RAS_SIGNALING_UP, P53_DN.V1_UP

#### C7 - Immunologic Signature Gene Sets (5,219 gene sets)
**Description:** Gene sets from immunological studies and cell type signatures
**Best for:** Immune cell analysis, immunology research
**Examples:** GSE15750_NAIVE_VS_MEMORY_CD8_TCELL_UP

#### C8 - Cell Type Signature Gene Sets (830 gene sets)
**Description:** Gene sets that contain cluster markers for cell types identified in single-cell sequencing studies
**Best for:** Cell type identification, single-cell analysis
**Examples:** EMBRYONIC_STEM_CELLS, CD8_T_CELLS

### MSigDB Query Parameters
- `collection`: Collection code (H, C1, C2, etc.) or list of collections
- `dbver`: MSigDB version (default: "2025.1.Hs")
- Use `--all` flag to query all collections simultaneously

### MSigDB Output Format
Returns gene set names and sizes:
```
HALLMARK_P53_PATHWAY (n=200 genes)
HALLMARK_DNA_REPAIR (n=150 genes)
```

### MSigDB Best Practices

**For Quick Overview:**
- Use collection **H** (Hallmark) - concise and well-defined

**For Comprehensive Analysis:**
- Use collections **H + C2 + C6** for pathway and oncogenic signatures
- Add **C5:BP** for detailed GO biological processes

**For Specific Domains:**
- **Cancer research**: H, C2, C6
- **Immunology**: H, C7
- **Single-cell studies**: H, C8
- **Regulatory networks**: C3
- **Functional annotation**: C5

**Performance Considerations:**
- Querying all collections can be slow (16,000+ gene sets)
- Start with H or C2 for faster results
- Use specific collections based on research question

### MSigDB Versions
The default version is `2025.1.Hs` (human). Other versions:
- `2025.1.Mm` - Mouse
- Older versions available but not recommended

---

## Database Comparison

| Feature | KEGG | Reactome | MSigDB |
|---------|------|----------|--------|
| **Focus** | Pathways across organisms | Human pathways + inference | Gene set collections |
| **Pathway Count** | ~500 pathways | ~2,500 pathways | 50-16,000 (varies by collection) |
| **Organism Support** | 7,000+ organisms | 20+ organisms | Human, Mouse |
| **Curation** | Manual | Manual | Mixed (manual + computational) |
| **Hierarchy** | Flat | Hierarchical | Organized by collection |
| **Updates** | Quarterly | Quarterly | ~Annual |
| **Best For** | Cross-species comparison | Detailed human pathways | GSEA, comprehensive signatures |
| **Limitations** | Limited human detail | Mainly human | Version-dependent, large |

---

## Choosing the Right Database

### Use KEGG when:
- Comparing pathways across different organisms
- Working with non-human model organisms
- Need standardized pathway IDs across species
- Interested in metabolic pathways

### Use Reactome when:
- Deep dive into human molecular mechanisms
- Need detailed reaction-level information
- Working with protein-protein interactions
- Analyzing signaling cascades

### Use MSigDB when:
- Performing gene set enrichment analysis (GSEA)
- Need comprehensive gene signatures
- Working with cancer or immunology data
- Comparing to published experimental signatures

### Use All Three when:
- Comprehensive pathway coverage is needed
- Cross-database validation is important
- Comparing different annotation systems
- Publication-quality analysis requiring multiple sources

---

## Tips for Effective Queries

1. **Start Broad**: Begin with unified query across all databases
2. **Refine by Collection**: Use MSigDB collection filtering for specific domains
3. **Check Organism Codes**: Verify correct organism/species codes before querying
4. **Consider Version**: MSigDB versions affect results; document which version you use
5. **Multi-gene Analysis**: Use UpSet plots to visualize pathway overlap across genes
6. **Export for Analysis**: Export to CSV/Excel for downstream statistical analysis

---

## Common Issues and Solutions

### Issue: No KEGG results for gene
**Solution:**
- Verify gene symbol spelling (case-sensitive)
- Try alternative gene names/symbols
- Check if organism code is correct
- Some genes may not be annotated in KEGG

### Issue: Limited Reactome results
**Solution:**
- Reactome focuses on human; other species have fewer pathways
- Try different identifier types (UniProt, NCBI, ENSEMBL)
- Check if gene is human or has human orthologs

### Issue: MSigDB query is slow
**Solution:**
- Limit to specific collections (H, C2, C6 instead of all)
- First download/cache of GMT files is slow; subsequent queries are fast
- Collections C5 and C7 are very large

### Issue: Gene name not found
**Solution:**
- Check gene symbol capitalization
- Try alternative gene names (e.g., "CDKN1A" vs "P21")
- Use UniProt or NCBI IDs for Reactome
- Some databases use different naming conventions

---

## API Rate Limits and Best Practices

### KEGG
- No strict rate limit but avoid excessive parallel requests
- Cache enabled for pathway names (repeated queries are fast)
- Retry logic handles temporary failures

### Reactome
- Generally fast and reliable
- No known rate limits for reasonable use
- Retry logic included

### MSigDB
- Downloaded via gseapy library (GMT files)
- First query per collection is slow (downloads ~1-10MB)
- Subsequent queries use cache (very fast)
- No API rate limits (local processing after download)

---

## Further Reading

**KEGG:**
- Website: https://www.kegg.jp/
- REST API: https://rest.kegg.jp/
- Documentation: https://www.kegg.jp/kegg/rest/keggapi.html

**Reactome:**
- Website: https://reactome.org/
- Content Service: https://reactome.org/ContentService/
- Documentation: https://reactome.org/dev/content-service

**MSigDB:**
- Website: https://www.gsea-msigdb.org/gsea/msigdb/
- Collections: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
- gseapy library: https://gseapy.readthedocs.io/
