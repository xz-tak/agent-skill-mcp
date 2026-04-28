#!/usr/bin/env python3
"""
DESeq2 RNA-seq Analysis Report Generator
Generates comprehensive markdown report from DESeq2 analysis outputs
Parameterized for use with any DESeq2 analysis

CODEX_XHIGH PROTOCOL:
This script generates draft reports with interpretation placeholders that include
CODEX_XHIGH directives for Codex to provide deep biological analysis. Each placeholder
requires Codex to:
1. Assume the role of a professional biologist (molecular, disease, computational biology)
2. Use extended thinking (codex xhigh) for mechanistic insights
3. Integrate findings with established biological knowledge
4. Use WebFetch for relevant data retrieval when needed
5. Generate specific, testable hypotheses
"""

import argparse
import os
import json
import pandas as pd
from datetime import datetime
from pathlib import Path


# ==============================================================================
# Argument Parsing
# ==============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate markdown report from DESeq2 analysis outputs'
    )
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing analysis outputs (figures/, deg/)')
    parser.add_argument('--output', required=True,
                        help='Output markdown report path')
    parser.add_argument('--prefix', default='RNAseq',
                        help='Prefix used in analysis output files')
    parser.add_argument('--title', default='RNA-seq DESeq2 Analysis Report',
                        help='Report title')
    parser.add_argument('--description', default='',
                        help='Study description for executive summary')
    parser.add_argument('--padj', type=float, default=0.05,
                        help='Adjusted p-value threshold used')
    parser.add_argument('--lfc', type=float, default=1.0,
                        help='Log2 fold change threshold used')
    parser.add_argument('--config', default=None,
                        help='JSON config file (overrides other arguments)')

    return parser.parse_args()


def load_config(args):
    """Load parameters from config file if provided."""
    params = vars(args)

    if args.config and os.path.exists(args.config):
        with open(args.config, 'r') as f:
            config = json.load(f)
        params.update(config)

    return params


# ==============================================================================
# Data Loading
# ==============================================================================

def load_data(params):
    """Load all analysis data files."""
    input_dir = Path(params['input_dir'])
    deg_dir = input_dir / 'deg'
    prefix = params['prefix']

    data = {}

    # Load DE counts summary
    de_counts_pattern = f"{prefix}_DE_lfc*_count.txt"
    de_counts_files = list(deg_dir.glob(de_counts_pattern))
    if de_counts_files:
        de_counts_file = de_counts_files[0]
        de_counts_df = pd.read_csv(de_counts_file, sep='\t')
        if 'comparison' in de_counts_df.columns:
            de_counts_df = de_counts_df.set_index('comparison')
        data['de_counts'] = de_counts_df
        print(f"  Loaded DE counts: {len(data['de_counts'])} comparisons")

    # Load combined summary stats
    summstats_file = deg_dir / f'{prefix}_summstats_all.txt'
    if summstats_file.exists():
        data['summstats'] = pd.read_csv(summstats_file, sep='\t')
        print(f"  Loaded summary stats: {len(data['summstats'])} entries")

    # Load GSEA results
    gsea_file = deg_dir / f'{prefix}_gseaGObp_all.txt'
    if gsea_file.exists() and gsea_file.stat().st_size > 0:
        try:
            data['gsea'] = pd.read_csv(gsea_file, sep='\t')
            print(f"  Loaded GSEA results: {len(data['gsea'])} terms")
        except Exception:
            data['gsea'] = None
    else:
        data['gsea'] = None

    # Load g:Profiler results
    gp_file = deg_dir / f'{prefix}_defisher_all.txt'
    if gp_file.exists() and gp_file.stat().st_size > 0:
        try:
            data['gprofiler'] = pd.read_csv(gp_file, sep='\t')
            print(f"  Loaded g:Profiler results: {len(data['gprofiler'])} terms")
        except Exception:
            data['gprofiler'] = None
    else:
        data['gprofiler'] = None

    return data


# ==============================================================================
# Helper Functions
# ==============================================================================

def get_de_count(de_counts, comparison):
    """Get DE gene counts for a comparison."""
    if de_counts is None or comparison not in de_counts.index:
        return 0, 0
    row = de_counts.loc[comparison]
    return int(row.get('up', 0)), int(row.get('dn', 0))


def format_comparison_name(name):
    """Format comparison name for display."""
    return name.replace('_', ' ').replace('  ', ' vs ')


def extract_top_pathways(data, comparison_name, source='GSEA', n_top=5):
    """Extract top pathways for a comparison."""
    if source == 'GSEA' and data.get('gsea') is not None:
        df = data['gsea']
        if 'comparison' in df.columns:
            df = df[df['comparison'] == comparison_name]
        if len(df) > 0 and 'NES' in df.columns:
            up = df[df['NES'] > 0].sort_values('NES', ascending=False)
            down = df[df['NES'] < 0].sort_values('NES', ascending=True)
            up_terms = up['Description'].head(n_top).tolist() if 'Description' in up.columns else []
            down_terms = down['Description'].head(n_top).tolist() if 'Description' in down.columns else []
            return {'up': up_terms, 'down': down_terms}

    elif source == 'gProfiler' and data.get('gprofiler') is not None:
        df = data['gprofiler']
        if 'comparison' in df.columns:
            df = df[df['comparison'] == comparison_name]
        if len(df) > 0:
            up_col = 'DE' if 'DE' in df.columns else None
            if up_col:
                up = df[df[up_col] == 'up'].sort_values('p_value' if 'p_value' in df.columns else df.columns[0])
                down = df[df[up_col] == 'dn'].sort_values('p_value' if 'p_value' in df.columns else df.columns[0])
                term_col = 'term_name' if 'term_name' in df.columns else 'Description' if 'Description' in df.columns else df.columns[0]
                return {'up': up[term_col].head(n_top).tolist(), 'down': down[term_col].head(n_top).tolist()}

    return {'up': [], 'down': []}


def format_pathway_list(pathways, max_show=5):
    """Format pathway list for display."""
    if not pathways:
        return "- None identified"
    return '\n'.join([f"- {p}" for p in pathways[:max_show]])


def get_top_genes(data, comparison_name, direction='up', n_top=10):
    """Extract top genes for a comparison."""
    if data.get('summstats') is None:
        return []

    df = data['summstats']
    if 'comparison' not in df.columns:
        return []

    comp_df = df[df['comparison'] == comparison_name].copy()
    if len(comp_df) == 0:
        return []

    # Filter for significant genes
    if 'padj' in comp_df.columns and 'log2FoldChange' in comp_df.columns:
        if direction == 'up':
            sig = comp_df[(comp_df['padj'] < 0.05) & (comp_df['log2FoldChange'] > 1)]
            sig = sig.sort_values('log2FoldChange', ascending=False)
        else:
            sig = comp_df[(comp_df['padj'] < 0.05) & (comp_df['log2FoldChange'] < -1)]
            sig = sig.sort_values('log2FoldChange', ascending=True)

        gene_col = 'symbol' if 'symbol' in sig.columns else 'gene' if 'gene' in sig.columns else None
        if gene_col:
            return sig[gene_col].head(n_top).tolist()

    return []


def generate_interpretation_placeholder(section_id, section_title, data_context):
    """Generate placeholder section for Codex to fill in biological interpretation.

    Includes CODEX_XHIGH directive for deep biological analysis with context-specific focus.
    """
    return f"""
<!-- CODEX_INTERPRETATION_START: {section_id} -->

#### {section_title}

**[INTERPRETATION REQUIRED]**

---
**⚠️ CODEX_XHIGH REQUIRED ⚠️**

**Codex's Role:** Professional Biologist with Broad Expertise
(Molecular Biology, Disease Biology across all therapeutic areas, Computational Biology)

**CONTEXT-SPECIFIC FOCUS:** Interpret findings in the context of THIS study's:
- Specific target/pathway under investigation
- Specific disease/indication being studied
- Cell type/model system used
- Treatments and their mechanisms of action

**CODEX_XHIGH Protocol:**
1. **Context Identification**: First identify the study-specific biological context
2. **Extended Analysis**: Deeply analyze the data in that specific context
3. **Mechanistic Insight**: Explain mechanisms SPECIFIC to this target/disease
4. **Adaptive Knowledge Retrieval**: Use WebFetch to retrieve from most relevant sources
   (PubMed, UniProt, KEGG, Reactome, OMIM, DrugBank, OpenTargets, or any other database)
5. **Clinical Relevance**: Identify implications for THIS therapeutic area
6. **Hypothesis Generation**: Propose hypotheses relevant to THIS study's context

---

*Codex should provide detailed biological interpretation based on the following data:*

{data_context}

<!-- CODEX_INTERPRETATION_END: {section_id} -->
"""


def format_pathway_table(data, comparison_name, direction='up', n_top=5):
    """Format pathway enrichment as a table with data for interpretation."""
    pathways = extract_top_pathways(data, comparison_name, 'GSEA', n_top)
    pathway_list = pathways.get(direction, [])

    if not pathway_list:
        return "No significantly enriched pathways identified."

    table = "| Pathway | Biological Significance |\n|---------|-------------------------|\n"
    for p in pathway_list[:n_top]:
        table += f"| {p} | [Interpretation needed] |\n"

    return table


# ==============================================================================
# Report Section Generators
# ==============================================================================

def generate_header(params):
    """Generate report header with CODEX_XHIGH notice."""
    return f"""# {params['title']}

**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

**Analysis Directory:** `{params['input_dir']}`

---

<!-- CODEX_XHIGH_INSTRUCTIONS_START
================================================================================
⚠️ CODEX_XHIGH REQUIRED FOR BIOLOGICAL INTERPRETATIONS ⚠️

Codex's Role: Professional Biologist with Broad Expertise
- Molecular Biology: Gene regulation, signaling pathways, protein function
- Disease Biology: Broad expertise across therapeutic areas (oncology, immunology,
  fibrosis, metabolic, neurological, cardiovascular, infectious disease, rare diseases)
- Computational Biology: Transcriptomics, pathway analysis, systems biology

CONTEXT-SPECIFIC FOCUS:
While having broad expertise, Codex MUST focus interpretations on:
- The SPECIFIC target/pathway being studied (infer from comparisons, top genes)
- The SPECIFIC disease/indication context relevant to this experiment
- The cell type/model system used and its biological relevance
- The mechanism of action of treatments tested in this study

CODEX_XHIGH Protocol for ALL interpretation sections:
1. First IDENTIFY the study context (target, disease, cell type, treatments)
2. Use extended thinking to deeply analyze data in that specific context
3. Provide MECHANISTIC explanations specific to this target/disease
4. Use WebFetch to ADAPTIVELY retrieve from most relevant sources (PubMed, UniProt,
   KEGG, Reactome, OMIM, DrugBank, OpenTargets, or ANY other relevant database)
5. Connect findings to the SPECIFIC disease biology under investigation
6. Generate testable hypotheses relevant to THIS therapeutic context

Quality Standards:
- Every interpretation must be SPECIFIC to this study's biological context
- Cite mechanisms relevant to the target/disease being studied
- Avoid generic interpretations - always tie back to study-specific biology
- Integrate multiple data types (DEGs + pathways + known biology for THIS context)
================================================================================
CODEX_XHIGH_INSTRUCTIONS_END -->

---
"""


def generate_toc():
    """Generate table of contents."""
    return """## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Input Data Description](#2-input-data-description)
3. [Parameters and Analysis Setup](#3-parameters-and-analysis-setup)
4. [Results](#4-results)
   - 4.1 [Quality Control and Sample Relationships](#41-quality-control-and-sample-relationships)
   - 4.2 [Differential Expression Overview](#42-differential-expression-overview)
   - 4.3 [Per-Comparison Results](#43-per-comparison-results)
   - 4.4 [Pathway Enrichment Summary](#44-pathway-enrichment-summary)
5. [Discussion](#5-discussion)
6. [Plot and Table Interpretation Guide](#6-plot-and-table-interpretation-guide)

---
"""


def generate_executive_summary(data, params):
    """Generate executive summary section."""
    de_counts = data.get('de_counts')

    if de_counts is not None:
        total_comparisons = len(de_counts)
        de_counts_copy = de_counts.copy()
        de_counts_copy['total'] = de_counts_copy['up'] + de_counts_copy['dn']
        top_comparison = de_counts_copy['total'].idxmax()
        max_degs = int(de_counts_copy['total'].max())
        comps_with_degs = int((de_counts_copy['total'] > 0).sum())
        total_up = int(de_counts_copy['up'].sum())
        total_dn = int(de_counts_copy['dn'].sum())
    else:
        total_comparisons = 0
        comps_with_degs = 0
        top_comparison = "N/A"
        max_degs = 0
        total_up = 0
        total_dn = 0

    description = params.get('description', '')
    if not description:
        description = "This report presents the results of differential gene expression (DEG) analysis from an RNA-seq experiment."

    return f"""## 1. Executive Summary

### Study Overview

{description}

**Analysis Method:** DESeq2 with VST normalization

### Key Statistics

| Metric | Value |
|--------|-------|
| Total Comparisons | {total_comparisons} |
| Comparisons with DEGs | {comps_with_degs} |
| Most DEGs | {format_comparison_name(top_comparison)} ({max_degs} genes) |
| Total Upregulated | {total_up} |
| Total Downregulated | {total_dn} |
| Significance Threshold | padj < {params['padj']}, |log2FC| > {params['lfc']} |

---
"""


def generate_input_data(params):
    """Generate input data description section."""
    return """## 2. Input Data Description

### Gene Expression Data

- **Format:** Raw counts (non-negative integers)
- **Rows:** Genes (ENSEMBL IDs or gene symbols)
- **Columns:** Samples

### Sample Metadata

- Sample identifiers matching counts columns
- Treatment/condition information
- Optional covariates (batch, etc.)

### Quality Control

- Low-count genes filtered (< 10 total counts across samples)
- Variance stabilizing transformation (VST) applied for visualization

---
"""


def generate_parameters(params):
    """Generate parameters and analysis setup section."""
    return f"""## 3. Parameters and Analysis Setup

### DESeq2 Configuration

```r
# Low count filter
dds <- dds[rowSums(counts(dds)) >= 10, ]

# DESeq2 analysis
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
```

### Significance Thresholds

| Parameter | Value | Description |
|-----------|-------|-------------|
| padj | < {params['padj']} | Benjamini-Hochberg adjusted p-value |
| log2FoldChange | > {params['lfc']} or < -{params['lfc']} | Minimum fold change |

### Enrichment Analysis Settings

**GSEA (Gene Set Enrichment Analysis):**
- Ontology: GO Biological Process (BP)
- Ranking metric: -log10(padj) x log2FoldChange
- P-value cutoff: 0.05
- Adjustment method: Benjamini-Hochberg

**g:Profiler:**
- Databases: GO:BP, Reactome, KEGG
- FDR correction method: g:SCS
- User threshold: 0.05

---
"""


def generate_results_qc(params):
    """Generate QC and sample relationships section."""
    prefix = params['prefix']
    return f"""## 4. Results

### 4.1 Quality Control and Sample Relationships

#### PCA Analysis

Principal Component Analysis (PCA) was performed on variance-stabilized transformed (VST) expression data to assess sample relationships and identify major sources of variation.

![PCA Basic](figures/{prefix}_PCA_basic_*.png)

**Figure 4.1.1:** Basic PCA plot colored by treatment group. Samples cluster primarily by condition, with clear separation between experimental groups.

![PCA Screeplot](figures/{prefix}_PCA_screeplot.png)

**Figure 4.1.2:** Scree plot showing variance explained by each principal component.

#### Key PCA Observations

1. **PC1** typically captures the major treatment effect
2. **PC2** often separates secondary experimental factors
3. Tight clustering within groups indicates low technical variability
4. Separation between groups indicates biological signal

---
"""


def generate_results_de_overview(data, params):
    """Generate DE overview section."""
    de_counts = data.get('de_counts')

    table_rows = []
    if de_counts is not None:
        for idx in de_counts.index:
            row = de_counts.loc[idx]
            up = int(row['up'])
            dn = int(row['dn'])
            total = up + dn
            table_rows.append(f"| {format_comparison_name(idx)} | {up} | {dn} | {total} |")

    table_str = '\n'.join(table_rows)

    return f"""### 4.2 Differential Expression Overview

#### DE Gene Counts Summary

The following table summarizes differentially expressed genes (DEGs) across all comparisons using thresholds of padj < {params['padj']} and |log2FC| > {params['lfc']}.

| Comparison | Upregulated | Downregulated | Total |
|------------|-------------|---------------|-------|
{table_str}

**Reference:** `deg/{params['prefix']}_DE_lfc{params['lfc']}padj{params['padj']}_count.txt`

![DE Summary Barplot](figures/DE_genes_summary_barplot.png)

**Figure 4.2.1:** Bar plot showing the number of upregulated (red) and downregulated (blue) genes per comparison.

#### Key Observations

1. **Largest effects:** Comparisons with the most DEGs indicate strongest biological signals
2. **Directional bias:** Note if certain treatments tend to upregulate or downregulate more genes
3. **Context dependency:** Effects may vary by experimental context

---
"""


def generate_results_per_comparison(data, params):
    """Generate per-comparison results section with interpretation placeholders."""
    de_counts = data.get('de_counts')

    if de_counts is None:
        return """### 4.3 Per-Comparison Results

No comparison results available.

---
"""

    sections = ["### 4.3 Per-Comparison Results\n"]

    for i, comp in enumerate(de_counts.index, 1):
        up, dn = get_de_count(de_counts, comp)
        pathways = extract_top_pathways(data, comp, 'GSEA', 5)

        # Get top genes for context
        top_up_genes = get_top_genes(data, comp, 'up', 10)
        top_dn_genes = get_top_genes(data, comp, 'down', 10)

        # Sanitize comparison name for file paths
        comp_file = comp.replace(' ', '_')
        formatted_name = format_comparison_name(comp)

        # Build data context for interpretation
        data_context = f"""
- **Comparison:** {formatted_name}
- **Upregulated genes:** {up}
- **Downregulated genes:** {dn}
- **Top upregulated genes:** {', '.join(top_up_genes[:10]) if top_up_genes else 'None identified'}
- **Top downregulated genes:** {', '.join(top_dn_genes[:10]) if top_dn_genes else 'None identified'}
- **Top enriched pathways (up):** {', '.join(pathways['up'][:5]) if pathways['up'] else 'None'}
- **Top enriched pathways (down):** {', '.join(pathways['down'][:5]) if pathways['down'] else 'None'}
"""

        sections.append(f"""
#### 4.3.{i} {formatted_name}

##### Differential Expression Summary

| Metric | Value |
|--------|-------|
| Upregulated genes | {up} |
| Downregulated genes | {dn} |
| Total DEGs | {up + dn} |

**Top 10 Upregulated Genes:** {', '.join(top_up_genes) if top_up_genes else 'None identified'}

**Top 10 Downregulated Genes:** {', '.join(top_dn_genes) if top_dn_genes else 'None identified'}

![Volcano {comp}](figures/volcano_{comp_file}.png)

**Figure 4.3.{i}:** Volcano plot showing differential expression. Red: significantly upregulated; Blue: significantly downregulated.

##### Pathway Enrichment

**Enriched Pathways (Upregulated genes):**
{format_pathway_list(pathways['up'])}

**Enriched Pathways (Downregulated genes):**
{format_pathway_list(pathways['down'])}

{generate_interpretation_placeholder(
    f"comparison_{i}_interpretation",
    "Biological Interpretation",
    data_context
)}

##### Reference Files
- DEG results: `deg/{params['prefix']}_{comp_file}.txt`
- GO enrichment: `deg/{params['prefix']}_{comp_file}_GOBP.txt`
- GSEA dotplot: `figures/GSEA_dotplot_{comp_file}.png`
""")

    sections.append("\n---\n")
    return '\n'.join(sections)


def generate_results_pathway_summary(params):
    """Generate pathway enrichment summary section."""
    return f"""### 4.4 Pathway Enrichment Summary

#### Overview of Enrichment Methods

Two complementary approaches were used for pathway analysis:

1. **GSEA (Gene Set Enrichment Analysis):** Rank-based method using all genes
2. **g:Profiler (Fisher's Exact Test):** Over-representation analysis of significant DEGs

#### Pathway Results Files

| Analysis | Combined File |
|----------|---------------|
| GSEA GO:BP | `deg/{params['prefix']}_gseaGObp_all.txt` |
| g:Profiler | `deg/{params['prefix']}_defisher_all.txt` |

Individual comparison results are available in the `deg/` directory.

---
"""


def generate_discussion(data, params):
    """Generate discussion section with interpretation placeholders."""
    de_counts = data.get('de_counts')

    # Build summary context for Claude
    if de_counts is not None:
        total_comparisons = len(de_counts)
        de_counts_copy = de_counts.copy()
        de_counts_copy['total'] = de_counts_copy['up'] + de_counts_copy['dn']
        comparisons_summary = []
        for idx in de_counts.index:
            row = de_counts.loc[idx]
            comparisons_summary.append(f"- {format_comparison_name(idx)}: {int(row['up'])} up, {int(row['dn'])} down")
        comparisons_text = '\n'.join(comparisons_summary)
    else:
        total_comparisons = 0
        comparisons_text = "No comparison data available"

    # Get pathway summary from GSEA data
    pathway_summary = ""
    if data.get('gsea') is not None:
        gsea_df = data['gsea']
        if 'Description' in gsea_df.columns and 'NES' in gsea_df.columns:
            top_up = gsea_df[gsea_df['NES'] > 0].groupby('Description').size().sort_values(ascending=False).head(10)
            top_dn = gsea_df[gsea_df['NES'] < 0].groupby('Description').size().sort_values(ascending=False).head(10)
            pathway_summary = f"""
**Most frequently enriched pathways (upregulated across comparisons):**
{chr(10).join([f'- {p} ({c} comparisons)' for p, c in top_up.items()]) if len(top_up) > 0 else '- None identified'}

**Most frequently enriched pathways (downregulated across comparisons):**
{chr(10).join([f'- {p} ({c} comparisons)' for p, c in top_dn.items()]) if len(top_dn) > 0 else '- None identified'}
"""

    discussion_context = f"""
**Analysis Overview:**
- Total comparisons analyzed: {total_comparisons}
- Significance thresholds: padj < {params['padj']}, |log2FC| > {params['lfc']}

**DEG Summary by Comparison:**
{comparisons_text}

{pathway_summary}
"""

    return f"""## 5. Discussion

{generate_interpretation_placeholder(
    "discussion_major_findings",
    "Summary of Major Findings",
    discussion_context
)}

{generate_interpretation_placeholder(
    "discussion_treatment_effects",
    "Treatment Effects Analysis",
    f'''**🔬 CODEX_XHIGH: Perform deep comparative analysis of treatment effects**

Based on the DEG counts and pathway enrichments above, provide detailed MECHANISTIC analysis of:
- Which treatments produced the strongest transcriptional response and WHY (molecular mechanism)
- What biological processes are most affected by each treatment (connect to known biology)
- How do the test compounds/antibodies modify the response compared to controls
- What does the balance of up vs down regulated genes suggest about cellular state changes
- Use WebFetch to adaptively retrieve from relevant sources (PubMed, databases) for this target'''
)}

{generate_interpretation_placeholder(
    "discussion_pathway_convergence",
    "Pathway Convergence and Divergence",
    f'''**🔬 CODEX_XHIGH: Analyze pathway patterns across all comparisons**

Perform deep analysis of which pathways are:
- Shared across multiple treatments (what does convergent biology reveal about core mechanisms?)
- Unique to specific treatments (what treatment-specific effects are observed?)
- Modified by test compounds vs controls (mechanism of action insights)
- Related to the target/pathway biology under study (connect to established signaling knowledge)
- Identify unexpected pathway connections and propose mechanistic explanations'''
)}

{generate_interpretation_placeholder(
    "discussion_biological_implications",
    "Biological Implications",
    f'''**🔬 CODEX_XHIGH: Synthesize findings into biological and clinical insights**

Provide expert-level synthesis answering:
- What do these results tell us about the target/pathway function in this cell system
- How do the signaling pathways and treatments interact to affect gene expression
- What is the mechanism of action suggested by the treatment effects
- How do these findings relate to the disease/indication under study
- What are the potential clinical or therapeutic applications of these findings'''
)}

{generate_interpretation_placeholder(
    "discussion_hypotheses",
    "Hypotheses for Future Investigation",
    f'''**🔬 CODEX_XHIGH: Generate specific, testable hypotheses**

Based on deep analysis of all findings, generate:
- Specific genes or pathways to validate experimentally (prioritized by evidence strength)
- Predicted functional outcomes of the treatment (what to expect in validation studies)
- Potential biomarkers of pathway activity (candidate markers with rationale)
- Concrete suggestions for follow-up experiments (experimental design recommendations)
- Novel mechanistic hypotheses that emerged from this analysis'''
)}

### Recommended Next Steps

#### From DEG Analysis

1. **Signature Development:**
   - Define activity signatures from treatment-responsive genes
   - Validate signatures in independent datasets

2. **Pathway Deep-Dive:**
   - Detailed analysis of top enriched pathways
   - Investigation of pathway crosstalk

3. **Network Analysis:**
   - Construct gene regulatory networks
   - Identify key transcription factors and hub genes

#### Validation

1. **Functional Assays:**
   - Validate key gene changes at protein level
   - Functional readouts for enriched pathways

2. **Independent Cohorts:**
   - Validate findings in independent samples
   - Cross-platform validation

### Limitations

- In vitro fibroblast model may not fully recapitulate in vivo biology
- Single timepoint analysis; temporal dynamics not captured
- Statistical thresholds may miss biologically relevant subtle changes

---
"""


def generate_interpretation_guide(params):
    """Generate plot and table interpretation guide."""
    return f"""## 6. Plot and Table Interpretation Guide

### PCA Plots

#### Basic PCA Plot
- **Purpose:** Visualize overall sample relationships based on gene expression
- **Interpretation:**
  - Each point represents one sample
  - Color indicates experimental group
  - Samples clustering together have similar expression profiles
  - Distance between points reflects transcriptional similarity

#### Scree Plot
- **Purpose:** Show variance captured by each principal component
- **Interpretation:**
  - Y-axis: Percentage of variance explained
  - X-axis: Principal component number
  - Elbow point indicates optimal number of components

### Volcano Plots

- **Purpose:** Visualize differential expression significance vs magnitude
- **Axes:**
  - X-axis: log2 Fold Change (effect size)
  - Y-axis: -log10(adjusted p-value) (significance)
- **Color coding:**
  - **Red:** Significantly upregulated (padj < {params['padj']}, log2FC > {params['lfc']})
  - **Blue:** Significantly downregulated (padj < {params['padj']}, log2FC < -{params['lfc']})
  - **Grey:** Not significant
- **Reference lines:**
  - Vertical grey lines: log2FC = +/-{params['lfc']}
  - Horizontal grey line: padj = {params['padj']}

### GSEA Dotplots

- **Purpose:** Show enriched GO Biological Process terms
- **Interpretation:**
  - Y-axis: GO term descriptions
  - X-axis: Gene ratio
  - Dot size: Number of genes
  - Dot color: Adjusted p-value
  - Facets: Activated (positive NES) vs Suppressed (negative NES)

### DE Summary Barplot

- **Purpose:** Compare DEG counts across all comparisons
- **Interpretation:**
  - X-axis: Comparison names
  - Y-axis: Number of DEGs
  - Red bars: Upregulated genes
  - Blue bars: Downregulated genes

### Output Tables

#### DEG Results (`deg/{params['prefix']}_*.txt`)

| Column | Description |
|--------|-------------|
| baseMean | Mean normalized count across all samples |
| log2FoldChange | Log2 ratio of expression (treatment/control) |
| lfcSE | Standard error of log2 fold change |
| stat | Wald test statistic |
| pvalue | Raw p-value |
| padj | Benjamini-Hochberg adjusted p-value |
| comparison | Comparison identifier |

#### GSEA Results (`deg/{params['prefix']}_*_GOBP.txt`)

| Column | Description |
|--------|-------------|
| ID | GO term identifier |
| Description | GO term name |
| setSize | Number of genes in the term |
| enrichmentScore | Raw enrichment score |
| NES | Normalized enrichment score |
| pvalue | Nominal p-value |
| p.adjust | Adjusted p-value |

---

## Appendix: File Inventory

### Figures Directory (`figures/`)

| Category | Files |
|----------|-------|
| PCA | `{params['prefix']}_PCA_*.png` |
| Volcano | `volcano_*.png` |
| GSEA | `GSEA_dotplot_*.png` |
| Summary | `DE_genes_summary_barplot.png` |

### DEG Results Directory (`deg/`)

| Category | Files |
|----------|-------|
| Per-comparison DEGs | `{params['prefix']}_*.txt` |
| Per-comparison GO | `{params['prefix']}_*_GOBP.txt` |
| Per-comparison g:Profiler | `{params['prefix']}_*_gp.txt` |
| Combined DEGs | `{params['prefix']}_summstats_all.txt` |
| Combined GSEA | `{params['prefix']}_gseaGObp_all.txt` |
| Combined g:Profiler | `{params['prefix']}_defisher_all.txt` |
| DE Counts Summary | `{params['prefix']}_DE_lfc*_count.txt` |

### R Objects

- `{params['prefix']}_analysis_data.rds` - All analysis objects
- `{params['prefix']}_analysis.RData` - Full R session

---

*Report generated by `generate_report.py`*
"""


# ==============================================================================
# Main Execution
# ==============================================================================

def main():
    print("DESeq2 RNA-seq Analysis Report Generator")
    print("=" * 50)

    # Parse arguments
    args = parse_args()
    params = load_config(args)

    print(f"\nInput directory: {params['input_dir']}")
    print(f"Output file: {params['output']}")
    print(f"Prefix: {params['prefix']}")

    # Load data
    print("\n=== Loading Analysis Data ===")
    data = load_data(params)

    print("\n=== Generating Report ===")

    # Assemble report sections
    report_sections = [
        generate_header(params),
        generate_toc(),
        generate_executive_summary(data, params),
        generate_input_data(params),
        generate_parameters(params),
        generate_results_qc(params),
        generate_results_de_overview(data, params),
        generate_results_per_comparison(data, params),
        generate_results_pathway_summary(params),
        generate_discussion(data, params),
        generate_interpretation_guide(params),
    ]

    # Combine sections
    report = '\n'.join(report_sections)

    # Write report
    output_path = Path(params['output'])
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        f.write(report)

    print(f"\n=== Report Generated ===")
    print(f"Output: {output_path}")
    print(f"Lines: {len(report.splitlines())}")

    print("\nReport sections:")
    print("  1. Executive Summary")
    print("  2. Input Data Description")
    print("  3. Parameters and Analysis Setup")
    print("  4. Results (QC, DE Overview, Per-Comparison, Pathways)")
    print("  5. Discussion")
    print("  6. Plot and Table Interpretation Guide")
    print("  + Appendix: File Inventory")


if __name__ == "__main__":
    main()
