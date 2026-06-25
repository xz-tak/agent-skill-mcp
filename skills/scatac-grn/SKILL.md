---
name: scatac-grn
description: >
  Query and retrieve base-GRN (gene regulatory network) prior networks from the
  Chen 2026 human scATAC atlas (1,028 tissue×celltype combos across 56 tissues,
  247 cell types, adult+fetal). Use when users want GRN data for specific tissues,
  cell types, or developmental stages, download gene×TF networks, or explore
  available combinations. Triggers on: "GRN for [tissue]", "gene regulatory network",
  "scATAC GRN", "base GRN", "download GRN", "what tissues have GRNs", "fetal/adult GRN",
  "gene×TF network", "Chen 2026", "CIRCE", regulatory network per tissue/cell type,
  or selecting combos for CellOracle/network analysis.
---

# scATAC-seq Base GRN Retrieval

Interactive discovery and download of tissue×cell-type gene regulatory networks from the Chen 2026 human scATAC-seq atlas.

## Overview

The Chen 2026 pipeline produced **1,028 base-GRN outputs** (645 fine-grained + 383 coarse) covering 56 tissues, 247 cell types, and 2 developmental stages (adult, fetal).

Each combo produces **3 parquet files**:

| File | Contents | Key Columns |
|------|----------|-------------|
| `gene_region` | Gene ↔ cis-regulatory region links | gene, region, score |
| `gene_tf` | Gene ↔ transcription factor regulation weights | gene, tf, importance |
| `region_tf` | Region ↔ TF binding predictions | region, tf, score |

### Granularity

- **fine** (645 combos): Specific cell types per tissue (e.g., "Brodmann area 10 × astrocyte")
- **coarse** (383 combos): Broader cell groupings per tissue (e.g., "brain × neuron")

### Consensus Type

- **multi_study**: Cross-validated across multiple studies — more robust
- **single_study**: Derived from one source — broader coverage but less validated

## S3 Configuration

```
Base URI:   s3://tec-rnd-sci-dev-gi2/gi2-xz/insights/project/Chen_2026/
Index:      grn_pipeline/base_grn_index.csv
Parquets:   Paths from index columns (gene_region_s3, gene_tf_s3, region_tf_s3)
```

All S3 paths in the index are **relative** to the base URI.

## Workflow

### Step 1: Load Index

Download the index CSV from S3 (cache at `/tmp/scatac_grn/base_grn_index.csv`, re-download if >24h old):

```python
import pandas as pd
import os, time

CACHE = "/tmp/scatac_grn/base_grn_index.csv"
S3_BASE = "s3://tec-rnd-sci-dev-gi2/gi2-xz/insights/project/Chen_2026/"
S3_INDEX = S3_BASE + "grn_pipeline/base_grn_index.csv"

os.makedirs("/tmp/scatac_grn", exist_ok=True)

# Re-download if missing or stale (>24h)
if not os.path.exists(CACHE) or (time.time() - os.path.getmtime(CACHE)) > 86400:
    os.system(f"aws s3 cp {S3_INDEX} {CACHE}")

idx = pd.read_csv(CACHE)
```

### Step 2: Filter

Apply **case-insensitive substring** matching on `tissue` and `cell_type`. Use **exact match** on `stage` (adult/fetal) and `consensus_type` (multi_study/single_study) when specified.

Always search BOTH fine and coarse granularities and present results from each.

```python
query_tissue = "brain"  # user's tissue query
query_cell = "astrocyte"  # user's cell type query (optional)

mask = idx["tissue"].str.lower().str.contains(query_tissue.lower(), na=False)
if query_cell:
    mask &= idx["cell_type"].str.lower().str.contains(query_cell.lower(), na=False)
# Optionally filter stage:
# mask &= idx["stage"] == "adult"

matches = idx[mask].sort_values(["granularity", "consensus_type", "tissue", "cell_type"])
```

### Step 3: MANDATORY — Confirm Combo with AskUserQuestion

**⚠️ BLOCKING: You MUST use AskUserQuestion to confirm the combo BEFORE any S3 data read or parquet loading. NEVER skip this step, even if the user's request seems unambiguous.**

After filtering the index, present matching combos via AskUserQuestion. Include the **best-matching option first** (marked Recommended) based on:
1. Prefer **fine** granularity (more cell-type specific)
2. Prefer **multi_study** consensus (cross-validated, more robust)
3. If tied, prefer higher `n_studies`

Each option label should include: `granularity | stage | consensus_type | n_studies`.
Each option description should include: the full combo name and total size.

```
AskUserQuestion:
  question: "Which {tissue} × {cell_type} combo should we use?"
  options:
    - label: "fine | adult | multi_study | 3 studies (Recommended)"
      description: "combo_name — 530 MB total"
    - label: "coarse | adult | single_study | 1 study"
      description: "combo_name — 480 MB total"
```

If >4 matches, show the top 4 most relevant (by the ranking above) and mention how many more exist.

**Only proceed to Step 4 after user confirms.**

### Step 4: Load Data from S3

Read parquets directly from S3 using pandas (no download required):

```python
import pandas as pd

S3_BASE = "s3://tec-rnd-sci-dev-gi2/gi2-xz/insights/project/Chen_2026/"
combo_row = matches[matches["combo"] == selected_combo].iloc[0]

gene_tf = pd.read_parquet(S3_BASE + combo_row["gene_tf_s3"])
gene_region = pd.read_parquet(S3_BASE + combo_row["gene_region_s3"])
```

### Step 5: Summarize

```python
print(f"gene_tf: {gene_tf.shape} — {gene_tf.columns.tolist()}")
print(f"gene_region: {gene_region.shape} — {gene_region.columns.tolist()}")

# Top TFs by frequency
if "tf" in gene_tf.columns:
    print("\nTop 10 TFs:")
    print(gene_tf["tf"].value_counts().head(10))
```

Leave DataFrames in the Python session for the user's downstream analysis.

After summary, cue the user:
```
print("\nTo visualize cis-regulatory interactions for a target gene, specify the gene and TFs of interest.")
```

## Index Columns Reference

| Column | Description |
|--------|-------------|
| `granularity` | "fine" or "coarse" |
| `combo` | Unique identifier: `{tissue}__{cell_type}` (underscores for spaces) |
| `tissue` | Human-readable tissue name |
| `cell_type` | Human-readable cell type |
| `stage` | "adult" or "fetal" |
| `consensus_type` | "multi_study" or "single_study" |
| `n_studies` | Number of contributing studies |
| `gene_region_s3` | Relative S3 path to gene↔region parquet |
| `gene_tf_s3` | Relative S3 path to gene↔TF parquet |
| `region_tf_s3` | Relative S3 path to region↔TF parquet |
| `gene_region_size_mb` | File size in MB |
| `gene_tf_size_mb` | File size in MB |
| `region_tf_size_mb` | File size in MB |
| `has_all_3_files` | Boolean — True if all 3 parquets exist |

## Usage Patterns

**Listing available tissues/cell types:**
```python
print(idx["tissue"].nunique(), "tissues")
print(idx["cell_type"].nunique(), "cell types")
print(idx.groupby("granularity").size())
```

**Fetal-only search:**
```python
fetal = idx[idx["stage"] == "fetal"]
```

**Multi-study only (recommended for robustness):**
```python
robust = idx[idx["consensus_type"] == "multi_study"]
```

## Recommendations

| Goal | Recommendation |
|------|---------------|
| Specific regulatory programs | Use **fine** granularity |
| Broader tissue-level patterns | Use **coarse** granularity |
| Higher confidence networks | Use **multi_study** consensus |
| Maximum tissue/cell coverage | Use **single_study** (more combos available) |
| CellOracle input | Download gene_tf parquet → convert to base GRN format |
| Network visualization | Use gene_tf for TF→gene edges, region_tf for TF→region edges |
| Cis-regulatory plot | See **Visualization** section below |

## Visualization: Cis-Regulatory Track Plot

Generate an interactive Plotly HTML showing TF binding to a target gene's regulatory landscape. Triggered **only** when user explicitly requests to visualize cis-regulatory interactions (e.g., "show regulation of NR2F2 by TWIST1", "plot cis-regulatory interactions", "visualize TF binding to IGFBP7").

### Prerequisites

- The script reads parquets directly from S3 (no local download required). If parquets are already local, pass `--data-dir` to skip S3.
- If multiple combos are relevant, use **AskUserQuestion** to confirm which combo to use

### Workflow

**⚠️ BLOCKING: Steps 1–2 use AskUserQuestion and MUST complete BEFORE any S3 read or plot generation.**

#### 1. Confirm combo

Load the index and filter for the user's tissue/cell type. Use **AskUserQuestion** to present matching combos (best match first, marked Recommended). NEVER assume a combo — always confirm.

#### 2. Get target gene and TFs from user

- User specifies target gene (e.g., "NR2F2") and TFs (e.g., "TWIST1, SOX6, PPARA::RXRA")
- If user says "top" or "all" TFs: load gene_tf from S3, find top 5 by max weight for that gene, then present them via **AskUserQuestion** for confirmation before plotting

#### 3. Validate

```python
# Check gene exists
gene_edges = gene_tf[gene_tf["gene"] == target_gene]
if len(gene_edges) == 0:
    print(f"ERROR: {target_gene} not found in gene_tf. Check spelling (case-sensitive).")
    # abort

# Check each TF
valid_tfs = []
missing_tfs = []
for tf in target_tfs:
    tf_edges = gene_edges[gene_edges["tf"] == tf]
    if len(tf_edges) > 0:
        valid_tfs.append(tf)
    else:
        missing_tfs.append(tf)
```

If any TFs are missing, use **AskUserQuestion** to show which are valid/invalid and let user confirm whether to proceed with remaining TFs. Abort only if ALL are missing.

#### 4. Generate plot

Run the standalone plot script:

```bash
python ~/.claude/skills/scatac-grn/generate_tracks_plot.py \
    --combo "{combo_name}" \
    --gene "{target_gene}" \
    --tfs "{tf1},{tf2},{tf3}" \
    --output-dir "{output_dir or cwd}"
```

Optional `--data-dir` uses local parquets instead of reading from S3 (faster if already downloaded).

The script produces an interactive Plotly HTML (`{combo}_{gene}_tracks.html`) with 4 tracks:

1. **Arc track** — co-accessibility arcs from distal CREs to promoter, colored per TF. Gray dotted arcs for distal CREs with no specified TF binding.
2. **TF binding track** — bar height = motif_score, one band per TF. Legend provides color key.
3. **CRE peaks** — blue=promoter, green=distal. Colored border + dot above = TF bound. "no match" below unbound peaks.
4. **Gene body** — blue rectangle with gene name in black + TSS direction arrow.

**Key features:** gap compression (auto-activates for sparse peaks), light theme, unified zoom across all tracks, ~4.7 MB self-contained HTML.

### Notes

- **Light theme** (white background) for paper/presentation readability
- **Gap compression** auto-activates when max inter-peak gap > 10× cluster span (e.g., GREM1). Adds "//" break, dashed lines, and distance label. Unified zoom works across all tracks.
- **Gray dotted arcs** show co-accessible distal CREs where none of the specified TFs bind — indicates regulation by other TFs
- **"no match"** label below unbound CRE peaks
- **Gene symbol always in black** for visibility on the blue gene body
- **No TF labels** on the TF binding panel (legend provides color key)
- Each TF gets a distinct color from the palette (up to 10 TFs)
- Colored dots above CRE peaks indicate which TFs bind at each site
- Output: `{combo}_{gene}_tracks.html` in working directory (~4.7 MB, fully self-contained, works offline)
