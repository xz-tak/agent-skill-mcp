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

### Step 3: Present Matches with AskUserQuestion

Group results by granularity. Present as multi-select options using AskUserQuestion. Include in each option label: `combo`, `stage`, `consensus_type`, `n_studies`.

**Recommendations to show the user:**
- Prefer **fine** granularity (more cell-type specific)
- Prefer **multi_study** consensus (cross-validated, more robust)
- Mention total file sizes from `gene_region_size_mb + gene_tf_size_mb + region_tf_size_mb`

If >4 matches per granularity, summarize the list and let the user narrow further or pick from a representative set.

### Step 4: Download Selected Combos

For each selected combo, download all 3 parquet files:

```python
import subprocess

combo_row = matches[matches["combo"] == selected_combo].iloc[0]
dest_dir = f"/tmp/scatac_grn/{combo_row['combo']}"
os.makedirs(dest_dir, exist_ok=True)

for col in ["gene_region_s3", "gene_tf_s3", "region_tf_s3"]:
    s3_path = S3_BASE + combo_row[col]
    local_path = os.path.join(dest_dir, os.path.basename(combo_row[col]))
    subprocess.run(["aws", "s3", "cp", s3_path, local_path], check=True)
```

### Step 5: Load and Summarize

```python
import pandas as pd

gene_tf = pd.read_parquet(f"{dest_dir}/{combo}_gene_tf.parquet")
gene_region = pd.read_parquet(f"{dest_dir}/{combo}_gene_region.parquet")
region_tf = pd.read_parquet(f"{dest_dir}/{combo}_region_tf.parquet")

print(f"gene_tf: {gene_tf.shape} — {gene_tf.columns.tolist()}")
print(f"gene_region: {gene_region.shape} — {gene_region.columns.tolist()}")
print(f"region_tf: {region_tf.shape} — {region_tf.columns.tolist()}")

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

- At least one combo's parquets must already be downloaded (gene_tf + gene_region in `/tmp/scatac_grn/{combo}/`)
- If multiple combos are downloaded, use **AskUserQuestion** to confirm which combo to use

### Workflow

#### 1. Confirm combo

If multiple combos have been downloaded this session, ask user which one to use for the plot.

#### 2. Get target gene and TFs from user

- User specifies target gene (e.g., "NR2F2") and TFs (e.g., "TWIST1, SOX6, PPARA::RXRA")
- If user says "top" or "all" TFs, pick top 5 by max weight for that gene

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

Run the script at `skill_test/generate_tracks_plot.py` as reference, adapting parameters. The plot has 4 tracks (top to bottom):

1. **Arc track** — co-accessibility arcs from distal CREs to promoter, colored per TF. Gray dotted arcs for distal CREs with no specified TF binding (shows the link exists but other TFs regulate there).
2. **TF binding track** — bar height = motif_score, one band per TF (stacked). No TF name labels on the panel (colors identify them via legend).
3. **CRE peaks** — blue=promoter, green=distal. Red/colored border + dot above = TF bound. "no match" annotation below unbound distal CREs.
4. **Gene body** — blue rectangle with gene name in **black** text + TSS direction arrow.

**Key features:**
- **Gap compression**: If peaks span vastly different scales (max gap > 10× cluster span), compress the large gap to 5% of original. Adds "//" break indicator and dashed vertical lines. All tracks share one x-axis (unified zoom/pan).
- **Light theme**: white background, dark text
- **No crop**: height=1000, margin.b=190, legend at y=-0.14

```python
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os

# Parameters (set from user input)
combo_name = "COMBO_NAME"       # from user selection
target_gene = "GENE_NAME"      # from user
target_tfs = ["TF1", "TF2"]   # validated list
output_dir = os.getcwd()       # or user-specified

# Color palette (up to 10 TFs)
TF_PALETTE = ['#e63946', '#f4a261', '#7c3aed', '#2a9d8f', '#264653',
              '#e76f51', '#219ebc', '#8338ec', '#fb5607', '#3a86ff']
TF_COLORS = {tf: TF_PALETTE[i % len(TF_PALETTE)] for i, tf in enumerate(target_tfs)}

BG = 'white'; TEXT = '#1a1a2e'; MUTED = '#555555'; BLUE = '#2196F3'; GREEN = '#4CAF50'
GAP_COMPRESS_RATIO = 0.05; SPLIT_THRESHOLD_RATIO = 10

# Load data
dest_dir = f"/tmp/scatac_grn/{combo_name}"
gene_tf = pd.read_parquet(f"{dest_dir}/{combo_name}_gene_tf.parquet")
gene_region = pd.read_parquet(f"{dest_dir}/{combo_name}_gene_region.parquet")

# TF binding info
tf_data = {}
for tf in target_tfs:
    tf_edges = gene_tf[(gene_tf["tf"] == tf) & (gene_tf["gene"] == target_gene)]
    tf_data[tf] = {"peaks": set(tf_edges["peak_id"].tolist()),
        "scores": tf_edges.set_index("peak_id")["motif_score"].to_dict(),
        "weights": tf_edges.set_index("peak_id")["weight"].to_dict()}

# Parse peaks
def parse_peak(pid):
    parts = pid.split("_"); return parts[0], int(parts[1]), int(parts[2])

gr = gene_region[gene_region["gene"] == target_gene].copy().reset_index(drop=True)
gr["chr"] = gr["peak_id"].apply(lambda p: parse_peak(p)[0])
gr["start"] = gr["peak_id"].apply(lambda p: parse_peak(p)[1])
gr["end"] = gr["peak_id"].apply(lambda p: parse_peak(p)[2])
main_chr = gr["chr"].value_counts().index[0]
gr = gr[gr["chr"] == main_chr].copy().sort_values("start").reset_index(drop=True)
promoter_peaks = gr[gr["link_type"] == "promoter"]
distal_peaks = gr[gr["link_type"] == "distal"]
tss_mid = int((promoter_peaks["end"].max() + promoter_peaks.loc[promoter_peaks["end"].idxmax(), "start"]) / 2)

# Gap compression (for widely-spaced peaks)
gaps = gr["start"].diff().dropna().values
max_gap = gaps.max(); cluster_span = (gr["end"].max() - gr["start"].min()) - max_gap
need_compress = max_gap > cluster_span * SPLIT_THRESHOLD_RATIO
gap_start = gap_end = None; compressed_gap_size = 0
if need_compress:
    gi = int(np.argmax(gaps))
    gap_start = int(gr.iloc[gi]["end"]) + 500; gap_end = int(gr.iloc[gi+1]["start"]) - 500
    compressed_gap_size = (gap_end - gap_start) * GAP_COMPRESS_RATIO

def cx(x):
    if not need_compress: return float(x)
    if x <= gap_start: return float(x)
    elif x >= gap_end: return float(gap_start + compressed_gap_size + (x - gap_end))
    else: return float(gap_start + ((x - gap_start) / (gap_end - gap_start)) * compressed_gap_size)
def cx_array(arr): return np.array([cx(x) for x in arr])

chrom_start_c = min(cx(s) for s in gr["start"]) - 3000
chrom_end_c = max(cx(e) for e in gr["end"]) + 3000
tss_mid_c = cx(tss_mid)

# Build figure (single axis, unified zoom)
fig = make_subplots(rows=4, cols=1, shared_xaxes=True,
    row_heights=[0.33, 0.25, 0.22, 0.10], vertical_spacing=0.05)

# Track 1: Arcs — TF-bound (colored) + unbound (gray dotted)
fig.add_trace(go.Scatter(x=[chrom_start_c, chrom_end_c], y=[0, 0],
    mode="lines", line=dict(width=0.5, color="rgba(0,0,0,0.1)"),
    hoverinfo="skip", showlegend=False), row=1, col=1)
for _, row in distal_peaks.iterrows():
    pm = (row["start"] + row["end"]) / 2; pmc = cx(pm); coa = row["coaccess_score"]
    btfs = [tf for tf in target_tfs if row["peak_id"] in tf_data[tf]["peaks"]]
    if btfs:
        for ti, tf in enumerate(btfs):
            off = 1.0 + ti * 0.12; xc = cx_array(np.linspace(pm, tss_mid, 60))
            sc = abs(pmc - tss_mid_c)
            yn = ((sc/2) * np.sin(np.linspace(0, np.pi, 60)) * off) / ((chrom_end_c - chrom_start_c)/2)
            ms = tf_data[tf]["scores"][row["peak_id"]]; wt = tf_data[tf]["weights"][row["peak_id"]]
            fig.add_trace(go.Scatter(x=xc, y=yn, mode="lines",
                line=dict(width=max(2, coa*4.5), color=TF_COLORS[tf], shape="spline"),
                opacity=0.5+coa*0.4, hoverinfo="text",
                text=f"<b>{tf}</b> → {row['peak_id']}<br>coaccess: {coa:.3f}<br>motif: {ms:.4f}, weight: {wt:.4f}",
                showlegend=False), row=1, col=1)
    else:
        xc = cx_array(np.linspace(pm, tss_mid, 60)); sc = abs(pmc - tss_mid_c)
        yn = ((sc/2) * np.sin(np.linspace(0, np.pi, 60))) / ((chrom_end_c - chrom_start_c)/2)
        fig.add_trace(go.Scatter(x=xc, y=yn, mode="lines",
            line=dict(width=max(1.5, coa*3), color="rgba(150,150,150,0.4)", shape="spline", dash="dot"),
            hoverinfo="text", text=f"<b>No {'/'.join(target_tfs)} binding</b><br>{row['peak_id']}<br>coaccess: {coa:.3f}",
            showlegend=False), row=1, col=1)

# Track 2: TF Binding (no TF name labels)
bh = 1.0 / len(target_tfs)
for ti, tf in enumerate(target_tfs):
    yb = ti * bh
    for _, row in gr.iterrows():
        if row["peak_id"] not in tf_data[tf]["peaks"]: continue
        ms = tf_data[tf]["scores"][row["peak_id"]]; wt = tf_data[tf]["weights"][row["peak_id"]]
        bt = yb + ms * bh * 0.9; sc, ec = cx(row["start"]), cx(row["end"])
        fig.add_trace(go.Scatter(x=[sc,sc,ec,ec,sc], y=[yb,bt,bt,yb,yb],
            fill="toself", fillcolor=TF_COLORS[tf], line=dict(color=TF_COLORS[tf], width=0.5),
            mode="lines", opacity=0.8, hoverinfo="text",
            text=f"<b>{tf}</b><br>{row['peak_id']}<br>motif: {ms:.4f}<br>weight: {wt:.4f}",
            showlegend=False), row=2, col=1)

# Track 3: CRE Peaks + dots + "no match" labels
for _, row in gr.iterrows():
    ip = row["link_type"] == "promoter"
    btfs = [tf for tf in target_tfs if row["peak_id"] in tf_data[tf]["peaks"]]
    ha = len(btfs) > 0
    co = "rgba(33,150,243,0.5)" if ip else "rgba(76,175,80,0.5)"
    bd = TF_COLORS[btfs[0]] if ha else ("rgba(33,150,243,1)" if ip else "rgba(76,175,80,1)")
    sc, ec = cx(row["start"]), cx(row["end"])
    fig.add_trace(go.Scatter(x=[sc,sc,ec,ec,sc], y=[0.1,0.9,0.9,0.1,0.1],
        fill="toself", fillcolor=co, line=dict(color=bd, width=2.5 if ha else 1),
        mode="lines", hoverinfo="text",
        text=f"<b>{row['peak_id']}</b><br>{row['link_type']}<br>coaccess: {row['coaccess_score']:.4f}<br>TFs: {', '.join(btfs) if btfs else 'none'}",
        showlegend=False), row=3, col=1)
    mc = (sc + ec) / 2
    if btfs:
        sp = (chrom_end_c - chrom_start_c) * 0.01; sx = mc - (len(btfs)-1)*sp/2
        for di, tf in enumerate(btfs):
            fig.add_trace(go.Scatter(x=[sx+di*sp], y=[1.1], mode="markers",
                marker=dict(size=6, color=TF_COLORS[tf], symbol="circle"), hoverinfo="text",
                text=f"{tf}: {tf_data[tf]['scores'][row['peak_id']]:.4f}", showlegend=False), row=3, col=1)
    else:
        fig.add_annotation(x=mc, y=-0.15, text="no match", showarrow=False,
            font=dict(size=7, color="#999"), xref="x3", yref="y3")

# Track 4: Gene body (black gene symbol)
gs = int(promoter_peaks["start"].min()); ge = int(promoter_peaks["end"].max())
gsc, gec = cx(gs), cx(ge)
fig.add_trace(go.Scatter(x=[gsc,gsc,gec,gec,gsc], y=[0.2,0.8,0.8,0.2,0.2],
    fill="toself", fillcolor="rgba(21,101,192,0.85)", line=dict(color="rgba(13,71,161,1)", width=2),
    mode="lines", hoverinfo="text", text=f"<b>{target_gene}</b><br>{main_chr}:{gs:,}–{ge:,}",
    showlegend=False), row=4, col=1)
fig.add_annotation(x=(gsc+gec)/2, y=0.5, text=f"<b>{target_gene}</b>", showarrow=False,
    font=dict(size=13, color="black"), xref="x4", yref="y4")
fig.add_annotation(x=gec+(chrom_end_c-chrom_start_c)*0.01, y=0.5,
    ax=gec+(chrom_end_c-chrom_start_c)*0.04, ay=0.5,
    xref="x4", yref="y4", axref="x4", ayref="y4",
    showarrow=True, arrowhead=3, arrowsize=1, arrowwidth=1.5, arrowcolor="#1565C0")

# Gap break indicator (if compressed)
if need_compress:
    gmc = cx((gap_start + gap_end) / 2)
    for ri in range(1, 5):
        fig.add_vline(x=cx(gap_start), row=ri, col=1, line=dict(width=1, color="rgba(0,0,0,0.15)", dash="dash"))
        fig.add_vline(x=cx(gap_end), row=ri, col=1, line=dict(width=1, color="rgba(0,0,0,0.15)", dash="dash"))
    fig.add_annotation(x=gmc, y=0.5, text="<b>//</b>", showarrow=False,
        font=dict(size=20, color="#bbb"), xref="x", yref="paper")
    fig.add_annotation(x=gmc, y=-0.03, text=f"{(gap_end-gap_start)/1000:.0f} kb",
        showarrow=False, font=dict(size=9, color=MUTED), xref="x4", yref="paper")

# Layout
fig.update_layout(height=1000, width=1000,
    title=dict(text=f"<b>{' + '.join(target_tfs)} → {target_gene}</b> Cis-Regulatory Landscape<br>"
        f"<span style='font-size:11px;color:{MUTED}'>{combo_name} | {main_chr}</span>",
        font=dict(size=15, color=TEXT), x=0.5),
    showlegend=False, hovermode="closest",
    plot_bgcolor=BG, paper_bgcolor=BG, font=dict(color=TEXT),
    margin=dict(l=90, r=30, t=80, b=190))
for i in range(1, 5):
    fig.update_xaxes(range=[chrom_start_c, chrom_end_c], showticklabels=(i==4), showgrid=False, zeroline=False, row=i, col=1)
# Custom ticks showing real genomic positions
tick_pos = sorted(set([int(gr["start"].min()), tss_mid, int(gr.iloc[-1]["start"]), int(gr.iloc[-1]["end"])]))
fig.update_xaxes(tickmode="array", tickvals=[cx(t) for t in tick_pos], ticktext=[f"{t:,}" for t in tick_pos],
    tickfont=dict(size=12, color=MUTED), title_text=f"Genomic Position ({main_chr})",
    title_font=dict(size=13, color=MUTED), showgrid=True, gridcolor="rgba(0,0,0,0.05)", row=4, col=1)
track_labels = ["Regulatory\nArcs", "TF Binding\n(motif score)", "CRE\nPeaks", "Gene"]
for i in range(1, 5):
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False,
        title_text=track_labels[i-1], title_font=dict(size=11, color=MUTED), title_standoff=5, row=i, col=1)

# Legend (well below plot, no overlap)
legend_items = ([(TF_COLORS[tf], tf) for tf in target_tfs] +
    [(BLUE, "Promoter CRE"), (GREEN, "Distal CRE"), ("#999", "No TF match (dotted)")])
for idx, (color, label) in enumerate(legend_items):
    rp = idx // 3; cp = idx % 3
    fig.add_annotation(x=0.02+cp*0.34, y=-0.14-rp*0.04,
        text=f"<span style='color:{color};font-size:16px'>■</span> {label}",
        showarrow=False, xref="paper", yref="paper", xanchor="left", font=dict(size=12, color=TEXT))

output_path = os.path.join(output_dir, f"{combo_name}_{target_gene}_tracks.html")
fig.write_html(output_path, include_plotlyjs=True)
print(f"Saved: {output_path}")
```

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
