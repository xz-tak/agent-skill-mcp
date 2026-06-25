#!/usr/bin/env python3
"""Generate cis-regulatory track plot for scATAC-seq GRN data.

Produces an interactive Plotly HTML with 4 tracks:
  1. Regulatory arcs (TF-bound colored, unbound gray dotted)
  2. TF binding bars (motif_score height)
  3. CRE peaks (promoter blue, distal green, dots + "no match")
  4. Gene body (blue rect, black gene symbol, TSS arrow)

Features: gap compression, light theme, unified zoom, self-contained HTML.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


S3_BASE = "s3://tec-rnd-sci-dev-gi2/gi2-xz/insights/project/Chen_2026/"
S3_INDEX = S3_BASE + "grn_pipeline/base_grn_index.csv"
INDEX_CACHE = "/tmp/scatac_grn/base_grn_index.csv"

TF_PALETTE = [
    '#e63946', '#f4a261', '#7c3aed', '#2a9d8f', '#264653',
    '#e76f51', '#219ebc', '#8338ec', '#fb5607', '#3a86ff',
]
BG = 'white'
TEXT = '#1a1a2e'
MUTED = '#555555'
BLUE = '#2196F3'
GREEN = '#4CAF50'
GAP_COMPRESS_RATIO = 0.05
SPLIT_THRESHOLD_RATIO = 10


def parse_peak(pid):
    parts = pid.split("_")
    return parts[0], int(parts[1]), int(parts[2])


def load_parquets(combo_name, data_dir=None):
    """Load gene_tf and gene_region parquets from local dir or S3."""
    if data_dir and os.path.isdir(data_dir):
        gene_tf = pd.read_parquet(os.path.join(data_dir, f"{combo_name}_gene_tf.parquet"))
        gene_region = pd.read_parquet(os.path.join(data_dir, f"{combo_name}_gene_region.parquet"))
        return gene_tf, gene_region

    import time
    os.makedirs(os.path.dirname(INDEX_CACHE), exist_ok=True)
    if not os.path.exists(INDEX_CACHE) or (time.time() - os.path.getmtime(INDEX_CACHE)) > 86400:
        idx = pd.read_csv(S3_INDEX)
        idx.to_csv(INDEX_CACHE, index=False)
    else:
        idx = pd.read_csv(INDEX_CACHE)

    rows = idx[idx["combo"] == combo_name]
    if len(rows) == 0:
        print(f"ERROR: combo '{combo_name}' not found in index.")
        print(f"Available combos (sample): {idx['combo'].head(10).tolist()}")
        sys.exit(1)
    row = rows.iloc[0]
    print(f"Reading from S3: {combo_name} ({row['granularity']}, {row['consensus_type']})")
    gene_tf = pd.read_parquet(S3_BASE + row["gene_tf_s3"])
    gene_region = pd.read_parquet(S3_BASE + row["gene_region_s3"])
    return gene_tf, gene_region


def build_figure(combo_name, target_gene, target_tfs, gene_tf, gene_region, output_dir):
    tf_colors = {tf: TF_PALETTE[i % len(TF_PALETTE)] for i, tf in enumerate(target_tfs)}

    gene_edges = gene_tf[gene_tf["gene"] == target_gene]
    if len(gene_edges) == 0:
        print(f"ERROR: '{target_gene}' not found in gene_tf. Available genes (sample):")
        print(gene_tf["gene"].drop_duplicates().sort_values().head(20).tolist())
        sys.exit(1)

    valid_tfs = []
    missing_tfs = []
    for tf in target_tfs:
        if len(gene_edges[gene_edges["tf"] == tf]) > 0:
            valid_tfs.append(tf)
        else:
            missing_tfs.append(tf)

    if missing_tfs:
        print(f"WARNING: TFs not found for {target_gene}: {missing_tfs}")
    if not valid_tfs:
        print(f"ERROR: None of the specified TFs regulate {target_gene}. Aborting.")
        sys.exit(1)

    target_tfs = valid_tfs
    tf_colors = {tf: TF_PALETTE[i % len(TF_PALETTE)] for i, tf in enumerate(target_tfs)}

    tf_data = {}
    for tf in target_tfs:
        tf_edges = gene_tf[(gene_tf["tf"] == tf) & (gene_tf["gene"] == target_gene)]
        tf_data[tf] = {
            "peaks": set(tf_edges["peak_id"].tolist()),
            "scores": tf_edges.set_index("peak_id")["motif_score"].to_dict(),
            "weights": tf_edges.set_index("peak_id")["weight"].to_dict(),
        }

    gr = gene_region[gene_region["gene"] == target_gene].copy().reset_index(drop=True)
    gr["chr"] = gr["peak_id"].apply(lambda p: parse_peak(p)[0])
    gr["start"] = gr["peak_id"].apply(lambda p: parse_peak(p)[1])
    gr["end"] = gr["peak_id"].apply(lambda p: parse_peak(p)[2])
    main_chr = gr["chr"].value_counts().index[0]
    gr = gr[gr["chr"] == main_chr].copy().sort_values("start").reset_index(drop=True)
    promoter_peaks = gr[gr["link_type"] == "promoter"]
    distal_peaks = gr[gr["link_type"] == "distal"]
    tss_mid = int(
        (promoter_peaks["end"].max() + promoter_peaks.loc[promoter_peaks["end"].idxmax(), "start"]) / 2
    )

    # Gap compression
    gaps = gr["start"].diff().dropna().values
    max_gap = gaps.max()
    cluster_span = (gr["end"].max() - gr["start"].min()) - max_gap
    need_compress = max_gap > cluster_span * SPLIT_THRESHOLD_RATIO
    gap_start = gap_end = None
    compressed_gap_size = 0
    if need_compress:
        gi = int(np.argmax(gaps))
        gap_start = int(gr.iloc[gi]["end"]) + 500
        gap_end = int(gr.iloc[gi + 1]["start"]) - 500
        compressed_gap_size = (gap_end - gap_start) * GAP_COMPRESS_RATIO

    def cx(x):
        if not need_compress:
            return float(x)
        if x <= gap_start:
            return float(x)
        elif x >= gap_end:
            return float(gap_start + compressed_gap_size + (x - gap_end))
        else:
            return float(gap_start + ((x - gap_start) / (gap_end - gap_start)) * compressed_gap_size)

    def cx_array(arr):
        return np.array([cx(x) for x in arr])

    chrom_start_c = min(cx(s) for s in gr["start"]) - 3000
    chrom_end_c = max(cx(e) for e in gr["end"]) + 3000
    tss_mid_c = cx(tss_mid)

    fig = make_subplots(
        rows=4, cols=1, shared_xaxes=True,
        row_heights=[0.33, 0.25, 0.22, 0.10], vertical_spacing=0.05,
    )

    # Track 1: Arcs
    fig.add_trace(go.Scatter(
        x=[chrom_start_c, chrom_end_c], y=[0, 0],
        mode="lines", line=dict(width=0.5, color="rgba(0,0,0,0.1)"),
        hoverinfo="skip", showlegend=False,
    ), row=1, col=1)

    for _, row in distal_peaks.iterrows():
        pm = (row["start"] + row["end"]) / 2
        pmc = cx(pm)
        coa = row["coaccess_score"]
        btfs = [tf for tf in target_tfs if row["peak_id"] in tf_data[tf]["peaks"]]
        if btfs:
            for ti, tf in enumerate(btfs):
                off = 1.0 + ti * 0.12
                xc = cx_array(np.linspace(pm, tss_mid, 60))
                sc = abs(pmc - tss_mid_c)
                yn = ((sc / 2) * np.sin(np.linspace(0, np.pi, 60)) * off) / ((chrom_end_c - chrom_start_c) / 2)
                ms = tf_data[tf]["scores"][row["peak_id"]]
                wt = tf_data[tf]["weights"][row["peak_id"]]
                fig.add_trace(go.Scatter(
                    x=xc, y=yn, mode="lines",
                    line=dict(width=max(2, coa * 4.5), color=tf_colors[tf], shape="spline"),
                    opacity=0.5 + coa * 0.4, hoverinfo="text",
                    text=f"<b>{tf}</b> → {row['peak_id']}<br>coaccess: {coa:.3f}<br>motif: {ms:.4f}, weight: {wt:.4f}",
                    showlegend=False,
                ), row=1, col=1)
        else:
            xc = cx_array(np.linspace(pm, tss_mid, 60))
            sc = abs(pmc - tss_mid_c)
            yn = ((sc / 2) * np.sin(np.linspace(0, np.pi, 60))) / ((chrom_end_c - chrom_start_c) / 2)
            fig.add_trace(go.Scatter(
                x=xc, y=yn, mode="lines",
                line=dict(width=max(1.5, coa * 3), color="rgba(150,150,150,0.4)", shape="spline", dash="dot"),
                hoverinfo="text",
                text=f"<b>No {'/'.join(target_tfs)} binding</b><br>{row['peak_id']}<br>coaccess: {coa:.3f}",
                showlegend=False,
            ), row=1, col=1)

    # Track 2: TF Binding
    bh = 1.0 / len(target_tfs)
    for ti, tf in enumerate(target_tfs):
        yb = ti * bh
        for _, row in gr.iterrows():
            if row["peak_id"] not in tf_data[tf]["peaks"]:
                continue
            ms = tf_data[tf]["scores"][row["peak_id"]]
            wt = tf_data[tf]["weights"][row["peak_id"]]
            bt = yb + ms * bh * 0.9
            sc, ec = cx(row["start"]), cx(row["end"])
            fig.add_trace(go.Scatter(
                x=[sc, sc, ec, ec, sc], y=[yb, bt, bt, yb, yb],
                fill="toself", fillcolor=tf_colors[tf],
                line=dict(color=tf_colors[tf], width=0.5),
                mode="lines", opacity=0.8, hoverinfo="text",
                text=f"<b>{tf}</b><br>{row['peak_id']}<br>motif: {ms:.4f}<br>weight: {wt:.4f}",
                showlegend=False,
            ), row=2, col=1)

    # Track 3: CRE Peaks
    for _, row in gr.iterrows():
        ip = row["link_type"] == "promoter"
        btfs = [tf for tf in target_tfs if row["peak_id"] in tf_data[tf]["peaks"]]
        ha = len(btfs) > 0
        co = "rgba(33,150,243,0.5)" if ip else "rgba(76,175,80,0.5)"
        bd = tf_colors[btfs[0]] if ha else ("rgba(33,150,243,1)" if ip else "rgba(76,175,80,1)")
        sc, ec = cx(row["start"]), cx(row["end"])
        fig.add_trace(go.Scatter(
            x=[sc, sc, ec, ec, sc], y=[0.1, 0.9, 0.9, 0.1, 0.1],
            fill="toself", fillcolor=co,
            line=dict(color=bd, width=2.5 if ha else 1),
            mode="lines", hoverinfo="text",
            text=f"<b>{row['peak_id']}</b><br>{row['link_type']}<br>coaccess: {row['coaccess_score']:.4f}<br>TFs: {', '.join(btfs) if btfs else 'none'}",
            showlegend=False,
        ), row=3, col=1)
        mc = (sc + ec) / 2
        if btfs:
            sp = (chrom_end_c - chrom_start_c) * 0.01
            sx = mc - (len(btfs) - 1) * sp / 2
            for di, tf in enumerate(btfs):
                fig.add_trace(go.Scatter(
                    x=[sx + di * sp], y=[1.1], mode="markers",
                    marker=dict(size=6, color=tf_colors[tf], symbol="circle"),
                    hoverinfo="text",
                    text=f"{tf}: {tf_data[tf]['scores'][row['peak_id']]:.4f}",
                    showlegend=False,
                ), row=3, col=1)
        else:
            fig.add_annotation(
                x=mc, y=-0.15, text="no match", showarrow=False,
                font=dict(size=7, color="#999"), xref="x3", yref="y3",
            )

    # Track 4: Gene body
    gs = int(promoter_peaks["start"].min())
    ge = int(promoter_peaks["end"].max())
    gsc, gec = cx(gs), cx(ge)
    fig.add_trace(go.Scatter(
        x=[gsc, gsc, gec, gec, gsc], y=[0.2, 0.8, 0.8, 0.2, 0.2],
        fill="toself", fillcolor="rgba(21,101,192,0.85)",
        line=dict(color="rgba(13,71,161,1)", width=2),
        mode="lines", hoverinfo="text",
        text=f"<b>{target_gene}</b><br>{main_chr}:{gs:,}–{ge:,}",
        showlegend=False,
    ), row=4, col=1)
    fig.add_annotation(
        x=(gsc + gec) / 2, y=0.5, text=f"<b>{target_gene}</b>", showarrow=False,
        font=dict(size=13, color="black"), xref="x4", yref="y4",
    )
    fig.add_annotation(
        x=gec + (chrom_end_c - chrom_start_c) * 0.01, y=0.5,
        ax=gec + (chrom_end_c - chrom_start_c) * 0.04, ay=0.5,
        xref="x4", yref="y4", axref="x4", ayref="y4",
        showarrow=True, arrowhead=3, arrowsize=1, arrowwidth=1.5, arrowcolor="#1565C0",
    )

    # Gap break indicator
    if need_compress:
        gmc = cx((gap_start + gap_end) / 2)
        for ri in range(1, 5):
            fig.add_vline(x=cx(gap_start), row=ri, col=1, line=dict(width=1, color="rgba(0,0,0,0.15)", dash="dash"))
            fig.add_vline(x=cx(gap_end), row=ri, col=1, line=dict(width=1, color="rgba(0,0,0,0.15)", dash="dash"))
        fig.add_annotation(
            x=gmc, y=0.5, text="<b>//</b>", showarrow=False,
            font=dict(size=20, color="#bbb"), xref="x", yref="paper",
        )
        fig.add_annotation(
            x=gmc, y=-0.03, text=f"{(gap_end - gap_start) / 1000:.0f} kb",
            showarrow=False, font=dict(size=9, color=MUTED), xref="x4", yref="paper",
        )

    # Layout
    fig.update_layout(
        height=1000, width=1000,
        title=dict(
            text=(
                f"<b>{' + '.join(target_tfs)} → {target_gene}</b> Cis-Regulatory Landscape<br>"
                f"<span style='font-size:11px;color:{MUTED}'>{combo_name} | {main_chr}</span>"
            ),
            font=dict(size=15, color=TEXT), x=0.5,
        ),
        showlegend=False, hovermode="closest",
        plot_bgcolor=BG, paper_bgcolor=BG, font=dict(color=TEXT),
        margin=dict(l=90, r=30, t=80, b=190),
    )
    for i in range(1, 5):
        fig.update_xaxes(
            range=[chrom_start_c, chrom_end_c],
            showticklabels=(i == 4), showgrid=False, zeroline=False, row=i, col=1,
        )

    tick_pos = sorted(set([
        int(gr["start"].min()), tss_mid, int(gr.iloc[-1]["start"]), int(gr.iloc[-1]["end"]),
    ]))
    fig.update_xaxes(
        tickmode="array", tickvals=[cx(t) for t in tick_pos],
        ticktext=[f"{t:,}" for t in tick_pos],
        tickfont=dict(size=12, color=MUTED),
        title_text=f"Genomic Position ({main_chr})",
        title_font=dict(size=13, color=MUTED),
        showgrid=True, gridcolor="rgba(0,0,0,0.05)", row=4, col=1,
    )

    track_labels = ["Regulatory\nArcs", "TF Binding\n(motif score)", "CRE\nPeaks", "Gene"]
    for i in range(1, 5):
        fig.update_yaxes(
            showticklabels=False, showgrid=False, zeroline=False,
            title_text=track_labels[i - 1],
            title_font=dict(size=11, color=MUTED), title_standoff=5, row=i, col=1,
        )

    # Legend
    legend_items = (
        [(tf_colors[tf], tf) for tf in target_tfs]
        + [(BLUE, "Promoter CRE"), (GREEN, "Distal CRE"), ("#999", "No TF match (dotted)")]
    )
    for idx, (color, label) in enumerate(legend_items):
        rp = idx // 3
        cp = idx % 3
        fig.add_annotation(
            x=0.02 + cp * 0.34, y=-0.14 - rp * 0.04,
            text=f"<span style='color:{color};font-size:16px'>■</span> {label}",
            showarrow=False, xref="paper", yref="paper", xanchor="left",
            font=dict(size=12, color=TEXT),
        )

    output_path = os.path.join(output_dir, f"{combo_name}_{target_gene}_tracks.html")
    fig.write_html(output_path, include_plotlyjs=True)
    print(f"Saved: {output_path}")
    return output_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate cis-regulatory track plot from scATAC-seq GRN parquets.",
    )
    parser.add_argument("--combo", required=True, help="Combo name (e.g., small_intestine__fibroblast)")
    parser.add_argument("--gene", required=True, help="Target gene symbol (e.g., NR2F2)")
    parser.add_argument("--tfs", required=True, help="Comma-separated TF list (e.g., TWIST1,SOX6,PPARA::RXRA)")
    parser.add_argument("--data-dir", default=None, help="Path to local parquets (if omitted, reads directly from S3)")
    parser.add_argument("--output-dir", default=None, help="Output directory for HTML (default: cwd)")
    args = parser.parse_args()

    combo_name = args.combo
    target_gene = args.gene
    target_tfs = [tf.strip() for tf in args.tfs.split(",")]
    output_dir = args.output_dir or os.getcwd()

    os.makedirs(output_dir, exist_ok=True)
    gene_tf, gene_region = load_parquets(combo_name, args.data_dir)
    build_figure(combo_name, target_gene, target_tfs, gene_tf, gene_region, output_dir)


if __name__ == "__main__":
    main()
