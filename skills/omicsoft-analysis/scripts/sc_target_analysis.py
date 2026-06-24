"""
Single-Cell Target Expression Analysis Script
Processes multiple scRNA-seq studies from S3, performing target gene expression
and signature enrichment analysis with static PDF reports and interactive HTML outputs.

Usage:
    conda activate spatial
    python sc_target_analysis.py --targets "IL11,IL11RA" --signatures "IL11_sig:IL13RA2,CEMIP,MMP3" --output-dir ./results
"""

import os
import sys
import argparse
import warnings
import tempfile
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import decoupler as dc
import boto3
import seaborn as sns
from scipy.stats import spearmanr, false_discovery_control, mannwhitneyu
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots

try:
    import xzsc as ci
    HAS_XZSC = True
except ImportError:
    HAS_XZSC = False

warnings.filterwarnings('ignore')


# ============================================================
# STUDY CONFIGURATIONS (hardcoded)
# ============================================================

STUDY_CONFIGS = {
    'jhu': {
        's3_path': 's3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/sc_ibd_mash_ssc_hs_internal_06102026/ibd/jhu_sc/jhu_patientsample_v2_scanpy.h5ad',
        'cell_type_col': 'cluster',
        'condition_col': 'comb_condition',
        'sample_col': 'Subject',
        'embedding_key': 'X_umap',
    },
    'rieder': {
        's3_path': 's3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/sc_ibd_mash_ssc_hs_internal_06102026/ibd/reider_sn/reider_bidmcanno_scanpy.h5ad',
        'cell_type_col': 'cluster',
        'condition_col': 'comb_condition',
        'sample_col': 'sample',
        'embedding_key': 'X_mde_scANVI_240923',
    },
    'umcg_fibroblast': {
        's3_path': 's3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/sc_ibd_mash_ssc_hs_internal_06102026/ibd/umcg_stromal_sc/umcg_fibroblasts_stromal_scanpy.h5ad',
        'cell_type_col': 'cluster',
        'condition_col': 'condition',
        'sample_col': 'sample',
        'embedding_key': 'X_umap',
    },
    'otar': {
        's3_path': 's3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/sc_ibd_mash_ssc_hs_internal_06102026/mash/otar2051/NASH_snucSeq_OTAR0251_scanpy.h5ad',
        'cell_type_col': 'Cell_type',
        'condition_col': 'condition',
        'sample_col': 'sample',
        'embedding_key': 'X_umap',
    },
    'tacolny': {
        's3_path': 's3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/sc_ibd_mash_ssc_hs_internal_06102026/mash/tacolny_schwabe/NASH_snucSeq_TaColNY_scanpy.h5ad',
        'cell_type_col': 'Cell_type',
        'condition_col': 'condition',
        'sample_col': 'sample',
        'embedding_key': 'X_umap',
    },
    'pf_lung': {
        's3_path': 's3://tec-rnd-sci-dev-gi2/gi2-xz/insights/project/TAK-999/PF_dataset/CellHint_lung_Adam2020Habermann2020Morse2019Reyfman2019_scanpy.h5ad',
        'cell_type_col': 'cluster',
        'condition_col': 'Disease',
        'sample_col': 'Donor',
        'embedding_key': 'X_umap',
    },
    'ssc_lung': {
        's3_path': 's3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/sc_ibd_mash_ssc_hs_internal_06102026/ssc/ssc_lung_atlas/integration_scanpy_11202025.h5ad',
        'cell_type_col': 'cluster',
        'condition_col': 'condition',
        'sample_col': 'orig.ident',
        'embedding_key': 'X_umap',
    },
    'hs_skin': {
        's3_path': 's3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/sc_ibd_mash_ssc_hs_internal_06102026/hs/rancho_hs_skin_atlas/SCDS_2025_Atlas-hidradenitis_suppurativa.h5ad',
        'cell_type_col': 'popv_majority_vote_prediction_ontology_name',
        'condition_col': 'condition',
        'sample_col': 'sample_id',
        'embedding_key': 'X_umap',
    },
    'ibd_atlas': {
        's3_path': 's3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/sc_ibd_mash_ssc_hs_internal_06102026/ibd/atlas042026_sc/integration_scanpy.h5ad',
        'cell_type_col': 'cluster',
        'condition_col': ['condition', 'disease', 'comb_condition'],
        'sample_col': 'sample',
        'embedding_key': 'X_umap',
        'comb_condition_cols': ['disease', 'condition'],
        'skip_aucell': True,
    },
    'umcg_change': {
        's3_path': 's3://tec-rnd-sci-dev-gi2/gi2-xz/omicsoft/sc_ibd_mash_ssc_hs_internal_06102026/ibd/umcg_change_sc/umcg_change_scanpy.h5ad',
        'cell_type_col': 'cluster',
        'condition_col': 'condition',
        'sample_col': 'sample',
        'embedding_key': 'X_umap',
    },
}


# ============================================================
# CLI ARGUMENT PARSING
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description='Single-Cell Target Expression Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Target-signature pairing format:
  Each --targets/--signatures pair is positionally matched (1st targets with 1st sig, etc.)
  Use semicolons to separate multiple sets within each argument.

  Example with 2 paired sets:
    --targets "IL11,IL11RA;TNFRSF25,TNFSF15"
    --signatures "IL11_sig:IL13RA2,CEMIP,MMP3;TL1A_sig:CHI3L1,MMP7,TCN1"

  This creates:
    - IL11_sig scored with: IL13RA2,CEMIP,MMP3 + IL11,IL11RA (paired targets)
    - TL1A_sig scored with: CHI3L1,MMP7,TCN1 + TNFRSF25,TNFSF15 (paired targets)

  If only 1 target set is given, it is unioned into ALL signatures (backward compat).
""")
    parser.add_argument('--targets', required=True,
                        help='Target genes. Multiple sets separated by semicolons: "IL11,IL11RA;TNFRSF25,TNFSF15"')
    parser.add_argument('--signatures', default=None,
                        help='Signatures paired with targets: "IL11_sig:IL13RA2,CEMIP;TL1A_sig:CHI3L1,MMP7"')
    parser.add_argument('--addon-signatures', default=None,
                        help='Additional signatures NOT paired with targets. Format: "Name1:Gene1,Gene2;Name2:GeneA,GeneB". '
                             'Scored independently via AUCell. Genes are NOT merged into target-paired signatures.')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for results')
    parser.add_argument('--sc-studies', default=None,
                        help='Comma-separated study names (default: all)')
    parser.add_argument('--max-workers', type=int, default=1,
                        help='Parallel workers (default: 1)')
    parser.add_argument('--temp-dir', default='/mnt/sagemaker-nvme/tmp',
                        help='Temp directory for S3 downloads')
    return parser.parse_args()


def build_signature_dict(targets_str, signatures_str=None, addon_signatures_str=None):
    """
    Parse CLI --targets, --signatures, and --addon-signatures into the signature_dict format.

    Pairing logic:
    - Target sets separated by semicolons: "set1_genes;set2_genes;set3_genes"
    - Signatures separated by semicolons: "sig1:genes;sig2:genes;sig3:genes"
    - Each target set is unioned into its positionally-matched signature
    - If only 1 target set is provided, it is unioned into ALL signatures (backward compat)
    - All targets (across all sets) are collected under signature_dict['target']
    - Addon signatures are added WITHOUT merging target genes
    """
    target_sets = [
        sorted(set(t.strip() for t in ts.split(',') if t.strip()))
        for ts in targets_str.split(';')
    ]
    all_targets = sorted(set(g for ts in target_sets for g in ts))
    signature_dict = {'target': all_targets}

    if signatures_str:
        sig_blocks = [s.strip() for s in signatures_str.split(';') if s.strip()]

        for i, sig_block in enumerate(sig_blocks):
            if ':' not in sig_block:
                continue
            name, genes_str = sig_block.split(':', 1)
            sig_genes = set(g.strip() for g in genes_str.split(',') if g.strip())

            if len(target_sets) == 1:
                sig_genes = sig_genes | set(target_sets[0])
            elif i < len(target_sets):
                sig_genes = sig_genes | set(target_sets[i])

            signature_dict[name.strip()] = sorted(sig_genes)

    if addon_signatures_str:
        for sig_block in (s.strip() for s in addon_signatures_str.split(';') if s.strip()):
            if ':' not in sig_block:
                continue
            name, genes_str = sig_block.split(':', 1)
            sig_genes = sorted(set(g.strip() for g in genes_str.split(',') if g.strip()))
            signature_dict[name.strip()] = sig_genes

    return signature_dict


# ============================================================
# CORE ANALYSIS FUNCTIONS
# ============================================================

def download_and_load_adata(s3_path, study_name, temp_dir):
    bucket_name = s3_path.split('/')[2]
    s3_key = '/'.join(s3_path.split('/')[3:])
    s3 = boto3.client('s3')

    os.makedirs(temp_dir, exist_ok=True)

    with tempfile.NamedTemporaryFile(suffix='.h5ad', dir=temp_dir, delete=True) as tmp:
        print(f"  [{study_name}] Downloading from S3...")
        s3.download_file(bucket_name, s3_key, tmp.name)
        print(f"  [{study_name}] Loading AnnData into memory...")
        adata = sc.read_h5ad(tmp.name)

    if 'log1p' in adata.uns:
        adata.uns['log1p']['base'] = None
    print(f"  [{study_name}] Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def resolve_embedding_key(adata, study_config, study_name):
    """Resolve the embedding key for a study, with fallback to X_umap."""
    embedding_key = study_config.get('embedding_key', 'X_umap')
    if embedding_key in adata.obsm:
        return embedding_key
    if embedding_key != 'X_umap' and 'X_umap' in adata.obsm:
        print(f"  [{study_name}] Embedding '{embedding_key}' not found, falling back to X_umap")
        return 'X_umap'
    if 'X_umap' in adata.obsm:
        return 'X_umap'
    print(f"  [{study_name}] No embedding found (tried '{embedding_key}' and 'X_umap'), skipping embedding plots")
    return None


def plot_corr_heatmap(adata, query_list, mode='auto', title='Target Correlation Heatmap (FDR Corrected)', figsize=(6, 4)):
    """
    Computes Spearman correlation for a list of genes or metadata (obs) columns
    and returns the matplotlib figure. Returns None if insufficient data.
    """
    if mode == 'auto':
        var_count = sum(1 for q in query_list if q in adata.var_names)
        obs_count = sum(1 for q in query_list if q in adata.obs.columns)
        chosen_mode = 'var' if var_count >= obs_count else 'obs'
    else:
        chosen_mode = mode

    if chosen_mode == 'var':
        valid_cols = [g for g in query_list if g in adata.var_names]
        if not valid_cols:
            return None
        sub_adata = adata[:, valid_cols]
        expr_matrix = sub_adata.X.toarray() if hasattr(sub_adata.X, "toarray") else sub_adata.X
        expr_df = pd.DataFrame(expr_matrix, columns=valid_cols)
    elif chosen_mode == 'obs':
        valid_cols = [c for c in query_list if c in adata.obs.columns]
        if not valid_cols:
            return None
        expr_df = adata.obs[valid_cols].copy()
        for col in expr_df.columns:
            if not pd.api.types.is_numeric_dtype(expr_df[col]):
                expr_df[col] = pd.Categorical(expr_df[col]).codes

    non_constant_cols = expr_df.columns[expr_df.nunique() > 1].tolist()
    if len(non_constant_cols) < 2:
        return None

    expr_df = expr_df[non_constant_cols]
    n_items = len(non_constant_cols)

    corr_matrix, p_matrix = spearmanr(expr_df, axis=0)
    corr_matrix = np.atleast_2d(corr_matrix)
    p_matrix = np.atleast_2d(p_matrix)

    if corr_matrix.shape == (1, 1) and n_items == 2:
        r_val = corr_matrix[0, 0]
        p_val = p_matrix[0, 0]
        corr_matrix = np.array([[1.0, r_val], [r_val, 1.0]])
        p_matrix = np.array([[0.0, p_val], [p_val, 0.0]])

    mask_off_diag = ~np.eye(n_items, dtype=bool)
    raw_p_flat = p_matrix[mask_off_diag]
    fdr_flat = false_discovery_control(raw_p_flat, method='bh')

    fdr_matrix = np.zeros_like(p_matrix)
    fdr_matrix[mask_off_diag] = fdr_flat

    annot_matrix = np.empty((n_items, n_items), dtype=object)
    for i in range(n_items):
        for j in range(n_items):
            if i == j:
                annot_matrix[i, j] = "1.0"
            else:
                r = corr_matrix[i, j]
                q_val = fdr_matrix[i, j]
                if q_val < 0.0001:
                    stars = "****"
                elif q_val < 0.001:
                    stars = "***"
                elif q_val < 0.01:
                    stars = "**"
                elif q_val < 0.05:
                    stars = "*"
                else:
                    stars = ""
                if np.isnan(r):
                    annot_matrix[i, j] = "NaN"
                else:
                    annot_matrix[i, j] = f"{r:.2f}{stars}"

    target_corr_df = pd.DataFrame(corr_matrix, index=non_constant_cols, columns=non_constant_cols)
    annot_df = pd.DataFrame(annot_matrix, index=non_constant_cols, columns=non_constant_cols)

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        target_corr_df, annot=annot_df, fmt="", cmap='bwr',
        center=0, vmin=-1, vmax=1, square=True, linewidths=.5, ax=ax
    )
    ax.set_title(title)
    fig.tight_layout()
    return fig


def run_aucell_scoring(adata, signature_dict):
    sig_df = pd.concat([
        pd.DataFrame({'target': v, 'collection': [k] * len(v), 'source': [k] * len(v)})
        for k, v in signature_dict.items() if k != 'target'
    ])

    if 'score_aucell' not in adata.obsm.keys():
        dc.mt.aucell(data=adata, net=sig_df, verbose=True)
        for c in adata.obsm["score_aucell"].columns:
            adata.obs[c] = adata.obsm["score_aucell"][c].tolist()

    return adata


def create_plotly_umap(adata, color_col, cell_type_col, condition_col, title, embedding_key='X_umap'):
    coords = adata.obsm[embedding_key]
    embed_label = embedding_key.replace('X_', '').upper() if embedding_key != 'X_umap' else 'UMAP'
    df = pd.DataFrame({
        f'{embed_label}1': coords[:, 0],
        f'{embed_label}2': coords[:, 1],
        'cell_type': adata.obs[cell_type_col].values,
        'condition': adata.obs[condition_col].values,
    })

    if color_col in adata.obs.columns:
        df['color_val'] = adata.obs[color_col].values
    elif color_col in adata.var_names:
        x = adata[:, color_col].X
        df['color_val'] = x.toarray().ravel() if hasattr(x, 'toarray') else np.asarray(x).ravel()
    else:
        df['color_val'] = 0

    is_numeric = pd.api.types.is_numeric_dtype(df['color_val'])

    fig = go.Figure()

    if is_numeric:
        fig.add_trace(go.Scattergl(
            x=df[f'{embed_label}1'], y=df[f'{embed_label}2'],
            mode='markers',
            marker=dict(
                color=df['color_val'], colorscale='GnBu', size=2,
                colorbar=dict(title=color_col), opacity=0.7
            ),
            customdata=np.column_stack([df['cell_type'], df['condition'], df['color_val']]),
            hovertemplate='Cell type: %{customdata[0]}<br>Condition: %{customdata[1]}<br>Value: %{customdata[2]:.3f}<extra></extra>',
        ))
    else:
        categories = df['color_val'].unique()
        for cat in categories:
            mask = df['color_val'] == cat
            sub = df[mask]
            fig.add_trace(go.Scattergl(
                x=sub[f'{embed_label}1'], y=sub[f'{embed_label}2'],
                mode='markers', marker=dict(size=2, opacity=0.7),
                name=str(cat),
                customdata=np.column_stack([sub['cell_type'], sub['condition']]),
                hovertemplate=f'{color_col}: {cat}<br>Cell type: %{{customdata[0]}}<br>Condition: %{{customdata[1]}}<extra></extra>',
            ))

    fig.update_layout(
        title=title, xaxis_title=f'{embed_label}1', yaxis_title=f'{embed_label}2',
        width=800, height=600, template='plotly_white'
    )
    return fig


def create_plotly_dotplot(adata, genes, groupby, title, mean_only_expressed=True):
    """Create interactive Plotly dotplot matching scanpy format with Reds colorscale."""
    obs_col = adata.obs[groupby]
    groups = list(obs_col.cat.categories) if hasattr(obs_col, "cat") else sorted(obs_col.unique())

    rows = []
    for group in groups:
        mask = adata.obs[groupby] == group
        subset = adata[mask]
        n_cells = int(mask.sum())
        for gene in genes:
            raw = subset[:, gene].X
            expr = raw.toarray().flatten() if hasattr(raw, "toarray") else np.array(raw).flatten()
            n_expressing = int((expr > 0).sum())
            pct_expressed = (n_expressing / n_cells * 100) if n_cells > 0 else 0.0
            if mean_only_expressed and n_expressing > 0:
                mean_expr = float(expr[expr > 0].mean())
            elif n_cells > 0:
                mean_expr = float(expr.mean())
            else:
                mean_expr = 0.0
            rows.append({
                "group": str(group), "gene": gene,
                "mean_expr": mean_expr, "pct_expressed": pct_expressed,
                "n_cells": n_cells, "n_expressing": n_expressing,
            })

    df = pd.DataFrame(rows)
    max_pct = max(df["pct_expressed"].max(), 1.0)
    dot_sizes = (4 + (df["pct_expressed"].clip(lower=0) / max_pct) * 14).tolist()

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["gene"].tolist(), y=df["group"].tolist(),
        mode="markers",
        marker=dict(
            size=dot_sizes,
            sizemode="diameter",
            color=df["mean_expr"].tolist(),
            colorscale="Reds",
            cmin=0,
            showscale=True,
            colorbar=dict(title="Mean Expression<br>(expressing cells)"),
            line=dict(width=0.5, color="black"),
        ),
        customdata=[[r["mean_expr"], r["pct_expressed"], r["n_cells"], r["n_expressing"]] for _, r in df.iterrows()],
        hovertemplate=(
            "<b>%{y}</b> | <b>%{x}</b><br>"
            "Mean Expression: %{customdata[0]:.3f}<br>"
            "Pct Expressed: %{customdata[1]:.1f}%<br>"
            "Total Cells: %{customdata[2]:.0f}<br>"
            "Expressing Cells: %{customdata[3]:.0f}<extra></extra>"
        ),
    ))
    fig.update_layout(
        title=title,
        xaxis_title="Gene",
        yaxis_title=groupby,
        yaxis=dict(categoryorder="array", categoryarray=[str(g) for g in reversed(groups)]),
        xaxis=dict(categoryorder="array", categoryarray=genes),
        template="plotly_white",
        height=max(400, len(groups) * 40 + 150),
        width=max(600, len(genes) * 80 + 200),
    )
    return fig


def create_plotly_score_dotplot(adata, score_names, groupby, title):
    """Dotplot for AUCell signature scores. X=signature names, Y=groups.
    Size = % cells with score > 0, Color = mean score (in active cells)."""
    obs_col = adata.obs[groupby]
    groups = list(obs_col.cat.categories) if hasattr(obs_col, "cat") else sorted(obs_col.unique())

    rows = []
    for group in groups:
        mask = adata.obs[groupby] == group
        n_cells = int(mask.sum())
        for score_name in score_names:
            scores = adata.obs.loc[mask, score_name].values
            n_active = int((scores > 0).sum())
            pct_active = (n_active / n_cells * 100) if n_cells > 0 else 0.0
            mean_score = float(scores[scores > 0].mean()) if n_active > 0 else 0.0
            rows.append({
                "group": str(group), "score": score_name,
                "mean_score": mean_score, "pct_active": pct_active,
                "n_cells": n_cells, "n_active": n_active,
            })

    df = pd.DataFrame(rows)
    max_pct = max(df["pct_active"].max(), 1.0)
    dot_sizes = (4 + (df["pct_active"].clip(lower=0) / max_pct) * 14).tolist()

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["score"].tolist(), y=df["group"].tolist(),
        mode="markers",
        marker=dict(
            size=dot_sizes, sizemode="diameter",
            color=df["mean_score"].tolist(), colorscale="Reds", cmin=0,
            showscale=True,
            colorbar=dict(title="Mean AUCell Score<br>(active cells)"),
            line=dict(width=0.5, color="black"),
        ),
        customdata=[[r["mean_score"], r["pct_active"], r["n_cells"], r["n_active"]] for _, r in df.iterrows()],
        hovertemplate=(
            "<b>%{y}</b> | <b>%{x}</b><br>"
            "Mean Score: %{customdata[0]:.4f}<br>"
            "Pct Active: %{customdata[1]:.1f}%<br>"
            "Total Cells: %{customdata[2]:.0f}<br>"
            "Active Cells: %{customdata[3]:.0f}<extra></extra>"
        ),
    ))
    fig.update_layout(
        title=title, xaxis_title="Signature", yaxis_title=groupby,
        yaxis=dict(categoryorder="array", categoryarray=[str(g) for g in reversed(groups)]),
        xaxis=dict(categoryorder="array", categoryarray=score_names),
        template="plotly_white",
        height=max(400, len(groups) * 40 + 150),
        width=max(600, len(score_names) * 120 + 200),
    )
    return fig


def _score_dotplot_static(adata, score_names, groupby, title):
    """Matplotlib dotplot for AUCell scores (for PDF/static)."""
    obs_col = adata.obs[groupby]
    groups = list(obs_col.cat.categories) if hasattr(obs_col, "cat") else sorted(obs_col.unique())

    rows = []
    for group in groups:
        mask = adata.obs[groupby] == group
        n_cells = int(mask.sum())
        for score_name in score_names:
            scores = adata.obs.loc[mask, score_name].values
            n_active = int((scores > 0).sum())
            pct_active = (n_active / n_cells * 100) if n_cells > 0 else 0.0
            mean_score = float(scores[scores > 0].mean()) if n_active > 0 else 0.0
            rows.append({"group": group, "score": score_name,
                         "mean_score": mean_score, "pct_active": pct_active})

    df = pd.DataFrame(rows)
    if df.empty:
        return None

    fig, ax = plt.subplots(figsize=(max(4, len(score_names) * 1.5 + 1), max(3, len(groups) * 0.4 + 1)))
    max_pct = max(df["pct_active"].max(), 1.0)
    sizes = (df["pct_active"] / max_pct * 200).clip(lower=5)

    x_map = {s: i for i, s in enumerate(score_names)}
    y_map = {g: i for i, g in enumerate(reversed(groups))}
    xs = df["score"].map(x_map)
    ys = df["group"].map(y_map)

    sc_plot = ax.scatter(xs, ys, s=sizes, c=df["mean_score"], cmap="Reds", edgecolors="black", linewidths=0.5)
    ax.set_xticks(range(len(score_names)))
    ax.set_xticklabels(score_names, rotation=45, ha='right')
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(list(reversed(groups)))
    ax.set_title(title)
    plt.colorbar(sc_plot, ax=ax, label="Mean Score")
    plt.tight_layout()
    return fig


def save_plotly_html_with_dropdown(figures_dict, output_path, title):
    if not figures_dict:
        return

    labels = list(figures_dict.keys())
    figs = list(figures_dict.values())

    combined_fig = go.Figure()
    buttons = []

    trace_counts = []
    for fig in figs:
        trace_counts.append(len(fig.data))
        for trace in fig.data:
            trace.visible = False
            combined_fig.add_trace(trace)

    total_traces = sum(trace_counts)
    offset = 0
    for i, (label, count) in enumerate(zip(labels, trace_counts)):
        visibility = [False] * total_traces
        for j in range(offset, offset + count):
            visibility[j] = True
        buttons.append(dict(
            label=label,
            method='update',
            args=[{'visible': visibility}]
        ))
        offset += count

    for j in range(trace_counts[0]):
        combined_fig.data[j].visible = True

    combined_fig.update_layout(
        title=title,
        updatemenus=[dict(
            type='dropdown', direction='down',
            x=1.0, xanchor='right', y=1.15, yanchor='top',
            buttons=buttons, active=0,
            showactive=True,
        )],
        width=900, height=650, template='plotly_white',
    )

    combined_fig.write_html(output_path, include_plotlyjs='cdn')


# ============================================================
# PER-STUDY ANALYSIS PIPELINE
# ============================================================

def process_study(study_name, study_config, signature_dict, output_dir, temp_dir):
    print(f"\n{'='*60}")
    print(f"PROCESSING: {study_name}")
    print(f"{'='*60}")

    cell_type_col = study_config['cell_type_col']
    condition_cols_raw = study_config['condition_col']
    condition_cols = [condition_cols_raw] if isinstance(condition_cols_raw, str) else list(condition_cols_raw)
    sample_col = study_config['sample_col']

    # 1. Download & Load
    adata = download_and_load_adata(study_config['s3_path'], study_name, temp_dir)

    # Create comb_condition if specified
    if 'comb_condition_cols' in study_config:
        cols = study_config['comb_condition_cols']
        adata.obs['comb_condition'] = (adata.obs[cols[0]].astype(str) + '_' + adata.obs[cols[1]].astype(str)).astype('category')

    # 2. Resolve embedding
    embedding_key = resolve_embedding_key(adata, study_config, study_name)

    # 3. Create output subdir
    study_out = os.path.join(output_dir, study_name)
    os.makedirs(study_out, exist_ok=True)

    # 4. AUCell scoring
    skip_aucell = study_config.get('skip_aucell', False)
    if not skip_aucell:
        print(f"  [{study_name}] Running AUCell scoring...")
        adata = run_aucell_scoring(adata, signature_dict)
    else:
        print(f"  [{study_name}] Skipping AUCell scoring (skip_aucell=True)")

    # Derived lists
    all_targets = signature_dict['target']
    sig_names = [k for k in signature_dict if k != 'target']
    target_genes_in_data = [g for g in all_targets if g in adata.var_names]
    sig_scores_in_data = [s for s in sig_names if s in adata.obs.columns] if not skip_aucell else []

    cell_types = adata.obs[cell_type_col].cat.categories.tolist() if hasattr(adata.obs[cell_type_col], 'cat') else sorted(adata.obs[cell_type_col].unique())

    # Collection dicts for interactive HTMLs
    umap_figures = {}
    heatmap_figures = {}
    all_pdf_plots = []
    dotplot_figures = []

    # 5. Generate all plots
    print(f"  [{study_name}] Generating plots...")

    pdf_path = os.path.join(study_out, f'{study_name}_report.pdf')
    with PdfPages(pdf_path) as pdf:

        # === UMAP Section ===
        if embedding_key:
            print(f"  [{study_name}] Embedding plots (using {embedding_key})...")
            embed_label = embedding_key.replace('X_', '').upper() if embedding_key != 'X_umap' else 'UMAP'

            for color_col in condition_cols + [cell_type_col]:
                fig_umap = create_plotly_umap(adata, color_col, cell_type_col, condition_cols[0],
                                              f'{study_name} - {embed_label} colored by {color_col}',
                                              embedding_key=embedding_key)
                umap_figures[color_col] = fig_umap

                with plt.rc_context({"figure.dpi": 300}):
                    sc.pl.embedding(adata, basis=embedding_key, color=color_col, show=False, return_fig=False)
                    fig_static = plt.gcf()
                    fig_static.suptitle(f'{study_name} - {color_col}')
                    pdf.savefig(fig_static, dpi=300, bbox_inches='tight')
                    all_pdf_plots.append((f'{embed_label} - {color_col}', _fig_to_base64(fig_static)))
                    plt.close(fig_static)

            for sig_name in sig_scores_in_data:
                fig_umap = create_plotly_umap(adata, sig_name, cell_type_col, condition_cols[0],
                                              f'{study_name} - {embed_label}: {sig_name} AUCell score',
                                              embedding_key=embedding_key)
                umap_figures[sig_name] = fig_umap

                with plt.rc_context({"figure.dpi": 300}):
                    sc.pl.embedding(adata, basis=embedding_key, color=sig_name, cmap='GnBu', show=False, return_fig=False)
                    fig_static = plt.gcf()
                    fig_static.suptitle(f'{study_name} - {sig_name}')
                    pdf.savefig(fig_static, dpi=300, bbox_inches='tight')
                    all_pdf_plots.append((f'{embed_label} - {sig_name}', _fig_to_base64(fig_static)))
                    plt.close(fig_static)

            for gene in target_genes_in_data:
                fig_umap = create_plotly_umap(adata, gene, cell_type_col, condition_cols[0],
                                              f'{study_name} - {embed_label}: {gene} expression',
                                              embedding_key=embedding_key)
                umap_figures[gene] = fig_umap

                with plt.rc_context({"figure.dpi": 300}):
                    sc.pl.embedding(adata, basis=embedding_key, color=gene, cmap='GnBu', show=False, return_fig=False)
                    fig_static = plt.gcf()
                    fig_static.suptitle(f'{study_name} - {gene}')
                    pdf.savefig(fig_static, dpi=300, bbox_inches='tight')
                    all_pdf_plots.append((f'{embed_label} - {gene}', _fig_to_base64(fig_static)))
                    plt.close(fig_static)
        else:
            print(f"  [{study_name}] Skipping embedding plots (no embedding available)")

        # === Dotplot Section ===
        print(f"  [{study_name}] Dotplots...")

        for sig_name, sig_genes in signature_dict.items():
            valid_genes = [g for g in sig_genes if g in adata.var_names]
            if not valid_genes:
                continue

            for condition_col in condition_cols:
                try:
                    sc.pl.dotplot(adata, valid_genes, groupby=condition_col,
                                  title=f'{sig_name} across {condition_col}',
                                  mean_only_expressed=True, show=False, return_fig=False)
                    fig_static = plt.gcf()
                    pdf.savefig(fig_static, dpi=300, bbox_inches='tight')
                    all_pdf_plots.append((f'Dotplot - {sig_name} by {condition_col}', _fig_to_base64(fig_static)))
                    plt.close(fig_static)
                    fig_plotly = create_plotly_dotplot(adata, valid_genes, groupby=condition_col,
                                                       title=f'{sig_name} across {condition_col}')
                    dotplot_figures.append((f'Dotplot - {sig_name} by {condition_col}', fig_plotly))
                except Exception:
                    plt.close('all')

            try:
                sc.pl.dotplot(adata, valid_genes, groupby=cell_type_col,
                              title=f'{sig_name} across {cell_type_col}',
                              mean_only_expressed=True, show=False, return_fig=False)
                fig_static = plt.gcf()
                pdf.savefig(fig_static, dpi=300, bbox_inches='tight')
                all_pdf_plots.append((f'Dotplot - {sig_name} by {cell_type_col}', _fig_to_base64(fig_static)))
                plt.close(fig_static)
                fig_plotly = create_plotly_dotplot(adata, valid_genes, groupby=cell_type_col,
                                                   title=f'{sig_name} across {cell_type_col}')
                dotplot_figures.append((f'Dotplot - {sig_name} by {cell_type_col}', fig_plotly))
            except Exception:
                plt.close('all')

        # === Signature Score Dotplots (Global) ===
        if sig_scores_in_data:
            for condition_col in condition_cols:
                try:
                    fig_plotly = create_plotly_score_dotplot(adata, sig_scores_in_data, groupby=condition_col,
                                                             title=f'Signature Scores by {condition_col}')
                    dotplot_figures.append((f'Score Dotplot - Signatures by {condition_col}', fig_plotly))
                    fig_static = _score_dotplot_static(adata, sig_scores_in_data, groupby=condition_col,
                                                        title=f'Signature Scores by {condition_col}')
                    if fig_static:
                        pdf.savefig(fig_static, dpi=300, bbox_inches='tight')
                        all_pdf_plots.append((f'Score Dotplot - Signatures by {condition_col}', _fig_to_base64(fig_static)))
                        plt.close(fig_static)
                except Exception:
                    plt.close('all')

            try:
                fig_plotly = create_plotly_score_dotplot(adata, sig_scores_in_data, groupby=cell_type_col,
                                                         title=f'Signature Scores by {cell_type_col}')
                dotplot_figures.append((f'Score Dotplot - Signatures by {cell_type_col}', fig_plotly))
                fig_static = _score_dotplot_static(adata, sig_scores_in_data, groupby=cell_type_col,
                                                    title=f'Signature Scores by {cell_type_col}')
                if fig_static:
                    pdf.savefig(fig_static, dpi=300, bbox_inches='tight')
                    all_pdf_plots.append((f'Score Dotplot - Signatures by {cell_type_col}', _fig_to_base64(fig_static)))
                    plt.close(fig_static)
            except Exception:
                plt.close('all')

        # === Stacked Violin Section (AUCell scores) ===
        print(f"  [{study_name}] Stacked violins...")

        if sig_scores_in_data:
            try:
                sc.pl.stacked_violin(adata, sig_scores_in_data, groupby=cell_type_col,
                                     title='AUCell scores across clusters', show=False, return_fig=False)
                fig_static = plt.gcf()
                pdf.savefig(fig_static, dpi=300, bbox_inches='tight')
                all_pdf_plots.append(('Stacked Violin - AUCell by cluster', _fig_to_base64(fig_static)))
                plt.close(fig_static)
            except Exception:
                plt.close('all')

            for condition_col in condition_cols:
                try:
                    sc.pl.stacked_violin(adata, sig_scores_in_data, groupby=condition_col,
                                         title=f'AUCell scores across {condition_col}', show=False, return_fig=False)
                    fig_static = plt.gcf()
                    pdf.savefig(fig_static, dpi=300, bbox_inches='tight')
                    all_pdf_plots.append((f'Stacked Violin - AUCell by {condition_col}', _fig_to_base64(fig_static)))
                    plt.close(fig_static)
                except Exception:
                    plt.close('all')

        # === Boxplot Section (Global) ===
        print(f"  [{study_name}] Global boxplots...")

        # === Correlation Heatmap Section (Global) ===
        print(f"  [{study_name}] Correlation heatmaps...")

        if len(target_genes_in_data) >= 2:
            fig_corr = plot_corr_heatmap(adata, target_genes_in_data,
                                          title=f'{study_name} - Target Gene Correlation (all cells)')
            if fig_corr:
                pdf.savefig(fig_corr, dpi=300, bbox_inches='tight')
                all_pdf_plots.append(('Correlation - Target genes (global)', _fig_to_base64(fig_corr)))
                heatmap_figures['Target genes (global)'] = fig_corr
                plt.close(fig_corr)
        else:
            print(f"  [{study_name}] Skipping target gene correlation: only {len(target_genes_in_data)} gene(s) in data")

        if len(sig_scores_in_data) >= 2:
            fig_corr = plot_corr_heatmap(adata, sig_scores_in_data, mode='obs',
                                          title=f'{study_name} - Signature Correlation (all cells)')
            if fig_corr:
                pdf.savefig(fig_corr, dpi=300, bbox_inches='tight')
                all_pdf_plots.append(('Correlation - Signature scores (global)', _fig_to_base64(fig_corr)))
                heatmap_figures['Signature scores (global)'] = fig_corr
                plt.close(fig_corr)
        else:
            print(f"  [{study_name}] Skipping signature correlation: only {len(sig_scores_in_data)} score(s)")

        # === Population Composition ===
        if HAS_XZSC:
            print(f"  [{study_name}] Population composition...")
            for condition_col in condition_cols:
                try:
                    pc = ci.ro.sc.PopulationComposition(adata, sample_col, condition_col, cell_type_col)
                    pc.boxplot(figsize=(18, 6))
                    fig_pc = plt.gcf()
                    fig_pc.suptitle(f'{study_name} - Population Composition ({condition_col})')
                    pdf.savefig(fig_pc, dpi=300, bbox_inches='tight')
                    all_pdf_plots.append((f'PopulationComposition ({condition_col})', _fig_to_base64(fig_pc)))
                    plt.close(fig_pc)

                    fractions = pc.fractions.copy()
                    if condition_col in fractions.columns and cell_type_col in fractions.columns:
                        fig_popcomp = go.Figure()
                        for ct in fractions[cell_type_col].unique():
                            ct_data = fractions[fractions[cell_type_col] == ct]
                            for cond in ct_data[condition_col].unique():
                                cond_data = ct_data[ct_data[condition_col] == cond]
                                fig_popcomp.add_trace(go.Box(
                                    y=cond_data['fraction'] if 'fraction' in cond_data.columns else cond_data.iloc[:, -1],
                                    name=f'{ct} | {cond}', boxpoints='all',
                                    hovertemplate=f'Cell type: {ct}<br>Condition: {cond}<br>Fraction: %{{y:.4f}}<extra></extra>',
                                ))
                        fig_popcomp.update_layout(title=f'{study_name} - Population Composition ({condition_col})',
                                                  width=1200, height=600, showlegend=False)
                except Exception as e:
                    print(f"  [{study_name}] PopComp ({condition_col}) failed: {e}")
                    plt.close('all')
        else:
            print(f"  [{study_name}] Skipping population composition (xzsc not installed)")

        # === Per-Cell-Type Section ===
        print(f"  [{study_name}] Per-cell-type analysis...")

        for ct in cell_types:
            ct_safe = str(ct).replace('/', '_').replace(' ', '_')
            ct_adata = adata[adata.obs[cell_type_col] == ct]

            if ct_adata.n_obs < 10:
                continue

            for sig_name, sig_genes in signature_dict.items():
                valid_genes = [g for g in sig_genes if g in adata.var_names]
                if not valid_genes:
                    continue
                for condition_col in condition_cols:
                    try:
                        sc.pl.dotplot(ct_adata, valid_genes, groupby=condition_col,
                                      title=f'{sig_name} in {ct} by {condition_col}',
                                      mean_only_expressed=True, show=False, return_fig=False)
                        fig_static = plt.gcf()
                        pdf.savefig(fig_static, dpi=300, bbox_inches='tight')
                        all_pdf_plots.append((f'Dotplot - {sig_name} in {ct} by {condition_col}', _fig_to_base64(fig_static)))
                        plt.close(fig_static)
                        fig_plotly = create_plotly_dotplot(ct_adata, valid_genes, groupby=condition_col,
                                                           title=f'{sig_name} in {ct} by {condition_col}')
                        dotplot_figures.append((f'Dotplot - {sig_name} in {ct} by {condition_col}', fig_plotly))
                    except Exception:
                        plt.close('all')

            # Signature score dotplot for this cell type
            if sig_scores_in_data:
                for condition_col in condition_cols:
                    try:
                        fig_plotly = create_plotly_score_dotplot(ct_adata, sig_scores_in_data, groupby=condition_col,
                                                                 title=f'Signature Scores in {ct} by {condition_col}')
                        dotplot_figures.append((f'Score Dotplot - Signatures in {ct} by {condition_col}', fig_plotly))
                        fig_static = _score_dotplot_static(ct_adata, sig_scores_in_data, groupby=condition_col,
                                                            title=f'Signature Scores in {ct} by {condition_col}')
                        if fig_static:
                            pdf.savefig(fig_static, dpi=300, bbox_inches='tight')
                            all_pdf_plots.append((f'Score Dotplot - Signatures in {ct} by {condition_col}', _fig_to_base64(fig_static)))
                            plt.close(fig_static)
                    except Exception:
                        plt.close('all')

            if len(target_genes_in_data) >= 2:
                try:
                    fig_corr = plot_corr_heatmap(ct_adata, target_genes_in_data,
                                                  title=f'Target Correlation in {ct}')
                    if fig_corr:
                        pdf.savefig(fig_corr, dpi=300, bbox_inches='tight')
                        all_pdf_plots.append((f'Correlation - Target genes ({ct})', _fig_to_base64(fig_corr)))
                        heatmap_figures[f'Target genes ({ct})'] = fig_corr
                        plt.close(fig_corr)
                except Exception:
                    plt.close('all')
            else:
                print(f"  [{study_name}] Skipping per-cell-type target correlation in {ct}: only {len(target_genes_in_data)} gene(s)")

            if len(sig_scores_in_data) >= 2:
                try:
                    fig_corr = plot_corr_heatmap(ct_adata, sig_scores_in_data, mode='obs',
                                                  title=f'Signature Correlation in {ct}')
                    if fig_corr:
                        pdf.savefig(fig_corr, dpi=300, bbox_inches='tight')
                        all_pdf_plots.append((f'Correlation - Signature scores ({ct})', _fig_to_base64(fig_corr)))
                        heatmap_figures[f'Signature scores ({ct})'] = fig_corr
                        plt.close(fig_corr)
                except Exception:
                    plt.close('all')
            else:
                print(f"  [{study_name}] Skipping per-cell-type signature correlation in {ct}: only {len(sig_scores_in_data)} score(s)")


    # 6. Save interactive HTML files
    print(f"  [{study_name}] Saving interactive HTML files...")

    save_plotly_html_with_dropdown(
        umap_figures,
        os.path.join(study_out, 'umap_interactive.html'),
        f'{study_name} - Interactive Embeddings'
    )

    if heatmap_figures:
        _save_heatmap_html(heatmap_figures, os.path.join(study_out, 'heatmap_interactive.html'),
                           f'{study_name} - Correlation Heatmaps')

    _save_all_plots_html(all_pdf_plots, os.path.join(study_out, 'all_plots_interactive.html'),
                         f'{study_name} - All Plots')

    if dotplot_figures:
        _save_dotplot_html(dotplot_figures, os.path.join(study_out, 'dotplot_interactive.html'),
                           f'{study_name} - Interactive Dotplots')

    print(f"  [{study_name}] COMPLETE. Output: {study_out}")

    del adata
    import gc
    gc.collect()


def _fig_to_base64(fig, dpi=300):
    """Convert a matplotlib figure to a base64-encoded PNG string."""
    import base64
    from io import BytesIO
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight')
    buf.seek(0)
    img_b64 = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    return img_b64


def _save_all_plots_html(plots_list, output_path, title):
    """Save all collected (name, content) tuples as HTML with searchable dropdown.
    Content can be a base64 PNG string or a plotly.graph_objects.Figure.
    Uses Plotly.react() with a single visible div for interactive figures."""
    if not plots_list:
        return

    has_plotly = any(isinstance(content, go.Figure) for _, content in plots_list)
    plotly_js_tag = '<script src="https://cdn.plot.ly/plotly-3.3.0.min.js"></script>' if has_plotly else ""

    labels = [name for name, _ in plots_list]
    fig_json_lines = []
    img_data_lines = []
    plot_types = []

    for name, content in plots_list:
        safe_name = name.replace("\\", "\\\\").replace('"', '\\"')
        if isinstance(content, go.Figure):
            fig_json_lines.append(f'figData["{safe_name}"] = {content.to_json()};')
            plot_types.append("plotly")
        else:
            img_data_lines.append(f'imgData["{safe_name}"] = "{content}";')
            plot_types.append("image")

    first_label = labels[0].replace("\\", "\\\\").replace('"', '\\"')

    html = f"""<!DOCTYPE html>
<html><head><title>{title}</title>
{plotly_js_tag}
<style>
body {{ font-family: Arial, sans-serif; margin: 20px; }}
#search-box {{ padding: 8px; width: 500px; font-size: 14px; border: 1px solid #ccc; border-radius: 4px; }}
#dropdown {{ max-height: 400px; overflow-y: auto; border: 1px solid #ccc; width: 500px; display: none; position: absolute; background: white; z-index: 1000; }}
#dropdown div {{ padding: 6px 10px; cursor: pointer; }}
#dropdown div:hover {{ background: #e0e0e0; }}
#dropdown div.highlight {{ background: #b3d4fc; }}
.controls {{ margin-bottom: 15px; position: relative; }}
#img-container img {{ max-width: 100%; border: 1px solid #ddd; }}
#plot-title {{ font-size: 16px; font-weight: bold; margin: 10px 0; color: #333; }}
</style>
</head><body>
<h2>{title}</h2>
<div class="controls">
<input id="search-box" type="text" placeholder="Type to search plots..." onfocus="showDropdown()" oninput="filterDropdown()">
<div id="dropdown"></div>
</div>
<div id="plot-title"></div>
<div id="plot-container" style="width:100%;"></div>
<div id="img-container"></div>
<script>
var figData = {{}};
var imgData = {{}};
{chr(10).join(fig_json_lines)}
{chr(10).join(img_data_lines)}
var labels = {labels};
var plotTypes = {plot_types};
var currentIdx = 0;
var highlightIdx = -1;
var filteredIndices = [];

function showDropdown() {{
    document.getElementById('dropdown').style.display = 'block';
    filterDropdown();
}}

function hideDropdown() {{
    setTimeout(function() {{ document.getElementById('dropdown').style.display = 'none'; highlightIdx = -1; }}, 200);
}}

function filterDropdown() {{
    var input = document.getElementById('search-box').value.toLowerCase();
    var dropdown = document.getElementById('dropdown');
    dropdown.innerHTML = '';
    filteredIndices = [];
    for (var i = 0; i < labels.length; i++) {{
        if (labels[i].toLowerCase().includes(input)) {{
            filteredIndices.push(i);
            (function(idx, fi) {{
                var div = document.createElement('div');
                div.textContent = labels[idx];
                div.onclick = function() {{ selectPlot(idx); }};
                dropdown.appendChild(div);
            }})(i, filteredIndices.length - 1);
        }}
    }}
    highlightIdx = -1;
    dropdown.style.display = filteredIndices.length > 0 ? 'block' : 'none';
    if (filteredIndices.length === 1 && labels[filteredIndices[0]].toLowerCase() === input) {{
        selectPlot(filteredIndices[0]);
    }}
}}

function handleKeydown(e) {{
    var dropdown = document.getElementById('dropdown');
    if (dropdown.style.display === 'none') return;
    var items = dropdown.children;
    if (e.key === 'ArrowDown') {{
        e.preventDefault();
        highlightIdx = Math.min(highlightIdx + 1, items.length - 1);
        updateHighlight(items);
    }} else if (e.key === 'ArrowUp') {{
        e.preventDefault();
        highlightIdx = Math.max(highlightIdx - 1, 0);
        updateHighlight(items);
    }} else if (e.key === 'Enter') {{
        e.preventDefault();
        if (highlightIdx >= 0 && highlightIdx < filteredIndices.length) {{
            selectPlot(filteredIndices[highlightIdx]);
        }}
    }} else if (e.key === 'Escape') {{
        dropdown.style.display = 'none';
        highlightIdx = -1;
    }}
}}

function updateHighlight(items) {{
    for (var i = 0; i < items.length; i++) {{
        items[i].classList.toggle('highlight', i === highlightIdx);
    }}
    if (highlightIdx >= 0 && items[highlightIdx]) {{
        items[highlightIdx].scrollIntoView({{block: 'nearest'}});
    }}
}}

function selectPlot(idx) {{
    var label = labels[idx];
    document.getElementById('search-box').value = label;
    document.getElementById('plot-title').textContent = label;
    document.getElementById('dropdown').style.display = 'none';
    currentIdx = idx;
    highlightIdx = -1;
    var plotDiv = document.getElementById('plot-container');
    var imgDiv = document.getElementById('img-container');
    if (plotTypes[idx] === 'plotly') {{
        imgDiv.style.display = 'none';
        plotDiv.style.display = 'block';
        var fig = figData[label];
        Plotly.react(plotDiv, fig.data, fig.layout, {{responsive: true}});
    }} else {{
        plotDiv.style.display = 'none';
        imgDiv.style.display = 'block';
        imgDiv.innerHTML = '<img src="data:image/png;base64,' + imgData[label] + '">';
    }}
}}

document.getElementById('search-box').onblur = hideDropdown;
document.getElementById('search-box').addEventListener('keydown', handleKeydown);
selectPlot(0);
</script></body></html>"""

    with open(output_path, 'w') as f:
        f.write(html)


def _save_dotplot_html(plots_list, output_path, title):
    """Save interactive Plotly dotplots as HTML with searchable dropdown and cell-type filter."""
    if not plots_list:
        return

    labels = [name for name, _ in plots_list]
    fig_json_lines = []
    for name, fig in plots_list:
        safe_name = name.replace("\\", "\\\\").replace('"', '\\"')
        fig_json_lines.append(f'figData["{safe_name}"] = {fig.to_json()};')

    cell_types = sorted(set(
        name.split(" in ")[-1] for name, _ in plots_list if " in " in name
    ))

    cell_type_options = "\n".join(
        f'<option value="{ct}">{ct}</option>' for ct in cell_types
    )

    first_label = labels[0].replace("\\", "\\\\").replace('"', '\\"')

    html = f"""<!DOCTYPE html>
<html><head><title>{title}</title>
<script src="https://cdn.plot.ly/plotly-3.3.0.min.js"></script>
<style>
body {{ font-family: Arial, sans-serif; margin: 20px; }}
#search-box {{ padding: 8px; width: 500px; font-size: 14px; border: 1px solid #ccc; border-radius: 4px; }}
#dropdown {{ max-height: 400px; overflow-y: auto; border: 1px solid #ccc; width: 500px; display: none; position: absolute; background: white; z-index: 1000; }}
#dropdown div {{ padding: 6px 10px; cursor: pointer; }}
#dropdown div:hover {{ background: #e0e0e0; }}
#dropdown div.highlight {{ background: #b3d4fc; }}
.controls {{ margin-bottom: 15px; position: relative; display: flex; gap: 15px; align-items: flex-start; }}
#ct-filter {{ min-width: 200px; max-width: 300px; height: 120px; font-size: 13px; }}
#plot-title {{ font-size: 16px; font-weight: bold; margin: 10px 0; color: #333; }}
.filter-group {{ display: flex; flex-direction: column; gap: 4px; }}
.filter-group label {{ font-size: 12px; color: #666; font-weight: bold; }}
.filter-group small {{ font-size: 11px; color: #999; }}
</style>
</head><body>
<h2>{title}</h2>
<div class="controls">
  <div class="filter-group">
    <label>Search plots:</label>
    <input id="search-box" type="text" placeholder="Type to search dotplots..." onfocus="showDropdown()" oninput="filterDropdown()">
    <div id="dropdown"></div>
  </div>
  <div class="filter-group">
    <label>Filter by cell type:</label>
    <select id="ct-filter" multiple onchange="filterDropdown()">
      <option value="" selected>Show All</option>
      {cell_type_options}
    </select>
    <small>Ctrl+click to multi-select</small>
  </div>
</div>
<div id="plot-title"></div>
<div id="plot-container" style="width:100%;"></div>
<script>
var figData = {{}};
{chr(10).join(fig_json_lines)}
var labels = {labels};
var currentIdx = 0;
var highlightIdx = -1;
var filteredIndices = [];

function getSelectedCellTypes() {{
    var sel = document.getElementById('ct-filter');
    var selected = [];
    for (var i = 0; i < sel.options.length; i++) {{
        if (sel.options[i].selected && sel.options[i].value !== '') {{
            selected.push(sel.options[i].value);
        }}
    }}
    return selected;
}}

function showDropdown() {{
    document.getElementById('dropdown').style.display = 'block';
    filterDropdown();
}}

function hideDropdown() {{
    setTimeout(function() {{ document.getElementById('dropdown').style.display = 'none'; highlightIdx = -1; }}, 200);
}}

function filterDropdown() {{
    var input = document.getElementById('search-box').value.toLowerCase();
    var ctFilter = getSelectedCellTypes();
    var dropdown = document.getElementById('dropdown');
    dropdown.innerHTML = '';
    filteredIndices = [];
    for (var i = 0; i < labels.length; i++) {{
        var matchesText = labels[i].toLowerCase().includes(input);
        var matchesCt = ctFilter.length === 0 || ctFilter.some(function(ct) {{
            return labels[i].includes(ct);
        }});
        if (matchesText && matchesCt) {{
            filteredIndices.push(i);
            (function(idx) {{
                var div = document.createElement('div');
                div.textContent = labels[idx];
                div.onclick = function() {{ selectPlot(idx); }};
                dropdown.appendChild(div);
            }})(i);
        }}
    }}
    highlightIdx = -1;
    dropdown.style.display = filteredIndices.length > 0 ? 'block' : 'none';
    if (filteredIndices.length === 1 && labels[filteredIndices[0]].toLowerCase() === input) {{
        selectPlot(filteredIndices[0]);
    }}
}}

function handleKeydown(e) {{
    var dropdown = document.getElementById('dropdown');
    if (dropdown.style.display === 'none') return;
    var items = dropdown.children;
    if (e.key === 'ArrowDown') {{
        e.preventDefault();
        highlightIdx = Math.min(highlightIdx + 1, items.length - 1);
        updateHighlight(items);
    }} else if (e.key === 'ArrowUp') {{
        e.preventDefault();
        highlightIdx = Math.max(highlightIdx - 1, 0);
        updateHighlight(items);
    }} else if (e.key === 'Enter') {{
        e.preventDefault();
        if (highlightIdx >= 0 && highlightIdx < filteredIndices.length) {{
            selectPlot(filteredIndices[highlightIdx]);
        }}
    }} else if (e.key === 'Escape') {{
        dropdown.style.display = 'none';
        highlightIdx = -1;
    }}
}}

function updateHighlight(items) {{
    for (var i = 0; i < items.length; i++) {{
        items[i].classList.toggle('highlight', i === highlightIdx);
    }}
    if (highlightIdx >= 0 && items[highlightIdx]) {{
        items[highlightIdx].scrollIntoView({{block: 'nearest'}});
    }}
}}

function selectPlot(idx) {{
    var label = labels[idx];
    document.getElementById('search-box').value = label;
    document.getElementById('plot-title').textContent = label;
    document.getElementById('dropdown').style.display = 'none';
    currentIdx = idx;
    highlightIdx = -1;
    var plotDiv = document.getElementById('plot-container');
    var fig = figData[label];
    Plotly.react(plotDiv, fig.data, fig.layout, {{responsive: true}});
}}

document.getElementById('search-box').onblur = hideDropdown;
document.getElementById('search-box').addEventListener('keydown', handleKeydown);
selectPlot(0);
</script></body></html>"""

    with open(output_path, 'w') as f:
        f.write(html)


def _save_heatmap_html(figures_dict, output_path, title):
    """Save matplotlib heatmap figures as base64-embedded HTML with dropdown."""
    import base64
    from io import BytesIO

    html_parts = [f"""<!DOCTYPE html>
<html><head><title>{title}</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 20px; }}
select {{ padding: 8px; font-size: 14px; margin-bottom: 15px; }}
img {{ max-width: 100%; border: 1px solid #ddd; }}
.plot-container {{ display: none; }}
.plot-container.active {{ display: block; }}
</style>
</head><body>
<h2>{title}</h2>
<select id="heatmap-select" onchange="switchHeatmap()">
"""]

    labels = list(figures_dict.keys())
    for i, label in enumerate(labels):
        selected = ' selected' if i == 0 else ''
        html_parts.append(f'<option value="plot-{i}"{selected}>{label}</option>\n')

    html_parts.append('</select>\n')

    for i, (label, fig) in enumerate(figures_dict.items()):
        buf = BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        active = ' active' if i == 0 else ''
        html_parts.append(f'<div class="plot-container{active}" id="plot-{i}"><img src="data:image/png;base64,{img_base64}"></div>\n')
        buf.close()

    html_parts.append("""
<script>
function switchHeatmap() {
    var sel = document.getElementById('heatmap-select').value;
    var plots = document.querySelectorAll('.plot-container');
    plots.forEach(function(p) { p.classList.remove('active'); });
    document.getElementById(sel).classList.add('active');
}
</script></body></html>""")

    with open(output_path, 'w') as f:
        f.write(''.join(html_parts))


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    args = parse_args()

    signature_dict = build_signature_dict(args.targets, args.signatures, args.addon_signatures)
    output_dir = args.output_dir
    temp_dir = args.temp_dir
    max_workers = args.max_workers

    # Filter studies
    if args.sc_studies:
        study_names = [s.strip() for s in args.sc_studies.split(',')]
        studies = {k: v for k, v in STUDY_CONFIGS.items() if k in study_names}
        skipped = [k for k in study_names if k not in STUDY_CONFIGS]
        if skipped:
            print(f"WARNING: Unknown study names (skipped): {skipped}")
            print(f"  Available: {list(STUDY_CONFIGS.keys())}")
    else:
        studies = STUDY_CONFIGS

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)

    print("=" * 60)
    print("ANALYSIS CONFIGURATION")
    print("=" * 60)
    print(f"Target genes: {signature_dict['target']}")
    print(f"Signatures to score: {[k for k in signature_dict if k != 'target']}")
    print(f"Studies: {list(studies.keys())}")
    print(f"Max parallel workers: {max_workers}")
    print(f"Output directory: {output_dir}")
    print("=" * 60)

    print(f"Running {len(studies)} studies: {list(studies.keys())}")

    if max_workers == 1:
        for name, cfg in studies.items():
            try:
                process_study(name, cfg, signature_dict, output_dir, temp_dir)
                print(f"[DONE] {name}")
            except Exception as e:
                print(f"[FAIL] {name}: {e}")
                traceback.print_exc()
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(process_study, name, cfg, signature_dict, output_dir, temp_dir): name
                for name, cfg in studies.items()
            }
            for future in as_completed(futures):
                study = futures[future]
                try:
                    future.result()
                    print(f"[DONE] {study}")
                except Exception as e:
                    print(f"[FAIL] {study}: {e}")
                    traceback.print_exc()

    print("\n" + "=" * 60)
    print("ALL STUDIES COMPLETE")
    print(f"Results in: {os.path.abspath(output_dir)}")
    print("=" * 60)


if __name__ == '__main__':
    main()
