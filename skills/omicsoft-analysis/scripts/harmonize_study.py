"""
Harmonize internal study data into S3 OmicSoft h5ad format.

Usage:
    conda run -n <env> python harmonize_study.py \
        --input <deg_results_dir_or_file> \
        --project-id <short_name> \
        --study-name <descriptive_name> \
        --organism human \
        --platform RNA-seq \
        --de-method DESeq2 \
        --output-dir <output_path> \
        [--gene-annotation <gene_annotation.txt>] \
        [--metadata <sample_metadata.csv>] \
        [--expression-matrix <normalized_counts.tsv>] \
        [--raw-counts <raw_counts.tsv>] \
        [--vocab <s3_vocabulary.json>] \
        [--validate-only]

This script:
1. Auto-detects gene identifier format
2. Maps gene IDs to HGNC symbols using annotation files
3. Proposes category mappings from s3_vocabulary.json
4. Generates h5ad with correct S3 schema
5. Runs concat-ready validation checks
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SKILL_DIR = Path(__file__).parent.parent
VOCAB_PATH = SKILL_DIR / "references" / "s3_vocabulary.json"

DEG_OBS_COLUMNS = [
    'database', 'comparison_id', 'comparison_contrast', 'comparison_category',
    'case_tissue', 'case_sample_material', 'case_disease_state', 'case_disease_subtype',
    'case_disease_group', 'case_disease_location', 'case_response',
    'control_tissue', 'control_sample_material', 'control_disease_state',
    'control_disease_subtype', 'control_disease_group', 'control_disease_location',
    'control_response', 'case_dosage', 'case_treatment', 'case_treatment_group',
    'case_treatment_status', 'case_treat_time', 'control_dosage', 'control_treatment',
    'control_treatment_group', 'control_treatment_status', 'control_treat_time',
    'case_age_category', 'case_gender', 'case_ethnicity',
    'control_age_category', 'control_gender', 'control_ethnicity',
    'case_sample_ids', 'control_sample_ids',
    'project_id', 'study', 'tissue', 'sample', 'disease', 'comparison'
]

EXPR_OBS_COLUMNS = [
    'database', 'tissue', 'project_id', 'disease_state', 'disease_stage',
    'ethnicity', 'gender', 'age_summary', 'treatment', 'sampling_time',
    'response', 'subject_id', 'cell_type', 'symptom', 'infection',
    'transfection', 'sample_integration_id', 'sample_pathology', 'sample_source',
    'sample_type', 'collection', 'data_source', 'organism', 'description',
    'title', 'platform_name', 'experiment_type', 'project_description',
    'project_title', 'comparison_group'
]


def load_vocabulary(vocab_path=None):
    """Load S3 controlled vocabulary."""
    path = Path(vocab_path) if vocab_path else VOCAB_PATH
    if not path.exists():
        print(f"WARNING: Vocabulary file not found at {path}")
        return {}
    with open(path) as f:
        return json.load(f)


def detect_gene_id_type(ids):
    """Auto-detect gene identifier type from a list of IDs."""
    sample = [str(x) for x in ids[:100]]

    ensembl_pattern = re.compile(r'^ENS[A-Z]*G\d+(\.\d+)?$')
    seqid_pattern = re.compile(r'^\d+-\d+$')
    probe_pattern = re.compile(r'^\d+_at$|^\d+_s_at$|^ILMN_\d+$')

    ensembl_count = sum(1 for s in sample if ensembl_pattern.match(s))
    seqid_count = sum(1 for s in sample if seqid_pattern.match(s))
    probe_count = sum(1 for s in sample if probe_pattern.match(s))
    all_upper = sum(1 for s in sample if s.isalpha() and s.isupper())

    total = len(sample)
    if ensembl_count / total > 0.8:
        has_version = any('.' in s for s in sample if ensembl_pattern.match(s))
        return 'ensembl_versioned' if has_version else 'ensembl'
    elif seqid_count / total > 0.8:
        return 'seqid'
    elif probe_count / total > 0.5:
        return 'probe'
    elif all_upper / total > 0.7:
        return 'hgnc_symbol'
    else:
        return 'unknown'


def strip_ensembl_version(gene_id):
    """Strip version suffix from Ensembl ID."""
    return gene_id.split('.')[0]


def map_genes(gene_ids, annotation_file, id_type):
    """Map gene IDs to HGNC symbols using annotation file."""
    if annotation_file is None or not Path(annotation_file).exists():
        if id_type == 'hgnc_symbol':
            return {g: g for g in gene_ids}, []
        print(f"ERROR: Gene annotation file required for id_type={id_type}")
        return {}, list(gene_ids)

    annot = pd.read_csv(annotation_file, sep='\t')

    mapped = {}
    unmapped = []

    if id_type in ('ensembl', 'ensembl_versioned'):
        stripped = {strip_ensembl_version(g): g for g in gene_ids}
        if 'gene_id' in annot.columns and 'gene_name' in annot.columns:
            annot_map = dict(zip(annot['gene_id'], annot['gene_name']))
        elif 'GeneID' in annot.columns and 'GeneName' in annot.columns:
            annot_map = dict(zip(annot['GeneID'], annot['GeneName']))
        else:
            print(f"WARNING: Cannot find gene_id/gene_name columns in annotation. Available: {annot.columns.tolist()}")
            return {g: g for g in gene_ids}, []

        for stripped_id, original_id in stripped.items():
            if stripped_id in annot_map and pd.notna(annot_map[stripped_id]):
                mapped[original_id] = annot_map[stripped_id]
            else:
                unmapped.append(original_id)

    elif id_type == 'seqid':
        if 'SeqId' in annot.columns and 'EntrezGeneSymbol' in annot.columns:
            annot_map = dict(zip(annot['SeqId'].astype(str), annot['EntrezGeneSymbol']))
        else:
            print(f"WARNING: Cannot find SeqId/EntrezGeneSymbol columns. Available: {annot.columns.tolist()}")
            return {g: g for g in gene_ids}, []

        for gid in gene_ids:
            if str(gid) in annot_map and pd.notna(annot_map[str(gid)]):
                mapped[gid] = annot_map[str(gid)]
            else:
                unmapped.append(gid)

    elif id_type == 'hgnc_symbol':
        mapped = {g: g for g in gene_ids}

    else:
        print(f"WARNING: Unsupported id_type={id_type}, using IDs as-is")
        mapped = {g: g for g in gene_ids}

    return mapped, unmapped


def find_best_vocab_match(value, vocab_list, threshold=0.6):
    """Find best matching term in vocabulary list using substring matching."""
    if not value or not vocab_list:
        return None

    value_lower = value.lower().strip()

    for term in vocab_list:
        if term.lower() == value_lower:
            return term

    for term in vocab_list:
        if value_lower in term.lower() or term.lower() in value_lower:
            return term

    return None


def propose_mappings(local_values, vocab_field, vocabulary):
    """Propose mappings from local values to S3 vocabulary."""
    if vocab_field not in vocabulary:
        return {v: v for v in local_values}

    vocab_list = vocabulary[vocab_field]
    proposals = {}

    for val in local_values:
        if pd.isna(val) or val in ('', 'NA', 'None', None):
            proposals[val] = 'NA'
            continue
        match = find_best_vocab_match(str(val), vocab_list)
        proposals[val] = match if match else str(val)

    return proposals


def compute_sig_score(log2fc_mat, padj_mat, pval_mat=None, cutoff_padj=0.05, cutoff_pval=None):
    """Compute sig_score with zero p-value floor."""
    f32_tiny = np.finfo(np.float32).tiny

    padj_work = padj_mat.copy()
    nonzero_padj = padj_work[padj_work > 0]
    if len(nonzero_padj) > 0:
        padj_floor = max(nonzero_padj.min() / 10.0, f32_tiny)
        padj_work[padj_work == 0] = padj_floor

    if pval_mat is not None:
        pval_work = pval_mat.copy()
        nonzero_pval = pval_work[pval_work > 0]
        if len(nonzero_pval) > 0:
            pval_floor = max(nonzero_pval.min() / 10.0, f32_tiny)
            pval_work[pval_work == 0] = pval_floor

    if cutoff_pval is not None and pval_mat is not None:
        is_significant = pval_work < cutoff_pval
    else:
        is_significant = padj_work < cutoff_padj

    sig_score = np.where(
        is_significant,
        log2fc_mat * -np.log10(padj_work),
        0.0
    )

    return sig_score.astype(np.float32)


def validate_concat_ready(adata, schema_type='deg'):
    """Validate h5ad is ready for concat with S3 OmicSoft data."""
    errors = []

    if adata.obs.index.name != 'sample_id':
        errors.append(f"obs.index.name = '{adata.obs.index.name}', expected 'sample_id'")
    if adata.var.index.name != 'gene_id':
        errors.append(f"var.index.name = '{adata.var.index.name}', expected 'gene_id'")

    if schema_type == 'deg':
        required_cols = DEG_OBS_COLUMNS
        required_layers = ['pval', 'padj', 'sig_score']
    else:
        required_cols = EXPR_OBS_COLUMNS
        required_layers = ['raw_counts']

    missing_cols = [c for c in required_cols if c not in adata.obs.columns]
    if missing_cols:
        errors.append(f"Missing obs columns: {missing_cols}")

    missing_layers = [l for l in required_layers if l not in adata.layers]
    if missing_layers:
        errors.append(f"Missing layers: {missing_layers}")

    if not adata.obs.index.is_unique:
        dups = adata.obs.index[adata.obs.index.duplicated()].tolist()
        errors.append(f"Duplicate obs indices: {dups[:5]}")

    numeric_vars = [v for v in adata.var_names[:100] if str(v).isdigit()]
    if numeric_vars:
        errors.append(f"Numeric var_names found: {numeric_vars[:5]}")

    if 'schema_version' not in adata.uns:
        errors.append("Missing uns['schema_version']")
    if 'schema_type' not in adata.uns:
        errors.append("Missing uns['schema_type']")

    if errors:
        print("VALIDATION FAILED:")
        for e in errors:
            print(f"  ✗ {e}")
        return False
    else:
        print("✓ Validation passed - concat-ready!")
        return True


def build_deg_h5ad(deg_results, gene_map, obs_metadata, project_id, study_name,
                   cutoff_padj=0.05, cutoff_pval=None, cutoff_log2fc=0.0):
    """Build DEG h5ad from harmonized data."""
    import anndata as ad
    from scipy.sparse import csr_matrix

    genes = sorted(set(gene_map.values()))
    comparisons = list(deg_results.keys())

    n_comp = len(comparisons)
    n_genes = len(genes)
    gene_idx = {g: i for i, g in enumerate(genes)}

    log2fc_mat = np.zeros((n_comp, n_genes), dtype=np.float32)
    pval_mat = np.ones((n_comp, n_genes), dtype=np.float32)
    padj_mat = np.ones((n_comp, n_genes), dtype=np.float32)

    for i, comp_name in enumerate(comparisons):
        df = deg_results[comp_name]
        for _, row in df.iterrows():
            raw_gene = row.get('gene_id', row.name)
            symbol = gene_map.get(raw_gene, raw_gene)
            if symbol in gene_idx:
                j = gene_idx[symbol]
                log2fc_mat[i, j] = row.get('log2FoldChange', row.get('logFC', 0.0))
                pval_mat[i, j] = row.get('pvalue', row.get('P.Value', 1.0))
                padj_mat[i, j] = row.get('padj', row.get('adj.P.Val', 1.0))

    sig_score_mat = compute_sig_score(log2fc_mat, padj_mat, pval_mat, cutoff_padj, cutoff_pval)

    obs_df = pd.DataFrame(index=comparisons)
    obs_df.index.name = 'sample_id'

    for col in DEG_OBS_COLUMNS:
        if col in obs_metadata.columns:
            obs_df[col] = obs_metadata[col].values if len(obs_metadata) == n_comp else 'NA'
        else:
            obs_df[col] = 'NA'

    obs_df['database'] = 'internal'
    obs_df['project_id'] = project_id
    obs_df['study'] = study_name
    obs_df['cutoff_padj'] = cutoff_padj
    obs_df['cutoff_pval'] = cutoff_pval if cutoff_pval else ''
    obs_df['cutoff_log2fc'] = cutoff_log2fc
    obs_df['source'] = 'internal'

    var_df = pd.DataFrame(index=genes)
    var_df.index.name = 'gene_id'

    adata = ad.AnnData(
        X=csr_matrix(log2fc_mat),
        obs=obs_df,
        var=var_df,
        layers={
            'pval': csr_matrix(pval_mat),
            'padj': csr_matrix(padj_mat),
            'sig_score': csr_matrix(sig_score_mat),
        }
    )

    adata.uns['schema_version'] = '1.0'
    adata.uns['schema_type'] = 'deg'
    adata.uns['available_layers'] = ['log2fc', 'pval', 'padj', 'sig_score']

    return adata


def build_expr_h5ad(expr_matrix, gene_map, obs_metadata, project_id,
                    raw_counts=None, platform_name='RNA-seq'):
    """Build EXPR h5ad from harmonized data."""
    import anndata as ad
    from scipy.sparse import csr_matrix

    mapped_cols = {}
    for col in expr_matrix.columns:
        symbol = gene_map.get(col, col)
        mapped_cols[col] = symbol
    expr_matrix = expr_matrix.rename(columns=mapped_cols)

    obs_df = obs_metadata.copy()
    obs_df.index.name = 'sample_id'

    for col in EXPR_OBS_COLUMNS:
        if col not in obs_df.columns:
            obs_df[col] = 'NA'

    obs_df['database'] = 'internal'
    obs_df['project_id'] = project_id
    obs_df['organism'] = 'Homo sapiens'
    obs_df['platform_name'] = platform_name
    obs_df['source'] = 'internal'

    var_df = pd.DataFrame(index=expr_matrix.columns)
    var_df.index.name = 'gene_id'

    X = csr_matrix(expr_matrix.values.astype(np.float32))

    layers = {}
    if raw_counts is not None:
        raw_mapped = raw_counts.rename(columns=mapped_cols)
        raw_aligned = raw_mapped.reindex(columns=expr_matrix.columns, fill_value=0)
        layers['raw_counts'] = csr_matrix(raw_aligned.values.astype(np.float32))
        raw_source = 'raw'
    else:
        layers['raw_counts'] = X.copy()
        raw_source = 'normalized_expression'

    adata = ad.AnnData(X=X, obs=obs_df, var=var_df, layers=layers)
    adata.uns['schema_version'] = '1.0'
    adata.uns['schema_type'] = 'expr'
    adata.uns['available_layers'] = ['normalized_expression', 'raw_counts']
    adata.uns['raw_counts_source'] = raw_source

    return adata


def main():
    parser = argparse.ArgumentParser(description='Harmonize internal study into S3 OmicSoft h5ad format')
    parser.add_argument('--input', required=True, help='Path to DEG results (file or directory)')
    parser.add_argument('--project-id', required=True, help='Short project identifier')
    parser.add_argument('--study-name', required=True, help='Descriptive study name')
    parser.add_argument('--organism', default='human', choices=['human', 'mouse'])
    parser.add_argument('--platform', default='RNA-seq', help='Platform name')
    parser.add_argument('--de-method', default='DESeq2', help='DE method (DESeq2, limma, edgeR)')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    parser.add_argument('--gene-annotation', help='Gene annotation file for ID mapping')
    parser.add_argument('--metadata', help='Sample metadata CSV/TSV')
    parser.add_argument('--expression-matrix', help='Normalized expression matrix')
    parser.add_argument('--raw-counts', help='Raw counts matrix')
    parser.add_argument('--vocab', help='Path to s3_vocabulary.json')
    parser.add_argument('--cutoff-padj', type=float, default=0.05)
    parser.add_argument('--cutoff-pval', type=float, default=None)
    parser.add_argument('--cutoff-log2fc', type=float, default=0.0)
    parser.add_argument('--validate-only', action='store_true', help='Only validate existing h5ad')

    args = parser.parse_args()

    vocab = load_vocabulary(args.vocab)

    if args.validate_only:
        import anndata as ad
        adata = ad.read_h5ad(args.input)
        schema_type = 'deg' if 'comparison_category' in adata.obs.columns else 'expr'
        validate_concat_ready(adata, schema_type)
        return

    print(f"=== Harmonizing study: {args.study_name} ===")
    print(f"  Project ID: {args.project_id}")
    print(f"  Organism: {args.organism}")
    print(f"  Platform: {args.platform}")
    print(f"  DE method: {args.de_method}")
    print()

    input_path = Path(args.input)
    if input_path.is_dir():
        deg_files = list(input_path.glob('*.tsv')) + list(input_path.glob('*.csv'))
        print(f"  Found {len(deg_files)} DEG result files")
    else:
        deg_files = [input_path]

    if not deg_files:
        print("ERROR: No DEG result files found")
        sys.exit(1)

    first_df = pd.read_csv(deg_files[0], sep='\t' if str(deg_files[0]).endswith('.tsv') else ',')
    gene_col = first_df.columns[0] if first_df.columns[0] not in ('log2FoldChange', 'logFC', 'pvalue') else first_df.index.name
    gene_ids = first_df.iloc[:, 0].tolist() if gene_col == first_df.columns[0] else first_df.index.tolist()

    id_type = detect_gene_id_type(gene_ids)
    print(f"  Detected gene ID type: {id_type}")

    gene_map, unmapped = map_genes(gene_ids, args.gene_annotation, id_type)
    print(f"  Mapped: {len(gene_map)} genes, Unmapped: {len(unmapped)} genes")

    gene_map_out = {
        'mapped': gene_map,
        'unmapped': unmapped,
        'id_type': id_type,
        'annotation_file': args.gene_annotation
    }
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / 'gene_map.json', 'w') as f:
        json.dump(gene_map_out, f, indent=2)
    print(f"  Saved gene_map.json")

    if vocab:
        print("\n  === Proposed Vocabulary Mappings ===")
        if args.metadata:
            meta_df = pd.read_csv(args.metadata, sep='\t' if args.metadata.endswith('.tsv') else ',')
            for col_name, vocab_field in [('tissue', 'tissue'), ('disease', 'disease_state'),
                                           ('treatment', 'treatment'), ('gender', 'gender')]:
                if col_name in meta_df.columns:
                    unique_vals = meta_df[col_name].dropna().unique()
                    proposals = propose_mappings(unique_vals, vocab_field, vocab)
                    print(f"\n  {col_name} -> {vocab_field}:")
                    for local, s3 in proposals.items():
                        marker = " (UNMAPPED)" if local == s3 and s3 not in vocab.get(vocab_field, []) else ""
                        print(f"    '{local}' -> '{s3}'{marker}")

    print(f"\n  === Loading DEG results ===")
    deg_results = {}
    for f in deg_files:
        comp_name = f.stem
        sep = '\t' if str(f).endswith('.tsv') else ','
        deg_results[comp_name] = pd.read_csv(f, sep=sep)
    print(f"  Loaded {len(deg_results)} comparisons")

    print(f"\n  === Building DEG h5ad ===")
    obs_meta = pd.DataFrame(index=list(deg_results.keys()))
    deg_adata = build_deg_h5ad(
        deg_results, gene_map, obs_meta,
        project_id=args.project_id,
        study_name=args.study_name,
        cutoff_padj=args.cutoff_padj,
        cutoff_pval=args.cutoff_pval,
        cutoff_log2fc=args.cutoff_log2fc
    )

    deg_path = out_dir / f"{args.project_id.lower()}_deg.h5ad"
    deg_adata.write_h5ad(deg_path)
    print(f"  Written: {deg_path} ({deg_adata.shape[0]} comparisons x {deg_adata.shape[1]} genes)")
    validate_concat_ready(deg_adata, 'deg')

    if args.expression_matrix:
        print(f"\n  === Building EXPR h5ad ===")
        sep = '\t' if args.expression_matrix.endswith('.tsv') else ','
        expr_df = pd.read_csv(args.expression_matrix, sep=sep, index_col=0)

        raw_df = None
        if args.raw_counts:
            raw_sep = '\t' if args.raw_counts.endswith('.tsv') else ','
            raw_df = pd.read_csv(args.raw_counts, sep=raw_sep, index_col=0)

        meta_df = pd.DataFrame(index=expr_df.index) if args.metadata is None else pd.read_csv(
            args.metadata, sep='\t' if args.metadata.endswith('.tsv') else ',', index_col=0
        )

        expr_adata = build_expr_h5ad(
            expr_df.T, gene_map, meta_df,
            project_id=args.project_id,
            raw_counts=raw_df.T if raw_df is not None else None,
            platform_name=args.platform
        )

        expr_path = out_dir / f"{args.project_id.lower()}_expr.h5ad"
        expr_adata.write_h5ad(expr_path)
        print(f"  Written: {expr_path} ({expr_adata.shape[0]} samples x {expr_adata.shape[1]} genes)")
        validate_concat_ready(expr_adata, 'expr')

    print("\n=== Harmonization complete ===")
    print(f"Output directory: {out_dir}")


if __name__ == '__main__':
    main()
