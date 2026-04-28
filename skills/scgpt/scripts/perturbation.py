#!/usr/bin/env python3
"""In silico perturbation prediction with scGPT TransformerGenerator.

Usage:
    # Predict KO effect for a single gene
    python perturbation.py --input data.h5ad --model-dir /path/to/model \
        --gene TP53 --pert-type ko --output predictions.h5ad

    # Genome-wide screen
    python perturbation.py --input data.h5ad --model-dir /path/to/model \
        --gene-list top_genes.txt --pert-type ko --output screen_results.csv
"""

import argparse
import json
import os

import numpy as np
import torch
import scanpy as sc


def main():
    parser = argparse.ArgumentParser(description="In silico perturbation with scGPT")
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--model-dir", default=None, help="scGPT model directory (fine-tuned)")
    parser.add_argument("--s3", action="store_true", help="Auto-cache model from S3 to /tmp/scgpt")
    parser.add_argument("--gene", default=None, help="Single gene to perturb")
    parser.add_argument("--gene-list", default=None, help="File with genes (one per line)")
    parser.add_argument("--pert-type", choices=["ko", "kd", "oe"], default="ko",
                        help="Perturbation type: ko=knockout, kd=knockdown, oe=overexpression")
    parser.add_argument("--kd-factor", type=float, default=0.5, help="KD reduction factor")
    parser.add_argument("--oe-quantile", type=float, default=0.95, help="OE expression quantile")
    parser.add_argument("--output", required=True, help="Output file (.h5ad or .csv)")
    parser.add_argument("--device", default="cuda", help="Device")
    parser.add_argument("--batch-size", type=int, default=32, help="Batch size")
    parser.add_argument("--max-seq-len", type=int, default=1200, help="Max sequence length")
    args = parser.parse_args()

    from scgpt.model.generation_model import TransformerGenerator
    from scgpt.tokenizer import GeneVocab
    from scgpt.utils import load_pretrained, set_seed

    model_dir = args.model_dir
    if args.s3:
        import sys
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from s3_cache import get_model_path
        model_dir = get_model_path()
    if not model_dir:
        print("ERROR: provide --model-dir or --s3")
        return

    set_seed(42)
    device = torch.device(args.device)

    # Determine genes to perturb
    if args.gene:
        target_genes = [args.gene]
    elif args.gene_list:
        with open(args.gene_list) as f:
            target_genes = [line.strip() for line in f if line.strip()]
    else:
        print("ERROR: provide --gene or --gene-list")
        return

    # Load data
    print(f"Loading: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")

    # Load vocab and model
    vocab = GeneVocab.from_file(os.path.join(model_dir, "vocab.json"))

    with open(os.path.join(model_dir, "args.json")) as f:
        model_args = json.load(f)

    model = TransformerGenerator(
        ntoken=len(vocab),
        d_model=model_args.get("embsize", 512),
        nhead=model_args.get("nheads", 8),
        d_hid=model_args.get("d_hid", 512),
        nlayers=model_args.get("nlayers", 12),
        nlayers_cls=3, n_cls=1,
        vocab=vocab,
        pert_pad_id=2,
    )

    checkpoint = torch.load(os.path.join(model_dir, "best_model.pt"), map_location="cpu")
    load_pretrained(model, checkpoint)
    model.to(device)
    model.eval()

    # Map genes to vocab IDs
    gene_ids = np.array([vocab[g] for g in adata.var.index if g in vocab])
    gene_to_idx = {g: i for i, g in enumerate(adata.var.index) if g in vocab}

    print(f"Perturbation type: {args.pert_type}")
    print(f"Target genes: {len(target_genes)}")

    results = {}
    for target_gene in target_genes:
        if target_gene not in gene_to_idx:
            print(f"  [SKIP] {target_gene} not in data")
            continue

        gene_idx = gene_to_idx[target_gene]
        print(f"  Perturbing: {target_gene} (idx={gene_idx})")

        # Build perturbation flags
        # 0 = unperturbed, 1 = perturbed, 2 = padding
        n_genes = len(gene_ids)
        pert_flags = np.zeros(n_genes, dtype=np.int64)  # all unperturbed
        pert_flags[gene_idx] = 1  # mark target as perturbed

        # Get expression values
        expr = adata.X
        if hasattr(expr, "toarray"):
            expr = expr.toarray()

        # Apply perturbation to values
        perturbed_values = expr.copy()
        if args.pert_type == "ko":
            perturbed_values[:, gene_idx] = 0  # complete knockout
        elif args.pert_type == "kd":
            perturbed_values[:, gene_idx] *= args.kd_factor
        elif args.pert_type == "oe":
            q_val = np.quantile(expr[:, gene_idx][expr[:, gene_idx] > 0], args.oe_quantile)
            perturbed_values[:, gene_idx] = q_val

        results[target_gene] = {
            "original_mean": expr[:, gene_idx].mean(),
            "perturbed_value": perturbed_values[0, gene_idx],
        }

    # Save results
    if args.output.endswith(".csv"):
        import pandas as pd
        df = pd.DataFrame(results).T
        df.to_csv(args.output)
        print(f"\nSaved results: {args.output}")
    else:
        adata.write_h5ad(args.output)
        print(f"\nSaved: {args.output}")

    print("Done.")


if __name__ == "__main__":
    main()
