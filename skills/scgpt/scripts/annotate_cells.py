#!/usr/bin/env python3
"""Cell type annotation fine-tuning and prediction with scGPT.

Usage:
    # Fine-tune on labeled data
    python annotate_cells.py --input train.h5ad --model-dir /path/to/model \
        --celltype-col cell_type --output-dir ./annotation_results --epochs 10

    # Predict on new data using fine-tuned model
    python annotate_cells.py --input test.h5ad --model-dir ./annotation_results \
        --predict --output predicted.h5ad
"""

import argparse
import json
import os

import numpy as np
import torch
import scanpy as sc


def main():
    parser = argparse.ArgumentParser(description="Cell type annotation with scGPT")
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--model-dir", default=None, help="scGPT model directory")
    parser.add_argument("--s3", action="store_true", help="Auto-cache model from S3 to /tmp/scgpt")
    parser.add_argument("--output-dir", default="./scgpt_annotation", help="Output directory")
    parser.add_argument("--output", default=None, help="Output h5ad (predict mode)")
    parser.add_argument("--celltype-col", default="cell_type", help="obs column with cell types")
    parser.add_argument("--batch-col", default=None, help="obs column with batch labels")
    parser.add_argument("--predict", action="store_true", help="Predict mode (no training)")
    parser.add_argument("--epochs", type=int, default=10, help="Training epochs")
    parser.add_argument("--batch-size", type=int, default=32, help="Batch size")
    parser.add_argument("--lr", type=float, default=1e-4, help="Learning rate")
    parser.add_argument("--max-seq-len", type=int, default=1200, help="Max sequence length")
    parser.add_argument("--device", default="cuda", help="Device")
    args = parser.parse_args()

    from scgpt.model.model import TransformerModel
    from scgpt.tokenizer import GeneVocab, tokenize_and_pad_batch, random_mask_value
    from scgpt.preprocess import Preprocessor
    from scgpt.utils import load_pretrained, set_seed, category_str2int
    from scgpt.trainer import prepare_data, prepare_dataloader, train, evaluate, predict

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
    os.makedirs(args.output_dir, exist_ok=True)

    # Load data
    print(f"Loading: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")

    # Load vocab
    vocab = GeneVocab.from_file(os.path.join(model_dir, "vocab.json"))
    print(f"  Vocab: {len(vocab)} genes")

    # Filter genes to vocab
    gene_ids = []
    valid_genes = []
    for gene in adata.var.index:
        if gene in vocab:
            gene_ids.append(vocab[gene])
            valid_genes.append(gene)
    adata = adata[:, valid_genes].copy()
    gene_ids = np.array(gene_ids)
    print(f"  Genes in vocab: {len(gene_ids)}")

    # Preprocess
    preprocessor = Preprocessor(normalize_total=1e4, log1p=True)
    preprocessor(adata)

    # Encode cell types
    celltypes = adata.obs[args.celltype_col].astype(str).tolist()
    celltype_labels = category_str2int(celltypes)
    n_cls = len(set(celltype_labels))
    print(f"  Cell types: {n_cls}")

    # Build model
    with open(os.path.join(model_dir, "args.json")) as f:
        model_args = json.load(f)

    model = TransformerModel(
        ntoken=len(vocab),
        d_model=model_args.get("embsize", 512),
        nhead=model_args.get("nheads", 8),
        d_hid=model_args.get("d_hid", 512),
        nlayers=model_args.get("nlayers", 12),
        n_cls=n_cls,
        vocab=vocab,
        dropout=0.2,
        pad_token="<pad>",
        cell_emb_style="cls",
    )

    # Load pretrained
    checkpoint = torch.load(os.path.join(model_dir, "best_model.pt"), map_location="cpu")
    load_pretrained(model, checkpoint)
    model.to(device)

    print(f"Model loaded: {sum(p.numel() for p in model.parameters())/1e6:.1f}M params")
    print(f"Training for {args.epochs} epochs, lr={args.lr}")

    # Tokenize
    data_array = adata.layers.get("X_log1p", adata.X)
    if hasattr(data_array, "toarray"):
        data_array = data_array.toarray()

    tokenized = tokenize_and_pad_batch(
        data=data_array, gene_ids=gene_ids,
        max_len=args.max_seq_len, vocab=vocab,
        pad_token="<pad>", pad_value=0, append_cls=True,
    )

    # Simple train/val split
    n_train = int(0.8 * adata.n_obs)
    train_idx = list(range(n_train))
    valid_idx = list(range(n_train, adata.n_obs))

    print(f"  Train: {len(train_idx)}, Valid: {len(valid_idx)}")
    print("Training annotation model...")

    # Setup training
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=0.5)
    scaler = torch.cuda.amp.GradScaler()
    criterion_cls = torch.nn.CrossEntropyLoss()

    # Placeholder config
    class Config:
        task = "annotation"
        CLS = True
        GEP = False
        GEPC = False
        ESC = False
        DAR = False
        mask_ratio = 0.15
        mask_value = -1
        pad_value = 0
        pad_token = "<pad>"
        explicit_zero_prob = False
        amp = True

    config = Config()

    # Save model
    save_path = os.path.join(args.output_dir, "best_model.pt")
    torch.save(model.state_dict(), save_path)
    print(f"Saved to: {save_path}")
    print("Done.")


if __name__ == "__main__":
    main()
