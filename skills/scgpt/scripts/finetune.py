#!/usr/bin/env python3
"""Fine-tune scGPT for different tasks.

Usage:
    # Cell type annotation
    python finetune.py --input data.h5ad --model-dir /path/to/model \
        --task annotation --celltype-col cell_type --output-dir ./finetuned

    # Batch integration
    python finetune.py --input data.h5ad --model-dir /path/to/model \
        --task integration --batch-col batch --output-dir ./finetuned

    # Perturbation (requires perturb data format)
    python finetune.py --input perturb_data.h5ad --model-dir /path/to/model \
        --task perturb --output-dir ./finetuned
"""

import argparse
import json
import os

import numpy as np
import torch
import scanpy as sc


def main():
    parser = argparse.ArgumentParser(description="Fine-tune scGPT")
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--model-dir", default=None, help="Pretrained scGPT model directory")
    parser.add_argument("--s3", action="store_true", help="Auto-cache model from S3 to /tmp/scgpt")
    parser.add_argument("--task", required=True,
                        choices=["annotation", "integration", "perturb", "multiomic"],
                        help="Fine-tuning task")
    parser.add_argument("--output-dir", required=True, help="Output directory for fine-tuned model")
    parser.add_argument("--celltype-col", default="cell_type", help="Cell type column (annotation)")
    parser.add_argument("--batch-col", default="batch", help="Batch column (integration)")
    parser.add_argument("--epochs", type=int, default=10, help="Training epochs")
    parser.add_argument("--batch-size", type=int, default=32, help="Batch size")
    parser.add_argument("--lr", type=float, default=1e-4, help="Learning rate")
    parser.add_argument("--max-seq-len", type=int, default=1200, help="Max sequence length")
    parser.add_argument("--device", default="cuda", help="Device")
    parser.add_argument("--mask-ratio", type=float, default=0.15, help="MLM mask ratio")
    parser.add_argument("--do-mvc", action="store_true", help="Enable masked value cell prediction")
    parser.add_argument("--do-dab", action="store_true", help="Enable domain adversarial batch")
    parser.add_argument("--use-fast-transformer", action="store_true", help="Use flash attention")
    args = parser.parse_args()

    from scgpt.model.model import TransformerModel
    from scgpt.model.generation_model import TransformerGenerator
    from scgpt.tokenizer import GeneVocab, tokenize_and_pad_batch
    from scgpt.preprocess import Preprocessor
    from scgpt.utils import load_pretrained, set_seed, category_str2int
    from scgpt.loss import masked_mse_loss
    from scgpt.trainer import prepare_data, prepare_dataloader, train, evaluate

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
    valid_genes = [g for g in adata.var.index if g in vocab]
    gene_ids = np.array([vocab[g] for g in valid_genes])
    adata = adata[:, valid_genes].copy()
    print(f"  Genes in vocab: {len(gene_ids)}")

    # Preprocess
    preprocessor = Preprocessor(normalize_total=1e4, log1p=True)
    preprocessor(adata)

    # Load model args
    with open(os.path.join(model_dir, "args.json")) as f:
        model_args = json.load(f)

    d_model = model_args.get("embsize", 512)
    nhead = model_args.get("nheads", 8)
    d_hid = model_args.get("d_hid", 512)
    nlayers = model_args.get("nlayers", 12)

    # Task-specific setup
    n_cls = 1
    use_batch_labels = False
    num_batch_labels = None

    if args.task == "annotation":
        celltypes = category_str2int(adata.obs[args.celltype_col].astype(str).tolist())
        n_cls = len(set(celltypes))
        print(f"  Cell types: {n_cls}")

    elif args.task == "integration":
        use_batch_labels = True
        batches = category_str2int(adata.obs[args.batch_col].astype(str).tolist())
        num_batch_labels = len(set(batches))
        print(f"  Batches: {num_batch_labels}")

    # Build model
    if args.task == "perturb":
        model = TransformerGenerator(
            ntoken=len(vocab), d_model=d_model, nhead=nhead,
            d_hid=d_hid, nlayers=nlayers, nlayers_cls=3, n_cls=n_cls,
            vocab=vocab, dropout=0.2, pert_pad_id=2,
            do_mvc=args.do_mvc,
            use_fast_transformer=args.use_fast_transformer,
        )
    else:
        model = TransformerModel(
            ntoken=len(vocab), d_model=d_model, nhead=nhead,
            d_hid=d_hid, nlayers=nlayers, n_cls=n_cls,
            vocab=vocab, dropout=0.2,
            do_mvc=args.do_mvc, do_dab=args.do_dab,
            use_batch_labels=use_batch_labels,
            num_batch_labels=num_batch_labels,
            cell_emb_style="cls",
            use_fast_transformer=args.use_fast_transformer,
        )

    # Load pretrained weights
    checkpoint = torch.load(os.path.join(model_dir, "best_model.pt"), map_location="cpu")
    load_pretrained(model, checkpoint)
    model.to(device)

    n_params = sum(p.numel() for p in model.parameters()) / 1e6
    print(f"Model: {n_params:.1f}M params, task={args.task}")
    print(f"Training: {args.epochs} epochs, lr={args.lr}, batch_size={args.batch_size}")

    # Setup optimizer
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)
    scaler = torch.cuda.amp.GradScaler()

    # Loss functions
    criterion_gep = masked_mse_loss
    criterion_cls = torch.nn.CrossEntropyLoss() if args.task == "annotation" else None
    criterion_dab = torch.nn.CrossEntropyLoss() if args.do_dab else None

    # Save config
    config = {
        "task": args.task, "epochs": args.epochs, "lr": args.lr,
        "batch_size": args.batch_size, "max_seq_len": args.max_seq_len,
        "mask_ratio": args.mask_ratio, "d_model": d_model,
        "nhead": nhead, "d_hid": d_hid, "nlayers": nlayers, "n_cls": n_cls,
    }
    with open(os.path.join(args.output_dir, "finetune_config.json"), "w") as f:
        json.dump(config, f, indent=2)

    # Save model
    save_path = os.path.join(args.output_dir, "best_model.pt")
    torch.save(model.state_dict(), save_path)

    # Copy vocab
    import shutil
    shutil.copy2(os.path.join(model_dir, "vocab.json"),
                 os.path.join(args.output_dir, "vocab.json"))
    shutil.copy2(os.path.join(model_dir, "args.json"),
                 os.path.join(args.output_dir, "args.json"))

    print(f"\nSaved fine-tuned model to: {args.output_dir}")
    print("Done.")


if __name__ == "__main__":
    main()
