#!/usr/bin/env python3
"""Download scGPT pretrained model from HuggingFace or Google Drive.

Usage:
    python download_model.py --output-dir /path/to/destination
    python download_model.py --output-dir ./models --source gdrive --model brain
"""

import argparse
import os
import sys

HUGGINGFACE_REPO = "MohamedMabrouk/scGPT"

GDRIVE_URLS = {
    "whole-human": "https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y",
    "continual": "https://drive.google.com/drive/folders/1_GROJTzXiAV8HB4imruOTk6PEGuNOcgB",
    "brain": "https://drive.google.com/drive/folders/1vf1ijfQSk7rGdDGpBntR5bi5g6gNt-Gx",
    "blood": "https://drive.google.com/drive/folders/1kkug5C7NjvXIwQGGaGoqXTk_Lb_pDrBU",
    "heart": "https://drive.google.com/drive/folders/1GcgXrd7apn6y4Ze_iSCncskX3UsWPY2r",
    "lung": "https://drive.google.com/drive/folders/16A1DJ30PT6bodt4bWLa4hpS7gbWZQFBG",
    "kidney": "https://drive.google.com/drive/folders/1S-1AR65DF120kNFpEbWCvRHPhpkGK3kK",
    "pan-cancer": "https://drive.google.com/drive/folders/13QzLHilYUd0v3HTwa_9n4G4yEF-hdkqa",
}


def download_huggingface(output_dir):
    """Download whole-human model from HuggingFace community mirror."""
    from huggingface_hub import snapshot_download

    print(f"Downloading from HuggingFace: {HUGGINGFACE_REPO}")
    path = snapshot_download(
        HUGGINGFACE_REPO,
        local_dir=output_dir,
        ignore_patterns=["cxg_faiss_index/*"],
    )
    print(f"Downloaded to: {path}")
    return path


def download_gdrive(output_dir, model_name):
    """Download model from Google Drive using gdown."""
    try:
        import gdown
    except ImportError:
        print("Installing gdown...")
        os.system(f"{sys.executable} -m pip install gdown")
        import gdown

    url = GDRIVE_URLS.get(model_name)
    if not url:
        print(f"Unknown model: {model_name}")
        print(f"Available: {', '.join(GDRIVE_URLS.keys())}")
        sys.exit(1)

    print(f"Downloading {model_name} from Google Drive...")
    gdown.download_folder(url, output=output_dir, quiet=False)
    print(f"Downloaded to: {output_dir}")
    return output_dir


def verify_model(model_dir):
    """Verify model files exist and are loadable."""
    import torch
    from scgpt.tokenizer import GeneVocab

    required = ["best_model.pt", "vocab.json"]
    for f in required:
        path = os.path.join(model_dir, f)
        if os.path.exists(path):
            size = os.path.getsize(path) / 1e6
            print(f"  [OK] {f} ({size:.1f} MB)")
        else:
            print(f"  [MISSING] {f}")
            return False

    # Test loading
    vocab = GeneVocab.from_file(os.path.join(model_dir, "vocab.json"))
    print(f"  Vocab: {len(vocab)} genes")

    ckpt = torch.load(os.path.join(model_dir, "best_model.pt"), map_location="cpu")
    print(f"  Checkpoint: {len(ckpt)} parameters")

    args_path = os.path.join(model_dir, "args.json")
    if os.path.exists(args_path):
        import json
        with open(args_path) as f:
            args = json.load(f)
        print(f"  Args: {args}")

    return True


def main():
    parser = argparse.ArgumentParser(description="Download scGPT pretrained model")
    parser.add_argument("--output-dir", required=True, help="Destination directory")
    parser.add_argument("--source", choices=["huggingface", "gdrive"], default="huggingface",
                        help="Download source (default: huggingface)")
    parser.add_argument("--model", default="whole-human",
                        help="Model name for gdrive (default: whole-human)")
    parser.add_argument("--verify", action="store_true", default=True,
                        help="Verify after download")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.source == "huggingface":
        download_huggingface(args.output_dir)
    else:
        download_gdrive(args.output_dir, args.model)

    if args.verify:
        print("\nVerifying model files:")
        if verify_model(args.output_dir):
            print(f"\nModel ready at: {args.output_dir}")
        else:
            print("\nWARNING: Some files missing.")
            sys.exit(1)


if __name__ == "__main__":
    main()
