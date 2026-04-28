#!/usr/bin/env python3
"""CLI wrapper for Geneformer TranscriptomeTokenizer.

Usage:
    python tokenize_data.py --input-dir data/ --output-dir tokenized/ --output-prefix my_data
    python tokenize_data.py --input-dir data/ --output-dir tokenized/ --output-prefix my_data --file-format h5ad
    python tokenize_data.py --input-dir data/ --output-dir tokenized/ --output-prefix my_data --custom-attrs '{"cell_type":"cell_type"}'
"""

import argparse
import json
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Tokenize scRNA-seq data using Geneformer TranscriptomeTokenizer"
    )
    parser.add_argument(
        "--input-dir", required=True, help="Directory containing .loom/.h5ad/.zarr files"
    )
    parser.add_argument(
        "--output-dir", required=True, help="Output directory for tokenized .dataset"
    )
    parser.add_argument(
        "--output-prefix", required=True, help="Output file prefix"
    )
    parser.add_argument(
        "--file-format",
        choices=["loom", "h5ad", "zarr"],
        default="loom",
        help="Input file format (default: loom)",
    )
    parser.add_argument(
        "--custom-attrs",
        default=None,
        help='JSON string for custom_attr_name_dict, e.g. \'{"cell_type":"cell_type"}\'',
    )
    parser.add_argument(
        "--nproc", type=int, default=4, help="Number of processes (default: 4)"
    )
    parser.add_argument(
        "--chunk-size", type=int, default=512, help="Chunk size for h5ad/zarr (default: 512)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    custom_attr_name_dict = None
    if args.custom_attrs is not None:
        try:
            custom_attr_name_dict = json.loads(args.custom_attrs)
        except json.JSONDecodeError as exc:
            print(f"Error parsing --custom-attrs JSON: {exc}", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    from geneformer import TranscriptomeTokenizer

    tokenizer = TranscriptomeTokenizer(
        custom_attr_name_dict=custom_attr_name_dict,
        nproc=args.nproc,
        chunk_size=args.chunk_size,
    )

    print(f"Tokenizing {args.file_format} files in: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Output prefix: {args.output_prefix}")

    tokenizer.tokenize_data(
        args.input_dir,
        args.output_dir,
        args.output_prefix,
        file_format=args.file_format,
    )

    output_path = os.path.join(args.output_dir, f"{args.output_prefix}.dataset")
    if os.path.exists(output_path):
        from datasets import load_from_disk

        ds = load_from_disk(output_path)
        print(f"\nTokenization complete:")
        print(f"  Output: {output_path}")
        print(f"  Cells:  {len(ds)}")
        print(f"  Columns: {ds.column_names}")
    else:
        print(f"\nTokenization complete. Check {args.output_dir} for output.")


if __name__ == "__main__":
    main()
