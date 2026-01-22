#!/usr/bin/env python3
"""
Convert large CSV pathway network file to parquet format.

This script uses PyArrow's streaming capabilities to handle very large files efficiently.
"""

import pyarrow as pa
import pyarrow.csv as csv
import pyarrow.parquet as pq
import sys
import os
from pathlib import Path

def convert_csv_to_parquet(csv_file: str, parquet_file: str = None):
    """
    Convert CSV to parquet using PyArrow streaming.

    Args:
        csv_file: Input CSV filename
        parquet_file: Output parquet filename (optional, defaults to same name with .parquet)
    """
    if not os.path.exists(csv_file):
        print(f"Error: File {csv_file} not found")
        sys.exit(1)

    # Generate output filename if not provided
    if parquet_file is None:
        parquet_file = str(Path(csv_file).with_suffix('.parquet'))

    print(f"Converting {csv_file} to {parquet_file}")

    # Get file size
    file_size_gb = os.path.getsize(csv_file) / (1024**3)
    print(f"Input file size: {file_size_gb:.2f} GB")

    try:
        print("Reading CSV with PyArrow (streaming)...")

        # Read CSV using PyArrow's streaming reader
        # This is much more memory efficient for large files
        table = csv.read_csv(csv_file)

        print(f"Read {table.num_rows:,} rows with {table.num_columns} columns")
        print(f"Columns: {', '.join(table.column_names)}")
        print("\nWriting to parquet with Snappy compression...")

        # Write to parquet with compression
        pq.write_table(
            table,
            parquet_file,
            compression='snappy',
            # Use row group size for better compression on large files
            row_group_size=1000000
        )

        print(f"\n{'='*80}")
        print(f"Conversion complete!")
        print(f"{'='*80}")
        print(f"Total rows: {table.num_rows:,}")

        # Get output file size
        output_size_gb = os.path.getsize(parquet_file) / (1024**3)
        print(f"Input CSV size: {file_size_gb:.2f} GB")
        print(f"Output Parquet size: {output_size_gb:.2f} GB")
        print(f"Compression ratio: {file_size_gb/output_size_gb:.2f}x")
        print(f"Space saved: {file_size_gb - output_size_gb:.2f} GB ({(1 - output_size_gb/file_size_gb)*100:.1f}%)")
        print(f"\nSaved to: {parquet_file}")
        print(f"{'='*80}")

    except Exception as e:
        print(f"\nError during conversion: {e}")
        import traceback
        traceback.print_exc()
        # Clean up partial file if error occurs
        if os.path.exists(parquet_file):
            os.remove(parquet_file)
            print(f"Cleaned up partial output file")
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python convert_csv_to_parquet.py <input.csv> [output.parquet]")
        print("\nExample:")
        print("  python convert_csv_to_parquet.py pathway_network.csv")
        print("  python convert_csv_to_parquet.py pathway_network.csv output.parquet")
        print("\nThis uses PyArrow's efficient streaming reader for large files.")
        sys.exit(1)

    csv_file = sys.argv[1]
    parquet_file = sys.argv[2] if len(sys.argv) > 2 else None

    convert_csv_to_parquet(csv_file, parquet_file)
