#!/usr/bin/env python3
"""
PrimeKG Data Download and Processing

Functions to download PrimeKG from Harvard Dataverse and split into train/test/valid.
Adapted from the dataset-primekg skill for use in MCP server context.
"""
import polars as pl
import os
import subprocess
import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)


def check_primekg_data(data_path: str) -> Dict:
    """
    Check if PrimeKG data files exist.

    Args:
        data_path: Base path for PrimeKG data (e.g., /path/to/data/primekg1)

    Returns:
        Dictionary with status information:
        - data_path: str - The checked path
        - exists: dict - Status of each required file
        - all_present: bool - True if all files exist
    """
    required_files = {
        "primekg.csv": os.path.join(os.path.dirname(data_path), "primekg.csv"),
        "train.txt": os.path.join(data_path, "raw", "train.txt"),
        "test.txt": os.path.join(data_path, "raw", "test.txt"),
        "valid.txt": os.path.join(data_path, "raw", "valid.txt"),
        "nodes.txt": os.path.join(data_path, "raw", "nodes.txt"),
    }

    status = {"data_path": data_path, "exists": {}, "all_present": True}

    for name, path in required_files.items():
        exists = os.path.exists(path)
        status["exists"][name] = exists
        if not exists:
            status["all_present"] = False

    return status


def download_primekg(data_path: str) -> str:
    """
    Download PrimeKG CSV from Harvard Dataverse.

    Args:
        data_path: Directory to save primekg.csv (typically parent of dataset directory)

    Returns:
        Path to downloaded primekg.csv

    Raises:
        RuntimeError: If download fails
    """
    os.makedirs(data_path, exist_ok=True)
    csv_path = os.path.join(data_path, "primekg.csv")

    # Skip if already exists
    if os.path.exists(csv_path):
        logger.info(f"PrimeKG CSV already exists at: {csv_path}")
        return csv_path

    logger.info(f"Downloading PrimeKG to: {data_path}")

    # Use wget to download
    download_cmd = [
        "wget",
        "-O", csv_path,
        "https://dataverse.harvard.edu/api/access/datafile/6180620",
        "--no-check-certificate",
        "-nc"  # no-clobber: don't overwrite existing files
    ]

    try:
        result = subprocess.run(
            download_cmd,
            check=True,
            capture_output=True,
            text=True
        )
        logger.info(f"Download successful: {csv_path}")
        return csv_path
    except subprocess.CalledProcessError as e:
        error_msg = f"Failed to download PrimeKG: {e.stderr}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)


def process_primekg(
    csv_path: str,
    output_path: str,
    train_frac: float = 0.8,
    test_frac: float = 0.1,
    valid_frac: float = 0.1,
    seed: int = 42,
) -> Dict:
    """
    Process PrimeKG CSV into train/test/valid splits.

    Args:
        csv_path: Path to primekg.csv
        output_path: Directory to save processed files (typically {dataset}/raw/)
        train_frac: Training set fraction (default: 0.8)
        test_frac: Test set fraction (default: 0.1)
        valid_frac: Validation set fraction (default: 0.1)
        seed: Random seed for reproducibility

    Returns:
        Dictionary with processing statistics:
        - nodes_count: int
        - train_edges: int
        - test_edges: int
        - valid_edges: int
        - output_path: str
    """
    logger.info(f"Loading PrimeKG from: {csv_path}")

    # Load CSV with proper schema
    df = pl.read_csv(
        csv_path,
        schema={
            "relation": pl.Categorical,
            "display_relation": pl.Categorical,
            "x_index": pl.Int64,
            "x_id": pl.String,
            "x_type": pl.Categorical,
            "x_name": pl.String,
            "x_source": pl.Categorical,
            "y_index": pl.Int64,
            "y_id": pl.String,
            "y_type": pl.Categorical,
            "y_name": pl.String,
            "y_source": pl.Categorical,
        },
    )

    logger.info(f"Loaded {df.shape[0]:,} edges")

    # Create nodes file
    logger.info("Creating nodes file...")
    nodes = (
        pl.concat(
            [
                df[["x_index", "x_id", "x_name", "x_type", "x_source"]].rename(
                    {
                        "x_index": "index",
                        "x_id": "source_id",
                        "x_name": "name",
                        "x_type": "type",
                        "x_source": "source",
                    }
                ),
                df[["y_index", "y_id", "y_name", "y_type", "y_source"]].rename(
                    {
                        "y_index": "index",
                        "y_id": "source_id",
                        "y_name": "name",
                        "y_type": "type",
                        "y_source": "source",
                    }
                ),
            ]
        )
        .with_columns(source_label=pl.col("source") + ":" + pl.col("source_id"))
        .unique()
    )

    logger.info(f"Created {nodes.shape[0]:,} unique nodes")

    # Create triples (edges)
    logger.info("Creating triples...")
    triples = df.with_columns(
        head_source_label=pl.col("x_source") + ":" + pl.col("x_id"),
        tail_source_label=pl.col("y_source") + ":" + pl.col("y_id"),
    )[["head_source_label", "display_relation", "tail_source_label"]].unique()

    # Split into train/test/valid
    logger.info(f"Splitting data: {train_frac}/{test_frac}/{valid_frac}")
    train = triples.sample(fraction=train_frac, seed=seed, with_replacement=False)
    remaining = triples.join(
        train,
        on=["head_source_label", "display_relation", "tail_source_label"],
        how="anti",
    )
    test = remaining.sample(
        fraction=test_frac / (test_frac + valid_frac), seed=seed, with_replacement=False
    )
    valid = remaining.join(
        test,
        on=["head_source_label", "display_relation", "tail_source_label"],
        how="anti",
    )

    logger.info(f"Train: {train.shape[0]:,} edges")
    logger.info(f"Test: {test.shape[0]:,} edges")
    logger.info(f"Valid: {valid.shape[0]:,} edges")

    # Export files
    os.makedirs(output_path, exist_ok=True)

    logger.info(f"Exporting to: {output_path}")
    train.write_csv(
        os.path.join(output_path, "train.txt"), separator="\t", include_header=False
    )
    test.write_csv(
        os.path.join(output_path, "test.txt"), separator="\t", include_header=False
    )
    valid.write_csv(
        os.path.join(output_path, "valid.txt"), separator="\t", include_header=False
    )
    nodes.drop("index").write_csv(
        os.path.join(output_path, "nodes.txt"), separator="\t", include_header=True
    )

    logger.info("Processing complete!")

    return {
        "nodes_count": nodes.shape[0],
        "train_edges": train.shape[0],
        "test_edges": test.shape[0],
        "valid_edges": valid.shape[0],
        "output_path": output_path,
    }


def setup_primekg(
    dataset_path: str,
    force_redownload: bool = False,
    train_frac: float = 0.8,
    test_frac: float = 0.1,
    valid_frac: float = 0.1,
    seed: int = 42,
) -> Dict:
    """
    Complete workflow: check, download, and process PrimeKG.

    Args:
        dataset_path: Path to dataset directory (e.g., /path/to/data/primekg1)
        force_redownload: If True, re-download even if files exist
        train_frac: Training set fraction (default: 0.8)
        test_frac: Test set fraction (default: 0.1)
        valid_frac: Validation set fraction (default: 0.1)
        seed: Random seed for reproducibility

    Returns:
        Dictionary with setup results:
        - status: str - "already_exists" or "completed"
        - nodes_count: int (if processed)
        - train_edges: int (if processed)
        - test_edges: int (if processed)
        - valid_edges: int (if processed)
        - output_path: str
    """
    # Check existing data
    status = check_primekg_data(dataset_path)

    if status["all_present"] and not force_redownload:
        logger.info("PrimeKG data already exists and is complete")
        logger.info(f"  Location: {dataset_path}")
        return {"status": "already_exists", "data_path": dataset_path, **status}

    # Get parent directory for CSV storage
    data_root = os.path.dirname(dataset_path)

    # Download if needed
    csv_path = os.path.join(data_root, "primekg.csv")
    if not os.path.exists(csv_path) or force_redownload:
        logger.info("Downloading PrimeKG...")
        download_primekg(data_root)
    else:
        logger.info(f"Using existing primekg.csv at: {csv_path}")

    # Process data
    logger.info("Processing PrimeKG...")
    output_path = os.path.join(dataset_path, "raw")
    results = process_primekg(
        csv_path, output_path, train_frac, test_frac, valid_frac, seed
    )

    logger.info("Setup complete!")
    logger.info(f"  Nodes: {results['nodes_count']:,}")
    logger.info(f"  Train edges: {results['train_edges']:,}")
    logger.info(f"  Test edges: {results['test_edges']:,}")
    logger.info(f"  Valid edges: {results['valid_edges']:,}")

    return {"status": "completed", **results}


if __name__ == "__main__":
    """Standalone execution for testing"""
    import argparse

    parser = argparse.ArgumentParser(description="PrimeKG ETL for ULTRA MCP Server")

    # Default to self-contained mcp/data directory
    mcp_dir = os.path.dirname(os.path.abspath(__file__))
    default_path = os.path.join(mcp_dir, "data", "primekg1")

    parser.add_argument(
        "--dataset_path",
        type=str,
        default=default_path,
        help="Path to dataset directory (default: mcp/data/primekg1)",
    )
    parser.add_argument(
        "--force_redownload",
        action="store_true",
        help="Force re-download of PrimeKG data even if it exists",
    )
    parser.add_argument("--train_frac", type=float, default=0.8)
    parser.add_argument("--test_frac", type=float, default=0.1)
    parser.add_argument("--valid_frac", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=42)

    args = parser.parse_args()

    # Setup logging for standalone execution
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    result = setup_primekg(
        dataset_path=args.dataset_path,
        force_redownload=args.force_redownload,
        train_frac=args.train_frac,
        test_frac=args.test_frac,
        valid_frac=args.valid_frac,
        seed=args.seed,
    )

    print("\nResults:")
    for key, value in result.items():
        print(f"  {key}: {value}")
