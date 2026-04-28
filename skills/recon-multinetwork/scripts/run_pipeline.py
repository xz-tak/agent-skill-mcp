#!/usr/bin/env python3
"""
ReCoN Multinetwork Pipeline Orchestrator

Chains modules M1-M7 sequentially with shared ReconConfig.
Supports partial runs via --start-from / --end-at flags.

Usage:
    # Full pipeline
    python run_pipeline.py --config config.json

    # Resume from module 4
    python run_pipeline.py --config config.json --start-from 4

    # Run only modules 4-7
    python run_pipeline.py --config config.json --start-from 4 --end-at 7

    # Run single module
    python run_pipeline.py --config config.json --start-from 3 --end-at 3
"""

import argparse
import sys
import traceback
from datetime import datetime
from pathlib import Path

from config import ReconConfig

# Lazy imports to avoid loading all modules at startup
_MODULE_REGISTRY = {
    1: "m1_data_prep",
    2: "m2_grn_pipeline",
    3: "m3_ccc_analysis",
    4: "m4_recon_analysis",
    5: "m5_multinetwork",
    6: "m6_differential",
    7: "m7_visualization",
    8: "m8_target_prediction",
}

_MODULE_NAMES = {
    1: "Data Preparation",
    2: "GRN Pipeline",
    3: "CCC Analysis",
    4: "ReCoN Analysis",
    5: "Multinetwork",
    6: "Differential Cascades",
    7: "Visualization",
    8: "Target Prediction",
}


def _import_module(module_num: int):
    """Lazy-import a pipeline module by number."""
    module_name = _MODULE_REGISTRY[module_num]
    return __import__(module_name)


def run_full_pipeline(
    config: ReconConfig,
    start_from: int = 1,
    end_at: int = 8,
) -> None:
    """
    Run the ReCoN pipeline from start_from to end_at (inclusive).

    Args:
        config: Pipeline configuration
        start_from: First module to run (1-7)
        end_at: Last module to run (1-7)
    """
    pipeline_start = datetime.now()

    print("=" * 60)
    print("ReCoN MULTINETWORK PIPELINE")
    print(f"Modules: {start_from} -> {end_at}")
    print(f"Started: {pipeline_start.isoformat()}")
    print("=" * 60)

    # Validate config
    errors = config.validate()
    if errors:
        print("\nConfig validation errors:")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)

    # Save initial config snapshot
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    config.to_json(output_dir / "config.json")

    completed = []
    failed = []

    for num in range(start_from, end_at + 1):
        name = _MODULE_NAMES[num]
        module_start = datetime.now()

        print(f"\n{'=' * 60}")
        print(f"MODULE {num}: {name}")
        print(f"{'=' * 60}")

        try:
            mod = _import_module(num)
            mod.main(config)
            elapsed = datetime.now() - module_start
            print(f"\nModule {num} completed in {elapsed}")
            completed.append(num)

            # Persist config after each module (modules may update it)
            config.to_json(output_dir / "config.json")

        except Exception as exc:
            elapsed = datetime.now() - module_start
            print(f"\nModule {num} FAILED after {elapsed}: {exc}")
            traceback.print_exc()
            failed.append(num)
            break

    # Summary
    total_elapsed = datetime.now() - pipeline_start
    print(f"\n{'=' * 60}")
    print("PIPELINE SUMMARY")
    print(f"{'=' * 60}")
    print(f"Completed modules: {completed}")
    if failed:
        print(f"Failed module: {failed}")
    print(f"Total time: {total_elapsed}")
    print(f"Config saved: {output_dir / 'config.json'}")


def main():
    parser = argparse.ArgumentParser(
        description="ReCoN Multinetwork Pipeline Orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_pipeline.py --config config.json
  python run_pipeline.py --config config.json --start-from 4
  python run_pipeline.py --config config.json --start-from 4 --end-at 7
        """,
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to JSON config file",
    )
    parser.add_argument(
        "--start-from",
        type=int,
        default=1,
        choices=range(1, 9),
        help="First module to run (1-8, default: 1)",
    )
    parser.add_argument(
        "--end-at",
        type=int,
        default=8,
        choices=range(1, 9),
        help="Last module to run (1-8, default: 8)",
    )

    args = parser.parse_args()

    if args.start_from > args.end_at:
        parser.error("--start-from must be <= --end-at")

    config = ReconConfig.from_json(args.config)
    run_full_pipeline(config, start_from=args.start_from, end_at=args.end_at)


if __name__ == "__main__":
    main()
