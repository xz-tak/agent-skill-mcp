#!/usr/bin/env python3
"""
Multi-Dataset DEG Analysis Adapter

Bridge script that adapts DESeq2 output for multi-dataset comparison analysis.
This enables comparison of in-house DESeq2 results with external datasets.

Usage:
    python multidata_adapter.py --deseq2_output /path/to/deseq2_output [options]

The adapter:
- Parses DESeq2 summstats_all.txt output
- Accepts external dataset configurations
- Generates unified config.json for multidata scripts
- Invokes generate_tables.py and generate_heatmaps.R (local scripts)

Modes:
- discovery: Find top N genes meeting reversal criteria
- gene_list: Analyze specific genes
- pathway_analysis: Run GSEA and compare pathway enrichments
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Optional


# =============================================================================
# CONFIGURATION
# =============================================================================

# Path to multidata scripts (same directory as this adapter)
SCRIPT_DIR = Path(__file__).parent.resolve()

# Default output directory name
DEFAULT_OUTPUT_DIR = "deg_multi_comparison"

# DESeq2 column mappings (DESeq2 output -> multidata expected)
DESEQ2_COLUMN_MAP = {
    "gene_col": "symbol",
    "log2fc_col": "log2FoldChange",
    "padj_col": "padj",
    "contrast_col": "comparison"
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def find_deseq2_summstats(deseq2_output: Path) -> Optional[Path]:
    """
    Find the summstats_all.txt file in DESeq2 output directory.

    Args:
        deseq2_output: Path to DESeq2 output directory

    Returns:
        Path to summstats file or None if not found
    """
    deg_dir = deseq2_output / "deg"

    # Look for *_summstats_all.txt pattern
    summstats_files = list(deg_dir.glob("*_summstats_all.txt"))

    if not summstats_files:
        return None

    # Return the first match (typically only one)
    return summstats_files[0]


def find_gsea_results(deseq2_output: Path) -> Optional[Path]:
    """
    Find the gsea_all.txt file in DESeq2 output directory.

    Args:
        deseq2_output: Path to DESeq2 output directory

    Returns:
        Path to GSEA file or None if not found
    """
    deg_dir = deseq2_output / "deg"

    gsea_files = list(deg_dir.glob("*_gsea_all.txt"))

    if not gsea_files:
        return None

    return gsea_files[0]


def create_deseq2_source_config(
    summstats_path: Path,
    comparisons: list,
    label_prefix: str = "InHouse"
) -> list:
    """
    Create data source configurations for DESeq2 comparisons.

    Args:
        summstats_path: Path to DESeq2 summstats_all.txt
        comparisons: List of comparison names to include
        label_prefix: Prefix for source labels

    Returns:
        List of data source configuration dictionaries
    """
    sources = []

    for comparison in comparisons:
        source = {
            "file": str(summstats_path),
            "gene_col": DESEQ2_COLUMN_MAP["gene_col"],
            "log2fc_col": DESEQ2_COLUMN_MAP["log2fc_col"],
            "padj_col": DESEQ2_COLUMN_MAP["padj_col"],
            "label": f"{label_prefix}_{comparison}",
            "contrast_filter": {
                "column": DESEQ2_COLUMN_MAP["contrast_col"],
                "value": comparison
            }
        }
        sources.append(source)

    return sources


def parse_external_config(config_path: Path) -> dict:
    """
    Parse external sources configuration file.

    Expected format:
    {
        "sources": [
            {
                "file": "/path/to/data.xlsx",
                "sheet": "Sheet1",
                "gene_col": "Gene",
                "log2fc_col": "log2FC",
                "padj_col": "FDR",
                "label": "External_Study1",
                "group": "External",
                "log2fc_cutoff": 0.5
            }
        ],
        "column_groups": {...},
        "score_logic": {...}
    }
    """
    with open(config_path) as f:
        return json.load(f)


def get_deseq2_comparisons(summstats_path: Path) -> list:
    """
    Extract available comparison names from DESeq2 summstats file.

    Args:
        summstats_path: Path to summstats_all.txt

    Returns:
        List of unique comparison names
    """
    import pandas as pd

    df = pd.read_csv(summstats_path, sep='\t')
    return df['comparison'].unique().tolist()


def build_config(
    mode: str,
    deseq2_sources: list,
    external_sources: list,
    column_groups: dict,
    score_logic: dict,
    output_dir: Path,
    prefix: str,
    genes: Optional[list] = None,
    top_n: int = 50,
    gsea_file: Optional[Path] = None,
    heatmap_config: Optional[dict] = None,
    column_mapping: Optional[dict] = None
) -> dict:
    """
    Build unified configuration for multidata analysis.

    Returns:
        Complete configuration dictionary
    """
    # Combine all data sources
    all_sources = deseq2_sources + external_sources

    config = {
        "mode": mode,
        "top_n": top_n,
        "data_sources": all_sources,
        "column_groups": column_groups,
        "score_logic": score_logic,
        "output": {
            "prefix": prefix,
            "directory": str(output_dir)
        }
    }

    if mode == "gene_list" and genes:
        config["genes"] = genes

    if mode == "pathway_analysis" and gsea_file:
        config["gsea_file"] = str(gsea_file)

    # Add heatmap configuration
    if heatmap_config:
        config["heatmap"] = heatmap_config
    else:
        # Default heatmap settings
        config["heatmap"] = {
            "color_scale": ["blue", "white", "red"],
            "row_annotation": "Score",
            "width": 14,
            "height": 10,
            "fontsize": 9,
            "star_fontsize": 6,
            "column_gap": 3
        }

    # Add column label mapping for display
    if column_mapping:
        config["column_mapping"] = column_mapping

    return config


def run_generate_tables(config_path: Path, output_dir: Path) -> bool:
    """
    Run generate_tables.py from deg-multidata-analysis.

    Returns:
        True if successful, False otherwise
    """
    script_path = SCRIPT_DIR / "generate_tables.py"

    if not script_path.exists():
        print(f"Error: generate_tables.py not found at {script_path}")
        return False

    cmd = [
        sys.executable,
        str(script_path),
        "--config", str(config_path)
    ]

    print(f"\nRunning: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)

    return result.returncode == 0


def run_generate_pathway_tables(config_path: Path, gsea_file: Path) -> bool:
    """
    Run generate_pathway_tables.py from deg-multidata-analysis.

    Returns:
        True if successful, False otherwise
    """
    script_path = SCRIPT_DIR / "generate_pathway_tables.py"

    if not script_path.exists():
        print(f"Error: generate_pathway_tables.py not found at {script_path}")
        return False

    cmd = [
        sys.executable,
        str(script_path),
        "--config", str(config_path),
        "--gsea_file", str(gsea_file)
    ]

    print(f"\nRunning: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)

    return result.returncode == 0


def run_generate_heatmaps(
    config_path: Path,
    conda_env: str = "r_env",
    pathway_mode: bool = False
) -> bool:
    """
    Run generate_heatmaps.R from deg-multidata-analysis.

    Returns:
        True if successful, False otherwise
    """
    script_path = SCRIPT_DIR / "generate_heatmaps.R"

    if not script_path.exists():
        print(f"Error: generate_heatmaps.R not found at {script_path}")
        return False

    cmd = [
        "conda", "run", "-n", conda_env,
        "Rscript", str(script_path),
        "--config", str(config_path)
    ]

    if pathway_mode:
        cmd.append("--pathway_mode")

    print(f"\nRunning: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)

    return result.returncode == 0


# =============================================================================
# MAIN FUNCTIONS
# =============================================================================

def run_discovery_mode(
    deseq2_output: Path,
    external_config: Optional[dict],
    comparisons: list,
    output_dir: Path,
    prefix: str,
    top_n: int,
    conda_env: str,
    heatmap_config: Optional[dict] = None
) -> bool:
    """Run discovery mode analysis."""
    print("\n=== Discovery Mode ===")
    print(f"Finding top {top_n} genes with reversal patterns")

    # Find DESeq2 summstats
    summstats_path = find_deseq2_summstats(deseq2_output)
    if not summstats_path:
        print(f"Error: Could not find summstats_all.txt in {deseq2_output}/deg/")
        return False

    print(f"DESeq2 summstats: {summstats_path}")

    # Create DESeq2 source configs
    deseq2_sources = create_deseq2_source_config(summstats_path, comparisons)

    # Get external sources
    external_sources = external_config.get("sources", []) if external_config else []

    # Get column groups, score logic, heatmap config, and column mapping
    if external_config:
        column_groups = external_config.get("column_groups", {})
        score_logic = external_config.get("score_logic", {})
        # Merge heatmap config from external config with command-line overrides
        ext_heatmap = external_config.get("heatmap", {})
        if heatmap_config:
            ext_heatmap.update(heatmap_config)
        heatmap_config = ext_heatmap if ext_heatmap else heatmap_config
        column_mapping = external_config.get("column_mapping", {})
    else:
        # Auto-generate column groups from DESeq2 comparisons
        column_groups = {
            "InHouse": [f"InHouse_{c}" for c in comparisons]
        }
        score_logic = {
            "UP": {"InHouse": "up"},
            "DOWN": {"InHouse": "down"}
        }
        column_mapping = None

    # Build config
    config = build_config(
        mode="discovery",
        deseq2_sources=deseq2_sources,
        external_sources=external_sources,
        column_groups=column_groups,
        score_logic=score_logic,
        output_dir=output_dir,
        prefix=prefix,
        top_n=top_n,
        heatmap_config=heatmap_config,
        column_mapping=column_mapping
    )

    # Save config
    output_dir.mkdir(parents=True, exist_ok=True)
    config_path = output_dir / "config.json"
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)
    print(f"Config saved: {config_path}")

    # Run table generation
    if not run_generate_tables(config_path, output_dir):
        print("Error: Table generation failed")
        return False

    # Run heatmap generation
    if not run_generate_heatmaps(config_path, conda_env):
        print("Warning: Heatmap generation failed")

    return True


def run_gene_list_mode(
    deseq2_output: Path,
    external_config: Optional[dict],
    comparisons: list,
    genes: list,
    output_dir: Path,
    prefix: str,
    conda_env: str,
    heatmap_config: Optional[dict] = None
) -> bool:
    """Run gene list mode analysis."""
    print("\n=== Gene List Mode ===")
    print(f"Analyzing {len(genes)} genes")

    # Find DESeq2 summstats
    summstats_path = find_deseq2_summstats(deseq2_output)
    if not summstats_path:
        print(f"Error: Could not find summstats_all.txt in {deseq2_output}/deg/")
        return False

    print(f"DESeq2 summstats: {summstats_path}")

    # Create DESeq2 source configs
    deseq2_sources = create_deseq2_source_config(summstats_path, comparisons)

    # Get external sources
    external_sources = external_config.get("sources", []) if external_config else []

    # Get column groups, score logic, heatmap config, and column mapping
    if external_config:
        column_groups = external_config.get("column_groups", {})
        score_logic = external_config.get("score_logic", {})
        ext_heatmap = external_config.get("heatmap", {})
        if heatmap_config:
            ext_heatmap.update(heatmap_config)
        heatmap_config = ext_heatmap if ext_heatmap else heatmap_config
        column_mapping = external_config.get("column_mapping", {})
    else:
        column_groups = {
            "InHouse": [f"InHouse_{c}" for c in comparisons]
        }
        score_logic = {
            "UP": {"InHouse": "up"},
            "DOWN": {"InHouse": "down"}
        }
        column_mapping = None

    # Build config
    config = build_config(
        mode="gene_list",
        deseq2_sources=deseq2_sources,
        external_sources=external_sources,
        column_groups=column_groups,
        score_logic=score_logic,
        output_dir=output_dir,
        prefix=prefix,
        genes=genes,
        heatmap_config=heatmap_config,
        column_mapping=column_mapping
    )

    # Save config
    output_dir.mkdir(parents=True, exist_ok=True)
    config_path = output_dir / "config.json"
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)
    print(f"Config saved: {config_path}")

    # Run table generation
    if not run_generate_tables(config_path, output_dir):
        print("Error: Table generation failed")
        return False

    # Run heatmap generation
    if not run_generate_heatmaps(config_path, conda_env):
        print("Warning: Heatmap generation failed")

    return True


def run_pathway_mode(
    deseq2_output: Path,
    external_config: Optional[dict],
    comparisons: list,
    output_dir: Path,
    prefix: str,
    top_n: int,
    conda_env: str,
    heatmap_config: Optional[dict] = None
) -> bool:
    """Run pathway analysis mode."""
    print("\n=== Pathway Analysis Mode ===")
    print(f"Finding top {top_n} pathways with reversal patterns")

    # Find DESeq2 GSEA results
    gsea_path = find_gsea_results(deseq2_output)
    if not gsea_path:
        print(f"Error: Could not find gsea_all.txt in {deseq2_output}/deg/")
        print("Run DESeq2 analysis with GSEA enabled first.")
        return False

    print(f"DESeq2 GSEA results: {gsea_path}")

    # For pathway mode, we use the GSEA results directly
    # Column groups should map to the 'source' column values in gsea_all.txt

    if external_config:
        column_groups = external_config.get("column_groups", {})
        score_logic = external_config.get("score_logic", {})
        ext_heatmap = external_config.get("heatmap", {})
        if heatmap_config:
            ext_heatmap.update(heatmap_config)
        heatmap_config = ext_heatmap if ext_heatmap else heatmap_config
        column_mapping = external_config.get("column_mapping", {})
    else:
        # Auto-generate from comparisons
        column_groups = {
            "InHouse": comparisons
        }
        score_logic = {
            "UP": {"InHouse": "up"},
            "DOWN": {"InHouse": "down"}
        }
        column_mapping = None

    # Default heatmap settings if not provided
    if not heatmap_config:
        heatmap_config = {
            "color_scale": ["blue", "white", "red"],
            "row_annotation": "Score",
            "width": 14,
            "height": 10,
            "fontsize": 9,
            "star_fontsize": 6,
            "column_gap": 3
        }

    # Build config
    config = {
        "mode": "pathway_analysis",
        "top_n": top_n,
        "column_groups": column_groups,
        "score_logic": score_logic,
        "heatmap": heatmap_config,
        "output": {
            "prefix": prefix,
            "directory": str(output_dir)
        }
    }

    if column_mapping:
        config["column_mapping"] = column_mapping

    # Save config
    output_dir.mkdir(parents=True, exist_ok=True)
    config_path = output_dir / "config.json"
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)
    print(f"Config saved: {config_path}")

    # Run pathway table generation
    if not run_generate_pathway_tables(config_path, gsea_path):
        print("Error: Pathway table generation failed")
        return False

    # Run heatmap generation (pathway mode)
    if not run_generate_heatmaps(config_path, conda_env, pathway_mode=True):
        print("Warning: Heatmap generation failed")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Multi-Dataset DEG Analysis Adapter",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Discovery mode (default output: deg_multi_comparison/)
    python multidata_adapter.py --deseq2_output /path/to/output --mode discovery --top_n 50

    # Gene list mode
    python multidata_adapter.py --deseq2_output /path/to/output --mode gene_list --genes "GREM1,IL11,NOG"

    # Pathway analysis mode
    python multidata_adapter.py --deseq2_output /path/to/output --mode pathway_analysis --top_n 100

    # With external datasets
    python multidata_adapter.py --deseq2_output /path/to/output --external_config external.json

    # Custom output directory
    python multidata_adapter.py --deseq2_output /path/to/output --output_dir /custom/path
        """
    )

    parser.add_argument(
        "--deseq2_output",
        required=True,
        help="Path to DESeq2 output directory"
    )
    parser.add_argument(
        "--external_config",
        help="Path to external sources config JSON"
    )
    parser.add_argument(
        "--mode",
        choices=["discovery", "gene_list", "pathway_analysis"],
        default="discovery",
        help="Analysis mode (default: discovery)"
    )
    parser.add_argument(
        "--genes",
        help="Comma-separated gene list (gene_list mode)"
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=50,
        help="Number of top genes/pathways (default: 50)"
    )
    parser.add_argument(
        "--comparisons",
        help="Comma-separated comparison names to include (default: all)"
    )
    parser.add_argument(
        "--output_dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})"
    )
    parser.add_argument(
        "--prefix",
        default="multidata",
        help="Output file prefix (default: multidata)"
    )
    parser.add_argument(
        "--conda_env",
        default="r_env",
        help="Conda environment for R scripts (default: r_env)"
    )

    # Heatmap configuration arguments
    parser.add_argument(
        "--heatmap_colors",
        default="blue,white,red",
        help="Comma-separated color scale (default: blue,white,red)"
    )
    parser.add_argument(
        "--heatmap_width",
        type=float,
        default=14,
        help="Heatmap width in inches (default: 14)"
    )
    parser.add_argument(
        "--heatmap_height",
        type=float,
        default=10,
        help="Heatmap height in inches (default: 10)"
    )
    parser.add_argument(
        "--heatmap_fontsize",
        type=int,
        default=9,
        help="Font size for labels (default: 9)"
    )
    parser.add_argument(
        "--heatmap_star_fontsize",
        type=int,
        default=6,
        help="Font size for significance stars (default: 6)"
    )
    parser.add_argument(
        "--heatmap_column_gap",
        type=int,
        default=3,
        help="Gap between column groups in mm (default: 3)"
    )

    args = parser.parse_args()

    # Convert paths
    deseq2_output = Path(args.deseq2_output)
    output_dir = Path(args.output_dir)

    if not deseq2_output.exists():
        print(f"Error: DESeq2 output directory not found: {deseq2_output}")
        sys.exit(1)

    # Load external config if provided
    external_config = None
    if args.external_config:
        external_config_path = Path(args.external_config)
        if not external_config_path.exists():
            print(f"Error: External config not found: {external_config_path}")
            sys.exit(1)
        external_config = parse_external_config(external_config_path)

    # Get comparisons
    summstats_path = find_deseq2_summstats(deseq2_output)
    if summstats_path:
        available_comparisons = get_deseq2_comparisons(summstats_path)
        print(f"Available comparisons: {available_comparisons}")
    else:
        available_comparisons = []

    if args.comparisons:
        comparisons = [c.strip() for c in args.comparisons.split(",")]
    else:
        comparisons = available_comparisons

    # Build heatmap configuration from command-line arguments
    heatmap_config = {
        "color_scale": [c.strip() for c in args.heatmap_colors.split(",")],
        "row_annotation": "Score",
        "width": args.heatmap_width,
        "height": args.heatmap_height,
        "fontsize": args.heatmap_fontsize,
        "star_fontsize": args.heatmap_star_fontsize,
        "column_gap": args.heatmap_column_gap
    }

    # Run appropriate mode
    success = False

    if args.mode == "discovery":
        success = run_discovery_mode(
            deseq2_output=deseq2_output,
            external_config=external_config,
            comparisons=comparisons,
            output_dir=output_dir,
            prefix=args.prefix,
            top_n=args.top_n,
            conda_env=args.conda_env,
            heatmap_config=heatmap_config
        )
    elif args.mode == "gene_list":
        if not args.genes:
            print("Error: --genes required for gene_list mode")
            sys.exit(1)
        genes = [g.strip() for g in args.genes.split(",")]
        success = run_gene_list_mode(
            deseq2_output=deseq2_output,
            external_config=external_config,
            comparisons=comparisons,
            genes=genes,
            output_dir=output_dir,
            prefix=args.prefix,
            conda_env=args.conda_env,
            heatmap_config=heatmap_config
        )
    elif args.mode == "pathway_analysis":
        success = run_pathway_mode(
            deseq2_output=deseq2_output,
            external_config=external_config,
            comparisons=comparisons,
            output_dir=output_dir,
            prefix=args.prefix,
            top_n=args.top_n,
            conda_env=args.conda_env,
            heatmap_config=heatmap_config
        )

    if success:
        print(f"\n=== Analysis Complete ===")
        print(f"Results saved to: {output_dir}")
    else:
        print("\n=== Analysis Failed ===")
        sys.exit(1)


if __name__ == "__main__":
    main()
