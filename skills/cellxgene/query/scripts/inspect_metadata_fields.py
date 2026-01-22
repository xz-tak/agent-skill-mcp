#!/usr/bin/env python3
"""
Quick script to inspect available metadata fields in CELLxGENE Census.

This script connects to CELLxGENE Census and displays all available
metadata fields (columns) without downloading large amounts of data.
"""

import cellxgene_census
import sys


def inspect_census_fields(organism="homo_sapiens", sample_size=100):
    """
    Inspect available metadata fields in CELLxGENE Census.

    Args:
        organism: Organism to inspect (default: homo_sapiens)
        sample_size: Number of rows to sample for inspection
    """
    print(f"\n{'='*80}")
    print(f"Inspecting CELLxGENE Census Metadata Fields")
    print(f"Organism: {organism}")
    print(f"{'='*80}\n")

    with cellxgene_census.open_soma(census_version="stable") as census:
        # Get experiment
        try:
            experiment = census["census_data"][organism]
        except KeyError:
            print(f"Error: Organism '{organism}' not found in Census")
            print(f"Available organisms: {list(census['census_data'].keys())}")
            return

        # Read a small sample of observation metadata
        print("Loading sample metadata...")
        obs_df = experiment.obs.read(coords=(slice(0, sample_size),)).concat().to_pandas()

        print(f"\n{'='*80}")
        print(f"OBSERVATION (CELL) METADATA FIELDS")
        print(f"{'='*80}\n")
        print(f"Total columns available: {len(obs_df.columns)}\n")

        # Display all column names with descriptions
        print(f"{'Field Name':<40} {'Data Type':<15} {'Example Values'}")
        print("-" * 80)

        for col in sorted(obs_df.columns):
            dtype = str(obs_df[col].dtype)

            # Get example values (non-null)
            example_vals = obs_df[col].dropna().unique()[:3]
            if len(example_vals) > 0:
                example_str = ", ".join([str(v)[:30] for v in example_vals])
                if len(example_str) > 40:
                    example_str = example_str[:37] + "..."
            else:
                example_str = "N/A"

            print(f"{col:<40} {dtype:<15} {example_str}")

        # Show value counts for key categorical fields
        print(f"\n{'='*80}")
        print(f"KEY FIELD STATISTICS (from sample of {sample_size} cells)")
        print(f"{'='*80}\n")

        key_fields = [
            'organism', 'tissue', 'tissue_general', 'cell_type',
            'assay', 'disease', 'sex', 'development_stage',
            'ethnicity', 'is_primary_data', 'organism_ontology_term_id'
        ]

        for field in key_fields:
            if field in obs_df.columns:
                print(f"\n{field.upper()}:")
                value_counts = obs_df[field].value_counts().head(10)
                for val, count in value_counts.items():
                    print(f"  {str(val):<50} {count:>5} cells")

                if len(obs_df[field].unique()) > 10:
                    print(f"  ... and {len(obs_df[field].unique()) - 10} more unique values")

        # Variable (gene) metadata
        print(f"\n{'='*80}")
        print(f"VARIABLE (GENE) METADATA FIELDS")
        print(f"{'='*80}\n")

        var_df = experiment.ms['RNA'].var.read(coords=(slice(0, 100),)).concat().to_pandas()

        print(f"Total columns available: {len(var_df.columns)}\n")

        print(f"{'Field Name':<40} {'Data Type':<15} {'Example Values'}")
        print("-" * 80)

        for col in sorted(var_df.columns):
            dtype = str(var_df[col].dtype)
            example_vals = var_df[col].dropna().unique()[:3]
            if len(example_vals) > 0:
                example_str = ", ".join([str(v)[:30] for v in example_vals])
                if len(example_str) > 40:
                    example_str = example_str[:37] + "..."
            else:
                example_str = "N/A"

            print(f"{col:<40} {dtype:<15} {example_str}")

        print(f"\n{'='*80}")
        print(f"SUMMARY")
        print(f"{'='*80}\n")
        print(f"Total observation fields: {len(obs_df.columns)}")
        print(f"Total variable fields: {len(var_df.columns)}")
        print(f"Sample size inspected: {len(obs_df)} cells, {len(var_df)} genes")
        print(f"\nThese fields are available in the metadata CSV file from queries.\n")


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Inspect CELLxGENE Census metadata fields'
    )
    parser.add_argument('--organism', type=str, default='homo_sapiens',
                        help='Organism to inspect (default: homo_sapiens, also available: mus_musculus)')
    parser.add_argument('--sample-size', type=int, default=100,
                        help='Number of cells to sample (default: 100)')

    args = parser.parse_args()

    try:
        inspect_census_fields(args.organism, args.sample_size)
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
