#!/usr/bin/env python3
"""
Example usage of CELLxGENE query script.

This script demonstrates how to use the CellxGeneQuery class programmatically.
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))

from query_cellxgene import CellxGeneQuery


def example_1_simple_tissue_query():
    """Example 1: Simple query for lung tissue."""
    print("\n" + "="*80)
    print("Example 1: Query human lung tissue (adult)")
    print("="*80 + "\n")

    with CellxGeneQuery() as querier:
        result = querier.query_data(
            tissue="lung",
        )

        if result:
            experiment, value_filter, var_filter, obs_df = result

            # Query the data
            with experiment.axis_query(
                measurement_name="RNA",
                obs_query=None if value_filter is None else {"value_filter": value_filter},
            ) as query:
                adata = query.to_anndata(X_name="raw")
                obs_df_filtered = query.obs().concat().to_pandas()

            # Generate and print summary
            summary = querier.generate_summary(adata, obs_df_filtered, value_filter or "None")
            querier.print_summary(summary)

            # Save results
            querier.save_results(
                adata,
                obs_df_filtered,
                summary,
                output_dir=Path("./output"),
                prefix="example1_lung"
            )


def example_2_cell_type_query():
    """Example 2: Query specific cell type across multiple tissues."""
    print("\n" + "="*80)
    print("Example 2: Query T cells from lung and intestine")
    print("="*80 + "\n")

    with CellxGeneQuery() as querier:
        result = querier.query_data(
            tissue=["lung", "intestine"],
            cell_type="T cell",
        )

        if result:
            experiment, value_filter, var_filter, obs_df = result

            with experiment.axis_query(
                measurement_name="RNA",
                obs_query=None if value_filter is None else {"value_filter": value_filter},
            ) as query:
                adata = query.to_anndata(X_name="raw")
                obs_df_filtered = query.obs().concat().to_pandas()

            summary = querier.generate_summary(adata, obs_df_filtered, value_filter or "None")
            querier.print_summary(summary)


def example_3_gene_set_query():
    """Example 3: Query with specific gene set."""
    print("\n" + "="*80)
    print("Example 3: Query liver with APOE gene family")
    print("="*80 + "\n")

    with CellxGeneQuery() as querier:
        result = querier.query_data(
            tissue="liver",
            gene_set=["APOE", "APOC1", "APOC2", "APOC3"],
        )

        if result:
            experiment, value_filter, var_filter, obs_df = result

            with experiment.axis_query(
                measurement_name="RNA",
                obs_query=None if value_filter is None else {"value_filter": value_filter},
                var_query=None if var_filter is None else {"value_filter": var_filter},
            ) as query:
                adata = query.to_anndata(X_name="raw")
                obs_df_filtered = query.obs().concat().to_pandas()

            summary = querier.generate_summary(adata, obs_df_filtered, value_filter or "None")
            querier.print_summary(summary)

            print(f"\nGenes in result: {list(adata.var_names)}")


def example_4_disease_query():
    """Example 4: Query disease samples."""
    print("\n" + "="*80)
    print("Example 4: Query COVID-19 lung samples")
    print("="*80 + "\n")

    with CellxGeneQuery() as querier:
        result = querier.query_data(
            tissue="lung",
            disease="COVID-19",
        )

        if result:
            experiment, value_filter, var_filter, obs_df = result

            with experiment.axis_query(
                measurement_name="RNA",
                obs_query=None if value_filter is None else {"value_filter": value_filter},
            ) as query:
                adata = query.to_anndata(X_name="raw")
                obs_df_filtered = query.obs().concat().to_pandas()

            summary = querier.generate_summary(adata, obs_df_filtered, value_filter or "None")
            querier.print_summary(summary)


def example_5_multi_species_query():
    """Example 5: Query multiple species."""
    print("\n" + "="*80)
    print("Example 5: Query brain neurons from human and mouse")
    print("="*80 + "\n")

    with CellxGeneQuery() as querier:
        result = querier.query_data(
            species=["human", "mouse"],
            tissue="brain",
            cell_type="neuron",
        )

        if result:
            experiment, value_filter, var_filter, obs_df = result

            with experiment.axis_query(
                measurement_name="RNA",
                obs_query=None if value_filter is None else {"value_filter": value_filter},
            ) as query:
                adata = query.to_anndata(X_name="raw")
                obs_df_filtered = query.obs().concat().to_pandas()

            summary = querier.generate_summary(adata, obs_df_filtered, value_filter or "None")
            querier.print_summary(summary)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run example queries")
    parser.add_argument("--example", type=int, default=None,
                        help="Run specific example (1-5), or all if not specified")

    args = parser.parse_args()

    examples = {
        1: example_1_simple_tissue_query,
        2: example_2_cell_type_query,
        3: example_3_gene_set_query,
        4: example_4_disease_query,
        5: example_5_multi_species_query,
    }

    if args.example:
        if args.example in examples:
            examples[args.example]()
        else:
            print(f"Example {args.example} not found. Available: 1-5")
    else:
        print("Running all examples...")
        print("WARNING: This will take a while and download significant data!")
        print("Press Ctrl+C to cancel, or Enter to continue...", end='')
        input()

        for i, func in examples.items():
            try:
                func()
            except KeyboardInterrupt:
                print("\nExamples interrupted by user")
                break
            except Exception as e:
                print(f"\nExample {i} failed with error: {e}")
                continue
