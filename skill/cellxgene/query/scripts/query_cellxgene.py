#!/usr/bin/env python3
"""
CELLxGENE Census Query Script

Query single-cell data from CELLxGENE Census with flexible filtering options.
Supports filtering by species, tissue, cell type, disease, sex, development stage,
drug treatment, and gene sets.
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Union, List, Set, Optional, Dict, Any

import cellxgene_census
import numpy as np
import pandas as pd
import tiledbsoma as soma


# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class CellxGeneQuery:
    """Query CELLxGENE Census with flexible filtering."""

    def __init__(self):
        self.census = None
        self.filters = {}
        self.query_log = []

    def __enter__(self):
        """Open census connection."""
        logger.info("Opening CELLxGENE Census connection...")
        self.census = cellxgene_census.open_soma(census_version="stable")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close census connection."""
        if self.census:
            self.census.close()
            logger.info("Closed CELLxGENE Census connection")

    @staticmethod
    def normalize_input(value: Union[str, List, Set, None]) -> Optional[Set[str]]:
        """Normalize input to a set of lowercase strings."""
        if value is None:
            return None
        if isinstance(value, str):
            return {value.lower()}
        if isinstance(value, (list, set)):
            return {str(v).lower() for v in value}
        return {str(value).lower()}

    @staticmethod
    def partial_match(query: str, target: str) -> bool:
        """Check if query partially matches target (case-insensitive)."""
        return query.lower() in target.lower()

    @staticmethod
    def any_match(query_set: Set[str], target: str) -> bool:
        """Check if any element in query_set matches target (case-insensitive)."""
        target_lower = target.lower()
        return any(q in target_lower for q in query_set)

    def interpret_development_stage(self, stage_input: Union[str, List, Set, None]) -> Optional[Set[str]]:
        """
        Interpret development stage input and map to appropriate terms.

        Common mappings:
        - adult -> adult, mature
        - embryonic -> embryo, embryonic
        - fetal -> fetus, fetal
        - postnatal -> postnatal, newborn, infant, child
        """
        if stage_input is None:
            return None

        stage_set = self.normalize_input(stage_input)
        expanded_stages = set()

        stage_mappings = {
            'adult': ['adult', 'mature'],
            'embryonic': ['embryo', 'embryonic'],
            'fetal': ['fetus', 'fetal'],
            'postnatal': ['postnatal', 'newborn', 'infant', 'child'],
            'prenatal': ['prenatal', 'embryo', 'fetus'],
        }

        for stage in stage_set:
            # Add the original term
            expanded_stages.add(stage)

            # Add expanded terms if mapping exists
            for key, values in stage_mappings.items():
                if stage in key or key in stage:
                    expanded_stages.update(values)

        return expanded_stages

    def build_value_filter(self, obs_df: pd.DataFrame) -> str:
        """
        Build SOMA value filter string based on provided criteria.

        Args:
            obs_df: Observation metadata DataFrame to check available values

        Returns:
            Filter string for SOMA query
        """
        filter_conditions = []

        # Species filter (default to human if not specified)
        species = self.filters.get('species')
        if species is None:
            species = {'homo sapiens'}
            logger.info("No species specified, defaulting to: Homo sapiens")

        if species:
            # Check for any match in organism
            matching_organisms = set()
            for org in obs_df['organism'].unique():
                if self.any_match(species, str(org)):
                    matching_organisms.add(org)

            if matching_organisms:
                organism_filter = ' or '.join([f'organism == "{org}"' for org in matching_organisms])
                filter_conditions.append(f'({organism_filter})')
                logger.info(f"Species filter matched: {matching_organisms}")
            else:
                logger.warning(f"No matching organisms found for: {species}")

        # Tissue filter
        tissue = self.filters.get('tissue')
        if tissue:
            matching_tissues = set()
            for tis in obs_df['tissue_general'].unique():
                if len(tissue) == 1 and isinstance(self.filters.get('tissue'), str):
                    # Partial match for string
                    if self.partial_match(list(tissue)[0], str(tis)):
                        matching_tissues.add(tis)
                else:
                    # Any match for list/set
                    if self.any_match(tissue, str(tis)):
                        matching_tissues.add(tis)

            if matching_tissues:
                tissue_filter = ' or '.join([f'tissue_general == "{t}"' for t in matching_tissues])
                filter_conditions.append(f'({tissue_filter})')
                logger.info(f"Tissue filter matched: {matching_tissues}")
            else:
                logger.warning(f"No matching tissues found for: {tissue}")

        # Cell type filter
        cell_type = self.filters.get('cell_type')
        if cell_type:
            matching_cell_types = set()
            for ct in obs_df['cell_type'].unique():
                if len(cell_type) == 1 and isinstance(self.filters.get('cell_type'), str):
                    # Partial match for string
                    if self.partial_match(list(cell_type)[0], str(ct)):
                        matching_cell_types.add(ct)
                else:
                    # Any match for list/set
                    if self.any_match(cell_type, str(ct)):
                        matching_cell_types.add(ct)

            if matching_cell_types:
                ct_filter = ' or '.join([f'cell_type == "{ct}"' for ct in matching_cell_types])
                filter_conditions.append(f'({ct_filter})')
                logger.info(f"Cell type filter matched: {matching_cell_types}")
            else:
                logger.warning(f"No matching cell types found for: {cell_type}")

        # Disease filter
        disease = self.filters.get('disease')
        if disease:
            matching_diseases = set()
            for dis in obs_df['disease'].unique():
                if len(disease) == 1 and isinstance(self.filters.get('disease'), str):
                    # Partial match for string
                    if self.partial_match(list(disease)[0], str(dis)):
                        matching_diseases.add(dis)
                else:
                    # Any match for list/set
                    if self.any_match(disease, str(dis)):
                        matching_diseases.add(dis)

            if matching_diseases:
                disease_filter = ' or '.join([f'disease == "{d}"' for d in matching_diseases])
                filter_conditions.append(f'({disease_filter})')
                logger.info(f"Disease filter matched: {matching_diseases}")
            else:
                logger.warning(f"No matching diseases found for: {disease}")

        # Sex filter
        sex = self.filters.get('sex')
        if sex:
            matching_sex = set()
            for s in obs_df['sex'].unique():
                if self.any_match(sex, str(s)):
                    matching_sex.add(s)

            if matching_sex:
                sex_filter = ' or '.join([f'sex == "{s}"' for s in matching_sex])
                filter_conditions.append(f'({sex_filter})')
                logger.info(f"Sex filter matched: {matching_sex}")
            else:
                logger.warning(f"No matching sex found for: {sex}")

        # Development stage filter (default to adult if not specified)
        dev_stage = self.filters.get('development_stage')
        if dev_stage is None:
            dev_stage = {'adult'}
            logger.info("No development stage specified, defaulting to: adult")

        if dev_stage:
            dev_stage_expanded = self.interpret_development_stage(list(dev_stage)[0] if len(dev_stage) == 1 else dev_stage)
            matching_stages = set()

            for stage in obs_df['development_stage'].unique():
                if self.any_match(dev_stage_expanded, str(stage)):
                    matching_stages.add(stage)

            if matching_stages:
                stage_filter = ' or '.join([f'development_stage == "{s}"' for s in matching_stages])
                filter_conditions.append(f'({stage_filter})')
                logger.info(f"Development stage filter matched: {matching_stages}")
            else:
                logger.warning(f"No matching development stages found for: {dev_stage}")

        # Drug treatment filter (using suspension_type as proxy, or self_reported_ethnicity_ontology_term_id)
        drug_treatment = self.filters.get('drug_treatment')
        if drug_treatment:
            # Note: CELLxGENE Census doesn't have a direct drug_treatment field
            # You might need to filter this from dataset metadata or use custom fields
            logger.warning("Drug treatment filtering is not directly supported in Census schema. Skipping.")

        if not filter_conditions:
            return None

        return ' and '.join(filter_conditions)

    def query_data(
        self,
        species: Union[str, List, Set, None] = None,
        tissue: Union[str, List, Set, None] = None,
        cell_type: Union[str, List, Set, None] = None,
        disease: Union[str, List, Set, None] = None,
        sex: Union[str, List, Set, None] = None,
        development_stage: Union[str, List, Set, None] = None,
        drug_treatment: Union[str, List, Set, None] = None,
        gene_set: Union[List, Set, None] = None,
    ) -> Optional[soma.Experiment]:
        """
        Query CELLxGENE Census with specified filters.

        Args:
            species: Species name(s) (default: human)
            tissue: Tissue name(s) (default: None)
            cell_type: Cell type name(s) (default: None)
            disease: Disease name(s) (default: None)
            sex: Sex (default: None)
            development_stage: Development stage(s) (default: adult)
            drug_treatment: Drug treatment(s) (default: None)
            gene_set: List of gene symbols to include (default: None, all genes)

        Returns:
            Query result or None if query fails
        """
        # Store filters
        self.filters = {
            'species': self.normalize_input(species),
            'tissue': self.normalize_input(tissue),
            'cell_type': self.normalize_input(cell_type),
            'disease': self.normalize_input(disease),
            'sex': self.normalize_input(sex),
            'development_stage': self.normalize_input(development_stage),
            'drug_treatment': self.normalize_input(drug_treatment),
            'gene_set': gene_set if gene_set is None else list(gene_set),
        }

        # Log filters
        self.query_log.append({
            'timestamp': datetime.now().isoformat(),
            'filters': {k: list(v) if isinstance(v, set) else v for k, v in self.filters.items()}
        })

        logger.info("=" * 80)
        logger.info("Query Filters:")
        for key, value in self.filters.items():
            if value is not None:
                logger.info(f"  {key}: {value}")
        logger.info("=" * 80)

        # Get organism (default to human)
        # Map species names to Census organism keys
        species_map = {
            'homo sapiens': 'homo_sapiens',
            'human': 'homo_sapiens',
            'mus musculus': 'mus_musculus',
            'mouse': 'mus_musculus',
        }

        organism_key = "homo_sapiens"  # Default
        organism_display = "Homo sapiens"

        if species:
            # Use first species if multiple
            organism_candidates = list(self.filters['species'])
            if organism_candidates:
                query_species = organism_candidates[0]
                organism_key = species_map.get(query_species, query_species)
                organism_display = query_species.replace('_', ' ').title()

        logger.info(f"Querying Census for organism: {organism_display} (key: {organism_key})")

        # Access the experiment
        try:
            experiment = self.census["census_data"][organism_key]
        except KeyError:
            logger.error(f"Organism '{organism_key}' not found in Census")
            logger.error(f"Available organisms: {list(self.census['census_data'].keys())}")
            return None

        # Get observation metadata to build filters
        logger.info("Loading observation metadata for filter building...")
        obs_df = experiment.obs.read().concat().to_pandas()
        logger.info(f"Total observations in Census: {len(obs_df)}")

        # Build value filter
        value_filter = self.build_value_filter(obs_df)

        if value_filter:
            logger.info(f"Applied filter: {value_filter}")
        else:
            logger.info("No filters applied, querying all data")

        # Build gene filter if gene_set is provided
        var_filter = None
        if gene_set:
            gene_symbols = [f'feature_name == "{gene}"' for gene in gene_set]
            var_filter = ' or '.join(gene_symbols)
            logger.info(f"Gene filter: {len(gene_set)} genes specified")

        return experiment, value_filter, var_filter, obs_df

    def generate_summary(
        self,
        adata,
        obs_df_filtered: pd.DataFrame,
        filters_applied: str
    ) -> Dict[str, Any]:
        """
        Generate summary statistics for queried data.

        Args:
            adata: AnnData object with query results
            obs_df_filtered: Filtered observation metadata
            filters_applied: String describing filters applied

        Returns:
            Dictionary with summary statistics
        """
        summary = {
            'timestamp': datetime.now().isoformat(),
            'filters_applied': filters_applied,
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'organisms': obs_df_filtered['organism'].value_counts().to_dict(),
            'tissues': obs_df_filtered['tissue_general'].value_counts().to_dict(),
            'cell_types': obs_df_filtered['cell_type'].value_counts().to_dict(),
            'diseases': obs_df_filtered['disease'].value_counts().to_dict(),
            'sexes': obs_df_filtered['sex'].value_counts().to_dict(),
            'development_stages': obs_df_filtered['development_stage'].value_counts().to_dict(),
            'assays': obs_df_filtered['assay'].value_counts().to_dict(),
            'donors': obs_df_filtered['donor_id'].nunique(),
        }

        return summary

    def print_summary(self, summary: Dict[str, Any]):
        """Print summary in a readable format."""
        logger.info("=" * 80)
        logger.info("QUERY RESULTS SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Timestamp: {summary['timestamp']}")
        logger.info(f"Filters Applied: {summary['filters_applied']}")
        logger.info(f"\nTotal Cells: {summary['n_cells']:,}")
        logger.info(f"Total Genes: {summary['n_genes']:,}")
        logger.info(f"Unique Donors: {summary['donors']:,}")

        logger.info(f"\nOrganisms:")
        for org, count in summary['organisms'].items():
            logger.info(f"  {org}: {count:,}")

        logger.info(f"\nTissues:")
        for tissue, count in sorted(summary['tissues'].items(), key=lambda x: x[1], reverse=True)[:10]:
            logger.info(f"  {tissue}: {count:,}")
        if len(summary['tissues']) > 10:
            logger.info(f"  ... and {len(summary['tissues']) - 10} more")

        logger.info(f"\nCell Types:")
        for ct, count in sorted(summary['cell_types'].items(), key=lambda x: x[1], reverse=True)[:10]:
            logger.info(f"  {ct}: {count:,}")
        if len(summary['cell_types']) > 10:
            logger.info(f"  ... and {len(summary['cell_types']) - 10} more")

        logger.info(f"\nDiseases:")
        for disease, count in sorted(summary['diseases'].items(), key=lambda x: x[1], reverse=True)[:10]:
            logger.info(f"  {disease}: {count:,}")

        logger.info(f"\nSexes:")
        for sex, count in summary['sexes'].items():
            logger.info(f"  {sex}: {count:,}")

        logger.info(f"\nDevelopment Stages:")
        for stage, count in summary['development_stages'].items():
            logger.info(f"  {stage}: {count:,}")

        logger.info(f"\nAssays:")
        for assay, count in sorted(summary['assays'].items(), key=lambda x: x[1], reverse=True)[:5]:
            logger.info(f"  {assay}: {count:,}")

        logger.info("=" * 80)

    def save_results(
        self,
        adata,
        obs_df_filtered: pd.DataFrame,
        summary: Dict[str, Any],
        output_dir: Path,
        prefix: str = "cellxgene_query"
    ):
        """
        Save query results to disk.

        Args:
            adata: AnnData object with query results
            obs_df_filtered: Filtered observation metadata
            summary: Summary statistics dictionary
            output_dir: Directory to save results
            prefix: Prefix for output files
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save AnnData object
        adata_path = output_dir / f"{prefix}_data.h5ad"
        logger.info(f"Saving AnnData object to: {adata_path}")
        adata.write_h5ad(adata_path)

        # Save metadata
        metadata_path = output_dir / f"{prefix}_metadata.csv"
        logger.info(f"Saving metadata to: {metadata_path}")
        obs_df_filtered.to_csv(metadata_path, index=False)

        # Save summary as JSON
        summary_path = output_dir / f"{prefix}_summary.json"
        logger.info(f"Saving summary to: {summary_path}")
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)

        # Save log
        log_path = output_dir / f"{prefix}_log.json"
        logger.info(f"Saving query log to: {log_path}")
        with open(log_path, 'w') as f:
            json.dump(self.query_log, f, indent=2)

        logger.info("All results saved successfully!")
        logger.info(f"Output directory: {output_dir}")


def main():
    """Main function to run query from command line."""
    parser = argparse.ArgumentParser(
        description='Query CELLxGENE Census with flexible filtering options',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Query human lung epithelial cells
  python query_cellxgene.py --tissue lung --cell-type epithelial

  # Query mouse and human brain neurons
  python query_cellxgene.py --species "mouse,human" --tissue brain --cell-type neuron

  # Query with specific genes
  python query_cellxgene.py --tissue liver --genes APOE,APOC1,APOC2 --output liver_apoe

  # Query with disease filter
  python query_cellxgene.py --tissue lung --disease "COVID-19" --development-stage adult
        """
    )

    # Filter arguments
    parser.add_argument('--species', type=str, default=None,
                        help='Species (comma-separated for multiple). Default: human')
    parser.add_argument('--tissue', type=str, default=None,
                        help='Tissue name (comma-separated for multiple)')
    parser.add_argument('--cell-type', type=str, default=None,
                        help='Cell type name (comma-separated for multiple)')
    parser.add_argument('--disease', type=str, default=None,
                        help='Disease name (comma-separated for multiple)')
    parser.add_argument('--sex', type=str, default=None,
                        help='Sex (comma-separated for multiple)')
    parser.add_argument('--development-stage', type=str, default=None,
                        help='Development stage (comma-separated for multiple). Default: adult')
    parser.add_argument('--drug-treatment', type=str, default=None,
                        help='Drug treatment (comma-separated for multiple)')
    parser.add_argument('--genes', type=str, default=None,
                        help='Gene symbols (comma-separated)')

    # Output arguments
    parser.add_argument('--output', type=str, default='cellxgene_query',
                        help='Output prefix for saved files')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Output directory')
    parser.add_argument('--no-interactive', action='store_true',
                        help='Skip interactive confirmation and download automatically')

    args = parser.parse_args()

    # Parse comma-separated arguments
    def parse_list_arg(arg):
        if arg is None:
            return None
        return [x.strip() for x in arg.split(',')]

    species = parse_list_arg(args.species)
    tissue = parse_list_arg(args.tissue)
    cell_type = parse_list_arg(args.cell_type)
    disease = parse_list_arg(args.disease)
    sex = parse_list_arg(args.sex)
    development_stage = parse_list_arg(args.development_stage)
    drug_treatment = parse_list_arg(args.drug_treatment)
    gene_set = parse_list_arg(args.genes)

    # Run query
    with CellxGeneQuery() as querier:
        result = querier.query_data(
            species=species,
            tissue=tissue,
            cell_type=cell_type,
            disease=disease,
            sex=sex,
            development_stage=development_stage,
            drug_treatment=drug_treatment,
            gene_set=gene_set,
        )

        if result is None:
            logger.error("Query failed")
            sys.exit(1)

        experiment, value_filter, var_filter, obs_df = result

        # Apply filters and get data
        logger.info("Fetching data from Census...")

        try:
            # Query the data
            with experiment.axis_query(
                measurement_name="RNA",
                obs_query=soma.AxisQuery(value_filter=value_filter) if value_filter else None,
                var_query=soma.AxisQuery(value_filter=var_filter) if var_filter else None,
            ) as query:
                # Read as AnnData
                logger.info("Reading data into AnnData format...")
                adata = query.to_anndata(X_name="raw")

                # Get filtered obs
                obs_df_filtered = query.obs().concat().to_pandas()

        except Exception as e:
            logger.error(f"Error during query execution: {e}")
            sys.exit(1)

        # Generate summary
        summary = querier.generate_summary(
            adata=adata,
            obs_df_filtered=obs_df_filtered,
            filters_applied=value_filter if value_filter else "None"
        )

        # Print summary
        querier.print_summary(summary)

        # Ask for download permission
        if not args.no_interactive:
            print("\n" + "=" * 80)
            print("Do you want to download and save this data? (yes/no): ", end='')
            response = input().strip().lower()

            if response not in ['yes', 'y']:
                logger.info("Download cancelled by user")
                sys.exit(0)

        # Save results
        querier.save_results(
            adata=adata,
            obs_df_filtered=obs_df_filtered,
            summary=summary,
            output_dir=Path(args.output_dir),
            prefix=args.output
        )


if __name__ == "__main__":
    main()
