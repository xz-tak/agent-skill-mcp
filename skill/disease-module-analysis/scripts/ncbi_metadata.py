#!/usr/bin/env python3
"""
NCBI Metadata Extraction Module

Extract metadata (species, tissue, cell type, disease/condition) from NCBI GEO datasets
"""

import re
import time
import logging
from typing import Dict, Optional
import requests
from xml.etree import ElementTree as ET


logger = logging.getLogger('ncbi_metadata')


class NCBIMetadataExtractor:
    """Extract metadata from NCBI GEO datasets."""

    def __init__(self, email='user@example.com', delay=0.4):
        """
        Initialize NCBI metadata extractor.

        Parameters
        ----------
        email : str
            Email for NCBI API (required by NCBI)
        delay : float
            Delay between API calls in seconds (NCBI limit: max 3 requests/sec)
        """
        self.email = email
        self.delay = delay
        self.base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        self.cache = {}

    def extract_geo_id(self, pert_name: str) -> Optional[str]:
        """
        Extract GEO series ID from perturbation name.

        Parameters
        ----------
        pert_name : str
            Perturbation name (e.g., 'GSE12345_GENE_UP')

        Returns
        -------
        geo_id : str or None
            GEO series ID (e.g., 'GSE12345') or None if not found
        """
        # Pattern for GEO series
        match = re.search(r'(GSE\d+)', pert_name)
        if match:
            return match.group(1)
        return None

    def query_geo_metadata(self, geo_id: str) -> Optional[Dict[str, str]]:
        """
        Query NCBI GEO for dataset metadata.

        Parameters
        ----------
        geo_id : str
            GEO series ID (e.g., 'GSE12345')

        Returns
        -------
        metadata : dict or None
            Dictionary with metadata fields or None if query failed
        """
        # Check cache
        if geo_id in self.cache:
            return self.cache[geo_id]

        try:
            # Search for GEO ID
            search_url = f"{self.base_url}esearch.fcgi"
            search_params = {
                'db': 'gds',
                'term': geo_id,
                'retmode': 'json',
                'email': self.email,
            }

            time.sleep(self.delay)
            response = requests.get(search_url, params=search_params, timeout=10)
            response.raise_for_status()
            search_result = response.json()

            # Get IDs
            id_list = search_result.get('esearchresult', {}).get('idlist', [])
            if not id_list:
                logger.warning(f"No results found for {geo_id}")
                return None

            # Fetch summary
            summary_url = f"{self.base_url}esummary.fcgi"
            summary_params = {
                'db': 'gds',
                'id': id_list[0],
                'retmode': 'json',
                'email': self.email,
            }

            time.sleep(self.delay)
            response = requests.get(summary_url, params=summary_params, timeout=10)
            response.raise_for_status()
            summary_result = response.json()

            # Extract metadata
            result = summary_result.get('result', {}).get(id_list[0], {})

            metadata = {
                'title': result.get('title', ''),
                'summary': result.get('summary', ''),
                'organism': result.get('taxon', ''),
                'sample_count': result.get('n_samples', ''),
                'platform': result.get('gpl', ''),
                'pubmed_id': result.get('pubmedids', [''])[0] if result.get('pubmedids') else '',
            }

            # Cache result
            self.cache[geo_id] = metadata
            return metadata

        except Exception as e:
            logger.warning(f"Error querying {geo_id}: {e}")
            return None

    def extract_metadata_fields(self, metadata: Dict[str, str]) -> Dict[str, str]:
        """
        Extract structured fields from metadata.

        Parameters
        ----------
        metadata : dict
            Raw metadata from NCBI

        Returns
        -------
        fields : dict
            Structured fields: species, tissue, cell_type, disease
        """
        if not metadata:
            return {
                'species': '',
                'tissue': '',
                'cell_type': '',
                'disease': '',
            }

        title = metadata.get('title', '').lower()
        summary = metadata.get('summary', '').lower()
        text = f"{title} {summary}"

        # Extract species
        species = self._extract_species(metadata.get('organism', ''), text)

        # Extract tissue/organ
        tissue = self._extract_tissue(text)

        # Extract cell type
        cell_type = self._extract_cell_type(text)

        # Extract disease/condition
        disease = self._extract_disease(text)

        return {
            'species': species,
            'tissue': tissue,
            'cell_type': cell_type,
            'disease': disease,
        }

    def _extract_species(self, organism: str, text: str) -> str:
        """Extract species from organism field and text."""
        # First check organism field
        if organism:
            organism_lower = organism.lower()
            if 'homo sapiens' in organism_lower or 'human' in organism_lower:
                return 'Human'
            elif 'mus musculus' in organism_lower or 'mouse' in organism_lower:
                return 'Mouse'
            elif 'rattus norvegicus' in organism_lower or 'rat' in organism_lower:
                return 'Rat'
            else:
                return organism

        # Check text
        if 'human' in text or 'homo sapiens' in text:
            return 'Human'
        elif 'mouse' in text or 'murine' in text or 'mus musculus' in text:
            return 'Mouse'
        elif 'rat' in text or 'rattus' in text:
            return 'Rat'

        return ''

    def _extract_tissue(self, text: str) -> str:
        """Extract tissue/organ from text."""
        tissues = {
            'liver': 'Liver',
            'hepat': 'Liver',
            'kidney': 'Kidney',
            'renal': 'Kidney',
            'heart': 'Heart',
            'cardiac': 'Heart',
            'brain': 'Brain',
            'cerebral': 'Brain',
            'lung': 'Lung',
            'pulmonary': 'Lung',
            'skin': 'Skin',
            'dermal': 'Skin',
            'blood': 'Blood',
            'muscle': 'Muscle',
            'bone': 'Bone',
            'intestin': 'Intestine',
            'colon': 'Colon',
            'stomach': 'Stomach',
            'pancrea': 'Pancreas',
            'spleen': 'Spleen',
            'thymus': 'Thymus',
            'adipose': 'Adipose',
            'breast': 'Breast',
            'prostate': 'Prostate',
            'ovary': 'Ovary',
            'testis': 'Testis',
            'uterus': 'Uterus',
            'placenta': 'Placenta',
            'eye': 'Eye',
            'retina': 'Retina',
        }

        found_tissues = []
        for keyword, tissue_name in tissues.items():
            if keyword in text:
                found_tissues.append(tissue_name)

        return '; '.join(list(set(found_tissues))) if found_tissues else ''

    def _extract_cell_type(self, text: str) -> str:
        """Extract cell type from text."""
        cell_types = {
            'fibroblast': 'Fibroblast',
            'macrophage': 'Macrophage',
            'monocyte': 'Monocyte',
            'neutrophil': 'Neutrophil',
            'lymphocyte': 'Lymphocyte',
            't cell': 'T cell',
            'b cell': 'B cell',
            'nk cell': 'NK cell',
            'dendritic': 'Dendritic cell',
            'endothel': 'Endothelial cell',
            'epithel': 'Epithelial cell',
            'keratinocyte': 'Keratinocyte',
            'hepatocyte': 'Hepatocyte',
            'neuron': 'Neuron',
            'astrocyte': 'Astrocyte',
            'microglia': 'Microglia',
            'stem cell': 'Stem cell',
            'cancer cell': 'Cancer cell',
            'tumor cell': 'Tumor cell',
            'cardiomyocyte': 'Cardiomyocyte',
            'adipocyte': 'Adipocyte',
            'osteoblast': 'Osteoblast',
            'osteoclast': 'Osteoclast',
            'myoblast': 'Myoblast',
            'pbmc': 'PBMC',
        }

        found_types = []
        for keyword, cell_name in cell_types.items():
            if keyword in text:
                found_types.append(cell_name)

        return '; '.join(list(set(found_types))) if found_types else ''

    def _extract_disease(self, text: str) -> str:
        """Extract disease/condition from text."""
        diseases = {
            'cancer': 'Cancer',
            'tumor': 'Cancer',
            'carcinoma': 'Carcinoma',
            'leukemia': 'Leukemia',
            'lymphoma': 'Lymphoma',
            'diabetes': 'Diabetes',
            'fibrosis': 'Fibrosis',
            'cirrhosis': 'Cirrhosis',
            'inflammation': 'Inflammation',
            'inflammatory': 'Inflammation',
            'infection': 'Infection',
            'viral': 'Viral infection',
            'bacterial': 'Bacterial infection',
            'alzheimer': 'Alzheimer',
            'parkinson': 'Parkinson',
            'atherosclerosis': 'Atherosclerosis',
            'hypertension': 'Hypertension',
            'obesity': 'Obesity',
            'asthma': 'Asthma',
            'copd': 'COPD',
            'arthritis': 'Arthritis',
            'lupus': 'Lupus',
            'psoriasis': 'Psoriasis',
            'wound': 'Wound',
            'injury': 'Injury',
            'ischemia': 'Ischemia',
            'hypoxia': 'Hypoxia',
            'sepsis': 'Sepsis',
            'autoimmune': 'Autoimmune',
            'metabolic': 'Metabolic disorder',
            'healthy': 'Healthy',
            'normal': 'Normal',
            'control': 'Control',
        }

        found_diseases = []
        for keyword, disease_name in diseases.items():
            if keyword in text:
                found_diseases.append(disease_name)

        return '; '.join(list(set(found_diseases))) if found_diseases else ''

    def enrich_dataframe(self, df, name_column='pert_GENE'):
        """
        Enrich dataframe with metadata columns.

        Parameters
        ----------
        df : DataFrame
            Input dataframe with perturbation names
        name_column : str
            Column containing perturbation names

        Returns
        -------
        df_enriched : DataFrame
            Dataframe with added metadata columns
        """
        import pandas as pd

        df = df.copy()

        # Initialize columns
        df['species'] = ''
        df['tissue'] = ''
        df['cell_type'] = ''
        df['disease'] = ''
        df['geo_id'] = ''

        # Process each row
        for idx, row in df.iterrows():
            pert_name = row[name_column]

            # Extract GEO ID
            geo_id = self.extract_geo_id(pert_name)
            if not geo_id:
                continue

            df.at[idx, 'geo_id'] = geo_id

            # Query metadata
            metadata = self.query_geo_metadata(geo_id)
            if not metadata:
                continue

            # Extract fields
            fields = self.extract_metadata_fields(metadata)
            df.at[idx, 'species'] = fields['species']
            df.at[idx, 'tissue'] = fields['tissue']
            df.at[idx, 'cell_type'] = fields['cell_type']
            df.at[idx, 'disease'] = fields['disease']

        return df


def enrich_perturbation_tables(df, pert_type='gene', email='user@example.com'):
    """
    Enrich perturbation tables with NCBI metadata.

    Parameters
    ----------
    df : DataFrame
        Perturbation results dataframe
    pert_type : str
        Type of perturbation: 'gene', 'tf', 'crispr', 'chem'
    email : str
        Email for NCBI API

    Returns
    -------
    df_enriched : DataFrame
        Enriched dataframe with metadata columns
    """
    # Determine column name based on perturbation type
    column_map = {
        'gene': 'pert_GENE',
        'tf': 'pert_TF',
        'crispr': 'pert_GENE',
        'chem': 'pert_CHEM',
    }

    name_column = column_map.get(pert_type.lower())
    if not name_column:
        raise ValueError(f"Unknown perturbation type: {pert_type}")

    if name_column not in df.columns:
        logger.warning(f"Column {name_column} not found in dataframe")
        return df

    # Only enrich if it's a GEO-based perturbation (not for CHEM)
    if pert_type.lower() == 'chem':
        logger.info("Skipping metadata enrichment for chemical perturbations")
        return df

    logger.info(f"Enriching {len(df)} {pert_type} perturbations with NCBI metadata...")

    extractor = NCBIMetadataExtractor(email=email)
    df_enriched = extractor.enrich_dataframe(df, name_column=name_column)

    # Count how many were enriched
    n_enriched = (df_enriched['geo_id'] != '').sum()
    logger.info(f"Successfully enriched {n_enriched}/{len(df)} perturbations with metadata")

    return df_enriched


if __name__ == '__main__':
    # Test the extractor
    import pandas as pd

    # Create test dataframe
    test_data = {
        'pert_GENE': [
            'GSE12345_GENE1_UP',
            'GSE67890_GENE2_DOWN',
            'MANUAL_GENE3',
        ],
        'p_value': [0.001, 0.01, 0.05],
    }
    df = pd.DataFrame(test_data)

    print("Original dataframe:")
    print(df)

    # Enrich
    df_enriched = enrich_perturbation_tables(df, pert_type='gene', email='test@example.com')

    print("\nEnriched dataframe:")
    print(df_enriched)
