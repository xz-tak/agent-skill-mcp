"""
NCBI Paper Query - Omics Extractor Module

Extracts omics dataset accessions from publications and fetches detailed metadata
from GEO, SRA, ArrayExpress, PRIDE, MetaboLights, and other databases.
"""

import re
import time
import logging
from typing import List, Dict, Optional, TYPE_CHECKING

import requests

from ..config import logger
from ..models import OmicsDataset

if TYPE_CHECKING:
    from .downloader import PaperDownloader

# PDF parsing
try:
    from PyPDF2 import PdfReader
except ImportError:
    PdfReader = None


class OmicsExtractor:
    """Extract omics dataset accessions from publications and fetch detailed metadata."""

    # Accession patterns for various databases
    ACCESSION_PATTERNS = {
        "GEO": [
            (r"\bGSE\d{3,8}\b", "series"),
            (r"\bGSM\d{3,8}\b", "sample"),
            (r"\bGPL\d{3,6}\b", "platform"),
        ],
        "SRA": [
            (r"\bSRP\d{5,9}\b", "study"),
            (r"\bSRR\d{6,10}\b", "run"),
            (r"\bSRX\d{5,9}\b", "experiment"),
            (r"\bSRS\d{5,9}\b", "sample"),
            (r"\bPRJNA\d{5,9}\b", "bioproject"),
            (r"\bPRJEB\d{5,9}\b", "bioproject"),
        ],
        "ArrayExpress": [
            (r"\bE-[A-Z]{4}-\d{3,6}\b", "experiment"),
        ],
        "PRIDE": [
            (r"\bPXD\d{5,8}\b", "project"),
        ],
        "MetaboLights": [
            (r"\bMTBLS\d{3,6}\b", "study"),
        ],
        "dbGaP": [
            (r"\bphs\d{6}\b", "study"),
        ],
        "ENA": [
            (r"\bERP\d{5,9}\b", "study"),
            (r"\bERR\d{6,10}\b", "run"),
        ],
        "SingleCellPortal": [
            (r"\bSCP\d{3,6}\b", "study"),
        ],
        "NGDC_GSA": [
            (r"\bCRA\d{6,9}\b", "project"),      # GSA project (e.g., CRA001234)
            (r"\bHRA\d{6,9}\b", "project"),      # GSA-Human project (e.g., HRA001234)
            (r"\bPRJCA\d{6,9}\b", "bioproject"), # BioProject China (e.g., PRJCA001234)
            (r"\bSAMC\d{6,9}\b", "sample"),      # Sample (e.g., SAMC001234)
            (r"\bCRR\d{6,10}\b", "run"),         # Run (e.g., CRR001234)
            (r"\bCRX\d{6,9}\b", "experiment"),   # Experiment (e.g., CRX001234)
        ],
    }

    # Database URLs
    DATABASE_URLS = {
        "GEO": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
        "SRA": "https://www.ncbi.nlm.nih.gov/sra/",
        "ArrayExpress": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/",
        "PRIDE": "https://www.ebi.ac.uk/pride/archive/projects/",
        "MetaboLights": "https://www.ebi.ac.uk/metabolights/",
        "dbGaP": "https://www.ncbi.nlm.nih.gov/gap/?term=",
        "ENA": "https://www.ebi.ac.uk/ena/browser/view/",
        "SingleCellPortal": "https://singlecell.broadinstitute.org/single_cell/study/",
        "NGDC_GSA": "https://ngdc.cncb.ac.cn/gsa/browse/",
    }

    def __init__(self):
        self.geo_api_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.sra_api_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.arrayexpress_api_url = "https://www.ebi.ac.uk/biostudies/api/v1/studies/"
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mozilla/5.0 (NCBI Paper Retriever)"
        })

    def extract_from_text(self, text: str, source: str = "unknown") -> List[OmicsDataset]:
        """Extract omics dataset accessions from text.

        Args:
            text: The text to search for accessions.
            source: Where the text came from ("data_availability", "key_resources",
                   "full_text", "abstract"). Used for tracking.

        Returns:
            List of OmicsDataset objects with source field populated.
        """
        datasets = []
        seen_accessions = set()

        for database, patterns in self.ACCESSION_PATTERNS.items():
            for pattern, data_type in patterns:
                matches = re.findall(pattern, text, re.IGNORECASE)
                for accession in matches:
                    accession_upper = accession.upper()
                    if accession_upper not in seen_accessions:
                        seen_accessions.add(accession_upper)
                        datasets.append(OmicsDataset(
                            accession=accession_upper,
                            database=database,
                            technology=data_type,
                            data_url=self.DATABASE_URLS.get(database, "") + accession_upper,
                            source=source
                        ))

        return datasets

    def _extract_data_availability_section(self, text: str) -> Optional[str]:
        """Extract the Data and Code Availability section from paper text.

        Common section names in various journals:
        - "Data Availability" (Nature)
        - "Data and Code Availability" (Cell)
        - "Data Availability Statement" (PLOS)
        - "Code and Data Availability"
        - "Availability of data and materials"
        - "Data Access"
        - "Accession Numbers"
        - "Data deposition"

        Args:
            text: Full text of the paper.

        Returns:
            The section text if found, None otherwise.
        """
        # Section header patterns (order by specificity)
        header_patterns = [
            r'Data\s+and\s+Code\s+Availability',
            r'Code\s+and\s+Data\s+Availability',
            r'Data\s+(?:and\s+Materials?\s+)?Availability(?:\s+Statement)?',
            r'Availability\s+of\s+(?:data|materials)(?:\s+and\s+materials)?',
            r'Data\s+Access(?:ibility)?',
            r'Accession\s+(?:Numbers?|Codes?)',
            r'Data\s+[Dd]eposition',
            r'Data\s+[Ss]haring',
            r'Resource\s+Availability',
        ]

        # Patterns for next section headers (to detect section boundaries)
        next_section_patterns = [
            r'(?:References|Acknowledgment|Author\s+Contributions?|Competing|Conflict|Declaration)',
            r'(?:Supplementar|Methods|Figure|Table|Supporting\s+Information)',
            r'(?:Funding|Ethics|Financial)',
        ]
        next_section = '|'.join(next_section_patterns)

        # Try each pattern
        for pattern in header_patterns:
            # Match section header to next section or end of document
            # Look for the header followed by content, stopping at next section
            regex = rf'(?:^|\n)\s*[\*\#]*\s*({pattern})[:\.\s\*]*\n?(.*?)(?=\n\s*(?:{next_section})|\Z)'
            match = re.search(regex, text, re.IGNORECASE | re.DOTALL)
            if match:
                section_text = match.group(2).strip()
                # Minimum meaningful content (avoid empty sections)
                if len(section_text) > 20:
                    # Limit section size to avoid capturing too much
                    if len(section_text) > 3000:
                        section_text = section_text[:3000]
                    logger.debug(f"Found Data Availability section ({len(section_text)} chars) via pattern: {pattern}")
                    return section_text

        # Alternative approach: look for paragraphs containing key phrases
        # This helps when the section header isn't clearly separated
        key_phrases = [
            r'(?:data|datasets?)\s+(?:have\s+been\s+)?deposited\s+(?:in|at|to)',
            r'(?:accession\s+(?:number|code|ID)|GEO\s+accession)',
            r'(?:available\s+(?:at|in|from|through)\s+(?:GEO|SRA|ArrayExpress|PRIDE))',
            r'(?:publicly\s+available\s+(?:at|in|from))',
        ]

        for phrase_pattern in key_phrases:
            # Find sentences containing the phrase
            regex = rf'([^.]*{phrase_pattern}[^.]*\.(?:[^.]*(?:GSE|SRP|PRJNA|E-[A-Z]{{4}}-)[^.]*\.)?)'
            matches = re.findall(regex, text, re.IGNORECASE)
            if matches:
                combined = ' '.join(matches)
                if len(combined) > 30:
                    logger.debug(f"Found Data Availability content via phrase search ({len(combined)} chars)")
                    return combined[:3000]

        return None

    def _extract_key_resources_table(self, text: str) -> Optional[str]:
        """Extract Key Resources Table section from paper text.

        Cell/Elsevier papers often have structured tables like:
        | REAGENT or RESOURCE | SOURCE | IDENTIFIER |
        | Deposited Data      |        |            |
        | RNA-seq data        | GEO    | GSE123456  |

        Args:
            text: Full text of the paper.

        Returns:
            Text containing database/accession info from table, or None.
        """
        # Patterns for Key Resources Table headers
        table_patterns = [
            r'KEY\s+RESOURCES?\s+TABLE',
            r'REAGENT\s+(?:or|and)\s+RESOURCE',
            r'STAR\s+METHODS',
        ]

        # Look for the section
        for pattern in table_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                # Extract text after the header
                start_pos = match.end()
                # Find reasonable end (next major section or limited length)
                section_text = text[start_pos:start_pos + 5000]

                # Look for "Deposited Data" subsection
                deposited_match = re.search(
                    r'Deposited\s+[Dd]ata.*?(?=\n\s*(?:Software|Experimental|Oligonucleotides|Chemicals|Other)|$)',
                    section_text,
                    re.IGNORECASE | re.DOTALL
                )
                if deposited_match:
                    deposited_text = deposited_match.group(0)
                    logger.debug(f"Found Deposited Data in Key Resources Table ({len(deposited_text)} chars)")
                    return deposited_text

                # If no specific "Deposited Data", search for accession patterns in table
                accession_patterns = r'(GSE\d+|SRP\d+|PRJNA\d+|E-[A-Z]{4}-\d+|PXD\d+)'
                accession_matches = re.findall(accession_patterns, section_text, re.IGNORECASE)
                if accession_matches:
                    # Return the section with accessions
                    logger.debug(f"Found {len(accession_matches)} accessions in Key Resources Table")
                    return section_text[:2000]

        return None

    def extract_from_pdf(self, pdf_path: str) -> List[OmicsDataset]:
        """Extract accessions from a PDF file, prioritizing Data Availability section.

        Extraction priority:
        1. Data Availability section (most reliable)
        2. Key Resources Table (common in Cell/Elsevier)
        3. Full text search (fallback)

        Args:
            pdf_path: Path to the PDF file.

        Returns:
            List of OmicsDataset objects with source field populated.
        """
        if PdfReader is None:
            logger.warning("PyPDF2 not available for PDF parsing")
            return []

        try:
            reader = PdfReader(pdf_path)
            full_text = ""
            for page in reader.pages:
                full_text += page.extract_text() + "\n"

            # Step 1: Try Data Availability section first (most reliable source)
            data_avail_section = self._extract_data_availability_section(full_text)
            if data_avail_section:
                datasets = self.extract_from_text(data_avail_section, source="data_availability")
                if datasets:
                    logger.info(f"Extracted {len(datasets)} accessions from Data Availability section")
                    return datasets
                else:
                    logger.debug("No accessions found in Data Availability section")

            # Step 2: Try Key Resources Table (common in Cell/Elsevier)
            key_resources_section = self._extract_key_resources_table(full_text)
            if key_resources_section:
                datasets = self.extract_from_text(key_resources_section, source="key_resources")
                if datasets:
                    logger.info(f"Extracted {len(datasets)} accessions from Key Resources Table")
                    return datasets
                else:
                    logger.debug("No accessions found in Key Resources Table")

            # Step 3: Fallback to full text search
            logger.debug("Falling back to full PDF text search for accessions")
            datasets = self.extract_from_text(full_text, source="full_text")
            if datasets:
                logger.info(f"Extracted {len(datasets)} accessions from full text search (fallback)")
            return datasets

        except Exception as e:
            logger.error(f"Failed to parse PDF: {e}")
            return []

    def extract_from_web(
        self,
        pmc_id: Optional[str],
        doi: Optional[str],
        downloader: 'PaperDownloader'
    ) -> List[OmicsDataset]:
        """
        Extract accessions from web page without downloading PDF.

        Uses Playwright to fetch HTML content (PMC first, then publisher),
        then applies existing extraction logic to find omics accessions.

        Args:
            pmc_id: PubMed Central ID (e.g., "PMC1234567")
            doi: Digital Object Identifier
            downloader: PaperDownloader instance for web fetching

        Returns:
            List of OmicsDataset objects with source field populated.
        """
        # Fetch HTML content (PMC first, then DOI)
        article_text = downloader.fetch_article_html(pmc_id, doi)

        if not article_text:
            logger.debug(f"Could not fetch article content (PMC: {pmc_id}, DOI: {doi})")
            return []

        logger.debug(f"Fetched article HTML ({len(article_text)} chars)")

        # Use shared extraction logic
        return self.extract_from_fetched_text(article_text)

    def extract_from_fetched_text(self, article_text: str) -> List[OmicsDataset]:
        """
        Extract accessions from pre-fetched article text without re-fetching.

        This method allows extraction from text that was already fetched,
        decoupling the web access check from accession extraction.

        Args:
            article_text: Pre-fetched article text content

        Returns:
            List of OmicsDataset objects with source field populated.
        """
        if not article_text:
            return []

        # Reuse existing extraction logic (same as extract_from_pdf)
        # Step 1: Try Data Availability section
        data_avail_section = self._extract_data_availability_section(article_text)
        if data_avail_section:
            datasets = self.extract_from_text(data_avail_section, source="data_availability")
            if datasets:
                logger.info(f"Extracted {len(datasets)} accessions from web Data Availability section")
                return datasets

        # Step 2: Try Key Resources Table
        key_resources_section = self._extract_key_resources_table(article_text)
        if key_resources_section:
            datasets = self.extract_from_text(key_resources_section, source="key_resources")
            if datasets:
                logger.info(f"Extracted {len(datasets)} accessions from web Key Resources Table")
                return datasets

        # Step 3: Fallback to full text
        datasets = self.extract_from_text(article_text, source="web_full_text")
        if datasets:
            logger.info(f"Extracted {len(datasets)} accessions from web full text (fallback)")
        return datasets

    def _parse_sample_conditions(self, accession: str) -> str:
        """Parse GEO sample characteristics to count samples per condition.

        Fetches sample-level metadata from GEO and groups samples by key
        characteristics like treatment, diagnosis, tissue, etc.

        Args:
            accession: GEO series accession (e.g., GSE12345)

        Returns:
            Formatted string like "Inflamed (6); Normal (12); Treated (4)"
        """
        try:
            # Fetch all samples for this series in text format
            url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=gsm&form=text&view=brief"
            response = self.session.get(url, timeout=30)
            if response.status_code != 200:
                return ""

            text = response.text

            # Priority fields for grouping samples (in order of preference)
            group_fields = [
                'treatment', 'treatment_induction', 'treatment_maintenance',
                'diagnosis', 'disease', 'condition', 'disease state', 'disease_state',
                'tissue', 'cell type', 'cell_type',
                'genotype', 'phenotype', 'group', 'sample_group', 'sample group',
                'status', 'state', 'response', 'responder'
            ]

            # Parse sample blocks and their characteristics
            sample_blocks = text.split('^SAMPLE = ')
            condition_counter = {}

            for block in sample_blocks[1:]:  # Skip first empty block
                lines = block.split('\n')

                # Collect all characteristics for this sample
                sample_chars = {}
                for line in lines:
                    if '!Sample_characteristics_ch' in line:
                        # Parse "!Sample_characteristics_ch1 = key: value"
                        match = re.search(r'=\s*([^:]+):\s*(.+)', line)
                        if match:
                            key = match.group(1).strip().lower()
                            value = match.group(2).strip()
                            # Skip empty or placeholder values
                            if value and value.lower() not in ('none', 'n/a', 'na', '.', '--'):
                                sample_chars[key] = value

                # Find the best grouping field for this sample
                for field in group_fields:
                    if field in sample_chars:
                        condition_value = sample_chars[field]
                        # Truncate very long condition names
                        if len(condition_value) > 50:
                            condition_value = condition_value[:47] + "..."
                        condition_counter[condition_value] = condition_counter.get(condition_value, 0) + 1
                        break

            # Format as "Condition (N); ..."
            if condition_counter:
                # Sort by count descending, then alphabetically
                sorted_conditions = sorted(condition_counter.items(), key=lambda x: (-x[1], x[0]))
                formatted = "; ".join(f"{cond} ({count})" for cond, count in sorted_conditions)
                return formatted

            return ""

        except Exception as e:
            logger.debug(f"Sample condition parsing failed for {accession}: {e}")
            return ""

    def fetch_geo_metadata(self, accession: str) -> Dict:
        """Fetch detailed metadata for a GEO accession."""
        try:
            # First try to get basic info via Entrez
            url = f"{self.geo_api_url}esearch.fcgi"
            params = {"db": "gds", "term": accession, "retmode": "json"}

            response = self.session.get(url, params=params, timeout=10)
            data = response.json()

            ids = data.get("esearchresult", {}).get("idlist", [])
            result = {}

            if ids:
                # Fetch summary
                summary_url = f"{self.geo_api_url}esummary.fcgi"
                summary_params = {"db": "gds", "id": ids[0], "retmode": "json"}
                summary_response = self.session.get(summary_url, params=summary_params, timeout=10)
                summary_data = summary_response.json()

                result_data = summary_data.get("result", {})
                if ids[0] in result_data:
                    record = result_data[ids[0]]
                    result = {
                        "platform": record.get("gpl", ""),
                        "sample_count": record.get("n_samples", 0),
                        "technology": record.get("gdstype", ""),
                        "title": record.get("title", ""),
                        "organism": record.get("taxon", ""),
                    }

            # Fetch detailed info from GEO directly for conditions
            geo_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=self&form=text&view=brief"
            geo_response = self.session.get(geo_url, timeout=15)

            if geo_response.status_code == 200:
                geo_text = geo_response.text

                # Extract summary/design for conditions
                design_match = re.search(r'!Series_overall_design\s*=\s*(.+?)(?=\n!|$)', geo_text, re.DOTALL)
                if design_match:
                    result["experiment_design"] = design_match.group(1).strip()[:500]

                # Extract sample characteristics for conditions
                samples = re.findall(r'!Series_sample_id\s*=\s*(\S+)', geo_text)
                if samples and not result.get("sample_count"):
                    result["sample_count"] = len(samples)

                # Extract type for assay description
                type_match = re.search(r'!Series_type\s*=\s*(.+?)(?=\n|$)', geo_text)
                if type_match:
                    result["assay_description"] = type_match.group(1).strip()

                # Extract conditions from summary
                summary_match = re.search(r'!Series_summary\s*=\s*(.+?)(?=\n!|$)', geo_text, re.DOTALL)
                if summary_match:
                    summary = summary_match.group(1).strip()
                    # Look for condition-related text
                    cond_patterns = [
                        r'(\d+)\s*(control|treated|patient|healthy|disease|sample)',
                        r'(control|treated|patient|healthy|disease)\s*(?:group|sample)?\s*(?:n\s*=\s*)?(\d+)',
                    ]
                    conditions = []
                    for pattern in cond_patterns:
                        matches = re.findall(pattern, summary, re.IGNORECASE)
                        for m in matches:
                            if isinstance(m, tuple):
                                conditions.append(" ".join(str(x) for x in m))
                    if conditions:
                        result["conditions"] = "; ".join(conditions[:5])

            # Parse sample-level conditions for detailed condition counts
            condition_counts = self._parse_sample_conditions(accession)
            if condition_counts:
                result["condition_counts"] = condition_counts

            return result

        except Exception as e:
            logger.debug(f"GEO metadata fetch failed for {accession}: {e}")
            return {}

    def fetch_sra_metadata(self, accession: str) -> Dict:
        """Fetch detailed metadata for an SRA accession."""
        try:
            # Determine database based on accession type
            if accession.startswith("PRJNA") or accession.startswith("PRJEB"):
                db = "bioproject"
            else:
                db = "sra"

            # Search for the accession
            url = f"{self.sra_api_url}esearch.fcgi"
            params = {"db": db, "term": accession, "retmode": "json"}

            response = self.session.get(url, params=params, timeout=10)
            data = response.json()

            ids = data.get("esearchresult", {}).get("idlist", [])
            result = {}

            if ids:
                # Fetch summary
                summary_url = f"{self.sra_api_url}esummary.fcgi"
                summary_params = {"db": db, "id": ids[0], "retmode": "json"}
                summary_response = self.session.get(summary_url, params=summary_params, timeout=10)
                summary_data = summary_response.json()

                result_data = summary_data.get("result", {})
                if ids[0] in result_data:
                    record = result_data[ids[0]]

                    if db == "bioproject":
                        result = {
                            "title": record.get("project_name", ""),
                            "organism": record.get("organism_name", ""),
                            "experiment_design": record.get("project_description", "")[:500] if record.get("project_description") else "",
                        }
                    else:
                        # SRA record
                        exp_xml = record.get("expxml", "")
                        runs = record.get("runs", "")

                        # Parse run count from runs field
                        run_count = len(re.findall(r'Run acc="[^"]+"', runs)) if runs else 0
                        result["sample_count"] = run_count

                        # Parse library strategy
                        lib_match = re.search(r'<Library_Strategy>([^<]+)</Library_Strategy>', exp_xml)
                        if lib_match:
                            result["assay_description"] = lib_match.group(1)

                        # Parse organism
                        org_match = re.search(r'<Organism>([^<]+)</Organism>', exp_xml)
                        if org_match:
                            result["organism"] = org_match.group(1)

                        # Parse title
                        title_match = re.search(r'<Title>([^<]+)</Title>', exp_xml)
                        if title_match:
                            result["title"] = title_match.group(1)

            return result

        except Exception as e:
            logger.debug(f"SRA metadata fetch failed for {accession}: {e}")
            return {}

    def fetch_arrayexpress_metadata(self, accession: str) -> Dict:
        """Fetch detailed metadata for an ArrayExpress accession."""
        try:
            # Try BioStudies API (ArrayExpress was migrated to BioStudies)
            url = f"{self.arrayexpress_api_url}{accession}"
            response = self.session.get(url, timeout=15)

            result = {}

            if response.status_code == 200:
                data = response.json()

                result = {
                    "title": data.get("title", ""),
                    "organism": "",
                    "sample_count": 0,
                    "experiment_design": "",
                    "assay_description": "",
                }

                # Extract attributes
                attributes = data.get("attributes", [])
                for attr in attributes:
                    name = attr.get("name", "").lower()
                    value = attr.get("value", "")

                    if "organism" in name:
                        result["organism"] = value
                    elif "assay" in name or "technology" in name:
                        result["assay_description"] = value
                    elif "sample" in name and "count" in name:
                        try:
                            result["sample_count"] = int(value)
                        except ValueError:
                            pass

                # Get section info
                sections = data.get("section", {})
                if sections:
                    subsections = sections.get("subsections", [])
                    for subsec in subsections:
                        if subsec.get("type") == "Study":
                            for attr in subsec.get("attributes", []):
                                if attr.get("name") == "Description":
                                    result["experiment_design"] = attr.get("value", "")[:500]

            return result

        except Exception as e:
            logger.debug(f"ArrayExpress metadata fetch failed for {accession}: {e}")
            return {}

    def enrich_datasets(self, datasets: List[OmicsDataset]) -> List[OmicsDataset]:
        """Enrich datasets with metadata from their databases."""
        for dataset in datasets:
            metadata = {}

            # Fetch metadata based on database type
            if dataset.database == "GEO" and dataset.accession.startswith("GSE"):
                metadata = self.fetch_geo_metadata(dataset.accession)
            elif dataset.database == "SRA":
                metadata = self.fetch_sra_metadata(dataset.accession)
            elif dataset.database == "ArrayExpress":
                metadata = self.fetch_arrayexpress_metadata(dataset.accession)

            # Apply metadata to dataset
            if metadata:
                if metadata.get("platform"):
                    dataset.platform = metadata["platform"]
                if metadata.get("sample_count"):
                    dataset.sample_count = metadata["sample_count"]
                if metadata.get("technology"):
                    dataset.technology = metadata["technology"]
                if metadata.get("title"):
                    dataset.title = metadata["title"]
                if metadata.get("organism"):
                    dataset.organism = metadata["organism"]
                if metadata.get("conditions"):
                    dataset.conditions = metadata["conditions"]
                if metadata.get("condition_counts"):
                    dataset.condition_counts = metadata["condition_counts"]
                if metadata.get("assay_description"):
                    dataset.assay_description = metadata["assay_description"]
                if metadata.get("experiment_design"):
                    dataset.experiment_design = metadata["experiment_design"]

            # Rate limit to avoid overloading APIs
            time.sleep(0.5)

        return datasets
