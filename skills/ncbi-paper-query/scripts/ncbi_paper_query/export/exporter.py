"""
NCBI Paper Query - Results Exporter Module

Exports results to CSV/Excel with proper formatting and generates
separate files for accessions and failed studies.
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import List, Optional

import pandas as pd

from ..config import logger
from ..models import Publication


class ResultsExporter:
    """Export results to CSV/Excel."""

    def __init__(self, output_dir: str = "output", study_name: str = None):
        """
        Initialize exporter with output directory structure.

        Args:
            output_dir: Base output directory (default: "output")
            study_name: If provided, creates output_dir/study_name/ subdirectory structure:
                       - output/{study_name}/results.csv
                       - output/{study_name}/downloads/
                       - output/{study_name}/accessions.csv
        """
        self.base_output_dir = Path(output_dir)
        self.study_name = study_name
        self._dirs_created = False

        if study_name:
            self.output_dir = self.base_output_dir / study_name
            self.downloads_dir = self.output_dir / "downloads"
        else:
            self.output_dir = self.base_output_dir
            self.downloads_dir = self.output_dir / "downloads"

    def _ensure_dirs(self):
        """Create output directories lazily (only when needed)."""
        if not self._dirs_created:
            self.output_dir.mkdir(parents=True, exist_ok=True)
            self.downloads_dir.mkdir(parents=True, exist_ok=True)
            self._dirs_created = True

    def publications_to_dataframe(self, publications: List[Publication]) -> pd.DataFrame:
        """Convert publications to a DataFrame."""
        rows = []

        for pub in publications:
            # Determine download status
            if pub.free_access:
                download_status = "Free"
            elif pub.download_path:
                download_status = "Downloaded"
            elif pub.web_scraped:
                download_status = "Full access without download"
            else:
                download_status = "Not Downloaded"

            # Base publication info
            row = {
                "PMID": pub.pmid,
                "DOI": pub.doi,
                "PMC_ID": pub.pmc_id,
                "Title": pub.title,
                "Abstract": pub.abstract[:500] + "..." if len(pub.abstract) > 500 else pub.abstract,
                "Authors": "; ".join(pub.authors[:5]) + ("..." if len(pub.authors) > 5 else ""),
                "Affiliations": "; ".join(pub.affiliations[:3]) + ("..." if len(pub.affiliations) > 3 else ""),
                "Journal": pub.journal,
                "Impact_Factor": pub.impact_factor,
                "Publication_Year": pub.publication_year,
                "Is_Preprint": pub.is_preprint,
                "Matched_Disease": pub.matched_disease,
                "Matched_Tissue": pub.matched_tissue,
                "Matched_Organism": pub.matched_organism,
                "Matched_Cell_Type": pub.matched_cell_type,
                "Relevance_Score": pub.relevance_score,
                "Free_Access": pub.free_access,
                "Requires_Subscription": pub.requires_subscription,
                "Access_URL": pub.access_url,
                "Download_Path": pub.download_path,
                "Download_Status": download_status,
            }

            # Omics datasets
            if pub.omics_datasets:
                accessions = [d.accession for d in pub.omics_datasets]
                databases = list(set(d.database for d in pub.omics_datasets))
                technologies = list(set(d.technology for d in pub.omics_datasets if d.technology))
                urls = [d.data_url for d in pub.omics_datasets if d.data_url]
                sources = [d.source for d in pub.omics_datasets]

                # Collect sample counts, condition_counts, and assay descriptions
                sample_counts = [str(d.sample_count) for d in pub.omics_datasets if d.sample_count]
                condition_counts = [d.condition_counts for d in pub.omics_datasets if d.condition_counts]
                assay_descriptions = [d.assay_description for d in pub.omics_datasets if d.assay_description]
                organisms = list(set(d.organism for d in pub.omics_datasets if d.organism))

                row["Omics_Accessions"] = "; ".join(accessions)
                row["Omics_Sources"] = "; ".join(sources)
                row["Omics_Databases"] = "; ".join(databases)
                row["Omics_Technologies"] = "; ".join(technologies)
                row["Omics_URLs"] = "; ".join(urls[:5])
                row["Omics_Count"] = len(pub.omics_datasets)
                row["Omics_Sample_Count"] = "; ".join(sample_counts) if sample_counts else ""
                row["Omics_Condition_Counts"] = " | ".join(condition_counts) if condition_counts else ""
                row["Omics_Assay_Description"] = "; ".join(list(set(assay_descriptions))[:3]) if assay_descriptions else ""
                row["Omics_Organisms"] = "; ".join(organisms) if organisms else ""
            else:
                row["Omics_Accessions"] = ""
                row["Omics_Sources"] = ""
                row["Omics_Databases"] = ""
                row["Omics_Technologies"] = ""
                row["Omics_URLs"] = ""
                row["Omics_Count"] = 0
                row["Omics_Sample_Count"] = ""
                row["Omics_Condition_Counts"] = ""
                row["Omics_Assay_Description"] = ""
                row["Omics_Organisms"] = ""

            rows.append(row)

        return pd.DataFrame(rows)

    def export_csv(self, publications: List[Publication], filename: str = None) -> str:
        """Export to CSV."""
        self._ensure_dirs()

        if filename is None:
            if self.study_name:
                filename = "results.csv"
            else:
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                filename = f"results_{timestamp}.csv"

        filepath = self.output_dir / filename
        df = self.publications_to_dataframe(publications)
        df.to_csv(filepath, index=False)

        logger.info(f"Exported {len(publications)} publications to {filepath}")
        return str(filepath)

    def export_excel(self, publications: List[Publication], filename: str = None) -> str:
        """Export to Excel with formatting."""
        self._ensure_dirs()

        if filename is None:
            if self.study_name:
                filename = "results.xlsx"
            else:
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                filename = f"results_{timestamp}.xlsx"

        filepath = self.output_dir / filename
        df = self.publications_to_dataframe(publications)

        with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Publications', index=False)

            # Auto-adjust column widths
            worksheet = writer.sheets['Publications']
            from openpyxl.utils import get_column_letter
            for idx, col in enumerate(df.columns):
                max_length = max(
                    df[col].astype(str).map(len).max(),
                    len(col)
                )
                col_letter = get_column_letter(idx + 1)  # openpyxl uses 1-based indexing
                worksheet.column_dimensions[col_letter].width = min(max_length + 2, 50)

        logger.info(f"Exported {len(publications)} publications to {filepath}")
        return str(filepath)

    def export_accessions(self, publications: List[Publication], filename: str = "accessions.csv") -> Optional[str]:
        """Export omics accessions to a separate CSV file.

        Creates a flat table with one row per accession, including:
        - PMID, Title, Journal of the source paper
        - Accession, Database, Technology, URL, Sample count
        """
        rows = []

        for pub in publications:
            if not pub.omics_datasets:
                continue

            for dataset in pub.omics_datasets:
                row = {
                    "PMID": pub.pmid,
                    "Title": pub.title[:100] + "..." if len(pub.title) > 100 else pub.title,
                    "Journal": pub.journal,
                    "Publication_Year": pub.publication_year,
                    "Accession": dataset.accession,
                    "Database": dataset.database,
                    "Technology": dataset.technology or "",
                    "Organism": dataset.organism or "",
                    "Sample_Count": dataset.sample_count or "",
                    "Condition_Counts": dataset.condition_counts or "",
                    "Assay_Description": dataset.assay_description or "",
                    "Data_URL": dataset.data_url or "",
                    "Source": dataset.source or "",
                }
                rows.append(row)

        if not rows:
            logger.info("No omics accessions to export")
            return None

        self._ensure_dirs()
        filepath = self.output_dir / filename
        df = pd.DataFrame(rows)
        df.to_csv(filepath, index=False)

        logger.info(f"Exported {len(rows)} omics accessions to {filepath}")
        return str(filepath)

    def export_failed_studies(
        self,
        publications: List[Publication],
        filename: str = "failed_studies.csv"
    ) -> Optional[str]:
        """
        Export papers that couldn't be accessed by any method to a separate CSV.

        A paper is considered "failed" if:
        - free_access = False
        - download_path = None
        - web_scraped = False

        This helps identify papers that may need manual retrieval or
        institutional access.

        Args:
            publications: List of Publication objects
            filename: Output filename (default: failed_studies.csv)

        Returns:
            Path to the exported file, or None if no failed studies.
        """
        rows = []

        for pub in publications:
            # Check if this is a failed study (couldn't be accessed)
            if not pub.free_access and not pub.download_path and not pub.web_scraped:
                row = {
                    "PMID": pub.pmid,
                    "DOI": pub.doi,
                    "PMC_ID": pub.pmc_id,
                    "Title": pub.title,
                    "Journal": pub.journal,
                    "Impact_Factor": pub.impact_factor,
                    "Publication_Year": pub.publication_year,
                    "Authors": "; ".join(pub.authors[:3]) + ("..." if len(pub.authors) > 3 else ""),
                    "Matched_Disease": pub.matched_disease,
                    "Matched_Tissue": pub.matched_tissue,
                    "Matched_Organism": pub.matched_organism,
                    "Relevance_Score": pub.relevance_score,
                    "Requires_Subscription": pub.requires_subscription,
                    "Access_URL": pub.access_url or f"https://pubmed.ncbi.nlm.nih.gov/{pub.pmid}/",
                }
                rows.append(row)

        if not rows:
            logger.info("No failed studies to export (all papers were accessible)")
            return None

        self._ensure_dirs()
        filepath = self.output_dir / filename
        df = pd.DataFrame(rows)
        df.to_csv(filepath, index=False)

        logger.info(f"Exported {len(rows)} failed studies to {filepath}")
        return str(filepath)
