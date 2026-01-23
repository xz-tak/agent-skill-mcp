"""
NCBI Paper Query - Publication Model

Represents a publication with all metadata.
"""

from dataclasses import dataclass, field
from typing import List, Optional

from .omics import OmicsDataset


@dataclass
class Publication:
    """Represents a publication with all metadata."""
    pmid: str
    title: str
    abstract: str
    authors: List[str]
    affiliations: List[str]
    journal: str
    publication_year: int
    doi: Optional[str] = None
    pmc_id: Optional[str] = None

    # Matched entities
    matched_disease: Optional[str] = None
    matched_tissue: Optional[str] = None
    matched_organism: Optional[str] = None
    matched_cell_type: Optional[str] = None
    relevance_score: float = 0.0

    # Journal metrics
    impact_factor: Optional[float] = None
    is_preprint: bool = False

    # Access information
    free_access: bool = False
    requires_subscription: bool = False
    access_url: Optional[str] = None
    download_path: Optional[str] = None
    web_scraped: bool = False  # True when content accessed via web-scraping without downloading

    # Omics datasets
    omics_datasets: List[OmicsDataset] = field(default_factory=list)
