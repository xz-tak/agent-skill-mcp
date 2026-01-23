"""
NCBI Paper Query - OmicsDataset Model

Represents an omics dataset extracted from a publication.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class OmicsDataset:
    """Represents an omics dataset extracted from a publication."""
    accession: str
    database: str  # GEO, SRA, ArrayExpress, PRIDE, MetaboLights
    technology: Optional[str] = None
    platform: Optional[str] = None
    sample_count: Optional[int] = None
    conditions: Optional[str] = None
    condition_counts: Optional[str] = None  # Formatted as "Condition1 (N1); Condition2 (N2)"
    data_url: Optional[str] = None
    assay_description: Optional[str] = None
    experiment_design: Optional[str] = None
    organism: Optional[str] = None
    title: Optional[str] = None
    source: str = "unknown"  # Where accession was found: "data_availability", "key_resources", "full_text", "abstract"
