"""
NCBI Paper Query Package

A modular toolkit for retrieving publications from NCBI/PubMed,
downloading papers, and extracting omics datasets.

Usage:
    from ncbi_paper_query import retrieve_publications, PubMedSearcher

    # Or run as module:
    python -m ncbi_paper_query --disease "Crohn's disease" --tissue intestine
"""

# Core configuration
from .config import (
    SKILL_DIR,
    NCBI_EMAIL,
    NCBI_API_KEY,
    CLAUDE_CLI_AVAILABLE,
    logger,
)

# Data models
from .models import OmicsDataset, Publication

# Components - now imported from modular files
from .retrieval import PubMedSearcher
from .matching import EntityMatcher
from .metadata import ImpactFactorLookup
from .download import PaperDownloader, OmicsExtractor
from .export import ResultsExporter
from .core import retrieve_publications, validate_results, compare_rounds

__version__ = "1.2.0"

__all__ = [
    # Config
    'SKILL_DIR',
    'NCBI_EMAIL',
    'NCBI_API_KEY',
    'CLAUDE_CLI_AVAILABLE',
    'logger',
    # Models
    'OmicsDataset',
    'Publication',
    # Components
    'PubMedSearcher',
    'EntityMatcher',
    'ImpactFactorLookup',
    'PaperDownloader',
    'OmicsExtractor',
    'ResultsExporter',
    # Functions
    'retrieve_publications',
    'validate_results',
    'compare_rounds',
]
