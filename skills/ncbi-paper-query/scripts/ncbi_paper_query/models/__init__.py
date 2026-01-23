"""
NCBI Paper Query - Data Models

Exports OmicsDataset and Publication dataclasses.
"""

from .omics import OmicsDataset
from .publication import Publication

__all__ = ['OmicsDataset', 'Publication']
