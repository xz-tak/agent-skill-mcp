"""
NCBI Paper Query - Core Module

Exports main retrieval and validation functions.
"""

from .retrieval import retrieve_publications
from .validation import validate_results, compare_rounds

__all__ = ['retrieve_publications', 'validate_results', 'compare_rounds']
