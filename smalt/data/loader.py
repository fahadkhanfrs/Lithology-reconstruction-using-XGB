"""
SMALT Data Loader module for lithology data ingestion, standardization,
QC, synthetic GR generation, and unified dataset export.
"""

from data.loader import (
    LithologLoader,
    DEFAULT_FACIES_MAPPING,
    DEFAULT_BASE_GR,
    FACIES_ALIASES,
    CRITICAL_COLUMNS,
    load_provenance_manifest,
)

__all__ = [
    "LithologLoader",
    "DEFAULT_FACIES_MAPPING",
    "DEFAULT_BASE_GR",
    "FACIES_ALIASES",
    "CRITICAL_COLUMNS",
    "load_provenance_manifest",
]
