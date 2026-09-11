"""
SMALT Data Ingestion and Standardization Subpackage.
"""

from data.loader import (
    LithologLoader,
    DEFAULT_FACIES_MAPPING,
    DEFAULT_BASE_GR,
    FACIES_ALIASES,
    CRITICAL_COLUMNS,
)

__all__ = [
    "LithologLoader",
    "DEFAULT_FACIES_MAPPING",
    "DEFAULT_BASE_GR",
    "FACIES_ALIASES",
    "CRITICAL_COLUMNS",
]
