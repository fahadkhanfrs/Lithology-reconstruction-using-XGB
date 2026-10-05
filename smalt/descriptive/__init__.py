"""
SMALT Descriptive Analysis and Stratigraphic Validation Package.
"""

from smalt.descriptive.analyzer import (
    LithologInspector,
    CANONICAL_FACIES_SCHEMA,
    PROVENANCE_METADATA,
    normalize_facies_label,
)

__all__ = [
    "LithologInspector",
    "CANONICAL_FACIES_SCHEMA",
    "PROVENANCE_METADATA",
    "normalize_facies_label",
]
