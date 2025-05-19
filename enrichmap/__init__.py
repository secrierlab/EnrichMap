"""Spatially aware gene set enrichment analysis for spatial transcriptomics data."""

from __future__ import annotations

from . import tools as tl
from . import plotting as pl

import sys
sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["tl", "pl"]})
del sys

__all__ = [
    "__version__",
    "pl",
    "tl"
]
