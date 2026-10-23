"""Parser backends.  BioTSR 1.0 ships the legacy-offset reader only."""

from . import pdb_legacy

__all__ = ["pdb_legacy"]
