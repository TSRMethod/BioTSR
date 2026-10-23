"""The structure layer -- the only layer that knows file formats."""

from .model import AtomRecord, Structure, TruncationEvent

__all__ = ["AtomRecord", "Structure", "TruncationEvent"]
