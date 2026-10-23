"""Composable modifiers.

One shared pipeline, three different internal hooks.  There are deliberately no
combinatorial per-combination classes -- no ``GroupedMirrorTSR``, no
``FilteredGroupedTSR`` -- because every combination is the same engine reading
different configuration.
"""

from .aa_grouping import apply_grouping, check_applicable
from .base import IncompatibleModifier, ModifierHook, ModifierName
from .mirror import (
    MIRROR_CATEGORY,
    MIRROR_DISABLED,
    MIRROR_NEGATE_WHEN_PREFOLD_GT_90,
    MirrorConvention,
)
from .size_filter import Comparison, SizeFilter

__all__ = [
    "Comparison",
    "IncompatibleModifier",
    "MIRROR_CATEGORY",
    "MIRROR_DISABLED",
    "MIRROR_NEGATE_WHEN_PREFOLD_GT_90",
    "MirrorConvention",
    "ModifierHook",
    "ModifierName",
    "SizeFilter",
    "apply_grouping",
    "check_applicable",
]
