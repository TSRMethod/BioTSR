"""The mirror modifier: a key sign convention, hooked at ``KEY_SIGN``.

The rule, validated against 84,727 real triplets across four families with zero
mismatches (ADR-018)::

    mirror disabled                     -> never negate, under any circumstance
    mirror enabled, theta_raw <= 90     -> positive key, category "regular"
    mirror enabled, theta_raw  > 90     -> negative key, category "mirror"

Two things the sign is keyed off, which are easy to confuse:

* the **pre-fold** angle ``theta_raw`` decides the sign;
* the **folded** angle ``theta`` decides the theta class.

So mirroring changes only the sign of a key whose magnitude is already fully
determined -- it is a pure sign multiplier, never a second key space of its own.

One edge case survives from legacy and is preserved rather than papered over:
``key == 0`` negates to itself, so a mirrored zero-key is indistinguishable from
an unmirrored one *by sign alone*.  The explicit ``category`` field closes that
hole without changing the key, which is why the result model reports category
alongside frequency rather than leaving the sign to carry the whole meaning.

This module is family-agnostic: it reads ``theta_raw`` and nothing else, so the
amino-acid atomic, drug, carbohydrate and nucleotide families reuse it unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..domain.result import DEFAULT_CATEGORY
from .base import ModifierHook, ModifierName

__all__ = [
    "HOOK",
    "NAME",
    "MIRROR_CATEGORY",
    "MirrorConvention",
    "MIRROR_DISABLED",
    "MIRROR_NEGATE_WHEN_PREFOLD_GT_90",
]

NAME = ModifierName.MIRROR
HOOK = ModifierHook.KEY_SIGN

#: The category recorded for a negated key.
MIRROR_CATEGORY = "mirror"


@dataclass(frozen=True)
class MirrorConvention:
    """Whether, and when, a key is negated.

    The convention enters ``KeySpace`` because it changes what an integer means:
    the same triangle yields ``k`` under one convention and ``-k`` under the
    other.
    """

    name: str
    enabled: bool
    threshold: float = 90.0

    def category_for(self, theta_raw: float) -> str:
        if self.enabled and theta_raw > self.threshold:
            return MIRROR_CATEGORY
        return DEFAULT_CATEGORY

    def apply(self, key: int, theta_raw: float) -> int:
        """Return the signed key.

        Strictly greater than the threshold, never `>=`: a pre-fold angle of
        exactly 90 is *not* mirrored.
        """
        if self.enabled and theta_raw > self.threshold:
            return -key
        return key


#: The default: never negate.
MIRROR_DISABLED = MirrorConvention(name="disabled_never_negate", enabled=False)

#: Mirror mode. Never inferred from method or family -- it is requested or it is
#: off.
MIRROR_NEGATE_WHEN_PREFOLD_GT_90 = MirrorConvention(
    name="negate_when_prefold_angle_gt_90", enabled=True, threshold=90.0
)
