"""Bin schemes: continuous value -> class index.

Schemes are **values**, not modules.  Each is explicitly named and content
fingerprinted, and the fingerprint hashes the scheme's actual semantics --
boundaries, declared class count, interval semantics, and out-of-range behaviour
-- together.  Hashing the boundary array alone is not enough: the two Ca distance
schemes share every cut point and differ *only* in what happens past the last
one, and they produce different keys for any triangle at or beyond 132.45 A
(ADR-014).

``n_classes`` is the **encoder radix**, which is a separate fact from the highest
class the binner can actually emit.  The legacy Ca binner declares 35 and can
only reach 34 -- a permanently unreachable radix slot that went unnoticed for
years.  ``LadderBinScheme`` therefore validates the two against each other and
requires an explicit, reasoned opt-out to express the legacy defect.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Tuple

from ..domain.keyspace import SchemeRef

__all__ = ["LadderBinScheme"]


@dataclass(frozen=True)
class LadderBinScheme:
    """A monotone ladder binner: ``class = 1 + index of the first boundary the value is below``.

    Values at or above every boundary fall into ``overflow_class``.  Intervals are
    lower-inclusive / upper-exclusive throughout, matching every legacy table in
    the family.
    """

    name: str
    boundaries: Tuple[float, ...]
    n_classes: int
    overflow_class: int
    interval_semantics: str = "lower_inclusive_upper_exclusive"
    allow_unreachable_classes: bool = False
    unreachable_reason: str = ""

    def __post_init__(self) -> None:
        if not self.boundaries:
            raise ValueError("a bin scheme needs at least one boundary")
        prev = None
        for b in self.boundaries:
            if prev is not None and not (b > prev):
                raise ValueError(
                    "boundaries must be strictly increasing: %r" % (self.boundaries,)
                )
            prev = b
        if self.overflow_class > self.n_classes:
            raise ValueError(
                "overflow_class %d exceeds declared n_classes %d"
                % (self.overflow_class, self.n_classes)
            )
        if self.max_emitted_class != self.n_classes and not self.allow_unreachable_classes:
            raise ValueError(
                "scheme %r declares %d classes but can only emit %d; if that is a "
                "legacy defect being reproduced, set allow_unreachable_classes=True "
                "with an unreachable_reason" % (self.name, self.n_classes, self.max_emitted_class)
            )
        if self.allow_unreachable_classes and not self.unreachable_reason:
            raise ValueError("allow_unreachable_classes requires an unreachable_reason")

    @property
    def max_emitted_class(self) -> int:
        """The highest class ``classify`` can actually return."""
        return max(len(self.boundaries), self.overflow_class)

    def classify(self, value: float) -> int:
        for i, boundary in enumerate(self.boundaries):
            if value < boundary:
                return i + 1
        return self.overflow_class

    def content_fingerprint(self) -> str:
        """Hash of the scheme's semantics.  Deliberately excludes ``name``.

        Names can lie -- that is the whole reason this exists -- so two schemes
        with identical semantics fingerprint identically regardless of naming,
        and two schemes with identical *boundaries* but different overflow
        behaviour fingerprint differently.
        """
        payload = {
            "boundaries": [repr(b) for b in self.boundaries],
            "n_classes": self.n_classes,
            "overflow_class": self.overflow_class,
            "max_emitted_class": self.max_emitted_class,
            "interval_semantics": self.interval_semantics,
        }
        blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(blob.encode("utf-8")).hexdigest()

    def as_ref(self) -> SchemeRef:
        return SchemeRef(
            name=self.name,
            content_fingerprint=self.content_fingerprint(),
            radix=self.n_classes,
        )
