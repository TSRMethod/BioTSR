"""The shared mixed-radix key encoder.

Deliberately narrow: it knows three radices and nothing else.  It does not know
how the PDB was parsed, how points were selected, how triangles were generated,
which method is running, or whether a future vectorized generator produced the
candidate.  That narrowness is why one encoder can serve every method in the
family.

    key = dLen*dTheta*L^2*(c1-1) + dLen*dTheta*L*(c2-1) + dLen*dTheta*(c3-1)
        + dTheta*(distance_class-1) + (theta_class-1)

``L``, ``dLen`` and ``dTheta`` are positional radices.  Changing any of them
re-maps *every* key -- the same integer denotes a different ``(labels, bins)``
tuple -- so there is no partial compatibility and no arithmetic migration between
radix triples.  That is why they live in ``KeySpace``.

Mirror signing is applied *after* this function, by
``modifiers/mirror.py``, because the sign depends on ``theta_raw`` and the
magnitude does not.  Keeping them apart is what makes ``abs(signed_key)`` provably
equal to the unsigned key for the same triangle -- a property the tests assert
directly rather than assume.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

__all__ = ["KeyFormulaSpec", "encode", "KEY_FORMULA_VERSION"]

#: Bumped only when the formula itself changes.  Enters ``KeySpace``.
KEY_FORMULA_VERSION = "mixed_radix_v1"


@dataclass(frozen=True)
class KeyFormulaSpec:
    """The three radices the formula multiplies by.

    ``label_range`` is the label scheme's ``(lowest, highest)`` code.  When it is
    given, :func:`encode` refuses a label outside it, and the range itself must
    fit in one radix window -- the condition under which no label digit can
    carry into its neighbour (ADR-107).  Every profile supplies it
    (``Profile.key_formula``); it is not part of ``KeySpace``, because it adds a
    check and changes no key.
    """

    n_labels: int
    n_theta: int
    n_distance: int
    version: str = KEY_FORMULA_VERSION
    label_range: Optional[Tuple[int, int]] = None

    def __post_init__(self) -> None:
        for field_name in ("n_labels", "n_theta", "n_distance"):
            value = getattr(self, field_name)
            if not isinstance(value, int) or value <= 0:
                raise ValueError("%s must be a positive int" % field_name)
        if self.label_range is not None:
            lo, hi = self.label_range
            if lo < 1 or hi < lo or hi - lo >= self.n_labels:
                raise ValueError(
                    "label range %d..%d does not fit in one radix-%d window of "
                    "codes >= 1; labels outside it would carry into the "
                    "neighbouring digit" % (lo, hi, self.n_labels))


def encode(
    ordered_labels: Sequence[int],
    theta_class: int,
    distance_class: int,
    spec: KeyFormulaSpec,
) -> int:
    """Encode one canonically ordered triangle into its integer key.

    ``ordered_labels`` is ``(l1, l2, l3)`` -- already canonically ordered.  The
    encoder does no ordering of its own; handing it unordered labels produces a
    wrong key silently, which is exactly why ordering is a separate, exhaustively
    tested component.
    """
    if len(ordered_labels) != 3:
        raise ValueError("expected exactly three ordered labels")
    c1, c2, c3 = ordered_labels
    if spec.label_range is not None:
        lo, hi = spec.label_range
        if not (lo <= c1 <= hi and lo <= c2 <= hi and lo <= c3 <= hi):
            raise ValueError(
                "label(s) %r outside the label scheme's codes %d..%d; a code "
                "outside that window would carry into the neighbouring label "
                "digit and alias another triple (ADR-107)"
                % ((c1, c2, c3), lo, hi))
    if theta_class < 1 or theta_class > spec.n_theta:
        raise ValueError(
            "theta_class %r outside 1..%d" % (theta_class, spec.n_theta)
        )
    if distance_class < 1 or distance_class > spec.n_distance:
        raise ValueError(
            "distance_class %r outside 1..%d" % (distance_class, spec.n_distance)
        )

    d_len = spec.n_distance
    d_theta = spec.n_theta
    n_labels = spec.n_labels

    return (
        d_len * d_theta * (n_labels ** 2) * (c1 - 1)
        + d_len * d_theta * n_labels * (c2 - 1)
        + d_len * d_theta * (c3 - 1)
        + d_theta * (distance_class - 1)
        + (theta_class - 1)
    )
