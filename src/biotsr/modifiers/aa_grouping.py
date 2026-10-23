"""The AA-grouping modifier: a label-scheme swap, hooked at ``LABELS``.

Grouping is **not** a flag on the engine and **not** a second selector.  It
replaces the label scheme before canonical ordering, and everything downstream --
ordering, geometry, binning, encoding, serialization -- is untouched code running
on different integers.

Two rules this module enforces:

1. **Cα-only.**  AA grouping requires the protein amino-acid label domain.  It is
   not exposed for amino-acid atomic, drug, carbohydrate or nucleotide TSR, *even
   though those points carry residue metadata* -- an explicit TSR Research
   Group decision.  The check is against the label scheme's declared ``domain``,
   not against which selector produced the points, so it stays correct as new
   selectors are added.
2. **The radix follows the codes.**  The standard grouped scheme
   (``AMINO_ACID_GROUP_23``) has 23 codes spanning 4..33 and declares radix 30,
   the span that keeps its encoding collision-free (ADR-107).  The legacy
   behaviour -- the 13 historical groups encoded with an unrelated 20-slot
   radix -- was a defect (ADR-015) and is not implemented.  Grouping keeps its
   intended effect: residues of one group share one code, so their triangles
   share keys on purpose.
"""

from __future__ import annotations

from ..labels.schemes import LabelDomain, LabelScheme
from .base import IncompatibleModifier, ModifierHook, ModifierName

__all__ = [
    "HOOK",
    "NAME",
    "REQUIRED_LABEL_DOMAIN",
    "check_applicable",
    "apply_grouping",
]

NAME = ModifierName.AA_GROUPING
HOOK = ModifierHook.LABELS

#: The one domain in which grouping is meaningful.
REQUIRED_LABEL_DOMAIN = LabelDomain.PROTEIN_AMINO_ACID


def check_applicable(base_scheme: LabelScheme, *, method: str = "") -> None:
    """Raise unless ``base_scheme`` is in the protein amino-acid label domain."""
    if base_scheme.domain != REQUIRED_LABEL_DOMAIN:
        raise IncompatibleModifier(
            "aa_grouping requires a Protein Ca label domain (%r); the %slabel "
            "scheme %r is in domain %r. Amino-acid grouping is defined only over "
            "amino-acid residue identities, so it is not available for atomic, "
            "drug, carbohydrate or nucleotide methods even though their points "
            "carry residue metadata."
            % (
                REQUIRED_LABEL_DOMAIN,
                ("%s method's " % method) if method else "",
                base_scheme.name,
                base_scheme.domain,
            )
        )


def apply_grouping(
    base_scheme: LabelScheme, grouped_scheme: LabelScheme, *, method: str = ""
) -> LabelScheme:
    """Swap ``base_scheme`` for ``grouped_scheme`` after checking applicability.

    ``grouped_scheme`` is supplied by the profile rather than hardcoded here,
    so this module never restates which grouped scheme a profile uses.
    """
    check_applicable(base_scheme, method=method)
    if grouped_scheme.domain != REQUIRED_LABEL_DOMAIN:
        raise IncompatibleModifier(
            "the grouped scheme %r must also be in domain %r, not %r"
            % (grouped_scheme.name, REQUIRED_LABEL_DOMAIN, grouped_scheme.domain)
        )
    missing = sorted(set(base_scheme.table) - set(grouped_scheme.table))
    if missing:
        raise IncompatibleModifier(
            "grouped scheme %r has no group for %d residue(s) the base scheme %r "
            "labels: %s" % (grouped_scheme.name, len(missing), base_scheme.name, ", ".join(missing))
        )
    return grouped_scheme
