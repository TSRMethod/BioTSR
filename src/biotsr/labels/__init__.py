"""Label schemes: residue name or element -> integer code, plus the encoder radix."""

from .elements import (
    ATOMIC_LABEL_RADIX,
    ELEMENT_LABELS_112,
    ELEMENT_RESOLVER,
    element_label_scheme,
)
from .schemes import (
    AMINO_ACID_30,
    AMINO_ACID_GROUP_23,
    MODIFIED_RESIDUE_CODES,
    STANDARD_AMINO_ACIDS,
    LabelDomain,
    LabelScheme,
    LabelSource,
    UnlabelledPointPolicy,
    label_triple_collisions,
    ordered_label_triples,
)

__all__ = [
    "AMINO_ACID_30",
    "ELEMENT_LABELS_112",
    "ATOMIC_LABEL_RADIX",
    "ELEMENT_RESOLVER",
    "element_label_scheme",
    "AMINO_ACID_GROUP_23",
    "MODIFIED_RESIDUE_CODES",
    "STANDARD_AMINO_ACIDS",
    "LabelDomain",
    "LabelScheme",
    "LabelSource",
    "UnlabelledPointPolicy",
    "label_triple_collisions",
    "ordered_label_triples",
]
