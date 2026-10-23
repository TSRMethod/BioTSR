"""Label schemes: residue/atom identity -> integer code.

The two standard schemes here are both in the **protein amino-acid** label domain:

* ``AMINO_ACID_30`` -- the default: the 20 regular amino acids and the
  supported modified residues, 30 codes, radix **30**;
* ``AMINO_ACID_GROUP_23`` -- the optional AA-grouping modifier's target: the
  same residues mapped onto 23 groups, radix **30**.

Both tables are the TSR Research Group's supplied tables
(``aminoAcidCode_lexicographic_new.txt`` and ``aminoAcidCode_grouping.txt``),
transcribed verbatim, and base 30 is the radix the TSR Research Group set for both
(ADR-107).  They share a radix but not a key space: the scheme's name and
content fingerprint enter ``KeySpace``, so a grouped key and an ungrouped key are
never compared even when they are the same integer.

A scheme declares a ``domain``.  That is what makes modifier compatibility a
property of the *points* rather than a hardcoded per-selector rule: AA grouping
requires the protein amino-acid domain and is rejected anywhere else, including
for atomic/ligand/carbohydrate/nucleotide schemes that will carry residue
metadata but are not in that domain (see ``modifiers/aa_grouping.py``).

**What makes an encoding collision-free.**  The encoder
(``encoding/key.py``) writes the three labels as positional digits
``c - 1`` in base ``n_labels``.  Codes do not start at 1 -- the residue codes run
4..33 -- so the digit of the smallest code is not 0.  What matters is the
**span**: when every code lies in one window ``lo .. lo + n_labels - 1``, the
three digits each stay inside one radix step, no digit can carry into its
neighbour, and distinct ordered label triples always give distinct label
components.  A code outside that window does carry: with the supplied table at
radix 20, ``(SEP, PHE, PHE) = (25, 4, 4)`` and ``(TPO, TPO, PHE) = (24, 24, 4)``
both give 9663.  The number of distinct codes is not the test -- the grouped
table has 23 codes but spans 4..33 (17..23 are unused), so radix 23 still
collides.  ``LabelScheme`` therefore refuses any scheme whose codes do not fit in
one radix window, and :func:`label_triple_collisions` enumerates every
canonically ordered triple to show that the shipped schemes have none.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Mapping, Optional, Tuple

from ..domain.keyspace import SchemeRef

__all__ = [
    "LabelDomain",
    "LabelSource",
    "LabelScheme",
    "AMINO_ACID_30",
    "AMINO_ACID_GROUP_23",
    "STANDARD_AMINO_ACIDS",
    "MODIFIED_RESIDUE_CODES",
    "UnlabelledPointPolicy",
    "label_triple_collisions",
    "ordered_label_triples",
]


class LabelDomain:
    """What a scheme's labels *mean*.

    Modifier compatibility is validated against this, not against which selector
    produced the points.
    """

    #: One code per amino-acid residue identity.  The Ca family, and the only
    #: domain in which AA grouping is meaningful.
    PROTEIN_AMINO_ACID = "protein_amino_acid"

    #: One code per chemical element, from the shared element table (radix
    #: 112).  Used by the atomic/ligand/carbohydrate/nucleotide families.
    ATOM_ELEMENT = "atom_element"


class LabelSource:
    """What a scheme looks a point up by."""

    #: The point's residue name (the Ca residue schemes).
    RESIDUE_NAME = "residue_name"
    #: The atom's element field, PDB columns 77-78, case-insensitive (the
    #: atomic scheme, ADR-113).  The atom name is never consulted.
    ELEMENT = "element"


class UnlabelledPointPolicy:
    """What happens to a point the scheme has no code for.

    Only ``ERROR`` is implemented, which is what the legacy protein Ca path
    does (an uncaught ``KeyError`` on an unknown residue name).
    ``DROP_AND_RECORD`` is named but not wired up until a method needs it, so
    nothing is silently dropped in the meantime.
    """

    ERROR = "error"
    DROP_AND_RECORD = "drop_and_record"


def normalize_element(symbol: str) -> str:
    """An element field in chemical case: ``" S" -> "S"``, ``"CL" -> "Cl"``, ``"cl" -> "Cl"``.

    The PDB writes the element right-justified and upper-case (``" S"``,
    ``"CL"``); the element table spells symbols in chemical case (``S``,
    ``Cl``).  Reading the field case-insensitively is the approved rule
    (ADR-108); nothing else is changed -- a field that is not an element symbol
    stays what it is and finds no label.
    """
    text = (symbol or "").strip()
    if not text:
        return ""
    return text[0].upper() + text[1:].lower()


@dataclass(frozen=True)
class LabelScheme:
    """A named lookup from a point's residue name or element to an integer code, plus its radix."""

    name: str
    table: Mapping[str, int]
    n_labels: int
    domain: str = LabelDomain.PROTEIN_AMINO_ACID
    description: str = ""
    #: What ``table`` is keyed by (:class:`LabelSource`).
    label_source: str = LabelSource.RESIDUE_NAME
    #: The named rule that turns a point into a table key, recorded in the
    #: content fingerprint; empty for the residue schemes, whose key is the
    #: residue name as written.
    resolver: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.n_labels, int) or isinstance(self.n_labels, bool) \
                or self.n_labels <= 0:
            raise ValueError("n_labels (the encoder radix) must be a positive int")
        codes = sorted(self.all_codes())
        if not codes:
            raise ValueError("a label scheme needs at least one code")
        bad = [c for c in codes if not isinstance(c, int) or isinstance(c, bool) or c < 1]
        if bad:
            raise ValueError(
                "scheme %r has code(s) %s; every code must be an integer >= 1, "
                "because the encoder writes code - 1 as a digit" % (self.name, bad))
        lo, hi = codes[0], codes[-1]
        if hi - lo + 1 > self.n_labels:
            overflow = [c for c in codes if c - lo >= self.n_labels]
            raise ValueError(
                "scheme %r has codes %d..%d (a span of %d) but declares radix %d: "
                "code(s) %s lie outside one radix window, so they would carry into "
                "the neighbouring label digit and different label triples could "
                "share a key (ADR-107). The radix must be at least the span of the "
                "codes; the number of distinct codes (%d) is not the test."
                % (self.name, lo, hi, hi - lo + 1, self.n_labels,
                   ", ".join(str(c) for c in overflow), len(codes)))

    def all_codes(self):
        """Every code the scheme can emit."""
        return set(self.table.values())

    @property
    def n_distinct_codes(self) -> int:
        """The size of the label alphabet, which is *not* necessarily the radix."""
        return len(self.all_codes())

    @property
    def code_range(self) -> Tuple[int, int]:
        codes = self.all_codes()
        return (min(codes), max(codes))

    @property
    def codes_are_consecutive(self) -> bool:
        lo, hi = self.code_range
        return hi - lo + 1 == self.n_distinct_codes

    @property
    def encoding_is_injective(self) -> bool:
        """Whether distinct ordered label triples always encode to distinct keys.

        True when every code lies in one window of ``n_labels`` consecutive
        integers: the positional digits then stay within one radix step, so no
        carry can alias two triples together.  Gaps inside the window are
        harmless.  A constructed scheme always satisfies it, because
        ``__post_init__`` refuses one that does not; it is kept as a property so
        the reason can be read off a scheme.
        """
        lo, hi = self.code_range
        return (hi - lo) < self.n_labels

    def label(self, name: str) -> Optional[int]:
        return self.table.get(name)

    def resolve_atom(self, atom_name: str, element: str = "") -> Optional[int]:
        """Resolve one point to its integer label, or ``None`` when it has none.

        An **element** scheme (the atomic scheme, ADR-113) looks up the atom's
        element field, read case-insensitively (:func:`normalize_element`):

        ==============================================  =========  =====
        record                                          element    label
        ==============================================  =========  =====
        ``CA`` of an amino acid, element ``C``          carbon     106
        ``CA`` calcium ion, element ``CA``              calcium     99
        ``SG`` of cysteine, element ``S``               sulfur     111
        ``PA`` of ATP, element ``P``                    phosphorus 110
        ``CL35`` of a ligand, element ``CL``            chlorine    93
        ``HO`` holmium ion, element ``HO``              holmium     37
        ==============================================  =========  =====

        ``atom_name`` is accepted so every scheme is called the same way, and is
        not read: the same element always gets the same label, whatever the
        atom is called.  A blank field, or one naming no element in the table,
        has no label; nothing is inferred from the name.

        A **residue** scheme (the Ca schemes) looks up the name it is given --
        the residue name -- exactly as written.
        """
        if self.label_source == LabelSource.ELEMENT:
            return self.table.get(normalize_element(element))
        return self.table.get(atom_name)

    def content_fingerprint(self) -> str:
        payload = {
            "table": {k: v for k, v in sorted(self.table.items())},
            "n_labels": self.n_labels,
            "domain": self.domain,
        }
        if self.label_source == LabelSource.RESIDUE_NAME:
            # The residue schemes' payload exactly as it has been since 1.0.0rc1,
            # including two fields only the retired atom-name lexicon used, so
            # their fingerprints -- and every Ca key space -- are unchanged.
            payload["element_table"] = {}
            payload["namespace_aware"] = False
        else:
            # How a point becomes a key into the table is part of what a label
            # means: the same table read through another rule labels atoms
            # differently, so the rule enters the fingerprint (ADR-113).
            payload["label_source"] = self.label_source
            payload["resolver"] = self.resolver
        blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(blob.encode("utf-8")).hexdigest()

    def as_ref(self) -> SchemeRef:
        return SchemeRef(
            name=self.name,
            content_fingerprint=self.content_fingerprint(),
            radix=self.n_labels,
        )


#: The 20 regular amino acids, in the order of the supplied tables' first rows.
STANDARD_AMINO_ACIDS = (
    "PHE", "CYS", "GLN", "GLU", "LEU", "HIS", "ILE", "GLY", "LYS", "MET",
    "ASP", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ALA", "ARG", "ASN",
)

#: The supported modified residues and their codes in the ungrouped table, in
#: the supplied table's order.  ``MSE`` (selenomethionine) shares ``MET``'s code
#: on purpose; every other one has a code of its own.  The Ca selectors read
#: these residues, including from ``HETATM`` records (ADR-107).
MODIFIED_RESIDUE_CODES = (
    ("TPO", 24),  # phosphothreonine
    ("SEP", 25),  # phosphoserine
    ("PTR", 26),  # phosphotyrosine
    ("MSE", 13),  # selenomethionine -> MET
    ("NIY", 27),  # 3-nitrotyrosine
    ("M3L", 28),  # N6,N6,N6-trimethyllysine
    ("TRF", 29),  # N1'-formyl-tryptophan
    ("CSO", 30),  # S-hydroxycysteine
    ("LLP", 31),  # lysine-pyridoxal-5'-phosphate
    ("CME", 32),  # S,S-(2-hydroxyethyl)thiocysteine
    ("CSD", 33),  # 3-sulfinoalanine
)

#: The default: the supplied ``aminoAcidCode_lexicographic_new.txt``
#: (sha256 ``c400c313bc67...25f4``), verbatim.  Its first 20 rows are the 20
#: regular amino-acid codes of ``TSR-Package/tsr_package/TSR.py``'s
#: ``aminoAcidLabel`` (4..23); the modified residues follow (24..33, MSE = 13).
#: 30 distinct codes, 4..33, at radix 30 (ADR-107).
AMINO_ACID_30 = LabelScheme(
    name="amino_acid_30",
    table={
        "PHE": 4, "CYS": 5, "GLN": 6, "GLU": 7, "LEU": 8, "HIS": 9,
        "ILE": 10, "GLY": 11, "LYS": 12, "MET": 13, "ASP": 14, "PRO": 15,
        "SER": 16, "THR": 17, "TRP": 18, "TYR": 19, "VAL": 20, "ALA": 21,
        "ARG": 22, "ASN": 23,
        "TPO": 24, "SEP": 25, "PTR": 26, "MSE": 13, "NIY": 27, "M3L": 28,
        "TRF": 29, "CSO": 30, "LLP": 31, "CME": 32, "CSD": 33,
    },
    n_labels=30,
    domain=LabelDomain.PROTEIN_AMINO_ACID,
    description=(
        "Residue labels: the 20 regular amino acids (codes 4..23) and the "
        "supported modified residues TPO 24, SEP 25, PTR 26, NIY 27, M3L 28, "
        "TRF 29, CSO 30, LLP 31, CME 32, CSD 33, with MSE sharing MET's 13; "
        "30 codes, 4..33, encoder radix 30."
    ),
)

#: The AA-grouping modifier's target: the supplied
#: ``aminoAcidCode_grouping.txt`` (sha256 ``f40b773fefe9...51b5``), verbatim.
#: Its first 20 rows are ``TSR.py``'s ``aminoAcidGroup`` (13 groups, 4..16); the
#: modified residues follow with the same codes as the ungrouped table, except
#: MSE, which joins MET's group (14).  23 distinct codes spanning 4..33 -- 17..23
#: are unused -- so the radix that keeps it collision-free is the span, 30, not
#: the code count (ADR-107).
AMINO_ACID_GROUP_23 = LabelScheme(
    name="amino_acid_group_23",
    table={
        "ASP": 4, "GLU": 4, "ASN": 5, "GLN": 5, "GLY": 6, "ALA": 7, "VAL": 7,
        "LEU": 8, "ILE": 8, "PRO": 9, "PHE": 10, "TRP": 10, "TYR": 11, "SER": 12,
        "THR": 12, "CYS": 13, "MET": 14, "HIS": 15, "LYS": 16, "ARG": 16,
        "TPO": 24, "SEP": 25, "PTR": 26, "MSE": 14, "NIY": 27, "M3L": 28,
        "TRF": 29, "CSO": 30, "LLP": 31, "CME": 32, "CSD": 33,
    },
    n_labels=30,
    domain=LabelDomain.PROTEIN_AMINO_ACID,
    description=(
        "Grouped residue labels: the 20 regular amino acids in 13 "
        "physicochemical groups (codes 4..16), MSE in MET's group, and each "
        "other supported modified residue in a group of its own (24..33); 23 "
        "codes spanning 4..33, encoder radix 30. A different key space from "
        "the ungrouped labels although the radix is the same."
    ),
)


def ordered_label_triples(scheme: LabelScheme):
    """Every label triple the canonical ordering can emit, ``c1 >= c2 >= c3``.

    The ordering sorts labels in descending order before the encoder sees them
    (``geometry/ordering.py``), so these are all the triples that can occur:
    4,960 for 30 codes, 2,300 for 23.
    """
    import itertools

    codes = sorted(scheme.all_codes(), reverse=True)
    return list(itertools.combinations_with_replacement(codes, 3))


def label_triple_collisions(scheme: LabelScheme):
    """``{label_component: [triples]}`` for every label component that two or
    more distinct ordered triples share under ``scheme``'s radix.

    Exhaustive, through the encoder itself: with one theta and one distance
    class the key *is* the label component.  Empty for every scheme BioTSR
    ships; used by the tests and by anyone checking a table.
    """
    from ..encoding.key import KeyFormulaSpec, encode

    spec = KeyFormulaSpec(n_labels=scheme.n_labels, n_theta=1, n_distance=1)
    seen = {}
    for triple in ordered_label_triples(scheme):
        seen.setdefault(encode(triple, 1, 1, spec), []).append(triple)
    return {k: v for k, v in seen.items() if len(v) > 1}
