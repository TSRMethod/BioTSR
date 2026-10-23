"""Identity value objects.

Every type here is a frozen, hashable, ordered value object with no I/O and no
algorithms.  ``domain`` is the dependency sink of the package: it
imports nothing from ``biotsr``.

Two requirements shape this module:

* **No four-character structure-id assumption.**  ``StructureId`` is a generic,
  non-empty identifier string.  PDB identifiers, four characters today, are
  being extended to longer forms (the wwPDB's extended PDB ID), so no length is validated, encoded, or assumed
  anywhere in the package -- including filename construction and parsing.
* **Full physical atom identity.**  ``AtomIdentity`` carries structure, model,
  chain, residue name, residue number, insertion code, atom name and altLoc.
  Atom name alone is never an identity.  This is the load-bearing requirement
  established by the cross atom-identity audits (ADR-016, ADR-017): two legacy
  cross implementations key their atom pools by ``(side, atom_name)`` and
  silently collapse distinct physical atoms.  The domain model must never force
  that defect on a cross method.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "StructureId",
    "ChainId",
    "ResidueIdentity",
    "AtomIdentity",
]


@dataclass(frozen=True, order=True)
class StructureId:
    """A structure identifier of arbitrary length.

    Deliberately *not* a "PDB ID": no four-character rule, no alphabet rule
    beyond excluding whitespace and path separators (which would make derived
    filenames ambiguous or unsafe).
    """

    value: str

    def __post_init__(self) -> None:
        v = self.value
        if not isinstance(v, str) or not v:
            raise ValueError("StructureId must be a non-empty string")
        if v != v.strip():
            raise ValueError("StructureId must not have leading/trailing whitespace: %r" % v)
        if any(ch.isspace() for ch in v):
            raise ValueError("StructureId must not contain whitespace: %r" % v)
        for bad in ("/", "\\", "\x00"):
            if bad in v:
                raise ValueError("StructureId must not contain %r: %r" % (bad, v))

    @classmethod
    def parse(cls, raw: str) -> "StructureId":
        """Normalize a user-supplied identifier (strip, upper-case).

        Upper-casing matches the legacy family's ``file_name.upper()`` habit and
        the historical dataset sheet's own convention; it applies to identifiers of any
        length.
        """
        return cls(str(raw).strip().upper())

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value


@dataclass(frozen=True, order=True)
class ChainId:
    """A chain identifier.

    The PDB format gives one column, but mmCIF allows longer labels, so no
    length is enforced.  The empty string is a legal chain id (a blank PDB
    chain column) and is distinct from "no chain requested".
    """

    value: str

    def __post_init__(self) -> None:
        if not isinstance(self.value, str):
            raise ValueError("ChainId value must be a string")
        if any(ch.isspace() for ch in self.value):
            raise ValueError("ChainId must not contain whitespace: %r" % self.value)
        if "/" in self.value or "\\" in self.value:
            raise ValueError("ChainId must not contain a path separator: %r" % self.value)

    @property
    def is_blank(self) -> bool:
        return self.value == ""

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value


@dataclass(frozen=True, order=True)
class ResidueIdentity:
    """Identity of one residue instance within a structure.

    ``insertion_code`` is a first-class separate field, never concatenated into
    the sequence number (the legacy ``line[22:27]`` habit).
    ``legacy_seq_label`` reconstructs the concatenated form only where a legacy
    output format demands it.
    """

    chain: ChainId
    name: str
    seq: int
    insertion_code: str = ""
    model: int = 1

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("residue name must be a non-empty string")
        if self.name != self.name.strip():
            raise ValueError("residue name must be stripped: %r" % self.name)
        if not isinstance(self.seq, int):
            raise ValueError("residue seq must be an int")
        if not isinstance(self.insertion_code, str) or len(self.insertion_code) > 1:
            raise ValueError("insertion_code must be '' or a single character")
        if self.insertion_code.isspace():
            raise ValueError("insertion_code must be '' rather than whitespace")

    @property
    def legacy_seq_label(self) -> str:
        """``line[22:27].strip()`` as the legacy formats write it: number+iCode."""
        return "%d%s" % (self.seq, self.insertion_code)


@dataclass(frozen=True, order=True)
class AtomIdentity:
    """Full physical identity of one atom.

    Never reduce this to ``atom_name``: identically-named atoms in different
    residues are different atoms, and collapsing them is a demonstrated legacy
    defect with real output impact (ADR-016).
    """

    structure: StructureId
    model: int
    chain: ChainId
    residue_name: str
    residue_seq: int
    insertion_code: str = ""
    atom_name: str = ""
    alt_loc: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.model, int):
            raise ValueError("model must be an int")
        if not isinstance(self.atom_name, str) or not self.atom_name:
            raise ValueError("atom_name must be a non-empty string")
        if self.atom_name != self.atom_name.strip():
            raise ValueError("atom_name must be stripped: %r" % self.atom_name)
        if not isinstance(self.alt_loc, str) or len(self.alt_loc) > 1:
            raise ValueError("alt_loc must be '' or a single character")
        if self.alt_loc.isspace():
            raise ValueError("alt_loc must be '' rather than whitespace")
        # Validate residue-shaped fields through the same rules as ResidueIdentity.
        ResidueIdentity(
            chain=self.chain,
            name=self.residue_name,
            seq=self.residue_seq,
            insertion_code=self.insertion_code,
            model=self.model,
        )

    @property
    def residue(self) -> ResidueIdentity:
        return ResidueIdentity(
            chain=self.chain,
            name=self.residue_name,
            seq=self.residue_seq,
            insertion_code=self.insertion_code,
            model=self.model,
        )

    def as_tuple(self):
        """Canonical, hashable, fully-discriminating tuple form."""
        return (
            self.structure.value,
            self.model,
            self.chain.value,
            self.residue_name,
            self.residue_seq,
            self.insertion_code,
            self.atom_name,
            self.alt_loc,
        )

    def describe(self) -> str:
        """Human-readable full identity, e.g. ``1UBQ/1/A/LEU/15/CA``."""
        icode = self.insertion_code or ""
        alt = ("/altloc=%s" % self.alt_loc) if self.alt_loc else ""
        return "%s/%d/%s/%s/%d%s/%s%s" % (
            self.structure.value,
            self.model,
            self.chain.value,
            self.residue_name,
            self.residue_seq,
            icode,
            self.atom_name,
            alt,
        )
