"""Domain value objects, with the two requirements they encode.

1. No four-character structure-id assumption anywhere.
2. Full physical atom identity -- never atom name alone.
"""

from __future__ import annotations

import pytest

from biotsr.domain.candidate import TriangleCandidate
from biotsr.domain.identity import AtomIdentity, ChainId, ResidueIdentity, StructureId
from biotsr.domain.point import LabeledPoint, Point, PointSet
from biotsr.domain.result import CategoryCounts, DEFAULT_CATEGORY


# ---------------------------------------------------------------------------
# StructureId -- length-agnostic by requirement
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value",
    [
        "1UBQ",
        "6G9M",
        "PDB_0000001UBQ",  # the 12-character extended PDB ID form
        "ABCDEFGHIJ",  # exactly ten
        "A_VERY_LONG_STRUCTURE_IDENTIFIER_INDEED",
        "X",
    ],
)
def test_structure_id_accepts_any_length(value):
    """PDB identifiers are growing beyond four characters (extended PDB IDs).
    No length is validated or assumed."""
    assert StructureId(value).value == value


def test_structure_id_parse_normalizes_without_truncating():
    assert StructureId.parse("  pdb_0000001ubq  ").value == "PDB_0000001UBQ"
    assert len(StructureId.parse("pdb_0000001ubq").value) == 14


@pytest.mark.parametrize("bad", ["", "1U BQ", "a/b", "a\\b", " 1UBQ", "1UBQ "])
def test_structure_id_rejects_unusable_values(bad):
    with pytest.raises(ValueError):
        StructureId(bad)


def test_structure_ids_are_comparable_and_hashable():
    ids = {StructureId("1UBQ"), StructureId("1UBQ"), StructureId("6G9M")}
    assert len(ids) == 2
    assert sorted(ids)[0] == StructureId("1UBQ")


# ---------------------------------------------------------------------------
# ChainId
# ---------------------------------------------------------------------------


def test_blank_chain_is_legal_and_distinct_from_a_named_one():
    assert ChainId("").is_blank
    assert ChainId("") != ChainId("A")


def test_chain_id_allows_multi_character_labels():
    """PDB gives one column; mmCIF does not, so no length is enforced."""
    assert ChainId("AAA").value == "AAA"


# ---------------------------------------------------------------------------
# ResidueIdentity
# ---------------------------------------------------------------------------


def test_insertion_code_is_a_separate_field():
    r = ResidueIdentity(chain=ChainId("A"), name="LEU", seq=15, insertion_code="B")
    assert r.seq == 15
    assert r.insertion_code == "B"
    assert r.legacy_seq_label == "15B"


def test_legacy_seq_label_without_an_insertion_code():
    r = ResidueIdentity(chain=ChainId("A"), name="LEU", seq=15)
    assert r.legacy_seq_label == "15"


def test_negative_residue_numbers_round_trip():
    r = ResidueIdentity(chain=ChainId("A"), name="MET", seq=-3)
    assert r.legacy_seq_label == "-3"


# ---------------------------------------------------------------------------
# AtomIdentity -- the load-bearing requirement
# ---------------------------------------------------------------------------


def _atom(**kwargs):
    base = dict(
        structure=StructureId("1KX5"),
        model=1,
        chain=ChainId("A"),
        residue_name="LYS",
        residue_seq=24,
        insertion_code="",
        atom_name="NZ",
        alt_loc="",
    )
    base.update(kwargs)
    return AtomIdentity(**base)


def test_same_atom_name_in_different_residues_is_a_different_atom():
    """The exact collision that under-counts triangles by 44% in the legacy cross
    pools: LYS 24 NZ and LYS 28 NZ are two physical atoms, not one."""
    a = _atom(residue_seq=24)
    b = _atom(residue_seq=28)
    assert a != b
    assert a.as_tuple() != b.as_tuple()
    assert len({a, b}) == 2


def test_every_identity_component_discriminates():
    base = _atom()
    variants = [
        _atom(structure=StructureId("1TUP")),
        _atom(model=2),
        _atom(chain=ChainId("B")),
        _atom(residue_name="ARG"),
        _atom(residue_seq=99),
        _atom(insertion_code="A"),
        _atom(atom_name="CA"),
        _atom(alt_loc="B"),
    ]
    for v in variants:
        assert v != base
        assert v.as_tuple() != base.as_tuple()
    assert len(set(variants) | {base}) == len(variants) + 1


def test_atom_name_alone_is_not_an_identity():
    a = _atom(residue_seq=24)
    b = _atom(residue_seq=28)
    assert a.atom_name == b.atom_name
    assert a != b, "identical names must not make identical atoms"


def test_identity_exposes_its_residue():
    a = _atom(residue_seq=24, insertion_code="C")
    r = a.residue
    assert (r.name, r.seq, r.insertion_code, r.chain) == ("LYS", 24, "C", ChainId("A"))


def test_describe_is_readable_and_complete():
    text = _atom(alt_loc="A").describe()
    for part in ("1KX5", "A", "LYS", "24", "NZ", "altloc=A"):
        assert part in text


@pytest.mark.parametrize(
    "kwargs",
    [
        {"atom_name": ""},
        {"atom_name": " CA "},
        {"alt_loc": "AB"},
        {"residue_name": " LEU"},
    ],
)
def test_malformed_identities_are_rejected(kwargs):
    with pytest.raises(ValueError):
        _atom(**kwargs)


# ---------------------------------------------------------------------------
# Point / PointSet
# ---------------------------------------------------------------------------


def test_distance_matches_the_legacy_expression():
    a, b = Point(1.0, 2.0, 3.0), Point(4.0, 6.0, 3.0)
    assert a.distance_to(b) == 5.0
    assert a.distance_to(b) == b.distance_to(a)


def test_midpoint():
    m = Point(0.0, 0.0, 0.0).midpoint_to(Point(2.0, 4.0, 6.0))
    assert m.as_tuple() == (1.0, 2.0, 3.0)


def test_pointset_preserves_order_and_is_immutable():
    pts = PointSet(
        [
            LabeledPoint(identity=_atom(residue_seq=i), coord=Point(float(i), 0, 0), label=8)
            for i in range(5)
        ]
    )
    assert len(pts) == 5
    assert [p.coord.x for p in pts] == [0.0, 1.0, 2.0, 3.0, 4.0]
    assert pts.labels == (8,) * 5
    with pytest.raises(TypeError):
        pts[0] = None  # type: ignore[index]


# ---------------------------------------------------------------------------
# Candidates and results
# ---------------------------------------------------------------------------


def test_candidate_indices_must_be_ascending():
    TriangleCandidate(0, 1, 2)
    for bad in [(1, 0, 2), (0, 2, 1), (0, 0, 1), (2, 1, 0)]:
        with pytest.raises(ValueError):
            TriangleCandidate(*bad)


def test_category_counts_default_bucket():
    c = CategoryCounts(by_category={DEFAULT_CATEGORY: 7})
    assert c.total == 7
    assert DEFAULT_CATEGORY == "regular"


def test_category_counts_support_multiple_buckets_without_a_schema_change():
    """A plain run writes one bucket; mirror and SSE need no new type."""
    c = CategoryCounts(by_category={"regular": 7, "mirror": 3})
    assert c.total == 10
