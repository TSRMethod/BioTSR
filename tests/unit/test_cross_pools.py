"""Cross contact pools: physical atom identity, and the cutoff's real meaning.

Two properties are load-bearing and both
are asserted with constructed, fully determined inputs rather than with a real
structure, so a failure points at the rule and not at a fixture:

* two different physical atoms that happen to share an atom name stay two atoms;
* the contact cutoff decides **pool membership** and is not an edge-length
  constraint on the resulting triangles.
"""

from __future__ import annotations

import pytest

from biotsr.cross.pools import (
    SIDE_A,
    SIDE_B,
    ContactFilter,
    ContactPoolPolicy,
    PoolIdentity,
    PoolOrder,
    build_contact_pools,
)
from biotsr.domain.identity import AtomIdentity, ChainId, StructureId
from biotsr.domain.point import LabeledPoint, Point


def atom(chain, resname, seq, name, *, alt_loc="", structure="TEST", model=1,
         icode=""):
    return AtomIdentity(
        structure=StructureId(structure),
        model=model,
        chain=ChainId(chain),
        residue_name=resname,
        residue_seq=seq,
        insertion_code=icode,
        atom_name=name,
        alt_loc=alt_loc,
    )


def point(identity, x, y, z, label=108):
    return LabeledPoint(identity=identity, coord=Point(x, y, z), label=label)


def standard(cutoff, **kw):
    return ContactPoolPolicy(name="test_standard", cutoff=cutoff, **kw)


# ===========================================================================
# Physical atom identity
# ===========================================================================


LYS24_NZ = atom("D", "LYS", 24, "NZ")
LYS28_NZ = atom("D", "LYS", 28, "NZ")


def two_lysines_and_one_entity_atom():
    """The exact 1KX5 shape: one entity atom in contact with LYS24 NZ and
    LYS28 NZ, two different atoms that share the name ``NZ``."""
    entity = [point(atom("J", "DT", 50, "O2"), 0.0, 0.0, 0.0, label=109)]
    protein = [
        point(LYS24_NZ, 2.4, 0.0, 0.0),
        point(LYS28_NZ, 0.0, 2.5, 0.0),
    ]
    return entity, protein


def test_same_named_atoms_in_different_residues_stay_distinct_in_standard():
    """Historical difference, recorded as a constant: the `cross_key.py` pools
    were keyed by ``(side, atom_name)``, so LYS24 NZ and LYS28 NZ shared one slot
    and the protein pool here was ``HISTORICAL_SIDE_AND_NAME_POOL_B`` = 1 (one
    atom lost to the identity, last writer kept).  The standard keeps both."""
    HISTORICAL_SIDE_AND_NAME_POOL_B = 1  # documentation only
    entity, protein = two_lysines_and_one_entity_atom()
    pools = build_contact_pools(entity, protein, standard(3.0))
    assert pools.side_a_pool == 1
    assert pools.side_b_pool == 2
    assert pools.n_identity_collisions == 0
    assert pools.n_atoms_lost_to_identity == 0
    kept = {p.identity for p in pools.points if p.side == SIDE_B}
    assert kept == {LYS24_NZ, LYS28_NZ}
    assert pools.side_b_pool != HISTORICAL_SIDE_AND_NAME_POOL_B


def test_the_same_pools_under_physical_identity_keep_both_atoms():
    """The 7LF9 ``1_7LF9_C_NAG`` shape: ASN297 OD1 and ASP265 OD1 both contact
    the glycan (three contact pairs: ASN, ASP, ASN again).  Under the standard
    identity nothing is displaced at all.

    Historical difference, recorded as a constant: under the `cross_key.py`
    ``(side, atom_name)`` identity the OD1 slot was written ASN, ASP, ASN and
    ended holding ASN alone -- ``HISTORICAL_SIDE_AND_NAME_POOL_B`` = 1 with two
    identity collisions.
    """
    HISTORICAL_SIDE_AND_NAME_POOL_B = 1  # documentation only
    asn = atom("A", "ASN", 297, "OD1")
    asp = atom("A", "ASP", 265, "OD1")
    entity = [
        point(atom("C", "NAG", 1, "C1"), 0.0, 1.0, 0.0, label=106),
        point(atom("C", "NAG", 1, "N2"), 10.0, 1.0, 0.0, label=108),
        point(atom("C", "NAG", 1, "O5"), 0.0, -1.0, 0.0, label=109),
    ]
    protein = [point(asn, 0.0, 0.0, 0.0), point(asp, 10.0, 0.0, 0.0)]
    pools = build_contact_pools(entity, protein, standard(2.0))
    assert pools.n_contact_pairs == 3
    assert pools.side_b_pool == 2
    assert pools.side_b_pool != HISTORICAL_SIDE_AND_NAME_POOL_B
    assert pools.n_identity_collisions == 0
    assert {p.identity for p in pools.points if p.side == SIDE_B} == {asn, asp}


def test_repeated_references_to_the_same_physical_atom_do_not_double_count():
    """Several contact rows for one atom collapse to one pool member."""
    nz = atom("D", "LYS", 24, "NZ")
    entity = [
        point(atom("J", "DT", 50, "O2"), 0.0, 0.0, 0.0),
        point(atom("J", "DT", 50, "O4"), 0.5, 0.0, 0.0),
        point(atom("J", "DT", 50, "OP1"), 1.0, 0.0, 0.0),
    ]
    protein = [point(nz, 0.5, 1.0, 0.0)]
    pools = build_contact_pools(entity, protein, standard(3.0))
    assert pools.n_contact_pairs == 3
    assert pools.side_b_pool == 1
    assert pools.n_identity_collisions == 0


def test_the_blocked_layout_records_no_duplicate_offers():
    """The blocked layout walks a deduplicated index list, so an atom in three
    contact pairs is offered once and there is no duplicate to count."""
    nz = atom("D", "LYS", 24, "NZ")
    entity = [
        point(atom("J", "DT", 50, "O2"), 0.0, 0.0, 0.0),
        point(atom("J", "DT", 50, "O4"), 0.5, 0.0, 0.0),
        point(atom("J", "DT", 50, "OP1"), 1.0, 0.0, 0.0),
    ]
    protein = [point(nz, 0.5, 1.0, 0.0)]
    assert build_contact_pools(entity, protein, standard(3.0)).n_duplicate_records == 0


@pytest.mark.parametrize(
    "other",
    [
        atom("E", "LYS", 24, "NZ"),                 # different chain
        atom("D", "LYS", 24, "NZ", alt_loc="B"),    # different conformer
        atom("D", "LYS", 24, "NZ", model=2),        # different model
        atom("D", "LYS", 24, "NZ", icode="A"),      # different insertion code
        atom("D", "ARG", 24, "NZ"),                 # different residue name
        atom("D", "LYS", 24, "NZ", structure="OTHER"),  # different structure
    ],
)
def test_every_identity_field_discriminates(other):
    """Full identity means full: each field alone keeps two atoms apart."""
    base = atom("D", "LYS", 24, "NZ")
    entity = [point(atom("J", "DT", 50, "O2"), 0.0, 0.0, 0.0, label=109)]
    protein = [point(base, 1.0, 0.0, 0.0), point(other, 0.0, 1.0, 0.0)]
    pools = build_contact_pools(entity, protein, standard(3.0))
    assert pools.side_b_pool == 2
    assert pools.n_identity_collisions == 0


def test_the_standard_identity_can_never_collide():
    """Not "does not happen to" -- cannot: the slot key IS the full identity."""
    entity = [point(atom("J", "DT", 50, "O2"), 0.0, 0.0, 0.0, label=109)]
    protein = [
        point(atom("D", "LYS", n, "NZ"), 0.1 * n, 1.0, 0.0) for n in range(1, 10)
    ]
    pools = build_contact_pools(entity, protein, standard(20.0))
    assert pools.side_b_pool == 9
    assert pools.n_identity_collisions == 0
    assert pools.skipped == ()


# ===========================================================================
# The contact cutoff is pool membership, not a triangle-edge constraint
# ===========================================================================


def test_a_triangle_may_contain_an_edge_far_longer_than_the_contact_cutoff():
    """One A-B contact satisfies the cutoff; another edge of the resulting
    mixed-side triangle does not.  The triangle is kept; edges are never
    re-checked against the cutoff."""
    a = [point(atom("J", "DT", 50, "O2"), 0.0, 0.0, 0.0, label=109)]
    b = [
        point(atom("D", "LYS", 24, "NZ"), 2.5, 0.0, 0.0),    # 2.5 A: in contact
        point(atom("D", "LYS", 28, "NZ"), 2.0, 2.0, 0.0),    # 2.83 A: in contact
    ]
    # ... but the B-B edge is 2.06 A; widen it so it exceeds the cutoff.
    b[1] = point(atom("D", "LYS", 28, "NZ"), -2.9, 0.4, 0.0)
    pools = build_contact_pools(a, b, standard(3.0))
    assert pools.side_a_pool == 1 and pools.side_b_pool == 2

    from math import dist

    coords = [p.coord.as_tuple() for p in pools.points]
    edges = [
        dist(coords[0], coords[1]), dist(coords[0], coords[2]),
        dist(coords[1], coords[2]),
    ]
    assert max(edges) > 3.0, "the test needs an edge longer than the cutoff"
    assert min(edges) <= 3.0

    from biotsr.cross.candidates import count, mixed_side

    assert len(list(mixed_side(pools.sides))) == count(1, 2) == 1


def test_an_atom_out_of_contact_never_enters_the_pool():
    """The other half of the same rule: membership really is decided by it."""
    a = [point(atom("J", "DT", 50, "O2"), 0.0, 0.0, 0.0, label=109)]
    b = [
        point(atom("D", "LYS", 24, "NZ"), 2.5, 0.0, 0.0),
        point(atom("D", "LYS", 28, "NZ"), 50.0, 0.0, 0.0),
    ]
    pools = build_contact_pools(a, b, standard(3.0))
    assert pools.side_b_pool == 1
    assert [p.identity for p in pools.points if p.side == SIDE_B] == [
        atom("D", "LYS", 24, "NZ")
    ]


def test_the_cutoff_boundary_is_inclusive():
    """The comparison is ``<=`` (as in every historical implementation)."""
    a = [point(atom("J", "DT", 50, "O2"), 0.0, 0.0, 0.0, label=109)]
    b = [point(atom("D", "LYS", 24, "NZ"), 3.0, 0.0, 0.0)]
    assert build_contact_pools(a, b, standard(3.0)).side_b_pool == 1
    assert build_contact_pools(a, b, standard(2.9999)).side_b_pool == 0


def test_contact_filtering_can_narrow_one_side_only():
    """An explicit one-sided override: the whole ligand takes part, the protein
    is narrowed."""
    entity = [
        point(atom("A", "ESW", 404, "C%d" % n), 0.0, 0.0, 3.0 * n, label=106)
        for n in range(1, 6)
    ]
    protein = [
        point(atom("A", "LYS", 207, "NZ"), 1.0, 0.0, 3.0),
        point(atom("A", "ASP", 30, "OD2"), 0.0, 0.0, 400.0),
    ]
    narrowed = build_contact_pools(
        protein, entity,
        ContactPoolPolicy(name="t", cutoff=3.5,
                          contact_filter=ContactFilter.SIDE_A_ONLY),
    )
    assert narrowed.side_a_pool == 1       # only the contacting protein atom
    assert narrowed.side_b_pool == 5       # every ligand atom, contact or not

    both = build_contact_pools(protein, entity, standard(3.5))
    assert both.side_a_pool == 1
    assert both.side_b_pool < 5


def test_no_contact_filter_uses_both_selections_whole():
    entity = [point(atom("A", "ESW", 404, "C1"), 0.0, 0.0, 0.0, label=106)]
    protein = [point(atom("A", "LYS", 207, "NZ"), 0.0, 0.0, 900.0)]
    pools = build_contact_pools(
        entity, protein,
        ContactPoolPolicy(name="none", cutoff=None,
                          contact_filter=ContactFilter.NONE),
    )
    assert pools.side_a_pool == 1 and pools.side_b_pool == 1
    assert pools.n_contact_pairs == 0


def test_a_cutoff_without_a_filter_is_refused_rather_than_ignored():
    with pytest.raises(ValueError) as excinfo:
        ContactPoolPolicy(name="x", cutoff=3.0, contact_filter=ContactFilter.NONE)
    assert "filters nothing" in str(excinfo.value)


def test_a_filter_without_a_cutoff_is_refused_rather_than_ignored():
    with pytest.raises(ValueError) as excinfo:
        ContactPoolPolicy(name="x", cutoff=None,
                          contact_filter=ContactFilter.BOTH_SIDES)
    assert "no cutoff" in str(excinfo.value)


@pytest.mark.parametrize("bad", [0, -1.0, -0.0])
def test_a_non_positive_cutoff_is_refused(bad):
    with pytest.raises(ValueError):
        ContactPoolPolicy(name="x", cutoff=bad)


# ===========================================================================
# Pool ordering
# ===========================================================================


def test_blocked_ordering_puts_side_a_first_then_side_b():
    a = [point(atom("J", "DT", 50, "O%d" % n), 0.0, 0.0, 0.0, label=109)
         for n in range(1, 4)]
    b = [point(atom("D", "LYS", 20 + n, "NZ"), 1.0, 0.0, 0.0) for n in range(3)]
    pools = build_contact_pools(a, b, standard(5.0))
    assert list(pools.sides) == [SIDE_A] * 3 + [SIDE_B] * 3


def test_pool_construction_is_deterministic():
    a = [point(atom("J", "DT", 50, "O%d" % n), 0.1 * n, 0.0, 0.0, label=109)
         for n in range(1, 6)]
    b = [point(atom("D", "LYS", 20 + n, "NZ"), 0.0, 1.0, 0.1 * n) for n in range(5)]
    first = build_contact_pools(a, b, standard(5.0))
    second = build_contact_pools(a, b, standard(5.0))
    assert [p.identity for p in first.points] == [p.identity for p in second.points]


# ===========================================================================
# Only the standard identity and order can be constructed
# ===========================================================================


def test_the_default_policy_is_the_standard_identity_and_order():
    policy = standard(3.5)
    assert policy.identity == PoolIdentity.FULL_ATOM_IDENTITY
    assert policy.order == PoolOrder.BLOCKED_BY_SIDE


@pytest.mark.parametrize(
    "identity",
    ["legacy_side_and_atom_name", "side_and_atom_name", "", "FULL_ATOM_IDENTITY"],
)
def test_a_policy_with_a_non_standard_identity_cannot_be_constructed(identity):
    with pytest.raises(ValueError) as excinfo:
        ContactPoolPolicy(name="x", cutoff=3.5, identity=identity)
    assert "pool identity" in str(excinfo.value)


@pytest.mark.parametrize(
    "order",
    ["interleaved_by_contact_scan", "blocked_by_contact_scan", ""],
)
def test_a_policy_with_a_non_standard_order_cannot_be_constructed(order):
    with pytest.raises(ValueError) as excinfo:
        ContactPoolPolicy(name="x", cutoff=3.5, order=order)
    assert "pool order" in str(excinfo.value)


@pytest.mark.parametrize(
    "removed",
    ["LEGACY_SIDE_AND_ATOM_NAME"],
)
def test_the_removed_pool_identity_is_not_offered(removed):
    assert not hasattr(PoolIdentity, removed)


@pytest.mark.parametrize(
    "removed",
    ["BLOCKED_BY_CONTACT_SCAN", "INTERLEAVED_BY_CONTACT_SCAN"],
)
def test_the_removed_pool_orders_are_not_offered(removed):
    assert not hasattr(PoolOrder, removed)
