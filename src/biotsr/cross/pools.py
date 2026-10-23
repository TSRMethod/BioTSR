"""Cross contact-pool construction.

    side-A selection ─┐
                      ├─> contact scan (cutoff) ─> pooled, side-tagged PointSet
    side-B selection ─┘

Three decisions live here and nowhere else, and each is a recorded fact rather
than a default:

**1. What a contact pool member *is* (``PoolIdentity``).**  Pool membership is
keyed by the full physical ``AtomIdentity``: structure, model, chain, residue
name, residue sequence number, insertion code, atom name, altLoc.  A physical
atom that appears in several contact pairs contributes **one** pool member;
two different physical atoms never share one.  (Two historical cross programs
keyed pools by side and atom name only, so ``LYS 24 NZ`` and ``LYS 28 NZ``
collapsed into one slot -- on 1KX5's ``50_1KX5_J_DT`` target, 48 triangles
instead of 70.  That defect is not reproduced.)

**2. Which side the cutoff narrows (``ContactFilter``).**  The standard narrows
both sides symmetrically; ``side_a_only`` / ``side_b_only`` narrow one.

**3. Pool ordering.**  Side A's members, then side B's, each in selection (file)
order (``BLOCKED_BY_SIDE``).  Order fixes candidate order and therefore the
triplet file's layout.

**The cutoff is pool membership, not a triangle-edge constraint.**  An atom
enters because it is within the cutoff of *one* atom on the other side.  The
triangle built from it may then have edges far longer than the cutoff, and that
is correct -- no implementation in the family re-checks edges.  Size
filtering is a separate, optional modifier that acts on ``max_distance`` after
geometry; the two are never conflated.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from ..domain.identity import AtomIdentity
from ..domain.point import LabeledPoint, PointSet
from ..domain.result import SkipRecord
from ..geometry.conventions import SQRT_MATH_SQRT, SqrtConvention

__all__ = [
    "SIDE_A",
    "SIDE_B",
    "PoolIdentity",
    "ContactFilter",
    "CONTACT_FILTER_RULES",
    "PoolOrder",
    "PoolProjection",
    "ContactPoolPolicy",
    "CrossPools",
    "build_contact_pools",
    "build_interface_calpha_pools",
    "STANDARD_SYMMETRIC_ATOM_POOLS",
    "STANDARD_INTERFACE_CALPHA_POOLS",
]

#: The two side tags.  One spelling, used in code, provenance and output.
SIDE_A = "A"
SIDE_B = "B"


class PoolIdentity:
    """What counts as "the same atom" when a contact pool is deduplicated."""

    #: structure / model / chain / residue name / residue seq / iCode / atom / altLoc.
    #: The only identity BioTSR implements. (Two historical cross programs keyed
    #: their pools by ``(side, atom_name)`` and so collapsed distinct atoms; that
    #: defect is not reproduced.)
    FULL_ATOM_IDENTITY = "full_atom_identity"


class ContactFilter:
    """Which sides the contact cutoff narrows.

    Named for the side that *is* narrowed, because that is the side whose pool
    shrinks and the one a reader is checking.
    """

    BOTH_SIDES = "both_sides"
    #: Only side A is narrowed by the cutoff -- side B enters whole.
    SIDE_A_ONLY = "side_a_only"
    #: The mirror image, for a family that lists its sides the other way round.
    SIDE_B_ONLY = "side_b_only"
    #: No contact filtering at all -- both selections are used whole.
    NONE = "none"


#: The scientific spelling of each mechanism value, for ``SamplingProfile``.
CONTACT_FILTER_RULES = {
    ContactFilter.BOTH_SIDES: "symmetric_atom_contact_filter",
    ContactFilter.SIDE_A_ONLY: "one_side_only_contact_filter_side_a",
    ContactFilter.SIDE_B_ONLY: "one_side_only_contact_filter_side_b",
    ContactFilter.NONE: "no_contact_filter",
}


class PoolProjection:
    """How a qualifying contact becomes a TSR vertex.

    The two are genuinely different science, not two spellings of one rule, so
    they are named rather than hidden behind a boolean.
    """

    #: The contacting heavy atom **is** the vertex. Every atomic cross method.
    #: An atom that makes no qualifying contact does not enter the pool, even
    #: when another atom of its own residue does.
    ATOM_LEVEL = "atom_level_contact_pool"
    #: A residue is classified as an interface residue when **any** pair of its
    #: heavy atoms and the opposite side's heavy atoms lies within the cutoff,
    #: and the residue is then represented by its Ca. The Ca itself need not be
    #: within the cutoff -- it usually is not, for a long side chain. Protein-Ca
    #: protein-protein cross only (ADR-052).
    HEAVY_ATOM_INTERFACE_TO_CALPHA = "heavy_atom_interface_to_ca_projection"


class PoolOrder:
    """How the two pools are laid out inside the pooled point set."""

    #: Side A's pool, then side B's, each in selection (file) order. The only
    #: layout implemented.
    BLOCKED_BY_SIDE = "blocked_by_side_selection_order"


@dataclass(frozen=True)
class ContactPoolPolicy:
    """Every pool-construction decision of one lineage, named.

    ``cutoff`` is in angstroms and is compared with ``<=``, matching all four
    legacy implementations.  ``None`` means no contact scan at all: both
    selections enter whole.  It changes which points were sampled, never what a
    key means, so it belongs in ``SamplingProfile`` -- and a test asserts that
    changing it leaves ``KeySpace`` untouched.
    """

    name: str
    cutoff: Optional[float]
    identity: str = PoolIdentity.FULL_ATOM_IDENTITY
    contact_filter: str = ContactFilter.BOTH_SIDES
    order: str = PoolOrder.BLOCKED_BY_SIDE
    #: How a qualifying contact becomes a vertex. See ``PoolProjection``.
    projection: str = PoolProjection.ATOM_LEVEL
    #: Which side is the OUTER loop of the contact scan. It decides nothing about
    #: which pairs qualify -- the relation is symmetric -- but it fixes the
    #: first-contact ordering of the pools and, for the interleaved layout, the
    #: order of the single pool. ``Drug_TSR.py`` loops the ligand outside the
    #: protein; ``Cross_TSR.py`` loops chain 1 outside chain 2; both contact
    #: generators loop the entity side outside the protein side.
    scan_outer_side: str = SIDE_A
    sqrt: SqrtConvention = SQRT_MATH_SQRT

    def __post_init__(self) -> None:
        if self.cutoff is not None:
            # The cutoff is a USER-SPECIFIED distance in angstroms, so it arrives
            # from wherever the caller got it -- a literal, a CSV cell, a config
            # file. Normalise it to float here, at the one place every
            # construction path passes through.
            #
            # This is not cosmetic. `SamplingProfile.contact_cutoff` is
            # serialized into a content fingerprint, and `3` and `3.0` serialize
            # differently: without this, two runs that produce byte-identical
            # keys and frequencies would be refused as "not sampled the same
            # way" purely because one caller wrote an int. A false incomparable
            # is exactly the failure this identity exists to prevent.
            if isinstance(self.cutoff, bool) or not isinstance(self.cutoff, (int, float)):
                raise TypeError(
                    "contact cutoff must be a number of angstroms or None, got "
                    "%r (%s)" % (self.cutoff, type(self.cutoff).__name__)
                )
            object.__setattr__(self, "cutoff", float(self.cutoff))
        if self.identity != PoolIdentity.FULL_ATOM_IDENTITY:
            raise ValueError("unsupported pool identity %r; BioTSR keys pools by "
                             "full physical atom identity" % (self.identity,))
        if self.order != PoolOrder.BLOCKED_BY_SIDE:
            raise ValueError("unsupported pool order %r; BioTSR lays pools out "
                             "by side in selection order" % (self.order,))
        if self.cutoff is not None and not (self.cutoff > 0):
            raise ValueError("contact cutoff must be positive or None, got %r" % (self.cutoff,))
        if self.cutoff is None and self.contact_filter != ContactFilter.NONE:
            raise ValueError(
                "policy %r has no cutoff but asks to filter %s; pass "
                "contact_filter=ContactFilter.NONE to use both selections whole"
                % (self.name, self.contact_filter)
            )
        if self.cutoff is not None and self.contact_filter == ContactFilter.NONE:
            raise ValueError(
                "policy %r sets a cutoff of %r but filters nothing; drop the "
                "cutoff or name a side to filter" % (self.name, self.cutoff)
            )

    @property
    def contact_filter_rule(self) -> str:
        """The scientific name of this policy's contact rule."""
        return CONTACT_FILTER_RULES[self.contact_filter]

    def to_dict(self):
        return {
            "name": self.name,
            "contact_cutoff": self.cutoff,
            "pool_identity": self.identity,
            "contact_filter": self.contact_filter,
            "contact_filter_rule": self.contact_filter_rule,
            "contact_projection": self.projection,
            "pool_order": self.order,
            "scan_outer_side": self.scan_outer_side,
        }


@dataclass(frozen=True)
class CrossPools:
    """The pooled, side-tagged points plus the ledger of how they were chosen."""

    points: PointSet
    policy: ContactPoolPolicy
    side_a_offered: int
    side_b_offered: int
    side_a_pool: int
    side_b_pool: int
    n_contact_pairs: int
    #: Overwriting events that displaced a *different* physical atom.
    n_identity_collisions: int = 0
    #: How many distinct physical atoms the pool lost: the sum over slots of
    #: (atoms that occupied it) - 1. This is the number that explains the pool
    #: size difference between the legacy and physical identities.
    n_atoms_lost_to_identity: int = 0
    n_duplicate_records: int = 0
    #: Ca projection only: interface residues that carried no usable Ca and so
    #: could not be represented. Recorded rather than silently substituted with
    #: some other atom (ADR-052).
    n_interface_residues_without_calpha: int = 0
    skipped: Tuple[SkipRecord, ...] = ()

    @property
    def sides(self) -> Tuple[str, ...]:
        return tuple(p.side for p in self.points)

    def to_dict(self):
        d = dict(self.policy.to_dict())
        d.update({
            "side_a_offered": self.side_a_offered,
            "side_b_offered": self.side_b_offered,
            "side_a_pool": self.side_a_pool,
            "side_b_pool": self.side_b_pool,
            "n_contact_pairs": self.n_contact_pairs,
            "n_identity_collisions": self.n_identity_collisions,
            "n_atoms_lost_to_identity": self.n_atoms_lost_to_identity,
            "n_duplicate_records": self.n_duplicate_records,
            "n_interface_residues_without_calpha":
                self.n_interface_residues_without_calpha,
        })
        return d


def _identity_key(policy: ContactPoolPolicy, side: str, identity: AtomIdentity):
    return (side,) + tuple(identity.as_tuple())


def _contact_pairs(points_a, points_b, policy):
    """Qualifying ``(index_a, index_b)`` pairs, in the policy's scan order.

    Which side is outer changes nothing about which pairs qualify -- "within the
    cutoff" is symmetric -- but it fixes the first-contact ordering of the pools,
    and that ordering is what the historical triplet files and first-seen key
    order depend on.  All four legacy implementations are transcribed by
    ``scan_outer_side``.
    """
    sqrt = policy.sqrt
    cutoff = policy.cutoff
    if policy.scan_outer_side == SIDE_A:
        for i, pa in enumerate(points_a):
            for j, pb in enumerate(points_b):
                # The legacy code compares the rooted distance, so root it; `<=`
                # is what all four implementations write.
                if sqrt(pa.coord.squared_distance_to(pb.coord)) <= cutoff:
                    yield i, j
    else:
        for j, pb in enumerate(points_b):
            for i, pa in enumerate(points_a):
                if sqrt(pa.coord.squared_distance_to(pb.coord)) <= cutoff:
                    yield i, j


def build_contact_pools(
    points_a: Sequence[LabeledPoint],
    points_b: Sequence[LabeledPoint],
    policy: ContactPoolPolicy,
) -> CrossPools:
    """Build the two cross contact pools and pool them into one tagged point set.

    ``points_a`` and ``points_b`` are whatever the two side selectors returned,
    already labelled.  Their ``side`` field is overwritten here, so a selector
    never has to know which side it is being used for -- which is what lets the
    drug, carbohydrate, nucleotide and protein selectors be reused unchanged.
    """
    a_in = [p for p in points_a]
    b_in = [p for p in points_b]

    skipped: List[SkipRecord] = []
    n_pairs = 0

    # --- which members each side contributes, and in what order ----------------
    if policy.cutoff is None:
        keep_a = list(range(len(a_in)))
        keep_b = list(range(len(b_in)))
        scan: List[Tuple[int, int]] = []
    else:
        scan = list(_contact_pairs(a_in, b_in, policy))
        n_pairs = len(scan)
        seen_a, seen_b = [], []
        in_a, in_b = set(), set()
        for i, j in scan:
            if i not in in_a:
                in_a.add(i)
                seen_a.append(i)
            if j not in in_b:
                in_b.add(j)
                seen_b.append(j)
        keep_a = (
            list(range(len(a_in)))
            if policy.contact_filter == ContactFilter.SIDE_B_ONLY
            else seen_a
        )
        keep_b = (
            list(range(len(b_in)))
            if policy.contact_filter == ContactFilter.SIDE_A_ONLY
            else seen_b
        )
        keep_a = sorted(keep_a)
        keep_b = sorted(keep_b)

    # --- dedup under the policy's identity, and count what it cost -------------
    slots: "Dict[object, int]" = {}
    members: List[LabeledPoint] = []
    # Every distinct physical atom that has occupied each slot. Its size tells how
    # many atoms a slot collapsed, which is the honest measure of identity loss;
    # `n_identity_collisions` counts overwriting *events*, which can exceed it
    # when a slot is written by A, then B, then A again.
    physical_of_slot: "Dict[object, set]" = {}
    n_collisions = 0
    n_duplicates = 0

    def admit(point: LabeledPoint, side: str) -> None:
        nonlocal n_collisions, n_duplicates
        key = _identity_key(policy, side, point.identity)
        physical = point.identity.as_tuple()
        member = LabeledPoint(
            identity=point.identity, coord=point.coord,
            label=point.label, side=side,
        )
        if key in slots:
            index = slots[key]
            previous = members[index]
            if previous.identity.as_tuple() == physical:
                # The same physical atom, offered again. Collapsing it is correct
                # under BOTH identities and is not the defect.
                n_duplicates += 1
            else:
                # Two DIFFERENT physical atoms in one slot. Only reachable under
                # the legacy identity; the loss is recorded, never silent.
                n_collisions += 1
                physical_of_slot.setdefault(key, set()).add(previous.identity.as_tuple())
                physical_of_slot[key].add(physical)
                skipped.append(
                    SkipRecord(
                        stage="cross_pool",
                        reason="atom_identity_overwritten",
                        detail=(
                            "%s displaced %s in pool slot %r under pool identity "
                            "%r; distinct physical atoms are being collapsed"
                            % (
                                point.identity.describe(),
                                previous.identity.describe(),
                                key,
                                policy.identity,
                            )
                        ),
                    )
                )
            # Last write wins, unconditionally -- `getAxisSeq[k] = v` in
            # `cross_key.py` re-assigns the slot on EVERY contact row, so a later
            # row restoring an earlier atom is itself an overwrite. Short-circuiting
            # on "we have seen this physical atom before" would keep the wrong
            # atom, which is a real difference: on 7LF9's 1_7LF9_C_NAG target the
            # OD1 slot is written by ASN297, then ASP265, then ASN297 again.
            members[index] = member
            return
        physical_of_slot.setdefault(key, set()).add(physical)
        slots[key] = len(members)
        members.append(member)

    for i in keep_a:
        admit(a_in[i], SIDE_A)
    for j in keep_b:
        admit(b_in[j], SIDE_B)

    pool_a = sum(1 for p in members if p.side == SIDE_A)
    pool_b = sum(1 for p in members if p.side == SIDE_B)
    lost = sum(len(v) - 1 for v in physical_of_slot.values() if len(v) > 1)

    return CrossPools(
        points=PointSet(members),
        policy=policy,
        side_a_offered=len(a_in),
        side_b_offered=len(b_in),
        side_a_pool=pool_a,
        side_b_pool=pool_b,
        n_contact_pairs=n_pairs,
        n_identity_collisions=n_collisions,
        n_atoms_lost_to_identity=lost,
        n_duplicate_records=n_duplicates,
        skipped=tuple(skipped),
    )


def _residue_label(residue) -> str:
    """``chain/NAME/seq[iCode]`` -- the residue-level counterpart of
    ``AtomIdentity.describe()``."""
    return "%s/%s/%s" % (
        residue.chain.value if hasattr(residue.chain, "value") else residue.chain,
        residue.name,
        residue.legacy_seq_label,
    )


def build_interface_calpha_pools(
    heavy_a: Sequence[LabeledPoint],
    heavy_b: Sequence[LabeledPoint],
    calpha_a: Sequence[LabeledPoint],
    calpha_b: Sequence[LabeledPoint],
    policy: ContactPoolPolicy,
) -> CrossPools:
    """Interface residues by heavy-atom contact, represented by their Ca.

    The standard rule for Ca protein-protein cross, in two steps that
    must not be collapsed into one (ADR-052):

    1. residue *R* on side A and residue *S* on side B are interface residues
       when **at least one** pair of their non-hydrogen atoms is within the
       cutoff -- so membership is decided over all heavy atoms, side chains
       included;
    2. each interface residue then contributes its **Ca**, and that Ca does not
       itself have to be within the cutoff.  For a long side chain it usually is
       not: a LYS NZ can reach the opposite chain while its Ca sits 8 A away.

    Testing the Ca distance instead would be a different, stricter method, and
    would silently drop exactly the long-side-chain contacts that define most
    protein-protein interfaces.  Which is why the two steps are separate code
    and separate names.

    ``heavy_*`` decide membership; ``calpha_*`` supply the points and their
    labels (the residue label space at radix 30 -- a Ca point *is* a
    residue).  An interface residue with no usable Ca is recorded, never
    replaced by another atom.
    """
    if policy.projection != PoolProjection.HEAVY_ATOM_INTERFACE_TO_CALPHA:
        raise ValueError(
            "build_interface_calpha_pools needs projection=%r, got %r"
            % (PoolProjection.HEAVY_ATOM_INTERFACE_TO_CALPHA, policy.projection)
        )
    if policy.identity != PoolIdentity.FULL_ATOM_IDENTITY:
        raise ValueError(
            "the Ca interface projection has no legacy counterpart and so no "
            "legacy identity: pool identity must be %r, got %r"
            % (PoolIdentity.FULL_ATOM_IDENTITY, policy.identity)
        )

    skipped: List[SkipRecord] = []
    ca_of: Dict[str, Dict[object, LabeledPoint]] = {SIDE_A: {}, SIDE_B: {}}
    order_of: Dict[str, Dict[object, int]] = {SIDE_A: {}, SIDE_B: {}}
    for side, points in ((SIDE_A, calpha_a), (SIDE_B, calpha_b)):
        for n, point in enumerate(points):
            key = point.identity.residue
            if key in ca_of[side]:
                # Two Ca records for one residue -- an altLoc the selector let
                # through. First wins, and the displacement is recorded.
                skipped.append(
                    SkipRecord(
                        stage="cross_pool",
                        reason="duplicate_calpha_for_residue",
                        detail="%s duplicates the Ca already held for %s"
                        % (point.identity.describe(), _residue_label(key)),
                    )
                )
                continue
            ca_of[side][key] = point
            order_of[side][key] = n

    # --- step 1: interface residues, by ANY heavy-atom pair within the cutoff --
    sqrt = policy.sqrt
    cutoff = policy.cutoff
    if cutoff is None:
        raise ValueError(
            "the Ca interface projection is defined by a contact cutoff; "
            "policy %r has none" % (policy.name,)
        )
    # Dicts, not sets: insertion order is the heavy-atom scan order, which is
    # deterministic, while set iteration order over string-hashed keys is not.
    interface_a: "Dict[object, None]" = {}
    interface_b: "Dict[object, None]" = {}
    n_pairs = 0
    for pa in heavy_a:
        ra = pa.identity.residue
        for pb in heavy_b:
            if sqrt(pa.coord.squared_distance_to(pb.coord)) <= cutoff:
                n_pairs += 1
                interface_a.setdefault(ra, None)
                interface_b.setdefault(pb.identity.residue, None)

    # --- step 2: project each interface residue onto its Ca -------------------
    members: List[LabeledPoint] = []
    n_without_ca = 0

    def project(side: str, residues) -> None:
        nonlocal n_without_ca
        available = ca_of[side]
        chosen = []
        for residue in residues:
            point = available.get(residue)
            if point is None:
                n_without_ca += 1
                skipped.append(
                    SkipRecord(
                        stage="cross_pool",
                        reason="interface_residue_without_calpha",
                        detail=(
                            "%s is an interface residue on side %s but the Ca "
                            "selection holds no Ca for it, so it contributes no "
                            "point; no other atom is substituted"
                            % (_residue_label(residue), side)
                        ),
                    )
                )
                continue
            chosen.append((order_of[side][residue], point))
        # Selection (file) order, which is what makes candidate enumeration and
        # therefore serialization deterministic.
        for _, point in sorted(chosen, key=lambda pair: pair[0]):
            members.append(
                LabeledPoint(
                    identity=point.identity, coord=point.coord,
                    label=point.label, side=side,
                )
            )

    project(SIDE_A, interface_a)
    project(SIDE_B, interface_b)

    pool_a = sum(1 for p in members if p.side == SIDE_A)
    pool_b = sum(1 for p in members if p.side == SIDE_B)
    return CrossPools(
        points=PointSet(members),
        policy=policy,
        side_a_offered=len(calpha_a),
        side_b_offered=len(calpha_b),
        side_a_pool=pool_a,
        side_b_pool=pool_b,
        n_contact_pairs=n_pairs,
        n_interface_residues_without_calpha=n_without_ca,
        skipped=tuple(skipped),
    )


# ---------------------------------------------------------------------------
# Named policies
# ---------------------------------------------------------------------------


def STANDARD_SYMMETRIC_ATOM_POOLS(cutoff: float) -> ContactPoolPolicy:
    """The standard rule for every **atomic** cross method.

    An atom of either side enters the pool iff it has at least one non-hydrogen
    atom on the opposite specified side within the cutoff.  Both sides are
    narrowed -- the rule is symmetric -- and it acts at the **atom** level: an
    atom is not admitted merely because another atom of its residue or of its
    ligand makes contact (ADR-050, ADR-051).
    """
    return ContactPoolPolicy(
        name="standard_symmetric_atom_contact_pool",
        cutoff=cutoff,
        identity=PoolIdentity.FULL_ATOM_IDENTITY,
        contact_filter=ContactFilter.BOTH_SIDES,
        order=PoolOrder.BLOCKED_BY_SIDE,
        projection=PoolProjection.ATOM_LEVEL,
    )


def STANDARD_INTERFACE_CALPHA_POOLS(cutoff: float) -> ContactPoolPolicy:
    """The standard rule for **Ca** protein-protein cross (ADR-052).

    Interface residues by heavy-atom contact, then one Ca per interface residue.
    See :func:`build_interface_calpha_pools`.
    """
    return ContactPoolPolicy(
        name="standard_interface_calpha_pool",
        cutoff=cutoff,
        identity=PoolIdentity.FULL_ATOM_IDENTITY,
        contact_filter=ContactFilter.BOTH_SIDES,
        order=PoolOrder.BLOCKED_BY_SIDE,
        projection=PoolProjection.HEAVY_ATOM_INTERFACE_TO_CALPHA,
    )
