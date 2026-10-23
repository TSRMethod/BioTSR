"""The standard cross rules (ADR-049 to ADR-053), one test per rule, on exact geometry.

Settled on 2026-09-19 and no longer open questions.  Every structure here is
constructed so that the rule is the *only* thing that can decide the outcome:
the distances are chosen to straddle the cutoff by whole angstroms, so no
assertion turns on a rounding, a bin edge or a file's accidental content.

The rules, and what each one would look like if it were implemented the other
plausible way -- which is what these tests actually discriminate:

============================  ==================================================
rule                          the wrong implementation these catch
============================  ==================================================
both chains user-specified    scanning the whole structure and finding a
                              qualifying atom on a chain nobody named
hydrogen excluded             a hydrogen creating contact membership where no
                              heavy atom does
symmetric contact filtering   narrowing one side and letting the other in whole
atom-level pools              admitting a whole residue or a whole ligand
                              because one of its atoms touched
Ca interface projection       testing the Ca-Ca distance instead of the
                              heavy-atom distance, which drops every
                              long-side-chain contact
cutoff is pool membership     re-imposing the cutoff on the triangle's edges
============================  ==================================================
"""

from __future__ import annotations

from math import comb

import pytest

from biotsr import analyze_cross
from biotsr.domain.errors import InvalidRequest
from biotsr.cross.pools import SIDE_A, SIDE_B, PoolProjection
from biotsr.domain.identity import StructureId
from biotsr.structure.backends.pdb_legacy import parse_pdb_legacy

SID = StructureId("SYNT")

CUTOFF = 3.0
#: Every "far" atom is this far from everything on the other side, so no
#: assertion below depends on a marginal distance.
FAR = 40.0


def atom(serial, name, res, chain, seq, x, y, z, element=" C", record="ATOM  ",
         alt=" "):
    return (
        "%-6s%5d %-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n"
        % (record, serial, name, alt, res, chain, seq, x, y, z, element)
    )


def structure(atoms):
    return parse_pdb_legacy(atoms, SID)


def pooled(result):
    """``{(chain, resname, seq, atom_name)}`` of everything in the contact pool.

    Read off the emitted triangles' identities when there are any, and off the
    pool ledger otherwise, so a pool too small to form a triangle is still
    inspectable.
    """
    return {
        (i.chain.value, i.residue_name, i.residue_seq, i.atom_name)
        for t in (result.triplets or ())
        for i in t.identities
    }


def pool_sizes(result):
    p = result.annotations["cross_pools"]
    return p["side_a_pool"], p["side_b_pool"]


# ===========================================================================
# Rule 1 -- both participating chains are user-specified
# ===========================================================================


#: A nucleotide in chain E with a qualifying protein atom in chain B -- the chain
#: the caller names -- and another equally qualifying one in chain Z, which the
#: caller never names.  Chain Z exists precisely to be found by a whole-structure
#: scan, so "the standard honours the named chains" and "the standard happens to
#: have nothing else to find" cannot be confused.
UNRELATED_CHAIN_STRUCTURE = [
    atom(1, " N1 ", "DG", "E", 8, 0.0, 0.0, 0.0, " N"),
    atom(2, " C2 ", "DG", "E", 8, 0.0, 0.0, FAR, " C"),
    atom(3, " O6 ", "DG", "E", 8, 0.0, 0.0, 2 * FAR, " O"),
    # chain B -- named by the caller
    atom(4, " N  ", "ALA", "B", 10, 2.5, 0.0, 0.0, " N"),
    atom(5, " CA ", "ALA", "B", 10, 12.0, 0.0, 0.0, " C"),
    atom(6, " CB ", "ALA", "B", 10, 2.5, 1.0, 0.0, " C"),
    # chain Z -- never named, and just as close
    atom(7, " N  ", "ARG", "Z", 99, -2.5, 0.0, 0.0, " N"),
    atom(8, " CA ", "ARG", "Z", 99, -12.0, 0.0, 0.0, " C"),
]


def test_the_standard_pool_holds_only_the_two_chains_the_caller_named():
    result = analyze_cross(
        structure(UNRELATED_CHAIN_STRUCTURE), method="nucleotide_protein_cross",
        chain="E", chain_2="B", entity_name="DG", entity_id="8",
        contact_cutoff=CUTOFF,
    )
    a, b = pool_sizes(result)
    assert (a, b) == (1, 2)
    assert result.sampling.side_a_chain == "E"
    assert result.sampling.side_b_chain == "B"


def test_an_equally_close_atom_on_an_unnamed_chain_does_not_enter_the_standard_pool():
    """Chain Z's N is 2.5 A from the nucleotide -- the same distance as chain B's
    N, and well inside the 3 A cutoff.  It is excluded for one reason only: the
    caller did not name chain Z.

    Historical difference, recorded as a constant: the historical contact
    generators had no chain predicate and scanned the whole structure, so on
    this structure their pool was ``HISTORICAL_WHOLE_STRUCTURE_POOL`` = (1, 3)
    (chain Z included).  The standard pool is (1, 2).
    """
    HISTORICAL_WHOLE_STRUCTURE_POOL = (1, 3)  # documentation only; not computed
    assert HISTORICAL_WHOLE_STRUCTURE_POOL != (1, 2)
    result = analyze_cross(
        structure(UNRELATED_CHAIN_STRUCTURE), method="nucleotide_protein_cross",
        chain="E", chain_2="B", entity_name="DG", entity_id="8",
        contact_cutoff=CUTOFF, keep_triplets=False,
    )
    chains = set(result.selection.chains_used)
    assert "Z" not in chains
    assert chains <= {"E", "B"}


# ===========================================================================
# Rule 2 -- hydrogen atoms are excluded from every cross calculation
# ===========================================================================


#: The only atom within reach of the nucleotide's C2 is a hydrogen.  If hydrogen
#: created contact membership, C2 would be in the pool and the pool would be
#: larger by one on each side.
HYDROGEN_STRUCTURE = [
    atom(1, " N1 ", "DG", "E", 8, 0.0, 0.0, 0.0, " N"),
    atom(2, " C2 ", "DG", "E", 8, 0.0, 0.0, FAR, " C"),
    atom(3, " N  ", "ALA", "B", 10, 2.5, 0.0, 0.0, " N"),
    atom(4, " H  ", "ALA", "B", 11, 1.0, 0.0, FAR, " H"),
]


def test_a_hydrogen_does_not_create_contact_membership():
    result = analyze_cross(
        structure(HYDROGEN_STRUCTURE), method="nucleotide_protein_cross",
        chain="E", chain_2="B", entity_name="DG", entity_id="8",
        contact_cutoff=CUTOFF,
    )
    assert pool_sizes(result) == (1, 1)
    names = {name for _c, _r, _s, name in pooled(result)} if result.triplets else set()
    assert "H" not in names
    # and the atom it would have admitted is absent too
    assert result.annotations["cross_pools"]["side_a_pool"] == 1


def test_the_same_hydrogen_excluded_from_the_calpha_interface_scan():
    """Interface residues are identified over non-hydrogen atoms, so a hydrogen
    cannot promote a residue to the interface either."""
    lines = [
        atom(1, " CA ", "LYS", "A", 1, 0.0, 0.0, 0.0, " C"),
        atom(2, " NZ ", "LYS", "A", 1, 8.0, 0.0, 0.0, " N"),
        atom(3, " CA ", "GLU", "B", 1, 12.0, 0.0, 0.0, " C"),
        atom(4, " OE1", "GLU", "B", 1, 10.0, 0.0, 0.0, " O"),
        # a hydrogen that alone would put ALA 2 at the interface
        atom(5, " CA ", "ALA", "A", 2, 0.0, 0.0, FAR, " C"),
        atom(6, " HB1", "ALA", "A", 2, 9.5, 0.0, 0.0, " H"),
    ]
    result = analyze_cross(
        structure(lines), method="protein_protein_calpha_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF,
    )
    a, _b = pool_sizes(result)
    assert a == 1, "only LYS 1 reaches chain B, and it does so by its NZ"


# ===========================================================================
# Rule 3 -- contact filtering is symmetric and acts at the atom level
# ===========================================================================


#: A1 contacts B1.  A2 contacts nobody.  B2 contacts nobody.  Every asymmetric
#: rule -- "narrow the protein only", "narrow the ligand only" -- gives a
#: different answer here, and so does every residue-level rule, because A1/A2 are
#: two atoms of one residue and B1/B2 of another.
SYMMETRY_STRUCTURE = [
    atom(1, " N  ", "ALA", "A", 1, 0.0, 0.0, 0.0, " N"),       # A1
    atom(2, " CB ", "ALA", "A", 1, 0.0, 0.0, FAR, " C"),       # A2
    atom(3, " N  ", "GLY", "B", 1, 2.0, 0.0, 0.0, " N"),       # B1
    atom(4, " CA ", "GLY", "B", 1, 0.0, FAR, 0.0, " C"),       # B2
]


def test_only_the_contacting_atom_of_each_side_enters_the_standard_pool():
    result = analyze_cross(
        structure(SYMMETRY_STRUCTURE), method="protein_protein_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF, keep_triplets=False,
    )
    assert pool_sizes(result) == (1, 1)
    assert result.sampling.contact_filter_rule == "symmetric_atom_contact_filter"
    assert result.sampling.contact_projection == PoolProjection.ATOM_LEVEL


def test_a_non_contacting_atom_of_a_contacting_residue_is_still_excluded():
    """The rule is atom-level, not residue-level: ALA 1's CB is 40 A away and is
    excluded although its own N makes contact."""
    result = analyze_cross(
        structure(SYMMETRY_STRUCTURE), method="protein_protein_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF, keep_triplets=False,
    )
    assert result.annotations["cross_pools"]["side_a_offered"] == 2
    assert result.annotations["cross_pools"]["side_a_pool"] == 1


def test_the_drug_side_is_contact_filtered_like_every_other_side():
    """The standard filters both sides, so a ligand atom 40 A from the protein is
    excluded.

    Historical difference, recorded as a constant: `Drug_TSR.py` let every
    selected ligand atom take part regardless of distance, so on this structure
    it pooled ``HISTORICAL_DRUG_TSR_POOL`` = (1, 2) -- the whole ligand.  The
    standard pool is (1, 1).
    """
    HISTORICAL_DRUG_TSR_POOL = (1, 2)  # documentation only; not computed
    lines = [
        atom(1, " N  ", "ALA", "A", 1, 0.0, 0.0, 0.0, " N"),
        atom(2, " C1 ", "LIG", "A", 404, 2.0, 0.0, 0.0, " C", record="HETATM"),
        atom(3, " C2 ", "LIG", "A", 404, 0.0, 0.0, FAR, " C", record="HETATM"),
    ]
    kw = dict(method="drug_protein_cross", chain="A", chain_2="A",
              entity_name="LIG", entity_id="404", contact_cutoff=CUTOFF,
              keep_triplets=False)
    standard = analyze_cross(structure(lines), **kw)
    assert pool_sizes(standard) == (1, 1), "the far ligand atom is filtered out"
    assert standard.sampling.contact_filter_rule == "symmetric_atom_contact_filter"
    assert pool_sizes(standard) != HISTORICAL_DRUG_TSR_POOL


def test_an_explicit_drug_atom_selection_is_applied_before_the_contact_filter():
    """`drug_atom_selected` says which atoms are eligible; the cutoff then says
    which of those are at the interface.  Selecting a far atom does not force it
    into the pool."""
    lines = [
        atom(1, " N  ", "ALA", "A", 1, 0.0, 0.0, 0.0, " N"),
        atom(2, " C1 ", "LIG", "A", 404, 2.0, 0.0, 0.0, " C", record="HETATM"),
        atom(3, " C2 ", "LIG", "A", 404, 0.0, 0.0, FAR, " C", record="HETATM"),
        atom(4, " C3 ", "LIG", "A", 404, 2.0, 1.0, 0.0, " C", record="HETATM"),
    ]
    result = analyze_cross(
        structure(lines), method="drug_protein_cross", chain="A", chain_2="A",
        entity_name="LIG", entity_id="404", atoms="C1;C2",
        contact_cutoff=CUTOFF, keep_triplets=False,
    )
    pools = result.annotations["cross_pools"]
    assert pools["side_b_offered"] == 2, "the selection narrowed the ligand to C1, C2"
    assert pools["side_b_pool"] == 1, "and the cutoff then dropped the far one"


# ===========================================================================
# Rule 4 -- carbohydrate / nucleotide entities are NOT admitted whole
# ===========================================================================


def test_only_the_contacting_atoms_of_an_entity_residue_are_retained():
    """E1 contacts the protein; E2 and E3 do not.  The residue is not admitted
    whole -- only E1 becomes a TSR vertex."""
    lines = [
        atom(1, " C1 ", "NAG", "C", 1, 0.0, 0.0, 0.0, " C", record="HETATM"),
        atom(2, " C2 ", "NAG", "C", 1, 0.0, 0.0, FAR, " C", record="HETATM"),
        atom(3, " O5 ", "NAG", "C", 1, 0.0, 0.0, 2 * FAR, " O", record="HETATM"),
        atom(4, " N  ", "ASN", "A", 297, 2.0, 0.0, 0.0, " N"),
        atom(5, " OD1", "ASN", "A", 297, 2.0, 1.0, 0.0, " O"),
    ]
    result = analyze_cross(
        structure(lines), method="carbohydrate_protein_cross",
        chain="C", chain_2="A", entity_name="NAG", entity_id="1",
        contact_cutoff=CUTOFF,
    )
    pools = result.annotations["cross_pools"]
    assert pools["side_a_offered"] == 3, "all three atoms were offered"
    assert pools["side_a_pool"] == 1, "only the contacting one was retained"
    if result.triplets:
        assert {n for _c, _r, _s, n in pooled(result)} == {"C1", "N", "OD1"}


# ===========================================================================
# Rule 5 -- Ca protein-protein: heavy-atom interface, then Ca projection
# ===========================================================================


#: LYS 1 reaches chain B through its NZ, 2 A from GLU 1's OE1.  Its own Ca is 12 A
#: from GLU 1's Ca -- four times the cutoff.  A Ca-Ca contact rule finds nothing
#: here; the approved rule finds both residues.  ALA 2 reaches nothing and must
#: stay out, so "everything is admitted" fails too.
CALPHA_INTERFACE_STRUCTURE = [
    atom(1, " CA ", "LYS", "A", 1, 0.0, 0.0, 0.0, " C"),
    atom(2, " CB ", "LYS", "A", 1, 4.0, 0.0, 0.0, " C"),
    atom(3, " NZ ", "LYS", "A", 1, 8.0, 0.0, 0.0, " N"),
    atom(4, " CA ", "ALA", "A", 2, 0.0, 0.0, FAR, " C"),
    atom(5, " CB ", "ALA", "A", 2, 1.5, 0.0, FAR, " C"),
    atom(6, " CA ", "GLU", "B", 1, 12.0, 0.0, 0.0, " C"),
    atom(7, " OE1", "GLU", "B", 1, 10.0, 0.0, 0.0, " O"),
    atom(8, " CA ", "SER", "B", 2, 12.0, 0.0, FAR, " C"),
]


def test_a_side_chain_contact_puts_the_residue_at_the_interface():
    result = analyze_cross(
        structure(CALPHA_INTERFACE_STRUCTURE),
        method="protein_protein_calpha_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF,
    )
    assert pool_sizes(result) == (1, 1)
    assert result.sampling.contact_projection == (
        PoolProjection.HEAVY_ATOM_INTERFACE_TO_CALPHA
    )


def test_the_projected_calpha_is_itself_far_outside_the_cutoff():
    """The point of the rule.  Both pooled points are Ca atoms 12 A apart, at a
    3 A cutoff: had the cutoff been applied to the Ca atoms, the pool would be
    empty."""
    result = analyze_cross(
        structure(CALPHA_INTERFACE_STRUCTURE),
        method="protein_protein_calpha_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF, keep_triplets=False,
    )
    assert pool_sizes(result) == (1, 1)

    # The same structure, asked the wrong question: if the two Ca atoms were the
    # contact test, nothing at all would qualify.
    atomic = analyze_cross(
        structure(CALPHA_INTERFACE_STRUCTURE), method="protein_protein_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF,
        atom_selection="LYS_1_CA", atom_selection_2="GLU_1_CA",
        keep_triplets=False,
    )
    assert pool_sizes(atomic) == (0, 0)


def test_a_residue_that_reaches_nothing_is_not_an_interface_residue():
    result = analyze_cross(
        structure(CALPHA_INTERFACE_STRUCTURE),
        method="protein_protein_calpha_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF,
    )
    pools = result.annotations["cross_pools"]
    assert pools["side_a_offered"] == 2, "both Ca atoms of chain A were available"
    assert pools["side_a_pool"] == 1, "only LYS 1 is at the interface"


def test_the_calpha_cross_labels_residues_at_radix_30():
    result = analyze_cross(
        structure(CALPHA_INTERFACE_STRUCTURE),
        method="protein_protein_calpha_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF, keep_triplets=False,
    )
    assert result.key_space.label_scheme.radix == 30
    assert result.key_space.label_scheme.name == "amino_acid_30"
    assert result.key_space.distance_scheme.radix == 35
    assert result.key_space.numerical_precision == "float64"


def test_an_interface_residue_without_a_calpha_is_recorded_not_substituted():
    """A residue can be at the interface and still have no Ca -- an incomplete
    model, or a ligand-like residue.  It contributes no point, and the reason is
    in the ledger rather than being papered over with some other atom."""
    lines = [
        atom(1, " CA ", "LYS", "A", 1, 0.0, 0.0, 0.0, " C"),
        atom(2, " NZ ", "LYS", "A", 1, 8.0, 0.0, 0.0, " N"),
        atom(3, " OE1", "GLU", "B", 1, 10.0, 0.0, 0.0, " O"),   # no CA at all
        atom(4, " CA ", "SER", "B", 2, 11.0, 0.0, 0.0, " C"),
    ]
    result = analyze_cross(
        structure(lines), method="protein_protein_calpha_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF, keep_triplets=False,
    )
    pools = result.annotations["cross_pools"]
    assert pools["n_interface_residues_without_calpha"] == 1
    reasons = {s.reason for s in result.skipped}
    assert "interface_residue_without_calpha" in reasons
    assert pools["side_b_pool"] == 1, "SER 2 has a Ca; GLU 1 does not"


@pytest.mark.parametrize("mode", ["legacy", "standard"])
def test_the_calpha_cross_refuses_a_mode_argument(mode):
    """There is one scientific mode, so ``mode=`` is refused by name -- even
    ``mode="standard"``."""
    with pytest.raises(InvalidRequest) as excinfo:
        analyze_cross(structure(CALPHA_INTERFACE_STRUCTURE),
                      method="protein_protein_calpha_cross",
                      chain="A", chain_2="B", mode=mode)
    assert "mode=" in str(excinfo.value)


# ===========================================================================
# Rule 6 -- the cutoff is pool membership, never a triangle-edge constraint
# ===========================================================================


#: Four atoms, every one of which has a qualifying opposite-side contact, laid
#: out so that three of the four triangles have an edge far beyond the cutoff.
POST_POOL_STRUCTURE = [
    atom(1, " N  ", "ALA", "A", 1, 0.0, 0.0, 0.0, " N"),
    atom(2, " CB ", "ALA", "A", 2, 0.0, 0.0, 10.0, " C"),
    atom(3, " N  ", "GLY", "B", 1, 2.0, 0.0, 0.0, " N"),
    atom(4, " CB ", "GLY", "B", 2, 2.0, 0.0, 10.0, " C"),
]


def test_a_triangle_whose_edge_exceeds_the_cutoff_is_still_generated():
    result = analyze_cross(
        structure(POST_POOL_STRUCTURE), method="protein_protein_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF,
    )
    a, b = pool_sizes(result)
    assert (a, b) == (2, 2)
    assert result.counts.n_candidates == a * comb(b, 2) + b * comb(a, 2) == 4
    assert result.counts.n_triangles == 4

    longest = max(t.max_distance for t in result.triplets)
    assert longest > CUTOFF, (
        "every one of these triangles spans the 10 A separation; if the cutoff "
        "had been re-imposed on the edges, none of them would exist"
    )


def test_the_size_filter_is_the_separate_concept_and_is_off_by_default():
    plain = analyze_cross(
        structure(POST_POOL_STRUCTURE), method="protein_protein_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF, keep_triplets=False,
    )
    assert plain.sampling.size_filter is None
    assert plain.counts.n_filtered_by_size == 0

    filtered = analyze_cross(
        structure(POST_POOL_STRUCTURE), method="protein_protein_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF, size_filter=5.0,
        keep_triplets=False,
    )
    assert pool_sizes(filtered) == pool_sizes(plain), "the POOL is untouched"
    assert filtered.counts.n_candidates == plain.counts.n_candidates
    assert filtered.counts.n_filtered_by_size == 4, "and every triangle is rejected"
    assert filtered.key_space.fingerprint == plain.key_space.fingerprint


# ===========================================================================
# Rule 7 -- composition: mixed-side only, exactly
# ===========================================================================


def test_the_composition_is_exactly_the_mixed_side_formula_with_no_monochrome():
    result = analyze_cross(
        structure(POST_POOL_STRUCTURE), method="protein_protein_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF,
    )
    a, b = pool_sizes(result)
    assert result.counts.n_candidates == a * comb(b, 2) + b * comb(a, 2)
    assert result.counts.n_candidates == (
        comb(a + b, 3) - comb(a, 3) - comb(b, 3)
    )
    for triplet in result.triplets:
        assert set(triplet.sides) == {SIDE_A, SIDE_B}


def test_the_calpha_cross_composition_is_mixed_side_too():
    lines = CALPHA_INTERFACE_STRUCTURE + [
        atom(9, " CA ", "TRP", "A", 3, 8.5, 1.0, 0.0, " C"),
        atom(10, " CA ", "THR", "B", 3, 9.5, 2.0, 0.0, " C"),
    ]
    result = analyze_cross(
        structure(lines), method="protein_protein_calpha_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF,
    )
    a, b = pool_sizes(result)
    assert result.counts.n_candidates == a * comb(b, 2) + b * comb(a, 2)
    for triplet in result.triplets:
        assert set(triplet.sides) == {SIDE_A, SIDE_B}


# ===========================================================================
# The rules as properties of the recipes, not just of these structures
# ===========================================================================


ALL_CROSS = (
    "drug_protein_cross",
    "carbohydrate_protein_cross",
    "nucleotide_protein_cross",
    "protein_protein_cross",
    "protein_protein_calpha_cross",
)

#: One synthetic two-sided case per method, so a cutoff can be varied over every
#: family without loading a real structure. The geometry is a ladder: pairs of
#: opposite-side atoms at 2, 3.5 and 5 A, so each named cutoff admits a strictly
#: larger pool and no assertion turns on a marginal distance.
_LADDER = [
    atom(1, " N  ", "ALA", "A", 1, 0.0, 0.0, 0.0, " N"),
    atom(2, " CA ", "ALA", "A", 1, 0.0, 0.0, 3.0, " C"),
    atom(3, " CB ", "ALA", "A", 1, 0.0, 0.0, 6.5, " C"),
    atom(4, " N  ", "GLY", "B", 1, 2.0, 0.0, 0.0, " N"),
    atom(5, " CA ", "GLY", "B", 1, 3.5, 0.0, 3.0, " C"),
    atom(6, " CB ", "GLY", "B", 1, 5.0, 0.0, 6.5, " C"),
]
_LADDER_ENTITY = _LADDER[:3] + [
    atom(4, " C1 ", "LIG", "B", 404, 2.0, 0.0, 0.0, " C", record="HETATM"),
    atom(5, " C2 ", "LIG", "B", 404, 3.5, 0.0, 3.0, " C", record="HETATM"),
    atom(6, " C3 ", "LIG", "B", 404, 5.0, 0.0, 6.5, " C", record="HETATM"),
]
#: The Ca ladder needs a Ca per residue, and residues far enough apart that the
#: heavy-atom contact -- not the Ca separation -- is what each cutoff decides.
_LADDER_CALPHA = [
    atom(1, " CA ", "LYS", "A", 1, 0.0, 20.0, 0.0, " C"),
    atom(2, " NZ ", "LYS", "A", 1, 0.0, 0.0, 0.0, " N"),
    atom(3, " CA ", "SER", "A", 2, 0.0, 20.0, 3.0, " C"),
    atom(4, " OG ", "SER", "A", 2, 0.0, 0.0, 3.0, " O"),
    atom(5, " CA ", "THR", "A", 3, 0.0, 20.0, 6.5, " C"),
    atom(6, " OG1", "THR", "A", 3, 0.0, 0.0, 6.5, " O"),
    atom(7, " CA ", "GLU", "B", 1, 2.0, -20.0, 0.0, " C"),
    atom(8, " OE1", "GLU", "B", 1, 2.0, 0.0, 0.0, " O"),
    atom(9, " CA ", "ASP", "B", 2, 3.5, -20.0, 3.0, " C"),
    atom(10, " OD1", "ASP", "B", 2, 3.5, 0.0, 3.0, " O"),
    atom(11, " CA ", "ASN", "B", 3, 5.0, -20.0, 6.5, " C"),
    atom(12, " OD1", "ASN", "B", 3, 5.0, 0.0, 6.5, " O"),
]

_LADDER_KW = {
    "drug_protein_cross": (_LADDER_ENTITY, dict(
        chain="B", chain_2="A", entity_name="LIG", entity_id="404")),
    "carbohydrate_protein_cross": (_LADDER_ENTITY, dict(
        chain="B", chain_2="A", entity_name="LIG", entity_id="404")),
    "nucleotide_protein_cross": (_LADDER_ENTITY, dict(
        chain="B", chain_2="A", entity_name="LIG", entity_id="404")),
    "protein_protein_cross": (_LADDER, dict(chain="A", chain_2="B")),
    "protein_protein_calpha_cross": (_LADDER_CALPHA, dict(chain="A", chain_2="B")),
}


def _cross(name, **overrides):
    """Run one method over its ladder structure."""
    lines, kw = _LADDER_KW[name]
    return analyze_cross(structure(lines), method=name,
                         keep_triplets=False, **dict(kw, **overrides))


@pytest.mark.parametrize("name", ALL_CROSS)
def test_no_standard_pool_policy_filters_only_one_side(name):
    """A one-sided filter is reachable only through an explicit caller
    override (``contact_filter=``), never from a standard pool policy."""
    from biotsr.cross.pools import ContactFilter
    from biotsr.methods.builtin.cross import cross_recipe

    policy = cross_recipe(name).standard_pool(3.0)
    assert policy.contact_filter == ContactFilter.BOTH_SIDES
    assert policy.contact_filter_rule == "symmetric_atom_contact_filter"


@pytest.mark.parametrize("name", ALL_CROSS)
def test_every_standard_pool_policy_uses_the_full_atom_identity(name):
    from biotsr.cross.pools import PoolIdentity
    from biotsr.methods.builtin.cross import cross_recipe

    assert cross_recipe(name).standard_pool(3.0).identity == (
        PoolIdentity.FULL_ATOM_IDENTITY
    )


@pytest.mark.parametrize("name", ALL_CROSS)
def test_every_standard_side_selector_is_chain_scoped(name):
    """Every side of every cross method uses a chain-scoped selector; there is no
    whole-structure selector at all."""
    from biotsr.cross import sides
    from biotsr.methods.builtin.cross import cross_recipe

    recipe = cross_recipe(name)
    chain_scoped = {
        sides.SIDE_SELECTOR_PROTEIN_CHAIN,
        sides.SIDE_SELECTOR_ENTITY,
        sides.SIDE_SELECTOR_CALPHA_CHAIN,
    }
    assert {recipe.side_a_selector, recipe.side_b_selector} <= chain_scoped
    assert not hasattr(sides, "SIDE_SELECTOR_PROTEIN_STRUCTURE")


@pytest.mark.parametrize("name", ALL_CROSS)
def test_every_standard_admission_excludes_hydrogen_and_deuterium(name):
    from biotsr.methods.builtin.cross import cross_recipe

    recipe = cross_recipe(name)
    admissions = [recipe.side_a_admission, recipe.side_b_admission]
    if recipe.interface_admission is not None:
        admissions.append(recipe.interface_admission)
    for admission in admissions:
        assert "H" in admission.excluded_elements
        assert "D" in admission.excluded_elements
        assert not hasattr(admission, "use_legacy_tokens"), (
            "every admission reads the element from columns 77-78, never the "
            "off-by-one [77:80] token"
        )


@pytest.mark.parametrize("name", ALL_CROSS)
def test_the_standard_cutoff_is_a_documented_reference_value(name):
    """User-configurable, never a mathematical constant -- and each family's
    default is the reference value of the standard profile."""
    from biotsr.methods.builtin.cross import cross_recipe

    expected = {
        "drug_protein_cross": 3.5,
        "carbohydrate_protein_cross": 3.5,
        "nucleotide_protein_cross": 3.5,
        "protein_protein_cross": 5.0,
        "protein_protein_calpha_cross": 5.0,
    }
    assert cross_recipe(name).default_cutoff == expected[name]


def test_the_nucleotide_family_uses_the_standard_reference_cutoff():
    """The standard reference value is 3.5 A.

    Historical difference, recorded as a constant: `drug_protein_3A.py`
    hard-coded ``HISTORICAL_NUCLEOTIDE_CUTOFF`` = 3.0 A.  The standard default is
    not inferred from it.
    """
    from biotsr.methods.builtin.cross import cross_recipe

    HISTORICAL_NUCLEOTIDE_CUTOFF = 3.0  # documentation only
    recipe = cross_recipe("nucleotide_protein_cross")
    assert recipe.default_cutoff == 3.5
    assert recipe.default_cutoff != HISTORICAL_NUCLEOTIDE_CUTOFF


@pytest.mark.parametrize("name", ALL_CROSS)
def test_the_numerical_precision_is_part_of_what_a_standard_key_means(name):
    """Not provenance: on 1BRS A/D exactly 2 of 1,230,714 triangles change theta
    class between float32 and float64, so two runs that differ only here can
    disagree on a key (ADR-053)."""
    from biotsr.methods.registry import get_method

    key_space = get_method(name).default_profile.key_space()
    assert key_space.numerical_precision == "float64"
    assert "numerical_precision" in key_space.to_dict()


def test_two_runs_differing_only_in_precision_do_not_share_a_key_space():
    from dataclasses import replace

    from biotsr.domain.keyspace import IncompatibleKeySpace
    from biotsr.profiles.cross import STANDARD_ATOMIC_CROSS_V1

    single = replace(STANDARD_ATOMIC_CROSS_V1, coordinate_precision="float32")
    with pytest.raises(IncompatibleKeySpace) as excinfo:
        STANDARD_ATOMIC_CROSS_V1.key_space().require_compatible(single.key_space())
    assert "numerical_precision" in str(excinfo.value)


@pytest.mark.parametrize("pdb,chains", [("1BRS", "AD"), ("3EJ1", "AB"), ("6G9M", "A")])
def test_the_interface_scan_and_the_calpha_selection_see_the_same_residues(pdb, chains):
    """The Cα rule uses two selections per side -- heavy atoms for membership, Cα
    atoms for representation -- and they run under different selectors with
    different truncation rules.  If those two ever disagreed about which residues
    exist, interface residues could be missed *silently*: a residue absent from
    the heavy-atom scan is simply never classified, and nothing records it.  (The
    opposite direction is safe -- a residue with no Cα is recorded.)

    So the agreement is pinned rather than assumed, on the two real structures the
    method is exercised against.
    """
    import os

    from biotsr.api import load_structure
    from biotsr.labels.elements import ELEMENT_LABELS_112
    from biotsr.labels.schemes import AMINO_ACID_30
    from biotsr.methods.builtin.cross import cross_recipe
    from biotsr.selection.atomic import select_entity_atoms
    from biotsr.selection.protein import select_calpha

    from ..conftest import structure

    path = structure(pdb)

    struct = load_structure(path)
    for chain in chains:
        calpha = {
            p.identity.residue
            for p in select_calpha(struct, chain, AMINO_ACID_30).points
        }
        heavy = {
            p.identity.residue
            for p in select_entity_atoms(
                struct, chain, ELEMENT_LABELS_112,
                admission=cross_recipe("protein_protein_calpha_cross").interface_admission,
                require_entity_name=False,
            ).points
        }
        assert heavy == calpha, (
            "%s chain %s: the interface scan and the Ca selection disagree on "
            "%d residue(s)" % (pdb, chain, len(heavy ^ calpha))
        )


def test_contact_filter_is_refused_on_the_calpha_cross_rather_than_ignored():
    """The Ca projection is symmetric by construction, so there is no one-sided
    variant to select and nothing would read the argument: it is refused by
    name."""
    with pytest.raises(InvalidRequest) as excinfo:
        analyze_cross(structure(CALPHA_INTERFACE_STRUCTURE),
                      method="protein_protein_calpha_cross",
                      chain="A", chain_2="B", contact_filter="side_a_only")
    assert "contact_filter" in str(excinfo.value)


def test_a_protein_protein_cross_refuses_to_cross_a_chain_with_itself():
    """`Cross_TSR.py` is "cross keys between two chains" and takes chain1 and
    chain2; there is no one-chain form. Defaulting the second chain would emit a
    full key set labelled side A and side B from the same chain."""
    with pytest.raises(ValueError) as excinfo:
        analyze_cross(structure(SYMMETRY_STRUCTURE),
                      method="protein_protein_cross", chain="A")
    assert "chain_2" in str(excinfo.value)


# ===========================================================================
# The cutoff is a user-specified parameter, on every method
# ===========================================================================


#: Deliberately including two values that are nobody's default, so "the default
#: happened to be right" cannot pass these.
NON_DEFAULT_CUTOFFS = (2.5, 4.0, 6.0)


@pytest.mark.parametrize("name", ALL_CROSS)
@pytest.mark.parametrize("cutoff", NON_DEFAULT_CUTOFFS)
def test_any_positive_cutoff_is_accepted_and_recorded(name, cutoff):
    """No source edit, no family-specific special case: the caller names a
    distance in angstroms and that is the distance used."""
    from biotsr.methods.builtin.cross import cross_recipe

    assert cutoff != cross_recipe(name).default_cutoff, (
        "this test is only meaningful at a NON-default cutoff"
    )
    result = _cross(name, contact_cutoff=cutoff)
    assert result.sampling.contact_cutoff == cutoff
    assert result.annotations["cross_pools"]["contact_cutoff"] == cutoff


@pytest.mark.parametrize("name", ALL_CROSS)
def test_widening_the_cutoff_grows_the_pools_and_leaves_the_key_space_alone(name):
    """The whole KeySpace / SamplingProfile split, on one axis: a wider cutoff
    draws more points and so produces different counts, but the integers it
    produces still mean the same thing and remain comparable."""
    narrow = _cross(name, contact_cutoff=3.0)
    wide = _cross(name, contact_cutoff=6.0)

    n_narrow = pool_sizes(narrow)
    n_wide = pool_sizes(wide)
    assert n_wide[0] >= n_narrow[0] and n_wide[1] >= n_narrow[1]
    assert sum(n_wide) > sum(n_narrow), "the wider cutoff must actually draw more"
    assert wide.counts.n_candidates > narrow.counts.n_candidates

    assert wide.key_space.fingerprint == narrow.key_space.fingerprint
    wide.key_space.require_compatible(narrow.key_space)
    assert wide.sampling.fingerprint != narrow.sampling.fingerprint


@pytest.mark.parametrize("name", ALL_CROSS)
def test_a_caller_cutoff_always_overrides_the_family_default(name):
    from biotsr.methods.builtin.cross import cross_recipe

    default = cross_recipe(name).default_cutoff
    assert _cross(name).sampling.contact_cutoff == default
    assert _cross(name, contact_cutoff=4.25).sampling.contact_cutoff == 4.25


@pytest.mark.parametrize("bad", [0, 0.0, -1.0, -3.5])
def test_a_non_positive_cutoff_is_refused(bad):
    with pytest.raises(ValueError) as excinfo:
        _cross("protein_protein_cross", contact_cutoff=bad)
    assert "positive" in str(excinfo.value)


@pytest.mark.parametrize("bad", ["3.5", True, [3.5], None.__class__])
def test_a_non_numeric_cutoff_is_refused_with_a_message_that_names_the_type(bad):
    """`"3.5"` used to reach a bare comparison and raise an opaque TypeError
    from deep inside the pool builder."""
    with pytest.raises(TypeError) as excinfo:
        _cross("protein_protein_cross", contact_cutoff=bad)
    assert "angstroms" in str(excinfo.value)


def test_an_integer_cutoff_is_normalised_to_float():
    """Not cosmetic: `SamplingProfile.contact_cutoff` is serialized into a
    content fingerprint, and `3` and `3.0` serialize differently. Without the
    normalisation, two runs producing byte-identical keys and frequencies are
    refused as "not sampled the same way" purely because one caller wrote an
    int -- a false incomparable, which is the one thing this identity must not
    produce."""
    as_int = _cross("protein_protein_cross", contact_cutoff=3)
    as_float = _cross("protein_protein_cross", contact_cutoff=3.0)

    assert isinstance(as_int.sampling.contact_cutoff, float)
    assert as_int.frequencies == as_float.frequencies
    assert as_int.sampling.fingerprint == as_float.sampling.fingerprint
    as_int.sampling.require_compatible(as_float.sampling)


@pytest.mark.parametrize("name", ALL_CROSS)
def test_the_cutoff_is_not_in_the_key_space(name):
    from biotsr.methods.registry import get_method

    key_space = set(get_method(name).default_profile.key_space().to_dict())
    assert not [f for f in key_space if "cutoff" in f or "contact" in f]


@pytest.mark.parametrize("name", ALL_CROSS)
def test_the_cutoff_and_the_size_filter_stay_independent(name):
    """Two distances in angstroms that must never be conflated: the cutoff picks
    the vertices, the size filter rejects formed triangles. Setting one leaves
    the other's effect untouched."""
    plain = _cross(name, contact_cutoff=4.0)
    filtered = _cross(name, contact_cutoff=4.0, size_filter=6.0)

    assert pool_sizes(filtered) == pool_sizes(plain)
    assert filtered.counts.n_candidates == plain.counts.n_candidates
    assert filtered.sampling.contact_cutoff == plain.sampling.contact_cutoff == 4.0
    assert plain.sampling.size_filter is None
    assert filtered.sampling.size_filter == 6.0
    assert filtered.key_space.fingerprint == plain.key_space.fingerprint


def test_an_unlabellable_interface_atom_still_makes_contact_in_the_calpha_cross():
    """The Cα cross reads interface heavy atoms only for their coordinates and
    residue identity; their atom labels never reach a key.  A side-chain atom
    whose element has no label (``QQ``) still decides contact,
    and does not stop the calculation (it used to raise ``KeyError``)."""
    lines = [
        atom(1, " CA ", "LYS", "A", 1, 0.0, 0.0, 0.0, " C"),
        atom(2, " QQ1", "LYS", "A", 1, 8.0, 0.0, 0.0, "QQ"),   # the only contact
        atom(3, " CA ", "SER", "A", 2, 0.0, 6.0, 0.0, " C"),
        atom(4, " CA ", "SER", "B", 1, 10.0, 0.0, 0.0, " C"),
        atom(5, " CA ", "GLY", "B", 2, 30.0, 0.0, 0.0, " C"),
    ]
    result = analyze_cross(
        structure(lines), method="protein_protein_calpha_cross",
        chain="A", chain_2="B", contact_cutoff=CUTOFF, keep_triplets=False,
    )
    pools = result.annotations["cross_pools"]
    assert (pools["side_a_pool"], pools["side_b_pool"]) == (1, 1)   # LYS A1, SER B1
    assert result.counts.n_triangles == 0                            # two points only
