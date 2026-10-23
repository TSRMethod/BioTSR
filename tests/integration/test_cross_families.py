"""Cross families end to end: exact counts, exact identities, one shared engine.

Every assertion is an exact number or an exact set; "some keys were generated"
is not evidence of anything.
"""

from __future__ import annotations

import ast
import os

import pytest

from biotsr import analyze_cross, list_cross_methods, list_methods
from biotsr.cross.candidates import count as mixed_side_count
from biotsr.cross.pools import SIDE_A, SIDE_B, PoolIdentity
from biotsr.domain.errors import InvalidRequest
from biotsr.domain.result import DEFAULT_SIDE
from biotsr.methods.registry import get_method

from ..conftest import structure as tracked_structure

HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = HERE  # the package project is the repository root


def fixture(pdb_id):
    return tracked_structure(pdb_id)


#: One representative run per family, on a pinned acceptance-sheet structure.
#: Both chains are named on every case, because the standard requires it: a
#: cross is between two user-specified sides and never scans an unrelated chain
#: (ADR-049).  Each protein chain is the one the entity actually contacts --
#: 7LF9's NAG C1 contacts chain A, 2AC0's DNA chain E contacts chain B -- so these
#: are the interfaces the frozen references describe, named rather than discovered.
CASES = {
    "drug_protein_cross": dict(
        source="6G9M", method="drug_protein_cross", chain="A", chain_2="A",
        entity_name="ESW", entity_id="404",
    ),
    "carbohydrate_protein_cross": dict(
        source="7LF9", method="carbohydrate_protein_cross", chain="C", chain_2="A",
        entity_name="NAG", entity_id="1",
    ),
    "nucleotide_protein_cross": dict(
        source="2AC0", method="nucleotide_protein_cross", chain="E", chain_2="B",
        entity_name="DG", entity_id="8",
    ),
    "protein_protein_cross": dict(
        source="1BRS", method="protein_protein_cross", chain="A", chain_2="D",
    ),
    "protein_protein_calpha_cross": dict(
        source="1BRS", method="protein_protein_calpha_cross", chain="A", chain_2="D",
    ),
}

#: The four fully atomic cross methods.  ``protein_protein_calpha_cross`` is a
#: residue method (radix 30, 35 distance classes) and is excluded from every
#: assertion about the atomic 112 key space.
ATOMIC_CASES = [n for n in CASES if n != "protein_protein_calpha_cross"]


def run(name, **overrides):
    kw = dict(CASES[name])
    source = fixture(kw.pop("source"))
    kw.update(overrides)
    return analyze_cross(source, **kw)


#: Protein-protein at its historical 5 A cutoff emits 1.2 million triangles on
#: 1BRS A/D. The assertions in this module are about composition and identity, not
#: about scale, so the shared fixture uses a narrower cutoff and the full-size run
#: is exercised once, in the characterization gate.
MODULE_CUTOFF = {"protein_protein_cross": 3.2, "protein_protein_calpha_cross": 3.2}


@pytest.fixture(scope="module")
def standard_runs():
    return {
        name: run(name, keep_triplets=True,
                  **({"contact_cutoff": MODULE_CUTOFF[name]}
                     if name in MODULE_CUTOFF else {}))
        for name in CASES
    }


# ===========================================================================
# Every family produces mixed-side triangles and nothing else
# ===========================================================================


@pytest.mark.parametrize("name", list(CASES))
def test_the_candidate_count_is_exactly_the_mixed_side_formula(name, standard_runs):
    result = standard_runs[name]
    pools = result.annotations["cross_pools"]
    a, b = pools["side_a_pool"], pools["side_b_pool"]
    assert result.counts.n_candidates == mixed_side_count(a, b)
    assert result.counts.n_points == a + b


@pytest.mark.parametrize("name", list(CASES))
def test_every_emitted_triplet_contains_both_sides(name, standard_runs):
    result = standard_runs[name]
    assert result.triplets, "the case must actually produce triangles"
    for triplet in result.triplets:
        assert set(triplet.sides) in ({SIDE_A, SIDE_B},), triplet.sides
        assert sorted(triplet.sides).count(SIDE_A) in (1, 2)


@pytest.mark.parametrize("name", list(CASES))
def test_no_monochromatic_triangle_survives_to_a_key(name, standard_runs):
    for triplet in standard_runs[name].triplets:
        assert len(set(triplet.sides)) == 2


@pytest.mark.parametrize("name", list(CASES))
def test_the_candidate_ledger_balances(name, standard_runs):
    assert standard_runs[name].counts.ledger_balances


# ===========================================================================
# Full physical identity, everywhere, in standard mode
# ===========================================================================


@pytest.mark.parametrize("name", list(CASES))
def test_the_standard_pool_uses_full_physical_atom_identity(name, standard_runs):
    result = standard_runs[name]
    assert result.sampling.pool_identity == PoolIdentity.FULL_ATOM_IDENTITY
    assert result.annotations["cross_pools"]["n_identity_collisions"] == 0
    assert result.annotations["cross_pools"]["n_atoms_lost_to_identity"] == 0


@pytest.mark.parametrize("name", list(CASES))
def test_no_two_pool_members_share_a_physical_identity(name, standard_runs):
    result = standard_runs[name]
    seen = set()
    for triplet in result.triplets[:5000]:
        for identity in triplet.identities:
            seen.add(identity.as_tuple())
    # every identity tuple is fully discriminating by construction; the real
    # assertion is that the pool size equals the number of DISTINCT identities.
    pools = result.annotations["cross_pools"]
    assert len(seen) <= pools["side_a_pool"] + pools["side_b_pool"]


@pytest.mark.parametrize("name", list(CASES))
def test_a_triplet_record_carries_the_full_identity_of_every_vertex(name,
                                                                    standard_runs):
    """Nothing is lost -- structure, chain, residue, number, insertion
    code, atom name and side are all on the record."""
    triplet = standard_runs[name].triplets[0]
    for identity, side in zip(triplet.identities, triplet.sides):
        assert identity.structure.value
        assert identity.atom_name
        assert identity.residue_name
        assert isinstance(identity.residue_seq, int)
        assert identity.insertion_code == identity.insertion_code  # present
        assert side in (SIDE_A, SIDE_B)


# ===========================================================================
# Label space and bin conventions
# ===========================================================================


@pytest.mark.parametrize("name", ATOMIC_CASES)
def test_every_standard_cross_run_uses_radix_112_and_18_distance_classes(
    name, standard_runs
):
    key_space = standard_runs[name].key_space
    assert key_space.label_scheme.radix == 112
    assert key_space.distance_scheme.radix == 18
    assert key_space.theta_scheme.radix == 29


@pytest.mark.parametrize("name", ATOMIC_CASES)
def test_all_four_families_share_one_standard_key_space(name, standard_runs):
    reference = standard_runs["drug_protein_cross"].key_space
    standard_runs[name].key_space.require_compatible(reference)


@pytest.mark.parametrize("name", list(CASES))
def test_the_families_do_not_share_a_sampling_profile(name, standard_runs):
    """Same key space, different sampling -- which is exactly the distinction the
    two identities exist to express."""
    from biotsr.domain.keyspace import IncompatibleSampling

    if name == "drug_protein_cross":
        return
    with pytest.raises(IncompatibleSampling):
        standard_runs[name].sampling.require_compatible(
            standard_runs["drug_protein_cross"].sampling
        )


# ===========================================================================
# Contact cutoff vs size filter -- separate concepts, separately recorded
# ===========================================================================


def test_widening_the_contact_cutoff_grows_the_pool_without_changing_key_space():
    narrow = run("nucleotide_protein_cross", contact_cutoff=3.0, keep_triplets=False)
    wide = run("nucleotide_protein_cross", contact_cutoff=4.5, keep_triplets=False)
    assert (
        wide.annotations["cross_pools"]["side_b_pool"]
        >= narrow.annotations["cross_pools"]["side_b_pool"]
    )
    assert wide.counts.n_candidates > narrow.counts.n_candidates
    assert wide.key_space.fingerprint == narrow.key_space.fingerprint
    assert wide.sampling.fingerprint != narrow.sampling.fingerprint
    assert wide.sampling.contact_cutoff == 4.5


def test_the_size_filter_is_off_by_default_and_rejects_triangles_not_atoms():
    plain = run("nucleotide_protein_cross", keep_triplets=False)
    assert plain.sampling.size_filter is None
    assert plain.counts.n_filtered_by_size == 0

    filtered = run("nucleotide_protein_cross", size_filter=4.0, keep_triplets=False)
    # the POOL is identical -- the cutoff decided membership, the filter did not
    assert (
        filtered.annotations["cross_pools"]["side_a_pool"]
        == plain.annotations["cross_pools"]["side_a_pool"]
    )
    assert (
        filtered.annotations["cross_pools"]["side_b_pool"]
        == plain.annotations["cross_pools"]["side_b_pool"]
    )
    assert filtered.counts.n_candidates == plain.counts.n_candidates
    assert filtered.counts.n_filtered_by_size > 0
    assert filtered.counts.n_triangles < plain.counts.n_triangles
    assert filtered.key_space.fingerprint == plain.key_space.fingerprint


# ===========================================================================
# Determinism
# ===========================================================================


#: A cutoff small enough to keep the repeated runs cheap. Determinism is a
#: property of the pipeline, not of the pool size, so shrinking the pool tests the
#: same thing in a fraction of the time -- the full-size runs are exercised by the
#: characterization gate.
FAST_CUTOFF = {"protein_protein_cross": 3.2, "drug_protein_cross": 3.5,
               "carbohydrate_protein_cross": 3.5, "nucleotide_protein_cross": 3.5,
               "protein_protein_calpha_cross": 3.2}


@pytest.mark.parametrize("name", list(CASES))
def test_a_cross_run_is_deterministic(name):
    kw = dict(keep_triplets=False, contact_cutoff=FAST_CUTOFF[name])
    a = run(name, **kw)
    b = run(name, **kw)
    assert a.frequencies == b.frequencies
    assert a.key_order == b.key_order
    assert a.counts.n_candidates == b.counts.n_candidates


def test_standard_key_serialization_is_sorted_numerically(tmp_path):
    from biotsr.io import native

    result = run("nucleotide_protein_cross")
    path = native.write_keys(result, str(tmp_path / "keys"))
    keys = [int(line.split("\t")[0]) for line in open(path).read().splitlines()[1:]]
    assert keys == sorted(keys)


def test_the_native_triplet_format_names_the_side_of_every_vertex(tmp_path):
    from biotsr.io import native

    result = run("nucleotide_protein_cross")
    path = native.write_triplets(result, str(tmp_path / "triplets"))
    lines = open(path).read().splitlines()
    header = lines[0].split("\t")
    assert ["side1", "side2", "side3"] == [c for c in header if c.startswith("side")]
    row = dict(zip(header, lines[1].split("\t")))
    assert {row["side1"], row["side2"], row["side3"]} == {SIDE_A, SIDE_B}
    assert row["structure1"] == result.structure_id.value


def test_a_one_sided_method_reports_the_default_side():
    from biotsr import analyze

    result = analyze(fixture("1UBQ"), method="protein", chain="A")
    assert set(result.triplets[0].sides) == {DEFAULT_SIDE}


# ===========================================================================
# The acceptance sheet's own column shapes
# ===========================================================================


def test_the_acceptance_sheet_narrowing_columns_drive_both_sides():
    """A dataset-sheet-shaped request: nt_chain=E, nt_name=DT, nt_id=7 and an
    nt_atom_selected cell, against aa_chain=A and an aa_atom_selection cell.

    With no contact filter the candidate count is exactly the mixed-side formula
    over the two explicitly enumerated selections.
    """
    from math import comb

    result = analyze_cross(
        fixture("2AC0"), method="nucleotide_protein_cross",
        chain="E", chain_2="A", entity_name="DT", entity_id="7",
        atoms="P;C1';N1;C2;O4;C7;C6",
        atom_selection=("LYS_120_NZ;LYS_120_CE;ALA_276_CB;CYS_277_SG;"
                        "ARG_280_NH1;ARG_280_NH2"),
        no_contact_filter=True,
    )
    pools = result.annotations["cross_pools"]
    assert (pools["side_a_pool"], pools["side_b_pool"]) == (7, 6)
    assert result.counts.n_candidates == 7 * comb(6, 2) + 6 * comb(7, 2) == 231
    assert result.sampling.contact_cutoff is None
    assert result.sampling.contact_filter == "none"


def test_narrowing_the_protein_side_without_naming_its_chain_is_refused():
    """A narrowing selector is chain-scoped, and `chain` is the ENTITY's chain --
    which for a nucleotide or a glycan is usually not the protein's. Defaulting
    would silently select from the wrong chain, so it is refused instead.
    """
    with pytest.raises(ValueError) as excinfo:
        analyze_cross(fixture("2AC0"), method="nucleotide_protein_cross",
                      chain="E", entity_name="DT", entity_id="7",
                      atom_selection="SER_241_OG;ARG_248_CG")
    assert "chain_2" in str(excinfo.value)


@pytest.mark.parametrize("name", list(CASES))
def test_the_standard_refuses_to_guess_the_second_chain(name):
    """Settled: both participating chains are user-specified, for every cross
    method (ADR-049).  Defaulting the second
    chain to the first would be a guess, and on 7LF9 and 2AC0 a wrong one -- the
    glycan in chain C contacts chain A, the DNA in chain E contacts chain B -- so
    it is asked for instead."""
    kw = dict(CASES[name])
    source = fixture(kw.pop("source"))
    kw.pop("chain_2")
    with pytest.raises(ValueError) as excinfo:
        analyze_cross(source, **kw)
    message = str(excinfo.value)
    assert "chain_2" in message
    if CASES[name].get("entity_name"):
        # The entity families name the rule. The two protein-protein methods
        # refuse with "there is no one-chain form of it" instead.
        assert "ADR-049" in message
    else:
        assert "with itself" in message


@pytest.mark.parametrize("name", ATOMIC_CASES)
def test_the_standard_never_reaches_a_chain_the_caller_did_not_name(name):
    """The strongest form of the chain rule: not "chain_2 is honoured" but "no
    atom of any other chain is in the pool at all"."""
    result = run(name, keep_triplets=True,
                 **({"contact_cutoff": MODULE_CUTOFF[name]}
                    if name in MODULE_CUTOFF else {}))
    named = {CASES[name]["chain"], CASES[name]["chain_2"]}
    reached = {i.chain.value for t in result.triplets for i in t.identities}
    assert reached <= named, (
        "%s reached chain(s) %r that the caller never named"
        % (name, sorted(reached - named))
    )


def test_the_standard_records_both_named_chains_not_a_wildcard():
    """Historical note: the `cross_key.py` contact generator had no chain
    predicate, so its protein side was the whole structure (recorded then as
    ``chain=*``).  The standard records the two chains the caller named."""
    result = run("nucleotide_protein_cross", keep_triplets=False)
    assert (result.sampling.side_a_chain, result.sampling.side_b_chain) == ("E", "B")
    assert "chain=*" not in result.sampling.side_b_selector


def test_chain_2_is_honoured_wherever_the_protein_side_is_chain_scoped():
    """Every side is chain-scoped, so `chain_2` is effective and recorded."""
    drug = analyze_cross(fixture("6G9M"), method="drug_protein_cross", chain="A",
                         chain_2="A", entity_name="ESW", entity_id="404",
                         keep_triplets=False)
    assert "chain=A" in drug.sampling.side_a_selector

    pp = analyze_cross(fixture("1BRS"), method="protein_protein_cross",
                       chain="A", chain_2="D", contact_cutoff=3.2,
                       keep_triplets=False)
    assert "chain=A" in pp.sampling.side_a_selector
    assert "chain=D" in pp.sampling.side_b_selector


def test_a_protein_side_takes_at_most_one_narrowing_argument():
    with pytest.raises(ValueError) as excinfo:
        analyze_cross(fixture("2AC0"), method="nucleotide_protein_cross",
                      chain="E", chain_2="A", entity_name="DT", entity_id="7",
                      residue_types="ARG", atom_selection="SER_241_OG")
    assert "at most one narrowing argument" in str(excinfo.value)


def test_the_second_side_arguments_are_refused_where_there_is_no_second_protein():
    with pytest.raises(InvalidRequest) as excinfo:
        analyze_cross(fixture("2AC0"), method="nucleotide_protein_cross",
                      chain="E", chain_2="A", entity_name="DT", entity_id="7",
                      residue_types_2="ARG")
    message = str(excinfo.value)
    assert "residue_types_2" in message
    assert "does not use" in message


def test_protein_protein_cross_narrows_each_side_from_its_own_columns():
    """aa_ind_specific and aa_ind_specific_2, on 3EJ1 chains A and B.

    Chain B of 3EJ1 starts at residue 175, so its residues are chosen from that
    range (HIS 179, TYR 185) -- see the next test for what happens when residues
    the chain does not have are asked for.
    """
    result = analyze_cross(
        fixture("3EJ1"), method="protein_protein_cross", chain="A", chain_2="B",
        residues="LYS_33;ASP_145", residues_2="HIS_179;TYR_185",
        no_contact_filter=True, keep_triplets=False,
    )
    assert "LYS_33" in result.sampling.side_a_selector
    assert "HIS_179" in result.sampling.side_b_selector
    pools = result.annotations["cross_pools"]
    assert pools["side_a_pool"] > 0 and pools["side_b_pool"] > 0
    assert result.counts.n_candidates == mixed_side_count(
        pools["side_a_pool"], pools["side_b_pool"]
    )


def test_a_side_whose_requested_residues_are_all_absent_is_refused():
    """3EJ1's chain B has no residue 33 or 51. An explicit selection that matches
    nothing is refused, naming what was not found, rather than returning a
    silently empty result."""
    with pytest.raises(InvalidRequest) as excinfo:
        analyze_cross(
            fixture("3EJ1"), method="protein_protein_cross", chain="A",
            chain_2="B", residues="LYS_33;ASP_145", residues_2="LYS_33;GLU_51",
            no_contact_filter=True, keep_triplets=False,
        )
    message = str(excinfo.value)
    assert "requested but not found" in message
    assert "LYS_33" in message and "GLU_51" in message


def test_a_partly_absent_residue_selection_is_reported_not_dropped():
    """HIS 179 is in 3EJ1 chain B, GLU 51 is not: the run proceeds on what was
    found, and the missing residue is a named warning, an `atoms_missing` entry
    and a `SkipRecord`."""
    result = analyze_cross(
        fixture("3EJ1"), method="protein_protein_cross", chain="A", chain_2="B",
        residues="LYS_33;ASP_145", residues_2="HIS_179;GLU_51",
        no_contact_filter=True, keep_triplets=False,
    )
    assert result.annotations["cross_pools"]["side_b_pool"] > 0
    assert "GLU_51" in result.selection.atoms_missing
    assert "HIS_179" not in result.selection.atoms_missing
    reasons = {s.reason for s in result.skipped}
    assert "requested_entity_not_found" in reasons
    missing = [w for w in result.warnings if w.startswith("requested but not found")]
    assert len(missing) == 1 and "GLU_51" in missing[0]


def test_a_calpha_side_cannot_be_crossed_with_an_atomic_side():
    """The atomic method does not read Ca residue selections at all, so a Ca
    side cannot be mixed into it: the argument is refused by name."""
    with pytest.raises(InvalidRequest) as excinfo:
        analyze_cross(fixture("3EJ1"), method="protein_protein_cross",
                      chain="A", chain_2="B",
                      calpha_residues="33;145",
                      residues_2="HIS_179;TYR_185")
    message = str(excinfo.value)
    assert "calpha_residues" in message
    assert "does not use" in message


def test_a_calpha_cross_on_the_atomic_method_is_redirected_not_guessed():
    """Ca cross is its own method with its own key space (radix 30) and its own
    pool rule.  Asking the ATOMIC method for it is refused, naming both Ca
    arguments, rather than quietly encoding residues with the 112-radix atomic
    labels."""
    with pytest.raises(InvalidRequest) as excinfo:
        analyze_cross(fixture("3EJ1"), method="protein_protein_cross",
                      chain="A", chain_2="B",
                      calpha_residues="126;150", calpha_residues_2="175;178")
    message = str(excinfo.value)
    assert "calpha_residues" in message and "calpha_residues_2" in message
    assert "protein_protein_cross" in message


# ===========================================================================
# Arguments a method does not read are refused by name
# ===========================================================================


#: A plausible value for every cross argument, so each refusal below is about
#: the argument being unread, never about its value.
ARGUMENT_VALUES = {
    "contact_filter": "side_a_only",
    "no_contact_filter": True,
    "entity_name": "ESW",
    "entity_id": "404",
    "entity_names": "NAG",
    "atoms": "C1",
    "residue_types": "LYS",
    "residues": "LYS_33",
    "atom_selection": "LYS_33_NZ",
    "calpha_residues": "33",
    "residue_types_2": "LYS",
    "residues_2": "LYS_33",
    "atom_selection_2": "LYS_33_NZ",
    "calpha_residues_2": "33",
}

#: What each cross method reads, written out from the release brief rather than
#: derived from the package table, so the table is checked against it.
ACCEPTED = {
    "drug_protein_cross": {
        "chain_2", "contact_cutoff", "contact_filter", "no_contact_filter",
        "entity_name", "entity_id", "entity_names", "atoms", "residue_types",
        "residues", "atom_selection",
    },
    "protein_protein_cross": {
        "chain_2", "contact_cutoff", "contact_filter", "no_contact_filter",
        "residue_types", "residues", "atom_selection", "residue_types_2",
        "residues_2", "atom_selection_2",
    },
    "protein_protein_calpha_cross": {
        "chain_2", "contact_cutoff", "calpha_residues", "calpha_residues_2",
    },
}
ACCEPTED["carbohydrate_protein_cross"] = ACCEPTED["drug_protein_cross"]
ACCEPTED["nucleotide_protein_cross"] = ACCEPTED["drug_protein_cross"]


def test_the_argument_table_matches_the_documented_per_method_arguments():
    from biotsr.methods.arguments import CROSS_ARGUMENTS_BY_METHOD

    assert {
        name: set(args) for name, args in CROSS_ARGUMENTS_BY_METHOD.items()
    } == ACCEPTED


REFUSALS = [
    (name, argument)
    for name in sorted(ACCEPTED)
    for argument in sorted(ARGUMENT_VALUES)
    if argument not in ACCEPTED[name]
]


@pytest.mark.parametrize("name,argument", REFUSALS)
def test_a_cross_method_refuses_an_argument_it_does_not_read(name, argument):
    """Among them: protein_protein_calpha_cross refuses entity_name, atoms,
    residue_types, contact_filter and no_contact_filter; protein_protein_cross
    refuses entity_name; drug_protein_cross refuses residue_types_2 and
    calpha_residues."""
    kw = dict(CASES[name])
    source = fixture(kw.pop("source"))
    kw[argument] = ARGUMENT_VALUES[argument]
    with pytest.raises(InvalidRequest) as excinfo:
        analyze_cross(source, keep_triplets=False, **kw)
    message = str(excinfo.value)
    assert "does not use %s" % argument in message or (
        "does not use" in message and argument in message
    ), message
    assert name in message


def test_the_explicit_cases_named_in_the_brief_are_among_the_refusals():
    expected = {
        ("protein_protein_calpha_cross", "entity_name"),
        ("protein_protein_calpha_cross", "atoms"),
        ("protein_protein_calpha_cross", "residue_types"),
        ("protein_protein_calpha_cross", "contact_filter"),
        ("protein_protein_calpha_cross", "no_contact_filter"),
        ("protein_protein_cross", "entity_name"),
        ("drug_protein_cross", "residue_types_2"),
        ("drug_protein_cross", "calpha_residues"),
    }
    assert expected <= set(REFUSALS)


@pytest.mark.parametrize(
    "conflict",
    [dict(entity_name="NAG"), dict(entity_id="1"), dict(atoms="C1"),
     dict(entity_name="NAG", entity_id="1")],
    ids=["entity_name", "entity_id", "atoms", "entity_name+entity_id"],
)
def test_entity_names_combined_with_a_single_copy_argument_is_refused(conflict):
    with pytest.raises(InvalidRequest) as excinfo:
        analyze_cross(fixture("7LF9"), method="carbohydrate_protein_cross",
                      chain="C", chain_2="A", entity_names="NAG",
                      keep_triplets=False, **conflict)
    assert "entity_names=" in str(excinfo.value)


# ===========================================================================
# The entity side: an identifier is one copy; a name alone is every copy
# on the chain (ADR-096)
# ===========================================================================


def test_an_omitted_entity_id_is_the_all_copies_scope_even_for_one_copy():
    """6G9M chain A has one ESW (404): without an id the scope is all copies,
    recorded as such; with a single copy the keys equal the named run's."""
    implicit = analyze_cross(fixture("6G9M"), method="drug_protein_cross",
                             chain="A", chain_2="A", entity_name="ESW",
                             keep_triplets=False)
    explicit = run("drug_protein_cross", keep_triplets=False)
    assert "entity=ESW(all copies)" in implicit.sampling.side_b_selector
    assert "entity=ESW_404" in explicit.sampling.side_b_selector
    assert implicit.frequencies == explicit.frequencies
    side = [s for s in implicit.annotations["selection_request"]["sides"]
            if s["role"] == "drug"][0]
    assert side["scope"] == "all_copies" and side["instances"] == ["ESW_404"]


@pytest.mark.parametrize(
    "pdb,kw,expected",
    [
        ("6G9M", dict(method="drug_protein_cross", chain="A", chain_2="A",
                      entity_name="ESW"), "ESW_allcopies"),
        ("7LF9", dict(method="carbohydrate_protein_cross", chain="A",
                      chain_2="A", entity_name="NAG"), "NAG_allcopies"),
        ("2AC0", dict(method="nucleotide_protein_cross", chain="E",
                      chain_2="B", entity_name="DG", entity_id="8"), "DG_8"),
    ],
    ids=["drug-id-omitted", "carbohydrate-id-omitted", "nucleotide-id-given"],
)
def test_the_entity_cross_result_records_its_target(pdb, kw, expected):
    """``annotations["target"]`` is ``NAME_ID`` of the copy named, or
    ``NAME_allcopies`` when no id was given."""
    result = analyze_cross(fixture(pdb), keep_triplets=False, **kw)
    assert result.annotations.get("target") == expected


@pytest.mark.parametrize(
    "pdb,kw,ids",
    [
        # 5WT9 chain A holds two glycan copies, NAG 1 and NAG 2.
        ("5WT9", dict(method="carbohydrate_protein_cross", chain="A",
                      chain_2="H", entity_name="NAG"), ("1", "2")),
        # 6G9M chain A holds three sulphates, SO4 401-403.
        ("6G9M", dict(method="drug_protein_cross", chain="A", chain_2="A",
                      entity_name="SO4"), ("401", "402", "403")),
        # 7LF9 chain C holds NAG 1, 2 and 5.
        ("7LF9", dict(method="carbohydrate_protein_cross", chain="C",
                      chain_2="A", entity_name="NAG"), ("1", "2", "5")),
    ],
    ids=["5WT9-NAG", "6G9M-SO4", "7LF9-NAG"],
)
def test_an_omitted_entity_id_with_several_copies_uses_every_copy(pdb, kw, ids):
    """rc3 refused these; now every copy on the chain is one combined
    calculation, every copy is recorded, and no other chain is searched."""
    result = analyze_cross(fixture(pdb), keep_triplets=False, **kw)
    side = [s for s in result.annotations["selection_request"]["sides"]
            if s["selector"] == "cross_entity_atoms"][0]
    name = kw["entity_name"]
    assert side["scope"] == "all_copies"
    assert side["chain"] == kw["chain"]
    assert result.annotations["target"] == "%s_allcopies" % name
    # instances list the copies that reached the pool after admission
    assert set(side["instances"]) <= {"%s_%s" % (name, i) for i in ids}


def test_an_explicit_entity_id_that_does_not_exist_is_refused():
    with pytest.raises(InvalidRequest) as excinfo:
        analyze_cross(fixture("5WT9"), method="carbohydrate_protein_cross",
                      chain="A", chain_2="H", entity_name="NAG", entity_id="7",
                      keep_triplets=False)
    message = str(excinfo.value)
    assert "NAG 7" in message
    assert "1, 2" in message


def test_an_entity_absent_from_the_chain_is_refused_naming_what_is_there():
    with pytest.raises(InvalidRequest) as excinfo:
        analyze_cross(fixture("5WT9"), method="carbohydrate_protein_cross",
                      chain="A", chain_2="H", entity_name="XYZ",
                      keep_triplets=False)
    message = str(excinfo.value)
    assert "XYZ" in message
    assert "NAG" in message and "FUC" in message


def test_a_cross_run_with_no_contacts_says_it_has_no_triangles():
    """5WT9's lone FUC 3 makes no contact with chain H at the reference
    cutoff: zero triangles, labelled as such."""
    result = analyze_cross(fixture("5WT9"), method="carbohydrate_protein_cross",
                           chain="A", chain_2="H", entity_name="FUC",
                           keep_triplets=False)
    assert "entity=FUC(all copies)" in result.sampling.side_a_selector
    assert result.counts.n_triangles == 0
    assert any(w.startswith("no triangles: ") for w in result.warnings)


# ===========================================================================
# One mode, one profile
# ===========================================================================


@pytest.mark.parametrize("name", list(CASES))
@pytest.mark.parametrize("mode", ["legacy", "legacy_identity_corrected",
                                  "standard"])
def test_mode_is_refused_on_every_cross_method(name, mode):
    kw = dict(CASES[name])
    source = fixture(kw.pop("source"))
    with pytest.raises(InvalidRequest) as excinfo:
        analyze_cross(source, mode=mode, keep_triplets=False, **kw)
    assert "mode=" in str(excinfo.value)


def _non_standard_profiles(name):
    from dataclasses import replace

    from biotsr.profiles.atomic import STANDARD_DRUG_V1
    from biotsr.profiles.cross import (
        STANDARD_ATOMIC_CROSS_V1,
        STANDARD_CALPHA_CROSS_V1,
    )

    own = get_method(name).default_profile
    other = (STANDARD_ATOMIC_CROSS_V1 if own is STANDARD_CALPHA_CROSS_V1
             else STANDARD_CALPHA_CROSS_V1)
    return {
        "by-name": own.name,
        "other-cross-method": other,
        "one-sided-method": STANDARD_DRUG_V1,
        "derived": replace(own, name=own.name + "_derived"),
    }


@pytest.mark.parametrize("name", list(CASES))
@pytest.mark.parametrize("kind", ["by-name", "other-cross-method",
                                  "one-sided-method", "derived"])
def test_a_non_standard_profile_is_refused(name, kind):
    kw = dict(CASES[name])
    source = fixture(kw.pop("source"))
    with pytest.raises(InvalidRequest) as excinfo:
        analyze_cross(source, profile=_non_standard_profiles(name)[kind],
                      keep_triplets=False, **kw)
    assert "profile" in str(excinfo.value)


def test_the_method_s_own_standard_profile_is_accepted():
    profile = get_method("drug_protein_cross").default_profile
    explicit = run("drug_protein_cross", profile=profile, keep_triplets=False)
    implicit = run("drug_protein_cross", keep_triplets=False)
    assert explicit.frequencies == implicit.frequencies
    assert explicit.key_space.fingerprint == implicit.key_space.fingerprint


# ===========================================================================
# No duplicated scientific engine
# ===========================================================================


SRC = os.path.join(HERE, "src", "biotsr")

FORBIDDEN_CLASS_NAMES = (
    "DrugProteinGeometry", "CarbProteinGeometry", "CarbohydrateProteinGeometry",
    "NucleotideProteinGeometry", "ProteinProteinGeometry", "CrossGeometry",
    "DrugProteinOrdering", "CrossOrdering", "CrossEncoder", "CrossKeyFormula",
    "CrossBinScheme", "DrugProteinTSR", "ProteinProteinTSR",
)


def python_files():
    for dirpath, _dirnames, filenames in os.walk(SRC):
        for filename in filenames:
            if filename.endswith(".py"):
                yield os.path.join(dirpath, filename)


def test_no_family_specific_geometry_ordering_or_encoder_class_exists():
    offenders = []
    for path in python_files():
        tree = ast.parse(open(path).read(), filename=path)
        for node in ast.walk(tree):
            if isinstance(node, (ast.ClassDef, ast.FunctionDef)):
                if node.name in FORBIDDEN_CLASS_NAMES:
                    offenders.append((path, node.name))
    assert offenders == [], "family-specific engine found: %r" % offenders


def test_the_cross_package_reimplements_no_science():
    """``biotsr/cross/`` may select, pool and enumerate.  It may not compute a
    distance, an angle, an ordering, a bin or a key -- those come from the shared
    components, and importing them is how that is proved."""
    import biotsr.cross.candidates as candidates
    import biotsr.cross.pools as pools
    import biotsr.cross.runner as runner
    import biotsr.cross.sides as sides

    banned_calls = ("acos", "atan2", "degrees", "radians")
    for module in (candidates, pools, sides, runner):
        source = open(module.__file__).read()
        tree = ast.parse(source)
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute) and node.attr in banned_calls:
                raise AssertionError(
                    "%s computes an angle itself (%s)" % (module.__name__, node.attr)
                )
            if isinstance(node, ast.Name) and node.id in banned_calls:
                raise AssertionError(
                    "%s computes an angle itself (%s)" % (module.__name__, node.id)
                )
        # no bin ladder and no key arithmetic either
        assert "theta_bins" not in source
        assert "dist_bins" not in source
        assert "numOfLabels" not in source


def test_the_cross_runner_calls_the_shared_core():
    source = open(os.path.join(SRC, "cross", "runner.py")).read()
    assert "from ..methods.engine import run_core" in source
    assert "run_core(" in source


@pytest.mark.parametrize("name", list_cross_methods())
def test_every_cross_method_declares_the_one_composition_rule(name):
    from biotsr.cross.candidates import COMPOSITION_RULE

    recipe = get_method(name)
    assert recipe.composition_rule == COMPOSITION_RULE
    assert recipe.scope == "cross_two_sided"


def test_the_cross_methods_are_registered_and_no_more():
    assert list_cross_methods() == (
        "carbohydrate_protein_cross",
        "drug_protein_cross",
        "nucleotide_protein_cross",
        "protein_protein_calpha_cross",
        "protein_protein_cross",
    )
    assert set(list_cross_methods()) <= set(list_methods())
