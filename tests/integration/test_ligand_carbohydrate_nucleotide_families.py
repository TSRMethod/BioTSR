"""Drug, carbohydrate and nucleotide end to end.

Real fixtures, real chemistry, the shared atomic core.  Nothing here constructs a
geometry engine, an encoder or a candidate iterator: every family is a method
name plus a target.
"""

from __future__ import annotations

import math

import pytest

import biotsr
from biotsr.methods.registry import get_method
from biotsr.selection.completeness import (
    CARBOHYDRATE_RESIDUE_NAMES,
    DNA_RESIDUE_NAMES,
    LEGACY_CARBOHYDRATE_ATOM_COUNTS,
    LEGACY_NUCLEOTIDE_ATOM_COUNTS,
    RNA_RESIDUE_NAMES,
)

ATOMIC_FAMILY_METHODS = (
    "drug", "drug_atom_selection",
    "carbohydrate", "carbohydrate_chain", "carbohydrate_atom_selection",
    "nucleotide", "nucleotide_chain", "nucleotide_atom_selection",
)


def _ledger_balances(result):
    c = result.counts
    return (
        c.n_candidates == c.n_valid_geometry + c.n_degenerate
        and c.n_triangles == c.n_valid_geometry - c.n_filtered_by_size
    )


# ===========================================================================
# Drug
# ===========================================================================


def test_drug_regular_uses_every_heavy_atom_of_the_ligand(mapk_path):
    result = biotsr.analyze(
        mapk_path, method="drug", chain="A", entity_name="ESW", entity_id="404"
    )
    assert result.counts.n_points == 37
    assert result.counts.n_candidates == math.comb(37, 3) == 7770
    assert result.counts.n_triangles == 7770
    assert result.counts.n_degenerate == 0
    assert _ledger_balances(result)
    # Hydrogens are excluded and the exclusion is recorded, not silent.
    assert sum(1 for s in result.skipped if s.reason == "hydrogen_or_deuterium") == 30


def test_drug_reads_the_wide_atom_name_slice(mapk_path):
    """6G9M's ESW carries a CL1 -- the narrow [13:16] slice turns it into L1.

    This is not hypothetical: under the narrow slice the historical lexicon had
    no entry for the atom and the ligand silently lost a vertex.
    """
    result = biotsr.analyze(
        mapk_path, method="drug", chain="A", entity_name="ESW", entity_id="404"
    )
    labels = {
        i.atom_name: l
        for t in result.triplets
        for i, l in zip(t.identities, t.labels)
    }
    assert "CL1" in labels
    assert "L1" not in labels
    assert labels["CL1"] == 93          # chlorine, not an unlabelled drop


def test_drug_atom_selection_narrows_to_the_named_atoms(mapk_path):
    """A drug_atom_selected cell: seven named atoms of ESW 404, CL1 among them."""
    result = biotsr.analyze(
        mapk_path, method="drug_atom_selection", chain="A",
        entity_name="ESW", entity_id="404", atoms="C3;C5;O6;N8;C9;C18;CL1",
    )
    assert result.counts.n_points == 7
    assert result.counts.n_triangles == math.comb(7, 3) == 35
    names = {i.atom_name for t in result.triplets for i in t.identities}
    assert names == {"C3", "C5", "O6", "N8", "C9", "C18", "CL1"}
    assert result.selection.atoms_missing == ()


def test_a_requested_drug_atom_that_is_absent_is_reported(mapk_path):
    result = biotsr.analyze(
        mapk_path, method="drug_atom_selection", chain="A",
        entity_name="ESW", entity_id="404", atoms="C3;C5;O6;N8;NOPE1",
    )
    assert result.selection.atoms_missing == ("NOPE1",)
    assert any(s.reason == "requested_atom_not_found" for s in result.skipped)


def test_drug_reads_hetatm_only(trypsin_path):
    """3PTB's BEN is a HETATM; a protein ATOM record must never enter a drug run."""
    result = biotsr.analyze(
        trypsin_path, method="drug", chain="A", entity_name="BEN", entity_id="1"
    )
    assert result.counts.n_points == 9
    assert result.counts.n_triangles == math.comb(9, 3) == 84
    assert {i.residue_name for t in result.triplets for i in t.identities} == {"BEN"}


def test_a_drug_run_does_not_stop_at_the_peptide_TER(trypsin_path):
    """A ligand always sits after the chain's TER; the peptide rule would see none."""
    from biotsr.selection.truncation import FIRST_MODEL_ONLY, LEGACY_SINGLE_PASS

    assert get_method("drug").default_profile.atom_admission.truncation == FIRST_MODEL_ONLY
    assert (
        get_method("protein").default_profile.parser_backend
        and LEGACY_SINGLE_PASS != FIRST_MODEL_ONLY
    )
    result = biotsr.analyze(
        trypsin_path, method="drug", chain="A", entity_name="BEN", entity_id="1"
    )
    assert result.counts.n_points == 9


def test_drug_has_no_completeness_rule(trypsin_path):
    """WebsiteTSR's ligand selector has none, so none is invented."""
    assert get_method("drug").completeness_policy is None
    assert get_method("drug").default_profile.completeness.name == "none"


# ===========================================================================
# Carbohydrate
# ===========================================================================


def test_carbohydrate_monomer_scope(glyco_path):
    """5WT9/A NAG 1: 14 atoms -> C(14,3) = 364, the historical baseline exactly."""
    result = biotsr.analyze(
        glyco_path, method="carbohydrate", chain="A",
        entity_name="NAG", entity_id="1",
    )
    assert result.counts.n_points == 14 == LEGACY_CARBOHYDRATE_ATOM_COUNTS["NAG"]
    assert result.counts.n_triangles == 364
    assert result.selection.residues_used == ("1_NAG",)
    assert _ledger_balances(result)


def test_carbohydrate_polymer_scope_pools_the_whole_chain(glyco_path):
    """5WT9/A: NAG + NAG + FUC -> 38 atoms -> C(38,3) = 8,436, historical baseline."""
    result = biotsr.analyze(glyco_path, method="carbohydrate_chain", chain="A")
    assert result.counts.n_points == 38
    assert result.counts.n_triangles == 8436
    assert len(result.selection.residues_used) == 3
    assert {r.split("_")[1] for r in result.selection.residues_used} <= set(
        CARBOHYDRATE_RESIDUE_NAMES
    )
    # Genuinely inter-residue: triangles span more than one monomer.
    spans = {
        len({i.residue_seq for i in t.identities}) for t in result.triplets
    }
    assert spans & {2, 3}


def test_the_polymer_scope_has_no_completeness_filter(glyco_path):
    """Characterized fact: the legacy grouped caller passes
    apply_ligand_atom_count_filter=False."""
    assert get_method("carbohydrate_chain").completeness_policy is None
    assert (
        get_method("carbohydrate_chain").default_profile.completeness.name == "none"
    )
    assert (
        get_method("carbohydrate").default_profile.completeness.name
        == "exact_legacy_carbohydrate_atom_count"
    )


def test_the_carbohydrate_vocabulary_is_declared_on_the_method_not_hidden():
    """Legacy scoping was JSON-driven and hidden; it is now visible and overridable."""
    recipe = get_method("carbohydrate_chain")
    assert recipe.default_entity_names == CARBOHYDRATE_RESIDUE_NAMES
    assert recipe.default_entity_names == tuple(
        sorted(LEGACY_CARBOHYDRATE_ATOM_COUNTS)
    )


def test_the_carbohydrate_vocabulary_can_be_narrowed_per_call(glyco_path):
    whole = biotsr.analyze(glyco_path, method="carbohydrate_chain", chain="A")
    nag_only = biotsr.analyze(
        glyco_path, method="carbohydrate_chain", chain="A", entity_names=("NAG",)
    )
    assert nag_only.counts.n_points < whole.counts.n_points
    assert {r.split("_")[1] for r in nag_only.selection.residues_used} == {"NAG"}


def test_carbohydrate_atom_selection_narrows_a_monomer(lf9_path):
    """A car_atom_selected cell: five named atoms of one NAG."""
    result = biotsr.analyze(
        lf9_path, method="carbohydrate_atom_selection", chain="D",
        entity_name="NAG", entity_id="1", atoms="C1;C2;N2;O5;O6",
    )
    assert result.counts.n_points == 5
    assert result.counts.n_triangles == math.comb(5, 3) == 10
    assert {i.atom_name for t in result.triplets for i in t.identities} == {
        "C1", "C2", "N2", "O5", "O6",
    }


def test_atom_narrowing_happens_after_completeness(lf9_path):
    """Asking for five atoms of a fourteen-atom NAG must not make it incomplete."""
    result = biotsr.analyze(
        lf9_path, method="carbohydrate_atom_selection", chain="D",
        entity_name="NAG", entity_id="1", atoms="C1;C2;N2;O5;O6",
    )
    assert result.counts.n_points == 5
    assert not [s for s in result.skipped if s.stage == "completeness"]


def test_an_incomplete_carbohydrate_monomer_is_rejected_with_a_reason(glyco_path):
    """A name outside the table is rejected whatever its atom count."""
    result = biotsr.analyze(
        glyco_path, method="carbohydrate", chain="A",
        entity_name="FUC", entity_id="3",
    )
    # 5WT9/A's FUC has ten atoms, which is what the table expects.
    if result.counts.n_points:
        assert result.counts.n_points == LEGACY_CARBOHYDRATE_ATOM_COUNTS["FUC"]
    else:
        assert any(s.stage == "completeness" for s in result.skipped)


def test_a_carbohydrate_run_does_not_stop_at_TER(lf9_path):
    """CarbohydratesTSR's own comment: glycan HETATMs appear after the TER."""
    result = biotsr.analyze(lf9_path, method="carbohydrate_chain", chain="D")
    assert result.counts.n_points > 0


# ===========================================================================
# Nucleotide
# ===========================================================================


def test_nucleotide_monomer_scope(bdna_path):
    result = biotsr.analyze(
        bdna_path, method="nucleotide", chain="A", entity_name="DC", entity_id="3"
    )
    assert result.counts.n_points == 19 == LEGACY_NUCLEOTIDE_ATOM_COUNTS["DC"]
    assert result.counts.n_triangles == math.comb(19, 3) == 969
    assert _ledger_balances(result)


def test_the_exact_completeness_rule_drops_the_terminal_nucleotide(bdna_path):
    """1BNA/A's 5'-terminal DC 1 is missing its phosphate: 12 seen, 11 admitted."""
    targets = biotsr.list_entity_targets(bdna_path, chain="A", method="nucleotide")
    assert len(targets) == 12
    admitted, rejected = [], []
    for seq, name in targets:
        result = biotsr.analyze(
            bdna_path, method="nucleotide", chain="A",
            entity_name=name, entity_id=seq, keep_triplets=False,
        )
        (admitted if result.counts.n_points else rejected).append((seq, name))
    assert len(admitted) == 11
    assert rejected == [("1", "DC")]


def test_a_rejected_nucleotide_records_its_reason(bdna_path):
    result = biotsr.analyze(
        bdna_path, method="nucleotide", chain="A", entity_name="DC", entity_id="1"
    )
    assert result.counts.n_points == 0
    reasons = [s for s in result.skipped if s.stage == "completeness"]
    assert reasons and "expected exactly 19" in reasons[0].detail


def test_dna_and_rna_are_distinguished_only_by_residue_name():
    """The only distinction the legacy source makes is the completeness table."""
    assert set(DNA_RESIDUE_NAMES) == {"DA", "DC", "DG", "DT"}
    assert set(RNA_RESIDUE_NAMES) == {"A", "C", "G", "U"}
    assert set(LEGACY_NUCLEOTIDE_ATOM_COUNTS) == set(DNA_RESIDUE_NAMES) | set(
        RNA_RESIDUE_NAMES
    )


def test_nucleotide_chain_scope_pools_the_strand(bdna_path):
    """STANDARD-ONLY: no legacy implementation pools nucleotide residues."""
    result = biotsr.analyze(
        bdna_path, method="nucleotide_chain", chain="A", keep_triplets=False
    )
    # All 12 nucleotides with the atoms they have (ADR-114): 227 points when
    # DC 1, which lacks its 5' phosphate, was dropped by its atom count; 243 now.
    assert result.counts.n_points == 243
    assert result.counts.n_triangles == math.comb(243, 3)
    assert len(result.selection.residues_used) == 12
    assert not result.annotations["legacy_compatibility_profile"]


def test_a_dna_chain_run_can_be_scoped_to_dna_names(bdna_path):
    whole = biotsr.analyze(
        bdna_path, method="nucleotide_chain", chain="A", keep_triplets=False
    )
    dna = biotsr.analyze(
        bdna_path, method="nucleotide_chain", chain="A",
        entity_names=DNA_RESIDUE_NAMES, keep_triplets=False,
    )
    assert dna.counts.n_points == whole.counts.n_points     # 1BNA is all DNA
    # an entity_names override that matches nothing is refused, not an empty result
    with pytest.raises(biotsr.InvalidRequest, match="none of the requested entities"):
        biotsr.analyze(
            bdna_path, method="nucleotide_chain", chain="A",
            entity_names=RNA_RESIDUE_NAMES, keep_triplets=False,
        )


def test_nucleotide_atom_selection(ac0_path):
    """An nt_atom_selected cell: seven named atoms of one DT, primes kept."""
    result = biotsr.analyze(
        ac0_path, method="nucleotide_atom_selection", chain="E",
        entity_name="DT", entity_id="7", atoms="P;C1';N1;C2;O4;C7;C6",
    )
    assert result.counts.n_points == 7
    assert result.counts.n_triangles == math.comb(7, 3) == 35
    assert {i.atom_name for t in result.triplets for i in t.identities} == {
        "P", "C1'", "N1", "C2", "O4", "C7", "C6",
    }


def test_nucleotide_reads_atom_records_only(ac0_path):
    result = biotsr.analyze(
        ac0_path, method="nucleotide", chain="E", entity_name="DT", entity_id="7"
    )
    assert {i.residue_name for t in result.triplets for i in t.identities} == {"DT"}


# ===========================================================================
# Modifiers, shared across the three families
# ===========================================================================

MODIFIER_CASES = [
    ("drug", "mapk_path", "A", dict(entity_name="ESW", entity_id="404")),
    ("carbohydrate", "glyco_path", "A", dict(entity_name="NAG", entity_id="1")),
    ("carbohydrate_chain", "glyco_path", "A", {}),
    ("nucleotide", "bdna_path", "A", dict(entity_name="DC", entity_id="3")),
]


@pytest.mark.parametrize("method,fixture_name,chain,kwargs", MODIFIER_CASES)
def test_mirror_never_changes_the_triangle_count(method, fixture_name, chain,
                                                 kwargs, request):
    path = request.getfixturevalue(fixture_name)
    plain = biotsr.analyze(path, method=method, chain=chain, **kwargs)
    mirrored = biotsr.analyze(path, method=method, chain=chain, mirror=True, **kwargs)
    assert plain.counts.n_triangles == mirrored.counts.n_triangles
    assert all(
        abs(m.key) == p.key for m, p in zip(mirrored.triplets, plain.triplets)
    )
    assert all(t.key > 0 for t in plain.triplets)
    assert any(t.key < 0 for t in mirrored.triplets)


@pytest.mark.parametrize("method,fixture_name,chain,kwargs", MODIFIER_CASES)
def test_mirror_keys_off_the_pre_fold_angle(method, fixture_name, chain,
                                            kwargs, request):
    path = request.getfixturevalue(fixture_name)
    mirrored = biotsr.analyze(path, method=method, chain=chain, mirror=True, **kwargs)
    for t in mirrored.triplets:
        assert (t.key < 0) == (t.theta_raw > 90)


@pytest.mark.parametrize("method,fixture_name,chain,kwargs", MODIFIER_CASES)
def test_size_filtering_is_strictly_below_the_cutoff(method, fixture_name, chain,
                                                     kwargs, request):
    path = request.getfixturevalue(fixture_name)
    cutoff = 4.0
    plain = biotsr.analyze(path, method=method, chain=chain, **kwargs)
    filtered = biotsr.analyze(
        path, method=method, chain=chain, size_filter=cutoff, **kwargs
    )
    assert filtered.counts.n_filtered_by_size > 0
    assert filtered.counts.n_degenerate == 0
    assert all(t.max_distance < cutoff for t in filtered.triplets)
    assert (
        filtered.counts.n_triangles
        == plain.counts.n_triangles - filtered.counts.n_filtered_by_size
    )


@pytest.mark.parametrize("method,fixture_name,chain,kwargs", MODIFIER_CASES)
def test_the_size_filter_boundary_excludes_a_triangle_exactly_at_the_cutoff(
    method, fixture_name, chain, kwargs, request
):
    path = request.getfixturevalue(fixture_name)
    plain = biotsr.analyze(path, method=method, chain=chain, **kwargs)
    distances = sorted({t.max_distance for t in plain.triplets})
    # A cutoff placed exactly on an observed max_distance, in the middle of the
    # range so both sides of the boundary have triangles.
    observed = distances[len(distances) // 2]
    n_at_boundary = sum(1 for t in plain.triplets if t.max_distance == observed)

    at = biotsr.analyze(path, method=method, chain=chain, size_filter=observed, **kwargs)
    above = biotsr.analyze(
        path, method=method, chain=chain,
        size_filter=math.nextafter(observed, math.inf), **kwargs,
    )
    # Strictly below: a triangle *at* the cutoff is excluded...
    assert all(t.max_distance < observed for t in at.triplets)
    # ...and one ulp above the cutoff admits exactly the boundary triangles.
    assert above.counts.n_triangles == at.counts.n_triangles + n_at_boundary


@pytest.mark.parametrize("method,fixture_name,chain,kwargs", MODIFIER_CASES)
def test_modifiers_compose_and_stay_deterministic(method, fixture_name, chain,
                                                  kwargs, request):
    path = request.getfixturevalue(fixture_name)
    seen = {}
    for mirror in (False, True):
        for cutoff in (None, 5.0):
            keys = [
                t.key
                for t in biotsr.analyze(
                    path, method=method, chain=chain,
                    mirror=mirror, size_filter=cutoff, **kwargs
                ).triplets
            ]
            again = [
                t.key
                for t in biotsr.analyze(
                    path, method=method, chain=chain,
                    mirror=mirror, size_filter=cutoff, **kwargs
                ).triplets
            ]
            assert keys == again
            seen[(mirror, cutoff)] = keys
    assert len(seen) == 4


# ===========================================================================
# Shared-architecture invariants
# ===========================================================================


@pytest.mark.parametrize("name", ATOMIC_FAMILY_METHODS)
def test_no_family_has_its_own_profile_class(name):
    from biotsr.profiles.standard import Profile

    assert type(get_method(name).default_profile) is Profile


def test_no_family_specific_geometry_order_or_encoder_class_exists():
    """ADR-027, enforced structurally rather than by convention."""
    import pathlib

    forbidden = (
        "DrugGeometry", "CarbohydrateGeometry", "NucleotideGeometry",
        "DrugOrdering", "CarbohydrateOrdering", "NucleotideOrdering",
        "DrugEncoder", "CarbohydrateEncoder", "NucleotideEncoder",
        "DrugTSR", "CarbohydrateTSR", "NucleotideTSR",
        "AminoAcidAtomicGeometry",
    )
    root = pathlib.Path(biotsr.__file__).parent
    offenders = []
    for path in root.rglob("*.py"):
        text = path.read_text()
        for symbol in forbidden:
            if ("class %s" % symbol) in text or ("def %s" % symbol) in text:
                offenders.append("%s: %s" % (path.name, symbol))
    assert offenders == []


def test_every_family_shares_the_one_geometry_and_key_formula():
    from biotsr.encoding.key import KEY_FORMULA_VERSION
    from biotsr.geometry.ordering import CANONICAL_ORDERING_VERSION

    spaces = [
        get_method(n).default_profile.key_space() for n in ATOMIC_FAMILY_METHODS
    ]
    assert {s.canonical_ordering_version for s in spaces} == {CANONICAL_ORDERING_VERSION}
    assert {s.key_formula_version for s in spaces} == {KEY_FORMULA_VERSION}
    assert {s.sqrt_convention for s in spaces} == {"math_sqrt"}
    assert {s.angle_convention for s in spaces} == {"pi:3.141592653589793"}


def test_every_atomic_family_produces_comparable_keys():
    """One key space, so drug, carbohydrate and nucleotide keys can be combined."""
    fingerprints = {
        get_method(n).default_profile.key_space().fingerprint
        for n in ATOMIC_FAMILY_METHODS
    }
    assert len(fingerprints) == 1


def test_an_atomic_result_cannot_be_combined_with_a_calpha_one(glyco_path,
                                                               mapk_path):
    from biotsr.domain.keyspace import IncompatibleKeySpace

    carb = biotsr.analyze(
        glyco_path, method="carbohydrate", chain="A",
        entity_name="NAG", entity_id="1", keep_triplets=False,
    )
    protein = biotsr.analyze(
        mapk_path, method="aa_selection", chain="A", residues="26-28;192;207",
        keep_triplets=False,
    )
    with pytest.raises(IncompatibleKeySpace):
        carb.key_space.require_compatible(protein.key_space)


@pytest.mark.parametrize("name", ATOMIC_FAMILY_METHODS)
def test_every_point_carries_full_physical_identity(name):
    """Never a bare atom name -- the demonstrated cross-family defect."""
    from biotsr.domain.identity import AtomIdentity

    assert AtomIdentity.__dataclass_fields__.keys() >= {
        "structure", "model", "chain", "residue_name", "residue_seq",
        "insertion_code", "atom_name", "alt_loc",
    }


def test_two_atoms_with_the_same_name_in_different_entities_stay_distinct(glyco_path):
    """Pooling across residues is exactly where a name-keyed pool collapses."""
    result = biotsr.analyze(glyco_path, method="carbohydrate_chain", chain="A")
    identities = [
        i for t in result.triplets for i in t.identities
    ]
    by_name = {}
    for identity in identities:
        by_name.setdefault(identity.atom_name, set()).add(identity.as_tuple())
    repeated = {n: v for n, v in by_name.items() if len(v) > 1}
    assert repeated, "the polymer fixture should contain repeated atom names"


def test_streaming_works_for_every_family(tmp_path, mapk_path, glyco_path,
                                          bdna_path):
    from biotsr.io.native import NativeTripletStreamWriter

    cases = [
        ("drug", mapk_path, "A", dict(entity_name="ESW", entity_id="404")),
        ("carbohydrate", glyco_path, "A", dict(entity_name="NAG", entity_id="1")),
        ("nucleotide", bdna_path, "A", dict(entity_name="DC", entity_id="3")),
    ]
    for method, path, chain, kwargs in cases:
        out = tmp_path / ("%s.txt" % method)
        with NativeTripletStreamWriter(str(out)) as writer:
            result = biotsr.analyze(
                path, method=method, chain=chain, keep_triplets=False,
                triplet_sink=writer.write, **kwargs
            )
        assert result.triplets is None
        assert writer.n_written == result.counts.n_triangles
