"""Comparison, alternate-location and molecule-instance rules, end to end on real structures.

* comparison (ADR-094): one key definition always; processing differences only
  with an explicit opt-in; what was selected never blocks a comparison;
* alternate locations (ADR-095): blank and A kept, B and later excluded, no
  fallback to another conformer;
* molecule instances (ADR-096): chain + name + identifier is one copy; chain +
  name is every copy on that chain, in one calculation; identifiers are PDB
  residue numbers with an optional insertion code, never positions.
"""

from __future__ import annotations

import os

import pytest

import biotsr
from biotsr.analysis import (
    PartialSamples,
    AnalysisSample,
    IncompatibleProcessing,
    clustermap_spec,
    common_key_intersection,
    distance_matrix,
    feature_matrix,
    key_occurrence_percent,
)
from biotsr.domain.errors import InvalidRequest
from biotsr.domain.keyspace import IncompatibleKeySpace
from biotsr.domain.selection_record import SCOPE_ALL_COPIES, SCOPE_ONE_INSTANCE

from ..conftest import STRUCTURES


def _pdb(name):
    return os.path.join(STRUCTURES, name + ".pdb")


def _sample(sample_id, **kwargs):
    api = biotsr.analyze_cross if "chain_2" in kwargs else biotsr.analyze
    return AnalysisSample.from_result(api(keep_triplets=False, **kwargs),
                                      sample_id=sample_id)


# ===========================================================================
# Comparison
# ===========================================================================


def test_lys_only_and_tyr_only_compare_without_any_opt_in():
    """The agreed case: different amino-acid types, one key
    definition, one processing -- compared, with no warning."""
    lys = _sample("lys", source=_pdb("1UBQ"), method="aa_ind", chain="A",
                  residue_types="LYS")
    tyr = _sample("tyr", source=_pdb("1UBQ"), method="aa_ind", chain="A",
                  residue_types="TYR")
    report = common_key_intersection([lys, tyr])
    assert report.warnings == ()
    assert report.comparison["processing_uniform"] is True
    assert [s["selection"][0]["arguments"] for s in report.samples] == [
        {"residue_types": ["LYS"]}, {"residue_types": ["TYR"]}]
    matrix = feature_matrix([lys, tyr])
    assert matrix.warnings == () and matrix.n_samples == 2


def test_an_empty_intersection_is_a_valid_result():
    a = _sample("ben", source=_pdb("3PTB"), method="drug", chain="A",
                entity_name="BEN", entity_id="1")
    b = _sample("esw", source=_pdb("6G9M"), method="drug_atom_selection",
                chain="A", entity_name="ESW", entity_id="404", atoms="CL1;C3;C5")
    report = common_key_intersection([a, b])
    assert report.is_empty
    assert [r.total_keys for r in report.rows] == [len(a.counts), len(b.counts)]


def test_different_methods_sharing_a_key_space_compare():
    """A ligand and a glycan, two methods: the same 112-code atom key space and
    the same processing, so comparable; the methods are recorded."""
    a = _sample("ben", source=_pdb("3PTB"), method="drug", chain="A",
                entity_name="BEN", entity_id="1")
    b = _sample("esw", source=_pdb("6G9M"), method="drug_atom_selection", chain="A",
                entity_name="ESW", entity_id="404", atoms="C3;C5;O6;N8")
    report = key_occurrence_percent([a, b])
    assert report.warnings == ()
    assert report.comparison["methods"] == ["drug", "drug_atom_selection"]


def test_the_same_method_with_another_key_definition_is_refused_even_with_the_opt_in():
    plain = _sample("plain", source=_pdb("1UBQ"), method="protein", chain="A")
    mirrored = _sample("mirror", source=_pdb("1UBQ"), method="protein", chain="A",
                       mirror=True)
    for allow in (False, True):
        with pytest.raises(IncompatibleKeySpace, match="mirror_convention"):
            feature_matrix([plain, mirrored], allow_mixed_processing=allow)


def test_different_label_spaces_are_refused_whatever_the_method_names():
    calpha = _sample("ca", source=_pdb("1UBQ"), method="protein", chain="A")
    atomic = _sample("at", source=_pdb("1UBQ"), method="aa_ind", chain="A",
                     residue_types="LYS")
    with pytest.raises(IncompatibleKeySpace, match="label_scheme"):
        common_key_intersection([calpha, atomic], allow_mixed_processing=True)


@pytest.mark.parametrize("first, second, field", [
    (dict(residue_types="LYS"), dict(residue_types="TYR", size_filter=10.0),
     "size_filter"),
])
def test_a_processing_difference_needs_acknowledgement(first, second, field):
    a = _sample("a", source=_pdb("1UBQ"), method="aa_ind", chain="A", **first)
    b = _sample("b", source=_pdb("1UBQ"), method="aa_ind", chain="A", **second)
    with pytest.raises(IncompatibleProcessing, match=field):
        common_key_intersection([a, b])
    report = common_key_intersection([a, b], allow_mixed_processing=True)
    assert report.comparison["processing_differences"] == [field]
    assert field in report.warnings[0]


def test_a_contact_cutoff_difference_needs_acknowledgement_and_chains_do_not():
    base = dict(source=_pdb("6G9M"), method="drug_protein_cross", chain="A",
                chain_2="A", entity_name="ESW", entity_id="404")
    a = _sample("a", **base)
    b = _sample("b", contact_cutoff=4.5, **base)
    with pytest.raises(IncompatibleProcessing, match="contact_cutoff"):
        feature_matrix([a, b])
    other = _sample("ben", source=_pdb("3PTB"), method="drug_protein_cross",
                    chain="A", chain_2="A", entity_name="BEN", entity_id="1")
    assert a.sampling.fingerprint != other.sampling.fingerprint
    # 6G9M's whole protein side holds CME 161, which no method reads, so ``a``
    # is partial and only that needs acknowledging (ADR-106); the different
    # structure, ligand, copy and chains need nothing.
    with pytest.raises(PartialSamples, match="CME_161"):
        feature_matrix([a, other])
    matrix = feature_matrix([a, other], allow_partial=True)
    assert matrix.comparison["processing_uniform"] is True
    assert list(matrix.comparison["partial_samples"]) == ["a"]
    assert len(matrix.warnings) == 1 and "explicitly allowed" in matrix.warnings[0]


def test_the_comparison_records_selection_settings_fingerprints_and_warnings():
    a = _sample("nag_all", source=_pdb("7LF9"), method="carbohydrate", chain="D",
                entity_name="NAG")
    b = _sample("man_all", source=_pdb("7LF9"), method="carbohydrate", chain="D",
                entity_name="MAN")
    matrix = feature_matrix([a, b])
    record = matrix.samples[0]
    assert record["method"] == "carbohydrate"
    assert record["selection"][0]["scope"] == SCOPE_ALL_COPIES
    assert record["selection"][0]["instances"] == ["NAG_1", "NAG_2", "NAG_5"]
    assert record["processing"]["altloc_policy"] == "blank_or_a"
    assert record["processing_fingerprint"] == matrix.samples[1]["processing_fingerprint"]
    assert record["key_space_fingerprint"] and record["sampling_fingerprint"]
    assert "result_warnings" in record


def test_the_default_tree_and_heat_map_need_no_cluster_count():
    samples = [_sample(t, source=_pdb("1UBQ"), method="aa_ind", chain="A",
                       residue_types=t) for t in ("LYS", "ARG", "GLU", "LEU", "ILE")]
    spec = clustermap_spec(distance_matrix(feature_matrix(samples)))
    assert spec.assignment is None and spec.adjusted_rand is None
    assert sorted(spec.ordered_sample_ids()) == sorted(s.sample_id for s in samples)
    assert len(spec.hierarchy.linkage.rows) == 4
    cut = clustermap_spec(distance_matrix(feature_matrix(samples)), height=0.9)
    assert cut.assignment.criterion == "distance"


# ===========================================================================
# Alternate locations
# ===========================================================================


def _altloc_structure(tmp_path):
    """LEU 2: side chain in A, B and C.  LEU 3: side chain only in B and C."""
    def atom(serial, name, alt, seq, xyz, element):
        return ("ATOM  %5d %-4s%1s%3s A%4d    %8.3f%8.3f%8.3f  1.00 10.00          %2s\n"
                % (serial, " " + name if len(name) < 4 else name, alt, "LEU", seq,
                   xyz[0], xyz[1], xyz[2], element))
    lines, n = [], 1
    for seq, alts in ((2, "ABC"), (3, "BC")):
        base = 6.0 * seq
        for name, el, d in (("N", "N", (0, 0, 0)), ("CA", "C", (1.46, 0, 0)),
                            ("C", "C", (2.0, 1.4, 0)), ("O", "O", (1.3, 2.4, 0))):
            lines.append(atom(n, name, " ", seq, (base + d[0], d[1], d[2]), el))
            n += 1
        for k, alt in enumerate(alts):
            for name, el, d in (("CB", "C", (2.0, -0.8, 1.2)), ("CG", "C", (3.4, -0.6, 1.8)),
                                ("CD1", "C", (3.8, -1.9, 2.5)), ("CD2", "C", (3.6, 0.6, 2.8))):
                lines.append(atom(n, name, alt, seq, (base + d[0] + 0.4 * k,
                                                      d[1] + 0.4 * k, d[2] - 0.4 * k), el))
                n += 1
    path = tmp_path / "ALTC.pdb"
    path.write_text("".join(lines) + "END\n")
    return str(path)


def test_blank_and_a_are_kept_and_b_and_later_are_excluded(tmp_path):
    """Before 1.0.0rc4 amino_acid_atomic kept C conformers too, so LEU 2 had 12
    atoms and was excluded as incomplete; now it is blank + A = 8 atoms."""
    path = _altloc_structure(tmp_path)
    result = biotsr.analyze(path, method="amino_acid_atomic", chain="A",
                            residue="2", residue_name="LEU", keep_triplets=False)
    assert result.counts.n_points == 8
    assert result.source.altloc_policy == "blank_or_a"
    assert {s.reason for s in result.skipped} == {"altloc_rejected"}
    assert any("excluded 8 atom record(s) with indicator B, C" in w
               for w in result.warnings)


def test_there_is_no_fallback_when_a_is_absent(tmp_path):
    """LEU 3's side chain exists only as B and C: none of it is used, and the
    completeness rule then excludes the residue, saying why. (1.0.0rc3 built it
    from conformer C.)"""
    path = _altloc_structure(tmp_path)
    result = biotsr.analyze(path, method="amino_acid_atomic", chain="A",
                            residue="3", residue_name="LEU", keep_triplets=False)
    assert result.counts.n_points == 0
    assert any("expected exactly 8 atoms, got 4" in s.detail for s in result.skipped
               if s.stage == "completeness")
    pooled = biotsr.analyze(path, method="aa_ind", chain="A", residue_types="LEU",
                            keep_triplets=False)
    assert list(pooled.selection.residues_used) == ["2_LEU"]
    assert any("3_LEU" in w for w in pooled.warnings)


@pytest.mark.parametrize("method", ["protein", "aa_selection", "amino_acid_atomic",
                                    "amino_acid_atom_selection", "aa_ind", "drug",
                                    "carbohydrate", "nucleotide_chain",
                                    "protein_protein_cross"])
def test_every_method_records_the_one_altloc_rule(method):
    kwargs = {
        "protein": dict(source=_pdb("1UBQ"), chain="A"),
        "aa_selection": dict(source=_pdb("1UBQ"), chain="A", residues="1-10"),
        "amino_acid_atomic": dict(source=_pdb("1UBQ"), chain="A", residue="15",
                                  residue_name="LEU"),
        "amino_acid_atom_selection": dict(source=_pdb("1UBQ"), chain="A",
                                          atoms="LYS_48_NZ;LEU_8_CD1;ILE_44_CD1"),
        "aa_ind": dict(source=_pdb("1UBQ"), chain="A", residue_types="LYS"),
        "drug": dict(source=_pdb("3PTB"), chain="A", entity_name="BEN", entity_id="1"),
        "carbohydrate": dict(source=_pdb("7LF9"), chain="D", entity_name="NAG",
                             entity_id="1"),
        "nucleotide_chain": dict(source=_pdb("1BNA"), chain="A"),
        "protein_protein_cross": dict(source=_pdb("1BRS"), chain="A", chain_2="D",
                                      residue_types="LYS", residue_types_2="ASP"),
    }[method]
    api = biotsr.analyze_cross if "chain_2" in kwargs else biotsr.analyze
    result = api(method=method, keep_triplets=False, **kwargs)
    assert result.source.altloc_policy == "blank_or_a"


# ===========================================================================
# Molecule instances and identifiers
# ===========================================================================


def test_name_and_identifier_select_that_one_copy():
    result = biotsr.analyze(_pdb("7LF9"), method="carbohydrate", chain="D",
                            entity_name="NAG", entity_id="5", keep_triplets=False)
    side = result.annotations["selection_request"]["sides"][0]
    assert side["scope"] == SCOPE_ONE_INSTANCE
    assert side["instances"] == ["NAG_5"]
    assert result.annotations["target"] == "NAG_5"
    assert list(result.selection.residues_used) == ["5_NAG"]


def test_a_name_without_identifier_is_every_copy_in_one_calculation():
    """NAG 1, 2 and 5 of 7LF9 chain D, pooled: 42 atoms, and triangles that span
    the copies -- more than the three single-copy calculations together."""
    every = biotsr.analyze(_pdb("7LF9"), method="carbohydrate", chain="D",
                           entity_name="NAG", keep_triplets=False)
    singles = [biotsr.analyze(_pdb("7LF9"), method="carbohydrate", chain="D",
                              entity_name="NAG", entity_id=i, keep_triplets=False)
               for i in ("1", "2", "5")]
    assert every.counts.n_points == sum(s.counts.n_points for s in singles) == 42
    assert every.counts.n_triangles > sum(s.counts.n_triangles for s in singles)
    side = every.annotations["selection_request"]["sides"][0]
    assert side["scope"] == SCOPE_ALL_COPIES
    assert side["instances"] == ["NAG_1", "NAG_2", "NAG_5"]
    assert every.annotations["target"] == "NAG_allcopies"


def test_all_copies_of_one_sugar_is_not_the_whole_glycan():
    """Every MAN of chain D (4 and 7) versus the entire glycan, NAG, BMA, MAN
    and GAL together."""
    man = biotsr.analyze(_pdb("7LF9"), method="carbohydrate", chain="D",
                         entity_name="MAN", keep_triplets=False)
    glycan = biotsr.analyze(_pdb("7LF9"), method="carbohydrate_chain", chain="D",
                            keep_triplets=False)
    assert list(man.selection.residues_used) == ["4_MAN", "7_MAN"]
    assert {r.split("_")[1] for r in glycan.selection.residues_used} == {
        "NAG", "BMA", "MAN", "GAL"}
    assert glycan.annotations["selection_request"]["sides"][0]["scope"] == (
        "chain_vocabulary")


def test_an_atom_subset_applies_to_each_copy_and_missing_atoms_are_named_per_copy():
    result = biotsr.analyze(_pdb("7LF9"), method="carbohydrate_atom_selection",
                            chain="D", entity_name="NAG", atoms="C1;C2;XX9",
                            keep_triplets=False)
    assert result.counts.n_points == 6                   # C1, C2 of three copies
    side = result.annotations["selection_request"]["sides"][0]
    assert side["atom_subset"]["applies_to"] == "each_instance"
    assert side["atom_subset"]["missing"] == {
        "NAG_1": ["XX9"], "NAG_2": ["XX9"], "NAG_5": ["XX9"]}
    assert any("XX9 in NAG_5" in w for w in result.warnings)


def test_a_cross_entity_side_follows_the_same_policy():
    one = biotsr.analyze_cross(_pdb("7LF9"), method="carbohydrate_protein_cross",
                               chain="D", chain_2="B", entity_name="NAG",
                               entity_id="1", keep_triplets=False)
    every = biotsr.analyze_cross(_pdb("7LF9"), method="carbohydrate_protein_cross",
                                 chain="D", chain_2="B", entity_name="NAG",
                                 keep_triplets=False)
    assert every.annotations["target"] == "NAG_allcopies"
    assert "entity=NAG(all copies)" in every.sampling.side_a_selector
    assert every.counts.n_points >= one.counts.n_points
    entity = [s for s in every.annotations["selection_request"]["sides"]
              if s["role"] == "carbohydrate"][0]
    assert entity["scope"] == SCOPE_ALL_COPIES


@pytest.mark.parametrize("entity_id, message", [
    ("3", "copies of NAG there have ids 1, 2, 5"),   # BMA 3 is not a NAG
    ("1A", "copies of NAG there have ids 1, 2, 5"),  # no insertion code there
    ("second", "not a PDB residue identifier"),      # never a position
    ("1 A", "not a PDB residue identifier"),
])
def test_an_identifier_is_the_residue_number_in_the_file(entity_id, message):
    with pytest.raises(InvalidRequest, match=message):
        biotsr.analyze(_pdb("7LF9"), method="carbohydrate", chain="D",
                       entity_name="NAG", entity_id=entity_id, keep_triplets=False)


def test_leading_zeros_name_the_same_residue():
    a = biotsr.analyze(_pdb("7LF9"), method="carbohydrate", chain="D",
                       entity_name="NAG", entity_id="005", keep_triplets=False)
    assert a.annotations["target"] == "NAG_5"


def _insertion_structure(tmp_path, split=False):
    """Two NAG copies numbered 52 and 52A on chain X (and, with ``split``, the
    identifier 52 used twice, non-contiguously)."""
    def atom(serial, name, seq, icode, x, el):
        return ("HETATM%5d %-4s %3s X%4d%1s   %8.3f%8.3f%8.3f  1.00 10.00          %2s\n"
                % (serial, " " + name, "NAG", seq, icode, x, 0.5 * serial, 0.1 * serial, el))
    names = [("C1", "C"), ("C2", "C"), ("O5", "O"), ("N2", "N")]
    lines, n = [], 1
    blocks = [(52, ""), (52, "A")] + ([(53, ""), (52, "")] if split else [])
    for seq, icode in blocks:
        for name, el in names:
            lines.append(atom(n, name, seq, icode, 1.3 * n, el))
            n += 1
    path = tmp_path / ("ICODE%s.pdb" % ("S" if split else ""))
    path.write_text("".join(lines) + "END\n")
    return str(path)


def test_an_insertion_code_is_part_of_the_identifier(tmp_path):
    path = _insertion_structure(tmp_path)
    plain = biotsr.analyze(path, method="drug", chain="X", entity_name="NAG",
                           entity_id="52", keep_triplets=False)
    coded = biotsr.analyze(path, method="drug", chain="X", entity_name="NAG",
                           entity_id="52A", keep_triplets=False)
    both = biotsr.analyze(path, method="drug", chain="X", entity_name="NAG",
                          keep_triplets=False)
    assert list(plain.selection.residues_used) == ["52_NAG"]
    assert list(coded.selection.residues_used) == ["52A_NAG"]
    assert both.annotations["selection_request"]["sides"][0]["instances"] == [
        "NAG_52", "NAG_52A"]


def test_an_identifier_used_for_two_separate_molecules_is_refused(tmp_path):
    path = _insertion_structure(tmp_path, split=True)
    for entity_id in ("52", None):
        with pytest.raises(InvalidRequest, match="two separate places"):
            biotsr.analyze(path, method="drug", chain="X", entity_name="NAG",
                           entity_id=entity_id, keep_triplets=False)


def test_the_chain_is_never_inferred():
    with pytest.raises(TypeError):
        biotsr.analyze(_pdb("7LF9"), method="carbohydrate", entity_name="NAG")
    with pytest.raises(InvalidRequest, match="chains in this structure"):
        biotsr.analyze(_pdb("7LF9"), method="carbohydrate", chain="Q",
                       entity_name="NAG")
    # Chain B has no NAG; the copies on chains A, C and D are never borrowed.
    with pytest.raises(InvalidRequest, match="not in chain B"):
        biotsr.analyze(_pdb("7LF9"), method="carbohydrate", chain="B",
                       entity_name="NAG", keep_triplets=False)
    # Chain A's single NAG is chain A's: all copies means on that chain only.
    only_a = biotsr.analyze(_pdb("7LF9"), method="carbohydrate", chain="A",
                            entity_name="NAG", keep_triplets=False)
    assert {r.split("_")[1] for r in only_a.selection.residues_used} == {"NAG"}
    assert only_a.selection.chains_used == ("A",)
    with pytest.raises(ValueError, match="chain_2"):
        biotsr.analyze_cross(_pdb("7LF9"), method="carbohydrate_protein_cross",
                             chain="D", entity_name="NAG")


def test_a_residue_type_without_an_occurrence_runs_each_residue_of_that_name():
    """amino_acid_atomic is per residue: LEU alone is one calculation per
    leucine; a single occurrence needs its identifier."""
    every = biotsr.analyze_entities(_pdb("1UBQ"), method="amino_acid_atomic",
                                    chain="A", residue_name="LEU",
                                    keep_triplets=False)
    assert list(every) == ["8_LEU", "15_LEU", "43_LEU", "50_LEU", "56_LEU",
                           "67_LEU", "69_LEU", "71_LEU", "73_LEU"]
    one = biotsr.analyze(_pdb("1UBQ"), method="amino_acid_atomic", chain="A",
                         residue="15", residue_name="LEU", keep_triplets=False)
    assert one.frequencies == every["15_LEU"].frequencies
    with pytest.raises(ValueError, match="needs a target residue"):
        biotsr.analyze(_pdb("1UBQ"), method="amino_acid_atomic", chain="A",
                       residue_name="LEU")
