"""The comparison rule, end to end: computed results, saved, loaded, compared.

ADR-094: one key definition always; processing differences only when
acknowledged; what was selected never blocks a comparison.  rc3 instead
compared a "sampling rule" that included residue-type sets, entity
vocabularies, atom-subset narrowing and the method, so ``aa_ind`` over LYS and
over TYR needed an opt-in.

Every case here goes the whole way a user's data goes: computed, written with
the native writer (or by the command line), read back with ``load_result``,
reduced to samples, and handed to the feature matrix, common keys, similarity
and clustering.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest

import biotsr
from biotsr.analysis import (
    AnalysisSample,
    IncompatibleProcessing,
    PartialSamples,
    cluster,
    common_key_intersection,
    distance_matrix,
    feature_matrix,
    key_occurrence_percent,
)
from biotsr.analysis.samples import processing_fingerprint
from biotsr.domain.keyspace import IncompatibleKeySpace
from biotsr.io.native import save

from ..conftest import SRC, structure


def _saved(tmp_path, result, name, label=None):
    """Write ``result`` natively, read it back, and reduce it to a sample."""
    out = tmp_path / name
    out.mkdir()
    written = save(result, str(out), triplets=False)
    loaded = biotsr.load_result(written["keys"])
    return AnalysisSample.from_result(loaded, sample_id=name, label=label)


def _all_operations(samples, **kwargs):
    """Every comparing operation, so a rule cannot hold in one and not another."""
    fm = feature_matrix(samples, **kwargs)
    common_key_intersection(samples, **kwargs)
    key_occurrence_percent(samples, **kwargs)
    dm = distance_matrix(fm)
    return fm, dm, cluster(dm, 1)


def _compared_freely(samples):
    fm, dm, assignment = _all_operations(samples)
    assert fm.warnings == () and dm.warnings == () and assignment.warnings == ()
    assert fm.comparison["processing_uniform"] is True
    return fm


def _refused_everywhere(samples, match):
    for op in (lambda: feature_matrix(samples),
               lambda: common_key_intersection(samples),
               lambda: key_occurrence_percent(samples)):
        with pytest.raises(IncompatibleProcessing, match=match):
            op()
    # The acknowledgement admits them, and every downstream product says so.
    fm, dm, assignment = _all_operations(samples, allow_mixed_processing=True)
    assert fm.comparison["mixed_processing_allowed"] is True
    assert fm.warnings and dm.warnings and assignment.warnings
    return fm


# ---------------------------------------------------------------------------
# Selections never block a comparison
# ---------------------------------------------------------------------------


def test_calpha_chains_of_different_structures_compare_without_any_warning(tmp_path):
    samples = [
        _saved(tmp_path, biotsr.analyze(structure(pdb), method="protein", chain=ch,
                                         keep_triplets=False), "%s_%s" % (pdb, ch))
        for pdb, ch in (("1UBQ", "A"), ("1BRS", "A"), ("1BRS", "D"))
    ]
    fm = _compared_freely(samples)
    assert fm.comparison["n_processing_policies"] == 1
    assert [p["selection"][0]["chain"] for p in fm.samples] == ["A", "A", "D"]


def test_lys_only_and_tyr_only_compare_without_an_opt_in(tmp_path):
    lys = _saved(tmp_path, biotsr.analyze(structure("1UBQ"), method="aa_ind", chain="A",
                                          residue_types="LYS", keep_triplets=False), "lys")
    tyr = _saved(tmp_path, biotsr.analyze(structure("1BRS"), method="aa_ind", chain="A",
                                          residue_types="TYR", keep_triplets=False), "tyr")
    fm = _compared_freely([lys, tyr])
    assert [p["selection"][0]["arguments"] for p in fm.samples] == [
        {"residue_types": ["LYS"]}, {"residue_types": ["TYR"]}]


def test_one_ligand_cross_across_structures_ligands_and_copies(tmp_path):
    """BEN in 3PTB and ESW 404 in 6G9M: different structures, ligands, copies
    and chains -- compared with no opt-in, each target reported."""
    ben = biotsr.analyze_cross(structure("3PTB"), method="drug_protein_cross",
                               chain="A", chain_2="A", entity_name="BEN",
                               entity_id="1", keep_triplets=False)
    esw = biotsr.analyze_cross(structure("6G9M"), method="drug_protein_cross",
                               chain="A", chain_2="A", entity_name="ESW",
                               entity_id="404", keep_triplets=False)
    a, b = _saved(tmp_path, ben, "3PTB_BEN"), _saved(tmp_path, esw, "6G9M_ESW")
    assert a.sampling.fingerprint != b.sampling.fingerprint
    assert processing_fingerprint(a) == processing_fingerprint(b)
    # 6G9M's whole protein side holds CME 161, which no method reads: the
    # sample is partial, and that alone needs acknowledging (ADR-106) -- the
    # different structures, ligands, copies and chains still need nothing.
    with pytest.raises(PartialSamples, match="side a: not supported or not represented: CME_161"):
        feature_matrix([a, b])
    fm, dm, assignment = _all_operations([a, b], allow_partial=True)
    assert fm.comparison["processing_uniform"] is True
    assert fm.comparison["mixed_processing_allowed"] is False
    assert list(fm.comparison["partial_samples"]) == ["6G9M_ESW"]
    assert all("partial" in w for w in fm.warnings + dm.warnings + assignment.warnings)
    drug_sides = [p["selection"][1] for p in fm.samples]
    assert [s["arguments"]["entity_name"] for s in drug_sides] == ["BEN", "ESW"]
    assert [s["arguments"]["entity_id"] for s in drug_sides] == ["1", "404"]


@pytest.mark.parametrize("first, second", [
    # a ligand's atom subset against the whole ligand
    (dict(), dict(atoms="C;N1")),
    # a narrowed protein side against the whole chain
    (dict(), dict(residue_types="ASP;SER;GLY")),
])
def test_a_narrower_cross_selection_is_a_selection_not_a_processing_change(
        tmp_path, first, second):
    base = dict(method="drug_protein_cross", chain="A", chain_2="A",
                entity_name="BEN", entity_id="1", keep_triplets=False)
    a = biotsr.analyze_cross(structure("3PTB"), **dict(base, **first))
    b = biotsr.analyze_cross(structure("3PTB"), **dict(base, **second))
    _compared_freely([_saved(tmp_path, a, "a"), _saved(tmp_path, b, "b")])


def test_a_different_entity_vocabulary_is_a_selection(tmp_path):
    default = biotsr.analyze(structure("5WT9"), method="carbohydrate_chain", chain="A",
                             keep_triplets=False)
    nag_only = biotsr.analyze(structure("5WT9"), method="carbohydrate_chain", chain="A",
                              entity_names="NAG", keep_triplets=False)
    _compared_freely([_saved(tmp_path, default, "default"),
                      _saved(tmp_path, nag_only, "nag")])


def test_a_whole_ligand_and_its_atom_selection_method_compare(tmp_path):
    """Two methods, one key space, the same processing (neither has a
    completeness rule): comparable, and the methods are recorded."""
    whole = biotsr.analyze(structure("3PTB"), method="drug", chain="A",
                           entity_name="BEN", keep_triplets=False)
    some = biotsr.analyze(structure("3PTB"), method="drug_atom_selection", chain="A",
                          entity_name="BEN", atoms="C1;C2;C3;C4", keep_triplets=False)
    fm = _compared_freely([_saved(tmp_path, whole, "whole"),
                           _saved(tmp_path, some, "some")])
    assert fm.comparison["methods"] == ["drug", "drug_atom_selection"]


# ---------------------------------------------------------------------------
# Processing differences: refused by default, named, recorded when acknowledged
# ---------------------------------------------------------------------------


def test_a_one_sided_contact_filter_is_a_processing_difference(tmp_path):
    both = biotsr.analyze_cross(structure("3PTB"), method="drug_protein_cross",
                                chain="A", chain_2="A", entity_name="BEN",
                                keep_triplets=False)
    one = biotsr.analyze_cross(structure("3PTB"), method="drug_protein_cross",
                               chain="A", chain_2="A", entity_name="BEN",
                               contact_filter="side_a_only", keep_triplets=False)
    fm = _refused_everywhere([_saved(tmp_path, both, "both"),
                              _saved(tmp_path, one, "one")], "contact_filter")
    assert "contact_filter" in fm.comparison["processing_differences"]


def test_a_completeness_rule_difference_is_a_processing_difference(tmp_path):
    ligand = biotsr.analyze(structure("3PTB"), method="drug", chain="A",
                            entity_name="BEN", keep_triplets=False)
    sugar = biotsr.analyze(structure("7LF9"), method="carbohydrate", chain="D",
                           entity_name="NAG", entity_id="1", keep_triplets=False)
    _refused_everywhere([_saved(tmp_path, ligand, "ligand"),
                         _saved(tmp_path, sugar, "sugar")], "completeness_policy")


def test_the_opt_in_never_relaxes_the_key_space(tmp_path):
    plain = biotsr.analyze(structure("1UBQ"), method="protein", chain="A",
                           keep_triplets=False)
    mirrored = biotsr.analyze(structure("1UBQ"), method="protein", chain="A",
                              mirror=True, keep_triplets=False)
    pair = [_saved(tmp_path, plain, "plain"), _saved(tmp_path, mirrored, "mirror")]
    for op in (feature_matrix, common_key_intersection, key_occurrence_percent):
        with pytest.raises(IncompatibleKeySpace):
            op(pair, allow_mixed_processing=True)


def test_the_report_records_each_samples_selection_settings_and_warnings(tmp_path):
    ok = biotsr.analyze(structure("1UBQ"), method="aa_ind_specific", chain="A",
                        residues="LYS_6;GLU_18", keep_triplets=False)
    partial = biotsr.analyze(structure("1BRS"), method="aa_ind_specific", chain="A",
                             residues="LYS_27;GLU_500", keep_triplets=False)
    samples = [_saved(tmp_path, ok, "ok"), _saved(tmp_path, partial, "partial")]
    with pytest.raises(PartialSamples, match="partial \\(not found: GLU_500\\)"):
        feature_matrix(samples)
    fm = feature_matrix(samples, allow_partial=True)
    assert len(fm.warnings) == 1 and "explicitly allowed" in fm.warnings[0]
    first, second = fm.samples
    assert (first["partial"], second["partial"]) == (False, True)
    assert second["omissions"]["primary"] == {"not_found": ["GLU_500"]}
    assert first["selection"][0]["arguments"]["residues"] == ["LYS_6", "GLU_18"]
    assert first["selection"][0]["scope"] == "named_residues"
    assert first["selection"][0]["instances"] == ["LYS_6", "GLU_18"]
    assert second["selection"][0]["excluded"] == ["GLU_500"]
    assert any("requested but not found" in w for w in second["result_warnings"])
    assert first["sampling"]["selector"] == "specific_residue_atoms"
    assert first["processing"]["altloc_policy"] == "blank_or_a"
    assert first["sampling_fingerprint"] == ok.sampling.fingerprint
    assert first["key_space_fingerprint"] == ok.key_space.fingerprint
    assert first["processing_fingerprint"] == second["processing_fingerprint"]


def test_a_result_without_a_selection_record_still_compares(tmp_path):
    """The selection record is provenance: nothing in the comparison reads it."""
    samples = []
    for name, rt in (("lys", "LYS"), ("tyr", "TYR")):
        result = biotsr.analyze(structure("1UBQ"), method="aa_ind", chain="A",
                                residue_types=rt, keep_triplets=False)
        written = save(result, str(tmp_path / name), triplets=False)
        with open(written["metadata"]) as fh:
            payload = json.load(fh)
        payload["provenance"].pop("selection_request", None)
        with open(written["metadata"], "w") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
        samples.append(AnalysisSample.from_result(
            biotsr.load_result(written["keys"]), sample_id=name))
    fm = feature_matrix(samples)
    assert fm.warnings == ()
    assert [p["selection"] for p in fm.samples] == [None, None]


# ---------------------------------------------------------------------------
# The command line: the same rule from saved results
# ---------------------------------------------------------------------------


def _cli(args):
    env = dict(os.environ, PYTHONPATH=SRC)
    return subprocess.run([sys.executable, "-m", "biotsr"] + args, env=env,
                          capture_output=True, text=True, timeout=300)


def test_the_command_line_applies_the_rule_to_saved_results(tmp_path):
    out = tmp_path / "out"
    for pdb, rt, extra in (("1UBQ", "LYS", []), ("1BRS", "TYR", []),
                           ("1BRS", "LYS", ["--size-filter", "12"])):
        done = _cli(["run", "--method", "aa_ind", "--pdb", structure(pdb), "--chain",
                     "A", "--residue-types", rt, "--output",
                     str(out / ("sized" if extra else "plain")), "--no-triplets"] + extra)
        assert done.returncode == 0, done.stderr
    same = _cli(["analyze", "common-keys", "--from-results", str(out / "plain"),
                 "--output", str(tmp_path / "a1")])
    assert same.returncode == 0, same.stderr
    mixed = _cli(["analyze", "common-keys", "--from-results", str(out),
                  "--output", str(tmp_path / "a2")])
    assert mixed.returncode == 2
    assert "size_filter" in mixed.stderr and "--allow-mixed-processing" in mixed.stderr
    renamed = _cli(["analyze", "common-keys", "--from-results", str(out),
                    "--output", str(tmp_path / "a3"), "--allow-mixed-sampling"])
    assert renamed.returncode == 2 and "renamed --allow-mixed-processing" in renamed.stderr
    allowed = _cli(["analyze", "common-keys", "--from-results", str(out),
                    "--output", str(tmp_path / "a4"), "--allow-mixed-processing"])
    assert allowed.returncode == 0, allowed.stderr
    report = json.loads((tmp_path / "a4" / "analysis" / "common_keys"
                         / "all_samples.json").read_text())
    assert report["comparison"]["n_processing_policies"] == 2
    assert report["comparison"]["processing_differences"] == ["size_filter"]
    assert any("size_filter" in w for w in report["warnings"])
    assert {s["sample_id"]: s["selection"][0]["arguments"]
            for s in report["samples"]} == {
        "1UBQ_A_aa_ind_LYS": {"residue_types": ["LYS"]},
        "1BRS_A_aa_ind_TYR": {"residue_types": ["TYR"]},
        "1BRS_A_aa_ind_LYS_size12": {"residue_types": ["LYS"]},
    }
