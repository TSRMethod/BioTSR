"""Processing identity on resume, and partial results through analysis.

**Processing identity (ADR-105).**  The same request -- structure, method,
chain, selection -- can produce different numbers when the rules that admit
atoms change.  ADR-101 made the amino-acid atomic methods read TPO, SEP and
PTR written as ``HETATM``; the key space and the version number stayed the
same.  ``tests/data/legacy_outputs/before_adr101_hetatm_modified_residues/`` is a real output of the
older implementation (its README gives the commit and command): a selection
of two ordinary residues and TPO 3, computed without TPO.  These tests check
that ``--if-exists skip`` never keeps it as if it were current, that analysis
never compares it with a current result without saying so, and that its key
space is still recognised as the same.

**Partial results (ADR-106).**  A result whose selection did not get
everything it asked for is kept and usable, but it is never presented as the
complete selection: the structured selection record says what is missing,
the batch says ``partial`` (also when a resumed run keeps it), the exit code
is 4, and every analysis refuses it unless ``allow_partial=True``
(``--allow-partial``) -- a separate acknowledgement from mixed processing,
which relaxes nothing about key definitions.
"""

from __future__ import annotations

import contextlib
import csv
import hashlib
import io
import json
import os
import shutil

import pytest

import biotsr
from biotsr.analysis import (
    AnalysisSample,
    IncompatibleProcessing,
    PartialSamples,
    common_key_intersection,
    feature_matrix,
)
from biotsr.cli import EXIT_INCOMPLETE, EXIT_OK, EXIT_PARTIAL, EXIT_USAGE, main
from biotsr.cli import output as cli_output
from biotsr.cli.errors import EXIT_CODES
from biotsr.domain.keyspace import IncompatibleKeySpace
from biotsr.io.native_read import load_result
from biotsr.api import processing_identity
from biotsr.selection.base import ModifiedResiduePolicy

from ..conftest import SYNTHETIC

HERE = os.path.dirname(os.path.abspath(__file__))
SYN = os.path.join(SYNTHETIC, "SYNPHOS.pdb")
LEGACY = os.path.join(HERE, "..", "data", "legacy_outputs", "before_adr101_hetatm_modified_residues")
PREFIX = "SYNPHOS_A_specific"
#: The request the older output answered, spelled exactly as it was run.
REQUEST = ["--method", "aa_ind_specific", "--chain", "A",
           "--residues", "SER_4;TPO_3;TYR_6", "--sample-id", PREFIX,
           "--no-triplets"]


def run(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main([str(a) for a in argv])
    return code, out.getvalue(), err.getvalue()


def summary(folder):
    with open(os.path.join(str(folder), "batch_summary.csv"), newline="") as fh:
        return list(csv.DictReader(fh))


def write(path, text):
    with open(str(path), "w", newline="") as fh:
        fh.write(text)
    return str(path)


def digest(folder):
    return {name: hashlib.sha256(open(os.path.join(folder, name), "rb").read()).hexdigest()
            for name in sorted(os.listdir(folder))}


@pytest.fixture
def legacy(tmp_path):
    """The older output in a folder of its own, beside a copy of its structure."""
    shutil.copy(SYN, str(tmp_path))
    out = tmp_path / "out"
    shutil.copytree(LEGACY, str(out), ignore=shutil.ignore_patterns("README.md"))
    return tmp_path


def sample(path):
    return AnalysisSample.from_result(load_result(str(path)))


# ===========================================================================
# The fixture is what its README says
# ===========================================================================


def test_the_older_output_is_unedited_and_was_computed_without_tpo():
    with open(os.path.join(LEGACY, "README.md")) as fh:
        readme = fh.read()
    for name, sha in digest(LEGACY).items():
        if name != "README.md":
            assert sha in readme, name
    with open(os.path.join(LEGACY, PREFIX + ".run.json")) as fh:
        record = json.load(fh)
    # The same version number as the build the processing rules were added in
    # (1.0.0rc4), which is why the version alone could not tell them apart.
    assert record["biotsr_version"] == "1.0.0rc4"
    assert "processing" not in record["request"]
    [entry] = record["results"]
    assert entry["counts"]["n_points"] == 18                 # SER 4 + TYR 6 only
    assert entry["warnings"] == ["requested but not found: TPO_3"]
    result = load_result(os.path.join(LEGACY, PREFIX + ".keys"))
    assert result.source.modified_residue_policy is None     # not recorded


def test_the_older_request_differs_from_the_current_one_only_in_its_processing_rules(legacy):
    """Without the processing rules, the old record would match exactly."""
    with open(str(legacy / "out" / (PREFIX + ".run.json"))) as fh:
        recorded = json.load(fh)["request"]
    current = biotsr.analyze(SYN, method="aa_ind_specific", chain="A",
                             residues="SER_4;TPO_3;TYR_6", keep_triplets=False)
    assert current.counts.n_points == 29                     # TPO 3's 11 atoms too
    from biotsr.cli.spec import build_spec

    spec = build_spec({"method": "aa_ind_specific", "pdb": SYN, "chain": "A",
                       "residues": "SER_4;TPO_3;TYR_6"})
    request = cli_output.request_record(spec, recorded["structure_sha256"], "SYNPHOS")
    differing = {k for k in set(recorded) | set(request) if recorded.get(k) != request.get(k)}
    assert differing == {"processing"}
    assert request["processing"] == processing_identity("aa_ind_specific")


# ===========================================================================
# Resume never keeps an output computed under other rules
# ===========================================================================


def test_run_with_skip_refuses_the_older_output_and_leaves_it_untouched(legacy):
    before = digest(str(legacy / "out"))
    code, _, err = run(["run", "--pdb", legacy / "SYNPHOS.pdb", "--output",
                        legacy / "out", "--if-exists", "skip"] + REQUEST)
    assert code == EXIT_USAGE, err
    assert "cannot be reused" in err and "processing rules" in err
    assert "modified_residue_policy=%s" % ModifiedResiduePolicy.ATOMIC in err
    assert "--if-exists overwrite" in err
    assert digest(str(legacy / "out")) == before


def test_batch_with_skip_reports_the_older_output_as_not_reusable(legacy):
    manifest = write(legacy / "m.csv", "pdb_id,chain,method,residues,sample_id\n"
                     "SYNPHOS,A,aa_ind_specific,SER_4;TPO_3;TYR_6,%s\n" % PREFIX)
    code, _, err = run(["batch", "--manifest", manifest, "--pdb-dir", legacy,
                        "--output", legacy / "out", "--if-exists", "skip",
                        "--no-triplets"])
    assert code == EXIT_PARTIAL, err
    [row] = summary(legacy / "out")
    assert row["status"] == "failed"
    assert "processing rules" in row["message"]


def test_overwrite_recomputes_under_the_current_rules_and_then_skip_keeps_it(legacy):
    args = ["run", "--pdb", legacy / "SYNPHOS.pdb", "--output", legacy / "out"] + REQUEST
    code, _, err = run(args + ["--if-exists", "overwrite"])
    assert code == EXIT_OK, err
    with open(str(legacy / "out" / (PREFIX + ".run.json"))) as fh:
        record = json.load(fh)
    assert record["request"]["processing"] == processing_identity("aa_ind_specific")
    assert record["results"][0]["counts"]["n_points"] == 29
    assert record["results"][0]["omissions"] == {}
    result = load_result(str(legacy / "out" / (PREFIX + ".keys")))
    assert result.source.modified_residue_policy == ModifiedResiduePolicy.ATOMIC
    code, out, err = run(args + ["--if-exists", "skip"])
    assert code == EXIT_OK, err
    assert "skipped" in out + err


def test_a_record_of_different_processing_rules_is_not_reused(tmp_path):
    args = ["run", "--pdb", SYN, "--output", tmp_path / "out"] + REQUEST
    assert run(args)[0] == EXIT_OK
    path = tmp_path / "out" / (PREFIX + ".run.json")
    with open(str(path)) as fh:
        record = json.load(fh)
    record["request"]["processing"]["modified_residue_policy"] = "hetatm_amino_acids:"
    with open(str(path), "w") as fh:
        json.dump(record, fh)
    code, _, err = run(args + ["--if-exists", "skip"])
    assert code == EXIT_USAGE
    assert "different processing rules" in err


# ===========================================================================
# Analysis: same key space, different processing -- and the old one is partial
# ===========================================================================


def test_analysis_keeps_the_key_space_but_names_the_processing_difference(legacy, tmp_path):
    old = sample(legacy / "out" / (PREFIX + ".keys"))
    code, _, err = run(["run", "--pdb", SYN, "--output", tmp_path / "new"] + REQUEST)
    assert code == EXIT_OK, err
    new = sample(tmp_path / "new" / (PREFIX + ".keys"))
    as_written = AnalysisSample(sample_id="older", counts=old.counts,
                                key_space=old.key_space, sampling=old.sampling,
                                annotations=old.annotations)
    # Since 1.0.0rc7 atoms are labelled by element (ADR-113): the older output's
    # key space is another one, and no opt-in compares it with a current result.
    assert as_written.key_space.fingerprint != new.key_space.fingerprint
    for kwargs in ({}, {"allow_mixed_processing": True, "allow_partial": True}):
        with pytest.raises(IncompatibleKeySpace):
            common_key_intersection([as_written, new], **kwargs)
    # The processing and partial checks below are the next two steps of the
    # rule; to reach them, the older counts are given the current key space
    # (the labels of SER, TYR and their atoms are the same codes in both).
    old = AnalysisSample(sample_id="older", counts=old.counts, key_space=new.key_space,
                         sampling=old.sampling, annotations=old.annotations)
    # The processing rules differ, and an unrecorded rule is never read as current.
    with pytest.raises(IncompatibleProcessing) as exc:
        common_key_intersection([old, new])
    assert "modified_residue_policy (not recorded vs '%s')" % ModifiedResiduePolicy.ATOMIC \
        in str(exc.value)
    # Allowing mixed processing does not acknowledge that the old one is partial.
    with pytest.raises(PartialSamples) as exc:
        common_key_intersection([old, new], allow_mixed_processing=True)
    assert "older" in str(exc.value) and "TPO_3" in str(exc.value)
    both = common_key_intersection([old, new], allow_mixed_processing=True,
                                   allow_partial=True)
    # The input rules (ADR-108) are not recorded by the older output either.
    assert both.comparison["processing_differences"] == [
        "input_rules", "modified_residue_policy"]
    assert both.comparison["partial_samples"] == {
        "older": {"messages": ["requested but not found: TPO_3"]}}


# ===========================================================================
# Partial results: same residues requested, one structure lacks one
# ===========================================================================


def _mix_manifest(folder):
    return write(folder / "m.csv", "pdb_id,chain,method,residue_types,sample_id\n"
                 "SYNPHOS,A,aa_ind_mix,TYR;LYS;SER,chainA\n"
                 "SYNPHOS,B,aa_ind_mix,TYR;LYS;SER,chainB\n")


@pytest.fixture
def mixed(tmp_path):
    shutil.copy(SYN, str(tmp_path))
    manifest = _mix_manifest(tmp_path)
    code, out, err = run(["batch", "--manifest", manifest, "--pdb-dir", tmp_path,
                          "--output", tmp_path / "out", "--no-triplets"])
    return tmp_path, manifest, code, err


def test_the_batch_marks_the_sample_lacking_a_residue_partial_and_exits_4(mixed):
    folder, _, code, err = mixed
    assert code == EXIT_INCOMPLETE, err
    rows = {r["sample_id"]: r for r in summary(folder / "out")}
    assert rows["chainA"]["status"] == "ok"
    assert rows["chainB"]["status"] == "partial"
    assert "requested but not found: residue type SER" in rows["chainB"]["message"]
    assert "1 partial" in err
    assert EXIT_CODES[EXIT_INCOMPLETE] == "completed, but some results are partial"


def test_the_omission_is_structured_and_survives_saving_and_loading(mixed):
    folder = mixed[0]
    with open(str(folder / "out" / "chainB.run.json")) as fh:
        entry = json.load(fh)["results"][0]
    assert entry["omissions"]["primary"] == {"not_found": ["SER"]}
    loaded = sample(folder / "out" / "chainB.keys")
    assert loaded.omissions["primary"] == {"not_found": ["SER"]}
    assert loaded.is_partial
    assert not sample(folder / "out" / "chainA.keys").is_partial


def test_a_resumed_batch_keeps_the_partial_result_and_still_says_so(mixed):
    folder, manifest, _, _ = mixed
    code, _, err = run(["batch", "--manifest", manifest, "--pdb-dir", folder,
                        "--output", folder / "out", "--no-triplets",
                        "--if-exists", "skip"])
    assert code == EXIT_INCOMPLETE, err
    rows = {r["sample_id"]: r for r in summary(folder / "out")}
    assert rows["chainA"]["status"] == rows["chainB"]["status"] == "skipped"
    assert rows["chainA"]["message"] == "verified existing output of the same request"
    assert "kept partial result: requested but not found: residue type SER" \
        in rows["chainB"]["message"]
    assert "2 skipped (1 of them partial)" in err


def test_run_exits_4_for_a_partial_result_and_for_keeping_one(tmp_path):
    args = ["run", "--pdb", SYN, "--method", "aa_ind_mix", "--chain", "B",
            "--residue-types", "TYR;LYS;SER", "--output", tmp_path, "--no-triplets"]
    code, _, err = run(args)
    assert code == EXIT_INCOMPLETE, err
    code, _, err = run(args + ["--if-exists", "skip"])
    assert code == EXIT_INCOMPLETE
    assert "kept partial result" in err


def test_analysis_refuses_a_partial_sample_until_it_is_acknowledged(mixed):
    folder = mixed[0]
    a, b = (sample(folder / "out" / n) for n in ("chainA.keys", "chainB.keys"))
    a = AnalysisSample(sample_id="chainA", counts=a.counts, key_space=a.key_space,
                       sampling=a.sampling, annotations=a.annotations)
    b = AnalysisSample(sample_id="chainB", counts=b.counts, key_space=b.key_space,
                       sampling=b.sampling, annotations=b.annotations)
    for build in (common_key_intersection, feature_matrix):
        with pytest.raises(PartialSamples) as exc:
            build([a, b])
        assert exc.value.partial == {"chainB": {"primary": {"not_found": ["SER"]},
                                                "messages": [
            "requested but not found: residue type SER"]}}
        assert "chainB (not found: SER)" in str(exc.value)
    matrix = feature_matrix([a, b], allow_partial=True)
    assert matrix.n_samples == 2                             # nothing dropped
    assert matrix.comparison["partial_allowed"] is True
    assert list(matrix.comparison["partial_samples"]) == ["chainB"]
    assert any("partial and were explicitly allowed" in w for w in matrix.warnings)
    rows = {p["sample_id"]: p for p in matrix.samples}
    assert rows["chainA"]["partial"] is False and rows["chainA"]["omissions"] == {}
    assert rows["chainB"]["omissions"]["primary"] == {"not_found": ["SER"]}
    from biotsr.analysis import similarity_matrix

    similarity = similarity_matrix(matrix, "generalized_jaccard")
    exported = similarity.provenance()
    assert exported["comparison"]["partial_samples"]["chainB"]["primary"] == {
        "not_found": ["SER"]}
    assert any("partial" in w for w in exported["warnings"])
    from biotsr.analysis import hierarchy

    tree = hierarchy(similarity.to_distance()).provenance()
    assert tree["comparison"]["partial_samples"]["chainB"]["primary"] == {
        "not_found": ["SER"]}
    assert {s["sample_id"]: s["partial"] for s in tree["samples"]} == {
        "chainA": False, "chainB": True}


def test_acknowledging_partial_samples_never_relaxes_key_definitions(mixed):
    folder = mixed[0]
    b = sample(folder / "out" / "chainB.keys")
    calpha = AnalysisSample.from_result(
        biotsr.analyze(SYN, method="protein", chain="A", keep_triplets=False),
        sample_id="calpha")
    with pytest.raises(IncompatibleKeySpace):
        feature_matrix([b, calpha], allow_partial=True, allow_mixed_processing=True)


def test_complete_lys_only_and_tyr_only_samples_compare_without_any_acknowledgement():
    lys = AnalysisSample.from_result(biotsr.analyze(
        SYN, method="aa_ind", chain="A", residue_types="LYS", keep_triplets=False),
        sample_id="lys")
    tyr = AnalysisSample.from_result(biotsr.analyze(
        SYN, method="aa_ind", chain="B", residue_types="TYR", keep_triplets=False),
        sample_id="tyr")
    matrix = feature_matrix([lys, tyr])
    assert matrix.comparison["partial_samples"] == {}
    assert matrix.comparison["processing_uniform"] is True
    assert matrix.warnings == ()


def test_the_cli_analysis_needs_allow_partial_and_then_exits_4(mixed):
    folder = mixed[0]
    code, _, err = run(["analyze", "common-keys", "--from-results", folder / "out",
                        "--output", folder / "a1"])
    assert code == EXIT_USAGE
    assert "chainB" in err and "--allow-partial" in err
    code, _, err = run(["analyze", "common-keys", "--from-results", folder / "out",
                        "--output", folder / "a2", "--allow-partial"])
    assert code == EXIT_INCOMPLETE, err
    with open(str(folder / "a2" / "analysis" / "common_keys" / "all_samples.json")) as fh:
        report = json.load(fh)
    assert report["comparison"]["partial_allowed"] is True
    assert report["comparison"]["partial_samples"]["chainB"]["primary"] == {
        "not_found": ["SER"]}
    by_id = {s["sample_id"]: s for s in report["samples"]}
    assert by_id["chainB"]["partial"] is True and by_id["chainA"]["partial"] is False


# ===========================================================================
# Each kind of omission is structured, not only worded
# ===========================================================================


def test_unsupported_and_completeness_rejected_residues_are_structured(tmp_path):
    # A copy of the synthetic structure without LYS 10's NZ: the exact-count
    # rule then rejects LYS 10.
    trimmed = tmp_path / "SYNPHOS.pdb"
    with open(SYN) as src, open(str(trimmed), "w") as dst:
        for line in src:
            if not (line.startswith("ATOM") and line[17:20] == "LYS"
                    and line[21] == "A" and line[12:16].strip() == "NZ"):
                dst.write(line)
    result = biotsr.analyze(str(trimmed), method="aa_ind_specific", chain="A",
                            residues="TYR_6;LYS_10;CME_8;VAL_99", keep_triplets=False)
    found = AnalysisSample.from_result(result).omissions["primary"]
    assert found == {"not_found": ["VAL_99"], "unsupported": ["CME_8"],
                     "declined": ["LYS_10"]}


def test_a_calpha_whole_chain_with_an_unlabelled_modified_residue_is_partial(tmp_path):
    """The Cα methods read the supported modified residues (ADR-107); one the
    residue labels do not cover -- CME 8 renamed to HYP here -- is named and
    makes the result partial."""
    renamed = tmp_path / "SYNPHOS.pdb"
    with open(SYN) as src, open(str(renamed), "w") as dst:
        for line in src:
            if line.startswith("HETATM") and line[17:20] == "CME":
                line = line[:17] + "HYP" + line[20:]
            dst.write(line)
    result = biotsr.analyze(str(renamed), method="protein", chain="A", keep_triplets=False)
    om = AnalysisSample.from_result(result).omissions["primary"]
    assert om == {"unsupported": ["HYP_8"]}
    complete = biotsr.analyze(SYN, method="protein", chain="A", keep_triplets=False)
    assert AnalysisSample.from_result(complete).omissions == {}


def test_a_selection_too_small_for_a_triangle_that_also_lacks_atoms_exits_4(tmp_path):
    """``no_triangles`` alone is a valid outcome (exit 0); with a named atom
    missing it is also incomplete, and the exit code says so."""
    shutil.copy(SYN, str(tmp_path))
    manifest = write(tmp_path / "m.csv",
                     "pdb_id,chain,method,entity_name,entity_id,atoms,sample_id\n"
                     "SYNPHOS,A,drug_atom_selection,LIG,101,C1;N1;XX9,lacking\n"
                     "SYNPHOS,A,drug_atom_selection,LIG,101,C1;N1,small\n")
    code, _, err = run(["batch", "--manifest", manifest, "--pdb-dir", tmp_path,
                        "--output", tmp_path / "out", "--no-triplets"])
    assert code == EXIT_INCOMPLETE, err
    rows = {r["sample_id"]: r for r in summary(tmp_path / "out")}
    assert rows["lacking"]["status"] == rows["small"]["status"] == "no_triangles"
    assert "XX9" in rows["lacking"]["message"]
    assert "2 with no triangles (1 of them partial)" in err
    small = write(tmp_path / "s.csv", "pdb_id,chain,method,entity_name,entity_id,atoms,sample_id\n"
                  "SYNPHOS,A,drug_atom_selection,LIG,101,C1;N1,small\n")
    code, _, err = run(["batch", "--manifest", small, "--pdb-dir", tmp_path,
                        "--output", tmp_path / "out2", "--no-triplets"])
    assert code == EXIT_OK, err
