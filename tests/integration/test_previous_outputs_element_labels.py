"""Outputs of the previous build, after ADR-113 and ADR-114.

``tests/data/legacy_outputs/before_adr113_adr114`` holds three real
calculations written by 1.0.0rc6 (see its README):

* an atomic one (``drug``) -- the same keys today, but in another key space,
  so it is never reused by ``--if-exists skip`` and never compared with a
  current result, whatever opt-in is given;
* a ``nucleotide_chain`` one -- computed under the old count rule (DC 1
  excluded), so its processing rules differ as well as its key space;
* a Calpha one -- nothing about it changed, so it is reused.

Refusing an output never touches it.
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
from biotsr.analysis import AnalysisSample, common_key_intersection, feature_matrix
from biotsr.cli import EXIT_OK, EXIT_USAGE, main
from biotsr.domain.keyspace import IncompatibleKeySpace
from biotsr.io.native_read import load_result

from ..conftest import STRUCTURES, SYNTHETIC

HERE = os.path.dirname(os.path.abspath(__file__))
PREVIOUS = os.path.join(HERE, "..", "data", "legacy_outputs", "before_adr113_adr114")
SYN = os.path.join(SYNTHETIC, "SYNPHOS.pdb")
EXCERPT_SHA256 = "f0888b7f91215e4317b724632cc98b8c049fc3e865377bfbf88e5816f49426de"

DRUG = ["--method", "drug", "--chain", "A", "--entity-name", "LIG", "--entity-id", "101",
        "--sample-id", "SYNPHOS_A_LIG101", "--no-triplets"]
CHAIN = ["--method", "nucleotide_chain", "--chain", "A",
         "--sample-id", "1BNA_A_1-3_chain", "--no-triplets"]
PROTEIN = ["--method", "protein", "--chain", "A", "--sample-id", "SYNPHOS_A_protein",
           "--no-triplets"]


def run(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main([str(a) for a in argv])
    return code, out.getvalue(), err.getvalue()


def digest(folder):
    return {name: hashlib.sha256(open(os.path.join(str(folder), name), "rb").read()).hexdigest()
            for name in sorted(os.listdir(str(folder)))}


def write_excerpt(folder):
    """Chain A residues 1-3 of 1BNA and a TER: the structure the old output read."""
    with open(os.path.join(STRUCTURES, "1BNA.pdb")) as fh:
        lines = [l for l in fh if l.startswith("ATOM  ") and l[21] == "A"
                 and l[22:26].strip() in ("1", "2", "3")]
    lines += ["%-6s%5d      %3s %1s%4d\n" % ("TER", len(lines) + 1, " DC", "A", 3), "END\n"]
    path = os.path.join(str(folder), "1BNA_A_1-3.pdb")
    with open(path, "w") as fh:
        fh.writelines(lines)
    return path


@pytest.fixture
def previous(tmp_path):
    shutil.copy(SYN, str(tmp_path))
    write_excerpt(tmp_path)
    shutil.copytree(PREVIOUS, str(tmp_path / "out"), ignore=shutil.ignore_patterns("README.md"))
    return tmp_path


def _run(previous, pdb, request, if_exists="skip"):
    return run(["run", "--pdb", previous / pdb, "--output", previous / "out",
                "--if-exists", if_exists] + request)


def test_the_previous_outputs_are_unedited_and_were_written_by_rc6(tmp_path):
    with open(os.path.join(PREVIOUS, "README.md")) as fh:
        readme = fh.read()
    for name, sha in digest(PREVIOUS).items():
        if name != "README.md":
            assert sha in readme, name
    with open(write_excerpt(tmp_path), "rb") as fh:
        assert hashlib.sha256(fh.read()).hexdigest() == EXCERPT_SHA256
    for prefix, points, space in (("SYNPHOS_A_LIG101", 8, "1187735d"),
                                  ("1BNA_A_1-3_chain", 41, "1187735d"),
                                  ("SYNPHOS_A_protein", 14, "da39d735")):
        with open(os.path.join(PREVIOUS, prefix + ".run.json")) as fh:
            record = json.load(fh)
        assert record["biotsr_version"] == "1.0.0rc6"
        assert "completeness_policy" not in record["request"]["processing"]
        [entry] = record["results"]
        assert entry["counts"]["n_points"] == points
        assert entry["key_space_fingerprint"].startswith(space)
    old = load_result(os.path.join(PREVIOUS, "SYNPHOS_A_LIG101.keys"))
    assert old.key_space.label_scheme.name == "atom_lexicon_standard_112"


def test_the_atomic_output_has_todays_keys_in_another_key_space(previous):
    old = load_result(str(previous / "out" / "SYNPHOS_A_LIG101.keys"))
    new = biotsr.analyze(SYN, method="drug", chain="A", entity_name="LIG", entity_id="101",
                         keep_triplets=False)
    assert dict(old.frequencies) == dict(new.frequencies)
    assert old.key_space.fingerprint != new.key_space.fingerprint
    assert new.key_space.label_scheme.name == "element_labels_112"


@pytest.mark.parametrize("pdb,request_args,words", [
    ("SYNPHOS.pdb", DRUG, ["different key space"]),
    ("1BNA_A_1-3.pdb", CHAIN, ["different key space",
                               "completeness_policy=available_nucleotide_atoms"]),
], ids=["drug", "nucleotide_chain"])
def test_skip_refuses_an_old_atomic_output_and_leaves_every_file_untouched(
        previous, pdb, request_args, words):
    before = digest(previous / "out")
    code, _, err = _run(previous, pdb, request_args)
    assert code == EXIT_USAGE, err
    assert "cannot be reused" in err
    for word in words:
        assert word in err, err
    assert digest(previous / "out") == before


def test_skip_reuses_the_unchanged_calpha_output(previous):
    before = digest(previous / "out")
    code, out, err = _run(previous, "SYNPHOS.pdb", PROTEIN)
    assert code == EXIT_OK, err
    assert digest(previous / "out") == before


def test_a_batch_resume_keeps_the_calpha_output_and_refuses_the_atomic_ones(previous):
    manifest = previous / "m.csv"
    manifest.write_text(
        "pdb_id,chain,method,entity_name,entity_id,sample_id\n"
        "SYNPHOS,A,drug,LIG,101,SYNPHOS_A_LIG101\n"
        "1BNA_A_1-3,A,nucleotide_chain,,,1BNA_A_1-3_chain\n"
        "SYNPHOS,A,protein,,,SYNPHOS_A_protein\n")
    before = digest(previous / "out")
    code, _, err = run(["batch", "--manifest", manifest, "--pdb-dir", previous,
                        "--output", previous / "out", "--if-exists", "skip",
                        "--no-triplets"])
    with open(str(previous / "out" / "batch_summary.csv"), newline="") as fh:
        rows = {r["sample_id"]: r for r in csv.DictReader(fh)}
    assert rows["SYNPHOS_A_protein"]["status"] == "skipped"
    for prefix in ("SYNPHOS_A_LIG101", "1BNA_A_1-3_chain"):
        assert rows[prefix]["status"] == "failed"
        assert "different key space" in rows[prefix]["message"]
    after = digest(previous / "out")
    assert {k: v for k, v in after.items() if k in before} == before
    assert code != EXIT_OK


def test_overwrite_recomputes_the_chain_under_the_current_rules(previous):
    code, _, err = _run(previous, "1BNA_A_1-3.pdb", CHAIN, if_exists="overwrite")
    assert code == EXIT_OK, err
    with open(str(previous / "out" / "1BNA_A_1-3_chain.run.json")) as fh:
        record = json.load(fh)
    [entry] = record["results"]
    assert entry["counts"]["n_points"] == 57
    assert entry["omissions"] == {}
    assert record["request"]["processing"]["completeness_policy"] == \
        "available_nucleotide_atoms"
    code, _, err = _run(previous, "1BNA_A_1-3.pdb", CHAIN)
    assert code == EXIT_OK, err


def test_no_opt_in_compares_an_old_atomic_output_with_a_current_one(previous):
    old = AnalysisSample.from_result(load_result(str(previous / "out" / "SYNPHOS_A_LIG101.keys")))
    old = AnalysisSample(sample_id="previous", counts=old.counts, key_space=old.key_space,
                         sampling=old.sampling, annotations=old.annotations)
    new = AnalysisSample.from_result(biotsr.analyze(SYN, method="drug", chain="A",
                                                    entity_name="LIG", entity_id="101"))
    for kwargs in ({}, {"allow_mixed_processing": True},
                   {"allow_mixed_processing": True, "allow_partial": True}):
        with pytest.raises(IncompatibleKeySpace):
            common_key_intersection([old, new], **kwargs)
        with pytest.raises(IncompatibleKeySpace):
            feature_matrix([old, new], **kwargs)


def test_the_analysis_command_refuses_them_too(previous, tmp_path):
    """An old and a current drug result in one folder: no opt-in compares them."""
    code, _, err = run(["run", "--pdb", previous / "SYNPHOS.pdb", "--output",
                        previous / "out", "--method", "drug", "--chain", "A",
                        "--entity-name", "LIG", "--entity-id", "101", "--sample-id", "now",
                        "--no-triplets"])
    assert code == EXIT_OK, err
    code, _, err = run(["analyze", "common-keys", "--from-results", previous / "out",
                        "--method", "drug", "--allow-mixed-processing", "--allow-partial",
                        "--output", tmp_path / "analysis"])
    assert code != EXIT_OK
    assert "key space" in err, err
