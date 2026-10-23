"""End-to-end command behaviour: batch, failure policy, analysis, describe.

These run real commands over real pinned fixtures.  What they assert is
orchestration -- that every input row is accounted for, that a failure stays
attached to the row that caused it, that provenance survives the round trip to
disk -- not science, which ``test_cli_equivalence.py`` owns.
"""

from __future__ import annotations

import contextlib
import csv
import io
import json
import os

import pytest

import biotsr
from biotsr.cli import EXIT_OK, EXIT_PARTIAL, EXIT_USAGE, main

from ..conftest import _require


def cli(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main(argv)
    return code, out.getvalue(), err.getvalue()


@pytest.fixture(scope="module")
def pdb_dir():
    tests_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    workspace = os.path.dirname(tests_dir)             # the repository root
    root = os.path.join(workspace, "tests", "data", "structures")
    _require(root)
    return root


def manifest(tmp_path, text, name="m.csv"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


GOOD = (
    "sample_id,pdb_id,method,chain,label\n"
    "ubq,1UBQ,protein,A,ubiquitin\n"
    "brs_a,1BRS,protein,A,barnase\n"
    "brs_b,1BRS,protein,B,barnase\n"
)

MIXED = GOOD + (
    "no_such_chain,1UBQ,protein,Q,ubiquitin\n"
    "no_such_method,1UBQ,wobble,A,ubiquitin\n"
    "no_such_file,9XYZ,protein,A,ubiquitin\n"
)


# ---------------------------------------------------------------------------
# Batch
# ---------------------------------------------------------------------------


def test_a_clean_batch_writes_one_calculation_per_row_into_one_folder(tmp_path, pdb_dir):
    out = tmp_path / "out"
    code, stdout, _err = cli(
        ["batch", "-M", manifest(tmp_path, GOOD), "--pdb-dir", pdb_dir,
         "-o", str(out), "--no-triplets", "-q"]
    )
    assert code == EXIT_OK
    assert stdout.strip() == str(out)
    for name in ("ubq", "brs_a", "brs_b"):
        for suffix in (".keys", ".meta.json", ".run.json"):
            assert os.path.isfile(os.path.join(str(out), name + suffix))
    assert not [n for n in os.listdir(str(out))
                if os.path.isdir(os.path.join(str(out), n))]
    assert os.path.isfile(os.path.join(str(out), "batch_summary.csv"))
    assert not os.path.exists(os.path.join(str(out), "batch_failures.csv"))


def test_every_input_row_appears_in_the_summary(tmp_path, pdb_dir):
    out = tmp_path / "out"
    code, _stdout, _err = cli(
        ["batch", "-M", manifest(tmp_path, MIXED), "--pdb-dir", pdb_dir,
         "-o", str(out), "--no-triplets", "-q"]
    )
    assert code == EXIT_PARTIAL
    rows = list(csv.DictReader(open(os.path.join(str(out), "batch_summary.csv"))))
    assert len(rows) == 6
    assert [r["rows"] for r in rows] == ["2", "3", "4", "5", "6", "7"]
    assert [r["status"] for r in rows] == (
        ["ok"] * 3 + ["failed"] * 3
    )


def test_every_failure_stays_attributable_to_its_own_row(tmp_path, pdb_dir):
    out = tmp_path / "out"
    cli(
        ["batch", "-M", manifest(tmp_path, MIXED), "--pdb-dir", pdb_dir,
         "-o", str(out), "--no-triplets", "-q"]
    )
    failures = list(
        csv.DictReader(open(os.path.join(str(out), "batch_failures.csv")))
    )
    assert len(failures) == 3
    by_row = {r["rows"]: r for r in failures}
    assert "chain 'Q'" in by_row["5"]["message"]
    assert "wobble" in by_row["6"]["message"]
    assert "9XYZ" in by_row["7"]["message"]
    for record in failures:
        assert record["message"], "a failure with no reason is a silent omission"


def test_a_failed_sample_leaves_no_output_directory(tmp_path, pdb_dir):
    out = tmp_path / "out"
    cli(
        ["batch", "-M", manifest(tmp_path, MIXED), "--pdb-dir", pdb_dir,
         "-o", str(out), "--no-triplets", "-q"]
    )
    assert not [n for n in os.listdir(str(out)) if n.startswith("no_such_chain")]
    assert not any(n.startswith(".biotsr-") for n in os.listdir(str(out)))


def test_stop_halts_at_the_first_failure_and_still_writes_the_report(
    tmp_path, pdb_dir
):
    out = tmp_path / "out"
    code, _stdout, _err = cli(
        ["batch", "-M", manifest(tmp_path, MIXED), "--pdb-dir", pdb_dir,
         "-o", str(out), "--no-triplets", "-q", "--on-error", "stop"]
    )
    assert code == EXIT_USAGE
    rows = list(csv.DictReader(open(os.path.join(str(out), "batch_summary.csv"))))
    # The unknown method and the missing structure are visible before anything
    # is computed, so `stop` halts there: nothing runs, and every row is listed.
    assert [r["status"] for r in rows] == [
        "not_started", "not_started", "not_started", "not_started", "failed", "failed"]
    assert os.path.isfile(os.path.join(str(out), "batch_failures.csv"))
    assert sorted(os.listdir(str(out))) == ["batch_failures.csv", "batch_summary.csv"]


def test_stop_halts_at_the_first_failure_found_while_computing(tmp_path, pdb_dir):
    out = tmp_path / "out"
    code, _stdout, _err = cli(
        ["batch", "-M", manifest(tmp_path, GOOD + "no_such_chain,1UBQ,protein,Q,x\n"
                                 "after,1BRS,protein,D,x\n"),
         "--pdb-dir", pdb_dir, "-o", str(out), "--no-triplets", "-q",
         "--on-error", "stop"]
    )
    # a row that fails while computing leaves the batch partial (exit 3)
    assert code == EXIT_PARTIAL
    rows = list(csv.DictReader(open(os.path.join(str(out), "batch_summary.csv"))))
    assert [r["status"] for r in rows] == ["ok", "ok", "ok", "failed", "not_started"]


def test_skip_resumes_an_interrupted_batch(tmp_path, pdb_dir):
    out = tmp_path / "out"
    cli(
        ["batch", "-M", manifest(tmp_path, GOOD), "--pdb-dir", pdb_dir,
         "-o", str(out), "--no-triplets", "-q"]
    )
    before = os.path.getmtime(os.path.join(str(out), "ubq.run.json"))
    code, _stdout, _err = cli(
        ["batch", "-M", manifest(tmp_path, GOOD), "--pdb-dir", pdb_dir,
         "-o", str(out), "--no-triplets", "-q", "--if-exists", "skip"]
    )
    assert code == EXIT_OK
    assert os.path.getmtime(os.path.join(str(out), "ubq.run.json")) == before
    rows = list(csv.DictReader(open(os.path.join(str(out), "batch_summary.csv"))))
    assert [r["status"] for r in rows] == ["skipped"] * 3


def test_a_second_batch_into_the_same_directory_is_refused_by_default(
    tmp_path, pdb_dir
):
    out = tmp_path / "out"
    cli(
        ["batch", "-M", manifest(tmp_path, GOOD), "--pdb-dir", pdb_dir,
         "-o", str(out), "--no-triplets", "-q"]
    )
    code, _stdout, _err = cli(
        ["batch", "-M", manifest(tmp_path, GOOD), "--pdb-dir", pdb_dir,
         "-o", str(out), "--no-triplets", "-q"]
    )
    assert code == EXIT_PARTIAL
    failures = list(
        csv.DictReader(open(os.path.join(str(out), "batch_failures.csv")))
    )
    assert len(failures) == 3
    assert all("already exists" in r["message"] for r in failures)


def test_command_line_values_are_row_defaults(tmp_path, pdb_dir):
    out = tmp_path / "out"
    text = "sample_id,pdb_id,chain\nubq,1UBQ,A\n"
    code, _stdout, _err = cli(
        ["batch", "-M", manifest(tmp_path, text), "--pdb-dir", pdb_dir,
         "-o", str(out), "--no-triplets", "-q", "-m", "protein", "--mirror"]
    )
    assert code == EXIT_OK
    payload = json.load(open(os.path.join(str(out), "ubq.run.json")))
    assert payload["method"] == "protein"
    assert payload["parameters"]["mirror"] is True


def test_a_row_method_and_a_command_line_method_conflict(tmp_path, pdb_dir):
    """Since 1.0.0rc4 --method chooses the methods for every row; a manifest
    that also names methods is refused before anything runs, naming both
    remedies, rather than one silently winning."""
    out = tmp_path / "out"
    text = (
        "sample_id,pdb_id,chain,method,residues\n"
        "ubq,1UBQ,A,aa_selection,10-30\n"
    )
    code, _stdout, err = cli(
        ["batch", "-M", manifest(tmp_path, text), "--pdb-dir", pdb_dir,
         "-o", str(out), "--no-triplets", "-q", "-m", "protein"]
    )
    assert code == EXIT_USAGE
    assert "--method also chooses methods" in err and "--ignore-column method" in err
    assert not out.exists()
    code, _stdout, _err = cli(
        ["batch", "-M", manifest(tmp_path, text), "--pdb-dir", pdb_dir,
         "-o", str(out), "--no-triplets", "-q"]
    )
    assert code == EXIT_OK
    assert json.load(open(os.path.join(str(out), "ubq.run.json")))["method"] == "aa_selection"


def test_a_batch_row_and_a_run_command_produce_the_same_files(tmp_path, pdb_dir):
    text = "sample_id,pdb_id,method,chain,sse\nubq,1UBQ,protein,A,true\n"
    batch_out = tmp_path / "batch"
    run_out = tmp_path / "run"
    assert cli(
        ["batch", "-M", manifest(tmp_path, text), "--pdb-dir", pdb_dir,
         "-o", str(batch_out), "--no-triplets", "-q"]
    )[0] == EXIT_OK
    assert cli(
        ["run", "-m", "protein", "--pdb", os.path.join(pdb_dir, "1UBQ.pdb"),
         "-c", "A", "--sse", "-o", str(run_out), "--no-triplets", "-q",
         "--sample-id", "ubq"]
    )[0] == EXIT_OK
    import filecmp

    names = sorted(os.listdir(str(run_out)))
    assert names == ["ubq.keys", "ubq.meta.json", "ubq.run.json"]
    for name in names:
        if name.endswith(".run.json"):
            continue
        assert filecmp.cmp(
            os.path.join(str(batch_out), name),
            os.path.join(str(run_out), name),
            shallow=False,
        ), name


# ---------------------------------------------------------------------------
# Provenance survives serialization
# ---------------------------------------------------------------------------


def test_provenance_survives_the_round_trip(tmp_path, pdb_dir):
    out = tmp_path / "out"
    cli(
        ["run", "-m", "protein", "--pdb", os.path.join(pdb_dir, "1BRS.pdb"),
         "-c", "A", "--sse", "--mirror", "--size-filter", "20", "-o", str(out),
         "--no-triplets", "-q", "--label", "barnase"]
    )
    prefix = "1BRS_A_protein_mirror_size20_sse"
    payload = json.load(open(os.path.join(str(out), prefix + ".run.json")))
    assert payload["prefix"] == prefix and payload["complete"] is True
    assert payload["method"] == "protein"
    assert "output_format" not in payload
    assert payload["profile"]["resolved"] == "standard_protein_calpha_v1"
    assert payload["profile"]["active"].startswith("standard_protein_calpha_v1+")
    assert payload["profile"]["sse_strand_identity"] == "sheet_id_and_serial"
    assert payload["parameters"]["mirror"] is True
    assert payload["parameters"]["size_filter"] == 20.0
    assert payload["parameters"]["sse"] is True
    assert payload["label"] == "barnase"

    sidecar = json.load(
        open(os.path.join(str(out), payload["results"][0]["files"]["metadata"]))
    )
    assert sidecar["key_space_fingerprint"] == payload["results"][0][
        "key_space_fingerprint"
    ]
    assert sidecar["sampling_fingerprint"] == payload["results"][0][
        "sampling_fingerprint"
    ]
    assert sidecar["sse"] is not None
    assert sidecar["sse_coverage"]["n_helix"] > 0


def test_the_completion_summary_names_the_standard_profile(tmp_path, pdb_dir):
    out = tmp_path / "out"
    _code, _stdout, err = cli(
        ["run", "-m", "protein", "--pdb", os.path.join(pdb_dir, "1UBQ.pdb"),
         "-c", "A", "-o", str(out), "--no-triplets"]
    )
    assert "standard_protein_calpha_v1" in err
    assert "compatibility" not in err


# ---------------------------------------------------------------------------
# describe
# ---------------------------------------------------------------------------


def test_describe_reads_a_completion_record(tmp_path, pdb_dir):
    out = tmp_path / "out"
    _code, stdout, _err = cli(
        ["run", "-m", "protein", "--pdb", os.path.join(pdb_dir, "1UBQ.pdb"),
         "-c", "A", "-o", str(out), "--no-triplets", "-q"]
    )
    sample = stdout.strip()
    code, text, _err = cli(["describe", sample])
    assert code == EXIT_OK
    assert "protein" in text
    assert "standard_protein_calpha_v1" in text


def test_describe_reads_a_sidecar(tmp_path, pdb_dir):
    out = tmp_path / "out"
    _code, stdout, _err = cli(
        ["run", "-m", "protein", "--pdb", os.path.join(pdb_dir, "1UBQ.pdb"),
         "-c", "A", "-o", str(out), "--no-triplets", "-q"]
    )
    record = stdout.strip()
    sidecar = record[: -len(".run.json")] + ".meta.json"
    code, text, _err = cli(["describe", sidecar])
    assert code == EXIT_OK
    assert "key space" in text
    assert "protein_calpha" in text


def test_describe_refuses_a_directory_that_is_not_a_sample(tmp_path):
    code, _out, err = cli(["describe", str(tmp_path)])
    assert code == EXIT_USAGE
    assert "run.json" in err


# ---------------------------------------------------------------------------
# analyze common-keys
# ---------------------------------------------------------------------------


def test_common_keys_matches_the_analysis_api(tmp_path, pdb_dir):
    out = tmp_path / "out"
    code, _stdout, _err = cli(
        ["analyze", "common-keys", "-M", manifest(tmp_path, GOOD),
         "--pdb-dir", pdb_dir, "-o", str(out), "-q"]
    )
    assert code == EXIT_OK
    payload = json.load(
        open(os.path.join(str(out), "analysis", "common_keys", "all_samples.json"))
    )

    from biotsr.analysis import AnalysisSample, common_key_intersection

    samples = [
        AnalysisSample.from_result(
            biotsr.analyze(
                os.path.join(pdb_dir, "%s.pdb" % pdb), method="protein",
                chain=chain, keep_triplets=False,
            ),
            sample_id=sample_id,
        )
        for sample_id, pdb, chain in
        (("ubq", "1UBQ", "A"), ("brs_a", "1BRS", "A"), ("brs_b", "1BRS", "B"))
    ]
    expected = common_key_intersection(samples)
    assert payload["n_common_keys"] == expected.n_common_keys
    assert payload["keys"] == list(expected.keys)
    assert {r["fileName"] for r in payload["rows"]} == {
        "ubq", "brs_a", "brs_b"
    }
    for row, reference in zip(
        sorted(payload["rows"], key=lambda r: r["fileName"]),
        sorted(expected.to_dict()["rows"], key=lambda r: r["fileName"]),
    ):
        assert row == reference


def test_common_keys_writes_a_csv_beside_the_json(tmp_path, pdb_dir):
    out = tmp_path / "out"
    cli(
        ["analyze", "common-keys", "-M", manifest(tmp_path, GOOD),
         "--pdb-dir", pdb_dir, "-o", str(out), "-q"]
    )
    rows = list(
        csv.DictReader(
            open(os.path.join(str(out), "analysis", "common_keys", "all_samples.csv"))
        )
    )
    assert len(rows) == 3
    assert set(rows[0]) == {
        "sample_id", "method", "total_keys", "common_keys", "sum_common_keys_freq",
        "key_space_fingerprint", "processing_fingerprint", "sampling_fingerprint",
    }


def test_common_keys_groups_by_label(tmp_path, pdb_dir):
    out = tmp_path / "out"
    code, _stdout, _err = cli(
        ["analyze", "common-keys", "-M", manifest(tmp_path, GOOD),
         "--pdb-dir", pdb_dir, "-o", str(out), "-q", "--group-by", "label"]
    )
    assert code == EXIT_OK
    root = os.path.join(str(out), "analysis", "common_keys")
    assert sorted(os.listdir(root)) == [
        "barnase.csv", "barnase.json", "ubiquitin.csv", "ubiquitin.json"
    ]
    barnase = json.load(open(os.path.join(root, "barnase.json")))
    assert barnase["n_samples"] == 2


def test_common_keys_group_by_needs_every_row_labelled(tmp_path, pdb_dir):
    out = tmp_path / "out"
    text = (
        "sample_id,pdb_id,method,chain,label\n"
        "ubq,1UBQ,protein,A,ubiquitin\n"
        "brs,1BRS,protein,A,\n"
    )
    code, _stdout, err = cli(
        ["analyze", "common-keys", "-M", manifest(tmp_path, text),
         "--pdb-dir", pdb_dir, "-o", str(out), "-q", "--group-by", "label"]
    )
    assert code == EXIT_USAGE
    assert "no label" in err


def test_common_keys_stops_rather_than_answering_over_a_subset(tmp_path, pdb_dir):
    """A silently reduced sample set is a wrong intersection, not a partial one."""
    out = tmp_path / "out"
    code, _stdout, err = cli(
        ["analyze", "common-keys", "-M", manifest(tmp_path, MIXED),
         "--pdb-dir", pdb_dir, "-o", str(out), "-q"]
    )
    assert code == EXIT_USAGE
    assert not os.path.exists(os.path.join(str(out), "analysis"))


def test_common_keys_needs_at_least_two_samples(tmp_path, pdb_dir):
    out = tmp_path / "out"
    text = "sample_id,pdb_id,method,chain\nubq,1UBQ,protein,A\n"
    code, _stdout, err = cli(
        ["analyze", "common-keys", "-M", manifest(tmp_path, text),
         "--pdb-dir", pdb_dir, "-o", str(out), "-q"]
    )
    assert code == EXIT_USAGE
    assert "at least two" in err


def test_common_keys_refuses_a_mixed_key_space(tmp_path, pdb_dir):
    """The analysis layer's own precondition, surfaced rather than worked around."""
    out = tmp_path / "out"
    text = (
        "sample_id,pdb_id,method,chain,entity_name,entity_id\n"
        "ubq,1UBQ,protein,A,,\n"
        "ben,3PTB,drug,A,BEN,1\n"
    )
    code, _stdout, err = cli(
        ["analyze", "common-keys", "-M", manifest(tmp_path, text),
         "--pdb-dir", pdb_dir, "-o", str(out), "-q"]
    )
    assert code == EXIT_USAGE
    assert "Traceback" not in err


# ---------------------------------------------------------------------------
# Working-directory independence
# ---------------------------------------------------------------------------


def test_the_cli_works_from_any_directory(tmp_path, pdb_dir, monkeypatch):
    elsewhere = tmp_path / "somewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    out = tmp_path / "out"
    code, _stdout, _err = cli(
        ["run", "-m", "protein", "--pdb", os.path.join(pdb_dir, "1UBQ.pdb"),
         "-c", "A", "-o", str(out), "--no-triplets", "-q"]
    )
    assert code == EXIT_OK
    assert os.listdir(elsewhere) == []


def test_relative_paths_resolve_against_the_working_directory(
    tmp_path, pdb_dir, monkeypatch
):
    import shutil

    work = tmp_path / "work"
    work.mkdir()
    shutil.copy(os.path.join(pdb_dir, "1UBQ.pdb"), str(work / "1UBQ.pdb"))
    monkeypatch.chdir(work)
    code, stdout, _err = cli(
        ["run", "-m", "protein", "--pdb", "1UBQ.pdb", "-c", "A", "-o", "results",
         "--no-triplets", "-q"]
    )
    assert code == EXIT_OK
    assert os.path.isfile(str(work / "results" / "1UBQ_A_protein.run.json"))
    assert stdout.strip() == os.path.join("results", "1UBQ_A_protein.run.json")
