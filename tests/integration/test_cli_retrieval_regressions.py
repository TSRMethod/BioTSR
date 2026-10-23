"""Regression tests for smaller command-line and retrieval defects.

* ``--jobs auto`` honoured a scheduler variable larger than the machine;
* ``--on-error stop`` exited 2 rather than the documented partial-batch 3;
* the ``invalid_cached`` remedy ("replace it yourself") could not work, because
  the download record still held the old checksum;
* an obsolete PDB entry downloaded with no notice;
* ``analyze common-keys`` accepted mixed sampling despite its help, recorded no
  provenance, and silently ignored computing flags with ``--from-results``;
* ``analyze common-keys --manifest`` was documented as writing triplet files.
"""

from __future__ import annotations

import csv
import json
import os

import pytest

import biotsr
from biotsr.cli import parallel
from biotsr.cli.main import main
from biotsr.retrieval import PROVENANCE_DIR

from ..conftest import STRUCTURES, WORKSPACE, structure
from ..support.stub_server import Stub, _entry, _gz


def test_jobs_auto_never_exceeds_the_machine(monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", str((os.cpu_count() or 1) * 50))
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    count, _source = parallel.available_cpus()
    assert 1 <= count <= (os.cpu_count() or 1)


def test_on_error_stop_exits_with_the_partial_batch_code(tmp_path, capsys):
    manifest = tmp_path / "m.csv"
    manifest.write_text(
        "sample_id,pdb_path,method,chain,residues\n"
        "good,%s,aa_selection,A,1-5\nbad,%s,aa_selection,A,600\n"
        "later,%s,aa_selection,A,6-9\n"
        % (structure("1UBQ"), structure("1UBQ"), structure("1UBQ")))
    code = main(["batch", "--manifest", str(manifest), "--output",
                 str(tmp_path / "out"), "--no-triplets", "--on-error", "stop"])
    assert code == 3
    err = capsys.readouterr().err
    assert err.count("600") >= 1
    assert "stopped at its first failed request" in err
    with open(tmp_path / "out" / "batch_summary.csv") as fh:
        statuses = {r["sample_id"]: r["status"] for r in csv.DictReader(fh)}
    assert statuses == {"good": "ok", "bad": "failed", "later": "not_started"}


@pytest.fixture
def archive():
    server = Stub()
    yield server
    server.close()


def test_the_invalid_cache_message_names_the_record_to_delete(tmp_path, archive):
    archive.ok("1ABC")
    report = biotsr.fetch_structures(["1ABC"], str(tmp_path), base_url=archive.base)
    assert report.records[0].status == "downloaded"
    (tmp_path / "1ABC.pdb").write_text("ATOM tampered\nEND\n")
    again = biotsr.fetch_structures(["1ABC"], str(tmp_path), base_url=archive.base)
    record = again.records[0]
    assert record.status == "invalid_cached"
    assert os.path.join(PROVENANCE_DIR, "1ABC.json") in record.message
    assert "--refresh" in record.message


def test_an_obsolete_entry_is_downloaded_with_a_notice(tmp_path, archive):
    text = _entry("1OLD").decode().replace(
        "END\n",
        "",
    )
    lines = text.splitlines()
    lines.insert(1, "OBSLTE     01-JAN-20 1OLD      2NEW")
    body = ("\n".join(lines) + "\nEND\n").encode()
    archive.serve("1OLD", {"status": 200, "body": _gz(body)})
    report = biotsr.fetch_structures(["1OLD"], str(tmp_path), base_url=archive.base)
    record = report.records[0]
    assert record.status == "downloaded"
    assert record.message.startswith("note: this entry is obsolete")
    assert "2NEW" in record.message


# ---------------------------------------------------------------------------
# analyze common-keys
# ---------------------------------------------------------------------------


def _mixed_manifest(tmp_path):
    manifest = tmp_path / "m.csv"
    manifest.write_text(
        "sample_id,pdb_path,method,chain,size_filter\n"
        "plain,%s,protein,A,\nfiltered,%s,protein,A,25\n"
        % (structure("1UBQ"), structure("1UBQ")))
    return str(manifest)


def test_mixed_processing_is_refused_by_default(tmp_path, capsys):
    code = main(["analyze", "common-keys", "--manifest", _mixed_manifest(tmp_path),
                 "--output", str(tmp_path / "out"), "-q"])
    assert code == 2
    err = capsys.readouterr().err
    assert "size_filter" in err and "--allow-mixed-processing" in err
    assert not (tmp_path / "out" / "analysis").exists()


def test_mixed_processing_can_be_allowed_and_is_recorded(tmp_path, capsys):
    out = tmp_path / "out"
    code = main(["analyze", "common-keys", "--manifest", _mixed_manifest(tmp_path),
                 "--output", str(out), "--allow-mixed-processing"])
    assert code == 0
    assert "not all processed with the same settings" in capsys.readouterr().err
    report = json.loads((out / "analysis" / "common_keys" / "all_samples.json")
                        .read_text())
    assert report["comparison"]["mixed_processing_allowed"] is True
    assert report["comparison"]["n_processing_policies"] == 2
    assert report["warnings"]
    by_id = {s["sample_id"]: s for s in report["samples"]}
    assert by_id["filtered"]["sampling"]["size_filter"] == 25.0
    assert "size_filter" not in by_id["plain"]["sampling"]
    assert by_id["plain"]["key_space_fingerprint"] == \
        by_id["filtered"]["key_space_fingerprint"]
    with open(out / "analysis" / "common_keys" / "all_samples.csv") as fh:
        header = next(csv.reader(fh))
    assert "sampling_fingerprint" in header and "key_space_fingerprint" in header
    assert "processing_fingerprint" in header


def test_different_structures_and_chains_are_not_mixed_processing(tmp_path):
    manifest = tmp_path / "m.csv"
    manifest.write_text(
        "sample_id,pdb_path,method,chain,chain_2\n"
        "ad,%s,protein_protein_calpha_cross,A,D\n"
        "be,%s,protein_protein_calpha_cross,B,E\n"
        % (structure("1BRS"), structure("1BRS")))
    out = tmp_path / "out"
    code = main(["analyze", "common-keys", "--manifest", str(manifest),
                 "--output", str(out), "-q"])
    assert code == 0
    report = json.loads((out / "analysis" / "common_keys" / "all_samples.json")
                        .read_text())
    assert report["comparison"]["processing_uniform"] is True
    fingerprints = {s["sampling_fingerprint"] for s in report["samples"]}
    assert len(fingerprints) == 2        # the chains differ ...
    assert len({s["processing_fingerprint"] for s in report["samples"]}) == 1


def test_different_key_spaces_are_refused_even_with_the_opt_in(tmp_path, capsys):
    manifest = tmp_path / "m.csv"
    manifest.write_text(
        "sample_id,pdb_path,method,chain,mirror\n"
        "plain,%s,protein,A,\nmirrored,%s,protein,A,true\n"
        % (structure("1UBQ"), structure("1UBQ")))
    code = main(["analyze", "common-keys", "--manifest", str(manifest),
                 "--output", str(tmp_path / "out"), "--allow-mixed-processing", "-q"])
    assert code == 2
    assert "cannot be compared" in capsys.readouterr().err


def test_common_keys_writes_keys_only_and_has_no_triplet_option(tmp_path):
    manifest = tmp_path / "m.csv"
    manifest.write_text(
        "sample_id,pdb_path,method,chain\na,%s,protein,A\nb,%s,protein,D\n"
        % (structure("1UBQ"), structure("1BRS")))
    out = tmp_path / "out"
    assert main(["analyze", "common-keys", "--manifest", str(manifest),
                 "--output", str(out), "-q"]) == 0
    assert not [n for n in os.listdir(out) if n.endswith(".triplets")]
    for sample in ("a", "b"):
        assert json.loads((out / (sample + ".run.json")).read_text())["triplets"] is False
    with pytest.raises(SystemExit):
        main(["analyze", "common-keys", "--manifest", str(manifest),
              "--output", str(out), "--no-triplets"])


def test_the_documentation_no_longer_claims_common_keys_writes_triplets():
    docs = os.path.join(WORKSPACE, "docs")
    if not os.path.isdir(docs):
        pytest.skip("repository-only evidence not present: docs")
    for name in ("KNOWN_LIMITATIONS.md", "CLI_GUIDE.md", "USER_GUIDE.md"):
        text = open(os.path.join(docs, name)).read()
        for line in text.splitlines():
            if "common-keys" in line and "--manifest" in line:
                assert "triplet" not in line or "no triplet" in line \
                    or "keys only" in line or "not" in line, (name, line)


def test_from_results_refuses_flags_it_would_ignore(tmp_path, capsys):
    out = tmp_path / "saved"
    manifest = tmp_path / "m.csv"
    manifest.write_text(
        "sample_id,pdb_path,method,chain\na,%s,protein,A\nb,%s,protein,D\n"
        % (structure("1UBQ"), structure("1BRS")))
    assert main(["batch", "--manifest", str(manifest), "--output", str(out),
                 "--no-triplets", "-q"]) == 0
    for flags in (["--mirror"], ["--jobs", "2"],
                  ["--pdb-dir", STRUCTURES], ["--if-exists", "skip"],
                  ["--timeout", "5"]):
        code = main(["analyze", "common-keys", "--from-results", str(out),
                     "--output", str(tmp_path / "a")] + flags)
        assert code == 2, flags
        assert "would have no effect" in capsys.readouterr().err
    assert main(["analyze", "common-keys", "--from-results", str(out),
                 "--output", str(tmp_path / "a"), "-q"]) == 0
