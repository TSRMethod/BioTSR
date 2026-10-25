"""A sample id is a prefix inside --output, and nothing outside it is touched.

Regression tests for the rc1 defect in which a manifest row with
``sample_id=../precious`` and ``--if-exists overwrite`` replaced -- and deleted the
contents of -- a sibling directory of ``--output``, and ``--sample-id /abs/path``
wrote outside ``--output``.  Every test plants a sentinel outside the output
root and checks it is byte-for-byte unchanged afterwards.
"""

from __future__ import annotations

import csv
import json
import os

import pytest

import biotsr
from biotsr import batch as biotsr_batch
from biotsr.cli.errors import CLIError
from biotsr.cli.main import main
from biotsr.cli.output import (
    record_name,
    validate_prefix,
    verify_existing,
)

from ..conftest import structure


@pytest.fixture
def sandbox(tmp_path):
    """``tmp/out`` is the output root; ``tmp/precious/file.txt`` is the sentinel."""
    precious = tmp_path / "precious"
    precious.mkdir()
    (precious / "file.txt").write_text("keep me\n")
    return tmp_path


def _sentinel_intact(sandbox):
    assert (sandbox / "precious" / "file.txt").read_text() == "keep me\n"
    assert sorted(os.listdir(sandbox / "precious")) == ["file.txt"]


def _manifest(path, rows, header=("sample_id", "pdb_path", "method", "chain")):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)
    return str(path)


def _summary(out):
    with open(os.path.join(out, "batch_summary.csv")) as fh:
        return list(csv.DictReader(fh))


@pytest.mark.parametrize("bad", [
    "../precious", "..", ".", "/abs/precious", "a/b", "a\\b", ".hidden",
    "tab\tname", " lead", "trail ", "", "analysis", "Batch_Summary.csv",
    "fetch_report.csv", ".biotsr-batch.lock", "x" * 201,
])
def test_an_unusable_sample_id_is_refused_by_name(bad):
    with pytest.raises(CLIError) as excinfo:
        validate_prefix(bad)
    assert "sample id" in excinfo.value.message


@pytest.mark.parametrize("good", ["ubq_regular", "1UBQ_A_protein", "a.b-c_d", "x" * 200])
def test_an_ordinary_sample_id_is_accepted(good):
    assert validate_prefix(good) == good


def test_run_refuses_a_traversing_sample_id_and_touches_nothing(sandbox):
    out = sandbox / "out"
    code = main(["run", "--method", "protein", "--pdb", structure("1UBQ"),
                 "--chain", "A", "--sample-id", "../precious", "--output", str(out),
                 "--if-exists", "overwrite", "--no-triplets", "-q"])
    assert code == 2
    _sentinel_intact(sandbox)
    assert not out.exists()


def test_run_refuses_an_absolute_sample_id(sandbox):
    target = sandbox / "elsewhere"
    code = main(["run", "--method", "protein", "--pdb", structure("1UBQ"),
                 "--chain", "A", "--sample-id", str(target),
                 "--output", str(sandbox / "out"), "--no-triplets", "-q"])
    assert code == 2
    assert not target.exists()


def test_batch_overwrite_cannot_replace_a_sibling_of_the_output(sandbox):
    """The exact rc1 reproduction: the row fails, the sibling survives."""
    manifest = _manifest(sandbox / "m.csv", [
        ["../precious", structure("1UBQ"), "protein", "A"],
        ["fine", structure("1UBQ"), "protein", "A"],
    ])
    out = sandbox / "out"
    code = main(["batch", "--manifest", manifest, "--output", str(out),
                 "--if-exists", "overwrite", "--no-triplets", "-q"])
    assert code == 3
    _sentinel_intact(sandbox)
    rows = {r["rows"]: r for r in _summary(out)}
    assert rows["2"]["status"] == "failed"
    assert "cannot be used as an output prefix" in rows["2"]["message"]
    assert rows["3"]["status"] == "ok"
    assert sorted(os.listdir(out)) == ["batch_failures.csv", "batch_summary.csv",
                                       "fine.keys", "fine.meta.json", "fine.run.json"]


def test_a_dry_run_reports_the_unusable_id_and_writes_nothing(sandbox):
    manifest = _manifest(sandbox / "m.csv", [
        ["/abs/escape", structure("1UBQ"), "protein", "A"],
        ["../precious", structure("1UBQ"), "protein", "A"],
    ])
    out = sandbox / "out"
    code = main(["batch", "--manifest", manifest, "--output", str(out),
                 "--if-exists", "overwrite", "--dry-run"])
    assert code == 3
    assert not out.exists()          # no directory, no lock, no summary
    _sentinel_intact(sandbox)


def test_the_python_batch_api_applies_the_same_rule(sandbox):
    manifest = _manifest(sandbox / "m.csv", [
        ["../precious", structure("1UBQ"), "protein", "A"],
    ])
    report = biotsr_batch.run(manifest, str(sandbox / "out"), if_exists="overwrite",
                              triplets=False)
    assert report.n_failed == 1
    assert "output prefix" in report.failures[0]["message"]
    _sentinel_intact(sandbox)


@pytest.mark.parametrize("linked", ["sample.run.json", "sample.keys"])
def test_a_symbolic_link_at_a_target_name_is_never_acted_through(sandbox, linked):
    """String checks cannot see a link; the file-system check does."""
    out = sandbox / "out"
    out.mkdir()
    os.symlink(sandbox / "precious" / "file.txt", out / linked)
    code = main(["run", "--method", "protein", "--pdb", structure("1UBQ"),
                 "--chain", "A", "--sample-id", "sample", "--output", str(out),
                 "--if-exists", "overwrite", "--no-triplets", "-q"])
    assert code == 2
    assert os.path.islink(out / linked)
    _sentinel_intact(sandbox)


def test_a_symbolic_link_is_refused_in_a_batch_too(sandbox):
    out = sandbox / "out"
    out.mkdir()
    os.symlink(sandbox / "precious" / "file.txt", out / "linked.run.json")
    manifest = _manifest(sandbox / "m.csv", [
        ["linked", structure("1UBQ"), "protein", "A"],
    ])
    code = main(["batch", "--manifest", manifest, "--output", str(out),
                 "--if-exists", "overwrite", "--no-triplets", "-q"])
    assert code == 3
    _sentinel_intact(sandbox)
    assert "symbolic link" in _summary(out)[0]["message"]


def test_a_generated_name_is_held_to_the_same_rule(sandbox, tmp_path):
    """A structure file named ``.hidden.pdb`` would generate a hidden name."""
    hidden = tmp_path / ".hidden.pdb"
    hidden.write_bytes(open(structure("1UBQ"), "rb").read())
    code = main(["run", "--method", "protein", "--pdb", str(hidden), "--chain", "A",
                 "--output", str(sandbox / "out"), "--no-triplets", "-q"])
    assert code == 2
    assert not (sandbox / "out").exists()


def test_skip_never_trusts_a_record_that_names_files_outside_the_folder(tmp_path):
    out = tmp_path / "out"
    code = main(["run", "--method", "protein", "--pdb", structure("1UBQ"),
                 "--chain", "A", "--output", str(out), "--no-triplets", "-q"])
    assert code == 0
    record = out / record_name("1UBQ_A_protein")
    payload = json.loads(record.read_text())
    request = payload["request"]
    assert verify_existing(str(out), "1UBQ_A_protein", request, False) is None
    payload["results"][0]["files"]["keys"] = "../../elsewhere"
    record.write_text(json.dumps(payload))
    assert "outside the folder" in verify_existing(str(out), "1UBQ_A_protein",
                                                   request, False)


def test_the_version_under_test_is_the_release():
    assert biotsr.__version__ == "1.0.0"
