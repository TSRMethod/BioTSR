"""The output-folder contract (ADR-098).

One shared folder; one prefix per calculation; files staged privately and
installed with the completion record written last; nothing replaced or deleted
that the calculation does not own; and what the completion record has to
contain for the calculation to be reproducible from disk alone.
"""

from __future__ import annotations

import json
import os

import pytest

from biotsr.cli import output as cli_output
from biotsr.cli.errors import CLIError


_KIND = {"keys": "keys", "meta.json": "metadata", "triplets": "triplets"}


def _write(writer, names=("p.keys", "p.meta.json"), content="x"):
    """Stage ``names`` and a minimal completion record for ``writer``."""
    for name in names:
        with open(os.path.join(writer.path, name), "w") as fh:
            fh.write(content)
    writer.payload = {
        "schema": cli_output.REQUEST_SCHEMA, "prefix": writer.prefix,
        "results": [{"name": writer.prefix,
                     "files": {_KIND[n.split(".", 1)[1]]: n for n in names},
                     "file_bytes": {n: len(content) for n in names}}],
    }


def _record(root, prefix="p"):
    return json.load(open(os.path.join(root, prefix + ".run.json")))


def test_a_calculation_is_staged_and_installed_with_its_record_last(tmp_path):
    root = str(tmp_path / "out")
    writer = cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_FAIL)
    with writer:
        assert os.path.basename(os.path.dirname(writer.path)).startswith(".biotsr-staging-")
        _write(writer)
        assert not os.path.exists(os.path.join(root, "p.keys"))
    assert sorted(os.listdir(root)) == ["p.keys", "p.meta.json", "p.run.json"]
    record = _record(root)
    assert record["complete"] is True and record["prefix"] == "p"
    assert cli_output.check_complete(root, "p") is None
    assert not os.path.exists(writer.staging)


def test_a_shared_staging_directory_is_left_to_the_invocation_that_owns_it(tmp_path):
    """Batch workers share one token and one staging directory.  A worker
    finishing while that directory is momentarily empty must not remove it:
    another worker may be creating its own staging inside it, which macOS
    refuses with EINVAL.  The batch removes the directory when it ends."""
    root = str(tmp_path / "out")
    token = cli_output.invocation_token()
    first = cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_FAIL,
                                         token=token, shared_staging=True)
    with first:
        _write(first)
    assert os.path.isdir(first.staging_root) and os.listdir(first.staging_root) == []
    second = cli_output.CalculationWriter(root, "q", cli_output.IF_EXISTS_FAIL,
                                          token=token, shared_staging=True)
    with second:
        _write(second, names=("q.keys", "q.meta.json"))
    assert cli_output.check_complete(root, "q") is None
    # A writer of its own invocation removes its staging directory.
    alone = cli_output.CalculationWriter(root, "r", cli_output.IF_EXISTS_FAIL)
    with alone:
        _write(alone, names=("r.keys", "r.meta.json"))
    assert not os.path.exists(alone.staging_root)


def test_an_interrupted_calculation_leaves_nothing_that_looks_finished(tmp_path):
    root = str(tmp_path / "out")
    writer = cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_FAIL)
    with pytest.raises(RuntimeError):
        with writer:
            _write(writer)
            raise RuntimeError("boom")
    assert not os.path.exists(root)        # the run created it, and leaves none


def test_an_interruption_during_installation_leaves_an_incomplete_record(tmp_path, monkeypatch):
    """Stopped between moving the data files and writing the final record:
    the record says incomplete, readers refuse it, and a resume replaces it."""
    root = str(tmp_path / "out")
    real = cli_output._write_json_atomically
    calls = []

    def flaky(root_, name, payload, token):
        calls.append(payload.get("complete"))
        if payload.get("complete"):
            raise KeyboardInterrupt
        return real(root_, name, payload, token)

    monkeypatch.setattr(cli_output, "_write_json_atomically", flaky)
    writer = cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_FAIL)
    with pytest.raises(KeyboardInterrupt):
        with writer:
            _write(writer)
    monkeypatch.setattr(cli_output, "_write_json_atomically", real)
    assert calls == [False, True]
    state, record = cli_output.existing_state(root, "p")
    assert state == "incomplete" and record["files_pending"] == ["p.keys", "p.meta.json"]
    assert "interrupted" in cli_output.check_complete(root, "p")
    with pytest.raises(CLIError, match="incomplete"):
        cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_FAIL).check()
    again = cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_OVERWRITE)
    with again:
        _write(again, content="yy")
    assert cli_output.check_complete(root, "p") is None
    assert open(os.path.join(root, "p.keys")).read() == "yy"


def test_an_existing_calculation_is_refused_by_default(tmp_path):
    root = str(tmp_path / "out")
    with cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_FAIL) as w:
        _write(w)
    writer = cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_FAIL)
    with pytest.raises(CLIError) as excinfo:
        with writer:
            pass
    assert "already exists" in excinfo.value.message
    assert "--if-exists" in excinfo.value.hint


def test_overwrite_replaces_owned_files_and_removes_owned_ones_no_longer_written(tmp_path):
    root = str(tmp_path / "out")
    with cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_FAIL) as w:
        _write(w, names=("p.keys", "p.meta.json", "p.triplets"))
    open(os.path.join(root, "notes.txt"), "w").write("mine")
    with cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_OVERWRITE) as w:
        _write(w, content="new")
    assert sorted(os.listdir(root)) == ["notes.txt", "p.keys", "p.meta.json", "p.run.json"]
    assert open(os.path.join(root, "p.keys")).read() == "new"
    assert open(os.path.join(root, "notes.txt")).read() == "mine"


def test_files_biotsr_did_not_write_are_never_replaced(tmp_path):
    """Not even with --if-exists overwrite: no record of this prefix owns them."""
    root = str(tmp_path / "out")
    os.makedirs(root)
    open(os.path.join(root, "p.keys"), "w").write("someone else's")
    for policy in cli_output.IF_EXISTS_CHOICES:
        with pytest.raises(CLIError, match="no BioTSR completion record"):
            cli_output.CalculationWriter(root, "p", policy).check()
    assert open(os.path.join(root, "p.keys")).read() == "someone else's"


def test_a_name_differing_only_in_letter_case_is_a_collision(tmp_path):
    root = str(tmp_path / "out")
    with cli_output.CalculationWriter(root, "P", cli_output.IF_EXISTS_FAIL) as w:
        _write(w, names=("P.keys", "P.meta.json"))
    with pytest.raises(CLIError, match="letter case"):
        cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_OVERWRITE).check()


def test_a_symbolic_link_at_a_target_name_is_never_written_through(tmp_path):
    root = str(tmp_path / "out")
    os.makedirs(root)
    target = tmp_path / "elsewhere.txt"
    target.write_text("precious")
    os.symlink(str(target), os.path.join(root, "p.run.json"))
    with pytest.raises(CLIError):
        cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_OVERWRITE).check()
    assert target.read_text() == "precious"


def test_skip_leaves_a_verified_calculation_untouched(tmp_path):
    root = str(tmp_path / "out")
    with cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_FAIL) as w:
        _write(w, content="old")
    writer = cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_SKIP,
                                          verify=lambda prefix: None)
    with writer:
        assert writer.skipped
    assert open(os.path.join(root, "p.keys")).read() == "old"


def test_a_stale_staging_directory_of_a_dead_process_is_cleared(tmp_path):
    root = str(tmp_path / "out")
    os.makedirs(root)
    dead = "999999-%s-0000abcd" % cli_output.host_tag()
    leftover = os.path.join(root, ".biotsr-staging-" + dead)
    os.makedirs(os.path.join(leftover, "p"))
    lock = os.path.join(root, ".biotsr-lock-p")
    open(lock, "w").write(dead)
    removed, foreign = cli_output.remove_stale_leftovers(root)
    assert sorted(removed) == sorted([leftover, lock])
    assert foreign == []


def test_nothing_is_written_outside_the_output_root(tmp_path, monkeypatch):
    """The CLI never writes to the working directory."""
    elsewhere = tmp_path / "cwd"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    root = str(tmp_path / "out")
    with cli_output.CalculationWriter(root, "p", cli_output.IF_EXISTS_FAIL) as w:
        _write(w)
    assert os.listdir(elsewhere) == []


@pytest.mark.parametrize("bad", ["", "../x", "a/b", ".hidden", "analysis",
                                 "BATCH_SUMMARY.csv", "x" * 201])
def test_a_prefix_is_a_name_inside_the_folder(bad):
    with pytest.raises(CLIError):
        cli_output.validate_prefix(bad)


# ---------------------------------------------------------------------------
# The completion record
# ---------------------------------------------------------------------------


def _run(tmp_path, ubq_pdb, argv_extra=(), **kwargs):
    from biotsr.cli.spec import build_spec, execute

    values = dict(
        {"method": "protein", "pdb": ubq_pdb, "chain": "A", "keep_triplets": False},
        **kwargs
    )
    spec = build_spec(values)
    outcome = execute(spec)
    out_dir = str(tmp_path / "s")
    os.makedirs(out_dir)
    payload = cli_output.write_outputs(
        outcome, out_dir, prefix="s", triplets=False, argv=list(argv_extra)
    )
    return payload, out_dir


def test_the_record_records_everything_needed_to_reproduce(tmp_path, ubq_pdb):
    payload, out_dir = _run(
        tmp_path, ubq_pdb, argv_extra=["run", "-m", "protein"], mirror=True,
        size_filter=25.0,
    )
    assert payload["prefix"] == "s"
    assert payload["schema"] == cli_output.REQUEST_SCHEMA
    assert payload["method"] == "protein"
    assert set(payload["profile"]) == {"resolved", "active", "sse_strand_identity"}
    assert "output_format" not in payload
    assert payload["profile"]["resolved"] == "standard_protein_calpha_v1"
    assert payload["profile"]["active"] == (
        "standard_protein_calpha_v1+mirror+size_filter(25)"
    )
    assert payload["parameters"]["mirror"] is True
    assert payload["parameters"]["size_filter"] == 25.0
    assert payload["argv"] == ["run", "-m", "protein"]
    assert payload["biotsr_version"]
    entry = payload["results"][0]
    assert entry["files"] == {"keys": "s.keys", "metadata": "s.meta.json"}
    assert entry["file_bytes"]["s.keys"] == os.path.getsize(os.path.join(out_dir, "s.keys"))
    assert entry["key_space_fingerprint"]
    assert entry["sampling_fingerprint"]
    assert entry["warnings"] == []


def test_the_record_records_only_what_the_user_asked_for(tmp_path, ubq_pdb):
    payload, _ = _run(tmp_path, ubq_pdb)
    assert "mirror" not in payload["parameters"]
    assert "size_filter" not in payload["parameters"]
    assert "sse" not in payload["parameters"]


def test_the_record_records_the_strand_identity_only_when_sse_is_on(tmp_path, ubq_pdb):
    payload, _ = _run(tmp_path, ubq_pdb)
    assert payload["profile"]["sse_strand_identity"] is None
    payload, _ = _run(tmp_path / "b", ubq_pdb, sse=True)
    assert payload["profile"]["sse_strand_identity"] == "sheet_id_and_serial"


def test_the_native_sidecar_is_written_and_is_the_packages_own(tmp_path, ubq_pdb):
    from biotsr import analyze

    payload, out_dir = _run(tmp_path, ubq_pdb)
    sidecar = [f for f in os.listdir(out_dir) if f.endswith(".meta.json")]
    assert len(sidecar) == 1
    written = json.load(open(os.path.join(out_dir, sidecar[0])))
    direct = analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False)
    assert written == json.loads(json.dumps(direct.metadata()))
