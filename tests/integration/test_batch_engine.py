"""The batch engine: column mapping, chains, output names, verified resume,
dry runs, downloads, interruption, memory, and serial/parallel equivalence.

Every test drives the real command line in-process (``biotsr.cli.main``) or the
Python entry point (``biotsr.batch``), over the small tracked structures.
"""

from __future__ import annotations

import builtins
import csv
import gc
import json
import os
import shutil

import pytest

import biotsr
from biotsr import batch as batch_api
from biotsr.cli import main
from biotsr.cli import batch as engine
from biotsr.cli import output as output_mod
from biotsr.cli import parallel
from biotsr.cli.errors import CLIError
from biotsr.cli.manifest import read_manifest

from ..conftest import STRUCTURES, structure
from ..support.stub_server import Stub, _entry

SCOP_HEADER = "protein,chain,group,group1\n"


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def summary(out):
    with open(os.path.join(str(out), "batch_summary.csv"), newline="") as fh:
        return list(csv.DictReader(fh))


def run_json(out, prefix):
    """The completion record ``PREFIX.run.json`` in the shared output folder."""
    with open(os.path.join(str(out), prefix + ".run.json")) as fh:
        return json.load(fh)


def calculation_files(out, prefix):
    """Every file of one calculation in the shared output folder, sorted."""
    return sorted(n for n in os.listdir(str(out)) if n.startswith(prefix + "."))


def cli(*argv):
    return main([str(a) for a in argv])


@pytest.fixture(scope="module")
def pdb_dir(tmp_path_factory):
    """A private copy of the structures, so a test may add or break files."""
    target = tmp_path_factory.mktemp("pdb")
    for name in ("1UBQ", "1BRS", "3PTB"):
        shutil.copy(structure(name), str(target))
    return str(target)


# ---------------------------------------------------------------------------
# A dataset sheet with extra columns
# ---------------------------------------------------------------------------


SCOP = SCOP_HEADER + "1UBQ,A,alpha,0\n1BRS,A,beta,1\n1BRS,D,beta,1\n"


@pytest.mark.parametrize("label,kept", [("group", "group1"), ("group1", "group")])
def test_a_dataset_sheet_runs_with_either_label_column(tmp_path, pdb_dir, label, kept):
    manifest = write(tmp_path, "scop.csv", SCOP)
    out = tmp_path / "out"
    code = cli("batch", "--manifest", manifest, "--method", "protein",
               "--column", "label=%s" % label, "--keep-column", kept,
               "--pdb-dir", pdb_dir, "-o", out, "--no-triplets", "-q")
    assert code == 0
    rows = summary(out)
    assert [r["sample_id"] for r in rows] == [
        "1UBQ_A_protein", "1BRS_A_protein", "1BRS_D_protein"]
    assert [r[kept] for r in rows] == (["0", "1", "1"] if kept == "group1"
                                       else ["alpha", "beta", "beta"])
    assert all(r["chain_source"] == "row" for r in rows)
    assert [(r["request"], r["rows"]) for r in rows] == [("1", "2"), ("2", "3"), ("3", "4")]
    payload = run_json(out, "1UBQ_A_protein")
    assert payload["complete"] is True and payload["prefix"] == "1UBQ_A_protein"
    assert payload["source"]["rows"] == [2]
    assert payload["label"] == ("alpha" if label == "group" else "0")
    assert payload["metadata"] == {kept: "0" if kept == "group1" else "alpha"}
    # labels and metadata are never scientific
    assert "aa_grouping" not in payload["parameters"]
    assert payload["request"]["modifiers"] == {}
    assert "label" not in payload["request"] and "metadata" not in payload["request"]


def test_the_sheet_is_refused_until_its_extra_columns_are_given_a_role(tmp_path, pdb_dir,
                                                                       capsys):
    manifest = write(tmp_path, "scop.csv", SCOP)
    code = cli("batch", "--manifest", manifest, "--method", "protein",
               "--pdb-dir", pdb_dir, "-o", tmp_path / "out")
    err = capsys.readouterr().err
    assert code == 2
    assert "unrecognised column(s): group, group1" in err
    for remedy in ("--column FIELD=HEADER", "--keep-column HEADER", "--ignore-column HEADER"):
        assert remedy in err
    assert "aa_grouping'?" not in err      # a class label is never steered to grouping
    assert not (tmp_path / "out").exists()


def test_chain_values_are_preserved_exactly(tmp_path):
    manifest = write(tmp_path, "m.csv", SCOP_HEADER + "1ABC,1,alpha,0\n1ABD,a,beta,1\n")
    rows = read_manifest(manifest, defaults={"method": "protein"},
                         columns={"label": "group"}, keep_columns=["group1"])
    assert [r["chain"] for r in rows] == ["1", "a"]


@pytest.mark.parametrize("kwargs,message", [
    ({"columns": {"label": "group", "sample_id": "group"}}, "mapped to more than one field"),
    ({"columns": {"label": "nope"}}, "has no column 'nope'"),
    ({"columns": {"lable": "group"}}, "did you mean 'label'"),
    ({"columns": {"pdb_id": "group"}}, "also has column 'protein'"),
    ({"columns": {"label": "group"}, "keep_columns": ["group"]}, "both mapped"),
    ({"keep_columns": ["group1"], "ignore_columns": ["group1", "group"]},
     "both kept and ignored"),
    ({"keep_columns": ["chain"], "ignore_columns": ["group", "group1"]},
     "already a manifest column"),
    ({"keep_columns": ["missing"]}, "has no column 'missing'"),
])
def test_every_column_ambiguity_is_refused_with_its_own_message(tmp_path, kwargs, message):
    manifest = write(tmp_path, "m.csv", SCOP_HEADER + "1UBQ,A,alpha,0\n")
    with pytest.raises(CLIError) as excinfo:
        read_manifest(manifest, **kwargs)
    assert message in excinfo.value.message + excinfo.value.hint


def test_a_mapping_that_resolves_a_conflict_by_ignoring_the_alias(tmp_path):
    manifest = write(tmp_path, "m.csv", "protein,PDB,chain\nxxxx,1UBQ,A\n")
    rows = read_manifest(manifest, columns={"pdb_id": "PDB"}, ignore_columns=["protein"])
    assert rows[0]["pdb_id"] == "1UBQ"


def test_the_column_option_must_be_field_equals_header(tmp_path, pdb_dir, capsys):
    manifest = write(tmp_path, "scop.csv", SCOP)
    assert cli("batch", "--manifest", manifest, "--column", "label",
               "--pdb-dir", pdb_dir, "-o", tmp_path / "o") == 2
    assert "is not FIELD=HEADER" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# Chains
# ---------------------------------------------------------------------------


def test_the_chain_default_fills_only_rows_without_one(tmp_path, pdb_dir):
    manifest = write(tmp_path, "m.csv", "pdb_id,chain,sample_id\n1BRS,,a\n1BRS,D,d\n")
    out = tmp_path / "out"
    assert cli("batch", "--manifest", manifest, "--method", "protein", "--chain", "A",
               "--pdb-dir", pdb_dir, "-o", out, "--no-triplets", "-q") == 0
    rows = summary(out)
    assert [(r["chain"], r["chain_source"]) for r in rows] == [("A", "default"), ("D", "row")]


def test_rows_without_a_chain_are_skipped_never_guessed(tmp_path, pdb_dir, capsys):
    """ADR-112: a row whose chain is blank, with no --chain default, is not run
    and is listed as chain_not_specified; the other rows run; exit 4."""
    manifest = write(tmp_path, "ids.csv", "pdb_id,chain\n1UBQ,\n1BRS,A\n")
    code = cli("batch", "--manifest", manifest, "--method", "protein",
               "--pdb-dir", pdb_dir, "-o", tmp_path / "out", "--no-triplets")
    err = capsys.readouterr().err
    assert code == 4
    rows = summary(tmp_path / "out")
    assert [r["status"] for r in rows] == ["chain_not_specified", "ok"]
    assert rows[0]["chain"] == "" and "never assumes a chain" in rows[0]["message"]
    assert "--chain" in rows[0]["message"] and "protein needs chain" in err
    assert not [f for f in os.listdir(str(tmp_path / "out")) if f.startswith("1UBQ")]


def test_a_chain_the_structure_lacks_fails_that_row_and_names_the_chains(tmp_path, pdb_dir):
    manifest = write(tmp_path, "m.csv", "pdb_id,chain\n1UBQ,Z\n1UBQ,A\n")
    out = tmp_path / "out"
    assert cli("batch", "--manifest", manifest, "--method", "protein",
               "--pdb-dir", pdb_dir, "-o", out, "--no-triplets", "-q") == 3
    first, second = summary(out)
    assert first["status"] == "failed" and "chains in this structure: A" in first["message"]
    assert second["status"] == "ok"


# ---------------------------------------------------------------------------
# Output names
# ---------------------------------------------------------------------------


def test_rows_that_would_share_an_output_are_refused_before_computing(tmp_path, pdb_dir,
                                                                      capsys):
    """Two different requests that would write one prefix are refused at planning,
    under every --if-exists policy, before anything is computed or created.

    Premise changed in rc4 (ADR-098): automatic prefixes now name the modifiers,
    so rows differing only in ``mirror`` no longer collide -- they plan distinct
    prefixes (checked below with a dry run).  A collision now needs two different
    requests given one prefix: a repeated ``sample_id``, also when the repeats
    differ only in letter case.
    """
    out = tmp_path / "out"
    variants = write(tmp_path, "variants.csv", "pdb_id,chain,mirror\n1UBQ,A,\n1UBQ,A,true\n")
    assert cli("batch", "--manifest", variants, "--method", "protein",
               "--pdb-dir", pdb_dir, "-o", out, "--dry-run") == 0
    planned = capsys.readouterr().out
    assert "1UBQ_A_protein  [would run]" in planned
    assert "1UBQ_A_protein_mirror  [would run]" in planned
    for repeated in ("same,same", "same,SAME"):
        first, second = repeated.split(",")
        manifest = write(tmp_path, "m.csv", "sample_id,pdb_id,chain,mirror\n"
                         "%s,1UBQ,A,\n%s,1UBQ,A,true\n" % (first, second))
        for policy in ("fail", "skip", "overwrite"):
            code = cli("batch", "--manifest", manifest, "--method", "protein",
                       "--pdb-dir", pdb_dir, "-o", out, "--if-exists", policy)
            assert code == 2
            err = capsys.readouterr().err
            assert "more than one different request" in err
            assert "%s (rows 2, 3)" % first in err and "sample_id" in err
    assert not out.exists()


def test_names_that_differ_only_in_case_collide(tmp_path, pdb_dir, capsys):
    manifest = write(tmp_path, "m.csv", "sample_id,pdb_id,chain\nubq,1UBQ,A\nUBQ,1UBQ,A\n")
    assert cli("batch", "--manifest", manifest, "--method", "protein",
               "--pdb-dir", pdb_dir, "-o", tmp_path / "out") == 2
    assert "letter case" in capsys.readouterr().err


def test_the_summary_is_readable_by_others_like_any_output(tmp_path, pdb_dir):
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets") == 0
    current = os.umask(0o022)
    os.umask(current)
    mode = os.stat(str(tmp_path / "out" / "batch_summary.csv")).st_mode & 0o777
    assert mode == 0o666 & ~current


def test_distinct_sample_ids_make_variants_runnable(tmp_path, pdb_dir):
    manifest = write(tmp_path, "m.csv",
                     "sample_id,pdb_id,chain,mirror\nplain,1UBQ,A,\nmirrored,1UBQ,A,true\n")
    out = tmp_path / "out"
    assert cli("batch", "--manifest", manifest, "--method", "protein",
               "--pdb-dir", pdb_dir, "-o", out, "--no-triplets", "-q") == 0
    assert run_json(out, "plain")["request"]["modifiers"] == {}
    assert run_json(out, "mirrored")["request"]["modifiers"] == {"mirror": True}
    for prefix in ("plain", "mirrored"):
        assert calculation_files(out, prefix) == [
            prefix + ".keys", prefix + ".meta.json", prefix + ".run.json"]


# ---------------------------------------------------------------------------
# Verified resume
# ---------------------------------------------------------------------------


def _batch(tmp_path, pdb_dir, text, *extra):
    manifest = write(tmp_path, "m.csv", text)
    return cli("batch", "--manifest", manifest, "--method", "protein",
               "--pdb-dir", pdb_dir, "-o", tmp_path / "out", "-q", *extra)


UBQ_ROW = "sample_id,pdb_id,chain,mirror\nubq,1UBQ,A,\n"


def test_skip_keeps_an_output_of_exactly_the_same_request(tmp_path, pdb_dir):
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets") == 0
    out = tmp_path / "out"
    files = calculation_files(out, "ubq")
    assert files == ["ubq.keys", "ubq.meta.json", "ubq.run.json"]
    before = {name: os.stat(str(out / name)).st_mtime_ns for name in files}
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets", "--if-exists", "skip") == 0
    row, = summary(out)
    assert row["status"] == "skipped" and row["n_triangles"] == "70300"
    assert {name: os.stat(str(out / name)).st_mtime_ns
            for name in calculation_files(out, "ubq")} == before


@pytest.mark.parametrize("second,differs", [
    ("sample_id,pdb_id,chain,mirror\nubq,1UBQ,A,true\n", "modifiers"),
])
def test_skip_refuses_an_output_of_a_different_request(tmp_path, pdb_dir, second, differs):
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets") == 0
    assert _batch(tmp_path, pdb_dir, second, "--no-triplets", "--if-exists", "skip") == 3
    row, = summary(tmp_path / "out")
    assert row["status"] == "failed"
    assert "different request" in row["message"] and differs in row["message"]


@pytest.mark.parametrize("first,second,differs", [
    ("sample_id,pdb_id,method,chain\ns,1BRS,protein,A\n",
     "sample_id,pdb_id,method,chain\ns,1BRS,protein,D\n", "chain"),
    ("sample_id,pdb_id,method,chain,residues\ns,1UBQ,aa_selection,A,10-12\n",
     "sample_id,pdb_id,method,chain,residues\ns,1UBQ,aa_selection,A,10-13\n", "selection"),
    ("sample_id,pdb_id,method,chain,chain_2,contact_cutoff\n"
     "s,1BRS,protein_protein_cross,A,D,4.0\n",
     "sample_id,pdb_id,method,chain,chain_2,contact_cutoff\n"
     "s,1BRS,protein_protein_cross,A,D,5.0\n", "cross"),
    ("sample_id,pdb_id,method,chain\ns,1UBQ,protein,A\n",
     "sample_id,pdb_id,method,chain,residues\ns,1UBQ,aa_selection,A,10-12\n", "method"),
])
def test_skip_refuses_an_output_whose_calculation_differs(tmp_path, first, second, differs):
    """Chains, selections, contact settings and the method are all part of the
    request a resume compares, not just its name."""
    def run(text, *extra):
        manifest = write(tmp_path, "m.csv", text)
        return cli("batch", "--manifest", manifest, "--pdb-dir", STRUCTURES,
                   "-o", tmp_path / "out", "-q", "--no-triplets", *extra)

    assert run(first) == 0
    assert run(second, "--if-exists", "skip") == 3
    row, = summary(tmp_path / "out")
    assert "different request" in row["message"] and differs in row["message"]


@pytest.mark.parametrize("damage", ["truncate", "append"])
def test_skip_refuses_a_file_changed_since_it_was_written(tmp_path, pdb_dir, damage):
    """A file whose size differs from the one its completion record holds makes
    skip refuse (rc4: sizes are recorded per file in ``PREFIX.run.json``)."""
    assert _batch(tmp_path, pdb_dir, UBQ_ROW) == 0
    triplets = tmp_path / "out" / "ubq.triplets"
    assert run_json(tmp_path / "out", "ubq")["results"][0]["files"]["triplets"] == triplets.name
    data = triplets.read_bytes()
    triplets.write_bytes(data[:-10] if damage == "truncate" else data + b"extra\n")
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--if-exists", "skip") == 3
    message = summary(tmp_path / "out")[0]["message"]
    assert "ubq.triplets" in message and "incomplete or modified" in message
    assert triplets.stat().st_size != len(data)       # refused, not repaired


def test_skip_refuses_an_output_of_a_different_structure_file(tmp_path, pdb_dir):
    local = tmp_path / "pdb"
    local.mkdir()
    shutil.copy(os.path.join(pdb_dir, "1UBQ.pdb"), str(local))
    assert _batch(tmp_path, str(local), UBQ_ROW, "--no-triplets") == 0
    with open(str(local / "1UBQ.pdb"), "a") as fh:
        fh.write("REMARK 999 A DIFFERENT FILE\n")
    assert _batch(tmp_path, str(local), UBQ_ROW, "--no-triplets", "--if-exists", "skip") == 3
    assert "structure_sha256" in summary(tmp_path / "out")[0]["message"]


def test_a_keys_only_output_does_not_satisfy_a_request_for_triplets(tmp_path, pdb_dir):
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets") == 0
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--if-exists", "skip") == 3
    assert "--no-triplets" in summary(tmp_path / "out")[0]["message"]


def test_a_superset_output_satisfies_a_narrower_request(tmp_path, pdb_dir):
    assert _batch(tmp_path, pdb_dir, UBQ_ROW) == 0                      # with triplets
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets", "--if-exists", "skip") == 0
    assert summary(tmp_path / "out")[0]["status"] == "skipped"


def test_skip_refuses_an_incomplete_output(tmp_path, pdb_dir):
    """A complete record whose listed key file is gone is refused, never reused."""
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets") == 0
    out = tmp_path / "out"
    keys = run_json(out, "ubq")["results"][0]["files"]["keys"]
    assert keys == "ubq.keys"
    os.remove(str(out / keys))
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets", "--if-exists", "skip") == 3
    row, = summary(out)
    assert row["status"] == "failed"
    assert "its file ubq.keys is missing" in row["message"]
    assert not (out / keys).exists()


def test_an_interrupted_installation_is_refused_by_fail_and_recomputed_by_skip(tmp_path,
                                                                               pdb_dir):
    """An intent record (``"complete": false``) left by an interrupted install.

    New in rc4 (ADR-098): such a record is never read as a result.  The default
    --if-exists fail refuses it and changes nothing; --if-exists skip (resume)
    recomputes it, ending with a complete record and intact files.
    """
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets") == 0
    out = tmp_path / "out"
    record = out / "ubq.run.json"
    payload = json.loads(record.read_text())
    payload.update(complete=False, files_pending=["ubq.keys", "ubq.meta.json"])
    record.write_text(json.dumps(payload))
    (out / "ubq.keys").write_text("cut short")
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets") == 3
    row, = summary(out)
    assert row["status"] == "failed" and "left incomplete" in row["message"]
    assert (out / "ubq.keys").read_text() == "cut short"
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets", "--if-exists", "skip") == 0
    row, = summary(out)
    assert row["status"] == "ok" and row["n_triangles"] == "70300"
    payload = run_json(out, "ubq")
    assert payload["complete"] is True and "files_pending" not in payload
    sizes = payload["results"][0]["file_bytes"]
    assert sizes["ubq.keys"] == (out / "ubq.keys").stat().st_size > len("cut short")
    assert calculation_files(out, "ubq") == ["ubq.keys", "ubq.meta.json", "ubq.run.json"]


@pytest.mark.parametrize("damage,message", [
    ("no_request", "its completion record does not hold the request"),
    ("old_schema", "written by a different BioTSR version (schema 3; this one writes 4)"),
])
def test_skip_refuses_an_output_that_does_not_record_its_request(tmp_path, pdb_dir,
                                                                 damage, message):
    """A record without the request cannot be verified, so it is never reused.

    rc4: a record of another schema than 4 (an earlier version's) is "foreign"
    and refused with its own message; a schema-4 record without ``request`` is
    refused as unverifiable.
    """
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets") == 0
    path = tmp_path / "out" / "ubq.run.json"
    payload = json.loads(path.read_text())
    if damage == "no_request":
        del payload["request"]
    else:
        payload["schema"] = 3
    path.write_text(json.dumps(payload))
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets", "--if-exists", "skip") == 3
    row, = summary(tmp_path / "out")
    assert row["status"] == "failed" and message in row["message"]
    assert json.loads(path.read_text()) == payload           # left as it was


def test_run_verifies_before_skipping_too(tmp_path, pdb_dir, capsys):
    """``biotsr run --if-exists skip`` verifies the request before keeping an output.

    Premise changed in rc4: the automatic prefix names ``mirror``, so a mirrored
    run no longer meets the plain one's output -- it gets its own prefix and
    runs.  A different request under the *same* prefix needs an explicit
    ``--sample-id``, and that is refused as not reusable.
    """
    out = tmp_path / "out"
    base = ("run", "-m", "protein", "--pdb", os.path.join(pdb_dir, "1UBQ.pdb"), "-c", "A",
            "-o", out, "--no-triplets", "-q")
    assert cli(*base) == 0
    record = out / "1UBQ_A_protein.run.json"
    written = record.read_bytes()
    assert cli(*(base + ("--if-exists", "skip"))) == 0
    assert record.read_bytes() == written
    named = base + ("--sample-id", "1UBQ_A_protein")
    assert cli(*(named + ("--mirror", "--if-exists", "skip"))) == 2
    assert "cannot be reused" in capsys.readouterr().err
    assert record.read_bytes() == written
    assert cli(*(base + ("--mirror", "--if-exists", "skip"))) == 0
    assert run_json(out, "1UBQ_A_protein_mirror")["request"]["modifiers"] == {"mirror": True}


# ---------------------------------------------------------------------------
# Dry run
# ---------------------------------------------------------------------------


def _tree(path):
    out = []
    for root, dirs, files in os.walk(str(path)):
        out.extend(os.path.relpath(os.path.join(root, n), str(path)) for n in dirs + files)
    return sorted(out)


def test_a_dry_run_writes_downloads_and_computes_nothing(tmp_path, pdb_dir, capsys):
    local = tmp_path / "pdb"
    local.mkdir()
    shutil.copy(os.path.join(pdb_dir, "1UBQ.pdb"), str(local))
    before = _tree(local)
    manifest = write(tmp_path, "scop.csv", SCOP_HEADER + "1UBQ,A,alpha,0\n9ZZZ,A,beta,1\n")
    code = cli("batch", "--manifest", manifest, "--method", "protein",
               "--column", "label=group", "--keep-column", "group1",
               "--pdb-dir", local, "--download-missing", "-o", tmp_path / "out",
               "--dry-run", "-v")
    out = capsys.readouterr().out
    assert code == 0
    assert "1 present, 1 distinct would be downloaded" in out
    assert "1 would run, 1 would be downloaded then run" in out
    assert "download" in out and "9ZZZ_A_protein" in out
    assert "deferred" in out
    assert not (tmp_path / "out").exists()
    assert _tree(local) == before


def test_a_dry_run_reports_rows_that_would_fail(tmp_path, pdb_dir, capsys):
    manifest = write(tmp_path, "m.csv", "pdb_id,chain\n1UBQ,A\n9ZZZ,A\n")
    code = cli("batch", "--manifest", manifest, "--method", "protein",
               "--pdb-dir", pdb_dir, "-o", tmp_path / "out", "--dry-run")
    assert code == 3
    assert "no structure file for '9ZZZ'" in capsys.readouterr().out
    assert not (tmp_path / "out").exists()


def test_a_dry_run_verifies_existing_outputs(tmp_path, pdb_dir, capsys):
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets") == 0
    before = _tree(tmp_path / "out")
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets", "--if-exists", "skip",
                  "--dry-run") == 0
    assert "1 would be skipped (verified)" in capsys.readouterr().out
    assert _tree(tmp_path / "out") == before


# ---------------------------------------------------------------------------
# Staging, interruption and failures
# ---------------------------------------------------------------------------


def test_another_runs_staging_directory_is_never_touched(tmp_path, pdb_dir):
    """rc4 staging names: ``.biotsr-staging-<token>/PREFIX/``, the token naming the
    owning process and machine."""
    out = tmp_path / "out"
    # another machine's invocation (its host tag is not this host's): its owner
    # cannot be checked from here, so its directory is left alone -- even when it
    # is staging the very prefix this run writes
    other_host = "0" * 8 if output_mod.host_tag() != "0" * 8 else "f" * 8
    foreign = out / (".biotsr-staging-99999-%s-deadbeef" % other_host)
    (foreign / "ubq").mkdir(parents=True)
    (foreign / "ubq" / "ubq.keys").write_text("another invocation's work")
    assert output_mod.stale_leftovers(str(out)) == ([], [str(foreign)])
    assert _batch(tmp_path, pdb_dir, UBQ_ROW, "--no-triplets") == 0
    assert (foreign / "ubq" / "ubq.keys").read_text() == "another invocation's work"
    assert sorted(os.listdir(str(foreign))) == ["ubq"]
    assert sorted(n for n in os.listdir(str(out)) if n.startswith(".")) == [foreign.name]
    assert summary(out)[0]["status"] == "ok"


THREE = "sample_id,pdb_id,chain\nfirst,1UBQ,A\nsecond,1BRS,A\nthird,1BRS,D\n"


def test_an_interrupted_serial_batch_reports_every_row(tmp_path, pdb_dir, monkeypatch):
    real = engine.execute
    calls = []

    def interrupt_on_second(spec):
        calls.append(spec.sample_id)
        if len(calls) == 2:
            raise KeyboardInterrupt()
        return real(spec)

    monkeypatch.setattr(engine, "execute", interrupt_on_second)
    assert _batch(tmp_path, pdb_dir, THREE, "--no-triplets") == 130
    rows = summary(tmp_path / "out")
    assert [r["status"] for r in rows] == ["ok", "interrupted", "not_started"]
    # only the finished request's files and the summary: no staging, lock or
    # partial file of the interrupted one
    assert sorted(os.listdir(str(tmp_path / "out"))) == [
        "batch_summary.csv", "first.keys", "first.meta.json", "first.run.json"]
    # resuming completes the set and keeps the finished row
    monkeypatch.setattr(engine, "execute", real)
    assert _batch(tmp_path, pdb_dir, THREE, "--no-triplets", "--if-exists", "skip") == 0
    assert [r["status"] for r in summary(tmp_path / "out")] == ["skipped", "ok", "ok"]


def test_on_error_stop_leaves_later_rows_unstarted(tmp_path, pdb_dir):
    text = "sample_id,pdb_id,chain\nfirst,1UBQ,A\nbroken,1UBQ,Z\nlast,1BRS,A\n"
    # a row that fails while computing makes the batch partial: exit 3, summary written
    assert _batch(tmp_path, pdb_dir, text, "--no-triplets", "--on-error", "stop") == 3
    assert [r["status"] for r in summary(tmp_path / "out")] == ["ok", "failed", "not_started"]


def test_on_error_stop_refuses_a_row_that_cannot_be_planned_before_computing(tmp_path, pdb_dir):
    text = "sample_id,pdb_id,chain\nfirst,1UBQ,A\nmissing,9ZZZ,A\n"
    assert _batch(tmp_path, pdb_dir, text, "--no-triplets", "--on-error", "stop") == 2
    assert [r["status"] for r in summary(tmp_path / "out")] == ["not_started", "failed"]


def test_unclaimed_rows_are_classified_from_the_worker_markers(tmp_path):
    """Markers are keyed by the request's number (``_Request.index``), which in
    rc4 differs from its source line(s)."""
    rows = []
    for index, line in ((1, 2), (2, 3), (3, 4)):
        row = engine._Request(index=index, lines=[line], where="m.csv line %d" % line,
                              method="protein", participant="", values={}, metadata={},
                              chain_source="row",
                              record={"status": "failed", "message": ""})
        rows.append(row)
    markers = tmp_path / "markers"
    markers.mkdir()
    (markers / "1.started").write_text("")
    engine._write_marker(str(markers / "1.done"), {"status": "ok", "n_triangles": 7})
    (markers / "2.started").write_text("")
    # markers named by source line would be misread: line 4 is request 3
    (markers / "4.started").write_text("")
    engine._classify_unreceived(rows, str(markers), "stopped while running",
                                engine.STATUS_INTERRUPTED)
    assert [r.record["status"] for r in rows] == ["ok", "interrupted", "not_started"]
    assert rows[0].record["n_triangles"] == 7
    assert rows[1].record["message"] == "stopped while running"


# ---------------------------------------------------------------------------
# Memory
# ---------------------------------------------------------------------------


def test_the_batch_keeps_no_result_after_writing_it(tmp_path, pdb_dir):
    text = "sample_id,pdb_id,chain\n" + "".join(
        "s%d,1UBQ,A\n" % i for i in range(4))
    gc.collect()
    report = batch_api.run(write(tmp_path, "m.csv", text), str(tmp_path / "out"),
                           pdb_dir=pdb_dir, defaults={"method": "protein"})
    assert report.n_ok == 4
    for record in report.records:
        assert all(isinstance(v, (str, int, float, type(None))) for v in record.values())
    gc.collect()
    alive = [o for o in gc.get_objects() if isinstance(o, biotsr.TSRResult)]
    assert alive == []


# ---------------------------------------------------------------------------
# Parallel execution
# ---------------------------------------------------------------------------


EQUIVALENCE = """sample_id,pdb_id,method,chain,chain_2,aa_grouping,mirror,size_filter,sse,residues,residue,residue_name,entity_name,entity_id
plain,1UBQ,protein,A,,,,,,,,,,
grouped,1UBQ,protein,A,,true,,,,,,,,
mirror,1UBQ,protein,A,,,true,,,,,,,
size,1UBQ,protein,A,,,,12.0,,,,,,
sse,1UBQ,protein,A,,,,,true,,,,,
selection,1UBQ,aa_selection,A,,,,,,1-10;20,,,,
residue,1UBQ,amino_acid_atomic,A,,,,,,,15,LEU,,
every_residue,1UBQ,amino_acid_atomic,A,,,,,,,,,,
drug,3PTB,drug,A,,,,,,,,,BEN,1
carb,5WT9,carbohydrate,A,,,,,,,,,NAG,1
dna,1BNA,nucleotide,A,,,,,,,,,DG,2
cross,3PTB,drug_protein_cross,A,A,,,,,,,,BEN,1
"""


_BATCH_FILES = ("batch_summary.csv", "batch_failures.csv", "fetch_report.csv")


def _files(root):
    """Every calculation file in the shared output folder: ``{name: bytes}``.

    rc4: all calculations sit side by side in ``--output``; the batch-level
    reports and hidden names are not calculation files.
    """
    out = {}
    for name in sorted(os.listdir(str(root))):
        path = os.path.join(str(root), name)
        if name.startswith(".") or name in _BATCH_FILES or not os.path.isfile(path):
            continue
        with open(path, "rb") as fh:
            out[name] = fh.read()
    return out


def test_parallel_output_is_identical_to_serial_output(tmp_path):
    manifest = write(tmp_path, "m.csv", EQUIVALENCE)
    outputs = {}
    for jobs in ("1", "3"):
        out = tmp_path / ("jobs%s" % jobs)
        code = cli("batch", "--manifest", manifest, "--pdb-dir", STRUCTURES, "-o", out,
                   "--jobs", jobs, "-q")
        assert code == 0
        outputs[jobs] = out
    serial, par = _files(outputs["1"]), _files(outputs["3"])
    assert sorted(serial) == sorted(par)
    # every request has its completion record among the compared files
    prefixes = [r["sample_id"] for r in summary(outputs["1"])]
    assert len(prefixes) == 12
    assert {p + ".run.json" for p in prefixes} == {n for n in serial if n.endswith(".run.json")}
    assert sum(n.endswith(".keys") for n in serial) > len(prefixes)    # + sweep members
    for key in serial:
        if key.endswith(".run.json"):
            a, b = json.loads(serial[key]), json.loads(par[key])
            assert a.pop("argv") != b.pop("argv")       # the only difference: provenance
            assert a == b, key
        else:
            assert serial[key] == par[key], key
    strip = lambda rows: [{k: v for k, v in r.items() if k != "output"} for r in rows]
    assert strip(summary(outputs["1"])) == strip(summary(outputs["3"]))


def test_parallel_failures_are_reported_against_their_rows(tmp_path, pdb_dir):
    text = "sample_id,pdb_id,chain\na,1UBQ,A\nb,1UBQ,Z\nc,1BRS,A\nd,1BRS,D\n"
    assert _batch(tmp_path, pdb_dir, text, "--no-triplets", "--jobs", "2") == 3
    assert [r["status"] for r in summary(tmp_path / "out")] == ["ok", "failed", "ok", "ok"]
    leftovers = [n for n in os.listdir(str(tmp_path / "out")) if n.startswith(".")]
    assert leftovers == []


def test_parallel_execution_without_joblib_says_how_to_install_it(tmp_path, pdb_dir,
                                                                   monkeypatch, capsys):
    real_import = builtins.__import__

    def no_joblib(name, *args, **kwargs):
        if name == "joblib" or name.startswith("joblib."):
            raise ImportError("No module named 'joblib'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_joblib)
    text = "sample_id,pdb_id,chain\na,1UBQ,A\nb,1BRS,A\n"
    assert _batch(tmp_path, pdb_dir, text, "--no-triplets", "--jobs", "2") == 2
    err = capsys.readouterr().err
    assert "pip install biotsr" in err and 'pip install "joblib>=1.3"' in err
    assert "[parallel]" not in err     # an extra is no longer needed
    assert _batch(tmp_path, pdb_dir, text, "--no-triplets", "--jobs", "1",
                  "--if-exists", "overwrite") == 0


def test_worker_counts_follow_the_scheduler_allocation(monkeypatch):
    for name in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "PBS_NUM_PPN", "NCPUS",
                 "LSB_DJOB_NUMPROC", "NSLOTS"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    count, source = parallel.available_cpus()
    assert count == 1 and "Slurm" in source
    assert parallel.resolve_jobs("auto") == 1
    warnings = []
    assert parallel.resolve_jobs("3", warn=warnings.append) == 3
    assert "more workers than the 1 CPU(s)" in warnings[0]
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "")
    monkeypatch.setenv("NSLOTS", "2")
    assert parallel.available_cpus()[0] <= 2


@pytest.mark.parametrize("value", ["0", "-2", "many"])
def test_a_bad_worker_count_is_refused(value):
    with pytest.raises(CLIError):
        parallel.resolve_jobs(value)


# ---------------------------------------------------------------------------
# Retrieval from the command line
# ---------------------------------------------------------------------------


@pytest.fixture
def archive(monkeypatch):
    server = Stub()
    monkeypatch.setenv("BIOTSR_PDB_BASE_URL", server.base)
    yield server
    server.close()


def test_fetch_reads_only_the_id_column_of_a_dataset_sheet(tmp_path, archive, capsys):
    for pdb_id in ("1ABC", "1ABD"):
        archive.ok(pdb_id)
    manifest = write(tmp_path, "scop.csv", SCOP_HEADER + "1abc,A,alpha,0\n1ABD,X,beta,1\n"
                     "1ABC,B,alpha,0\n")
    pdb = tmp_path / "pdb"
    assert cli("fetch", "--manifest", manifest, "--pdb-dir", pdb, "-q") == 0
    assert sorted(n for n in os.listdir(str(pdb)) if n.endswith(".pdb")) == [
        "1ABC.pdb", "1ABD.pdb"]
    rows = list(csv.DictReader(open(str(pdb / "fetch_report.csv"))))
    assert [(r["pdb_id"], r["status"]) for r in rows] == [
        ("1ABC", "downloaded"), ("1ABD", "downloaded")]


def test_fetch_takes_a_mapped_id_column_and_reports_failures(tmp_path, archive):
    archive.ok("1ABC")
    manifest = write(tmp_path, "other.csv", "PDB,notes\n1ABC,x\n9ZZZ,y\n")
    pdb = tmp_path / "pdb"
    assert cli("fetch", "--manifest", manifest, "--column", "pdb_id=PDB",
               "--pdb-dir", pdb, "-q") == 3
    rows = list(csv.DictReader(open(str(pdb / "fetch_report.csv"))))
    assert [(r["pdb_id"], r["status"]) for r in rows] == [
        ("1ABC", "downloaded"), ("9ZZZ", "not_found")]


def test_fetch_dry_run_writes_nothing(tmp_path, archive, capsys):
    pdb = tmp_path / "pdb"
    assert cli("fetch", "1ABC", "--pdb-dir", pdb, "--dry-run") == 0
    assert not pdb.exists() and archive.hits == {}


def test_batch_downloads_what_is_missing_and_fails_what_cannot_be_fetched(tmp_path, archive,
                                                                           pdb_dir):
    with open(os.path.join(pdb_dir, "1UBQ.pdb"), "rb") as fh:
        archive.ok("2UBQ", fh.read().replace(b"1UBQ", b"2UBQ"))
    local = tmp_path / "pdb"
    local.mkdir()
    shutil.copy(os.path.join(pdb_dir, "1BRS.pdb"), str(local))
    manifest = write(tmp_path, "scop.csv",
                     SCOP_HEADER + "2UBQ,A,alpha,0\n1BRS,D,beta,1\n9ZZZ,A,beta,1\n")
    out = tmp_path / "out"
    code = cli("batch", "--manifest", manifest, "--method", "protein",
               "--column", "label=group", "--ignore-column", "group1",
               "--pdb-dir", local, "--download-missing", "--download-workers", "2",
               "-o", out, "--no-triplets", "-q")
    assert code == 3
    rows = summary(out)
    assert [r["status"] for r in rows] == ["ok", "ok", "failed"]
    assert rows[0]["n_triangles"] == "70300"
    assert "could not be downloaded" in rows[2]["message"]
    fetched = list(csv.DictReader(open(str(out / "fetch_report.csv"))))
    assert [(r["pdb_id"], r["status"]) for r in fetched] == [
        ("2UBQ", "downloaded"), ("9ZZZ", "not_found")]
    assert archive.hits.get("/download/1BRS.pdb.gz") is None     # present: not fetched


def test_download_missing_refuses_a_truncated_cached_file(tmp_path, archive):
    local = tmp_path / "pdb"
    local.mkdir()
    (local / "1ABC.pdb").write_bytes(_entry("1ABC").rsplit(b"END", 1)[0])
    manifest = write(tmp_path, "m.csv", "pdb_id,chain\n1ABC,A\n")
    out = tmp_path / "out"
    assert cli("batch", "--manifest", manifest, "--method", "protein", "--pdb-dir", local,
               "--download-missing", "-o", out, "-q") == 3
    row, = summary(out)
    assert "END record" in row["message"] and "--refresh" in row["message"]


# ---------------------------------------------------------------------------
# The Python entry point
# ---------------------------------------------------------------------------


def test_the_python_entry_point_matches_the_command_line(tmp_path, pdb_dir):
    manifest = write(tmp_path, "scop.csv", SCOP)
    rows = batch_api.read_manifest(manifest, columns={"label": "group"},
                                   keep_columns=["group1"], defaults={"method": "protein"})
    assert rows[0]["label"] == "alpha" and rows[0]["metadata"] == {"group1": "0"}
    assert rows[0]["line"] == 2 and rows[0]["chain_source"] == "row"
    report = batch_api.run(manifest, str(tmp_path / "py"), pdb_dir=pdb_dir,
                           defaults={"method": "protein"}, columns={"label": "group"},
                           keep_columns=["group1"], triplets=False)
    assert (report.n_ok, report.n_failed) == (3, 0)
    assert cli("batch", "--manifest", manifest, "--method", "protein",
               "--column", "label=group", "--keep-column", "group1", "--pdb-dir", pdb_dir,
               "-o", tmp_path / "cli", "--no-triplets", "-q") == 0
    assert _files(tmp_path / "py").keys() == _files(tmp_path / "cli").keys() == {
        "%s_protein.%s" % (sample, suffix)
        for sample in ("1UBQ_A", "1BRS_A", "1BRS_D")
        for suffix in ("keys", "meta.json", "run.json")}
    with pytest.raises(batch_api.BatchError):
        batch_api.run(manifest, str(tmp_path / "x"), pdb_dir=pdb_dir)   # group unmapped


def test_analyze_reads_samples_back_including_verified_skips(tmp_path, pdb_dir):
    manifest = write(tmp_path, "m.csv", "sample_id,pdb_id,chain,label\n"
                     "ubq,1UBQ,A,x\nbarnase,1BRS,A,y\nbarstar,1BRS,D,y\n")
    out = tmp_path / "out"
    base = ("analyze", "common-keys", "--manifest", manifest, "--method", "protein",
            "--pdb-dir", pdb_dir, "-o", out, "-q")
    assert cli(*base) == 0
    first = (out / "analysis" / "common_keys" / "all_samples.csv").read_text()
    assert cli(*(base + ("--if-exists", "skip"))) == 0
    assert [r["status"] for r in summary(out)] == ["skipped"] * 3
    assert (out / "analysis" / "common_keys" / "all_samples.csv").read_text() == first


# ---------------------------------------------------------------------------
# Streaming triplets
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("argv,call", [
    (("-m", "protein", "--pdb", "1UBQ", "-c", "A", "--mirror"),
     lambda p: biotsr.analyze(p("1UBQ"), method="protein", chain="A", mirror=True)),
    (("-m", "drug_protein_cross", "--pdb", "3PTB", "-c", "A", "--chain-2", "A",
      "--entity-name", "BEN", "--entity-id", "1"),
     lambda p: biotsr.analyze_cross(p("3PTB"), method="drug_protein_cross", chain="A",
                                    chain_2="A", entity_name="BEN", entity_id="1")),
])
def test_streamed_triplets_are_the_bytes_the_serializer_writes(tmp_path, argv, call):
    """The CLI streams native triplets to disk as they are computed; the files
    must be exactly those of saving an in-memory result."""
    from biotsr.io import native

    argv = [a if a not in ("1UBQ", "3PTB") else structure(a) for a in argv]
    out = tmp_path / "cli"
    assert cli("run", *argv, "-o", out, "-q") == 0
    reference = tmp_path / "ref"
    native.save(call(structure), str(reference))
    # rc4: one shared folder; the API's default prefix is the command line's
    prefix = native.default_prefix(call(structure))
    expected = sorted(prefix + s for s in (".keys", ".meta.json", ".triplets"))
    assert sorted(os.listdir(str(reference))) == expected
    assert sorted(os.listdir(str(out))) == sorted(expected + [prefix + ".run.json"])
    for name in expected:
        with open(os.path.join(str(out), name), "rb") as a, \
                open(str(reference / name), "rb") as b:
            assert a.read() == b.read(), name
    assert run_json(out, prefix)["results"][0]["files"]["triplets"] == prefix + ".triplets"


def test_a_failed_streaming_run_leaves_nothing_behind(tmp_path, pdb_dir):
    out = tmp_path / "never"
    assert cli("run", "-m", "protein", "--pdb", os.path.join(pdb_dir, "1UBQ.pdb"),
               "-c", "Z", "-o", out, "-q") == 2
    assert not out.exists()
