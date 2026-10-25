"""Worker death, signals and leftovers: every row ends with a truthful status.

Regression tests for two rc1 defects:

* a worker process killed during ``--jobs N`` (the out-of-memory killer, a
  signal) lost the whole batch summary -- exit 1, no ``batch_summary.csv``,
  finished rows unreported;
* SIGTERM/SIGINT outside the compute phase (while structures downloaded, or
  during ``biotsr run``) exited with no summary and could leave a stale lock or
  a staging directory, which ``analyze --from-results`` then read as a sample.

Since rc4 every calculation is installed into one shared ``--output`` folder
under its prefix (``PREFIX.keys``, ``PREFIX.meta.json``, ``PREFIX.triplets``,
``PREFIX.run.json``), staged first in ``.biotsr-staging-<token>/PREFIX/``.  The
other hidden names an invocation can leave are ``.biotsr-tmp-<token>-*``,
``.biotsr-batch-<token>-*`` (worker markers) and ``.biotsr-lock-PREFIX``, whose
content is the owner's token.

Every scenario runs ``biotsr`` in its own subprocess, so a real signal reaches a
real process.  ``SIGKILL`` cannot be handled by any program; what it leaves
behind is recognised and cleaned up by the next run on the same machine.
"""

from __future__ import annotations

import csv
import json
import os
import signal
import subprocess
import sys
import textwrap
import time

import pytest

from biotsr.cli.output import host_tag, stale_leftovers

from ..conftest import SRC, structure
from ..support.stub_server import Stub, _entry, _gz

pytestmark = pytest.mark.skipif(
    sys.platform == "win32", reason="POSIX signals"
)

#: Installed through PYTHONPATH as ``sitecustomize`` in the subprocess and in its
#: joblib workers: the worker that picks up the named sample kills itself with
#: SIGKILL, exactly as the out-of-memory killer would.
_KILLER = textwrap.dedent('''
    import importlib.abc, importlib.util, os, signal, sys

    _TARGET = os.environ.get("BIOTSR_TEST_KILL_SAMPLE")

    class _Hook(importlib.abc.MetaPathFinder):
        def find_spec(self, name, path, target=None):
            if name != "biotsr.cli.batch":
                return None
            sys.meta_path.remove(self)
            spec = importlib.util.find_spec(name)
            original = spec.loader.exec_module

            def exec_module(module):
                original(module)
                real = module.run_task

                def run_task(task):
                    # Only worker processes run tasks when --jobs > 1. Die the
                    # way the out-of-memory killer would: after the row has
                    # started, before it finishes. (rc4: a task names its
                    # output prefix, which is the row's sample_id.)
                    if task["prefix"] == _TARGET:
                        marker = os.path.join(task["markers"],
                                              "%d.started" % task["index"])
                        open(marker, "w").close()
                        os.kill(os.getpid(), signal.SIGKILL)
                    return real(task)

                module.run_task = run_task

            spec.loader.exec_module = exec_module
            return spec

    if _TARGET:
        sys.meta_path.insert(0, _Hook())
''')


def _env(extra=None, pythonpath=()):
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(list(pythonpath) + [SRC])
    env.pop("BIOTSR_PDB_BASE_URL", None)
    env.update(extra or {})
    return env


def _biotsr(args, **kwargs):
    return [sys.executable, "-m", "biotsr"] + list(args)


def _manifest(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "pdb_path", "method", "chain"])
        writer.writerows(rows)
    return str(path)


def _summary(out):
    with open(os.path.join(out, "batch_summary.csv")) as fh:
        return {r["sample_id"]: r for r in csv.DictReader(fh)}


_SUFFIXES = (".keys", ".meta.json", ".triplets", ".run.json")


def _calculation_files(out, prefix):
    """The files of one calculation present in the shared output folder."""
    return sorted(n for n in os.listdir(out) if n in {prefix + s for s in _SUFFIXES})


def _wait_for(predicate, timeout=60.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if predicate():
            return True
        time.sleep(0.05)
    return False


# ---------------------------------------------------------------------------
# A worker process dies
# ---------------------------------------------------------------------------


def _default_sigint():
    """In the child: take Ctrl-C normally even if this test process ignores it.

    A process started in the background by a non-interactive shell -- a suite
    run under ``nohup``/``&`` or by a scheduler -- inherits SIGINT as ignored,
    and Python then installs no KeyboardInterrupt handler in its children.
    """
    signal.signal(signal.SIGINT, signal.SIG_DFL)


def test_a_killed_worker_keeps_every_finished_row_and_reports_the_rest(tmp_path):
    hook = tmp_path / "hook"
    hook.mkdir()
    (hook / "sitecustomize.py").write_text(_KILLER)
    rows = [["s%d" % i, structure("1UBQ"), "protein", "A"] for i in range(1, 7)]
    manifest = _manifest(tmp_path / "m.csv", rows)
    out = tmp_path / "out"
    proc = subprocess.Popen(
        _biotsr(["batch", "--manifest", manifest, "--output", str(out),
                 "--no-triplets", "--jobs", "2", "-q"]),
        env=_env({"BIOTSR_TEST_KILL_SAMPLE": "s3"}, pythonpath=[str(hook)]),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    stdout, stderr = proc.communicate(timeout=300)
    assert proc.returncode == 3, stderr

    summary = _summary(out)
    assert set(summary) == {"s1", "s2", "s3", "s4", "s5", "s6"}
    assert summary["s3"]["status"] == "failed"
    assert "worker pool failed" in summary["s3"]["message"]
    statuses = {r["status"] for r in summary.values()}
    assert statuses <= {"ok", "failed", "not_started"}
    for sample_id, row in summary.items():
        if row["status"] == "ok":
            # Nothing is guessed: an ok row has a complete, verifiable output.
            assert _calculation_files(str(out), sample_id) == [
                sample_id + ".keys", sample_id + ".meta.json", sample_id + ".run.json"]
            with open(os.path.join(out, sample_id + ".run.json")) as fh:
                assert json.load(fh)["complete"] is True
        else:
            assert _calculation_files(str(out), sample_id) == []
    # The killed worker's staging directory belonged to this invocation and was
    # removed by the coordinator; nothing hidden is left behind.
    assert [n for n in os.listdir(out) if n.startswith(".")] == []

    # Resume: finished rows are verified and kept, the rest are computed.
    resumed = subprocess.run(
        _biotsr(["batch", "--manifest", manifest, "--output", str(out),
                 "--no-triplets", "--jobs", "2", "--if-exists", "skip", "-q"]),
        env=_env(), capture_output=True, text=True, timeout=300,
    )
    assert resumed.returncode == 0, resumed.stderr
    after = _summary(out)
    assert {r["status"] for r in after.values()} <= {"ok", "skipped"}
    for sample_id, row in summary.items():
        if row["status"] == "ok":
            assert after[sample_id]["status"] == "skipped"


def test_more_workers_than_rows_starts_only_as_many_as_needed(tmp_path):
    manifest = _manifest(tmp_path / "m.csv", [
        ["a", structure("1UBQ"), "protein", "A"],
        ["b", structure("1UBQ"), "protein", "A"],
    ])
    done = subprocess.run(
        _biotsr(["batch", "--manifest", manifest, "--output", str(tmp_path / "out"),
                 "--no-triplets", "--jobs", "8"]),
        env=_env(), capture_output=True, text=True, timeout=300,
    )
    assert done.returncode == 0, done.stderr
    assert "workers    2 processes" in done.stderr


# ---------------------------------------------------------------------------
# Signals outside the compute phase
# ---------------------------------------------------------------------------


@pytest.fixture
def slow_archive():
    server = Stub()
    yield server
    server.close()


@pytest.mark.parametrize("signum", [signal.SIGTERM, signal.SIGINT])
def test_a_signal_while_downloading_writes_both_reports(tmp_path, slow_archive, signum):
    for pdb_id in ("1AAA", "1AAB", "1AAC"):
        slow_archive.serve(pdb_id, {"status": 200, "body": _gz(_entry(pdb_id)),
                                    "delay": 30})
    manifest = tmp_path / "m.csv"
    manifest.write_text("pdb_id,method,chain\n1AAA,protein,A\n1AAB,protein,A\n"
                        "1AAC,protein,A\n")
    out, pdb_dir = tmp_path / "out", tmp_path / "pdb"
    proc = subprocess.Popen(
        _biotsr(["batch", "--manifest", str(manifest), "--pdb-dir", str(pdb_dir),
                 "--download-missing", "--download-workers", "1", "--output",
                 str(out), "--no-triplets", "--timeout", "30"]),
        env=_env({"BIOTSR_PDB_BASE_URL": slow_archive.base}),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        preexec_fn=_default_sigint,
    )
    assert _wait_for(lambda: slow_archive.active > 0), "no download started"
    started = time.time()
    proc.send_signal(signum)
    stdout, stderr = proc.communicate(timeout=60)
    assert proc.returncode == 130, stderr
    # No retry starts after the signal, queued ids never start, and the request
    # in flight is stopped at once -- not at its 30-second --timeout, nor at
    # the server's 30-second delay (test_download_cancellation.py).
    assert time.time() - started < 10
    summary = _summary(out)
    assert {r["status"] for r in summary.values()} == {"not_started"}
    assert all("downloading" in r["message"] for r in summary.values())
    with open(out / "fetch_report.csv") as fh:
        fetched = {r["pdb_id"]: r["status"] for r in csv.DictReader(fh)}
    # One download worker: the first id is in flight, the others are queued
    # and never start (test_download_cancellation.py).
    assert fetched == {"1AAA": "interrupted", "1AAB": "not_started",
                       "1AAC": "not_started"}
    assert not (out / ".biotsr-batch.lock").exists()
    assert not pdb_dir.exists() or not [
        n for n in os.listdir(pdb_dir) if n.endswith(".pdb")]


def test_fetch_interrupted_writes_its_report(tmp_path, slow_archive):
    slow_archive.serve("1AAA", {"status": 200, "body": _gz(_entry("1AAA")),
                                "delay": 30})
    pdb_dir = tmp_path / "pdb"
    proc = subprocess.Popen(
        _biotsr(["fetch", "1AAA", "1AAB", "--pdb-dir", str(pdb_dir),
                 "--download-workers", "1", "--timeout", "30"]),
        env=_env({"BIOTSR_PDB_BASE_URL": slow_archive.base}),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        preexec_fn=_default_sigint,
    )
    assert _wait_for(lambda: slow_archive.active > 0)
    started = time.time()
    proc.send_signal(signal.SIGTERM)
    _out, stderr = proc.communicate(timeout=60)
    assert time.time() - started < 10
    assert proc.returncode == 130, stderr
    with open(pdb_dir / "fetch_report.csv") as fh:
        rows = {r["pdb_id"]: r["status"] for r in csv.DictReader(fh)}
    assert rows == {"1AAA": "interrupted", "1AAB": "not_started"}


def test_sigterm_during_run_leaves_no_staging_directory(tmp_path):
    out = tmp_path / "out"
    proc = subprocess.Popen(
        _biotsr(["run", "--method", "protein", "--pdb", structure("1TUP"),
                 "--chain", "A", "--output", str(out)]),
        env=_env(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    assert _wait_for(lambda: out.exists() and any(
        n.startswith(".") for n in os.listdir(out))), "no staging directory appeared"
    proc.send_signal(signal.SIGTERM)
    _out, stderr = proc.communicate(timeout=120)
    assert proc.returncode == 130, stderr
    assert not out.exists() or os.listdir(out) == []


# ---------------------------------------------------------------------------
# What an uncatchable termination leaves behind
# ---------------------------------------------------------------------------


def _dead_pid():
    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    proc.wait()
    return proc.pid


def _token(pid, host, tail):
    """An invocation token as rc4 writes it: ``<pid>-<host tag>-<8 hex>``."""
    return "%d-%s-%s" % (pid, host, tail)


def test_leftovers_of_a_dead_process_are_removed_and_live_ones_are_not(tmp_path):
    """Every kind of hidden name an invocation creates is recognised by its owner.

    rc4 names: a staging directory ``.biotsr-staging-<token>``, a temporary file
    ``.biotsr-tmp-<token>-*``, a worker-marker directory ``.biotsr-batch-<token>-*``
    and a per-prefix lock ``.biotsr-lock-PREFIX`` whose content is the token.
    """
    out = tmp_path / "out"
    out.mkdir()
    other = "ffffffff" if host_tag() != "ffffffff" else "00000000"
    dead_token = _token(_dead_pid(), host_tag(), "0000abcd")
    live_token = _token(os.getpid(), host_tag(), "0000abce")
    foreign_token = _token(_dead_pid(), other, "0000abcf")
    dead = [".biotsr-staging-" + dead_token, ".biotsr-batch-%s-q1w2e3" % dead_token]
    dead_files = [".biotsr-tmp-%s-old.run.json" % dead_token, ".biotsr-lock-old"]
    live = [".biotsr-staging-" + live_token]
    live_files = [".biotsr-lock-current"]
    foreign = [".biotsr-staging-" + foreign_token]
    foreign_files = [".biotsr-lock-far"]
    unrelated = [".mine", ".biotsr-staging-not-a-token"]
    for name in dead + live + foreign + unrelated:
        (out / name / "old").mkdir(parents=True)
        (out / name / "old" / "old.keys").write_text("x")
    (out / dead_files[0]).write_text("{}")
    for name, token in ((dead_files[1], dead_token), (live_files[0], live_token),
                        (foreign_files[0], foreign_token)):
        (out / name).write_text(token)
    stale, foreign_found = stale_leftovers(str(out))
    assert [os.path.basename(p) for p in stale] == sorted(dead + dead_files)
    assert [os.path.basename(p) for p in foreign_found] == sorted(foreign + foreign_files)

    manifest = _manifest(tmp_path / "m.csv", [["a", structure("1UBQ"), "protein", "A"]])
    done = subprocess.run(
        _biotsr(["batch", "--manifest", manifest, "--output", str(out),
                 "--no-triplets"]),
        env=_env(), capture_output=True, text=True, timeout=300,
    )
    assert done.returncode == 0, done.stderr
    names = set(os.listdir(out))
    assert not names & set(dead + dead_files)
    assert set(live + live_files + foreign + foreign_files + unrelated) <= names
    for name in live + foreign + unrelated:
        assert (out / name / "old" / "old.keys").read_text() == "x"
    assert (out / live_files[0]).read_text() == live_token
    assert (out / foreign_files[0]).read_text() == foreign_token
    assert "another machine" in done.stderr
    assert _calculation_files(str(out), "a") == ["a.keys", "a.meta.json", "a.run.json"]


def test_a_dry_run_only_reports_leftovers(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    dead = ".biotsr-staging-" + _token(_dead_pid(), host_tag(), "0000abcd")
    (out / dead / "a").mkdir(parents=True)
    manifest = _manifest(tmp_path / "m.csv", [["a", structure("1UBQ"), "protein", "A"]])
    done = subprocess.run(
        _biotsr(["batch", "--manifest", manifest, "--output", str(out), "--dry-run"]),
        env=_env(), capture_output=True, text=True, timeout=300,
    )
    assert done.returncode == 0, done.stderr
    assert "1 staging name(s) of ended runs would be removed" in done.stdout
    assert (out / dead / "a").is_dir()
    assert sorted(os.listdir(out)) == [dead]


def test_analysis_never_reads_a_staging_directory_or_an_incomplete_result(tmp_path):
    out = tmp_path / "out"
    manifest = _manifest(tmp_path / "m.csv", [
        ["a", structure("1UBQ"), "protein", "A"],
        ["b", structure("1BRS"), "protein", "D"],
    ])
    done = subprocess.run(
        _biotsr(["batch", "--manifest", manifest, "--output", str(out),
                 "--no-triplets", "-q"]),
        env=_env(), capture_output=True, text=True, timeout=300,
    )
    assert done.returncode == 0, done.stderr
    # A staging copy of a calculation, as an interrupted run would leave it
    # (rc4: .biotsr-staging-<token>/PREFIX/), including its completion record.
    staging = out / (".biotsr-staging-" + _token(_dead_pid(), host_tag(), "0000abcd")) / "a"
    staging.mkdir(parents=True)
    for name in _calculation_files(str(out), "a"):
        (staging / name).write_bytes((out / name).read_bytes())
    assert sorted(os.listdir(staging)) == ["a.keys", "a.meta.json", "a.run.json"]
    analysed = subprocess.run(
        _biotsr(["analyze", "common-keys", "--from-results", str(out),
                 "--output", str(tmp_path / "analysis")]),
        env=_env(), capture_output=True, text=True, timeout=300,
    )
    assert analysed.returncode == 0, analysed.stderr
    report = json.loads((tmp_path / "analysis" / "analysis" / "common_keys"
                         / "all_samples.json").read_text())
    assert [s["sample_id"] for s in report["samples"]] == ["a", "b"]

    def analyse(name):
        return subprocess.run(
            _biotsr(["analyze", "common-keys", "--from-results", str(out),
                     "--output", str(tmp_path / name)]),
            env=_env(), capture_output=True, text=True, timeout=300,
        )

    # A key file cut short after its completion record was written.
    keys = out / "b.keys"
    intact = keys.read_bytes()
    keys.write_bytes(intact[:-40])
    refused = analyse("analysis2")
    assert refused.returncode == 2
    assert "not complete" in refused.stderr and "b.keys" in refused.stderr
    assert not (tmp_path / "analysis2" / "analysis").exists()

    # rc4: an intent record left by an interrupted installation (complete:false)
    # is refused too, even with every file intact.
    keys.write_bytes(intact)
    record = out / "b.run.json"
    payload = json.loads(record.read_text())
    payload["complete"] = False
    record.write_text(json.dumps(payload))
    refused = analyse("analysis3")
    assert refused.returncode == 2
    assert "not complete" in refused.stderr and "b.run.json" in refused.stderr
    assert "interrupted" in refused.stderr

    # A key file that no completion record lists is never read as a sample.
    payload["complete"] = True
    record.write_text(json.dumps(payload))
    (out / "c.keys").write_bytes(intact)
    refused = analyse("analysis4")
    assert refused.returncode == 2
    assert "c.keys" in refused.stderr and "no completion record lists it" in refused.stderr
