"""Several invocations writing into one output folder at the same time.

Each calculation is staged in a directory unique to its invocation and
installed under a per-prefix lock that refuses to replace another run's
calculation, and an interrupted invocation removes only its own staging
directory.  These tests run real processes side by side and interrupt one of
them.
"""

from __future__ import annotations

import csv
import json
import os
import signal
import subprocess
import sys
import time

import pytest

from biotsr.cli import batch as engine
from biotsr.cli import output as _output
from biotsr.cli.errors import CLIError

from ..conftest import STRUCTURES, WORKSPACE

pytestmark = pytest.mark.skipif(not hasattr(signal, "SIGINT") or os.name != "posix",
                                reason="POSIX signals")

#: 3EJ1 chain A: about four million triangles, long enough to interrupt.
SLOW = os.path.join(STRUCTURES, "3EJ1.pdb")
FAST = os.path.join(STRUCTURES, "1UBQ.pdb")


def _env():
    env = dict(os.environ)
    env["PYTHONPATH"] = os.path.join(WORKSPACE, "src") + os.pathsep + env.get("PYTHONPATH", "")
    return env


def _biotsr(*argv, **kwargs):
    return [sys.executable, "-m", "biotsr"] + [str(a) for a in argv]


def _default_sigint():
    """In the child: take Ctrl-C normally even if this test process ignores it.

    A process started in the background by a non-interactive shell inherits
    SIGINT as ignored, and Python then installs no KeyboardInterrupt handler.
    """
    signal.signal(signal.SIGINT, signal.SIG_DFL)


def _wait_for_staging(out, prefix, timeout=60.0):
    """Until ``out/.biotsr-staging-<token>/<prefix>`` exists."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        if os.path.isdir(out):
            for name in os.listdir(out):
                if name.startswith(".biotsr-staging-") and os.path.isdir(
                        os.path.join(out, name, prefix)):
                    return
        time.sleep(0.05)
    raise AssertionError("no staging directory for %s appeared in %s" % (prefix, out))


def _files(folder, prefix):
    """Every file of one calculation, by name; its record without the command line."""
    out = {}
    for name in sorted(os.listdir(folder)):
        if not name.startswith(prefix + "."):
            continue
        with open(os.path.join(folder, name), "rb") as fh:
            data = fh.read()
        if name.endswith(".run.json"):
            payload = json.loads(data)
            payload.pop("argv")
            data = json.dumps(payload, sort_keys=True).encode()
        out[name] = data
    return out


def _calculation(prefix, triplets=True):
    names = [prefix + ".keys", prefix + ".meta.json", prefix + ".run.json"]
    return sorted(names + ([prefix + ".triplets"] if triplets else []))


def _interrupt(proc):
    proc.send_signal(signal.SIGINT)
    return proc.wait(timeout=120)


def test_an_interrupted_run_leaves_a_run_sharing_its_output_intact(tmp_path):
    out = str(tmp_path / "out")
    slow = subprocess.Popen(_biotsr("run", "-m", "protein", "--pdb", SLOW, "-c", "A",
                                    "-o", out, "-q"),
                            env=_env(), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            preexec_fn=_default_sigint)
    try:
        _wait_for_staging(out, "3EJ1_A_protein")
        fast = subprocess.run(_biotsr("run", "-m", "protein", "--pdb", FAST, "-c", "A",
                                      "-o", out, "-q"),
                              env=_env(), capture_output=True, text=True, timeout=300)
        assert fast.returncode == 0, fast.stderr
        assert slow.poll() is None, "the slow run finished before it could be interrupted"
    finally:
        code = _interrupt(slow)
    assert code == 130
    # Only the finished calculation is left: no staging, no lock, nothing partial.
    assert sorted(os.listdir(out)) == _calculation("1UBQ_A_protein")
    # The surviving sample is exactly what the same run writes on its own ...
    ref = str(tmp_path / "ref")
    alone = subprocess.run(_biotsr("run", "-m", "protein", "--pdb", FAST, "-c", "A",
                                   "-o", ref, "-q"),
                           env=_env(), capture_output=True, text=True, timeout=300)
    assert alone.returncode == 0, alone.stderr
    assert _files(out, "1UBQ_A_protein") == _files(ref, "1UBQ_A_protein")
    # ... and it passes the complete-and-same-request check a resume applies.
    again = subprocess.run(_biotsr("run", "-m", "protein", "--pdb", FAST, "-c", "A",
                                   "-o", out, "--if-exists", "skip", "-q"),
                           env=_env(), capture_output=True, text=True, timeout=300)
    assert again.returncode == 0, again.stderr


def test_an_interrupted_batch_removes_only_its_own_files(tmp_path):
    out = str(tmp_path / "out")
    manifest = tmp_path / "m.csv"
    manifest.write_text("sample_id,pdb_path,method,chain\n"
                        "slow,%s,protein,A\nlater,%s,protein,A\n" % (SLOW, FAST))
    batch = subprocess.Popen(_biotsr("batch", "--manifest", manifest, "-o", out, "-q"),
                             env=_env(), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             preexec_fn=_default_sigint)
    try:
        _wait_for_staging(out, "slow")
        other = subprocess.run(_biotsr("run", "-m", "protein", "--pdb", FAST, "-c", "A",
                                       "-o", out, "--sample-id", "other", "-q"),
                               env=_env(), capture_output=True, text=True, timeout=300)
        assert other.returncode == 0, other.stderr
        assert batch.poll() is None
    finally:
        code = _interrupt(batch)
    assert code == 130
    # lock released, staging removed, nothing of "slow" left
    assert sorted(os.listdir(out)) == ["batch_summary.csv"] + _calculation("other")
    with open(os.path.join(out, "batch_summary.csv"), newline="") as fh:
        statuses = [(r["sample_id"], r["status"]) for r in csv.DictReader(fh)]
    assert statuses == [("slow", "interrupted"), ("later", "not_started")]
    payload = json.load(open(os.path.join(out, "other.run.json")))
    assert payload["results"][0]["counts"]["n_triangles"] == 70300


# ---------------------------------------------------------------------------
# One batch per output directory
# ---------------------------------------------------------------------------


def _manifest(tmp_path):
    path = tmp_path / "m.csv"
    path.write_text("sample_id,pdb_path,method,chain\nubq,%s,protein,A\n" % FAST)
    return str(path)


def test_a_second_batch_into_a_busy_output_is_refused(tmp_path, capsys):
    out = tmp_path / "out"
    out.mkdir()
    lock = out / engine._OutputLock.NAME
    lock.write_text(json.dumps({"pid": os.getpid(), "host": engine._host(),
                                "token": "someone-else"}))
    code = engine_cli("batch", "--manifest", _manifest(tmp_path), "-o", out,
                      "--no-triplets", "-q")
    assert code == 2
    err = capsys.readouterr().err
    assert "another batch is writing into" in err and "--output" in err
    assert lock.exists() and not (out / "ubq.keys").exists()
    # A dry run takes no lock and leaves the other batch's alone.
    assert engine_cli("batch", "--manifest", _manifest(tmp_path), "-o", out,
                      "--no-triplets", "--dry-run", "-q") == 0
    assert json.loads(lock.read_text())["token"] == "someone-else"


def test_a_lock_left_by_a_dead_process_is_taken_over(tmp_path):
    finished = subprocess.run([sys.executable, "-c", "import os; print(os.getpid())"],
                              capture_output=True, text=True)
    dead_pid = int(finished.stdout)
    out = tmp_path / "out"
    out.mkdir()
    (out / engine._OutputLock.NAME).write_text(
        json.dumps({"pid": dead_pid, "host": engine._host(), "token": "stale"}))
    assert engine_cli("batch", "--manifest", _manifest(tmp_path), "-o", out,
                      "--no-triplets", "-q") == 0
    assert sorted(os.listdir(str(out))) == ["batch_summary.csv"] + _calculation(
        "ubq", triplets=False)


def engine_cli(*argv):
    from biotsr.cli import main

    return main([str(a) for a in argv])


# ---------------------------------------------------------------------------
# Two invocations finishing the same sample
# ---------------------------------------------------------------------------


def _token(suffix):
    """An invocation token in the ``<pid>-<hosttag8>-<hex8>`` form.

    Both "invocations" carry this process's pid and host tag, so each sees the
    other's staging directory as owned by a live process and leaves it alone.
    """
    return "%d-%s-%s" % (os.getpid(), _output.host_tag(), suffix * 8)


def _staged(root, token, text, if_exists=_output.IF_EXISTS_FAIL):
    writer = _output.CalculationWriter(str(root), "s", if_exists, token=token)
    writer.__enter__()
    with open(os.path.join(writer.path, "s.keys"), "w") as fh:
        fh.write(text)
    writer.payload = {"schema": _output.REQUEST_SCHEMA, "prefix": "s",
                      "results": [{"name": "s", "files": {"keys": "s.keys"},
                                   "file_bytes": {"s.keys": len(text)}}]}
    return writer


def test_the_second_run_to_finish_a_calculation_never_replaces_the_first(tmp_path):
    first = _staged(tmp_path, _token("a"), "first")
    second = _staged(tmp_path, _token("b"), "second")
    first.__exit__(None, None, None)
    with pytest.raises(CLIError, match="appeared in .* while it was being computed"):
        second.__exit__(None, None, None)
    assert sorted(os.listdir(str(tmp_path))) == ["s.keys", "s.run.json"]
    assert (tmp_path / "s.keys").read_text() == "first"


def test_an_installation_in_progress_is_never_interleaved_with_another(tmp_path):
    """The per-prefix lock: a second installer meeting a live owner's lock
    refuses rather than mixing its files with the first one's."""
    second = _staged(tmp_path, _token("b"), "second")
    lock = tmp_path / ".biotsr-lock-s"
    lock.write_text(_token("a"))                     # a live invocation installing s
    with pytest.raises(CLIError, match="another run is installing s"):
        second.__exit__(None, None, None)
    assert lock.read_text() == _token("a")
    assert not (tmp_path / "s.keys").exists()


def test_only_an_explicit_overwrite_replaces_a_finished_calculation(tmp_path):
    first = _staged(tmp_path, _token("a"), "first")
    second = _staged(tmp_path, _token("b"), "second", if_exists=_output.IF_EXISTS_OVERWRITE)
    first.__exit__(None, None, None)
    second.__exit__(None, None, None)
    assert sorted(os.listdir(str(tmp_path))) == ["s.keys", "s.run.json"]
    assert (tmp_path / "s.keys").read_text() == "second"
