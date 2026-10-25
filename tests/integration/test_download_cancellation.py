"""What cancelling a download actually guarantees, phase by phase.

rc2 documented that "an in-flight read still ends at --timeout".  That is not
a bound: ``--timeout`` limits each wait, so a response that keeps delivering
small pieces never times out, and the whole response was read in one blocking
call that a cancellation could not interrupt.  These tests pin what is true now
(``biotsr.retrieval``, "Timeouts and cancellation"):

* a connection being established is abandoned within a fraction of a second;
* a stalled response, and one still trickling in, stop at once;
* a backoff wait ends at once, and no retry follows;
* queued identifiers never start;
* the report names every identifier truthfully, and nothing partial, and
  nothing at all after the cancellation, reaches the structure directory.

Each scenario runs ``biotsr fetch`` in its own process and sends it a real
signal.  ``--timeout 30`` throughout, so finishing in a few seconds cannot be
the timeout doing the work.
"""

from __future__ import annotations

import concurrent.futures
import csv
import hashlib
import os
import signal
import socket
import subprocess
import sys
import time

import pytest

import biotsr

from ..conftest import SRC
from ..support.stub_server import Stub, _entry, _gz

pytestmark = pytest.mark.skipif(sys.platform == "win32", reason="POSIX signals")

#: Seconds a cancelled fetch may take to exit: interpreter teardown and report
#: writing, never a network wait (every scenario's network wait is >= 30 s).
PROMPT = 5.0


@pytest.fixture
def archive():
    server = Stub()
    yield server
    server.close()


def _default_sigint():
    """In the child: take Ctrl-C normally even if this test process ignores it.

    A process started in the background by a non-interactive shell -- a suite
    run under ``nohup``/``&`` or by a scheduler -- inherits SIGINT as ignored,
    and Python then installs no KeyboardInterrupt handler in its children.
    """
    signal.signal(signal.SIGINT, signal.SIG_DFL)


def test_the_signal_setup_changes_only_the_child():
    """``_default_sigint`` runs in the child before it starts; the test process
    keeps its own disposition -- here deliberately ``SIG_IGN``, as a suite run
    in the background inherits -- and gets back exactly what it had."""
    before = signal.getsignal(signal.SIGINT)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        probe = subprocess.run(
            [sys.executable, "-c", "import signal; print(signal.getsignal("
             "signal.SIGINT) is signal.default_int_handler)"],
            stdout=subprocess.PIPE, text=True, preexec_fn=_default_sigint)
        assert probe.stdout.strip() == "True"
        assert signal.getsignal(signal.SIGINT) is signal.SIG_IGN
    finally:
        signal.signal(signal.SIGINT, before)
    assert signal.getsignal(signal.SIGINT) is before


def _start(args, base):
    env = dict(os.environ, PYTHONPATH=SRC, BIOTSR_PDB_BASE_URL=base)
    return subprocess.Popen([sys.executable, "-m", "biotsr", "fetch"] + args,
                            env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, preexec_fn=_default_sigint)


def _wait_for(predicate, timeout=20.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if predicate():
            return True
        time.sleep(0.05)
    return False


def _cancel(proc, signum=signal.SIGINT):
    started = time.time()
    proc.send_signal(signum)
    _out, err = proc.communicate(timeout=60)
    return time.time() - started, err


def _report(pdb_dir):
    with open(os.path.join(str(pdb_dir), "fetch_report.csv")) as fh:
        return {r["pdb_id"]: r for r in csv.DictReader(fh)}


def _leftovers(pdb_dir):
    """Anything in the structure directory besides the report and finished files."""
    out = []
    for root, _dirs, names in os.walk(str(pdb_dir)):
        out.extend(n for n in names if ".part-" in n or n.startswith(".fetch_report."))
    return out


@pytest.mark.parametrize("signum", [signal.SIGINT, signal.SIGTERM])
def test_a_stalled_response_stops_at_once(tmp_path, archive, signum):
    archive.serve("1AAA", {"status": 200, "body": _gz(_entry("1AAA")), "delay": 60})
    pdb = tmp_path / "pdb"
    proc = _start(["1AAA", "--pdb-dir", str(pdb), "--timeout", "30"], archive.base)
    assert _wait_for(lambda: archive.active > 0)
    elapsed, err = _cancel(proc, signum)
    assert proc.returncode == 130, err
    assert elapsed < PROMPT
    row = _report(pdb)["1AAA"]
    assert row["status"] == "interrupted" and row["attempts"] == "1"
    assert not (pdb / "1AAA.pdb").exists() and _leftovers(pdb) == []


def test_a_response_still_trickling_in_stops_at_once(tmp_path, archive):
    """Ten bytes every 0.2 s: never idle for --timeout, so the timeout never
    ends it; the cancellation does."""
    body = _gz(_entry("1AAA")) + b""
    archive.serve("1AAA", {"status": 200, "body": body, "trickle": (10, 0.2)})
    pdb = tmp_path / "pdb"
    proc = _start(["1AAA", "--pdb-dir", str(pdb), "--timeout", "30"], archive.base)
    assert _wait_for(lambda: archive.active > 0)
    time.sleep(1.0)  # well into the body
    elapsed, err = _cancel(proc)
    assert proc.returncode == 130, err
    assert elapsed < PROMPT
    assert _report(pdb)["1AAA"]["status"] == "interrupted"
    assert not (pdb / "1AAA.pdb").exists() and _leftovers(pdb) == []


def test_the_timeout_bounds_each_wait_not_a_whole_download(tmp_path, archive):
    """The documented limitation, pinned: a response that takes longer than
    --timeout in total but never pauses that long is completed, not cut off."""
    body = _gz(_entry("1AAA"))
    pieces = -(-len(body) // 40)
    archive.serve("1AAA", {"status": 200, "body": body, "trickle": (40, 0.1)})
    started = time.time()
    report = biotsr.fetch_structures(["1AAA"], str(tmp_path), base_url=archive.base,
                                     timeout=0.5, retries=0)
    assert report.records[0].status == "downloaded", report.records[0].message
    assert time.time() - started > 0.5 and pieces * 0.1 > 0.5


def test_a_backoff_wait_ends_at_once_and_no_retry_follows(tmp_path, archive):
    archive.serve("1AAA", {"status": 503, "body": b"",
                           "headers": {"Retry-After": "45"}})
    pdb = tmp_path / "pdb"
    proc = _start(["1AAA", "--pdb-dir", str(pdb), "--timeout", "30",
                   "--retries", "3"], archive.base)
    assert _wait_for(lambda: archive.hits.get("/download/1AAA.pdb.gz", 0) == 1)
    time.sleep(0.5)  # now inside the 45-second wait
    elapsed, err = _cancel(proc)
    assert proc.returncode == 130, err
    assert elapsed < PROMPT
    time.sleep(0.5)
    assert archive.hits["/download/1AAA.pdb.gz"] == 1
    row = _report(pdb)["1AAA"]
    assert row["status"] == "interrupted" and row["attempts"] == "1"


def test_queued_identifiers_never_start(tmp_path, archive):
    archive.serve("1AAA", {"status": 200, "body": _gz(_entry("1AAA")), "delay": 60})
    for pdb_id in ("1AAB", "1AAC"):
        archive.ok(pdb_id)
    pdb = tmp_path / "pdb"
    proc = _start(["1AAA", "1AAB", "1AAC", "--pdb-dir", str(pdb),
                   "--download-workers", "1", "--timeout", "30"], archive.base)
    assert _wait_for(lambda: archive.active > 0)
    elapsed, err = _cancel(proc)
    assert proc.returncode == 130, err
    assert elapsed < PROMPT
    rows = _report(pdb)
    assert {k: r["status"] for k, r in rows.items()} == {
        "1AAA": "interrupted", "1AAB": "not_started", "1AAC": "not_started"}
    time.sleep(0.3)
    assert set(archive.hits) == {"/download/1AAA.pdb.gz"}


class _SlowShutdownPool(concurrent.futures.ThreadPoolExecutor):
    """A pool whose shutdown starts late, as on a busy machine: long enough for
    a worker freed by the cancellation to take the next queued id first."""

    def shutdown(self, wait=True, **kwargs):
        time.sleep(0.5)
        super().shutdown(wait=wait, **kwargs)


def _interrupted_once_downloading(archive):
    """Stands in for ``as_completed``: the main thread is interrupted (as by
    Ctrl-C or SIGTERM) once a download is in flight."""
    def as_completed(futures, timeout=None):
        assert _wait_for(lambda: archive.active > 0), "no download started"
        raise KeyboardInterrupt
        yield  # pragma: no cover -- a generator, like the original
    return as_completed


def test_a_worker_freed_by_the_cancellation_starts_no_queued_id(
        tmp_path, archive, monkeypatch):
    """Regression: the cancellation stopped the request in flight before it
    cancelled the queued ids, so the worker it freed could take the next id
    first.  That id made no request but was reported ``interrupted``.  Here the
    pool's shutdown is held back, so a freed worker always gets that chance:
    the queued id must still not start -- no request, no file -- and be
    reported ``not_started``, and the one in flight ``interrupted``."""
    from biotsr import retrieval

    archive.serve("1AAA", {"status": 200, "body": _gz(_entry("1AAA")), "delay": 60})
    archive.ok("1AAB")
    monkeypatch.setattr(retrieval.concurrent.futures, "ThreadPoolExecutor",
                        _SlowShutdownPool)
    monkeypatch.setattr(retrieval.concurrent.futures, "as_completed",
                        _interrupted_once_downloading(archive))
    pdb = tmp_path / "pdb"
    started = time.time()
    with pytest.raises(retrieval.FetchInterrupted) as caught:
        biotsr.fetch_structures(["1AAA", "1AAB"], str(pdb), workers=1, timeout=30,
                                base_url=archive.base)
    assert time.time() - started < PROMPT
    records = caught.value.report.by_id()
    assert {k: r.status for k, r in records.items()} == {
        "1AAA": "interrupted", "1AAB": "not_started"}
    assert records["1AAB"].message == "not fetched: the run was interrupted first"
    time.sleep(0.3)
    assert set(archive.hits) == {"/download/1AAA.pdb.gz"}
    assert os.listdir(str(pdb)) == []


def test_an_id_taken_up_after_the_cancellation_does_nothing(tmp_path, archive):
    """However a worker comes to take an id after the cancellation -- one
    whose own download finished just as it arrived, say -- that id does not
    start: no request, no look at the cache, ``not_started``."""
    from biotsr import retrieval

    archive.ok("1AAB")
    pdb = tmp_path / "pdb"
    pdb.mkdir()
    (pdb / "1AAB.pdb").write_text("not a valid cached file\n")
    fetcher = retrieval._Fetcher(str(pdb), timeout=30, retries=0, refresh=False,
                                 base_url=archive.base, backoff=1.0)
    fetcher.stop()
    record = fetcher.fetch("1aab", "1AAB")
    assert (record.pdb_id, record.requested, record.status) == (
        "1AAB", "1aab", "not_started")
    assert record.message == "not fetched: the run was interrupted first"
    assert (record.path, record.url, record.attempts) == ("", "", 0)
    assert archive.hits == {}


def test_a_finished_download_is_reported_as_finished(tmp_path, archive):
    archive.ok("1AAB")
    archive.serve("1AAA", {"status": 200, "body": _gz(_entry("1AAA")), "delay": 60})
    pdb = tmp_path / "pdb"
    proc = _start(["1AAB", "1AAA", "--pdb-dir", str(pdb), "--download-workers", "2",
                   "--timeout", "30"], archive.base)
    assert _wait_for(lambda: (pdb / "1AAB.pdb").exists() and archive.active > 0)
    elapsed, err = _cancel(proc)
    assert proc.returncode == 130, err
    rows = _report(pdb)
    assert rows["1AAB"]["status"] == "downloaded"
    assert rows["1AAA"]["status"] == "interrupted"
    assert rows["1AAB"]["sha256"] == hashlib.sha256(
        (pdb / "1AAB.pdb").read_bytes()).hexdigest()


def test_an_interrupted_refresh_keeps_the_previous_file_and_its_record(tmp_path, archive):
    pdb = tmp_path / "pdb"
    archive.ok("1AAA")
    assert biotsr.fetch_structures(["1AAA"], str(pdb), base_url=archive.base
                                   ).records[0].status == "downloaded"
    before = (pdb / "1AAA.pdb").read_bytes()
    record = (pdb / ".biotsr-fetch" / "1AAA.json").read_bytes()
    archive.serve("1AAA", {"status": 200, "body": _gz(_entry("1AAA") + b"\n"),
                           "trickle": (10, 0.2)})
    proc = _start(["1AAA", "--pdb-dir", str(pdb), "--refresh", "--timeout", "30"],
                  archive.base)
    assert _wait_for(lambda: archive.active > 0)
    time.sleep(0.6)
    elapsed, err = _cancel(proc)
    assert proc.returncode == 130, err
    assert elapsed < PROMPT
    assert (pdb / "1AAA.pdb").read_bytes() == before
    assert (pdb / ".biotsr-fetch" / "1AAA.json").read_bytes() == record
    row = _report(pdb)["1AAA"]
    assert row["status"] == "interrupted" and "unchanged" in row["message"]
    assert _leftovers(pdb) == []


def test_an_id_cancelled_while_its_other_formats_are_checked_is_interrupted(tmp_path, archive):
    """The PDB-format file is missing (404) and the mmCIF check hangs: the id
    did not finish, so it is interrupted, not "not found"."""
    archive.serve("1AAA", {"status": 404, "body": b""})
    archive.serve("1AAA", {"status": 200, "body": b"", "delay": 60}, suffix=".cif")
    pdb = tmp_path / "pdb"
    proc = _start(["1AAA", "--pdb-dir", str(pdb), "--timeout", "30"], archive.base)
    assert _wait_for(lambda: archive.hits.get("/download/1AAA.cif", 0) == 1)
    elapsed, err = _cancel(proc)
    assert proc.returncode == 130, err
    assert elapsed < PROMPT
    row = _report(pdb)["1AAA"]
    assert row["status"] == "interrupted" and "another format" in row["message"]


def _stalled_listener():
    """A loopback port whose accept queue is full, so a new connect() hangs."""
    server = socket.socket()
    server.bind(("127.0.0.1", 0))
    server.listen(0)
    port = server.getsockname()[1]
    held = []
    for _ in range(1100):
        client = socket.socket()
        client.settimeout(0.3)
        try:
            client.connect(("127.0.0.1", port))
        except OSError:
            client.close()
            return server, port, held
        held.append(client)
    for c in held:
        c.close()
    server.close()
    return None, None, []


def test_a_connection_being_established_is_abandoned_promptly(tmp_path):
    server, port, held = _stalled_listener()
    if server is None:
        pytest.skip("this system never stalls a loopback connect on a full backlog")
    try:
        pdb = tmp_path / "pdb"
        proc = _start(["1AAA", "--pdb-dir", str(pdb), "--timeout", "30"],
                      "http://127.0.0.1:%d/download" % port)
        time.sleep(1.5)  # started, and blocked connecting
        elapsed, err = _cancel(proc)
        assert proc.returncode == 130, err
        assert elapsed < PROMPT
        row = _report(pdb)["1AAA"]
        assert row["status"] == "interrupted" and row["attempts"] == "1"
    finally:
        for c in held:
            c.close()
        server.close()
