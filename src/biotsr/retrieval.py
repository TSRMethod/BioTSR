"""Opt-in retrieval of PDB-format structures from the RCSB PDB archive.

BioTSR never touches the network unless asked: ``biotsr fetch``,
``biotsr batch --download-missing`` and :func:`fetch_structures` are the only
ways in, so a compute node with no internet access runs every calculation
exactly as before.

What a download has to survive before it is used
------------------------------------------------

A structure file that parses is not proof of a complete download: a transfer cut
off at a line boundary still parses, and silently computes keys over part of a
protein.  So every download is checked in layers, and only a file that passes
all of them is installed:

1. **Transfer.** HTTPS only, with certificate verification; a redirect that
   would leave HTTPS is refused.  The body must match ``Content-Length`` when the
   server sends one.
2. **Compression.** Files are fetched as ``.pdb.gz`` and decompressed in full.
   The gzip stream must end properly, its CRC-32 and length trailer must match
   (zlib checks both), and nothing may follow it.
3. **Content.** The decompressed text must contain ``ATOM``/``HETATM`` records,
   end with an ``END`` record (every RCSB entry does), name the requested entry
   in its ``HEADER`` when it has one, and be read without error by the same
   reader every calculation uses.
4. **Installation.** The file is written to a temporary name in the target
   directory, flushed to disk, validated, and moved into place with one
   ``os.replace``.  An interrupted download never leaves a partial file under
   the final name.

Each installed file gets a small provenance record (URL, UTC retrieval time,
SHA-256 of the compressed and decompressed bytes) in ``.biotsr-fetch/`` beside
it.  RCSB remediates entries over time, so the record says which bytes were
used; every result's sidecar also records the SHA-256 of the file it was
computed from.

The cache
---------

The structure directory is the cache.  A file already there is reused only when
it is valid: when it has a provenance record, its SHA-256 must still match it;
when it has none (a file you supplied, or one from another tool), it must pass
the content checks above.  An invalid cached file is reported and never reused
or overwritten silently.

``refresh=True`` downloads again even when a valid file is cached, and replaces
it only after the new download passes every check.  If the refresh fails, the
previous file is kept untouched and the failure is reported.

Only the PDB format is read.  RCSB serves some large entries only as mmCIF;
those are reported as ``pdb_format_unavailable`` rather than as missing.

Timeouts and cancellation
-------------------------

``timeout`` bounds **each wait**, not a whole download: establishing a TCP
connection (per address tried), each step of the TLS handshake, the wait for
the response headers, and each read of the body.  A server that keeps sending
data, however slowly, therefore keeps a download alive; the only bound on one
attempt's total duration is the ``MAX_DOWNLOAD_BYTES`` size limit.  Each
transient failure is retried ``retries`` times with exponential backoff
(at most ``MAX_BACKOFF`` seconds per wait, or longer when ``Retry-After``
asks for less than that).

Cancellation (Ctrl-C, or SIGTERM through the command line) is prompt:

* no new attempt, retry or queued identifier starts after it;
* a backoff wait ends at once;
* a TCP connection being established is abandoned within ``CONNECT_POLL``
  seconds;
* a request waiting for its response headers, or reading its body, is stopped
  at once: its socket is shut down;
* every identifier is reported -- finished ones as they finished, the one in
  flight as ``interrupted``, the rest as ``not_started`` -- and no partial or
  unvalidated file is ever installed (a file is validated in full before the one
  ``os.replace`` that installs it).

Two steps cannot be interrupted: resolving the host name (the system resolver
has its own timeouts) and a TLS handshake step in progress (bounded by
``timeout``).  Cancellation waits for them to end.

Standard library only.  This is the only module in the package that opens a
network connection.
"""

from __future__ import annotations

import concurrent.futures
import csv
import datetime
import errno
import hashlib
import json
import os
import random
import re
import selectors
import socket
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zlib
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

__all__ = [
    "DEFAULT_BASE_URL",
    "DEFAULT_WORKERS",
    "MAX_WORKERS",
    "DEFAULT_TIMEOUT",
    "DEFAULT_RETRIES",
    "DOWNLOADED",
    "REFRESHED",
    "CACHED",
    "WOULD_DOWNLOAD",
    "INVALID_CACHED",
    "NOT_FOUND",
    "PDB_FORMAT_UNAVAILABLE",
    "INVALID_ID",
    "REFRESH_FAILED",
    "FAILED",
    "USABLE_STATUSES",
    "FetchRecord",
    "FetchReport",
    "fetch_structures",
    "validate_download_settings",
    "FetchInterrupted",
    "INTERRUPTED",
    "NOT_STARTED",
    "normalise_pdb_id",
    "cached_structure_path",
    "check_cached_file",
    "file_sha256",
    "publish_file_mode",
    "REPORT_COLUMNS",
]

#: Where entries are downloaded from.  ``BIOTSR_PDB_BASE_URL`` overrides it, for
#: a site mirror; it must be HTTPS as well.
DEFAULT_BASE_URL = "https://files.rcsb.org/download"
#: Concurrent downloads.  Network I/O, not computation, so threads, bounded, and
#: independent of how many processes compute.
DEFAULT_WORKERS = 4
MAX_WORKERS = 16
#: Seconds to wait for a connection, or for the next piece of a response.  A
#: limit on each wait, not on a whole download: a server that keeps sending, however
#: slowly, is never timed out (see "Timeouts and cancellation" above).
DEFAULT_TIMEOUT = 30.0
#: Bytes read per step of a response body; cancellation is checked between steps.
READ_CHUNK = 64 * 1024
#: How often a connection attempt in progress looks for a cancellation.
CONNECT_POLL = 0.1
#: Further attempts after the first, for errors that can be transient.
DEFAULT_RETRIES = 3
#: Upper bound on one wait between attempts, whatever ``Retry-After`` says.
MAX_BACKOFF = 60.0
#: No single PDB-format entry is anywhere near this; it only bounds memory if a
#: server misbehaves.
MAX_DOWNLOAD_BYTES = 512 * 1024 * 1024

PROVENANCE_DIR = ".biotsr-fetch"

DOWNLOADED = "downloaded"
REFRESHED = "refreshed"
CACHED = "cached"
WOULD_DOWNLOAD = "would_download"
INVALID_CACHED = "invalid_cached"
NOT_FOUND = "not_found"
PDB_FORMAT_UNAVAILABLE = "pdb_format_unavailable"
INVALID_ID = "invalid_id"
REFRESH_FAILED = "refresh_failed"
FAILED = "failed"
#: Ctrl-C or SIGTERM arrived while this id was being fetched: its work had
#: begun and was stopped, and nothing of it was installed.
INTERRUPTED = "interrupted"
#: Ctrl-C or SIGTERM arrived before this id's fetch began: nothing was done for
#: it -- no request, and no look at the cache.
NOT_STARTED = "not_started"
_NOT_STARTED_MESSAGE = "not fetched: the run was interrupted first"

#: After these, a valid file is at ``record.path``.  ``refresh_failed`` keeps the
#: previous valid file, so it is usable, but the run is still reported as failed.
USABLE_STATUSES = (DOWNLOADED, REFRESHED, CACHED, REFRESH_FAILED)

REPORT_COLUMNS = (
    "pdb_id",
    "requested",
    "status",
    "path",
    "bytes",
    "sha256",
    "attempts",
    "url",
    "retrieved_at",
    "message",
)

_ID = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
_LOOPBACK = ("127.0.0.1", "::1", "localhost")


@dataclass
class FetchRecord:
    """What happened to one identifier."""

    pdb_id: str
    requested: str
    status: str
    path: str = ""
    bytes: int = 0
    sha256: str = ""
    attempts: int = 0
    url: str = ""
    retrieved_at: str = ""
    message: str = ""

    @property
    def usable(self) -> bool:
        return self.status in USABLE_STATUSES

    @property
    def failed(self) -> bool:
        return self.status not in (DOWNLOADED, REFRESHED, CACHED, WOULD_DOWNLOAD)

    def row(self) -> Dict[str, object]:
        return {name: getattr(self, name) for name in REPORT_COLUMNS}


@dataclass
class FetchReport:
    """Every identifier's outcome, in the order they were first requested."""

    records: List[FetchRecord] = field(default_factory=list)

    @property
    def failures(self) -> List[FetchRecord]:
        return [r for r in self.records if r.failed]

    @property
    def ok(self) -> bool:
        return not self.failures

    def by_id(self) -> Dict[str, FetchRecord]:
        return {r.pdb_id: r for r in self.records if r.pdb_id}

    def counts(self) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for record in self.records:
            out[record.status] = out.get(record.status, 0) + 1
        return out

    def write_csv(self, path: str) -> str:
        """Write the report atomically; returns ``path``."""
        directory = os.path.dirname(os.path.abspath(path))
        os.makedirs(directory, exist_ok=True)
        fd, tmp = tempfile.mkstemp(prefix=".fetch_report.", dir=directory)
        try:
            with os.fdopen(fd, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=list(REPORT_COLUMNS))
                writer.writeheader()
                for record in self.records:
                    writer.writerow(record.row())
            publish_file_mode(tmp)
            os.replace(tmp, path)
        except BaseException:
            _remove_quietly(tmp)
            raise
        return path


# ---------------------------------------------------------------------------
# Identifiers and the cache
# ---------------------------------------------------------------------------


def normalise_pdb_id(text: str) -> Optional[str]:
    """``"1ubq"`` -> ``"1UBQ"``; ``None`` for anything that is not a PDB id."""
    candidate = str(text or "").strip()
    if not _ID.match(candidate):
        return None
    return candidate.upper()


def cached_structure_path(pdb_dir: str, pdb_id: str) -> Optional[str]:
    """The existing file for ``pdb_id`` in ``pdb_dir``, looked up the way a batch
    looks it up (``.pdb`` or ``.ent``; as written, upper and lower case)."""
    for stem in (pdb_id, pdb_id.upper(), pdb_id.lower()):
        for ext in (".pdb", ".ent"):
            candidate = os.path.join(pdb_dir, stem + ext)
            if os.path.isfile(candidate):
                return candidate
    return None


def file_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _provenance_path(pdb_dir: str, pdb_id: str) -> str:
    return os.path.join(pdb_dir, PROVENANCE_DIR, "%s.json" % pdb_id)


def _read_provenance(pdb_dir: str, pdb_id: str) -> Optional[dict]:
    path = _provenance_path(pdb_dir, pdb_id)
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def check_cached_file(path: str, pdb_id: str, pdb_dir: str) -> Optional[str]:
    """``None`` when the cached file is valid, otherwise the reason it is not.

    Read-only.  A file BioTSR downloaded is checked against its recorded
    SHA-256; any other file against the content rules (ATOM/HETATM records, a
    final ``END`` record, a matching ``HEADER`` id).  A file without an ``END``
    record is refused even though it may parse: that is exactly what a download
    cut off at a line boundary looks like.
    """
    provenance = _read_provenance(pdb_dir, pdb_id)
    if provenance and provenance.get("file") == os.path.basename(path):
        actual = file_sha256(path)
        if actual != provenance.get("sha256"):
            return (
                "it no longer matches the SHA-256 recorded when it was downloaded "
                "(%s...), so it was changed or truncated afterwards"
                % str(provenance.get("sha256"))[:12]
            )
        return None
    with open(path, "rb") as fh:
        data = fh.read()
    return _content_problem(data, pdb_id)


def _content_problem(data: bytes, pdb_id: str) -> Optional[str]:
    if not data.strip():
        return "the file is empty"
    if b"\x00" in data:
        return "the file contains NUL bytes, so it is not PDB-format text"
    text = data.decode("latin-1")
    records = [line[:6].rstrip() for line in text.splitlines() if line.strip()]
    if not any(r in ("ATOM", "HETATM") for r in records):
        return "it has no ATOM or HETATM records"
    if not records or records[-1] != "END":
        return (
            "it does not end with an END record; every RCSB entry does, so the "
            "file is probably truncated"
        )
    for line in text.splitlines():
        if line.startswith("HEADER"):
            named = line[62:66].strip().upper()
            if named and named != pdb_id.upper():
                return "its HEADER names entry %s, not %s" % (named, pdb_id)
            break
    return None


def _obsolete_notice(data: bytes) -> str:
    """A notice when the entry declares itself obsolete (an ``OBSLTE`` record).

    RCSB still serves obsolete entries, and they are valid structures, so the
    file is used -- but a dataset built on one should know it was superseded.
    """
    for line in data.decode("latin-1").splitlines():
        if line.startswith("OBSLTE"):
            replaced_by = " ".join(line[31:].split())
            return (
                "note: this entry is obsolete (OBSLTE record)%s"
                % ("; superseded by %s" % replaced_by if replaced_by else "")
            )
        if line.startswith(("ATOM", "HETATM", "MODEL")):
            break
    return ""


# ---------------------------------------------------------------------------
# Transport
# ---------------------------------------------------------------------------


class _Cancelled(Exception):
    """The run was interrupted; this id stops without another attempt."""


class _Transient(Exception):
    """Worth another attempt: the next one may succeed."""

    def __init__(self, message: str, retry_after: Optional[float] = None):
        super().__init__(message)
        self.retry_after = retry_after


class _Permanent(Exception):
    """Another attempt would get the same answer."""

    def __init__(self, message: str, status: str = FAILED):
        super().__init__(message)
        self.status = status


def _check_url(url: str) -> None:
    parts = urllib.parse.urlsplit(url)
    if parts.scheme == "https":
        return
    if parts.scheme == "http" and parts.hostname in _LOOPBACK:
        return  # a local mirror or a test server; never a remote host
    raise ValueError(
        "refusing %r: structures are downloaded over HTTPS only (plain HTTP is "
        "accepted for a loopback address and nothing else)" % url
    )


class _HTTPSOnlyRedirects(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        _check_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


class _Transfers(object):
    """The sockets of the requests in flight, one per worker thread.

    Registered as soon as a connection is established, so a cancellation can
    shut each of them down: a thread blocked waiting for response headers, or
    in the middle of a body, then returns at once instead of at the next
    timeout.  A connection made after the cancellation is refused on
    registration.
    """

    def __init__(self, cancel: threading.Event):
        self.cancel = cancel
        self._lock = threading.Lock()
        self._sockets: Dict[int, socket.socket] = {}

    def add(self, sock: socket.socket) -> None:
        with self._lock:
            if self.cancel.is_set():
                _shutdown(sock)
                raise _Cancelled()
            self._sockets[threading.get_ident()] = sock

    def release(self) -> None:
        with self._lock:
            self._sockets.pop(threading.get_ident(), None)

    def abort(self) -> None:
        with self._lock:
            sockets = list(self._sockets.values())
        for sock in sockets:
            _shutdown(sock)


def _shutdown(sock) -> None:
    """End every blocked read on ``sock``; a socket already closed is left alone.

    ``socket.socket.shutdown`` is called on the object itself (for a TLS socket,
    on its file descriptor without touching the TLS state another thread may
    be using).  A closed socket object holds no descriptor, so this can never
    reach a descriptor the process has since reused.
    """
    try:
        socket.socket.shutdown(sock, socket.SHUT_RDWR)
    except (OSError, ValueError):
        pass


def _cancellable_connect(address, timeout, source_address=None, *, cancel):
    """``socket.create_connection``, abandoning the attempt when ``cancel`` is set.

    The connection is made non-blocking and polled every ``CONNECT_POLL``
    seconds, so a cancellation does not wait for ``timeout``.  Name resolution
    happens first and is not interruptible.
    """
    host, port = address
    if not isinstance(timeout, (int, float)) or timeout <= 0:
        timeout = DEFAULT_TIMEOUT
    error: Optional[BaseException] = None
    for family, kind, proto, _name, where in socket.getaddrinfo(
            host, port, 0, socket.SOCK_STREAM):
        if cancel.is_set():
            raise _Cancelled()
        sock = socket.socket(family, kind, proto)
        try:
            if source_address:
                sock.bind(source_address)
            sock.setblocking(False)
            deadline = time.monotonic() + timeout
            code = sock.connect_ex(where)
            # EINTR from a non-blocking connect means "still in progress".
            if code not in (0, errno.EINPROGRESS, errno.EWOULDBLOCK, errno.EALREADY,
                            errno.EAGAIN, errno.EINTR):
                raise OSError(code, os.strerror(code))
            if code != 0:
                with selectors.DefaultSelector() as selector:
                    selector.register(sock, selectors.EVENT_WRITE)
                    while True:
                        if cancel.is_set():
                            raise _Cancelled()
                        remaining = deadline - time.monotonic()
                        if remaining <= 0:
                            raise socket.timeout("timed out")
                        if selector.select(min(CONNECT_POLL, remaining)):
                            break
                code = sock.getsockopt(socket.SOL_SOCKET, socket.SO_ERROR)
                if code:
                    raise OSError(code, os.strerror(code))
            sock.settimeout(timeout)
            return sock
        except _Cancelled:
            sock.close()
            raise
        except OSError as exc:
            sock.close()
            error = exc
    if error is not None:
        raise error
    raise OSError("no address found for %s" % host)


def _connection_classes(transfers: "_Transfers"):
    """HTTP(S) connections that connect cancellably and register their socket."""
    import functools
    import http.client

    connect = functools.partial(_cancellable_connect, cancel=transfers.cancel)

    class _HTTPConnection(http.client.HTTPConnection):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._create_connection = connect

        def connect(self):
            super().connect()
            transfers.add(self.sock)

    class _HTTPSConnection(http.client.HTTPSConnection):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._create_connection = connect

        def connect(self):
            super().connect()
            transfers.add(self.sock)

    return _HTTPConnection, _HTTPSConnection


def _opener(transfers: Optional["_Transfers"] = None):
    # Imported here, not at module scope: a Python built without OpenSSL can
    # still import biotsr and compute; only a download needs TLS.
    import ssl

    context = ssl.create_default_context()
    if transfers is None:
        return urllib.request.build_opener(
            _HTTPSOnlyRedirects(),
            urllib.request.HTTPSHandler(context=context),
        )
    http_class, https_class = _connection_classes(transfers)

    class _HTTP(urllib.request.HTTPHandler):
        def http_open(self, req):
            return self.do_open(http_class, req)

    class _HTTPS(urllib.request.HTTPSHandler):
        def https_open(self, req):
            return self.do_open(https_class, req, context=context)

    return urllib.request.build_opener(
        _HTTPSOnlyRedirects(), _HTTP(), _HTTPS(context=context),
    )


def _user_agent() -> str:
    from . import __version__

    return "biotsr/%s (Python urllib)" % __version__


def _retry_after(headers) -> Optional[float]:
    value = (headers or {}).get("Retry-After") if headers is not None else None
    if not value:
        return None
    try:
        return max(0.0, float(value))
    except ValueError:
        pass
    try:
        from email.utils import parsedate_to_datetime

        when = parsedate_to_datetime(value)
        now = datetime.datetime.now(when.tzinfo)
        return max(0.0, (when - now).total_seconds())
    except (TypeError, ValueError):
        return None


def _read_body(response, cancel: Optional[threading.Event]) -> bytes:
    """The body in ``READ_CHUNK`` steps, at most ``MAX_DOWNLOAD_BYTES + 1`` bytes."""
    chunks = []
    total = 0
    while total <= MAX_DOWNLOAD_BYTES:
        if cancel is not None and cancel.is_set():
            raise _Cancelled()
        block = response.read(min(READ_CHUNK, MAX_DOWNLOAD_BYTES + 1 - total))
        if not block:
            break
        chunks.append(block)
        total += len(block)
    return b"".join(chunks)


def _get(opener, url: str, timeout: float,
         cancel: Optional[threading.Event] = None) -> bytes:
    request = urllib.request.Request(
        url, headers={"User-Agent": _user_agent(), "Accept-Encoding": "identity"}
    )

    def stopped() -> bool:
        return cancel is not None and cancel.is_set()

    try:
        with opener.open(request, timeout=timeout) as response:
            expected = response.headers.get("Content-Length")
            body = _read_body(response, cancel)
    except _Cancelled:
        raise
    except urllib.error.HTTPError as exc:
        if stopped():
            raise _Cancelled()
        if exc.code in (404, 410):
            raise _Permanent("HTTP %d" % exc.code, status=NOT_FOUND)
        if exc.code == 429 or exc.code >= 500:
            raise _Transient("HTTP %d" % exc.code, _retry_after(exc.headers))
        raise _Permanent("HTTP %d %s" % (exc.code, exc.reason))
    except urllib.error.URLError as exc:
        import ssl

        if stopped():
            raise _Cancelled()
        reason = exc.reason
        if isinstance(reason, ssl.SSLCertVerificationError):
            raise _Permanent("TLS certificate verification failed: %s" % reason)
        raise _Transient("connection failed: %s" % reason)
    except (socket.timeout, TimeoutError):
        if stopped():
            raise _Cancelled()
        raise _Transient("timed out after %gs" % timeout)
    except (ConnectionError, OSError) as exc:
        # http.client.IncompleteRead and RemoteDisconnected land here too.
        if stopped():
            raise _Cancelled()
        raise _Transient("transfer interrupted: %s" % exc)
    except Exception as exc:  # noqa: BLE001 -- http.client.HTTPException et al.
        if stopped():
            raise _Cancelled()
        raise _Transient("transfer failed: %s: %s" % (type(exc).__name__, exc))
    if stopped():
        # A socket shut down by the cancellation ends the body early, which
        # would otherwise look like an ordinary short transfer.
        raise _Cancelled()
    if len(body) > MAX_DOWNLOAD_BYTES:
        raise _Permanent("response larger than %d bytes" % MAX_DOWNLOAD_BYTES)
    if expected is not None:
        try:
            wanted = int(expected)
        except ValueError:
            wanted = None
        if wanted is not None and wanted != len(body):
            raise _Transient(
                "incomplete transfer: received %d of %d bytes" % (len(body), wanted)
            )
    return body


def _exists(opener, url: str, timeout: float) -> Optional[bool]:
    """``True``/``False`` for a HEAD request, ``None`` when it cannot be told."""
    request = urllib.request.Request(
        url, method="HEAD", headers={"User-Agent": _user_agent()}
    )
    try:
        with opener.open(request, timeout=timeout):
            return True
    except _Cancelled:
        raise
    except urllib.error.HTTPError as exc:
        return False if exc.code in (404, 410) else None
    except Exception:  # noqa: BLE001 -- a diagnostic, never the error itself
        return None


def _gunzip_strict(data: bytes) -> bytes:
    """Decompress one complete gzip member, or raise ``_Transient``.

    zlib in gzip mode verifies the CRC-32 and the length trailer; ``eof`` proves
    the stream ended, and ``unused_data`` that nothing follows it.
    """
    decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
    try:
        out = decompressor.decompress(data) + decompressor.flush()
    except zlib.error as exc:
        raise _Transient("corrupt gzip data: %s" % exc)
    if not decompressor.eof:
        raise _Transient("truncated gzip stream")
    if decompressor.unused_data:
        raise _Transient("unexpected data after the gzip stream")
    return out


# ---------------------------------------------------------------------------
# One identifier
# ---------------------------------------------------------------------------


def _utc_now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


_FILE_MODE: Optional[int] = None
_MODE_LOCK = threading.Lock()


def publish_file_mode(path: str) -> None:
    """Give a file made by ``tempfile.mkstemp`` ordinary permissions.

    ``mkstemp`` creates files readable by their owner only, and ``os.replace``
    keeps that mode, so a structure directory shared by a group on a cluster
    would hold files nobody else can read.  The final file gets what ``open``
    would have given it: ``0o666`` less the process umask.
    """
    global _FILE_MODE
    with _MODE_LOCK:
        if _FILE_MODE is None:
            # os.umask is the only portable way to read the umask; it is read
            # once, under a lock, and restored immediately.
            current = os.umask(0o022)
            os.umask(current)
            _FILE_MODE = 0o666 & ~current
    try:
        os.chmod(path, _FILE_MODE)
    except OSError:
        pass


def _remove_quietly(path: str) -> None:
    try:
        os.remove(path)
    except OSError:
        pass


def _fsync_directory(directory: str) -> None:
    try:
        fd = os.open(directory, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(fd)
    except OSError:
        pass
    finally:
        os.close(fd)


def _write_atomically(directory: str, name: str, data: bytes, prefix: str,
                      validate: Optional[Callable[[str], None]] = None) -> str:
    """Write ``data`` as ``directory/name`` via a private temporary file.

    The temporary name carries this process's id and a random suffix, so two
    concurrent invocations never touch each other's files, and only the file
    this call created is ever removed.
    """
    fd, tmp = tempfile.mkstemp(prefix="%s.part-%d-" % (prefix, os.getpid()),
                               dir=directory)
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        if validate is not None:
            validate(tmp)
        publish_file_mode(tmp)
        final = os.path.join(directory, name)
        os.replace(tmp, final)
    except BaseException:
        _remove_quietly(tmp)
        raise
    _fsync_directory(directory)
    return final


def _parse_check(pdb_id: str) -> Callable[[str], None]:
    def check(path: str) -> None:
        from .domain import StructureId
        from .structure.backends.pdb_legacy import read_pdb_legacy

        try:
            structure = read_pdb_legacy(path, StructureId.parse(pdb_id))
        except Exception as exc:  # noqa: BLE001 -- reported, never swallowed
            raise _Permanent("the structure reader refused it: %s" % exc)
        if not list(structure.chains()):
            raise _Permanent("the structure reader found no chains in it")

    return check


class _Fetcher(object):
    def __init__(self, pdb_dir, *, timeout, retries, refresh, base_url, backoff,
                 sleep=time.sleep):
        self.pdb_dir = pdb_dir
        self.timeout = float(timeout)
        self.retries = int(retries)
        self.refresh = bool(refresh)
        self.base_url = base_url.rstrip("/")
        self.backoff = float(backoff)
        self.sleep = sleep
        self.cancel = threading.Event()
        self.transfers = _Transfers(self.cancel)
        self._local = threading.local()

    def opener(self):
        if getattr(self._local, "opener", None) is None:
            self._local.opener = _opener(self.transfers)
        return self._local.opener

    def stop(self) -> None:
        """Cancel: no new attempt starts, and every request in flight ends now."""
        self.cancel.set()
        self.transfers.abort()

    def url(self, pdb_id: str, suffix: str) -> str:
        return "%s/%s%s" % (self.base_url, pdb_id, suffix)

    def fetch(self, requested: str, pdb_id: str) -> FetchRecord:
        if self.cancel.is_set():
            # Taken from the queue after the cancellation (a worker freed just
            # as it arrived): this id does not start.
            return FetchRecord(pdb_id, requested, NOT_STARTED,
                               message=_NOT_STARTED_MESSAGE)
        record = FetchRecord(pdb_id=pdb_id, requested=requested, status=FAILED)
        existing = cached_structure_path(self.pdb_dir, pdb_id)
        if existing and not self.refresh:
            problem = check_cached_file(existing, pdb_id, self.pdb_dir)
            record.path = existing
            if problem is None:
                record.status = CACHED
                record.bytes = os.path.getsize(existing)
                record.sha256 = file_sha256(existing)
                provenance = _read_provenance(self.pdb_dir, pdb_id) or {}
                record.url = provenance.get("url", "")
                record.retrieved_at = provenance.get("retrieved_at", "")
                return record
            record.status = INVALID_CACHED
            record.message = (
                "%s is not reused: %s. Pass --refresh (refresh=True) to download "
                "it again; to use a file of your own instead, replace it and also "
                "delete its download record %s, which still holds the old "
                "checksum" % (existing, problem,
                              os.path.join(self.pdb_dir, PROVENANCE_DIR,
                                           "%s.json" % pdb_id))
            )
            return record

        url = self.url(pdb_id, ".pdb.gz")
        record.url = url
        try:
            data, text = self._download(url, record)
            problem = _content_problem(text, pdb_id)
            if problem:
                raise _Permanent("the downloaded file is not usable: %s" % problem)
            if self.cancel.is_set():
                # Complete and valid, but the run was cancelled before it was
                # installed: nothing is written after a cancellation.
                raise _Cancelled()
            # A refresh replaces the very file a batch would read, whatever its
            # spelling (1ubq.ent); a new entry is written as <ID>.pdb.
            name = os.path.basename(existing) if existing else "%s.pdb" % pdb_id
            parse_check = _parse_check(pdb_id)

            def validate(path: str) -> None:
                parse_check(path)
                # Checked again here, the last moment before os.replace: a
                # cancellation that arrived while the file was parsed still
                # leaves nothing installed.
                if self.cancel.is_set():
                    raise _Cancelled()

            final = _write_atomically(
                self.pdb_dir, name, text, prefix=".%s.pdb" % pdb_id,
                validate=validate,
            )
            record.path = final
            record.bytes = len(text)
            record.sha256 = hashlib.sha256(text).hexdigest()
            record.retrieved_at = _utc_now()
            record.status = REFRESHED if existing else DOWNLOADED
            record.message = _obsolete_notice(text)
            self._record_provenance(pdb_id, name, url, data, record)
            return record
        except _Permanent as exc:
            if exc.status == NOT_FOUND:
                try:
                    record.status, record.message = self._not_found(pdb_id)
                except _Cancelled:
                    # The 404 was real, but whether an mmCIF entry exists
                    # was not established: the id did not finish.
                    record.status = INTERRUPTED
                    record.message = ("stopped: the run was interrupted while "
                                      "checking whether %s exists in another "
                                      "format" % pdb_id)
                    return record
            else:
                record.status = FAILED
                record.message = str(exc)
        except _Cancelled:
            record.status = INTERRUPTED
            record.message = "stopped: the run was interrupted"
            if existing:
                record.message += "; the previous file %s is unchanged" % existing
                record.path = existing
            return record
        except _Transient as exc:
            record.status = FAILED
            record.message = "%s (after %d attempt%s)" % (
                exc, record.attempts, "" if record.attempts == 1 else "s")
        except Exception as exc:  # noqa: BLE001 -- e.g. a full or read-only disk
            record.status = FAILED
            record.message = "%s: %s" % (type(exc).__name__, exc)
        if existing:
            record.message = (
                "refresh failed, the previous file %s is kept unchanged: %s"
                % (existing, record.message)
            )
            record.status = REFRESH_FAILED
            record.path = existing
        return record

    def _download(self, url: str, record: FetchRecord) -> Tuple[bytes, bytes]:
        """``(compressed, decompressed)``, retrying what can be transient."""
        opener = self.opener()
        attempt = 0
        while True:
            if self.cancel.is_set():
                raise _Cancelled()
            attempt += 1
            record.attempts = attempt
            try:
                data = _get(opener, url, self.timeout, self.cancel)
                # A corrupt body is a transfer problem: retry it like one.
                return data, _gunzip_strict(data)
            except _Transient as exc:
                if attempt > self.retries:
                    raise
                wait = self.backoff * (2 ** (attempt - 1))
                wait += random.uniform(0, self.backoff)
                if exc.retry_after is not None:
                    wait = max(wait, exc.retry_after)
                wait = min(wait, MAX_BACKOFF)
                if self.sleep is time.sleep:
                    if self.cancel.wait(wait):
                        raise _Cancelled()
                else:
                    self.sleep(wait)
            finally:
                self.transfers.release()

    def _not_found(self, pdb_id: str) -> Tuple[str, str]:
        """Classify a 404; raises ``_Cancelled`` if the run is cancelled meanwhile."""
        if self.cancel.is_set():
            raise _Cancelled()
        try:
            cif = _exists(self.opener(), self.url(pdb_id, ".cif"), self.timeout)
        finally:
            self.transfers.release()
        if self.cancel.is_set():
            # A HEAD request cut short by the cancellation answers nothing.
            raise _Cancelled()
        if cif:
            return (
                PDB_FORMAT_UNAVAILABLE,
                "RCSB has %s only in mmCIF format; BioTSR reads PDB format, so "
                "convert the entry or supply the file yourself" % pdb_id,
            )
        if cif is None:
            return (NOT_FOUND, "no PDB-format file for %s (whether an mmCIF "
                               "version exists could not be checked)" % pdb_id)
        return NOT_FOUND, "no entry %s at %s" % (pdb_id, self.base_url)

    def _record_provenance(self, pdb_id, name, url, data, record) -> None:
        directory = os.path.join(self.pdb_dir, PROVENANCE_DIR)
        os.makedirs(directory, exist_ok=True)
        payload = {
            "pdb_id": pdb_id,
            "file": name,
            "url": url,
            "retrieved_at": record.retrieved_at,
            "sha256": record.sha256,
            "bytes": record.bytes,
            "compressed_sha256": hashlib.sha256(data).hexdigest(),
            "compressed_bytes": len(data),
            "retrieved_by": _user_agent(),
        }
        blob = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
        _write_atomically(directory, "%s.json" % pdb_id, blob, prefix=".%s" % pdb_id)


# ---------------------------------------------------------------------------
# The public entry point
# ---------------------------------------------------------------------------


def validate_download_settings(*, workers, timeout, retries, backoff=1.0) -> None:
    """Refuse download settings that cannot be honoured, naming the setting.

    Shared by :func:`fetch_structures`, ``biotsr fetch``, ``biotsr batch
    --download-missing`` and :func:`biotsr.batch.run`, so a value one of them
    refuses is never silently replaced by a default in another.  Nothing is
    clamped: a worker count above the maximum is refused rather than lowered.
    """
    def whole(name, value):
        if isinstance(value, bool) or not isinstance(value, int):
            if isinstance(value, float) and value.is_integer():
                return int(value)
            raise ValueError("%s must be a whole number (got %r)" % (name, value))
        return value

    workers = whole("download workers", workers)
    if not 1 <= workers <= MAX_WORKERS:
        raise ValueError(
            "download workers must be from 1 to %d (got %d)" % (MAX_WORKERS, workers)
        )
    retries = whole("retries", retries)
    if retries < 0:
        raise ValueError("retries cannot be negative (got %d)" % retries)
    for name, value, allow_zero in (("timeout", timeout, False),
                                    ("backoff", backoff, True)):
        if isinstance(value, bool):
            raise ValueError("%s must be a number of seconds (got %r)" % (name, value))
        try:
            number = float(value)
        except (TypeError, ValueError):
            raise ValueError("%s must be a number of seconds (got %r)" % (name, value))
        if number != number or number in (float("inf"), float("-inf")):
            raise ValueError("%s must be a finite number of seconds (got %r)"
                             % (name, value))
        if number < 0 or (number == 0 and not allow_zero):
            raise ValueError("%s must be greater than 0 seconds (got %r)"
                             % (name, value))


class FetchInterrupted(KeyboardInterrupt):
    """Ctrl-C or SIGTERM during :func:`fetch_structures`.

    ``report`` holds every id: the ones that finished, as they finished, and the
    rest as ``interrupted`` (in progress) or ``not_started``.
    """

    def __init__(self, report: "FetchReport"):
        super().__init__("interrupted")
        self.report = report


def _unique_ids(ids: Iterable[str]) -> Tuple[List[Tuple[str, str]], List[FetchRecord]]:
    """``([(requested, pdb_id), ...], invalid_records)``, first spelling wins."""
    seen: Dict[str, str] = {}
    wanted: List[Tuple[str, str]] = []
    invalid: List[FetchRecord] = []
    for requested in ids:
        text = str(requested or "").strip()
        if not text:
            continue
        pdb_id = normalise_pdb_id(text)
        if pdb_id is None:
            invalid.append(FetchRecord(
                pdb_id="", requested=text, status=INVALID_ID,
                message="%r is not a PDB id (four characters, the first a digit)"
                        % text,
            ))
            continue
        if pdb_id not in seen:
            seen[pdb_id] = text
            wanted.append((text, pdb_id))
    return wanted, invalid


def fetch_structures(
    ids: Sequence[str],
    pdb_dir: str,
    *,
    workers: int = DEFAULT_WORKERS,
    timeout: float = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    refresh: bool = False,
    dry_run: bool = False,
    base_url: Optional[str] = None,
    backoff: float = 1.0,
    progress: Optional[Callable[[FetchRecord], None]] = None,
) -> FetchReport:
    """Make sure ``pdb_dir`` holds a valid PDB-format file for every id.

    ``ids`` may repeat and mix case; each entry is handled once.  Returns a
    :class:`FetchReport` with one record per distinct id (plus one per string
    that is not an id), in first-seen order.  Nothing is raised for an id that
    fails -- its record says why -- so one bad id never costs the others.

    ``dry_run=True`` makes no request and writes nothing: cached files are
    checked, and everything else is reported as ``would_download``.

    Example::

        report = biotsr.fetch_structures(["1UBQ", "4hhb"], "pdb/", workers=4)
        for record in report.failures:
            print(record.pdb_id, record.status, record.message)
    """
    base_url = base_url or os.environ.get("BIOTSR_PDB_BASE_URL") or DEFAULT_BASE_URL
    _check_url(base_url)
    validate_download_settings(workers=workers, timeout=timeout, retries=retries,
                               backoff=backoff)
    workers = int(workers)
    wanted, invalid = _unique_ids(ids)

    if dry_run:
        records = []
        for requested, pdb_id in wanted:
            existing = cached_structure_path(pdb_dir, pdb_id) if os.path.isdir(pdb_dir) else None
            if existing and not refresh:
                problem = check_cached_file(existing, pdb_id, pdb_dir)
                if problem is None:
                    records.append(FetchRecord(pdb_id, requested, CACHED, path=existing,
                                               bytes=os.path.getsize(existing)))
                else:
                    records.append(FetchRecord(pdb_id, requested, INVALID_CACHED,
                                               path=existing, message=problem))
            else:
                records.append(FetchRecord(
                    pdb_id, requested, WOULD_DOWNLOAD,
                    path=os.path.join(pdb_dir, "%s.pdb" % pdb_id),
                    url="%s/%s.pdb.gz" % (base_url.rstrip("/"), pdb_id),
                    message="would replace the cached file" if existing else "",
                ))
        return FetchReport(records + invalid)

    os.makedirs(pdb_dir, exist_ok=True)
    fetcher = _Fetcher(pdb_dir, timeout=timeout, retries=retries, refresh=refresh,
                       base_url=base_url, backoff=backoff)
    results: Dict[str, FetchRecord] = {}
    if wanted:
        pool = concurrent.futures.ThreadPoolExecutor(
            max_workers=min(workers, len(wanted)), thread_name_prefix="biotsr-fetch"
        )
        futures = {}
        try:
            futures = {
                pool.submit(fetcher.fetch, requested, pdb_id): pdb_id
                for requested, pdb_id in wanted
            }
            for future in concurrent.futures.as_completed(futures):
                record = future.result()
                results[futures[future]] = record
                if progress is not None:
                    progress(record)
        except KeyboardInterrupt:
            # Stop promptly: no new attempt starts, a backoff wait ends at once,
            # queued ids never start, a connection being made is abandoned and a
            # request in flight has its socket shut down. Only name resolution
            # and a TLS handshake step are waited for (module docstring).
            # The queued ids are cancelled first, while every worker is still
            # busy: stopping the requests in flight frees their workers, which
            # would otherwise take the next queued id before it is cancelled.
            pool.shutdown(wait=False, cancel_futures=True)
            fetcher.stop()
            pool.shutdown(wait=True)
            records = []
            for requested, pdb_id in wanted:
                if pdb_id in results:
                    records.append(results[pdb_id])
                    continue
                future = next((f for f, i in futures.items() if i == pdb_id), None)
                done = None
                if future is not None and future.done() and not future.cancelled():
                    try:
                        done = future.result()
                    except BaseException:  # noqa: BLE001
                        done = None
                if done is not None:
                    records.append(done)
                else:
                    records.append(FetchRecord(
                        pdb_id, requested, NOT_STARTED,
                        message=_NOT_STARTED_MESSAGE))
            raise FetchInterrupted(FetchReport(records + invalid)) from None
        except BaseException:
            fetcher.stop()
            pool.shutdown(wait=True, cancel_futures=True)
            raise
        pool.shutdown(wait=True)
    return FetchReport([results[pdb_id] for _requested, pdb_id in wanted] + invalid)
