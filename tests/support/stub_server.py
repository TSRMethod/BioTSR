"""A scripted local HTTP server standing in for the structure archive.

Used by the retrieval and batch tests so that no test ever reaches the real
archive.  It listens on a loopback address, the one place the retriever accepts
plain HTTP.
"""

from __future__ import annotations

import gzip
import http.server
import socketserver
import threading
import time


def _gz(data: bytes) -> bytes:
    return gzip.compress(data, mtime=0)


def _entry(pdb_id: str) -> bytes:
    """A small but complete PDB-format entry for ``pdb_id``."""
    lines = [
        "HEADER    TEST PROTEIN                            01-JAN-00   %s              " % pdb_id,
    ]
    atoms = [
        (1, "N", "MET", 1, 27.340, 24.430, 2.614),
        (2, "CA", "MET", 1, 26.266, 25.413, 2.842),
        (3, "C", "MET", 1, 26.913, 26.639, 3.531),
        (9, "CA", "GLN", 2, 26.335, 27.770, 7.258),
        (18, "CA", "ILE", 3, 23.970, 29.004, 8.940),
    ]
    for serial, name, res, seq, x, y, z in atoms:
        lines.append(
            "ATOM  %5d  %-3s %3s A%4d    %8.3f%8.3f%8.3f  1.00 10.00           %s"
            % (serial, name, res, seq, x, y, z, name[0])
        )
    lines.append("END")
    return ("\n".join(lines) + "\n").encode("ascii")


class _LoopbackServer(http.server.ThreadingHTTPServer):
    def server_bind(self):
        # HTTPServer.server_bind resolves the host's FQDN, which can stall for
        # tens of seconds on a machine without reverse DNS; loopback needs none.
        socketserver.TCPServer.server_bind(self)
        self.server_name, self.server_port = self.server_address[:2]


class Stub(object):
    """A scripted HTTP server: ``routes[path]`` is a list of responses, used in
    order (the last one repeats).  A response is a dict with ``status``,
    ``body``, optional ``headers``, ``delay`` (seconds before answering),
    ``length`` (a Content-Length that may lie) and ``trickle`` -- ``(bytes,
    seconds)``: send the body that many bytes at a time, pausing between
    pieces, so a response keeps progressing without ever finishing quickly."""

    def __init__(self):
        self.routes = {}
        self.hits = {}
        self.active = 0
        self.max_active = 0
        self.lock = threading.Lock()
        stub = self

        class Handler(http.server.BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def _answer(self, head_only):
                with stub.lock:
                    stub.active += 1
                    stub.max_active = max(stub.max_active, stub.active)
                    stub.hits[self.path] = stub.hits.get(self.path, 0) + 1
                    script = stub.routes.get(self.path, [{"status": 404, "body": b""}])
                    response = script[min(stub.hits[self.path], len(script)) - 1]
                try:
                    time.sleep(response.get("delay", 0))
                    body = response.get("body", b"")
                    self.send_response(response["status"])
                    for key, value in response.get("headers", {}).items():
                        self.send_header(key, value)
                    self.send_header("Content-Length",
                                     str(response.get("length", len(body))))
                    self.end_headers()
                    if not head_only:
                        trickle = response.get("trickle")
                        if trickle:
                            size, pause = trickle
                            for start in range(0, len(body), size):
                                self.wfile.write(body[start:start + size])
                                self.wfile.flush()
                                time.sleep(pause)
                        else:
                            self.wfile.write(body)
                            self.wfile.flush()
                    if response.get("length", len(body)) != len(body):
                        self.close_connection = True
                except (BrokenPipeError, ConnectionResetError):
                    pass
                finally:
                    with stub.lock:
                        stub.active -= 1

            def do_GET(self):
                self._answer(head_only=False)

            def do_HEAD(self):
                self._answer(head_only=True)

        self.server = _LoopbackServer(("127.0.0.1", 0), Handler)
        self.server.daemon_threads = True
        self.thread = threading.Thread(target=self.server.serve_forever,
                                       kwargs={"poll_interval": 0.02}, daemon=True)
        self.thread.start()
        self.base = "http://127.0.0.1:%d/download" % self.server.server_address[1]

    def serve(self, pdb_id, *responses, suffix=".pdb.gz"):
        self.routes["/download/%s%s" % (pdb_id, suffix)] = list(responses)

    def ok(self, pdb_id, data=None):
        self.serve(pdb_id, {"status": 200, "body": _gz(data or _entry(pdb_id))})

    def close(self):
        self.server.shutdown()
        self.server.server_close()
