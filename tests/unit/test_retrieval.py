"""Structure retrieval, against a local stub server -- never the real archive.

Every behaviour the retrieval contract promises is exercised here: validated
downloads, cache reuse and cache refusal, refresh semantics, not-found versus
mmCIF-only, retries on transient failures, rejection of incomplete or corrupt
transfers without leaving a file behind, bounded concurrency and HTTPS-only
transport.  The server listens on a loopback address, which is the one place
plain HTTP is accepted.
"""

from __future__ import annotations

import hashlib
import json
import os
import time

import pytest

import biotsr
from biotsr import retrieval as R

from ..conftest import STRUCTURES
from ..support.stub_server import Stub, _entry, _gz

with open(os.path.join(STRUCTURES, "1UBQ.pdb"), "rb") as _fh:
    UBQ = _fh.read()


@pytest.fixture
def stub():
    server = Stub()
    yield server
    server.close()


def fetch(stub, ids, pdb_dir, **kwargs):
    kwargs.setdefault("backoff", 0.01)
    kwargs.setdefault("timeout", 5.0)
    return R.fetch_structures(ids, str(pdb_dir), base_url=stub.base, **kwargs)


def _parts(pdb_dir):
    return [n for n in os.listdir(str(pdb_dir)) if ".part-" in n]


# ---------------------------------------------------------------------------
# Success, the cache and provenance
# ---------------------------------------------------------------------------


def test_a_download_is_validated_installed_and_recorded(stub, tmp_path):
    stub.ok("1ABC")
    report = fetch(stub, ["1abc"], tmp_path)
    record = report.records[0]
    assert (record.pdb_id, record.requested, record.status) == ("1ABC", "1abc", "downloaded")
    path = tmp_path / "1ABC.pdb"
    assert record.path == str(path)
    assert path.read_bytes() == _entry("1ABC")
    assert record.sha256 == hashlib.sha256(_entry("1ABC")).hexdigest()
    assert record.attempts == 1
    assert record.url.endswith("/download/1ABC.pdb.gz")
    provenance = json.loads((tmp_path / ".biotsr-fetch" / "1ABC.json").read_text())
    assert provenance["sha256"] == record.sha256
    assert provenance["compressed_sha256"] == hashlib.sha256(_gz(_entry("1ABC"))).hexdigest()
    assert provenance["url"] == record.url and provenance["retrieved_at"]
    assert _parts(tmp_path) == []
    # and the file is exactly what a calculation reads
    assert biotsr.load_structure(str(path)).chains()


def test_a_real_structure_round_trips_byte_for_byte(stub, tmp_path):
    stub.ok("1UBQ", UBQ)
    record = fetch(stub, ["1UBQ"], tmp_path).records[0]
    assert record.status == "downloaded"
    assert (tmp_path / "1UBQ.pdb").read_bytes() == UBQ


def test_a_valid_cached_file_is_reused_without_a_request(stub, tmp_path):
    stub.ok("1ABC")
    fetch(stub, ["1ABC"], tmp_path)
    before = (tmp_path / "1ABC.pdb").stat().st_mtime_ns
    report = fetch(stub, ["1ABC"], tmp_path)
    assert report.records[0].status == "cached"
    assert stub.hits["/download/1ABC.pdb.gz"] == 1
    assert (tmp_path / "1ABC.pdb").stat().st_mtime_ns == before


def test_a_file_supplied_by_hand_is_reused_when_it_is_complete(stub, tmp_path):
    (tmp_path / "1abc.pdb").write_bytes(_entry("1ABC"))
    record = fetch(stub, ["1ABC"], tmp_path).records[0]
    assert record.status == "cached"
    assert os.path.samefile(record.path, str(tmp_path / "1abc.pdb"))
    assert stub.hits == {}


def test_a_cached_file_without_an_end_record_is_refused_not_reused(stub, tmp_path):
    truncated = _entry("1ABC").rsplit(b"END", 1)[0]
    (tmp_path / "1ABC.pdb").write_bytes(truncated)
    record = fetch(stub, ["1ABC"], tmp_path).records[0]
    assert record.status == "invalid_cached"
    assert "END record" in record.message and "--refresh" in record.message
    assert (tmp_path / "1ABC.pdb").read_bytes() == truncated   # untouched
    assert stub.hits == {}


def test_a_downloaded_file_that_changed_afterwards_is_refused(stub, tmp_path):
    stub.ok("1ABC")
    fetch(stub, ["1ABC"], tmp_path)
    with open(str(tmp_path / "1ABC.pdb"), "ab") as fh:
        fh.write(b"REMARK edited\n")
    record = fetch(stub, ["1ABC"], tmp_path).records[0]
    assert record.status == "invalid_cached"
    assert "SHA-256 recorded when it was downloaded" in record.message


# ---------------------------------------------------------------------------
# Refresh
# ---------------------------------------------------------------------------


def test_refresh_replaces_a_valid_file_only_after_the_new_one_validates(stub, tmp_path):
    stub.ok("1ABC")
    fetch(stub, ["1ABC"], tmp_path)
    newer = _entry("1ABC").replace(b"TEST PROTEIN", b"TEST REMEDIE")
    stub.ok("1ABC", newer)
    record = fetch(stub, ["1ABC"], tmp_path, refresh=True).records[0]
    assert record.status == "refreshed"
    assert (tmp_path / "1ABC.pdb").read_bytes() == newer
    provenance = json.loads((tmp_path / ".biotsr-fetch" / "1ABC.json").read_text())
    assert provenance["sha256"] == hashlib.sha256(newer).hexdigest()


def test_a_failed_refresh_keeps_the_previous_file(stub, tmp_path):
    stub.ok("1ABC")
    fetch(stub, ["1ABC"], tmp_path)
    original = (tmp_path / "1ABC.pdb").read_bytes()
    stub.serve("1ABC", {"status": 200, "body": b"not gzip at all"})
    record = fetch(stub, ["1ABC"], tmp_path, refresh=True, retries=1).records[0]
    assert record.status == "refresh_failed"
    assert record.usable and record.failed
    assert "previous file" in record.message and "kept unchanged" in record.message
    assert (tmp_path / "1ABC.pdb").read_bytes() == original
    assert _parts(tmp_path) == []


def test_refresh_repairs_an_invalid_cached_file(stub, tmp_path):
    (tmp_path / "1ABC.pdb").write_bytes(b"ATOM  truncated\n")
    stub.ok("1ABC")
    record = fetch(stub, ["1ABC"], tmp_path, refresh=True).records[0]
    assert record.status == "refreshed"
    assert (tmp_path / "1ABC.pdb").read_bytes() == _entry("1ABC")


# ---------------------------------------------------------------------------
# Not found, and PDB format unavailable
# ---------------------------------------------------------------------------


def test_a_missing_entry_is_not_found_and_not_retried(stub, tmp_path):
    record = fetch(stub, ["9ZZZ"], tmp_path, retries=3).records[0]
    assert record.status == "not_found"
    assert stub.hits["/download/9ZZZ.pdb.gz"] == 1
    assert not os.path.exists(str(tmp_path / "9ZZZ.pdb"))


def test_an_mmcif_only_entry_is_reported_as_such(stub, tmp_path):
    stub.serve("7BIG", {"status": 200, "body": b"data_7BIG\n"}, suffix=".cif")
    record = fetch(stub, ["7BIG"], tmp_path).records[0]
    assert record.status == "pdb_format_unavailable"
    assert "mmCIF" in record.message


# ---------------------------------------------------------------------------
# Retries and transfer integrity
# ---------------------------------------------------------------------------


def test_a_server_error_is_retried_until_it_succeeds(stub, tmp_path):
    stub.serve("1ABC", {"status": 500, "body": b"busy"},
               {"status": 503, "body": b"busy"},
               {"status": 200, "body": _gz(_entry("1ABC"))})
    record = fetch(stub, ["1ABC"], tmp_path, retries=3).records[0]
    assert record.status == "downloaded"
    assert record.attempts == 3


def test_retries_are_bounded(stub, tmp_path):
    stub.serve("1ABC", {"status": 500, "body": b"busy"})
    record = fetch(stub, ["1ABC"], tmp_path, retries=2).records[0]
    assert record.status == "failed"
    assert stub.hits["/download/1ABC.pdb.gz"] == 3
    assert "after 3 attempts" in record.message


def test_429_honours_retry_after(stub, tmp_path):
    stub.serve("1ABC", {"status": 429, "body": b"", "headers": {"Retry-After": "0.3"}},
               {"status": 200, "body": _gz(_entry("1ABC"))})
    started = time.time()
    record = fetch(stub, ["1ABC"], tmp_path).records[0]
    assert record.status == "downloaded" and record.attempts == 2
    assert time.time() - started >= 0.3


def test_a_timeout_is_retried_then_reported(stub, tmp_path):
    stub.serve("1ABC", {"status": 200, "body": _gz(_entry("1ABC")), "delay": 1.0})
    record = fetch(stub, ["1ABC"], tmp_path, timeout=0.2, retries=1).records[0]
    assert record.status == "failed"
    assert record.attempts == 2
    assert not os.path.exists(str(tmp_path / "1ABC.pdb"))


def test_a_truncated_transfer_never_becomes_a_file(stub, tmp_path):
    body = _gz(_entry("1ABC"))
    stub.serve("1ABC", {"status": 200, "body": body[: len(body) // 2], "length": len(body)})
    record = fetch(stub, ["1ABC"], tmp_path, retries=1).records[0]
    assert record.status == "failed"
    assert not os.path.exists(str(tmp_path / "1ABC.pdb"))
    assert _parts(tmp_path) == []


def test_a_gzip_stream_cut_short_is_rejected_even_with_a_matching_length(stub, tmp_path):
    body = _gz(_entry("1ABC"))[:-6]   # the length header agrees; the trailer is gone
    stub.serve("1ABC", {"status": 200, "body": body})
    record = fetch(stub, ["1ABC"], tmp_path, retries=0).records[0]
    assert record.status == "failed"
    assert "gzip" in record.message
    assert not os.path.exists(str(tmp_path / "1ABC.pdb"))


def test_a_gzip_crc_mismatch_is_rejected(stub, tmp_path):
    body = bytearray(_gz(_entry("1ABC")))
    body[-8] ^= 0xFF                    # corrupt the CRC-32 in the trailer
    stub.serve("1ABC", {"status": 200, "body": bytes(body)})
    record = fetch(stub, ["1ABC"], tmp_path, retries=0).records[0]
    assert record.status == "failed" and "gzip" in record.message


def test_data_after_the_gzip_stream_is_rejected(stub, tmp_path):
    stub.serve("1ABC", {"status": 200, "body": _gz(_entry("1ABC")) + b"garbage"})
    record = fetch(stub, ["1ABC"], tmp_path, retries=0).records[0]
    assert record.status == "failed" and "after the gzip stream" in record.message


def test_a_decompressed_file_that_ends_early_is_rejected(stub, tmp_path):
    cut = _entry("1ABC").rsplit(b"END", 1)[0]
    stub.ok("1ABC", cut)
    record = fetch(stub, ["1ABC"], tmp_path).records[0]
    assert record.status == "failed" and "END record" in record.message
    assert not os.path.exists(str(tmp_path / "1ABC.pdb"))


def test_a_file_for_another_entry_is_rejected(stub, tmp_path):
    stub.ok("1ABC", _entry("2XYZ"))
    record = fetch(stub, ["1ABC"], tmp_path).records[0]
    assert record.status == "failed" and "HEADER names entry 2XYZ" in record.message


# ---------------------------------------------------------------------------
# Identifiers, concurrency, transport and dry runs
# ---------------------------------------------------------------------------


def test_ids_are_normalised_deduplicated_and_invalid_ones_reported(stub, tmp_path):
    stub.ok("1ABC")
    report = fetch(stub, ["1abc", "1ABC", " 1aBc ", "notanid", "", "12"], tmp_path)
    statuses = [(r.pdb_id, r.requested, r.status) for r in report.records]
    assert statuses == [
        ("1ABC", "1abc", "downloaded"),
        ("", "notanid", "invalid_id"),
        ("", "12", "invalid_id"),
    ]
    assert stub.hits["/download/1ABC.pdb.gz"] == 1
    assert not report.ok and len(report.failures) == 2


def test_concurrent_downloads_are_bounded_by_workers(stub, tmp_path):
    ids = ["1A%02d" % i for i in range(12)]
    for pdb_id in ids:
        stub.serve(pdb_id, {"status": 200, "body": _gz(_entry(pdb_id)), "delay": 0.15})
    report = fetch(stub, ids, tmp_path, workers=3)
    assert all(r.status == "downloaded" for r in report.records)
    assert [r.pdb_id for r in report.records] == ids      # reported in request order
    assert 1 < stub.max_active <= 3


def test_plain_http_to_a_remote_host_is_refused(tmp_path):
    with pytest.raises(ValueError, match="HTTPS only"):
        R.fetch_structures(["1ABC"], str(tmp_path), base_url="http://example.org/download")


def test_a_redirect_away_from_https_is_refused(stub, tmp_path):
    stub.serve("1ABC", {"status": 302, "body": b"",
                        "headers": {"Location": "http://example.org/1ABC.pdb.gz"}})
    record = fetch(stub, ["1ABC"], tmp_path, retries=0).records[0]
    assert record.status == "failed" and "HTTPS" in record.message


def test_a_dry_run_makes_no_request_and_writes_nothing(stub, tmp_path):
    (tmp_path / "2DEF.pdb").write_bytes(_entry("2DEF"))
    before = sorted(os.listdir(str(tmp_path)))
    report = fetch(stub, ["1ABC", "2DEF"], tmp_path, dry_run=True)
    assert [(r.pdb_id, r.status) for r in report.records] == [
        ("1ABC", "would_download"), ("2DEF", "cached")]
    assert stub.hits == {}
    assert sorted(os.listdir(str(tmp_path))) == before
    missing = tmp_path / "not_created"
    fetch(stub, ["1ABC"], missing, dry_run=True)
    assert not missing.exists()


def test_the_report_csv_has_one_row_per_id(stub, tmp_path):
    stub.ok("1ABC")
    report = fetch(stub, ["1ABC", "9ZZZ"], tmp_path)
    path = report.write_csv(str(tmp_path / "fetch_report.csv"))
    lines = open(path).read().splitlines()
    assert lines[0].split(",") == list(R.REPORT_COLUMNS)
    assert [line.split(",")[2] for line in lines[1:]] == ["downloaded", "not_found"]


def test_bad_arguments_are_refused(tmp_path):
    for kwargs in ({"workers": 0}, {"retries": -1}, {"timeout": 0}):
        with pytest.raises(ValueError):
            R.fetch_structures(["1ABC"], str(tmp_path), dry_run=True, **kwargs)


def test_the_package_exports_the_entry_point():
    assert biotsr.fetch_structures is R.fetch_structures


def test_installed_files_get_ordinary_permissions_not_mkstemp_ones(stub, tmp_path):
    stub.ok("1ABC")
    report = fetch(stub, ["1ABC"], tmp_path)
    report.write_csv(str(tmp_path / "fetch_report.csv"))
    current = os.umask(0o022)
    os.umask(current)
    expected = 0o666 & ~current
    for path in (tmp_path / "1ABC.pdb", tmp_path / ".biotsr-fetch" / "1ABC.json",
                 tmp_path / "fetch_report.csv"):
        assert (path.stat().st_mode & 0o777) == expected, path
