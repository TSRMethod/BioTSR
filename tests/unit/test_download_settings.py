"""Download settings are validated identically by every entry point.

Regression tests for the rc1 defect in which ``biotsr batch --download-missing``
silently replaced ``--timeout 0`` with 30 s and ``--download-workers 0`` with 4,
while ``biotsr fetch`` refused them -- and in which a non-finite timeout
(``nan``) was accepted everywhere.
"""

from __future__ import annotations

import math

import pytest

import biotsr
from biotsr import batch as biotsr_batch
from biotsr.cli.main import main
from biotsr.retrieval import MAX_WORKERS, validate_download_settings

from ..conftest import structure


@pytest.mark.parametrize("workers", [0, -1, MAX_WORKERS + 1, True, 2.5, "4"])
def test_invalid_worker_counts_are_refused(workers):
    with pytest.raises(ValueError, match="download workers"):
        validate_download_settings(workers=workers, timeout=30, retries=3)


@pytest.mark.parametrize("timeout", [0, -1, math.nan, math.inf, "x", True])
def test_invalid_timeouts_are_refused(timeout):
    with pytest.raises(ValueError, match="timeout"):
        validate_download_settings(workers=4, timeout=timeout, retries=3)


@pytest.mark.parametrize("retries", [-1, 1.5, True, "3"])
def test_invalid_retry_counts_are_refused(retries):
    with pytest.raises(ValueError, match="retries"):
        validate_download_settings(workers=4, timeout=30, retries=retries)


@pytest.mark.parametrize("workers,timeout,retries", [
    (1, 0.5, 0), (MAX_WORKERS, 30, 3), (4, 30.0, 10), (4.0, 1, 2.0),
])
def test_valid_settings_pass(workers, timeout, retries):
    validate_download_settings(workers=workers, timeout=timeout, retries=retries)


def test_fetch_structures_applies_the_same_rules(tmp_path):
    for kwargs in ({"workers": 0}, {"workers": MAX_WORKERS + 1}, {"timeout": 0},
                   {"timeout": math.nan}, {"retries": -1}):
        with pytest.raises(ValueError):
            biotsr.fetch_structures(["1UBQ"], str(tmp_path), dry_run=True, **kwargs)


@pytest.mark.parametrize("flag,value,named", [
    ("--timeout", "0", "--timeout"),
    ("--timeout", "nan", "--timeout"),
    ("--timeout", "-5", "--timeout"),
    ("--download-workers", "0", "--download-workers"),
    ("--download-workers", "17", "--download-workers"),
    ("--retries", "-1", "--retries"),
])
@pytest.mark.parametrize("command", ["batch", "fetch", "analyze"])
def test_every_command_refuses_the_same_values(tmp_path, capsys, flag, value, named,
                                               command):
    manifest = tmp_path / "m.csv"
    manifest.write_text("pdb_id,method,chain\n1UBQ,protein,A\n")
    pdb_dir = tmp_path / "pdb"
    pdb_dir.mkdir()
    pdb_dir.joinpath("1UBQ.pdb").write_bytes(open(structure("1UBQ"), "rb").read())
    if command == "batch":
        argv = ["batch", "--manifest", str(manifest), "--pdb-dir", str(pdb_dir),
                "--download-missing", "--output", str(tmp_path / "out")]
    elif command == "analyze":
        argv = ["analyze", "common-keys", "--manifest", str(manifest), "--pdb-dir",
                str(pdb_dir), "--download-missing", "--output", str(tmp_path / "out")]
    else:
        argv = ["fetch", "1UBQ", "--pdb-dir", str(pdb_dir), "--dry-run"]
    code = main(argv + [flag, value])
    assert code == 2
    err = capsys.readouterr().err
    assert named in err
    assert not (tmp_path / "out").exists()


def test_the_python_batch_api_refuses_them_too(tmp_path):
    manifest = tmp_path / "m.csv"
    manifest.write_text("pdb_id,method,chain\n1UBQ,protein,A\n")
    for kwargs in ({"timeout": 0}, {"download_workers": 0}, {"retries": -1},
                   {"if_exists": "replace"}, {"on_error": "ignore"}):
        with pytest.raises(biotsr_batch.BatchError):
            biotsr_batch.run(str(manifest), str(tmp_path / "out"),
                             pdb_dir=str(tmp_path), **kwargs)
    assert not (tmp_path / "out").exists()


def test_the_timeout_help_states_its_actual_scope(capsys):
    with pytest.raises(SystemExit):
        main(["fetch", "--help"])
    text = " ".join(capsys.readouterr().out.split())
    assert "It limits each wait, not a whole download" in text
    assert "however slowly, is never timed out" in text
