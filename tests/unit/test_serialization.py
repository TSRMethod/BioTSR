"""Serializer determinism and naming; only the native writer ships."""

from __future__ import annotations

import json
import os

import pytest

import biotsr
from biotsr.io import native

from ..conftest import FIXTURE_PDB  # noqa: F401  (path constants live there)

from ..conftest import _require


@pytest.fixture(scope="module")
def standard_result():
    path = os.path.join(FIXTURE_PDB, "1UBQ.pdb")
    _require(path)
    return biotsr.analyze(path, method="protein", chain="A")


# ---------------------------------------------------------------------------
# Native naming
# ---------------------------------------------------------------------------


def test_native_files_share_one_prefix_and_carry_no_bin_counts(standard_result, tmp_path):
    """``PREFIX.keys``, ``PREFIX.meta.json``, ``PREFIX.triplets`` (ADR-098); the
    bin schemes are in the metadata, not in the name."""
    assert native.default_prefix(standard_result) == "1UBQ_A_protein"
    assert native.keys_filename("p") == "p.keys"
    assert native.triplets_filename("p") == "p.triplets"
    assert native.metadata_filename("p") == "p.meta.json"
    written = native.save(standard_result, str(tmp_path))
    assert sorted(os.path.basename(p) for p in written.values()) == [
        "1UBQ_A_protein.keys", "1UBQ_A_protein.meta.json", "1UBQ_A_protein.triplets"]
    meta = json.loads(open(written["metadata"]).read())
    assert meta["key_space"]["theta_scheme"]["radix"] == 29
    assert meta["key_space"]["distance_scheme"]["radix"] == 35


def test_native_prefixes_make_no_length_assumption(standard_result, tmp_path):
    written = native.save(standard_result, str(tmp_path), prefix="PDB_0000001UBQ_A")
    assert os.path.basename(written["keys"]) == "PDB_0000001UBQ_A.keys"
    assert native.prefix_of(written["keys"]) == "PDB_0000001UBQ_A"
    assert native.prefix_of("A_VERY_LONG_ID_B.meta.json") == "A_VERY_LONG_ID_B"


def test_prefix_of_rejects_nonsense():
    for bad in ("nostem", ".keys", "1UBQ_A.keys_29_35", "1UBQ_A.other"):
        with pytest.raises(ValueError):
            native.prefix_of(bad)


# ---------------------------------------------------------------------------
# Native determinism
# ---------------------------------------------------------------------------


def test_native_keys_are_sorted_numerically(standard_result, tmp_path):
    path = str(tmp_path / "keys")
    native.write_keys(standard_result, path)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh]
    assert header == list(native.KEY_COLUMNS) == ["key", "category", "frequency"]
    keys = [int(r[0]) for r in rows]
    assert keys == sorted(keys)
    assert len(set(keys)) == standard_result.counts.n_unique_keys


def test_native_keys_report_category_and_frequency_explicitly(standard_result, tmp_path):
    """The requirement: category of each key AND its frequency."""
    path = str(tmp_path / "keys")
    native.write_keys(standard_result, path)
    with open(path) as fh:
        fh.readline()
        rows = [line.rstrip("\n").split("\t") for line in fh]
    assert all(len(r) == 3 for r in rows)
    assert {r[1] for r in rows} == {"regular"}, "a run without mirror has one category"
    assert sum(int(r[2]) for r in rows) == standard_result.counts.n_triangles


def test_native_output_is_byte_stable_across_repeated_writes(standard_result, tmp_path):
    a, b = str(tmp_path / "a"), str(tmp_path / "b")
    native.write_keys(standard_result, a)
    native.write_keys(standard_result, b)
    assert open(a, "rb").read() == open(b, "rb").read()

    ta, tb = str(tmp_path / "ta"), str(tmp_path / "tb")
    native.write_triplets(standard_result, ta)
    native.write_triplets(standard_result, tb)
    assert open(ta, "rb").read() == open(tb, "rb").read()


def test_native_triplets_carry_a_header_and_both_angles(standard_result, tmp_path):
    path = str(tmp_path / "t")
    native.write_triplets(standard_result, path)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        first = fh.readline().rstrip("\n").split("\t")
    assert header == list(native.TRIPLET_COLUMNS)
    assert "theta_raw" in header and "theta" in header
    assert len(first) == len(header)
    theta = float(first[header.index("theta")])
    theta_raw = float(first[header.index("theta_raw")])
    assert 0 <= theta <= 90
    assert theta_raw >= theta - 1e-9


def test_native_triplets_include_full_atom_identity(standard_result, tmp_path):
    path = str(tmp_path / "t")
    native.write_triplets(standard_result, path)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        row = fh.readline().rstrip("\n").split("\t")
    for n in (1, 2, 3):
        assert row[header.index("chain%d" % n)] == "A"
        assert row[header.index("atom%d" % n)] == "CA"
        assert row[header.index("residue%d" % n)].isalpha()


def test_streaming_and_in_memory_native_writers_agree(standard_result, tmp_path):
    direct = str(tmp_path / "direct")
    native.write_triplets(standard_result, direct)
    streamed = str(tmp_path / "streamed")
    with native.NativeTripletStreamWriter(streamed) as w:
        for record in standard_result.triplets:
            w.write(record)
    assert open(direct, "rb").read() == open(streamed, "rb").read()


def test_metadata_sidecar_carries_both_identities(standard_result, tmp_path):
    path = str(tmp_path / "meta.json")
    native.write_metadata(standard_result, path)
    meta = json.load(open(path))
    assert meta["key_space_fingerprint"]
    assert meta["sampling_fingerprint"]
    assert meta["key_space"]["angle_convention"].startswith("pi")
    assert meta["counts"]["n_triangles"] == standard_result.counts.n_triangles
    assert meta["source"]["sha256"]
    assert meta["structure_id"] == "1UBQ"


def test_save_writes_keys_metadata_and_triplets(standard_result, tmp_path):
    written = native.save(standard_result, str(tmp_path))
    assert set(written) == {"keys", "metadata", "triplets"}
    for path in written.values():
        assert os.path.getsize(path) > 0


# ---------------------------------------------------------------------------
# No legacy writer
# ---------------------------------------------------------------------------


def test_no_legacy_writer_ships(standard_result):
    """The historical ``.3Dkeys_*`` writers were removed; native is the only format."""
    import importlib

    import biotsr.io as io_pkg

    for name in ("legacy_text", "legacy_sse", "legacy_atomic", "legacy_cross"):
        assert not hasattr(io_pkg, name), name
        with pytest.raises(ImportError):
            importlib.import_module("biotsr.io." + name)
    assert "3Dkeys" not in native.default_prefix(standard_result)
