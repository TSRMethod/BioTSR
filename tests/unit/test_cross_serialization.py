"""Cross serialization: a standard result carries enough provenance to identify
both sides of every triangle, and its key file is sorted numerically.
"""

from __future__ import annotations

import os

import pytest

from biotsr import analyze_cross
from biotsr.cross.pools import SIDE_A, SIDE_B
from biotsr.io import native

from ..conftest import structure as tracked_structure

HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = HERE  # the package project is the repository root


def fixture(pdb_id):
    return tracked_structure(pdb_id)


@pytest.fixture(scope="module")
def standard_result():
    return analyze_cross(fixture("2AC0"), method="nucleotide_protein_cross",
                         chain="E", chain_2="B", entity_name="DG", entity_id="8")


# ===========================================================================
# Provenance survives serialization
# ===========================================================================


def test_the_native_triplet_line_identifies_both_sides_of_every_triangle(
    standard_result, tmp_path
):
    path = native.write_triplets(standard_result, str(tmp_path / "t"))
    lines = open(path).read().splitlines()
    header = lines[0].split("\t")
    for row in (dict(zip(header, line.split("\t"))) for line in lines[1:]):
        sides = {row["side1"], row["side2"], row["side3"]}
        assert sides == {SIDE_A, SIDE_B}
        for n in (1, 2, 3):
            assert row["structure%d" % n] == "2AC0"
            assert row["residue%d" % n]
            assert row["chain%d" % n]
            assert row["seq%d" % n]
            assert row["atom%d" % n]


def test_nothing_is_lost_between_the_record_and_the_line(standard_result, tmp_path):
    path = native.write_triplets(standard_result, str(tmp_path / "t"))
    lines = open(path).read().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    assert len(rows) == len(standard_result.triplets)
    for row, record in zip(rows, standard_result.triplets):
        assert int(row["key"]) == record.key
        for n, (identity, side) in enumerate(
            zip(record.identities, record.sides), start=1
        ):
            assert row["side%d" % n] == side
            assert row["structure%d" % n] == identity.structure.value
            assert row["chain%d" % n] == identity.chain.value
            assert row["residue%d" % n] == identity.residue_name
            assert row["seq%d" % n] == identity.residue.legacy_seq_label
            assert row["atom%d" % n] == identity.atom_name


def test_the_standard_key_file_is_sorted_numerically(standard_result, tmp_path):
    path = native.write_keys(standard_result, str(tmp_path / "k"))
    keys = [int(line.split("\t")[0]) for line in open(path).read().splitlines()[1:]]
    assert keys == sorted(keys)
    assert len(keys) == len(set(keys))
