"""Manifest reading: nothing a manifest says is silently dropped.

Regression tests for two rc1 defects: cells beyond the header were silently
dropped (``1UBQ,A,1,2,3,4`` under a ``residues`` column used only ``1``), and a
UTF-8 byte-order mark -- what spreadsheet programs write -- made the first header
"unrecognised" with an invisible character.
"""

from __future__ import annotations

import pytest

from biotsr.cli.errors import CLIError
from biotsr.cli.manifest import read_id_column, read_manifest


def _write(tmp_path, text, name="m.csv", bom=False):
    path = tmp_path / name
    data = text.encode("utf-8")
    path.write_bytes((b"\xef\xbb\xbf" if bom else b"") + data)
    return str(path)


def test_a_byte_order_mark_is_accepted(tmp_path):
    path = _write(tmp_path, "pdb_path,method,chain\nx.pdb,protein,A\n", bom=True)
    rows = read_manifest(path)
    assert rows[0]["pdb_path"] == "x.pdb"
    assert rows[0]["chain"] == "A"


def test_a_byte_order_mark_is_accepted_when_fetch_reads_ids(tmp_path):
    path = _write(tmp_path, "protein,anything\n1UBQ,x\n", bom=True)
    assert read_id_column(path) == [(2, "1UBQ")]


def test_surplus_cells_are_refused_with_the_row_and_the_remedy(tmp_path):
    path = _write(tmp_path,
                  "pdb_path,method,chain,residues\nx.pdb,aa_selection,A,1,2,3,4\n")
    with pytest.raises(CLIError) as excinfo:
        read_manifest(path)
    assert "line 2" in excinfo.value.message
    assert "'2'" in excinfo.value.message
    assert "no header" in excinfo.value.message
    assert "';'" in excinfo.value.hint or "quote" in excinfo.value.hint


def test_a_quoted_list_in_one_cell_is_one_value(tmp_path):
    path = _write(tmp_path,
                  'pdb_path,method,chain,residues\nx.pdb,aa_selection,A,"1,2,3"\n')
    assert read_manifest(path)[0]["residues"] == "1,2,3"


def test_trailing_empty_cells_are_harmless(tmp_path):
    path = _write(tmp_path, "pdb_path,method,chain\nx.pdb,protein,A,,\n")
    assert read_manifest(path)[0]["chain"] == "A"


def test_a_value_under_a_blank_header_is_refused(tmp_path):
    path = _write(tmp_path, "pdb_path,,chain\nx.pdb,protein,A\n")
    with pytest.raises(CLIError, match="no header"):
        read_manifest(path, defaults={"method": "protein"})


def test_a_duplicated_header_is_refused(tmp_path):
    path = _write(tmp_path, "pdb_path,chain,chain\nx.pdb,A,\n")
    with pytest.raises(CLIError, match="more than one column named 'chain'"):
        read_manifest(path)


def test_the_profile_column_is_refused_with_a_reason(tmp_path):
    path = _write(tmp_path, "pdb_path,method,chain,profile\nx.pdb,protein,A,standard\n")
    with pytest.raises(CLIError) as excinfo:
        read_manifest(path)
    assert "'profile' column" in excinfo.value.message
    assert "standard profile" in excinfo.value.hint
    # ... and it can be dropped explicitly.
    rows = read_manifest(path, ignore_columns=["profile"])
    assert "profile" not in rows[0]


def test_aliases_mapping_and_metadata_behave_as_documented(tmp_path):
    path = _write(tmp_path,
                  "protein,aa_chain,group,group1,comment_x\n1UBQ,A,g1,meta,ignored\n")
    rows = read_manifest(path, columns={"label": "group"}, keep_columns=["group1"],
                         defaults={"method": "protein"})
    row = rows[0]
    assert row["pdb_id"] == "1UBQ"
    assert row["chain"] == "A"
    assert row["label"] == "g1"
    assert row["_metadata"] == {"group1": "meta"}
    assert row["_chain_source"] == "row"


def test_two_aliases_of_one_field_in_a_row_are_refused(tmp_path):
    path = _write(tmp_path, "protein,aa_chain,drug_chain\n1UBQ,A,B\n")
    with pytest.raises(CLIError, match="both give 'chain'"):
        read_manifest(path, defaults={"method": "drug"})


def test_a_row_without_a_chain_is_never_given_one(tmp_path):
    path = _write(tmp_path, "protein,aa_chain\n1UBQ,\n")
    rows = read_manifest(path, defaults={"method": "protein"})
    assert "chain" not in rows[0]
    assert rows[0]["_chain_source"] == ""


def test_an_unrecognised_column_is_still_an_error(tmp_path):
    path = _write(tmp_path, "pdb_path,method,chain,size_fliter\nx.pdb,protein,A,25\n")
    with pytest.raises(CLIError) as excinfo:
        read_manifest(path)
    assert "size_fliter" in excinfo.value.message
    assert "size_filter" in excinfo.value.hint
