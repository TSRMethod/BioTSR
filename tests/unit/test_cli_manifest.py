"""The canonical dataset manifest.

One row is one run.  A blank cell is "not set", never a default.  An unknown
header is an error, never an ignored column -- silently dropping a column a user
believed was driving the calculation is the single worst thing a batch driver can
do, and it is what the legacy ``--*_column`` override flags made easy.
"""

from __future__ import annotations

import os

import pytest

from biotsr.cli.errors import CLIError
from biotsr.cli.manifest import (
    CANONICAL_COLUMNS,
    COLUMN_ALIASES,
    read_manifest,
    resolve_structure_path,
)
from biotsr.cli.spec import build_spec


def write(tmp_path, text, name="m.csv"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


# ---------------------------------------------------------------------------
# Shape
# ---------------------------------------------------------------------------


def test_one_row_is_one_run(tmp_path):
    path = write(tmp_path, "pdb_id,method,chain\n1UBQ,protein,A\n1TUP,protein,A\n")
    rows = read_manifest(path)
    assert len(rows) == 2
    assert rows[0]["pdb_id"] == "1UBQ"
    assert rows[1]["pdb_id"] == "1TUP"


def test_a_blank_cell_is_absent_not_a_default(tmp_path):
    path = write(tmp_path, "pdb_id,method,chain,mirror,size_filter\n1UBQ,protein,A,,\n")
    row = read_manifest(path)[0]
    assert "mirror" not in row
    assert "size_filter" not in row
    spec = build_spec(dict(row, pdb="x.pdb"))
    assert spec.modifiers == {}


def test_blank_lines_are_skipped(tmp_path):
    path = write(tmp_path, "pdb_id,method,chain\n1UBQ,protein,A\n\n,,\n1TUP,protein,A\n")
    assert len(read_manifest(path)) == 2


def test_comment_columns_are_ignored(tmp_path):
    """Historical dataset sheets carry several (``comment_1`` ...)."""
    path = write(
        tmp_path,
        "pdb_id,method,chain,comment_1,comment_2\n1UBQ,protein,A,hello,world\n",
    )
    row = read_manifest(path)[0]
    assert "comment_1" not in row


def test_an_unknown_column_is_an_error_naming_the_known_ones(tmp_path):
    path = write(tmp_path, "pdb_id,method,chain,wobble\n1UBQ,protein,A,x\n")
    with pytest.raises(CLIError) as excinfo:
        read_manifest(path)
    assert "wobble" in excinfo.value.message
    assert "pdb_id" in excinfo.value.hint


def test_an_empty_manifest_is_an_error(tmp_path):
    with pytest.raises(CLIError):
        read_manifest(write(tmp_path, ""))
    with pytest.raises(CLIError):
        read_manifest(write(tmp_path, "pdb_id,method,chain\n"))


def test_a_missing_manifest_is_an_error(tmp_path):
    with pytest.raises(CLIError) as excinfo:
        read_manifest(str(tmp_path / "nope.csv"))
    assert "not found" in excinfo.value.message


# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("text", ["true", "TRUE", "True", "yes", "y", "1", "on"])
def test_truthy_spellings(tmp_path, text):
    path = write(tmp_path, "pdb_id,method,chain,mirror\n1UBQ,protein,A,%s\n" % text)
    assert read_manifest(path)[0]["mirror"] is True


@pytest.mark.parametrize("text", ["false", "FALSE", "no", "n", "0", "off"])
def test_falsy_spellings(tmp_path, text):
    path = write(tmp_path, "pdb_id,method,chain,mirror\n1UBQ,protein,A,%s\n" % text)
    assert read_manifest(path)[0]["mirror"] is False


def test_a_non_truth_value_is_an_error_naming_the_row(tmp_path):
    path = write(tmp_path, "pdb_id,method,chain,mirror\n1UBQ,protein,A,maybe\n")
    with pytest.raises(CLIError) as excinfo:
        read_manifest(path)
    assert "line 2" in excinfo.value.message
    assert "maybe" in excinfo.value.message


def test_numbers_are_parsed_as_numbers(tmp_path):
    path = write(
        tmp_path, "pdb_id,method,chain,size_filter\n1UBQ,protein,A,25.0\n"
    )
    assert read_manifest(path)[0]["size_filter"] == 25.0


def test_a_non_number_is_an_error_naming_the_row(tmp_path):
    path = write(
        tmp_path, "pdb_id,method,chain,size_filter\n1UBQ,protein,A,big\n"
    )
    with pytest.raises(CLIError) as excinfo:
        read_manifest(path)
    assert "line 2" in excinfo.value.message
    assert "not a number" in excinfo.value.message


# ---------------------------------------------------------------------------
# Aliases
# ---------------------------------------------------------------------------


def test_dataset_sheet_spellings_resolve(tmp_path):
    path = write(
        tmp_path,
        "protein,aa_chain,method,aa_group\n1UBQ,A,protein,ubiquitin\n",
    )
    row = read_manifest(path)[0]
    assert row["pdb_id"] == "1UBQ"
    assert row["chain"] == "A"
    assert row["label"] == "ubiquitin"


def test_the_known_header_mismatch_resolves_without_renaming_the_data(tmp_path):
    """``car_atom_selected`` in the sheets, ``carb_atom_selected`` in their definitions."""
    assert COLUMN_ALIASES["car_atom_selected"] == "atoms"
    assert COLUMN_ALIASES["carb_atom_selected"] == "atoms"
    path = write(
        tmp_path,
        "pdb_id,method,chain,entity_name,entity_id,car_atom_selected\n"
        "5WT9,carbohydrate_atom_selection,A,NAG,1,C1;C2;O5\n",
    )
    assert read_manifest(path)[0]["atoms"] == "C1;C2;O5"


def test_two_aliases_of_the_same_column_in_one_row_is_an_error(tmp_path):
    """A wide acceptance sheet is not a manifest, and this is how it says so."""
    path = write(
        tmp_path,
        "protein,aa_chain,drug_chain,method\n1UBQ,A,A,protein\n",
    )
    with pytest.raises(CLIError) as excinfo:
        read_manifest(path)
    assert "aa_chain" in excinfo.value.message
    assert "drug_chain" in excinfo.value.message
    assert "one run" in excinfo.value.hint


def test_an_empty_alias_cell_does_not_collide(tmp_path):
    """The wide sheet reshaped by hand: the irrelevant chain column is blank."""
    path = write(
        tmp_path,
        "protein,aa_chain,drug_chain,method\n1UBQ,A,,protein\n",
    )
    assert read_manifest(path)[0]["chain"] == "A"


def test_every_alias_targets_a_real_column():
    for alias, canonical in COLUMN_ALIASES.items():
        assert canonical in CANONICAL_COLUMNS, (alias, canonical)


def test_canonical_columns_cover_the_whole_flag_vocabulary():
    for name in ("method", "chain", "chain_2", "sse", "mirror",
                 "aa_grouping", "size_filter", "contact_cutoff", "contact_filter",
                 "entity_name", "entity_id", "entity_names", "atoms", "residues",
                 "residue_types", "residue", "residue_name", "label", "sample_id"):
        assert name in CANONICAL_COLUMNS, name
    assert "profile" not in CANONICAL_COLUMNS


# ---------------------------------------------------------------------------
# Defaults from the command line
# ---------------------------------------------------------------------------


def test_command_line_defaults_apply_where_a_row_is_silent(tmp_path):
    path = write(tmp_path, "pdb_id,chain\n1UBQ,A\n")
    row = read_manifest(path, defaults={"method": "protein"})[0]
    assert row["method"] == "protein"


def test_a_row_overrides_a_command_line_default(tmp_path):
    path = write(tmp_path, "pdb_id,chain,method\n1UBQ,A,aa_selection\n")
    row = read_manifest(path, defaults={"method": "protein"})[0]
    assert row["method"] == "aa_selection"


# ---------------------------------------------------------------------------
# Structure resolution
# ---------------------------------------------------------------------------


def test_pdb_path_is_used_as_given(tmp_path):
    path, structure_id = resolve_structure_path(
        {"pdb_path": "/data/1UBQ.pdb"}, "", "row"
    )
    assert path == "/data/1UBQ.pdb"
    assert structure_id is None


def test_pdb_id_is_resolved_against_the_directory(tmp_path):
    (tmp_path / "1UBQ.pdb").write_text("")
    path, structure_id = resolve_structure_path(
        {"pdb_id": "1UBQ"}, str(tmp_path), "row"
    )
    assert path == str(tmp_path / "1UBQ.pdb")
    assert structure_id == "1UBQ"


def test_a_lower_case_pdb_id_finds_an_upper_case_file(tmp_path):
    """The acceptance-sheet CSV writes ``1ubq``; RCSB downloads write ``1UBQ.pdb``.

    ``os.path.samefile`` rather than a string comparison, because a
    case-insensitive filesystem resolves the first candidate it is handed and the
    test is about finding the right *file*, not about which spelling won.
    """
    (tmp_path / "1UBQ.pdb").write_text("")
    path, _structure_id = resolve_structure_path(
        {"pdb_id": "1ubq"}, str(tmp_path), "row"
    )
    assert os.path.samefile(path, str(tmp_path / "1UBQ.pdb"))


def test_an_ent_file_is_found_too(tmp_path):
    (tmp_path / "1UBQ.ent").write_text("")
    path, _structure_id = resolve_structure_path(
        {"pdb_id": "1UBQ"}, str(tmp_path), "row"
    )
    assert os.path.samefile(path, str(tmp_path / "1UBQ.ent"))


def test_pdb_path_wins_but_pdb_id_still_names_the_structure():
    path, structure_id = resolve_structure_path(
        {"pdb_path": "/data/download_4f2q.pdb", "pdb_id": "4F2Q"}, "", "row"
    )
    assert path == "/data/download_4f2q.pdb"
    assert structure_id == "4F2Q"


def test_a_pdb_id_with_no_directory_says_so():
    with pytest.raises(CLIError) as excinfo:
        resolve_structure_path({"pdb_id": "1UBQ"}, "", "m.csv line 2")
    assert "--pdb-dir" in excinfo.value.hint


def test_a_missing_structure_names_what_was_looked_for(tmp_path):
    with pytest.raises(CLIError) as excinfo:
        resolve_structure_path({"pdb_id": "9XYZ"}, str(tmp_path), "m.csv line 2")
    assert "9XYZ" in excinfo.value.message
    assert "m.csv line 2" in excinfo.value.message


def test_a_row_with_no_structure_at_all_says_so():
    with pytest.raises(CLIError) as excinfo:
        resolve_structure_path({"chain": "A"}, "/data", "m.csv line 2")
    assert "no structure" in excinfo.value.message


# ---------------------------------------------------------------------------
# A manifest row and a command line are validated identically
# ---------------------------------------------------------------------------


def test_a_manifest_row_goes_through_the_same_validation(tmp_path):
    path = write(tmp_path, "pdb_id,method,chain,sse\n3PTB,drug,A,true\n")
    row = read_manifest(path)[0]
    with pytest.raises(CLIError) as excinfo:
        build_spec(dict(row, pdb="x.pdb"), where=row["_where"])
    assert "does not support 'sse'" in excinfo.value.message
    assert "line 2" in excinfo.value.message
