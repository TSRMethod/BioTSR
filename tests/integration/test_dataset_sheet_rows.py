"""Dataset-sheet rows: the historical one-entry-per-row format, end to end.

The rows come from ``tests/data/manifests/dataset_sheet.csv``, a synthetic sheet
written for these tests from the public structures in ``tests/data/structures``
(its README says how).  It keeps the format's real quirks: the
``car_atom_selected`` header spelling (documented as ``carb_atom_selected``),
``comment_*`` columns, and both residue spellings in one cell.

Nothing here invents a meaning for a blank or undocumented field, and the
header mismatch is resolved through an explicit alias, never by rewriting data.

A real dataset sheet kept outside the repository can be checked against the
same vocabulary by pointing ``BIOTSR_EXTERNAL_DATASET_SHEET`` at it; without that
variable the check is skipped as an opt-in check that was not requested.
"""

from __future__ import annotations

import csv
import math
import os

import pytest

import biotsr
from biotsr.io.input_columns import (
    DATASET_COLUMNS,
    HEADER_ALIASES,
    canonical_column,
    method_for_column,
)
from biotsr.methods.registry import list_methods

from ..conftest import STRUCTURES, TESTS_DIR

SHEET = os.path.join(TESTS_DIR, "data", "manifests", "dataset_sheet.csv")
EXTERNAL_SHEET_ENV = "BIOTSR_EXTERNAL_DATASET_SHEET"
#: Skip reason for the opt-in check over a sheet kept outside the repository.
OPT_IN_NOT_REQUESTED = "opt-in check not requested"


def _read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


ROWS = _read(SHEET)
DRUG_ROWS = [r for r in ROWS if r["drug_name"].strip()]
CARB_ROWS = [r for r in ROWS if r["carb_name"].strip()]
NT_ROWS = [r for r in ROWS if r["nt_name"].strip()]


def _structure(pdb_id):
    return os.path.join(STRUCTURES, "%s.pdb" % pdb_id)


def _row_id(*columns):
    return lambda r: "_".join(r[c] for c in columns)


def _atoms(cell):
    return [a for a in cell.split(";") if a.strip()]


def _undocumented(headers):
    unknown = []
    for header in headers:
        if header.startswith("comment_"):
            continue
        try:
            canonical_column(header)
        except KeyError:
            unknown.append(header)
    return unknown


# ===========================================================================
# The column vocabulary
# ===========================================================================


def test_the_sheet_exercises_every_documented_column():
    """The fixture covers the whole vocabulary, so the checks below cover it too."""
    canonical = {canonical_column(h) for h in ROWS[0] if not h.startswith("comment_")}
    assert canonical == set(DATASET_COLUMNS)


def test_every_header_is_documented_or_an_explicit_alias():
    """No column is silently ignored and none is silently renamed."""
    assert _undocumented(ROWS[0].keys()) == []


def test_an_external_dataset_sheet_uses_only_documented_columns():
    """The same check over a real sheet kept outside the repository."""
    path = os.environ.get(EXTERNAL_SHEET_ENV)
    if not path:
        pytest.skip("%s: an external dataset sheet (set %s)"
                    % (OPT_IN_NOT_REQUESTED, EXTERNAL_SHEET_ENV))
    rows = _read(path)
    assert rows, path
    assert _undocumented(rows[0].keys()) == []
    from biotsr.domain.identity import StructureId

    for row in rows:
        assert StructureId.parse(row["protein"]).value == row["protein"].strip().upper()


def test_the_carb_atom_selected_header_mismatch_is_an_alias_not_a_rename():
    headers = set(ROWS[0])
    assert "car_atom_selected" in headers          # the spelling sheets carry
    assert "carb_atom_selected" not in headers     # the documented name
    assert HEADER_ALIASES["car_atom_selected"] == "carb_atom_selected"
    assert canonical_column("car_atom_selected") == "carb_atom_selected"
    assert method_for_column("car_atom_selected") == "carbohydrate_atom_selection"


def test_an_unknown_header_raises_rather_than_being_guessed_at():
    with pytest.raises(KeyError, match="unknown dataset column"):
        canonical_column("some_column_nobody_documented")


def test_aa_group_is_not_bound_to_the_aa_grouping_modifier():
    spec = DATASET_COLUMNS["aa_group"]
    assert spec.method is None
    assert "annotation" in spec.meaning
    assert spec.unimplemented_reason


def test_every_method_a_column_names_is_registered():
    registered = set(list_methods())
    for spec in DATASET_COLUMNS.values():
        if spec.method is not None:
            assert spec.method in registered, spec.column


def test_no_column_maps_to_a_cross_method():
    """Cross keys pair two sides; no single sheet column describes one."""
    for spec in DATASET_COLUMNS.values():
        assert spec.method is None or "cross" not in spec.method


# ===========================================================================
# The rows themselves
# ===========================================================================


def test_the_sheet_describes_all_three_entity_families():
    assert (len(DRUG_ROWS), len(CARB_ROWS), len(NT_ROWS)) == (2, 2, 1)
    for row in ROWS:
        assert os.path.exists(_structure(row["protein"])), row["protein"]


@pytest.mark.parametrize("row", DRUG_ROWS, ids=_row_id("protein", "drug_name"))
def test_a_drug_row_runs_whole_and_with_its_atom_selection(row):
    path = _structure(row["protein"])
    kwargs = dict(chain=row["drug_chain"], entity_name=row["drug_name"],
                  entity_id=row["drug_id"], keep_triplets=False)
    whole = biotsr.analyze(path, method="drug", **kwargs)
    assert whole.counts.n_points > 3
    wanted = _atoms(row["drug_atom_selected"])
    if not wanted:
        return
    selected = biotsr.analyze(path, method=method_for_column("drug_atom_selected"),
                              atoms=row["drug_atom_selected"], **kwargs)
    assert selected.counts.n_points == len(wanted)
    assert selected.counts.n_triangles == math.comb(len(wanted), 3)
    assert selected.selection.atoms_missing == ()
    assert selected.key_space.label_scheme.radix == 112
    assert selected.key_space.distance_scheme.radix == 18


@pytest.mark.parametrize("row", CARB_ROWS, ids=_row_id("protein", "carb_name", "carb_id"))
def test_a_carbohydrate_row_runs_monomer_and_selected_atoms(row):
    path = _structure(row["protein"])
    kwargs = dict(chain=row["carb_chain"], entity_name=row["carb_name"],
                  entity_id=row["carb_id"], keep_triplets=False)
    monomer = biotsr.analyze(path, method="carbohydrate", **kwargs)
    assert monomer.counts.n_points > 3
    assert monomer.selection.residues_used == (
        "%s_%s" % (row["carb_id"], row["carb_name"]),
    )
    # The header is `car_atom_selected`; the documented name is
    # `carb_atom_selected`. The alias resolves the *meaning*; the data keeps its
    # own spelling, which is why the row is still read by its real header.
    selected = biotsr.analyze(path, method=method_for_column("car_atom_selected"),
                              atoms=row["car_atom_selected"], **kwargs)
    assert selected.counts.n_points == len(_atoms(row["car_atom_selected"]))
    assert selected.selection.atoms_missing == ()


@pytest.mark.parametrize("row", NT_ROWS, ids=_row_id("protein", "nt_chain", "nt_name", "nt_id"))
def test_a_nucleotide_row_runs_residue_and_selected_atoms(row):
    path = _structure(row["protein"])
    kwargs = dict(chain=row["nt_chain"], entity_name=row["nt_name"],
                  entity_id=row["nt_id"], keep_triplets=False)
    residue = biotsr.analyze(path, method="nucleotide", **kwargs)
    assert residue.counts.n_points > 3
    selected = biotsr.analyze(path, method=method_for_column("nt_atom_selected"),
                              atoms=row["nt_atom_selected"], **kwargs)
    assert selected.counts.n_points == len(_atoms(row["nt_atom_selected"]))
    assert selected.selection.atoms_missing == ()


def _kinase_row():
    return next(r for r in ROWS if r["protein"] == "6G9M")


def test_the_aa_selection_cell_stays_in_the_calpha_key_space():
    row = _kinase_row()
    result = biotsr.analyze(
        _structure("6G9M"), method=method_for_column("aa_selection"),
        chain=row["aa_chain"], residues=row["aa_selection"], keep_triplets=False,
    )
    assert result.counts.n_points == 5            # 26, 27, 28, 192, 207
    assert result.key_space.label_scheme.radix == 30
    assert result.key_space.distance_scheme.radix == 35


def test_the_aa_atom_selection_cell_stays_in_the_atomic_key_space():
    row = _kinase_row()
    result = biotsr.analyze(
        _structure("6G9M"), method=method_for_column("aa_atom_selection"),
        chain=row["aa_chain"], atoms=row["aa_atom_selection"], keep_triplets=False,
    )
    assert result.key_space.label_scheme.radix == 112
    assert result.key_space.distance_scheme.radix == 18
    assert result.counts.n_points == len(_atoms(row["aa_atom_selection"])) == 5
    assert result.selection.atoms_missing == ()


def test_a_cell_mixing_both_residue_spellings_selects_both_residues():
    row = _kinase_row()
    result = biotsr.analyze(
        _structure("6G9M"), method=method_for_column("aa_ind_specific"),
        chain=row["aa_chain"], residues=row["aa_ind_specific"], keep_triplets=False,
    )
    assert set(result.selection.residues_used) == {"207_LYS", "266_SER"}


def test_structure_ids_are_never_assumed_to_be_four_characters():
    from biotsr.domain.identity import StructureId

    long_id = StructureId.parse("SOME-VERY-LONG-STRUCTURE-ID-2026")
    assert long_id.value == "SOME-VERY-LONG-STRUCTURE-ID-2026"
    for row in ROWS:
        assert StructureId.parse(row["protein"]).value == row["protein"].upper()
