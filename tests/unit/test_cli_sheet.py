"""The built-in wide-sheet layout agrees with the documented column vocabulary."""

from __future__ import annotations

import pytest

from biotsr.cli.errors import CLIError
from biotsr.cli.sheet import builtin_layout, parse_map_options
from biotsr.io.input_columns import DATASET_COLUMNS, HEADER_ALIASES
from biotsr.methods.arguments import accepted_arguments
from biotsr.methods.registry import get_method, list_methods


def test_the_tsr_layout_reads_only_documented_columns():
    shared, methods = builtin_layout("tsr")
    documented = set(DATASET_COLUMNS) | set(HEADER_ALIASES)
    used = set(shared.values())
    for participants in methods.values():
        for fields in participants.values():
            used.update(fields.values())
    assert used <= documented, sorted(used - documented)


def test_every_documented_method_column_is_read_by_its_method():
    """aa_ind_2 is aa_ind on aa_chain_2, carb_atom_selected feeds
    carbohydrate_atom_selection, and so on: the layout and the dictionary agree
    on which method and which chain column each column belongs to."""
    _shared, methods = builtin_layout("tsr")
    for spec in DATASET_COLUMNS.values():
        if spec.method is None:
            continue
        readers = [fields for fields in methods[spec.method].values()
                   if spec.column in fields.values()]
        assert readers, spec.column
        for fields in readers:
            assert fields["chain"] == spec.chain_column, spec.column


def test_every_layout_field_is_an_argument_the_method_reads():
    _shared, methods = builtin_layout("tsr")
    assert set(methods) == set(list_methods())
    for method, participants in methods.items():
        allowed = {"chain"} | set(accepted_arguments(get_method(method)))
        for fields in participants.values():
            assert set(fields) <= allowed, (method, set(fields) - allowed)


def test_a_cross_method_reads_no_narrowing_column_unless_mapped():
    _shared, methods = builtin_layout("tsr")
    for method in ("drug_protein_cross", "carbohydrate_protein_cross",
                   "nucleotide_protein_cross", "protein_protein_cross",
                   "protein_protein_calpha_cross"):
        assert list(methods[method]) == ["1"]
        fields = set(methods[method]["1"])
        assert fields <= {"chain", "chain_2", "entity_name", "entity_id"}, method


def test_aa_group_is_an_annotation_not_the_grouping_modifier():
    shared, methods = builtin_layout("tsr")
    assert shared["label"] == "aa_group"
    for participants in methods.values():
        for fields in participants.values():
            assert "aa_grouping" not in fields and "aa_group" not in fields.values()


def test_the_custom_layout_starts_empty():
    assert builtin_layout("custom") == ({}, {})


@pytest.mark.parametrize("spelling, parsed", [
    ("carbohydrate.chain=glycan", ("carbohydrate", "1", "chain", "glycan")),
    ("aa_ind.2.residue_types=t2", ("aa_ind", "2", "residue_types", "t2")),
])
def test_map_options_parse(spelling, parsed):
    assert parse_map_options([spelling])[0][:4] == parsed


@pytest.mark.parametrize("bad", ["carbohydrate", "carbohydrate.chain", "=x",
                                 "a.b.c.d=x", ".chain=x"])
def test_malformed_map_options_are_refused(bad):
    with pytest.raises(CLIError):
        parse_map_options([bad])


def test_every_cross_side_candidate_is_a_documented_column_its_method_reads():
    """--cross-side chooses among the dictionary's columns, and each feeds an
    argument the cross method accepts; every cross method has two sides."""
    from biotsr.cli.sheet import CROSS_SIDE_COLUMNS
    from biotsr.api import list_cross_methods

    assert set(CROSS_SIDE_COLUMNS) == set(list_cross_methods())
    documented = set(DATASET_COLUMNS)
    for method, sides in CROSS_SIDE_COLUMNS.items():
        assert len(sides) == 2, method
        accepted = set(accepted_arguments(get_method(method)))
        for columns in sides.values():
            assert set(columns) <= documented, (method, set(columns) - documented)
            assert set(columns.values()) <= accepted, (method, columns)


@pytest.mark.parametrize("spelling, parsed", [
    ("drug_protein_cross.protein=aa_ind", {("drug_protein_cross", "protein"): "aa_ind"}),
    ("carbohydrate_protein_cross.carbohydrate=car_atom_selected",
     {("carbohydrate_protein_cross", "carbohydrate"): "carb_atom_selected"}),
    ("protein_protein_cross.2=whole", {("protein_protein_cross", "2"): "whole"}),
])
def test_cross_side_options_parse(spelling, parsed):
    from biotsr.cli.sheet import parse_cross_side_options

    assert parse_cross_side_options([spelling]) == parsed


@pytest.mark.parametrize("bad", ["drug_protein_cross", "drug_protein_cross.protein",
                                 "drug_protein_cross=aa_ind", ".protein=aa_ind",
                                 "drug_protein_cross.protein=aa_selection"])
def test_malformed_cross_side_options_are_refused(bad):
    from biotsr.cli.sheet import parse_cross_side_options

    with pytest.raises(CLIError):
        parse_cross_side_options([bad])
