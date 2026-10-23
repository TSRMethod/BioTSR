"""Residue-level assignment: the conventions, and what they do to real records."""

from __future__ import annotations

import pytest

from biotsr.api import load_structure
from biotsr.domain.errors import InvalidRequest
from biotsr.domain.identity import AtomIdentity, ChainId, StructureId
from biotsr.domain.point import Point
from biotsr.sse import (
    STRAND_IDENTITIES,
    STRAND_IDENTITY_SHEET_ID_AND_SERIAL,
    build_assignment_map,
)
from biotsr.sse.assignment import SSEAnnotator
from biotsr.sse.states import HELIX, NONE, SHEET, SSEElementId
from biotsr.structure.backends.pdb_legacy import parse_pdb_legacy
from biotsr.structure.model import AtomRecord, SecondaryStructureRecord, Structure

#: Chain A's residues, in file order, unless a test lays out its own.
_RESIDUES = [(n, "") for n in range(1, 41)]


def _identity(seq, icode="", chain="A", name="ALA"):
    return AtomIdentity(
        structure=StructureId("TEST"),
        model=1,
        chain=ChainId(chain),
        residue_name=name,
        residue_seq=seq,
        insertion_code=icode,
        atom_name="CA",
    )


def _structure(records, residues=_RESIDUES, chain="A"):
    atoms = tuple(
        AtomRecord(identity=_identity(seq, icode, chain=chain),
                   coord=Point(float(n), 0.0, 0.0), hetero=False, element="C",
                   file_order=n, calpha_name="CA")
        for n, (seq, icode) in enumerate(residues))
    return Structure(
        structure_id=StructureId("TEST"),
        atoms=atoms,
        secondary_structure=tuple(records),
        parser_backend="test",
    )


def _record(kind, serial, start, end, chain="A", group_id="X", end_chain=None,
            start_icode="", end_icode=""):
    return SecondaryStructureRecord(
        kind=kind,
        serial=serial,
        group_id=group_id,
        chain=ChainId(chain),
        end_chain=ChainId(end_chain if end_chain is not None else chain),
        start_seq=start,
        end_seq=end,
        start_icode=start_icode,
        end_icode=end_icode,
    )


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def test_helix_and_sheet_records_are_parsed_at_their_spec_offsets(ubq_pdb):
    st = load_structure(ubq_pdb)
    kinds = [r.kind for r in st.secondary_structure]
    assert kinds == ["HELIX", "HELIX"] + ["SHEET"] * 5
    helix1 = st.secondary_structure[0]
    assert (helix1.serial, helix1.group_id) == ("1", "H1")
    assert (helix1.start_seq, helix1.end_seq) == (23, 34)
    assert helix1.chain == ChainId("A") == helix1.end_chain


def test_the_sheet_serial_is_the_strand_number_not_the_sheet_id(ubq_pdb):
    """Both halves of a β-strand's identity are parsed and kept.

    1UBQ has a single physical sheet, ``BET``, with five strands. Five distinct
    strands is the *correct* reading -- categories distinguish individual
    β-strands, so a triangle spanning three of them is ``9_3b3`` and that is
    right, inside one β-sheet (ADR-063). With one sheet, ``serNum`` alone happens
    to be enough; it stops being enough the moment a second sheet exists.
    """
    st = load_structure(ubq_pdb)
    sheets = [r for r in st.secondary_structure if r.kind == "SHEET"]
    assert {r.group_id for r in sheets} == {"BET"}
    assert [r.serial for r in sheets] == ["1", "2", "3", "4", "5"]


def test_an_unreadable_endpoint_is_kept_and_refuses_that_chains_sse():
    """A hybrid-36 or otherwise unreadable endpoint is never dropped silently:
    the record is kept, and SSE annotation of its chain is refused (ADR-110)."""
    lines = ["HELIX    1  H1 ILE A A000  GLU A A002  1\n"]
    st = parse_pdb_legacy(lines, StructureId("X"))
    (record,) = st.secondary_structure
    assert (record.start_seq, record.start_text) == (None, "A000")
    assert any("cannot read" in w for w in st.warnings)
    with pytest.raises(InvalidRequest, match="cannot read"):
        build_assignment_map(_structure(st.secondary_structure), "A")


# ---------------------------------------------------------------------------
# The conventions
# ---------------------------------------------------------------------------


def test_a_records_range_is_inclusive_of_both_endpoints():
    m = build_assignment_map(_structure([_record(HELIX, "1", 10, 12)]), "A")
    assert [m.assignment_for(_identity(s)).state for s in (9, 10, 11, 12, 13)] == [
        NONE, HELIX, HELIX, HELIX, NONE,
    ]


def test_records_on_other_chains_are_ignored():
    m = build_assignment_map(_structure([_record(HELIX, "1", 10, 12, chain="B")]), "A")
    assert m.assignment_for(_identity(11)).state == NONE
    assert not m.has_records


def test_helix_beats_sheet_when_a_residue_is_in_both():
    st = _structure([_record(SHEET, "4", 10, 20), _record(HELIX, "1", 15, 25)])
    m = build_assignment_map(st, "A")
    assert m.assignment_for(_identity(17)).token == "HELIX_1"
    assert m.assignment_for(_identity(12)).token == "SHEET_X/4"
    assert any("both a HELIX and a SHEET" in w for w in m.warnings)


def test_the_first_record_wins_when_two_of_a_kind_overlap():
    st = _structure([_record(HELIX, "1", 10, 20), _record(HELIX, "2", 15, 25)])
    m = build_assignment_map(st, "A")
    assert m.assignment_for(_identity(17)).group == SSEElementId(serial="1")
    assert m.assignment_for(_identity(22)).group == SSEElementId(serial="2")


_ICODES = [(26, ""), (27, ""), (27, "A"), (27, "B"), (28, ""),
           (184, "A"), (184, ""), (185, "")]


def test_a_range_covers_insertion_coded_residues_between_its_endpoints():
    """Endpoints are matched with their insertion codes and the range is the
    chain's own order between them (ADR-110): ``26..28`` covers 27A and 27B."""
    m = build_assignment_map(_structure([_record(HELIX, "1", 26, 28)], _ICODES), "A")
    assert [m.assignment_for(_identity(s, i)).state for s, i in _ICODES[:5]] == [HELIX] * 5
    assert m.assignment_for(_identity(184, "A")).state == NONE


def test_an_insertion_coded_endpoint_covers_only_its_own_residues():
    """``27A..27B`` is those two residues, never residue 27 (the numeric
    reading covered 27 and missed both)."""
    m = build_assignment_map(_structure(
        [_record(HELIX, "1", 27, 27, start_icode="A", end_icode="B")], _ICODES), "A")
    assert [m.assignment_for(_identity(s, i)).state for s, i in _ICODES[:5]] == [
        NONE, NONE, HELIX, HELIX, NONE]


def test_the_chain_order_decides_when_an_insertion_code_precedes_its_number():
    """3PTB writes 184A before 184: ``184A..185`` covers 184 too."""
    m = build_assignment_map(_structure(
        [_record(SHEET, "1", 184, 185, start_icode="A")], _ICODES), "A")
    assert [m.assignment_for(_identity(s, i)).state for s, i in _ICODES[5:]] == [SHEET] * 3


@pytest.mark.parametrize("record, reason", [
    (_record(HELIX, "1", 10, 12, chain="A", end_chain="B"), "ends on another chain"),
    (_record(HELIX, "1", 12, 10), "ends before it starts"),
    (_record(HELIX, "1", 10, 99), "end residue 99, which is not in the chain"),
    (_record(HELIX, "1", 27, 28, start_icode="C"), "start residue 27C, which is not in the chain"),
])
def test_an_endpoint_that_cannot_be_resolved_refuses_the_chains_sse(record, reason):
    """No numeric fallback and no silent coil (ADR-110)."""
    residues = _RESIDUES if record.end_seq != 28 else _ICODES
    with pytest.raises(InvalidRequest, match=reason):
        build_assignment_map(_structure([record], residues), "A")


def test_an_endpoint_written_in_two_separate_places_refuses():
    residues = [(10, ""), (11, ""), (12, ""), (10, "")]
    with pytest.raises(InvalidRequest, match="two separate places"):
        build_assignment_map(_structure([_record(HELIX, "1", 10, 12)], residues), "A")


def test_a_problem_on_another_chain_does_not_stop_this_one():
    st = _structure([_record(HELIX, "1", 10, 12),
                     _record(HELIX, "2", 500, 501, chain="B")])
    assert build_assignment_map(st, "A").assignment_for(_identity(11)).state == HELIX


def test_grouping_by_sheet_id_alone_is_refused():
    """The third option that is not an option.

    A β-sheet contains one or more strands, so grouping by ``sheetID`` would
    merge strands the categories exist to tell apart -- it would turn 1UBQ's five
    strands into one group and every 1UBQ β-triangle into ``7_3b1``. It is
    refused by name, and the message says why (ADR-063).
    """
    with pytest.raises(ValueError, match="unknown strand identity convention"):
        build_assignment_map(_structure([]), "A", strand_identity="by_sheet_id")


def test_the_standard_is_the_compound_strand_identity():
    m = build_assignment_map(_structure([]), "A")
    assert m.strand_identity == STRAND_IDENTITY_SHEET_ID_AND_SERIAL
    assert m.to_dict()["strand_identity"] == "sheet_id_and_serial"


def test_the_serial_only_strand_identity_is_refused():
    """The historical serial-only identity was removed; only one convention exists."""
    assert STRAND_IDENTITIES == (STRAND_IDENTITY_SHEET_ID_AND_SERIAL,)
    with pytest.raises(ValueError, match="unknown strand identity convention"):
        build_assignment_map(_structure([]), "A", strand_identity="legacy_serial_only")


def test_a_strand_carries_its_sheet_and_a_helix_does_not():
    """A strand is (sheetID, serNum); a helix is named by its serial alone."""
    strand = build_assignment_map(
        _structure([_record(SHEET, "1", 10, 12, group_id="A")]), "A"
    )
    assert strand.assignment_for(_identity(11)).group == SSEElementId("1", "A")
    assert strand.assignment_for(_identity(11)).token == "SHEET_A/1"
    helix = build_assignment_map(
        _structure([_record(HELIX, "3", 10, 12, group_id="HB")]), "A"
    )
    assert helix.assignment_for(_identity(11)).group == SSEElementId(serial="3")
    assert helix.assignment_for(_identity(11)).token == "HELIX_3"


def test_two_sheets_numbering_their_strands_from_one_stay_apart():
    """The collision the historical identity cannot see."""
    records = [
        _record(SHEET, "1", 10, 12, group_id="A"),
        _record(SHEET, "1", 30, 32, group_id="B"),
    ]
    std = build_assignment_map(_structure(records), "A")
    assert (
        std.assignment_for(_identity(11)).group
        != std.assignment_for(_identity(31)).group
    )
    assert std.assignment_for(_identity(11)).token == "SHEET_A/1"
    assert std.assignment_for(_identity(31)).token == "SHEET_B/1"


# ---------------------------------------------------------------------------
# The missing-annotation distinction
# ---------------------------------------------------------------------------


def test_a_structure_with_no_records_warns_rather_than_reporting_all_coil_silently():
    m = build_assignment_map(_structure([]), "A")
    assert not m.has_records
    assert len(m) == 0
    assert any("no HELIX or SHEET record" in w for w in m.warnings)


def test_that_warning_reaches_the_result(bdna_path):
    """1BNA is B-DNA: no HELIX, no SHEET, so every triangle would be ``17_3c``."""
    import biotsr

    st = load_structure(bdna_path)
    assert st.secondary_structure == ()


def test_coverage_counts_the_residues_actually_in_play(ubq_pdb):
    import biotsr

    result = biotsr.analyze(ubq_pdb, method="protein", chain="A", sse=True)
    coverage = result.annotations["sse_coverage"]
    assert sum(coverage.values()) == result.selection.n_points == 76
    assert coverage == {"n_helix": 16, "n_sheet": 33, "n_coil": 27}


def test_the_assignment_map_matches_the_pdb_records_on_1ubq(ubq_pdb):
    """Rebuilt from 1UBQ's HELIX/SHEET records: inclusive ranges, helix first,
    the first record winning an overlap, and every strand named by its sheet."""
    st = load_structure(ubq_pdb)
    m = build_assignment_map(st, "A")
    expected = {}
    for seq in range(23, 35):
        expected[seq] = "HELIX_1"
    for seq in range(56, 60):
        expected[seq] = "HELIX_2"
    for serial, (lo, hi) in zip(
        ("1", "2", "3", "4", "5"),
        ((10, 17), (1, 7), (64, 72), (40, 45), (48, 50)),
    ):
        for seq in range(lo, hi + 1):
            expected.setdefault(seq, "SHEET_BET/%s" % serial)
    assert {k[1]: v.token for k, v in m.by_residue.items()} == expected


# ---------------------------------------------------------------------------
# The annotator
# ---------------------------------------------------------------------------


def test_the_annotator_reads_identities_and_nothing_else(ubq_pdb):
    st = load_structure(ubq_pdb)
    annotator = SSEAnnotator(assignments=build_assignment_map(st, "A"))
    triple, category = annotator.annotate(
        [_identity(1), _identity(2), _identity(3)]
    )
    assert [a.token for a in triple] == ["SHEET_BET/2"] * 3
    assert category == "7_3b1"
