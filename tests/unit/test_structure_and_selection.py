"""PDB parsing, legacy truncation, and the protein C-alpha selector."""

from __future__ import annotations

import pytest

from biotsr.domain.errors import InvalidRequest
from biotsr.domain.identity import ChainId, StructureId
from biotsr.labels.schemes import AMINO_ACID_30, UnlabelledPointPolicy
from biotsr.selection.protein import select_calpha
from biotsr.selection.truncation import legacy_truncation_line
from biotsr.structure.backends.pdb_legacy import parse_pdb_legacy, read_pdb_legacy

SID = StructureId("TEST")


def atom_line(
    serial=1,
    name=" CA ",
    alt=" ",
    res="LEU",
    chain="A",
    seq=15,
    icode=" ",
    x=None,
    y=2.0,
    z=3.0,
    record="ATOM  ",
    element=" C",
):
    # Each residue and record its own position: two selected atoms at one
    # position are refused (ADR-108).  Residue 15, record 1 sits at x = 1.0.
    if x is None:
        x = 1.0 + 3.8 * (seq - 15) + 0.1 * (serial - 1) + (0.5 if icode.strip() else 0.0)
    return (
        "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f  1.00 20.00          %2s\n"
        % (record, serial, name, alt, res, chain, seq, icode, x, y, z, element)
    )


def parse(lines):
    return parse_pdb_legacy(lines, SID)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def test_fields_are_parsed_into_full_identity():
    s = parse([atom_line(seq=15, icode="B", alt="A", chain="C", res="MET", x=1.0)])
    assert len(s.atoms) == 1
    a = s.atoms[0].identity
    assert a.structure == SID
    assert a.chain == ChainId("C")
    assert a.residue_name == "MET"
    assert a.residue_seq == 15
    assert a.insertion_code == "B"
    assert a.atom_name == "CA"
    assert a.alt_loc == "A"
    assert s.atoms[0].coord.as_tuple() == (1.0, 2.0, 3.0)


def test_hetatm_is_flagged_not_dropped():
    s = parse([atom_line(record="HETATM", res="HOH", name=" O  ")])
    assert len(s.atoms) == 1
    assert s.atoms[0].hetero is True


def test_blank_insertion_code_and_altloc_become_empty_strings():
    a = parse([atom_line()]).atoms[0].identity
    assert a.insertion_code == ""
    assert a.alt_loc == ""


def test_file_order_is_recorded():
    s = parse([atom_line(seq=1), "REMARK dummy\n", atom_line(seq=2)])
    assert [a.file_order for a in s.atoms] == [0, 2]


def test_model_number_is_tracked_onto_atoms():
    s = parse(
        ["MODEL        1\n", atom_line(seq=1), "ENDMDL\n", "MODEL        2\n", atom_line(seq=2)]
    )
    assert [a.identity.model for a in s.atoms] == [1, 2]


def test_truncation_events_are_recorded_with_their_context():
    s = parse(
        [
            atom_line(seq=1),
            "TER     123      LEU A  15\n",
            "MODEL        2\n",
            "ENDMDL\n",
        ]
    )
    kinds = [(e.kind, e.file_order) for e in s.truncation_events]
    assert kinds == [("TER", 1), ("MODEL", 2), ("ENDMDL", 3)]
    ter = s.truncation_events[0]
    assert ter.chain == ChainId("A")
    assert s.truncation_events[1].model_number == 2


def test_parsing_is_chain_agnostic():
    """Truncation is a selection decision, so the parser keeps everything."""
    s = parse([atom_line(chain="A"), atom_line(chain="B"), atom_line(chain="C")])
    assert len(s.atoms) == 3


def test_short_ter_line_warns_instead_of_raising():
    """The legacy reader indexes line[21] unconditionally and would raise."""
    s = parse([atom_line(), "TER\n"])
    assert any("TER record shorter" in w for w in s.warnings)
    assert s.truncation_events[0].chain == ChainId("")


def test_unparseable_coordinates_are_skipped_with_a_warning():
    bad = "ATOM      1  CA  LEU A  15     ABCDEFGH   2.000   3.000  1.00 20.00           C\n"
    s = parse([bad])
    assert s.atoms == ()
    assert any("coordinates" in w for w in s.warnings)


def test_the_two_legacy_atom_name_slices_are_checked_not_assumed():
    """line[13:15] (the legacy CA predicate) vs line[13:16] (the name field).

    They agree for every standard protein atom.  A record where they disagree
    produces a loud warning rather than a silent divergence.
    """
    ordinary = parse([atom_line(name=" CA ")])
    assert ordinary.warnings == ()
    # A contrived four-character name whose 2-char slice reads "CA" but whose
    # 3-char slice does not.
    contrived = parse([atom_line(name="1CA2")])
    assert any("slices disagree" in w for w in contrived.warnings)


def test_leading_zero_sequence_number_is_reported():
    line = "ATOM      1  CA  LEU A0015       1.000   2.000   3.000  1.00 20.00           C\n"
    s = parse([line])
    assert any("round-trip" in w for w in s.warnings)


def test_structure_id_from_filename_makes_no_length_assumption(tmp_path):
    path = tmp_path / "PDB_0000001UBQ.pdb"
    path.write_text(atom_line())
    s = read_pdb_legacy(str(path))
    assert s.structure_id.value == "PDB_0000001UBQ"
    assert s.source_sha256 is not None


# ---------------------------------------------------------------------------
# Legacy truncation rule
# ---------------------------------------------------------------------------


def test_ter_truncates_only_for_the_matching_chain():
    s = parse([atom_line(chain="A"), "TER     123      LEU B  15\n", atom_line(chain="A", seq=16)])
    assert legacy_truncation_line(s, ChainId("B")) == 1
    assert legacy_truncation_line(s, ChainId("A")) is None


def test_endmdl_truncates_for_every_chain():
    s = parse([atom_line(), "ENDMDL\n", atom_line(seq=16)])
    assert legacy_truncation_line(s, ChainId("A")) == 1
    assert legacy_truncation_line(s, ChainId("Z")) == 1


def test_model_one_does_not_truncate_but_model_two_does():
    s = parse(["MODEL        1\n", atom_line(), "MODEL        2\n", atom_line(seq=16)])
    assert legacy_truncation_line(s, ChainId("A")) == 2


def test_the_earliest_applicable_event_wins():
    s = parse(
        [
            atom_line(chain="A"),
            "TER     123      LEU A  15\n",
            "ENDMDL\n",
        ]
    )
    assert legacy_truncation_line(s, ChainId("A")) == 1


# ---------------------------------------------------------------------------
# Selector
# ---------------------------------------------------------------------------


def select(lines, chain="A", scheme=AMINO_ACID_30, **kwargs):
    return select_calpha(parse(lines), chain, scheme, **kwargs)


def test_selects_only_calpha_atoms_of_the_requested_chain():
    lines = [
        atom_line(chain="A", seq=1, res="MET", name=" CA "),
        atom_line(chain="A", seq=1, res="MET", name=" CB "),
        atom_line(chain="B", seq=1, res="MET", name=" CA "),
        atom_line(chain="A", seq=2, res="LEU", name=" CA "),
    ]
    r = select(lines)
    assert len(r.points) == 2
    assert [p.identity.residue_seq for p in r.points] == [1, 2]
    assert r.chains_used == ("A",)


def hetero_residue(res, seq, chain="A", alt=" ", ca_element=" C", backbone=True,
                   serial=100):
    """The backbone of one ``HETATM`` residue: N, CA and C (or the CA alone)."""
    names = ((" N  ", " N"), (" CA ", ca_element), (" C  ", " C")) if backbone \
        else ((" CA ", ca_element),)
    return [atom_line(serial=serial + n, record="HETATM", res=res, seq=seq, chain=chain,
                      alt=alt, name=name, element=element, x=float(seq), y=float(n))
            for n, (name, element) in enumerate(names)]


def test_a_supported_modified_residue_written_as_hetatm_is_selected():
    """ADR-107: the Ca selector reads a HETATM Ca of a supported modified
    residue whose residue carries N, a carbon CA and C, labelled from the table."""
    lines = ([atom_line(seq=1, res="LEU")] + hetero_residue("MSE", 2)
             + hetero_residue("SEP", 3) + hetero_residue("CSD", 4)
             + [atom_line(seq=5, res="GLY")])
    r = select(lines)
    assert [p.identity.residue_name for p in r.points] == ["LEU", "MSE", "SEP", "CSD", "GLY"]
    assert [p.label for p in r.points] == [8, 13, 25, 33, 11]
    assert r.residues_used == ("1", "2", "3", "4", "5")


def test_a_supported_modified_residue_written_as_atom_is_selected():
    """An ATOM record of a supported modified residue is an ordinary residue
    record; until 1.0.0rc5 it had no code and stopped the calculation."""
    r = select([atom_line(seq=1, res="MSE"), atom_line(seq=2, res="TPO")])
    assert [p.label for p in r.points] == [13, 24]


@pytest.mark.parametrize("lines, why", [
    (hetero_residue("MSE", 2, backbone=False), "a CA without N and C"),
    (hetero_residue("SEP", 2)[1:], "a CA and C without N"),
    (hetero_residue("SEP", 2)[:2], "an N and CA without C"),
    (hetero_residue("SEP", 2, ca_element="  "), "a CA with a blank element field"),
])
def test_a_supported_hetatm_ca_needs_no_other_atom(lines, why):
    """ADR-107 (corrected): the Ca is read on its own -- carbon or blank element,
    supported name, requested chain; N and C are not required."""
    r = select([atom_line(seq=1, res="LEU")] + lines)
    assert [p.identity.residue_name for p in r.points][1:] == [lines[0][17:20]], why
    assert not r.skipped, why


@pytest.mark.parametrize("lines, why", [
    (hetero_residue("HYP", 2), "a modified residue the labels do not cover"),
    (hetero_residue("LIG", 2), "a ligand with N, CA and C atoms"),
    (hetero_residue("CA", 2, ca_element="CA", backbone=False), "a calcium ion"),
    (hetero_residue("SEP", 2, chain="B"), "a supported residue of another chain"),
])
def test_no_other_hetatm_ca_is_selected(lines, why):
    r = select([atom_line(seq=1, res="LEU")] + lines)
    assert [p.identity.residue_name for p in r.points] == ["LEU"], why
    assert not r.skipped, why


def test_a_supported_hetatm_ca_with_another_element_is_skipped_and_recorded():
    r = select([atom_line(seq=1, res="LEU")] + hetero_residue("SEP", 2, ca_element="CA"))
    assert [p.identity.residue_name for p in r.points] == ["LEU"]
    assert [s.reason for s in r.skipped] == ["calpha_element_not_carbon"]
    assert "'CA'" in r.skipped[0].detail


def test_the_element_rule_is_not_applied_to_atom_records():
    """A remaining asymmetry, pinned so it is visible: ``ATOM`` records were never
    element-checked (no record type is, for the regular residues), so an
    ``ATOM`` SEP whose CA declares calcium is read while the same ``HETATM`` CA
    is not.  ADR-107 records it; no real file seen so far writes it."""
    atom = select([atom_line(seq=1, res="SEP", element="CA")])
    hetero = select(hetero_residue("SEP", 1, ca_element="CA"))
    assert [p.identity.residue_name for p in atom.points] == ["SEP"]
    assert hetero.points.labels == ()


def test_the_altloc_rule_applies_to_modified_residues():
    lines = (hetero_residue("SEP", 2, alt="A") + hetero_residue("SEP", 2, alt="B", serial=200)
             + hetero_residue("TPO", 3, alt="B", serial=300))
    r = select([atom_line(seq=1, res="LEU")] + lines)
    assert [(p.identity.residue_name, p.identity.alt_loc) for p in r.points] == [
        ("LEU", ""), ("SEP", "A")]
    assert [s.reason for s in r.skipped] == ["altloc_rejected", "altloc_rejected"]


def test_a_modified_residue_after_the_chains_ter_is_not_read():
    lines = ([atom_line(seq=1, res="LEU"), "TER       2      LEU A   1\n"]
             + hetero_residue("SEP", 900))
    r = select(lines)
    assert [p.identity.residue_name for p in r.points] == ["LEU"]


def test_altloc_policy_accepts_blank_and_A_and_records_the_rest():
    lines = [
        atom_line(seq=1, alt=" "),
        atom_line(seq=2, alt="A"),
        atom_line(seq=3, alt="B"),
    ]
    r = select(lines)
    assert len(r.points) == 2
    assert [s.reason for s in r.skipped] == ["altloc_rejected"]


def test_a_unk_residue_in_the_chain_refuses_the_selection():
    """No calculation is done for UNK, and none for the rest without it
    (ADR-109)."""
    with pytest.raises(InvalidRequest, match="UNK_15"):
        select([atom_line(res="UNK"), atom_line(seq=2, res="LEU")])


def test_unk_can_still_be_reported_without_refusing():
    r = select([atom_line(res="UNK"), atom_line(seq=2, res="LEU")], refuse_unk=False)
    assert len(r.points) == 1
    assert r.skipped[0].reason == "excluded_residue_name"


def test_points_come_out_in_file_order():
    lines = [atom_line(seq=s) for s in (30, 10, 20)]
    r = select(lines)
    assert [p.identity.residue_seq for p in r.points] == [30, 10, 20]


def test_labels_come_from_the_scheme():
    r = select([atom_line(res="MET"), atom_line(seq=2, res="ILE")])
    assert [p.label for p in r.points] == [13, 10]


def test_unknown_residue_refuses_naming_it():
    """A residue the labels cannot represent refuses the calculation, naming
    the residue and what to do (ADR-108)."""
    with pytest.raises(InvalidRequest, match="residue XYZ 15 of chain A cannot be labelled"):
        select([atom_line(res="XYZ")])


def test_unlabelled_drop_policy_records_instead_of_raising():
    r = select([atom_line(res="XYZ"), atom_line(seq=2, res="LEU")],
               unlabelled_policy=UnlabelledPointPolicy.DROP_AND_RECORD)
    assert len(r.points) == 1
    assert r.skipped[0].reason == "unlabelled_point"


def test_selection_stops_at_the_legacy_truncation_point():
    lines = [
        atom_line(chain="A", seq=1),
        "TER     123      LEU A   1\n",
        atom_line(chain="A", seq=2),
    ]
    r = select(lines, chain="A")
    assert len(r.points) == 1
    assert r.truncated_at_line == 1


def test_unsupported_policies_are_refused_rather_than_ignored():
    with pytest.raises(ValueError, match="altloc"):
        select([atom_line()], altloc_policy="all")
    with pytest.raises(ValueError, match="truncation"):
        select([atom_line()], truncation="none")
