"""``nucleotide_chain`` takes every nucleotide with the atoms it has (ADR-114).

The chain-only rule for terminal residues (ADR-111) is applied without a
terminal-residue detector: no atom count is required of a nucleotide in a chain
calculation, wherever it sits, and every nucleotide whose atoms differ from its
standard set is named in a warning.  That warning describes the structure file;
the calculation used every atom present, so the result is not partial.  The
residue-specific methods (``nucleotide``, ``nucleotide_atom_selection``) keep
the exact atom count.

1BNA is the unmodified RCSB entry.  The other structures are 1BNA with a few
records removed, renamed or added, written by the tests below.
"""

from __future__ import annotations

import os

import pytest

import biotsr
from biotsr.domain.errors import InvalidRequest
from biotsr.domain.selection_record import selection_omissions
from biotsr.selection.atomic import ATOM_SET_NOTE
from biotsr.selection.completeness import (
    AVAILABLE_NUCLEOTIDE_ATOMS,
    LEGACY_NUCLEOTIDE_ATOM_COUNTS,
    STANDARD_NUCLEOTIDE_ATOM_NAMES,
)

from ..conftest import STRUCTURES

BNA = os.path.join(STRUCTURES, "1BNA.pdb")
NOTE = ATOM_SET_NOTE.split(":")[0]


def _lines():
    with open(BNA) as fh:
        return fh.readlines()


def _is(line, chain, seq, atom=None):
    return (line.startswith("ATOM  ") and line[21] == chain
            and line[22:26].strip() == str(seq)
            and (atom is None or line[12:16].strip() == atom))


def _write(tmp_path, name, lines):
    path = tmp_path / name
    path.write_text("".join(lines))
    return str(path)


def _excerpt(tmp_path, residues=(1, 2, 3), lines=None, name="excerpt.pdb"):
    """Chain A residues ``residues`` of 1BNA (or of ``lines``), with their own TER."""
    out = [l for l in (lines or _lines()) if l.startswith("ATOM  ") and l[21] == "A"
           and l[22:26].strip() in {str(r) for r in residues}]
    out += ["%-6s%5d      %3s %1s%4d\n" % ("TER", len(out) + 1, out[-1][17:20], "A",
                                          int(out[-1][22:26])), "END\n"]
    return _write(tmp_path, name, out)


def _points_by_residue(result):
    out = {}
    for t in result.triplets:
        for i in t.identities:
            out.setdefault("%s_%s" % (i.residue_name, i.residue.legacy_seq_label),
                           set()).add(i.atom_name)
    return out


# ===========================================================================
# The standard atom sets used for the report
# ===========================================================================


def test_the_standard_sets_have_exactly_the_historical_counts():
    for name, atoms in STANDARD_NUCLEOTIDE_ATOM_NAMES.items():
        assert len(atoms) == len(set(atoms)) == LEGACY_NUCLEOTIDE_ATOM_COUNTS[name]
        assert "OP3" not in atoms and atoms[:3] == ("P", "OP1", "OP2")


def test_the_policy_admits_every_listed_nucleotide_whatever_its_count():
    for count in (1, 16, 19, 25):
        assert AVAILABLE_NUCLEOTIDE_ATOMS.admit("DC", count).admitted
    refused = AVAILABLE_NUCLEOTIDE_ATOMS.admit("PSU", 20)
    assert not refused.admitted and refused.reason == "unlisted_entity_type"


# ===========================================================================
# 1BNA chain A: the 5' end
# ===========================================================================


@pytest.fixture(scope="module")
def chain_a():
    return biotsr.analyze(BNA, method="nucleotide_chain", chain="A", keep_triplets=False)


def test_1bna_chain_a_gains_exactly_dc_1s_sixteen_atoms(chain_a):
    """227 points before (DC 1 dropped by its count, 16 of 19); 243 now."""
    assert chain_a.counts.n_points == 243
    assert chain_a.counts.n_triangles == 243 * 242 * 241 // 6 == 2362041
    assert list(chain_a.selection.residues_used) == [
        "1_DC", "2_DG", "3_DC", "4_DG", "5_DA", "6_DA", "7_DT", "8_DT",
        "9_DC", "10_DG", "11_DC", "12_DG"]
    dc1 = [l for l in _lines() if _is(l, "A", 1)]
    assert len(dc1) == 16
    assert 243 - 227 == len(dc1)


def test_dc_1_is_named_with_what_it_lacks_and_nothing_else_is(chain_a):
    notes = [w for w in chain_a.warnings if w.startswith(NOTE)]
    assert notes == [ATOM_SET_NOTE % "DC_1 lacks P, OP1, OP2"]
    assert not any("completeness rule" in w for w in chain_a.warnings)


def test_the_note_does_not_make_the_result_partial(chain_a):
    sides = chain_a.annotations["selection_request"]["sides"]
    assert selection_omissions(sides, chain_a.warnings) == {}
    assert sides[0].get("excluded") == [] and "declined" not in sides[0]


def test_both_strands_are_treated_alike():
    chain_b = biotsr.analyze(BNA, method="nucleotide_chain", chain="B", keep_triplets=False)
    assert chain_b.counts.n_points == 243
    assert chain_b.warnings[-1] == ATOM_SET_NOTE % "DC_13 lacks P, OP1, OP2"


# ===========================================================================
# Inside the chain
# ===========================================================================


def test_an_incomplete_internal_nucleotide_takes_part_with_the_atoms_it_has(tmp_path):
    """DG 4 without C8 and N7, between DC 3 and DA 5: its other 20 atoms are
    points, and it is named."""
    lines = [l for l in _lines() if not (_is(l, "A", 4, "C8") or _is(l, "A", 4, "N7"))]
    path = _excerpt(tmp_path, (3, 4, 5), lines, "gap.pdb")
    result = biotsr.analyze(path, method="nucleotide_chain", chain="A")
    assert result.counts.n_points == 19 + 20 + 21
    by_residue = _points_by_residue(result)
    assert by_residue["DG_4"] == set(STANDARD_NUCLEOTIDE_ATOM_NAMES["DG"]) - {"C8", "N7"}
    assert by_residue["DC_3"] == set(STANDARD_NUCLEOTIDE_ATOM_NAMES["DC"])
    assert by_residue["DA_5"] == set(STANDARD_NUCLEOTIDE_ATOM_NAMES["DA"])
    assert result.warnings == (ATOM_SET_NOTE % "DG_4 lacks C8, N7",)
    sides = result.annotations["selection_request"]["sides"]
    assert selection_omissions(sides, result.warnings) == {}


def test_a_lone_phosphate_is_still_the_nucleotide_it_belongs_to(tmp_path):
    """Only P, OP1 and OP2 of DA 5 present -- as in some real chain breaks."""
    keep = {"P", "OP1", "OP2"}
    lines = [l for l in _lines() if not (_is(l, "A", 5) and l[12:16].strip() not in keep)]
    path = _excerpt(tmp_path, (4, 5, 6), lines, "lone.pdb")
    result = biotsr.analyze(path, method="nucleotide_chain", chain="A", keep_triplets=False)
    assert result.counts.n_points == 22 + 3 + 21
    assert list(result.selection.residues_used) == ["4_DG", "5_DA", "6_DA"]
    assert result.warnings[-1].startswith(ATOM_SET_NOTE % "DA_5 lacks O5', C5', C4'")


def test_an_extra_or_differently_named_atom_is_used_and_named(tmp_path):
    """OP1 of DG 4 written O1P (the pre-2007 spelling), and an OP3 added to DG 2:
    every record is a point; the note says what differs."""
    lines = []
    for line in _lines():
        if _is(line, "A", 4, "OP1"):
            line = line[:12] + " O1P" + line[16:]
        lines.append(line)
        if _is(line, "A", 2, "OP1"):
            x = float(line[30:38]) + 1.5
            lines.append(line[:12] + " OP3" + line[16:30] + "%8.3f" % x + line[38:])
    path = _excerpt(tmp_path, (1, 2, 3, 4), lines, "names.pdb")
    result = biotsr.analyze(path, method="nucleotide_chain", chain="A", keep_triplets=False)
    assert result.counts.n_points == 16 + 23 + 19 + 22
    assert result.warnings[-1] == ATOM_SET_NOTE % (
        "DC_1 lacks P, OP1, OP2; DG_2 has OP3 in addition; DG_4 lacks OP1; "
        "has O1P in addition")


# ===========================================================================
# What does not change
# ===========================================================================


def test_the_residue_specific_methods_keep_the_exact_count():
    one = biotsr.analyze(BNA, method="nucleotide", chain="A", entity_name="DC",
                         entity_id="1")
    assert one.counts.n_points == 0
    assert "expected exactly 19 atoms, got 16" in one.warnings[-1]
    narrowed = biotsr.analyze(BNA, method="nucleotide_atom_selection", chain="A",
                              entity_name="DC", entity_id="1", atoms="C1';N1;C2")
    assert narrowed.counts.n_points == 0
    two = biotsr.analyze(BNA, method="nucleotide", chain="A", entity_name="DG",
                         entity_id="2")
    assert two.counts.n_points == 22


def test_the_carbohydrate_rules_are_untouched():
    assert biotsr.get_method("carbohydrate").completeness_policy == \
        "exact_legacy_carbohydrate_atom_count"
    assert biotsr.get_method("carbohydrate_chain").completeness_policy is None
    assert biotsr.get_method("nucleotide_protein_cross").completeness_policy is None


def test_a_nucleotide_outside_the_vocabulary_is_still_not_read(tmp_path):
    """Naming a non-standard nucleotide in the vocabulary does not admit it."""
    lines = [l[:17] + "5CM" + l[20:] if _is(l, "A", 3) else l for l in _lines()]
    path = _excerpt(tmp_path, (1, 2, 3, 4), lines, "mod.pdb")
    result = biotsr.analyze(path, method="nucleotide_chain", chain="A",
                            entity_names="DC;DG;5CM", keep_triplets=False)
    assert "3_5CM" not in result.selection.residues_used
    assert any("5CM is not one of the 8 residue types" in w for w in result.warnings)


def test_alternate_locations_models_and_duplicates_follow_the_input_rules(tmp_path):
    lines = []
    for line in _lines():
        if _is(line, "A", 2, "C8"):
            lines.append(line[:16] + "A" + line[17:])
            x = float(line[30:38]) + 0.4
            lines.append(line[:16] + "B" + line[17:30] + "%8.3f" % x + line[38:])
            continue
        lines.append(line)
    lines = [l for l in lines if l.startswith("ATOM  ") and l[21] == "A"
             and l[22:26].strip() in ("1", "2", "3")]
    model2 = [l for l in _lines() if _is(l, "A", 4)]
    lines = (["MODEL        1\n"] + lines + ["TER\n", "ENDMDL\n", "MODEL        2\n"]
             + model2 + ["ENDMDL\n", "END\n"])
    result = biotsr.analyze(_write(tmp_path, "alt.pdb", lines), method="nucleotide_chain",
                            chain="A", keep_triplets=False)
    assert result.counts.n_points == 57
    assert list(result.selection.residues_used) == ["1_DC", "2_DG", "3_DC"]
    assert any(w.startswith("alternate locations: kept blank and A; excluded 1")
               for w in result.warnings)
    duplicate = [l for l in _lines()]
    at = next(i for i, l in enumerate(duplicate) if _is(l, "A", 2, "N3"))
    duplicate.insert(at, duplicate[at])
    with pytest.raises(InvalidRequest, match="N3"):
        biotsr.analyze(_excerpt(tmp_path, (1, 2, 3), duplicate, "dup.pdb"),
                       method="nucleotide_chain", chain="A")


def test_the_processing_identity_records_the_chain_rule(tmp_path):
    rules = biotsr.api.processing_identity("nucleotide_chain")
    assert rules["completeness_policy"] == "available_nucleotide_atoms"
    assert biotsr.api.processing_identity("nucleotide")["completeness_policy"] == \
        "exact_legacy_nucleotide_atom_count"
    for calpha in ("protein", "aa_selection", "protein_protein_calpha_cross"):
        assert "completeness_policy" not in biotsr.api.processing_identity(calpha)
    result = biotsr.analyze(_excerpt(tmp_path), method="nucleotide_chain", chain="A",
                            keep_triplets=False)
    assert result.sampling.completeness_policy == "available_nucleotide_atoms"


def test_the_three_nucleotide_excerpt_has_all_57_atoms(tmp_path):
    result = biotsr.analyze(_excerpt(tmp_path), method="nucleotide_chain", chain="A",
                            keep_triplets=False)
    assert result.counts.n_points == 16 + 22 + 19
    assert result.counts.n_triangles == 57 * 56 * 55 // 6
