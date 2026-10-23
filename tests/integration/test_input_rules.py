"""The 1.0.0rc6 reading and selection rules, one decision at a time.

ADR-108 (reading atom records), ADR-109 (selections), ADR-110 (SSE endpoints),
ADR-111 (terminal residues) and ADR-112 (required chains).  Small synthetic
structures isolate each rule; the real fixtures show that what the rules do not
touch is unchanged.
"""

from __future__ import annotations

import csv
import os

import pytest

import biotsr
from biotsr import InvalidRequest
from biotsr.cli.main import main

from ..conftest import structure


# ---------------------------------------------------------------------------
# A tiny PDB writer: every column where the PDB format puts it
# ---------------------------------------------------------------------------


def atom(name, res, chain, seq, xyz, element, record="ATOM", alt="", icode="",
         charge="", serial=1):
    """One ``ATOM``/``HETATM`` line; ``name`` is the 4-column field as written."""
    x, y, z = xyz
    return ("%-6s%5d %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f  1.00 20.00          %2s%-2s"
            % (record, serial, name, alt, res, chain, seq, icode, x, y, z,
               element, charge))


def residue(res, chain, seq, x0, icode="", record="ATOM", extra=()):
    """A glycine-shaped residue: N, CA, C, O (plus ``extra`` lines)."""
    return [
        atom(" N  ", res, chain, seq, (x0, 0.0, 0.0), " N", record, icode=icode),
        atom(" CA ", res, chain, seq, (x0 + 1.2, 0.9, 0.0), " C", record, icode=icode),
        atom(" C  ", res, chain, seq, (x0 + 2.4, 0.1, 0.2), " C", record, icode=icode),
        atom(" O  ", res, chain, seq, (x0 + 2.5, -1.1, 0.4), " O", record, icode=icode),
    ] + list(extra)


def chain_of(names, chain="A", start=1, spacing=3.8, offset=(0.0, 0.0)):
    lines = []
    for i, name in enumerate(names):
        lines += residue(name, chain, start + i, i * spacing)
    if offset != (0.0, 0.0):
        lines = [_shift(line, offset) for line in lines]
    return lines


def _shift(line, offset):
    y = float(line[38:46]) + offset[0]
    z = float(line[46:54]) + offset[1]
    return line[:38] + "%8.3f%8.3f" % (y, z) + line[54:]


def ter(chain):
    return "TER   %5d      %3s %1s%4d" % (99, "GLY", chain, 2)


def write(tmp_path, lines, name="SYN.pdb"):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\nEND\n")
    return str(path)


# ---------------------------------------------------------------------------
# A. Atom names: columns 13-16 in atomic paths; the Calpha test unchanged
# ---------------------------------------------------------------------------


def test_an_atom_record_metal_keeps_its_name_and_element_label(tmp_path):
    """An ATOM-record zinc (``ZN`` from column 13) used to be read as ``N`` and
    labelled nitrogen (108); it is now ``ZN``, labelled zinc."""
    zinc = atom("ZN  ", "ZN", "A", 3, (4.0, 6.0, 1.0), "ZN")
    path = write(tmp_path, residue("GLY", "A", 1, 0.0) + residue("ALA", "A", 2, 3.8) + [zinc])
    result = biotsr.analyze(path, method="amino_acid_atom_selection", chain="A",
                            atoms="GLY_1_CA;ALA_2_CA;ZN_3_ZN")
    labels = {t.identities[i].atom_name: t.labels[i]
              for t in result.triplets for i in range(3)}
    assert labels == {"CA": 106, "ZN": 105}


def test_calcium_and_ca1_are_never_calpha_points(tmp_path):
    lines = chain_of(["GLY", "ALA", "SER"]) + [
        atom("CA  ", "CA", "A", 9, (20.0, 5.0, 5.0), "CA", record="HETATM"),
        atom(" CA1", "LIG", "A", 10, (25.0, 5.0, 5.0), " C", record="HETATM"),
    ]
    result = biotsr.analyze(write(tmp_path, lines), method="protein", chain="A")
    assert result.selection.residues_used == ("1", "2", "3")


def test_ordinary_atoms_and_frozen_results_are_unchanged():
    """Standard atoms start in column 14, so the two slices agree; the frozen
    references (test_standard_references.py) prove it at scale."""
    ubq = biotsr.load_structure(structure("1UBQ"))
    protein_atoms = [a for a in ubq.atoms if not a.hetero]
    assert all(a.atom_name == a.calpha_name for a in protein_atoms)


# ---------------------------------------------------------------------------
# B. Narrowed Calpha cross: both selections narrow the interface
# ---------------------------------------------------------------------------


def _heavy(st, chain, seqs=None):
    out = []
    for a in st.atoms:
        if a.chain.value != chain or a.hetero or a.element.strip().upper() in ("H", "D"):
            continue
        if a.identity.alt_loc not in ("", "A"):
            continue
        if seqs is not None and a.identity.residue.legacy_seq_label not in seqs:
            continue
        out.append(a)
    return out


def test_a_narrowed_side_narrows_the_partner_interface():
    """1BRS, chain A narrowed to residues 27, 83, 87 and 102 against chain D.
    Every partner residue in the pool touches (heavy atoms within 5 A) one of
    the selected residues; none is kept for touching only the rest of A."""
    path = structure("1BRS")
    st = biotsr.load_structure(path)
    selected = ("27", "83", "87", "102")
    result = biotsr.analyze_cross(path, method="protein_protein_calpha_cross",
                                  chain="A", chain_2="D",
                                  calpha_residues=";".join(selected))
    a_atoms = _heavy(st, "A", selected)
    expected_d = set()
    for b in _heavy(st, "D"):
        if any(sum((p - q) ** 2 for p, q in zip(a.coord.as_tuple(), b.coord.as_tuple()))
               <= 25.0 for a in a_atoms):
            expected_d.add(b.identity.residue.legacy_seq_label)
    pool_d = {i.residue.legacy_seq_label for t in result.triplets
              for i in t.identities if i.chain.value == "D"}
    assert pool_d and pool_d <= expected_d
    whole = biotsr.analyze_cross(path, method="protein_protein_calpha_cross",
                                 chain="A", chain_2="D")
    whole_d = {i.residue.legacy_seq_label for t in whole.triplets
               for i in t.identities if i.chain.value == "D"}
    assert pool_d < whole_d
    # Measured: 798 triangles under 1.0.0rc5, which kept 7 partner residues
    # that touch only unselected residues of A.
    assert result.counts.n_triangles == 336


def test_a_single_narrowed_residue_matches_the_independent_prediction():
    """A27 against D: 3 partner residues and 3 triangles, as an independent
    contact calculation predicted before the change."""
    result = biotsr.analyze_cross(structure("1BRS"), method="protein_protein_calpha_cross",
                                  chain="A", chain_2="D", calpha_residues="27")
    assert (result.counts.n_points, result.counts.n_triangles) == (4, 3)


def test_a_whole_chain_calpha_cross_is_unchanged():
    """The whole-chain interface is decided as before (counts measured with
    1.0.0rc5 on 1BRS A/D)."""
    result = biotsr.analyze_cross(structure("1BRS"), method="protein_protein_calpha_cross",
                                  chain="A", chain_2="D")
    assert (result.counts.n_points, result.counts.n_triangles) == (41, 8151)


# ---------------------------------------------------------------------------
# C. Terminal residues
# ---------------------------------------------------------------------------


def test_a_terminal_residue_with_oxt_is_excluded_and_named_by_residue_methods():
    path = structure("1UBQ")
    result = biotsr.analyze(path, method="aa_ind", chain="A", residue_types="GLY")
    assert result.counts.n_points == 20
    assert any("76_GLY: GLY: expected exactly 4 atoms, got 5" in w for w in result.warnings)
    single = biotsr.analyze(path, method="amino_acid_atomic", chain="A",
                            residue="76", residue_name="GLY")
    assert single.counts.n_points == 0


def test_cross_methods_keep_terminal_atoms_including_oxt(tmp_path):
    """No completeness rule on a cross side: a C-terminal OXT in contact is a
    vertex like any other atom."""
    oxt = atom(" OXT", "ALA", "A", 2, (6.3, -1.0, -0.8), " O")
    lines = residue("GLY", "A", 1, 0.0) + residue("ALA", "A", 2, 3.8, extra=[oxt]) \
        + [ter("A")] + [_shift(l, (0.0, 3.0)) for l in residue("GLY", "B", 1, 3.0)]
    result = biotsr.analyze_cross(write(tmp_path, lines), method="protein_protein_cross",
                                  chain="A", chain_2="B")
    names = {i.atom_name for t in result.triplets for i in t.identities}
    assert "OXT" in names


def test_calpha_methods_keep_terminal_points():
    result = biotsr.analyze(structure("1UBQ"), method="protein", chain="A")
    assert result.selection.residues_used[-1] == "76"


# ---------------------------------------------------------------------------
# D. Cross sides: distinct protein chains, no shared atoms
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["protein_protein_cross", "protein_protein_calpha_cross"])
def test_protein_protein_cross_needs_two_different_chains(method):
    with pytest.raises(InvalidRequest, match="name the same chain"):
        biotsr.analyze_cross(structure("1BRS"), method=method, chain="A", chain_2="A")


def test_a_molecule_side_may_share_the_protein_chain_identifier(tmp_path):
    nucleotide = [atom(name, "DC", "A", 50, xyz, el) for name, xyz, el in (
        (" P  ", (1.0, 2.5, 0.5), " P"), (" OP1", (2.0, 3.0, 0.8), " O"),
        (" C5'", (3.0, 3.4, 1.0), " C"))]
    lines = chain_of(["GLY", "ALA"]) + nucleotide
    result = biotsr.analyze_cross(write(tmp_path, lines), method="nucleotide_protein_cross",
                                  chain="A", chain_2="A", entity_name="DC")
    assert result.counts.n_triangles > 0


def test_two_sides_sharing_an_atom_are_refused_and_named(tmp_path):
    """A drug side naming TPO selects the very atoms the protein side reads
    (TPO is a protein residue, ADR-101): refused, nothing subtracted."""
    tpo = residue("TPO", "A", 3, 7.6, record="HETATM")
    lines = chain_of(["GLY", "ALA"]) + tpo
    with pytest.raises(InvalidRequest, match="select 4 of the same atom"):
        biotsr.analyze_cross(write(tmp_path, lines), method="drug_protein_cross",
                             chain="A", chain_2="A", entity_name="TPO")


# ---------------------------------------------------------------------------
# E. SSE: endpoints with insertion codes, in chain order
# ---------------------------------------------------------------------------


def _sse_file(tmp_path, helix):
    lines = [helix]
    x = 0.0
    for seq, icode in ((13, ""), (14, ""), (14, "A"), (14, "B"), (15, ""), (16, "")):
        lines += residue("ALA", "A", seq, x, icode=icode)
        x += 3.8
    return write(tmp_path, lines)


def test_insertion_coded_residues_inside_a_helix_are_helix(tmp_path):
    helix = "HELIX    1   1 ALA A   14A ALA A   14B 5                                   2"
    path = _sse_file(tmp_path, helix)
    plain = biotsr.analyze(path, method="protein", chain="A")
    annotated = biotsr.analyze(path, method="protein", chain="A", sse=True)
    assert annotated.annotations["sse_coverage"] == {"n_helix": 2, "n_sheet": 0, "n_coil": 4}
    # An annotation: the keys and both fingerprints are unchanged.
    assert annotated.key_space.fingerprint == plain.key_space.fingerprint
    assert annotated.counts.n_triangles == plain.counts.n_triangles


def test_an_unresolvable_endpoint_refuses_sse_for_that_chain(tmp_path):
    helix = "HELIX    1   1 ALA A   14C ALA A   15  5                                   2"
    path = _sse_file(tmp_path, helix)
    with pytest.raises(InvalidRequest, match="start residue 14C, which is not in the chain"):
        biotsr.analyze(path, method="protein", chain="A", sse=True)
    # The calculation without SSE is not affected.
    assert biotsr.analyze(path, method="protein", chain="A").counts.n_points == 6


# ---------------------------------------------------------------------------
# F. Element fields
# ---------------------------------------------------------------------------


def _ala_with_hydrogens(charge_text, element=" H"):
    extra = [atom(" CB ", "ALA", "A", 1, (1.0, 2.3, 0.4), " C"),
             atom(" H  ", "ALA", "A", 1, (-0.6, 0.8, 0.0), element, charge=charge_text),
             atom(" HA ", "ALA", "A", 1, (1.5, 1.2, -1.0), element, charge=charge_text)]
    return residue("ALA", "A", 1, 0.0, extra=extra)


@pytest.mark.parametrize("charge", ["", "0", "1+"])
def test_amino_acid_atomic_reads_hydrogen_from_columns_77_78(tmp_path, charge):
    """Charge text in columns 79-80 no longer hides a hydrogen (it used to make
    ALA count 7 atoms against 5 and exclude the residue)."""
    path = write(tmp_path, _ala_with_hydrogens(charge))
    result = biotsr.analyze(path, method="amino_acid_atomic", chain="A",
                            residue="1", residue_name="ALA")
    assert result.counts.n_points == 5


def test_a_lower_case_element_is_read_without_regard_to_case(tmp_path):
    path = write(tmp_path, _ala_with_hydrogens("", element=" h"))
    result = biotsr.analyze(path, method="amino_acid_atomic", chain="A",
                            residue="1", residue_name="ALA")
    assert result.counts.n_points == 5


def test_a_selected_atom_with_a_blank_element_is_refused(tmp_path):
    path = write(tmp_path, _ala_with_hydrogens("", element="  "))
    with pytest.raises(InvalidRequest, match="blank element field"):
        biotsr.analyze(path, method="aa_ind", chain="A", residue_types="ALA")


def test_a_blank_element_elsewhere_does_not_stop_the_calculation(tmp_path):
    blank = [_shift(l, (12.0, 0.0)) for l in _ala_with_hydrogens("", element="  ")]
    blank = [l[:21] + "B" + l[22:] for l in blank]
    path = write(tmp_path, chain_of(["GLY", "ALA", "SER"]) + [ter("A")] + blank)
    result = biotsr.analyze(path, method="aa_ind", chain="A", residue_types="GLY")
    assert result.counts.n_points == 4


def test_a_blank_element_calpha_of_a_hetatm_residue_is_still_read(tmp_path):
    """The approved Calpha exception (ADR-107) is kept."""
    mse = [atom(" CA ", "MSE", "A", 3, (8.8, 0.9, 0.0), "  ", record="HETATM")]
    path = write(tmp_path, chain_of(["GLY", "ALA"]) + mse)
    result = biotsr.analyze(path, method="protein", chain="A")
    assert result.selection.residues_used == ("1", "2", "3")


# ---------------------------------------------------------------------------
# G. Repeated, conflicting and coincident atoms
# ---------------------------------------------------------------------------


def test_a_repeated_record_is_refused(tmp_path):
    lines = chain_of(["GLY", "ALA"])
    path = write(tmp_path, lines + [lines[5]])
    with pytest.raises(InvalidRequest, match="repeated records of one atom: line 6 and line 9"):
        biotsr.analyze(path, method="aa_ind_specific", chain="A", residues="ALA_2")


def test_conflicting_records_are_refused_naming_what_differs(tmp_path):
    lines = chain_of(["GLY", "ALA"])
    moved = lines[5][:30] + "%8.3f" % 9.9 + lines[5][38:]
    with pytest.raises(InvalidRequest, match="differ in coordinates"):
        biotsr.analyze(write(tmp_path, lines + [moved]), method="protein", chain="A")


def test_two_atoms_at_one_position_are_refused(tmp_path):
    lines = chain_of(["GLY", "ALA"])
    clash = lines[4][:12] + " CB " + lines[4][16:]   # ALA CB on top of ALA N
    clash = clash[:76] + " C" + clash[78:]
    with pytest.raises(InvalidRequest, match="exactly the same coordinates"):
        biotsr.analyze(write(tmp_path, lines + [clash]), method="aa_ind",
                       chain="A", residue_types="ALA")


# ---------------------------------------------------------------------------
# H. TER: protein sides stop at their chain's TER
# ---------------------------------------------------------------------------


def _glycan_after_ter():
    nag = [atom(name, "NAG", "A", 101, xyz, el, record="HETATM") for name, xyz, el in (
        (" C1 ", (4.0, 2.0, 1.0), " C"), (" O5 ", (5.0, 2.4, 1.2), " O"),
        (" C2 ", (3.0, 2.6, 1.5), " C"))]
    free_ala = residue("ALA", "A", 50, 5.2, record="ATOM")
    free_ala = [_shift(l, (0.5, 1.0)) for l in free_ala]
    return chain_of(["GLY", "ALA"]) + [ter("A")] + nag + free_ala


def test_carbohydrate_cross_protein_side_stops_at_the_chains_ter(tmp_path):
    path = write(tmp_path, _glycan_after_ter())
    result = biotsr.analyze_cross(path, method="carbohydrate_protein_cross",
                                  chain="A", chain_2="A", entity_name="NAG")
    protein = {i.residue_seq for t in result.triplets for i in t.identities
               if i.residue_name != "NAG"}
    assert protein and 50 not in protein


# ---------------------------------------------------------------------------
# I. UNK
# ---------------------------------------------------------------------------


def test_unk_in_a_selected_whole_chain_is_refused(tmp_path):
    path = write(tmp_path, chain_of(["GLY", "UNK", "SER"]))
    with pytest.raises(InvalidRequest, match=r"includes residue\(s\) named UNK \(UNK_2\)"):
        biotsr.analyze(path, method="protein", chain="A")


def test_a_selection_leaving_unk_out_is_computed(tmp_path):
    path = write(tmp_path, chain_of(["GLY", "UNK", "SER", "ALA"]))
    result = biotsr.analyze(path, method="aa_selection", chain="A", residues="1;3;4")
    assert result.counts.n_triangles == 1


def test_unk_on_an_unselected_chain_does_not_matter(tmp_path):
    lines = chain_of(["GLY", "ALA", "SER"]) + [ter("A")] \
        + chain_of(["UNK"], chain="B", offset=(20.0, 0.0))
    assert biotsr.analyze(write(tmp_path, lines), method="protein", chain="A").counts.n_points == 3


@pytest.mark.parametrize("method", ["drug_protein_cross", "carbohydrate_protein_cross"])
def test_unk_on_a_whole_protein_cross_side_is_refused(tmp_path, method):
    entity = "NAG" if method.startswith("carb") else "LIG"
    lig = [atom(name, entity, "B", 1, xyz, el, record="HETATM") for name, xyz, el in (
        (" C1 ", (4.0, 2.0, 1.0), " C"), (" O5 ", (5.0, 2.4, 1.2), " O"),
        (" C2 ", (3.0, 2.6, 1.5), " C"))]
    path = write(tmp_path, chain_of(["GLY", "UNK"]) + [ter("A")] + lig)
    with pytest.raises(InvalidRequest, match="UNK_2"):
        biotsr.analyze_cross(path, method=method, chain="B", chain_2="A",
                             entity_name=entity)


def test_unk_named_by_an_atomic_selection_is_refused():
    with pytest.raises(InvalidRequest, match="UNK"):
        biotsr.analyze(structure("1UBQ"), method="aa_ind", chain="A", residue_types="UNK")


# ---------------------------------------------------------------------------
# J. An atom that cannot be labelled
# ---------------------------------------------------------------------------


def test_an_unlabelable_selected_atom_is_refused_naming_it(tmp_path):
    odd = atom(" QQ ", "LIG", "B", 1, (1.0, 1.0, 1.0), "XX", record="HETATM")
    other = [atom(" C%d " % i, "LIG", "B", 1, (1.0 + i, 2.0, 0.0), " C", record="HETATM")
             for i in (1, 2)]
    path = write(tmp_path, [odd] + other)
    with pytest.raises(InvalidRequest, match=r"atom name 'QQ', element 'XX'\) cannot be labelled"):
        biotsr.analyze(path, method="drug", chain="B", entity_name="LIG")
    # Not selected, not labelled: narrowing to the other atoms is computed.
    result = biotsr.analyze(path, method="drug_atom_selection", chain="B",
                            entity_name="LIG", atoms="C1;C2")
    assert result.counts.n_points == 2


# ---------------------------------------------------------------------------
# K. Confirmed smaller rules
# ---------------------------------------------------------------------------


def test_the_nucleotide_cross_default_cutoff_is_3_5():
    assert biotsr.api.default_contact_cutoff("nucleotide_protein_cross") == 3.5


# ---------------------------------------------------------------------------
# L. Hybrid-36 residue numbers
# ---------------------------------------------------------------------------


def _hybrid(line, text):
    return line[:22] + "%4s" % text + line[26:]


def test_an_unreadable_residue_number_in_the_selection_is_refused(tmp_path):
    lines = chain_of(["GLY", "ALA", "SER"])
    lines = lines[:8] + [_hybrid(l, "A000") for l in lines[8:]]
    path = write(tmp_path, lines)
    st = biotsr.load_structure(path)
    assert len(st.unreadable_numbers) == 4
    with pytest.raises(InvalidRequest, match="residue number BioTSR cannot read"):
        biotsr.analyze(path, method="protein", chain="A")
    with pytest.raises(InvalidRequest, match="cannot read"):
        biotsr.analyze(path, method="aa_ind", chain="A", residue_types="SER")


def test_an_unreadable_number_outside_the_selection_does_not_matter(tmp_path):
    lines = chain_of(["GLY", "ALA", "SER"])
    lines = lines[:8] + [_hybrid(l, "A000") for l in lines[8:]]
    path = write(tmp_path, lines + [ter("A")] + chain_of(["GLY", "ALA", "SER"], chain="B",
                                                         offset=(20.0, 0.0)))
    assert biotsr.analyze(path, method="protein", chain="B").counts.n_points == 3
    # A residue-type selection that names other residues is not affected either.
    assert biotsr.analyze(path, method="aa_ind", chain="A",
                          residue_types="GLY").counts.n_points == 4


# ---------------------------------------------------------------------------
# Processing identity: rc5 outputs are not reused, nor mixed silently
# ---------------------------------------------------------------------------


def test_results_and_processing_identity_record_the_input_rules():
    result = biotsr.analyze(structure("1UBQ"), method="protein", chain="A")
    assert result.source.input_rules == "input_rules_v2"
    assert biotsr.api.processing_identity("protein")["input_rules"] == "input_rules_v2"


# ---------------------------------------------------------------------------
# ADR-112: a request without its chain is skipped, never guessed
# ---------------------------------------------------------------------------


def test_a_manifest_row_without_a_chain_is_skipped_as_chain_not_specified(tmp_path):
    manifest = tmp_path / "rows.csv"
    with open(str(manifest), "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["pdb_path", "chain", "method"])
        writer.writerow([structure("1UBQ"), "A", "protein"])
        writer.writerow([structure("1UBQ"), "", "protein"])
        writer.writerow([structure("1BRS"), "A", "protein_protein_calpha_cross"])
    out = str(tmp_path / "out")
    code = main(["batch", "--manifest", str(manifest), "--output", out])
    assert code == 4
    with open(os.path.join(out, "batch_summary.csv")) as fh:
        rows = list(csv.DictReader(fh))
    statuses = [r["status"] for r in rows]
    assert statuses == ["ok", "chain_not_specified", "chain_not_specified"]
    assert "protein needs chain" in rows[1]["message"]
    assert "chain_2" in rows[2]["message"]


def test_an_explicit_chain_default_still_fills_blank_rows(tmp_path):
    manifest = tmp_path / "rows.csv"
    manifest.write_text("pdb_path,chain,method\n%s,,protein\n" % structure("1UBQ"))
    out = str(tmp_path / "out")
    assert main(["batch", "--manifest", str(manifest), "--output", out,
                 "--chain", "A"]) == 0


# ---------------------------------------------------------------------------
# Closing checks (1.0.0rc7): what each check reads, and what no opt-in relaxes
# ---------------------------------------------------------------------------


def _alternates(first, second):
    """GLY 1, ALA 2 with GLY 1's CA written twice, as ``first`` and ``second``."""
    lines = chain_of(["GLY", "ALA"])
    ca = lines[1]
    moved = ca[:30] + "%8.3f" % (float(ca[30:38]) + 0.3) + ca[38:]
    return lines[:1] + [ca[:16] + first + ca[17:], moved[:16] + second + moved[17:]] \
        + lines[2:]


@pytest.mark.parametrize("method,arguments", [
    ("protein", {}),
    ("aa_ind_specific", {"residues": "GLY_1"}),
    ("aa_ind", {"residue_types": "GLY"}),
])
def test_duplicates_are_judged_after_the_conformation_is_chosen(tmp_path, method, arguments):
    """CA of GLY 1 as altLoc A and B: B is set aside first, so one record of
    the atom remains and nothing is refused.  Blank and A are both kept by the
    blank-or-A rule, so they are two records of one atom and are refused."""
    kept = write(tmp_path, _alternates("A", "B"), "ab.pdb")
    result = biotsr.analyze(kept, method=method, chain="A", **arguments)
    assert result.counts.n_points > 0
    assert any(w.startswith("alternate locations: kept blank and A; excluded 1")
               for w in result.warnings)
    both = write(tmp_path, _alternates(" ", "A"), "blank_a.pdb")
    with pytest.raises(InvalidRequest, match="conflicting records of one atom"):
        biotsr.analyze(both, method=method, chain="A", **arguments)


def test_an_unreadable_number_is_refused_only_where_the_selection_reads(tmp_path):
    """A hybrid-36 water after the chain's TER, and a hybrid-36 ligand on the
    same chain: the Calpha chain, another ligand and the protein's residues are
    computed; a selection that would read the unreadable records is refused."""
    water = atom(" O  ", "HOH", "A", "A00A", (30.0, 0.0, 0.0), " O", record="HETATM")
    ligand = [atom(" C%d " % i, "LIG", "A", "A00B", (20.0 + i, 3.0, 0.0), " C",
                   record="HETATM") for i in (1, 2, 3)]
    other = [atom(" C%d " % i, "BEN", "A", "  90", (20.0 + i, -3.0, 0.0), " C",
                  record="HETATM") for i in (1, 2, 3)]
    path = write(tmp_path, chain_of(["GLY", "ALA", "SER"]) + [ter("A")] + [water]
                 + ligand + other)
    assert len(biotsr.load_structure(path).unreadable_numbers) == 4
    assert biotsr.analyze(path, method="protein", chain="A").counts.n_points == 3
    assert biotsr.analyze(path, method="aa_ind", chain="A",
                          residue_types="GLY").counts.n_points == 4
    assert biotsr.analyze(path, method="drug", chain="A",
                          entity_name="BEN").counts.n_points == 3
    with pytest.raises(InvalidRequest, match="cannot read"):
        biotsr.analyze(path, method="drug", chain="A", entity_name="LIG")


@pytest.mark.parametrize("method,arguments", [
    ("protein", {}),
    ("aa_selection", {"residues": "1-3"}),
])
def test_the_calpha_blank_element_exception_is_unchanged(tmp_path, method, arguments):
    """A supported modified residue whose CA has a blank element field is a
    Calpha point (ADR-107); the atomic methods, which label by element, refuse
    a blank element they would select (ADR-108)."""
    mse = [atom(" N  ", "MSE", "A", 3, (7.6, 0.0, 0.0), " N", record="HETATM"),
           atom(" CA ", "MSE", "A", 3, (8.8, 0.9, 0.0), "  ", record="HETATM"),
           atom(" C  ", "MSE", "A", 3, (10.0, 0.1, 0.2), " C", record="HETATM")]
    path = write(tmp_path, chain_of(["GLY", "ALA"]) + mse)
    result = biotsr.analyze(path, method=method, chain="A", **arguments)
    assert result.selection.residues_used == ("1", "2", "3")
    assert result.key_space.label_scheme.name == "amino_acid_30"
    with pytest.raises(InvalidRequest, match="blank element field"):
        biotsr.analyze(path, method="drug", chain="A", entity_name="MSE")


def test_no_opt_in_compares_atomic_and_calpha_key_spaces():
    from biotsr.analysis import AnalysisSample, common_key_intersection
    from biotsr.domain.keyspace import IncompatibleKeySpace

    atomic = AnalysisSample.from_result(biotsr.analyze(
        structure("1UBQ"), method="aa_ind", chain="A", residue_types="LEU",
        keep_triplets=False))
    calpha = AnalysisSample.from_result(biotsr.analyze(
        structure("1UBQ"), method="protein", chain="A", keep_triplets=False))
    for kwargs in ({}, {"allow_mixed_processing": True, "allow_partial": True}):
        with pytest.raises(IncompatibleKeySpace):
            common_key_intersection([atomic, calpha], **kwargs)
