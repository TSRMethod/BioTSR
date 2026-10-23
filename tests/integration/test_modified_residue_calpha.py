"""Modified residues in the protein Cα methods, at radix 30 (ADR-107).

* Hand-calculated synthetic triangles, through the whole API: the pair that
  shared a key at radix 20 -- SEP-PHE-PHE and TPO-TPO-PHE in the same angle and
  distance classes -- now has two keys, each computed by hand.
* Intentional grouping still merges what it should, and only that.
* Real structures: CME 161 of 6G9M chain A, and the synthetic SYNPHOS chains in
  the Cα methods and the Cα cross.
* A modified residue the residue labels do not cover is named, never read.
* An output of the previous build (radix 20, no modified residues) is never
  reused, and never compared with a current result.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import os
import shutil

import pytest

import biotsr
from biotsr.analysis import AnalysisSample, common_key_intersection
from biotsr.cli import EXIT_OK, EXIT_USAGE, main
from biotsr.domain.errors import InvalidRequest
from biotsr.domain.keyspace import IncompatibleKeySpace
from biotsr.io.native_read import load_result
from biotsr.selection.base import ModifiedResiduePolicy

from ..conftest import STRUCTURES, SYNTHETIC

HERE = os.path.dirname(os.path.abspath(__file__))
SYN = os.path.join(SYNTHETIC, "SYNPHOS.pdb")
PREVIOUS = os.path.join(HERE, "..", "data", "legacy_outputs", "before_adr107_calpha_radix20")
PREFIX = "SYNPHOS_A_protein"
REQUEST = ["--method", "protein", "--chain", "A", "--sample-id", PREFIX, "--no-triplets"]

#: The synthetic geometry: base from (0,0,0) to (8,0,0), apex (3.5,4,0). The
#: median from the base midpoint (4,0,0) to the apex makes 82.875 degrees with
#: the base -- theta class 26 (81.36 <= theta < 83.51) -- and the longest edge is
#: 8.0 A, distance class 3 (7.00 <= d < 9.00).
GEOMETRY = ((0.0, 0.0, 0.0), (8.0, 0.0, 0.0), (3.5, 4.0, 0.0))


def residue_lines(chain, seq, name, xyz, serial, hetero):
    """N, CA and C of one residue, the CA at ``xyz``."""
    record = "HETATM" if hetero else "ATOM  "
    x, y, z = xyz
    out = []
    for n, (atom, dx, element) in enumerate(((" N  ", -1.2, "N"), (" CA ", 0.0, "C"),
                                             (" C  ", 1.2, "C"))):
        out.append("%-6s%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n"
                   % (record, serial + n, atom, name, chain, seq, x + dx, y, z + 1.5 if
                      atom != " CA " else z, element))
    return out


def write_pdb(path, chains):
    """``chains``: {chain: [(name, xyz), ...]}; modified residues as HETATM."""
    lines, serial = [], 1
    for chain, residues in chains.items():
        for seq, (name, xyz) in enumerate(residues, 1):
            hetero = name not in biotsr.labels.STANDARD_AMINO_ACIDS
            lines += residue_lines(chain, seq, name, xyz, serial, hetero)
            serial += 3
        lines.append("TER   %5d      %3s %1s%4d\n" % (serial, residues[-1][0], chain,
                                                      len(residues)))
        serial += 1
    lines.append("END\n")
    with open(str(path), "w") as fh:
        fh.writelines(lines)
    return str(path)


def key_by_hand(c1, c2, c3, theta_class, distance_class, radix=30):
    return (35 * 29 * (radix ** 2 * (c1 - 1) + radix * (c2 - 1) + (c3 - 1))
            + 29 * (distance_class - 1) + (theta_class - 1))


@pytest.fixture
def pair(tmp_path):
    """Chain A: SEP-PHE-PHE; chain B: TPO-TPO-PHE; the same three positions."""
    p0, p1, p2 = GEOMETRY
    return write_pdb(tmp_path / "PAIR.pdb", {
        "A": [("SEP", p0), ("PHE", p1), ("PHE", p2)],
        "B": [("TPO", p0), ("TPO", p1), ("PHE", p2)],
    })


# ===========================================================================
# Hand-calculated keys, through the whole API
# ===========================================================================


def test_the_pair_that_collided_at_radix_20_has_two_hand_calculated_keys(pair):
    a = biotsr.analyze(pair, method="protein", chain="A")
    b = biotsr.analyze(pair, method="protein", chain="B")
    [ta], [tb] = a.triplets, b.triplets
    assert ta.labels == (25, 4, 4) and tb.labels == (24, 24, 4)
    assert (ta.theta_class, ta.distance_class) == (tb.theta_class, tb.distance_class) == (26, 3)
    assert ta.theta == pytest.approx(82.875, abs=1e-3)
    assert ta.key == key_by_hand(25, 4, 4, 26, 3) == 22018478
    assert tb.key == key_by_hand(24, 24, 4, 26, 3) == 21713978
    # At radix 20 the same labels and classes give one number (9663 * 1015 + 83).
    assert key_by_hand(25, 4, 4, 26, 3, radix=20) == key_by_hand(24, 24, 4, 26, 3, radix=20)
    assert a.key_space.fingerprint == b.key_space.fingerprint
    assert set(a.frequencies) & set(b.frequencies) == set()


def test_the_two_triangles_in_one_chain_keep_their_own_counts(tmp_path):
    p0, p1, p2 = GEOMETRY
    far = [(x + 200.0, y, z) for x, y, z in GEOMETRY]
    path = write_pdb(tmp_path / "ONE.pdb", {"A": [
        ("SEP", p0), ("PHE", p1), ("PHE", p2),
        ("TPO", far[0]), ("TPO", far[1]), ("PHE", far[2])]})
    result = biotsr.analyze(path, method="protein", chain="A", keep_triplets=False)
    assert result.counts.n_points == 6 and result.n_triangles == 20
    assert result.frequencies[22018478] == 1
    assert result.frequencies[21713978] == 1


def test_grouping_separates_the_modified_residues_and_merges_a_group(tmp_path):
    """Grouped: SEP (25) and TPO (24) stay apart; SER and THR (12) share a code
    on purpose, so SER-PHE-PHE and THR-PHE-PHE share a key."""
    p0, p1, p2 = GEOMETRY
    path = write_pdb(tmp_path / "GROUP.pdb", {
        "A": [("SEP", p0), ("PHE", p1), ("PHE", p2)],
        "B": [("TPO", p0), ("TPO", p1), ("PHE", p2)],
        "C": [("SER", p0), ("PHE", p1), ("PHE", p2)],
        "D": [("THR", p0), ("PHE", p1), ("PHE", p2)],
    })
    grouped = {c: biotsr.analyze(path, method="protein", chain=c, aa_grouping=True)
               for c in "ABCD"}
    keys = {c: r.triplets[0].key for c, r in grouped.items()}
    # PHE is group 10 in the grouped table.
    assert keys["A"] == key_by_hand(25, 10, 10, 26, 3)
    assert keys["B"] == key_by_hand(24, 24, 10, 26, 3)
    assert keys["C"] == keys["D"] == key_by_hand(12, 10, 10, 26, 3)
    assert len({keys["A"], keys["B"], keys["C"]}) == 3
    plain = {c: biotsr.analyze(path, method="protein", chain=c).triplets[0].key
             for c in "CD"}
    assert plain["C"] != plain["D"]
    assert grouped["A"].key_space.fingerprint != biotsr.analyze(
        path, method="protein", chain="A").key_space.fingerprint


# ===========================================================================
# Real and synthetic structures
# ===========================================================================


def test_cme_161_of_6g9m_is_a_calpha_point_labelled_32():
    result = biotsr.analyze(os.path.join(STRUCTURES, "6G9M.pdb"), method="aa_selection",
                            chain="A", residues="159-163")
    assert result.selection.residues_used == ("159", "160", "161", "162", "163")
    labels = {i.residue_name: l for t in result.triplets
              for i, l in zip(t.identities, t.labels)}
    assert labels["CME"] == 32
    assert result.warnings == ()


def test_every_modified_residue_of_synphos_chain_a_is_a_point():
    result = biotsr.analyze(SYN, method="protein", chain="A")
    assert result.counts.n_points == 14 and result.n_triangles == 364
    labels = {i.residue_name: l for t in result.triplets
              for i, l in zip(t.identities, t.labels)}
    assert {k: labels[k] for k in ("TPO", "SEP", "PTR", "CME")} == {
        "TPO": 24, "SEP": 25, "PTR": 26, "CME": 32}
    # The ligand and the calcium ion after TER, whose atom is also named CA,
    # are not points.
    assert "LIG" not in labels and "CA" not in labels
    assert result.warnings == ()
    assert result.source.modified_residue_policy == ModifiedResiduePolicy.CALPHA
    record = result.annotations["selection_request"]["sides"][0]
    assert "unsupported" not in record


def test_the_calpha_cross_reads_modified_residues_on_both_sides():
    result = biotsr.analyze_cross(SYN, method="protein_protein_calpha_cross",
                                  chain="A", chain_2="B")
    points = {(i.chain.value, i.residue_name, l) for t in result.triplets
              for i, l in zip(t.identities, t.labels)}
    for expected in (("A", "TPO", 24), ("A", "SEP", 25), ("A", "PTR", 26),
                     ("B", "TPO", 24)):
        assert expected in points
    assert result.warnings == ()
    assert result.source.modified_residue_policy == ModifiedResiduePolicy.CALPHA


def test_the_atomic_methods_still_read_only_tpo_sep_and_ptr():
    mixed = biotsr.analyze(SYN, method="aa_ind_mix", chain="A",
                           residue_types="TYR;TPO;CME", keep_triplets=False)
    assert mixed.annotations["selection_request"]["sides"][0]["unsupported"] == ["CME"]
    assert mixed.source.modified_residue_policy == ModifiedResiduePolicy.ATOMIC


def test_a_modified_residue_without_a_label_is_named_not_read(tmp_path):
    p0, p1, p2 = GEOMETRY
    path = write_pdb(tmp_path / "HYP.pdb", {"A": [
        ("GLY", p0), ("HYP", p1), ("ALA", p2), ("LEU", (10.0, 0.0, 0.0))]})
    result = biotsr.analyze(path, method="protein", chain="A")
    assert result.selection.residues_used == ("1", "3", "4")
    assert result.warnings[0].startswith("not represented: HYP_2 (")
    record = result.annotations["selection_request"]["sides"][0]
    assert record["unsupported"] == ["HYP_2"]


# ===========================================================================
# An output of the previous build is never reused nor compared
# ===========================================================================


def run(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main([str(a) for a in argv])
    return code, out.getvalue(), err.getvalue()


def digest(folder):
    return {name: hashlib.sha256(open(os.path.join(folder, name), "rb").read()).hexdigest()
            for name in sorted(os.listdir(folder))}


@pytest.fixture
def previous(tmp_path):
    shutil.copy(SYN, str(tmp_path))
    shutil.copytree(PREVIOUS, str(tmp_path / "out"),
                    ignore=shutil.ignore_patterns("README.md"))
    return tmp_path


def test_the_previous_output_is_unedited_and_is_radix_20():
    with open(os.path.join(PREVIOUS, "README.md")) as fh:
        readme = fh.read()
    for name, sha in digest(PREVIOUS).items():
        if name != "README.md":
            assert sha in readme, name
    result = load_result(os.path.join(PREVIOUS, PREFIX + ".keys"))
    assert result.key_space.label_scheme.name == "amino_acid_20"
    assert result.key_space.label_scheme.radix == 20
    assert result.counts.n_points == 10 and result.n_triangles == 120
    assert result.source.modified_residue_policy == ModifiedResiduePolicy.ATOMIC


def test_skip_refuses_the_previous_output_and_leaves_it_untouched(previous):
    before = digest(str(previous / "out"))
    code, _, err = run(["run", "--pdb", previous / "SYNPHOS.pdb", "--output",
                        previous / "out", "--if-exists", "skip"] + REQUEST)
    assert code == EXIT_USAGE, err
    assert "different processing rules" in err
    assert ModifiedResiduePolicy.CALPHA in err
    assert digest(str(previous / "out")) == before


def test_skip_refuses_a_calpha_output_of_another_key_space_even_under_current_rules(previous):
    """The processing rule edited to the current one: the key space still
    tells the output apart."""
    path = previous / "out" / (PREFIX + ".run.json")
    with open(str(path)) as fh:
        record = json.load(fh)
    record["request"]["processing"]["modified_residue_policy"] = ModifiedResiduePolicy.CALPHA
    record["request"]["processing"]["input_rules"] = biotsr.api.processing_identity(
        "protein")["input_rules"]
    with open(str(path), "w") as fh:
        json.dump(record, fh)
    code, _, err = run(["run", "--pdb", previous / "SYNPHOS.pdb", "--output",
                        previous / "out", "--if-exists", "skip"] + REQUEST)
    assert code == EXIT_USAGE, err
    assert "different key space" in err


def test_overwrite_recomputes_it_and_then_skip_keeps_it(previous):
    args = ["run", "--pdb", previous / "SYNPHOS.pdb", "--output", previous / "out"] + REQUEST
    code, _, err = run(args + ["--if-exists", "overwrite"])
    assert code == EXIT_OK, err
    result = load_result(str(previous / "out" / (PREFIX + ".keys")))
    assert result.counts.n_points == 14 and result.key_space.label_scheme.radix == 30
    code, out, err = run(args + ["--if-exists", "skip"])
    assert code == EXIT_OK, err


def test_no_analysis_compares_the_previous_output_with_a_current_one(tmp_path):
    old = AnalysisSample.from_result(load_result(os.path.join(PREVIOUS, PREFIX + ".keys")))
    old = AnalysisSample(sample_id="previous", counts=old.counts, key_space=old.key_space,
                         sampling=old.sampling, annotations=old.annotations)
    new = AnalysisSample.from_result(biotsr.analyze(SYN, method="protein", chain="A"))
    for kwargs in ({}, {"allow_mixed_processing": True, "allow_partial": True}):
        with pytest.raises(IncompatibleKeySpace):
            common_key_intersection([old, new], **kwargs)


# ===========================================================================
# The Ca of a supported residue is read on its own (ADR-107, corrected): every
# Calpha method, and every residue left out is named with its reason
# ===========================================================================

_OFFSET = {"N": (-1.2, 0.5), "CA": (0.0, 0.0), "C": (1.2, 0.5), "O": (1.4, 1.6),
           "CB": (0.0, -1.5), "OG": (0.5, -2.6), "P": (1.0, -3.8)}


def _residue(record, name, chain, seq, x, y, serial, atoms, elements=None, alts=None):
    out = []
    for n, atom in enumerate(atoms):
        dx, dy = _OFFSET[atom]
        element = (elements or {}).get(atom, atom[0])
        out.append("%-6s%5d  %-3s%1s%3s %1s%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n"
                   % (record, serial + n, atom, (alts or {}).get(atom, " "), name, chain,
                      seq, x + dx, y + dy, 0.0, element))
    return out


def _fixture(path, record, name, atoms, elements=None, alts=None):
    """Chain A: GLY 1, ALA 2, <the case residue 3>, LEU 4, VAL 5; chain B: three
    residues 4.5 A away, so every chain-A residue is at the interface."""
    lines, serial = [], 1
    for seq, res in ((1, "GLY"), (2, "ALA"), (3, None), (4, "LEU"), (5, "VAL")):
        if res is None:
            lines += _residue(record, name, "A", 3, 7.6, 0.0, serial, atoms, elements, alts)
        else:
            lines += _residue("ATOM", res, "A", seq, 3.8 * (seq - 1), 0.0, serial,
                              ("N", "CA", "C"))
        serial += 10
    lines.append("TER   %5d      VAL A   5\n" % serial)
    for seq, res in ((1, "LYS"), (2, "GLU"), (3, "ASP")):
        lines += _residue("ATOM", res, "B", seq, 3.8 * seq, 4.5, serial + seq * 10,
                          ("N", "CA", "C"))
    lines.append("END\n")
    with open(str(path), "w") as fh:
        fh.writelines(lines)
    return str(path)


_SIDE = ("CB", "OG", "P")
#: case -> (record, name, atoms, elements, alts, expected reason or None if read)
BACKBONE_CASES = {
    "complete": ("HETATM", "SEP", ("N", "CA", "C", "O") + _SIDE, None, None, None),
    "ca_only": ("HETATM", "SEP", ("CA",), None, None, None),
    "missing_n": ("HETATM", "SEP", ("CA", "C", "O") + _SIDE, None, None, None),
    "missing_c": ("HETATM", "SEP", ("N", "CA") + _SIDE, None, None, None),
    "blank_elements": ("HETATM", "SEP", ("N", "CA", "C", "O"),
                       {a: "" for a in ("N", "CA", "C", "O")}, None, None),
    "n_only_as_altloc_b": ("HETATM", "SEP", ("N", "CA", "C"), None, {"N": "B"}, None),
    "atom_record_ca_only": ("ATOM", "SEP", ("CA",), None, None, None),
    "ca_declared_calcium": ("HETATM", "SEP", ("N", "CA", "C"), {"CA": "CA"}, None,
                            "its CA atom's element field is 'CA', not carbon"),
    "ca_only_as_altloc_b": ("HETATM", "SEP", ("N", "CA", "C"), None, {"CA": "B"},
                            "its CA atom has only alternate-location indicator(s) B"),
    "no_ca": ("HETATM", "SEP", ("N", "C", "O") + _SIDE, None, None, "no CA atom"),
}


@pytest.mark.parametrize("case", sorted(BACKBONE_CASES))
def test_every_calpha_method_reads_or_names_the_supported_residue(case, tmp_path):
    record, name, atoms, elements, alts, reason = BACKBONE_CASES[case]
    path = _fixture(tmp_path / "CASE.pdb", record, name, atoms, elements, alts)
    chain = biotsr.analyze(path, method="protein", chain="A", keep_triplets=False)
    listed = biotsr.analyze(path, method="aa_selection", chain="A", residues="1-5",
                            keep_triplets=False)
    if elements and any(not e for e in elements.values()):
        # A blank-element HETATM CA is still a Calpha point (ADR-107), but the
        # Calpha cross also reads the residue's heavy atoms for contact, and a
        # blank element there cannot be told from a hydrogen (ADR-108).
        with pytest.raises(InvalidRequest, match="blank element field"):
            biotsr.analyze_cross(path, method="protein_protein_calpha_cross",
                                 chain="A", chain_2="B")
        assert "3" in chain.selection.residues_used
        assert "3" in listed.selection.residues_used
        return
    cross = biotsr.analyze_cross(path, method="protein_protein_calpha_cross",
                                 chain="A", chain_2="B")
    cross_points = {i.residue.legacy_seq_label for t in cross.triplets
                    for i in t.identities if i.chain.value == "A"}
    for result, side in ((chain, chain.annotations["selection_request"]["sides"][0]),
                         (listed, listed.annotations["selection_request"]["sides"][0]),
                         (cross, cross.annotations["selection_request"]["sides"][0])):
        assert "not_found" not in side, case          # present: never "not found"
        if reason is None:
            assert "unsupported" not in side, case
            assert not any(w.startswith("not represented") for w in result.warnings), case
        else:
            assert side["unsupported"] == ["SEP_3"], case
            assert any(w.startswith("not represented: SEP_3 (%s" % reason)
                       for w in result.warnings), (case, result.warnings)
    used = "3" in chain.selection.residues_used
    assert used == (reason is None), case
    assert ("3" in listed.selection.residues_used) == (reason is None), case
    assert ("3" in cross_points) == (reason is None), case
    if reason is not None:
        # The interface rule still sees the residue: it is recorded as an
        # interface residue without a Ca, never replaced by another atom.
        assert cross.annotations["cross_pools"]["n_interface_residues_without_calpha"] == 1
        assert AnalysisSample.from_result(chain).omissions["primary"] == {
            "unsupported": ["SEP_3"]}


def test_a_calcium_ion_is_never_a_point_and_a_requested_one_is_not_called_missing(tmp_path):
    path = _fixture(tmp_path / "ION.pdb", "HETATM", "CA", ("CA",), {"CA": "CA"})
    chain = biotsr.analyze(path, method="protein", chain="A", keep_triplets=False)
    assert chain.selection.residues_used == ("1", "2", "4", "5")
    assert chain.warnings == ()                    # an ion is not part of the chain's residues
    listed = biotsr.analyze(path, method="aa_selection", chain="A", residues="1-5",
                            keep_triplets=False)
    side = listed.annotations["selection_request"]["sides"][0]
    assert side["unsupported"] == ["CA_3"] and "not_found" not in side
    assert listed.warnings[0].startswith(
        "not represented: CA_3 (a HETATM record that is not an amino-acid residue")
    absent = biotsr.analyze(path, method="aa_selection", chain="A", residues="1-2;9",
                            keep_triplets=False)
    assert absent.annotations["selection_request"]["sides"][0]["not_found"] == ["9"]


def test_an_output_of_the_backbone_rule_is_not_reused(tmp_path):
    """An earlier 1.0.0rc5 build also required N and C and recorded the residue
    names alone; such an output is refused on resume, never kept."""
    args = ["run", "--pdb", SYN, "--output", tmp_path / "out"] + REQUEST
    assert run(args)[0] == EXIT_OK
    path = tmp_path / "out" / (PREFIX + ".run.json")
    with open(str(path)) as fh:
        record = json.load(fh)
    current = record["request"]["processing"]["modified_residue_policy"]
    assert current == ModifiedResiduePolicy.CALPHA
    assert current.endswith(";ca_element:C_or_blank")
    record["request"]["processing"]["modified_residue_policy"] = current.split(";")[0]
    with open(str(path), "w") as fh:
        json.dump(record, fh)
    before = digest(str(tmp_path / "out"))
    code, _, err = run(args + ["--if-exists", "skip"])
    assert code == EXIT_USAGE
    assert "different processing rules" in err and "choose another --output" in err
    assert digest(str(tmp_path / "out")) == before
