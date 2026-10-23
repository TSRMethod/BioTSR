"""Atomic labels come from the element field (ADR-113).

The expected labels below are written out by element -- carbon 106, nitrogen
108, oxygen 109, phosphorus 110, sulfur 111, chlorine 93, calcium 99,
manganese 100, holmium 37 -- and are **not** read from the package's table, so
the table is checked against them rather than with itself.

The records are real: each line is copied byte for byte from the RCSB entry
named beside it (model 1).  They are written into one test file, each as its
own residue, which is the only thing synthetic about it.  Before 1.0.0rc7 the
label came from the atom name, through the historical atom lexicon, and the
comment on each record gives the label it got then.
"""

from __future__ import annotations

import ast
import hashlib
import pathlib

import pytest

import biotsr
from biotsr.labels import elements
from biotsr.labels.elements import (
    ATOMIC_LABEL_RADIX,
    ELEMENT_LABELS_112,
    ELEMENT_LABELS_RESOURCE,
    ELEMENT_RESOLVER,
    ElementTableError,
    parse_element_labels,
)
from biotsr.labels.schemes import LabelDomain, LabelSource
from biotsr.methods.registry import get_method
from biotsr.profiles.atomic import ATOMIC_PROFILES
from biotsr.profiles.cross import STANDARD_ATOMIC_CROSS_V1
from biotsr.selection.atomic import STANDARD_HETERO_ADMISSION, select_entity_atoms

SRC_ROOT = pathlib.Path(biotsr.__file__).parent

#: The packaged table, by content.  A change to it is a change of every atomic
#: key space and must come with a new scheme and a changelog entry.
RESOURCE_SHA256 = "d7d8e489ffbb5dab7a433a61d55dd1a1b1680ef56cbfd67764a7b8971c335d9f"

#: Independent expectations, by element symbol as the PDB writes it.
BY_ELEMENT = {"C": 106, "N": 108, "O": 109, "P": 110, "S": 111, "CL": 93, "CA": 99,
              "MN": 100, "HO": 37}

#: (entry, record, label before 1.0.0rc7 -- None: could not be labelled)
RECORDS = [
    ("1EVZ", "HETATM 2604  PA  NAD A 400      36.226  47.383  22.318  1.00 61.39           P  ", 55),
    ("1EVZ", "HETATM 2622  N1A NAD A 400      33.143  54.146  30.032  1.00 75.42           N  ", 107),
    ("1UFJ", "HETATM 1276  CB1 CZM A 180      -1.249   2.719   1.899  1.00 24.03           C  ", 107),
    ("4X6J", "HETATM 1775 CL2  3Y2 A 307      47.301  22.874  21.171  1.00 27.26          CL  ", 106),
    ("3FSF", "HETATM 2785 CL35 FSS A 450      19.426  14.165  29.412  1.00 43.93          CL  ", None),
    ("1RER", "HETATM 9171 HO    HO A 411      60.702 -26.581  68.384  1.00 18.63          HO  ", 37),
    ("1ZHM", "HETATM 1929  CB2 GSH A1002       7.682  28.957  35.665  1.00 31.40           C  ", 107),
    ("1ZHM", "HETATM 1930  SG2 GSH A1002       8.790  27.564  35.213  1.00 31.45           S  ", 111),
    ("2AXT", "HETATM41130 MN1  OEC A 565      26.316  38.265  68.065  1.00 63.84          MN  ", 100),
    ("2AXT", "HETATM41134 CA1  OEC A 565      29.091  39.753  66.828  1.00 62.59          CA  ", 106),
    ("1BXQ", "HETATM 2440  CR  PP8 A 324      18.950   7.001  19.711  1.00  8.81           C  ", 22),
    ("1BXQ", "HETATM 2455  CL1 PP8 A 324      24.401   7.533  17.885  1.00  4.69           C  ", 93),
    ("1BXQ", "HETATM 2456  CL2 PP8 A 324      24.164   7.013  20.300  1.00  6.45           C  ", 106),
]

#: ATP 1388 of chain A of 1OL6, every heavy atom.  Its alpha and beta
#: phosphorus atoms, PA and PB, were labelled 55 (protactinium) and 56 (lead)
#: by name before 1.0.0rc7.
ATP_1OL6 = """\
HETATM 1847  PG  ATP A1388     209.746  28.040-102.392  1.00 92.20           P
HETATM 1848  O1G ATP A1388     210.725  28.862-103.148  1.00 87.86           O
HETATM 1849  O2G ATP A1388     208.564  28.830-101.874  1.00 86.46           O
HETATM 1850  O3G ATP A1388     210.406  27.473-101.247  1.00 88.27           O
HETATM 1851  PB  ATP A1388     209.841  25.357-103.825  1.00 85.16           P
HETATM 1852  O1B ATP A1388     210.262  25.368-105.232  1.00 84.82           O
HETATM 1853  O2B ATP A1388     210.819  24.932-102.833  1.00 83.82           O
HETATM 1854  O3B ATP A1388     209.260  26.795-103.376  1.00 85.30           O
HETATM 1855  PA  ATP A1388     207.157  24.685-103.046  1.00 73.06           P
HETATM 1856  O1A ATP A1388     207.212  24.983-101.607  1.00 72.92           O
HETATM 1857  O2A ATP A1388     206.477  25.701-103.825  1.00 71.11           O
HETATM 1858  O3A ATP A1388     208.620  24.422-103.627  1.00 77.47           O
HETATM 1859  O5' ATP A1388     206.500  23.262-103.300  1.00 68.55           O
HETATM 1860  C5' ATP A1388     206.362  22.909-104.714  1.00 57.30           C
HETATM 1861  C4' ATP A1388     204.929  22.513-105.058  1.00 52.44           C
HETATM 1862  O4' ATP A1388     204.617  21.307-104.295  1.00 52.05           O
HETATM 1863  C3' ATP A1388     203.836  23.484-104.623  1.00 49.17           C
HETATM 1864  O3' ATP A1388     203.605  24.474-105.549  1.00 50.48           O
HETATM 1865  C2' ATP A1388     202.635  22.557-104.339  1.00 47.17           C
HETATM 1866  O2' ATP A1388     201.860  22.192-105.451  1.00 43.99           O
HETATM 1867  C1' ATP A1388     203.292  21.348-103.790  1.00 48.11           C
HETATM 1868  N9  ATP A1388     203.307  21.308-102.317  1.00 47.58           N
HETATM 1869  C8  ATP A1388     204.261  22.009-101.585  1.00 46.08           C
HETATM 1870  N7  ATP A1388     203.974  22.131-100.298  1.00 47.49           N
HETATM 1871  C5  ATP A1388     202.790  21.409-100.150  1.00 47.28           C
HETATM 1872  C6  ATP A1388     201.951  21.042 -99.006  1.00 48.17           C
HETATM 1873  N6  ATP A1388     202.236  21.280 -97.730  1.00 49.18           N
HETATM 1874  N1  ATP A1388     200.767  20.435 -99.293  1.00 47.30           N
HETATM 1875  C2  ATP A1388     200.468  20.145-100.587  1.00 48.26           C
HETATM 1876  N3  ATP A1388     201.298  20.127-101.620  1.00 48.71           N
HETATM 1877  C4  ATP A1388     202.333  20.984-101.412  1.00 47.49           C
"""


def _write(path, lines):
    path.write_text("".join(l.rstrip("\n").ljust(80) + "\n" for l in lines) + "END\n")
    return str(path)


def _label_of(path, record):
    structure = biotsr.load_structure(path)
    selection = select_entity_atoms(
        structure, "A", ELEMENT_LABELS_112, admission=STANDARD_HETERO_ADMISSION,
        entity_names=(record[17:20].strip(),), entity_seqs=(record[22:27].strip(),),
        atom_names=(record[12:16].strip(),))
    (point,) = selection.points
    return point.label


# ===========================================================================
# The table
# ===========================================================================


def test_the_table_has_110_elements_and_109_labels_in_one_radix_window():
    table = ELEMENT_LABELS_112.table
    assert len(table) == 110
    assert len(set(table.values())) == 109
    assert sorted(set(table.values())) == list(range(4, 113))
    assert ELEMENT_LABELS_112.n_labels == ATOMIC_LABEL_RADIX == 112
    assert ELEMENT_LABELS_112.encoding_is_injective


def test_the_existing_element_codes_are_kept():
    """No code was renumbered when names stopped choosing labels."""
    table = ELEMENT_LABELS_112.table
    assert (table["C"], table["N"], table["O"], table["P"], table["S"]) == (106, 108, 109, 110, 111)
    assert (table["Cl"], table["Ho"], table["Ca"], table["Sg"], table["Zn"]) == (93, 37, 99, 72, 105)
    assert table["H"] == table["D"] == 107


def test_hydrogen_and_deuterium_are_the_only_shared_label():
    by_label = {}
    for element, label in ELEMENT_LABELS_112.table.items():
        by_label.setdefault(label, []).append(element)
    assert {k: sorted(v) for k, v in by_label.items() if len(v) > 1} == {107: ["D", "H"]}
    assert elements.SHARED_LABELS == (("H", "D"),)


def test_the_superheavy_elements_the_lexicon_never_listed_have_no_label():
    for symbol in ("Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"):
        assert symbol not in ELEMENT_LABELS_112.table


def test_the_packaged_table_is_pinned_by_content():
    text = elements._read_resource(ELEMENT_LABELS_RESOURCE)
    assert hashlib.sha256(text.encode("utf-8")).hexdigest() == RESOURCE_SHA256
    assert text.splitlines()[0] == "element,label"
    assert RESOURCE_SHA256[:12] in ELEMENT_LABELS_112.description


def test_only_the_element_table_ships_and_the_lexicon_is_gone():
    resources = sorted(p.name for p in (SRC_ROOT / "labels" / "resources").iterdir()
                       if not p.name.startswith("."))
    assert resources == [ELEMENT_LABELS_RESOURCE]
    import importlib

    with pytest.raises(ImportError):
        importlib.import_module("biotsr.labels.lexicon")
    assert not hasattr(biotsr, "ATOM_LEXICON_STANDARD_112")
    assert biotsr.ELEMENT_LABELS_112 is ELEMENT_LABELS_112


def test_the_scheme_records_its_source_and_resolver_in_the_fingerprint():
    assert ELEMENT_LABELS_112.domain == LabelDomain.ATOM_ELEMENT
    assert ELEMENT_LABELS_112.label_source == LabelSource.ELEMENT
    assert ELEMENT_LABELS_112.resolver == ELEMENT_RESOLVER == "element_field_v1"
    import dataclasses

    other = dataclasses.replace(ELEMENT_LABELS_112, resolver="element_field_v2")
    assert other.content_fingerprint() != ELEMENT_LABELS_112.content_fingerprint()


# ===========================================================================
# The loader refuses a malformed table
# ===========================================================================

_GOOD = "element,label\nC,106\nN,108\nH,107\nD,107\n"


def test_a_well_formed_table_parses():
    assert parse_element_labels(_GOOD) == {"C": 106, "N": 108, "H": 107, "D": 107}


@pytest.mark.parametrize("text,words", [
    ("atom,seq,ATOM\nC,106,C\n", "first line must be exactly"),
    ("element,label,ATOM\nC,106,C\n", "first line must be exactly"),
    ("", "first line must be exactly"),
    ("element,label\n", "lists no elements"),
    ("element,label\nC,106,CA\n", "expected 2 fields"),
    ("element,label\nC\n", "expected 2 fields"),
    ("element,label\nCL,93\n", "not an element symbol in chemical case"),
    ("element,label\ncl,93\n", "not an element symbol in chemical case"),
    ("element,label\nC1,106\n", "not an element symbol in chemical case"),
    ("element,label\n,106\n", "not an element symbol in chemical case"),
    ("element,label\nC,x\n", "not a positive integer"),
    ("element,label\nC,0\n", "not a positive integer"),
    ("element,label\nC,-3\n", "not a positive integer"),
    ("element,label\nC, 106\n", "not a positive integer"),
    ("element,label\nC,106\nC,108\n", "a conflicting duplicate"),
    ("element,label\nC,106\nC,106\n", "listed again"),
    ("element,label\nC,106\nN,106\n", "share label 106"),
    ("element,label\nH,107\nD,108\n", "H/D must share one label"),
])
def test_a_malformed_table_is_refused_naming_the_problem(text, words):
    with pytest.raises(ElementTableError, match=words.replace("(", r"\(")):
        parse_element_labels(text)


def test_a_conflicting_duplicate_names_both_lines():
    with pytest.raises(ElementTableError) as info:
        parse_element_labels("element,label\nC,106\nN,108\nC,93\n")
    assert "line 4" in str(info.value) and "first on line 2" in str(info.value)


# ===========================================================================
# The resolver reads the element, never the name
# ===========================================================================


@pytest.mark.parametrize("element", ["CL", "Cl", "cl", " CL", "CL "])
def test_the_element_field_is_read_case_insensitively(element):
    assert ELEMENT_LABELS_112.resolve_atom("CL1", element) == 93


@pytest.mark.parametrize("name", ["C1", "CA", "CL1", "PA", "SG", "N1A", "XYZ", ""])
def test_the_atom_name_never_changes_the_label(name):
    assert ELEMENT_LABELS_112.resolve_atom(name, "C") == 106
    assert ELEMENT_LABELS_112.resolve_atom(name, "P") == 110


@pytest.mark.parametrize("element", ["", " ", "X", "Q", "C1", "Uuo", "Ds"])
def test_a_blank_unknown_or_unsupported_element_has_no_label(element):
    """Nothing is inferred from the name: a CA atom with an unreadable element
    is not assumed to be carbon (or calcium)."""
    assert ELEMENT_LABELS_112.resolve_atom("CA", element) is None


def test_no_cys_sg_special_case_exists_anywhere_in_the_package():
    """No executable line mentions CYS and SG together: a label never depends on
    a residue-and-atom-name rule.  Docstrings and comments are excluded."""
    offenders = []
    for path in sorted(SRC_ROOT.rglob("*.py")):
        source = path.read_text()
        tree = ast.parse(source)
        blanked = set()
        for node in ast.walk(tree):
            body = getattr(node, "body", None)
            if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                                 ast.AsyncFunctionDef)) and body:
                first = body[0]
                if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) \
                        and isinstance(first.value.value, str):
                    blanked.update(range(first.lineno, (first.end_lineno or first.lineno) + 1))
        for n, line in enumerate(source.splitlines(), 1):
            code = "" if n in blanked else line.split("#", 1)[0]
            if ('"CYS"' in code or "'CYS'" in code) and ('"SG"' in code or "'SG'" in code):
                offenders.append("%s:%d" % (path.name, n))
    assert offenders == []


# ===========================================================================
# Real records
# ===========================================================================


@pytest.mark.parametrize("entry,record,before", RECORDS, ids=[
    "%s_%s_%s" % (e, r[17:20].strip(), r[12:16].strip()) for e, r, _ in RECORDS])
def test_a_real_record_is_labelled_by_its_element(entry, record, before, tmp_path):
    path = _write(tmp_path / "records.pdb", [r for _, r, _ in RECORDS])
    expected = BY_ELEMENT[record[76:78].strip()]
    assert _label_of(path, record) == expected
    # The record's element field decides it, whatever it was labelled before.
    assert before != expected or record[12:16].strip() in ("HO", "SG2", "MN1", "CL2")


def test_the_atp_phosphorus_atoms_are_phosphorus_in_a_drug_calculation(tmp_path):
    path = _write(tmp_path / "1OL6_ATP.pdb", ATP_1OL6.splitlines())
    result = biotsr.analyze(path, method="drug", chain="A", entity_name="ATP",
                            entity_id="1388")
    labels = {}
    for t in result.triplets:
        for identity, label in zip(t.identities, t.labels):
            labels.setdefault(identity.atom_name, set()).add(label)
    assert labels["PA"] == labels["PB"] == labels["PG"] == {110}
    for name, found in labels.items():
        assert found == {BY_ELEMENT[name[0]]}, name
    assert result.counts.n_points == 31 and result.counts.n_triangles == 4495
    names = {t.label_names for t in result.triplets}
    assert ("P", "P", "P") in names


def test_renaming_atoms_does_not_change_a_single_key(tmp_path):
    """The same atoms under other names: identical labels, identical keys."""
    original = _write(tmp_path / "a.pdb", ATP_1OL6.splitlines())
    renamed_lines = []
    for n, line in enumerate(ATP_1OL6.splitlines()):
        # PA -> CL1, PB -> CA, O1G -> N1A ...: names that, read as names, once
        # meant other elements.  The element field is untouched.
        name = ("CL%-2d" % n if n % 3 == 0 else "CA%-2d" % n if n % 3 == 1 else "N%-3d" % n)
        renamed_lines.append(line[:12] + name + line[16:])
    renamed = _write(tmp_path / "b.pdb", renamed_lines)
    a = biotsr.analyze(original, method="drug", chain="A", entity_name="ATP",
                       entity_id="1388", keep_triplets=False)
    b = biotsr.analyze(renamed, method="drug", chain="A", entity_name="ATP",
                       entity_id="1388", keep_triplets=False)
    assert dict(a.frequencies) == dict(b.frequencies)
    assert a.key_space.fingerprint == b.key_space.fingerprint


def test_hydrogen_and_deuterium_are_excluded_before_labelling(tmp_path):
    """They carry 107 in the table and never reach it: the admission drops them."""
    lines = ATP_1OL6.splitlines()[:5] + [
        "HETATM 1900  H1  ATP A1388     209.000  27.000-102.000  1.00 90.00           H",
        "HETATM 1901  D2  ATP A1388     208.000  26.000-101.000  1.00 90.00           D",
        "HETATM 1902  HX  ATP A1388     207.000  25.000-100.000  1.00 90.00           h",
    ]
    path = _write(tmp_path / "h.pdb", lines)
    result = biotsr.analyze(path, method="drug", chain="A", entity_name="ATP",
                            entity_id="1388")
    assert result.counts.n_points == 5
    assert all(107 not in t.labels for t in result.triplets)


# ===========================================================================
# One resolver for every atomic method
# ===========================================================================


def test_every_atomic_profile_and_method_shares_the_one_scheme_object():
    profiles = list(ATOMIC_PROFILES) + [STANDARD_ATOMIC_CROSS_V1]
    assert all(p.label_scheme is ELEMENT_LABELS_112 for p in profiles)
    atomic = [m for m in biotsr.list_methods()
              if get_method(m).default_profile.label_scheme.domain == LabelDomain.ATOM_ELEMENT]
    assert len(atomic) == 17
    for name in atomic:
        assert get_method(name).default_profile.label_scheme is ELEMENT_LABELS_112, name


def test_the_calpha_schemes_are_untouched():
    from biotsr.labels import AMINO_ACID_30, AMINO_ACID_GROUP_23

    assert AMINO_ACID_30.label_source == LabelSource.RESIDUE_NAME
    assert AMINO_ACID_30.content_fingerprint() == \
        "d6fbfd306d14ec59720c1fe189cfac3fa30c12d6bdf2139b4a56f5445eb81d92"
    assert AMINO_ACID_GROUP_23.content_fingerprint() == \
        "79686ef72c0030b64e1b4d7edf7c93718a6fc5acb1a36e873a286a4b8314ffef"
    for method in ("protein", "aa_selection", "protein_protein_calpha_cross"):
        assert get_method(method).default_profile.key_space().fingerprint == \
            "da39d73510c54a4c07529e17f34449789584fad4ea79b020b1badd3b5550d30b"
