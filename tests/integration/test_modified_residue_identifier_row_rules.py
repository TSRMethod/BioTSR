"""Modified residues, residue identifiers, duplicate rows and cross sides, end to end.

* modified amino acids (ADR-101): TPO, SEP and PTR -- the modified residues
  with a defined atom count, labelled by element -- are read from ``HETATM``
  records by the amino-acid atomic methods and the atomic cross sides; any
  other modified residue is named as not supported, never silently dropped,
  and a calculation missing something it named is reported ``partial``.  The
  Calpha methods read the eleven modified residues of the residue labels
  (ADR-107), CME included;
* residue identifiers and ranges (ADR-102): one grammar for the API, the
  command line and a manifest; insertion-code ranges, descending ranges and a
  range spanning an insertion-code residue are refused with the remedy;
* duplicate rows (ADR-103): one calculation of one sample is computed once
  with every source row kept; distinct sample ids stay distinct samples; rows
  that would merge with different labels are refused;
* cross sides from a wide sheet (ADR-104): a side the row could narrow in
  several ways must be chosen (``--cross-side``), and the preview shows both
  sides of every request.

``tests/data/synthetic/SYNPHOS.pdb`` is synthetic (its README says how it was
written): chain A holds HETATM TPO 3, SEP 5, PTR 7 and CME 8 and the
insertion-code residues 10A and 10B; chain B holds a TPO; a ligand LIG 101
and a calcium ion follow TER.
"""

from __future__ import annotations

import contextlib
import csv
import io
import json
import os
import shutil

import pytest

import biotsr
import biotsr.api
import biotsr.batch
from biotsr.cli import EXIT_INCOMPLETE, EXIT_OK, EXIT_USAGE, main
from biotsr.domain.errors import InvalidRequest
from biotsr.selection.atomic import hetero_amino_acid_residues

from ..conftest import STRUCTURES, SYNTHETIC

SYN = os.path.join(SYNTHETIC, "SYNPHOS.pdb")


def run(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main([str(a) for a in argv])
    return code, out.getvalue(), err.getvalue()


def summary(folder):
    with open(os.path.join(str(folder), "batch_summary.csv"), newline="") as fh:
        return list(csv.DictReader(fh))


def side(result, index=0):
    return result.annotations["selection_request"]["sides"][index]


def analyze(**kwargs):
    return biotsr.analyze(SYN, keep_triplets=False, **kwargs)


def write(path, text):
    with open(str(path), "w", newline="") as fh:
        fh.write(text)
    return str(path)


# ===========================================================================
# Modified amino acids
# ===========================================================================


def test_the_synthetic_structure_is_exactly_what_its_generator_writes():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "make_synthetic", os.path.join(SYNTHETIC, "make_synthetic.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with open(SYN, newline="") as fh:
        assert fh.read() == module.build()


def test_the_structure_finder_names_modified_residues_not_ions_or_ligands():
    structure = biotsr.load_structure(SYN)
    assert hetero_amino_acid_residues(structure, "A") == (
        ("3", "TPO"), ("5", "SEP"), ("7", "PTR"), ("8", "CME"))
    assert hetero_amino_acid_residues(structure, "B") == (("3", "TPO"),)


@pytest.mark.parametrize("name, n_atoms", [("TPO", 11), ("SEP", 10), ("PTR", 16)])
def test_the_phosphorylated_residues_are_read_with_their_defined_atom_counts(name, n_atoms):
    result = analyze(method="aa_ind", chain="A", residue_types=name)
    assert result.counts.n_points == n_atoms
    assert result.warnings == ()
    assert side(result).get("unsupported") is None


def test_tpo_is_a_complete_residue_of_the_per_residue_method_and_the_sweep():
    one = analyze(method="amino_acid_atomic", chain="A", residue="3", residue_name="TPO")
    assert one.counts.n_points == 11 and one.counts.n_triangles == 165
    swept = biotsr.analyze_entities(SYN, chain="A", keep_triplets=False)
    assert {"3_TPO", "5_SEP", "7_PTR"} <= set(swept)
    assert not any(key.endswith("_CME") for key in swept)


def test_a_mixed_selection_names_the_unsupported_residue_and_is_not_complete():
    result = analyze(method="aa_ind_mix", chain="A", residue_types="TYR;TPO;CME")
    assert set(result.selection.residues_used) == {"3_TPO", "6_TYR"}
    warning = [w for w in result.warnings if w.startswith("not supported:")]
    assert warning and "residue type CME" in warning[0]
    assert "reads only TPO, SEP, PTR" in warning[0]
    assert side(result)["unsupported"] == ["CME"]


def test_a_request_of_unsupported_residues_only_is_refused_with_the_reason():
    with pytest.raises(InvalidRequest, match="not supported: residue type CME"):
        analyze(method="aa_ind", chain="A", residue_types="CME")
    with pytest.raises(InvalidRequest, match="residue 8 CME .* is not supported"):
        analyze(method="amino_acid_atomic", chain="A", residue="8", residue_name="CME")
    with pytest.raises(InvalidRequest, match="CME_8 .*not supported"):
        biotsr.analyze_entities(SYN, chain="A", residue_name="CME", keep_triplets=False)
    # A ligand is not an amino acid, modified or not: simply absent.
    with pytest.raises(InvalidRequest, match="no residue of type LIG"):
        analyze(method="aa_ind", chain="A", residue_types="LIG")


def test_named_residues_and_atoms_separate_not_found_from_not_supported():
    result = analyze(method="aa_ind_specific", chain="A", residues="TPO_3;CME_8;LYS_99")
    record = side(result)
    assert record["instances"] == ["TPO_3"]
    assert record["not_found"] == ["LYS_99"] and record["unsupported"] == ["CME_8"]
    assert "requested but not found: LYS_99" in result.warnings
    atoms = analyze(method="amino_acid_atom_selection", chain="A",
                    atoms="TPO_3_P;PTR_7_P;CME_8_SG;TYR_6_OH")
    assert atoms.counts.n_points == 3
    assert side(atoms)["unsupported"] == ["CME_8_SG"]
    assert not any(w.startswith("requested but not found") for w in atoms.warnings)


def test_the_calpha_methods_read_the_supported_modified_residues():
    """ADR-107 (until 1.0.0rc5 every one of them was named as not represented)."""
    chain = analyze(method="protein", chain="A")
    assert {"3", "5", "7", "8"} <= set(chain.selection.residues_used)
    assert "unsupported" not in side(chain)
    assert not any(w.startswith("not represented") for w in chain.warnings)
    listed = analyze(method="aa_selection", chain="A", residues="1-4")
    assert side(listed)["instances"] == ["1", "2", "3", "4"]
    assert "unsupported" not in side(listed)
    assert side(listed)["excluded"] == []


def test_the_atomic_cross_sides_read_tpo_and_name_cme():
    result = biotsr.analyze_cross(SYN, method="protein_protein_cross", chain="A",
                                  chain_2="B", keep_triplets=False)
    used = set(result.selection.residues_used)
    assert {"3_TPO", "5_SEP", "7_PTR"} <= used and "8_CME" not in used
    assert side(result, 0)["unsupported"] == ["CME_8"]
    assert "unsupported" not in side(result, 1)
    calpha = biotsr.analyze_cross(SYN, method="protein_protein_calpha_cross",
                                  chain="A", chain_2="B", keep_triplets=False)
    assert "unsupported" not in side(calpha, 0)
    assert "unsupported" not in side(calpha, 1)


def test_the_alternate_location_rule_applies_to_modified_residues(tmp_path):
    """Blank and A are kept, B is excluded, and TPO's exact atom count then
    declines a TPO that lost an atom to it -- recorded, never relaxed."""
    lines = open(SYN).read().splitlines(True)
    edited = []
    for line in lines:
        if line.startswith("HETATM") and line[17:20] == "TPO" and line[21] == "A" \
                and line[12:16].strip() == "O3P":
            line = line[:16] + "B" + line[17:]
        edited.append(line)
    path = write(tmp_path / "ALTB.pdb", "".join(edited))
    result = biotsr.analyze(path, method="aa_ind_mix", chain="A",
                            residue_types="TPO;TYR", keep_triplets=False)
    assert "3_TPO" not in result.selection.residues_used
    assert any(s.reason == "altloc_rejected" for s in result.skipped)
    assert any(s.stage == "completeness" and "TPO" in s.detail for s in result.skipped)
    named = biotsr.analyze(path, method="aa_ind_specific", chain="A",
                           residues="TPO_3;TYR_6", keep_triplets=False)
    assert any(w.startswith("requested but excluded by the method's completeness rule: "
                            "TPO_3") for w in named.warnings)


def test_a_batch_reports_a_calculation_that_lacks_a_named_item_as_partial(tmp_path):
    pdb = tmp_path / "pdb"
    pdb.mkdir()
    shutil.copy(SYN, str(pdb))
    manifest = write(tmp_path / "m.csv",
                     "pdb_id,chain,method,residues,sample_id\n"
                     "SYNPHOS,A,aa_ind_specific,TPO_3;TYR_6,complete\n"
                     "SYNPHOS,A,aa_ind_specific,TPO_3;CME_8,partly\n")
    code, out, err = run(["batch", "--manifest", manifest, "--pdb-dir", pdb,
                          "--output", tmp_path / "out", "--no-triplets"])
    assert code == EXIT_INCOMPLETE, err  # partial results: exit 4 (ADR-106)
    rows = {r["sample_id"]: r for r in summary(tmp_path / "out")}
    assert rows["complete"]["status"] == "ok"
    assert rows["partly"]["status"] == "partial"
    assert "not supported: CME_8" in rows["partly"]["message"]
    assert "1 partial" in err


def test_a_per_residue_sweep_names_the_modified_residues_it_leaves_out(tmp_path):
    assert biotsr.api.unsupported_residues(SYN, chain="A",
                                       method="amino_acid_atomic") == ("CME_8",)
    assert biotsr.api.unsupported_residues(SYN, chain="A", method="protein") == ()
    code, _out, err = run(["run", "--method", "amino_acid_atomic", "--pdb", SYN,
                           "--chain", "A", "--output", tmp_path / "sweep",
                           "--no-triplets"])
    assert code == EXIT_INCOMPLETE, err  # partial results: exit 4 (ADR-106)
    assert "not supported: CME_8" in err and "not part of this per-residue sweep" in err
    record = json.load(open(str(tmp_path / "sweep" /
                                "SYNPHOS_A_amino_acid_atomic_perresidue.run.json")))
    assert [r["name"] for r in record["results"]][:3] == [
        "SYNPHOS_A_amino_acid_atomic_perresidue_GLY_1",
        "SYNPHOS_A_amino_acid_atomic_perresidue_ALA_2",
        "SYNPHOS_A_amino_acid_atomic_perresidue_TPO_3"]
    assert record["warnings"][0].startswith("not supported: CME_8")


# ===========================================================================
# Residue identifiers and ranges
# ===========================================================================


@pytest.mark.parametrize("residues, message", [
    ("9-11", "span residue\\(s\\) 10A, 10B"),
    ("10A-11", "'10A-11' is not supported"),
    ("11-9", "'11-9' is descending"),
    ("9-", "'9-' is not supported"),
    ("ten", "'ten' is not a PDB residue identifier"),
])
def test_the_range_grammar_refuses_what_it_does_not_define(residues, message):
    with pytest.raises(InvalidRequest, match=message):
        analyze(method="aa_selection", chain="A", residues=residues)


def test_listing_the_insertion_code_residues_includes_them_and_order_never_matters():
    listed = analyze(method="aa_selection", chain="A", residues="9-11;10A;10B")
    shuffled = analyze(method="aa_selection", chain="A", residues="10B;11;10A;9-10")
    one_by_one = analyze(method="aa_selection", chain="A", residues="9;10;11")
    assert listed.selection.residues_used == ("9", "10", "10A", "10B", "11")
    assert shuffled.keys == listed.keys
    assert one_by_one.selection.residues_used == ("9", "10", "11")
    assert analyze(method="aa_selection", chain="A",
                   residues="009;10").selection.residues_used == ("9", "10")


def test_insertion_codes_in_every_kind_of_identifier():
    assert analyze(method="aa_ind_specific", chain="A",
                   residues="GLY_10A;ALA10B").selection.residues_used == (
        "10A_GLY", "10B_ALA")
    one = analyze(method="amino_acid_atomic", chain="A", residue="10A",
                  residue_name="GLY")
    assert one.counts.n_points == 4
    atoms = analyze(method="amino_acid_atom_selection", chain="A",
                    atoms="GLY_10A_CA;ALA_10B_CB;LYS_10_NZ")
    assert atoms.counts.n_points == 3


def test_the_api_the_command_line_and_a_manifest_agree(tmp_path):
    api = analyze(method="aa_selection", chain="A", residues="9-11;10A;10B")
    code, _out, err = run(["run", "--method", "aa_selection", "--pdb", SYN,
                           "--chain", "A", "--residues", "9-11;10A;10B",
                           "--output", tmp_path / "cli", "--no-triplets"])
    assert code == EXIT_OK, err
    cli = biotsr.load_result(str(tmp_path / "cli" /
                                 "SYNPHOS_A_aa_selection_res-9-11+10A+10B.keys"))
    assert dict(cli.keys) == dict(api.keys)
    pdb = tmp_path / "pdb"
    pdb.mkdir()
    shutil.copy(SYN, str(pdb))
    manifest = write(tmp_path / "m.csv", "pdb_id,chain,method,residues\n"
                                         'SYNPHOS,A,aa_selection,"9-11;10A;10B"\n'
                                         "SYNPHOS,A,aa_selection,9-11\n")
    code, _out, err = run(["batch", "--manifest", manifest, "--pdb-dir", pdb,
                           "--output", tmp_path / "csv", "--no-triplets"])
    rows = summary(tmp_path / "csv")
    assert rows[0]["status"] == "ok"
    assert rows[1]["status"] == "failed" and "span residue(s) 10A, 10B" in rows[1]["message"]
    batch = biotsr.load_result(str(tmp_path / "csv" /
                                   "SYNPHOS_A_aa_selection_res-9-11+10A+10B.keys"))
    assert dict(batch.keys) == dict(api.keys)
    code, _out, err = run(["run", "--method", "aa_selection", "--pdb", SYN,
                           "--chain", "A", "--residues", "11-9",
                           "--output", tmp_path / "bad", "--no-triplets"])
    assert code != EXIT_OK and "descending" in err


def test_the_calpha_cross_applies_the_same_range_grammar():
    with pytest.raises(InvalidRequest, match="span residue\\(s\\) 10A, 10B"):
        biotsr.analyze_cross(SYN, method="protein_protein_calpha_cross", chain="A",
                             chain_2="B", calpha_residues="8-11", keep_triplets=False)


# ===========================================================================
# Duplicate rows: computations versus samples
# ===========================================================================


def _batch(tmp_path, manifest_text, *extra, output="out"):
    manifest = write(tmp_path / ("%s.csv" % output), manifest_text)
    return run(["batch", "--manifest", manifest, "--pdb-dir", STRUCTURES,
                "--output", tmp_path / output, "--no-triplets"] + list(extra))


def _samples(tmp_path, output="out", *extra):
    code, _out, err = run(["analyze", "common-keys", "--from-results", tmp_path / output,
                           "--output", tmp_path / (output + "_an")] + list(extra))
    assert code == EXIT_OK, err
    report = json.load(open(str(tmp_path / (output + "_an") / "analysis" /
                                "common_keys" / "all_samples.json")))
    return [r["fileName"] for r in report["rows"]]


def test_distinct_sample_ids_are_distinct_samples_of_one_calculation(tmp_path):
    code, _out, err = _batch(tmp_path,
                             "pdb_id,chain,method,residue_types,sample_id,label\n"
                             "1UBQ,A,aa_ind,LYS,ubq_first,group1\n"
                             "1UBQ,A,aa_ind,LYS,ubq_second,group2\n")
    assert code == EXIT_OK, err
    rows = summary(tmp_path / "out")
    assert [(r["rows"], r["sample_id"], r["status"]) for r in rows] == [
        ("2", "ubq_first", "ok"), ("3", "ubq_second", "ok")]
    assert _samples(tmp_path) == ["ubq_first", "ubq_second"]
    first = biotsr.load_result(str(tmp_path / "out" / "ubq_first.keys"))
    second = biotsr.load_result(str(tmp_path / "out" / "ubq_second.keys"))
    assert dict(first.keys) == dict(second.keys)
    records = [json.load(open(str(tmp_path / "out" / ("%s.run.json" % name))))
               for name in ("ubq_first", "ubq_second")]
    assert [r["label"] for r in records] == ["group1", "group2"]


def test_one_calculation_of_one_sample_is_computed_once_with_every_row(tmp_path):
    code, _out, err = _batch(tmp_path,
                             "pdb_id,chain,method,residue_types,label,note\n"
                             "1UBQ,A,aa_ind,LYS,group1,first entry\n"
                             "1UBQ,A,aa_ind,LYS,group1,second entry\n"
                             "1UBQ,A,aa_ind,TYR,group1,another sample\n",
                             "--keep-column", "note")
    assert code == EXIT_OK, err
    rows = summary(tmp_path / "out")
    assert [r["rows"] for r in rows] == ["2;3", "4"]
    assert "one sample from 2 source requests (line 2, 3)" in rows[0]["message"]
    assert rows[0]["note"] == "line 2=first entry; line 3=second entry"
    assert _samples(tmp_path) == ["1UBQ_A_aa_ind_LYS", "1UBQ_A_aa_ind_TYR"]


def test_rows_that_would_merge_with_different_labels_are_refused(tmp_path):
    code, _out, err = _batch(tmp_path,
                             "pdb_id,chain,method,residue_types,label\n"
                             "1UBQ,A,aa_ind,LYS,group1\n"
                             "1UBQ,A,aa_ind,LYS,group2\n")
    assert code == EXIT_USAGE
    assert "requested by rows with different labels" in err
    assert "label 'group1' on line 2; label 'group2' on line 3" in err
    assert not (tmp_path / "out").exists()


def test_the_collision_check_still_refuses_two_requests_with_one_sample_id(tmp_path):
    code, _out, err = _batch(tmp_path,
                             "pdb_id,chain,method,residue_types,sample_id\n"
                             "1UBQ,A,aa_ind,LYS,same\n"
                             "1UBQ,A,aa_ind,TYR,same\n")
    assert code == EXIT_USAGE and "more than one different request" in err


def test_two_participants_asking_for_one_calculation_are_one_sample(tmp_path):
    sheet = write(tmp_path / "sheet.csv",
                  "protein,aa_chain,aa_chain_2,aa_group,aa_ind,aa_ind_2\n"
                  "1UBQ,A,A,ubiquitin,LYS,LYS\n")
    code, _out, err = run(["batch", "--manifest", sheet, "--sheet", "tsr",
                           "--method", "aa_ind", "--participant", "all",
                           "--pdb-dir", STRUCTURES, "--output", tmp_path / "out",
                           "--no-triplets"])
    assert code == EXIT_OK, err
    rows = summary(tmp_path / "out")
    assert len(rows) == 1 and rows[0]["participant"] == "1;2"
    assert "one sample from 2 source requests" in rows[0]["message"]


# ===========================================================================
# Cross sides from a wide sheet
# ===========================================================================

WIDE = ("protein,aa_chain,drug_chain,carb_chain,aa_chain_2,aa_group,aa_ind,aa_ind_mix,"
        "aa_ind_specific,aa_ind_2,aa_ind_specific_2,drug_name,drug_id,carb_name,"
        "carb_id,aa_selection,aa_selection_2,aa_atom_selection,drug_atom_selected,"
        "car_atom_selected\n"
        "6G9M,A,A,,,kinase,TYR,TYR;ARG,ARG_15;GLU12,,,ESW,404,,,12-15,,TYR_30_OH,"
        "C3;C5;O6,\n"
        "3PTB,A,A,,,protease,,,,,,BEN,1,,,,,,,\n"
        "7LF9,B,,D,,glycoprotein,TYR,,,,,,,NAG,1,,,,,C1;C2;N2\n"
        "1BRS,A,,,D,barnase,LYS,,LYS_19;TYR24,ASP,LYS1;ARG_11,,,,,1-30,1-20,,,\n")


def _sheet(tmp_path, *extra, output="out"):
    sheet = write(tmp_path / "wide.csv", WIDE)
    return run(["batch", "--manifest", sheet, "--sheet", "tsr", "--pdb-dir", STRUCTURES,
                "--output", tmp_path / output, "--no-triplets",
                "--inapplicable", "skip"] + list(extra))


def test_a_populated_narrowing_column_needs_an_explicit_side_choice(tmp_path):
    code, _out, err = _sheet(tmp_path, "--method", "drug_protein_cross", "--dry-run")
    assert code == EXIT_USAGE
    assert "does not choose among them" in err
    assert "the drug side of drug_protein_cross: drug_atom_selected on line 2" in err
    assert "the protein side of drug_protein_cross: aa_atom_selection on line 2; " \
           "aa_ind on line 2; aa_ind_mix on line 2; aa_ind_specific on line 2" in err
    assert "--cross-side drug_protein_cross.drug=whole" in err


def test_a_chosen_column_narrows_the_side_and_the_others_are_named(tmp_path):
    code, out, err = _sheet(tmp_path, "--method", "drug_protein_cross",
                            "--cross-side", "drug_protein_cross.protein=aa_ind",
                            "--cross-side", "drug_protein_cross.drug=drug_atom_selected",
                            "--dry-run")
    assert code == EXIT_OK, err
    lines = out.splitlines()
    first = [i for i, l in enumerate(lines) if "6G9M_A_A_drug_protein_cross" in l][0]
    assert lines[first + 1].strip() == "drug, chain A: ESW copy 404, atoms=C3;C5;O6"
    assert lines[first + 2].strip() == "protein, chain A: residue_types=TYR"
    assert "not applied: aa_ind_mix (the protein side uses aa_ind)" in lines[first + 3]
    # Row 3 names a drug but leaves both chosen columns blank: it does not
    # describe the narrowed calculation, and is never computed whole.
    assert any("inapplicable: the row leaves aa_ind blank" in l or
               "inapplicable: the row leaves drug_atom_selected blank" in l
               for l in lines if l.split()[1:2] == ["3"])


def test_whole_sides_are_chosen_and_said(tmp_path):
    code, out, err = _sheet(tmp_path, "--method", "carbohydrate_protein_cross",
                            "--cross-side", "carbohydrate_protein_cross.carbohydrate=whole",
                            "--cross-side", "carbohydrate_protein_cross.protein=whole",
                            "--dry-run")
    assert code == EXIT_OK, err
    assert "carbohydrate, chain D: NAG copy 1, all atoms" in out
    assert "protein, chain B: whole chain" in out
    assert "not applied: car_atom_selected (the carbohydrate side is used whole)" in out


def test_protein_protein_sides_read_their_own_participant_columns(tmp_path):
    code, out, err = _sheet(tmp_path, "--method", "protein_protein_cross",
                            "--cross-side", "protein_protein_cross.1=aa_ind_specific",
                            "--cross-side", "protein_protein_cross.2=aa_ind_specific_2",
                            "--dry-run")
    assert code == EXIT_OK, err
    assert "side 1, chain A: residues=LYS_19;TYR24" in out
    assert "side 2, chain D: residues=LYS1;ARG_11" in out
    code, out, err = _sheet(tmp_path, "--method", "protein_protein_calpha_cross",
                            "--cross-side", "protein_protein_calpha_cross.1=aa_selection",
                            "--cross-side", "protein_protein_calpha_cross.2=whole",
                            output="ca")
    assert code == EXIT_INCOMPLETE, err  # partial results: exit 4 (ADR-106)
    row = [r for r in summary(tmp_path / "ca") if r["rows"] == "5"][0]
    assert row["selection"] == ("side 1, chain A: calpha_residues=1-30 | "
                                "side 2, chain D: whole chain")
    # Barnase chain A starts at residue 3: 1-30 asks for two residues it
    # does not have, so the calculation is partial, and says which.
    assert row["status"] == "partial"
    assert "requested but not found: 1, 2 (protein side, chain A)" in row["message"]
    assert "not applied: aa_selection_2 (side 2 is used whole)" in row["message"]


@pytest.mark.parametrize("extra, message", [
    (["--cross-side", "drug_protein_cross.protein=aa_selection"],
     "does not narrow the protein side"),
    (["--cross-side", "drug_protein_cross.ligand=whole"], "has no side 'ligand'"),
    (["--cross-side", "drug.protein=whole"], "not a cross method"),
    (["--cross-side", "drug_protein_cross.protein=aa_ind",
      "--cross-side", "drug_protein_cross.protein=whole"], "twice"),
    (["--cross-side", "drug_protein_cross.protein=whole",
      "--map", "drug_protein_cross.residue_types=aa_ind"], "chosen twice"),
    (["--cross-side", "carbohydrate_protein_cross.protein=whole"],
     "not among the --method choices"),
])
def test_unclear_side_choices_are_refused(tmp_path, extra, message):
    code, _out, err = _sheet(tmp_path, "--method", "drug_protein_cross", "--dry-run",
                             *extra)
    assert code == EXIT_USAGE and message in err, err


def test_cross_side_outside_a_sheet_is_refused(tmp_path):
    manifest = write(tmp_path / "m.csv", "pdb_id,chain,method\n1UBQ,A,protein\n")
    code, _out, err = run(["batch", "--manifest", manifest, "--pdb-dir", STRUCTURES,
                           "--output", tmp_path / "out", "--cross-side",
                           "drug_protein_cross.protein=whole"])
    assert code == EXIT_USAGE and "--cross-side chooses among a wide sheet" in err


@pytest.mark.parametrize("participant, expected", [
    ("1", {("aa_ind", "A")}),
    ("2", {("aa_ind", "D")}),
    ("all", {("aa_ind", "A"), ("aa_ind", "D")}),
])
def test_participants_choose_the_chain_and_selection_columns(tmp_path, participant,
                                                             expected):
    code, out, err = _sheet(tmp_path, "--method", "aa_ind", "--participant", participant,
                            "--dry-run")
    assert code == EXIT_OK, err
    planned = {(line.split()[2], line.split()[4]) for line in out.splitlines()
               if "1BRS" in line and "[would run]" in line}
    assert planned == expected
    if participant != "1":
        assert "chain D: residue_types=ASP" in out


def test_the_python_api_takes_the_same_side_choices(tmp_path):
    sheet = write(tmp_path / "wide.csv", WIDE)
    report = biotsr.batch.run(
        sheet, str(tmp_path / "out"), pdb_dir=STRUCTURES, sheet="tsr",
        methods=["drug_protein_cross"], inapplicable="skip", dry_run=True,
        cross_sides={"drug_protein_cross.protein": "whole",
                     "drug_protein_cross.drug": "whole"})
    planned = [r for r in report.records if r["status"] == "would_run"]
    assert [r["selection"] for r in planned] == [
        "drug, chain A: ESW copy 404, all atoms | protein, chain A: whole chain",
        "drug, chain A: BEN copy 1, all atoms | protein, chain A: whole chain"]
