"""Wide dataset sheets, methods chosen on the command line, one shared output folder.

``tests/data/manifests/wide_sheet.csv`` is synthetic (its README says how it was
made): one row per structure, several participant chains, molecules with and
without identifiers, the second protein chain's columns, and both residue
spellings.  These tests run it through ``biotsr batch --sheet`` (ADR-097) and
check the output folder contract (ADR-098).
"""

from __future__ import annotations

import csv
import io
import json
import os
import shutil

import pytest

from biotsr.cli import EXIT_OK, EXIT_PARTIAL, EXIT_USAGE, main

from ..conftest import STRUCTURES, TESTS_DIR

SHEET = os.path.join(TESTS_DIR, "data", "manifests", "wide_sheet.csv")


def run(argv):
    import contextlib

    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main([str(a) for a in argv])
    return code, out.getvalue(), err.getvalue()


def batch(tmp_path, *extra, sheet=SHEET, output="out"):
    return run(["batch", "--manifest", sheet, "--pdb-dir", STRUCTURES,
                "--output", tmp_path / output, "--no-triplets"] + list(extra))


def summary(tmp_path, output="out"):
    with open(os.path.join(str(tmp_path / output), "batch_summary.csv"), newline="") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------------------
# Planning: methods outside the sheet, method-aware columns
# ---------------------------------------------------------------------------


def test_the_preview_lists_every_request_with_its_source_row_and_prefix(tmp_path):
    code, out, err = batch(tmp_path, "--sheet", "tsr", "--method", "carbohydrate",
                           "--method", "aa_ind", "--inapplicable", "skip", "--dry-run")
    assert code == EXIT_OK, err
    assert not (tmp_path / "out").exists()            # a dry run writes nothing
    lines = out.splitlines()
    header = [l for l in lines if "prefix / outcome" in l][0]
    assert "rows" in header and "method" in header and "chains" in header
    assert any("7LF9_D_carbohydrate_NAG_5" in l and "[would run]" in l for l in lines)
    assert any("7LF9_D_carbohydrate_MAN_allcopies" in l for l in lines)
    # The protein side of three glycan rows is one calculation, serving all three.
    assert any("4;5;6" in l and "7LF9_B_aa_ind_TYR" in l for l in lines)
    assert any("inapplicable: the row has no carb_name / carb_id" in l for l in lines)
    unused = [l for l in lines if l.startswith("unused")][0]
    for column in ("drug_name", "aa_ind_2", "aa_selection", "car_atom_selected"):
        assert column in unused
    assert "comment_1" not in unused


def test_several_methods_expand_one_row_and_every_output_names_its_source(tmp_path):
    code, out, err = batch(tmp_path, "--sheet", "tsr", "--method", "drug",
                           "--method", "aa_ind_specific", "--inapplicable", "skip")
    assert code == EXIT_OK, err
    rows = summary(tmp_path)
    by_status = {}
    for r in rows:
        by_status.setdefault(r["status"], []).append(r)
    ok = {(r["rows"], r["method"]): r["sample_id"] for r in by_status["ok"]}
    assert ok[("2", "drug")] == "6G9M_A_drug_ESW_404"
    assert ok[("2", "aa_ind_specific")] == "6G9M_A_aa_ind_specific_ARG15+GLU12"
    assert ok[("3", "drug")] == "3PTB_A_drug_BEN_1"
    assert ok[("8", "aa_ind_specific")] == "1BRS_A_aa_ind_specific_LYS19+TYR24"
    assert {r["status"] for r in rows} == {"ok", "inapplicable"}
    record = json.load(open(str(tmp_path / "out" / "6G9M_A_drug_ESW_404.run.json")))
    assert record["source"] == {"manifest": "wide_sheet.csv", "rows": [2],
                                "participant": "1",
                                "columns": {"chain": "drug_chain",
                                            "entity_name": "drug_name",
                                            "entity_id": "drug_id"}}
    assert record["complete"] is True


def test_an_inapplicable_row_fails_unless_skipping_is_asked_for(tmp_path):
    code, _out, err = batch(tmp_path, "--sheet", "tsr", "--method", "nucleotide")
    assert code == EXIT_PARTIAL
    rows = summary(tmp_path)
    failed = [r for r in rows if r["status"] == "failed"]
    assert len(failed) == 7 and len([r for r in rows if r["status"] == "ok"]) == 1
    assert "the row has no nt_name / nt_id for nucleotide" in failed[0]["message"]
    assert "--inapplicable skip" in failed[0]["message"]


def test_the_second_participant_runs_only_when_asked_for(tmp_path):
    code, _out, err = batch(tmp_path, "--sheet", "tsr", "--method", "aa_ind",
                            "--participant", "2", "--inapplicable", "skip")
    assert code == EXIT_OK, err
    ok = [r for r in summary(tmp_path) if r["status"] == "ok"]
    assert [(r["rows"], r["participant"], r["sample_id"]) for r in ok] == [
        ("8", "2", "1BRS_D_aa_ind_ASP")]


def test_a_cross_method_resolves_both_sides_from_their_own_columns(tmp_path):
    """Both sides of the glycan rows are chosen whole explicitly: the rows also
    populate car_atom_selected and aa_ind, which could narrow them."""
    code, _out, err = batch(tmp_path, "--sheet", "tsr",
                            "--method", "carbohydrate_protein_cross",
                            "--method", "protein_protein_cross",
                            "--cross-side", "carbohydrate_protein_cross.carbohydrate=whole",
                            "--cross-side", "carbohydrate_protein_cross.protein=whole",
                            "--map", "protein_protein_cross.residue_types=aa_ind",
                            "--map", "protein_protein_cross.residue_types_2=aa_ind_2",
                            "--inapplicable", "skip")
    assert code == EXIT_OK, err
    ok = {r["sample_id"]: r for r in summary(tmp_path) if r["status"] == "ok"}
    assert "7LF9_D_B_carbohydrate_protein_cross_NAG_1" in ok
    assert "7LF9_D_B_carbohydrate_protein_cross_MAN_allcopies" in ok
    pp = ok["1BRS_A_D_protein_protein_cross_A-LYS_D-ASP"]
    assert (pp["chain"], pp["chain_2"]) == ("A", "D")


def test_one_sheet_two_methods_two_different_chain_columns(tmp_path):
    """carbohydrate reads carb_chain (D), aa_ind reads aa_chain (B), from the same rows."""
    code, _out, err = batch(tmp_path, "--sheet", "tsr", "--method", "carbohydrate",
                            "--method", "aa_ind", "--inapplicable", "skip")
    assert code == EXIT_OK, err
    rows = [r for r in summary(tmp_path) if r["status"] == "ok"]
    chains = {(r["method"], r["rows"]): r["chain"] for r in rows}
    assert chains[("carbohydrate", "4")] == "D"
    assert chains[("aa_ind", "4;5;6;7")] == "B"


@pytest.mark.parametrize("args, message", [
    (["--map", "nonsense.chain=aa_chain"], "unknown method"),
    (["--map", "carbohydrate.colour=aa_chain"], "does not read 'colour'"),
    (["--map", "carbohydrate.chain=no_such_column"], "no column 'no_such_column'"),
    (["--map", "carbohydrate.chain=aa_chain", "--map", "carbohydrate.chain=nt_chain"],
     "both map carbohydrate"),
    (["--map", "carbohydrate.3.chain=aa_chain"], "participant '3'"),
    (["--map", "carbohydrate_protein_cross.2.chain=aa_chain"], "one participant"),
    (["--map", "carbohydrate"], "METHOD[.PARTICIPANT].FIELD=HEADER"),
])
def test_unknown_or_conflicting_mappings_are_refused_with_the_remedy(tmp_path, args, message):
    code, _out, err = batch(tmp_path, "--sheet", "tsr", "--method", "carbohydrate",
                            "--dry-run", *args)
    assert code == EXIT_USAGE
    assert message in err


def test_a_custom_layout_reads_arbitrary_headers(tmp_path):
    sheet = tmp_path / "mine.csv"
    sheet.write_text("entry,glycan_chain,sugar,sugar_id,notes\n"
                     "7LF9,D,NAG,5,first\n7LF9,D,MAN,,second\n")
    code, _out, err = batch(tmp_path, "--sheet", "custom", "--method", "carbohydrate",
                            "--column", "pdb_id=entry",
                            "--map", "carbohydrate.chain=glycan_chain",
                            "--map", "carbohydrate.entity_name=sugar",
                            "--map", "carbohydrate.entity_id=sugar_id",
                            "--keep-column", "notes", sheet=str(sheet))
    assert code == EXIT_OK, err
    rows = summary(tmp_path)
    assert [(r["sample_id"], r["notes"]) for r in rows] == [
        ("7LF9_D_carbohydrate_NAG_5", "first"),
        ("7LF9_D_carbohydrate_MAN_allcopies", "second")]


def test_a_map_file_is_the_same_as_map_options(tmp_path):
    maps = tmp_path / "map.csv"
    maps.write_text("method,participant,field,column\n"
                    "carbohydrate,,chain,glycan_chain\n"
                    "carbohydrate,1,entity_name,sugar\n")
    sheet = tmp_path / "mine.csv"
    sheet.write_text("pdb_id,glycan_chain,sugar\n7LF9,D,BMA\n")
    code, out, err = batch(tmp_path, "--sheet", "custom", "--method", "carbohydrate",
                           "--map-file", maps, "--dry-run", sheet=str(sheet))
    assert code == EXIT_OK, err
    assert "7LF9_D_carbohydrate_BMA_allcopies" in out


def test_both_spellings_of_the_carbohydrate_atom_column_are_refused_together(tmp_path):
    sheet = tmp_path / "both.csv"
    sheet.write_text("protein,carb_chain,carb_name,carb_id,car_atom_selected,"
                     "carb_atom_selected\n7LF9,D,NAG,1,C1;C2,C1;C2\n")
    code, _out, err = batch(tmp_path, "--sheet", "tsr",
                            "--method", "carbohydrate_atom_selection", "--dry-run",
                            sheet=str(sheet))
    assert code == EXIT_USAGE and "two columns for 'carb_atom_selected'" in err


def test_the_car_atom_selected_alias_feeds_the_atom_selection_method(tmp_path):
    code, _out, err = batch(tmp_path, "--sheet", "tsr",
                            "--method", "carbohydrate_atom_selection",
                            "--inapplicable", "skip")
    assert code == EXIT_OK, err
    ok = [r["sample_id"] for r in summary(tmp_path) if r["status"] == "ok"]
    assert ok == ["7LF9_D_carbohydrate_atom_selection_NAG_%s_atoms-C1+C2+N2" % i
                  for i in (1, 2, 5)]


def test_a_sheet_takes_no_method_column_and_no_default_chain(tmp_path):
    sheet = tmp_path / "with_method.csv"
    sheet.write_text("protein,aa_chain,method\n1UBQ,A,protein\n")
    code, _out, err = batch(tmp_path, "--sheet", "tsr", "--method", "protein",
                            "--dry-run", sheet=str(sheet))
    assert code == EXIT_USAGE and "'method' column" in err
    code, _out, err = batch(tmp_path, "--sheet", "tsr", "--method", "protein",
                            "--chain", "A", "--dry-run")
    assert code == EXIT_USAGE and "--chain does not apply to a sheet" in err
    code, _out, err = batch(tmp_path, "--sheet", "tsr", "--dry-run")
    assert code == EXIT_USAGE and "--method" in err


def test_a_blank_participant_chain_is_reported_never_inferred(tmp_path):
    sheet = tmp_path / "no_chain2.csv"
    sheet.write_text("protein,aa_chain,aa_chain_2,aa_ind,aa_ind_2\n1BRS,A,,LYS,ASP\n")
    code, _out, err = batch(tmp_path, "--sheet", "tsr", "--method", "aa_ind",
                            "--participant", "all", sheet=str(sheet))
    assert code == 4   # ADR-112: not a failure, and not a complete batch
    rows = summary(tmp_path)
    assert [(r["participant"], r["status"]) for r in rows] == [
        ("1", "ok"), ("2", "chain_not_specified")]
    assert "aa_chain_2 blank" in rows[1]["message"]
    assert "never assumes a chain" in rows[1]["message"]


# ---------------------------------------------------------------------------
# One-request-per-row manifests and command-line methods
# ---------------------------------------------------------------------------


def test_a_manifest_method_column_and_method_option_conflict(tmp_path):
    manifest = tmp_path / "m.csv"
    manifest.write_text("pdb_id,method,chain\n1UBQ,protein,A\n")
    code, _out, err = batch(tmp_path, "--method", "aa_selection", "--dry-run",
                            sheet=str(manifest))
    assert code == EXIT_USAGE
    assert "--method also chooses methods" in err and "--ignore-column method" in err


def test_methods_from_the_command_line_run_on_every_manifest_row(tmp_path):
    manifest = tmp_path / "m.csv"
    manifest.write_text("pdb_id,chain,residue_types\n1UBQ,A,LYS\n1UBQ,A,TYR\n")
    code, _out, err = batch(tmp_path, "--method", "aa_ind", "--method", "protein",
                            sheet=str(manifest))
    # protein takes no residue_types, so its two requests fail by name
    assert code == EXIT_PARTIAL
    rows = summary(tmp_path)
    assert [(r["rows"], r["method"], r["status"]) for r in rows] == [
        ("2", "aa_ind", "ok"), ("2", "protein", "failed"),
        ("3", "aa_ind", "ok"), ("3", "protein", "failed")]
    assert "does not use --residue-types" in rows[1]["message"]


# ---------------------------------------------------------------------------
# The shared output folder
# ---------------------------------------------------------------------------


def _folder(tmp_path, name):
    root = tmp_path / name
    return sorted(n for n in os.listdir(str(root))
                  if not n.startswith(".") and n != "batch_summary.csv")


def test_one_shared_folder_holds_every_calculation_side_by_side(tmp_path):
    code, _out, err = batch(tmp_path, "--sheet", "tsr", "--method", "carbohydrate",
                            "--inapplicable", "skip")
    assert code == EXIT_OK, err
    names = _folder(tmp_path, "out")
    assert "7LF9_D_carbohydrate_NAG_1.keys" in names
    assert "7LF9_D_carbohydrate_NAG_1.meta.json" in names
    assert "7LF9_D_carbohydrate_NAG_1.run.json" in names
    assert not any(n.endswith(".triplets") for n in names)          # keys only
    assert not any(os.path.isdir(str(tmp_path / "out" / n)) for n in names)
    assert not [n for n in os.listdir(str(tmp_path / "out")) if n.startswith(".biotsr-")]


def test_serial_and_parallel_runs_write_identical_calculations(tmp_path):
    args = ("--sheet", "tsr", "--method", "carbohydrate", "--method", "aa_ind",
            "--inapplicable", "skip")
    assert batch(tmp_path, *args, output="serial")[0] == EXIT_OK
    code, _out, err = batch(tmp_path, *args, "--jobs", "3", output="parallel")
    assert code == EXIT_OK, err
    assert _folder(tmp_path, "serial") == _folder(tmp_path, "parallel")
    for name in _folder(tmp_path, "serial"):
        a = open(str(tmp_path / "serial" / name), "rb").read()
        b = open(str(tmp_path / "parallel" / name), "rb").read()
        if name.endswith(".run.json"):
            ja, jb = json.loads(a), json.loads(b)
            for volatile in ("argv",):
                ja.pop(volatile), jb.pop(volatile)
            assert ja == jb, name
        else:
            assert a == b, name


def test_resume_keeps_verified_calculations_and_recomputes_missing_ones(tmp_path):
    args = ("--sheet", "tsr", "--method", "carbohydrate", "--inapplicable", "skip")
    assert batch(tmp_path, *args)[0] == EXIT_OK
    root = tmp_path / "out"
    os.remove(str(root / "7LF9_D_carbohydrate_NAG_2.keys"))
    record = json.load(open(str(root / "7LF9_D_carbohydrate_NAG_5.run.json")))
    record["complete"] = False                      # as an interrupted install leaves it
    (root / "7LF9_D_carbohydrate_NAG_5.run.json").write_text(json.dumps(record))
    code, _out, err = batch(tmp_path, *args)
    assert code == EXIT_PARTIAL                     # --if-exists fail: refuses both
    code, _out, err = batch(tmp_path, *args, "--if-exists", "skip")
    assert code == EXIT_PARTIAL                     # NAG_2's record no longer matches
    rows = {r["sample_id"]: r for r in summary(tmp_path)}
    assert rows["7LF9_D_carbohydrate_NAG_1"]["status"] == "skipped"
    assert rows["7LF9_D_carbohydrate_NAG_5"]["status"] == "ok"   # incomplete: redone
    assert "is missing" in rows["7LF9_D_carbohydrate_NAG_2"]["message"]
    code, _out, err = batch(tmp_path, *args, "--if-exists", "overwrite")
    assert code == EXIT_OK, err


def test_a_keys_only_calculation_does_not_satisfy_a_request_for_triplets(tmp_path):
    args = ["batch", "--manifest", SHEET, "--pdb-dir", STRUCTURES, "--sheet", "tsr",
            "--method", "drug", "--inapplicable", "skip", "--output", tmp_path / "out"]
    assert run(args + ["--no-triplets"])[0] == EXIT_OK
    code, _out, _err = run(args + ["--if-exists", "skip"])
    assert code == EXIT_PARTIAL
    assert all("--no-triplets" in r["message"] for r in summary(tmp_path)
               if r["status"] == "failed")


def test_unrelated_files_in_the_folder_are_never_touched(tmp_path):
    root = tmp_path / "out"
    root.mkdir()
    (root / "my_notes.txt").write_text("keep me")
    (root / "6G9M_A_drug_ESW_404.keys").write_text("not BioTSR's")
    code, _out, _err = batch(tmp_path, "--sheet", "tsr", "--method", "drug",
                             "--inapplicable", "skip", "--if-exists", "overwrite")
    assert code == EXIT_PARTIAL
    assert (root / "my_notes.txt").read_text() == "keep me"
    assert (root / "6G9M_A_drug_ESW_404.keys").read_text() == "not BioTSR's"
    failed = [r for r in summary(tmp_path) if r["status"] == "failed"]
    assert [r["sample_id"] for r in failed] == ["6G9M_A_drug_ESW_404"]
    assert "no BioTSR completion record lists them" in failed[0]["message"]


def test_two_requests_with_one_prefix_are_refused_before_anything_runs(tmp_path):
    manifest = tmp_path / "m.csv"
    manifest.write_text("sample_id,pdb_id,method,chain\nsame,1UBQ,protein,A\n"
                        "SAME,3PTB,protein,A\n")
    code, _out, err = batch(tmp_path, sheet=str(manifest))
    assert code == EXIT_USAGE
    assert "would each be written by more than one different request" in err
    assert not (tmp_path / "out").exists()


def test_a_per_residue_sweep_reserves_its_prefix_for_its_residues(tmp_path):
    manifest = tmp_path / "m.csv"
    manifest.write_text(
        "sample_id,pdb_id,method,chain,residue,residue_name\n"
        ",1UBQ,amino_acid_atomic,A,,LEU\n"
        "1UBQ_A_amino_acid_atomic_perresidue-LEU_LEU_15,1UBQ,protein,A,,\n")
    code, _out, err = batch(tmp_path, sheet=str(manifest))
    assert code == EXIT_USAGE and "more than one different request" in err


def test_the_analysis_reads_the_folder_by_its_completion_records(tmp_path):
    assert batch(tmp_path, "--sheet", "tsr", "--method", "carbohydrate",
                 "--inapplicable", "skip")[0] == EXIT_OK
    code, _out, err = run(["analyze", "common-keys", "--from-results", tmp_path / "out",
                           "--output", tmp_path / "an"])
    assert code == EXIT_OK, err
    report = json.load(open(str(tmp_path / "an" / "analysis" / "common_keys" /
                                "all_samples.json")))
    names = [r["fileName"] for r in report["rows"]]
    assert names == sorted(names)                     # by record path: reproducible
    assert names == ["7LF9_D_carbohydrate_MAN_allcopies", "7LF9_D_carbohydrate_NAG_1",
                     "7LF9_D_carbohydrate_NAG_2", "7LF9_D_carbohydrate_NAG_5"]
    # A key file no completion record lists is refused, not analysed.
    shutil.copy(str(tmp_path / "out" / "7LF9_D_carbohydrate_NAG_1.keys"),
                str(tmp_path / "out" / "stray.keys"))
    code, _out, err = run(["analyze", "common-keys", "--from-results", tmp_path / "out",
                           "--output", tmp_path / "an2"])
    assert code == EXIT_USAGE and "not part of a completed calculation" in err


def test_a_folder_of_the_candidate_layout_is_refused_by_name(tmp_path):
    old = tmp_path / "old" / "1UBQ_A_protein"
    old.mkdir(parents=True)
    (old / "run.json").write_text("{}")
    (old / "1UBQ_A.keys_29_35").write_text("key\tcategory\tfrequency\n")
    code, _out, err = run(["analyze", "common-keys", "--from-results", tmp_path / "old",
                           "--output", tmp_path / "an"])
    assert code == EXIT_USAGE and "release candidate" in err
    code, _out, err = run(["describe", old])
    assert code == EXIT_USAGE and "release candidate" in err


def test_describe_lists_a_folder_and_reads_a_prefix(tmp_path):
    assert batch(tmp_path, "--sheet", "tsr", "--method", "drug",
                 "--inapplicable", "skip")[0] == EXIT_OK
    code, out, err = run(["describe", tmp_path / "out"])
    assert code == EXIT_OK, err
    assert "6G9M_A_drug_ESW_404" in out and "complete" in out and "keys only" in out
    code, out, err = run(["describe", tmp_path / "out" / "3PTB_A_drug_BEN_1"])
    assert code == EXIT_OK, err
    assert "source rows      3 of wide_sheet.csv (participant 1)" in out


def test_method_selects_calculations_of_a_shared_folder_for_analysis(tmp_path):
    """One folder, two key spaces: --method analyses one subset at a time;
    calculations not selected are still checked for completeness."""
    # One aa_selection row names residues 1 and 2, which its chain lacks:
    # that row is partial, so the batch exits 4 (ADR-106).
    assert batch(tmp_path, "--sheet", "tsr", "--method", "carbohydrate",
                 "--method", "aa_selection", "--inapplicable", "skip")[0] == 4
    code, _out, err = run(["analyze", "common-keys", "--from-results", tmp_path / "out",
                           "--output", tmp_path / "an0"])
    assert code == EXIT_USAGE and "do not share a key space" in err
    code, _out, err = run(["analyze", "common-keys", "--from-results", tmp_path / "out",
                           "--method", "carbohydrate", "--output", tmp_path / "an1"])
    assert code == EXIT_OK, err
    report = json.load(open(str(tmp_path / "an1" / "analysis" / "common_keys" /
                                "all_samples.json")))
    assert {s["method"] for s in report["samples"]} == {"carbohydrate"}
    code, _out, err = run(["analyze", "common-keys", "--from-results", tmp_path / "out",
                           "--method", "drug", "--output", tmp_path / "an2"])
    assert code == EXIT_USAGE and "no completed calculation of drug" in err
