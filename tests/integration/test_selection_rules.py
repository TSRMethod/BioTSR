"""The standard-method selection rules, as public behaviour.

* **One molecule per single-entity method.**  ``drug``, ``carbohydrate``,
  ``nucleotide`` and their ``*_atom_selection`` variants analyse one copy: an
  omitted ``entity_id`` is allowed only when exactly one copy is in the chain.
  rc1 silently pooled every copy (5WT9 NAG without an id: two monomers, 3,276
  triangles, no warning).
* **Pooled methods stay pooled** (``carbohydrate_chain``, ``nucleotide_chain``).
* **Nothing requested is silently missing.**  An explicit selection that matches
  nothing is refused; a partially matched one says what was not found; a
  selection too small for a triangle is labelled as such.
* **Every method runs under its standard profile**, through Python as on the
  command line.
* **Arguments a method does not read are refused by name** (rc1 accepted them;
  ``protein_protein_calpha_cross --entity-name FOO`` put FOO in the output name).
* **Metal ions are labelled as the elements they are** (the historical defect
  recorded as ADR-082 is not reachable).
"""

from __future__ import annotations

import json
import math
from dataclasses import replace

import pytest

import biotsr
from biotsr import InvalidRequest
from biotsr import batch as biotsr_batch
from biotsr.cli.main import main

from ..conftest import structure


# ---------------------------------------------------------------------------
# Molecule instances (ADR-096): an identifier is one copy; a name alone is
# every copy on the chain, in one calculation
# ---------------------------------------------------------------------------


def test_an_omitted_entity_id_is_every_copy_even_when_there_is_one():
    implicit = biotsr.analyze(structure("3PTB"), method="drug", chain="A",
                              entity_name="BEN", keep_triplets=False)
    explicit = biotsr.analyze(structure("3PTB"), method="drug", chain="A",
                              entity_name="BEN", entity_id="1", keep_triplets=False)
    assert implicit.annotations["target"] == "BEN_allcopies"
    assert implicit.metadata()["provenance"]["target"] == "BEN_allcopies"
    assert explicit.annotations["target"] == "BEN_1"
    assert implicit.selection.residues_used == ("1_BEN",)
    assert implicit.frequencies == explicit.frequencies
    assert implicit.key_space.fingerprint == explicit.key_space.fingerprint
    scopes = [r.annotations["selection_request"]["sides"][0]["scope"]
              for r in (implicit, explicit)]
    assert scopes == ["all_copies", "one_instance"]


@pytest.mark.parametrize("method,extra,points", [
    ("carbohydrate", {}, 28),
    ("carbohydrate_atom_selection", {"atoms": "C1;C2;O5"}, 6),
])
def test_several_copies_without_an_id_are_pooled_and_recorded(method, extra, points):
    """rc3 refused these; both NAG copies of 5WT9 chain A now form one
    calculation, with triangles spanning the two."""
    result = biotsr.analyze(structure("5WT9"), method=method, chain="A",
                            entity_name="NAG", keep_triplets=False, **extra)
    side = result.annotations["selection_request"]["sides"][0]
    assert side["instances"] == ["NAG_1", "NAG_2"]
    assert result.counts.n_points == points
    assert result.counts.n_triangles == math.comb(points, 3)


def test_a_nucleotide_named_without_an_id_is_every_copy_on_the_strand():
    """1BNA chain A: DC 1, 3, 9 and 11; the incomplete DC 1 is excluded by the
    completeness rule, named in a warning, and the others pooled."""
    result = biotsr.analyze(structure("1BNA"), method="nucleotide", chain="A",
                            entity_name="DC", keep_triplets=False)
    side = result.annotations["selection_request"]["sides"][0]
    assert side["scope"] == "all_copies"
    assert "DC_1" in side["excluded"] and "DC_1" not in side["instances"]
    assert len(side["instances"]) >= 2
    assert any("completeness rule" in w and "1_DC" in w for w in result.warnings)


def test_an_explicit_id_that_does_not_exist_is_refused():
    with pytest.raises(InvalidRequest, match="BEN 2 is not in chain A.*ids 1"):
        biotsr.analyze(structure("3PTB"), method="drug", chain="A",
                       entity_name="BEN", entity_id="2")


def test_an_absent_entity_is_refused_naming_what_is_there():
    with pytest.raises(InvalidRequest) as excinfo:
        biotsr.analyze(structure("3PTB"), method="drug", chain="A", entity_name="XYZ")
    assert "XYZ is not in chain A" in str(excinfo.value)
    assert "BEN" in str(excinfo.value)


def test_the_cli_names_the_all_copies_scope_in_the_output(tmp_path, capsys):
    out = tmp_path / "out"
    code = main(["run", "--method", "carbohydrate", "--pdb", structure("5WT9"),
                 "--chain", "A", "--entity-name", "NAG", "--output", str(out),
                 "--no-triplets", "-q"])
    assert code == 0
    record = json.loads((out / "5WT9_A_carbohydrate_NAG_allcopies.run.json").read_text())
    assert record["results"][0]["target"] == "NAG_allcopies"


def test_pooled_methods_still_pool_every_copy():
    chain = biotsr.analyze(structure("5WT9"), method="carbohydrate_chain", chain="A",
                           keep_triplets=False)
    assert len(chain.selection.residues_used) >= 2
    assert chain.annotations["target"] == "A"
    two_nag = biotsr.analyze(structure("5WT9"), method="carbohydrate_chain",
                             chain="A", entity_names="NAG", keep_triplets=False)
    assert two_nag.selection.residues_used == ("1_NAG", "2_NAG")
    assert two_nag.counts.n_points == 28
    assert two_nag.counts.n_triangles == math.comb(28, 3) == 3276


# ---------------------------------------------------------------------------
# Empty and partially matched selections (decision D)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method,kwargs", [
    ("aa_selection", {"residues": "600;601"}),
    ("aa_ind", {"residue_types": "TRP"}),
    ("aa_ind_specific", {"residues": "LYS_600;GLU601"}),
    ("amino_acid_atom_selection", {"atoms": "LYS_600_NZ;GLU_601_OE1"}),
    ("amino_acid_atomic", {"residue": "600", "residue_name": "LYS"}),
    ("amino_acid_atomic", {"residue": "15", "residue_name": "GLY"}),
    ("carbohydrate_chain", {"entity_names": "MAN"}),
])
def test_an_explicit_selection_that_matches_nothing_is_refused(method, kwargs):
    with pytest.raises(InvalidRequest):
        biotsr.analyze(structure("1UBQ") if method != "carbohydrate_chain"
                       else structure("5WT9"),
                       method=method, chain="A", **kwargs)


def test_named_atoms_absent_from_the_entity_are_refused():
    with pytest.raises(InvalidRequest, match="none of the requested atoms"):
        biotsr.analyze(structure("3PTB"), method="drug_atom_selection", chain="A",
                       entity_name="BEN", entity_id="1", atoms="ZZ1;ZZ2")


def test_partly_missing_items_are_named_in_a_warning():
    result = biotsr.analyze(structure("1UBQ"), method="aa_selection", chain="A",
                            residues="10-14;500", keep_triplets=False)
    assert result.counts.n_points == 5
    assert "requested but not found: 500" in result.warnings
    mixed = biotsr.analyze(structure("1UBQ"), method="aa_ind_mix", chain="A",
                           residue_types="LYS;TRP", keep_triplets=False)
    assert any(w.startswith("requested but not found: residue type TRP")
               for w in mixed.warnings)
    atoms = biotsr.analyze(structure("3PTB"), method="drug_atom_selection",
                           chain="A", entity_name="BEN", entity_id="1",
                           atoms="C1;C2;C3;ZZ9", keep_triplets=False)
    assert atoms.counts.n_points == 3
    assert "requested but not found: ZZ9" in atoms.warnings


def test_a_selection_too_small_for_a_triangle_is_labelled():
    result = biotsr.analyze(structure("1UBQ"), method="aa_selection", chain="A",
                            residues="1-2", keep_triplets=False)
    assert result.counts.n_triangles == 0
    assert any(w.startswith("no triangles: the selection has 2 points")
               for w in result.warnings)


def test_the_cli_shows_the_warning_and_records_it(tmp_path, capsys):
    out = tmp_path / "out"
    code = main(["run", "--method", "aa_selection", "--pdb", structure("1UBQ"),
                 "--chain", "A", "--residues", "1-2;500", "--output", str(out),
                 "--no-triplets"])
    assert code == 4  # written, but residue 500 is missing: partial (ADR-106)
    err = capsys.readouterr().err
    assert "requested but not found: 500" in err
    assert "no triangles" in err
    run = json.loads((out / "1UBQ_A_aa_selection_res-1-2+500.run.json").read_text())
    assert any(w.startswith("no triangles") for w in run["results"][0]["warnings"])


def test_a_batch_labels_a_zero_triangle_row(tmp_path):
    manifest = tmp_path / "m.csv"
    manifest.write_text("sample_id,pdb_path,method,chain,residues\n"
                        "tiny,%s,aa_selection,A,1-2\nfine,%s,aa_selection,A,1-5\n"
                        % (structure("1UBQ"), structure("1UBQ")))
    report = biotsr_batch.run(str(manifest), str(tmp_path / "out"), triplets=False)
    by_id = {r["sample_id"]: r for r in report.records}
    assert by_id["tiny"]["status"] == "no_triangles"
    assert by_id["tiny"]["message"].startswith("no triangles")
    assert by_id["fine"]["status"] == "ok"
    assert report.n_failed == 0


# ---------------------------------------------------------------------------
# Standard profiles only (decision B)
# ---------------------------------------------------------------------------


def test_the_methods_own_standard_profile_is_accepted():
    recipe = biotsr.get_method("protein")
    a = biotsr.analyze(structure("1UBQ"), method="protein", chain="A",
                       profile=recipe.default_profile, keep_triplets=False)
    b = biotsr.analyze(structure("1UBQ"), method="protein", chain="A",
                       keep_triplets=False)
    assert a.frequencies == b.frequencies


@pytest.mark.parametrize("profile", [
    "standard",
    "legacy_tsr_package_v1",
    biotsr.STANDARD_AMINOACID_ATOMIC_V1,
    replace(biotsr.STANDARD_PROTEIN_CALPHA_V1, name="my_variant"),
])
def test_any_other_profile_is_refused(profile):
    with pytest.raises(InvalidRequest, match="standard profile"):
        biotsr.analyze(structure("1UBQ"), method="protein", chain="A",
                       profile=profile)


def test_a_cross_method_refuses_a_foreign_profile_and_a_mode():
    with pytest.raises(InvalidRequest):
        biotsr.analyze_cross(structure("1BRS"), method="protein_protein_cross",
                             chain="A", chain_2="D",
                             profile=biotsr.STANDARD_PROTEIN_CALPHA_V1)
    with pytest.raises(InvalidRequest, match="mode="):
        biotsr.analyze_cross(structure("1BRS"), method="protein_protein_cross",
                             chain="A", chain_2="D", mode="legacy")


def test_the_cli_has_no_profile_or_format_option(capsys):
    for option in ("--profile", "--format"):
        with pytest.raises(SystemExit) as excinfo:
            main(["run", "--method", "protein", "--pdb", structure("1UBQ"),
                  "--chain", "A", "--output", "out", option, "x"])
        assert excinfo.value.code == 2


# ---------------------------------------------------------------------------
# Arguments a method does not read (finding 3)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method,kwargs,named", [
    ("protein", {"residues": "1-5"}, "residues"),
    ("protein", {"entity_name": "BEN"}, "entity_name"),
    ("protein", {"chain_2": "B"}, "chain_2"),
    ("aa_ind", {"residue_name": "TYR"}, "residue_name"),
    ("drug", {"residue_types": "LYS", "entity_name": "BEN"}, "residue_types"),
    ("carbohydrate_chain", {"entity_name": "NAG"}, "entity_name"),
])
def test_python_refuses_an_unused_argument_by_name(method, kwargs, named):
    pdb = structure("3PTB") if method in ("drug",) else structure("5WT9") \
        if method == "carbohydrate_chain" else structure("1UBQ")
    with pytest.raises(InvalidRequest) as excinfo:
        biotsr.analyze(pdb, method=method, chain="A", **kwargs)
    assert named in str(excinfo.value)


def test_the_cli_refuses_unused_cross_arguments_and_writes_nothing(tmp_path, capsys):
    out = tmp_path / "out"
    code = main(["run", "--method", "protein_protein_calpha_cross", "--pdb",
                 structure("1BRS"), "--chain", "A", "--chain-2", "D",
                 "--entity-name", "FOO", "--atoms", "XX", "--output", str(out),
                 "--no-triplets"])
    assert code == 2
    err = capsys.readouterr().err
    assert "--entity-name" in err and "--atoms" in err and "--calpha-residues" in err
    assert not out.exists()


def test_residue_name_without_residue_is_a_named_per_residue_sweep(tmp_path, capsys):
    """Not silently a sweep of everything: --residue-name LEU alone runs each
    leucine, says so in the prefix, and --residue alone is refused."""
    out = tmp_path / "out"
    code = main(["run", "--method", "amino_acid_atomic", "--pdb", structure("1UBQ"),
                 "--chain", "A", "--residue-name", "LEU", "--output", str(out),
                 "--no-triplets", "-q"])
    assert code == 0
    record = json.loads(
        (out / "1UBQ_A_amino_acid_atomic_perresidue-LEU.run.json").read_text())
    assert record["per_residue"] is True
    assert {e["entity"] for e in record["results"]} == {
        "8_LEU", "15_LEU", "43_LEU", "50_LEU", "56_LEU", "67_LEU", "69_LEU",
        "71_LEU", "73_LEU"}
    code = main(["run", "--method", "amino_acid_atomic", "--pdb", structure("1UBQ"),
                 "--chain", "A", "--residue", "15", "--output", str(tmp_path / "x")])
    assert code == 2
    assert "--residue without --residue-name" in capsys.readouterr().err


def test_methods_lists_each_methods_own_arguments(capsys):
    assert main(["methods", "--json"]) == 0
    records = {r["name"]: r for r in json.loads(capsys.readouterr().out)}
    assert records["protein_protein_calpha_cross"]["target_arguments"] == [
        "calpha_residues", "calpha_residues_2", "chain_2", "contact_cutoff"]
    assert "entity_name" not in records["protein_protein_cross"]["target_arguments"]
    assert records["protein"]["target_arguments"] == []
    assert all("profiles" not in r for r in records.values())
    assert records["protein"]["profile"] == "standard_protein_calpha_v1"


# ---------------------------------------------------------------------------
# Metal ions (decision E)
# ---------------------------------------------------------------------------


def _labels(result):
    return {(t.identities[i].residue_name, t.labels[i])
            for t in result.triplets for i in range(3)}


def test_manganese_is_labelled_manganese():
    """1KX5 chain I: seven Mn ions, every vertex labelled 100. The historical
    narrow atom-name slice spelled them N and labelled them nitrogen (108)."""
    result = biotsr.analyze(structure("1KX5"), method="carbohydrate_chain",
                            chain="I", entity_names="MN")
    assert result.counts.n_points == 7
    assert result.counts.n_triangles == 35
    assert _labels(result) == {("MN", 100)}


def test_zinc_is_labelled_zinc():
    """One Zn per chain in 1TUP, so no triangle -- but the point's label is the
    element's code (105), never nitrogen."""
    from biotsr.selection.atomic import STANDARD_HETERO_ADMISSION, select_entity_atoms

    structure_ = biotsr.load_structure(structure("1TUP"))
    selection = select_entity_atoms(
        structure_, "A", biotsr.ELEMENT_LABELS_112,
        admission=STANDARD_HETERO_ADMISSION, entity_names=("ZN",),
    )
    assert [p.label for p in selection.points] == [105]
    result = biotsr.analyze(structure("1TUP"), method="carbohydrate_chain",
                            chain="A", entity_names="ZN")
    assert result.counts.n_points == 1
    assert any(w.startswith("no triangles") for w in result.warnings)


def test_the_parser_reads_the_full_atom_name_field_and_keeps_the_calpha_token():
    """The structural reason the metal-ion defect cannot recur (ADR-108):
    every atom's identity carries PDB columns 13-16; only the Calpha test reads
    columns 14-16, so calcium never becomes a Calpha."""
    structure_ = biotsr.load_structure(structure("1TUP"))
    zinc = [a for a in structure_.atoms if a.residue_name == "ZN"]
    assert zinc and all(a.atom_name == "ZN" and a.calpha_name == "N" for a in zinc)
    calpha = [a for a in structure_.atoms if a.calpha_name == "CA" and not a.hetero]
    assert calpha and all(a.atom_name == "CA" for a in calpha)
