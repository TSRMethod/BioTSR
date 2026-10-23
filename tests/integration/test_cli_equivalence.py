"""**The CLI gate: the CLI is never a second scientific implementation.**

For every representative workflow, the result the CLI serializes must be
identical -- key for key, count for count, fingerprint for fingerprint -- to the
result a direct ``analyze()`` call returns for the same input, method, profile
and parameters.  Not "close", not "the same number of keys": the same mapping.

The comparison is made two ways on purpose:

1. **In process** -- run the spec builder and executor and compare ``TSRResult``
   objects field by field.  This catches an argument the CLI routes to the wrong
   keyword, which a file comparison could hide if both sides happened to be wrong
   in the same way.
2. **On disk** -- run ``main()`` end to end and compare the files it wrote with
   the files the API's own serializer writes for the direct result.  This catches
   a CLI that computes the right thing and then serializes it differently.

If these two ever disagree with each other, the CLI has grown an opinion.
"""

from __future__ import annotations

import filecmp
import io
import json
import os
import contextlib

import pytest

import biotsr
from biotsr.cli import EXIT_OK, main
from biotsr.cli.spec import build_spec, execute
from biotsr.io import native

# ---------------------------------------------------------------------------
# The workflow matrix.
#
# (tag, CLI values, direct-API callable). Each covers one representative
# workflow shape. Only small pinned fixtures are used, and no new scientific
# reference is created -- the API's own answer is the reference.
# ---------------------------------------------------------------------------


FIXTURES = (
    "1UBQ", "1TUP", "1BRS", "3PTB", "5WT9", "1BNA", "1KX5", "2AC0", "7LF9",
)


def _cases(paths):
    ubq, brs, ptb, wt9, bna, kx5, ac0, lf9 = (
        paths["1UBQ"], paths["1BRS"], paths["3PTB"],
        paths["5WT9"], paths["1BNA"], paths["1KX5"], paths["2AC0"],
        paths["7LF9"],
    )
    return [
        (
            "protein_calpha",
            {"method": "protein", "pdb": ubq, "chain": "A"},
            lambda: biotsr.analyze(ubq, method="protein", chain="A"),
        ),
        (
            "protein_calpha_grouped_mirror_size",
            {"method": "protein", "pdb": ubq, "chain": "A", "aa_grouping": True,
             "mirror": True, "size_filter": 25.0},
            lambda: biotsr.analyze(
                ubq, method="protein", chain="A", aa_grouping=True, mirror=True,
                size_filter=25.0,
            ),
        ),
        (
            "protein_calpha_sse",
            {"method": "protein", "pdb": brs, "chain": "A", "sse": True},
            lambda: biotsr.analyze(brs, method="protein", chain="A", sse=True),
        ),
        (
            "protein_calpha_sse_mirror",
            {"method": "protein", "pdb": brs, "chain": "A", "sse": True,
             "mirror": True},
            lambda: biotsr.analyze(
                brs, method="protein", chain="A", sse=True, mirror=True
            ),
        ),
        (
            "aa_selection_sse",
            {"method": "aa_selection", "pdb": brs, "chain": "A", "sse": True,
             "residues": "10-40;60-70"},
            lambda: biotsr.analyze(
                brs, method="aa_selection", chain="A", sse=True,
                residues="10-40;60-70",
            ),
        ),
        (
            "amino_acid_atomic_one_residue",
            {"method": "amino_acid_atomic", "pdb": ubq, "chain": "A",
             "residue": "15", "residue_name": "LEU"},
            lambda: biotsr.analyze(
                ubq, method="amino_acid_atomic", chain="A", residue="15",
                residue_name="LEU",
            ),
        ),
        (
            "amino_acid_atom_selection",
            {"method": "amino_acid_atom_selection", "pdb": ubq, "chain": "A",
             "atoms": "LEU_15_CB;ILE_13_CG1;VAL_17_CG2;LYS_11_NZ"},
            lambda: biotsr.analyze(
                ubq, method="amino_acid_atom_selection", chain="A",
                atoms="LEU_15_CB;ILE_13_CG1;VAL_17_CG2;LYS_11_NZ",
            ),
        ),
        (
            "aa_ind_mix",
            {"method": "aa_ind_mix", "pdb": ubq, "chain": "A",
             "residue_types": "LEU;ILE"},
            lambda: biotsr.analyze(
                ubq, method="aa_ind_mix", chain="A", residue_types="LEU;ILE"
            ),
        ),
        (
            "drug",
            {"method": "drug", "pdb": ptb, "chain": "A", "entity_name": "BEN",
             "entity_id": "1"},
            lambda: biotsr.analyze(
                ptb, method="drug", chain="A", entity_name="BEN", entity_id="1"
            ),
        ),
        (
            "carbohydrate_chain",
            {"method": "carbohydrate_chain", "pdb": wt9, "chain": "A"},
            lambda: biotsr.analyze(wt9, method="carbohydrate_chain", chain="A"),
        ),
        (
            "nucleotide",
            {"method": "nucleotide", "pdb": bna, "chain": "A",
             "entity_name": "DC", "entity_id": "3"},
            lambda: biotsr.analyze(
                bna, method="nucleotide", chain="A", entity_name="DC",
                entity_id="3",
            ),
        ),
        (
            "protein_protein_cross_all_atom",
            {"method": "protein_protein_cross", "pdb": brs, "chain": "A",
             "chain_2": "D", "contact_cutoff": 4.0},
            lambda: biotsr.analyze_cross(
                brs, method="protein_protein_cross", chain="A", chain_2="D",
                contact_cutoff=4.0,
            ),
        ),
        (
            "protein_protein_calpha_cross",
            {"method": "protein_protein_calpha_cross", "pdb": brs, "chain": "A",
             "chain_2": "D", "contact_cutoff": 8.0},
            lambda: biotsr.analyze_cross(
                brs, method="protein_protein_calpha_cross", chain="A",
                chain_2="D", contact_cutoff=8.0,
            ),
        ),
        (
            "drug_protein_cross",
            {"method": "drug_protein_cross", "pdb": ptb, "chain": "A",
             "chain_2": "A", "entity_name": "BEN", "entity_id": "1",
             "contact_cutoff": 4.0},
            lambda: biotsr.analyze_cross(
                ptb, method="drug_protein_cross", chain="A", chain_2="A",
                entity_name="BEN", entity_id="1", contact_cutoff=4.0,
            ),
        ),
        (
            "cross_size_filtered",
            {"method": "protein_protein_cross", "pdb": brs, "chain": "A",
             "chain_2": "D", "contact_cutoff": 4.0, "size_filter": 12.0},
            lambda: biotsr.analyze_cross(
                brs, method="protein_protein_cross", chain="A", chain_2="D",
                contact_cutoff=4.0, size_filter=12.0,
            ),
        ),
        (
            "nucleotide_protein_cross",
            {"method": "nucleotide_protein_cross", "pdb": ac0, "chain": "E",
             "chain_2": "B", "entity_name": "DG", "entity_id": "8"},
            lambda: biotsr.analyze_cross(
                ac0, method="nucleotide_protein_cross", chain="E", chain_2="B",
                entity_name="DG", entity_id="8",
            ),
        ),
        (
            "carbohydrate_protein_cross",
            {"method": "carbohydrate_protein_cross", "pdb": lf9, "chain": "C",
             "chain_2": "A", "entity_name": "NAG", "entity_id": "1"},
            lambda: biotsr.analyze_cross(
                lf9, method="carbohydrate_protein_cross", chain="C", chain_2="A",
                entity_name="NAG", entity_id="1",
            ),
        ),
        (
            "protein_calpha_on_a_helix_rich_chain",
            {"method": "protein", "pdb": kx5, "chain": "A", "sse": True,
             "size_filter": 20.0},
            lambda: biotsr.analyze(
                kx5, method="protein", chain="A", sse=True, size_filter=20.0
            ),
        ),
    ]


@pytest.fixture(scope="module")
def paths():
    from ..conftest import structure

    return {name: structure(name) for name in FIXTURES}


@pytest.fixture(scope="module")
def cases(paths):
    return {tag: (values, direct) for tag, values, direct in _cases(paths)}


TAGS = [tag for tag, _v, _d in _cases({n: n for n in FIXTURES})]


# ---------------------------------------------------------------------------
# 1. In process
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("tag", TAGS)
def test_cli_result_is_identical_to_the_api_result(tag, cases):
    values, direct = cases[tag]
    expected = direct()
    outcome = execute(build_spec(dict(values)))
    assert len(outcome.results) == 1
    produced = list(outcome.results.values())[0]

    assert produced.method == expected.method
    assert produced.structure_id == expected.structure_id
    assert produced.chain == expected.chain
    # The numbers themselves.
    assert dict(produced.frequencies) == dict(expected.frequencies)
    assert produced.key_order == expected.key_order
    assert dict(produced.category_totals()) == dict(expected.category_totals())
    # The candidate ledger, including the two ways a candidate can be dropped.
    assert produced.counts.to_dict() == expected.counts.to_dict()
    # Both identities.
    assert produced.key_space.fingerprint == expected.key_space.fingerprint
    assert produced.sampling.fingerprint == expected.sampling.fingerprint
    # Selection and provenance.
    assert produced.selection.to_dict() == expected.selection.to_dict()
    assert produced.annotations["profile"] == expected.annotations["profile"]
    assert (
        produced.annotations["legacy_compatibility_profile"]
        == expected.annotations["legacy_compatibility_profile"]
    )
    # The full self-describing sidecar payload, which is what lands on disk.
    assert produced.metadata() == expected.metadata()


@pytest.mark.parametrize("tag", [t for t in TAGS if "sse" in t])
def test_sse_stratification_is_identical(tag, cases):
    values, direct = cases[tag]
    expected = direct()
    produced = list(execute(build_spec(dict(values))).results.values())[0]
    assert produced.annotations["sse"] == expected.annotations["sse"]
    assert produced.annotations["sse_coverage"] == expected.annotations["sse_coverage"]
    # Every key's per-category breakdown, not merely the totals.
    assert {k: dict(v.by_category) for k, v in produced.keys.items()} == {
        k: dict(v.by_category) for k, v in expected.keys.items()
    }


def test_sse_annotation_changes_no_key_through_the_cli(cases):
    """The SSE key-invariance guarantee, re-asserted through the command line."""
    plain = execute(build_spec(
        {"method": "protein", "pdb": cases["protein_calpha_sse"][0]["pdb"],
         "chain": "A"}
    ))
    annotated = execute(build_spec(dict(cases["protein_calpha_sse"][0])))
    a = list(plain.results.values())[0]
    b = list(annotated.results.values())[0]
    assert dict(a.frequencies) == dict(b.frequencies)
    assert a.key_space.fingerprint == b.key_space.fingerprint
    assert a.sampling.fingerprint == b.sampling.fingerprint


# ---------------------------------------------------------------------------
# 2. On disk
# ---------------------------------------------------------------------------


def _cli(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main(argv)
    return code, out.getvalue(), err.getvalue()


def _flags(values):
    argv = []
    for key, value in sorted(values.items()):
        flag = "--" + key.replace("_", "-")
        if value is True:
            argv.append(flag)
        else:
            argv.extend([flag, str(value)])
    return argv


@pytest.mark.parametrize("tag", TAGS)
def test_the_files_the_cli_writes_are_the_serializers_own(tag, cases, tmp_path):
    values, direct = cases[tag]
    out = tmp_path / "cli"
    code, stdout, _err = _cli(["run"] + _flags(values) + ["-o", str(out), "-q"])
    assert code == EXIT_OK

    record = stdout.strip()
    assert record.endswith(".run.json") and os.path.isfile(record)
    folder = os.path.dirname(record)

    # The API names the files by the same grammar the CLI used for the request.
    reference = tmp_path / "api"
    reference.mkdir()
    native.save(direct(), str(reference))

    produced = sorted(f for f in os.listdir(folder) if not f.endswith(".run.json"))
    expected = sorted(os.listdir(str(reference)))
    assert produced == expected
    assert os.path.basename(record)[: -len(".run.json")] == native.prefix_of(expected[0])
    for name in expected:
        assert filecmp.cmp(
            os.path.join(folder, name),
            os.path.join(str(reference), name),
            shallow=False,
        ), "%s differs between the CLI and the API" % name


def test_the_per_entity_sweep_matches_analyze_entities(tmp_path, paths):
    expected = biotsr.analyze_entities(
        paths["1UBQ"], method="amino_acid_atomic", chain="A", keep_triplets=False
    )
    outcome = execute(
        build_spec(
            {"method": "amino_acid_atomic", "pdb": paths["1UBQ"], "chain": "A",
             "keep_triplets": False}
        )
    )
    assert sorted(outcome.results) == sorted(expected)
    for name, result in outcome.results.items():
        assert dict(result.frequencies) == dict(expected[name].frequencies)
        assert result.metadata() == expected[name].metadata()


def test_the_standard_strand_identity_survives_the_cli_on_the_multi_sheet_chain(
    tmp_path, paths
):
    """1TUP chain A is the corpus's only multi-sheet chain, so it is the one
    place where how a beta-strand is identified decides a category (ADR-063).
    The CLI must carry the standard convention through unchanged -- the same
    keys and the same per-key category breakdown as the API.
    """
    standard = list(
        execute(
            build_spec(
                {"method": "protein", "pdb": paths["1TUP"], "chain": "A",
                 "sse": True, "keep_triplets": False}
            )
        ).results.values()
    )[0]
    expected_standard = biotsr.analyze(
        paths["1TUP"], method="protein", chain="A", sse=True, keep_triplets=False
    )
    assert standard.metadata() == expected_standard.metadata()
    assert {k: dict(v.by_category) for k, v in standard.keys.items()} == {
        k: dict(v.by_category) for k, v in expected_standard.keys.items()
    }
    assert standard.counts.n_triangles == expected_standard.counts.n_triangles
    assert sum(standard.category_totals().values()) == standard.counts.n_triangles


def test_the_cli_writes_one_file_set_per_entity(tmp_path, paths):
    out = tmp_path / "cli"
    code, stdout, _err = _cli(
        ["run", "-m", "amino_acid_atomic", "--pdb", paths["1UBQ"], "-c", "A",
         "-o", str(out), "-q", "--no-triplets"]
    )
    assert code == EXIT_OK
    record = stdout.strip()
    folder = os.path.dirname(record)
    keys = [f for f in os.listdir(folder) if f.endswith(".keys")]
    assert len(keys) == 75
    assert all(k.startswith("1UBQ_A_amino_acid_atomic_perresidue_") for k in keys)
    payload = json.load(open(record))
    assert len(payload["results"]) == 75 and payload["per_residue"] is True
    assert sorted(e["files"]["keys"] for e in payload["results"]) == sorted(keys)
