#!/usr/bin/env python3
"""Command-line gates: the command line, checked as a contract rather than a feature.

Four categories, each answering a different question:

``equivalence``
    Is what the CLI produces the same science the Python API produces?  This is
    the central gate: the CLI must never become a second implementation.
``contract``
    Does the CLI keep the promises a script depends on -- deterministic output
    prefixes, deterministic exit codes, one shared output folder in which a
    calculation is installed with its completion record written last (so it is
    never half-written and never mistaken for complete), a failure that stays
    attached to its input row, and the comparison rule applied to saved results?
``discovery``
    Can a user find out what the package does without reading its source, and is
    what they are told generated from the registry rather than restated?
``safety``
    Does the CLI refuse the things it must refuse -- a scientific default of its
    own, a profile or output-format choice, a deferred method, an accidental
    overwrite, a traceback in a user's face?

Usage::

    python3 validation/scripts/gates_cli.py [--json PATH]
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import csv
import filecmp
import io
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
WORKSPACE = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(WORKSPACE, "src")
CLI_DIR = os.path.join(SRC, "biotsr", "cli")
FIXTURES = os.path.join(WORKSPACE, "tests", "data", "structures")
WIDE_SHEET = os.path.join(WORKSPACE, "tests", "data", "manifests", "wide_sheet.csv")
EVIDENCE = os.path.join(WORKSPACE, "validation", "evidence", "cli")

sys.path.insert(0, SRC)

import biotsr  # noqa: E402
from biotsr.cli import EXIT_OK, EXIT_PARTIAL, EXIT_USAGE, main as cli_main  # noqa: E402
from biotsr.cli import output as cli_output  # noqa: E402
from biotsr.cli import spec as cli_spec  # noqa: E402
from biotsr.cli.errors import CLIError  # noqa: E402
from biotsr.cli.manifest import CANONICAL_COLUMNS, COLUMN_ALIASES, read_manifest  # noqa: E402
from biotsr.cli.parser import build_parser  # noqa: E402
from biotsr.io import native  # noqa: E402
from biotsr.io import naming  # noqa: E402
from biotsr.methods.arguments import accepted_arguments  # noqa: E402
from biotsr.methods.registry import get_method  # noqa: E402

RESULTS = []


#: The 1.0 methods (README, "Method families"): exactly these are registered
#: and listed, so no deferred candidate strategy is reachable.
RELEASED_METHODS = frozenset((
    "aa_ind", "aa_ind_mix", "aa_ind_specific", "aa_selection",
    "amino_acid_atom_selection", "amino_acid_atomic", "carbohydrate",
    "carbohydrate_atom_selection", "carbohydrate_chain", "carbohydrate_protein_cross",
    "drug", "drug_atom_selection", "drug_protein_cross", "nucleotide",
    "nucleotide_atom_selection", "nucleotide_chain", "nucleotide_protein_cross",
    "protein", "protein_protein_calpha_cross", "protein_protein_cross",
))


def gate(category, name, detail=""):
    def record(passed, note=""):
        RESULTS.append(
            {
                "category": category,
                "name": name,
                "passed": bool(passed),
                "detail": note or detail,
            }
        )
    return record


def check(category, name, detail=""):
    """Decorator: run the body, record pass/fail, never let one gate stop the rest."""
    def decorate(fn):
        try:
            note = fn()
            RESULTS.append(
                {"category": category, "name": name, "passed": True,
                 "detail": note or detail}
            )
        except AssertionError as exc:
            RESULTS.append(
                {"category": category, "name": name, "passed": False,
                 "detail": str(exc)[:300]}
            )
        except Exception as exc:  # noqa: BLE001
            RESULTS.append(
                {"category": category, "name": name, "passed": False,
                 "detail": "%s: %s" % (type(exc).__name__, str(exc)[:260])}
            )
        return fn
    return decorate


def cli(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = cli_main(argv)
    return code, out.getvalue(), err.getvalue()


def fixture(name):
    path = os.path.join(FIXTURES, "%s.pdb" % name)
    if os.path.exists(path):
        return path
    raise AssertionError("pinned fixture %s.pdb is missing" % name)


TMP = tempfile.mkdtemp(prefix="tsr-cli-gates-")


def workdir(name):
    path = os.path.join(TMP, name)
    os.makedirs(path, exist_ok=True)
    return path


# ===========================================================================
# equivalence -- the CLI is never a second implementation
# ===========================================================================

EQUIVALENCE_CASES = [
    (
        "protein_calpha",
        ["-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A"],
        lambda: biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                            keep_triplets=False),
    ),
    (
        "protein_calpha_grouped",
        ["-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A", "--aa-grouping"],
        lambda: biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                            aa_grouping=True, keep_triplets=False),
    ),
    (
        "protein_calpha_mirror",
        ["-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A", "--mirror"],
        lambda: biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                            mirror=True, keep_triplets=False),
    ),
    (
        "protein_calpha_size_filter",
        ["-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
         "--size-filter", "25"],
        lambda: biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                            size_filter=25.0, keep_triplets=False),
    ),
    (
        "protein_calpha_all_three_modifiers",
        ["-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A", "--aa-grouping",
         "--mirror", "--size-filter", "25"],
        lambda: biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                            aa_grouping=True, mirror=True, size_filter=25.0,
                            keep_triplets=False),
    ),
    (
        "protein_calpha_sse",
        ["-m", "protein", "--pdb", fixture("1BRS"), "-c", "A", "--sse"],
        lambda: biotsr.analyze(fixture("1BRS"), method="protein", chain="A",
                            sse=True, keep_triplets=False),
    ),
    (
        "protein_calpha_sse_with_modifiers",
        ["-m", "protein", "--pdb", fixture("1BRS"), "-c", "A", "--sse",
         "--mirror", "--aa-grouping"],
        lambda: biotsr.analyze(fixture("1BRS"), method="protein", chain="A",
                            sse=True, mirror=True, aa_grouping=True,
                            keep_triplets=False),
    ),
    (
        "aa_selection",
        ["-m", "aa_selection", "--pdb", fixture("1UBQ"), "-c", "A",
         "--residues", "10-30;54;70-72"],
        lambda: biotsr.analyze(fixture("1UBQ"), method="aa_selection", chain="A",
                            residues="10-30;54;70-72", keep_triplets=False),
    ),
    (
        "aa_selection_sse",
        ["-m", "aa_selection", "--pdb", fixture("1BRS"), "-c", "A",
         "--residues", "10-40", "--sse"],
        lambda: biotsr.analyze(fixture("1BRS"), method="aa_selection", chain="A",
                            residues="10-40", sse=True, keep_triplets=False),
    ),
    (
        "amino_acid_atomic",
        ["-m", "amino_acid_atomic", "--pdb", fixture("1UBQ"), "-c", "A",
         "--residue", "15", "--residue-name", "LEU"],
        lambda: biotsr.analyze(fixture("1UBQ"), method="amino_acid_atomic",
                            chain="A", residue="15", residue_name="LEU",
                            keep_triplets=False),
    ),
    (
        "amino_acid_atom_selection",
        ["-m", "amino_acid_atom_selection", "--pdb", fixture("1UBQ"), "-c", "A",
         "--atoms", "LEU_15_CB;ILE_13_CG1;VAL_17_CG2;LYS_11_NZ"],
        lambda: biotsr.analyze(fixture("1UBQ"),
                            method="amino_acid_atom_selection", chain="A",
                            atoms="LEU_15_CB;ILE_13_CG1;VAL_17_CG2;LYS_11_NZ",
                            keep_triplets=False),
    ),
    (
        "aa_ind",
        ["-m", "aa_ind", "--pdb", fixture("1UBQ"), "-c", "A",
         "--residue-types", "LEU"],
        lambda: biotsr.analyze(fixture("1UBQ"), method="aa_ind", chain="A",
                            residue_types="LEU", keep_triplets=False),
    ),
    (
        "aa_ind_mix",
        ["-m", "aa_ind_mix", "--pdb", fixture("1UBQ"), "-c", "A",
         "--residue-types", "LEU;ILE"],
        lambda: biotsr.analyze(fixture("1UBQ"), method="aa_ind_mix", chain="A",
                            residue_types="LEU;ILE", keep_triplets=False),
    ),
    (
        "aa_ind_specific",
        ["-m", "aa_ind_specific", "--pdb", fixture("1UBQ"), "-c", "A",
         "--residues", "LEU_15;ILE_13;VAL_17"],
        lambda: biotsr.analyze(fixture("1UBQ"), method="aa_ind_specific", chain="A",
                            residues="LEU_15;ILE_13;VAL_17", keep_triplets=False),
    ),
    (
        "drug",
        ["-m", "drug", "--pdb", fixture("3PTB"), "-c", "A",
         "--entity-name", "BEN", "--entity-id", "1"],
        lambda: biotsr.analyze(fixture("3PTB"), method="drug", chain="A",
                            entity_name="BEN", entity_id="1", keep_triplets=False),
    ),
    (
        "carbohydrate",
        ["-m", "carbohydrate", "--pdb", fixture("5WT9"), "-c", "A",
         "--entity-name", "NAG", "--entity-id", "1"],
        lambda: biotsr.analyze(fixture("5WT9"), method="carbohydrate", chain="A",
                            entity_name="NAG", entity_id="1", keep_triplets=False),
    ),
    (
        "carbohydrate_chain",
        ["-m", "carbohydrate_chain", "--pdb", fixture("5WT9"), "-c", "A"],
        lambda: biotsr.analyze(fixture("5WT9"), method="carbohydrate_chain",
                            chain="A", keep_triplets=False),
    ),
    (
        "nucleotide",
        ["-m", "nucleotide", "--pdb", fixture("1BNA"), "-c", "A",
         "--entity-name", "DC", "--entity-id", "3"],
        lambda: biotsr.analyze(fixture("1BNA"), method="nucleotide", chain="A",
                            entity_name="DC", entity_id="3", keep_triplets=False),
    ),
    (
        "nucleotide_chain",
        ["-m", "nucleotide_chain", "--pdb", fixture("1BNA"), "-c", "A"],
        lambda: biotsr.analyze(fixture("1BNA"), method="nucleotide_chain",
                            chain="A", keep_triplets=False),
    ),
    (
        "protein_protein_cross_all_atom",
        ["-m", "protein_protein_cross", "--pdb", fixture("1BRS"), "-c", "A",
         "--chain-2", "D", "--contact-cutoff", "4.0"],
        lambda: biotsr.analyze_cross(fixture("1BRS"), method="protein_protein_cross",
                                  chain="A", chain_2="D", contact_cutoff=4.0,
                                  keep_triplets=False),
    ),
    (
        "protein_protein_calpha_cross",
        ["-m", "protein_protein_calpha_cross", "--pdb", fixture("1BRS"),
         "-c", "A", "--chain-2", "D", "--contact-cutoff", "8.0"],
        lambda: biotsr.analyze_cross(fixture("1BRS"),
                                  method="protein_protein_calpha_cross",
                                  chain="A", chain_2="D", contact_cutoff=8.0,
                                  keep_triplets=False),
    ),
    (
        "drug_protein_cross",
        ["-m", "drug_protein_cross", "--pdb", fixture("3PTB"), "-c", "A",
         "--chain-2", "A", "--entity-name", "BEN", "--entity-id", "1",
         "--contact-cutoff", "4.0"],
        lambda: biotsr.analyze_cross(fixture("3PTB"), method="drug_protein_cross",
                                  chain="A", chain_2="A", entity_name="BEN",
                                  entity_id="1", contact_cutoff=4.0,
                                  keep_triplets=False),
    ),
    (
        "cross_arbitrary_cutoff",
        ["-m", "protein_protein_cross", "--pdb", fixture("1BRS"), "-c", "A",
         "--chain-2", "D", "--contact-cutoff", "3.7"],
        lambda: biotsr.analyze_cross(fixture("1BRS"), method="protein_protein_cross",
                                  chain="A", chain_2="D", contact_cutoff=3.7,
                                  keep_triplets=False),
    ),
    (
        "cross_size_filtered",
        ["-m", "protein_protein_cross", "--pdb", fixture("1BRS"), "-c", "A",
         "--chain-2", "D", "--contact-cutoff", "4.0", "--size-filter", "12"],
        lambda: biotsr.analyze_cross(fixture("1BRS"), method="protein_protein_cross",
                                  chain="A", chain_2="D", contact_cutoff=4.0,
                                  size_filter=12.0, keep_triplets=False),
    ),
]


def _run_equivalence():
    for tag, flags, direct in EQUIVALENCE_CASES:
        out = workdir("eq_%s" % tag)

        def body(tag=tag, flags=flags, direct=direct, out=out):
            code, stdout, err = cli(
                ["run"] + flags + ["-o", out, "-q", "--no-triplets",
                                   "--if-exists", "overwrite"]
            )
            assert code == EXIT_OK, "exit %d: %s" % (code, err.strip()[:200])
            # `biotsr run` prints the completion record it installed, which
            # sits directly in --output beside the calculation's data files.
            record_path = stdout.strip()
            assert os.path.dirname(record_path) == out, record_path
            assert record_path.endswith(cli_output.RECORD_SUFFIX), record_path
            prefix = os.path.basename(record_path)[: -len(cli_output.RECORD_SUFFIX)]
            expected = direct()

            # The API's own serializer, under its own default prefix: the file
            # names agreeing means the CLI named the request exactly as
            # native.default_prefix names the result from its provenance.
            reference = os.path.join(TMP, "eq_reference_%s" % tag)
            shutil.rmtree(reference, ignore_errors=True)
            os.makedirs(reference)
            native.save(expected, reference, triplets=False)
            assert native.default_prefix(expected) == prefix, (
                native.default_prefix(expected), prefix)

            produced_files = sorted(
                f for f in os.listdir(out) if f != os.path.basename(record_path)
            )
            expected_files = sorted(os.listdir(reference))
            assert produced_files == expected_files, (
                "%r vs %r" % (produced_files, expected_files)
            )
            for name in expected_files:
                assert filecmp.cmp(
                    os.path.join(out, name),
                    os.path.join(reference, name),
                    shallow=False,
                ), "%s differs" % name
            payload = json.load(open(record_path))
            assert payload["complete"] is True and payload["prefix"] == prefix
            assert len(payload["results"]) == 1, len(payload["results"])
            entry = payload["results"][0]
            assert entry["files"] == {"keys": native.keys_filename(prefix),
                                      "metadata": native.metadata_filename(prefix)}
            assert entry["key_space_fingerprint"] == expected.key_space.fingerprint
            assert entry["sampling_fingerprint"] == expected.sampling.fingerprint
            assert entry["counts"]["n_triangles"] == expected.counts.n_triangles
            assert entry["counts"]["n_unique_keys"] == expected.counts.n_unique_keys
            return "%s: %d triangles, %d keys, %d files" % (
                prefix, expected.counts.n_triangles, expected.counts.n_unique_keys,
                len(expected_files),
            )

        check("equivalence", "run/%s" % tag)(body)


def _run_more_equivalence():
    @check("equivalence", "per_entity_sweep_matches_analyze_entities")
    def _():
        expected = biotsr.analyze_entities(
            fixture("1UBQ"), method="amino_acid_atomic", chain="A",
            keep_triplets=False,
        )
        outcome = cli_spec.execute(
            cli_spec.build_spec(
                {"method": "amino_acid_atomic", "pdb": fixture("1UBQ"),
                 "chain": "A", "keep_triplets": False}
            )
        )
        assert sorted(outcome.results) == sorted(expected)
        for name, produced in outcome.results.items():
            assert produced.metadata() == expected[name].metadata(), name
        return "%d residues, metadata identical" % len(expected)

    @check("equivalence", "sse_annotation_changes_no_key_through_the_cli")
    def _():
        plain = cli_spec.execute(cli_spec.build_spec(
            {"method": "protein", "pdb": fixture("1BRS"), "chain": "A",
             "keep_triplets": False}))
        annotated = cli_spec.execute(cli_spec.build_spec(
            {"method": "protein", "pdb": fixture("1BRS"), "chain": "A",
             "sse": True, "keep_triplets": False}))
        a = list(plain.results.values())[0]
        b = list(annotated.results.values())[0]
        assert dict(a.frequencies) == dict(b.frequencies)
        assert a.key_space.fingerprint == b.key_space.fingerprint
        assert a.sampling.fingerprint == b.sampling.fingerprint
        return "%d keys, identical either way" % len(a.frequencies)

    @check("equivalence", "sse_stratification_is_identical")
    def _():
        expected = biotsr.analyze(fixture("1BRS"), method="protein", chain="A",
                               sse=True, keep_triplets=False)
        produced = list(cli_spec.execute(cli_spec.build_spec(
            {"method": "protein", "pdb": fixture("1BRS"), "chain": "A",
             "sse": True, "keep_triplets": False})).results.values())[0]
        assert {k: dict(v.by_category) for k, v in produced.keys.items()} == {
            k: dict(v.by_category) for k, v in expected.keys.items()
        }
        assert produced.annotations["sse"] == expected.annotations["sse"]
        return "%d categories" % len(expected.category_totals())


# ===========================================================================
# contract
# ===========================================================================


def _run_contract():
    @check("contract", "output_layout_is_one_shared_folder_of_prefixed_files")
    def _():
        """Every calculation's files sit directly in --output under one prefix
        (ADR-098): two calculations into one folder give eight files side by
        side, and no per-sample directory."""
        out = workdir("contract_layout")
        code, stdout, err = cli(
            ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
             "-o", out, "-q", "--if-exists", "overwrite"]
        )
        assert code == EXIT_OK, err
        assert stdout.strip() == os.path.join(out, "1UBQ_A_protein.run.json"), stdout
        code, stdout, err = cli(
            ["run", "-m", "aa_ind", "--pdb", fixture("1UBQ"), "-c", "A",
             "--residue-types", "LEU", "-o", out, "-q", "--if-exists", "overwrite"]
        )
        assert code == EXIT_OK, err
        assert stdout.strip() == os.path.join(out, "1UBQ_A_aa_ind_LEU.run.json"), stdout
        names = sorted(os.listdir(out))
        expected = sorted(
            prefix + suffix
            for prefix in ("1UBQ_A_protein", "1UBQ_A_aa_ind_LEU")
            for suffix in (native.KEYS_SUFFIX, native.METADATA_SUFFIX,
                           native.TRIPLETS_SUFFIX, cli_output.RECORD_SUFFIX)
        )
        assert names == expected, names
        assert not any(os.path.isdir(os.path.join(out, n)) for n in names)
        return "2 calculations side by side, 4 files each, no subdirectory"

    @check("contract", "filenames_come_from_the_native_serializer")
    def _():
        result = biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                                keep_triplets=False)
        prefix = native.default_prefix(result)
        assert prefix == "1UBQ_A_protein", prefix
        spec = cli_spec.build_spec(
            {"method": "protein", "pdb": fixture("1UBQ"), "chain": "A"})
        assert cli_spec.output_prefix(spec) == prefix, cli_spec.output_prefix(spec)
        names = (native.keys_filename(prefix), native.metadata_filename(prefix),
                 native.triplets_filename(prefix))
        assert names == ("1UBQ_A_protein.keys", "1UBQ_A_protein.meta.json",
                         "1UBQ_A_protein.triplets"), names
        assert cli_output.record_name(prefix) == "1UBQ_A_protein.run.json"
        out = workdir("contract_layout")
        for name in names:
            assert os.path.isfile(os.path.join(out, name)), name
        return ", ".join(names)

    @check("contract", "output_prefixes_are_deterministic")
    def _():
        request = {"structure_id": "1UBQ", "method": "protein", "chain": "A"}
        first = naming.request_prefix(request)
        assert first == naming.request_prefix(dict(request))
        assert first != naming.request_prefix(dict(request, chain="B"))
        values = {"method": "aa_ind", "pdb": fixture("1UBQ"), "chain": "A",
                  "residue_types": "TYR;ARG"}
        a = cli_spec.output_prefix(cli_spec.build_spec(dict(values)))
        b = cli_spec.output_prefix(cli_spec.build_spec(dict(values)))
        assert a == b == "1UBQ_A_aa_ind_TYR+ARG", (a, b)
        return "%s, %s" % (first, a)

    @check("contract", "output_prefixes_are_portable")
    def _():
        import re

        name = cli_spec.output_prefix(cli_spec.build_spec(
            {"method": "aa_selection", "pdb": fixture("1UBQ"), "chain": "A",
             "residues": "10-12;48;70-71"}))
        assert set(name) <= naming.SAFE_CHARACTERS, name
        assert name == "1UBQ_A_aa_selection_res-10-12+48+70-71", name
        # Where the readable form is inexact, the request fingerprint is
        # appended, so two requests never share a name -- not even on a
        # case-insensitive file system.
        odd = naming.request_prefix(
            {"structure_id": "1 UBQ/x", "method": "protein", "chain": "A"})
        upper = naming.request_prefix(
            {"structure_id": "1UBQ", "method": "protein", "chain": "A"})
        lower = naming.request_prefix(
            {"structure_id": "1UBQ", "method": "protein", "chain": "a"})
        for prefix in (odd, lower):
            assert set(prefix) <= naming.SAFE_CHARACTERS, prefix
            assert re.search(r"_h[0-9a-f]{10}$", prefix), prefix
        assert lower.casefold() != upper.casefold(), (lower, upper)
        return "%s; %s" % (name, lower)

    @check("contract", "a_long_structure_id_survives")
    def _():
        import re

        name = naming.request_prefix(
            {"structure_id": "PDB_0123456789", "method": "protein", "chain": "A"})
        assert "PDB_0123456789" in name, name
        # Past MAX_READABLE the name is shortened and fingerprinted, and two ids
        # that differ only beyond the cut still get different names.
        a = naming.request_prefix(
            {"structure_id": "X" * 120 + "1", "method": "protein", "chain": "A"})
        b = naming.request_prefix(
            {"structure_id": "X" * 120 + "2", "method": "protein", "chain": "A"})
        for long_name in (a, b):
            assert len(long_name) <= naming.MAX_READABLE + 12, len(long_name)
            assert re.search(r"_h[0-9a-f]{10}$", long_name), long_name
        assert a != b
        return "%s; %d-character ids shortened and fingerprinted" % (name, 121)

    @check("contract", "a_second_run_is_refused_by_default")
    def _():
        out = workdir("contract_collision")
        argv = ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
                "-o", out, "-q", "--no-triplets"]
        assert cli(argv)[0] == EXIT_OK
        code, _stdout, err = cli(argv)
        assert code == EXIT_USAGE
        assert "already exists" in err
        return "exit %d" % code

    @check("contract", "overwrite_replaces_owned_files_and_keeps_every_other_file")
    def _():
        """--if-exists overwrite replaces the calculation's own files, removes the
        ones its new version no longer writes (here the triplets), and leaves a
        file it never wrote alone."""
        out = workdir("contract_overwrite")
        base = ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
                "-o", out, "-q"]
        assert cli(base + ["--if-exists", "overwrite"])[0] == EXIT_OK
        assert os.path.isfile(os.path.join(out, "1UBQ_A_protein.triplets"))
        notes = os.path.join(out, "notes.txt")
        open(notes, "w").write("mine")
        assert cli(base + ["--no-triplets", "--if-exists", "overwrite"])[0] == EXIT_OK
        names = sorted(os.listdir(out))
        assert names == ["1UBQ_A_protein.keys", "1UBQ_A_protein.meta.json",
                         "1UBQ_A_protein.run.json", "notes.txt"], names
        assert open(notes).read() == "mine"
        record = json.load(open(os.path.join(out, "1UBQ_A_protein.run.json")))
        assert record["complete"] is True and record["triplets"] is False
        return "stale triplets removed, an unrelated file kept"

    @check("contract", "a_file_biotsr_did_not_write_is_never_replaced")
    def _():
        out = workdir("contract_foreign")
        target = os.path.join(out, "1UBQ_A_protein.keys")
        open(target, "w").write("someone else's")
        for policy in cli_output.IF_EXISTS_CHOICES:
            code, _stdout, err = cli(
                ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
                 "-o", out, "-q", "--no-triplets", "--if-exists", policy])
            assert code == EXIT_USAGE, (policy, code)
            assert "no BioTSR completion record" in err, err[:200]
        assert open(target).read() == "someone else's"
        assert os.listdir(out) == ["1UBQ_A_protein.keys"], os.listdir(out)
        return "refused under all %d --if-exists policies" % len(
            cli_output.IF_EXISTS_CHOICES)

    @check("contract", "skip_leaves_the_existing_output_untouched")
    def _():
        out = workdir("contract_skip")
        cli(["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
             "-o", out, "-q", "--no-triplets", "--if-exists", "overwrite"])
        paths = [os.path.join(out, "1UBQ_A_protein" + suffix)
                 for suffix in (cli_output.RECORD_SUFFIX, native.KEYS_SUFFIX)]
        before = [(os.path.getmtime(p), open(p, "rb").read()) for p in paths]
        code, _stdout, _err = cli(
            ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
             "-o", out, "-q", "--no-triplets", "--if-exists", "skip"]
        )
        assert code == EXIT_OK
        assert [(os.path.getmtime(p), open(p, "rb").read()) for p in paths] == before
        return "record and keys unchanged"

    @check("contract", "writes_are_staged_and_moved_atomically")
    def _():
        source = ast.parse(open(os.path.join(CLI_DIR, "output.py")).read())
        names = {
            n.attr for n in ast.walk(source) if isinstance(n, ast.Attribute)
        }
        assert "replace" in names, "the calculation's files are not moved atomically"
        out = workdir("contract_staging")
        writer = cli_output.CalculationWriter(out, "p", cli_output.IF_EXISTS_FAIL)
        assert os.path.dirname(writer.staging).startswith(
            os.path.join(out, ".biotsr-staging-")), writer.staging
        assert os.path.basename(writer.staging) == "p"
        cli(["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
             "-o", out, "-q", "--no-triplets", "--if-exists", "overwrite"])
        left = [n for n in os.listdir(out) if n.startswith(".")]
        assert left == [], left
        return "staged in .biotsr-staging-<token>/PREFIX, os.replace, nothing left behind"

    @check("contract", "the_completion_record_is_written_last")
    def _():
        """An intent record (complete: false) goes in before any data file, and
        the complete record is written only once every data file is in place."""
        out = workdir("contract_record_last")
        real = cli_output._write_json_atomically
        seen = []

        def spy(root, name, payload, token):
            present = sorted(n for n in os.listdir(root) if not n.startswith("."))
            seen.append((payload.get("complete"), list(payload.get("files_pending") or ()),
                         present, name))
            return real(root, name, payload, token)

        cli_output._write_json_atomically = spy
        try:
            code, stdout, err = cli(
                ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
                 "-o", out, "-q", "--no-triplets"])
        finally:
            cli_output._write_json_atomically = real
        assert code == EXIT_OK, err
        record_name = "1UBQ_A_protein.run.json"
        data = ["1UBQ_A_protein.keys", "1UBQ_A_protein.meta.json"]
        assert [s[0] for s in seen] == [False, True], [s[0] for s in seen]
        assert [s[3] for s in seen] == [record_name, record_name]
        intent, final = seen
        assert intent[1] == data, intent[1]
        assert intent[2] == [], "data files were in place before the intent record"
        assert final[2] == sorted(data + [record_name]), final[2]
        record = json.load(open(os.path.join(out, record_name)))
        assert record["complete"] is True and "files_pending" not in record
        sizes = record["results"][0]["file_bytes"]
        assert sizes == {n: os.path.getsize(os.path.join(out, n)) for n in data}, sizes
        return "intent record, then %d data files, then the complete record" % len(data)

    @check("contract", "an_incomplete_calculation_is_refused_until_recomputed")
    def _():
        """Interrupted between installing the data files and the final record:
        the record says incomplete, and run, analyze and describe all say so
        rather than read it; --if-exists skip recomputes it."""
        out = workdir("contract_incomplete")
        prefix = "1UBQ_A_protein"
        argv = ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
                "-o", out, "-q", "--no-triplets"]
        real = cli_output._write_json_atomically

        def interrupted(root, name, payload, token):
            if payload.get("complete"):
                raise KeyboardInterrupt
            return real(root, name, payload, token)

        cli_output._write_json_atomically = interrupted
        try:
            code, _stdout, _err = cli(argv)
        finally:
            cli_output._write_json_atomically = real
        assert code == 130, code
        assert [n for n in os.listdir(out) if n.startswith(".")] == [], os.listdir(out)
        state, record = cli_output.existing_state(out, prefix)
        assert state == "incomplete", state
        assert record["files_pending"] == [prefix + ".keys", prefix + ".meta.json"]
        assert "interrupted" in cli_output.check_complete(out, prefix)

        code, _stdout, err = cli(argv)
        assert code == EXIT_USAGE and "incomplete" in err, (code, err[:200])
        analysis = os.path.join(TMP, "contract_incomplete_analysis")
        code, _stdout, err = cli(["analyze", "common-keys", "--from-results", out,
                                  "-o", analysis, "-q"])
        assert code == EXIT_USAGE and "not complete" in err, (code, err[:200])
        assert not os.path.exists(os.path.join(analysis, "analysis"))
        code, text, _err = cli(["describe", out])
        assert code == EXIT_OK and prefix in text and "incomplete" in text, text

        code, _stdout, err = cli(argv + ["--if-exists", "skip"])
        assert code == EXIT_OK, err
        assert cli_output.check_complete(out, prefix) is None
        assert json.load(open(os.path.join(out, prefix + ".run.json")))["complete"] is True
        return "exit 130; refused by run, analyze and describe; skip recomputed it"

    @check("contract", "a_failed_run_creates_no_output_directory")
    def _():
        out = os.path.join(TMP, "contract_never")
        code, _stdout, _err = cli(
            ["run", "-m", "protein", "--pdb", "/no/such/file.pdb", "-c", "A",
             "-o", out, "-q"]
        )
        assert code == EXIT_USAGE
        assert not os.path.exists(out), "an empty output root was left behind"
        return "nothing created"

    @check("contract", "nothing_is_written_outside_the_output_root")
    def _():
        work = workdir("contract_cwd")
        out = workdir("contract_cwd_out")
        here = os.getcwd()
        try:
            os.chdir(work)
            cli(["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
                 "-o", out, "-q", "--no-triplets", "--if-exists", "overwrite"])
        finally:
            os.chdir(here)
        assert os.listdir(work) == [], os.listdir(work)
        return "working directory untouched"

    @check("contract", "main_returns_a_code_and_never_raises_for_a_malformed_command")
    def _():
        """Every exit leaves through ``return``, not ``SystemExit``.

        ``argparse``'s own ``--help`` exits 0. ``biotsr analyze`` with no operation
        once routed through it, so a command that did nothing reported success
        and raised out of ``main()`` instead of returning a code an in-process
        caller could read.
        """
        import io as _io

        for argv in ([], ["analyze"]):
            out, err = _io.StringIO(), _io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = cli_main(argv)     # must not raise
            assert code == EXIT_USAGE, (argv, code)
        return "no command and `analyze` with no operation both return 2"

    @check("contract", "exit_codes_are_the_documented_six")
    def _():
        from biotsr.cli.errors import EXIT_CODES

        # 4 = completed, but some results are partial (ADR-106).
        assert sorted(EXIT_CODES) == [0, 1, 2, 3, 4, 130]
        return ", ".join("%d=%s" % (k, v) for k, v in sorted(EXIT_CODES.items()))

    @check("contract", "batch_records_every_input_row")
    def _():
        out = workdir("contract_batch")
        manifest = os.path.join(TMP, "contract_batch.csv")
        open(manifest, "w").write(
            "sample_id,pdb_path,method,chain\n"
            "a,%s,protein,A\n"
            "b,%s,protein,Q\n"
            "c,%s,wobble,A\n" % (fixture("1UBQ"), fixture("1UBQ"), fixture("1UBQ"))
        )
        code, _stdout, _err = cli(
            ["batch", "-M", manifest, "-o", out, "-q", "--no-triplets",
             "--if-exists", "overwrite"]
        )
        assert code == EXIT_PARTIAL, code
        rows = list(csv.DictReader(open(os.path.join(out, "batch_summary.csv"))))
        assert len(rows) == 3
        assert [r["status"] for r in rows] == ["ok", "failed", "failed"]
        # The one calculation that ran is in the shared folder; the failed
        # rows left nothing.
        assert sorted(os.listdir(out)) == [
            "a.keys", "a.meta.json", "a.run.json", "batch_failures.csv",
            "batch_summary.csv"], sorted(os.listdir(out))
        return "3 rows, 1 ok, 2 failed"

    @check("contract", "every_batch_failure_carries_its_reason")
    def _():
        out = workdir("contract_batch")
        failures = list(
            csv.DictReader(open(os.path.join(out, "batch_failures.csv")))
        )
        assert len(failures) == 2
        assert all(r["message"] for r in failures)
        # ``rows`` names the manifest line(s) each request came from.
        assert [r["rows"] for r in failures] == ["3", "4"], [r["rows"] for r in failures]
        assert [r["sample_id"] for r in failures] == ["b", "c"]
        return "; ".join(r["message"][:40] for r in failures)

    @check("contract", "batch_exits_3_when_and_only_when_something_failed")
    def _():
        out = workdir("contract_batch_clean")
        manifest = os.path.join(TMP, "contract_batch_clean.csv")
        open(manifest, "w").write(
            "sample_id,pdb_path,method,chain\na,%s,protein,A\nb,%s,protein,A\n"
            % (fixture("1UBQ"), fixture("1BRS"))
        )
        code, _stdout, _err = cli(
            ["batch", "-M", manifest, "-o", out, "-q", "--no-triplets",
             "--if-exists", "overwrite"]
        )
        assert code == EXIT_OK
        assert not os.path.exists(os.path.join(out, "batch_failures.csv"))
        return "exit 0, no failure report"

    @check("contract", "batch_stop_halts_at_the_first_failure")
    def _():
        out = workdir("contract_batch_stop")
        manifest = os.path.join(TMP, "contract_batch.csv")
        code, _stdout, _err = cli(
            ["batch", "-M", manifest, "-o", out, "-q", "--no-triplets",
             "--on-error", "stop", "--if-exists", "overwrite"]
        )
        assert code == EXIT_USAGE
        rows = list(csv.DictReader(open(os.path.join(out, "batch_summary.csv"))))
        # Row c names an unknown method, which planning sees before anything is
        # computed, so nothing runs; every row is still recorded.
        assert [r["status"] for r in rows] == ["not_started", "not_started", "failed"]
        computed = [n for n in os.listdir(out)
                    if n.endswith((cli_output.RECORD_SUFFIX, native.KEYS_SUFFIX))]
        assert computed == [], computed
        return "stopped before computing; 3 rows recorded"

    @check("contract", "a_batch_row_and_a_run_command_agree_byte_for_byte")
    def _():
        manifest = os.path.join(TMP, "contract_same.csv")
        open(manifest, "w").write(
            "sample_id,pdb_path,method,chain,sse\ns,%s,protein,A,true\n"
            % fixture("1BRS")
        )
        batch_out = workdir("contract_same_batch")
        run_out = workdir("contract_same_run")
        assert cli(["batch", "-M", manifest, "-o", batch_out, "-q",
                    "--no-triplets", "--if-exists", "overwrite"])[0] == EXIT_OK
        assert cli(["run", "-m", "protein", "--pdb", fixture("1BRS"), "-c", "A",
                    "--sse", "-o", run_out, "-q", "--no-triplets",
                    "--sample-id", "s", "--if-exists", "overwrite"])[0] == EXIT_OK
        names = ("s.keys", "s.meta.json")
        for name in names:
            assert filecmp.cmp(os.path.join(batch_out, name),
                               os.path.join(run_out, name), shallow=False), name
        a = json.load(open(os.path.join(batch_out, "s.run.json")))
        b = json.load(open(os.path.join(run_out, "s.run.json")))
        assert a["request"] == b["request"], (a["request"], b["request"])
        assert a["request_sha256"] == b["request_sha256"]
        return "identical files and identical normalised requests"

    @check("contract", "provenance_survives_serialization")
    def _():
        import hashlib

        out = workdir("contract_provenance")
        code, stdout, _err = cli(
            ["run", "-m", "protein", "--pdb", fixture("1BRS"), "-c", "A",
             "--sse", "--mirror", "--size-filter", "20", "-o", out, "-q",
             "--no-triplets", "--label", "barnase", "--if-exists", "overwrite"]
        )
        assert code == EXIT_OK
        payload = json.load(open(stdout.strip()))
        for field in ("method", "profile", "parameters", "results", "biotsr_version",
                      "prefix", "request", "request_sha256", "argv", "command"):
            assert field in payload, field
        assert "output_format" not in payload
        assert payload["schema"] == cli_output.REQUEST_SCHEMA
        assert payload["complete"] is True
        assert payload["prefix"] == "1BRS_A_protein_mirror_size20_sse", payload["prefix"]
        assert sorted(payload["profile"]) == [
            "active", "resolved", "sse_strand_identity"], sorted(payload["profile"])
        assert payload["profile"]["resolved"]
        assert payload["profile"]["sse_strand_identity"] == "sheet_id_and_serial"
        assert isinstance(payload["results"][0]["warnings"], list)
        assert payload["parameters"]["mirror"] is True
        assert payload["parameters"]["size_filter"] == 20.0
        assert payload["label"] == "barnase"
        assert payload["results"][0]["key_space_fingerprint"]
        assert payload["results"][0]["sampling_fingerprint"]
        assert payload["request_sha256"] == cli_output.request_sha256(payload["request"])
        with open(fixture("1BRS"), "rb") as fh:
            assert payload["request"]["structure_sha256"] == hashlib.sha256(
                fh.read()).hexdigest()
        return "method, profile, both fingerprints, strand identity, parameters, request"

    @check("contract", "the_completion_record_records_only_what_the_user_specified")
    def _():
        out = workdir("contract_minimal")
        code, stdout, _err = cli(
            ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
             "-o", out, "-q", "--no-triplets", "--if-exists", "overwrite"]
        )
        assert code == EXIT_OK
        payload = json.load(open(stdout.strip()))
        for absent in ("mirror", "aa_grouping", "sse", "size_filter"):
            assert absent not in payload["parameters"], absent
        assert payload["request"]["modifiers"] == {}, payload["request"]["modifiers"]
        return "no unrequested parameter echoed back"

    @check("contract", "describe_lists_the_calculations_of_a_folder")
    def _():
        out = workdir("contract_minimal")
        code, text, _err = cli(["describe", out])
        assert code == EXIT_OK
        lines = text.splitlines()
        assert len(lines) == 1, lines
        for part in ("1UBQ_A_protein", "protein", "70300 triangles", "complete",
                     "(keys only)"):
            assert part in lines[0], (part, lines[0])
        return lines[0].split()[0]

    @check("contract", "describe_reads_back_a_completion_record_by_its_prefix")
    def _():
        out = workdir("contract_minimal")
        code, text, _err = cli(["describe", os.path.join(out, "1UBQ_A_protein")])
        assert code == EXIT_OK
        assert "protein" in text and "standard_protein_calpha_v1" in text
        assert "1UBQ_A_protein" in text and "complete         True" in text
        return "%d lines" % len(text.splitlines())

    @check("contract", "describe_reads_back_a_native_sidecar")
    def _():
        sidecar = os.path.join(workdir("contract_minimal"),
                               native.metadata_filename("1UBQ_A_protein"))
        code, text, _err = cli(["describe", sidecar])
        assert code == EXIT_OK
        assert "key space" in text and "sampling" in text
        return os.path.basename(sidecar)

    @check("contract", "the_manifest_is_long_not_wide")
    def _():
        manifest = os.path.join(TMP, "wide.csv")
        open(manifest, "w").write(
            "protein,aa_chain,drug_chain,method\n1UBQ,A,A,protein\n"
        )
        try:
            read_manifest(manifest)
        except CLIError as exc:
            assert "aa_chain" in exc.message and "drug_chain" in exc.message
            return "a wide row is refused, naming both columns"
        raise AssertionError("a wide row was accepted")

    @check("contract", "a_blank_manifest_cell_is_not_a_default")
    def _():
        manifest = os.path.join(TMP, "blank.csv")
        open(manifest, "w").write(
            "pdb_id,method,chain,mirror,size_filter\n1UBQ,protein,A,,\n"
        )
        row = read_manifest(manifest)[0]
        assert "mirror" not in row and "size_filter" not in row
        return "absent, not False"

    @check("contract", "an_unknown_manifest_column_is_refused")
    def _():
        manifest = os.path.join(TMP, "unknown.csv")
        open(manifest, "w").write("pdb_id,method,chain,wobble\n1UBQ,protein,A,x\n")
        try:
            read_manifest(manifest)
        except CLIError as exc:
            assert "wobble" in exc.message
            return "named, not ignored"
        raise AssertionError("an unknown column was silently dropped")

    @check("contract", "every_manifest_alias_targets_a_real_column")
    def _():
        for alias, canonical in COLUMN_ALIASES.items():
            assert canonical in CANONICAL_COLUMNS, (alias, canonical)
        return "%d aliases, %d canonical columns" % (
            len(COLUMN_ALIASES), len(CANONICAL_COLUMNS)
        )

    @check("contract", "common_keys_matches_the_analysis_api")
    def _():
        from biotsr.analysis import AnalysisSample, common_key_intersection

        manifest = os.path.join(TMP, "analysis.csv")
        open(manifest, "w").write(
            "sample_id,pdb_path,method,chain,label\n"
            "ubq,%s,protein,A,ubiquitin\n"
            "brs_a,%s,protein,A,barnase\n"
            "brs_b,%s,protein,B,barnase\n"
            % (fixture("1UBQ"), fixture("1BRS"), fixture("1BRS"))
        )
        out = workdir("contract_analysis")
        code, _stdout, err = cli(
            ["analyze", "common-keys", "-M", manifest, "-o", out, "-q",
             "--if-exists", "overwrite"]
        )
        assert code == EXIT_OK, err.strip()[:200]
        payload = json.load(
            open(os.path.join(out, "analysis", "common_keys", "all_samples.json"))
        )
        samples = [
            AnalysisSample.from_result(
                biotsr.analyze(path, method="protein", chain=chain,
                            keep_triplets=False),
                sample_id=sid,
            )
            for sid, path, chain in (
                ("ubq", fixture("1UBQ"), "A"),
                ("brs_a", fixture("1BRS"), "A"),
                ("brs_b", fixture("1BRS"), "B"),
            )
        ]
        expected = common_key_intersection(samples)
        assert payload["keys"] == list(expected.keys)
        assert payload["n_common_keys"] == expected.n_common_keys
        return "%d common keys over 3 samples" % expected.n_common_keys

    @check("contract", "common_keys_groups_by_label")
    def _():
        manifest = os.path.join(TMP, "analysis.csv")
        out = workdir("contract_analysis_grouped")
        code, _stdout, _err = cli(
            ["analyze", "common-keys", "-M", manifest, "-o", out, "-q",
             "--group-by", "label", "--if-exists", "overwrite"]
        )
        assert code == EXIT_OK
        root = os.path.join(out, "analysis", "common_keys")
        assert sorted(os.listdir(root)) == [
            "barnase.csv", "barnase.json", "ubiquitin.csv", "ubiquitin.json"
        ]
        return "2 groups"

    @check("contract", "common_keys_refuses_to_answer_over_a_subset")
    def _():
        manifest = os.path.join(TMP, "analysis_bad.csv")
        open(manifest, "w").write(
            "sample_id,pdb_path,method,chain\n"
            "ubq,%s,protein,A\nbad,%s,protein,Q\n"
            % (fixture("1UBQ"), fixture("1UBQ"))
        )
        out = workdir("contract_analysis_bad")
        code, _stdout, _err = cli(
            ["analyze", "common-keys", "-M", manifest, "-o", out, "-q",
             "--if-exists", "overwrite"]
        )
        assert code != EXIT_OK
        assert not os.path.exists(os.path.join(out, "analysis"))
        return "exit %d, no analysis written" % code

    @check("contract", "a_wide_sheet_dry_run_lists_every_request")
    def _():
        """`--sheet tsr` with methods chosen on the command line (ADR-097): the
        preview lists every request with its source rows and prefix, and writes
        nothing."""
        out = os.path.join(TMP, "contract_sheet")
        code, text, err = cli(
            ["batch", "--manifest", WIDE_SHEET, "--sheet", "tsr",
             "--method", "carbohydrate", "--method", "aa_ind",
             "--inapplicable", "skip", "--pdb-dir", FIXTURES, "--output", out,
             "--dry-run"])
        assert code == EXIT_OK, err[:200]
        assert not os.path.exists(out), "a dry run wrote %s" % out
        lines = text.splitlines()
        assert any("prefix / outcome" in l for l in lines), "no request table"
        requests = [l for l in lines if l.rstrip().endswith("]")
                    and ("[would run]" in l or "[inapplicable" in l)]
        assert "(8 source rows, 13 requests)" in text + err, (text + err)[:300]
        assert len(requests) == 13, len(requests)
        assert any("7LF9_D_carbohydrate_NAG_5" in l and "[would run]" in l
                   for l in lines)
        assert any("7LF9_D_carbohydrate_MAN_allcopies" in l for l in lines)
        # One protein-side calculation serves three glycan rows.
        assert any("4;5;6" in l and "7LF9_B_aa_ind_TYR" in l for l in lines)
        assert any("inapplicable: the row has no carb_name / carb_id" in l
                   for l in lines)
        # Methods come only from --method: a sheet without one is refused.
        code, _text, err = cli(["batch", "--manifest", WIDE_SHEET, "--sheet", "tsr",
                                "--pdb-dir", FIXTURES, "--output", out, "--dry-run"])
        assert code == EXIT_USAGE and "--method" in err, (code, err[:200])
        assert not os.path.exists(out)
        return "13 requests from 8 rows listed, nothing written; no --method refused"

    def _saved(folder, method, pdb, extra):
        code, _stdout, err = cli(
            ["run", "-m", method, "--pdb", fixture(pdb), "-c", "A", "-o", folder,
             "-q", "--no-triplets", "--if-exists", "overwrite"] + extra)
        assert code == EXIT_OK, err[:200]

    @check("contract", "common_keys_compares_lys_only_with_tyr_only_without_an_opt_in")
    def _():
        """What was selected never blocks a comparison (ADR-094): LYS-only and
        TYR-only samples of two structures compare with no flag and no warning."""
        root = workdir("contract_rule")
        plain = os.path.join(root, "plain")
        _saved(plain, "aa_ind", "1UBQ", ["--residue-types", "LYS"])
        _saved(plain, "aa_ind", "1BRS", ["--residue-types", "TYR"])
        out = os.path.join(TMP, "contract_rule_same")
        code, _stdout, err = cli(["analyze", "common-keys", "--from-results", plain,
                                  "-o", out, "-q"])
        assert code == EXIT_OK, err[:300]
        report = json.load(open(os.path.join(out, "analysis", "common_keys",
                                             "all_samples.json")))
        assert report["warnings"] == [], report["warnings"]
        comparison = report["comparison"]
        assert comparison["processing_uniform"] is True
        assert comparison["mixed_processing_allowed"] is False
        assert {s["sample_id"]: s["selection"][0]["arguments"]
                for s in report["samples"]} == {
            "1UBQ_A_aa_ind_LYS": {"residue_types": ["LYS"]},
            "1BRS_A_aa_ind_TYR": {"residue_types": ["TYR"]},
        }, [s["sample_id"] for s in report["samples"]]
        return "%d common keys, no warning" % report["n_common_keys"]

    @check("contract", "a_processing_difference_is_refused_unless_acknowledged")
    def _():
        root = workdir("contract_rule")
        _saved(os.path.join(root, "sized"), "aa_ind", "1BRS",
               ["--residue-types", "LYS", "--size-filter", "12"])
        code, _stdout, err = cli(["analyze", "common-keys", "--from-results", root,
                                  "-o", os.path.join(TMP, "contract_rule_mixed"), "-q"])
        assert code == EXIT_USAGE, code
        assert "size_filter" in err and "--allow-mixed-processing" in err, err[:300]
        assert not os.path.exists(os.path.join(TMP, "contract_rule_mixed", "analysis"))
        # The rc3 spelling is refused, naming the new one.
        code, _stdout, err = cli(["analyze", "common-keys", "--from-results", root,
                                  "-o", os.path.join(TMP, "contract_rule_old"), "-q",
                                  "--allow-mixed-sampling"])
        assert code == EXIT_USAGE and "renamed --allow-mixed-processing" in err, err[:300]
        out = os.path.join(TMP, "contract_rule_allowed")
        code, _stdout, err = cli(["analyze", "common-keys", "--from-results", root,
                                  "-o", out, "-q", "--allow-mixed-processing"])
        assert code == EXIT_OK, err[:300]
        report = json.load(open(os.path.join(out, "analysis", "common_keys",
                                             "all_samples.json")))
        comparison = report["comparison"]
        assert comparison["mixed_processing_allowed"] is True
        assert comparison["n_processing_policies"] == 2
        assert comparison["processing_differences"] == ["size_filter"]
        assert any("size_filter" in w for w in report["warnings"]), report["warnings"]
        return "refused naming size_filter; old flag refused; acknowledged and recorded"

    @check("contract", "the_key_space_refusal_survives_allow_mixed_processing")
    def _():
        folder = workdir("contract_keyspace")
        _saved(folder, "protein", "1UBQ", [])
        _saved(folder, "protein", "1UBQ", ["--mirror"])
        out = os.path.join(TMP, "contract_keyspace_out")
        code, _stdout, err = cli(["analyze", "common-keys", "--from-results", folder,
                                  "-o", out, "-q", "--allow-mixed-processing"])
        assert code == EXIT_USAGE, code
        assert "key space" in err and "mirror_convention" in err, err[:300]
        assert "No option relaxes this" in err, err[:300]
        assert not os.path.exists(os.path.join(out, "analysis"))
        return "plain vs mirror refused even with --allow-mixed-processing"

    @check("contract", "the_analysis_refuses_a_keys_file_no_record_lists")
    def _():
        folder = workdir("contract_stray")
        _saved(folder, "protein", "1UBQ", [])
        _saved(folder, "protein", "1BRS", [])
        shutil.copy(os.path.join(folder, "1UBQ_A_protein.keys"),
                    os.path.join(folder, "stray.keys"))
        out = os.path.join(TMP, "contract_stray_out")
        code, _stdout, err = cli(["analyze", "common-keys", "--from-results", folder,
                                  "-o", out, "-q"])
        assert code == EXIT_USAGE, code
        assert "stray.keys" in err and "no completion record lists it" in err, err[:300]
        # A listed file whose size changed is refused too.
        os.remove(os.path.join(folder, "stray.keys"))
        with open(os.path.join(folder, "1BRS_A_protein.keys"), "a") as fh:
            fh.write("\n")
        code, _stdout, err = cli(["analyze", "common-keys", "--from-results", folder,
                                  "-o", out, "-q"])
        assert code == EXIT_USAGE and "not complete" in err, (code, err[:300])
        assert not os.path.exists(os.path.join(out, "analysis"))
        return "an unlisted keys file and a modified one both refused"

    @check("contract", "the_old_per_sample_layout_is_refused_by_name")
    def _():
        """Outputs of 1.0.0rc3 and earlier (a directory per sample with run.json,
        files like 1UBQ_A.keys_29_35) are never read as current."""
        root = workdir("contract_old_layout")
        sample = os.path.join(root, "1UBQ_A_protein")
        os.makedirs(sample, exist_ok=True)
        result = biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                                keep_triplets=False)
        written = native.save(result, sample, triplets=False)
        os.rename(written["keys"], os.path.join(sample, "1UBQ_A.keys_29_35"))
        os.rename(written["metadata"],
                  os.path.join(sample, "1UBQ_A.keys_29_35.meta.json"))
        open(os.path.join(sample, "run.json"), "w").write("{}")
        out = os.path.join(TMP, "contract_old_layout_out")
        code, _stdout, err = cli(["analyze", "common-keys", "--from-results", root,
                                  "-o", out, "-q"])
        assert code == EXIT_USAGE and "release candidate" in err, (code, err[:300])
        code, _stdout, err = cli(["describe", sample])
        assert code == EXIT_USAGE and "release candidate" in err, (code, err[:300])
        assert not os.path.exists(out)
        return "analyze and describe both refuse it, naming the release candidate"


# ===========================================================================
# discovery
# ===========================================================================


def _run_discovery():
    @check("discovery", "tsr_methods_lists_every_registered_method")
    def _():
        code, text, _err = cli(["methods"])
        assert code == EXIT_OK
        for name in biotsr.list_methods():
            assert name in text, name
        return "%d methods" % len(biotsr.list_methods())

    @check("discovery", "tsr_methods_is_generated_from_the_registry")
    def _():
        code, text, _err = cli(["methods", "--json"])
        records = {r["name"]: r for r in json.loads(text)}
        assert set(records) == set(biotsr.list_methods())
        for name, record in records.items():
            recipe = get_method(name)
            assert record["modifiers"] == sorted(recipe.supported_modifiers)
            assert record["description"] == recipe.description
            assert record["profile"] == recipe.default_profile.name
            assert record["profile_description"] == recipe.default_profile.description
            assert record["target_arguments"] == sorted(accepted_arguments(recipe))
        return "registry fields match one for one"

    @check("discovery", "each_method_describes_its_standard_profile")
    def _():
        code, text, _err = cli(["methods", "protein"])
        assert code == EXIT_OK
        profile = get_method("protein").default_profile.name
        assert profile in text, profile
        assert "legacy_" not in text, "a legacy profile name is still listed"
        return profile

    @check("discovery", "help_exists_for_every_command")
    def _():
        parser = build_parser()
        for argv in (["--help"], ["methods", "--help"], ["run", "--help"],
                     ["batch", "--help"], ["analyze", "--help"],
                     ["analyze", "common-keys", "--help"],
                     ["describe", "--help"]):
            out = io.StringIO()
            try:
                with contextlib.redirect_stdout(out):
                    parser.parse_args(argv)
            except SystemExit as exc:
                assert exc.code == 0, argv
            assert "usage:" in out.getvalue(), argv
            assert len(out.getvalue()) > 200, argv
        return "7 help screens"

    @check("discovery", "help_marks_required_arguments_and_defaults")
    def _():
        parser = build_parser()
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                parser.parse_args(["run", "--help"])
        except SystemExit:
            pass
        text = out.getvalue()
        assert "Required" in text
        assert "default" in text
        assert "standard" in text
        return "required, defaults and the standard profile all stated"

    @check("discovery", "help_uses_no_internal_terminology")
    def _():
        parser = build_parser()
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                parser.parse_args(["run", "--help"])
        except SystemExit:
            pass
        text = out.getvalue().lower()
        for internal in ("keyspace", "samplingprofile", "tsrresult", "recipe",
                         "triplet_sink", "candidate strategy"):
            assert internal not in text, internal
        return "no internal names leak"

    @check("discovery", "every_example_in_the_main_help_parses")
    def _():
        import shlex

        parser = build_parser()
        examples = _examples_from_help(parser)
        assert examples, "no examples found in the help text"
        checked = 0
        for command in examples:
            argv = shlex.split(command)[1:]  # drop the leading "biotsr"
            try:
                parser.parse_args(argv)
            except SystemExit as exc:
                raise AssertionError("%r does not parse (exit %s)" % (command, exc.code))
            checked += 1
        return "%d examples parse" % checked

    @check("discovery", "every_tsr_command_in_the_cli_guide_parses")
    def _():
        """Documentation that has drifted is worse than none.

        Every ``biotsr ...`` line in ``docs/CLI_GUIDE.md`` -- inside fenced blocks,
        in prose, and in the SLURM script -- is handed to the real parser. A flag
        that is renamed or removed breaks this gate, which is the point.
        """
        import re
        import shlex

        guide = os.path.join(WORKSPACE, "docs", "CLI_GUIDE.md")
        assert os.path.isfile(guide), guide
        text = open(guide).read()
        # Join shell line continuations so a wrapped example is one command.
        text = re.sub(r"\\\s*\n\s*", " ", text)
        parser = build_parser()
        # A documented command may be embedded in a pipeline or a $(...)
        # substitution; what has to parse is the biotsr invocation itself.
        # ``biotsr(?=[ \t])`` rather than ``\btsr\b`` so an error message quoted in the
        # guide ("biotsr: error: ...") is not mistaken for a command.
        pattern = re.compile(r"(?:^|\$\()\s*(?:python -m )?biotsr(?=[ \t])(.*)$")
        checked = 0
        in_shell = False
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith("```"):
                in_shell = stripped.startswith("```bash")
                continue
            if not in_shell:
                # Only fenced shell blocks are commands. Prose and captured
                # output legitimately contain the word "biotsr".
                continue
            line = stripped
            if line.startswith("#") or not line:
                continue
            if line.startswith("$ "):
                line = line[2:].strip()
            if line.startswith("SAMPLE="):
                line = line[len("SAMPLE="):].strip()
            match = pattern.search(line)
            if match is None:
                continue
            try:
                argv = shlex.split(match.group(1))
            except ValueError:
                continue
            for stop in ("|", ")", "#"):
                if stop in argv:
                    argv = argv[: argv.index(stop)]
            argv = [a for a in argv if not a.endswith(")") or a == ")"]
            if not argv:
                continue
            try:
                with contextlib.redirect_stdout(io.StringIO()), \
                        contextlib.redirect_stderr(io.StringIO()):
                    parser.parse_args(argv)
            except SystemExit as exc:
                assert exc.code == 0, "%r does not parse" % line
            checked += 1
        assert checked >= 25, "only %d commands found in the guide" % checked
        return "%d documented commands parse" % checked

    @check("discovery", "the_cli_guide_documents_every_command_and_exit_code")
    def _():
        guide = open(os.path.join(WORKSPACE, "docs", "CLI_GUIDE.md")).read()
        for command in ("biotsr run", "biotsr batch", "biotsr methods", "biotsr describe",
                        "biotsr analyze common-keys", "python -m biotsr"):
            assert command in guide, command
        for topic in ("exit code", "manifest", "standard",
                      "SLURM", "--if-exists", "--group-by"):
            assert topic in guide, topic
        return "all commands, exit codes, manifest, the standard profile and HPC covered"

    @check("discovery", "exit_codes_are_documented_in_the_help")
    def _():
        parser = build_parser()
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                parser.parse_args(["--help"])
        except SystemExit:
            pass
        text = out.getvalue()
        for code in ("0", "1", "2", "3"):
            assert code in text
        assert "exit codes" in text
        return "all four documented"


def _examples_from_help(parser):
    """Every ``biotsr ...`` line from every command's epilog, unwrapped."""
    import re

    def subparsers(p):
        for action in p._actions:
            choices = getattr(action, "choices", None)
            if isinstance(choices, dict):   # a subparsers action, not a choice flag
                for sub in choices.values():
                    yield sub

    texts = [parser.epilog or ""]
    for sub in subparsers(parser):
        texts.append(sub.epilog or "")
        for deeper in subparsers(sub):
            texts.append(deeper.epilog or "")
    commands = []
    for text in texts:
        joined = re.sub(r"\\\s*\n\s*", " ", text)
        for line in joined.splitlines():
            line = line.strip()
            if line.startswith("biotsr "):
                commands.append(line)
    return commands


# ===========================================================================
# safety
# ===========================================================================


def _run_safety():
    @check("safety", "the_standard_profile_is_the_default_everywhere")
    def _():
        for name in biotsr.list_methods():
            profile = get_method(name).default_profile
            assert profile.name.startswith("standard_"), (name, profile.name)
            assert not profile.is_legacy_compatibility, name
        return "%d methods default to standard" % len(biotsr.list_methods())

    @check("safety", "the_cli_declares_no_scientific_default")
    def _():
        parser = build_parser()
        run = None
        for action in parser._actions:
            if getattr(action, "choices", None) and "run" in action.choices:
                run = action.choices["run"]
        scientific = {"--mirror", "--aa-grouping", "--sse", "--size-filter",
                      "--contact-cutoff", "--contact-filter",
                      "--no-contact-filter", "--residues", "--residue-types",
                      "--atoms", "--entity-names"}
        for action in run._actions:
            for option in action.option_strings:
                if option in scientific:
                    assert action.default is None, option
        return "%d scientific flags, none with a CLI default" % len(scientific)

    @check("safety", "the_cli_holds_no_numeric_scientific_constant")
    def _():
        # 130 is the exit code for an interrupted run (128 + SIGINT).
        allowed = {0, 1, 2, 3, -1, 4, 5, 8, 10, 12, 16, 24, 32, 36, 76, 88, 130}
        # Named, non-scientific integers, each allowed in one file only (the
        # same table as tests/unit/test_cli_architecture.py): in output.py, 127
        # is the ASCII DEL control character an output prefix may not contain
        # and 200 the longest prefix (in bytes) that stays a portable file
        # name; in batch.py, 7, 28 and 34 are column widths of the dry-run
        # preview and 31 the truncation of its selection column.
        allowed_ints = {"output.py": {127, 200}, "batch.py": {7, 28, 31, 34}}
        # The one non-scientific float: the similarity above which a misspelt
        # manifest column gets a "did you mean" suggestion.
        allowed_floats = {"manifest.py": {0.8}}
        offenders = []
        for name in sorted(os.listdir(CLI_DIR)):
            if not name.endswith(".py"):
                continue
            tree = ast.parse(open(os.path.join(CLI_DIR, name)).read())
            docs = _docstring_ids(tree)
            for node in ast.walk(tree):
                if id(node) in docs:
                    continue
                if isinstance(node, ast.Constant):
                    if isinstance(node.value, bool):
                        continue
                    if isinstance(node.value, float):
                        if node.value not in allowed_floats.get(name, ()):
                            offenders.append((name, node.value))
                    elif isinstance(node.value, int) and node.value not in allowed \
                            and node.value not in allowed_ints.get(name, ()):
                        offenders.append((name, node.value))
        assert offenders == [], str(offenders)
        return "no distance, angle, radix or bin count in the CLI"

    @check("safety", "the_cli_imports_no_scientific_internal")
    def _():
        forbidden = {"geometry", "binning", "encoding", "candidates", "cross",
                     "selection", "structure", "sse"}
        offenders = []
        for name in sorted(os.listdir(CLI_DIR)):
            if not name.endswith(".py"):
                continue
            tree = ast.parse(open(os.path.join(CLI_DIR, name)).read())
            for node in ast.walk(tree):
                module = ""
                if isinstance(node, ast.ImportFrom):
                    module = node.module or ""
                elif isinstance(node, ast.Import):
                    module = ",".join(a.name for a in node.names)
                for part in module.replace("biotsr.", "").split("."):
                    if part in forbidden:
                        offenders.append((name, module))
        assert offenders == [], str(offenders)
        return "api, io, profiles and the registry only"

    @check("safety", "the_cli_reimplements_no_serializer_filename")
    def _():
        markers = (".keys_", ".triplets_", ".3Dkeys", "theta29", "dist18")
        offenders = []
        for name in sorted(os.listdir(CLI_DIR)):
            if not name.endswith(".py"):
                continue
            tree = ast.parse(open(os.path.join(CLI_DIR, name)).read())
            docs = _docstring_ids(tree)
            for node in ast.walk(tree):
                if id(node) in docs:
                    continue
                if isinstance(node, ast.Constant) and isinstance(node.value, str):
                    for marker in markers:
                        if marker in node.value:
                            offenders.append((name, node.value))
        assert offenders == [], str(offenders)
        return "every filename comes from biotsr.io"

    @check("safety", "the_cli_offers_no_profile_or_format_choice")
    def _():
        """``--profile`` and ``--format`` no longer exist: every method runs under
        its standard profile and writes the native format, so either flag is an
        unrecognised argument, refused with exit 2 before anything is computed."""
        refused = 0
        for extra in (["--profile", "standard"], ["--format", "native"]):
            out = os.path.join(TMP, "safety_no_choice")
            argv = ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
                    "-o", out, "-q"] + extra
            err = io.StringIO()
            try:
                with contextlib.redirect_stdout(io.StringIO()), \
                        contextlib.redirect_stderr(err):
                    code = cli_main(argv)
            except SystemExit as exc:
                code = exc.code
            assert code == EXIT_USAGE, (extra, code)
            assert extra[0] in err.getvalue(), (extra, err.getvalue()[:120])
            assert not os.path.exists(out), "output written despite %s" % extra[0]
            refused += 1
        return "%d flags refused, nothing written" % refused

    @check("safety", "sse_is_refused_outside_the_two_calpha_methods")
    def _():
        refused = 0
        for name in biotsr.list_methods():
            if "sse" in get_method(name).supported_modifiers:
                continue
            try:
                cli_spec.build_spec(
                    {"method": name, "pdb": fixture("1UBQ"), "chain": "A",
                     "sse": True, "entity_name": "X", "residues": "1-2",
                     "residue_types": "LEU", "atoms": "A_1_CA", "residue": "1",
                     "residue_name": "LEU", "chain_2": "B"}
                )
            except CLIError as exc:
                assert "sse" in exc.message
                refused += 1
        assert refused == 18, refused
        return "18 of 20 methods refuse --sse by name"

    @check("safety", "aa_grouping_is_refused_outside_the_calpha_methods")
    def _():
        refused = 0
        for name in biotsr.list_methods():
            if "aa_grouping" in get_method(name).supported_modifiers:
                continue
            try:
                cli_spec.build_spec(
                    {"method": name, "pdb": fixture("1UBQ"), "chain": "A",
                     "aa_grouping": True}
                )
            except CLIError:
                refused += 1
        assert refused == 18, refused
        return "18 of 20 methods refuse --aa-grouping"

    @check("safety", "a_nonpositive_size_filter_is_refused")
    def _():
        for value in (0.0, -1.0, -25.0):
            try:
                cli_spec.build_spec(
                    {"method": "protein", "pdb": fixture("1UBQ"), "chain": "A",
                     "size_filter": value}
                )
            except CLIError as exc:
                assert "greater than 0" in exc.message
                continue
            raise AssertionError("accepted size_filter=%r" % value)
        return "0 and negatives refused"

    @check("safety", "a_nonpositive_contact_cutoff_is_refused")
    def _():
        for value in (0.0, -3.0):
            try:
                cli_spec.build_spec(
                    {"method": "protein_protein_cross", "pdb": fixture("1BRS"),
                     "chain": "A", "chain_2": "D", "contact_cutoff": value}
                )
            except CLIError as exc:
                assert "greater than 0" in exc.message
                continue
            raise AssertionError("accepted contact_cutoff=%r" % value)
        return "0 and negatives refused"

    @check("safety", "any_positive_contact_cutoff_is_accepted")
    def _():
        for value in (0.1, 2.5, 3.0, 3.7, 4.0, 5.0, 8.0, 15.0, 100.0):
            spec = cli_spec.build_spec(
                {"method": "protein_protein_cross", "pdb": fixture("1BRS"),
                 "chain": "A", "chain_2": "D", "contact_cutoff": value}
            )
            assert spec.cross["contact_cutoff"] == value
        return "9 arbitrary cutoffs accepted"

    @check("safety", "the_contact_cutoff_is_sampling_not_keyspace")
    def _():
        a = biotsr.analyze_cross(
            fixture("1BRS"), method="protein_protein_cross", chain="A",
            chain_2="D", contact_cutoff=4.0, keep_triplets=False,
        )
        b = biotsr.analyze_cross(
            fixture("1BRS"), method="protein_protein_cross", chain="A",
            chain_2="D", contact_cutoff=6.0, keep_triplets=False,
        )
        assert a.key_space.fingerprint == b.key_space.fingerprint
        assert a.sampling.fingerprint != b.sampling.fingerprint
        return "same KeySpace, different SamplingProfile"

    @check("safety", "the_size_filter_is_distinct_from_the_contact_cutoff")
    def _():
        spec = cli_spec.build_spec(
            {"method": "protein_protein_cross", "pdb": fixture("1BRS"),
             "chain": "A", "chain_2": "D", "contact_cutoff": 4.0,
             "size_filter": 12.0}
        )
        assert spec.cross["contact_cutoff"] == 4.0
        assert spec.modifiers["size_filter"] == 12.0
        assert "size_filter" not in spec.cross
        return "routed to different arguments"

    @check("safety", "an_empty_selection_is_refused_rather_than_written")
    def _():
        out = workdir("safety_empty")
        code, _stdout, err = cli(
            ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "Q",
             "-o", out, "-q"]
        )
        assert code == EXIT_USAGE
        assert "no atoms" in err and "chains in this structure" in err
        written = [n for n in os.listdir(out) if n.startswith("1UBQ_Q_protein")]
        assert written == [], written
        return "refused, and the available chains named"

    @check("safety", "no_user_error_prints_a_traceback")
    def _():
        cases = [
            ["run", "-m", "protein", "--pdb", "/no/such.pdb", "-c", "A"],
            ["run", "-m", "wobble", "--pdb", fixture("1UBQ"), "-c", "A"],
            ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "Z"],
            ["run", "-m", "drug", "--pdb", fixture("3PTB"), "-c", "A"],
            ["run", "-m", "drug", "--pdb", fixture("3PTB"), "-c", "A",
             "--entity-name", "XXX", "--entity-id", "1"],
            ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
             "--entity-name", "BEN"],
            ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
             "--size-filter", "-1"],
            ["run", "-m", "aa_selection", "--pdb", fixture("1UBQ"), "-c", "A",
             "--residues", "not-a-selection"],
            ["run", "-m", "amino_acid_atom_selection", "--pdb", fixture("1UBQ"),
             "-c", "A", "--atoms", "NOPE_999_ZZ"],
            ["batch", "-M", "/no/such.csv"],
            ["batch", "-M", WIDE_SHEET, "--sheet", "tsr", "--pdb-dir", FIXTURES],
            ["analyze", "common-keys", "--from-results", "/no/such/dir"],
            ["analyze", "common-keys", "--from-results", FIXTURES,
             "--allow-mixed-sampling"],
            ["describe", "/no/such/path"],
        ]
        out = workdir("safety_errors")
        for argv in cases:
            code, _stdout, err = cli(argv + (["-o", out] if argv[0] != "describe" else []))
            assert code == EXIT_USAGE, (argv, code)
            assert "Traceback" not in err, argv
            assert err.startswith("biotsr: error:"), (argv, err[:80])
        return "%d error paths, all one line" % len(cases)

    @check("safety", "no_deferred_candidate_strategy_is_registered")
    def _():
        names = set(biotsr.list_methods())
        assert names == RELEASED_METHODS, sorted(names ^ RELEASED_METHODS)
        code, listed, _err = cli(["methods", "--json"])
        assert {r["name"] for r in json.loads(listed)} == RELEASED_METHODS
        return "20 methods, none deferred"

    @check("safety", "no_atomic_or_cross_sse_is_reachable")
    def _():
        code, text, _err = cli(["methods", "--json"])
        records = json.loads(text)
        with_sse = {r["name"] for r in records if "sse" in r["modifiers"]}
        assert with_sse == {"protein", "aa_selection"}, with_sse
        return "sse on protein and aa_selection only"

    @check("safety", "the_cli_adds_no_dependency")
    def _():
        """The CLI imports only what the package's own allow-list permits.

        ``STDLIB_ALLOWED`` is **read out of the hygiene test** rather than
        restated here. An earlier version of this gate carried its own copy, the
        two drifted, and this gate passed while
        ``test_core_imports_nothing_outside_the_standard_library`` failed on the
        same files -- which is the whole failure mode a second list creates.
        One list, one place, parsed at run time.
        """
        base = set(_hygiene_set("STDLIB_ALLOWED")) | {"contextlib", "biotsr"}
        # The hygiene test's per-module carve-outs, read from the same file:
        # joblib, the parallel backend (ADR-086), in parallel.py only; socket,
        # for the host name in a staging token, in output.py only.
        extra = {"parallel.py": _hygiene_set("PARALLEL_ALLOWED"),
                 "output.py": _hygiene_set("CLI_OUTPUT_ALLOWED")}
        offenders = []
        for name in sorted(os.listdir(CLI_DIR)):
            if not name.endswith(".py"):
                continue
            allowed = base | extra.get(name, set())
            tree = ast.parse(open(os.path.join(CLI_DIR, name)).read())
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        top = alias.name.split(".")[0]
                        if top not in allowed:
                            offenders.append((name, alias.name))
                elif isinstance(node, ast.ImportFrom) and not node.level:
                    top = (node.module or "").split(".")[0]
                    if top not in allowed:
                        offenders.append((name, node.module))
        assert offenders == [], str(offenders)
        return "standard library only, joblib in parallel.py, socket in output.py"

    @check("safety", "the_package_declares_only_joblib_as_a_dependency")
    def _():
        text = open(os.path.join(WORKSPACE, "pyproject.toml")).read()
        assert '\ndependencies = ["joblib>=1.3"]\n' in text
        assert 'biotsr = "biotsr.cli:main"' in text
        return "dependencies = joblib>=1.3 (ADR-086), entry point declared"


def _hygiene_set(wanted):
    """A literal set (``STDLIB_ALLOWED`` and the per-module carve-outs) from
    ``tests/unit/test_candidates_and_hygiene.py``.

    Parsed rather than imported: the test module is not importable without
    pytest's fixtures, and the set is a plain literal, so reading it out of the
    AST is both sufficient and immune to anything else in the file.
    """
    path = os.path.join(
        WORKSPACE, "tests", "unit",
        "test_candidates_and_hygiene.py",
    )
    tree = ast.parse(open(path).read())
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == wanted:
                    return set(ast.literal_eval(node.value))
    raise AssertionError("%s not found in %s" % (wanted, path))


def _docstring_ids(tree):
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            if ast.get_docstring(node, clean=False) is not None and node.body:
                first = node.body[0]
                if isinstance(first, ast.Expr):
                    out.add(id(first.value))
    return out


# ===========================================================================


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", dest="json_path", default="")
    args = parser.parse_args()

    try:
        _run_equivalence()
        _run_more_equivalence()
        _run_contract()
        _run_discovery()
        _run_safety()
    finally:
        shutil.rmtree(TMP, ignore_errors=True)

    by_category = {}
    for record in RESULTS:
        bucket = by_category.setdefault(
            record["category"], {"total": 0, "failed": 0}
        )
        bucket["total"] += 1
        bucket["failed"] += not record["passed"]

    current = None
    for record in RESULTS:
        if record["category"] != current:
            current = record["category"]
            print("\n== %s ==" % current)
        print(
            "  %s %-58s %s"
            % ("PASS" if record["passed"] else "FAIL", record["name"],
               record["detail"])
        )

    failed = sum(b["failed"] for b in by_category.values())
    print()
    print("%d gates, %d failed  %r" % (len(RESULTS), failed, by_category))

    payload = {
        "n_gates": len(RESULTS),
        "n_failed": failed,
        "by_category": by_category,
        "gates": RESULTS,
    }
    path = args.json_path or os.path.join(EVIDENCE, "GATE_RESULTS.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print("-> %s" % os.path.relpath(path, WORKSPACE))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
