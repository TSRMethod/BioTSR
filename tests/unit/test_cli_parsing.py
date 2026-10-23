"""Argument parsing, validation, resolution and exit codes.

These tests never run a calculation.  What they check is the layer in front of
one: that a request is understood, that a request that cannot be honoured is
refused before anything is written, and that the refusal says what to do instead.
"""

from __future__ import annotations

import io

import pytest

from biotsr.cli import EXIT_OK, EXIT_USAGE, CLIError, main
from biotsr.cli.parser import build_parser
from biotsr.cli.spec import build_spec, output_prefix
from biotsr.io.naming import request_prefix


def default_sample_id(method, structure_id, chain, chain_2=None, target=None):
    """The automatic prefix of a request with no selection beyond ``target``."""
    selection = {"entity_name": target} if target else {}
    return request_prefix({"structure_id": structure_id, "method": method,
                           "chain": chain, "chain_2": chain_2,
                           "selection": selection})


def run(argv):
    """``(exit_code, stdout, stderr)`` for one invocation."""
    out, err = io.StringIO(), io.StringIO()
    import contextlib

    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main(argv)
    return code, out.getvalue(), err.getvalue()


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def test_flags_are_the_api_keywords_with_dashes():
    """One vocabulary, not two: ``size_filter`` is ``--size-filter``, always."""
    parser = build_parser()
    args = parser.parse_args(
        ["run", "-m", "protein", "--pdb", "x.pdb", "-c", "A", "-o", "out",
         "--size-filter", "25.0", "--aa-grouping", "--mirror", "--sse"]
    )
    assert args.size_filter == 25.0
    assert args.aa_grouping is True
    assert args.mirror is True
    assert args.sse is True


def test_unset_modifiers_stay_unset():
    parser = build_parser()
    args = parser.parse_args(
        ["run", "-m", "protein", "--pdb", "x.pdb", "-c", "A", "-o", "out"]
    )
    assert args.mirror is None
    assert args.aa_grouping is None
    assert args.sse is None
    assert args.size_filter is None
    spec = build_spec({"method": "protein", "pdb": "x.pdb", "chain": "A"})
    assert spec.modifiers == {}
    assert "mirror" not in spec.user_parameters()


def test_output_is_required():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["run", "-m", "protein", "--pdb", "x.pdb", "-c", "A"])


def test_defaults_are_the_documented_ones():
    parser = build_parser()
    args = parser.parse_args(
        ["run", "-m", "protein", "--pdb", "x.pdb", "-c", "A", "-o", "out"]
    )
    assert args.if_exists == "fail"
    assert not hasattr(args, "format")
    assert not hasattr(args, "profile")
    args = parser.parse_args(["batch", "-M", "m.csv", "-o", "out"])
    assert args.on_error == "continue"


def test_analysis_defaults_to_stopping_on_the_first_failure():
    """A key intersection over a silently reduced sample set is a wrong answer."""
    parser = build_parser()
    args = parser.parse_args(
        ["analyze", "common-keys", "-M", "m.csv", "-o", "out"]
    )
    assert args.on_error == "stop"


# ---------------------------------------------------------------------------
# Method resolution
# ---------------------------------------------------------------------------


def test_method_is_resolved_from_the_registry():
    spec = build_spec({"method": "protein", "pdb": "x.pdb", "chain": "A"})
    assert spec.method == "protein"
    assert not spec.is_cross
    spec = build_spec(
        {"method": "protein_protein_cross", "pdb": "x.pdb", "chain": "A",
         "chain_2": "B"}
    )
    assert spec.is_cross


def test_an_unknown_method_names_the_known_ones():
    with pytest.raises(CLIError) as excinfo:
        build_spec({"method": "protien", "pdb": "x.pdb", "chain": "A"})
    assert "unknown method" in excinfo.value.message
    assert "protein" in excinfo.value.hint


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_a_modifier_the_method_does_not_support_is_refused_by_name():
    with pytest.raises(CLIError) as excinfo:
        build_spec(
            {"method": "drug", "pdb": "x.pdb", "chain": "A",
             "entity_name": "BEN", "sse": True}
        )
    assert "does not support 'sse'" in excinfo.value.message
    assert "--mirror" in excinfo.value.hint


def test_aa_grouping_is_refused_outside_the_calpha_methods():
    with pytest.raises(CLIError):
        build_spec(
            {"method": "nucleotide_chain", "pdb": "x.pdb", "chain": "A",
             "aa_grouping": True}
        )


def test_a_target_argument_the_method_never_reads_is_refused():
    with pytest.raises(CLIError) as excinfo:
        build_spec(
            {"method": "protein", "pdb": "x.pdb", "chain": "A",
             "entity_name": "BEN"}
        )
    assert "--entity-name" in excinfo.value.message


def test_a_missing_required_target_names_an_example():
    with pytest.raises(CLIError) as excinfo:
        build_spec({"method": "drug", "pdb": "x.pdb", "chain": "A"})
    assert "--entity-name" in excinfo.value.message
    assert "ESW" in excinfo.value.hint


def test_a_chain_is_always_required():
    with pytest.raises(CLIError) as excinfo:
        build_spec({"method": "protein", "pdb": "x.pdb"})
    assert "chain" in excinfo.value.message


def test_a_nonpositive_size_filter_is_refused():
    for value in (0.0, -1.0):
        with pytest.raises(CLIError) as excinfo:
            build_spec(
                {"method": "protein", "pdb": "x.pdb", "chain": "A",
                 "size_filter": value}
            )
        assert "greater than 0" in excinfo.value.message


def test_a_nonpositive_contact_cutoff_is_refused_and_any_positive_one_is_not():
    with pytest.raises(CLIError):
        build_spec(
            {"method": "protein_protein_cross", "pdb": "x.pdb", "chain": "A",
             "chain_2": "B", "contact_cutoff": 0.0}
        )
    for value in (0.5, 3.0, 5.0, 7.5, 12.0, 100.0):
        spec = build_spec(
            {"method": "protein_protein_cross", "pdb": "x.pdb", "chain": "A",
             "chain_2": "B", "contact_cutoff": value}
        )
        assert spec.cross["contact_cutoff"] == value


def test_the_contact_cutoff_is_sampling_and_the_size_filter_is_a_modifier():
    """They are different parameters and the CLI keeps them apart.

    ``--contact-cutoff`` decides which atoms enter the pool; ``--size-filter``
    decides which triangles survive.  An atom admitted at 3 A can still form a
    20 A triangle, and every implementation keeps it.
    """
    spec = build_spec(
        {"method": "protein_protein_cross", "pdb": "x.pdb", "chain": "A",
         "chain_2": "B", "contact_cutoff": 3.0, "size_filter": 20.0}
    )
    assert spec.cross["contact_cutoff"] == 3.0
    assert spec.modifiers["size_filter"] == 20.0
    assert "size_filter" not in spec.cross


def test_no_contact_filter_and_a_cutoff_contradict_each_other():
    with pytest.raises(CLIError) as excinfo:
        build_spec(
            {"method": "protein_protein_cross", "pdb": "x.pdb", "chain": "A",
             "chain_2": "B", "contact_cutoff": 5.0, "no_contact_filter": True}
        )
    assert "contradict" in excinfo.value.message


# ---------------------------------------------------------------------------
# Output naming
# ---------------------------------------------------------------------------


def test_the_default_sample_id_is_deterministic():
    first = default_sample_id("protein", "1UBQ", "A")
    second = default_sample_id("protein", "1UBQ", "A")
    assert first == second == "1UBQ_A_protein"


def test_the_default_sample_id_separates_runs_that_differ():
    names = {
        default_sample_id("protein", "1UBQ", "A"),
        default_sample_id("protein", "1UBQ", "B"),
        default_sample_id("aa_selection", "1UBQ", "A"),
        default_sample_id("protein", "1TUP", "A"),
        default_sample_id("drug", "3PTB", "A", target="BEN"),
        default_sample_id("protein_protein_cross", "1BRS", "A", chain_2="D"),
    }
    assert len(names) == 6


def test_the_sample_id_is_shell_safe():
    name = request_prefix({"structure_id": "1UBQ", "method": "aa_selection",
                           "chain": "A", "selection": {"residues": "10-12;48;70-71"}})
    assert name == "1UBQ_A_aa_selection_res-10-12+48+70-71"
    assert all(c.isalnum() or c in "._+-" for c in name), name
    awkward = request_prefix({"structure_id": "1UBQ", "method": "amino_acid_atom_selection",
                              "chain": "A", "selection": {"atoms": "LYS_48 NZ;$HOME"}})
    assert all(c.isalnum() or c in "._+-" for c in awkward), awkward
    assert " " not in awkward


def test_a_long_identifier_survives():
    """No four-character PDB assumption exists anywhere in this package."""
    name = default_sample_id("protein", "PDB_0123456789", "A")
    assert "PDB_0123456789" in name


# ---------------------------------------------------------------------------
# Exit codes and error presentation
# ---------------------------------------------------------------------------


def test_a_bad_request_exits_two_and_prints_no_traceback(tmp_path):
    code, out, err = run(
        ["run", "-m", "protien", "--pdb", "x.pdb", "-c", "A",
         "-o", str(tmp_path / "out")]
    )
    assert code == EXIT_USAGE
    assert "Traceback" not in err
    assert err.startswith("biotsr: error:")
    assert out == ""


def test_a_missing_file_exits_two(tmp_path):
    code, _out, err = run(
        ["run", "-m", "protein", "--pdb", str(tmp_path / "nope.pdb"), "-c", "A",
         "-o", str(tmp_path / "out")]
    )
    assert code == EXIT_USAGE
    assert "Traceback" not in err


def test_methods_exits_zero_and_writes_to_stdout():
    code, out, err = run(["methods"])
    assert code == EXIT_OK
    assert "protein" in out
    assert err == ""


def test_methods_json_is_machine_readable():
    import json

    code, out, _err = run(["methods", "--json"])
    assert code == EXIT_OK
    records = json.loads(out)
    names = {r["name"] for r in records}
    assert "protein" in names and "drug_protein_cross" in names
    protein = [r for r in records if r["name"] == "protein"][0]
    assert "sse" in protein["modifiers"]
    assert set(protein) == {
        "name", "family", "scope", "description", "modifiers",
        "target_arguments", "profile", "profile_description",
    }
    assert protein["profile"] == "standard_protein_calpha_v1"
    assert protein["profile_description"]


def test_no_command_prints_help_and_exits_two():
    code, _out, _err = run([])
    assert code == EXIT_USAGE


def test_analyze_without_an_operation_exits_two_by_returning():
    """A malformed command, not a request for help.

    ``argparse``'s own ``--help`` exits 0 through ``SystemExit``. Routing
    ``biotsr analyze`` with no operation through it would report success for a
    command that did nothing, and would raise out of ``main()`` instead of
    returning, so an in-process caller could not read a code at all.
    """
    code, out, err = run(["analyze"])
    assert code == EXIT_USAGE
    assert "common-keys" in out or "common-keys" in err
    assert "biotsr: error:" in err


def test_main_returns_a_code_rather_than_raising_for_every_malformed_command():
    for argv in ([], ["analyze"]):
        out, err = io.StringIO(), io.StringIO()
        import contextlib

        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = main(argv)          # must not raise SystemExit
        assert code == EXIT_USAGE, argv


def test_quiet_suppresses_the_summary_but_not_errors(tmp_path):
    code, _out, err = run(
        ["run", "-m", "protien", "--pdb", "x.pdb", "-c", "A",
         "-o", str(tmp_path / "out"), "-q"]
    )
    assert code == EXIT_USAGE
    assert "biotsr: error:" in err


# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "argv",
    [
        ["--help"],
        ["methods", "--help"],
        ["run", "--help"],
        ["batch", "--help"],
        ["analyze", "--help"],
        ["analyze", "common-keys", "--help"],
        ["describe", "--help"],
    ],
)
def test_every_command_has_help(argv):
    parser = build_parser()
    out = io.StringIO()
    import contextlib

    with contextlib.redirect_stdout(out), pytest.raises(SystemExit) as excinfo:
        parser.parse_args(argv)
    assert excinfo.value.code == 0
    text = out.getvalue()
    assert "usage:" in text
    assert len(text) > 200


def test_help_distinguishes_required_from_optional():
    parser = build_parser()
    out = io.StringIO()
    import contextlib

    with contextlib.redirect_stdout(out), pytest.raises(SystemExit):
        parser.parse_args(["run", "--help"])
    text = out.getvalue()
    assert "--output" in text and "Required" in text
    assert "default" in text


def test_help_states_the_standard_profile_and_offers_no_other():
    parser = build_parser()
    out = io.StringIO()
    import contextlib

    with contextlib.redirect_stdout(out), pytest.raises(SystemExit):
        parser.parse_args(["run", "--help"])
    text = out.getvalue()
    assert "standard" in text
    assert "--profile" not in text
    assert "--format" not in text


def test_help_does_not_leak_internal_terminology():
    """A user should not have to know what a selector or a recipe is."""
    parser = build_parser()
    out = io.StringIO()
    import contextlib

    with contextlib.redirect_stdout(out), pytest.raises(SystemExit):
        parser.parse_args(["run", "--help"])
    text = out.getvalue().lower()
    for internal in ("keyspace", "samplingprofile", "tsrresult", "recipe",
                     "candidate strategy", "triplet_sink"):
        assert internal not in text, "help text exposes %r" % internal
