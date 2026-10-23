"""The commands.  Each one orchestrates; none of them computes.

If a line of this file ever looks like science -- a distance, a bin, an ordering,
a label, a category -- it is in the wrong file.  What lives here is: read the
user's request, hand it to ``biotsr.api``, hand the answer to ``biotsr.io``, say what
happened.
"""

from __future__ import annotations

import csv
import json
import os
from typing import Any, Dict, List

from .. import api as _api
from ..domain.selection_record import STRUCTURE_NOTES
from ..methods.arguments import accepted_arguments
from ..methods.registry import get_method
from ..retrieval import (
    MAX_WORKERS,
    FetchInterrupted,
    fetch_structures,
    file_sha256,
    validate_download_settings,
)
from . import output as _output
from .batch import (
    FAILURES_NAME,
    FETCH_REPORT_NAME,
    ON_ERROR_CONTINUE,
    ON_ERROR_STOP,
    STATUS_FAILED,
    STATUS_NO_TRIANGLES,
    PARTIAL_WARNINGS,
    STATUS_PARTIAL,
    STATUS_OK,
    STATUS_SKIPPED,
    SUMMARY_COLUMNS,
    SUMMARY_NAME,
    BatchInterrupted,
    BatchOptions,
    BatchReport,
    _skipped_record,
    compute_into,
    run_batch,
    structure_id_for,
)
from .errors import EXIT_INCOMPLETE, EXIT_INTERRUPTED, EXIT_OK, EXIT_PARTIAL, CLIError
from .manifest import parse_column_options, read_id_column
from .parallel import resolve_jobs
from .spec import build_spec, output_prefix, sweeps_residues

__all__ = [
    "cmd_methods",
    "cmd_run",
    "cmd_batch",
    "cmd_fetch",
    "cmd_analyze",
    "cmd_describe",
    "ON_ERROR_CONTINUE",
    "ON_ERROR_STOP",
    "SUMMARY_COLUMNS",
]


# ---------------------------------------------------------------------------
# biotsr methods
# ---------------------------------------------------------------------------


def _method_record(name: str) -> Dict[str, Any]:
    recipe = get_method(name)
    return {
        "name": recipe.name,
        "family": recipe.family,
        "scope": recipe.scope,
        "description": recipe.description,
        "modifiers": sorted(recipe.supported_modifiers),
        "target_arguments": sorted(accepted_arguments(recipe)),
        "profile": recipe.default_profile.name,
        "profile_description": recipe.default_profile.description,
    }


def cmd_methods(args, console) -> int:
    names = _api.list_methods()
    if args.method:
        if args.method not in names:
            raise CLIError(
                "unknown method %r" % args.method,
                "supported methods: %s" % ", ".join(names),
            )
        names = (args.method,)

    if args.json:
        console.out(
            json.dumps([_method_record(n) for n in names], indent=2, sort_keys=True)
        )
        return EXIT_OK

    if len(names) == 1:
        record = _method_record(names[0])
        console.out(record["name"])
        console.out("  family     %s" % record["family"])
        console.out("  scope      %s" % record["scope"])
        console.out("  " + _wrap(record["description"], "  "))
        console.out("  modifiers  %s" % (", ".join(
            "--" + m.replace("_", "-") for m in record["modifiers"]) or "none"))
        if record["target_arguments"]:
            console.out("  targets    %s" % ", ".join(
                "--" + a.replace("_", "-") for a in record["target_arguments"]))
        console.out("  profile    %s" % record["profile"])
        console.out("      " + _wrap(record["profile_description"], "      "))
        return EXIT_OK

    console.out("%-30s %-8s %s" % ("METHOD", "FAMILY", "MODIFIERS"))
    for name in names:
        recipe = get_method(name)
        console.out(
            "%-30s %-8s %s"
            % (
                name,
                recipe.family,
                ", ".join(sorted(recipe.supported_modifiers)) or "-",
            )
        )
    console.out("")
    console.out(
        "%d methods. `biotsr methods NAME` describes one, including its "
        "arguments and the standard profile it runs under."
        % len(names)
    )
    console.out(
        "Candidate-strategy extensions (for example CA-kNN) are not part of "
        "this release and are absent on purpose."
    )
    return EXIT_OK


def _wrap(text: str, indent: str, width: int = 76) -> str:
    words = (text or "").split()
    lines: List[str] = []
    current: List[str] = []
    for word in words:
        if current and len(" ".join(current + [word])) > width - len(indent):
            lines.append(" ".join(current))
            current = [word]
        else:
            current.append(word)
    if current:
        lines.append(" ".join(current))
    return ("\n" + indent).join(lines)


# ---------------------------------------------------------------------------
# biotsr run
# ---------------------------------------------------------------------------


def _values_from_args(args) -> Dict[str, Any]:
    """The flat mapping ``build_spec`` validates, from an argparse namespace."""
    values = dict(vars(args))
    for drop in (
        "command", "analysis", "func", "quiet", "verbose", "debug", "output",
        "manifest", "pdb_dir", "on_error", "if_exists", "no_triplets",
        "json", "path", "group_by", "version", "column", "keep_column",
        "ignore_column", "jobs", "download_missing", "download_workers",
        "timeout", "retries", "dry_run", "from_results", "allow_mixed_processing",
        "allow_partial",
        "sheet", "map", "map_file", "participant", "inapplicable", "cross_side",
        "renamed_allow_mixed_sampling",
    ):
        values.pop(drop, None)
    return {k: v for k, v in values.items() if v is not None and v is not False}


def _summarise(outcome, console, payload) -> None:
    spec = outcome.spec
    total_triangles = sum(r.counts.n_triangles for r in outcome.results.values())
    total_keys = sum(r.counts.n_unique_keys for r in outcome.results.values())
    console.say("method     %s" % spec.method)
    console.say("profile    %s" % (outcome.base_profile or outcome.profile))
    if spec.modifiers.get("sse"):
        console.say("sse        on, strand identity %s"
                    % payload["profile"]["sse_strand_identity"])
    if spec.modifiers:
        console.say("modifiers  %s" % ", ".join(
            "%s=%s" % (k, v) for k, v in sorted(spec.modifiers.items())))
    if len(outcome.results) > 1:
        console.say("entities   %d" % len(outcome.results))
    console.say("triangles  %d" % total_triangles)
    console.say("keys       %d distinct" % total_keys)
    for warning in outcome.warnings:
        console.warn(warning)
    for name, result in outcome.results.items():
        for warning in result.warnings:
            if _is_selection_warning(warning):
                console.warn("%s%s" % ("%s: " % name if name else "", warning))
            else:
                console.detail("  warning: %s" % warning)


def _is_selection_warning(text: str) -> bool:
    """The warnings a user must see without --verbose: what a selection lacked,
    and what the structure file lacks (ADR-114)."""
    return text.startswith(("no triangles:",) + PARTIAL_WARNINGS + STRUCTURE_NOTES)


def _run_one(spec, args, console, argv=None, command="run"):
    triplets = not getattr(args, "no_triplets", False)
    prefix = output_prefix(spec)
    try:
        sha = file_sha256(spec.pdb_path)
    except OSError as exc:
        raise CLIError("cannot read structure %s: %s"
                       % (spec.pdb_path, exc.strerror or exc))
    request = _output.request_record(spec, sha, structure_id_for(spec),
                                     per_residue=sweeps_residues(spec))
    writer = _output.CalculationWriter(
        args.output, prefix, args.if_exists,
        verify=lambda name: _output.verify_existing(args.output, name, request,
                                                    triplets),
    )
    # Policy first, then the calculation, then the files: a refused request
    # costs nothing, and a failed one leaves nothing behind.
    if writer.check():
        return prefix, None, True
    if args.if_exists == _output.IF_EXISTS_SKIP:
        state, _old = _output.existing_state(args.output, prefix)
        if state == "incomplete":
            writer.if_exists = _output.IF_EXISTS_OVERWRITE
    outcome, payload = compute_into(
        spec, writer, prefix=prefix, triplets=triplets,
        command=command, argv=argv, request=request,
    )
    return prefix, (outcome, payload), False


def cmd_run(args, console) -> int:
    spec = build_spec(dict(_values_from_args(args), keep_triplets=not args.no_triplets))
    prefix, produced, skipped = _run_one(
        spec, args, console, argv=args._argv, command="run"
    )
    if skipped:
        kept = _skipped_record(_output.read_record(args.output, prefix))
        console.say("skipped %s (verified existing output of the same request)"
                    % prefix)
        if kept["kept_partial"]:
            console.warn(kept["message"])
            return EXIT_INCOMPLETE
        return EXIT_OK
    outcome, payload = produced
    _summarise(outcome, console, payload)
    names = [name for entry in payload["results"] for name in entry["files"].values()]
    console.say("output     %s (%d file%s, completion record %s)" % (
        os.path.join(args.output, prefix), len(names) + 1,
        "" if not names else "s", _output.record_name(prefix)))
    console.out(os.path.join(args.output, _output.record_name(prefix)))
    # A partial result is written and usable, but it is not the selection that
    # was asked for, and the exit code says so (ADR-106).
    if any(entry.get("omissions") for entry in payload["results"]) or any(
            w.startswith(PARTIAL_WARNINGS) for w in payload.get("warnings") or ()):
        return EXIT_INCOMPLETE
    return EXIT_OK


# ---------------------------------------------------------------------------
# biotsr batch
# ---------------------------------------------------------------------------


def _batch_options(args, command="batch") -> BatchOptions:
    from .sheet import parse_cross_side_options, parse_map_options, read_map_file

    defaults: Dict[str, Any] = {}
    chain = getattr(args, "chain", None)
    if chain:
        defaults["chain"] = chain
    for name in ("mirror", "aa_grouping", "sse"):
        if getattr(args, name, False):
            defaults[name] = True
    for name in ("size_filter", "contact_cutoff", "contact_filter"):
        value = getattr(args, name, None)
        if value is not None:
            defaults[name] = value
    _check_download_options(args)
    sheet = getattr(args, "sheet", None)
    maps = list(read_map_file(args.map_file)) if getattr(args, "map_file", None) else []
    maps.extend(parse_map_options(getattr(args, "map", None) or ()))
    if maps and not sheet:
        raise CLIError(
            "--map and --map-file describe a wide sheet; pass --sheet tsr (the TSR "
            "dataset-sheet columns) or --sheet custom (only your mappings)",
            "a one-request-per-row manifest maps columns with --column FIELD=HEADER",
        )
    if sheet and chain:
        raise CLIError(
            "--chain does not apply to a sheet: each method reads its own chain "
            "column, and BioTSR never supplies a chain a row leaves blank",
            "map the chain column with --map METHOD.chain=HEADER",
        )
    participant = getattr(args, "participant", None)
    if participant and not sheet:
        raise CLIError("--participant applies to a sheet (--sheet)")
    cross_sides = parse_cross_side_options(getattr(args, "cross_side", None) or ())
    if cross_sides and not sheet:
        raise CLIError(
            "--cross-side chooses among a wide sheet's columns (--sheet tsr)",
            "a one-request-per-row manifest names each argument in its own column")
    return BatchOptions(
        manifest=args.manifest,
        output=args.output,
        pdb_dir=args.pdb_dir or "",
        methods=tuple(getattr(args, "method", None) or ()),
        defaults=defaults,
        columns=parse_column_options(getattr(args, "column", None) or ()),
        keep_columns=tuple(getattr(args, "keep_column", None) or ()),
        ignore_columns=tuple(getattr(args, "ignore_column", None) or ()),
        sheet=sheet,
        maps=tuple(maps),
        cross_sides=cross_sides,
        participant=participant or "1",
        inapplicable=getattr(args, "inapplicable", None) or "fail",
        triplets=not args.no_triplets,
        if_exists=args.if_exists,
        on_error=args.on_error,
        jobs=resolve_jobs(getattr(args, "jobs", None), warn=console_warning(args)),
        download_missing=bool(getattr(args, "download_missing", False)),
        download_workers=args.download_workers,
        timeout=args.timeout,
        retries=args.retries,
        dry_run=bool(getattr(args, "dry_run", False)),
        command=command,
        argv=args._argv,
    )


def _check_download_options(args) -> None:
    """The download settings, refused by flag name when they cannot be honoured.

    ``fetch``, ``batch`` and ``analyze`` validate them identically, through
    :func:`biotsr.retrieval.validate_download_settings`; nothing is replaced by
    a default.
    """
    try:
        validate_download_settings(
            workers=args.download_workers, timeout=args.timeout, retries=args.retries,
        )
    except ValueError as exc:
        text = str(exc)
        for name, option in (("download workers", "--download-workers"),
                             ("timeout", "--timeout"), ("retries", "--retries")):
            if text.startswith(name):
                text = option + text[len(name):]
                break
        raise CLIError(
            text,
            "--download-workers is a whole number from 1 to %d; --timeout is the "
            "number of seconds to wait for a connection or for the next piece of "
            "a response (not for a whole download) and must be greater than 0; "
            "--retries is a whole number of further attempts, 0 or more"
            % MAX_WORKERS,
        )


def console_warning(args):
    """Warnings from option resolution go to stderr, before the console exists."""
    import sys

    def warn(text):
        if not getattr(args, "quiet", False):
            sys.stderr.write("biotsr: warning: %s\n" % text)

    return warn


def _batch(args, console, command="batch") -> BatchReport:
    try:
        return run_batch(_batch_options(args, command), console)
    except BatchInterrupted as exc:
        _report_batch(exc.report, args, console)
        raise


def cmd_batch(args, console) -> int:
    report = _batch(args, console)
    if getattr(args, "dry_run", False):
        return EXIT_PARTIAL if any(
            r["status"] == STATUS_FAILED for r in report.records) else EXIT_OK
    return _report_batch(report, args, console)


def _report_batch(report: BatchReport, args, console) -> int:
    console.say("")
    for path in report.removed_leftovers:
        console.detail("removed    %s (left by a process that no longer exists)"
                       % path)
    line = "requests   %d ok, %d skipped, %d failed" % (
        report.n_ok, report.n_skipped, report.n_failed)
    if report.n_skipped_partial:
        line = line.replace(" skipped,", " skipped (%d of them partial),"
                            % report.n_skipped_partial)
    if report.n_partial:
        line += ", %d partial" % report.n_partial
    if report.n_no_triangles:
        line += ", %d with no triangles" % report.n_no_triangles
        if report.n_no_triangles_partial:
            line += " (%d of them partial)" % report.n_no_triangles_partial
    if report.n_inapplicable:
        line += ", %d inapplicable" % report.n_inapplicable
    if report.n_chain_not_specified:
        line += ", %d chain not specified" % report.n_chain_not_specified
    if report.n_interrupted or report.n_not_started:
        line += ", %d interrupted, %d not started" % (
            report.n_interrupted, report.n_not_started)
    console.say(line)
    console.say("summary    %s" % os.path.join(args.output, SUMMARY_NAME))
    if report.n_failed:
        console.say("failures   %s" % os.path.join(args.output, FAILURES_NAME))
    if report.fetch_report_path:
        console.say("downloads  %s" % report.fetch_report_path)
    if report.pool_error:
        console.say("workers    the pool failed: %s" % report.pool_error)
    if report.interrupted:
        return EXIT_INTERRUPTED
    console.out(args.output)
    incomplete = report.n_failed or report.n_interrupted or report.n_not_started
    if incomplete or report.pool_error:
        return EXIT_PARTIAL
    if report.n_incomplete:
        return EXIT_INCOMPLETE
    return EXIT_OK


# ---------------------------------------------------------------------------
# biotsr fetch
# ---------------------------------------------------------------------------


def cmd_fetch(args, console) -> int:
    ids: List[str] = list(args.ids or ())
    if args.id_file:
        if not os.path.isfile(args.id_file):
            raise CLIError("--id-file not found: %s" % args.id_file)
        with open(args.id_file) as fh:
            for line in fh:
                text = line.split("#", 1)[0].strip()
                ids.extend(t for t in text.replace(",", " ").split() if t)
    if args.manifest:
        columns = parse_column_options(args.column or ())
        unknown = sorted(set(columns) - {"pdb_id"})
        if unknown:
            raise CLIError(
                "fetch reads only the PDB id column; --column %s=... does not "
                "apply here" % unknown[0],
                "use --column pdb_id=HEADER to name the id column",
            )
        ids.extend(pdb_id for _line, pdb_id in
                   read_id_column(args.manifest, columns.get("pdb_id")))
    elif args.column:
        raise CLIError("--column needs --manifest")
    if not ids:
        raise CLIError(
            "no PDB ids given",
            "name them (biotsr fetch 1UBQ 4HHB ...), or pass --id-file FILE or "
            "--manifest CSV",
        )

    _check_download_options(args)
    try:
        report = fetch_structures(
            ids, args.pdb_dir,
            workers=args.download_workers,
            timeout=args.timeout,
            retries=args.retries,
            refresh=args.refresh,
            dry_run=args.dry_run,
            progress=lambda r: console.detail("%-6s %s" % (r.pdb_id, r.status)),
        )
    except FetchInterrupted as exc:
        path = args.report or os.path.join(args.pdb_dir, FETCH_REPORT_NAME)
        exc.report.write_csv(path)
        console.error("interrupted; the report lists every id: %s" % path)
        return EXIT_INTERRUPTED
    counts = report.counts()
    console.say("ids        %d requested, %d distinct" % (
        len(ids), len(report.records)))
    console.say("result     %s" % ", ".join(
        "%d %s" % (n, s) for s, n in sorted(counts.items())))
    for record in report.failures:
        console.error("%s: %s" % (record.pdb_id or record.requested, record.message))
    for record in report.records:
        if not record.failed and record.message.startswith("note:"):
            console.warn("%s: %s" % (record.pdb_id, record.message[len("note: "):]))
    if args.dry_run:
        console.say("dry run    nothing was downloaded or written")
        return EXIT_PARTIAL if any(
            r.status not in ("cached", "would_download") for r in report.records
        ) else EXIT_OK
    path = args.report or os.path.join(args.pdb_dir, FETCH_REPORT_NAME)
    report.write_csv(path)
    console.say("report     %s" % path)
    console.out(args.pdb_dir)
    return EXIT_OK if report.ok else EXIT_PARTIAL


# ---------------------------------------------------------------------------
# biotsr analyze common-keys
# ---------------------------------------------------------------------------

ANALYSIS_DIR = "analysis"


def _candidate_layout(root: str) -> str:
    """A file of the release-candidate layout under ``root``, or ``""``."""
    import re

    pattern = re.compile(r"\.(keys|triplets)_\d+_\d+(\.meta\.json)?$")
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if not d.startswith(".")]
        for name in filenames:
            if pattern.search(name) or name == "run.json":
                return os.path.join(dirpath, name)
    return ""


def _saved_results(root: str, console, methods=()) -> List[Any]:
    """Every completed calculation under ``root``: ``[(sample_id, result, label)]``.

    A calculation is read through its completion record (``PREFIX.run.json``),
    which lists its files with their sizes; a record that is not complete, a
    listed file that is missing or of another size, and a key file that no
    completion record lists are each refused by name, so a partial output is
    never analysed.  Sorted by record path, so the report order does not depend
    on the file system's listing order.  The sample id is the calculation's
    prefix (one per residue for a per-residue sweep).  ``methods``, when given,
    keeps only the calculations of those methods (``--method``); the others
    are still checked for completeness, so a partial output is never hidden.
    """
    from ..io.native import KEYS_SUFFIX
    from ..io.native_read import CANDIDATE_LAYOUT_MESSAGE, load_result

    if not os.path.isdir(root):
        raise CLIError(
            "--from-results %s is not a directory" % root,
            "point it at an --output folder written by `biotsr run` or `biotsr batch`",
        )
    old = _candidate_layout(root)
    if old:
        raise CLIError(CANDIDATE_LAYOUT_MESSAGE % old,
                       "recompute the dataset with this version into a new --output")
    records = _output.find_records(root)
    if not records:
        raise CLIError(
            "no completed calculations under %s" % root,
            "a calculation is complete when its PREFIX%s completion record is; "
            "`biotsr run` and `biotsr batch` write one per calculation"
            % _output.RECORD_SUFFIX,
        )
    listed = set()
    out = []
    for path in records:
        folder = os.path.dirname(path)
        prefix = os.path.basename(path)[: -len(_output.RECORD_SUFFIX)]
        problem = _output.check_complete(folder, prefix)
        if problem:
            raise CLIError(
                "the calculation %s is not complete: %s" % (path, problem),
                "recompute it (biotsr run or batch with --if-exists skip or "
                "overwrite), or remove its files before analysing",
            )
        record = _output.read_record(folder, prefix)
        for entry in record.get("results") or ():
            keys = os.path.join(folder, entry["files"]["keys"])
            listed.add(os.path.abspath(keys))
            if methods and record.get("method") not in methods:
                continue
            out.append((entry.get("name") or prefix, load_result(keys),
                        record.get("label"), keys))
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if not d.startswith(".") and d != "analysis"]
        for name in filenames:
            if name.endswith(KEYS_SUFFIX) and not name.startswith("."):
                full = os.path.abspath(os.path.join(dirpath, name))
                if full not in listed:
                    raise CLIError(
                        "%s is not part of a completed calculation: no completion "
                        "record lists it" % full,
                        "an interrupted run can leave such a file; recompute it, or "
                        "move it out of %s before analysing" % root,
                    )
    if methods and not out:
        raise CLIError(
            "no completed calculation of %s under %s" % (", ".join(methods), root),
            "`biotsr describe %s` lists the calculations there" % root,
        )
    console.say("loaded     %d saved result(s) from %d calculation(s) in %s%s"
                % (len(out), len(records), root,
                   " (method %s)" % ", ".join(methods) if methods else ""))
    return out


def _samples_from_saved(args, console):
    """`--from-results`: load, then reduce to analysis samples.

    The id is the calculation's prefix, which is what the user sees in the
    output folder; two calculations with one prefix in different folders of
    the tree are refused, because a report with two rows under one name cannot
    be matched back to what produced them.
    """
    from ..analysis import AnalysisSample

    root = args.from_results
    samples = []
    seen = {}
    methods = tuple(getattr(args, "method", None) or ())
    for sample_id, result, label, keys_path in _saved_results(root, console, methods):
        if sample_id in seen:
            raise CLIError(
                "two saved results under %s would both be reported as %r"
                % (root, sample_id),
                "they are %s and %s. Every sample needs a distinguishable id; "
                "analyse those folders separately" % (seen[sample_id], keys_path),
            )
        seen[sample_id] = keys_path
        samples.append(AnalysisSample.from_result(result, sample_id=sample_id,
                                                  label=label))
    return samples


def _samples_from_report(report: BatchReport, output: str):
    """The analysis samples of a batch, read back from the outputs it wrote.

    Each request's result was written and released during the batch, so the
    batch never held the dataset's triplets.  The native key files are read
    back here, in plan order, and each sample's key counts (one entry per
    distinct key) are kept for the analysis, so memory grows with samples x
    distinct keys.  A request skipped because its existing output was verified
    to be the same request contributes exactly as a freshly computed one does.
    """
    from ..analysis import AnalysisSample
    from ..io.native_read import load_result

    samples = []
    for record in report.records:
        if record["status"] not in (STATUS_OK, STATUS_PARTIAL, STATUS_SKIPPED,
                                    STATUS_NO_TRIANGLES):
            continue
        completion = _output.read_record(output, record["sample_id"]) or {}
        for entry in completion.get("results", []):
            result = load_result(os.path.join(output, entry["files"]["keys"]))
            samples.append(AnalysisSample.from_result(
                result, sample_id=entry.get("name") or record["sample_id"],
                label=record.get("label")))
    return samples


#: Flags that change how results are computed, written or downloaded. With
#: ``--from-results`` nothing is computed, written per sample or downloaded, so
#: each of them would be silently ignored; they are refused instead.
_NOT_FOR_SAVED_RESULTS = (
    "chain", "aa_grouping", "mirror", "size_filter", "sse",
    "pdb_dir", "contact_cutoff", "contact_filter", "column", "keep_column",
    "ignore_column", "download_missing", "download_workers", "timeout",
    "retries", "jobs", "if_exists", "on_error", "sheet", "map", "map_file",
    "participant", "inapplicable", "cross_side",
)


def cmd_analyze(args, console) -> int:
    from ..analysis import (
        common_key_intersection,
        common_key_intersection_by_group,
    )
    from ..analysis import PartialSamples
    from ..domain.keyspace import IncompatibleKeySpace, IncompatibleSampling

    from_results = getattr(args, "from_results", None)
    if from_results and args.manifest:
        raise CLIError(
            "--manifest and --from-results both name a dataset; pass one",
            "--manifest computes the results now, --from-results reads results "
            "that were computed earlier",
        )
    if not from_results and not args.manifest:
        raise CLIError(
            "no dataset given",
            "pass --manifest CSV to compute the samples, or --from-results DIR "
            "to analyse results a previous run already saved",
        )
    if getattr(args, "renamed_allow_mixed_sampling", False):
        raise CLIError(
            "--allow-mixed-sampling was renamed --allow-mixed-processing in "
            "1.0.0rc4, and its meaning narrowed: selections (structures, chains, "
            "residue types, molecules, residues, atoms) never need it any more; "
            "only processing settings that differ do",
            "pass --allow-mixed-processing",
        )
    allow_mixed = bool(getattr(args, "allow_mixed_processing", False))
    allow_partial = bool(getattr(args, "allow_partial", False))

    if from_results:
        explicit = getattr(args, "_explicit", set())
        meaningless = [d for d in _NOT_FOR_SAVED_RESULTS if d in explicit]
        if meaningless:
            raise CLIError(
                "--from-results reads results that were already computed, so %s "
                "would have no effect"
                % ", ".join("--" + d.replace("_", "-") for d in meaningless),
                "every saved sample keeps the method, chain, modifiers and "
                "settings it was computed with; drop these flags, or use "
                "--manifest to compute the samples now",
            )
        # Nothing is computed, so there is no batch, no failure report and no
        # partial exit: a saved result either loads or is refused by name.
        samples = _samples_from_saved(args, console)
        report = None
    else:
        # Key frequencies are all a common-key analysis reads, so the samples
        # are always written keys-only.
        args.no_triplets = True
        report = _batch(args, console, command="analyze")
        samples = _samples_from_report(report, args.output)

    if len(samples) < 2:
        raise CLIError(
            "common keys needs at least two samples; %d produced results"
            % len(samples),
            (
                "check --from-results %s" % from_results
                if from_results
                else "check the manifest, or the failure report at %s"
                % os.path.join(args.output, FAILURES_NAME)
            ),
        )

    root = os.path.join(args.output, ANALYSIS_DIR, "common_keys")
    written = []
    try:
        if args.group_by:
            missing = [s.sample_id for s in samples if not s.label]
            if missing:
                raise CLIError(
                    "--group-by %s but %d sample(s) have no %s"
                    % (args.group_by, len(missing), args.group_by),
                    "every row must carry a value in that column for the groups "
                    "to partition the set",
                )
            groups = common_key_intersection_by_group(
                samples, lambda s: s.label or "", allow_mixed_processing=allow_mixed,
                allow_partial=allow_partial,
            )
            os.makedirs(root, exist_ok=True)
            for name, intersection in groups.items():
                written.extend(_write_intersection(root, name, intersection))
                console.say(
                    "%-24s %d samples, %d common keys"
                    % (name, intersection.n_samples, intersection.n_common_keys)
                )
                for warning in intersection.warnings:
                    console.warn("%s: %s" % (name, warning))
        else:
            intersection = common_key_intersection(
                samples, allow_mixed_processing=allow_mixed,
                allow_partial=allow_partial)
            os.makedirs(root, exist_ok=True)
            written.extend(_write_intersection(root, "all_samples", intersection))
            console.say(
                "common keys %d over %d samples"
                % (intersection.n_common_keys, intersection.n_samples)
            )
            for warning in intersection.warnings:
                console.warn(warning)
    except IncompatibleKeySpace as exc:
        raise CLIError(
            "the samples cannot be compared: %s" % exc,
            "their keys are defined differently (label scheme, bins, mirror "
            "signing or another key convention), so the same number means "
            "different triangles; analyse each key space separately. No option "
            "relaxes this",
        )
    except IncompatibleSampling as exc:
        raise CLIError(
            str(exc),
            "pass --allow-mixed-processing to compare them anyway; the report then "
            "records each sample's settings and a warning",
        )
    except PartialSamples as exc:
        raise CLIError(
            str(exc),
            "pass --allow-partial to analyse them as computed; the report then "
            "records what each one lacks. Or fix the requests (or the structures) "
            "so every selection is complete. Partial results stay in %s either "
            "way" % (from_results or args.output),
        )

    console.say("analysis   %s" % root)
    # Partial samples reach this point only when acknowledged; the exit code
    # still says the analysis rests on incomplete selections (ADR-106).
    partial_used = any(s.is_partial for s in samples)
    if from_results:
        # No batch ran, so nothing failed. `--group-by label` reads the label
        # each completion record carries (the manifest's label column when it
        # was computed); a calculation without one is refused above, never
        # given a label derived from its file name.
        console.out(args.output)
        return EXIT_INCOMPLETE if partial_used else EXIT_OK
    code = _report_batch(report, args, console)
    if code == EXIT_OK and partial_used:
        return EXIT_INCOMPLETE
    return code


def _write_intersection(root: str, name: str, intersection) -> List[str]:
    payload = dict(intersection.to_dict())
    payload["group"] = name
    json_path = os.path.join(root, "%s.json" % name)
    with open(json_path, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    csv_path = os.path.join(root, "%s.csv" % name)
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["sample_id", "method", "total_keys", "common_keys",
             "sum_common_keys_freq", "key_space_fingerprint",
             "processing_fingerprint", "sampling_fingerprint"]
        )
        for row, sample in zip(intersection.rows, intersection.samples):
            writer.writerow(
                [
                    row.sample_id,
                    sample.get("method") or "",
                    row.total_keys,
                    row.n_common_keys,
                    row.sum_common_key_frequency,
                    sample["key_space_fingerprint"],
                    sample["processing_fingerprint"],
                    sample["sampling_fingerprint"],
                ]
            )
    return [json_path, csv_path]


# ---------------------------------------------------------------------------
# biotsr describe
# ---------------------------------------------------------------------------


def cmd_describe(args, console) -> int:
    path = args.path
    if os.path.isdir(path):
        return _describe_folder(path, console, args.json)
    if not os.path.isfile(path):
        # A prefix names a calculation: PREFIX -> PREFIX.run.json.
        candidate = path + _output.RECORD_SUFFIX
        if os.path.isfile(candidate):
            path = candidate
        else:
            raise CLIError(
                "no such file or calculation: %s" % path,
                "point describe at an output folder, a calculation's prefix "
                "(e.g. results/1UBQ_A_protein), its .run.json or its .meta.json",
            )
    if path.endswith(".keys") or path.endswith(".triplets"):
        from ..io.native import METADATA_SUFFIX, prefix_of

        path = os.path.join(os.path.dirname(path), prefix_of(path) + METADATA_SUFFIX)
    try:
        with open(path, "r") as fh:
            payload = json.load(fh)
    except (OSError, ValueError) as exc:
        raise CLIError("%s is not readable as JSON: %s" % (path, exc))

    if args.json:
        console.out(json.dumps(payload, indent=2, sort_keys=True))
        return EXIT_OK

    if "results" in payload and "profile" in payload:
        _describe_run(payload, console)
    else:
        _describe_sidecar(payload, console)
    return EXIT_OK


def _describe_folder(path: str, console, as_json: bool) -> int:
    """One line per calculation in an output folder, from its completion record."""
    from ..io.native_read import CANDIDATE_LAYOUT_MESSAGE

    old = os.path.join(path, "run.json")
    if os.path.isfile(old):
        raise CLIError(CANDIDATE_LAYOUT_MESSAGE % old)
    rows = []
    for record_path in sorted(
            os.path.join(path, n) for n in os.listdir(path)
            if n.endswith(_output.RECORD_SUFFIX) and not n.startswith(".")):
        prefix = os.path.basename(record_path)[: -len(_output.RECORD_SUFFIX)]
        record = _output.read_record(path, prefix) or {}
        problem = _output.check_complete(path, prefix)
        entries = record.get("results") or []
        rows.append({
            "prefix": prefix,
            "status": "complete" if problem is None else "incomplete: %s" % problem,
            "method": record.get("method"),
            "rows": (record.get("source") or {}).get("rows"),
            "results": len(entries),
            "n_triangles": sum(e.get("counts", {}).get("n_triangles", 0)
                               for e in entries),
            "triplets": record.get("triplets"),
        })
    if as_json:
        console.out(json.dumps(rows, indent=2, sort_keys=True))
        return EXIT_OK
    if not rows:
        raise CLIError("%s holds no calculation (no PREFIX%s completion record)"
                       % (path, _output.RECORD_SUFFIX))
    for row in rows:
        console.out("%-56s %-28s %8d triangles  %s%s" % (
            row["prefix"], row["method"] or "?", row["n_triangles"], row["status"],
            "" if row["triplets"] else "  (keys only)"))
    return EXIT_OK


def _describe_run(payload, console) -> None:
    profile = payload.get("profile", {})
    console.out("command          biotsr %s" % payload.get("command", "?"))
    console.out("biotsr version   %s" % payload.get("biotsr_version", "?"))
    console.out("method           %s" % payload.get("method", "?"))
    console.out("prefix           %s%s" % (payload.get("prefix", "?"),
                                          " (one per residue)" if payload.get(
                                              "per_residue") else ""))
    console.out("complete         %s" % payload.get("complete"))
    source = payload.get("source") or {}
    if source.get("rows"):
        console.out("source rows      %s of %s%s" % (
            ", ".join(str(r) for r in source["rows"]), source.get("manifest"),
            " (participant %s)" % source["participant"]
            if source.get("participant") else ""))
    console.out("profile          %s" % profile.get("resolved"))
    if profile.get("active") and profile.get("active") != profile.get("resolved"):
        console.out("active profile   %s" % profile["active"])
    if profile.get("sse_strand_identity"):
        console.out("sse strand id    %s" % profile["sse_strand_identity"])
    console.out("parameters")
    for key in sorted(payload.get("parameters", {})):
        console.out("  %-16s %s" % (key, payload["parameters"][key]))
    console.out("results          %d" % len(payload.get("results", [])))
    for entry in payload.get("results", []):
        counts = entry.get("counts", {})
        console.out(
            "  %-24s %s  keyspace %s  sampling %s"
            % (
                entry.get("name") or entry.get("target") or "-",
                "%d triangles, %d keys"
                % (counts.get("n_triangles", 0), counts.get("n_unique_keys", 0)),
                entry.get("key_space_fingerprint", "")[:12],
                entry.get("sampling_fingerprint", "")[:12],
            )
        )
        for warning in entry.get("warnings") or ():
            console.out("    warning: %s" % warning)


def _describe_sidecar(payload, console) -> None:
    console.out("method           %s" % payload.get("method"))
    console.out("structure        %s" % payload.get("structure_id"))
    console.out("chain            %s" % payload.get("chain"))
    console.out("key space        %s" % payload.get("key_space_fingerprint"))
    for key in sorted(payload.get("key_space", {})):
        console.out("  %-16s %s" % (key, payload["key_space"][key]))
    console.out("sampling         %s" % payload.get("sampling_fingerprint"))
    for key in sorted(payload.get("sampling", {})):
        console.out("  %-16s %s" % (key, payload["sampling"][key]))
    counts = payload.get("counts", {})
    console.out("counts")
    for key in sorted(counts):
        console.out("  %-16s %s" % (key, counts[key]))
    if payload.get("sse_coverage"):
        console.out("sse coverage     %s" % payload["sse_coverage"])
    warnings = payload.get("warnings") or []
    if warnings:
        console.out("warnings         %d" % len(warnings))
        for warning in warnings:
            console.out("  %s" % warning)
