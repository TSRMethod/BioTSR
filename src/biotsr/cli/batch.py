"""Planning and running a dataset: the engine behind ``biotsr batch``.

A batch is run in two passes.

**Planning** reads the whole dataset and decides everything that can be
decided without computing: which column is which field, which method each
request runs (from the manifest's ``method`` column, or from ``--method`` on
the command line, one request per source row and method), each request's
chains and where they came from, every run specification, every output
prefix, which structures are present, which would have to be downloaded, and
what happens to outputs that already exist.  Problems that make the dataset as
a whole unusable -- an ambiguous method choice, two different requests that
would write the same output prefix -- are refused here, before a single
structure is read, so a thousand-request job cannot die at request 900 over
something visible at request 1.  Two rows asking for exactly the same
calculation (a wide sheet repeats a protein chain on every ligand row) become
one request, and the summary lists every source row it serves.

**Execution** computes the requests that are left, one after another in this
process (``jobs=1``) or in worker processes (``jobs>1``, joblib).  Either way a
request is computed by the same function with the same arguments, into a
staging directory private to this invocation, and installed into the shared
output folder with a completion record written last (``biotsr.cli.output``).
What comes back to the coordinating process is a one-line record, never a
result: memory does not grow with the number of requests.  The summary lists
every request in plan order with one of these statuses:

``ok``            computed and written
``no_triangles``  computed and written, but the selection has fewer than three
                  points, so the result has no keys; ``message`` says why
``skipped``       an existing output was verified to be this exact request
``inapplicable``  the source row does not describe this method's calculation,
                  and ``--inapplicable skip`` was given; ``message`` says why
``failed``        the request could not be computed; ``message`` says why
``interrupted``   it was being computed when the run was stopped
``not_started``   the run stopped before it began
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import signal
import tempfile
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from .. import TSRError
from ..domain.selection_record import (
    INCOMPLETE_WARNINGS,
    STRUCTURE_NOTES,
    selection_omissions,
    selection_sides,
)
from ..io import native
from ..methods.registry import get_method
from ..retrieval import (
    DEFAULT_RETRIES,
    DEFAULT_TIMEOUT,
    DEFAULT_WORKERS,
    cached_structure_path,
    check_cached_file,
    fetch_structures,
    file_sha256,
    normalise_pdb_id,
    publish_file_mode,
    FetchInterrupted,
    MAX_WORKERS,
    validate_download_settings,
)
from . import output as _output
from .errors import EXIT_PARTIAL, CLIError
from .manifest import CHAIN_FROM_DEFAULT, CHAIN_FROM_ROW, read_manifest, resolve_structure_path
from .spec import (
    build_spec,
    can_stream_triplets,
    execute,
    execute_streaming,
    output_prefix,
    sweeps_residues,
)

__all__ = [
    "BatchOptions",
    "BatchPlan",
    "BatchReport",
    "INAPPLICABLE_FAIL",
    "INAPPLICABLE_SKIP",
    "ON_ERROR_CONTINUE",
    "ON_ERROR_STOP",
    "STATUS_OK",
    "STATUS_NO_TRIANGLES",
    "STATUS_PARTIAL",
    "STATUS_SKIPPED",
    "STATUS_INAPPLICABLE",
    "STATUS_CHAIN_NOT_SPECIFIED",
    "STATUS_FAILED",
    "STATUS_INTERRUPTED",
    "STATUS_NOT_STARTED",
    "SUMMARY_NAME",
    "FAILURES_NAME",
    "FETCH_REPORT_NAME",
    "SUMMARY_COLUMNS",
    "compute_into",
    "plan_batch",
    "run_batch",
    "structure_id_for",
]

ON_ERROR_STOP = "stop"
ON_ERROR_CONTINUE = "continue"

#: What a source row that does not describe a method's calculation becomes.
INAPPLICABLE_FAIL = "fail"
INAPPLICABLE_SKIP = "skip"
INAPPLICABLE_CHOICES = (INAPPLICABLE_FAIL, INAPPLICABLE_SKIP)

STATUS_OK = "ok"
#: Written like ``ok``, but the calculation produced no triangle: not an ordinary
#: result, and labelled so a summary is never read as if it were one.
STATUS_NO_TRIANGLES = "no_triangles"
#: Written like ``ok``, but something the request named explicitly -- a
#: residue, residue type, atom or molecule -- was not used: not found, not
#: supported (a modified residue BioTSR does not read), or excluded by the
#: method's completeness rule.  The message names each one, so a partial
#: selection is never read as the complete one that was asked for.
STATUS_PARTIAL = "partial"
#: The result warnings that say a calculation is partial (``api`` writes
#: them).  The status itself is decided from the structured selection record
#: (``selection_omissions``); these words also mark a result written before
#: that record existed, and a per-residue sweep's own omissions.
PARTIAL_WARNINGS = INCOMPLETE_WARNINGS
STATUS_SKIPPED = "skipped"
STATUS_INAPPLICABLE = "inapplicable"
#: The request needs a chain (or, for a cross method, a second chain) that
#: neither its row nor an explicit ``--chain`` default supplies.  Nothing is
#: computed and no chain is assumed; the batch exits 4, as for a partial result
#: (ADR-112).
STATUS_CHAIN_NOT_SPECIFIED = "chain_not_specified"
STATUS_FAILED = "failed"
STATUS_INTERRUPTED = "interrupted"
STATUS_NOT_STARTED = "not_started"

SUMMARY_NAME = "batch_summary.csv"
FAILURES_NAME = "batch_failures.csv"
FETCH_REPORT_NAME = "fetch_report.csv"

SUMMARY_COLUMNS = (
    "request",
    "rows",
    "participant",
    "sample_id",
    "status",
    "method",
    "profile",
    "structure_id",
    "chain",
    "chain_2",
    "selection",
    "n_results",
    "n_triangles",
    "n_unique_keys",
    "output",
    "message",
    "chain_source",
)

#: How a planned request will be handled.
ACTION_RUN = "run"
ACTION_SKIP = "skip"
ACTION_FAIL = "fail"
ACTION_INAPPLICABLE = "inapplicable"
ACTION_NO_CHAIN = "chain_not_specified"
#: Dry run only: the structure is not present yet, so it would be downloaded.
ACTION_DOWNLOAD = "download"


@dataclass
class BatchOptions:
    """Everything one batch needs besides the dataset rows themselves."""

    manifest: str
    output: str
    pdb_dir: str = ""
    #: The methods chosen on the command line (``--method``, repeatable).  Each
    #: source row becomes one request per method.  Empty: each row names its
    #: own (a ``method`` column).
    methods: Sequence[str] = ()
    #: Row defaults: a row's own value always wins (``{"chain": "A"}``).
    defaults: Mapping[str, Any] = field(default_factory=dict)
    #: ``{"label": "group"}`` -- read field ``label`` from column ``group``.
    columns: Mapping[str, str] = field(default_factory=dict)
    keep_columns: Sequence[str] = ()
    ignore_columns: Sequence[str] = ()
    #: A wide-sheet layout (``"tsr"`` or ``"custom"``), or ``None`` for a
    #: one-request-per-row manifest.
    sheet: Optional[str] = None
    #: ``--map`` / ``--map-file`` entries (``biotsr.cli.sheet``).
    maps: Sequence[Tuple[str, str, str, str, str]] = ()
    #: ``--cross-side``: ``{(method, side): column or "whole"}``.
    cross_sides: Mapping[Tuple[str, str], str] = field(default_factory=dict)
    participant: str = "1"
    inapplicable: str = INAPPLICABLE_FAIL
    triplets: bool = True
    if_exists: str = _output.IF_EXISTS_FAIL
    on_error: str = ON_ERROR_CONTINUE
    jobs: int = 1
    download_missing: bool = False
    download_workers: int = DEFAULT_WORKERS
    timeout: float = DEFAULT_TIMEOUT
    retries: int = DEFAULT_RETRIES
    dry_run: bool = False
    command: str = "batch"
    argv: Optional[List[str]] = None


@dataclass
class _Request:
    index: int
    lines: List[int]
    where: str
    method: str
    participant: str
    values: Dict[str, Any]
    metadata: Dict[str, str]
    chain_source: str
    record: Dict[str, Any]
    columns: Dict[str, str] = field(default_factory=dict)
    spec: Any = None
    prefix: str = ""
    fetch_id: Optional[str] = None
    request: Optional[Dict[str, Any]] = None
    action: str = ACTION_RUN
    #: Said in the summary whatever the outcome: that this request is one
    #: sample serving several source requests (``_merge_identical``).
    note: str = ""


@dataclass
class BatchPlan:
    options: BatchOptions
    rows: List[_Request]
    metadata_columns: Tuple[str, ...]
    #: Columns of a wide sheet that no planned request reads.
    unused_columns: Tuple[str, ...] = ()
    #: The ``--output`` lock while the batch runs (never taken by a dry run).
    lock: Any = None

    def count(self, action: str) -> int:
        return sum(1 for row in self.rows if row.action == action)


@dataclass
class BatchReport:
    """What a batch did.  ``records`` are the summary rows, in plan order."""

    records: List[Dict[str, Any]]
    summary_path: str = ""
    failures_path: str = ""
    fetch_report_path: str = ""
    interrupted: bool = False
    #: Set when the worker pool itself failed (a worker process died, or the
    #: executor broke), describing what happened.  The requests it affected say
    #: so individually; this says why.
    pool_error: str = ""
    #: Staging leftovers of earlier, abruptly terminated invocations on this
    #: machine that were removed before this batch ran, and ones belonging to
    #: another machine that were left alone because they cannot be checked.
    removed_leftovers: Tuple[str, ...] = ()
    foreign_leftovers: Tuple[str, ...] = ()
    unused_columns: Tuple[str, ...] = ()

    def _count(self, status: str) -> int:
        return sum(1 for r in self.records if r["status"] == status)

    @property
    def n_ok(self) -> int:
        return self._count(STATUS_OK)

    @property
    def n_no_triangles(self) -> int:
        return self._count(STATUS_NO_TRIANGLES)

    @property
    def n_partial(self) -> int:
        return self._count(STATUS_PARTIAL)

    @property
    def n_skipped(self) -> int:
        return self._count(STATUS_SKIPPED)

    @property
    def n_skipped_partial(self) -> int:
        """Skipped requests whose kept output is itself partial."""
        return sum(1 for r in self.records
                   if r["status"] == STATUS_SKIPPED and r.get("kept_partial"))

    @property
    def n_no_triangles_partial(self) -> int:
        """Requests with no triangles whose selection also fell short."""
        return sum(1 for r in self.records
                   if r["status"] == STATUS_NO_TRIANGLES and r.get("partial_selection"))

    @property
    def n_incomplete(self) -> int:
        """Every request whose result does not hold the whole requested selection.

        ``partial`` rows, kept outputs that are partial, and rows with no
        triangles whose selection also lacked something it named.  Any of them
        makes the batch exit 4 (ADR-106), and so does a request whose chain
        was not specified (ADR-112).
        """
        return (self.n_partial + self.n_skipped_partial + self.n_no_triangles_partial
                + self.n_chain_not_specified)

    @property
    def n_chain_not_specified(self) -> int:
        return self._count(STATUS_CHAIN_NOT_SPECIFIED)

    @property
    def n_inapplicable(self) -> int:
        return self._count(STATUS_INAPPLICABLE)

    @property
    def n_failed(self) -> int:
        return self._count(STATUS_FAILED)

    @property
    def n_interrupted(self) -> int:
        return self._count(STATUS_INTERRUPTED)

    @property
    def n_not_started(self) -> int:
        return self._count(STATUS_NOT_STARTED)

    @property
    def failures(self) -> List[Dict[str, Any]]:
        return [r for r in self.records if r["status"] == STATUS_FAILED]


class BatchInterrupted(KeyboardInterrupt):
    """Ctrl-C or SIGTERM during a batch, after its reports were written."""

    def __init__(self, report: BatchReport):
        super().__init__("interrupted")
        self.report = report


def structure_id_for(spec) -> str:
    if spec.structure_id:
        return spec.structure_id
    base = os.path.basename(spec.pdb_path)
    return base[:-4] if base.lower().endswith(".pdb") else base


# ---------------------------------------------------------------------------
# Planning
# ---------------------------------------------------------------------------


def _lines(rows, limit: int = 10) -> str:
    numbers = [str(r) for r in rows[:limit]]
    more = len(rows) - len(numbers)
    return ", ".join(numbers) + (" and %d more" % more if more > 0 else "")


def _locate(row: Mapping[str, Any], options: BatchOptions, where: str):
    """``(path, structure_id, id_to_download)`` for one request.

    With ``download_missing`` a request naming a ``pdb_id`` that is not in the
    structure directory is planned against the path its download will have.
    """
    if row.get("pdb_path") or not options.download_missing:
        path, structure_id = resolve_structure_path(row, options.pdb_dir or "", where)
        return path, structure_id, None
    pdb_id = row.get("pdb_id")
    if not pdb_id:
        path, structure_id = resolve_structure_path(row, options.pdb_dir or "", where)
        return path, structure_id, None
    if not options.pdb_dir:
        raise CLIError(
            "%s: --download-missing needs --pdb-dir, the directory downloads "
            "are kept in" % where
        )
    structure_id = row.get("structure_id") or pdb_id
    existing = cached_structure_path(options.pdb_dir, str(pdb_id))
    normalised = normalise_pdb_id(str(pdb_id))
    if existing:
        if normalised:
            problem = check_cached_file(existing, normalised, options.pdb_dir)
            if problem:
                raise CLIError(
                    "%s: the cached structure %s is not used: %s; `biotsr fetch %s "
                    "--pdb-dir %s --refresh` downloads it again"
                    % (where, existing, problem, normalised, options.pdb_dir)
                )
        return existing, structure_id, None
    if normalised is None:
        raise CLIError(
            "%s: %r is not in %s and is not a PDB id that can be downloaded"
            % (where, pdb_id, options.pdb_dir)
        )
    return os.path.join(options.pdb_dir, normalised + ".pdb"), structure_id, normalised


def _selection_text(values: Mapping[str, Any]) -> str:
    keys = ("residue", "residue_name", "residues", "residue_types", "atoms",
            "entity_name", "entity_id", "entity_names", "atom_selection",
            "calpha_residues", "residue_types_2", "residues_2", "atom_selection_2",
            "calpha_residues_2", "contact_cutoff", "contact_filter",
            "no_contact_filter", "aa_grouping", "mirror", "size_filter", "sse")
    return "; ".join("%s=%s" % (k, values[k]) for k in keys
                     if values.get(k) not in (None, "", False))


def _side_texts(spec) -> List[str]:
    """Each side's resolved selection, for the preview and the summary.

    One line per participant of the calculation: a one-sided method has one;
    a cross method two, each naming its role, chain and narrowing ("whole
    chain" when nothing narrows it), so a side computed whole is never taken
    for a narrowed one.
    """
    selection, cross = spec.selection, spec.cross

    def narrowing(values, suffix=""):
        parts = ["%s=%s" % (k, values[k + suffix])
                 for k in ("residue_types", "residues", "atom_selection",
                           "calpha_residues") if values.get(k + suffix)]
        return "; ".join(parts) or "whole chain"

    if not spec.is_cross:
        return ["chain %s: %s" % (spec.chain, _selection_text(selection)
                                  or "whole chain")]
    from ..api import cross_entity_role

    role = cross_entity_role(spec.method)
    if role is None:
        return ["side 1, chain %s: %s" % (spec.chain, narrowing(cross)),
                "side 2, chain %s: %s" % (cross.get("chain_2"), narrowing(cross, "_2"))]
    if cross.get("entity_names"):
        molecule = "every %s of the chain" % cross["entity_names"]
    else:
        molecule = "%s %s" % (cross.get("entity_name"),
                              "copy %s" % cross["entity_id"] if cross.get("entity_id")
                              else "(all copies)")
    atoms = "atoms=%s" % cross["atoms"] if cross.get("atoms") else "all atoms"
    return ["%s, chain %s: %s, %s" % (role, spec.chain, molecule, atoms),
            "protein, chain %s: %s" % (cross.get("chain_2"), narrowing(cross))]


def _base_record(values: Mapping[str, Any], lines: Sequence[int], method: str,
                 participant: str, chain_source: str) -> Dict[str, Any]:
    return {
        "request": 0,
        "rows": ";".join(str(n) for n in lines),
        "participant": participant,
        "sample_id": values.get("sample_id", ""),
        "status": STATUS_FAILED,
        "method": method,
        "profile": "",
        "structure_id": "",
        "chain": values.get("chain", ""),
        "chain_2": values.get("chain_2", ""),
        "selection": _selection_text(values),
        "n_results": 0,
        "n_triangles": 0,
        "n_unique_keys": 0,
        "output": "",
        "message": "",
        "chain_source": chain_source,
    }


def _raw_requests(options: BatchOptions):
    """``(requests, metadata_columns, unused_columns)`` before any spec is built."""
    from .sheet import read_sheet

    methods = list(dict.fromkeys(options.methods or ()))
    if options.sheet:
        sheet_requests, _map, unused = read_sheet(
            options.manifest, layout=options.sheet, methods=methods,
            participant=options.participant, columns=options.columns,
            maps=options.maps, keep_columns=options.keep_columns,
            ignore_columns=options.ignore_columns, cross_sides=options.cross_sides,
        )
        out = []
        for item in sheet_requests:
            values = dict(options.defaults)
            values.update(item.values)
            values["method"] = item.method
            request = _Request(
                index=0, lines=[item.line], where=item.where, method=item.method,
                participant=item.participant, values=values,
                metadata=dict(item.metadata),
                chain_source=CHAIN_FROM_ROW if item.values.get("chain") else (
                    CHAIN_FROM_DEFAULT if values.get("chain") else ""),
                record={}, columns=dict(item.columns),
            )
            if item.not_applied:
                request.note = "not applied: %s" % "; ".join(
                    "%s (%s)" % pair for pair in item.not_applied)
            missing = {f: c for f, c in item.chain_missing.items()
                       if not values.get(f)}
            if item.inapplicable:
                request.action = ACTION_INAPPLICABLE
                request.values["_inapplicable"] = item.inapplicable
            elif item.problem:
                request.values["_problem"] = item.problem
            elif missing:
                request.values["_no_chain"] = _no_chain_message(
                    item.method, list(missing), "the row leaves %s blank" % " and ".join(
                        missing.values()), item.participant)
            out.append(request)
        return out, tuple(options.keep_columns), tuple(unused)

    manifest_rows = read_manifest(
        options.manifest,
        defaults=options.defaults,
        columns=options.columns,
        keep_columns=options.keep_columns,
        ignore_columns=options.ignore_columns,
    )
    named = [r["_row"] for r in manifest_rows if r.get("method")]
    if methods and named:
        raise CLIError(
            "%s names a method in %d row(s) (line %s) and --method also chooses "
            "methods" % (options.manifest, len(named), _lines(named)),
            "use one: drop --method to run each row's own method, or pass "
            "--ignore-column method to run the --method list on every row",
        )
    if not methods:
        missing = [r["_row"] for r in manifest_rows if not r.get("method")]
        if missing:
            raise CLIError(
                "%d row(s) of %s name no method: line %s"
                % (len(missing), options.manifest, _lines(missing)),
                "choose the methods with --method NAME (repeatable), which runs "
                "each of them on every row, or give the manifest a method column",
            )
    out = []
    for raw in manifest_rows:
        values = {k: v for k, v in raw.items() if not k.startswith("_")}
        for method in (methods or [values["method"]]):
            row_values = dict(values, method=method)
            request = _Request(
                index=0, lines=[raw["_row"]], where=raw["_where"], method=method,
                participant="", values=row_values, metadata=dict(raw["_metadata"]),
                chain_source=raw["_chain_source"], record={},
            )
            if methods:
                reason = _inapplicable_reason(method, row_values)
                if reason:
                    request.action = ACTION_INAPPLICABLE
                    request.values["_inapplicable"] = reason
            if request.action != ACTION_INAPPLICABLE:
                missing = _missing_chains(method, row_values)
                if missing:
                    request.values["_no_chain"] = _no_chain_message(
                        method, missing, "the row gives no %s and no --chain "
                        "default was given" % " or ".join(missing), "")
            out.append(request)
    metadata_columns = tuple(manifest_rows[0]["_metadata"]) if manifest_rows else ()
    return out, metadata_columns, ()


def _missing_chains(method: str, values: Mapping[str, Any]) -> List[str]:
    """The chain fields ``method`` needs that ``values`` leaves empty."""
    try:
        recipe = get_method(method)
    except KeyError:
        return []  # an unknown method is reported, by name, when the row is planned
    needed = ["chain", "chain_2"] if recipe.family == "cross" else ["chain"]
    return [f for f in needed if not values.get(f)]


def _no_chain_message(method: str, missing, where: str, participant: str) -> str:
    return ("%s%s needs %s, and %s; BioTSR never assumes a chain: give it in the "
            "row, map a column to it (--column chain=HEADER or --map), or pass "
            "--chain for rows that leave it blank"
            % (method, " (participant 2)" if participant == "2" else "",
               " and ".join(missing), where))


def _inapplicable_reason(method: str, values: Mapping[str, Any]) -> str:
    """Why a long-manifest row does not describe ``method``'s calculation, or ``""``.

    A row that sets none of the arguments the method needs beyond the chain
    does not describe it (a ``drug`` row for ``carbohydrate``).
    """
    from ..methods.arguments import required_arguments

    recipe = get_method(method)
    if recipe.family == "cross":
        needed = ("entity_name", "entity_names") if recipe.name in (
            "drug_protein_cross", "carbohydrate_protein_cross",
            "nucleotide_protein_cross") else ("chain_2",)
    elif recipe.scope == "per_entity":
        needed = ()        # no residue: one calculation per residue of the chain
    else:
        needed = tuple(required_arguments(recipe))
    if not needed or any(values.get(n) for n in needed):
        return ""
    return "the row sets no %s for %s" % (" / ".join(needed), method)


def plan_batch(options: BatchOptions) -> BatchPlan:
    """Read the dataset and decide every request, without computing or writing.

    Raises :class:`CLIError` for a dataset that cannot be run as a whole:
    unreadable or ambiguous columns, a method choice made twice, rows without a
    chain, or two different requests that would write the same output prefix.
    A request that is wrong on its own is planned as a failure and reported
    with its source row(s).
    """
    if options.inapplicable not in INAPPLICABLE_CHOICES:
        raise CLIError("inapplicable=%r is not one of %s"
                       % (options.inapplicable, ", ".join(INAPPLICABLE_CHOICES)))
    raw, metadata_columns, unused = _raw_requests(options)
    clash = sorted(set(metadata_columns) & set(SUMMARY_COLUMNS))
    if clash:
        raise CLIError(
            "--keep-column %s would collide with the batch summary's own column"
            % ", ".join(clash),
            "rename that column in the manifest, or drop it with --ignore-column",
        )
    expanded = len(raw) > len({tuple(r.lines) for r in raw})

    for request in raw:
        values = request.values
        request.record = _base_record(values, request.lines, request.method,
                                      request.participant, request.chain_source)
        if request.action == ACTION_INAPPLICABLE:
            reason = values.pop("_inapplicable")
            if options.inapplicable == INAPPLICABLE_SKIP:
                request.record["status"] = STATUS_INAPPLICABLE
                request.record["message"] = reason
            else:
                request.action = ACTION_FAIL
                request.record["message"] = (
                    "%s (pass --inapplicable skip to record such rows as "
                    "inapplicable instead of failed)" % reason)
            continue
        problem = values.pop("_problem", "")
        if problem:
            _fail(request, problem)
            continue
        no_chain = values.pop("_no_chain", "")
        if no_chain:
            request.action = ACTION_NO_CHAIN
            request.record["status"] = STATUS_CHAIN_NOT_SPECIFIED
            request.record["message"] = no_chain
            continue
        try:
            path, structure_id, fetch_id = _locate(values, options, request.where)
            spec_values = {k: v for k, v in values.items()}
            spec_values["pdb"] = path
            if structure_id:
                spec_values["structure_id"] = structure_id
            spec_values["keep_triplets"] = options.triplets
            if expanded and spec_values.get("sample_id"):
                # One source row, several requests: the row's sample_id names
                # them all, each followed by its method (and participant).
                spec_values["sample_id"] = "%s_%s%s" % (
                    spec_values["sample_id"], request.method,
                    "_p2" if request.participant == "2" else "")
            request.spec = build_spec(spec_values, source_row=request.lines[0],
                                      where=request.where)
            request.fetch_id = fetch_id
            request.prefix = output_prefix(request.spec)
            request.record.update(
                label=request.spec.label,
                method=request.spec.method,
                profile=get_method(request.spec.method).default_profile.name,
                chain=request.spec.chain,
                chain_2=request.spec.cross.get("chain_2", ""),
                sample_id=request.prefix,
                output=os.path.join(options.output, request.prefix),
            )
            if request.spec.is_cross:
                request.record["selection"] = " | ".join(_side_texts(request.spec))
            if fetch_id and options.dry_run:
                request.action = ACTION_DOWNLOAD
        except CLIError as exc:
            _fail(request, exc.message)
        except TSRError as exc:
            _fail(request, str(exc))

    rows = _merge_identical(raw)
    _refuse_colliding_outputs(rows, options)
    for number, request in enumerate(rows, 1):
        request.index = number
        request.record["request"] = number
    return BatchPlan(options=options, rows=rows, metadata_columns=metadata_columns,
                     unused_columns=tuple(unused))


def _spec_key(spec) -> str:
    return json.dumps({
        "method": spec.method, "pdb": spec.pdb_path, "structure_id": spec.structure_id,
        "chain": spec.chain, "selection": dict(spec.selection),
        "cross": dict(spec.cross), "modifiers": dict(spec.modifiers),
        "sample_id": spec.sample_id,
    }, sort_keys=True, default=str)


def _merge_identical(requests: List[_Request]) -> List[_Request]:
    """One request per distinct calculation *and* sample; every source row kept.

    Two requests merge only when they are the same calculation of the same
    sample: the same method, structure, chains, selection, modifiers and
    sample id (an explicit ``sample_id``, or the automatic one, which names
    exactly those).  Rows with distinct sample ids are never merged -- each is
    its own sample, computed and written under its own prefix.  The merged
    request lists every source row and participant.

    A merge never chooses between rows: requests of one calculation and sample
    whose labels differ are refused (one output cannot carry two labels; give
    the rows distinct sample ids), and a kept metadata column whose values
    differ records each row's value (``line 11=x; line 12=y``).
    """
    out: List[_Request] = []
    by_key: Dict[str, _Request] = {}
    members: Dict[int, List[_Request]] = {}
    for request in requests:
        if request.spec is None or request.action != ACTION_RUN and \
                request.action != ACTION_DOWNLOAD:
            out.append(request)
            continue
        key = _spec_key(request.spec)
        first = by_key.get(key)
        if first is None:
            by_key[key] = request
            members[id(request)] = [request]
            out.append(request)
            continue
        members[id(first)].append(request)

    conflicts = []
    for first in out:
        group = members.get(id(first))
        if not group or len(group) == 1:
            continue
        labels = {}
        for member in group:
            labels.setdefault(member.spec.label or "", []).extend(member.lines)
        if len(labels) > 1:
            conflicts.append("%s: %s" % (first.prefix, "; ".join(
                "label %r on line %s" % (label, _lines(lines))
                for label, lines in labels.items())))
            continue
        for name in first.metadata:
            values = {}
            for member in group:
                for line in member.lines:
                    values.setdefault(line, member.metadata.get(name, ""))
            if len(set(values.values())) > 1:
                first.metadata[name] = "; ".join(
                    "line %d=%s" % (line, value) for line, value in values.items())
        for member in group[1:]:
            first.lines.extend(n for n in member.lines if n not in first.lines)
        first.record["rows"] = ";".join(str(n) for n in first.lines)
        parts = list(dict.fromkeys(m.participant for m in group if m.participant))
        first.participant = ";".join(parts)
        first.record["participant"] = first.participant
        first.note = "; ".join(filter(None, (first.note, (
            "one sample from %d source requests (line %s): the same calculation "
            "and sample id, computed once" % (len(group), _lines(first.lines))))))
    if conflicts:
        raise CLIError(
            "%d calculation(s) are requested by rows with different labels: %s%s"
            % (len(conflicts), "; ".join(conflicts[:5]),
               " and %d more" % (len(conflicts) - 5) if len(conflicts) > 5 else ""),
            "one output is one sample and carries one label; give those rows "
            "distinct sample_id values to keep them as separate samples, or "
            "correct the label",
        )
    return out


def _fail(row: _Request, message: str) -> None:
    row.action = ACTION_FAIL
    row.record["status"] = STATUS_FAILED
    row.record["message"] = _unprefixed(message, row.where)


def _unprefixed(message: str, where: str) -> str:
    """The message without a leading ``<where>: ``; the report adds the row itself."""
    prefix = where + ": "
    return message[len(prefix):] if message.startswith(prefix) else message


def _refuse_colliding_outputs(rows: List[_Request], options: BatchOptions) -> None:
    """Two different requests, one output prefix: refuse before anything runs.

    Compared case-insensitively, because on the default macOS and Windows file
    systems ``1UBQ_A_protein.keys`` and ``1ubq_A_protein.keys`` are one file.
    A per-residue sweep writes ``<prefix>_<NAME>_<id>`` for each residue, so a
    prefix extending a sweep's prefix collides with it too.
    """
    by_name: Dict[str, List[_Request]] = {}
    for row in rows:
        if row.prefix and row.action in (ACTION_RUN, ACTION_DOWNLOAD):
            by_name.setdefault(row.prefix.casefold(), []).append(row)
    groups = [group for group in by_name.values() if len(group) > 1]
    sweeps = [row for row in rows if row.prefix and row.spec is not None
              and sweeps_residues(row.spec)
              and row.action in (ACTION_RUN, ACTION_DOWNLOAD)]
    for sweep in sweeps:
        stem = (sweep.prefix + "_").casefold()
        clashing = [r for key, group in by_name.items() if key.startswith(stem)
                    for r in group]
        if clashing:
            groups.append([sweep] + clashing)
    if not groups:
        return
    described = "; ".join(
        "%s (rows %s)" % (group[0].prefix,
                          _lines([n for r in group for n in r.lines]))
        for group in groups[:5]
    )
    raise CLIError(
        "%d output prefix(es) in %s would each be written by more than one "
        "different request: %s%s"
        % (len(groups), options.manifest, described,
           " and %d more" % (len(groups) - 5) if len(groups) > 5 else ""),
        "give those rows distinct sample_id values. Automatic prefixes name the "
        "structure, chains, method, target and modifiers, so this happens when "
        "two rows give one structure id to different files, or when sample ids "
        "repeat; names differing only in letter case also collide on "
        "case-insensitive file systems",
    )


def _resolve_existing(row: _Request, options: BatchOptions) -> None:
    """Hash the structure, build the request, and apply ``--if-exists``."""
    try:
        sha = file_sha256(row.spec.pdb_path)
    except OSError as exc:
        _fail(row, "%s: cannot read structure %s: %s"
              % (row.where, row.spec.pdb_path, exc.strerror or exc))
        return
    row.request = _output.request_record(row.spec, sha, structure_id_for(row.spec),
                                         per_residue=sweeps_residues(row.spec))
    state, _record = _output.existing_state(options.output, row.prefix)
    if state == "absent":
        return
    if state == "foreign":
        _fail(row, "output %s cannot be written: %s" % (
            row.prefix, _output.describe_foreign(options.output, row.prefix, _record)))
        return
    if state == "incomplete":
        if options.if_exists == _output.IF_EXISTS_FAIL:
            _fail(row, "output %s was left incomplete by an interrupted run (pass "
                       "--if-exists skip to resume, or --if-exists overwrite)"
                  % row.prefix)
        return
    if options.if_exists == _output.IF_EXISTS_FAIL:
        _fail(row, "output %s already exists (pass --if-exists overwrite to "
                   "replace it, or --if-exists skip to keep it)" % row.prefix)
    elif options.if_exists == _output.IF_EXISTS_SKIP:
        problem = _output.verify_existing(options.output, row.prefix, row.request,
                                          options.triplets)
        if problem:
            _fail(row, "output %s exists but cannot be reused: %s (it is left "
                       "untouched; pass --if-exists overwrite to recompute it, or "
                       "choose a new --output)" % (row.prefix, problem))
        else:
            row.action = ACTION_SKIP
            row.record.update(_skipped_record(_record))
            row.record.update(_counts_from_record(_record))


def _kept_omissions(record: Optional[Mapping[str, Any]]) -> List[str]:
    """What a kept calculation's record says it did not get, in words.

    Read from the completion record: each result's structured ``omissions``
    and its warnings, and the calculation's own warnings (a per-residue
    sweep's).  A record without ``omissions`` predates them and is read from
    its warnings alone.
    """
    record = record or {}
    notes: List[str] = []
    entries = record.get("results") or []
    for warning in list(record.get("warnings") or ()) + [
            w for e in entries for w in (e.get("warnings") or ())]:
        if str(warning).startswith(INCOMPLETE_WARNINGS) and warning not in notes:
            notes.append(warning)
    if not notes and any(e.get("omissions") for e in entries):
        notes.append("the kept result's selection record lists omissions")
    return notes


def _skipped_record(record: Optional[Mapping[str, Any]]) -> Dict[str, Any]:
    """Status and message of a request whose verified existing output is kept.

    A kept output that is itself partial says so: resuming never turns an
    incomplete selection into one that reads as complete.
    """
    notes = _kept_omissions(record)
    message = "verified existing output of the same request"
    if notes:
        message += "; kept partial result: " + "; ".join(notes)
    return {"status": STATUS_SKIPPED, "message": message,
            "kept_partial": bool(notes)}


def _counts_from_record(record: Optional[Mapping[str, Any]]) -> Dict[str, Any]:
    """The summary numbers of a verified existing calculation, from its record."""
    entries = (record or {}).get("results") or []
    if not entries:
        return {}
    return {
        "structure_id": entries[0].get("structure_id") or "",
        "n_results": len(entries),
        "n_triangles": sum(e.get("counts", {}).get("n_triangles", 0) for e in entries),
        "n_unique_keys": sum(e.get("counts", {}).get("n_unique_keys", 0)
                             for e in entries),
    }


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------

#: The name a streamed triplet file has in the staging directory until it is
#: known to be complete.
_STREAMING_NAME = ".streaming-triplets.part"


def compute_into(spec, writer, *, prefix: str, triplets: bool, command: str,
                 argv=None, request=None, metadata=None, source=None):
    """Compute ``spec`` into ``writer``'s staging directory; return ``(outcome, record)``.

    When the triplet file is wanted and the method produces one result, the
    records are streamed into the staging directory as they are computed,
    instead of being held on the result until the end: the file is
    byte-for-byte the one ``native.save`` writes, and a large chain no longer
    needs its whole triplet list in memory.  Anything that fails removes the
    staging directory with everything in it, and the output folder is unchanged.
    """
    stream = triplets and spec.keep_triplets and can_stream_triplets(spec)
    per_residue = sweeps_residues(spec)
    with writer:
        if writer.skipped:
            return None, None
        streamed = None
        if stream:
            partial = os.path.join(writer.path, _STREAMING_NAME)
            sink = native.NativeTripletStreamWriter(partial)
            try:
                outcome = execute_streaming(spec, sink.write)
            finally:
                sink.close()
            name = native.triplets_filename(prefix)
            os.replace(partial, os.path.join(writer.path, name))
            streamed = {"": name}
        else:
            outcome = execute(spec)
        writer.payload = _output.write_outputs(
            outcome,
            writer.path,
            prefix=prefix,
            per_residue=per_residue,
            triplets=triplets,
            command=command,
            argv=argv,
            request=request,
            metadata=metadata,
            source=source,
            streamed_triplets=streamed,
        )
    return outcome, writer.payload


def _task(row: _Request, options: BatchOptions, token: str, markers: str) -> Dict[str, Any]:
    return {
        "index": row.index,
        "spec": row.spec,
        "prefix": row.prefix,
        "request": row.request,
        "metadata": row.metadata,
        "source": {
            "manifest": os.path.basename(options.manifest),
            "rows": list(row.lines),
            "participant": row.participant or None,
            "columns": dict(row.columns) or None,
        },
        "output": options.output,
        "triplets": options.triplets,
        # Planning already applied --if-exists to what was on disk then. If the
        # calculation appears while this request computes, it belongs to
        # another run, so only an explicit overwrite may replace it.
        "if_exists": (_output.IF_EXISTS_OVERWRITE
                      if options.if_exists == _output.IF_EXISTS_OVERWRITE
                      else _output.IF_EXISTS_SKIP
                      if options.if_exists == _output.IF_EXISTS_SKIP
                      else _output.IF_EXISTS_FAIL),
        "command": options.command,
        "argv": options.argv,
        "token": token,
        "markers": markers,
    }


def run_task(task: Mapping[str, Any]) -> Dict[str, Any]:
    """Compute and install one request; return its small record.  Never raises.

    Runs in this process for ``jobs=1`` and in a worker process otherwise.  The
    result object lives only inside this call: it is written, summarised to a
    few numbers, and dropped.
    """
    markers = task.get("markers")
    if markers:
        _touch(os.path.join(markers, "%d.started" % task["index"]))
    record: Dict[str, Any] = {"status": STATUS_FAILED, "message": ""}
    try:
        request = task["request"]
        writer = _output.CalculationWriter(
            task["output"], task["prefix"], task["if_exists"], token=task["token"],
            verify=lambda prefix: _output.verify_existing(
                task["output"], prefix, request, task["triplets"]),
            shared_staging=True,
        )
        # An incomplete record found at planning is replaced, never kept.
        if task["if_exists"] == _output.IF_EXISTS_SKIP:
            state, _old = _output.existing_state(task["output"], task["prefix"])
            if state == "incomplete":
                writer.if_exists = _output.IF_EXISTS_OVERWRITE
        outcome, _payload = compute_into(
            task["spec"], writer,
            prefix=task["prefix"],
            triplets=task["triplets"],
            command=task["command"],
            argv=task["argv"],
            request=request,
            metadata=task["metadata"],
            source=task["source"],
        )
        if outcome is None:
            record.update(_skipped_record(
                _output.read_record(task["output"], task["prefix"])))
            return record
        first = next(iter(outcome.results.values()))
        n_triangles = sum(r.counts.n_triangles for r in outcome.results.values())
        # Partial is decided from each result's structured selection record
        # (and, for a sweep, from the residues it names as left out); the
        # warnings only supply the words.
        partial = any(
            selection_omissions(
                selection_sides(r.annotations.get("selection_request")), r.warnings)
            for r in outcome.results.values()
        ) or any(w.startswith(PARTIAL_WARNINGS) for w in outcome.warnings)
        notes = []
        for warnings in [outcome.warnings] + [r.warnings for r in outcome.results.values()]:
            for warning in warnings:
                if warning.startswith(("no triangles:",) + PARTIAL_WARNINGS
                                      + STRUCTURE_NOTES) \
                        and warning not in notes:
                    notes.append(warning)
        record.update(
            status=(STATUS_NO_TRIANGLES if not n_triangles
                    else STATUS_PARTIAL if partial else STATUS_OK),
            partial_selection=bool(partial),
            structure_id=first.structure_id.value,
            n_results=len(outcome.results),
            n_triangles=n_triangles,
            n_unique_keys=sum(r.counts.n_unique_keys for r in outcome.results.values()),
            message="; ".join(notes),
        )
        del outcome
    except CLIError as exc:
        record["message"] = exc.message
    except TSRError as exc:
        # A package refusal is the request being wrong, exactly as a CLIError
        # is, so it is reported verbatim with no exception type in front of it.
        # The type prefix below is reserved for "this is a defect in biotsr".
        record["message"] = str(exc)
    except Exception as exc:  # noqa: BLE001 -- attributed, never swallowed
        record["message"] = "%s: %s" % (type(exc).__name__, exc)
    if markers:
        _write_marker(os.path.join(markers, "%d.done" % task["index"]), record)
    return record


def _touch(path: str) -> None:
    try:
        with open(path, "w"):
            pass
    except OSError:
        pass


def _write_marker(path: str, record: Mapping[str, Any]) -> None:
    try:
        with open(path + ".tmp", "w") as fh:
            json.dump(dict(record), fh)
        os.replace(path + ".tmp", path)
    except OSError:
        pass


def _read_marker(path: str) -> Optional[Dict[str, Any]]:
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


class _SigtermAsInterrupt(object):
    """Treat SIGTERM -- what a scheduler sends at a time limit -- like Ctrl-C."""

    def __enter__(self):
        self.previous = None
        try:
            self.previous = signal.signal(signal.SIGTERM, self._handler)
        except (ValueError, OSError, AttributeError):  # not the main thread
            self.previous = None
        return self

    @staticmethod
    def _handler(signum, frame):
        raise KeyboardInterrupt()

    def __exit__(self, *exc):
        if self.previous is not None:
            signal.signal(signal.SIGTERM, self.previous)
        return False


def run_batch(options: BatchOptions, console=None) -> BatchReport:
    """Plan and run a dataset.  Writes the summary unless ``dry_run``.

    Raises :class:`CLIError` for a dataset that cannot be run at all, and
    :class:`BatchInterrupted` (after writing the summary) on Ctrl-C or SIGTERM
    at any point after planning: while structures download, while existing
    outputs are checked, or while requests compute.  With ``on_error="stop"``
    the first failure stops the batch and is raised as a :class:`CLIError` with
    the partial-batch exit code, also after the summary is written.
    """
    console = console or _SilentConsole()
    validate_batch_options(options)
    with _SigtermAsInterrupt():
        plan = plan_batch(options)
        sources = len({n for r in plan.rows for n in r.lines})
        console.say("manifest   %s (%d source row%s, %d request%s)" % (
            options.manifest, sources, "" if sources == 1 else "s",
            len(plan.rows), "" if len(plan.rows) == 1 else "s"))
        if plan.unused_columns and not options.dry_run:
            console.say("unused     %s (read by none of the selected methods)"
                        % ", ".join(plan.unused_columns))
        if options.dry_run:
            return _dry_run(plan, console)
        if options.jobs > 1:
            # A missing dependency is the installation, not a request: refuse
            # it before anything is downloaded, computed or locked.
            from .parallel import require_joblib

            require_joblib()
        with _OutputLock(options.output) as lock:
            plan.lock = lock
            return _run_planned(plan, console)


def validate_batch_options(options: BatchOptions) -> None:
    """Refuse option values that cannot be honoured, before anything is read.

    The same rules ``biotsr fetch`` and :func:`biotsr.fetch_structures` apply to
    the download settings, so a value one entry point refuses is never silently
    replaced by a default in another.
    """
    if options.if_exists not in _output.IF_EXISTS_CHOICES:
        raise CLIError(
            "if_exists=%r is not one of %s"
            % (options.if_exists, ", ".join(_output.IF_EXISTS_CHOICES))
        )
    if options.on_error not in (ON_ERROR_CONTINUE, ON_ERROR_STOP):
        raise CLIError(
            "on_error=%r is not one of %s, %s"
            % (options.on_error, ON_ERROR_CONTINUE, ON_ERROR_STOP)
        )
    if options.inapplicable not in INAPPLICABLE_CHOICES:
        raise CLIError("inapplicable=%r is not one of %s"
                       % (options.inapplicable, ", ".join(INAPPLICABLE_CHOICES)))
    if isinstance(options.jobs, bool) or not isinstance(options.jobs, int) \
            or options.jobs < 1:
        raise CLIError("jobs must be a whole number of at least 1 (got %r)"
                       % (options.jobs,))
    if options.sheet and options.defaults.get("method"):
        raise CLIError("a sheet takes its methods from --method")
    for method in options.methods or ():
        try:
            get_method(method)
        except KeyError:
            from ..api import list_methods

            raise CLIError("unknown method %r" % method,
                           "supported methods: %s" % ", ".join(list_methods()))
    try:
        validate_download_settings(
            workers=options.download_workers,
            timeout=options.timeout,
            retries=options.retries,
        )
    except ValueError as exc:
        raise CLIError(str(exc), _DOWNLOAD_HINT)


_DOWNLOAD_HINT = (
    "--download-workers is a whole number from 1 to %d, --timeout a number of "
    "seconds greater than 0, and --retries a whole number of 0 or more"
    % MAX_WORKERS
)


def _run_planned(plan: BatchPlan, console) -> BatchReport:
    options = plan.options
    token = plan.lock.token
    removed, foreign = _output.remove_stale_leftovers(options.output)
    for path in removed:
        console.detail("removed    %s (left by a process that no longer exists)"
                       % path)
    if foreign:
        console.error(
            "%d staging name%s in %s belong to another machine and were left "
            "alone: %s" % (len(foreign), "" if len(foreign) == 1 else "s",
                           options.output, ", ".join(os.path.basename(f)
                                                     for f in foreign[:5]))
        )

    fetch_report_path = ""
    try:
        if options.on_error == ON_ERROR_STOP:
            _stop_if_failed(plan, console)

        fetch_report_path = _download(plan, console)
        if options.on_error == ON_ERROR_STOP:
            _stop_if_failed(plan, console)

        for row in plan.rows:
            if row.action == ACTION_RUN:
                _resolve_existing(row, options)
        if options.on_error == ON_ERROR_STOP:
            _stop_if_failed(plan, console)
    except _DownloadInterrupted as exc:
        _interrupted_before_compute(plan, console, exc.report_path,
                                    "stopped while structures were downloading")
    except KeyboardInterrupt:
        _interrupted_before_compute(plan, console, fetch_report_path,
                                    "stopped before computing began")

    for row in plan.rows:
        if row.action == ACTION_FAIL:
            console.error("%s: %s" % (_where(row), row.record["message"]))
        elif row.action == ACTION_SKIP:
            console.detail("skipped %s" % row.prefix)
            if row.record.get("kept_partial"):
                console.warn("%s: %s" % (_where(row), row.record["message"]))
        elif row.action == ACTION_INAPPLICABLE:
            console.detail("inapplicable %s %s: %s"
                           % (_where(row), row.method, row.record["message"]))
        elif row.action == ACTION_NO_CHAIN:
            console.warn("%s: %s" % (_where(row), row.record["message"]))

    pending = [row for row in plan.rows if row.action == ACTION_RUN]
    interrupted = False
    pool_error = ""
    try:
        if options.jobs > 1 and len(pending) > 1:
            interrupted, pool_error = _execute_parallel(pending, plan, token, console)
        else:
            interrupted = _execute_serial(pending, plan, token, console)
    finally:
        _clean_staging(plan, token)
    report = _write_reports(plan, console)
    report.fetch_report_path = fetch_report_path
    report.pool_error = pool_error
    report.removed_leftovers = tuple(removed)
    report.foreign_leftovers = tuple(foreign)
    report.unused_columns = plan.unused_columns
    if interrupted:
        report.interrupted = True
        raise BatchInterrupted(report)
    if options.on_error == ON_ERROR_STOP and report.failures:
        first = report.failures[0]
        raise CLIError(
            "the batch stopped at its first failed request (request %s, rows %s); "
            "the summary lists every request" % (first["request"], first["rows"]),
            "rerun with --if-exists skip to keep the calculations already written",
            exit_code=EXIT_PARTIAL,
        )
    return report


def _where(row: _Request) -> str:
    """``sheet.csv line 5 (carbohydrate)`` -- where a request came from."""
    base = row.where if len(row.lines) == 1 else "%s (lines %s)" % (
        row.where.rsplit(" line ", 1)[0], _lines(row.lines))
    if row.participant == "2":
        return "%s (%s, participant 2)" % (base, row.method)
    return "%s (%s)" % (base, row.method)


class _DownloadInterrupted(KeyboardInterrupt):
    """Ctrl-C or SIGTERM while structures were downloading."""

    def __init__(self, report_path: str):
        super().__init__("interrupted")
        self.report_path = report_path


def _interrupted_before_compute(plan: BatchPlan, console, fetch_report_path: str,
                                message: str) -> None:
    """Record an interruption that came before any request computed, then raise."""
    for row in plan.rows:
        if row.action in (ACTION_RUN, ACTION_DOWNLOAD):
            row.record["status"] = STATUS_NOT_STARTED
            row.record["message"] = message
    report = _write_reports(plan, console)
    report.fetch_report_path = fetch_report_path
    report.interrupted = True
    raise BatchInterrupted(report)


class _OutputLock(object):
    """One batch at a time per ``--output``, so batch reports are never clobbered.

    Calculations are safe to share between runs -- each is staged under a name
    unique to its invocation and installed under its own lock, refusing to
    replace another run's work -- but ``batch_summary.csv`` and
    ``batch_failures.csv`` are one file per output folder.  A second batch into
    a folder a live batch is writing is refused before it computes anything.  A
    lock left by a process that no longer exists on this machine is taken over;
    one from another machine cannot be checked, so the message says how to
    clear it.
    """

    NAME = ".biotsr-batch.lock"

    def __init__(self, output: str):
        self.path = os.path.join(output, self.NAME)
        self.token = _output.invocation_token()
        self._held = False

    def __enter__(self) -> "_OutputLock":
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        record = {"pid": os.getpid(), "host": _host(), "token": self.token}
        for _attempt in range(2):
            try:
                fh = open(self.path, "x")       # exclusive: fails if it exists
            except FileExistsError:
                holder = self._holder()
                if holder is not None and not _alive_here(holder):
                    try:
                        os.remove(self.path)
                    except FileNotFoundError:
                        pass
                    continue
                raise CLIError(
                    "another batch is writing into %s (%s)"
                    % (os.path.dirname(self.path), _describe(holder)),
                    "its batch_summary.csv would overwrite this one's; give each "
                    "batch its own --output, or delete %s if no batch is running"
                    % self.path,
                )
            try:
                with fh:
                    json.dump(record, fh)
            except BaseException:
                # Interrupted between creating the lock and recording its
                # owner: remove it, or the next batch would find a lock it
                # cannot attribute to anyone.
                try:
                    os.remove(self.path)
                except OSError:
                    pass
                raise
            self._held = True
            return self
        raise CLIError("could not take the batch lock %s" % self.path,
                       "delete it if no batch is running")

    def __exit__(self, *exc):
        if self._held and (self._holder() or {}).get("token") == self.token:
            try:
                os.remove(self.path)
            except FileNotFoundError:
                pass
        return False

    def _holder(self) -> Optional[Dict[str, Any]]:
        try:
            with open(self.path) as fh:
                record = json.load(fh)
        except (OSError, ValueError):
            return None
        return record if isinstance(record, dict) else None


def _host() -> str:
    return os.uname().nodename if hasattr(os, "uname") else ""


def _alive_here(holder: Mapping[str, Any]) -> bool:
    """Whether the process that took a lock may still be running."""
    if holder.get("host") != _host():
        return True                         # another machine: cannot tell
    try:
        os.kill(int(holder.get("pid")), 0)
    except ProcessLookupError:
        return False
    except (PermissionError, TypeError, ValueError, OverflowError):
        return True
    return True


def _describe(holder: Optional[Mapping[str, Any]]) -> str:
    if not holder:
        return "lock file unreadable"
    return "pid %s on %s" % (holder.get("pid"), holder.get("host") or "unknown host")


def _stop_if_failed(plan: BatchPlan, console) -> None:
    failed = [row for row in plan.rows if row.action == ACTION_FAIL]
    if not failed:
        return
    for row in plan.rows:
        if row.action not in (ACTION_FAIL, ACTION_INAPPLICABLE, ACTION_NO_CHAIN):
            row.record["status"] = STATUS_NOT_STARTED
            row.record["message"] = ""
    _write_reports(plan, console)
    raise CLIError(
        "%s: %s" % (_where(failed[0]), failed[0].record["message"]),
        "--on-error stop halts before computing anything when a request cannot "
        "be run; the summary lists every request",
    )


def _record_result(row: _Request, record: Mapping[str, Any], console) -> None:
    row.record.update(record)
    if row.record["status"] in (STATUS_OK, STATUS_PARTIAL):
        console.detail("%-7s %-44s %8d triangles"
                       % (row.record["status"], row.prefix, row.record["n_triangles"]))
        if row.record.get("message"):
            console.warn("%s: %s" % (_where(row), row.record["message"]))
    elif row.record["status"] == STATUS_NO_TRIANGLES:
        console.warn("%s: %s" % (_where(row), row.record["message"]))
    elif row.record["status"] == STATUS_SKIPPED:
        console.detail("skipped %s" % row.prefix)
        if row.record.get("kept_partial"):
            console.warn("%s: %s" % (_where(row), row.record["message"]))
    else:
        console.error("%s: %s" % (_where(row), row.record["message"]))


def _execute_serial(pending: List[_Request], plan: BatchPlan, token: str,
                    console) -> bool:
    options = plan.options
    index = 0
    try:
        for index, row in enumerate(pending):
            record = run_task(_task(row, options, token, ""))
            _record_result(row, record, console)
            if record["status"] == STATUS_FAILED and options.on_error == ON_ERROR_STOP:
                for later in pending[index + 1:]:
                    later.record["status"] = STATUS_NOT_STARTED
                break
    except KeyboardInterrupt:
        pending[index].record["status"] = STATUS_INTERRUPTED
        pending[index].record["message"] = "stopped while running"
        for later in pending[index + 1:]:
            later.record["status"] = STATUS_NOT_STARTED
        return True
    return False


def _execute_parallel(pending: List[_Request], plan: BatchPlan, token: str,
                      console) -> Tuple[bool, str]:
    """Run ``pending`` in worker processes.  Returns ``(interrupted, pool_error)``.

    Whatever happens to the pool, every request ends with a truthful status: a
    request whose record came back is reported as it was; one whose worker
    installed its result but whose record was not received yet is reported
    from the marker the worker left; one that had started and never finished is
    reported as interrupted (Ctrl-C, SIGTERM) or failed (the pool broke); one
    that never started is ``not_started``.  Nothing is guessed to have
    succeeded.
    """
    from .parallel import parallel_map

    options = plan.options
    os.makedirs(options.output, exist_ok=True)
    markers = tempfile.mkdtemp(prefix=".biotsr-batch-%s-" % token, dir=options.output)
    # More workers than requests would only start processes that never get work.
    workers = min(options.jobs, len(pending))
    console.say("workers    %d processes" % workers)
    received = 0
    results = None
    try:
        tasks = (_task(row, options, token, markers) for row in pending)
        results = parallel_map(run_task, tasks, workers)
        for row, record in zip(pending, results):
            received += 1
            _record_result(row, record, console)
            if record["status"] == STATUS_FAILED and options.on_error == ON_ERROR_STOP:
                break
        if received < len(pending):
            _close(results)
            _classify_unreceived(pending[received:], markers,
                                 "stopped after an earlier request failed",
                                 STATUS_INTERRUPTED)
        return False, ""
    except KeyboardInterrupt:
        _close(results)
        _classify_unreceived(pending[received:], markers, "stopped while running",
                             STATUS_INTERRUPTED)
        return True, ""
    except Exception as exc:  # noqa: BLE001 -- the pool broke; requests say how
        # joblib raises TerminatedWorkerError when a worker process dies -- the
        # operating system's out-of-memory killer, a signal, a crash in native
        # code -- and other exceptions when the executor itself fails. Either
        # way the requests already finished are kept and reported, and the
        # ones that were running when it happened are failed, not guessed.
        _close(results)
        reason = "%s: %s" % (type(exc).__name__, _first_line(str(exc)))
        _classify_unreceived(
            pending[received:], markers,
            "the worker pool failed while this request was running (%s); this "
            "request, or another one running at the same time, caused it. Rerun "
            "with --if-exists skip to keep every finished calculation" % reason,
            STATUS_FAILED,
        )
        console.error("worker pool failed: %s" % reason)
        return False, reason
    finally:
        shutil.rmtree(markers, ignore_errors=True)


def _first_line(text: str) -> str:
    text = (text or "").strip()
    return text.splitlines()[0] if text else ""


def _close(results) -> None:
    close = getattr(results, "close", None)
    if close is not None:
        try:
            close()
        except Exception:  # noqa: BLE001 -- shutting down; the report matters more
            pass


def _classify_unreceived(rows: List[_Request], markers: str, running_message: str,
                         running_status: str) -> None:
    """Requests whose record never reached this process: done, running or unstarted.

    Results come back in plan order, so a request can finish before an earlier
    one is reported. The worker's marker file says what actually happened.
    """
    for row in rows:
        done = _read_marker(os.path.join(markers, "%d.done" % row.index))
        if done is not None:
            row.record.update(done)
        elif os.path.exists(os.path.join(markers, "%d.started" % row.index)):
            row.record["status"] = running_status
            row.record["message"] = running_message
        else:
            row.record["status"] = STATUS_NOT_STARTED
            row.record["message"] = ""


def _clean_staging(plan: BatchPlan, token: str) -> None:
    """Remove this invocation's own staging directory, and nothing else."""
    path = os.path.join(plan.options.output, ".biotsr-staging-" + token)
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path, ignore_errors=True)


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------


def _download(plan: BatchPlan, console) -> str:
    options = plan.options
    wanted = [row for row in plan.rows if row.action == ACTION_RUN and row.fetch_id]
    if not wanted:
        return ""
    ids = []
    for row in wanted:
        if row.fetch_id not in ids:
            ids.append(row.fetch_id)
    console.say("download   %d structure(s) into %s (%d at a time)"
                % (len(ids), options.pdb_dir, min(options.download_workers, len(ids))))
    try:
        report = fetch_structures(
            ids, options.pdb_dir,
            workers=options.download_workers,
            timeout=options.timeout,
            retries=options.retries,
            progress=lambda r: console.detail("fetch   %-6s %s" % (r.pdb_id, r.status)),
        )
    except FetchInterrupted as exc:
        os.makedirs(options.output, exist_ok=True)
        path = exc.report.write_csv(os.path.join(options.output, FETCH_REPORT_NAME))
        raise _DownloadInterrupted(path) from None
    os.makedirs(options.output, exist_ok=True)
    path = report.write_csv(os.path.join(options.output, FETCH_REPORT_NAME))
    by_id = report.by_id()
    for row in wanted:
        fetched = by_id.get(row.fetch_id)
        if fetched is None or not fetched.usable:
            _fail(row, "%s: structure %s could not be downloaded: %s"
                  % (row.where, row.fetch_id,
                     fetched.message if fetched else "no record"))
    counts = report.counts()
    console.say("download   %s; report %s" % (
        ", ".join("%d %s" % (n, s) for s, n in sorted(counts.items())), path))
    return path


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------


def _write_csv_atomically(path: str, fieldnames, rows) -> None:
    directory = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(prefix="." + os.path.basename(path) + ".", dir=directory)
    try:
        with os.fdopen(fd, "w", newline="") as fh:
            # Records also carry fields the summary does not show (the label).
            writer = csv.DictWriter(fh, fieldnames=list(fieldnames),
                                    extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        publish_file_mode(tmp)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise


def _write_reports(plan: BatchPlan, console) -> BatchReport:
    options = plan.options
    fieldnames = list(SUMMARY_COLUMNS) + list(plan.metadata_columns)
    records = []
    for row in plan.rows:
        record = dict(row.record)
        if row.note and row.note not in record.get("message", ""):
            record["message"] = "; ".join(filter(None, (record.get("message"),
                                                        row.note)))
        for name in plan.metadata_columns:
            record[name] = row.metadata.get(name, "")
        records.append(record)
    os.makedirs(options.output, exist_ok=True)
    summary = os.path.join(options.output, SUMMARY_NAME)
    _write_csv_atomically(summary, fieldnames, records)
    failed = [r for r in records if r["status"] == STATUS_FAILED]
    failures_path = os.path.join(options.output, FAILURES_NAME)
    if failed:
        _write_csv_atomically(failures_path, fieldnames, failed)
    elif os.path.exists(failures_path):
        os.remove(failures_path)
    return BatchReport(records=records, summary_path=summary,
                       failures_path=failures_path if failed else "",
                       unused_columns=plan.unused_columns)


# ---------------------------------------------------------------------------
# Dry run
# ---------------------------------------------------------------------------

_PREVIEW = "  %-4s %-8s %-28s %-3s %-7s %-34s %s"


def _dry_run(plan: BatchPlan, console) -> BatchReport:
    """Report what a run would do.  Downloads nothing, creates nothing.

    Structures already present are hashed (read-only) so existing outputs can be
    checked against ``--if-exists``; for structures that would be downloaded,
    that check waits for the real run.  Every planned request is listed: its
    source row(s), method, participant, chains, selection, output prefix and
    what would happen to it.
    """
    options = plan.options
    for row in plan.rows:
        if row.action == ACTION_RUN:
            _resolve_existing(row, options)
    rows = plan.rows
    from_row = sum(1 for r in rows if r.chain_source == CHAIN_FROM_ROW)
    from_default = sum(1 for r in rows if r.chain_source == CHAIN_FROM_DEFAULT)
    present = sum(1 for r in rows if r.spec is not None and not r.fetch_id)
    to_download = sorted({r.fetch_id for r in rows if r.action == ACTION_DOWNLOAD})
    console.out("dry run    nothing will be computed, downloaded or written")
    console.out("requests   %d from %d source row(s)"
                % (len(rows), len({n for r in rows for n in r.lines})))
    console.out("chains     %d from the row, %d from --chain" % (from_row, from_default))
    if options.columns:
        console.out("columns    %s" % ", ".join(
            "%s <- %s" % (k, v) for k, v in sorted(options.columns.items())))
    if plan.metadata_columns:
        console.out("metadata   %s (kept, not scientific)" % ", ".join(plan.metadata_columns))
    if plan.unused_columns:
        console.out("unused     %s (read by none of the selected methods)"
                    % ", ".join(plan.unused_columns))
    console.out("structures %d present, %d distinct would be downloaded"
                % (present, len(to_download)))
    stale, foreign = _output.stale_leftovers(options.output)
    if stale or foreign:
        console.out("leftovers  %d staging name(s) of ended runs would be removed; "
                    "%d from another machine would be left alone"
                    % (len(stale), len(foreign)))
    console.out("outcome    %d would run, %d would be downloaded then run, "
                "%d would be skipped (verified), %d inapplicable, %d chain not "
                "specified, %d would fail"
                % (plan.count(ACTION_RUN), plan.count(ACTION_DOWNLOAD),
                   plan.count(ACTION_SKIP), plan.count(ACTION_INAPPLICABLE),
                   plan.count(ACTION_NO_CHAIN), plan.count(ACTION_FAIL)))
    console.out(_PREVIEW % ("req", "rows", "method", "pt", "chains", "selection",
                            "prefix / outcome"))
    outcome_text = {ACTION_RUN: "would run", ACTION_DOWNLOAD: "would download, run",
                    ACTION_SKIP: "would skip (verified)",
                    ACTION_INAPPLICABLE: "inapplicable", ACTION_FAIL: "would fail",
                    ACTION_NO_CHAIN: "chain not specified"}
    for row in rows:
        chains = row.record.get("chain") or "-"
        if row.record.get("chain_2"):
            chains += "+" + row.record["chain_2"]
        selection = row.record.get("selection") or "-"
        if len(selection) > 34:
            selection = selection[:31] + "..."
        tail = row.prefix or "-"
        status = outcome_text.get(row.action, row.action)
        if row.action in (ACTION_FAIL, ACTION_INAPPLICABLE, ACTION_NO_CHAIN):
            status = "%s: %s" % (status, row.record["message"])
        console.out(_PREVIEW % (row.index, ";".join(str(n) for n in row.lines)[:8],
                                row.method[:28], row.participant or "-", chains[:7],
                                selection, "%s  [%s]" % (tail, status)))
        if row.spec is not None:
            for text in _side_texts(row.spec):
                console.out("         %s" % text)
        if row.note:
            console.out("         %s" % row.note)
    console.out(
        "deferred   checked when each structure is read: that the chains exist in "
        "it and that every selection matches its residues, molecules or atoms%s"
        % ("; for structures not yet downloaded, also their content and whether "
           "existing outputs can be reused" if to_download else "")
    )
    records = []
    for row in rows:
        record = dict(row.record)
        if row.note:
            record["message"] = "; ".join(filter(None, (record.get("message"),
                                                        row.note)))
        record["status"] = {ACTION_RUN: "would_run", ACTION_DOWNLOAD: "would_download",
                            ACTION_SKIP: "would_skip",
                            ACTION_INAPPLICABLE: STATUS_INAPPLICABLE,
                            ACTION_NO_CHAIN: STATUS_CHAIN_NOT_SPECIFIED,
                            }.get(row.action, STATUS_FAILED)
        records.append(record)
    return BatchReport(records=records, unused_columns=plan.unused_columns)


class _SilentConsole(object):
    def say(self, text):
        pass

    detail = error = hint = out = warn = say
