"""Run a dataset manifest or wide sheet from Python: the same engine ``biotsr batch`` uses.

::

    from biotsr import batch

    report = batch.run(
        "sample_details.csv",
        "results/",
        pdb_dir="pdb/",
        methods=["protein"],               # every row, once per method
        columns={"label": "group"},        # read the label from column "group"
        keep_columns=["group1"],           # carried into the records, never scientific
        triplets=False,
        jobs=4,                            # default 1: one request at a time
    )
    print(report.n_ok, report.n_failed, report.summary_path)

    # a wide TSR dataset sheet: each method reads its own columns
    report = batch.run("sheet.csv", "results/", pdb_dir="pdb/", sheet="tsr",
                       methods=["carbohydrate", "aa_ind"], inapplicable="skip")

Every rule of the command line applies unchanged: a blank cell is "not set",
an unrecognised manifest column is an error unless it is mapped (``columns``),
kept (``keep_columns``) or dropped (``ignore_columns``), every request must
resolve its chain(s) -- a chain is never inferred -- and two different requests
that would write the same output prefix are refused before anything is
computed.  ``methods`` runs each listed method on every row (not together with
a manifest ``method`` column); ``sheet`` (``"tsr"`` or ``"custom"``) reads a
wide sheet method-aware, with ``maps`` entries ``(method, participant, field,
column, description)`` as ``--map`` makes them, ``participant`` (``"1"``,
``"2"``, ``"all"``) and ``inapplicable`` (``"fail"`` or ``"skip"``).  See
``docs/CLI_GUIDE.md`` for both formats.

Problems with the request as a whole raise :class:`BatchError` (the command
line's exit code 2); a row that fails on its own is recorded in the report and
in ``batch_failures.csv``, and does not stop the others unless
``on_error="stop"``.
"""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Union

from .cli.batch import (
    FAILURES_NAME,
    FETCH_REPORT_NAME,
    ON_ERROR_CONTINUE,
    ON_ERROR_STOP,
    STATUS_FAILED,
    STATUS_INAPPLICABLE,
    STATUS_INTERRUPTED,
    STATUS_NOT_STARTED,
    STATUS_NO_TRIANGLES,
    STATUS_PARTIAL,
    STATUS_OK,
    STATUS_SKIPPED,
    SUMMARY_COLUMNS,
    SUMMARY_NAME,
    BatchInterrupted,
    BatchOptions,
    BatchReport,
    run_batch,
)
from .cli.errors import CLIError as BatchError
from .cli.manifest import read_manifest as _read_manifest
from .cli.output import IF_EXISTS_FAIL
from .cli.parallel import available_cpus, resolve_jobs
from .retrieval import DEFAULT_RETRIES, DEFAULT_TIMEOUT, DEFAULT_WORKERS

__all__ = [
    "run",
    "read_manifest",
    "available_cpus",
    "BatchError",
    "BatchInterrupted",
    "BatchReport",
    "SUMMARY_COLUMNS",
    "SUMMARY_NAME",
    "FAILURES_NAME",
    "FETCH_REPORT_NAME",
    "STATUS_OK",
    "STATUS_NO_TRIANGLES",
    "STATUS_PARTIAL",
    "STATUS_SKIPPED",
    "STATUS_INAPPLICABLE",
    "STATUS_FAILED",
    "STATUS_INTERRUPTED",
    "STATUS_NOT_STARTED",
    "ON_ERROR_CONTINUE",
    "ON_ERROR_STOP",
]


def read_manifest(
    path: str,
    *,
    columns: Optional[Mapping[str, str]] = None,
    keep_columns: Sequence[str] = (),
    ignore_columns: Sequence[str] = (),
    defaults: Optional[Mapping[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """Every row of a manifest, resolved but not run.

    Each row is a mapping of canonical field names to values, plus ``line`` (its
    line number), ``metadata`` (the kept columns, verbatim) and ``chain_source``
    (``"row"`` or ``"default"``).
    """
    rows = _read_manifest(path, defaults=defaults, columns=columns,
                          keep_columns=keep_columns, ignore_columns=ignore_columns)
    out = []
    for row in rows:
        public = {k: v for k, v in row.items() if not k.startswith("_")}
        public["line"] = row["_row"]
        public["metadata"] = dict(row["_metadata"])
        public["chain_source"] = row["_chain_source"]
        out.append(public)
    return out


def run(
    manifest: str,
    output: str,
    *,
    pdb_dir: str = "",
    methods: Sequence[str] = (),
    defaults: Optional[Mapping[str, Any]] = None,
    columns: Optional[Mapping[str, str]] = None,
    keep_columns: Sequence[str] = (),
    ignore_columns: Sequence[str] = (),
    sheet: Optional[str] = None,
    maps: Sequence[Any] = (),
    cross_sides: Union[Sequence[str], Mapping[str, str], None] = None,
    participant: str = "1",
    inapplicable: str = "fail",
    triplets: bool = True,
    if_exists: str = IF_EXISTS_FAIL,
    on_error: str = ON_ERROR_CONTINUE,
    jobs: Union[int, str] = 1,
    download_missing: bool = False,
    download_workers: int = DEFAULT_WORKERS,
    timeout: float = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    dry_run: bool = False,
) -> BatchReport:
    """Plan and run ``manifest`` into ``output``; returns a :class:`BatchReport`.

    ``jobs`` is a worker count or ``"auto"`` (the CPUs this process may use,
    honouring scheduler allocations); anything above 1 runs rows in worker
    processes through joblib, which BioTSR installs.
    ``download_missing`` fetches absent ``pdb_id`` structures into ``pdb_dir``
    first (see :func:`biotsr.fetch_structures`); ``download_workers``,
    ``timeout`` and ``retries`` are validated exactly as ``biotsr fetch``
    validates them, and an invalid value raises :class:`BatchError` rather than
    being replaced by a default.  ``dry_run`` plans without
    downloading, computing or writing anything; its report's statuses are
    ``would_run``, ``would_download``, ``would_skip``, ``inapplicable`` and
    ``failed``.  ``maps`` may also be strings in the ``--map`` spelling
    (``"carbohydrate.chain=glycan_chain"``).  ``cross_sides`` chooses how a
    cross method's side is selected from a ``tsr`` sheet, as ``--cross-side``
    does: ``{"drug_protein_cross.protein": "aa_ind"}`` or the strings
    ``["drug_protein_cross.protein=aa_ind"]``; ``"whole"`` uses the whole side.
    """
    from .cli.sheet import parse_cross_side_options, parse_map_options

    entries = []
    for item in maps or ():
        entries.extend(parse_map_options([item]) if isinstance(item, str) else [tuple(item)])
    options = BatchOptions(
        manifest=manifest,
        output=output,
        pdb_dir=pdb_dir or "",
        methods=tuple(methods or ()),
        defaults=dict(defaults or {}),
        columns=dict(columns or {}),
        keep_columns=tuple(keep_columns or ()),
        ignore_columns=tuple(ignore_columns or ()),
        sheet=sheet,
        maps=tuple(entries),
        cross_sides=parse_cross_side_options(
            ["%s=%s" % item for item in cross_sides.items()]
            if isinstance(cross_sides, Mapping) else list(cross_sides or ())),
        participant=participant,
        inapplicable=inapplicable,
        triplets=bool(triplets),
        if_exists=if_exists,
        on_error=on_error,
        jobs=resolve_jobs(jobs),
        download_missing=bool(download_missing),
        download_workers=download_workers,
        timeout=timeout,
        retries=retries,
        dry_run=bool(dry_run),
        command="python",
        argv=None,
    )
    return run_batch(options)
