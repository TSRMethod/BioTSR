"""The output-folder contract: one shared folder, one prefix per calculation.

One rule: **a run writes only inside ``--output``, and only files it owns.**
No file is ever written to the working directory, next to the input PDB, or
relative to the installation.  Every calculation's files sit side by side in
the output folder under one prefix (ADR-098)::

    <output>/
        7LF9_D_carbohydrate_NAG_1.keys
        7LF9_D_carbohydrate_NAG_1.meta.json
        7LF9_D_carbohydrate_NAG_1.triplets
        7LF9_D_carbohydrate_NAG_1.run.json
        1UBQ_A_amino_acid_atomic_perresidue_LEU_15.keys     (a per-residue sweep:
        ...                                                  one prefix per residue,
        1UBQ_A_amino_acid_atomic_perresidue.run.json         one record for all)
        batch_summary.csv

The prefix comes from the request (``biotsr.io.naming``), or from
``sample_id`` when the caller gives one.  ``PREFIX.run.json`` is the
calculation's **completion record**: the command, the arguments, the source
row(s), the full normalised request, and every file with its size.  It is the
one file the command line adds to the native result files, and the last one
written.

**Installation.**  A calculation is computed into a staging directory private
to the invocation (``.biotsr-staging-<token>/``).  Its files are then installed
into the folder under a per-prefix lock:

1. an *intent* record is written (``"complete": false``, listing the files),
2. each data file is moved into place with one ``os.replace``,
3. files of an earlier version of the same calculation that the new one does
   not have are removed -- only files that calculation's own record lists,
4. the completion record is rewritten with ``"complete": true`` and every
   file's size.

A run interrupted at any point leaves either no record, a staging directory
that the next run removes, or an intent record that every reader refuses;
never a partial result that looks complete.  A file at a target name that no
record of this calculation lists -- another program's file, or another
calculation's name differing only in letter case -- is never replaced or
deleted.

**Skipping is verified.**  ``--if-exists skip`` keeps an existing calculation
only when its completion record is complete, records exactly the new request,
was computed under the same processing rules, covers the requested files
(keys-only does not satisfy a request for triplets), and every file is present
at its recorded size.  An incomplete record left by an interrupted run is
recomputed.

**The key space is checked too** (ADR-107).  Every result entry of a completion record
holds its ``KeySpace`` fingerprint, and skip requires it to equal the one this
version computes for the request: a calculation written by a build with other
label tables or radices -- protein C-alpha results from before 1.0.0rc5, at
radix 20 -- is refused and left untouched; only ``--if-exists overwrite``
recomputes it.

**Processing rules are part of the request** (ADR-105).  The same request --
same structure, method, chain and selection -- can produce different numbers
when the rules that admit atoms change: the modified amino acids TPO, SEP and
PTR written as ``HETATM`` were excluded before ADR-101 and are read now, with
the key space and the package version unchanged.  The request record therefore
carries :func:`biotsr.api.processing_identity`, and a record
without it was written under rules this version cannot vouch for: it is never
skipped as if it were current, only recomputed (``--if-exists overwrite``).
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import secrets
import shutil
import socket
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Tuple

from ..api import processing_identity
from ..domain.selection_record import selection_omissions, selection_sides
from ..io import native
from .errors import CLIError

__all__ = [
    "IF_EXISTS_FAIL",
    "IF_EXISTS_OVERWRITE",
    "IF_EXISTS_SKIP",
    "IF_EXISTS_CHOICES",
    "RECORD_SUFFIX",
    "REQUEST_SCHEMA",
    "RESERVED_NAMES",
    "CalculationWriter",
    "check_complete",
    "existing_state",
    "find_records",
    "invocation_token",
    "read_record",
    "record_name",
    "request_record",
    "request_sha256",
    "stale_leftovers",
    "remove_stale_leftovers",
    "validate_prefix",
    "verify_existing",
    "write_outputs",
]

#: Refuse and change nothing.  The default: a research output folder that
#: silently absorbs a second, differently-parameterised run is how
#: irreproducible results happen.
IF_EXISTS_FAIL = "fail"
#: Replace the calculation's files.  Only files its own completion record lists
#: are replaced or removed.
IF_EXISTS_OVERWRITE = "overwrite"
#: Keep a verified existing calculation and report it as skipped; recompute one
#: an interrupted run left incomplete.  For resuming.
IF_EXISTS_SKIP = "skip"
IF_EXISTS_CHOICES = (IF_EXISTS_FAIL, IF_EXISTS_OVERWRITE, IF_EXISTS_SKIP)

#: ``PREFIX.run.json``: the completion record.
RECORD_SUFFIX = ".run.json"

#: Bumped when the meaning of a completion record changes.
#: 4 (1.0.0rc4): one shared folder with a completion record per calculation
#: (``PREFIX.run.json``); earlier versions wrote one directory per sample with
#: a ``run.json`` inside, which this version does not read.
REQUEST_SCHEMA = 4

#: Names an output folder uses for itself.  A prefix may not take one of them.
#: Compared case-insensitively, as the default macOS and Windows file systems do.
RESERVED_NAMES = (
    "analysis",
    "batch_summary.csv",
    "batch_failures.csv",
    "fetch_report.csv",
    ".biotsr-batch.lock",
)

#: The longest prefix accepted, in bytes of UTF-8: common file systems allow
#: 255 per name, and the longest suffix plus temporary-name decoration must fit.
MAX_PREFIX_BYTES = 200

#: Hidden names this package creates in an output folder, each carrying the
#: token of the invocation that owns it.
_STAGING = ".biotsr-staging-"
_TEMPORARY = ".biotsr-tmp-"
_LOCK = ".biotsr-lock-"


def validate_prefix(prefix: str, where: str = "") -> str:
    """Return ``prefix`` if it can name a calculation's files in the output folder.

    A prefix is a file-name **stem**, never a path: every separator, ``.`` and
    ``..``, a leading dot (hidden and staging names), control characters, the
    names the folder uses for itself and names too long for a file system are
    refused.
    """
    prefix_text = "%s: " % where if where else ""
    text = str(prefix)
    problem = None
    if not text.strip():
        problem = "it is empty"
    elif text != text.strip():
        problem = "it starts or ends with whitespace"
    elif os.path.isabs(text) or (len(text) > 1 and text[1] == ":"):
        problem = "it is an absolute path"
    elif "/" in text or "\\" in text:
        problem = "it contains a path separator"
    elif text in (".", ".."):
        problem = "it names the output folder or its parent"
    elif text.startswith("."):
        problem = "it starts with '.', which is reserved for hidden and staging names"
    elif any(ord(ch) < 32 or ord(ch) == 127 for ch in text):
        problem = "it contains a control character"
    elif text.casefold() in {r.casefold() for r in RESERVED_NAMES}:
        problem = "the output folder uses that name for itself"
    elif len(text.encode("utf-8")) > MAX_PREFIX_BYTES:
        problem = "it is longer than %d bytes" % MAX_PREFIX_BYTES
    if problem:
        raise CLIError(
            "%ssample id %r cannot be used as an output prefix: %s"
            % (prefix_text, text, problem),
            "a sample id becomes the shared prefix of one calculation's files "
            "in --output: use letters, digits, '.', '_', '+' and '-', not a path",
        )
    return text


def record_name(prefix: str) -> str:
    return prefix + RECORD_SUFFIX


# ---------------------------------------------------------------------------
# Tokens and leftovers
# ---------------------------------------------------------------------------


def host_tag() -> str:
    """Eight hex digits naming this machine, for the names of staging files.

    A hash rather than the host name, so a staging name says which machine owns
    it without spelling the machine's name into a shared tree.
    """
    name = os.uname().nodename if hasattr(os, "uname") else socket.gethostname()
    return hashlib.sha256(name.encode("utf-8")).hexdigest()[:8]


def invocation_token() -> str:
    """``<pid>-<host tag>-<random>``, unique to this invocation.

    The process and machine that own a staging directory are readable from its
    name, so a later run can tell a leftover of a process that no longer exists
    from the work of one that is still running.
    """
    return "%d-%s-%s" % (os.getpid(), host_tag(), secrets.token_hex(4))


_TOKEN = r"(?P<pid>\d+)-(?P<host>[0-9a-f]{8})-[0-9a-f]{8}"
_LEFTOVER_PATTERNS = (
    re.compile(r"^%s%s$" % (re.escape(_STAGING), _TOKEN)),
    re.compile(r"^%s%s-.+$" % (re.escape(_TEMPORARY), _TOKEN)),
    re.compile(r"^\.biotsr-batch-%s-[^/]*$" % _TOKEN),
)


def _owner_state(pid: int, tag: str) -> str:
    """``"dead"``, ``"alive"`` or ``"unknown"`` (another machine)."""
    if tag != host_tag():
        return "unknown"
    if pid == os.getpid():
        return "alive"
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return "dead"
    except (PermissionError, OverflowError, ValueError):
        return "alive"
    return "alive"


def _token_state(token: str) -> str:
    match = re.match("^%s$" % _TOKEN, str(token or ""))
    if not match:
        return "unknown"
    return _owner_state(int(match.group("pid")), match.group("host"))


def stale_leftovers(root: str):
    """``(stale, foreign)``: staging, temporary and marker names directly in ``root``.

    ``stale`` belong to an invocation that ran on this machine and whose process
    no longer exists -- an abrupt termination that gave it no chance to clean
    up.  ``foreign`` belong to another machine and cannot be checked.  Names of
    a live invocation are in neither list.  Only names this package creates are
    considered.  Read-only.
    """
    stale, foreign = [], []
    try:
        names = os.listdir(root)
    except OSError:
        return stale, foreign
    for name in sorted(names):
        if name.startswith(_LOCK):
            # A per-calculation lock names its owner inside, not in its name.
            path = os.path.join(root, name)
            try:
                with open(path) as fh:
                    state = _token_state(fh.read().strip())
            except OSError:
                continue
            if state == "dead":
                stale.append(path)
            elif state == "unknown":
                foreign.append(path)
            continue
        for pattern in _LEFTOVER_PATTERNS:
            match = pattern.match(name)
            if not match:
                continue
            path = os.path.join(root, name)
            if os.path.islink(path):
                break
            state = _owner_state(int(match.group("pid")), match.group("host"))
            if state == "dead":
                stale.append(path)
            elif state == "unknown":
                foreign.append(path)
            break
    return stale, foreign


def remove_stale_leftovers(root: str):
    """Remove :func:`stale_leftovers`; return ``(removed, foreign)``."""
    stale, foreign = stale_leftovers(root)
    removed = []
    for path in stale:
        if os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
        else:
            try:
                os.remove(path)
            except OSError:
                pass
        if not os.path.exists(path):
            removed.append(path)
    return removed, foreign


# ---------------------------------------------------------------------------
# Completion records
# ---------------------------------------------------------------------------


def read_record(root: str, prefix: str) -> Optional[Dict[str, Any]]:
    """The completion record of ``prefix``, or ``None`` when there is none."""
    path = os.path.join(root, record_name(prefix))
    if not os.path.lexists(path):
        return None
    if os.path.islink(path):
        return {"_problem": "is a symbolic link"}
    try:
        with open(path, encoding="utf-8") as fh:
            payload = json.load(fh)
    except (OSError, ValueError):
        return {"_problem": "cannot be read as JSON"}
    if not isinstance(payload, dict):
        return {"_problem": "is not a JSON object"}
    return payload


def _record_files(record: Mapping[str, Any]) -> List[str]:
    """Every file name a completion record owns, including itself."""
    names = set()
    for entry in record.get("results") or ():
        for name in (entry.get("files") or {}).values():
            names.add(str(name))
    for name in record.get("files_pending") or ():
        names.add(str(name))
    prefix = record.get("prefix")
    if prefix:
        names.add(record_name(str(prefix)))
    return sorted(n for n in names if n and os.path.basename(n) == n
                  and not n.startswith("."))


def existing_state(root: str, prefix: str) -> Tuple[str, Optional[Dict[str, Any]]]:
    """What is already at ``prefix`` in ``root``: ``(state, record)``.

    ``"absent"``      nothing of this calculation is there;
    ``"complete"``    a complete completion record of this exact prefix;
    ``"incomplete"``  an intent record: an interrupted earlier run;
    ``"foreign"``     files at this prefix's names that no record of this
                      prefix owns, a record of another prefix (a name
                      differing only in letter case), an unreadable record, or
                      one written by a different version.
    """
    record = read_record(root, prefix)
    if record is not None:
        if record.get("_problem") or record.get("prefix") != prefix:
            return "foreign", record
        if record.get("schema") != REQUEST_SCHEMA:
            return "foreign", record
        return ("complete" if record.get("complete") is True else "incomplete"), record
    for suffix in (native.KEYS_SUFFIX, native.METADATA_SUFFIX, native.TRIPLETS_SUFFIX):
        if os.path.lexists(os.path.join(root, prefix + suffix)):
            return "foreign", None
    return "absent", None


def describe_foreign(root: str, prefix: str, record: Optional[Mapping[str, Any]]) -> str:
    """Why a ``"foreign"`` state cannot be written over, for a message."""
    if record is None:
        return "files with its names exist but no BioTSR completion record lists them"
    if record.get("_problem"):
        return "its completion record %s %s" % (record_name(prefix), record["_problem"])
    if record.get("prefix") != prefix:
        return ("the completion record there belongs to %r, a name differing "
                "only in letter case" % (record.get("prefix"),))
    return ("its completion record was written by a different BioTSR version "
            "(schema %r; this one writes %d)" % (record.get("schema"), REQUEST_SCHEMA))


def find_records(root: str) -> List[str]:
    """Every completion record under ``root``, sorted by path.

    Recursive, so a tree of output folders can be read as one dataset; hidden
    directories (staging) and ``analysis`` folders are not entered.
    """
    out = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames
                             if not d.startswith(".") and d != "analysis")
        for name in filenames:
            if name.endswith(RECORD_SUFFIX) and not name.startswith("."):
                out.append(os.path.join(dirpath, name))
    return sorted(out)


def request_record(spec, structure_sha256: str, structure_id: str,
                   per_residue: bool = False) -> Dict[str, Any]:
    """The normalised calculation request, as the completion record holds it.

    Everything that decides which numbers a calculation contains, or which
    files hold them: the structure by content (its SHA-256, not its path, which
    can move), the identifier, the method, both chains, every selection, every
    modifier, every contact setting, and the processing rules that admit atoms
    (``processing``: the alternate-location and modified-residue rules, which
    can change while everything else stays the same).  Labels, source rows and kept columns
    are not in it -- they annotate a result without changing it -- and neither
    is whether triplets were written, which is compared separately because a
    calculation with triplets also satisfies a keys-only request.
    """
    return {
        "schema": REQUEST_SCHEMA,
        "method": spec.method,
        "structure_sha256": structure_sha256,
        "structure_id": structure_id,
        "chain": spec.chain,
        "per_residue": bool(per_residue),
        "modifiers": {k: spec.modifiers[k] for k in sorted(spec.modifiers)},
        "selection": {k: spec.selection[k] for k in sorted(spec.selection)},
        "cross": {k: spec.cross[k] for k in sorted(spec.cross)},
        "processing": processing_identity(spec.method),
    }


def _rules_text(rules: Mapping[str, Any]) -> str:
    return ", ".join("%s=%s" % (k, rules[k]) for k in sorted(rules))


def request_sha256(request: Mapping[str, Any]) -> str:
    blob = json.dumps(request, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def verify_existing(root: str, prefix: str, request: Mapping[str, Any],
                    triplets: bool) -> Optional[str]:
    """``None`` when ``prefix`` in ``root`` is a complete result of exactly ``request``.

    Otherwise the reason it is not, written for the person reading the batch
    report.  Read-only.
    """
    state, record = existing_state(root, prefix)
    if state == "absent":
        return "it does not exist"
    if state == "foreign":
        return describe_foreign(root, prefix, record)
    if state == "incomplete":
        return "an earlier run was interrupted before completing it"
    recorded = record.get("request")
    if not isinstance(recorded, Mapping):
        return "its completion record does not hold the request"
    # The processing rules first, with their own reason: the request can be
    # identical while the rules that admitted its atoms are not (ADR-105).
    # When the keys are in another key space as well, that is said first: it
    # is the more basic difference (ADR-113).
    wanted = request.get("processing")
    if wanted is not None:
        rules = recorded.get("processing")
        if not isinstance(rules, Mapping):
            return _also_key_space(record, request) + (
                "its completion record does not state the processing rules it "
                "was computed under (it was written before they were recorded, "
                "and outputs that old may have excluded the modified amino "
                "acids TPO, SEP and PTR written as HETATM), so it is not taken "
                "to match the current rules (%s)" % _rules_text(wanted))
        if dict(rules) != dict(wanted):
            return _also_key_space(record, request) + (
                "it was computed under different processing rules (%s; this "
                "version applies %s)" % (_rules_text(rules), _rules_text(wanted)))
    if dict(recorded) != json.loads(json.dumps(dict(request), default=str)):
        differing = sorted(
            key for key in set(recorded) | set(request)
            if json.dumps(recorded.get(key), sort_keys=True, default=str)
            != json.dumps(request.get(key), sort_keys=True, default=str)
        )
        return "it was produced by a different request (it differs in: %s)" % (
            ", ".join(differing))
    if triplets and not record.get("triplets"):
        return "it was written with --no-triplets, and triplets were requested"
    entries = record.get("results") or []
    if not entries:
        return "its completion record lists no results"
    # The key space last: the same request under the same processing rules can
    # still have been computed with different label tables or radices by
    # another build (ADR-107). Its keys would mean something else.
    problem = _key_space_problem(entries, request)
    if problem:
        return problem
    for entry in entries:
        files = entry.get("files") or {}
        if "keys" not in files or "metadata" not in files:
            return "its completion record lists no key or metadata file"
        if triplets and "triplets" not in files:
            return "it has no triplet file, and triplets were requested"
        problem = _check_files(root, files.values(), entry.get("file_bytes") or {})
        if problem:
            return problem
    return None


def _key_space_problem(entries, request: Mapping[str, Any]) -> Optional[str]:
    """Why ``entries``' keys are not in the key space this version gives
    ``request``, or ``None`` when they are."""
    expected = expected_key_space_fingerprint(request)
    for entry in entries:
        recorded_space = entry.get("key_space_fingerprint")
        if expected is not None and recorded_space != expected:
            return ("its keys were computed in a different key space (fingerprint "
                    "%s; this version computes %s for this request), so they do "
                    "not mean what this version's keys mean"
                    % (str(recorded_space)[:12] if recorded_space else "not recorded",
                       expected[:12]))
    return None


def _also_key_space(record: Mapping[str, Any], request: Mapping[str, Any]) -> str:
    """The key-space reason followed by ``"; and "``, or ``""``."""
    problem = _key_space_problem(record.get("results") or [], request)
    return problem + "; and " if problem else ""


def expected_key_space_fingerprint(request: Mapping[str, Any]) -> Optional[str]:
    """The ``KeySpace`` fingerprint this version gives ``request``'s results.

    Derived from the method's standard profile and the request's modifiers,
    exactly as a calculation resolves them; ``None`` for a request naming a
    method this version does not have.
    """
    from ..methods.registry import get_method

    try:
        recipe = get_method(str(request.get("method")))
    except (KeyError, ValueError):
        return None
    modifiers = request.get("modifiers") or {}
    profile = recipe.default_profile.with_modifiers(
        aa_grouping=bool(modifiers.get("aa_grouping")),
        mirror=bool(modifiers.get("mirror")),
        size_filter=modifiers.get("size_filter"),
        method=recipe.name,
    )
    return profile.key_space().fingerprint


def _check_files(root: str, names: Iterable[str], sizes: Mapping[str, int]) -> Optional[str]:
    for name in names:
        if os.path.basename(name) != name or name.startswith("."):
            return "its completion record names a file outside the folder: %r" % name
        full = os.path.join(root, name)
        if os.path.islink(full) or not os.path.isfile(full):
            return "its file %s is missing" % name
        if name not in sizes:
            return "its file %s has no recorded size" % name
        if os.path.getsize(full) != sizes[name]:
            return ("its file %s is %d bytes, not the %d recorded when it was "
                    "written (incomplete or modified)"
                    % (name, os.path.getsize(full), sizes[name]))
    return None


def check_complete(root: str, prefix: str) -> Optional[str]:
    """``None`` when ``prefix`` is a complete calculation whose files are intact."""
    state, record = existing_state(root, prefix)
    if state == "foreign":
        return describe_foreign(root, prefix, record)
    if state != "complete":
        return {"absent": "it has no completion record",
                "incomplete": "its completion record says an earlier run was "
                              "interrupted before it was complete"}[state]
    for entry in record.get("results") or ():
        problem = _check_files(root, (entry.get("files") or {}).values(),
                               entry.get("file_bytes") or {})
        if problem:
            return problem
    return None


def _write_json_atomically(root: str, name: str, payload: Mapping[str, Any],
                           token: str) -> None:
    from ..retrieval import publish_file_mode

    tmp = os.path.join(root, "%s%s-%s" % (_TEMPORARY, token, name))
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
            fh.write("\n")
        publish_file_mode(tmp)
        os.replace(tmp, os.path.join(root, name))
    except BaseException:
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------


class _PrefixLock(object):
    """``.biotsr-lock-<prefix>``: one installer per calculation at a time.

    Taken with an exclusive create.  A lock whose owner no longer exists on
    this machine is taken over; one held by a live process, or by another
    machine, is reported.
    """

    def __init__(self, root: str, prefix: str, token: str):
        self.path = os.path.join(root, _LOCK + prefix)
        self.token = token
        self.held = False

    def __enter__(self):
        for _attempt in range(3):
            try:
                fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except FileExistsError:
                try:
                    with open(self.path) as fh:
                        owner = fh.read().strip()
                except OSError:
                    owner = ""
                if owner and _token_state(owner) == "dead":
                    try:
                        os.remove(self.path)
                    except FileNotFoundError:
                        pass
                    continue
                raise CLIError(
                    "another run is installing %s in %s (lock %s)"
                    % (os.path.basename(self.path)[len(_LOCK):],
                       os.path.dirname(self.path), os.path.basename(self.path)),
                    "give each run its own --output, or distinct sample ids; delete "
                    "the lock only if no BioTSR run is writing into that folder",
                )
            with os.fdopen(fd, "w") as fh:
                fh.write(self.token)
            self.held = True
            return self
        raise CLIError("could not take the lock %s" % self.path)

    def __exit__(self, *exc):
        if self.held:
            try:
                os.remove(self.path)
            except FileNotFoundError:
                pass
        return False


class CalculationWriter(object):
    """Stages one calculation privately and installs it into the shared folder.

    ``check()`` applies the collision policy before anything is computed;
    ``with writer:`` creates the staging directory (``writer.path``) and, on a
    clean exit, installs whatever was written there with ``writer.payload`` as
    its completion record.  On an exception the staging directory is removed
    and nothing in the folder changes.
    """

    def __init__(self, root: str, prefix: str, if_exists: str,
                 token: Optional[str] = None,
                 verify: Optional[Callable[[str], Optional[str]]] = None,
                 shared_staging: bool = False):
        self.root = root
        self.prefix = validate_prefix(prefix)
        self.token = token or invocation_token()
        # The workers of a batch share one token and so one staging directory,
        # which only the batch removes, when it ends: a worker removing it
        # while it happens to be empty races another worker creating its own
        # staging inside it (macOS answers EINVAL).
        self._owns_staging_root = not shared_staging
        self.if_exists = if_exists
        self.verify = verify
        self.skipped = False
        self.payload: Optional[Dict[str, Any]] = None
        self.staging_root = os.path.join(root, _STAGING + self.token)
        self.staging = os.path.join(self.staging_root, self.prefix)
        self._created_root = False

    def check(self) -> bool:
        """Apply the collision policy **before** anything is computed or created.

        Returns ``True`` when the caller should skip this calculation.  A run
        that is refused, or fails for any other reason, leaves no output folder
        behind: ``biotsr run --pdb /nope.pdb -o results/`` creates nothing.
        """
        if os.path.lexists(self.root) and not os.path.isdir(self.root):
            raise CLIError("--output %s exists and is not a directory" % self.root)
        self._refuse_case_collisions()
        state, record = existing_state(self.root, self.prefix)
        if state == "foreign":
            raise CLIError(
                "output %s in %s cannot be written: %s"
                % (self.prefix, self.root, describe_foreign(self.root, self.prefix, record)),
                "BioTSR never replaces or deletes a file it did not write as "
                "this calculation; remove those files or choose another sample id",
            )
        if state == "complete":
            if self.if_exists == IF_EXISTS_FAIL:
                raise CLIError(
                    "output %s already exists in %s" % (self.prefix, self.root),
                    "pass --if-exists overwrite to replace it, or --if-exists "
                    "skip to keep it and continue",
                )
            if self.if_exists == IF_EXISTS_SKIP:
                problem = self.verify(self.prefix) if self.verify else None
                if problem:
                    raise CLIError(
                        "output %s exists but cannot be reused: %s"
                        % (self.prefix, problem),
                        "--if-exists skip only keeps a calculation produced by "
                        "exactly this request; pass --if-exists overwrite to "
                        "recompute it, choose another --output, or give the row "
                        "its own sample_id",
                    )
                self.skipped = True
        elif state == "incomplete" and self.if_exists == IF_EXISTS_FAIL:
            raise CLIError(
                "output %s in %s was left incomplete by an interrupted run"
                % (self.prefix, self.root),
                "pass --if-exists skip (resume) or --if-exists overwrite to "
                "compute it again",
            )
        return self.skipped

    def _refuse_case_collisions(self) -> None:
        """A file whose name differs only in letter case is a collision too."""
        try:
            names = os.listdir(self.root)
        except OSError:
            return
        wanted = {(self.prefix + s).casefold(): self.prefix + s for s in (
            native.KEYS_SUFFIX, native.METADATA_SUFFIX, native.TRIPLETS_SUFFIX,
            RECORD_SUFFIX)}
        for name in names:
            exact = wanted.get(name.casefold())
            if exact is not None and exact != name:
                raise CLIError(
                    "output %s in %s collides with the existing file %s, whose "
                    "name differs only in letter case" % (self.prefix, self.root, name),
                    "on macOS and Windows these are one file; give the row its "
                    "own sample id",
                )

    def __enter__(self):
        if self.check():
            return self
        self._created_root = not os.path.exists(self.root)
        os.makedirs(self.root, exist_ok=True)
        os.makedirs(self.staging_root, exist_ok=True)
        os.makedirs(self.staging)
        return self

    def __exit__(self, exc_type, exc, tb):
        if self.skipped:
            return False
        try:
            if exc_type is None:
                if self.payload is None:
                    raise RuntimeError("nothing was written for %s" % self.prefix)
                self._install()
        except BaseException:
            self._discard(remove_root=True)
            raise
        self._discard(remove_root=exc_type is not None)
        return False

    def _discard(self, remove_root: bool) -> None:
        shutil.rmtree(self.staging, ignore_errors=True)
        if self._owns_staging_root:
            try:
                os.rmdir(self.staging_root)
            except OSError:
                pass
        if remove_root and self._created_root:
            # A failed run leaves no output folder behind, not even an empty
            # one -- unless something else has been written into it.
            try:
                os.rmdir(self.root)
            except OSError:
                pass

    def _install(self) -> None:
        new_names = sorted(os.listdir(self.staging))
        record = record_name(self.prefix)
        with _PrefixLock(self.root, self.prefix, self.token):
            state, old = existing_state(self.root, self.prefix)
            if state == "foreign":
                raise CLIError(
                    "output %s in %s cannot be written: %s"
                    % (self.prefix, self.root,
                       describe_foreign(self.root, self.prefix, old)),
                    "remove those files or choose another sample id",
                )
            if state == "complete" and self.if_exists != IF_EXISTS_OVERWRITE:
                raise CLIError(
                    "output %s appeared in %s while it was being computed"
                    % (self.prefix, self.root),
                    "another run is writing the same calculation into the same "
                    "--output; give each run its own output folder or distinct "
                    "sample ids",
                )
            if state == "incomplete" and self.if_exists == IF_EXISTS_FAIL:
                raise CLIError(
                    "output %s in %s was left incomplete by another run while "
                    "this one was computing" % (self.prefix, self.root),
                    "pass --if-exists skip or --if-exists overwrite to replace it",
                )
            owned = set(_record_files(old)) if old else set()
            listing = {name.casefold(): name for name in os.listdir(self.root)}
            for name in new_names + [record]:
                present = listing.get(name.casefold())
                if present is None:
                    continue
                full = os.path.join(self.root, present)
                if os.path.islink(full):
                    raise CLIError(
                        "%s in %s is a symbolic link" % (present, self.root),
                        "BioTSR never writes through, replaces or deletes a link",
                    )
                if present != name or name not in owned:
                    raise CLIError(
                        "%s already exists in %s and is not part of an earlier "
                        "version of this calculation" % (present, self.root),
                        "BioTSR never replaces or deletes a file it did not write "
                        "as this calculation; remove it or choose another sample id",
                    )
            intent = dict(self.payload)
            intent.update(complete=False, files_pending=new_names, writer=self.token)
            _write_json_atomically(self.root, record, intent, self.token)
            for name in new_names:
                os.replace(os.path.join(self.staging, name),
                           os.path.join(self.root, name))
            for name in sorted(owned - set(new_names) - {record}):
                full = os.path.join(self.root, name)
                if os.path.isfile(full) and not os.path.islink(full):
                    os.remove(full)
            final = dict(self.payload)
            final.update(complete=True)
            _write_json_atomically(self.root, record, final, self.token)

    @property
    def path(self) -> str:
        return self.staging


def write_outputs(
    outcome,
    out_dir: str,
    *,
    prefix: str,
    per_residue: bool = False,
    triplets: bool = True,
    command: str = "run",
    argv=None,
    request: Optional[Mapping[str, Any]] = None,
    metadata: Optional[Mapping[str, str]] = None,
    source: Optional[Mapping[str, Any]] = None,
    streamed_triplets: Optional[Mapping[str, str]] = None,
) -> Dict[str, Any]:
    """Serialize one outcome into a staging directory and return its completion record.

    One result is written under ``prefix``; the per-residue sweep writes each
    residue under ``<prefix>_<NAME>_<id>``.  The record lists every file with
    its size; :class:`CalculationWriter` writes it into the folder last.
    ``streamed_triplets`` maps an entity to the triplet file already written
    while it was computed.
    """
    entries = []
    for name, result in outcome.results.items():
        if per_residue:
            seq, _, residue = str(name).partition("_")
            member = "%s_%s_%s" % (prefix, residue, seq)
        else:
            member = prefix
        streamed = (streamed_triplets or {}).get(name)
        written = native.save(result, out_dir, prefix=member,
                              triplets=triplets and streamed is None)
        files = {kind: os.path.basename(path) for kind, path in written.items()}
        if streamed:
            files["triplets"] = streamed
        entries.append(
            {
                "name": member,
                "entity": name,
                "files": files,
                # Each file's size once fully written, so a later reader can
                # tell a complete file from one cut short or appended to since.
                "file_bytes": {
                    file_name: os.path.getsize(os.path.join(out_dir, file_name))
                    for file_name in files.values()
                },
                "structure_id": result.structure_id.value,
                "chain": result.chain,
                "target": result.annotations.get("target"),
                "key_space_fingerprint": result.key_space.fingerprint,
                "sampling_fingerprint": result.sampling.fingerprint,
                "counts": dict(result.counts.to_dict()),
                "category_totals": dict(result.category_totals()),
                "n_warnings": len(result.warnings),
                "warnings": list(result.warnings),
                # What the selection asked for and did not get, from its
                # structured record; ``{}`` when it got everything (ADR-106).
                "omissions": selection_omissions(
                    selection_sides(result.annotations.get("selection_request")),
                    result.warnings),
            }
        )

    payload: Dict[str, Any] = {
        "schema": REQUEST_SCHEMA,
        "prefix": prefix,
        "biotsr_version": _version(),
        "command": command,
        "argv": list(argv) if argv is not None else None,
        "method": outcome.spec.method,
        "per_residue": bool(per_residue),
        "profile": {
            # ``resolved`` is the method's standard profile; ``active`` is what
            # ``Profile.with_modifiers`` produced from it, which carries the
            # modifiers in its name.
            "resolved": outcome.base_profile or outcome.profile,
            "active": outcome.profile,
            "sse_strand_identity": _strand_identity(outcome),
        },
        "parameters": dict(outcome.spec.user_parameters()),
        "label": outcome.spec.label,
        "triplets": triplets,
        "results": entries,
    }
    if source:
        payload["source"] = dict(source)
    if request is not None:
        payload["request"] = dict(request)
        payload["request_sha256"] = request_sha256(request)
    if metadata:
        payload["metadata"] = dict(metadata)
    if getattr(outcome, "warnings", ()):
        payload["warnings"] = list(outcome.warnings)
    return payload


def _strand_identity(outcome):
    """Which convention named a beta-strand, or ``None`` when SSE was not used.

    Recorded because it is the one annotation convention that is *user-visible
    and profile-dependent* while appearing in neither identity fingerprint: two
    runs with the same KeySpace and the same SamplingProfile can still put their
    frequencies in different SSE categories (ADR-063).
    """
    if not outcome.spec.modifiers.get("sse"):
        return None
    from ..methods.registry import get_method

    return get_method(outcome.spec.method).default_profile.sse_strand_identity


def _version():
    from .. import __version__

    return __version__
