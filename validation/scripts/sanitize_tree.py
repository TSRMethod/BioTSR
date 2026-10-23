#!/usr/bin/env python3
"""Remove workstation paths and personal account identifiers from HEAD.

Generated evidence can pick up identifying material that is not science: the
absolute path of the machine a gate ran on, a temporary directory, an account
identifier in a home or cluster directory, an e-mail address copied from historical
packaging metadata, or Office metadata inside a tracked spreadsheet.

This script sanitizes the **current tree only**.  It deliberately does not touch
git history (a separate, human decision), and it never removes a person's
*name* where the name is authorship or provenance: only paths, account
identifiers and e-mail addresses are replaced, each with a placeholder that says
what was there.  Every rule is a pattern, so this file itself contains none of
the values it removes.

What is left alone, on purpose:

* every scientific value, hash, count and fingerprint;
* the pinned reference data under ``tests/data/`` beyond the patterns below;
* ``source.sha256`` -- the identity of an input structure, which is not a path;
* the cell contents of any tracked spreadsheet (only package metadata is
  removed, so every value it records is unchanged).

Idempotent.  ``--check`` exits 1 if anything would still change.

Usage::

    python3 validation/scripts/sanitize_tree.py [--check]
"""

from __future__ import annotations

import argparse
import io
import os
import re
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
WORKSPACE = os.path.dirname(os.path.dirname(HERE))

#: Only these trees are rewritten; code and current documentation are fixed at
#: the source (the writers now record relative paths) rather than by pattern.
SCOPES = ("validation/evidence/", "tests/data/")
TEXT_SUFFIXES = (".json", ".md", ".txt", ".csv")

#: (pattern, replacement, what it removes)
RULES = (
    # this checkout, as an absolute path on whichever machine produced the file
    (re.compile(re.escape(WORKSPACE) + r"(?=[/\s\"']|$)"), "<workspace>",
     "absolute path of the checkout"),
    # temporary directories (macOS per-user, and any other under /private/tmp)
    (re.compile(r"(/private)?/var/folders/[A-Za-z0-9_]{2}/[A-Za-z0-9_]+/T"), "<tmp>",
     "per-user temporary directory"),
    (re.compile(r"/private/tmp/[^\s\"'`)]+"), "<tmp>", "temporary directory"),
    # home and cluster directories named by an account
    (re.compile(r"/Users/(?!<)[^/\s\"']+"), "/Users/<user>", "home directory named by a user"),
    (re.compile(r"/home/(?!<)[A-Za-z][A-Za-z0-9_.-]*(/[A-Za-z]?\d{6,})?"), "/home/<user>",
     "home directory named by a user"),
    (re.compile(r"/(work|project)/[a-z]+\d+"), r"/\1/<user>",
     "cluster work/project directory named by an account"),
    # e-mail addresses (names are kept; the address is not provenance)
    (re.compile(r"[A-Za-z0-9._%+-]+@(?!users\.noreply\.github\.com)[A-Za-z0-9-]+(\.[A-Za-z0-9-]+)*"
                r"\.[A-Za-z]{2,}"), "<email redacted>", "e-mail address"),
)


def tracked_files():
    out = subprocess.run(["git", "-C", WORKSPACE, "ls-files"], stdout=subprocess.PIPE,
                         universal_newlines=True, check=True).stdout.split("\n")
    return [p for p in out if p.startswith(SCOPES)]


def sanitize_text(text):
    counts = {}
    for pattern, replacement, what in RULES:
        text, n = pattern.subn(replacement, text)
        if n:
            counts[what] = counts.get(what, 0) + n
    return text, counts


# --- spreadsheets (.xlsx) ----------------------------------------------------

_ABS_PATH = re.compile(
    rb"<mc:AlternateContent[^>]*>\s*<mc:Choice[^>]*>\s*<x15ac:absPath[^>]*/>\s*"
    rb"</mc:Choice>\s*</mc:AlternateContent>"
)
_CUSTOM_REL = re.compile(rb'<Relationship [^>]*Target="docProps/custom\.xml"[^>]*/>')
_CUSTOM_TYPE = re.compile(rb'<Override [^>]*PartName="/docProps/custom\.xml"[^>]*/>')


def sanitize_xlsx(path):
    """Drop the cloud-drive path and the sensitivity-label part; keep every cell."""
    with zipfile.ZipFile(path) as zin:
        infos = zin.infolist()
        parts = {i.filename: zin.read(i.filename) for i in infos}
    changes = {}
    if "docProps/custom.xml" in parts:
        changes["document sensitivity-label properties (docProps/custom.xml)"] = 1
    new = {}
    for name, data in parts.items():
        if name == "docProps/custom.xml":
            continue
        if name == "xl/workbook.xml":
            data, n = _ABS_PATH.subn(b"", data)
            if n:
                changes["cloud-drive path with account identifier (x15ac:absPath)"] = n
        elif name == "_rels/.rels":
            data = _CUSTOM_REL.sub(b"", data)
        elif name == "[Content_Types].xml":
            data = _CUSTOM_TYPE.sub(b"", data)
        new[name] = data
    if not changes:
        return None, {}
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zout:
        for info in infos:
            if info.filename not in new:
                continue
            out = zipfile.ZipInfo(info.filename, date_time=info.date_time)
            out.compress_type = info.compress_type
            out.external_attr = info.external_attr
            zout.writestr(out, new[info.filename])
    return buf.getvalue(), changes


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    pending = 0
    for rel in tracked_files():
        path = os.path.join(WORKSPACE, rel)
        if rel.endswith(".xlsx"):
            data, counts = sanitize_xlsx(path)
        elif rel.endswith(TEXT_SUFFIXES):
            with open(path, encoding="utf-8", newline="") as fh:
                original = fh.read()
            text, counts = sanitize_text(original)
            data = text.encode("utf-8") if counts else None
        else:
            continue
        if not counts:
            continue
        pending += 1
        print("%s %s: %s" % ("would sanitize" if args.check else "sanitized", rel,
                             "; ".join("%d x %s" % (n, w) for w, n in sorted(counts.items()))))
        if not args.check:
            with open(path, "wb") as fh:
                fh.write(data)
    print("%d file(s) %s" % (pending, "need sanitizing" if args.check else "sanitized"))
    return 1 if (args.check and pending) else 0


if __name__ == "__main__":
    sys.exit(main())
