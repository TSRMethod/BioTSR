"""Every command and code block in the README and the tutorials, executed.

The README and ``docs/tutorials/`` are the first thing a new user runs, and the
examples on them decay silently unless something runs them.  This module parses
each page into its fenced blocks and executes them, in order, from a fresh
working directory holding the pinned fixture structures the page uses:

* ``bash`` -- each command is run; ``biotsr ...`` is run as
  ``python -m biotsr ...`` against this source tree, so the page is tested
  whether or not the console script is on ``PATH``.  A command must exit 0.
* ``text`` directly after a ``bash`` block -- every line must appear, with
  whitespace normalised, in what the preceding commands printed.
* ``python`` -- executed in one namespace per page; the pages state their
  claims as ``assert`` statements.
* a block announced by ``<!-- file: NAME -->`` is written to ``NAME`` first.
* a block announced by ``<!-- not-executed -->`` is parsed with ``shlex`` and
  not run: installation commands, which would reach the network, and commands
  whose documented exit status is not 0 (a partial result exits 4, a refused
  analysis 2) -- each of those is followed on its page by a ``python`` block
  that runs the same command through ``biotsr.cli.main`` and asserts the exit
  status and output.
* a block announced by ``<!-- requires: plot -->`` draws a figure: it runs when
  matplotlib (the ``plot`` extra) is installed and is skipped when it is not.

The pages use only structures pinned under ``tests/data/structures/``; in a
checkout without them the module skips, naming the missing path.
"""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import sys

import pytest

import biotsr

from ..conftest import FIXTURE_PDB, SRC, SYNTHETIC, WORKSPACE, _require

TUTORIALS = os.path.join(WORKSPACE, "docs", "tutorials")
README = os.path.join(WORKSPACE, "README.md")
USER_GUIDE = os.path.join(WORKSPACE, "docs", "USER_GUIDE.md")

#: The pinned structures a page may use, copied into its working directory and
#: into ``pdb/`` beneath it (the directory the manifest examples resolve against).
STRUCTURES = ("1UBQ", "1BRS", "3PTB", "1BNA", "7LF9", "6G9M", "2AC0")
#: Synthetic structures from ``tests/data/synthetic`` (its README).
SYNTHETIC_STRUCTURES = ("SYNPHOS",)

#: The user guide also runs one example of every atomic and cross family.
USER_GUIDE_STRUCTURES = STRUCTURES + ("1BNA", "2AC0", "5WT9", "7LF9")

_FENCE = re.compile(r"^```([A-Za-z0-9_-]*)\s*$")
_DIRECTIVE = re.compile(
    r"^<!--\s*(file:\s*(?P<file>\S+)|(?P<skip>not-executed)|(?P<plot>requires: plot))\s*-->\s*$")


def plotting_available():
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        return False
    return True


def parse_blocks(path):
    """``[(lang, text, directive)]`` for every fenced block on the page."""
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    blocks, directive, i = [], None, 0
    while i < len(lines):
        m = _DIRECTIVE.match(lines[i].strip())
        if m:
            directive = (("file", m.group("file")) if m.group("file")
                         else ("plot", None) if m.group("plot") else ("skip", None))
            i += 1
            continue
        fence = _FENCE.match(lines[i])
        if fence:
            lang = fence.group(1)
            body = []
            i += 1
            while i < len(lines) and not lines[i].startswith("```"):
                body.append(lines[i])
                i += 1
            blocks.append((lang, "\n".join(body) + "\n", directive))
            directive = None
        elif lines[i].strip():
            directive = None
        i += 1
    return blocks


def shell_commands(text):
    """Commands in a ``bash`` block: continuations joined, comments dropped."""
    joined = text.replace("\\\n", " ")
    commands = []
    for line in joined.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        commands.append(shlex.split(line, comments=True))
    return commands


def _normalise(text):
    return [" ".join(line.split()) for line in text.splitlines() if line.strip()]


def _environment():
    env = dict(os.environ)
    env["PYTHONPATH"] = SRC + os.pathsep + env.get("PYTHONPATH", "")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return env


def _run(argv, cwd):
    if argv[0] == "biotsr":
        argv = [sys.executable, "-m", "biotsr"] + argv[1:]
    elif argv[0] in ("python", "python3"):
        argv = [sys.executable] + argv[1:]
    else:
        raise AssertionError("a documented command must be biotsr or python: %r" % argv)
    proc = subprocess.run(argv, cwd=cwd, env=_environment(), stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, universal_newlines=True)
    assert proc.returncode == 0, "%s exited %d:\n%s%s" % (
        " ".join(argv), proc.returncode, proc.stdout, proc.stderr)
    return proc.stdout + proc.stderr


def _workdir(tmp_path, structures=STRUCTURES):
    sources = [os.path.join(FIXTURE_PDB, s + ".pdb") for s in structures]
    sources += [os.path.join(SYNTHETIC, s + ".pdb") for s in SYNTHETIC_STRUCTURES]
    for src in sources:
        _require(src)
    os.makedirs(os.path.join(str(tmp_path), "pdb"))
    for src in sources:
        shutil.copy(src, str(tmp_path))
        shutil.copy(src, os.path.join(str(tmp_path), "pdb"))
    return str(tmp_path)


def execute_page(path, cwd):
    """Run every block of one page; return the number of blocks executed."""
    namespace = {"__name__": "__tutorial__"}
    last_output = None
    executed = 0
    previous_cwd = os.getcwd()
    os.chdir(cwd)
    try:
        for lang, text, directive in parse_blocks(path):
            if directive and directive[0] == "skip":
                if lang in ("bash", "sh", "shell"):
                    shell_commands(text)       # must at least parse
                continue
            if directive and directive[0] == "plot" and not plotting_available():
                continue       # the figure needs the plot extra; the numbers do not
            if directive and directive[0] == "file":
                with open(os.path.join(cwd, directive[1]), "w", encoding="utf-8") as fh:
                    fh.write(text)
                continue
            if lang in ("bash", "sh", "shell"):
                last_output = "".join(_run(c, cwd) for c in shell_commands(text))
                executed += 1
            elif lang == "text":
                if last_output is None:
                    continue       # a layout illustration, not command output
                produced = "\n".join(_normalise(last_output))
                for expected in _normalise(text):
                    assert expected in produced, (
                        "%s: expected output line %r not printed; got:\n%s"
                        % (os.path.basename(path), expected, last_output))
                last_output = None
                executed += 1
            elif lang == "python":
                exec(compile(text, path, "exec"), namespace)
                executed += 1
            else:
                last_output = None
    finally:
        os.chdir(previous_cwd)
    return executed


def _pages():
    if not os.path.isdir(TUTORIALS):
        return []
    return sorted(f for f in os.listdir(TUTORIALS) if re.match(r"\d\d_.*\.md$", f))


def test_the_tutorials_exist_and_cover_the_release_workflows():
    _require(TUTORIALS)
    pages = _pages()
    assert len(pages) >= 7, pages
    joined = " ".join(pages)
    for topic in ("protein_calpha", "atomic", "sse", "cross", "batch", "load", "analysis"):
        assert topic in joined, topic


@pytest.mark.parametrize("page", _pages() or ["<none>"])
def test_every_tutorial_block_runs(page, tmp_path):
    if page == "<none>":
        _require(TUTORIALS)
    executed = execute_page(os.path.join(TUTORIALS, page), _workdir(tmp_path))
    assert executed >= 2, "%s executed only %d block(s)" % (page, executed)


def test_every_user_guide_block_runs(tmp_path):
    _require(USER_GUIDE)
    executed = execute_page(USER_GUIDE, _workdir(tmp_path, USER_GUIDE_STRUCTURES))
    assert executed >= 40, executed


def test_the_user_guide_indexes_every_released_method():
    """The method table at the top of the guide names all twenty methods."""
    with open(_require(USER_GUIDE), encoding="utf-8") as fh:
        text = fh.read()
    table = text.split("## The methods at a glance", 1)[1].split("\n## ", 1)[0]
    indexed = set(re.findall(r"^\| `([a-z_]+)` \|", table, re.M))
    assert indexed == set(biotsr.list_methods())


def test_every_readme_block_runs(tmp_path):
    executed = execute_page(README, _workdir(tmp_path))
    assert executed >= 3
