"""Every relative Markdown link in the current documentation resolves.

The user-facing pages link to one another and to files in the repository, and
a moved or renamed page breaks those links silently: GitHub renders a dead
link without complaint.  This module parses each current page, skips fenced
code blocks and external URLs, and requires every relative target to exist.

Every page under ``docs/`` is current documentation, so every page is checked.
"""

from __future__ import annotations

import glob
import os
import re

import pytest

from ..conftest import REPOSITORY_ONLY, WORKSPACE, evidence_available

#: The documentation stays in the repository; the source distribution carries
#: only README.md and CHANGELOG.md, whose links point into it.
pytestmark = pytest.mark.skipif(
    not evidence_available(os.path.join(WORKSPACE, "docs", "README.md")),
    reason="%s: docs/" % REPOSITORY_ONLY,
)

#: Current pages: the root documents and every page under ``docs/``.
CURRENT_PAGES = sorted(
    {os.path.relpath(p, WORKSPACE) for p in (
        [os.path.join(WORKSPACE, name) for name in (
            "README.md", "CHANGELOG.md", "RELEASE_CHECKLIST.md")]
        + glob.glob(os.path.join(WORKSPACE, "docs", "*.md"))
        + glob.glob(os.path.join(WORKSPACE, "docs", "tutorials", "*.md"))
        + glob.glob(os.path.join(WORKSPACE, "docs", "architecture", "*.md"))
        + glob.glob(os.path.join(WORKSPACE, "docs", "methods", "*.md"))
    ) if os.path.exists(p)}
)

_LINK = re.compile(r"(?<!!)\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
_FENCE = re.compile(r"^\s*(```|~~~)")


def relative_links(path):
    """``[(line_number, target)]`` for every relative link outside code blocks."""
    links, in_code = [], False
    with open(path, encoding="utf-8") as fh:
        for number, line in enumerate(fh, 1):
            if _FENCE.match(line):
                in_code = not in_code
                continue
            if in_code:
                continue
            # Inline code spans are not links.
            text = re.sub(r"`[^`]*`", "", line)
            for target in _LINK.findall(text):
                if re.match(r"^[a-z][a-z0-9+.-]*:", target, re.I) or target.startswith("#"):
                    continue
                links.append((number, target))
    return links


def broken_links(rel):
    path = os.path.join(WORKSPACE, rel)
    base = os.path.dirname(path)
    broken = []
    for number, target in relative_links(path):
        file_part = target.split("#", 1)[0]
        if not os.path.exists(os.path.normpath(os.path.join(base, file_part))):
            broken.append("%s:%d -> %s" % (rel, number, target))
    return broken


def test_the_current_pages_are_found():
    assert "README.md" in CURRENT_PAGES
    assert os.path.join("docs", "README.md") in CURRENT_PAGES
    assert len(CURRENT_PAGES) >= 20, CURRENT_PAGES


@pytest.mark.parametrize("rel", CURRENT_PAGES)
def test_every_relative_link_resolves(rel):
    assert broken_links(rel) == []


#: README.md is the PyPI long description and CHANGELOG.md ships in the sdist;
#: neither renders relative links outside GitHub, so they link absolutely.
SHIPPED_PAGES = ("README.md", "CHANGELOG.md")

_REPOSITORY_URL = re.compile(
    r"https://github\.com/TSRMethod/BioTSR/(?:blob|tree)/main/([^)\s#>`\"]+)")


@pytest.mark.parametrize("rel", SHIPPED_PAGES)
def test_the_shipped_pages_use_no_relative_links(rel):
    assert relative_links(os.path.join(WORKSPACE, rel)) == []


@pytest.mark.parametrize("rel", CURRENT_PAGES)
def test_every_link_into_the_repository_resolves(rel):
    with open(os.path.join(WORKSPACE, rel), encoding="utf-8") as fh:
        targets = _REPOSITORY_URL.findall(fh.read())
    missing = [t for t in targets if not os.path.exists(os.path.join(WORKSPACE, t))]
    assert missing == []


def test_every_repository_link_in_the_package_resolves():
    """Help text, messages and the project URLs point users at the docs."""
    files = [os.path.join(WORKSPACE, "pyproject.toml")] + glob.glob(
        os.path.join(WORKSPACE, "src", "biotsr", "**", "*.py"), recursive=True)
    missing = []
    for path in files:
        with open(path, encoding="utf-8") as fh:
            for target in _REPOSITORY_URL.findall(fh.read()):
                target = target.rstrip(".,;:")
                if not os.path.exists(os.path.join(WORKSPACE, target)):
                    missing.append("%s -> %s" % (os.path.relpath(path, WORKSPACE), target))
    assert missing == []
    with open(os.path.join(WORKSPACE, "pyproject.toml"), encoding="utf-8") as fh:
        assert "/blob/main/docs/USER_GUIDE.md" in fh.read()
