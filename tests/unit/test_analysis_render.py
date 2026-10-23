"""The optional figure: its dependency, its absence, and what it draws.

rc2 declared ``plot = ["matplotlib", "seaborn"]`` while the renderer also
imported pandas, and ``seaborn.clustermap`` needs scipy at run time, which
neither seaborn nor the extra installs -- so a clean ``biotsr[plot]`` could not
draw.  It also switched the caller's matplotlib backend and closed every open
figure.  The renderer now draws with matplotlib alone on a standalone
``Figure``.  The tests that draw need the ``plot`` extra and skip without it;
the others always run.
"""

from __future__ import annotations

import ast
import os
import re
import subprocess
import sys
import textwrap

import pytest

from biotsr.analysis import clustermap_spec, distance_matrix, feature_matrix

from ..conftest import OPTIONAL_DEPENDENCY, PACKAGE_ROOT, SRC
from .test_analysis_features import sample

RENDER = os.path.join(SRC, "biotsr", "analysis", "render.py")
PYPROJECT = os.path.join(PACKAGE_ROOT, "pyproject.toml")


def _render_imports():
    with open(RENDER) as fh:
        tree = ast.parse(fh.read())
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and not node.level:
            names.add(node.module.split(".")[0])
    stdlib = {"__future__", "typing"}
    return sorted(names - stdlib)


def _plot_extra():
    with open(PYPROJECT) as fh:
        text = fh.read()
    m = re.search(r"^plot\s*=\s*\[(.*?)\]", text, re.M | re.S)
    assert m, "no plot extra in pyproject.toml"
    return sorted(re.split(r"[<>=!~ ]", s.strip().strip('"'))[0]
                  for s in m.group(1).split(",") if s.strip())


def test_the_plot_extra_declares_exactly_what_the_renderer_imports():
    assert _render_imports() == ["matplotlib"]
    assert _plot_extra() == _render_imports()


def test_the_renderer_never_touches_pyplot_or_the_backend():
    """No pyplot import, no ``use(...)`` and no ``close(...)`` call anywhere."""
    with open(RENDER) as fh:
        tree = ast.parse(fh.read())
    imported = set()
    calls = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.add(node.module or "")
            imported.update("%s.%s" % (node.module, a.name) for a in node.names)
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            calls.add(node.func.attr)
    assert not any("pyplot" in name for name in imported)
    assert not {"use", "close", "switch_backend"} & calls


def _blocks():
    return [
        sample("q1", {10: 5, 11: 5}, label="Q"), sample("p1", {1: 5, 2: 5, 3: 1}, label="P"),
        sample("p2", {1: 5, 2: 4, 3: 1}, label="P"), sample("q2", {10: 5, 11: 4}, label="Q"),
        sample("p3", {1: 4, 2: 5, 3: 2}, label="P"),
    ]


_BLOCK = textwrap.dedent('''
    import importlib.abc, sys

    class _NoPlotting(importlib.abc.MetaPathFinder):
        def find_spec(self, name, path, target=None):
            if name.split(".")[0] == "matplotlib":
                raise ImportError("blocked for the test", name=name.split(".")[0])
            return None

    sys.meta_path.insert(0, _NoPlotting())
''')


def test_every_number_is_computed_without_the_plotting_stack_and_a_figure_says_how(tmp_path):
    (tmp_path / "sitecustomize.py").write_text(_BLOCK)
    script = textwrap.dedent('''
        import sys
        sys.path.insert(0, %r)
        from biotsr.analysis import AnalysisSample, feature_matrix, distance_matrix, clustermap_spec
        from biotsr.analysis.render import MissingRenderDependency, render_clustermap
        from tests.unit.test_analysis_render import _blocks
        s = _blocks()
        dm = distance_matrix(feature_matrix(s))
        spec = clustermap_spec(dm, [x.label for x in s], n_clusters=2)
        assert spec.adjusted_rand == 1.0
        print(dm.to_csv().splitlines()[0])
        print(spec.to_assignments_csv().splitlines()[0])
        try:
            render_clustermap(spec, "x.png")
        except MissingRenderDependency as exc:
            print("REFUSED", exc)
        assert "matplotlib" not in sys.modules
    ''' % PACKAGE_ROOT)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(tmp_path), SRC, PACKAGE_ROOT]))
    done = subprocess.run([sys.executable, "-c", script], env=env, cwd=str(tmp_path),
                          capture_output=True, text=True, timeout=120)
    assert done.returncode == 0, done.stderr
    assert "generalized_jaccard_distance,q1,p1,p2,q2,p3" in done.stdout
    assert "sample_id,cluster,leaf_position,class_label" in done.stdout
    assert 'REFUSED' in done.stdout and 'pip install "biotsr[plot]"' in done.stdout
    assert not (tmp_path / "x.png").exists()


# ---------------------------------------------------------------------------
# Drawing (needs the plot extra)
# ---------------------------------------------------------------------------


def _plotting():
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        pytest.skip("%s: matplotlib (the plot extra)" % OPTIONAL_DEPENDENCY)


def _axes(fig):
    """(heat map, stripes, legend) axes, identified by what they show."""
    heat = next(a for a in fig.axes if a.images and a.get_yticklabels()
                and a.yaxis.get_ticks_position() == "right")
    stripe = next(a for a in fig.axes if a.images and a is not heat
                  and [t.get_text() for t in a.get_xticklabels()][-1:] == ["cluster"])
    (legend,) = fig.legends
    return heat, stripe, legend


def test_a_drawn_clustermap_keeps_rows_stripes_and_tree_aligned(tmp_path):
    _plotting()
    import matplotlib

    from biotsr.analysis.render import _palette, render_clustermap

    backend = matplotlib.get_backend()
    samples = _blocks()
    dm = distance_matrix(feature_matrix(samples))
    spec = clustermap_spec(dm, [s.label for s in samples], n_clusters=2)
    fig = render_clustermap(spec, str(tmp_path / "map.png"))
    assert (tmp_path / "map.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert matplotlib.get_backend() == backend
    heat, stripe, legend = _axes(fig)
    order = list(spec.ordered_sample_ids())
    # tick labels on both axes are the dendrogram order
    assert [t.get_text() for t in heat.get_yticklabels()] == order
    assert [t.get_text() for t in heat.get_xticklabels()] == order
    # every drawn cell is the matrix value for that (row id, column id)
    drawn = heat.images[0].get_array()
    for r, a in enumerate(order):
        for c, b in enumerate(order):
            i, j = dm.sample_ids.index(a), dm.sample_ids.index(b)
            assert drawn[r][c] == dm.value(i, j)
    assert heat.images[0].get_clim() == (0.0, 1.0)
    # each row's stripe colours are its own class's and cluster's
    classes = _palette(matplotlib, spec.class_palette, len(spec.unique_classes))
    clusters = _palette(matplotlib, spec.cluster_palette, len(spec.unique_clusters))
    rgb = stripe.images[0].get_array()
    for r, sid in enumerate(order):
        i = dm.sample_ids.index(sid)
        assert tuple(rgb[r][0]) == pytest.approx(classes[spec.row_class_slots()[i]])
        assert tuple(rgb[r][1]) == pytest.approx(clusters[spec.row_cluster_slots()[i]])
    # Observed classes and algorithm clusters are separate, distinctly named stripes.
    assert [t.get_text() for t in stripe.get_xticklabels()] == ["observed class", "cluster"]
    # the legend names every class and cluster, under its heading
    assert [t.get_text() for t in legend.get_texts()] == [
        "observed class", "P", "Q", "cluster (maxclust 2)", "1", "2"]
    labels = [a.get_ylabel() for a in fig.axes]
    assert "generalized Jaccard distance (0-1)" in labels


def test_drawing_leaves_the_callers_figures_alone(tmp_path):
    _plotting()
    import matplotlib.pyplot as plt

    from biotsr.analysis.render import render_clustermap

    mine = plt.figure()
    before = plt.get_fignums()
    samples = _blocks()
    spec = clustermap_spec(distance_matrix(feature_matrix(samples)), n_clusters=2)
    render_clustermap(spec, str(tmp_path / "x.png"))
    assert plt.get_fignums() == before and plt.gcf() is mine
    plt.close(mine)


def test_an_unlabelled_dataset_draws_with_a_cluster_stripe_only(tmp_path):
    _plotting()
    from biotsr.analysis.render import render_clustermap

    samples = [sample("a", {1: 3}), sample("b", {1: 2, 2: 1}), sample("c", {5: 1})]
    spec = clustermap_spec(distance_matrix(feature_matrix(samples)), n_clusters=2)
    fig = render_clustermap(spec, str(tmp_path / "small.svg"))
    assert (tmp_path / "small.svg").stat().st_size > 0
    _heat, stripe, legend = _axes(fig)
    assert [t.get_text() for t in stripe.get_xticklabels()] == ["cluster"]
    assert [t.get_text() for t in legend.get_texts()][0] == "cluster (maxclust 2)"


def test_a_percent_matrix_is_drawn_and_labelled_as_percent(tmp_path):
    _plotting()
    from biotsr.analysis.render import render_clustermap

    samples = _blocks()
    pct = distance_matrix(feature_matrix(samples)).scaled(100)
    spec = clustermap_spec(pct, [s.label for s in samples], n_clusters=2)
    fig = render_clustermap(spec)
    heat, _stripe, _legend = _axes(fig)
    assert heat.images[0].get_clim() == (0.0, 100.0)
    assert "generalized Jaccard distance (percent, 0-100)" in [
        a.get_ylabel() for a in fig.axes]


def test_more_classes_than_a_palette_has_never_share_a_colour():
    _plotting()
    import matplotlib

    from biotsr.analysis.render import _palette

    colours = _palette(matplotlib, "tab10", 23)
    assert len(set(colours)) == 23


def test_the_default_figure_needs_no_cut_and_draws_no_cluster_stripe(tmp_path):
    """ADR-099: tree and ordered heat map without any cluster count; the
    observed classes are the only stripe, and no partition is claimed."""
    _plotting()
    from biotsr.analysis.render import render_clustermap

    samples = _blocks()
    spec = clustermap_spec(distance_matrix(feature_matrix(samples)),
                           [s.label for s in samples])
    assert spec.assignment is None and spec.adjusted_rand is None
    fig = render_clustermap(spec, str(tmp_path / "uncut.png"))
    stripe = [a for a in fig.axes
              if [t.get_text() for t in a.get_xticklabels()][-1:] == ["observed class"]]
    assert len(stripe) == 1
    legend = fig.legends[0]
    assert not any(t.get_text().startswith("cluster") for t in legend.get_texts())
    assert "cluster" not in fig._suptitle.get_text()
    assert (tmp_path / "uncut.png").stat().st_size > 0
