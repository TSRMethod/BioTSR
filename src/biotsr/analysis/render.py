"""Draw a :class:`ClustermapSpec` -- the one module in ``biotsr`` with a third-party import.

Everything else in the package -- key generation and the whole numerical
analysis layer -- runs on the standard library alone, and
``tests/unit/test_candidates_and_hygiene.py`` enforces it.  Drawing a figure
cannot, so it is optional:

* install it with ``pip install "biotsr[plot]"``, which installs exactly what
  this module imports: ``matplotlib`` (a test compares the two lists);
* nothing else in ``biotsr`` imports this module, so ``import biotsr`` and
  every analysis computation work without matplotlib; asking for a figure
  without it raises :class:`MissingRenderDependency` with the command that
  installs it.

The figure is drawn with matplotlib alone -- the heat map, both dendrograms,
the annotation stripes, the legends and the colour bar -- from the spec's own
linkage and leaf order, so the plotted order is the computed order by
construction.  (seaborn's ``clustermap`` would also need scipy at run time,
which seaborn does not install.)  It is drawn on a standalone
``matplotlib.figure.Figure``: pyplot is never imported, so the caller's backend,
current figure and open figures are never touched, and nothing needs closing.

Nothing here computes a number.  The matrix, the linkage, the leaf order, any
cluster labels, the palette slots and any ARI all arrive on the spec.
"""

from __future__ import annotations

from typing import Any, List, Optional, Sequence, Tuple

from .visualization import ClustermapSpec

__all__ = ["render_clustermap", "MissingRenderDependency", "INSTALL_HINT"]

#: What to run when the plotting stack is missing.
INSTALL_HINT = 'pip install "biotsr[plot]"'

_WHITE = (1.0, 1.0, 1.0)


class MissingRenderDependency(ImportError):
    """Raised when a figure is asked for and matplotlib is absent."""


def _imports():
    try:
        import matplotlib
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
        from matplotlib.patches import Patch
    except ImportError as exc:
        raise MissingRenderDependency(
            "drawing a clustermap needs matplotlib, which BioTSR installs only "
            "with its 'plot' extra. Install it with:  %s  (every number -- the "
            "matrix, linkage, leaf order, cluster labels and ARI -- is already on "
            "the ClustermapSpec without it). Missing: %s"
            % (INSTALL_HINT, getattr(exc, "name", None) or exc)
        ) from exc
    return matplotlib, Figure, Patch, FigureCanvasAgg


def _palette(matplotlib, name: str, count: int) -> List[Tuple[float, float, float]]:
    """``count`` distinct colours from a named matplotlib colormap.

    A qualitative map with enough entries gives its first ``count``; otherwise
    the map is sampled evenly.
    """
    cmap = matplotlib.colormaps[name]
    if getattr(cmap, "colors", None) is not None:
        if count <= cmap.N:
            return [tuple(cmap(i)[:3]) for i in range(count)]
        # More values than the qualitative map has colours: never repeat one.
        cmap = matplotlib.colormaps["turbo"]
    if count <= 1:
        return [tuple(cmap(0.0)[:3])]
    return [tuple(cmap(i / (count - 1))[:3]) for i in range(count)]


def _wrapped(text: str, width: int) -> str:
    """``text`` broken at spaces into lines of at most ``width`` characters."""
    lines, current = [], ""
    for word in text.split():
        if current and len(current) + 1 + len(word) > width:
            lines.append(current)
            current = word
        else:
            current = "%s %s" % (current, word) if current else word
    return "\n".join(lines + [current])


def _dendrogram_segments(rows: Sequence[Sequence[float]], leaves: Sequence[int],
                         n: int):
    """``[(positions, heights)]`` polylines of the tree; leaf ``leaves[r]`` sits at ``r + 0.5``."""
    position = {}
    height = {}
    for rank, leaf in enumerate(leaves):
        position[leaf] = rank + 0.5
        height[leaf] = 0.0
    segments = []
    for k, (left, right, merge, _size) in enumerate(rows):
        left, right = int(left), int(right)
        node = n + k
        position[node] = (position[left] + position[right]) / 2.0
        height[node] = merge
        segments.append((
            (position[left], position[left], position[right], position[right]),
            (height[left], merge, merge, height[right]),
        ))
    return segments


def render_clustermap(spec: ClustermapSpec, output_path: Optional[str] = None,
                      dpi: Optional[int] = None) -> Any:
    """Draw ``spec`` and, when ``output_path`` is given, save it there.

    The heat map is ``spec.matrix`` as given, rows and columns in the
    dendrogram's leaf order (``spec.ordered_sample_ids()``); its colour bar is
    labelled ``spec.matrix.value_label`` and spans the matrix's full range, 0 to
    its scale.  Row stripes appear only when there is something to show: the
    observed class labels (when some sample has one; an unlabelled sample is
    white) and the flat-cluster labels (only when the spec was cut by an
    explicit rule), each with its own legend.  Without either, the figure is
    the tree and the ordered heat map.  The file format follows the extension
    of ``output_path``.

    Returns the ``matplotlib.figure.Figure``, which belongs to the caller.
    """
    matplotlib, Figure, Patch, FigureCanvasAgg = _imports()

    ids = list(spec.sample_ids)
    n = len(ids)
    order = list(spec.hierarchy.leaves)
    ordered_ids = [ids[i] for i in order]
    square = spec.matrix.square()
    data = [[square[i][j] for j in order] for i in order]

    stripes = []   # (name, colour per sample in input order)
    handles = []
    if spec.has_classes:
        class_colours = _palette(matplotlib, spec.class_palette, len(spec.unique_classes))
        stripes.append(("observed class", [_WHITE if s is None else class_colours[s]
                                           for s in spec.row_class_slots()]))
        handles.append(Patch(facecolor="none", edgecolor="none", label="observed class"))
        handles.extend(Patch(facecolor=class_colours[i], label=name)
                       for name, i in spec.class_palette_index().items())
        if any(c is None for c in spec.class_labels):
            handles.append(Patch(facecolor=_WHITE, edgecolor="0.6", label="(no label)"))
    if spec.has_clusters:
        cluster_colours = _palette(matplotlib, spec.cluster_palette,
                                   len(spec.unique_clusters))
        stripes.append(("cluster", [cluster_colours[s] for s in spec.row_cluster_slots()]))
        handles.append(Patch(facecolor="none", edgecolor="none",
                             label="cluster (%s)" % spec.cut_description()))
        handles.extend(Patch(facecolor=cluster_colours[i], label=str(label))
                       for label, i in spec.cluster_palette_index().items())

    width, height = spec.figsize
    fig = Figure(figsize=(width, height))
    grid = fig.add_gridspec(
        2, 3,
        width_ratios=(0.15, 0.03 * max(len(stripes), 1), 1.0),
        height_ratios=(0.15, 1.0),
        wspace=0.02, hspace=0.02,
    )
    ax_heat = fig.add_subplot(grid[1, 2])
    ax_rows = fig.add_subplot(grid[1, 0], sharey=ax_heat)
    ax_stripe = fig.add_subplot(grid[1, 1], sharey=ax_heat) if stripes else None
    ax_cols = fig.add_subplot(grid[0, 2], sharex=ax_heat)

    image = ax_heat.imshow(
        data, cmap=spec.cmap, vmin=0.0, vmax=float(spec.matrix.scale),
        aspect="auto", interpolation="nearest", extent=(0, n, n, 0),
    )
    ticks = [r + 0.5 for r in range(n)]
    size = max(4.0, min(10.0, 400.0 / max(n, 1)))
    ax_heat.set_xticks(ticks)
    ax_heat.set_xticklabels(ordered_ids, rotation=90, fontsize=size)
    ax_heat.yaxis.tick_right()
    ax_heat.set_yticks(ticks)
    ax_heat.set_yticklabels(ordered_ids, fontsize=size)
    ax_heat.set_xlim(0, n)
    ax_heat.set_ylim(n, 0)

    if stripes:
        stripe_rgb = [[stripes[s][1][i] for s in range(len(stripes))] for i in order]
        ax_stripe.imshow(stripe_rgb, aspect="auto", interpolation="nearest",
                         extent=(0, len(stripes), n, 0))
        ax_stripe.set_xticks([s + 0.5 for s in range(len(stripes))])
        ax_stripe.set_xticklabels([name for name, _ in stripes], rotation=90,
                                  fontsize=8)
        ax_stripe.tick_params(left=False, labelleft=False)

    rows = spec.hierarchy.linkage.rows
    top = max([row[2] for row in rows] + [0.0]) or 1.0
    for xs, ys in _dendrogram_segments(rows, order, n):
        ax_cols.plot(xs, ys, color="0.25", linewidth=0.8)
        ax_rows.plot(ys, xs, color="0.25", linewidth=0.8)
    ax_cols.set_ylim(0, top * 1.05)
    ax_rows.set_xlim(top * 1.05, 0)
    for ax in (ax_cols, ax_rows):
        ax.set_axis_off()

    cax = fig.add_axes((0.02, 0.86, 0.02, 0.11))
    bar = fig.colorbar(image, cax=cax)
    bar.set_label(spec.matrix.value_label, fontsize=8)
    bar.ax.tick_params(labelsize=7)

    # Outside the heat map's tick labels, at the right edge; saving with
    # bbox_inches="tight" widens the image to keep it.
    if handles:
        renderer = FigureCanvasAgg(fig).get_renderer()
        labels_right = max((t.get_window_extent(renderer).x1
                            for t in ax_heat.get_yticklabels()), default=0.0)
        anchor = max(1.0, labels_right / fig.bbox.width + 0.01)
        fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(anchor, 0.86),
                   frameon=False, fontsize=8)

    if spec.title:
        title = spec.title
    elif spec.has_clusters:
        formed = spec.assignment.n_clusters
        title = "%s, average linkage; flat cut %s: %d cluster%s%s" % (
            spec.matrix.value_label, spec.cut_description(), formed,
            "" if formed == 1 else "s",
            "" if spec.adjusted_rand is None
            else "; ARI vs observed class %.3f" % spec.adjusted_rand,
        )
    else:
        title = "%s, average linkage" % spec.matrix.value_label
    # Above the plotting area, wrapped, so a long title (a cut with its ARI)
    # never runs into the colour bar at the top left.
    fig.suptitle(_wrapped(title, 72), fontsize=11, y=1.0, va="bottom")

    if output_path is not None:
        if dpi is None:
            fig.savefig(output_path, bbox_inches="tight")
        else:
            fig.savefig(output_path, bbox_inches="tight", dpi=dpi)
    return fig
