"""Modifier vocabulary and the pipeline hooks modifiers attach to.

At the user-facing layer every modifier is the same kind of thing: an optional
keyword on one shared ``analyze(...)`` call.  Internally they are **not** one
abstraction, and collapsing them would misrepresent the science.  Each acts at a
different, named point in the pipeline, and that is exactly what keeps them from
becoming one leaky hook:

``` text
Structure
  -> Selection
  -> Candidate strategy
  -> LABELS            <- aa_grouping (a label-scheme swap)
  -> Canonical ordering
  -> Geometry
  -> CANDIDATE_FILTER  <- size_filter (a predicate on max_distance)
  -> Binning
  -> Encoding
  -> KEY_SIGN          <- mirror (a sign convention keyed off theta_raw)
  -> ANNOTATION        <- sse (a residue-level annotation naming the bucket)
  -> Result
```

The ordering is forced by data dependencies, not by taste:

* AA grouping must act **before** canonical ordering, because ordering compares
  labels.  A consequence worth stating loudly: grouping can therefore change
  which vertex becomes the apex, which changes ``theta``, which changes the
  mirror sign.  Grouping and mirror are not independent -- composing them is not
  the same as applying each to the other's output.
* Size filtering needs ``max_distance``, so it cannot run before geometry.
* Mirror needs ``theta_raw``, so it can only apply once geometry is done; it
  signs a key whose magnitude is already determined.
* SSE is last on purpose, and it is the only entry that is **downstream of the
  key rather than upstream of it**. Everything before ANNOTATION can change a
  key; nothing at ANNOTATION can, because by then the key -- sign included -- is
  final and the annotator is handed vertex identities and nothing else.

Which side of the two identities each lands on:

| Modifier | Changes key meaning? | Changes which triangles? | Identity |
|---|---|---|---|
| ``aa_grouping`` | **yes** (label scheme + radix) | no | ``KeySpace`` |
| ``mirror`` | **yes** (sign convention) | no | ``KeySpace`` |
| ``size_filter`` | no | **yes** | ``SamplingProfile`` |
| ``sse`` | **no** | **no** | *neither* |

``sse`` is the only entry in that table that lands on neither identity, and that
is the whole scientific claim of SSE support: annotating a run does not change what
a key means and does not change which triangles were drawn, so an SSE result and
a plain result over the same structure share a ``KeySpace`` fingerprint, share a
``SamplingProfile`` fingerprint, and are directly comparable. It is therefore not
a modifier in the sense the three above are -- it is listed in
``supported_modifiers`` purely so that requesting it on a family where the 18
categories are undefined fails loudly, naming what *is* supported (ADR-064).
"""

from __future__ import annotations

from ..domain.errors import IncompatibleRepresentation

__all__ = ["ModifierHook", "ModifierName", "IncompatibleModifier"]


class IncompatibleModifier(IncompatibleRepresentation):
    """Raised when a modifier cannot apply to the points in play.

    Always names what *is* supported.  Never a silent no-op -- a modifier that
    quietly does nothing is worse than one that refuses.
    """


class ModifierHook:
    """The named pipeline point a modifier attaches to."""

    LABELS = "labels"
    CANDIDATE_FILTER = "candidate_filter"
    KEY_SIGN = "key_sign"
    #: Downstream of the key: names the frequency bucket, moves nothing.
    ANNOTATION = "annotation"


class ModifierName:
    """The user-facing modifier names, as they appear in ``analyze(...)``."""

    AA_GROUPING = "aa_grouping"
    MIRROR = "mirror"
    SIZE_FILTER = "size_filter"
    SSE = "sse"
