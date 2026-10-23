"""Serialization: the native, self-describing result format.

It is read as well as written: :func:`load_result` restores a saved result
with both scientific identities and their exact fingerprints.

``input_columns`` is deliberately absent from ``__all__``: it documents the
columns of the historical dataset sheets and is not a serialization format the
package offers.
"""

from . import native, native_read
from .native import save, write_keys, write_metadata, write_triplets
from .native_read import load_result, read_keys, read_metadata

__all__ = [
    "native",
    "native_read",
    "save",
    "write_keys",
    "write_metadata",
    "write_triplets",
    "load_result",
    "read_keys",
    "read_metadata",
]
