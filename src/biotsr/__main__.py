"""``python -m biotsr`` -- the same entry point as the ``biotsr`` console script.

Offered because a console script is only on ``PATH`` when the environment's
``bin`` directory is, and on a cluster it often is not.  ``python -m biotsr`` works
from any interpreter that can import the package, which is the invariant an HPC
job script can actually rely on.
"""

from __future__ import annotations

import sys

from .cli.main import main

if __name__ == "__main__":
    sys.exit(main())
