"""coords-orchestrator: one reference frame for sequences, per-residue tables and protein structures.

Build the maps once from a ``coords.toml`` with ``python -m coords build coords.toml``, then read through
:class:`Frame`. See ``docs/`` and the README.
"""
__version__ = "0.1.0"

from .frame import Frame
from .types import HoleError, ResKey, Status

__all__ = ["Frame", "ResKey", "Status", "HoleError", "__version__"]
