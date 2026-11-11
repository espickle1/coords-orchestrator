"""Shared types: ``ResKey`` (a residue address) and the ``Status`` vocabulary.

A bare integer is never a residue address. Every structure residue that crosses a function
boundary is a ``ResKey``; every lookup that cannot produce a value returns a ``Status`` instead.
"""
from __future__ import annotations

from typing import NamedTuple


class ResKey(NamedTuple):
    """Address of one residue in a structure: author numbering, as deposited.

    ``chain`` is the author chain id (``auth_asym_id`` / ``pdb_strand_id``), ``auth_seq_id`` the
    author residue number, ``ins_code`` the insertion code ("" when absent).
    """

    model: int
    chain: str
    auth_seq_id: int
    ins_code: str = ""

    def __str__(self) -> str:  # 'A:195' or 'A:52A'
        return f"{self.chain}:{self.auth_seq_id}{self.ins_code}"


class Status:
    """Why a lookup produced no value. ``OK`` is the only status that carries a value."""

    OK = "ok"
    UNRESOLVED_SCHEME = "unresolved_scheme"      # a site whose numbering scheme is unknown
    NOT_ON_STRUCTURE = "not_on_structure"        # reference position has no residue in this chain
    NOT_OBSERVED = "not_observed"                # in SEQRES, no coordinates
    GAP_IN_SEQUENCE = "gap_in_sequence"          # f_i(p) undefined: deletion in the sequence
    OTHER_FAMILY = "other_family"                # neighbor on a chain belonging to another family
    UNREGISTERED_CHAIN = "unregistered_chain"    # neighbor on a chain with no map
    # A structure residue the reference has no position for: an expression tag, a fusion partner, an
    # insertion relative to the anchor. g is built from SEQRES, which can carry extras.
    NO_REF_POSITION = "no_ref_position"
    # The aligned residue is X (an unsequenced stretch). A per-residue table may have a row for it,
    # but that row describes an unknown token, not residue p.
    UNKNOWN_RESIDUE = "unknown_residue"

    ALL = (OK, UNRESOLVED_SCHEME, NOT_ON_STRUCTURE, NOT_OBSERVED, GAP_IN_SEQUENCE,
           OTHER_FAMILY, UNREGISTERED_CHAIN, NO_REF_POSITION, UNKNOWN_RESIDUE)


class HoleError(LookupError):
    """Raised instead of returning a hole when ``on_hole='raise'``."""
