"""A tour of coords-orchestrator on synthetic data: build a toy project, then read through the maps.

    python examples/walkthrough.py

Nothing here is real data. ``coords.testing.make_toy_project`` writes a project (references for two
groups, nine sequences with known indels, four structures, a site table) into a temporary directory.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

try:
    import coords  # noqa: F401
except ImportError:                                  # running from a checkout without installing
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from coords import Frame, ResKey
from coords.build import build_all
from coords.config import load_config
from coords.testing import make_toy_project

pd.set_option("display.width", 160)
pd.set_option("display.max_columns", 20)


def show(title: str) -> None:
    print(f"\n=== {title}")


with tempfile.TemporaryDirectory() as tmp:
    toy = make_toy_project(Path(tmp) / "toy")
    build_all(load_config(toy.config), verbose=False)       # same as: python -m coords build coords.toml
    fr = Frame.load(toy.maps)

    show("1. a deletion upstream moves every later row, and the deleted position is a hole")
    print("reference positions 29, 30, 31 in 'a_kin_del' (position 30 was deleted):")
    print(fr.rows("kin", "a_kin_del", [29, 30, 31]).tolist(), " <- the middle one is NA, not row 30")
    print(fr.row_table("kin", "a_kin_del", [29, 30, 31])[["ref_pos", "row", "ref_aa", "seq_aa", "status"]])

    show("2. an insertion upstream shifts rows the other way")
    print("positions 19, 20, 21 in 'a_kin_ins' (two residues inserted before 20):",
          fr.rows("kin", "a_kin_ins", [19, 20, 21]).tolist(), " (bare p - 1 would be [18, 19, 20])")

    show("3. an unsequenced run is a hole, but its flanks keep their rows")
    print(fr.rows("kin", "a_kin_x", [59, 60, 69, 70]).tolist(), " (positions 60-69 are X)")

    show("4. author numbers are not reference positions (offset changes inside the chain, plus an insertion code)")
    for p in (40, 41, 42):
        print(f"  reference position {p}: structure toy_a chain A residue", fr.struct("kin", "toy_a", p, chain="A")[0])
    print("  and back:", fr.ref_pos("toy_a", ResKey(1, "A", 49, "A")), "<- 49A is position 41, 49 is position 40")

    show("5. the same region numbered two ways in two structures")
    print("  toy_a :", fr.struct("kin", "toy_a", 20, chain="A")[0], "   toy_a2:", fr.struct("kin", "toy_a2", 20)[0])

    show("6. a missing loop is a hole on that chain only")
    print("  chain A:", fr.struct("kin", "toy_a", 62, chain="A"), "  all copies:", fr.struct("kin", "toy_a", 62))

    show("7. neighbours of a site, each with the row to read in every sequence")
    nb = fr.neighbors("kin", "toy_a", 30, 4.5, center="all", seq_ids=["a_kin_ref", "a_kin_ins"])
    print(nb[nb["is_query"] | (nb["neighbor_ref_pos"] == 31)]
          [["neighbor_reskey", "neighbor_ref_pos", "distance", "seq_id", "row", "seq_aa"]].to_string(index=False))

    show("8. sites in different numbering schemes land on one reference position")
    s = fr.sites("kin")
    print(s[["position", "scheme", "ref_pos", "ref_aa", "status", "phenotype"]].to_string(index=False))

    show("9. every hole was recorded, with its reason")
    print(fr.holes["status"].value_counts().to_string())
