"""Command line (``python -m coords``) and the plumbing scripts share to read through a Frame.

    python -m coords build coords.toml     # write the maps named by the config's out_dir
    python -m coords info maps             # summarise what was built and what was flagged

For a script that reads through a Frame:

    ap = argparse.ArgumentParser()
    coords.cli.add_frame_args(ap)
    args = ap.parse_args()
    fr = coords.cli.frame_from_args(args)
    ...
    coords.cli.write_holes(fr, out_dir / "stem_holes.csv")

There is no pass-through mode: ``--on-hole`` is ``skip`` (record and continue, the default) or ``raise``.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional, Union

import pandas as pd

from .config import ConfigError
from .frame import Frame


def add_frame_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    g = parser.add_argument_group("coordinate maps (coords)")
    g.add_argument("--coords-dir", type=Path, default=Path("maps"),
                   help="maps built by `python -m coords build` (the config's out_dir; default: ./maps)")
    g.add_argument("--on-hole", choices=("skip", "raise"), default="skip",
                   help="a position with no residue in a sequence/structure: skip the row and record it "
                        "(default) or raise. Holes are never filled with a neighbouring row.")
    return parser


def frame_from_args(args: argparse.Namespace) -> Frame:
    return Frame.load(coords_dir=args.coords_dir, on_hole=args.on_hole)


def write_holes(fr: Frame, path: Union[str, Path]) -> pd.DataFrame:
    """Write ``Frame.holes`` (deduplicated, with counts) next to a script's outputs and return it."""
    holes = fr.holes
    path = Path(path)
    if len(holes):
        key = list(holes.columns)
        out = holes.astype(str).groupby(key, dropna=False).size().reset_index(name="n")
        out = out.sort_values(["status", "n"], ascending=[True, False])
    else:
        out = pd.DataFrame(columns=["status", "n"])
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    return out


def _info(maps_dir: Path) -> int:
    mf = maps_dir / "manifest.json"
    if not mf.is_file():
        print(f"{maps_dir}: no manifest.json; run `python -m coords build` first", file=sys.stderr)
        return 2
    m = json.loads(mf.read_text(encoding="utf-8"))
    print(f"project   {m['project']}   anchor {m['anchor']}   groups {', '.join(m['groups'])}")
    print(f"families  {', '.join(m['families'])}")
    print(f"built     {m['built_utc']}   coords {m['coords_version']}   biopython {m['aligner']['biopython']}")
    for k, v in m["counts"].items():
        print(f"  {k:24s} {v}")
    flagged = [c for c in m["per_chain"] if c["quality"] != "ok"]
    if flagged:
        print("flagged chains (quality != ok):")
        for c in flagged:
            print(f"  {c['structure']}:{c['chain']}  family={c['family'] or '-'}  quality={c['quality']}  "
                  f"anchor_identity={c['anchor_identity']}")
    seq = maps_dir / "sequences.csv"
    if seq.is_file():
        s = pd.read_csv(seq, keep_default_na=False, dtype=str)
        bad = s[s["quality"] != "ok"] if "quality" in s.columns else s.iloc[0:0]
        if len(bad):
            print("flagged sequences (quality != ok):")
            for r in bad.itertuples():
                print(f"  {r.seq_id}  quality={r.quality}  identity={r.identity}  coverage={r.coverage_of_ref}")
    return 0


def main(argv: Optional[List[str]] = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m coords", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="build the maps from a coords.toml")
    b.add_argument("config", type=Path)
    b.add_argument("-q", "--quiet", action="store_true")
    i = sub.add_parser("info", help="summarise a built maps directory and what was flagged")
    i.add_argument("maps_dir", type=Path)
    args = ap.parse_args(argv)
    if args.cmd == "build":
        from .build import build_all
        from .config import load_config
        try:
            build_all(load_config(args.config), verbose=not args.quiet)
        except ConfigError as e:
            print(f"error: {e}", file=sys.stderr)
            sys.exit(2)
    elif args.cmd == "info":
        sys.exit(_info(args.maps_dir))


if __name__ == "__main__":
    main()
