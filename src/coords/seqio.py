"""Sequence text helpers: FASTA reading and the sequence hash.

The hash is the sha256 of the normalised sequence string (upper-case, no whitespace, no trailing
``*``). A map is a function of that string, so it is also what a per-residue table must be checked
against before it is read through a map (:meth:`coords.Frame.check_hash`).
"""
from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Dict, Union


def normalize(seq: str) -> str:
    """Upper-case, strip whitespace and a trailing stop."""
    return re.sub(r"\s+", "", str(seq)).upper().rstrip("*")


def sha256_seq(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii", errors="replace")).hexdigest()


def read_fasta(path: Union[str, Path]) -> Dict[str, str]:
    """``{first header token: normalised sequence}``. A duplicate id keeps its first record."""
    out: Dict[str, str] = {}
    name, buf = None, []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None and name not in out:
                out[name] = "".join(buf)
            tok = line[1:].split()
            name, buf = (tok[0] if tok else ""), []
        elif name is not None:
            buf.append(line.strip())
    if name is not None and name not in out:
        out[name] = "".join(buf)
    return {k: normalize(v) for k, v in out.items()}


def write_fasta(path: Union[str, Path], records: Dict[str, str], width: int = 80) -> None:
    lines = []
    for name, seq in records.items():
        lines.append(f">{name}")
        lines.extend(seq[i:i + width] for i in range(0, len(seq), width))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
