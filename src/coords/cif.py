"""Minimal mmCIF reader for what the coordinate maps need.

Reads five categories: ``_entry``, ``_entity``, ``_entity_poly``, ``_pdbx_poly_seq_scheme`` and
``_atom_site``. Pure Python, no gemmi/biotite dependency. Protein chains only: nucleic-acid
chains, ligands and waters are read past and dropped.

Two schema variants of ``_pdbx_poly_seq_scheme`` occur in RCSB mmCIF files:

* 12-column (``asym_id entity_id seq_id mon_id ndb_seq_num pdb_seq_num auth_seq_num pdb_mon_id
  auth_mon_id pdb_strand_id pdb_ins_code hetero``): ``pdb_mon_id`` is ``?`` for a residue with no
  coordinates, so observation is read from the file.
* 7-column (``entity_id asym_id mon_id seq_id pdb_strand_id pdb_seq_num pdb_ins_code``): no
  observed flag, so observation comes from intersecting with ``_atom_site``.

Both are normalised to the same frame by :func:`read_structure`.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O", "HYP": "P", "SEP": "S", "TPO": "T", "PTR": "Y",
    "CSO": "C", "CME": "C", "KCX": "K", "LLP": "K", "MLY": "K", "PCA": "Q",
}
_TOKEN = re.compile(r"""'(?:[^']|'(?=\S))*'(?=\s|$)|"(?:[^"]|"(?=\S))*"(?=\s|$)|\S+""")
_NULL = {".", "?"}
WANTED = ("_entry", "_entity", "_entity_poly", "_pdbx_poly_seq_scheme", "_atom_site")


def _unquote(tok: str) -> str:
    if len(tok) >= 2 and tok[0] == tok[-1] and tok[0] in "'\"":
        return tok[1:-1]
    return tok


def _blocks(text: str, wanted: Sequence[str]) -> Iterator[Tuple[str, List[str], List[List[str]]]]:
    """Yield ``(category, item_names, rows)`` for each wanted category in the first data block.

    Key-value categories come back as a single row. Semicolon text fields become one token.
    """
    lines = text.split("\n")
    n = len(lines)
    i = 0
    wanted_set = set(wanted)
    kv: Dict[str, Dict[str, str]] = {}

    def read_text_field(i: int) -> Tuple[str, int]:
        """lines[i] starts with ';'. Return (text, index after the closing ';' line)."""
        buf = [lines[i][1:]]
        i += 1
        while i < n and not lines[i].startswith(";"):
            buf.append(lines[i])
            i += 1
        return "\n".join(buf).strip(), i + 1

    while i < n:
        line = lines[i]
        if line.startswith("data_") and kv:
            break                                    # first block only
        if line.startswith("loop_"):
            i += 1
            tags: List[str] = []
            while i < n and lines[i].startswith("_"):
                tags.append(lines[i].split()[0])
                i += 1
            cat = tags[0].split(".")[0] if tags else ""
            if cat not in wanted_set:
                while i < n and not (lines[i].startswith("loop_") or lines[i].startswith("_")
                                     or lines[i].startswith("data_") or lines[i].startswith("#")):
                    if lines[i].startswith(";"):
                        _, i = read_text_field(i)
                    else:
                        i += 1
                continue
            toks: List[str] = []
            while i < n:
                l = lines[i]
                if l.startswith("loop_") or l.startswith("_") or l.startswith("data_") or l.startswith("#"):
                    break
                if l.startswith(";"):
                    t, i = read_text_field(i)
                    toks.append(t)
                    continue
                toks.extend(_unquote(t) for t in _TOKEN.findall(l))
                i += 1
            k = len(tags)
            if k and len(toks) % k:
                raise ValueError(f"truncated loop for {cat}")
            yield cat, [t.split(".", 1)[1] for t in tags], [toks[j:j + k] for j in range(0, len(toks), k)]
            continue
        if line.startswith("_"):
            parts = line.split(None, 1)
            tag = parts[0]
            cat, _, item = tag.partition(".")
            if len(parts) == 2:
                val = _unquote(parts[1].strip())
                i += 1
            else:
                i += 1
                if i < n and lines[i].startswith(";"):
                    val, i = read_text_field(i)
                elif i < n:
                    val = _unquote(lines[i].strip())
                    i += 1
                else:
                    val = ""
            if cat in wanted_set:
                kv.setdefault(cat, {})[item] = val
            continue
        i += 1
    for cat, d in kv.items():
        yield cat, list(d), [list(d.values())]


def _frame(items: List[str], rows: List[List[str]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=items, dtype=object)


def _int_or_na(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.where(~s.isin(_NULL)), errors="coerce").astype("Int64")


@dataclass
class Structure:
    """What a map needs from one CIF."""

    structure: str                         # the key callers use: the file stem unless overridden
    path: Path
    entities: pd.DataFrame                 # entity_id, type, description, strand_ids, seq_can
    scheme: pd.DataFrame                   # model-independent SEQRES table, see read_structure
    atoms: pd.DataFrame                    # protein atoms, first model, one altloc
    scheme_variant: str = ""               # "12col" | "7col"
    notes: List[str] = field(default_factory=list)
    entry_id: str = ""                     # ``_entry.id`` (the PDB id for RCSB files), "" if absent


def read_structure(path: Union[str, Path], structure: Optional[str] = None) -> Structure:
    """Parse one mmCIF into entities, a normalised SEQRES scheme and a protein atom table.

    ``scheme`` columns: entity_id, asym_id (label), chain (author), seq_id (1-based SEQRES index),
    mon_id, aa (one letter, "X" if non-standard), auth_seq_id (Int64, NA if none), ins_code,
    observed (bool). ``atoms`` columns: model, chain (author), asym_id, auth_seq_id, ins_code,
    comp_id, aa, atom, x, y, z, element.
    """
    path = Path(path)
    text = path.read_text()
    cats = {c: _frame(items, rows) for c, items, rows in _blocks(text, WANTED)}
    structure = structure or path.stem

    # ---- entities ---------------------------------------------------------------------
    ent = cats.get("_entity", pd.DataFrame(columns=["id", "type", "pdbx_description"]))
    ep = cats.get("_entity_poly", pd.DataFrame(columns=["entity_id", "pdbx_strand_id"]))
    seq_col = "pdbx_seq_one_letter_code_can" if "pdbx_seq_one_letter_code_can" in ep.columns else None
    entities = ent.rename(columns={"id": "entity_id", "pdbx_description": "description"})
    entities = entities[[c for c in ("entity_id", "type", "description") if c in entities.columns]]
    ep2 = ep.rename(columns={"pdbx_strand_id": "strand_ids"})
    if "type" in ep2.columns:
        ep2 = ep2.rename(columns={"type": "poly_type"})
    keep = ["entity_id"] + [c for c in ("strand_ids", "poly_type") if c in ep2.columns]
    if seq_col:
        ep2["seq_can"] = ep[seq_col].astype(str).str.replace(r"\s+", "", regex=True)
        keep.append("seq_can")
    entities = entities.merge(ep2[keep], on="entity_id", how="left")

    # ---- atoms (protein, first model, one altloc) -------------------------------------
    a = cats.get("_atom_site")
    if a is None or a.empty:
        raise ValueError(f"{path}: no _atom_site loop")
    col = lambda *names: next((c for c in names if c in a.columns), None)  # noqa: E731
    c_group = col("group_PDB")
    c_chain = col("auth_asym_id", "label_asym_id")
    c_seq = col("auth_seq_id", "label_seq_id")
    keep_mask = np.ones(len(a), dtype=bool)
    if c_group:
        keep_mask &= a[c_group].isin(["ATOM", "HETATM"]).to_numpy()
    atoms = pd.DataFrame({
        "model": pd.to_numeric(a[col("pdbx_PDB_model_num")] if col("pdbx_PDB_model_num") else "1",
                               errors="coerce").fillna(1).astype(int),
        "chain": a[c_chain].astype(str),
        "asym_id": a[col("label_asym_id")].astype(str) if col("label_asym_id") else a[c_chain].astype(str),
        "auth_seq_id": _int_or_na(a[c_seq]),
        "ins_code": a[col("pdbx_PDB_ins_code")].where(~a[col("pdbx_PDB_ins_code")].isin(_NULL), "")
        if col("pdbx_PDB_ins_code") else "",
        "comp_id": a[col("auth_comp_id", "label_comp_id")].astype(str).str.upper(),
        "atom": a[col("auth_atom_id", "label_atom_id")].astype(str),
        "altloc": a[col("label_alt_id")].where(~a[col("label_alt_id")].isin(_NULL), "")
        if col("label_alt_id") else "",
        "x": pd.to_numeric(a[col("Cartn_x")]), "y": pd.to_numeric(a[col("Cartn_y")]),
        "z": pd.to_numeric(a[col("Cartn_z")]),
        "element": a[col("type_symbol")].astype(str) if col("type_symbol") else "",
        "group": a[c_group] if c_group else "ATOM",
    })[keep_mask].reset_index(drop=True)
    atoms["aa"] = atoms["comp_id"].map(AA3)
    # protein residues only: standard residue names, or any ATOM-record polymer residue
    atoms = atoms[atoms["aa"].notna() | (atoms["group"] == "ATOM")].copy()
    atoms["aa"] = atoms["aa"].fillna("X")
    atoms = atoms[atoms["auth_seq_id"].notna()]
    # keep polypeptide chains only: nucleic acids carry ATOM records too
    if "poly_type" in entities.columns and "_pdbx_poly_seq_scheme" in cats:
        sch = cats["_pdbx_poly_seq_scheme"]
        etype = dict(zip(entities["entity_id"], entities["poly_type"].fillna("")))
        prot_chains = {str(c) for c, e in zip(sch["pdb_strand_id"], sch["entity_id"])
                       if str(etype.get(e, "")).startswith("polypeptide")}
        atoms = atoms[atoms["chain"].isin(prot_chains)]
    first = int(atoms["model"].min())
    atoms = atoms[atoms["model"] == first]
    rank = atoms["altloc"].map(lambda v: 0 if v in ("", ".") else (1 if v == "A" else 2))
    atoms = atoms.assign(_rank=rank, _order=np.arange(len(atoms)))
    # one altloc per (residue, atom): blank, else A, else the first remaining
    atoms = (atoms.sort_values(["_rank", "_order"])
             .drop_duplicates(["model", "chain", "auth_seq_id", "ins_code", "atom"], keep="first")
             .sort_values("_order").drop(columns=["_rank", "_order", "altloc", "group"]))
    atoms["auth_seq_id"] = atoms["auth_seq_id"].astype(int)
    atoms = atoms.reset_index(drop=True)

    # ---- SEQRES scheme ----------------------------------------------------------------
    s = cats.get("_pdbx_poly_seq_scheme")
    if s is None or s.empty:
        raise ValueError(f"{path}: no _pdbx_poly_seq_scheme loop")
    variant = "12col" if "pdb_mon_id" in s.columns else "7col"
    scheme = pd.DataFrame({
        "entity_id": s["entity_id"],
        "asym_id": s["asym_id"],
        "chain": s["pdb_strand_id"],
        "seq_id": _int_or_na(s["seq_id"]).astype(int),
        "mon_id": s["mon_id"].str.upper(),
        "auth_seq_id": _int_or_na(s["pdb_seq_num"]),
        "ins_code": s["pdb_ins_code"].where(~s["pdb_ins_code"].isin(_NULL), ""),
    })
    scheme["aa"] = scheme["mon_id"].map(AA3).fillna("X")
    if "poly_type" in entities.columns:
        etype = dict(zip(entities["entity_id"], entities["poly_type"].fillna("")))
        scheme["poly_type"] = scheme["entity_id"].map(etype).fillna("")
    else:
        scheme["poly_type"] = ""
    if variant == "12col":
        scheme["observed"] = ~s["pdb_mon_id"].isin(_NULL).to_numpy()
        # auth_seq_num is "?" for unobserved residues; pdb_seq_num is still the inferred number.
    else:
        scheme["observed"] = False                       # filled from atoms below
    first_atoms = atoms[["chain", "auth_seq_id", "ins_code"]].drop_duplicates()
    first_atoms["_in_atoms"] = True
    scheme = scheme.merge(first_atoms, left_on=["chain", "auth_seq_id", "ins_code"],
                          right_on=["chain", "auth_seq_id", "ins_code"], how="left")
    in_atoms = scheme.pop("_in_atoms").eq(True)
    notes: List[str] = []
    if variant == "7col":
        scheme["observed"] = in_atoms
    else:
        disagree = int((scheme["observed"] != in_atoms).sum())
        if disagree:
            notes.append(f"{disagree} scheme rows disagree between pdb_mon_id and _atom_site")
    entry = cats.get("_entry")
    entry_id = str(entry["id"].iloc[0]) if entry is not None and "id" in entry.columns and len(entry) else ""
    return Structure(structure, path, entities, scheme, atoms, variant, notes, entry_id)


def seqres_by_chain(st: Structure) -> Dict[str, Tuple[str, pd.DataFrame]]:
    """``{chain: (SEQRES one-letter string, scheme rows in seq_id order)}`` for author chain ids."""
    out: Dict[str, Tuple[str, pd.DataFrame]] = {}
    for chain, g in st.scheme.groupby("chain", sort=False):
        g = g.sort_values("seq_id")
        out[str(chain)] = ("".join(g["aa"]), g.reset_index(drop=True))
    return out
