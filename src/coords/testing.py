"""Synthetic data for tests and examples. Nothing here reads real data.

* :func:`write_mmcif`: a small mmCIF writer (both ``_pdbx_poly_seq_scheme`` variants) with the cases that
  break integer-passing code: unobserved residues, insertion codes, offsets that change inside a chain,
  alternate locations, nucleic acids, ligands and waters.
* :func:`make_toy_project`: a complete project (``coords.toml``, references, sequences, structures, sites)
  with known ground truth, used by this package's tests and by ``examples/``.
* :func:`synthetic_frame`: a :class:`~coords.Frame` built straight from sequences, for a downstream
  package's own tests (no files, no structures).
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from . import align as al
from .frame import SCHEMA_VERSION, Frame
from .seqio import sha256_seq, write_fasta
from .types import Status

AA = "ACDEFGHIKLMNPQRSTVWY"
THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY",
         "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
         "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}


# -----------------------------------------------------------------------------------------
# sequences
# -----------------------------------------------------------------------------------------
def random_protein(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(AA) for _ in range(n))


def homolog(seq: str, seed: int, substitution: float = 0.3,
            deletions: Sequence[Tuple[int, int]] = (), insertions: Sequence[Tuple[int, str]] = (),
            flank: int = 6) -> str:
    """A diverged copy: random substitutions, then deletions ``(start0, length)`` and insertions
    ``(index0, text)`` (both addressed in the *original* coordinates, applied back to front).

    Residues within ``flank`` of an indel are not substituted, so the indel has one best placement and
    a test can assert exact positions (with diverged flanks the aligner may tie between placements)."""
    rng = random.Random(seed)
    s = list(seq)
    edit_at = [i for i, _ in deletions] + [i + n - 1 for i, n in deletions] + [i for i, _ in insertions]
    for i, a in enumerate(s):
        if rng.random() < substitution and not any(abs(i - e) <= flank for e in edit_at):
            s[i] = rng.choice([x for x in AA if x != a])
    edits = [(i, "del", n) for i, n in deletions] + [(i, "ins", t) for i, t in insertions]
    for i, kind, arg in sorted(edits, key=lambda e: -e[0]):
        if kind == "del":
            del s[i:i + arg]
        else:
            s[i:i] = list(arg)
    return "".join(s)


# -----------------------------------------------------------------------------------------
# mmCIF writer
# -----------------------------------------------------------------------------------------
@dataclass
class ChainSpec:
    """One chain to write. ``auth[i]`` is the author number of SEQRES residue ``i`` (``None`` = no
    coordinates); ``ins[i]`` its insertion code. Residue ``i`` has CA at ``origin + (3.8 i, 0, 0)`` and a
    CB at ``CA + sidechain``, so distances between chains are set by ``origin``."""

    chain: str
    entity: str
    seq: str
    auth: Sequence[Optional[int]]
    ins: Sequence[str] = ()
    description: str = ""
    poly_type: str = "polypeptide(L)"
    origin: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    sidechain: Tuple[float, float, float] = (0.0, 2.0, 0.0)
    altloc_residues: Sequence[int] = ()            # SEQRES indices that get a second location (+50 A in x)
    ligand: bool = False                           # add one HETATM ligand and one water on this chain

    def ins_at(self, i: int) -> str:
        return self.ins[i] if i < len(self.ins) else ""


def _q(v: str) -> str:
    return "." if v == "" else v


def write_mmcif(path: Path, entry_id: str, chains: Sequence[ChainSpec], variant: str = "7col") -> Path:
    """Write ``chains`` as a minimal mmCIF. ``variant`` is the ``_pdbx_poly_seq_scheme`` schema:
    ``"7col"`` (no observed flag; unobserved residues have ``?`` for the number) or ``"12col"`` (an
    observed flag; unobserved residues keep an inferred number)."""
    if variant not in ("7col", "12col"):
        raise ValueError("variant must be '7col' or '12col'")
    lines = [f"data_{entry_id}", "#", f"_entry.id {entry_id}", "#"]
    entities: Dict[str, ChainSpec] = {}
    strands: Dict[str, List[str]] = {}
    for c in chains:
        entities.setdefault(c.entity, c)
        strands.setdefault(c.entity, []).append(c.chain)
    lines += ["loop_", "_entity.id", "_entity.type", "_entity.pdbx_description"]
    for eid, c in entities.items():
        lines.append(f"{eid} polymer '{c.description or 'chain ' + c.chain}'")
    if any(c.ligand for c in chains):
        lines += ["900 non-polymer NAG", "901 water water"]
    lines.append("#")
    lines += ["loop_", "_entity_poly.entity_id", "_entity_poly.type", "_entity_poly.pdbx_strand_id",
              "_entity_poly.pdbx_seq_one_letter_code_can"]
    for eid, c in entities.items():
        lines.append(f"{eid} '{c.poly_type}' {','.join(strands[eid])} {c.seq}")
    lines.append("#")

    if variant == "7col":
        lines += ["loop_", "_pdbx_poly_seq_scheme.entity_id", "_pdbx_poly_seq_scheme.asym_id",
                  "_pdbx_poly_seq_scheme.mon_id", "_pdbx_poly_seq_scheme.seq_id",
                  "_pdbx_poly_seq_scheme.pdb_strand_id", "_pdbx_poly_seq_scheme.pdb_seq_num",
                  "_pdbx_poly_seq_scheme.pdb_ins_code"]
    else:
        lines += ["loop_"] + [f"_pdbx_poly_seq_scheme.{k}" for k in (
            "asym_id", "entity_id", "seq_id", "mon_id", "ndb_seq_num", "pdb_seq_num", "auth_seq_num",
            "pdb_mon_id", "auth_mon_id", "pdb_strand_id", "pdb_ins_code", "hetero")]
    for c in chains:
        nucleic = not c.poly_type.startswith("polypeptide")
        for i, aa in enumerate(c.seq):
            mon = ("D" + aa) if nucleic else THREE[aa]
            a = c.auth[i]
            ins = _q(c.ins_at(i))
            asym = c.chain.lower()
            if variant == "7col":
                lines.append(f"{c.entity} {asym} {mon} {i + 1} {c.chain} {'?' if a is None else a} {ins}")
            else:
                inferred = i + 1 if a is None else a
                lines.append(f"{asym} {c.entity} {i + 1} {mon} {i + 1} {inferred} {'?' if a is None else a} "
                             f"{'?' if a is None else mon} {'?' if a is None else mon} {c.chain} {ins} n")
    lines.append("#")

    lines += ["loop_"] + [f"_atom_site.{k}" for k in (
        "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id", "label_asym_id",
        "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z",
        "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num")]
    serial = 0

    def atom(group, name, alt, comp, asym, ent, seqid, ins, xyz, auth, chain, elem):
        nonlocal serial
        serial += 1
        lines.append(f"{group} {serial} {elem} {name} {alt} {comp} {asym} {ent} {seqid} {ins} "
                     f"{xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f} {auth} {comp} {chain} {name} 1")

    for c in chains:
        nucleic = not c.poly_type.startswith("polypeptide")
        asym = c.chain.lower()
        for i, aa in enumerate(c.seq):
            a = c.auth[i]
            if a is None:
                continue
            ca = (c.origin[0] + 3.8 * i, c.origin[1], c.origin[2])
            ins = _q(c.ins_at(i))
            if nucleic:
                atom("ATOM", "P", ".", "D" + aa, asym, c.entity, i + 1, ins, ca, a, c.chain, "P")
                continue
            comp = THREE[aa]
            if i in c.altloc_residues:
                atom("ATOM", "CA", "A", comp, asym, c.entity, i + 1, ins, ca, a, c.chain, "C")
                atom("ATOM", "CA", "B", comp, asym, c.entity, i + 1, ins, (ca[0] + 50.0, ca[1], ca[2]), a, c.chain, "C")
            else:
                atom("ATOM", "CA", ".", comp, asym, c.entity, i + 1, ins, ca, a, c.chain, "C")
            if aa != "G":
                cb = (ca[0] + c.sidechain[0], ca[1] + c.sidechain[1], ca[2] + c.sidechain[2])
                atom("ATOM", "CB", ".", comp, asym, c.entity, i + 1, ins, cb, a, c.chain, "C")
        if c.ligand:
            atom("HETATM", "C1", ".", "NAG", asym, 900, ".", "?", (c.origin[0], c.origin[1] + 9.0, c.origin[2]),
                 301, c.chain, "C")
            atom("HETATM", "O", ".", "HOH", asym, 901, ".", "?", (c.origin[0] + 1.0, c.origin[1] + 9.0, c.origin[2]),
                 401, c.chain, "O")
    lines.append("#")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# -----------------------------------------------------------------------------------------
# the toy project
# -----------------------------------------------------------------------------------------
def toy_a_auth(p: int) -> Tuple[int, str]:
    """Author number of reference position ``p`` in structure ``toy_a``: offset 9, then an insertion
    code at p = 41 (49A), then offset 8. The offset changes inside the chain."""
    if p <= 40:
        return p + 9, ""
    if p == 41:
        return 49, "A"
    return p + 8, ""


@dataclass
class Toy:
    """Paths and ground truth of the toy project written by :func:`make_toy_project`."""

    root: Path
    config: Path
    maps: Path
    anchor_refs: Dict[str, str]          # family -> anchor (group A) reference
    group_b_refs: Dict[str, str]         # family -> group B reference
    seqs: Dict[str, str]                 # seq_id -> sequence in the sequence table
    far_kin: str                         # a ~40 % identical homologue of the anchor kin, in structure toy_far


KIN_LEN, BIND_LEN = 150, 90
MISSING_LOOP = range(60, 66)             # positions 60..65 unobserved in toy_a chain A
N_TERM_UNOBSERVED = range(1, 10)         # positions 1..9 unobserved in toy_a chain A


def _unambiguous_kin(start_seed: int = 11) -> str:
    """The anchor kinase. Deleting one of two identical neighbours gives the same string, so an aligner may
    place the gap on either; the toy's deletion (position 30) and insertion (before position 20) must
    have one best placement, so the first seed whose residues there are all distinct is used."""
    for seed in range(start_seed, start_seed + 200):
        k = random_protein(KIN_LEN, seed)
        if len({k[28], k[29], k[30]}) == 3 and k[18] != "W" and k[19] != "W" and k[17] != k[18]:
            return k
    raise RuntimeError("no unambiguous toy sequence found")


def make_toy_project(root: Path) -> Toy:
    """Write a complete, small project under ``root``. Group A is the anchor; group B a 70 %-identical
    homologue with a 3-residue deletion (positions 41-43) and a 2-residue insertion after position 101."""
    root = Path(root)
    (root / "refs").mkdir(parents=True, exist_ok=True)
    kin_a, bind_a = _unambiguous_kin(), random_protein(BIND_LEN, 1012)   # seeds must not collide: same seed = same prefix
    kin_b = homolog(kin_a, 13, 0.30, deletions=[(40, 3)], insertions=[(101, "GS")])
    bind_b = homolog(bind_a, 14, 0.25)
    far = homolog(kin_a, 15, 0.60)
    write_fasta(root / "refs" / "A.fasta", {"kin": kin_a, "bind": bind_a})
    write_fasta(root / "refs" / "B.fasta", {"kin": kin_b, "bind": bind_b})

    mutated_b = list(kin_b)
    for i in (20, 70, 120):
        mutated_b[i] = "W" if mutated_b[i] != "W" else "F"
    x_run = kin_a[:59] + "X" * 10 + kin_a[69:]
    seqs = {
        "a_kin_ref": kin_a,
        "a_kin_del": kin_a[:29] + kin_a[30:],                  # position 30 deleted
        "a_kin_ins": kin_a[:19] + "WW" + kin_a[19:],           # two residues inserted before position 20
        "a_kin_x": x_run,                                      # positions 60..69 unsequenced
        "a_kin_frag": kin_a[39:120],                           # positions 40..120 only
        "a_bind": bind_a,
        "b_kin_ref": kin_b,
        "b_kin_mut": "".join(mutated_b),
        "b_bind": bind_b,
    }
    write_fasta(root / "seqs.fasta", {"a_kin_del": seqs["a_kin_del"], "insertion_variant": seqs["a_kin_ins"]})
    groups = {"a_kin_ref": ("A", "kin"), "a_kin_del": ("A", "kin"), "a_kin_ins": ("A", "kin"),
              "a_kin_x": ("A", "kin"), "a_kin_frag": ("A", "kin"), "a_bind": ("A", "bind"),
              "b_kin_ref": ("B", "kin"), "b_kin_mut": ("B", "kin"), "b_bind": ("B", "bind")}
    rows = []
    for sid, (g, f) in groups.items():
        r = dict(seq_id=sid, group=g, family=f, sequence="", fasta="", record="", seq_hash="", source="toy")
        if sid == "a_kin_del":
            r.update(fasta="seqs.fasta", seq_hash=sha256_seq(seqs[sid]))
        elif sid == "a_kin_ins":
            r.update(fasta="seqs.fasta", record="insertion_variant")
        else:
            r.update(sequence=seqs[sid])
        rows.append(r)
    pd.DataFrame(rows).to_csv(root / "sequences.csv", index=False)

    pd.DataFrame([
        dict(family="kin", position=25, residue=kin_a[24], scheme="full", phenotype="x"),
        dict(family="kin", position=20, residue=kin_a[24], scheme="mature", phenotype="y"),   # 20 + 5 = 25
        dict(family="kin", position=10, residue="W", scheme="mature", phenotype="z"),
        dict(family="kin", position=10, residue="", scheme="", phenotype="blank"),
        dict(family="kin", position=5, residue="", scheme="tag", phenotype="t"),
        dict(family="kin", position=5, residue="", scheme="not_a_scheme", phenotype="u"),
        dict(family="bind", position=7, residue=bind_a[6], scheme="full", phenotype="b"),
        dict(family="bind", position=1000, residue="", scheme="full", phenotype="far"),
        dict(family="zzz", position=3, residue="", scheme="full", phenotype="nofam"),
    ]).to_csv(root / "sites.csv", index=False)

    # ---- structures
    sd = root / "structures"
    n = KIN_LEN
    a_auth = [None if (p in N_TERM_UNOBSERVED or p in MISSING_LOOP) else toy_a_auth(p)[0] for p in range(1, n + 1)]
    b_auth = [toy_a_auth(p)[0] for p in range(1, n + 1)]
    ins = [toy_a_auth(p)[1] for p in range(1, n + 1)]
    unrelated = random_protein(60, 99)
    write_mmcif(sd / "toy_a.cif", "TOYA", [
        ChainSpec("A", "1", kin_a, a_auth, ins, "kinase", origin=(0, 0, 0), sidechain=(0, 2.0, 0),
                  altloc_residues=[19], ligand=True),
        ChainSpec("B", "1", kin_a, b_auth, ins, "kinase", origin=(0, 5.0, 0), sidechain=(0, -2.0, 0)),
        ChainSpec("C", "2", bind_a, list(range(1, BIND_LEN + 1)), (), "binder", origin=(0, 0, 6.0)),
        ChainSpec("D", "3", unrelated, list(range(1, 61)), (), "antibody-like", origin=(0, -5.0, 0)),
        ChainSpec("E", "4", "ACGTACGT", list(range(1, 9)), (), "dna", poly_type="polydeoxyribonucleotide",
                  origin=(0, 20.0, 0)),
    ], variant="7col")
    write_mmcif(sd / "toy_a2.cif", "TOYB", [
        ChainSpec("A", "1", kin_a, [None, None] + [500 + p for p in range(3, n + 1)], (), "kinase"),
    ], variant="12col")
    write_mmcif(sd / "toy_b.cif", "TOYC", [
        ChainSpec("A", "1", kin_b, [100 + i + 1 for i in range(len(kin_b))], (), "kinase B"),
    ], variant="7col")
    write_mmcif(sd / "toy_far.cif", "TOYD", [
        ChainSpec("A", "1", far, list(range(1, len(far) + 1)), (), "far homologue"),
    ], variant="12col")

    cfg = root / "coords.toml"
    cfg.write_text("""[project]
name = "toy"
anchor = "A"
out_dir = "maps"

[references.A]
fasta = "refs/A.fasta"
label = "anchor group"

[references.B]
fasta = "refs/B.fasta"
label = "diverged group"

[sequences]
table = "sequences.csv"

[structures]
dir = "structures"

[sites]
table = "sites.csv"

[sites.schemes]
full = 0
mature = 5
tag = "unresolved"
""", encoding="utf-8")
    return Toy(root, cfg, root / "maps", {"kin": kin_a, "bind": bind_a}, {"kin": kin_b, "bind": bind_b}, seqs, far)


# -----------------------------------------------------------------------------------------
# a Frame without files
# -----------------------------------------------------------------------------------------
def synthetic_frame(anchor_refs: Dict[str, str], group_refs: Optional[Dict[str, Dict[str, str]]] = None,
                    sequences: Optional[Dict[str, Tuple[str, str, str]]] = None,
                    sites: Optional[pd.DataFrame] = None, anchor: str = "anchor",
                    on_hole: str = "skip") -> Frame:
    """A :class:`Frame` built from sequences alone (no structures, nothing on disk).

    anchor_refs   ``{family: sequence}`` of the anchor group.
    group_refs    ``{group: {family: sequence}}`` for other groups.
    sequences     ``{seq_id: (group, family, sequence)}``; the group is ``anchor`` or a key of ``group_refs``.
    """
    refs: Dict[Tuple[str, str], str] = {(anchor, f): s for f, s in anchor_refs.items()}
    for g, fams in (group_refs or {}).items():
        for f, s in fams.items():
            refs[(g, f)] = s
    ref2anchor = {}
    for (g, f), s in refs.items():
        ref2anchor[(g, f)] = np.arange(len(s)) if g == anchor else al.align_pair(refs[(anchor, f)], s).q_to_ref
    seq_maps, unknown, info = {}, {}, []
    for sid, (g, f, seq) in (sequences or {}).items():
        h = sha256_seq(seq)
        pa = al.align_pair(refs[(g, f)], seq)
        seq_maps[(g, f, h[:12])] = pa.ref_to_q
        unknown[(g, f, h[:12])] = np.array([j >= 0 and seq[j] == "X" for j in pa.ref_to_q])
        info.append(dict(seq_id=sid, group=g, family=f, seq_hash=h, length=len(seq), sequence=seq))
    seq_df = pd.DataFrame(info, columns=["seq_id", "group", "family", "seq_hash", "length", "sequence"])
    chains = pd.DataFrame({"structure": pd.Series(dtype=str), "chain": pd.Series(dtype=str),
                           "family": pd.Series(dtype=str), "group": pd.Series(dtype=str)})
    smaps = pd.DataFrame({"structure": pd.Series(dtype=str), "chain": pd.Series(dtype=str),
                          "auth_seq_id": pd.Series(dtype="Int64"), "ins_code": pd.Series(dtype=str),
                          "family": pd.Series(dtype=str), "ref_pos": pd.Series(dtype="Int64")})
    if sites is None:
        sites = pd.DataFrame(columns=["family", "position", "residue", "scheme", "ref_pos", "ref_aa",
                                      "status", "matches_ref"])
    sites = sites.copy()
    sites["ref_pos"] = sites["ref_pos"].astype("Int64")
    manifest = {"schema_version": SCHEMA_VERSION, "anchor": anchor,
                "groups": [anchor] + list(group_refs or {}), "families": list(anchor_refs)}
    return Frame(Path("."), anchor, refs, seq_df, seq_maps, ref2anchor, chains, smaps, sites, manifest,
                 on_hole=on_hole, seq_unknown=unknown)
