"""``Frame``: the one place that holds the reference integer and applies the maps before every read.

    fr   = Frame.load("maps")                              # the out_dir of coords.toml
    rows = fr.rows("kinase", seq_id, [141, 150])           # 0-based table rows, holes as pd.NA
    keys = fr.struct("kinase", "1abc", 141)                # [ResKey, ...] one per copy
    nb   = fr.neighbors("kinase", "1abc", 141, 8.0, seq_ids=[...])

Rules:

1. Every position that enters or leaves is a reference position ``p`` on the anchor group's reference of a
   family. ``q_i`` (string position) and ``a`` (author number) exist only inside a read.
2. Maps are the CSVs written by ``python -m coords build``; the Frame only merges against them.
3. A hole is never filled. A lookup that cannot return a value records a :class:`Status` in
   ``Frame.holes``; with ``on_hole="raise"`` it raises :class:`HoleError` instead. There is no
   pass-through mode.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from . import align as al
from .seqio import read_fasta, sha256_seq
from .types import HoleError, ResKey, Status

Number = Union[int, np.integer]
SCHEMA_VERSION = 1
_SEQ_NUMERIC_INT = ("length", "ref_length", "n_aligned", "n_unknown", "n_inserted", "n_gap")
_SEQ_NUMERIC_FLOAT = ("identity", "coverage_of_ref")


class Frame:
    def __init__(self, coords_dir: Path, anchor: str, refs: Dict[Tuple[str, str], str],
                 sequences: pd.DataFrame, seq_maps: Dict[tuple, np.ndarray],
                 ref2anchor: Dict[Tuple[str, str], np.ndarray], chains: pd.DataFrame, smaps: pd.DataFrame,
                 sites: pd.DataFrame, manifest: dict, on_hole: str = "skip",
                 seq_unknown: Optional[Dict[tuple, np.ndarray]] = None,
                 structures: Optional[pd.DataFrame] = None):
        if on_hole not in ("skip", "raise"):
            raise ValueError("on_hole must be 'skip' or 'raise'")
        self.coords_dir = Path(coords_dir)
        self.anchor = anchor
        self.refs = refs
        self.sequences = sequences.set_index("seq_id", drop=False)
        self._seq_maps = seq_maps
        self._seq_unknown = seq_unknown or {}
        self._ref2anchor = ref2anchor
        self._anchor2own = {k: al.invert(v, len(refs[(anchor, k[1])])) for k, v in ref2anchor.items()}
        self.chains = chains
        self.smaps = smaps
        self.sites_table = sites
        self.manifest = manifest
        self.structures = structures if structures is not None else pd.DataFrame(
            columns=["structure", "entry_id", "path", "sha256"])
        self.on_hole = on_hole
        self._holes: List[dict] = []
        self._atoms: Dict[str, pd.DataFrame] = {}
        self._seq_cache: Dict[str, Optional[str]] = {}
        self._build_struct_index()

    # ------------------------------------------------------------------ construction
    @classmethod
    def load(cls, coords_dir: Union[str, Path] = "maps", on_hole: str = "skip",
             anchor: Optional[str] = None) -> "Frame":
        """Read the maps written by ``python -m coords build`` (the config's ``out_dir``)."""
        d = Path(coords_dir)
        need = ["manifest.json", "sequences.csv", "sequence_maps.csv", "reference_maps.csv",
                "structure_chains.csv", "structure_maps.csv", "structures.csv", "sites_resolved.csv"]
        missing = [n for n in need if not (d / n).exists()]
        if missing:
            raise FileNotFoundError(f"{d}: missing {missing}. Run `python -m coords build <config>` first.")
        manifest = json.loads((d / "manifest.json").read_text(encoding="utf-8"))
        if manifest.get("schema_version") != SCHEMA_VERSION:
            raise ValueError(f"{d}: maps have schema_version {manifest.get('schema_version')!r}, this "
                             f"coords expects {SCHEMA_VERSION}. Rebuild.")
        if anchor is not None and manifest["anchor"] != anchor:
            raise ValueError(f"maps were built with anchor {manifest['anchor']!r}, not {anchor!r}")
        anchor = manifest["anchor"]
        refs: Dict[Tuple[str, str], str] = {}
        for f in sorted((d / "reference").glob("*/*.fasta")):
            refs[(f.parent.name, f.stem)] = next(iter(read_fasta(f).values()))

        sequences = pd.read_csv(d / "sequences.csv", keep_default_na=False, na_values=[""], dtype=str)
        for c in _SEQ_NUMERIC_INT:
            if c in sequences.columns:
                sequences[c] = pd.to_numeric(sequences[c]).astype("Int64")
        for c in _SEQ_NUMERIC_FLOAT:
            if c in sequences.columns:
                sequences[c] = pd.to_numeric(sequences[c])
        sm = pd.read_csv(d / "sequence_maps.csv", keep_default_na=False, na_values=[""],
                         usecols=["group", "family", "hash12", "ref_pos", "seq_pos", "status"],
                         dtype={"seq_pos": "Int64", "ref_pos": "int64", "hash12": str})
        seq_maps: Dict[tuple, np.ndarray] = {}
        seq_unknown: Dict[tuple, np.ndarray] = {}
        for (group, family, h), g in sm.groupby(["group", "family", "hash12"], sort=False):
            arr = np.full(len(refs[(group, family)]), -1, dtype=np.int64)
            unk = np.zeros(len(arr), dtype=bool)
            ok = g["seq_pos"].notna().to_numpy()
            arr[g["ref_pos"].to_numpy()[ok] - 1] = g["seq_pos"].to_numpy(dtype="float")[ok].astype(np.int64) - 1
            unk[g["ref_pos"].to_numpy() - 1] = (g["status"] == Status.UNKNOWN_RESIDUE).to_numpy()
            seq_maps[(group, family, h)] = arr
            seq_unknown[(group, family, h)] = unk

        rm = pd.read_csv(d / "reference_maps.csv", keep_default_na=False, na_values=[""], dtype={"to_pos": "Int64"})
        ref2anchor: Dict[Tuple[str, str], np.ndarray] = {}
        for (g, f), seq in refs.items():
            if g == anchor:
                ref2anchor[(g, f)] = np.arange(len(seq))
        for (family, group), g in rm.groupby(["family", "from_group"], sort=False):
            arr = np.full(len(refs[(group, family)]), -1, dtype=np.int64)
            ok = g["to_pos"].notna().to_numpy()
            arr[g["from_pos"].to_numpy()[ok] - 1] = g["to_pos"].to_numpy(dtype="float")[ok].astype(np.int64) - 1
            ref2anchor[(group, family)] = arr

        chains = pd.read_csv(d / "structure_chains.csv", keep_default_na=False, na_values=[""],
                             dtype={"chain": str, "family": str, "group": str, "structure": str})
        chains["family"] = chains["family"].fillna("")
        smaps = pd.read_csv(d / "structure_maps.csv", keep_default_na=False, na_values=[""],
                            dtype={"auth_seq_id": "Int64", "ref_pos": "Int64", "ins_code": str, "chain": str,
                                   "family": str, "group_ref_pos": "Int64", "structure": str})
        smaps["ins_code"] = smaps["ins_code"].fillna("")
        smaps["family"] = smaps["family"].fillna("")
        sites = pd.read_csv(d / "sites_resolved.csv", keep_default_na=False, na_values=[""],
                            dtype={"ref_pos": "Int64"})
        for c in ("residue", "scheme", "ref_aa"):            # text columns: an empty cell is "", not NaN
            if c in sites.columns:
                sites[c] = sites[c].fillna("").astype(str)
        structures = pd.read_csv(d / "structures.csv", keep_default_na=False, na_values=[""], dtype=str)
        return cls(d, anchor, refs, sequences, seq_maps, ref2anchor, chains, smaps, sites, manifest, on_hole,
                   seq_unknown, structures)

    def _build_struct_index(self) -> None:
        s = self.smaps
        self._by_reskey: Dict[tuple, int] = {}
        self._by_pos: Dict[tuple, List[int]] = defaultdict(list)
        ok_auth = s["auth_seq_id"].notna().to_numpy()
        cols = (s["structure"].to_numpy(), s["chain"].to_numpy(), s["auth_seq_id"].to_numpy(dtype="float"),
                s["ins_code"].to_numpy(), s["family"].to_numpy(), s["ref_pos"].to_numpy(dtype="float"))
        for i in range(len(s)):
            st, ch, auth, ins, fam, rp = (c[i] for c in cols)
            if ok_auth[i]:
                self._by_reskey[(st, ch, int(auth), ins)] = i
            if fam and rp == rp:
                self._by_pos[(st, fam, int(rp))].append(i)
        self._chain_order = {(r.structure, str(r.chain)): k for k, r in enumerate(self.chains.itertuples())}
        self._chain_family = {(r.structure, str(r.chain)): (r.family or "") for r in self.chains.itertuples()}

    # ------------------------------------------------------------------ introspection
    @property
    def groups(self) -> List[str]:
        return list(self.manifest.get("groups", sorted({g for g, _ in self.refs})))

    @property
    def families(self) -> List[str]:
        return list(self.manifest.get("families", sorted({f for g, f in self.refs if g == self.anchor})))

    # ------------------------------------------------------------------ holes
    def record_hole(self, status: str, **info) -> None:
        """Record a hole (or raise, under ``on_hole='raise'``). Callers that find a hole of their own
        use this so it lands in the same log as the Frame's."""
        rec = {"status": status, **info}
        if self.on_hole == "raise":
            raise HoleError(f"{status}: {info}")
        self._holes.append(rec)

    _hole = record_hole

    @property
    def holes(self) -> pd.DataFrame:
        """Every hole recorded since the Frame was created (or since :meth:`clear_holes`)."""
        return pd.DataFrame(self._holes) if self._holes else pd.DataFrame(columns=["status"])

    def clear_holes(self) -> None:
        self._holes.clear()

    # ------------------------------------------------------------------ references
    def ref_len(self, family: str) -> int:
        return len(self.refs[(self.anchor, family)])

    def ref_aa(self, family: str, p: int) -> str:
        return self.refs[(self.anchor, family)][int(p) - 1]

    def _check_p(self, family: str, p: Number) -> int:
        p = int(p)
        n = self.ref_len(family)
        if not 1 <= p <= n:
            raise ValueError(f"p={p} is outside the {family} reference (1..{n}); resolve sites through "
                             "Frame.sites() first")
        return p

    # ------------------------------------------------------------------ f_i : p -> q_i
    def _seq_entry(self, family: str, seq_id: str):
        if seq_id not in self.sequences.index:
            raise KeyError(f"unknown seq_id {seq_id!r}; not in {self.coords_dir / 'sequences.csv'}")
        r = self.sequences.loc[seq_id]
        if r["family"] != family:
            raise ValueError(f"{seq_id} is a {r['family']} sequence, not {family}")
        return r["group"], str(r["seq_hash"])[:12]

    def _lookup(self, family: str, seq_id: str, p: int) -> Tuple[Optional[int], str]:
        """``(0-based string index, Status.OK)`` or ``(None, why)``. An aligned ``X`` is a hole of kind
        ``unknown_residue``: a per-residue table's row for it describes an unknown token, not residue p."""
        group, h = self._seq_entry(family, seq_id)
        own = self._anchor2own[(group, family)][p - 1] if group != self.anchor else p - 1
        if own < 0:
            return None, Status.GAP_IN_SEQUENCE
        q = self._seq_maps[(group, family, h)][own]
        if q < 0:
            return None, Status.GAP_IN_SEQUENCE
        unk = self._seq_unknown.get((group, family, h))
        if unk is not None and unk.size and unk[own]:
            return None, Status.UNKNOWN_RESIDUE
        return int(q), Status.OK

    def _q_index(self, family: str, seq_id: str, p: int) -> Optional[int]:
        """0-based string index for reference position ``p``, or None (hole of any kind)."""
        return self._lookup(family, seq_id, p)[0]

    def seq_pos(self, family: str, seq_id: str, p: Number) -> Optional[int]:
        """1-based position of ``p`` in sequence ``seq_id``; None (and a recorded hole) if undefined."""
        p = self._check_p(family, p)
        q, why = self._lookup(family, seq_id, p)
        if q is None:
            self.record_hole(why, family=family, seq_id=seq_id, ref_pos=p)
            return None
        return q + 1

    def inverse_map(self, family: str, seq_id: str, include_unknown: bool = False) -> np.ndarray:
        """``q_i -> p`` for a whole sequence: an int array over 0-based string indices holding the
        1-based anchor position, ``0`` where the residue has no anchor position (an insertion) or, unless
        ``include_unknown``, where the residue is ``X``. No per-position hole is recorded."""
        group, h = self._seq_entry(family, seq_id)
        n_anchor = self.ref_len(family)
        own = np.arange(n_anchor) if group == self.anchor else self._anchor2own[(group, family)]
        smap = self._seq_maps[(group, family, h)]
        q = np.where(own >= 0, smap[np.clip(own, 0, len(smap) - 1)], -1)
        unk = self._seq_unknown.get((group, family, h))
        if not include_unknown and unk is not None and unk.size:
            q = np.where((own >= 0) & unk[np.clip(own, 0, len(unk) - 1)], -1, q)
        if "length" in self.sequences.columns and pd.notna(self.sequences.at[seq_id, "length"]):
            n_seq = int(self.sequences.at[seq_id, "length"])
        else:
            text = self.sequence(seq_id)
            n_seq = len(text) if text is not None else int(q.max()) + 1
        out = np.zeros(n_seq, dtype=np.int64)
        ok = q >= 0
        out[q[ok]] = np.flatnonzero(ok) + 1
        return out

    def to_anchor(self, group: str, family: str, own_pos: Number) -> Optional[int]:
        """Anchor position of ``own_pos`` on ``group``'s reference, or None (a hole recorded as
        ``no_ref_position``: an insertion relative to the anchor, or off the end)."""
        arr = self._ref2anchor[(group, family)]
        own_pos = int(own_pos)
        if not 1 <= own_pos <= len(arr) or arr[own_pos - 1] < 0:
            self.record_hole(Status.NO_REF_POSITION, group=group, family=family, own_pos=own_pos)
            return None
        return int(arr[own_pos - 1]) + 1

    def own_to_anchor_array(self, group: str, family: str) -> np.ndarray:
        """Whole-reference map own position -> anchor position (1-based, 0 = none), as a copy."""
        arr = self._ref2anchor[(group, family)]
        return np.where(arr >= 0, arr + 1, 0)

    def rows(self, family: str, seq_id: str, ps: Iterable[Number]) -> pd.arrays.IntegerArray:
        """0-based row index for each reference position; ``pd.NA`` where the sequence has a hole."""
        out = []
        for p in ps:
            q = self.seq_pos(family, seq_id, p)
            out.append(pd.NA if q is None else q - 1)
        return pd.array(out, dtype="Int64")

    def row_table(self, family: str, seq_id: str, ps: Iterable[Number]) -> pd.DataFrame:
        """Long table ``ref_pos, seq_pos, row, ref_aa, seq_aa, status`` for one sequence."""
        self._seq_entry(family, seq_id)
        seq = self.sequence(seq_id)
        rec = []
        for p in ps:
            p = self._check_p(family, p)
            q, why = self._lookup(family, seq_id, p)
            if q is None:
                self.record_hole(why, family=family, seq_id=seq_id, ref_pos=p)
            rec.append({"seq_id": seq_id, "ref_pos": p, "seq_pos": pd.NA if q is None else q + 1,
                        "row": pd.NA if q is None else q, "ref_aa": self.ref_aa(family, p),
                        "seq_aa": "" if q is None or seq is None else seq[q], "status": why})
        df = pd.DataFrame(rec, columns=["seq_id", "ref_pos", "seq_pos", "row", "ref_aa", "seq_aa", "status"])
        df["seq_pos"] = df["seq_pos"].astype("Int64")
        df["row"] = df["row"].astype("Int64")
        return df

    def sequence(self, seq_id: str) -> Optional[str]:
        """Sequence text for ``seq_id``, verified against the stored hash; None if it cannot be read
        (the FASTA moved, or the text no longer matches)."""
        if seq_id in self._seq_cache:
            return self._seq_cache[seq_id]
        r = self.sequences.loc[seq_id]
        text = None
        if "sequence" in r.index and isinstance(r["sequence"], str) and r["sequence"]:
            text = r["sequence"]
        elif "fasta" in r.index and isinstance(r["fasta"], str) and r["fasta"]:
            fa = (self.coords_dir / r["fasta"]).resolve()
            if fa.is_file():
                recs = read_fasta(fa)
                rec = r["record"] if "record" in r.index and isinstance(r["record"], str) else seq_id
                text = recs.get(rec) or recs.get(seq_id) or recs.get(r["family"])
        if text is not None and sha256_seq(text) != r["seq_hash"]:
            text = None
        self._seq_cache[seq_id] = text
        return text

    _sequence = sequence

    def check_hash(self, seq_id: str, sequence: Optional[str] = None, seq_hash: Optional[str] = None) -> None:
        """Raise if the sequence (or a stored ``seq_hash``) is not the string the map was built from.

        A per-residue table's rows are only sound for the string it was computed on.
        """
        want = self.sequences.loc[seq_id, "seq_hash"]
        got = seq_hash if seq_hash is not None else (sha256_seq(sequence) if sequence is not None else None)
        if got is None:
            raise ValueError("pass sequence= or seq_hash=")
        if got != want:
            raise ValueError(f"{seq_id}: sequence hash {got[:12]} != map hash {want[:12]}; the coordinate maps "
                             "describe a different string than the table being read. Rebuild with "
                             "`python -m coords build <config>`.")

    # ------------------------------------------------------------------ sites
    def sites(self, family: Optional[str] = None) -> pd.DataFrame:
        """The site table in the anchor frame, with ``scheme`` and ``status``. Unresolved rows stay in
        the table (status ``unresolved_scheme``) and are recorded as holes; callers filter on status."""
        df = self.sites_table if family is None else self.sites_table[self.sites_table["family"] == family]
        for r in df[df["status"] != Status.OK].itertuples():
            self.record_hole(r.status, family=r.family, position=r.position, scheme=getattr(r, "scheme", ""))
        return df.copy()

    # ------------------------------------------------------------------ g : p <-> ResKey
    def _row_key(self, i: int) -> ResKey:
        r = self.smaps.iloc[i]
        return ResKey(int(r["model"]), str(r["chain"]), int(r["auth_seq_id"]), str(r["ins_code"]))

    def struct(self, family: str, structure: str, p: Number, chain: Optional[str] = None) -> List[ResKey]:
        """Observed residues at reference position ``p``, one per copy of ``family`` in ``structure``
        (or only on ``chain``; holes are then recorded for that chain only)."""
        p = self._check_p(family, p)
        idx = self._by_pos.get((structure, family, p), [])
        if chain is not None:
            idx = [i for i in idx if str(self.smaps.at[i, "chain"]) == str(chain)]
        if not idx:
            self.record_hole(Status.NOT_ON_STRUCTURE, family=family, structure=structure, ref_pos=p,
                             **({"chain": chain} if chain is not None else {}))
            return []
        keys = []
        for i in sorted(idx, key=lambda i: self._chain_order[(structure, str(self.smaps.at[i, "chain"]))]):
            r = self.smaps.iloc[i]
            if not r["observed"] or pd.isna(r["auth_seq_id"]):
                self.record_hole(Status.NOT_OBSERVED, family=family, structure=structure, ref_pos=p,
                                 chain=r["chain"])
                continue
            keys.append(self._row_key(i))
        return keys

    def ref_pos(self, structure: str, reskey: ResKey) -> Optional[int]:
        """Reference position of a structure residue, or None (hole recorded)."""
        i = self._by_reskey.get((structure, reskey.chain, int(reskey.auth_seq_id), reskey.ins_code))
        if i is None:
            self.record_hole(Status.UNREGISTERED_CHAIN, structure=structure, reskey=str(reskey))
            return None
        r = self.smaps.iloc[i]
        if pd.isna(r["ref_pos"]):
            self.record_hole(r["status"], structure=structure, reskey=str(reskey))
            return None
        return int(r["ref_pos"])

    def chain_family(self, structure: str, chain: str) -> str:
        return self._chain_family.get((structure, str(chain)), "")

    # ------------------------------------------------------------------ geometry
    def structure_path(self, structure: str) -> Path:
        hit = self.structures[self.structures["structure"] == structure]
        if hit.empty:
            raise FileNotFoundError(f"unknown structure {structure!r}; known: {sorted(self.structures['structure'])}")
        path = (self.coords_dir / hit.iloc[0]["path"]).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"{path} (recorded in structures.csv) does not exist; the maps were built "
                                    "for a file that has moved")
        return path

    def atoms(self, structure: str) -> pd.DataFrame:
        """Protein atoms of ``structure`` (first model, one altloc), read lazily from its CIF."""
        if structure not in self._atoms:
            from .cif import read_structure
            self._atoms[structure] = read_structure(self.structure_path(structure), structure).atoms
        return self._atoms[structure]

    def neighbors(self, family: str, structure: str, p: Number, radius: float, center: str = "CA",
                  copy: Union[str, int] = "first", seq_ids: Optional[Sequence[str]] = None,
                  cross_family: bool = False, include_holes: bool = False) -> pd.DataFrame:
        """Residues within ``radius`` Å of reference position ``p``, in long format.

        center   ``"CA"``: CA–CA distance. ``"all"``: smallest heavy-atom distance between the two
                 residues (a CA-only centre misses side-chain contacts).
        copy     ``"first"`` (default), ``"all"``, or a chain id: which copies of the query residue to use.
        seq_ids  if given, each neighbor row is crossed with every sequence of the neighbor's family in this
                 list, adding ``seq_id, seq_pos, row, seq_aa`` (the row to read).
        cross_family  neighbors on another family's chain (a heteromer) carry ``status=other_family`` and
                 are dropped by default, because their ``neighbor_ref_pos`` is in another family's frame.
                 With ``cross_family=True`` they are kept as ``ok`` and crossed with the sequences of
                 *their own* family.
        Rows with ``status != ok`` are recorded in ``Frame.holes`` and dropped unless ``include_holes``.
        """
        if center not in ("CA", "all"):
            raise ValueError("center must be 'CA' or 'all'")
        p = self._check_p(family, p)
        keys = self.struct(family, structure, p)
        if copy == "first":
            keys = keys[:1]
        elif copy != "all":
            keys = [k for k in keys if k.chain == str(copy)]
        atoms = self.atoms(structure)
        if center == "CA":
            atoms = atoms[atoms["atom"] == "CA"]
        coords = atoms[["x", "y", "z"]].to_numpy()
        a_chain = atoms["chain"].to_numpy()
        a_auth = atoms["auth_seq_id"].to_numpy()
        a_ins = atoms["ins_code"].to_numpy()
        from scipy.spatial import cKDTree
        tree = cKDTree(coords) if len(coords) else None
        rows: List[dict] = []
        for qk in keys:
            sel = (a_chain == qk.chain) & (a_auth == qk.auth_seq_id) & (a_ins == qk.ins_code)
            if not sel.any() or tree is None:
                self.record_hole(Status.NOT_OBSERVED, family=family, structure=structure, ref_pos=p, chain=qk.chain)
                continue
            best: Dict[tuple, float] = {}
            for xyz in coords[sel]:
                for j in tree.query_ball_point(xyz, radius):
                    d = float(np.linalg.norm(coords[j] - xyz))
                    k = (a_chain[j], int(a_auth[j]), a_ins[j])
                    if d < best.get(k, 1e9):
                        best[k] = d
            for (ch, auth, ins), d in sorted(best.items(), key=lambda kv: (kv[1], kv[0])):
                nk = ResKey(qk.model, ch, auth, ins)
                i = self._by_reskey.get((structure, ch, auth, ins))
                if i is None:
                    nfam, nref, status, naa = self.chain_family(structure, ch), None, Status.UNREGISTERED_CHAIN, ""
                else:
                    r = self.smaps.iloc[i]
                    nfam, naa = r["family"], r["struct_aa"]
                    nref = None if pd.isna(r["ref_pos"]) else int(r["ref_pos"])
                    if not nfam:
                        status = Status.UNREGISTERED_CHAIN
                    elif nref is None:
                        status = Status.NO_REF_POSITION
                    elif nfam != family and not cross_family:
                        status = Status.OTHER_FAMILY
                    else:
                        status = Status.OK
                rows.append({"structure": structure, "family": family, "query_p": p, "query_reskey": qk,
                             "neighbor_reskey": nk, "neighbor_family": nfam, "neighbor_ref_pos": nref,
                             "neighbor_aa": naa, "ref_aa": self.ref_aa(nfam, nref) if nref and nfam else "",
                             "distance": round(d, 3), "is_query": nk == qk, "status": status})
        df = pd.DataFrame(rows, columns=["structure", "family", "query_p", "query_reskey", "neighbor_reskey",
                                         "neighbor_family", "neighbor_ref_pos", "neighbor_aa", "ref_aa",
                                         "distance", "is_query", "status"])
        df["neighbor_ref_pos"] = df["neighbor_ref_pos"].astype("Int64")
        for r in df[df["status"] != Status.OK].itertuples():
            self.record_hole(r.status, family=family, structure=structure, ref_pos=p,
                             neighbor=str(r.neighbor_reskey), neighbor_family=r.neighbor_family)
        if seq_ids is not None:
            df = self._cross(df, seq_ids)
        if not include_holes:
            df = df[df["status"] == Status.OK].reset_index(drop=True)
        return df

    def _cross(self, nb: pd.DataFrame, seq_ids: Sequence[str]) -> pd.DataFrame:
        by_family: Dict[str, List[str]] = defaultdict(list)
        for sid in seq_ids:
            if sid not in self.sequences.index:
                raise KeyError(f"unknown seq_id {sid!r}")
            by_family[self.sequences.at[sid, "family"]].append(sid)
        out: List[dict] = []
        for r in nb.to_dict("records"):
            if r["status"] != Status.OK:
                out.append({**r, "seq_id": None, "seq_pos": pd.NA, "row": pd.NA, "seq_aa": ""})
                continue
            for sid in by_family.get(r["neighbor_family"], []):
                q, why = self._lookup(r["neighbor_family"], sid, int(r["neighbor_ref_pos"]))
                seq = self.sequence(sid)
                if q is None:
                    self.record_hole(why, family=r["neighbor_family"], seq_id=sid, ref_pos=int(r["neighbor_ref_pos"]))
                    out.append({**r, "seq_id": sid, "seq_pos": pd.NA, "row": pd.NA, "seq_aa": "", "status": why})
                else:
                    out.append({**r, "seq_id": sid, "seq_pos": q + 1, "row": q,
                                "seq_aa": seq[q] if seq else "", "status": Status.OK})
        df = pd.DataFrame(out, columns=list(nb.columns) + ["seq_id", "seq_pos", "row", "seq_aa"])
        df["seq_pos"] = df["seq_pos"].astype("Int64")
        df["row"] = df["row"].astype("Int64")
        return df
