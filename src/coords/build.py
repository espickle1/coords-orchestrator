"""Build the coordinate maps from a ``coords.toml`` and write them to its ``out_dir``.

    python -m coords build coords.toml

Artifacts (CSV, human-readable, consumed by merge):

    reference/<group>/<family>.fasta   the designated reference S_0, one record each
    reference_maps.csv                 every non-anchor group's reference -> the anchor's, per residue
    sequences.csv                      one row per input sequence: hash, length, identity, quality
    sequence_maps.csv                  f_i per distinct sequence (hash12 = first 12 hex of the sha256 in
                                       sequences.csv): ref_pos -> seq_pos on the sequence's own group
    structures.csv                     one row per structure file: key (file stem), entry id, path, sha256
    structure_chains.csv               one row per protein chain: family, group, identities, quality
    structure_maps.csv                 g per SEQRES residue: author number <-> anchor ref_pos
    sites_resolved.csv                 the site table in the anchor frame, with scheme and status
    manifest.json                      input hashes, aligner + version, counts, per-chain identity

``ref_pos`` in ``sequence_maps.csv`` counts residues of the sequence's own-group reference;
``structure_maps.csv`` is already in the anchor frame. Reference-to-reference composition happens at
read time, in :class:`coords.Frame`.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from . import __version__
from . import align as al
from . import cif
from .config import UNRESOLVED, Config, ConfigError, Thresholds, load_config, load_references
from .seqio import normalize, read_fasta, sha256_seq, write_fasta
from .types import Status

SCHEMA_VERSION = 1
SEQ_REQUIRED = ("seq_id", "group", "family")
SEQ_COMPUTED = ("fasta", "record", "sequence", "seq_hash", "length", "ref_length", "n_aligned", "n_unknown",
                "identity", "coverage_of_ref", "n_inserted", "n_gap", "quality")


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _rel(path: Path, start: Path) -> str:
    """Posix path of ``path`` relative to ``start``; absolute if they share no root (other drive)."""
    try:
        return Path(os.path.relpath(path, start)).as_posix()
    except ValueError:
        return Path(path).as_posix()


# -----------------------------------------------------------------------------------------
# references and reference-to-reference maps
# -----------------------------------------------------------------------------------------
def build_reference_maps(refs: Dict[Tuple[str, str], str], cfg: Config):
    """Return ``(DataFrame, to_anchor, stats)``; ``to_anchor[(group, family)]`` maps own idx -> anchor idx."""
    rows: List[pd.DataFrame] = []
    to_anchor: Dict[Tuple[str, str], np.ndarray] = {}
    stats: Dict[str, dict] = {}
    anchor = cfg.anchor
    for (group, family), own in refs.items():
        if group == anchor:
            to_anchor[(group, family)] = np.arange(len(own))
            continue
        anc = refs[(anchor, family)]
        pa = al.align_pair(anc, own, cfg.aligner)               # ref = anchor, query = own
        to_anchor[(group, family)] = pa.q_to_ref
        ok = pa.q_to_ref >= 0
        df = pd.DataFrame({
            "family": family, "from_group": group, "to_group": anchor,
            "from_pos": np.arange(1, len(own) + 1),
            "to_pos": pd.array(np.where(ok, pa.q_to_ref + 1, 0), dtype="Int64"),
            "from_aa": list(own),
            "to_aa": [anc[j] if j >= 0 else "" for j in pa.q_to_ref],
        })
        df.loc[~ok, "to_pos"] = pd.NA
        df["identical"] = (df["from_aa"] == df["to_aa"]) & ok
        df["status"] = np.where(ok, Status.OK, Status.NO_REF_POSITION)
        rows.append(df)
        stats[f"{group}->{anchor}.{family}"] = {
            "identity": round(pa.identity, 4), "n_aligned": pa.n_aligned,
            "from_len": len(own), "anchor_len": len(anc)}
    maps = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["family", "from_group", "to_group", "from_pos", "to_pos", "from_aa", "to_aa", "identical", "status"])
    return maps, to_anchor, stats


# -----------------------------------------------------------------------------------------
# sequences
# -----------------------------------------------------------------------------------------
def load_sequences(cfg: Config, refs) -> pd.DataFrame:
    """The input sequence table with the sequence text resolved and its hash verified.

    Required columns: ``seq_id, group, family`` and ``sequence`` or ``fasta`` (path relative to the
    config file; ``record`` names the FASTA record, default ``seq_id`` then ``family``). An optional
    ``seq_hash`` is verified. Any other column is carried through to ``sequences.csv``.
    """
    spec = cfg.sequences
    empty = pd.DataFrame(columns=list(SEQ_REQUIRED) + ["sequence"])
    if spec is None:
        return empty
    if not spec.table.is_file():
        raise ConfigError(f"[sequences] table not found: {spec.table}")
    t = pd.read_csv(spec.table, keep_default_na=False, na_values=[""], dtype=str)
    missing = [c for c in SEQ_REQUIRED if c not in t.columns]
    if missing:
        raise ConfigError(f"{spec.table}: missing column(s) {missing}")
    if "sequence" not in t.columns and "fasta" not in t.columns:
        raise ConfigError(f"{spec.table}: needs a `sequence` column or a `fasta` column")
    clash = [c for c in t.columns if c in SEQ_COMPUTED and c not in ("sequence", "fasta", "record", "seq_hash")]
    if clash:
        raise ConfigError(f"{spec.table}: column(s) {clash} are written by the build; rename them")
    if t["seq_id"].duplicated().any():
        raise ConfigError(f"{spec.table}: duplicate seq_id {sorted(t['seq_id'][t['seq_id'].duplicated()])[:5]}")
    for c in ("sequence", "fasta", "record", "seq_hash"):
        if c not in t.columns:
            t[c] = pd.NA
    cache: Dict[Path, Dict[str, str]] = {}
    seqs, fastas, records, bad = [], [], [], []
    for r in t.itertuples():
        if (r.group, r.family) not in refs:
            bad.append(f"{r.seq_id}: no reference for group {r.group!r} family {r.family!r}")
            seqs.append(""), fastas.append(""), records.append("")
            continue
        if isinstance(r.sequence, str) and r.sequence:
            seq, fa, rec = normalize(r.sequence), "", ""
        elif isinstance(r.fasta, str) and r.fasta:
            path = (cfg.base / r.fasta).resolve()
            if path not in cache:
                if not path.is_file():
                    raise ConfigError(f"{r.seq_id}: fasta not found: {path}")
                cache[path] = read_fasta(path)
            recs = cache[path]
            rec = next((k for k in (r.record, r.seq_id, r.family) if isinstance(k, str) and k in recs), None)
            if rec is None:
                bad.append(f"{r.seq_id}: no record {r.record if isinstance(r.record, str) else r.seq_id!r} "
                           f"(or {r.family!r}) in {path.name}")
                seqs.append(""), fastas.append(""), records.append("")
                continue
            seq, fa = recs[rec], _rel(path, cfg.out_dir)
        else:
            bad.append(f"{r.seq_id}: neither `sequence` nor `fasta` given")
            seqs.append(""), fastas.append(""), records.append("")
            continue
        if not seq:
            bad.append(f"{r.seq_id}: empty sequence")
        if isinstance(r.seq_hash, str) and r.seq_hash and sha256_seq(seq) != r.seq_hash:
            bad.append(f"{r.seq_id}: sha256 of the sequence differs from the table's seq_hash")
        seqs.append(seq), fastas.append(fa), records.append(rec)
    if bad:
        raise ConfigError("sequence table problems:\n  " + "\n  ".join(bad[:20])
                          + (f"\n  ... and {len(bad) - 20} more" if len(bad) > 20 else ""))
    t["sequence"], t["fasta"], t["record"] = seqs, fastas, records
    t["seq_hash"] = [sha256_seq(s) for s in seqs]
    return t


def build_sequence_maps(refs, corpus: pd.DataFrame, cfg: Config):
    """Align every distinct sequence to its own-group S_0. Returns (sequences, maps)."""
    thr: Thresholds = cfg.thresholds
    maps: List[pd.DataFrame] = []
    info: List[dict] = []
    done: Dict[Tuple[str, str, str], dict] = {}
    seen12: Dict[Tuple[str, str, str], str] = {}
    passthrough = [c for c in corpus.columns if c not in list(SEQ_REQUIRED) + ["sequence", "fasta", "record", "seq_hash"]]
    for r in corpus.to_dict("records"):
        group, family, seq, h = r["group"], r["family"], r["sequence"], r["seq_hash"]
        ref = refs[(group, family)]
        key = (group, family, h)
        h12 = (group, family, h[:12])
        if seen12.setdefault(h12, h) != h:
            raise ConfigError(f"hash12 collision between two distinct sequences ({h[:12]}); report this")
        if key not in done:
            pa = al.align_pair(ref, seq, cfg.aligner)
            ok = pa.ref_to_q >= 0
            df = pd.DataFrame({
                "group": group, "family": family, "hash12": h[:12],
                "ref_pos": np.arange(1, len(ref) + 1),
                "seq_pos": pd.array(np.where(ok, pa.ref_to_q + 1, 0), dtype="Int64"),
                "ref_aa": list(ref),
                "seq_aa": [seq[j] if j >= 0 else "" for j in pa.ref_to_q],
            })
            df.loc[~ok, "seq_pos"] = pd.NA
            df["identical"] = (df["ref_aa"] == df["seq_aa"]) & ok
            unk = ok & (df["seq_aa"] == "X").to_numpy()
            df["status"] = np.where(~ok, Status.GAP_IN_SEQUENCE, np.where(unk, Status.UNKNOWN_RESIDUE, Status.OK))
            maps.append(df)
            n_unknown = int(unk.sum())
            cov = pa.coverage_of_ref(len(ref))
            if not (pa.identity >= thr.min_seq_identity and pa.n_known >= thr.min_seq_aligned):
                q = "unreliable"
            elif cov - n_unknown / len(ref) < thr.min_seq_coverage:
                q = "partial"                  # trustworthy where defined, but most positions are holes
            else:
                q = "ok"
            done[key] = dict(
                length=len(seq), ref_length=len(ref), n_aligned=pa.n_aligned, n_unknown=n_unknown,
                identity=round(float(pa.identity), 4) if pa.n_known else float("nan"),
                coverage_of_ref=round(cov, 4), n_inserted=int((pa.q_to_ref < 0).sum()),
                n_gap=int((pa.ref_to_q < 0).sum()), quality=q)
        row = {"seq_id": r["seq_id"], "group": group, "family": family, **{c: r[c] for c in passthrough},
               "fasta": r["fasta"], "record": r["record"], "seq_hash": h, **done[key]}
        if not r["fasta"]:
            row["sequence"] = seq                   # inline input: keep the text so a Frame can read it
        info.append(row)
    cols = (["seq_id", "group", "family"] + passthrough + ["fasta", "record", "seq_hash", "length", "ref_length",
            "n_aligned", "n_unknown", "identity", "coverage_of_ref", "n_inserted", "n_gap", "quality"])
    seq_df = pd.DataFrame(info)
    if "sequence" in seq_df.columns:
        cols.append("sequence")
    seq_df = seq_df.reindex(columns=cols) if len(seq_df) else pd.DataFrame(columns=cols)
    map_df = pd.concat(maps, ignore_index=True) if maps else pd.DataFrame(
        columns=["group", "family", "hash12", "ref_pos", "seq_pos", "ref_aa", "seq_aa", "identical", "status"])
    return seq_df, map_df


# -----------------------------------------------------------------------------------------
# structures
# -----------------------------------------------------------------------------------------
def _assign_chain(seqres: str, refs, cfg: Config):
    """Best designated reference for a chain, by number of identical aligned residues.

    Candidates are tried anchor group first, so a tie goes to the anchor. A chain is registered only if
    its best alignment clears three gates: identity, coverage of the chain, and a z-score of the
    alignment score against shuffled copies of the chain (free end gaps let an unrelated short chain
    reach a deceptive identity; its z-score stays near 0, a real homologue's is large). Returns
    ``((group, family, PairAlignment, z) | None, scores)``; ``scores`` lets a reviewer see what the
    chain was weighed against.
    """
    thr = cfg.thresholds
    best, scores = None, {}
    order = [(g, f) for g in cfg.groups for (gg, f) in refs if gg == g]
    for (group, family) in order:
        pa = al.align_pair(refs[(group, family)], seqres, cfg.aligner)
        scores[f"{group}.{family}"] = pa.n_identical
        if best is None or pa.n_identical > best[2].n_identical:
            best = (group, family, pa)
    if best is None or best[2].identity < thr.min_chain_identity or \
            best[2].coverage_of_query(len(seqres)) < thr.min_chain_coverage:
        return None, scores
    z = _shuffle_z(refs[(best[0], best[1])], seqres, best[2].score, thr.zscore_shuffles, cfg)
    if not z >= thr.min_chain_zscore:
        return None, scores
    return (*best, z), scores


def _shuffle_z(ref: str, query: str, score: float, n: int, cfg: Config) -> float:
    """z-score of an alignment score against ``n`` alignments of shuffled copies of ``query``."""
    import random
    rng = random.Random(sha256_seq(query))              # deterministic: the build is reproducible
    vals = []
    for _ in range(max(n, 2)):
        s = list(query)
        rng.shuffle(s)
        vals.append(al.align_pair(ref, "".join(s), cfg.aligner).score)
    sd = float(np.std(vals))
    return (score - float(np.mean(vals))) / sd if sd > 0 else float("inf")


def find_structures(cfg: Config) -> List[Path]:
    if cfg.structures is None:
        return []
    if not cfg.structures.dir.is_dir():
        raise ConfigError(f"[structures] dir not found: {cfg.structures.dir}")
    files = sorted(p for p in cfg.structures.dir.glob(cfg.structures.glob) if p.is_file())
    stems: Dict[str, Path] = {}
    for p in files:
        if p.stem in stems:
            raise ConfigError(f"two structure files share the key {p.stem!r}: {stems[p.stem]} and {p}. "
                              "The file stem is the structure key; rename one.")
        stems[p.stem] = p
    return files


def build_structure_maps(refs, to_anchor, cfg: Config, files: List[Path]):
    anchor = cfg.anchor
    thr = cfg.thresholds
    chain_rows: List[dict] = []
    res_rows: List[pd.DataFrame] = []
    struct_rows: List[dict] = []
    notes: Dict[str, List[str]] = {}
    assigned: Dict[str, tuple] = {}        # seqres string -> assignment (homo-oligomers share it)
    for path in files:
        stem = path.stem
        st = cif.read_structure(path, structure=stem)
        notes[stem] = st.notes
        struct_rows.append(dict(structure=stem, entry_id=st.entry_id, path=_rel(path, cfg.out_dir),
                                sha256=_sha256_file(path)))
        model = int(st.atoms["model"].min()) if len(st.atoms) else 1
        ent = st.entities.set_index("entity_id") if len(st.entities) else pd.DataFrame()
        for chain, (seqres, g) in cif.seqres_by_chain(st).items():
            if not g["poly_type"].iloc[0].startswith("polypeptide"):
                continue
            entity = str(g["entity_id"].iloc[0])
            desc = str(ent.loc[entity, "description"]) if len(ent) and entity in ent.index else ""
            if seqres not in assigned:
                assigned[seqres] = _assign_chain(seqres, refs, cfg)
            res, _scores = assigned[seqres]
            n = len(seqres)
            if res is None:
                family = group = ""
                q_to_own = np.full(n, -1)
                to_anc = np.full(n, -1)
                hopa = float("nan")
                coverage = 0.0
                own_ref = None
                zscore = float("nan")
            else:
                group, family, pa, zscore = res
                own_ref = refs[(group, family)]
                q_to_own = pa.q_to_ref
                ta = to_anchor[(group, family)]
                to_anc = np.where(q_to_own >= 0, ta[np.clip(q_to_own, 0, len(ta) - 1)], -1)
                hopa = pa.identity
                coverage = pa.coverage_of_query(n)
            anchor_ref = refs[(anchor, family)] if family else ""
            ok = to_anc >= 0
            aa = np.array(list(seqres))
            ref_aa = np.array([anchor_ref[j] if j >= 0 else "" for j in to_anc])
            own_aa = np.array([own_ref[j] if (own_ref and j >= 0) else "" for j in q_to_own])
            n_ok = int(ok.sum())
            ident_anchor = float((aa[ok] == ref_aa[ok]).sum()) / n_ok if n_ok else float("nan")
            observed = g["observed"].to_numpy()
            status = np.where(~observed, Status.NOT_OBSERVED, np.where(ok, Status.OK, Status.NO_REF_POSITION))
            if family == "":
                status = np.full(n, Status.UNREGISTERED_CHAIN)
            quality = ("unregistered" if family == "" else
                       "ok" if ident_anchor >= thr.low_chain_identity else "low_identity")
            # direct alignment to the anchor, for the cross-group agreement check
            direct_ident = direct_agree = float("nan")
            if family and group != anchor:
                pd_ = al.align_pair(anchor_ref, seqres, cfg.aligner)
                direct_ident = pd_.identity
                both = ok & (pd_.q_to_ref >= 0)
                direct_agree = float((to_anc[both] == pd_.q_to_ref[both]).mean()) if both.any() else float("nan")
            chain_rows.append(dict(
                structure=stem, entry_id=st.entry_id, model=model, chain=chain,
                asym_id=str(g["asym_id"].iloc[0]), entity=entity, description=desc,
                family=family, group=group, seqres_len=n, n_observed=int(observed.sum()),
                n_mapped=n_ok, hop_a_identity=round(hopa, 4) if hopa == hopa else np.nan,
                anchor_identity=round(ident_anchor, 4) if ident_anchor == ident_anchor else np.nan,
                direct_identity=round(direct_ident, 4) if direct_ident == direct_ident else np.nan,
                direct_agreement=round(direct_agree, 4) if direct_agree == direct_agree else np.nan,
                coverage=round(coverage, 4), score_z=round(zscore, 2) if zscore == zscore else np.nan,
                quality=quality, n_ins_codes=int((g["ins_code"] != "").sum()), scheme_variant=st.scheme_variant))
            df = pd.DataFrame({
                "structure": stem, "model": model, "chain": chain, "entity": entity, "family": family,
                "seqres_idx": g["seq_id"].to_numpy(),
                "auth_seq_id": g["auth_seq_id"].astype("Int64").to_numpy(),
                "ins_code": g["ins_code"].to_numpy(),
                "ref_pos": pd.array(np.where(ok, to_anc + 1, 0), dtype="Int64"),
                "struct_aa": aa, "ref_aa": ref_aa,
                "identical": (aa == ref_aa) & ok,
                "observed": observed, "status": status,
                "group_ref_pos": pd.array(np.where(q_to_own >= 0, q_to_own + 1, 0), dtype="Int64"),
                "group_ref_aa": own_aa,
            })
            df.loc[~ok, "ref_pos"] = pd.NA
            df.loc[df["group_ref_pos"] == 0, "group_ref_pos"] = pd.NA
            res_rows.append(df)
    chains = pd.DataFrame(chain_rows, columns=_CHAIN_COLS)
    smaps = pd.concat(res_rows, ignore_index=True) if res_rows else pd.DataFrame(columns=_SMAP_COLS)
    structs = pd.DataFrame(struct_rows, columns=["structure", "entry_id", "path", "sha256"])
    return chains, smaps, structs, notes


_CHAIN_COLS = ["structure", "entry_id", "model", "chain", "asym_id", "entity", "description", "family", "group",
               "seqres_len", "n_observed", "n_mapped", "hop_a_identity", "anchor_identity", "direct_identity",
               "direct_agreement", "coverage", "score_z", "quality", "n_ins_codes", "scheme_variant"]
_SMAP_COLS = ["structure", "model", "chain", "entity", "family", "seqres_idx", "auth_seq_id", "ins_code",
              "ref_pos", "struct_aa", "ref_aa", "identical", "observed", "status", "group_ref_pos", "group_ref_aa"]


# -----------------------------------------------------------------------------------------
# sites
# -----------------------------------------------------------------------------------------
SITE_COLS = ["family", "position", "residue", "scheme", "ref_pos", "ref_aa", "status", "matches_ref"]


def resolve_sites(cfg: Config, refs) -> pd.DataFrame:
    """Put each site on a reference position of the anchor, naming the scheme it came from.

    ``ref_pos = position + offset`` of the row's scheme. A row whose scheme is missing, ``unresolved``
    or not listed in ``[sites.schemes]`` is excluded (``status=unresolved_scheme``), never guessed. A
    resolved position outside the anchor reference, or a family without a reference, is
    ``no_ref_position``. All other input columns are kept.
    """
    spec = cfg.sites
    if spec is None:
        return pd.DataFrame(columns=SITE_COLS)
    if not spec.table.is_file():
        raise ConfigError(f"[sites] table not found: {spec.table}")
    t = pd.read_csv(spec.table, keep_default_na=False, na_values=[""], dtype=str)
    need = {"family": spec.family_column, "position": spec.position_column, "scheme": spec.scheme_column}
    missing = [c for c in need.values() if c not in t.columns]
    if missing:
        raise ConfigError(f"{spec.table}: missing column(s) {missing}")
    rename = {spec.family_column: "family", spec.position_column: "position", spec.scheme_column: "scheme"}
    if spec.residue_column in t.columns:
        rename[spec.residue_column] = "residue"
    t = t.rename(columns=rename)
    if "residue" not in t.columns:
        t["residue"] = ""
    out = []
    for r in t.to_dict("records"):
        family, scheme = r["family"], r["scheme"]
        try:
            pos = int(r["position"])
        except (TypeError, ValueError):
            raise ConfigError(f"{spec.table}: position {r['position']!r} for family {family!r} is not an integer")
        base = {**r, "position": pos, "ref_pos": pd.NA, "ref_aa": "", "matches_ref": pd.NA}
        off = spec.schemes.get(scheme) if isinstance(scheme, str) else None
        if not isinstance(scheme, str) or scheme == UNRESOLVED or scheme not in spec.schemes or off is None:
            out.append({**base, "status": Status.UNRESOLVED_SCHEME})
            continue
        ref = refs.get((cfg.anchor, family))
        p = pos + off
        if ref is None or not 1 <= p <= len(ref):
            out.append({**base, "status": Status.NO_REF_POSITION})
            continue
        resid = r.get("residue") or ""
        out.append({**base, "ref_pos": p, "ref_aa": ref[p - 1], "status": Status.OK,
                    "matches_ref": (ref[p - 1] == resid) if resid else pd.NA})
    df = pd.DataFrame(out)
    extras = [c for c in df.columns if c not in SITE_COLS]
    df = df[["family", "position", "residue", "scheme", "ref_pos", "ref_aa", "status", "matches_ref"] + extras]
    df["ref_pos"] = df["ref_pos"].astype("Int64")
    df["matches_ref"] = df["matches_ref"].astype("boolean")
    return df


# -----------------------------------------------------------------------------------------
# driver
# -----------------------------------------------------------------------------------------
def build_all(cfg: Config, verbose: bool = True) -> dict:
    out_dir = cfg.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    log = (lambda *a: print(*a, file=sys.stderr)) if verbose else (lambda *a: None)

    refs = load_references(cfg)
    for (group, family), seq in refs.items():
        d = out_dir / "reference" / group
        d.mkdir(parents=True, exist_ok=True)
        write_fasta(d / f"{family}.fasta", {f"{group}.{family} | {cfg.references[group].label} | {len(seq)} aa": seq})
    log(f"references: {len(refs)} ({len(cfg.references)} groups)")

    ref_maps, to_anchor, ref_stats = build_reference_maps(refs, cfg)
    ref_maps.to_csv(out_dir / "reference_maps.csv", index=False)

    corpus = load_sequences(cfg, refs)
    sequences, seq_maps = build_sequence_maps(refs, corpus, cfg)
    sequences.to_csv(out_dir / "sequences.csv", index=False)
    seq_maps.to_csv(out_dir / "sequence_maps.csv", index=False)
    log(f"sequences: {len(sequences)}  sequence_maps rows: {len(seq_maps)}")

    files = find_structures(cfg)
    chains, struct_maps, structs, notes = build_structure_maps(refs, to_anchor, cfg, files)
    chains.to_csv(out_dir / "structure_chains.csv", index=False)
    struct_maps.to_csv(out_dir / "structure_maps.csv", index=False)
    structs.to_csv(out_dir / "structures.csv", index=False)
    log(f"structures: {len(structs)}  chains: {len(chains)}  structure_maps rows: {len(struct_maps)}")

    sites = resolve_sites(cfg, refs)
    sites.to_csv(out_dir / "sites_resolved.csv", index=False)

    inputs = {_rel(spec.fasta, cfg.base): _sha256_file(spec.fasta) for spec in cfg.references.values()}
    if cfg.sequences is not None:
        inputs[_rel(cfg.sequences.table, cfg.base)] = _sha256_file(cfg.sequences.table)
    if cfg.sites is not None:
        inputs[_rel(cfg.sites.table, cfg.base)] = _sha256_file(cfg.sites.table)
    for p, row in zip(files, structs.to_dict("records")):
        inputs[_rel(p, cfg.base)] = row["sha256"]
    anchor_families = sorted(f for (g, f) in refs if g == cfg.anchor)
    chain_cols = ["structure", "chain", "family", "group", "hop_a_identity", "anchor_identity",
                  "direct_identity", "direct_agreement", "coverage", "score_z", "quality"]
    per_chain = chains[chain_cols].astype(object).where(chains[chain_cols].notna(), None).to_dict("records")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "coords_version": __version__,
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "project": cfg.name,
        "anchor": cfg.anchor,
        "groups": cfg.groups,
        "families": anchor_families,
        "references": {f"{g}.{f}": {"label": cfg.references[g].label, "length": len(q)} for (g, f), q in refs.items()},
        "aligner": {**cfg.aligner.as_dict(), "biopython": al.biopython_version(),
                    "note": "a different aligner can place divergent pairs differently"},
        "thresholds": {k: getattr(cfg.thresholds, k) for k in cfg.thresholds.__dataclass_fields__},
        "reference_to_reference": ref_stats,
        "counts": {
            "sequences": int(len(sequences)),
            "unique_sequence_maps": int(seq_maps[["group", "family", "hash12"]].drop_duplicates().shape[0]),
            "sequences_unreliable": int((sequences["quality"] == "unreliable").sum()) if len(sequences) else 0,
            "sequences_partial": int((sequences["quality"] == "partial").sum()) if len(sequences) else 0,
            "structures": int(len(structs)),
            "protein_chains": int(len(chains)),
            "chains_unregistered": int((chains["quality"] == "unregistered").sum()),
            "chains_low_identity": int((chains["quality"] == "low_identity").sum()),
            "sites": int(len(sites)),
            "sites_unresolved": int((sites["status"] != Status.OK).sum()) if len(sites) else 0,
        },
        "per_chain": per_chain,
        "structure_notes": {k: v for k, v in notes.items() if v},
        "inputs_sha256": inputs,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n", encoding="utf-8")
    log(f"manifest written to {out_dir}")
    return manifest


def build(config_path) -> dict:
    """Load ``coords.toml`` and build its maps. Returns the manifest."""
    return build_all(load_config(config_path))
